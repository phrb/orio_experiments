import sys, time
import math
import random
import numpy
import orio.main.tuner.search.search
from orio.main.util.globals import *
import copy
import json
import dataset
import os

import rpy2.rinterface as ri
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import DataFrame, IntVector, FloatVector, StrVector, BoolVector, Formula, r

class Doptanova(orio.main.tuner.search.search.Search):
    '''
    The search engine that uses a random search approach, enhanced with a local search that finds
    the best neighboring coordinate.

    Below is a list of algorithm-specific arguments used to steer the search algorithm.
      local_distance            the distance number used in the local search to find the best
                                neighboring coordinate located within the specified distance
    '''

    # algorithm-specific argument names
    __LOCAL_DIST = 'local_distance'  # default: 0

    def __init__(self, params):
        '''To instantiate a random search engine'''
        numpy.random.seed(39920)

        self.base      = importr("base")
        self.utils     = importr("utils")
        self.stats     = importr("stats")
        self.algdesign = importr("AlgDesign")
        self.car       = importr("car")

        self.total_runs = 20
        orio.main.tuner.search.search.Search.__init__(self, params)

        self.parameter_ranges = {}

        for i in range(len(self.params["axis_val_ranges"])):
            self.parameter_ranges[self.params["axis_names"][i]] = [0,
                    len(self.params["axis_val_ranges"][i])]

        info("Parameters: " + str(self.parameter_ranges))

        # set all algorithm-specific arguments to their default values
        self.local_distance = 0

        # read all algorithm-specific arguments
        self.__readAlgoArgs()

        self.init_samp = 2 * self.total_runs  # BN: used to be hard-coded to 10,000

        # complain if both the search time limit and the total number of search runs are undefined
        if self.time_limit <= 0 and self.total_runs <= 0:
            err((
                '%s search requires search time limit or '
                + 'total number of search runs to be defined') %
                self.__class__.__name__)

    def isclose(self, a, b, rel_tol = 1e-09, abs_tol = 0.0):
        return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

    def opt_federov(self, design_formula, trials, constraint, data):
        info("Starting \"optMonteCarlo\" run")
        info(str(data))

        self.base.set_seed(77126)

        candidate_multiplier = 20
        repetitions          = 8

        output = self.algdesign.optMonteCarlo(frml        = Formula(design_formula),
                                              data        = data,
                                              nCand       = candidate_multiplier * trials,
                                              nRepeats    = repetitions,
                                              constraints = constraint,
                                              nTrials     = trials)

        return output

    def transform_lm(self, design, lm_formula):
        info("Power Transform Step:")
        response  = lm_formula.split("~")[0].strip()
        variables = lm_formula.split("~")[1].strip()

        r_snippet = """boxcox_t <- powerTransform(%s, data = %s)
        regression <- lm(bcPower(%s, boxcox_t$lambda) ~ %s, data = %s)
        regression""" %(lm_formula, design.r_repr(), response, variables, design.r_repr())

        transformed_lm = robjects.r(r_snippet)
        return transformed_lm

    def anova(self, design, formula):
        regression = self.stats.lm(Formula(formula), data = design)
        heteroscedasticity_test = self.car.ncvTest(regression)
        info("Heteroscedasticity Test p-value: " + str(heteroscedasticity_test.rx("p")[0][0]))

        heteroscedasticity_threshold = 0.05
        if heteroscedasticity_test.rx("p")[0][0] < heteroscedasticity_threshold:
            regression = self.transform_lm(design, formula)
            heteroscedasticity_test = self.car.ncvTest(regression)
            info("Heteroscedasticity Test p-value: " + str(heteroscedasticity_test.rx("p")[0][0]))

        summary_regression = self.stats.summary_aov(regression)
        info("Regression Step:" + str(summary_regression))

        prf_values = {}
        for k, v in zip(self.base.rownames(summary_regression[0]), summary_regression[0][4]):
            if k.strip() != "Residuals":
                prf_values[k.strip()] = v

        return regression, prf_values

    def predict_best(self, regression, data):
        info("Predicting Best")
        predicted = self.stats.predict(regression, data)

        predicted_min = min(predicted)
        identical_predictions = 0

        for k in range(len(predicted)):
            if self.isclose(predicted[k], predicted_min, rel_tol = 1e-5):
                identical_predictions += 1

        info("Identical predictions (tol = 1e-5): {0}".format(identical_predictions))
        return data.rx(predicted.ro == self.base.min(predicted), True)

    def prune_data(self, data, predicted_best, fixed_variables):
        info("Pruning Data")
        conditions = []

        for k, v in fixed_variables.items():
            info("Predicted best column " + str(k) + ": " + str(predicted_best.rx2(str(k))))
            if conditions == []:
                conditions = data.rx2(str(k)).ro == predicted_best.rx2(str(k))
            else:
                conditions = conditions.ro & (data.rx2(str(k)).ro == predicted_best.rx2(str(k)))

        pruned_data = data.rx(conditions, True)

        info("Dimensions of Pruned Data: " + str(self.base.dim(pruned_data)).strip())
        info("Pruned data names: " + str(self.base.names(pruned_data)))
        info(str(self.utils.str(pruned_data)))
        return pruned_data

    def get_ordered_fixed_variables(self, ordered_keys, prf_values, threshold = 5, prf_threshold = 0.1):
        unique_variables = []
        info("prf_values in function call: " + str(prf_values))
        for k in ordered_keys:
            if k not in unique_variables and prf_values[str(k)] < prf_threshold:
                unique_variables.append(k)
            if len(unique_variables) > threshold:
                break

        return unique_variables

    def get_fixed_variables(self, predicted_best, ordered_prf_keys,
                            prf_values, fixed_factors):
        info("Getting fixed variables")
        variables = ordered_prf_keys
        #variables = [v.strip("I)(/1 ") for v in variables]

        unique_variables = self.get_ordered_fixed_variables(ordered_prf_keys, prf_values)

        fixed_variables = fixed_factors
        for v in unique_variables:
            fixed_variables[v] = predicted_best.rx2(1, str(v))[0]

        info("Fixed Variables: " + str(fixed_variables))
        return fixed_variables

    def prune_model(self, factors, inverse_factors, ordered_prf_keys,
                    prf_values):
        info("Pruning Model")
        variables = ordered_prf_keys
        #variables = [v.strip("I)(/1 ") for v in variables]

        info("prf_values in prune_data: " + str(prf_values))
        unique_variables = self.get_ordered_fixed_variables(ordered_prf_keys, prf_values)
        pruned_factors = [f for f in factors if not f in unique_variables]
        pruned_inverse_factors = [f for f in inverse_factors if not f in unique_variables]

        return pruned_factors, pruned_inverse_factors

    def get_federov_data(self, factors):
        low_level_limits  = IntVector([self.parameter_ranges[f][0] for f in factors])
        high_level_limits = IntVector([self.parameter_ranges[f][1] - 1 for f in factors])
        factor_centers    = IntVector([0 for f in factors])
        factor_levels     = IntVector([self.parameter_ranges[f][1] for f in factors])
        factor_round      = IntVector([0 for f in factors])
        # TODO Find a way to keep track of factor information
        is_factor         = BoolVector([False for f in factors])
        mix               = BoolVector([False for f in factors])

        opt_federov_data = {
                             "var": StrVector(factors),
                             "low": low_level_limits,
                             "high": high_level_limits,
                             "center": factor_centers,
                             "nLevels": factor_levels,
                             "round": factor_round,
                             "factor": is_factor,
                             "mix": mix
                           }

        opt_federov_dataframe = DataFrame(opt_federov_data)
        opt_federov_dataframe = opt_federov_dataframe.rx(StrVector(["var",
                                                                   "low",
                                                                   "high",
                                                                   "center",
                                                                   "nLevels",
                                                                   "round",
                                                                   "factor",
                                                                   "mix"]))
        return opt_federov_dataframe

    def get_updated_constraints(self, factors, fixed_variables):
        info("Updating Constraints")
        constraint_text = self.constraint
        variable_ranges = self.axis_val_ranges
        variable_names = self.params["axis_names"]

        for k, v in fixed_variables.items():
            current_value = str(variable_ranges[variable_names.index(k)][int(v)])
            constraint_text = constraint_text.replace(k, current_value)

        for i in range(len(factors)):
            current_value = ("variable_ranges[variable_names.index"
                             "(factors[{0}])][int(round(x[{0}]))]".format(i))

            old_value = variable_names[variable_names.index(factors[i])]
            constraint_text = constraint_text.replace(old_value, current_value)

        @ri.rternalize
        def constraint(x):
            type(variable_ranges)
            type(variable_names)
            type(factors)
            result = eval(constraint_text)
            return result

        info("Updated Constraint: " + str(self.constraint))

        return constraint

    def measure_design(self, design, response, fixed_factors):
        info("Measuring design of size " + str(len(design[0])))

        design_names    = [str(n) for n in self.base.names(design)]
        initial_factors = self.params["axis_names"]
        measurements    = []

        info("Current Design Names: " + str(design_names))
        info("Initial Factors: " + str(initial_factors))

        for line in range(len(design[0])):
            design_line = [int(v[0]) for v in design.rx(line + 1, True)]

            candidate = [0] * len(initial_factors)

            for k, v in fixed_factors.items():
                candidate[initial_factors.index(k)] = int(v)

            for i in range(len(design_names)):
                candidate[initial_factors.index(design_names[i])] = design_line[i]

            info("Initial Design Line: " + str(design_line))
            info("Fixed Factors: " + str(fixed_factors))
            info("Testing candidate " + str(line + 1) + ": " + str(candidate))

            measurement = self.getPerfCosts([candidate])
            if measurement != {}:
                measurements.append(float(numpy.mean(measurement[str(candidate)][0])))
            else:
                measurements.append(float('inf'))

        info("Measurements: " + str(measurements))

        design = self.base.cbind(design, DataFrame({response[0]: FloatVector(measurements)}))

        design = design.rx(self.base.is_finite(design.rx2(response[0])), True)

        info(str(design))

        return design

    def dopt_anova_step(self, response, factors, inverse_factors, step_space,
                        search_space, fixed_factors, budget):
        full_model     = "".join([" ~ ",
                                  " + ".join(factors)])

        # Leaving inverses out for now, since they do not work well with
        # Federov
        #
        # if len(inverse_factors) > 0:
        #     full_model += " + " + " + ".join(["I(1 / {0})".format(f) for f in
        #         inverse_factors])

        info(str(full_model))

        design_formula = full_model
        lm_formula     = response[0] + full_model
        trials         = int(2 * len(factors))

        fixed_variables = fixed_factors
        info("Fixed Factors: " + str(fixed_factors))

        if budget - len(step_space[0]) < 0 and trials < len(step_space[0]):
            info("Full data does not fit on budget")
            info("Computing D-Optimal Design")
            constraint = self.get_updated_constraints(factors, fixed_variables)

            info("Computing D-Optimal Design with " + str(trials) +
                 " experiments")
            info("Design Formula: " + str(design_formula))

            opt_federov_dataframe = self.get_federov_data(factors)

            output = self.opt_federov(design_formula, trials, constraint,
                                      opt_federov_dataframe)

            design = output.rx("design")[0]

            info(str(design))
            info("D-Efficiency Approximation: " + str(output.rx("Dea")[0]))

            design = self.measure_design(design, response, fixed_factors)

            info("Step Space Names: " + str(self.base.names(step_space)))

            used_experiments       = len(design[0])
            regression, prf_values = self.anova(design, lm_formula)
            ordered_prf_keys       = sorted(prf_values, key = prf_values.get)
            predicted_best         = self.predict_best(regression, step_space)
            info("prf_values: " + str(prf_values))
            fixed_variables        = self.get_fixed_variables(predicted_best, ordered_prf_keys,
                                                              prf_values, fixed_factors)

            info("prf_values after get vars: " + str(prf_values))
            pruned_space = self.prune_data(step_space, predicted_best, fixed_variables)
            pruned_factors, pruned_inverse_factors = self.prune_model(factors, inverse_factors,
                                                                      ordered_prf_keys, prf_values)
        else:
            info(("Full data fits on budget, or too few data points"
                  " for a D-Optimal design. Picking best value."))

            used_experiments = len(step_space[0])
            prf_values = []
            ordered_prf_keys = []
            pruned_space = []
            pruned_factors = []
            pruned_inverse_factors = []

            step_data = self.measure_design(step_space, response, fixed_factors)
            predicted_best = step_data.rx((step_data.rx2(response[0]).ro == min(step_data.rx(response[0])[0])),
                                      True)

        info("Best Predicted: " + str(predicted_best))
        info("Pruned Factors: " + str(pruned_factors))
        info("Fixed Factors: " + str(fixed_variables))

        return {"prf_values": prf_values,
                "ordered_prf_keys": ordered_prf_keys,
                "predicted_best": predicted_best,
                "pruned_space": pruned_space,
                "pruned_factors": pruned_factors,
                "pruned_inverse_factors": pruned_inverse_factors,
                "fixed_factors": fixed_variables,
                "used_experiments": used_experiments}


    def dopt_anova(self, initial_factors, initial_inverse_factors, search_space):
        response = ["cost_mean"]

        step_factors = initial_factors
        step_inverse_factors = initial_inverse_factors
        step_space = search_space

        fixed_factors = {}

        initial_budget = 180
        budget = initial_budget
        used_experiments = 0
        iterations = 10

        for i in range(iterations):
            info("Step {0}".format(i))
            if step_space == []:
                break

            step_data = self.dopt_anova_step(response,
                                             step_factors,
                                             step_inverse_factors,
                                             step_space,
                                             search_space,
                                             fixed_factors,
                                             budget)

            step_space = step_data["pruned_space"]
            step_factors = step_data["pruned_factors"]
            step_inverse_factors = step_data["pruned_inverse_factors"]
            budget -= step_data["used_experiments"]
            used_experiments += step_data["used_experiments"]
            fixed_factors = step_data["fixed_factors"]

            info("Fixed Factors: " + str(fixed_factors))

            starting_point = numpy.mean((self.getPerfCosts([[0] * self.total_dims]).values()[0])[0])
            info(str(starting_point))

            predicted_best = [int(v[0]) for v in step_data["predicted_best"].rx(1, True)]
            info(str(predicted_best))
            predicted_best_value = numpy.mean((self.getPerfCosts([predicted_best]).values()[0])[0])

            info("Slowdown: " + str(predicted_best_value / starting_point))
            info("Budget: {0}/{1}".format(used_experiments, initial_budget))

        return True

    def searchBestCoord(self, startCoord = None):
        '''
        To explore the search space and retun the coordinate that yields the best performance
        (i.e. minimum performance cost).
        '''
        info('\n----- begin random search -----')

        info(str(self.params["axis_names"]))
        info(str(self.total_dims))

        initial_factors = self.params["axis_names"]
        initial_inverse_factors = initial_factors

        best_coord = None
        best_perf_cost = self.MAXFLOAT
        old_perf_cost = best_perf_cost

        # record the number of runs
        runs = 0
        sruns = 0
        fruns = 0
        start_time = time.time()

        full_candidate_set = {}
        search_space = []

        self.seed_space_size = 400000

        info("Building seed search space (does not spend evaluations)")
        if not os.path.isfile("search_space_{0}.db".format(self.seed_space_size)):
            info("Generating new search space for this size")
            search_space_database = dataset.connect("sqlite:///search_space_{0}.db".format(self.seed_space_size))
            experiments = search_space_database['experiments']

            while len(search_space) < self.seed_space_size:
                if len(full_candidate_set) % 20000 == 0:
                    info("Evaluated coordinates: " + str(len(full_candidate_set)))

                candidate_point = self.getRandomCoord()
                candidate_point_key = str(candidate_point)

                if candidate_point_key not in full_candidate_set:
                    full_candidate_set[candidate_point_key] = candidate_point
                    try:
                        perf_params = self.coordToPerfParams(candidate_point)
                        is_valid = eval(self.constraint, copy.copy(perf_params),
                                        dict(self.input_params))
                    except Exception, e:
                        err('failed to evaluate the constraint expression: "%s"\n%s %s'
                            % (self.constraint, e.__class__.__name__, e))

                    if is_valid:
                        experiments.insert({"value": str(candidate_point)})
                        search_space.append(candidate_point)
                        if len(search_space) % 5000 == 0:
                            info("Valid coordinates: " + str(len(search_space)))

            info("Valid/Tested configurations: " + str(len(search_space)) + "/" +
                str(len(full_candidate_set)))
        else:
            info("Using pre-generated space for this size")
            search_space_database = dataset.connect("sqlite:///search_space_{0}.db".format(self.seed_space_size))
            for experiment in search_space_database['experiments']:
                search_space.append(eval(experiment["value"]))

        info("Starting DOPT-anova")

        r_search_space = {}
        for i in range(len(search_space[0])):
            r_row = [self.dim_uplimits[i] - 1, 0]
            for col in search_space:
                r_row.append(col[i])

            r_search_space[initial_factors[i]] = IntVector(r_row)

        data = DataFrame(r_search_space)
        data = data.rx(StrVector(initial_factors))

        self.dopt_anova(initial_factors, initial_inverse_factors, data)

        sys.exit()

        perf_cost, mean_perf_cost = self.MAXFLOAT, self.MAXFLOAT

        params = self.coordToPerfParams(coord)
        end_time = time.time()
        search_time = start_time - end_time
        speedup = float(eval_cost[0]) / float(best_perf_cost)
        search_time = time.time() - start_time

        info('----- end random search -----')

        info('----- begin random search summary -----')
        info(' total completed runs: %s' % runs)
        info(' total successful runs: %s' % sruns)
        info(' total failed runs: %s' % fruns)
        info(' speedup: %s' % speedup)
        info(' fount at: %s' % num_eval_best)
        info('----- end random search summary -----')

        # return the best coordinate
        return best_coord, best_perf_cost, search_time, sruns

    # Private methods

    def __readAlgoArgs(self):
        '''To read all algorithm-specific arguments'''

        # check for algorithm-specific arguments
        for vname, rhs in self.search_opts.iteritems():
            print vname, rhs
            # local search distance
            if vname == self.__LOCAL_DIST:
                if not isinstance(rhs, int) or rhs < 0:
                    err('orio.main.tuner.search.randomsearch: %s argument "%s" must be a positive integer or zero'
                        % (self.__class__.__name__, vname))
                self.local_distance = rhs

            # unrecognized algorithm-specific argument
            elif vname == 'total_runs':
                self.total_runs = rhs
            else:
                err('orio.main.tuner.search.randomsearch: unrecognized %s algorithm-specific argument: "%s"'
                    % (self.__class__.__name__, vname))
