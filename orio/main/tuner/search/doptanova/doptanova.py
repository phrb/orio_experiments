import sys, time
import math
import random
import numpy
import orio.main.tuner.search.search
from orio.main.util.globals import *
import copy
import json

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
        numpy.random.seed(7777)

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

        self.base.set_seed(1337)
        output = self.algdesign.optMonteCarlo(frml = Formula(design_formula),
                                              data = data,
                                              nCand = 10 * trials,
                                              constraints = constraint,
                                              nTrials = trials)

        return output

    def transform_lm(self, design, lm_formula):
        info("Power Transform Step:")
        response = lm_formula.split("~")[0].strip()
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

        if heteroscedasticity_test.rx("p")[0][0] < 0.05:
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
            if conditions == []:
                conditions = data.rx2(str(k)).ro == predicted_best.rx2(str(k))
            else:
                conditions = conditions.ro & (data.rx2(str(k)).ro == predicted_best.rx2(str(k)))

        pruned_data = data.rx(conditions, True)

        info("Dimensions of Pruned Data: " + str(self.base.dim(pruned_data)).strip())
        return pruned_data

    def get_fixed_variables(self, predicted_best, ordered_prf_keys,
                            fixed_factors, threshold = 2):
        info("Getting fixed variables")
        variables = ordered_prf_keys
        variables = [v.strip("I)(/1 ") for v in variables]

        unique_variables = []

        for v in variables:
            if v not in unique_variables:
                unique_variables.append(v)
            if len(unique_variables) >= threshold:
                break

        fixed_variables = fixed_factors
        for v in unique_variables:
            fixed_variables[v] = predicted_best.rx2(str(v))[0]

        info("Fixed Variables: " + str(fixed_variables))
        return fixed_variables

    def prune_model(self, factors, inverse_factors, ordered_prf_keys,
                    threshold = 2):
        info("Pruning Model")
        variables = ordered_prf_keys
        variables = [v.strip("I)(/1 ") for v in variables]

        unique_variables = []

        for v in variables:
            if v not in unique_variables:
                unique_variables.append(v)
            if len(unique_variables) >= threshold:
                break

        pruned_factors = [f for f in factors if not f in unique_variables]
        pruned_inverse_factors = [f for f in inverse_factors if not f in unique_variables]

        return pruned_factors, pruned_inverse_factors

    def dopt_anova_step(self, response, factors, inverse_factors, # data, step_data,
                        fixed_factors, budget):
        full_model     = "".join([" ~ ",
                                  " + ".join(factors)])

        #if len(inverse_factors) > 0:
        #    full_model += " + " + " + ".join(["I(1 / {0})".format(f) for f in
        #        inverse_factors])

        #info(str(full_model))

        low_level_limits = IntVector([self.parameter_ranges[f][0] for f in factors])
        high_level_limits = IntVector([self.parameter_ranges[f][1] - 1 for f in factors])
        factor_centers = IntVector([0 for f in factors])
        factor_levels = IntVector([self.parameter_ranges[f][1] for f in factors])
        factor_round = IntVector([0 for f in factors])
        # Find a way to keep track of factor information
        is_factor = BoolVector([False for f in factors])
        mix = BoolVector([False for f in factors])

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

        #info(str(opt_federov_dataframe))

        design_formula = full_model
        lm_formula     = response[0] + full_model
        trials         = int(round(2 * (len(factors) + len(inverse_factors) + 1)))
        trials         = int(2 * len(factors))

        fixed_variables = fixed_factors
        info("Fixed Factors: " + str(fixed_factors))

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
            # Using variables here so they get "rternalized"
            # There must be a better way
            type(variable_ranges)
            type(variable_names)
            type(factors)
            result = eval(constraint_text)
            return result

        info("Updated Constraint: " + str(self.constraint))

        info("Computing D-Optimal Design with " + str(trials) +
             " experiments")
        info("Design Formula: " + str(design_formula))

        output = self.opt_federov(design_formula, trials, constraint,
                                  opt_federov_dataframe)

        design = output.rx("design")[0]

        info(str(design))
        info("D-Efficiency Approximation: " + str(output.rx("Dea")[0]))
        info("Measuring design of size " + str(len(design)))

        measurements = []

        for line in range(len(design[0])):
            #TODO Update design line with fixed parameters!
            candidate = [int(v[0]) for v in design.rx(line + 1, True)]
            info("Testing candidate " + str(line + 1) + ": " + str(candidate))

            measurement = self.getPerfCosts([candidate])
            if measurement != {}:
                measurements.append(float(numpy.mean(measurement[str(candidate)][0])))
            else:
                measurements.append(float('inf'))

        info("Measurements: " + str(measurements))

        design = self.base.cbind(design, DataFrame({response[0]: FloatVector(measurements)}))

        info(str(design))
        used_experiments = len(design[0])
        regression, prf_values = self.anova(design, lm_formula)
        ordered_prf_keys       = sorted(prf_values, key = prf_values.get)
        predicted_best         = self.predict_best(regression, step_data)
        fixed_variables        = self.get_fixed_variables(predicted_best, ordered_prf_keys,
                                                          fixed_factors)

        pruned_factors, pruned_inverse_factors = self.prune_model(factors, inverse_factors,
                                                                  ordered_prf_keys)
        result = {"prf_values": prf_values,
                  "ordered_prf_keys": ordered_prf_keys,
                  "predicted_best": predicted_best,
                  # "pruned_data": pruned_data,
                  "pruned_factors": pruned_factors,
                  "pruned_inverse_factors": pruned_inverse_factors,
                  "fixed_factors": fixed_variables,
                  "used_experiments": used_experiments}

        info(str(result))

        sys.exit()

        return result

    def dopt_anova(self, initial_factors, initial_inverse_factors):
        response = ["cost_mean"]

        #info(str(initial_factors))

        #data = data.rx(StrVector(initial_factors))
        #data = data.rx(StrVector(initial_factors + response))
        #data_best = data.rx((data.rx2(response[0]).ro == min(data.rx(response[0])[0])),
        #                    True)

        step_factors = initial_factors
        step_inverse_factors = initial_inverse_factors
        # step_space = data

        fixed_factors = {}

        initial_budget = 1
        budget = initial_budget
        used_experiments = 0
        iterations = 1

        for i in range(iterations):
            info("Step {0}".format(i))
            # if step_space == []:
            #     break

            step_data = self.dopt_anova_step(response,
                                             step_factors,
                                             step_inverse_factors,
            #                                 data,
            #                                 step_space,
                                             fixed_factors,
                                             budget)

            # step_space = step_data["pruned_data"]
            step_factors = step_data["pruned_factors"]
            step_inverse_factors = step_data["pruned_inverse_factors"]
            budget -= step_data["used_experiments"]
            used_experiments += step_data["used_experiments"]
            fixed_factors = step_data["fixed_factors"]

            info("Fixed Factors: " + str(fixed_factors))

            # if step_space != []:
            #     step_best = step_space.rx((step_space.rx2(response[0]).ro ==
            #         min(step_space.rx(response[0])[0])), True)

            #     info("Best Step Slowdown: " +
            #             str(step_best.rx(response[0])[0][0] /
            #                 data_best.rx(response[0])[0][0]))

            info("Slowdown: " +
                    str(step_data["predicted_best"].rx(response[0])[0][0] /
                        data_best.rx(response[0])[0][0]))
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

        # d = {"constraints": self.constraint}
        # d.update(perf_params)
        # d.update(dict(self.input_params))
        # exec "def f(): return eval(constraints)" in d

        best_coord = None
        best_perf_cost = self.MAXFLOAT
        old_perf_cost = best_perf_cost

        # record the number of runs
        runs = 0
        sruns = 0
        fruns = 0
        start_time = time.time()

        # Trying to generate experiments at random:
        #
        # full_candidate_set = {}
        # search_space = []

        # self.seed_space_size = 1000

        # info("Building seed search space (does not spend evaluations)")
        # while len(search_space) < self.seed_space_size:
        #     if len(full_candidate_set) % 20000 == 0:
        #         info("Evaluated coordinates: " + str(len(full_candidate_set)))

        #     candidate_point = self.getRandomCoord()
        #     candidate_point_key = str(candidate_point)

        #     if candidate_point_key not in full_candidate_set:
        #         full_candidate_set[candidate_point_key] = candidate_point
        #         try:
        #             perf_params = self.coordToPerfParams(candidate_point)
        #             is_valid = eval(self.constraint, copy.copy(perf_params),
        #                             dict(self.input_params))
        #         except Exception, e:
        #             err('failed to evaluate the constraint expression: "%s"\n%s %s'
        #                 % (self.constraint, e.__class__.__name__, e))

        #         if not is_valid:
        #             continue

        #         search_space.append(candidate_point)
        #         if len(search_space) % 10000 == 0:
        #             info("Valid coordinates: " + str(len(search_space)))

        # info("Valid/Tested configurations: " + str(len(search_space)) + "/" +
        #     str(len(full_candidate_set)))

        # info("Starting DOPT-anova")

        # r_search_space = {}
        # for i in range(len(search_space[0])):
        #     r_row = [self.dim_uplimits[i] - 1, 0]
        #     for col in search_space:
        #         r_row.append(col[i])

        #     r_search_space[initial_factors[i]] = IntVector(r_row)

        # data = DataFrame(r_search_space)

        # info(str(search_space))
        # info(str(self.utils.str(data.rx(StrVector(initial_factors)))))

        self.dopt_anova(initial_factors, initial_inverse_factors)

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
