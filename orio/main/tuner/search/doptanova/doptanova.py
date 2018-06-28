import sys, time
import math
import random
import numpy
import scipy.linalg
import orio.main.tuner.search.search
from orio.main.util.globals import *
import copy
import json
import dataset
import os

import rpy2.rinterface as ri
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import DataFrame, IntVector, FloatVector, StrVector, BoolVector, Formula, NULL, r

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
        #numpy.random.seed(39920)
        #numpy.random.seed(99291)
        numpy.random.seed(22010)

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

        self.parameter_values = {}

        for i in range(len(self.params["axis_val_ranges"])):
            self.parameter_values[self.params["axis_names"][i]] = self.params["axis_val_ranges"][i]

        info("Parameter Range Values: " + str(self.parameter_values))

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


    def clean_search_space(self, federov_search_space, factors, inverse_factors, fixed_factors):
        data = {}
        r_snippet = """data <- %s
        clean_data <- Filter(function(x)(length(unique(x)) > 1), data)
        removed_factors <- names(Filter(function(x)(length(unique(x)) == 1), data))
        removed_inverse_factors <- names(Filter(function(x)(length(unique(x)) == 2), data))

        list("clean_data" = clean_data, "removed_factors" = removed_factors, "removed_inverse_factors" = removed_inverse_factors)
        """ % (federov_search_space.r_repr())

        output = robjects.r(r_snippet)

        removed_factors = output[1]
        removed_inverse_factors = output[2]

        info("Clean Info:")
        info("Removed Factors: " + str(removed_factors))
        info("Removed Inverse Factors: " + str(removed_inverse_factors))

        factors = [f for f in factors if f not in removed_factors]
        inverse_factors = [f for f in inverse_factors if f not in removed_inverse_factors + removed_factors]

        for f in removed_factors:
            fixed_factors[f] = int(federov_search_space.rx(1, f)[0])

        info("New Factors: " + str(factors))
        info("New Inverse Factors: " + str(inverse_factors))
        info("New Fixed Factors: " + str(fixed_factors))

        data["search_space"]    = output[0]
        data["factors"]         = factors
        data["inverse_factors"] = inverse_factors
        data["fixed_factors"]   = fixed_factors
        return data

    def generate_valid_sample(self, sample_size, fixed_variables):
        search_space_dataframe = {}

        for n in self.axis_names:
            search_space_dataframe[n] = []

        search_space = {}
        evaluated = 0

        info("Generating valid search space of size {0} (does not spend evaluations)".format(sample_size))

        while len(search_space) < sample_size:
            candidate_point      = self.getRandomCoord()
            candidate_point_key  = str(candidate_point)
            evaluated           += 1

            if candidate_point_key not in search_space:
                perf_params = self.coordToPerfParams(candidate_point)

                for k, v in fixed_variables.items():
                    perf_params[k] = self.parameter_values[k][int(v)]

                is_valid = eval(self.constraint, copy.copy(perf_params),
                                dict(self.input_params))

                if is_valid:
                    search_space[candidate_point_key] = candidate_point

                    for n in perf_params:
                        candidate_value = self.parameter_values[n].index(perf_params[n])
                        search_space_dataframe[n].append(candidate_value)

                    if len(search_space) % int(sample_size / 5) == 0:
                        info("Valid coordinates: " + str(len(search_space)))

                if evaluated % 10000 == 0:
                    info("Tested coordinates: " + str(evaluated))

        info("Valid/Tested configurations: " + str(len(search_space)) + "/" +
             str(evaluated))

        for k in search_space_dataframe:
            search_space_dataframe[k] = IntVector(search_space_dataframe[k])

        search_space_dataframe_r = DataFrame(search_space_dataframe)
        search_space_dataframe_r = search_space_dataframe_r.rx(StrVector(self.axis_names))

        info("Generated Search Space:")
        info(str(self.utils.str(search_space_dataframe_r)))

        return search_space_dataframe_r

    def opt_monte(self, design_formula, trials, constraint, data):
        info("Starting \"optMonteCarlo\" run")
        info(str(data))

        #self.base.set_seed(77126)
        #self.base.set_seed(66182)
        self.base.set_seed(22331)

        candidate_multiplier = 10
        repetitions          = 1

        output = self.algdesign.optMonteCarlo(frml        = Formula(design_formula),
                                              data        = data,
                                              constraints = constraint,
                                              nCand       = candidate_multiplier * trials,
                                              nRepeats    = repetitions,
                                              nTrials     = trials)
        return output

    def opt_federov(self, design_formula, trials, data):
        info("Starting \"optFederov\" run")
        info("Using Search Space:")
        info(str(self.utils.str(data)))

        #self.base.set_seed(77126)
        #self.base.set_seed(66182)
        self.base.set_seed(22331)

        #info("Correlation between variables in the dataset:")
        #info(str(self.stats.cor(data)))

        info("Data Dimensions: " + str(self.base.dim(data)))

        output = self.algdesign.optFederov(frml         = Formula(design_formula),
                                           data         = data,
                                           nTrials      = trials,
                                           maxIteration = 100000,
                                           center       = True,
                                           nullify      = 2)

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

    def predict_best(self, regression, size, fixed_variables):
        info("Predicting Best")
        data = self.generate_valid_sample(size, fixed_variables)
        predicted = self.stats.predict(regression, data)

        predicted_min = min(predicted)
        identical_predictions = 0

        for k in range(len(predicted)):
            if self.isclose(predicted[k], predicted_min, rel_tol = 1e-5):
                identical_predictions += 1

        info("Identical predictions (tol = 1e-5): {0}".format(identical_predictions))
        return data.rx(predicted.ro == self.base.min(predicted), True)

    def get_design_best(self, design, response):
        info("Getting Best from Design")
        info("Response: " + str(response))
        info("Design Names: " + str(self.base.names(design)))
        info("Design Response: " + str(design.rx(str(response[0]))))
        return design.rx(design.rx2(str(response[0])).ro == self.base.min(design.rx(str(response[0]))), True)

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

    def get_ordered_fixed_variables(self, ordered_keys, prf_values, threshold = 7, prf_threshold = 0.1):
        ordered_keys     = [k.replace("I(1/(1 + ", "").strip(") ") for k in ordered_keys]
        unique_variables = []
        for k in ordered_keys:
            if k not in unique_variables and prf_values[str(k)] < prf_threshold:
                unique_variables.append(k)
            if len(unique_variables) >= threshold:
                break

        if unique_variables == []:
            unique_variables.append(ordered_keys[0])

        return unique_variables

    def get_fixed_variables(self, predicted_best, ordered_prf_keys,
                            prf_values, fixed_factors):
        info("Getting fixed variables")
        unique_variables = self.get_ordered_fixed_variables(ordered_prf_keys, prf_values)
        fixed_variables  = fixed_factors
        for v in unique_variables:
            fixed_variables[v] = predicted_best.rx2(1, str(v))[0]

        return fixed_variables

    def prune_model(self, factors, inverse_factors, ordered_prf_keys,
                    prf_values):
        info("Pruning Model")
        unique_variables       = self.get_ordered_fixed_variables(ordered_prf_keys, prf_values)
        pruned_factors         = [f for f in factors if not f in unique_variables]
        pruned_inverse_factors = [f for f in inverse_factors if not f in unique_variables]
        return pruned_factors, pruned_inverse_factors

    def get_federov_data(self, factors):
        low_level_limits  = IntVector([self.parameter_ranges[f][0] for f in factors])
        high_level_limits = IntVector([self.parameter_ranges[f][1] - 1 for f in factors])
        factor_centers    = IntVector([0 for f in factors])
        factor_levels     = IntVector([self.parameter_ranges[f][1] for f in factors])
        factor_round      = IntVector([0 for f in factors])
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
        variable_ranges = self.parameter_values

        for k, v in fixed_variables.items():
            current_value = str(variable_ranges[k][int(v)])
            constraint_text = constraint_text.replace(k, current_value)

        for i in range(len(factors)):
            current_value = ("variable_ranges[\"{0}\"][int(x[{1}])]".format(factors[i], i))

            old_value = factors[i]
            constraint_text = constraint_text.replace(old_value, current_value)

        info("Variable Ranges: " + str(variable_ranges))
        info("Current Factors: " + str(factors))
        info("Constraint Text: " + str(self.constraint))

        @ri.rternalize
        def constraint(x):
            type(variable_ranges)
            result = eval(constraint_text)
            return result

        info("Updated Constraint: " + str(constraint_text))
        return constraint

    def measure_design(self, design, response, fixed_factors):
        info("Measuring design of size " + str(len(design[0])))

        design_names    = [str(n) for n in self.base.names(design)]
        initial_factors = self.params["axis_names"]
        measurements    = []

        info("Current Design Names: " + str(design_names))
        info("Initial Factors: " + str(initial_factors))

        for line in range(1, len(design[0]) + 1):
            design_line = [int(v[0]) for v in design.rx(line, True)]

            candidate = [0] * len(initial_factors)

            for k, v in fixed_factors.items():
                candidate[initial_factors.index(k)] = int(v)

            for i in range(len(design_names)):
                candidate[initial_factors.index(design_names[i])] = design_line[i]

            # info("Initial Design Line: " + str(design_line))
            # info("Fixed Factors: " + str(fixed_factors))
            # info("Testing candidate " + str(line + 1) + ": " + str(candidate))

            measurement = self.getPerfCosts([candidate])
            if measurement != {}:
                measurements.append(float(numpy.mean(measurement[str(candidate)][0])))
            else:
                measurements.append(float('inf'))

        # info("Measurements: " + str(measurements))

        design = self.base.cbind(design, DataFrame({response[0]: FloatVector(measurements)}))

        design = design.rx(self.base.is_finite(design.rx2(response[0])), True)

        info("Complete design, with measurements:")
        info(str(design))

        return design

    def dopt_anova_step(self, response, factors, inverse_factors,
                        fixed_factors, budget, step_number):
        trials = int(2 * (len(factors) + len(inverse_factors)))

        federov_samples = 20 * trials
        prediction_samples = 4 * federov_samples

        federov_search_space = self.generate_valid_sample(federov_samples, fixed_factors)
        federov_search_space = federov_search_space.rx(StrVector(factors))

        clean_search_space_data = self.clean_search_space(federov_search_space,
                                                          factors,
                                                          inverse_factors,
                                                          fixed_factors)

        federov_search_space = clean_search_space_data["search_space"]
        factors              = clean_search_space_data["factors"]
        inverse_factors      = clean_search_space_data["inverse_factors"]
        fixed_factors        = clean_search_space_data["fixed_factors"]

        full_model     = "".join([" ~ ",
                                  " + ".join(factors)])

        if len(inverse_factors) > 0:
            full_model += " + " + " + ".join(["I(1 / (1 + {0}))".format(f) for f in
                inverse_factors])

        info("Full Model: " + str(full_model))

        design_formula = full_model
        lm_formula     = response[0] + full_model

        fixed_variables = fixed_factors
        info("Fixed Factors: " + str(fixed_factors))
        info("Computing D-Optimal Design")
        #constraint = self.get_updated_constraints(factors, fixed_variables)

        info("Computing D-Optimal Design with " + str(trials) +
             " experiments")
        info("Design Formula: " + str(design_formula))

        # opt_federov_dataframe = self.get_federov_data(factors)
        # output = self.opt_monte(design_formula, trials, constraint,
        #                         opt_federov_dataframe)

        # X = self.stats.model_matrix(Formula(design_formula), federov_search_space)
        # info("Model Matrix Names: " + str(self.base.names(X)))
        # X_t = self.base.t(X)
        # X_I = numpy.matmul(X_t, X)
        # numpy.set_printoptions(precision = 2, linewidth = 120)

        # X_I_P, X_I_L, X_I_U = scipy.linalg.lu(X_I)
        # determinant = numpy.linalg.det(X_I)

        # info("Model Matrix Det: " + str(determinant))
        # info("LU Decomposition L:")
        # info(str(X_I_L))
        # info("LU Decomposition U:")
        # info(str(X_I_U))
        # info("LU Decomposition P:")
        # info(str(X_I_P))

        # X_I = numpy.array(X_I)

        # info("Computing Dependency Between Columns")
        # for i in range(X_I.shape[1]):
        #     for j in range(X_I.shape[1]):
        #         if i != j:
        #             inner_product = numpy.inner(
        #                     X_I[:,i],
        #                     X_I[:,j]
        #                     )
        #             norm_i = numpy.linalg.norm(X_I[:,i])
        #             norm_j = numpy.linalg.norm(X_I[:,j])

        #             if numpy.abs(inner_product - norm_j * norm_i) < 1E-5:
        #                 info("(i: " + str(i) + ", j: " + str(j) + ")")
        #                 info('Dependent')
        #                 info('I: ' + str(X_I[:,i]))
        #                 info('J: ' + str(X_I[:,j]))
        #                 info('Prod: ' + str(inner_product))
        #                 info('Norm i: ' + str(norm_i))
        #                 info('Norm j: ' + str(norm_j))

        output = self.opt_federov(design_formula, trials, federov_search_space)

        design = output.rx("design")[0]

        info(str(design))
        info("D-Efficiency Approximation: " + str(output.rx("Dea")[0]))

        design = self.measure_design(design, response, fixed_factors)

        used_experiments       = len(design[0])
        regression, prf_values = self.anova(design, lm_formula)
        ordered_prf_keys       = sorted(prf_values, key = prf_values.get)
        predicted_best         = self.predict_best(regression, prediction_samples, fixed_variables)
        design_best            = self.get_design_best(design, response)
        fixed_variables        = self.get_fixed_variables(predicted_best, ordered_prf_keys,
                                                          prf_values, fixed_factors)

        pruned_factors, pruned_inverse_factors = self.prune_model(factors, inverse_factors,
                                                                  ordered_prf_keys, prf_values)

        info("Best Predicted: " + str(predicted_best))
        info("Best From Design: " + str(design_best))
        info("Pruned Factors: " + str(pruned_factors))
        info("Fixed Factors: " + str(fixed_variables))

        return {"prf_values": prf_values,
                "ordered_prf_keys": ordered_prf_keys,
                "design_best": design_best,
                "predicted_best": predicted_best,
                "pruned_factors": pruned_factors,
                "pruned_inverse_factors": pruned_inverse_factors,
                "fixed_factors": fixed_variables,
                "used_experiments": used_experiments}


    def dopt_anova(self, initial_factors, initial_inverse_factors):
        response = ["cost_mean"]

        step_factors = initial_factors
        step_inverse_factors = initial_inverse_factors

        # iterations = int((len(initial_factors) + len(initial_inverse_factors)) / 2)
        iterations = 4

        fixed_factors = {}

        initial_budget = 180
        budget = initial_budget
        used_experiments = 0
        best_value = float("inf")
        best_point = []

        for i in range(iterations):
            info("Step {0}".format(i))

            step_data = self.dopt_anova_step(response,
                                             step_factors,
                                             step_inverse_factors,
                                             fixed_factors,
                                             budget,
                                             i)

            step_factors          = step_data["pruned_factors"]
            step_inverse_factors  = step_data["pruned_inverse_factors"]
            budget               -= step_data["used_experiments"]
            used_experiments     += step_data["used_experiments"]
            fixed_factors         = step_data["fixed_factors"]

            starting_point = numpy.mean((self.getPerfCosts([[0] * self.total_dims]).values()[0])[0])
            info("Starting Point (-O3):")
            info(str(starting_point))

            predicted_best = [int(v[0]) for v in step_data["predicted_best"].rx(1, True)]
            info("Predicted Best Point:")
            info(str(predicted_best))
            predicted_best_value = numpy.mean((self.getPerfCosts([predicted_best]).values()[0])[0])

            design_best = step_data["design_best"]
            info("Design Best Point:")
            info(str(design_best))

            info("Fixed Factors: " + str(fixed_factors))

            design_best_slowdown    = float(design_best.rx(1, response[0])[0]) / starting_point
            predicted_best_slowdown = predicted_best_value / starting_point

            info("Slowdown (Design Best): " + str(float(design_best.rx(1, response[0])[0]) / starting_point))
            info("Slowdown (Predicted Best): " + str(predicted_best_value / starting_point))
            info("Budget: {0}/{1}".format(used_experiments, initial_budget))

            if design_best_slowdown < predicted_best_slowdown:
                current_best = design_best
                current_best_value = float(design_best.rx(1, response[0])[0])
            else:
                current_best = predicted_best
                current_best_value = predicted_best_value

            if current_best_value < best_value or best_point == []:
                info("Updating Global Best Slowdown: " + str(current_best_value))
                best_point = current_best
                best_value = current_best_value

        target_names = [n for n in self.base.names(best_point) if n != response[0]]
        info("Target Names: " + str(target_names))

        design_line = [int(v[0]) for v in best_point.rx(1, StrVector(target_names))]
        candidate = [0] * len(initial_factors)

        for k, v in fixed_factors.items():
            candidate[initial_factors.index(k)] = int(v)

        for i in range(len(target_names)):
            candidate[initial_factors.index(target_names[i])] = design_line[i]

        return candidate

    def searchBestCoord(self, startCoord = None):
        '''
        To explore the search space and retun the coordinate that yields the best performance
        (i.e. minimum performance cost).
        '''
        info('\n----- begin random search -----')

        initial_factors = self.params["axis_names"]
        initial_inverse_factors = [f for f in initial_factors if self.parameter_ranges[f][1] > 2]
        #initial_inverse_factors = []

        info("Initial Factors: " + str(initial_factors))
        info("Initial Inverse Factors: " + str(initial_inverse_factors))

        best_coord = None
        best_perf_cost = self.MAXFLOAT
        old_perf_cost = best_perf_cost

        # record the number of runs
        runs = 0
        sruns = 0
        fruns = 0
        start_time = time.time()

        info("Starting DOPT-anova")

        best_point = self.dopt_anova(initial_factors, initial_inverse_factors)

        info("Ending DOPT-ANOVA")
        info("Best Point: " + str(best_point))

        predicted_best_value = numpy.mean((self.getPerfCosts([best_point]).values()[0])[0])
        starting_point       = numpy.mean((self.getPerfCosts([[0] * self.total_dims]).values()[0])[0])
        speedup              = starting_point / predicted_best_value

        end_time = time.time()
        search_time = start_time - end_time

        info("Speedup: " + str(speedup))

        info('----- end random search -----')

        info('----- begin random search summary -----')
        info(' total completed runs: %s' % runs)
        info(' total successful runs: %s' % sruns)
        info(' total failed runs: %s' % fruns)
        info(' speedup: %s' % speedup)
        info('----- end random search summary -----')

        # return the best coordinate
        return best_point, predicted_best_value, search_time, sruns

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
