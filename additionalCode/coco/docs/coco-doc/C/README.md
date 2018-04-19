The COCO/NumBBO experiments interface
=====================================

<a href="https://github.com/numbbo/coco">COCO (COmparing Continuous Optimisers)</a> is a platform 
for systematic and sound comparisons of real-parameter global optimizers mainly developed within the 
<a href="http://numbbo.gforge.inria.fr/doku.php">NumBBO project</a>. COCO provides benchmark function 
testbeds, experimentation templates which are easy to parallelize, and tools for processing and 
visualizing data generated by one or several optimizers.

For a getting started guide see [here](https://github.com/numbbo/coco/blob/master/README.md#getting-started). 

__Reimplementation of COCO in ANSI C__

In order to allow for easier maintenance and further extensions of the COCO platform, it was rewritten
entirely from 2014 till 2016. Now, a single implementation in ANSI C (aka C89) is used and called from
the other languages to conduct the experiments. This documentation of the COCO C code serves therefore
as the basic reference for:
- [How to conduct benchmarking experiments in C](#benchmarking)
- [How to write new test problems and combine them into test suites](#new-suites)
- [How to write additional performance indicators and logging functionality](#new-indicators)

__Pointers to the source code and other documentation__

The actual source code of COCO can be found at http://github.com/numbbo/coco

More information about the biobjective test suite (bbob-biobj) can be found at
http://numbbo.github.io/coco-doc/bbob-biobj/functions/

The experimental procedure is described in
http://numbbo.github.io/coco-doc/experimental-setup/

## How to conduct benchmarking experiments in C <a name="benchmarking"></a>

The best way to create a benchmark experiment is to copy the example experiment and
make the required changes to include the chosen optimizer. 

In order to simplify the interface between the optimizers and the COCO platform, a static pointer
to a COCO problem and a function type for evaluation functions are used:

    static coco_problem_t *PROBLEM;
    typedef void (*evaluate_function_t)(const double *x, double *y);

A simplified version of benchmarking a single run of the  algorithm ``my_optimizer`` on the ``bbob-biobj`` suite with default 
parameters is invoked in the following way (see below for explanation of the 
[suite parameters](#suite-parameters) and [observer parameters](#observer-parameters)):

    coco_suite_t *suite;
    coco_observer_t *observer;

    suite = coco_suite("bbob-biobj", "", "");
    observer = coco_observer("bbob-biobj", "");

    while ((PROBLEM = coco_suite_get_next_problem(suite, observer)) != NULL) {
      size_t dimension = coco_problem_get_dimension(PROBLEM);

      my_optimizer(evaluate_function, 
                   dimension,
                   coco_problem_get_number_of_objectives(PROBLEM),
                   coco_problem_get_smallest_values_of_interest(PROBLEM),
                   coco_problem_get_largest_values_of_interest(PROBLEM),
                   dimension * BUDGET_MULTIPLIER,
                   random_generator);
    }  

    coco_observer_free(observer);
    coco_suite_free(suite);

The ``coco_suite_t`` object is a collection of (in this case biobjective) optimization problems of 
type ``coco_problem_t``. The while loop iterates through all problems of the suite and optimizes 
each of them with ``my_optimizer`` (a simple random search is used in the ``example_experiment``). 
The ``coco_observer_t`` object takes care of logging the performance of the optimizer. The interface
to ``my_optimizer`` includes the following parameters:
- the function that evaluates solutions on the optimization problem in question,
- the number of variables (dimension),
- the number of objectives,
- the smallest and largest values of interest, which define the region of interest in the decision space,
- the maximal budget of evaluations and
- the random generator.

The optimizer should be run until ``dimension * BUDGET_MULTIPLIER`` number of evaluations have 
been reached. In the ``example_experiment``, the ``BUDGET_MULTIPLIER`` is conservatively set using

    static const size_t BUDGET_MULTIPLIER = 2;

so that the experiment runs quickly. The budget needs to be increased for real benchmarking
experiments, but this should be done gradually (it might be sensible to test ``BUDGET_MULTIPLIER = 1e2`` before any larger values are used) to see how it effects the running time of the benchmark. 

The actual ``example_experiment`` contains an additional loop that supports __independent restarts__
by ``my_optimizer`` and takes care of breaking the loop when the target has been hit or the 
budget of function evaluations has been exhausted. While the simple random search used in the 
example does not trigger restarts by itself, a more sophisticated optimizer should (in order to avoid
being stuck in a local optimum). When restarting the algorithm the optimizer should not be 
doing the exactly same thing in every run. 

The ``example_experiment`` records the time needed for optimizing a problem and can therefore 
serve also as a __timing experiment__ for an algorithm. 

Note that the benchmarking procedure remains the same whether we are dealing with single- 
or multi-objective problems and algorithms. To perform benchmarking on a different suite and with a 
different observer, it is enough to replace ``"bbob-biobj"`` with the name of the desired suite and observer. 

In the above example, the suite and observer are called without additional parameters (the empty 
strings ``""`` are used), which means that their default values apply. These can be changed by 
calling:

    suite = coco_suite("bbob-biobj", suite_instance, suite_options);
    observer = coco_observer("bbob-biobj", observer_options);

where ``suite_instance``, ``suite_options`` and ``observer_options`` are strings with parameters 
encoded  as pairs ``"key: value"``. When the value consists of one or more integers, it can be 
encoded using the syntax ``m-n`` (meaning all integer values from m to n), ``-n`` (meaning all 
values up to n), ``n-`` (meaning all values from n on) and even ``-`` (meaning all available 
values); or by simply listing the values separated by commas (as in ``2,3,5``). No spaces are 
allowed in the definition of a range or list of values. 

### Suite parameters <a name="suite-parameters"></a>

The suite contains a collection of problems constructed by a Cartesian product of the suite's 
optimization functions, dimensions and instances. The functions and dimensions are defined by the 
suite name, while the instances are defined with the ``suite_instance`` parameter. The suite can be 
filtered by specifying functions, dimensions and instances through the ``suite_options`` parameter. 

Possible keys and values for ``suite_instance`` are:
- either ``"year: YEAR"``, where ``YEAR`` is usually the year of the corresponding [BBOB 
workshop](http://numbbo.github.io/workshops) defining the instances used in that year's benchmark,
- or ``"instances: VALUES"``, where ``VALUES`` is a list or a range ``m-n`` of instances to be included 
in the suite (starting from 1).

If both ``year`` and ``instances`` appear in the ``suite_instance`` string, only the first one is 
taken into account. If no ``suite_instance`` is given, it defaults to the year of the current BBOB 
workshop. 

Possible keys and values for ``suite_options`` are:
- ``dimensions: LIST``, where ``LIST`` is the list of dimensions to keep in the suite (range-style
syntax is not allowed here), 
- ``dimension_indices: VALUES``, where ``VALUES`` is a list or a range of dimension indices (starting 
from 1) to keep in the suite, and
- ``function_indices: VALUES``, where ``VALUES`` is a list or a range of function indices (starting 
from 1) to keep in the suite, and
- ``instance_indices: VALUES``, where ``VALUES`` is a list or a range of instance indices (starting 
from 1) to keep in the suite. 

If both ``dimensions`` and ``dimension_indices`` appear in the ``suite_options`` string, only the first 
one is taken into account. If no ``suite_options`` is given, no filtering by functions, dimensions and
instances is performed, i.e. the experiment will be run on the entire benchmark suite. 

For example, the call:

    suite = coco_suite("bbob-biobj", 
                       "instances: 10-20", 
                       "dimensions: 2,3,5,10,20 instance_indices:1-5");

first creates the biobjective suite with instances 10 to 20, but then uses only the first five 
dimensions (skipping dimension 40) and the first five instances (i.e. instances 10 to 14) of the suite. 

This kind of filtering can be helpful when parallelizing the benchmark.

See [biobjective test suite](http://numbbo.github.io/coco-doc/bbob-biobj/functions/) and 
[bbob test sute](http://coco.lri.fr/downloads/download15.03/bbobdocfunctions.pdf) for more detailed information on the two 
currently supported suites.

### Observer parameters <a name="observer-parameters"></a>

The observer controls the logging that is performed within the benchmark. Some observer parameters are 
general, while others are specific to the chosen observer. 

Possible keys and values for the general ``observer_options`` are:
- ``result_folder: NAME``, determines the folder within the "exdata" folder into which the results will 
be output. If the folder with the given name already exists, first NAME_001 will be tried, then NAME_002 
and so on. The default value is "default".
- ``algorithm_name: NAME``, where ``NAME`` is a short name of the algorithm that will be used in plots 
(no spaces are allowed). The default value is "ALG".
- ``algorithm_info: STRING`` stores the description of the algorithm. If it contains spaces, it must be 
surrounded by double quotes. The default value is "" (no description).
- ``number_target_triggers: VALUE`` defines the number of targets between each 10**i and 10**(i+1)
(equally spaced in the logarithmic scale) that trigger logging. The default value is 100.
- ``target_precision: VALUE`` defines the precision used for targets (there are no targets for
abs(values) < target_precision). The default value is 1e-8.
- ``number_evaluation_triggers: VALUE`` defines the number of evaluations to be logged between each 10**i
and 10**(i+1). The default value is 20.
- ``base_evaluation_triggers: VALUES`` defines the base evaluations used to produce an additional
evaluation-based logging. The numbers of evaluations that trigger logging are every
base_evaluation * dimension * (10**i). For example, if base_evaluation_triggers = "1,2,5", the logger will
be triggered by evaluations dim*1, dim*2, dim*5, 10*dim*1, 10*dim*2, 10*dim*5, 100*dim*1, 100*dim*2,
100*dim*5, ... The default value is "1,2,5". 
- ``precision_x: VALUE`` defines the precision used when outputting variables and corresponds to the 
number of digits to be printed after the decimal point. The default value is 8.
- ``precision_f: VALUE`` defines the precision used when outputting f values and corresponds to the 
number of digits to be printed after the decimal point. The default value is 15.

Possible keys and values for the ``observer_options`` of the ``bbob-biobj`` observer are:
- ``log_nondominated: STRING`` determines how the nondominated solutions are handled. ``STRING`` can take 
on the values ``none`` (don't log nondominated solutions), ``final`` (log only the final nondominated 
solutions), ``all`` (log every solution that is nondominated at creation time) and ``read`` (the nondominated 
solutions are not logged, but are passed to the logger as input - this is a functionality needed in 
pre-processing of the data). The default value is all.
- ``log_decision_variables: STRING`` determines whether the decision variables are to be logged
in addition to the objective variables in the output of nondominated solutions. ``STRING`` can take 
on the values ``none`` (don't output decision variables), ``low_dim``(output decision variables only 
for dimensions lower or equal to 5) and ``all`` (output all decision variables). The default value is 
log_dim. 
- ``compute_indicators: VALUE`` determines whether to compute and output performance indicators 
(``1``) or not (``0``). The default value is 1.
- ``produce_all_data: VALUE`` determines whether to produce all data required for the workshop. If 
set to ``1``, it overwrites some other options and is equivalent to setting ``log_nondominated`` to 
``all``, ``log_decision_variables`` to ``low_dim`` and ``compute_indicators`` to ``1``. If set to 
``0``, it does not change the values of the other options. The default value is 0.

The benchmark can also be run without any observer, which produces no output, by invoking either 
``""`` or ``"no_observer"`` in place of the observer name. 

### Problem evaluation <a name="problem-evaluation"></a>

In order to evaluate the problem, the following method needs to be invoked:

    void coco_evaluate_function(coco_problem_t *problem, const double *x, double *y);

It will evaluate the problem function in point ``x`` and save the result in ``y``.

In order to evaluate the constraints of the problem, the following method needs to be invoked:

    void coco_evaluate_function(coco_problem_t *problem, const double *x, double *y);

It will evaluate the problem constraints in point ``x`` and save the result in ``y``. Note: while this
functionality is provided, the framework does not yet include problems with constraints.

### Problem properties <a name="problem-properties"></a>

Problem properties can be accessed in the following way:

    /* Returns the number of variables i.e. dimension of the problem */
    size_t coco_problem_get_dimension(const coco_problem_t *problem);

    /* Returns a vector of size 'dimension' with lower bounds of the region of interest in the decision space. */
    const double *coco_problem_get_smallest_values_of_interest(const coco_problem_t *problem);

    /* Returns a vector of size 'dimension' with upper bounds of the region of interest in  the decision space. */
    const double *coco_problem_get_largest_values_of_interest(const coco_problem_t *problem);

    /* Returns the number of objectives of the problem */
    size_t coco_problem_get_number_of_objectives(const coco_problem_t *problem);

    /* Returns the number of evaluations done on the problem */
    size_t coco_problem_get_evaluations(coco_problem_t *problem);

See the ``coco.h`` file for more information on these and other functions that can be used to interface 
COCO problem and other COCO structures. 

## How to write new test problems and combine them into test suites <a name="new-suites"></a>

A test suite is a collection of test problems to be solved during the same benchmarking experiment. Examples of suites in COCO are the single-objective ``bbob`` suite with 24 functions, 6 dimensions and 15 instances (i.e., 2160 problem instances in total) and the biobjective ``bbob-biobj`` suite with 55 functions, 6 dimensions and 10 instances (i.e., 3300 problem instances in total). Note that although the terms *function* and *problem* are sometimes used ambiguously, their meaning should be clear from the context. 

Writing a new test suite entails:

1. Implementing a new set of functions (or reusing existing functions),
2. Defining problem instances, and
3. Collecting problems into a suite.

### Implementing new test functions

Test functions are implemented in the C files starting with ``f_``. Let us use the sphere function from ``f_sphere.c`` to illustrate how this is done. Each function is implemented using three methods:

    /* Returns the square of x */
    static double f_sphere_raw(const double *x, const size_t dimension);

    /* Uses the f_sphere_raw method to compute the function value and store it in y, and the */
    /* problem properties to access other data, for example, the dimension */
    static void f_sphere_evaluate(coco_problem_t *problem, const double *x, double *y);

    /* Creates the sphere problem as a function of the dimension */
    static coco_problem_t *f_sphere_allocate(const size_t dimension);

Implementing a new test function called *blue* would mean defining three new methods ``f_blue_raw``, ``f_blue_evaluate`` and ``f_blue_allocate``. The actual function would be defined in ``f_blue_raw``, while the the other two methods would be very similar to ``f_sphere_evaluate`` and ``f_sphere_allocate`` and would require little effort to implement.

### Defining problem instances

In order to make the optimization problems more challenging, **transformations** such as shifts, oscillations, conditioning and others can be *wrapped around* the basic function or other transformations. For example, the BBOB sphere problem \f[y = \sum_{i=1}^D (x_i-x_i^{\mathrm{opt}})^2 + f^{\mathrm{opt}},\f] was created from the basic sphere function \f[y = \sum_{i=1}^D x_i^2,\f] using two transformations, a shift in the decision space by \f$x^{\mathrm{opt}}\f$ and a shift of the function value by \f$f^{\mathrm{opt}}\f$.    

Transformations take a problem (a ``coco_problem_t`` object often referred to as the inner problem) with some parameters and return the transformed problem (again a ``coco_problem_t`` object). Depending on whether they act on decision variables or the objective value, they are implemented in C files starting with ``f_transform_vars`` or ``f_transform_obj``, respectively. 

For example, the BBOB sphere problem is implemented as:

    problem = f_sphere_allocate(dimension);
    problem = transform_vars_shift(problem, xopt, 0);
    problem = transform_obj_shift(problem, fopt);

Note that transformations of the decision variables first perform the transformation and only then evaluate the inner problem using the new transformed variables, while the transformations of the objective variable first evaluate the inner problem and then transform its output. This is why the BBOB Rastrigin problem  \f[y = 10 \left( D - \sum_{i=1}^D \cos{(2 \pi z_i)} \right) + ||z||^2 + f^{\mathrm{opt}}, \quad \mathrm{where} \quad \mathbf{z} = \Delta^{10}T^{0.2}_{\mathrm{asy}}(T_{\mathrm{osz}}(\mathbf{x} - \mathbf{x}^{\mathrm{opt}}))\f] is defined using the following  order of transformations:

    problem = f_rastrigin_allocate(dimension);
    problem = transform_vars_conditioning(problem, 10.0);
    problem = transform_vars_asymmetric(problem, 0.2);
    problem = transform_vars_oscillate(problem);
    problem = transform_vars_shift(problem, xopt, 0);
    problem = transform_obj_shift(problem, fopt);
 
Varying the values of \f$x^{\mathrm{opt}}\f$ and \f$f^{\mathrm{opt}}\f$ yields **different instances** of the same BBOB problem. Each problem instance is therefore defined by dimension and instance number and implemented in a method such as:

    static coco_problem_t *f_sphere_bbob_problem_allocate(const size_t dimension,
                                                          const size_t instance,
                                                          ...);

The ``f_blue_problem_allocate`` method implementing the *blue* problem would therefore contain a call to ``f_blue_allocate``, possibly some transformations, and finally calls to methods ``coco_problem_set_id``, ``coco_problem_set_name`` and ``coco_problem_set_type`` to set these problem properties. Note that problem allocation methods need to be deterministic (return the same object given the same argument values). 

### Collecting problems into a suite

Once all the required problems are given, they need to be combined into a suite. A suite called *red* would have to implement the following methods stored in the ``suite_red.c`` file:

    static coco_suite_t *suite_red_initialize(void);
    static coco_problem_t *suite_red_get_problem(coco_suite_t *suite,
                                                 const size_t function_idx,
                                                 const size_t dimension_idx,
                                                 const size_t instance_idx);

The initialization is very simple, it requires a call to the suite allocation method, where the number of functions, available dimensions and default instances are set:

    static coco_suite_t *coco_suite_allocate(const char *suite_name,
                                             const size_t number_of_functions,
                                             const size_t number_of_dimensions,
                                             const size_t *dimensions,
                                             const char *default_instances);

The ``suite_red_get_problem`` method has to return the right problem given the suite and function, dimension and instance indices. 

In case the suites instances depend on the year (as is the case with the ``bbob`` and ``bbob-biobj`` suites), a  method that returns a string of instances for the given year can be defined as:

    static const char *suite_red_get_instances_by_year(const int year);

In order for the newly-implemented suite to be included in COCO, the following methods from ``coco_suite.c`` need to be updated (two lines per suite need to be added to each of these methods):

    static coco_suite_t *coco_suite_intialize(const char *suite_name); 
    static coco_problem_t *coco_suite_get_problem_from_indices(coco_suite_t *suite,
                                                               const size_t function_idx,
                                                               const size_t dimension_idx,
                                                               const size_t instance_idx);
    static const char *coco_suite_get_instances_by_year(const coco_suite_t *suite, const int year);

## How to write additional performance indicators and logging functionality <a name="new-indicators"></a>

Here we provide guidelines for implementing a new observer/logger and adding a new performance indicator to the ``bbob-biobj`` logger.  

### Implementing a new observer/logger

First, let us clarify the difference between observers and loggers. An **observer** (a ``coco_observer_t`` object) is a stand-alone entity that can exist independently from problems and test suites. It is defined only by its name and options (see [observer parameters](#observer-parameters)). On the other hand, a **logger** is a COCO problem (a ``coco_problem_t`` object) that is wrapped around the problem to be observed. It exists only for the time span in which its underlying problem exists. The logger needs information from the observer to be able to create log files with some continuity. The observer's task is therefore keeping track of the logging performed on the whole suite, while the actual logging is done by the logger.  

When a new logging functionality is required, both a new observer and a new logger need to be defined. Observers are 'derived' from the ``coco_observer_t`` object, which already contains various information that can be used by an observer. Creating a *yellow* observer means implementing the following structure and methods in the ``observer_yellow.c`` file:

    /* Structure containing data specific to the yellow observer */
    typedef struct {...} observer_yellow_data_t;

    /* Method for freeing the data contained in observer_yellow_data_t (if needed) */
    static void observer_yellow_free(void *data);

    /* Observer constructor that initializes observer_yellow_data_t and connects the observer with the yellow logger */
    /* (set the observer's logger_allocate_function and logger_free_function fields) */
    static void observer_yellow(coco_observer_t *observer, const char *options, coco_option_keys_t **option_keys);

The yellow logger in the ``logger_yellow.c`` file needs to implement the following structure and methods:

    /* Structure containing data specific to the yellow logger */  
    /* (typically pointers to data files etc.) */   
    typedef struct {...} logger_yellow_data_t;

    /* Method for freeing the data contained in logger_yellow_data_t (if needed) */
    static void logger_yellow_free(void *data);

    /* Method that evaluates the inner problem and performs logging */
    static void logger_yellow_evaluate(coco_problem_t *problem, const double *x, double *y);

	/* Logger constructor that initializes logger_yellow_data_t */
    static coco_problem_t *logger_yellow(coco_observer_t *observer, coco_problem_t *inner_problem);

In order to add the observer to the existing observers in COCO, the method ``coco_observer(...)`` in ``coco_observer.c`` needs to be updated (two lines per observer have to be added). Moreover, information on the new observer and its parameters should also be added to this documentation file (see [observer parameters](#observer-parameters)).  

### Adding a new performance indicator to the ``bbob-biobj`` logger

So far, the ``bbob-biobj`` logger contains a single performance indicator - the hypervolume of all nondominated solutions. We now present how other indicators can be added to this logger. 

Indicators of the ``bbob-biobj`` logger are of type ``logger_biobj_indicator_t`` and are stored in the array ``indicators`` in the ``logger_biobj_data_t`` data structure. Currently, ``indicators`` contains a single indicator, but can be easily extended to contain more. For example, to add a *green* indicator, the global counter ``LOGGER_BIOBJ_NUMBER_OF_INDICATORS`` in ``logger_biobj.c`` needs to be increased, the global variable ``logger_biobj_indicators`` needs to be extended with the string ``"green"``, and the array ``suite_biobj_best_values_green`` containing best green indicator values for each problem instance in the ``bbob-biobj`` test suite needs to be created and invoked within the ``suite_biobj_get_best_value(...)`` function in ``suite_biobj.c``. 

Computing/updating the indicator value is done within the following two methods:

    /* Updates the AVL tree containing nondominated solutions */
    static int logger_biobj_tree_update(logger_biobj_data_t *logger,
                                        logger_biobj_avl_item_t *node_item);

    /* Outputs data to the .dat and .tdat files */
    static void logger_biobj_output(logger_biobj_data_t *logger,
                                    const int update_performed,
                                    const logger_biobj_avl_item_t *node_item)

The ``bbob-biobj`` logger keeps an archive of all nondominated solutions in the form of an AVL tree. The ``logger_biobj_tree_update`` method checks for domination of the given solution (``node_item``) and updates the archive and the values of the indicators if the given node is not weakly dominated by existing nodes in the archive. Here the green indicator value needs to be updated any time a solution is added to or removed from the archive.

When the archive changes, the overall indicator values are further updated in the ``logger_biobj_output`` method just before being output. An update of the overall green indicator value is required at this point.

Other functionalities, such as initializing, freeing and outputting to the indicator-specific files should 'work out of the box' without requiring additional tweaking for individual indicators. 