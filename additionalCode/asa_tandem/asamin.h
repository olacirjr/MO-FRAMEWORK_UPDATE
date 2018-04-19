/***********************************************************************
* MATLAB Gateway Routine for Lester Ingber's Adaptive Simulated
* Annealing (ASA)
* 
* Copyright (c) 1999-2004  Shinichi Sakata.  All Rights Reserved.
***********************************************************************/
/* $Id: asamin.h,v 1.33 2006/10/23 00:37:57 ssakata Exp ssakata $ */
#define NLHS 5
#define NRHS 28
#define NRHS_ERR_MSG MKSTR2(Error: NRHS arguments are required.)
#define NLHS_ERR_MSG MKSTR2(Error: NLHS variables are required on the righ-hand side.)
#define MAXLEN_ASA_OUT_FILE 80
#define MAXLEN_COST_FUNC_NAME 80
#define MAXLEN_CMD 80

#define SHUFFLE 256		/* size of random array */

static void
set_real_to_option (const mxArray *value,
		     double *var);

static void
set_int_to_option (const mxArray *value,
		     int *var);

static void
set_longint_to_option (const mxArray *value,
		       long int *var);

static void
set_string_to_option (const mxArray *value,
		     char *var);

static double
myrand (LONG_INT * rand_seed);

static double
randflt (LONG_INT *rand_seed);

static double
resettable_randflt (LONG_INT * rand_seed, int reset);

static double
cost_function_without_test (double *x,
			    double *parameter_lower_bound,
			    double *parameter_upper_bound,
			    double *cost_tangents,
			    double *cost_curvature,
			    ALLOC_INT *parameter_dimension,
			    int *parameter_int_real,
			    int *cost_flag,
			    int *exit_code,
			    USER_DEFINES * USER_OPTIONS);

static double
cost_function_with_test (double *x,
			 double *parameter_lower_bound,
			 double *parameter_upper_bound,
			 double *cost_tangents,
			 double *cost_curvature,
			 ALLOC_INT *parameter_dimension,
			 int *parameter_int_real,
			 int *cost_flag,
			 int *exit_code,
			 USER_DEFINES * USER_OPTIONS);

static void
user_acceptance_test (double current_cost,
		      ALLOC_INT * parameter_dimension,
		      USER_DEFINES * USER_OPTIONS);

static void
reset_options (LONG_INT * rand_seed,
	       int * test_in_cost_func,
	       int * use_rejected_cost,
	       USER_DEFINES * USER_OPTIONS);

static void
exit_function(void);
