#include <cstdio>
#include "include/data_structure.h"

struct Options {
	
}

struct PreVars {
	int n_sample;
	double opt_R_max;
	int opt_ts_print;
	int opt_par;
	int p_micro_step;	
	// current path : in/sample_inputs (copied from in)
	// opt_R_max // 0 : fixed (default) // 1 : uniform random
	// opt_ts_print // 0 : all (default) // 1 : the first & the last // 2 : where it failed?
	// opt_par // 0 : Human-readble inputs, *_bc.json, *_load.json, timesteps.json
	// 1 : Parallelized inputs, *.bcv, *.load, *_p.in.st, *_t.in.st
}

struct PostVars {
	int np_micro;
	int n_sample;
	int n_micro_step;
	int p_micro_step;
	int n_micro_elem;
	double model_accu;
	double fully_damaged;
}


int pre_processing_initialization(PreVars *prv)
{
  int err = 0;

  prv->n_sample       = 0; 
  prv->opt_R_max      = 0.0;   
  prv->opt_ts_print   = 0;
  prv->opt_par        = 0; 
  prv->p_micro_step   = 0;

  return err;
}

/// destruct time stepping variable
/// free memory spaces for member arrays and structs
///
/// \param[in, out] ts an object for time stepping
/// \return non-zero on internal error
int destruct_pre_processing(PreVars *prv)
{
  int err = 0;

  if(NULL != prv->n_sample)       free(prv->n_sample);
  if(NULL != prv->opt_R_max) 	  free(prv->opt_R_max);
  if(NULL != prv->opt_ts_print)   free(prv->opt_ts_print);
  if(NULL != prv->opt_par)        free(prv->opt_par);
  if(NULL != prv->p_micro_step)   free(prv->p_micro_step);

  err += pre_processing_initialization(prv);

  return err;
}