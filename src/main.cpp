#include <cstdlib>
#include <cstdarg>
#include <resource>

#include <cassert>
#include <vector>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "include/utils.cpp"

int finalize(const Multiphysics& mp,
      			 double total_time,
      			 double out_time,
      			 double restart_time,
      			 double startup_time,
      			 double ms_solve_time,
      			 double ms_comm_time,
      			 int myrank)
{
  int err = 0;
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  PGFEM_printf("\n");
  PGFEM_printf("Time of analysis on processor [%d] - "
               " System %ld.%ld, User %ld.%ld.\n\n",
               myrank, usage.ru_stime.tv_sec, usage.ru_stime.tv_usec,
               usage.ru_utime.tv_sec, usage.ru_utime.tv_usec);

  PGFEM_printf("Total time (no Network Init)  = %f\n\n", total_time);
  PGFEM_printf("Startup time                  = %f\n", startup_time);
  PGFEM_printf("Multiscale iteration time     = %f\n", ms_solve_time);
  PGFEM_printf("Multiscale communication time = %f\n", ms_comm_time);
  PGFEM_printf("Output write time             = %f\n", out_time);
  PGFEM_printf("Restart write time            = %f\n\n", restart_time);

  return err;
}


int main()
{
	int err = 0;

  err += read_inputs();

  
  err += initialze_samples();

  
  err += determine_accept_reject_samples();


  err += generate_2d_layer();


  err += write_offline_library();


  err += finalize();  

	return err;
} 

