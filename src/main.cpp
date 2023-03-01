#include <cstdlib>
#include <cstdarg>
#include <resource>

#include <cassert>
#include <vector>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "include/utils.h"
#include "include/enumerations.h"


int finalize (double total_time)
{
  int err = 0;
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  PGFEM_printf("\n");
  PGFEM_printf("Time of analysis - "
               " System %ld.%ld, User %ld.%ld.\n\n",
               usage.ru_stime.tv_sec, usage.ru_stime.tv_usec,
               usage.ru_utime.tv_sec, usage.ru_utime.tv_usec);

  PGFEM_printf("Total time = %f\n\n", total_time);
  // PGFEM_printf("Startup time                  = %f\n", startup_time);
  // PGFEM_printf("Output write time             = %f\n", out_time);
  // PGFEM_printf("Restart write time            = %f\n\n", restart_time);

  return err;
}


int main(int argc, char* argv[])
{    
	int err = 0;
  int process_type;
  double total_time = 0.0;

  // read input commands
  parse_read_options(&process_type,argc,argv);

  try {  
    // run pre-/post-processing routine
    switch (process_type) {

    case PRE_PROCESSING:

      err += pre_1();

      err += pre_2();

      err += pre_3();

      break;

    case POST_PROCESSING:

      err += post_read_inputs();


      err += post_initialze_samples();


      err += post_determine_accept_reject_samples();


      err += post_generate_2d_layer();


      err += post_write_offline_library();

      break;

    default:
      PGFEM_printerr("ERROR, unrecognised type in %s\n",__func__);
      PGFEM_Abort();
      abort();
    }

  total_time += CLOCK(); // measure time spent  

  err += finalize(total_time);
    
} catch (const std::exception& ex) {

  cout << "Something went wrong: " << ex.what() << endl;
  abort();
} catch (...) {} 

	return err;  
} 

