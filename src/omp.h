#ifdef _OPENMP
#include <pthread.h>
#include <omp.h>
#else

int omp_get_num_procs(){
  return 1;
}

int omp_get_thread_limit(){
  return 1;
}

int omp_get_max_threads(){
  return 1;
}

int omp_get_thread_num(){
  return 0;
}
#endif
