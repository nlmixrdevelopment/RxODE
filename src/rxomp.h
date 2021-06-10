#ifdef _OPENMP
#include <pthread.h>
#include <omp.h>
#else

static inline int omp_get_num_procs(){
  return 1;
}

static inline int omp_get_thread_limit(){
  return 1;
}

static inline int omp_get_max_threads(){
  return 1;
}

static inline int omp_get_thread_num(){
  return 0;
}
#endif
