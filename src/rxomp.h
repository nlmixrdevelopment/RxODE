#ifdef _OPENMP
#include <pthread.h>
#include <omp.h>

static inline int rx_get_thread(int mx) {
  int tn = omp_get_thread_num();
  if (tn < 0) return 0;
  if (tn <= mx) return tn;
  return 0;
}

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

static inline int omp_get_thread_num() {
  return 0;
}

static inline int rx_get_thread(int mx) {
  return 0;
}
 
#endif
