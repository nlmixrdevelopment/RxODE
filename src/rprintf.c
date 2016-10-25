#include <R.h>
#include <Rinternals.h>

void F77_SUB(rprintf)(char* msg) {
   Rprintf(msg);
   Rprintf("\n");
}

// may be redundant
void F77_SUB(rprintf2)(char* msg) {
   Rprintf(msg);
   Rprintf("\n");
}

void F77_SUB(rprintfid)(char* msg, int *i, double *d) {
   Rprintf(msg, *i, *d);
   Rprintf("\n");
}

void F77_SUB(rprintfdi)(char* msg, double *d, int *i) {
   Rprintf(msg, *d, *i);
   Rprintf("\n");
}

void F77_SUB(rprintfdid)(char* msg, double *d1, int *i, double *d2) {
   Rprintf(msg, *d1, *i, *d2);
   Rprintf("\n");
}

void F77_SUB(rprintfd1)(char* msg, double *d) {
   Rprintf(msg, *d);
   Rprintf("\n");
}

void F77_SUB(rprintfd2)(char* msg, double *d1, double *d2) {
   Rprintf(msg, *d1, *d2);
   Rprintf("\n");
}

void F77_SUB(rprintfi1)(char* msg, int *i) {
   Rprintf(msg, *i);
   Rprintf("\n");
}

void F77_SUB(rprintfi2)(char* msg, int *i1, int *i2) {
   Rprintf(msg, *i1, *i2);
   Rprintf("\n");
}

void F77_SUB(rprintfi3)(char* msg, int *i1, int *i2, int* i3) {
   Rprintf(msg, *i1, *i2, *i3);
   Rprintf("\n");
}
