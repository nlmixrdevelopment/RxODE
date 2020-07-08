#include <R.h>
#include <Rversion.h>
#include <Rinternals.h>
#include <stdint.h>    // for uint64_t rather than unsigned long long
#include "../inst/include/RxODE.h"
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

// For RxODE sortType = 1
static int *TMP=NULL;               // UINT16_MAX*sizeof(int) for each thread; used by counting sort in radix_r()
static uint8_t *UGRP=NULL;          // 256 bytes for each thread; used by counting sort in radix_r() when sortType==0 (byte appearance order)

// FIXME key per thread
static uint8_t **key = NULL;

static int dround=0; // No rounding by default
static uint64_t dmask=0;
// Original comes from:
// https://github.com/Rdatatable/data.table/blob/588e0725320eacc5d8fc296ee9da4967cee198af/src/forder.c#L396-L414
// for signed integers it's easy: flip sign bit to swap positives and negatives; the resulting unsigned is in the right order with INT_MIN ending up as 0
// for floating point finite you have to flip the other bits too if it was signed: http://stereopsis.com/radix.html
// CHANGES for RxODE:
// By definition in RxODE, the value is always finite, drop the NA, NaN, and +-Inf
uint64_t dtwiddle(const void *p, int i)
{
  union {
    double d;
    uint64_t u64;
  } u;  // local for thread safety
  u.d = ((double *)p)[i];
  if (u.d==0) u.d=0; // changes -0.0 to 0.0,  issue #743
  u.u64 ^= (u.u64 & 0x8000000000000000) ? 0xffffffffffffffff : 0x8000000000000000; // always flip sign bit and if negative (sign bit was set) flip other bits too
  u.u64 += (u.u64 & dmask) << 1/*is this shift really correct. No need to shift*/  ;   // when dround==1|2, if 8th|16th bit is set, round up before chopping last 1|2 bytes
  return u.u64 >> (dround*8);
}

// Note range_d will be calculated in rxData while loading into the parallel R interface instead.




