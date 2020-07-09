#include <R.h>
#include <Rversion.h>
#include <Rinternals.h>
#include <stdint.h>    // for uint64_t rather than unsigned long long
#include <stdbool.h>
#include "../inst/include/RxODE.h"
#define min2( a , b )  ( (a) < (b) ? (a) : (b) )
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// For RxODE sortType = 1
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

// Adapted from:
// https://github.com/Rdatatable/data.table/blob/588e0725320eacc5d8fc296ee9da4967cee198af/src/forder.c#L824-L1258
static bool sort_ugrp(uint8_t *x, const int n)
// x contains n unique bytes; sort them in-place using insert sort
// always ascending. desc and nalast are done in WRITE_KEY because columns may cross byte boundaries
// maximum value for n is 256 which is one too big for uint8_t, hence int n
{
  bool skip = true;            // was x already sorted? if so, the caller can skip reordering
  for (int i=1; i<n; i++) {
    uint8_t tmp = x[i];
    if (tmp>x[i-1]) continue;  // x[i-1]==x[i] doesn't happen because x is unique
    skip = false;
    int j = i-1;
    do {
      x[j+1] = x[j];
    } while (--j>=0 && tmp<x[j]);
    x[j+1] = tmp;
  }
  return skip;
}
// Modified because:
//  - sortType=1 always for RxODE (ascending)
//  - retgrp = 0 (all pushs are meaningless)
//  - Keys are stored based on core
//  - modified so that radix_r is 0 order not 1 order like R (since we are using it in C)
void radix_r(const int from, const int to, const int radix,
	     rx_solving_options_ind *ind, rx_solve *rx) {
  uint8_t **key = rx->keys[omp_get_thread_num()];
  int *anso = ind->ix;
  const int my_n = to-from+1;
  if (my_n==1) {  // minor TODO: batch up the 1's instead in caller (and that's only needed when retgrp anyway)
    return;
  }
  else if (my_n<=256) {
    // if nth==1
    // Rprintf(_("insert clause: radix=%d, my_n=%d, from=%d, to=%d\n"), radix, my_n, from, to);
    // insert sort with some twists:
    // i) detects if grouped; if sortType==0 can then skip
    // ii) keeps group appearance order at byte level to minimize movement
    uint8_t *restrict my_key = key[radix]+from;  // safe to write as we don't use this radix again
    uint8_t o[my_n];
    // if last key (i.e. radix+1==nradix) there are no more keys to reorder so we could reorder osub by reference directly and save allocating and populating o just
    // to use it once. However, o's type is uint8_t so many moves within this max-256 vector should be faster than many moves in osub (4 byte or 8 byte ints) [1 byte
    // type is always aligned]
    bool skip = true;
    int start = 1;
    while (start<my_n && my_key[start]>=my_key[start-1]) start++;
    if (start<my_n) {
      skip = false;  // finding start is really just to take skip out of the loop below
      for (int i=0; i<start; i++) o[i]=i;  // always at least sets o[0]=0
      for (int i=start; i<my_n; i++) {
	uint8_t ktmp = my_key[i];
	int j=i-1;
	while (j>=0 && ktmp<my_key[j]) {
	  my_key[j+1] = my_key[j];
	  o[j+1] = o[j];
	  j--;
	}
	my_key[j+1] = ktmp;  // redundant write when while() did nothing, but that's unlikely given the common case of pre-ordered is handled by skip==true
	o[j+1] = i;          // important to initialize o[] even when while() did nothing.
      }
    }
    if (!skip) {
      // reorder osub and each remaining ksub
      int TMP[my_n];  // on stack fine since my_n is very small (<=256)
      const int *restrict osub = anso+from;
      for (int i=0; i<my_n; i++) TMP[i] = osub[o[i]];
      memcpy((int *restrict)(anso+from), TMP, my_n*sizeof(int));
      for (int r=radix+1; r< rx->nradix; r++) {
        const uint8_t *restrict ksub = key[r]+from;
        for (int i=0; i<my_n; i++) ((uint8_t *)TMP)[i] = ksub[o[i]];
        memcpy((uint8_t *restrict)(key[r]+from), (uint8_t *)TMP, my_n);
      }
    }
    // my_key is now grouped (and sorted by group too if sort!=0)
    // all we have left to do is find the group sizes and either recurse or push
    if (radix+1==rx->nradix) {
      return;
    }
    int ngrp=0, my_gs[my_n];  //minor TODO: could know number of groups with certainty up above
    my_gs[ngrp]=1;
    for (int i=1; i<my_n; i++) {
      if (my_key[i]!=my_key[i-1]) my_gs[++ngrp] = 1;
      else my_gs[ngrp]++;
    }
    ngrp++;
    if (radix+1==rx->nradix || ngrp==my_n) {  // ngrp==my_n => unique groups all size 1 and we can stop recursing now
    } else {
      for (int i=0, f=from; i<ngrp; i++) {
        radix_r(f, f+my_gs[i]-1, radix+1, ind, rx);
        f+=my_gs[i];
      }
    }
    return;
  }
  else if (my_n<=UINT16_MAX) {    // UINT16_MAX==65535 (important not 65536)
    // if (nth==1) Rprintf(_("counting clause: radix=%d, my_n=%d\n"), radix, my_n);
    uint16_t my_counts[256] = {0};  // Needs to be all-0 on entry. This ={0} initialization should be fast as it's on stack. Otherwise, we have to manage
                                    // a stack of counts anyway since this is called recursively and these counts are needed to make the recursive calls.
                                    // This thread-private stack alloc has no chance of false sharing and gives omp and compiler best chance.
    uint8_t *restrict my_ugrp = rx->UGRP + omp_get_thread_num()*256;  // uninitialized is fine; will use the first ngrp items. Only used if sortType==0
    // TODO: ensure my_counts, my_grp and my_tmp below are cache line aligned on both Linux and Windows.
    const uint8_t *restrict my_key = key[radix]+from;
    int ngrp = 0;          // number of groups (items in ugrp[]). Max value 256 but could be uint8_t later perhaps if 0 is understood as 1.
    bool skip = true;      // i) if already _grouped_ and sortType==0 then caller can skip, ii) if already _grouped and sorted__ when sort!=0 then caller can skip too
    for (int i=0; i<my_n; i++) {
      uint8_t elem = my_key[i];
      if (++my_counts[elem]==1) {
	// first time seen this value.  i==0 always does this branch
	my_ugrp[ngrp++]=elem;
      } else if (skip && elem!=my_key[i-1]) {   // does not happen for i==0
	// seen this value before and it isn't the previous value, so data is not grouped
	// including "skip &&" first is to avoid the != comparison
	skip=false;
      }
    }
    if (!skip) {
      // reorder anso and remaining radix keys

      // avoid allocating and populating order vector (my_n long); we just use counts several times to push rather than pull
      // with contiguous-read from osub and ksub, 256 write live cache-lines is worst case. However, often there are many fewer ugrp and only that number of
      // write cache lines will be active. These write-cache lines will be constrained within the UINT16_MAX width, so should be close by in cache, too.
      // If there is a good degree of grouping, there contiguous-read/write both ways happens automatically in this approach.

      // cumulate; for forwards-assign to give cpu prefetch best chance (cpu may not support prefetch backwards).
      uint16_t my_starts[256], my_starts_copy[256];
      // TODO: could be allocated up front (like my_TMP below), or are they better on stack like this? TODO: allocating up front would provide to cache-align them.
      for (int i=0, sum=0; i<ngrp; i++) { uint8_t w=my_ugrp[i]; int tmp=my_counts[w]; my_starts[w]=my_starts_copy[w]=sum; sum+=tmp; }  // cumulate in ugrp appearance order

      int *restrict my_TMP = rx->TMP + omp_get_thread_num()*UINT16_MAX; // Allocated up front to save malloc calls which i) block internally and ii) could fail
      if (radix==0) {
        // anso contains 1:n so skip reading and copying it. Only happens when nrow<65535. Saving worth the branch (untested) when user repeatedly calls a small-n small-cardinality order.
        for (int i=0; i<my_n; i++) anso[my_starts[my_key[i]]++] = i;  // +1 as R is 1-based.
        // The loop counter could be uint_fast16_t since max i here will be UINT16_MAX-1 (65534), hence ++ after last iteration won't overflow 16bits. However, have chosen signed
        // integer for counters for now, as signed probably very slightly faster than unsigned on most platforms from what I can gather.
      }

      // reorder remaining key columns (radix+1 onwards).   This could be done in one-step too (a single pass through x[],  with a larger TMP
      //    that's how its done in the batched approach below.  Which is better?  The way here is multiple (but contiguous) passes through (one-byte) my_key
      if (radix+1<rx->nradix) {
        for (int r=radix+1; r<rx->nradix; r++) {
          memcpy(my_starts, my_starts_copy, 256*sizeof(uint16_t));  // restore starting offsets
          //for (int i=0,last=0; i<256; i++) { int tmp=my_counts[i]; if (tmp==0) continue; my_counts[i]=last; last=tmp; }  // rewind ++'s to offsets
          const uint8_t *restrict ksub = key[r]+from;
          for (int i=0; i<my_n; i++) ((uint8_t *)my_TMP)[my_starts[my_key[i]]++] = ksub[i];
          memcpy(key[r]+from, my_TMP, my_n);
        }
      }
    }

    if (radix+1== rx->nradix) {
      return;  // we're done. avoid allocating and populating very last group sizes for last key
    }
    int my_gs[ngrp==0 ? 256 : ngrp];  // ngrp==0 when sort and skip==true; we didn't count the non-zeros in my_counts yet in that case
    ngrp=0;
    for (int i=0; i<256; i++) if (my_counts[i]) my_gs[ngrp++]=my_counts[i];  // this casts from uint16_t to int32, too
    if (radix+1==rx->nradix) {
      // aside: cannot be all size 1 (a saving used in my_n<=256 case above) because my_n>256 and ngrp<=256
    } else {
      // this single thread will now descend and resolve all groups, now that the groups are close in cache
      for (int i=0, my_from=from; i<ngrp; i++) {
        radix_r(my_from, my_from+my_gs[i]-1, radix+1, ind, rx);
        my_from+=my_gs[i];
      }
    }
    return;
  }
  // else parallel batches. This is called recursively but only once or maybe twice before resolving to UINT16_MAX branch above

  int batchSize = min2(UINT16_MAX, 1+my_n);  // (my_n-1)/nBatch + 1;   //UINT16_MAX == 65535
  int nBatch = (my_n-1)/batchSize + 1;   // TODO: make nBatch a multiple of nThreads?
  int lastBatchSize = my_n - (nBatch-1)*batchSize;
  uint16_t *counts = calloc(nBatch*256,sizeof(uint16_t));
  uint8_t  *ugrps =  malloc(nBatch*256*sizeof(uint8_t));
  int      *ngrps =  calloc(nBatch    ,sizeof(int));
  /* if (!counts || !ugrps || !ngrps) STOP(_("Failed to allocate parallel counts. my_n=%d, nBatch=%d"), my_n, nBatch); */

  bool skip=true;
  const int n_rem = rx->nradix-radix-1;   // how many radix are remaining after this one
  {
    int     *my_otmp = malloc(batchSize * sizeof(int)); // thread-private write
    uint8_t *my_ktmp = malloc(batchSize * sizeof(uint8_t) * n_rem);
    // TODO: move these up above and point restrict[me] to them. Easier to Error that way if failed to alloc.
    for (int batch=0; batch<nBatch; batch++) {
      const int my_n = (batch==nBatch-1) ? lastBatchSize : batchSize;  // lastBatchSize == batchSize when my_n is a multiple of batchSize
      const int my_from = from + batch*batchSize;
      uint16_t *restrict      my_counts = counts + batch*256;
      uint8_t  *restrict      my_ugrp   = ugrps  + batch*256;
      int                     my_ngrp   = 0;
      bool                    my_skip   = true;
      const uint8_t *restrict my_key    = key[radix] + my_from;
      const uint8_t *restrict byte = my_key;
      for (int i=0; i<my_n; i++, byte++) {
        if (++my_counts[*byte]==1) {   // always true first time when i==0
          my_ugrp[my_ngrp++] = *byte;
        } else if (my_skip && byte[0]!=byte[-1]) {   // include 'my_skip &&' to save != comparison after it's realized this batch is not grouped
          my_skip=false;
        }
      }
      ngrps[batch] = my_ngrp;  // write once to this shared cache line
      if (!my_skip) {
        skip = false;          // naked write to this shared byte is ok because false is only value written
        // gather this batch's anso and remaining keys. If we sorting too, urgrp is sorted later for that. Here we want to benefit from skip within batch
        // as much as possible which is a good chance since batchSize is relatively small (65535)
        for (int i=0, sum=0; i<my_ngrp; i++) { int tmp = my_counts[my_ugrp[i]]; my_counts[my_ugrp[i]]=sum; sum+=tmp; } // cumulate counts of this batch
        const int *restrict osub = anso+my_from;
        byte = my_key;
        for (int i=0; i<my_n; i++, byte++) {
          int dest = my_counts[*byte]++;
          my_otmp[dest] = *osub++;  // wastefully copies out 1:n when radix==0, but do not optimize as unlikely worth code complexity. my_otmp is not large, for example. Use first TEND() to decide.
          for (int r=0; r<n_rem; r++) my_ktmp[r*my_n + dest] = key[radix+1+r][my_from+i];   // reorder remaining keys
        }
        // or could do multiple passes through my_key like in the my_n<=65535 approach above. Test which is better depending on if TEND() points here.

        // we haven't completed all batches, so we don't know where these groups should place yet
        // So for now we write the thread-private small now-grouped buffers back in-place. The counts and groups across all batches will be used below to move these blocks.
        memcpy(anso+my_from, my_otmp, my_n*sizeof(int));
        for (int r=0; r<n_rem; r++) memcpy(key[radix+1+r]+my_from, my_ktmp+r*my_n, my_n*sizeof(uint8_t));

        // revert cumulate back to counts ready for vertical cumulate
        for (int i=0, last=0; i<my_ngrp; i++) { int tmp = my_counts[my_ugrp[i]]; my_counts[my_ugrp[i]]-=last; last=tmp; }
      }
    }
    free(my_otmp);
    free(my_ktmp);
  }

  // If my_n input is grouped and ugrp is sorted too (to illustrate), status now would be :
  // counts:                 ugrps:   ngrps:
  // 1: 20 18  2  0  0  0    0 1 2    3
  // 2:  0  0 17 21  5  0    2 3 4    3
  // 3:  0  0  0  0 15 19    4 5      2
  // If the keys within each and every batch were grouped, skip will be true.
  // Now we test if groups occurred in order across batches (like illustration above), and if not set skip=false

  uint8_t ugrp[256];  // head(ugrp,ngrp) will contain the unique values in appearance order
  bool    seen[256];  // is the value present in ugrp already
  int     ngrp=0;     // max value 256 so not uint8_t
  uint8_t last_seen=0;  // the last grp seen in the previous batch.  initialized 0 is not used
  for (int i=0; i<256; i++) seen[i]=false;
  for (int batch=0; batch<nBatch; batch++) {
    const uint8_t *restrict my_ugrp = ugrps + batch*256;
    if (ngrp==256 && !skip) break;  // no need to carry on
    for (int i=0; i<ngrps[batch]; i++) {
      if (!seen[my_ugrp[i]]) {
        seen[my_ugrp[i]] = true;
        ugrp[ngrp++] = last_seen = my_ugrp[i];
      } else if (skip && my_ugrp[i]!=last_seen) {   // ==last_seen would occur accross batch boundaries, like 2=>17 and 5=>15 in illustration above
        skip=false;
      }
    }
  }

  // If skip==true (my_n was pre-grouped) and
  // i) sortType==0 (not sorting groups) then osub and ksub don't need reordering. ugrp may happen to be sorted too but nothing special to do in that case.
  // ii) sortType==1|-1 and ugrp is already sorted then osub and ksub don't need reordering either.

  if (!sort_ugrp(ugrp, ngrp))
    skip=false;

  // now cumulate counts vertically to see where the blocks in the batches should be placed in the result across all batches
  // the counts are uint16_t but the cumulate needs to be int32_t (or int64_t in future) to hold the offsets
  // If skip==true and we're already done, we still need the first row of this cummulate (diff to get total group sizes) to push() or recurse below

  int *starts = calloc(nBatch*256, sizeof(int));  // keep starts the same shape and ugrp order as counts
  for (int j=0, sum=0; j<ngrp; j++) {  // iterate through columns (ngrp bytes)
    uint16_t *tmp1 = counts+ugrp[j];
    int      *tmp2 = starts+ugrp[j];
    for (int batch=0; batch<nBatch; batch++) {
      *tmp2 = sum;
      tmp2 += 256;
      sum += *tmp1;
      tmp1 += 256;
    }
  }
  // the first row now (when diff'd) now contains the size of each group across all batches

  if (!skip) {
    int *TMP = malloc(my_n * sizeof(int));
    /* if (!TMP) STOP(_("Unable to allocate TMP for my_n=%d items in parallel batch counting"), my_n); */
    for (int batch=0; batch<nBatch; batch++) {
      const int *restrict      my_starts = starts + batch*256;
      const uint16_t *restrict my_counts = counts + batch*256;
      const int *restrict      osub = anso + from + batch*batchSize;  // the groups sit here contiguously
      const uint8_t *restrict  byte = ugrps + batch*256;              // in appearance order always logged here in ugrps
      const int                my_ngrp = ngrps[batch];
      for (int i=0; i<my_ngrp; i++, byte++) {
        const uint16_t len = my_counts[*byte];
        memcpy(TMP+my_starts[*byte], osub, len*sizeof(int));
        osub += len;
      }
    }
    memcpy(anso+from, TMP, my_n*sizeof(int));

    for (int r=0; r<n_rem; r++) {    // TODO: groups of sizeof(anso)  4 byte int currently  (in future 8).  To save team startup cost (but unlikely significant anyway)
      for (int batch=0; batch<nBatch; batch++) {
        const int *restrict      my_starts = starts + batch*256;
        const uint16_t *restrict my_counts = counts + batch*256;
        const uint8_t *restrict  ksub = key[radix+1+r] + from + batch*batchSize;  // the groups sit here contiguosly
        const uint8_t *restrict  byte = ugrps + batch*256;                        // in appearance order always logged here in ugrps
        const int                my_ngrp = ngrps[batch];
        for (int i=0; i<my_ngrp; i++, byte++) {
          const uint16_t len = my_counts[*byte];
          memcpy((uint8_t *)TMP + my_starts[*byte], ksub, len);
          ksub += len;
        }
      }
      memcpy(key[radix+1+r]+from, (uint8_t *)TMP, my_n);
    }
    free(TMP);
  }

  int my_gs[ngrp];
  for (int i=1; i<ngrp; i++) my_gs[i-1] = starts[ugrp[i]] - starts[ugrp[i-1]];   // use the first row of starts to get totals
  my_gs[ngrp-1] = my_n - starts[ugrp[ngrp-1]];

  if (radix+1==rx->nradix) {
    // aside: ngrp==my_n (all size 1 groups) isn't a possible short-circuit here similar to my_n>256 case above, my_n>65535 but ngrp<=256
  }
  else {
    // TODO: explicitly repeat parallel batch for any skew bins
    for (int i=0; i<ngrp; i++) {
      int start = from + starts[ugrp[i]];
      radix_r(start, start+my_gs[i]-1, radix+1, ind, rx);
    }
  }
  free(counts);
  free(starts);
  free(ugrps);
  free(ngrps);
}
