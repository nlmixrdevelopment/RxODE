/*
  Copyright 2002-2004 John Plevyak, All Rights Reserved
*/
#ifndef _d_H_
#define _d_H_

#define __USE_MINGW_ANSI_STDIO 1
#ifdef MEMWATCH
#define MEMWATCH_STDIO 1
#include "../../src/memwatch-2.67/memwatch.h"
#define MEM_GROW_MACRO
#endif
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#if !defined(__FreeBSD__) || (__FreeBSD_version >= 500000)
#include <inttypes.h>
#endif
#include <limits.h>
#include <sys/types.h>
#ifndef __MINGW32__
#include <sys/mman.h>
#include <sys/uio.h>
#endif
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <dirent.h>
#include <ctype.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#define D_MAJOR_VERSION 1
#define D_MINOR_VERSION 30
#define D_BUILD_VERSION "R-62630b2ba313e6d2f0ae08f452efc8a1e702f731"

#ifdef LEAK_DETECT
#define GC_DEBUG
#include "gc.h"
#define CHECK_LEAKS() GC_gcollect()
#else
#ifdef USE_GC
#include "gc.h"
#define malloc dont_use_malloc_use_MALLOC_instead
#define relloc dont_use_realloc_use_R_chk_realloc_instead
#define free dont_use_free_use_Free_instead
#else
#endif
#endif

// enough already with the signed/unsiged char issues
#define isspace_(_c) isspace((unsigned char)(_c))
#define isdigit_(_c) isdigit((unsigned char)(_c))
#define isxdigit_(_c) isxdigit((unsigned char)(_c))
#define isprint_(_c) isprint((unsigned char)(_c))

#define D_VERSION			(\
(D_MAJOR_VERSION << 24) + (D_MINOR_VERSION << 16) + \
D_BUILD_VERSION)
                         
/* Compilation Options 
*/

#define round2(_x,_n) ((_x + ((_n)-1)) & ~((_n)-1))
#define tohex1(_x) \
((((_x)&15) > 9) ? (((_x)&15) - 10 + 'A') : (((_x)&15) + '0'))
#define tohex2(_x) \
((((_x)>>4) > 9) ? (((_x)>>4) - 10 + 'A') : (((_x)>>4) + '0'))
#define numberof(_x) ((sizeof(_x))/(sizeof((_x)[0])))

typedef int8_t int8;
typedef uint8_t uint8;
typedef int32_t int32;
typedef uint32_t uint32;
typedef int64_t int64;
typedef uint64_t  uint64;
typedef int16_t int16;
typedef uint16_t uint16;
typedef unsigned int uint;

#include "dparse.h"
#include "util.h"
#include "gram.h"
#include "lr.h"
#include "lex.h"
#include "scan.h"
#include "parse.h"
#include "write_tables.h"
#include "read_binary.h"

#ifdef D_DEBUG
#define DBG(_x) if (d_debug_level>1) { _x; }
#else
#define DBG(_x)
#endif

void d_version(char *);

#define USE_SCANNER		1

#endif