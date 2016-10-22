/*
  Copyright 2002-2004 John Plevyak, All Rights Reserved
*/
#ifndef _d_H_
#define _d_H_

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
#define D_BUILD_VERSION "R-6a201e22e57e77c297d3c43fdb245125ed3b7d64"

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

typedef char int8;
typedef unsigned char uint8;
typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;
typedef short int16;
typedef unsigned short uint16;
#ifdef __MINGW32__
/* already part of most systems */
typedef unsigned long ulong;
typedef uint32 uint; 
#endif

#include "dparse.h"
#include "arg.h"
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
