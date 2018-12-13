#include <sys/stat.h> 
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>   /* dj: import intptr_t */
#include "ode.h"
#include <dparser.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include "tran.g.d_parser.c"
#define max(a,b) (a)>(b) ? (a):(b)
#define min(a,b) (a)<(b) ? (a):(b)
#define MXSYM 50000
#define MXDER 5000
#define MXLEN 12000
#define MXBUF 48000
#define SBPTR sb.s+sb.o
#define SBTPTR sbt.s+sbt.o

#define NOASSIGN "'<-' not supported, use '=' instead or set 'options(RxODE.syntax.assign = TRUE)'."
#define NEEDSEMI "Lines need to end with ';' or to match R's handling of line endings set 'options(RxODE.syntax.require.semicolon = FALSE)'."
#define NEEDPOW "'**' not supported, use '^' instead or set 'options(RxODE.syntax.star.pow = TRUE)'."
#define NODOT "'.' in variables and states not supported, use '_' instead or set 'options(RxODE.syntax.allow.dots = TRUE)'."
#define NOINI0 "'%s(0)' for initialization not allowed.  To allow set 'options(RxODE.syntax.allow.ini0 = TRUE)'."
#define NOSTATE "Defined 'df(%s)/dy(%s)', but '%s' is not a state!"
#define NOSTATEVAR "Defined 'df(%s)/dy(%s)', but '%s' is not a state or variable!"

#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#if (__STDC_VERSION__ >= 199901L)
#include <stdint.h>
#endif

char *repl_str(const char *str, const char *from, const char *to) {
  // From http://creativeandcritical.net/str-replace-c by Laird Shaw
  /* Adjust each of the below values to suit your needs. */

  /* Increment positions cache size initially by this number. */
  size_t cache_sz_inc = 16;
  /* Thereafter, each time capacity needs to be increased,
   * multiply the increment by this factor. */
  const size_t cache_sz_inc_factor = 3;
  /* But never increment capacity by more than this number. */
  const size_t cache_sz_inc_max = 1048576;

  char *pret, *ret = NULL;
  const char *pstr2, *pstr = str;
  size_t i, count = 0;
#if (__STDC_VERSION__ >= 199901L)
  uintptr_t *pos_cache_tmp, *pos_cache = NULL;
#else
  ptrdiff_t *pos_cache_tmp, *pos_cache = NULL;
#endif
  size_t cache_sz = 0;
  size_t cpylen, orglen, retlen, tolen, fromlen = strlen(from);

  /* Find all matches and cache their positions. */
  while ((pstr2 = strstr(pstr, from)) != NULL) {
    count++;

    /* Increase the cache size when necessary. */
    if (cache_sz < count) {
      cache_sz += cache_sz_inc;
      pos_cache_tmp = R_chk_realloc(pos_cache, sizeof(*pos_cache) * cache_sz);
      if (pos_cache_tmp == NULL) {
        goto end_repl_str;
      } else pos_cache = pos_cache_tmp;
      cache_sz_inc *= cache_sz_inc_factor;
      if (cache_sz_inc > cache_sz_inc_max) {
        cache_sz_inc = cache_sz_inc_max;
      }
    }

    pos_cache[count-1] = pstr2 - str;
    pstr = pstr2 + fromlen;
  }

  orglen = pstr - str + strlen(pstr);

  /* Allocate memory for the post-replacement string. */
  if (count > 0) {
    tolen = strlen(to);
    retlen = orglen + (tolen - fromlen) * count;
  } else        retlen = orglen;
  ret = R_chk_calloc(1,retlen + 1);
  if (ret == NULL) {
    goto end_repl_str;
  }

  if (count == 0) {
    /* If no matches, then just duplicate the string. */
    strcpy(ret, str);
  } else {
    /* Otherwise, duplicate the string whilst performing
     * the replacements using the position cache. */
    pret = ret;
    memcpy(pret, str, pos_cache[0]);
    pret += pos_cache[0];
    for (i = 0; i < count; i++) {
      memcpy(pret, to, tolen);
      pret += tolen;
      pstr = str + pos_cache[i] + fromlen;
      cpylen = (i == count-1 ? orglen : pos_cache[i+1]) - pos_cache[i] - fromlen;
      memcpy(pret, pstr, cpylen);
      pret += cpylen;
    }
    ret[retlen] = '\0';
  }

 end_repl_str:
  /* Free the cache and return the post-replacement string,
   * which will be NULL in the event of an error. */
  Free(pos_cache);
  return ret;
}

// from mkdparse_tree.h
typedef void (print_node_fn_t)(int depth, char *token_name, char *token_value, void *client_data);

int R_get_option(const char *option, int def){
  SEXP s, t;
  int ret, pro=0;
  PROTECT(t = s = allocList(3));pro++;
  SET_TYPEOF(s, LANGSXP);
  SETCAR(t, install("getOption")); t = CDR(t);
  SETCAR(t, mkString(option)); t = CDR(t);
  if (def){
    SETCAR(t, ScalarLogical(1));
  } else {
    SETCAR(t, ScalarLogical(0));
  }
  ret = INTEGER(eval(s,R_GlobalEnv))[0];
  UNPROTECT(pro);
  return ret;
}

// Taken from dparser and changed to use R_alloc
int r_buf_read(const char *pathname, char **buf, int *len) {
  struct stat sb;
  int fd;

  *buf = 0;
  *len = 0;
  fd = open(pathname, O_RDONLY);
  if (fd <= 0) 
    return -1;
  memset(&sb, 0, sizeof(sb));
  fstat(fd, &sb);
  *len = sb.st_size;
  *buf = (char*)R_alloc(*len + 2,sizeof(char));
  // MINGW likes to convert cr lf => lf which messes with the size
  size_t real_size = read(fd, *buf, *len);
  (*buf)[real_size] = 0;
  (*buf)[real_size + 1] = 0;
  *len = real_size;
  close(fd);
  return *len;
}

// Taken from dparser and changed to use R_alloc
char * r_sbuf_read(const char *pathname) {
  char *buf;
  int len;
  if (r_buf_read(pathname, &buf, &len) < 0)
    return NULL;
  return buf;
}


// Taken from dparser and changed to use Calloc
char * rc_dup_str(const char *s, const char *e) {
  int l = e ? e-s : (int)strlen(s);
  char *ss = Calloc(l+1,char);
  memcpy(ss, s, l);
  ss[l] = 0;
  return ss;
}

// Taken from dparser and changed to use R_alloc
char * r_dup_str(const char *s, const char *e) {
  int l = e ? e-s : (int)strlen(s);
  char *ss = (char*)R_alloc(l+1,sizeof(char));
  memcpy(ss, s, l);
  ss[l] = 0;
  return ss;
}

int rx_syntax_error = 0, rx_suppress_syntax_info=0, rx_podo = 0;
static void trans_syntax_error_report_fn(char *err) {
  if (!rx_suppress_syntax_info)
    Rprintf("%s\n",err);
  rx_syntax_error = 1;
}


extern D_ParserTables parser_tables_RxODE;

unsigned int found_jac = 0, found_print = 0;
int rx_syntax_assign = 0, rx_syntax_star_pow = 0,
  rx_syntax_require_semicolon = 0, rx_syntax_allow_dots = 0,
  rx_syntax_allow_ini0 = 1, rx_syntax_allow_ini = 1, rx_syntax_allow_assign_state = 0,
  maxSumProdN = 0, SumProdLD = 0;

char s_aux_info[64*MXSYM];


typedef struct symtab {
  char ss[64*MXSYM];                     /* symbol string: all vars*/
  char de[64*MXSYM];             /* symbol string: all Des*/
  char ddt[MXSYM];
  int deo[MXSYM];        /* offest of des */
  int vo[MXSYM];        /* offset of symbols */
  int lh[MXSYM];        /* lhs symbols? =9 if a state var*/
  int ini[MXSYM];        /* initial variable assignment =2 if there are two assignments */
  int ini0[MXSYM];        /* state initial variable assignment =2 if there are two assignments */
  int di[MXDER];        /* ith of state vars */
  int idi[MXDER];       /* should ith state variable be ignored 0/1 */
  int fdi[MXDER];        /* Functional initialization of state variable */
  int nv;                       /* nbr of symbols */
  int ix;                       /* ith of curr symbol */
  int id;                       /* ith of curr symbol */
  int fn;                       /* curr symbol a fn?*/
  int nd;                       /* nbr of dydt */
  int pos;
  int pos_de;
  int ini_i; // #ini
  int statei; // # states
  int fdn; // # conditional states
  int sensi;
  int li; // # lhs
  int pi; // # param
  // Save Jacobian information
  int df[MXSYM];
  int dy[MXSYM];
  int sdfdy[MXSYM];
  int cdf;
  int ndfdy;
  int maxtheta;
  int maxeta;
} symtab;
symtab tb;

typedef struct sbuf {
  char *s;        /* curr print buffer */
  int sN;
  int o;                        /* offset of print buffer */
} sbuf;

sbuf sb;                        /* buffer w/ current parsed & translated line */
                                /* to be stored in a temp file */
sbuf sbt;

char *sgets(char * str, int num, sbuf *sbb){
  if (sbb->sN == 0) return NULL;
  // Adapted from
  // https://stackoverflow.com/questions/2068975/can-cs-fgets-be-coaxed-to-work-with-a-string-not-from-a-file
  // By Mark Williams
  // To work with sbufs
  char *next = sbb->s + sbb->o;
  int  numread = 0;

  while ( numread + 1 < num && *next ) {
    int isnewline = ( *next == '\n' );
    *str++ = *next++;
    numread++;
    // newline terminates the line but is included
    if ( isnewline )
      break;
  }

  if ( numread == 0 )
    return NULL;  // "eof"

  // must have hit the null terminator or*end of line
  *str = '\0';  // null terminate this string
  // set up input for next call
  sbb->o += numread;
  return str;
}

void sIni(sbuf *sbb){
  if (sbb->sN <= 0){
    sbb->s = Calloc(MXBUF, char);
    sbb->sN = MXBUF;
    sbb->s[0]='\0';
    sbb->o=0;
  } else {
    sbb->s[0]='\0';
    sbb->o=0;
  }
}

void sAppendN(sbuf *sbb, const char *what, int n){
  if (sbb->sN <= n + sbb->o){
    int mx = sbb->o + n + MXBUF;
    sbb->s = Realloc(sbb->s, mx, char);
    sbb->sN = mx;
  }
  sprintf(sbb->s+sbb->o, "%s", what);
  sbb->o +=n;
}

static void sPut(sbuf *sbb, char what){
  if (sbb->sN <= 1 + sbb->o){
    int mx = sbb->o + 1 + MXBUF;
    sbb->s = Realloc(sbb->s, mx, char);
    sbb->sN = mx;
  }
  sprintf(sbb->s+sbb->o, "%c", what);
  sbb->o++;
}

void sAppend(sbuf *sbb, const char *format, ...){
  char what[MXBUF*2];
  int n = 0;
  va_list argptr, copy;
  va_start(argptr, format);
  va_copy(copy, argptr);
  // Try first.
  n = vsnprintf(what, MXBUF*2, format, argptr);
  va_end(argptr);
  char *what2;
  int use2=0;
  if (n >= MXBUF*2){
    // Its too big;  Allocate it.
    what2 = Calloc(n+1, char);
    vsnprintf(what2, n, format, copy);
    use2=1;
  }
  va_end(copy);
  
  if (sbb->sN <= n + 1 + sbb->o){
    int mx = sbb->o + n + 1 + MXBUF;
    sbb->s = Realloc(sbb->s, mx, char);
    sbb->sN = mx;
  }
  if (use2){
    sprintf(sbb->s+sbb->o, "%s", what2);
    Free(what2);
  } else {
    sprintf(sbb->s+sbb->o, "%s", what);
  }
  sbb->o +=n;
}

sbuf sbPm, sbNrm;
char *extra_buf, *model_prefix, *md5;

sbuf sbOut[5];

static FILE *fpIO;

/* new symbol? if no, find it's ith */
int new_or_ith(const char *s) {
  int i, len, len_s=(int)strlen(s);

  if (tb.fn) return 0;
  if (!strcmp("t", s)) return 0;
  if (!strcmp("time", s)) return 0;
  if (!strcmp("podo", s)) return 0;
  if (!strcmp("rx__PTR__", s)) return 0;
  if (!strcmp("tlast", s)) return 0;
  // Ignore M_ constants
  if (!strcmp("M_E", s)) return 0;
  if (!strcmp("M_LOG2E", s)) return 0;
  if (!strcmp("M_LOG10E", s)) return 0;
  if (!strcmp("M_LN2", s)) return 0;
  if (!strcmp("M_LN10", s)) return 0;
  if (!strcmp("M_PI", s)) return 0;
  if (!strcmp("M_PI_2", s)) return 0;
  if (!strcmp("M_PI_4", s)) return 0;
  if (!strcmp("M_1_PI", s)) return 0;
  if (!strcmp("M_2_PI", s)) return 0;
  if (!strcmp("M_2_SQRTPI", s)) return 0;
  if (!strcmp("M_SQRT2", s)) return 0;
  if (!strcmp("M_SQRT1_2", s)) return 0;
  if (!strcmp("M_SQRT_3", s)) return 0;
  if (!strcmp("M_SQRT_32", s)) return 0;
  if (!strcmp("M_LOG10_2", s)) return 0;
  if (!strcmp("M_2PI", s)) return 0;
  if (!strcmp("M_SQRT_PI", s)) return 0;
  if (!strcmp("M_1_SQRT_2PI", s)) return 0;
  if (!strcmp("M_SQRT_2dPI", s)) return 0;
  if (!strcmp("M_LN_SQRT_PI", s)) return 0;
  if (!strcmp("M_LN_SQRT_2PI", s)) return 0;
  if (!strcmp("M_LN_SQRT_PId2", s)) return 0;
  // Ignore THETA[] and ETA
  if (strstr("[", s) != NULL) return 0;
  if (!tb.nv) return 1;

  for (i=0; i<tb.nv; i++) {
    len = tb.vo[i+1] - tb.vo[i] - 1;  /* -1 for added ',' */
    if (!strncmp(tb.ss+tb.vo[i], s, max(len, len_s))) { /* note we need take the max in order not to match a sub-string */
      tb.ix = i;
      return 0;
    }
  }
  return 1;
}

int new_de(const char *s){
  int i, len, len_s=(int)strlen(s);
  for (i=0; i<tb.nd; i++) {
    len = tb.deo[i+1] - tb.deo[i] - 1;
    if (!strncmp(tb.de+tb.deo[i], s, max(len, len_s))) { /* note we need take the max in order not to match a sub-string */
      tb.id = i;
      return 0;
    }
  }
  return 1;
}

void wprint_node(int depth, char *name, char *value, void *client_data) {
  int i;
  if (!strcmp("time",value)){
    sAppendN(&sb, "t", 1);
    sAppendN(&sbt, "t", 1);
  } else if (!strcmp("podo",value)){
    sAppendN(&sb, "_solveData->subjects[_cSub].podo", 32);
    sAppendN(&sbt, "podo", 4);
    rx_podo = 1;
  } else if (!strcmp("tlast",value)){
    sAppendN(&sb, "_solveData->subjects[_cSub].tlast", 33);
    sAppendN(&sbt, "tlast", 5);
  } else if (!strcmp("rx__PTR__",value)){
    sAppendN(&sb, "_solveData, _cSub", 17);
    sAppendN(&sbt, "rx__PTR__", 9);
  } else if (!strcmp("identifier",name) && !strcmp("gamma",value)){
    sAppendN(&sb, "lgammafn", 8);
    sAppendN(&sbt, "lgammafn", 8);
  } else if (!strcmp("identifier",name) && !strcmp("lfactorial",value)){
    sAppendN(&sb, "lgamma1p", 8);
    sAppendN(&sbt, "lgamma1p", 8);
  } else if (!strcmp("identifier",name) && !strcmp("log",value)){
    sAppendN(&sb, "_safe_log", 9);
    sAppendN(&sbt, "log", 3);
  } else if (!strcmp("identifier",name) && !strcmp("abs",value)){
    sAppendN(&sb, "fabs", 4);
    sAppendN(&sbt,"abs", 3);
  } else {
    // Apply fix for dot.syntax
    for (i = 0; i < (int)strlen(value); i++){
      if (value[i] == '.' && !strcmp("identifier_r",name)){
	sAppendN(&sb, "_DoT_", 5);
	sAppendN(&sbt, ".", 1);
        if (!rx_syntax_allow_dots){
          trans_syntax_error_report_fn(NODOT);
        }
      } else {
	sPut(&sb, value[i]);
	sPut(&sbt, value[i]);
      }
    }
  }
}

void wprint_parsetree(D_ParserTables pt, D_ParseNode *pn, int depth, print_node_fn_t fn, void *client_data) {
  char *name = (char*)pt.symbols[pn->symbol].name;
  int nch = d_get_number_of_children(pn), i, k, ii, found, safe_zero = 0;
  char *value = (char*)rc_dup_str(pn->start_loc.s, pn->end);
  char buf[1024];
  if ((!strcmp("identifier", name) || !strcmp("identifier_r", name) ||
       !strcmp("identifier_r_no_output",name)  ||
       !strcmp("theta0_noout", name) || 
       !strcmp("theta0", name)) &&
      new_or_ith(value)) {
    /* printf("[%d]->%s\n",tb.nv,value); */
    sprintf(tb.ss+tb.pos, "%s,", value);
    tb.pos += (int)strlen(value)+1;
    // Ignored variables
    if (!strcmp("rx_lambda_", value) || !strcmp("rx_yj_", value)){
      tb.lh[tb.nv] = 11; // Suppress param printout.
    }
    tb.vo[++tb.nv] = tb.pos;
    
  }
  if (!strcmp("(", name) ||
      !strcmp(")", name) ||
      !strcmp(",", name)
      ) {
    sPut(&sb, name[0]);
    if (!(strcmp(",", name)) && depth == 1){
      sAppendN(&sb, "(double)", 8);
    }
    sPut(&sbt, name[0]);
  }
  if (!strcmp("identifier", name) ||
      !strcmp("identifier_r", name) ||
      !strcmp("constant", name) ||
      !strcmp("theta0", name) ||
      !strcmp("+", name) ||
      !strcmp("-", name) ||
      !strcmp("*", name) ||
      !strcmp("/", name) ||

      !strcmp("&&", name) ||
      !strcmp("||", name) ||
      !strcmp("!=", name) ||
      !strcmp("==", name) ||
      !strcmp("<=", name) ||
      !strcmp(">=", name) ||
      !strcmp("!", name) ||
      !strcmp("<", name) ||
      !strcmp(">", name) ||

      !strcmp("=", name)
      )
    fn(depth, name, value, client_data);

  // Operator synonyms  
  if (!strcmp("<-",name)){
    sAppendN(&sb, " =", 2);
    sAppendN(&sbt, "=", 1);
  }

  // Suppress LHS calculation with ~
  if (!strcmp("~",name)){
    sAppendN(&sb, " =", 2);
    sAppendN(&sbt, "~", 1);
    tb.lh[tb.ix] = 10; // Suppress LHS printout.
  }
  
  if (!strcmp("|",name)){
    sAppendN(&sb, " ||", 3);
    sAppendN(&sbt, "||", 2);
  }

  if (!strcmp("&",name)){
    sAppendN(&sb, " &&", 3);
    sAppendN(&sbt, "&&", 2);
  }

  if (!strcmp("<>",name) ||
      !strcmp("~=",name) ||
      !strcmp("/=",name) 
      ){
    sAppendN(&sb, " !=", 3);
    sAppendN(&sbt, "!=", 2);
  }
  Free(value);

  //depth++;
  if (nch != 0) {
    if (!strcmp("power_expression", name)) {
      sAppendN(&sb, " R_pow(", 7);
    }
    for (i = 0; i < nch; i++) {
      if (!rx_syntax_assign  &&
          ((i == 4 && !strcmp("derivative", name)) ||
           (i == 6 && !strcmp("jac", name)) ||
           (i == 6 && !strcmp("dfdy", name)))) {
        D_ParseNode *xpn = d_get_child(pn,i);
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (!strcmp("<-",v)){
          trans_syntax_error_report_fn(NOASSIGN);
        }
        Free(v);
        continue;
      }
      if (!strcmp("derivative", name) && i< 2) continue;
      if (!strcmp("der_rhs", name)    && i< 2) continue;
      if (!strcmp("inf_rhs", name)    && i< 2) continue;
      if (!strcmp("derivative", name) && i==3) continue;
      if (!strcmp("der_rhs", name)    && i==3) continue;
      if (!strcmp("inf_rhs", name)    && i==3) continue;
      if (!strcmp("derivative", name) && i==4){
	D_ParseNode *xpn = d_get_child(pn,i);
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (!strcmp("~",v)){
        }
        Free(v);
	continue;
      }
      
      if (!strcmp("jac", name)     && i< 2)   continue;
      if (!strcmp("jac_rhs", name) && i< 2)   continue;
      if (!strcmp("jac", name)     && i == 3) continue;
      if (!strcmp("jac_rhs", name) && i == 3) continue;
      if (!strcmp("jac", name)     && i == 5) continue;
      if (!strcmp("jac_rhs", name) && i == 5) continue;
      if (!strcmp("jac", name)     && i == 6) continue;

      if (!strcmp("dfdy", name)     && i< 2)   continue;
      if (!strcmp("dfdy_rhs", name) && i< 2)   continue;
      if (!strcmp("dfdy", name)     && i == 3) continue;
      if (!strcmp("dfdy_rhs", name) && i == 3) continue;
      if (!strcmp("dfdy", name)     && i == 5) continue;
      if (!strcmp("dfdy_rhs", name) && i == 5) continue;
      if (!strcmp("dfdy", name)     && i == 6) continue;
      if (!strcmp("ini0", name)     && i == 1) continue;

      if (!strcmp("transit2", name) && i == 1) continue;
      if (!strcmp("transit3", name) && i == 1) continue;

      if (!strcmp("lfactorial",name) && i != 1) continue;
      if (!strcmp("factorial",name) && i != 0) continue;

      if ((!strcmp("theta",name) || !strcmp("eta",name)) && i != 2) continue;
      
      tb.fn = (!strcmp("function", name) && i==0) ? 1 : 0;

      if (tb.fn) depth = 0;

      D_ParseNode *xpn = d_get_child(pn,i);
      
      if (tb.fn){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (!strcmp("prod",v) || !strcmp("sum",v) || !strcmp("sign",v) ||
	    !strcmp("max",v) || !strcmp("min",v)){
	  ii = d_get_number_of_children(d_get_child(pn,3))+1;
	  if (!strcmp("prod", v)){
            sAppend(&sb, "_prod(_p, _input, _prodType(), %d, (double) ", ii);
            if (maxSumProdN < ii){
              maxSumProdN = ii;
            }
          } else if (!strcmp("sum", v)){
	    sAppend(&sb, "_sum(_p, _pld, -__MAX_PROD__, _sumType(), %d, (double) ", ii);
            if (SumProdLD < ii){
              SumProdLD = ii;
            }
	  } else {
	    sAppend(&sb, "_%s(%d, (double) ", v, ii);
	  }
	  sAppend(&sbt, "%s(", v);
          Free(v);
          i = 1;// Parse next arguments
	  depth=1;
	  continue;
        }
        Free(v);
      }
      
      if (!strcmp("theta",name)){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        sprintf(buf,"_THETA_%s_",v);
	ii = strtoimax(v,NULL,10);
	if (ii > tb.maxtheta){
	  tb.maxtheta =ii;
	}
	if (new_or_ith(buf)){
          sprintf(tb.ss+tb.pos, "%s,", buf);
          tb.pos += (int)strlen(buf)+1;
          tb.vo[++tb.nv] = tb.pos;
        }
        sAppend(&sb,"_THETA_%s_",v);
        sAppend(&sbt,"THETA[%s]",v);
        Free(v);
        continue;
      }

      if (!strcmp("eta",name)){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	ii = strtoimax(v,NULL,10);
        if (ii > tb.maxeta){
          tb.maxeta =ii;
        }
        sprintf(buf,"_ETA_%s_",v);
        if (new_or_ith(buf)){
	  sprintf(tb.ss+tb.pos, "%s,", buf);
          tb.pos += (int)strlen(buf)+1;
          tb.vo[++tb.nv] = tb.pos;
        }
        sAppend(&sb, "_ETA_%s_",v);
        sAppend(&sbt,"ETA[%s]",v);
        Free(v);
        continue;
      }
      wprint_parsetree(pt, xpn, depth, fn, client_data);
      if (rx_syntax_require_semicolon && !strcmp("end_statement",name) && i == 0){
        if (xpn->start_loc.s ==  xpn->end){
          trans_syntax_error_report_fn(NEEDSEMI);
        } 
      }
      
      if (!strcmp("mult_part",name)){
	char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	if (i == 0){
	  if (!strcmp("/",v)){
	    sAppendN(&sb,"safe_zero(", 10);
            safe_zero = 1;
	  } else {
	    safe_zero = 0;
	  }
	}
	if (i == 1){
	  if (safe_zero){
	    sAppendN(&sb, ")", 1);
	  }
	  safe_zero = 0;
	}
	Free(v);
      }
      if (!strcmp("print_command",name)){
        found_print = 1;
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if  (!strncmp(v,"print",5)){
	  sAppendN(&sbPm, "full_print;\n", 12);
	  sAppendN(&sbNrm, "print;\n", 7);
        } else {
	  sAppend(&sbPm, "%s;\n", v);
	  sAppend(&sbNrm, "%s;\n", v);
        }
        Free(v);
      }
      if (!strcmp("printf_statement",name)){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (i == 0){
          if (!strncmp(v,"ode0",4)){
	    sb.o =0;
	    sbt.o=0;
	    sAppendN(&sb, "ODE0_Rprintf(", 12);
            sAppendN(&sbt,"ode0_printf(", 12);
          } else if (!strncmp(v,"jac0",4)) {
	    sb.o =0;
	    sbt.o=0;
            sAppendN(&sb,"JAC0_Rprintf(", 12);
            sAppendN(&sbt,"jac0_printf(", 12);
          } else if (!strncmp(v,"ode",3)){
	    sb.o =0;
	    sbt.o=0;
            sAppendN(&sb,"ODE_Rprintf(", 11);
            sAppendN(&sbt,"ode_printf(", 11);
          } else if (!strncmp(v,"jac",3)){
	    sb.o =0;
	    sbt.o=0;
            sAppendN(&sb,"JAC_Rprintf(", 11);
            sAppendN(&sbt,"jac_printf(", 11);
          } else if (!strncmp(v,"lhs",3)){
	    sb.o =0;
	    sbt.o=0;
            sAppendN(&sb,"LHS_Rprintf(", 11);
            sAppendN(&sbt,"lhs_printf(", 11);
          } else {
	    sb.o =0;
	    sbt.o=0;
            sAppendN(&sb,"Rprintf(", 7);
            sAppendN(&sbt,"printf(", 7);
          }
        }
        if (i == 2){
          sAppend(&sb,"%s",v);
	  sAppend(&sbt,"%s",v);
        }
        if (i == 4){
	  sAppend(&sbPm, "%s;\n", sb.s);
	  sAppend(&sbNrm, "%s;\n", sbt.s);
        }
        Free(v);
        continue;
      } 

      if ( (!strcmp("jac",name) || !strcmp("jac_rhs",name) ||
            !strcmp("dfdy",name) || !strcmp("dfdy_rhs",name)) && i == 2){
        found_jac = 1;
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (!strcmp("jac_rhs",name) || !strcmp("dfdy_rhs",name)){
          // Continuation statement
          sAppend(&sb, "__PDStateVar__[[%s,",v);
          sAppend(&sbt,"df(%s)/dy(",v);
        } else {
          // New statment
          sb.o = 0;
          sbt.o = 0;
	  sAppend(&sb,"__PDStateVar__[[%s,",v);
          sAppend(&sbt,"df(%s)/dy(",v);
	  new_or_ith(v);
	  tb.cdf = tb.ix;
        }
        Free(v);
        continue;
      }
      if (!strcmp("factorial_exp",name) && i == 0){
        sb.o--;
        sAppendN(&sb, "exp(lgamma1p(", 13);
        continue;
      }
      if (!strcmp("lfactorial_exp",name) && i == 0){
        sAppendN(&sb, "lgamma1p(", 9);
        sAppendN(&sbt, "log((", 5);
        continue;
      }
      if (!strcmp("lfactorial_exp",name) && i == 2){
        sAppendN(&sb, ")", 1);
        sAppendN(&sbt, ")!)", 3);
        continue;
      }
      if (!strcmp("factorial_exp",name) && i == 3) {
        sb.o--;
        sAppendN(&sb, ")", 1);
        continue;
      }      
      if (!strcmp("factorial",name)){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        sAppend(&sb, "exp(lgamma1p(%s))",v);
        sAppend(&sbt, "%s!",v);
        Free(v);
        continue;
      }
      if (!strcmp("lfactorial",name)){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        sAppend(&sb, "lgamma1p(%s)",v);
        sAppend(&sbt, "log(%s!)",v);
        Free(v);
        continue;
      }
      if ((!strcmp("jac",name)  || !strcmp("jac_rhs",name) ||
           !strcmp("dfdy",name) || !strcmp("dfdy_rhs",name)) && i == 4){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	ii = 0;
	if (strstr(v,"THETA[") != NULL){
	  sprintf(buf,"_THETA_%.*s_",(int)(strlen(v))-7,v+6);
	  sAppend(&sbt, "%s)",v);
	  sAppend(&sb, "%s]]",buf);
	  ii = 1;
	} else if (strstr(v,"ETA[") != NULL) {
	  sprintf(buf,"_ETA_%.*s_",(int)(strlen(v))-5,v+4);
          sAppend(&sbt, "%s)",v);
          sAppend(&sb, "%s]]",buf);
          ii = 1;
        } else {
	  sAppend(&sb, "%s]]",v);
          sAppend(&sbt, "%s)",v);
        }
        if (!strcmp("jac",name) ||
            strcmp("dfdy",name) == 0){
          sAppendN(&sb ," = ", 3);
          sAppendN(&sbt ,"=", 1);
	  if (ii == 1){
	    new_or_ith(buf);
          } else {
	    new_or_ith(v);
          }
	  found = -1;
	  for (ii = 0; ii < tb.ndfdy; ii++){
            if (tb.df[ii] == tb.cdf && tb.dy[ii] == tb.ix){
	      found = ii;
	      break;
	    }
	  }
	  if (found < 0){
            tb.df[tb.ndfdy] = tb.cdf;
	    tb.dy[tb.ndfdy] = tb.ix;
	    tb.ndfdy = tb.ndfdy+1;
	    tb.cdf = -1;
          }
        }
        Free(v);
        continue;
      }
      
      //inits
      if (!strcmp("selection_statement", name) && i==1) {
	sb.o = 0;
	sbt.o = 0;
        sAppendN(&sb, "if (", 4);
        sAppendN(&sbt,"if (", 4);
        continue;
      }
      if (!strcmp("selection_statement", name) && i==3) {
        sAppend(&sb, " {", 2);
        sAppend(&sbt, "{", 1);
	sAppend(&sbPm, "%s\n", sb.s);
	sAppend(&sbNrm, "%s\n", sbt.s);
        continue;
      }
      if (!strcmp("selection_statement__8", name) && i==0) {
	sAppendN(&sbPm, "}\nelse {\n", 9);
	sAppendN(&sbNrm, "}\nelse {\n", 9);
        continue;
      }

      if (!strcmp("power_expression", name) && i==0) {
        sAppendN(&sb, ",", 1);
        sAppendN(&sbt, "^", 1);
      }
      if (!rx_syntax_star_pow && i == 1 &&!strcmp("power_expression", name)){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (!strcmp("**",v)){
          trans_syntax_error_report_fn(NEEDPOW);
        }
        Free(v);
      }
      if (!strcmp("transit2", name) && i == 0){
        sAppendN(&sb, "_transit3P(t, _cSub, ", 21);
        sAppendN(&sbt,"transit(", 8);
        rx_podo = 1;
      }
      if (!strcmp("transit3", name) && i == 0){
        sAppendN(&sb, "_transit4P(t, _cSub, ", 21);
        sAppendN(&sbt,"transit(", 8);
        rx_podo = 1;
      }
      if (!strcmp("derivative", name) && i==5) {
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	if (!strcmp("+", v) || 
	    !strcmp("-", v)){
          // = + is output  or = InfusionRate + is outupt.
        } else {
	  // = + is output  or = InfusionRate + is outupt.
          sAppendN(&sb, "+ ", 2);
        }
	Free(v);
	continue;
      }
      if (!strcmp("derivative", name) && i==2) {
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        sprintf(tb.ddt, "%s",v);
        if (new_de(v)){
	  tb.statei++;
	  if (strstr(v, "rx__sens_")){
	    tb.sensi++;
	  }
	  sb.o=0;
          sAppend(&sb, "__DDtStateVar__[%d] = _IR[%d] ", tb.nd, tb.nd);
	  sbt.o=0;
          sAppend(&sbt, "d/dt(%s)", v);
	  new_or_ith(v);
          /* Rprintf("%s; tb.ini = %d; tb.ini0 = %d; tb.lh = %d\n",v,tb.ini[tb.ix],tb.ini0[tb.ix],tb.lh[tb.ix]); */
          if (!rx_syntax_allow_assign_state && ((tb.ini[tb.ix] == 1 && tb.ini0[tb.ix] == 0) || tb.lh[tb.ix] == 1)){
            sprintf(buf,"Cannot assign state variable %s; For initial condition assignment use '%s(0) = #'.\n  Changing states can break sensitivity analysis (for nlmixr glmm/focei).\n  To override this behavior set 'options(RxODE.syntax.assign.state = TRUE)'.\n",v,v);
            trans_syntax_error_report_fn(buf);
          }
	  tb.lh[tb.ix] = 9;
          tb.di[tb.nd] = tb.ix;
          sprintf(tb.de+tb.pos_de, "%s,", v);
          tb.pos_de += (int)strlen(v)+1;
	  Free(v);
	  xpn = d_get_child(pn,4);
          v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
          if (!strcmp("~",v)){
            tb.idi[tb.nd] = 1;
	    sAppendN(&sbt, "~", 1);
          } else {
	    tb.idi[tb.nd] = 0;
	    sAppendN(&sbt, "=", 1);
	  }
          tb.deo[++tb.nd] = tb.pos_de;
        } else {
	  new_or_ith(v);
	  /* printf("de[%d]->%s[%d]\n",tb.id,v,tb.ix); */
	  sb.o=0;
          sAppend(&sb, "__DDtStateVar__[%d] = ", tb.id);
	  sbt.o=0;
          sAppend(&sbt, "d/dt(%s)=", v);
        }
        Free(v);
        continue;
      }
      if (!strcmp("der_rhs", name)) {
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (new_de(v)){
          sprintf(buf,"Tried to use d/dt(%s) before it was defined",v);
          trans_syntax_error_report_fn(buf);
        } else {
	  sb.o=0;
          sAppend(&sb, "__DDtStateVar__[%d]", tb.id);
	  sbt.o = 0;
          sAppend(&sbt, "d/dt(%s)", v);
        }
        Free(v);
        continue;
      }

      if (!strcmp("inf_rhs", name)) {
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (new_de(v)){
          sprintf(buf,"Tried to use rxRate(%s) before d/dt(%s) was defined",v,v);
          trans_syntax_error_report_fn(buf);
        } else {
	  if (strcmp(tb.ddt, v)){
	    sAppend(&sb, "_InfusionRate[%d]", tb.id);
            sAppend(&sbt, "rxRate(%s)", v);
          } else {
	    sAppendN(&sb, "0.0", 3);
            sAppend(&sbt, "rxRate(%s)", v);
	  }
        }
        Free(v);
        continue;
      }

      if (!strcmp("ini0f", name) && rx_syntax_allow_ini && i == 0){
	sprintf(sb.s,"(__0f__)");
	sb.o = 8;
	sbt.o = 0;
	char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	sAppend(&sb,"%s",v);
	sAppend(&sbt,"%s(0)",v);
      }

      if ((!strcmp("assignment", name) || !strcmp("ini", name) || !strcmp("ini0", name)) && i==0) {
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	tb.ddt[0]='\0';
        if ((rx_syntax_allow_ini && !strcmp("ini", name)) || !strcmp("ini0", name)){
	  sb.o=0;
          sAppendN(&sb,"(__0__)", 7);
	  tb.ini_i++;
          for (k = 0; k < (int)strlen(v); k++){
            if (v[k] == '.'){
              if (rx_syntax_allow_dots){
                sAppendN(&sb,"_DoT_", 5);
              } else {
                trans_syntax_error_report_fn(NODOT);
              }
            } else {
              sPut(&sb, v[k]);
            }
          }
          if (!strcmp("ini",name) && !new_de(v)){
            sprintf(buf,"Cannot assign state variable %s; For initial condition assigment use '%s(0) ='.\n",v,v);
            trans_syntax_error_report_fn(buf);
          }
          if (!rx_syntax_allow_ini0 && !strcmp("ini0",name)){
            sprintf(buf,NOINI0,v);
            trans_syntax_error_report_fn(buf);
          }
        } else {
          sb.o = 0;
          for (k = 0; k < (int)strlen(v); k++){
            if (v[k] == '.'){
              if (rx_syntax_allow_dots){              
                sAppendN(&sb,"_DoT_", 5);
              } else {
                trans_syntax_error_report_fn(NODOT);
              }
            } else {
              sPut(&sb, v[k]);
            }
          }
          if (!new_de(v)){
            sprintf(buf,"Cannot assign state variable %s; For initial condition assigment use '%s(0) ='.\n",v,v);
            trans_syntax_error_report_fn(buf);
          }
        }
	sbt.o=0;
        sAppend(&sbt, "%s", v);
	if (!strcmp("ini0",name)){
	  sbt.o=0;
	  sAppend(&sbt,"%s(0)",v);
	}
	new_or_ith(v);
	if (!strcmp("assignment", name)  || (!rx_syntax_allow_ini && !strcmp("ini", name))){
          tb.lh[tb.ix] = 1;
        } else if (!strcmp("ini", name) || !strcmp("ini0",name)){
          if (tb.ini[tb.ix] == 0){
            // If there is only one initialzation call, then assume
            // this is a parameter with an initial value.
            tb.ini[tb.ix] = 1;
            if (!strcmp("ini0",name)){
	      tb.ini0[tb.ix] = 1;
            } else {
	      tb.ini0[tb.ix] = 0;
              if (strstr(v,"rx_") != NULL){
                tb.lh[tb.ix] = 1;
              }
            }
          } else {
            // There is more than one call to this variable, it is a
            // conditional variabile
            tb.lh[tb.ix] = 1;
            if (!strcmp("ini0", name) && tb.ini0[tb.ix] == 1){
              sprintf(buf,"Cannot have conditional initial conditions for %s",v);
              trans_syntax_error_report_fn(buf);
            }
	    tb.ini0[tb.ix] = 0;
          }
        }
        Free(v);
      }
    }

    if (!strcmp("assignment", name) || !strcmp("ini", name) || !strcmp("derivative", name) || !strcmp("jac",name) || !strcmp("dfdy",name) ||
        !strcmp("ini0",name) || !strcmp("ini0f",name)){
      sAppend(&sbPm, "%s;\n", sb.s);
      sAppend(&sbNrm, "%s;\n", sbt.s);
    }

    if (!rx_syntax_assign && (!strcmp("assignment", name) || !strcmp("ini", name) || !strcmp("ini0",name) || !strcmp("ini0f",name))){
      if (!strcmp("ini0",name)){
        i = 2;
      } else {
        i = 1;
      }
      D_ParseNode *xpn = d_get_child(pn,i);
      char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
      if (!strcmp("<-",v)){
        trans_syntax_error_report_fn(NOASSIGN);
      }
      Free(v);
    }
    
    if (!strcmp("selection_statement", name)){
      sAppendN(&sbPm, "}\n", 2);
      sAppendN(&sbNrm, "}\n", 2);
    }
    
    if (!strcmp("power_expression", name)) {
      sAppendN(&sb, ")", 1);
    }

  }
}

void retieve_var(int i, char *buf) {
  int len;

  len = tb.vo[i+1] - tb.vo[i] - 1;
  strncpy(buf, tb.ss+tb.vo[i], len);
  buf[len] = 0;
}

void err_msg(int chk, const char *msg, int code)
{
  if(!chk) {
    error("%s",msg);
  }
}

/* when prnt_vars() is called, user defines the behavior in "case" */
void prnt_vars(int scenario, int lhs, const char *pre_str, const char *post_str, int show_ode) {
  int i, j, k;
  char buf[64], buf1[64],buf2[64];
  sAppend(&sbOut[show_ode], "%s", pre_str);
  if (scenario == 0 || scenario == 2){
    // show_ode = 1 dydt
    // show_ode = 2 Jacobian
    // show_ode = 3 Ini statement
    // show_ode = 0 LHS
    if (show_ode == 2 || show_ode == 0){
      //__DDtStateVar_#__
      for (i = 0; i < tb.nd; i++){
	if (scenario == 0){
	  sAppend(&sbOut[show_ode],"  __DDtStateVar_%d__,\n",i);
	} else {
	  sAppend(&sbOut[show_ode],"  (void)__DDtStateVar_%d__;\n",i);
	}
      }
    }
    // Now get Jacobain information  __PDStateVar_df_dy__ if needed
    if (show_ode != 3){
      for (i = 0; i < tb.ndfdy; i++){
        retieve_var(tb.df[i], buf1);
        retieve_var(tb.dy[i], buf2);
        // This is for dydt/ LHS/ or jacobian for df(state)/dy(parameter)
        if (show_ode == 1 || show_ode == 0 || tb.sdfdy[i] == 1){
	  if (scenario == 0){
	    sAppend(&sbOut[show_ode],"  __PDStateVar_%s_SeP_%s__,\n",buf1,buf2);
          } else {
	    sAppend(&sbOut[show_ode],"  (void)__PDStateVar_%s_SeP_%s__;\n",buf1,buf2);
	  }
        }
      }
    }
  }
  for (i=0, j=0; i<tb.nv; i++) {
    if (lhs && tb.lh[i]>0) continue;
    retieve_var(i, buf);
    switch(scenario) {
    case 0:   // Case 0 is for declaring the variables
      sAppendN(&sbOut[show_ode],"  ", 2);
      for (k = 0; k < (int)strlen(buf); k++){
        if (buf[k] == '.'){
          sAppend(&sbOut[show_ode],"_DoT_");
          if (!rx_syntax_allow_dots){
            trans_syntax_error_report_fn(NODOT);
          }
        } else {
          sPut(&sbOut[show_ode],buf[k]);
        }
      }
      if (!strcmp("rx_lambda_", buf) || !strcmp("rx_yj_", buf)){
	sAppendN(&sbOut[show_ode], "__", 2);
      }
      if (i <tb.nv-1)
        sAppendN(&sbOut[show_ode], ",\n", 2);
      else
        sAppendN(&sbOut[show_ode], ";\n", 2);
      break;
    case 2: // Case 2 is for suppressing all the warnings for the variables by using (void)var;
      // See https://stackoverflow.com/questions/1486904/how-do-i-best-silence-a-warning-about-unused-variables
      sAppend(&sbOut[show_ode],"  ");
      sAppend(&sbOut[show_ode],"(void)");
      for (k = 0; k < (int)strlen(buf); k++){
        if (buf[k] == '.'){
          sAppendN(&sbOut[show_ode],"_DoT_", 5);
          if (!rx_syntax_allow_dots){
            trans_syntax_error_report_fn(NODOT);
          }
        } else {
          sPut(&sbOut[show_ode],buf[k]);
        }
      }
      if (!strcmp("rx_lambda_", buf) || !strcmp("rx_yj_", buf)){
        sAppendN(&sbOut[show_ode], "__", 2);
      }
      sAppendN(&sbOut[show_ode], ";\n", 2);
      break;
    case 1:
      // Case 1 is for declaring the par_ptr.
      sAppendN(&sbOut[show_ode],"  ", 2);
      for (k = 0; k < (int)strlen(buf); k++){
        if (buf[k] == '.'){
          sAppendN(&sbOut[show_ode],"_DoT_", 5);
          if (!rx_syntax_allow_dots){
            trans_syntax_error_report_fn(NODOT);
          }
        } else {
          sPut(&sbOut[show_ode], buf[k]);
        }
      }
      sAppend(&sbOut[show_ode], " = _PP[%d];\n", j++);
      break;
    default: break;
    }
  }
  sAppend(&sbOut[show_ode], "%s", post_str);
}

void print_aux_info(char *model, const char *prefix, const char *libname, const char *pMd5, const char *timeId){
  int i, j, k, islhs,pi = 0,li = 0, o=0, o2=0, statei = 0, ini_i = 0, sensi=0, normi=0,fdi=0,
    in_str=0;
  char *s2;
  char sLine[MXLEN+1];
  char buf[512], buf2[512];
  for (i=0; i<tb.nv; i++) {
    islhs = tb.lh[i];
    if (islhs>1) continue;      /* is a state var */
    retieve_var(i, buf);
    if (islhs == 1){
      sprintf(s_aux_info+o, "  SET_STRING_ELT(lhs,%d,mkChar(\"%s\"));\n", li++, buf);
    } else {
      for (j = 1; j <= tb.maxtheta;j++){
        sprintf(buf2,"_THETA_%d_",j);
        if (!strcmp(buf,buf2)){
          sprintf(buf,"THETA[%d]",j);
        }
      }
      for (j = 1; j <= tb.maxeta;j++){
        sprintf(buf2,"_ETA_%d_",j);
        if (!strcmp(buf,buf2)){
          sprintf(buf,"ETA[%d]",j);
        }
      }
      sprintf(s_aux_info+o, "    SET_STRING_ELT(params,%d,mkChar(\"%s\"));\n", pi++, buf);
    }
    o = (int)strlen(s_aux_info);
  }
  for (i=0; i<tb.nd; i++) {                     /* name state vars */
    retieve_var(tb.di[i], buf);
    if (strstr(buf, "rx__sens_")){
      sprintf(s_aux_info+o, "    SET_STRING_ELT(sens,%d,mkChar(\"%s\"));\n", sensi++, buf);
      o = (int)strlen(s_aux_info);
      sprintf(s_aux_info+o, "    SET_STRING_ELT(state,%d,mkChar(\"%s\"));\n", statei++, buf);
      o = (int)strlen(s_aux_info);
      sprintf(s_aux_info+o, "    _SR[%d] = %d;\n", statei-1, tb.idi[i]);
    } else {
      sprintf(s_aux_info+o, "    SET_STRING_ELT(state,%d,mkChar(\"%s\"));\n", statei++, buf);
      o = (int)strlen(s_aux_info);
      sprintf(s_aux_info+o, "    SET_STRING_ELT(normState,%d,mkChar(\"%s\"));\n", normi++, buf);
      o = (int)strlen(s_aux_info);
      sprintf(s_aux_info+o, "    _SR[%d] = %d;\n", statei-1, tb.idi[i]);
    }
    if (tb.fdi[i]){
      o = (int)strlen(s_aux_info);
      sprintf(s_aux_info+o, "    SET_STRING_ELT(fn_ini,%d,mkChar(\"%s\"));\n", fdi++, buf);
    }
    o = (int)strlen(s_aux_info);
  }
  for (i=0; i<tb.ndfdy; i++) {                     /* name state vars */
    retieve_var(tb.df[i], buf);
    sprintf(s_aux_info+o, "    SET_STRING_ELT(dfdy,%d,mkChar(\"df(%s)/dy(", i, buf);
    o = (int)strlen(s_aux_info);
    retieve_var(tb.dy[i], buf);
    for (j = 1; j <= tb.maxtheta;j++){
      sprintf(buf2,"_THETA_%d_",j);
      if (!strcmp(buf,buf2)){
        sprintf(buf,"THETA[%d]",j);
      }
    }
    for (j = 1; j <= tb.maxeta;j++){
      sprintf(buf2,"_ETA_%d_",j);
      if (!strcmp(buf,buf2)){
        sprintf(buf,"ETA[%d]",j);
      }
    }
    sprintf(s_aux_info+o, "%s)\"));\n",buf);
    o = (int)strlen(s_aux_info);
  }
  sAppend(&sbOut[4], "extern SEXP %smodel_vars(){\n  int pro=0;\n", prefix);
  sAppend(&sbOut[4], "  SEXP _mv = PROTECT(_rxGetModelLib(\"%smodel_vars\"));pro++;\n", prefix);
  sAppendN(&sbOut[4], "  if (!_rxIsCurrentC(_mv)){\n", 28);
  sAppendN(&sbOut[4], "    SEXP lst      = PROTECT(allocVector(VECSXP, 16));pro++;\n", 60);
  sAppendN(&sbOut[4], "    SEXP names    = PROTECT(allocVector(STRSXP, 16));pro++;\n", 60);
  sAppend(&sbOut[4], "    SEXP params   = PROTECT(allocVector(STRSXP, %d));pro++;\n",pi);
  sAppend(&sbOut[4], "    SEXP lhs      = PROTECT(allocVector(STRSXP, %d));pro++;\n",li);
  sAppend(&sbOut[4], "    SEXP state    = PROTECT(allocVector(STRSXP, %d));pro++;\n",statei);
  sAppend(&sbOut[4], "    SEXP stateRmS = PROTECT(allocVector(INTSXP, %d));pro++;\n",statei);
  sAppendN(&sbOut[4], "    SEXP timeInt = PROTECT(allocVector(INTSXP, 1));pro++;\n", 58);
  sAppend(&sbOut[4], "    INTEGER(timeInt)[0] = %s;\n", timeId);
  sAppend(&sbOut[4], "    SEXP sens     = PROTECT(allocVector(STRSXP, %d));pro++;\n",sensi);
  sAppend(&sbOut[4], "    SEXP normState= PROTECT(allocVector(STRSXP, %d));pro++;\n",statei-sensi);
  sAppend(&sbOut[4], "    SEXP fn_ini   = PROTECT(allocVector(STRSXP, %d));pro++;\n",fdi);
  sAppend(&sbOut[4], "    SEXP dfdy     = PROTECT(allocVector(STRSXP, %d));pro++;\n",tb.ndfdy);
  sAppendN(&sbOut[4], "    SEXP tran     = PROTECT(allocVector(STRSXP, 14));pro++;\n", 60);
  sAppendN(&sbOut[4], "    SEXP trann    = PROTECT(allocVector(STRSXP, 14));pro++;\n", 60);
  sAppendN(&sbOut[4], "    SEXP mmd5     = PROTECT(allocVector(STRSXP, 2));pro++;\n", 59);
  sAppendN(&sbOut[4], "    SEXP mmd5n    = PROTECT(allocVector(STRSXP, 2));pro++;\n", 59);
  sAppendN(&sbOut[4], "    SEXP model    = PROTECT(allocVector(STRSXP, 1));pro++;\n", 59);
  sAppendN(&sbOut[4], "    SEXP modeln   = PROTECT(allocVector(STRSXP, 1));pro++;\n", 59);
  sAppendN(&sbOut[4], "    SEXP solve    = PROTECT(allocVector(VECSXP, 4));pro++;\n", 59);
  sAppendN(&sbOut[4], "    SEXP solven   = PROTECT(allocVector(STRSXP, 4));pro++;\n", 59);
  sAppend(&sbOut[4], "    SEXP initsr   = PROTECT(allocVector(REALSXP, %d));pro++;\n",statei);
  sAppend(&sbOut[4], "    SEXP scaler   = PROTECT(allocVector(REALSXP, %d));pro++;\n",statei);
  sAppend(&sbOut[4], "    SEXP infusionr= PROTECT(allocVector(REALSXP, %d));pro++;\n",statei);
  sAppend(&sbOut[4], "    SEXP badDosei = PROTECT(allocVector(INTSXP, %d));pro++;\n",statei);
  sAppendN(&sbOut[4], "    SEXP version    = PROTECT(allocVector(STRSXP, 3));pro++;\n", 61);
  sAppendN(&sbOut[4], "    SEXP versionn   = PROTECT(allocVector(STRSXP, 3));pro++;\n", 61);
  
  sAppend(&sbOut[4],  __VER_0__);
  sAppend(&sbOut[4],  __VER_1__);
  sAppend(&sbOut[4],  __VER_2__);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(versionn,0,mkChar(\"version\"));\n", 50);
  sAppendN(&sbOut[4], "    SET_STRING_ELT(versionn,1,mkChar(\"repo\"));\n", 47);
  sAppendN(&sbOut[4], "    SET_STRING_ELT(versionn,2,mkChar(\"md5\"));\n", 46);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(solven,0,mkChar(\"inits\"));\n", 46);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(solve,  0,initsr);\n", 38);
  sAppendN(&sbOut[4], "    SET_STRING_ELT(solven,1,mkChar(\"scale\"));\n", 46);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(solve,  1,scaler);\n", 38);
  sAppendN(&sbOut[4], "    SET_STRING_ELT(solven,2,mkChar(\"infusion\"));\n", 49);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(solve,  2,infusionr);\n", 41);
  sAppendN(&sbOut[4], "    SET_STRING_ELT(solven,3,mkChar(\"badDose\"));\n", 48);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(solve,  3,badDosei);\n", 40);
  sAppendN(&sbOut[4], "    setAttrib(solve, R_NamesSymbol, solven);\n", 45);
  sAppend(&sbOut[4], "%s",s_aux_info);
  // Save for outputting in trans
  tb.fdn = fdi;
  tb.pi = pi;
  tb.li = li;
  tb.sensi  = sensi;
  sAppendN(&sbOut[4], "    SET_STRING_ELT(modeln,0,mkChar(\"normModel\"));\n", 50);
  sAppendN(&sbOut[4], "    SET_STRING_ELT(model,0,mkChar(\"", 35);
  in_str=0;
  for (i = 0; i < sbNrm.o; i++){
    if (sbNrm.s[i] == '"'){
      if (in_str==1){
	in_str=0;
      } else {
	in_str=1;
      }
      sAppendN(&sbOut[4], "\\\"", 2);
    } else if (sbNrm.s[i] == '\''){
      if (in_str==1){
	in_str=0;
      } else {
	in_str=1;
      }
      sAppendN(&sbOut[4], "'", 1);
    } else if (sbNrm.s[i] == ' '){
      if (in_str==1){
	sAppendN(&sbOut[4], " ", 1);
      }
    } else if (sbNrm.s[i] == '\n'){
      sAppendN(&sbOut[4], "\\n", 2);
    } else if (sbNrm.s[i] == '\t'){
      sAppendN(&sbOut[4], "\\t", 2);
    } else if (sbNrm.s[i] == '\\'){
      sAppendN(&sbOut[4], "\\\\", 2);
    } else if (sbNrm.s[i] >= 33  && sbNrm.s[i] <= 126){ // ASCII only
      sPut(&sbOut[4], sbNrm.s[i]);
    }
  }
  sAppendN(&sbOut[4], "\"));\n", 5);
  
  s_aux_info[0] = '\0';
  o    = 0;
  
  sbPm.o=0;
  while(sgets(sLine, MXLEN, &sbPm)) {
    s2 = strstr(sLine,"(__0__)");
    if (s2){
      // See if this is a reclaimed initilization variable.
      for (i=0; i<tb.nv; i++) {
        if (tb.ini[i] == 1 && tb.lh[i] != 1){
          //(__0__)V2 =
          retieve_var(i, buf);
          sprintf(buf2,"(__0__)");
          o2 = 7;
          for (k = 0; k < (int)strlen(buf); k++){
            if (buf[k] == '.'){
              sprintf(buf2+o2,"_DoT_");
              if (!rx_syntax_allow_dots){
                trans_syntax_error_report_fn(NODOT);
              }
              o2+=5;
            } else {
              sprintf(buf2+o2,"%c",buf[k]);
              o2++;
            }
          }
          sprintf(buf2+o2,"=");
          s2 = strstr(sLine,buf2);
          if (s2){
            sprintf(s_aux_info+o,"    SET_STRING_ELT(inin,%d,mkChar(\"%s\"));\n",ini_i, buf);
            o = (int)strlen(s_aux_info);
            sprintf(s_aux_info+o,"    REAL(ini)[%d] = %.*s;\n",(int)(ini_i++), (int)(strlen(sLine))-(int)(strlen(buf2))-2,sLine + (int)(strlen(buf2)));
            o = (int)strlen(s_aux_info);
            continue;
          }
        }
      }
      continue;
    }
  }
  // putin constants
  for (i=0; i<tb.nv; i++) {
    if (tb.ini[i] == 0 && tb.lh[i] != 1) {
      retieve_var(i, buf);
      // Put in constants
      if  (!strcmp("pi",buf)){
        sprintf(s_aux_info+o,"    SET_STRING_ELT(inin,%d,mkChar(\"pi\"));\n",ini_i);
        o = (int)strlen(s_aux_info);
        // Use well more digits than double supports
        sprintf(s_aux_info+o,"    REAL(ini)[%d] = M_PI;\n",ini_i++);
        o = (int)strlen(s_aux_info);
      }
    }
  }
  tb.ini_i = ini_i;
  sAppend(&sbOut[4], "    SEXP ini    = PROTECT(allocVector(REALSXP,%d));pro++;\n",ini_i);
  sAppend(&sbOut[4], "    SEXP inin   = PROTECT(allocVector(STRSXP, %d));pro++;\n",ini_i);
  sAppend(&sbOut[4], "%s",s_aux_info);
  // Vector Names
  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,0,mkChar(\"params\"));\n", 46);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  0,params);\n", 36);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,1,mkChar(\"lhs\"));\n", 43);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  1,lhs);\n", 33);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,2,mkChar(\"state\"));\n", 45);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  2,state);\n", 35);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,3,mkChar(\"trans\"));\n", 45);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  3,tran);\n", 34);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,5,mkChar(\"model\"));\n", 45);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  5,model);\n", 35);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,4,mkChar(\"ini\"));\n", 43);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  4,ini);\n", 33);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,6,mkChar(\"md5\"));\n", 43);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  6,mmd5);\n", 34);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,7,mkChar(\"podo\"));\n", 44);
  sAppend(&sbOut[4], "    SET_VECTOR_ELT(lst,  7,ScalarLogical(%d));\n",rx_podo);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,8,mkChar(\"dfdy\"));\n", 44);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  8,dfdy);\n", 34);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,9,mkChar(\"sens\"));\n", 44);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  9,sens);\n", 34);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,10,mkChar(\"fn.ini\"));\n", 47);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  10,fn_ini);\n", 37);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,11,mkChar(\"state.ignore\"));\n", 53);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  11,stateRmS);\n", 39);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,12,mkChar(\"solve\"));\n", 46);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  12,solve);\n", 36);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,13,mkChar(\"version\"));\n", 48);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  13,version);\n", 38);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,14,mkChar(\"normal.state\"));\n", 53);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  14,normState);\n", 40);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(names,15,mkChar(\"timeId\"));\n", 47);
  sAppendN(&sbOut[4], "    SET_VECTOR_ELT(lst,  15,timeInt);\n", 38);

  // const char *rxVersion(const char *what)
  
  // md5 values
  sAppendN(&sbOut[4], "    SET_STRING_ELT(mmd5n,0,mkChar(\"file_md5\"));\n", 48);
  sAppend(&sbOut[4], "    SET_STRING_ELT(mmd5,0,mkChar(\"%s\"));\n",md5);
  sAppendN(&sbOut[4], "    SET_STRING_ELT(mmd5n,1,mkChar(\"parsed_md5\"));\n", 50);
  sAppend(&sbOut[4], "    SET_STRING_ELT(mmd5,1,mkChar(\"%s\"));\n", pMd5);
  
  // now trans output
  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,0,mkChar(\"lib.name\"));\n", 48);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 0,mkChar(\"%s\"));\n", libname);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,1,mkChar(\"jac\"));\n", 43);
  if (found_jac == 1){
    sAppendN(&sbOut[4], "    SET_STRING_ELT(tran,1,mkChar(\"fulluser\"));\n", 47); // Full User Matrix
  } else {
    sAppendN(&sbOut[4], "    SET_STRING_ELT(tran,1,mkChar(\"fullint\"));\n", 46); // Full Internal Matrix
  }
  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,2,mkChar(\"prefix\"));\n", 46);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 2,mkChar(\"%s\"));\n", prefix);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,3,mkChar(\"dydt\"));\n", 44);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 3,mkChar(\"%sdydt\"));\n", prefix);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,4,mkChar(\"calc_jac\"));\n", 48);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 4,mkChar(\"%scalc_jac\"));\n", prefix);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,5,mkChar(\"calc_lhs\"));\n", 48);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 5,mkChar(\"%scalc_lhs\"));\n", prefix);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,6,mkChar(\"model_vars\"));\n", 50);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 6,mkChar(\"%smodel_vars\"));\n", prefix);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,7,mkChar(\"ode_solver\"));\n", 50);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 7,mkChar(\"%sode_solver\"));\n", prefix);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,8,mkChar(\"inis\"));\n", 44);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 8,mkChar(\"%sinis\"));\n", prefix);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,  9,mkChar(\"dydt_lsoda\"));\n", 52);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran,   9,mkChar(\"%sdydt_lsoda\"));\n", prefix);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,10,mkChar(\"calc_jac_lsoda\"));\n", 55);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 10,mkChar(\"%scalc_jac_lsoda\"));\n", prefix);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,11,mkChar(\"ode_solver_solvedata\"));\n", 61);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 11,mkChar(\"%sode_solver_solvedata\"));\n", prefix);
  
  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,12,mkChar(\"ode_solver_get_solvedata\"));\n", 65);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 12,mkChar(\"%sode_solver_get_solvedata\"));\n", prefix);

  sAppendN(&sbOut[4], "    SET_STRING_ELT(trann,13,mkChar(\"dydt_liblsoda\"));\n", 54);
  sAppend(&sbOut[4], "    SET_STRING_ELT(tran, 13,mkChar(\"%sdydt_liblsoda\"));\n", prefix);
  
  sAppendN(&sbOut[4], "    setAttrib(tran, R_NamesSymbol, trann);\n", 43);
  sAppendN(&sbOut[4], "    setAttrib(mmd5, R_NamesSymbol, mmd5n);\n", 43);
  sAppendN(&sbOut[4], "    setAttrib(model, R_NamesSymbol, modeln);\n", 45);
  sAppendN(&sbOut[4], "    setAttrib(ini, R_NamesSymbol, inin);\n", 41);
  sAppendN(&sbOut[4], "    setAttrib(version, R_NamesSymbol, versionn);\n", 49);
  sAppendN(&sbOut[4], "    setAttrib(lst, R_NamesSymbol, names);\n", 42);
  sAppendN(&sbOut[4], "    SEXP cls = PROTECT(allocVector(STRSXP, 1));pro++;\n", 54);
  sAppendN(&sbOut[4], "    SET_STRING_ELT(cls, 0, mkChar(\"rxModelVars\"));\n", 51);
  sAppendN(&sbOut[4], "    classgets(lst, cls);\n", 25);
  sAppendN(&sbOut[4], "    _assign_ptr(lst);\n", 22);
  sAppendN(&sbOut[4], "    UNPROTECT(pro);\n", 20);
  
  sAppendN(&sbOut[4], "    return lst;\n", 16);
  sAppendN(&sbOut[4], "  } else {\n", 11);
  sAppendN(&sbOut[4], "    UNPROTECT(pro);\n", 20);
  sAppendN(&sbOut[4], "    return _mv;\n", 16);
  sAppendN(&sbOut[4], "  }\n", 4);
  sAppendN(&sbOut[4], "}\n", 2);

  sAppend(&sbOut[4],"extern void %sdydt_lsoda(int *neq, double *t, double *A, double *DADT)\n{\n  %sdydt(neq, *t, A, DADT);\n}\n", prefix, prefix);
  sAppend(&sbOut[4], "extern int %sdydt_liblsoda(double t, double *y, double *ydot, void *data)\n{\n  int *neq = (int*)(data);\n  %sdydt(neq, t, y, ydot);\n  return(0);\n}\n",
	  prefix,prefix);
  sAppend(&sbOut[4],"extern void %scalc_jac_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){\n  // Update all covariate parameters\n  %scalc_jac(neq, *t, A, JAC, *nrowpd);\n}\n",
	  prefix, prefix);
  sAppend(&sbOut[4],"\n//Initilize the dll to match RxODE's calls\nvoid R_init_%s(DllInfo *info){\n  // Get C callables on load; Otherwise it isn't thread safe\n\n", libname);
  sAppendN(&sbOut[4],"  _assign_ptr=(RxODE_assign_ptr) R_GetCCallable(\"RxODE\",\"RxODE_assign_fn_pointers\");\n  _rxRmModelLib=(_rxRmModelLibType) R_GetCCallable(\"RxODE\",\"rxRmModelLib\");\n  _rxGetModelLib=(_rxGetModelLibType) R_GetCCallable(\"RxODE\",\"rxGetModelLib\");\n  _old_c = (RxODE_ode_solver_old_c) R_GetCCallable(\"RxODE\",\"rxSolveOldC\");\n  _RxODE_rxAssignPtr=(_rx_asgn)R_GetCCallable(\"RxODE\",\"_RxODE_rxAssignPtr\");\n  _rxIsCurrentC = (_rxIsCurrentC_type)R_GetCCallable(\"RxODE\",\"rxIsCurrentC\");\n  _sumPS  = (_rxSumType) R_GetCCallable(\"PreciseSums\",\"PreciseSums_sum_r\");\n  _prodPS = (_rxProdType) R_GetCCallable(\"PreciseSums\",\"PreciseSums_prod_r\");\n  _prodType=(RxODE_fn0i)R_GetCCallable(\"PreciseSums\", \"PreciseSums_prod_get\");\n  _sumType=(RxODE_fn0i)R_GetCCallable(\"PreciseSums\", \"PreciseSums_sum_get\");\n  _ptrid=(RxODE_fn0i)R_GetCCallable(\"RxODE\", \"RxODE_current_fn_pointer_id\");\n  _powerD=(RxODE_fn3i)R_GetCCallable(\"RxODE\", \"powerD\");\n  _powerDi=(RxODE_fn3i)R_GetCCallable(\"RxODE\", \"powerDi\");\n  _powerDD=(RxODE_fn3i)R_GetCCallable(\"RxODE\", \"powerDD\");\n  _powerDDD=(RxODE_fn3i)R_GetCCallable(\"RxODE\", \"powerDDD\");\n  solveLinB=(solveLinB_p)R_GetCCallable(\"RxODE\", \"solveLinB\");\n  _update_par_ptr=(_update_par_ptr_p) R_GetCCallable(\"RxODE\",\"_update_par_ptr\");\n  // Register the outside functions", 1272);
  sAppend(&sbOut[4], "  R_RegisterCCallable(\"%s\",\"%sode_solver\", (DL_FUNC) %sode_solver);\n", libname, prefix, prefix);
  sAppend(&sbOut[4], "  R_RegisterCCallable(\"%s\",\"%sinis\",(DL_FUNC) %sinis);\n", libname, prefix, prefix);
  sAppend(&sbOut[4], "  R_RegisterCCallable(\"%s\",\"%sdydt\",(DL_FUNC) %sdydt);\n", libname, prefix, prefix);
  sAppend(&sbOut[4], "  R_RegisterCCallable(\"%s\",\"%scalc_lhs\",(DL_FUNC) %scalc_lhs);\n", libname, prefix, prefix);
  sAppend(&sbOut[4], "  R_RegisterCCallable(\"%s\",\"%scalc_jac\",(DL_FUNC) %scalc_jac);\n", libname, prefix, prefix);
  sAppend(&sbOut[4], "  R_RegisterCCallable(\"%s\",\"%sdydt_lsoda\", (DL_FUNC) %sdydt_lsoda);\n", libname, prefix, prefix);
  sAppend(&sbOut[4], "  R_RegisterCCallable(\"%s\",\"%scalc_jac_lsoda\", (DL_FUNC) %scalc_jac_lsoda);\n", libname, prefix, prefix);
  sAppend(&sbOut[4], "  R_RegisterCCallable(\"%s\",\"%sode_solver_solvedata\", (DL_FUNC) %sode_solver_solvedata);\n", libname, prefix, prefix);
  sAppend(&sbOut[4], "  R_RegisterCCallable(\"%s\",\"%sode_solver_get_solvedata\", (DL_FUNC) %sode_solver_get_solvedata);\n", libname, prefix, prefix);
  sAppend(&sbOut[4], "  R_RegisterCCallable(\"%s\",\"%sdydt_liblsoda\", (DL_FUNC) %sdydt_liblsoda);\n", libname, prefix, prefix);
  sAppend(&sbOut[4], "\n  static const R_CMethodDef cMethods[] = {\n    {\"%sode_solver\", (DL_FUNC) &%sode_solver, 15, %sode_solverrx_t},\n    {NULL, NULL, 0, NULL}\n  };\n",
	  prefix, prefix, prefix);
  sAppend(&sbOut[4], "\n  R_CallMethodDef callMethods[]  = {\n    {\"%smodel_vars\", (DL_FUNC) &%smodel_vars, 0},\n    {NULL, NULL, 0}\n  };\n",
	  prefix, prefix);
  sAppendN(&sbOut[4], "\nR_registerRoutines(info, cMethods, callMethods, NULL, NULL);\n  R_useDynamicSymbols(info,FALSE);\n}\n", 99);
  sAppend(&sbOut[4], "\n  void R_unload_%s (DllInfo *info){\n  // Free resources required for single subject solve.\n  SEXP _mv = PROTECT(_rxGetModelLib(\"%smodel_vars\"));\n",
	  libname, prefix);
  sAppend(&sbOut[4], "  if (!isNull(_mv)){\n    _rxRmModelLib(\"%smodel_vars\");\n  }\n  UNPROTECT(1);\n}\n", prefix);  
}


void codegen(char *model, int show_ode, const char *prefix, const char *libname, const char *pMd5, const char *timeId, const char *fixInis) {
  if (show_ode == 4) {
    print_aux_info(model, prefix, libname, pMd5, timeId);
  } else {
    int i, j, k, print_ode=0, print_vars = 0, print_parm = 0, print_jac=0, o;
    char sLine[MXLEN+1];
    char buf[64];
    char from[512], to[512], df[128], dy[128], state[128];
    char *s2;
    if (show_ode == 1){
      sAppendN(&sbOut[show_ode],"#include <RxODE_model.h>\n",25);
      sAppend(&sbOut[show_ode], "extern void  %sode_solver_solvedata (rx_solve *solve){\n  _solveData = solve;\n}\n",prefix);
      sAppend(&sbOut[show_ode], "extern rx_solve *%sode_solver_get_solvedata(){\n  return _solveData;\n}\n", prefix);
      sAppend(&sbOut[show_ode], "SEXP %smodel_vars();\n", prefix);
      sAppend(&sbOut[show_ode], "extern void %sode_solver(int *neq, double *theta, double *time, int *evid, int *ntime, double *inits, double *dose, double *ret, double *atol, double *rtol, int *stiff, int *transit_abs, int *nlhs, double *lhs, int *rc){\n  %s\n  _old_c(neq, _theta, time, evid, ntime, inits, dose, ret, atol, rtol, stiff, transit_abs, nlhs, lhs, rc);\n}\n", prefix, fixInis);
      sAppend(&sbOut[show_ode], "static R_NativePrimitiveArgType %sode_solverrx_t[] = {\n  INTSXP,REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP\n};\n", prefix);
      sAppendN(&sbOut[show_ode],"\n", 1);
      for (i = 0; i < (int)strlen(extra_buf); i++){
	if (extra_buf[i] == '"'){
	  sAppendN(&sbOut[show_ode],"\"", 1);
	} else if (extra_buf[i] == '\n'){
	  sAppendN(&sbOut[show_ode],"\n", 1);
	} else if (extra_buf[i] == '\t'){
	  sAppendN(&sbOut[show_ode],"\t", 1);
	} else if (extra_buf[i] >= 32  && extra_buf[i] <= 126){ // ASCII only
	  sPut(&sbOut[show_ode], extra_buf[i]);//fprintf(outpt,"%c",extra_buf[i]);
	}
      }
      sAppendN(&sbOut[show_ode], "\n// prj-specific differential eqns\nvoid ", 40);
      sAppend(&sbOut[show_ode], "%sdydt(int *_neq, double t, double *__zzStateVar__, double *__DDtStateVar__)\n{\n  int _cSub = _neq[1];\n", prefix);
    } else if (show_ode == 2){
      sAppend(&sbOut[show_ode], "// Jacobian derived vars\nvoid %scalc_jac(int *_neq, double t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {\n  int _cSub=_neq[1];\n", prefix);
    } else if (show_ode == 3){
      sAppend(&sbOut[show_ode],  "// Functional based initial conditions.\nvoid %sinis(int _cSub, double *__zzStateVar__){\n  double t=0;\n", prefix);
    } else {
      sAppend(&sbOut[show_ode],  "// prj-specific derived vars\nvoid %scalc_lhs(int _cSub, double t, double *__zzStateVar__, double *_lhs) {\n", prefix);
    }
    if (found_print){
      sAppendN(&sbOut[show_ode], "\n  int __print_ode__ = 0, __print_vars__ = 0,__print_parm__ = 0,__print_jac__ = 0;\n", 83);
    }
    if ((show_ode == 2 && found_jac == 1) || show_ode != 2){
      prnt_vars(0, 0, "  double ", "\n",show_ode);     /* declare all used vars */
      if (maxSumProdN > 0 || SumProdLD > 0){
	int mx = maxSumProdN;
	if (SumProdLD > mx) mx = SumProdLD;
	sAppend(&sbOut[show_ode],  "  double _p[%d], _input[%d];\n", mx, mx);
	sAppend(&sbOut[show_ode],  "  double _pld[%d];\n", mx);
      }
      prnt_vars(2, 0, "  (void)t;\n", "\n",show_ode);     /* declare all used vars */
      if (maxSumProdN){
	sAppendN(&sbOut[show_ode],  "  (void)_p;\n  (void)_input;\n", 28);
	if (SumProdLD){
	  sAppendN(&sbOut[show_ode],  "  (void)_pld;\n", 14);
	}
      }
      if (show_ode == 3){
	sAppendN(&sbOut[show_ode], "  _update_par_ptr(0.0, _cSub, _solveData, _idx);\n", 49);
      } else {
	sAppendN(&sbOut[show_ode], "  _update_par_ptr(t, _cSub, _solveData, _idx);\n", 47);
      }
      prnt_vars(1, 1, "", "\n",show_ode);                   /* pass system pars */
      for (i=0; i<tb.nd; i++) {                   /* name state vars */
	retieve_var(tb.di[i], buf);
	sAppendN(&sbOut[show_ode], "  ", 2);
	for (k = 0; k < (int)strlen(buf); k++){
	  if (buf[k] == '.'){
	    sAppendN(&sbOut[show_ode], "_DoT_", 5);
	    if (!rx_syntax_allow_dots){
	      trans_syntax_error_report_fn(NODOT);
	    }
	  } else {
	    sPut(&sbOut[show_ode], buf[k]);
	  }
	}
	sAppend(&sbOut[show_ode], " = __zzStateVar__[%d];\n", i);
      }
      sAppendN(&sbOut[show_ode], "\n", 1);
    
      sbPm.o=0;
      char *s;
      while(sgets(sLine, MXLEN, &sbPm)) {  /* parsed eqns */
	s = strstr(sLine,"(__0__)");
	if (s){
	  // See if this is a reclaimed initilization variable.
	  for (i=0; i<tb.nv; i++) {
	    if (tb.ini[i] == 1 && tb.lh[i] == 1){
	      //(__0__)V2=
	      retieve_var(i, buf);
	      s = strstr(sLine,buf);
	      if (s){
		sAppend(&sbOut[show_ode], "  %s\n",sLine + 7);
		continue;
	      }
	    }
	  }
	  continue;
	}
	if (show_ode == 3 && strstr(sLine,"full_print;")){
	  continue;
	}
	s = strstr(sLine,"(__0f__)");
	if (s){
	  if (show_ode == 3){
	    // FIXME
	    for (i = 0; i < tb.nd; i++){
	      retieve_var(tb.di[i], buf);
	      sprintf(to,"(__0f__)%s =",buf);
	      if (strstr(sLine,to)){
		if (!tb.fdi[i]){
		  tb.fdi[i] = 1;
		  tb.fdn++;
		}
		sAppend(&sbOut[show_ode],  "  //if (ISNA(%s)){ // Always apply; Othersiwse, Since this updates vector, it may cause this to only be updated once.\n    %s  //}\n",buf,sLine+8);
	      } else {
		sprintf(to,"(__0f__)%s=",buf);
		if (strstr(sLine,to)){
		  if (!tb.fdi[i]){
		    tb.fdi[i] = 1;
		    tb.fdn++;
		  }
		  sAppend(&sbOut[show_ode],  "  //if (ISNA(%s)){ // Always apply; Othersiwse, Since this updates vector, it may cause this to only be updated once.\n    %s  //}\n",buf,sLine+8);
		}
	      }
	    }
	    continue;
	  } else {
	    continue;
	  }
	}
	s = strstr(sLine,"ode_print;");
	if (show_ode == 1 && !s) s = strstr(sLine,"full_print;");
	if (show_ode != 1 && s) continue;
	else if (s) {
	  sAppendN(&sbOut[show_ode], "  Rprintf(\"================================================================================\\n\");\n", 97);
	  sAppendN(&sbOut[show_ode], "  Rprintf(\"ODE Count: %%d\\tTime (t): %%f\\n\", (&_solveData->subjects[_cSub])->dadt_counter[0], t);\n", 98);
	  sAppendN(&sbOut[show_ode], "  Rprintf(\"================================================================================\\n\");\n", 97);
	  sAppendN(&sbOut[show_ode], "  __print_ode__ = 1;\n", 21);
	  sAppendN(&sbOut[show_ode], "  __print_vars__ = 1;\n", 22);
	  sAppendN(&sbOut[show_ode], "  __print_parm__ = 1;\n", 22);
	  print_ode  = 1;
	  print_vars = 1;
	  print_parm = 1;
	  continue;
	}      
	s = strstr(sLine,"ODE_Rprintf");
	if ((show_ode != 1) && s) continue;
      
	s = strstr(sLine,"ODE0_Rprintf");
	if ((show_ode != 1) && s) continue;
      
	s = strstr(sLine, "__DDtStateVar__");
	if (s){
	  if (show_ode == 3){
	    continue;
	  }
	  if (show_ode!= 1){
	    for (i = 0; i < tb.nd; i++){
	      // Replace __DDtStateVar__[#] -> __DDtStateVar_#__
	      sprintf(to,"__DDtStateVar_%d__",i);
	      sprintf(from,"__DDtStateVar__[%d]",i);
	      s2 = repl_str(sLine,from,to);
	      strcpy(sLine, s2);
	      Free(s2);
	      s2=NULL;
	    }
	  }
	}
      
	s = strstr(sLine,"JAC_Rprintf");
	if ((show_ode != 2) && s) continue;

	s = strstr(sLine,"jac_print;");
	if (show_ode == 2 && !s) s = strstr(sLine,"full_print;");
	if (show_ode != 2 && s) continue;
	else if (s) {
	  sAppendN(&sbOut[show_ode], "  Rprintf(\"================================================================================\\n\");\n", 97);
	  sAppendN(&sbOut[show_ode], "  Rprintf(\"JAC Count: %%d\\tTime (t): %%f\\n\",(&_solveData->subjects[_cSub])->jac_counter[0], t);\n", 96);
	  sAppendN(&sbOut[show_ode], "  Rprintf(\"================================================================================\\n\");\n", 97);
	  sAppendN(&sbOut[show_ode], "  __print_ode__ = 1;\n", 21);
	  sAppendN(&sbOut[show_ode], "  __print_jac__ = 1;\n", 21);
	  sAppendN(&sbOut[show_ode], "  __print_vars__ = 1;\n", 22);
	  sAppendN(&sbOut[show_ode], "  __print_parm__ = 1;\n", 22);
	  print_ode  = 1;
	  print_vars = 1;
	  print_parm = 1;
	  print_jac = 1;
	  continue;
	}

	s = strstr(sLine,"JAC0_Rprintf");
	if ((show_ode != 2) && s) continue;

	s = strstr(sLine,"LHS_Rprintf");
	if ((show_ode != 0) && s) continue;

	s = strstr(sLine,"lhs_print;");
	if (show_ode == 0 && !s) s = strstr(sLine,"full_print;");
	if (show_ode != 0 && s) continue;
	else if (s) {
	  sAppendN(&sbOut[show_ode], "  Rprintf(\"================================================================================\\n\");\n", 97);
	  sAppendN(&sbOut[show_ode], "  Rprintf(\"LHS Time (t): %%f\\n\",t);\n", 36);
	  sAppendN(&sbOut[show_ode], "  Rprintf(\"================================================================================\\n\");\n", 97);
	  sAppendN(&sbOut[show_ode], "  __print_vars__ = 1;\n", 22);
	  sAppendN(&sbOut[show_ode], "  __print_parm__ = 1;\n", 22);
	  print_vars = 1;
	  print_parm = 1;
	  continue;
	}
      
	s = strstr(sLine,"__PDStateVar__");
	if (s){
	  if (show_ode == 3){
	    continue;
	  }
	  for (i = 0; i < tb.ndfdy; i++){
	    retieve_var(tb.df[i], df);
	    retieve_var(tb.dy[i], dy);
	    sprintf(from,"__PDStateVar__[[%s,%s]]",df,dy);
	    if (show_ode == 2 && tb.sdfdy[i] == 0){
	      // __PDStateVar__[__CMT_NUM_y__*(__NROWPD__)+__CMT_NUM_dy__]
	      sprintf(to,"__PDStateVar__[");
	      o = (int)strlen(to);
	      for (j=0; j<tb.nd; j++) {                     /* name state vars */
		retieve_var(tb.di[j], state);
		if (!strcmp(state, df)){
		  sprintf(to+o,"%d*(__NROWPD__)+",j);
		  o = (int)strlen(to);
		  break;
		}
	      }
	      for (j=0; j<tb.nd; j++){
		retieve_var(tb.di[j], state);
		if (!strcmp(state, dy)){
		  sprintf(to+o,"%d]",j);
		  o = (int)strlen(to);
		  break;
		}
	      }
	    } else {
	      sprintf(to,"__PDStateVar_%s_SeP_%s__",df,dy);
	    }
	    s2 = repl_str(sLine,from,to);
	    strcpy(sLine, s2);
	    Free(s2);
	    s2=NULL;
	  }        
	}
	// Replace __DDtStateVar__[#] -> __DDtStateVar_#__
	/* sprintf(to,""); */
	to[0]='\0';
	sprintf(from," ");
	s2 = repl_str(sLine,from,to);
	strcpy(sLine, s2);
	Free(s2);
	s2=NULL;
	int mx = maxSumProdN;
	if (SumProdLD > mx) mx = SumProdLD;
	sprintf(to,"%d", mx);
	sprintf(from, "__MAX_PROD__");
	s2 = repl_str(sLine,from,to);
	strcpy(sLine, s2);
	Free(s2);
	s2=NULL;
	sAppend(&sbOut[show_ode],  "  %s", sLine);
      }
    }
    if (print_ode && show_ode != 0){
      sAppendN(&sbOut[show_ode], "  if (__print_ode__ == 1){\n", 27);
      for (i=0; i<tb.nd; i++) {                   /* name state vars */
	retieve_var(tb.di[i], buf);
	sAppend(&sbOut[show_ode],  "    Rprintf(\"d/dt(%s)[%d]:\\t%%f\\t%s:\\t%%f\\n\", __DDtStateVar__[%d],%s);\n", buf, i,buf,i,buf);
      }
      sAppendN(&sbOut[show_ode], "  }\n", 4);
    }
    if (print_jac && show_ode == 2){
      sAppendN(&sbOut[show_ode], "  if (__print_jac__ == 1){\n", 27);
      sAppendN(&sbOut[show_ode], "  Rprintf(\"Fixme\\n\");", 21);
      sAppendN(&sbOut[show_ode], "  }\n", 4);
    }
    if (print_vars){
      sAppendN(&sbOut[show_ode], "  if (__print_vars__ == 1){\n", 28);
      sAppendN(&sbOut[show_ode], "    Rprintf(\".Left Handed Variables:.........................................................\\n\");\n", 99);
      for (i=0, j=0; i<tb.nv; i++) {
	if (tb.lh[i] != 1) continue;
	retieve_var(i, buf);
	sAppend(&sbOut[show_ode],  "    Rprintf(\"%s = %%f\\n\", %s);\n", buf, buf);
      }
      sAppendN(&sbOut[show_ode], "  }\n", 4);
    }
    if (print_parm){
      sAppendN(&sbOut[show_ode], "  if (__print_parm__ == 1){\n", 28);
      sAppendN(&sbOut[show_ode], "    Rprintf(\".User Supplied Variables:.......................................................\\n\");\n", 99);
      for (i=0, j=0; i<tb.nv; i++) {
	if (tb.lh[i]>0) continue;
	j++;
	retieve_var(i, buf);
	sAppend(&sbOut[show_ode],  "    Rprintf(\"%s=%%f\\t_par_ptr[%d]=%%f\\n\",%s,_PP[%d]);\n", buf, j-1, buf,j-1);
      }
      sAppendN(&sbOut[show_ode], "  }\n", 4);
    }
    if (print_jac || print_vars || print_ode || print_parm){
      sAppendN(&sbOut[show_ode], "  if (__print_jac__ || __print_vars__ || __print_ode__ || __print_parm__){\n", 75);
      sAppendN(&sbOut[show_ode], "    Rprintf(\"================================================================================\\n\\n\\n\");\n  }\n", 107);
    }
    if (show_ode == 1){
      sAppendN(&sbOut[show_ode],  "  (&_solveData->subjects[_cSub])->dadt_counter[0]++;\n}\n\n", 56);
    } else if (show_ode == 2){
      //sAppendN(&sbOut[show_ode], "  free(__ld_DDtStateVar__);\n");
      sAppendN(&sbOut[show_ode],  "  (&_solveData->subjects[_cSub])->jac_counter[0]++;\n", 52);
      sAppendN(&sbOut[show_ode],  "}\n", 2);
    } else if (show_ode == 3){
      for (i = 0; i < tb.nd; i++){
	retieve_var(tb.di[i], buf);
	sAppend(&sbOut[show_ode], "  __zzStateVar__[%d]=",i);
	for (k = 0; k < (int)strlen(buf); k++){
	  if (buf[k] == '.'){
	    sAppendN(&sbOut[show_ode], "_DoT_", 5);
	    if (!rx_syntax_allow_dots){
	      trans_syntax_error_report_fn(NODOT);
	    }
	  } else {
	    sPut(&sbOut[show_ode], buf[k]);
	  }
	}
	sAppendN(&sbOut[show_ode],  ";\n", 2);
      }
      sAppendN(&sbOut[show_ode],  "}\n", 2);
    } else {
      sAppendN(&sbOut[show_ode],  "\n", 1);
      for (i=0, j=0; i<tb.nv; i++) {
	if (tb.lh[i] != 1) continue;
	retieve_var(i, buf);
	sAppend(&sbOut[show_ode],  "  _lhs[%d]=", j);
	for (k = 0; k < (int)strlen(buf); k++){
	  if (buf[k] == '.'){
	    sAppendN(&sbOut[show_ode], "_DoT_", 5);
	    if (!rx_syntax_allow_dots){
	      trans_syntax_error_report_fn(NODOT);
	    }
	  } else {
	    sPut(&sbOut[show_ode], buf[k]);
	  }
	}
	sAppendN(&sbOut[show_ode],  ";\n", 2);
	j++;
      }
      sAppendN(&sbOut[show_ode],  "}\n", 2);
    }
  }
}
  
void reset (){
  // Reset sb/sbt string buffers
  sIni(&sb);
  sIni(&sbt);

  // Reset Arrays
  memset(tb.ss,		0, 64*MXSYM*sizeof(char));
  memset(tb.de,		0, 64*MXSYM*sizeof(char));
  memset(tb.deo,	0, MXSYM*sizeof(int));
  memset(tb.vo,		0, MXSYM*sizeof(int));
  memset(tb.lh,		0, MXSYM*sizeof(int));
  memset(tb.ini,	0, MXSYM*sizeof(int));
  memset(tb.di,		0, MXDER*sizeof(int));
  memset(tb.fdi,        0, MXDER*sizeof(int));
  memset(tb.dy,		0, MXSYM*sizeof(int));
  memset(tb.sdfdy,	0, MXSYM*sizeof(int));
  // Reset integers
  tb.nv		= 0;
  tb.ix		= 0;
  tb.id		= 0;
  tb.fn		= 0;
  tb.nd		= 0;
  tb.pos	= 0;
  tb.pos_de	= 0;
  tb.ini_i	= 0;
  tb.statei	= 0;
  tb.sensi	= 0;
  tb.li		= 0;
  tb.pi		= 0;
  tb.cdf	= 0;
  tb.ndfdy	= 0;
  tb.maxtheta   = 0;
  tb.maxeta     = 0;
  tb.fdn        = 0;
  // reset globals
  found_print = 0;
  found_jac = 0;
  rx_syntax_error = 0;
  rx_suppress_syntax_info=0;
  rx_podo = 0;
  memset(s_aux_info,         0, 64*MXSYM);
  rx_syntax_assign = 0;
  rx_syntax_star_pow = 0;
  rx_syntax_require_semicolon = 0;
  rx_syntax_allow_dots = 0;
  rx_syntax_allow_ini0 = 1;
  rx_syntax_allow_ini = 1;

  maxSumProdN = 0;
  SumProdLD = 0;  
}

void writeSb(sbuf *sbb, FILE *fp){
  // Adapted from ideas by Christian H
  // http://forums.codeguru.com/showthread.php?77477-What-is-the-fastest-way-to-write-data-to-a-file
  unsigned totalWritten=0;
  const unsigned OS_PAGESIZE = 4*1024;
  while( totalWritten < sbb->o) {
    register unsigned toWrite = min( OS_PAGESIZE, sbb->o - totalWritten);
    register unsigned written = fwrite(sbb->s + totalWritten, 1, toWrite, fp);
    if( toWrite != written){
      fclose(fp);
      error("IO error writing parsed C file.");
    } else{
      totalWritten += written; // add the written bytes
    }
  }
  if (totalWritten != sbb->o) {
    fclose(fp);
    error("IO error writing parsed C file.");
  }
}

char *gBuf;
void trans_internal(char* parse_file, int isStr){
  char buf1[512], buf2[512], bufe[2048];
  int i,j,found,islhs;
  D_ParseNode *pn;
  /* any number greater than sizeof(D_ParseNode_User) will do;
     below 1024 is used */
  D_Parser *p = new_D_Parser(&parser_tables_RxODE, 1024);
  p->save_parse_tree = 1;
  if (isStr){
    gBuf = parse_file;
  } else {
      gBuf = r_sbuf_read(parse_file);
      err_msg((intptr_t) gBuf, "error: empty buf for FILE_to_parse\n", -2);
  }
  sIni(&sbNrm);
  sIni(&sbPm);  
  if ((pn=dparse(p, gBuf, (int)strlen(gBuf))) && !p->syntax_errors) {
    wprint_parsetree(parser_tables_RxODE, pn, 0, wprint_node, NULL);
    // Determine Jacobian vs df/dvar
    for (i=0; i<tb.ndfdy; i++) {                     /* name state vars */
      retieve_var(tb.df[i], buf1);
      found=0;
      for (j=0; j<tb.nd; j++) {                     /* name state vars */
        retieve_var(tb.di[j], buf2);
	if (!strcmp(buf1, buf2)){
	  found=1;
          break;
	}
      }
      if (!found){
	retieve_var(tb.dy[i], buf2);
	sprintf(bufe,NOSTATE,buf1,buf2,buf1);
	trans_syntax_error_report_fn(bufe);
      }
      // Now the dy()
      retieve_var(tb.dy[i], buf1);
      found=0;
      for (j=0; j<tb.nd; j++) {                     /* name state vars */
        retieve_var(tb.di[j], buf2);
        if (!strcmp(buf1, buf2)){
          found=1;
          break;
        }
      }
      if (!found){
	for (j=0; j<tb.nv; j++) {
          islhs = tb.lh[j];
	  retieve_var(j, buf2);
          if (islhs>1) continue; /* is a state var */
          retieve_var(j, buf2);
          if ((islhs != 1 || tb.ini[j] == 1) &&!strcmp(buf1, buf2)){
	    found=1;
	    // This is a df(State)/dy(Parameter)
	    tb.sdfdy[i] = 1;
	    break;
	  }
        }
      }
      if (!found){
        retieve_var(tb.df[i], buf1);
      	retieve_var(tb.dy[i], buf2);
      	sprintf(bufe,NOSTATEVAR,buf1,buf2,buf2);
        trans_syntax_error_report_fn(bufe);
      }
    }
  } else {
    rx_syntax_error = 1;
  }
  free_D_Parser(p);
}

SEXP _RxODE_trans(SEXP parse_file, SEXP extra_c, SEXP prefix, SEXP model_md5, SEXP parseStr){
  char *in;
  char buf[1024], buf2[512], df[128], dy[128];
  char snum[512];
  char *s2;
  char sLine[MXLEN+1];
  int i, j, islhs, pi=0, li=0, ini_i = 0,o2=0,k=0, l=0;

  double d;
  int isStr =INTEGER(parseStr)[0];
  reset(); 
  rx_syntax_assign = R_get_option("RxODE.syntax.assign",1);
  rx_syntax_star_pow = R_get_option("RxODE.syntax.star.pow",1);
  rx_syntax_require_semicolon = R_get_option("RxODE.syntax.require.semicolon",0);
  rx_syntax_allow_dots = R_get_option("RxODE.syntax.allow.dots",1);
  rx_suppress_syntax_info = R_get_option("RxODE.suppress.syntax.info",0);
  rx_syntax_allow_ini0 = R_get_option("RxODE.syntax.allow.ini0",1);
  rx_syntax_allow_ini  = R_get_option("RxODE.syntax.allow.ini",1);
  rx_syntax_allow_assign_state = R_get_option("RxODE.syntax.assign.state",0);
  rx_syntax_error = 0;
  set_d_use_r_headers(0);
  set_d_rdebug_grammar_level(0);
  set_d_verbose_level(0);
  rx_podo = 0;
  if (isString(extra_c) && length(extra_c) == 1){
    in = r_dup_str(CHAR(STRING_ELT(extra_c,0)),0);
    extra_buf = r_sbuf_read(in);
    if (!((intptr_t) extra_buf)){ 
      extra_buf = (char *) R_alloc(1,sizeof(char));
      extra_buf[0]='\0';
    }
  } else {
    extra_buf = (char *) R_alloc(1,sizeof(char));
    extra_buf[0] = '\0';
  }

  /* orig = r_dup_str(CHAR(STRING_ELT(orig_file,0)),0); */
  in = r_dup_str(CHAR(STRING_ELT(parse_file,0)),0);
  
  if (isString(prefix) && length(prefix) == 1){
    model_prefix = r_dup_str(CHAR(STRING_ELT(prefix,0)),0);
  } else {
    error("model prefix must be specified");
  }

  if (isString(model_md5) && length(model_md5) == 1){
    md5 = r_dup_str(CHAR(STRING_ELT(model_md5,0)),0);
  } else {
    md5 = R_alloc(1,sizeof(char));
    md5[0] = '\0';
  }
  
  trans_internal(in, isStr);
  
  for (i=0; i<tb.nv; i++) {
    islhs = tb.lh[i];
    if (islhs>1) continue;      /* is a state var */
    if (islhs == 1){
      li++;
    } else {
      pi++;
    }
  }
  tb.pi=pi;
  tb.li=li;
  
  int pro = 0;
  SEXP lst   = PROTECT(allocVector(VECSXP, 11));pro++;
  SEXP names = PROTECT(allocVector(STRSXP, 11));pro++;
  
  SEXP tran  = PROTECT(allocVector(STRSXP, 19));pro++;
  SEXP trann = PROTECT(allocVector(STRSXP, 19));pro++;

  SEXP state    = PROTECT(allocVector(STRSXP,tb.statei));pro++;
  SEXP stateRmS = PROTECT(allocVector(INTSXP,tb.statei));pro++;
  int *stateRm  = INTEGER(stateRmS);
  SEXP sens     = PROTECT(allocVector(STRSXP,tb.sensi));pro++;
  SEXP fn_ini   = PROTECT(allocVector(STRSXP, tb.fdn));pro++;

  SEXP dfdy = PROTECT(allocVector(STRSXP,tb.ndfdy));pro++;
  
  SEXP params = PROTECT(allocVector(STRSXP, tb.pi));pro++;
  SEXP lhs    = PROTECT(allocVector(STRSXP, tb.li));pro++;

  SEXP inin   = PROTECT(allocVector(STRSXP, tb.ini_i));pro++;
  SEXP ini    = PROTECT(allocVector(REALSXP, tb.ini_i));pro++;

  SEXP model  = PROTECT(allocVector(STRSXP,1));pro++;
  SEXP modeln = PROTECT(allocVector(STRSXP,1));pro++;
  k=0;j=0;l=0;
  for (i=0; i<tb.nd; i++) {                     /* name state vars */
    retieve_var(tb.di[i], buf);
    if (strstr(buf,"rx__sens_")){
      SET_STRING_ELT(sens,j++,mkChar(buf));
      SET_STRING_ELT(state,k++,mkChar(buf));
      stateRm[k-1]=tb.idi[i];
    } else {
      SET_STRING_ELT(state,k++,mkChar(buf));
      stateRm[k-1]=tb.idi[i];
    }
    if (tb.fdi[i]){
      SET_STRING_ELT(fn_ini,l++,mkChar(buf));
    }
  }
  for (i=0; i<tb.ndfdy; i++) {                     /* name state vars */
    retieve_var(tb.df[i], df);
    retieve_var(tb.dy[i], dy);
    for (j = 1; j <= tb.maxtheta;j++){
      sprintf(buf,"_THETA_%d_",j);
      if (!strcmp(dy,buf2)){
        sprintf(dy,"THETA[%d]",j);
      }
    }
    for (j = 1; j <= tb.maxeta;j++){
      sprintf(buf,"_ETA_%d_",j);
      if (!strcmp(dy,buf)){
        sprintf(dy,"ETA[%d]",j);
      }
    }
    sprintf(buf,"df(%s)/dy(%s)",df,dy);
    SET_STRING_ELT(dfdy,i,mkChar(buf));
  }
  li=0, pi=0;
  for (i=0; i<tb.nv; i++) {
    islhs = tb.lh[i];
    if (islhs>1) continue;      /* is a state var */
    retieve_var(i, buf);
    if (islhs == 1){
      SET_STRING_ELT(lhs,li++,mkChar(buf));
    } else {
      for (j = 1; j <= tb.maxtheta;j++){
	sprintf(buf2,"_THETA_%d_",j);
	if (!strcmp(buf, buf2)){
	  sprintf(buf,"THETA[%d]",j);
	}
      }
      for (j = 1; j <= tb.maxeta;j++){
        sprintf(buf2,"_ETA_%d_",j);
        if (!strcmp(buf, buf2)){
          sprintf(buf,"ETA[%d]",j);
        }
      }
      SET_STRING_ELT(params,pi++,mkChar(buf));
    }
  }
  SET_STRING_ELT(names,0,mkChar("trans"));
  SET_VECTOR_ELT(lst,  0,tran);
  
  SET_STRING_ELT(names,3,mkChar("params"));
  SET_VECTOR_ELT(lst,  3,params);

  SET_STRING_ELT(names,1,mkChar("lhs"));
  SET_VECTOR_ELT(lst,  1,lhs);

  SET_STRING_ELT(names,2,mkChar("state"));
  SET_VECTOR_ELT(lst,  2,state);
  
  SET_STRING_ELT(names,4,mkChar("ini"));
  SET_VECTOR_ELT(lst,  4,ini);

  SET_STRING_ELT(names,5,mkChar("model"));
  SET_VECTOR_ELT(lst,  5,model);

  SET_STRING_ELT(names,6,mkChar("podo"));
  SET_VECTOR_ELT(lst,  6,ScalarLogical(rx_podo));

  SET_STRING_ELT(names,7,mkChar("dfdy"));
  SET_VECTOR_ELT(lst,  7,dfdy);

  SET_STRING_ELT(names,8,mkChar("sens"));
  SET_VECTOR_ELT(lst,  8,sens);
  
  SET_STRING_ELT(names,8,mkChar("fn.ini"));
  SET_VECTOR_ELT(lst,  8,fn_ini);

  SET_STRING_ELT(names,9,mkChar("state.ignore"));
  SET_VECTOR_ELT(lst,  9,stateRmS);
  
  SET_STRING_ELT(trann,0,mkChar("jac"));
  if (found_jac == 1){
    SET_STRING_ELT(tran,0,mkChar("fulluser")); // Full User Matrix
  } else {
    SET_STRING_ELT(tran,0,mkChar("fullint")); // Full Internal Matrix
  }
  SET_STRING_ELT(trann,1,mkChar("prefix"));
  SET_STRING_ELT(tran,1,mkChar(buf));

  sprintf(buf,"%sdydt",model_prefix);
  SET_STRING_ELT(trann,2,mkChar("dydt"));
  SET_STRING_ELT(tran,2,mkChar(buf)) ;

  sprintf(buf,"%scalc_jac",model_prefix);
  SET_STRING_ELT(trann,3,mkChar("calc_jac"));
  SET_STRING_ELT(tran, 3,mkChar(buf));

  sprintf(buf,"%scalc_lhs",model_prefix);
  SET_STRING_ELT(trann,4,mkChar("calc_lhs"));
  SET_STRING_ELT(tran, 4,mkChar(buf));

  sprintf(buf,"%smodel_vars",model_prefix);
  SET_STRING_ELT(trann,5,mkChar("model_vars"));
  SET_STRING_ELT(tran, 5,mkChar(buf));

  sprintf(buf,"%sode_solver",model_prefix);
  SET_STRING_ELT(trann,6,mkChar("ode_solver"));
  SET_STRING_ELT(tran, 6,mkChar(buf));

  sprintf(buf,"%sode_solver_sexp",model_prefix);
  SET_STRING_ELT(trann,7,mkChar("ode_solver_sexp"));
  SET_STRING_ELT(tran, 7,mkChar(buf));
  
  sprintf(buf,"%sode_solver_focei_eta",model_prefix);
  SET_STRING_ELT(trann,8,mkChar("ode_solver_focei_eta"));
  SET_STRING_ELT(tran, 8,mkChar(buf));

  sprintf(buf,"%sode_solver_ptr",model_prefix);
  SET_STRING_ELT(trann,9,mkChar("ode_solver_ptr"));
  SET_STRING_ELT(tran, 9,mkChar(buf));

  sprintf(buf,"%sinis",model_prefix);
  SET_STRING_ELT(trann,10,mkChar("inis"));
  SET_STRING_ELT(tran, 10,mkChar(buf));

  sprintf(buf,"%sode_solver_xptr",model_prefix);
  SET_STRING_ELT(trann,11,mkChar("ode_solver_xptr"));
  SET_STRING_ELT(tran, 11,mkChar(buf));

  sprintf(buf,"%sdydt_lsoda",model_prefix);
  SET_STRING_ELT(trann,12,mkChar("dydt_lsoda"));
  SET_STRING_ELT(tran, 12,mkChar(buf));

  sprintf(buf,"%scalc_jac_lsoda",model_prefix);
  SET_STRING_ELT(trann,13,mkChar("calc_jac_lsoda"));
  SET_STRING_ELT(tran, 13,mkChar(buf));

  sprintf(buf,"%sode_solver_solvedata",model_prefix);
  SET_STRING_ELT(trann,14,mkChar("ode_solver_solvedata"));
  SET_STRING_ELT(tran, 14,mkChar(buf));

  sprintf(buf,"%sode_solver_get_solvedata",model_prefix);
  SET_STRING_ELT(trann,15,mkChar("ode_solver_get_solvedata"));
  SET_STRING_ELT(tran, 15,mkChar(buf));

  sprintf(buf,"%d",tb.statei);
  SET_STRING_ELT(trann,16,mkChar("neq"));
  SET_STRING_ELT(tran, 16,mkChar(buf));

  sprintf(buf,"%d", tb.li);
  SET_STRING_ELT(trann,17,mkChar("nlhs"));
  SET_STRING_ELT(tran, 17,mkChar(buf));

  sprintf(buf,"%sdydt_liblsoda",model_prefix);
  SET_STRING_ELT(trann,18,mkChar("dydt_liblsoda"));
  SET_STRING_ELT(tran, 18,mkChar(buf));

  sbPm.o=0;
  while(sgets(sLine, MXLEN, &sbPm)) {
    s2 = strstr(sLine,"(__0__)");
    if (s2){
      // See if this is a reclaimed initilization variable.
      for (i=0; i<tb.nv; i++) {
        if (tb.ini[i] == 1 && tb.lh[i] != 1){
          //(__0__)V2=
          retieve_var(i, buf);
          sprintf(buf2,"(__0__)");
          o2 = 7;
          for (k = 0; k < (int)strlen(buf); k++){
            if (buf[k] == '.'){
              sprintf(buf2+o2,"_DoT_");
              if (!rx_syntax_allow_dots){
                trans_syntax_error_report_fn(NODOT);
              }
              o2+=5;
            } else {
              sprintf(buf2+o2,"%c",buf[k]);
              o2++;
            }
          }
          sprintf(buf2+o2,"=");
          s2 = strstr(sLine,buf2);
          if (s2){
            /* Rprintf("%s[%d]->\n",buf,ini_i++); */
            SET_STRING_ELT(inin,ini_i,mkChar(buf));
            sprintf(snum,"%.*s",(int)(strlen(sLine)-strlen(buf2) - 2), sLine + strlen(buf2));
            sscanf(snum, "%lf", &d);
            REAL(ini)[ini_i++] = d;
            continue;
          }
        }
      }
      continue;
    }
  }
  // putin constants
  for (i=0; i<tb.nv; i++) {
    if (tb.ini[i] == 0 && tb.lh[i] != 1) {
      retieve_var(i, buf);
      // Put in constants
      if  (!strcmp("pi",buf)){
        SET_STRING_ELT(inin,ini_i,mkChar("pi"));
        REAL(ini)[ini_i++] = M_PI;
      }
    }
  }
  
  SET_STRING_ELT(modeln,0,mkChar("normModel"));
  SET_STRING_ELT(model,0,mkChar(sbNrm.s));
  
  setAttrib(ini,   R_NamesSymbol, inin);
  setAttrib(tran,  R_NamesSymbol, trann);
  setAttrib(lst,   R_NamesSymbol, names);
  setAttrib(model, R_NamesSymbol, modeln);
  SEXP cls = PROTECT(allocVector(STRSXP, 1));pro++;
  SET_STRING_ELT(cls, 0, mkChar("rxModelVars"));
  classgets(lst, cls);
  
  UNPROTECT(pro);
  if (rx_syntax_error){
    error("Syntax Errors (see above)");
  }
  /* Free(sbPm); Free(sbNrm); */
  return lst;
}

SEXP _RxODE_codegen(SEXP c_file, SEXP prefix, SEXP libname,
		    SEXP pMd5, SEXP timeId, SEXP fixInis){
  sbuf *sbb;
  if (!sbPm.o || !sbNrm.o){
    error("Nothing in output queue to write.");
  }
  if (!isString(c_file) || length(c_file) != 1){
    error("c_file should only be 1 file");
  }
  fpIO = fopen(CHAR(STRING_ELT(c_file,0)), "wb");
  err_msg((intptr_t) fpIO, "error opening output c file\n", -2);
  for (int i = 0; i < 5; i++){
    codegen(gBuf, i, CHAR(STRING_ELT(prefix,0)),
	    CHAR(STRING_ELT(libname, 0)),
	    CHAR(STRING_ELT(pMd5,0)),
	    CHAR(STRING_ELT(timeId, 0)),
	    CHAR(STRING_ELT(fixInis, 0)));
  }  
  writeSb(&sbOut[1], fpIO);
  writeSb(&sbOut[2], fpIO);
  writeSb(&sbOut[3], fpIO);
  writeSb(&sbOut[0], fpIO);
  writeSb(&sbOut[4], fpIO);
  fclose(fpIO);
  for (int i = 0; i < 5; i++){
    sbb = &sbOut[i];
    sbb->o = 0;
  }
  sbPm.o =0;
  sbNrm.o=0;
  return R_NilValue;
}


