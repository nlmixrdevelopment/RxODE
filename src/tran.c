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

#define STRINGIFY(...) STRINGIFY_AUX(__VA_ARGS__)
#define STRINGIFY_AUX(...) #__VA_ARGS__

#define gCode(i) (&sbOut)->s[0]='\0';		\
  (&sbOut)->o=0;				\
  codegen(gBuf, i, CHAR(STRING_ELT(prefix,0)),	\
	  CHAR(STRING_ELT(libname, 0)),		\
	  CHAR(STRING_ELT(pMd5,0)),		\
	  CHAR(STRING_ELT(timeId, 0)),		\
	  CHAR(STRING_ELT(fixInis, 0)),         \
	  CHAR(STRING_ELT(fixInis, 1)),         \
	  CHAR(STRING_ELT(libname, 1)));					\
  writeSb(&sbOut, fpIO);

#define aAppendN(str, len) sAppendN(&sb, str, len); sAppendN(&sbDt, str, len);
#define aProp(prop) curLineProp(&sbPm, prop); curLineProp(&sbPmDt, prop);
#define aType(type) curLineType(&sbPm, type); curLineType(&sbPmDt, type); 

#define FBIO 1
#define ALAG 2
#define RATE 3
#define DUR 4
#define TINI 5
#define TLOGIC 6
#define PODE0 7
#define PJAC 8
#define PJAC0 9
#define PODE 10
#define PPRN 11
#define TDDT 12
#define TJAC 13
#define TF0 14
#define PLHS 15
#define PFPRN 16
#define TASSIGN 17
#define TMTIME 18

#define NOASSIGN "'<-' not supported, use '=' instead or set 'options(RxODE.syntax.assign = TRUE)'"
#define NEEDSEMI "Lines need to end with ';'\n     To match R's handling of line endings set 'options(RxODE.syntax.require.semicolon = FALSE)'"
#define NEEDPOW "'**' not supported, use '^' instead or set 'options(RxODE.syntax.star.pow = TRUE)'"
#define NODOT "'.' in variables and states not supported, use '_' instead or set 'options(RxODE.syntax.allow.dots = TRUE)'"
#define NOINI0 "'%s(0)' for initialization not allowed.  To allow set 'options(RxODE.syntax.allow.ini0 = TRUE)'"
#define NOSTATE "Defined 'df(%s)/dy(%s)', but '%s' is not a state"
#define NOSTATEVAR "Defined 'df(%s)/dy(%s)', but '%s' is not a state or variable"
#define ODEFIRST "ODEs compartment 'd/dt(%s)' must be defined before changing its properties (f/alag/rate/dur).\nIf you want to change this set 'options(RxODE.syntax.require.ode.first = FALSE).\nBe warned this will RxODE numbers compartments based on first occurance of property or ODE"

#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#if (__STDC_VERSION__ >= 199901L)
#include <stdint.h>
#endif

void setInits(SEXP init);
int getInits(char *s_aux_info, int *o);

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

int syntaxErrorExtra = 0;
int isEsc=0;
const char *lastStr;
int lastStrLoc=0;
// Taken from dparser and changed to use Calloc
char * rc_dup_str(const char *s, const char *e) {
  lastStr=s;
  int l = e ? e-s : (int)strlen(s);
  syntaxErrorExtra=min(l-1, 40);
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

int rx_syntax_error = 0, rx_suppress_syntax_info=0, rx_podo = 0, rx_syntax_require_ode_first = 1;

extern D_ParserTables parser_tables_RxODE;

unsigned int found_jac = 0, nmtime=0;
int rx_syntax_assign = 0, rx_syntax_star_pow = 0,
  rx_syntax_require_semicolon = 0, rx_syntax_allow_dots = 0,
  rx_syntax_allow_ini0 = 1, rx_syntax_allow_ini = 1, rx_syntax_allow_assign_state = 0,
  maxSumProdN = 0, SumProdLD = 0, good_jac=1, extraCmt=0;

char s_aux_info[64*MXSYM*4];


typedef struct symtab {
  char ss[64*MXSYM];                     /* symbol string: all vars*/
  char de[64*MXSYM];             /* symbol string: all Des*/
  char ddt[MXSYM];
  int deo[MXSYM];        /* offest of des */
  int vo[MXSYM];        /* offset of symbols */
  int lh[MXSYM];        /* lhs symbols? =9 if a state var*/
  int ini[MXSYM];        /* initial variable assignment =2 if there are two assignments */
  int mtime[MXSYM];
  double iniv[MXSYM];        /* Initial values */
  int ini0[MXSYM];        /* state initial variable assignment =2 if there are two assignments */
  int di[MXDER];        /* ith of state vars */
  int idi[MXDER];       /* should ith state variable be ignored 0/1 */
  int idu[MXDER];       /* Has the ith state been used in a derivative expression? */
  int fdi[MXDER];        /* Functional initialization of state variable */
  int dvid[MXDER];
  int dvidn;
  int nv;                       /* nbr of symbols */
  int ix;                       /* ith of curr symbol */
  int id;                       /* ith of curr symbol */
  int fn;                       /* curr symbol a fn?*/
  int nd;                       /* nbr of dydt */
  int pos;
  int pos_de;
  int ini_i; // #ini
  int statei; // # states
  int nExtra;
  int fdn; // # conditional states
  int sensi;
  int li; // # lhs
  int pi; // # param
  int isPi; // # pi?
  int linCmt; // Unparsed linear compartment
  // Save Jacobian information
  int df[MXSYM];
  int dy[MXSYM];
  int sdfdy[MXSYM];
  int cdf;
  int ndfdy;
  int maxtheta;
  int maxeta;
  int hasDepot;
  int hasCentral;
  int hasDepotCmt;
  int hasCentralCmt;
  int hasKa;
} symtab;
symtab tb;

typedef struct sbuf {
  char *s;        /* curr print buffer */
  int sN;
  int o;                        /* offset of print buffer */
} sbuf;

sbuf sb, sbDt;                        /* buffer w/ current parsed & translated line */
sbuf sbt;

void sIniTo(sbuf *sbb, int to){
  if (sbb->sN <= to) {
    if (sbb->sN <= 0){
      sbb->s = Calloc(to, char);
      sbb->sN = to;
      sbb->s[0]='\0';
      sbb->o=0;
    } else {
      Free(sbb->s);
      sbb->s = Calloc(to, char);
      sbb->sN = to;
      sbb->s[0]='\0';
      sbb->o=0;
    }
  } else {
    sbb->s[0]='\0';
    sbb->o=0;
  }
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

void sFree(sbuf *sbb){
  Free(sbb->s);
  sbb->sN=0;
  sbb->o=0;
}

void sFreeIni(sbuf *sbb){
  sFree(sbb);
  sIni(sbb);
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

typedef struct vLines {
  char *s;
  int sN;
  int o;
  int n;
  int nL;
  char **line;
  int *lProp;
  int *lType;
} vLines;

void lineIni(vLines *sbb){
  if (sbb->sN <= 0){
    sbb->s = Calloc(MXBUF, char);
    sbb->sN = MXBUF;
    sbb->s[0]='\0';
    sbb->o=0;
  } else {
    sbb->s[0]='\0';
    sbb->o=0;
  }
  if (sbb->nL < 1000){
    Free(sbb->lProp);
    Free(sbb->line);
    Free(sbb->lType);
    sbb->lProp = Calloc(1000, int);
    sbb->lType = Calloc(1000, int);
    sbb->line = Calloc(1000, char*);
    sbb->nL=1000;
  }
  sbb->lProp[0] = -1;
  sbb->lType[0] = 0;
  sbb->n = 0;
  sbb->o=0;
}

void lineFree(vLines *sbb){
  Free(sbb->s);
  Free(sbb->lProp);
  Free(sbb->lType);
  Free(sbb->line);
  sbb->sN=0;
  sbb->nL=0;
}

void addLine(vLines *sbb, const char *format, ...){
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
  if (sbb->n + 1 >= sbb->nL){
    int mx = sbb->n + 1000;
    sbb->lProp = Realloc(sbb->lProp, mx, int);
    sbb->lType = Realloc(sbb->lType, mx, int);
    sbb->line = Realloc(sbb->line, mx, char*);
  }
  sbb->line[sbb->n]=&(sbb->s[sbb->o]);
  sbb->o +=n+1; // Add the \0 at the end.
  sbb->n = sbb->n+1;
  sbb->lProp[sbb->n] = -1;
  sbb->lType[sbb->n] = 0;  

}

void curLineProp(vLines *sbb, int propId){
  sbb->lProp[sbb->n] = propId;
}

void curLineType(vLines *sbb, int propId){
  sbb->lType[sbb->n] = propId;
}


vLines sbPm, sbPmDt;
sbuf sbNrm;

char *extra_buf, *model_prefix, *md5 = NULL;
int foundF=0,foundLag=0, foundRate=0, foundDur=0, foundF0=0, needSort=0;

sbuf sbOut;

static FILE *fpIO;

int lastSyntaxErrorLine=0;
static void trans_syntax_error_report_fn(char *err);
static void trans_syntax_error_report_fn0(char *err);
char *getLine (char *src, int line, int *lloc);
void updateSyntaxCol();

/* new symbol? if no, find it's ith */
int new_or_ith(const char *s) {
  int i, len, len_s=(int)strlen(s);

  if (tb.fn) return 0;
  if (!strcmp("t", s)) return 0;
  if (!strcmp("rate", s)){
    updateSyntaxCol();
    trans_syntax_error_report_fn("'rate' cannot be a variable in an RxODE model");
    return 0;
  }
  if (!strcmp("dur", s)){
    updateSyntaxCol();
    trans_syntax_error_report_fn("'dur' cannot be a variable in an RxODE model");
    return 0;
  }
  if (!strcmp("amt", s)){
    updateSyntaxCol();
    trans_syntax_error_report_fn("'amt' cannot be a variable in an RxODE model");
    return 0;
  }
  if (!strcmp("ss", s)){
    updateSyntaxCol();
    trans_syntax_error_report_fn("'ss' cannot be a variable in an RxODE model");
    return 0;
  }
  if (!strcmp("addl", s)){
    updateSyntaxCol();
    trans_syntax_error_report_fn("'addl' cannot be a variable in an RxODE model");
    return 0;
  }
  if (!strcmp("evid", s)){
    updateSyntaxCol();
    trans_syntax_error_report_fn("'evid' cannot be a variable in an RxODE model");
    return 0;
  }
  if (!strcmp("ii", s)){
    updateSyntaxCol();
    trans_syntax_error_report_fn("'ii' cannot be a variable in an RxODE model");
    return 0;
  }
  if (!strcmp("dvid", s)){
    updateSyntaxCol();
    trans_syntax_error_report_fn("'dvid' cannot be a variable in an RxODE model");
    return 0;
  }
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
  if (!strcmp("pi", s)) tb.isPi=1;
  if (!tb.hasKa && !strcmp("ka", s)) tb.hasKa=1;
  if (!tb.hasKa && !strcmp("Ka", s)) tb.hasKa=1;
  if (!tb.hasKa && !strcmp("KA", s)) tb.hasKa=1;
  if (!tb.hasKa && !strcmp("kA", s)) tb.hasKa=1;
  if (!strcmp("newind", s)) return 0;
  if (!strcmp("NEWIND", s)) return 0;
  // Ignore THETA[] and ETA
  if (strstr("[", s) != NULL) return 0;

  for (i=0; i<tb.nv; i++) {
    len = tb.vo[i+1] - tb.vo[i] - 1;  /* -1 for added ',' */
    if (!strncmp(tb.ss+tb.vo[i], s, max(len, len_s))) { /* note we need take the max in order not to match a sub-string */
      tb.ix = i;
      return 0;
    }
  }
  return 1;
}

typedef struct nodeInfo {
  int alag;
  int assignment;
  int constant;
  int der_rhs;
  int derivative;
  int dfdy;
  int dfdy_rhs;
  int dur;
  int end_statement;
  int eta;
  int factorial;
  int factorial_exp;
  int fbio;
  int function;
  int identifier;
  int identifier_r;
  int identifier_r_no_output;
  int ini0;
  int ini0f;
  int ini;
  int jac;
  int jac_rhs;
  int lfactorial;
  int lfactorial_exp;
  int max;
  int min;
  int mtime;
  int mult_part;
  int power_expression;
  /* int print_command; */
  int printf_statement;
  int prod;
  int rate;
  int selection_statement;
  int selection_statement__8;
  int sign;
  int sum;
  int theta0;
  int theta0_noout;
  int theta;
  int transit2;
  int transit3;
  int cmt_statement;
  int dvid_statementI;
} nodeInfo;

#define NIB(what) ni.what
#define nodeHas(what) (NIB(what) == -1 ? (NIB(what) = !strcmp(STRINGIFY(what), name)) : NIB(what))
//#define nodeHas(what) (!strcmp(STRINGIFY(what), name))

void niReset(nodeInfo *ni){
  ni->mtime = -1;
  ni->alag = -1;
  ni->assignment = -1;
  ni->constant = -1;
  ni->der_rhs = -1;
  ni->derivative = -1;
  ni->dfdy = -1;
  ni->dfdy_rhs = -1;
  ni->dur = -1;
  ni->end_statement = -1;
  ni->eta = -1;
  ni->factorial = -1;
  ni->factorial_exp = -1;
  ni->fbio = -1;
  ni->function = -1;
  ni->identifier = -1;
  ni->identifier_r = -1;
  ni->identifier_r_no_output = -1;
  ni->ini = -1;
  ni->ini0 = -1;
  ni->ini0f = -1;
  ni->jac = -1;
  ni->jac_rhs = -1;
  ni->lfactorial = -1;
  ni->lfactorial_exp = -1;
  ni->max = -1;
  ni->min = -1;
  ni->mult_part = -1;
  ni->power_expression = -1;
  /* ni->print_command = -1; */
  ni->printf_statement = -1;
  ni->prod = -1;
  ni->rate = -1;
  ni->rate = -1;
  ni->selection_statement = -1;
  ni->selection_statement__8 = -1;
  ni->sign = -1;
  ni->sum = -1;
  ni->theta = -1;
  ni->theta0 = -1;
  ni->theta0_noout = -1;
  ni->transit2 = -1;
  ni->transit3 = -1;
  ni->cmt_statement = -1;
  ni->dvid_statementI = -1;
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
  nodeInfo ni;
  niReset(&ni);
  if (!strcmp("time",value)){
    aAppendN("t", 1);
    sAppendN(&sbt, "t", 1);
  } else if (!strcmp("podo",value)){
    aAppendN("_solveData->subjects[_cSub].podo", 32);
    sAppendN(&sbt, "podo", 4);
    rx_podo = 1;
  } else if (!strcmp("CMT",value)){
    aAppendN("_CMT", 4);
    sAppendN(&sbt, "CMT", 3);
  } else if (!strcmp("tlast",value)){
    aAppendN("_solveData->subjects[_cSub].tlast", 33);
    sAppendN(&sbt, "tlast", 5);
  } else if (!strcmp("rx__PTR__",value)){
    aAppendN("_solveData, _cSub", 17);
    sAppendN(&sbt, "rx__PTR__", 9);
  } else if (nodeHas(identifier) && !strcmp("gamma",value)){
    aAppendN("lgammafn", 8);
    sAppendN(&sbt, "lgammafn", 8);
  } else if (nodeHas(identifier) && !strcmp("lfactorial",value)){
    aAppendN("lgamma1p", 8);
    sAppendN(&sbt, "lgamma1p", 8);
  } else if (nodeHas(identifier) && !strcmp("log",value)){
    aAppendN("_safe_log", 9);
    sAppendN(&sbt, "log", 3);
  } else if (nodeHas(identifier) && !strcmp("abs",value)){
    aAppendN("fabs", 4);
    sAppendN(&sbt,"abs", 3);
  } else if (nodeHas(identifier) && !strcmp("linCmt",value)) {
    aAppendN("linCmt", 6);
    sAppendN(&sbt,"linCmt", 6);
    tb.linCmt=1;
  } else if (nodeHas(identifier) && !strcmp("solveLinB",value)){
    aAppendN("solveLinB", 9);
    sAppendN(&sbt,"solveLinB", 9);
    tb.linCmt=2;
  } else {
    // Apply fix for dot.syntax
    for (i = 0; i < (int)strlen(value); i++){
      if (value[i] == '.' && nodeHas(identifier_r)){
	aAppendN("_DoT_", 5);
	sAppendN(&sbt, ".", 1);
        if (rx_syntax_allow_dots == 0){
	  updateSyntaxCol();
          trans_syntax_error_report_fn(NODOT);
        }
      } else {
	sPut(&sb, value[i]);
	sPut(&sbDt, value[i]);
	sPut(&sbt, value[i]);
      }
    }
  }
}
char *gBuf;
int gBufLast;
D_Parser *curP;
void wprint_parsetree(D_ParserTables pt, D_ParseNode *pn, int depth, print_node_fn_t fn, void *client_data) {
  char *name = (char*)pt.symbols[pn->symbol].name;
  nodeInfo ni;
  niReset(&ni);
  int nch = d_get_number_of_children(pn), i, k, ii, found, safe_zero = 0;
  char *value = (char*)rc_dup_str(pn->start_loc.s, pn->end);
  char buf[1024];
  char buf2[1024];
  buf[0]='\0';
  buf2[0]='\0';
  double d;
  if ((nodeHas(identifier) || nodeHas(identifier_r) ||
       nodeHas(identifier_r_no_output)  ||
       nodeHas(theta0_noout) || 
       nodeHas(theta0)) &&
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
    sPut(&sbDt, name[0]);
    if (!(strcmp(",", name)) && depth == 1){
      aAppendN("(double)", 8);
    }
    sPut(&sbt, name[0]);
  }
  if (nodeHas(identifier) ||
      nodeHas(identifier_r) ||
      nodeHas(constant) ||
      nodeHas(theta0) ||
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
    aAppendN(" =", 2);
    sAppendN(&sbt, "=", 1);
  } else if (!strcmp("~",name)){
    // Suppress LHS calculation with ~
    aAppendN(" =", 2);
    sAppendN(&sbt, "~", 1);
    tb.lh[tb.ix] = 10; // Suppress LHS printout.
  } else if (!strcmp("|",name)){
    aAppendN(" ||", 3);
    sAppendN(&sbt, "||", 2);
  } else if (!strcmp("&",name)){
    aAppendN(" &&", 3);
    sAppendN(&sbt, "&&", 2);
  }

  Free(value);

  //depth++;
  if (nch != 0) {
    if (nodeHas(power_expression)) {
      aAppendN(" R_pow(", 7);
    }
    for (i = 0; i < nch; i++) {
      if (!rx_syntax_assign  &&
          ((i == 4 && nodeHas(derivative)) ||
           (i == 6 && nodeHas(dfdy)))) {
        D_ParseNode *xpn = d_get_child(pn,i);
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	if (!strcmp("<-",v)){
	  updateSyntaxCol();
          trans_syntax_error_report_fn(NOASSIGN);
        }
        Free(v);
        continue;
      }
      if ((i == 3 || i == 4 || i < 2) &&
	  (nodeHas(derivative) ||nodeHas(fbio) || nodeHas(alag) ||
	   nodeHas(rate) || nodeHas(dur))) continue;
      
      if ((i == 3 || i < 2) && nodeHas(der_rhs)) continue;
      

      if (nodeHas(dfdy)     && i< 2)   continue;
      if (nodeHas(dfdy_rhs) && i< 2)   continue;
      if (nodeHas(dfdy)     && i == 3) continue;
      if (nodeHas(dfdy_rhs) && i == 3) continue;
      if (nodeHas(dfdy)     && i == 5) continue;
      if (nodeHas(dfdy_rhs) && i == 5) continue;
      
      if (nodeHas(dfdy)     && i == 6) continue;
      if (nodeHas(ini0)     && i == 1) continue;

      if (nodeHas(transit2) && i == 1) continue;
      if (nodeHas(transit3) && i == 1) continue;


      if ((nodeHas(theta) || nodeHas(eta)) && i != 2) continue;
      if (nodeHas(mtime) && (i == 0 || i == 1 || i == 3)) continue;
      if (nodeHas(cmt_statement) && (i == 0 || i == 1 || i == 3)) continue;
      if (i != 0 && nodeHas(dvid_statementI)) continue;
      tb.fn = (nodeHas(function) && i==0) ? 1 : 0;

      if (tb.fn) depth = 0;

      D_ParseNode *xpn = d_get_child(pn,i);

      if (nodeHas(dvid_statementI)){
	if (tb.dvidn == 0){
	  // dvid->cmt translation
	  sb.o=0;sbDt.o=0; sbt.o=0;
	  xpn = d_get_child(pn,2);
	  char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	  tb.dvid[0]=atoi(v);
	  if (tb.dvid[0] == 0) error("dvid() cannot have zeros in it");
	  sAppend(&sbt, "dvid(%d", tb.dvid[0]);
	  xpn = d_get_child(pn,3);
	  tb.dvidn = d_get_number_of_children(xpn)+1;
	  D_ParseNode *xpn2;
	  for (i = 0; i < tb.dvidn-1; i++){
	    xpn2 = d_get_child(xpn, i);
	    v = (char*)rc_dup_str(xpn2->start_loc.s, xpn2->end);
	    tb.dvid[i+1]=atoi(v+1);
	    if (tb.dvid[i+1] == 0) error("dvid() cannot have zeros in it");
	    sAppend(&sbt, ",%d", tb.dvid[i+1]);
	    Free(v);
	  }
	  sAppend(&sbNrm, "%s);\n", sbt.s);
	  Free(v);
	  continue;
	} else {
	  error("RxODE only supports one dvid() statement per model");
	}
	continue;
      }
      if (tb.fn){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (!strcmp("prod",v) || !strcmp("sum",v) || !strcmp("sign",v) ||
	    !strcmp("max",v) || !strcmp("min",v)){
	  ii = d_get_number_of_children(d_get_child(pn,3))+1;
	  if (!strcmp("prod", v)){
            sAppend(&sb, "_prod(_p, _input, _prodType(), %d, (double) ", ii);
	    sAppend(&sbDt, "_prod(_p, _input, _prodType(), %d, (double) ", ii);
            if (maxSumProdN < ii){
              maxSumProdN = ii;
            }
          } else if (!strcmp("sum", v)){
	    sAppend(&sb, "_sum(_p, _pld, -__MAX_PROD__, _sumType(), %d, (double) ", ii);
	    sAppend(&sbDt, "_sum(_p, _pld, -__MAX_PROD__, _sumType(), %d, (double) ", ii);
            if (SumProdLD < ii){
              SumProdLD = ii;
            }
	  } else {
	    sAppend(&sb, "_%s(%d, (double) ", v, ii);
	    sAppend(&sbDt, "_%s(%d, (double) ", v, ii);
	  }
	  sAppend(&sbt, "%s(", v);
          Free(v);
          i = 1;// Parse next arguments
	  depth=1;
	  continue;
        }
        Free(v);
      }
      
      if (nodeHas(theta)){
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
	sAppend(&sbDt,"_THETA_%s_",v);
        sAppend(&sbt,"THETA[%s]",v);
        Free(v);
        continue;
      }

      if (nodeHas(eta)){
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
	sAppend(&sbDt, "_ETA_%s_",v);
        sAppend(&sbt,"ETA[%s]",v);
        Free(v);
        continue;
      }
      wprint_parsetree(pt, xpn, depth, fn, client_data);
      if (rx_syntax_require_semicolon && nodeHas(end_statement) && i == 0){
        if (xpn->start_loc.s ==  xpn->end){
	  updateSyntaxCol();
          trans_syntax_error_report_fn(NEEDSEMI);
        } 
      }
      
      if (nodeHas(mult_part)){
	char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	if (i == 0){
	  if (!strcmp("/",v)){
	    aAppendN("safe_zero(", 10);
            safe_zero = 1;
	  } else {
	    safe_zero = 0;
	  }
	}
	if (i == 1){
	  if (safe_zero){
	    aAppendN(")", 1);
	  }
	  safe_zero = 0;
	}
	Free(v);
      }
      if (nodeHas(printf_statement)){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (i == 0){
	  sb.o =0; sbDt.o =0;
	  sbt.o=0;
	  aType(PPRN);
	  aAppendN("Rprintf(", 8);
	  sAppendN(&sbt,"printf(", 7);
	  sb.o--;sbDt.o--;sbt.o--;
        }
        if (i == 2){
          sAppend(&sb,"%s",v);
	  sAppend(&sbDt,"%s",v);
	  sAppend(&sbt,"%s",v);
        }
        if (i == 4){
	  addLine(&sbPm, "%s;\n", sb.s);
	  addLine(&sbPmDt, "%s;\n", sbDt.s);
	  sAppend(&sbNrm, "%s;\n", sbt.s);
        }
        Free(v);
        continue;
      } 

      if ((nodeHas(dfdy) || nodeHas(dfdy_rhs)) && i == 2){
        found_jac = 1;
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (nodeHas(dfdy_rhs)){
          // Continuation statement
	  switch(sbPm.lType[sbPm.n]){
	  case FBIO:
	    updateSyntaxCol();
	    trans_syntax_error_report_fn("Bioavailability cannot depend on Jacobian values");
	    break;
	  case ALAG:
	    updateSyntaxCol();
	    trans_syntax_error_report_fn("Absorption Lag-time cannot depend on Jacobian values");
	    break;
	  case RATE:
	    updateSyntaxCol();
	    trans_syntax_error_report_fn("Model-based rate cannot depend on Jacobian values");
	    break;
	  case DUR:
	    updateSyntaxCol();
	    trans_syntax_error_report_fn("Model-based duration cannot depend on Jacobian values");
	    break;
	  default: {
	    aType(TJAC);
	    sAppend(&sbDt, "__PDStateVar_%s_SeP_",v);
	    sAppend(&sbt,"df(%s)/dy(",v);
	    if (new_de(v)){
	      updateSyntaxCol();
	      sprintf(buf,"d/dt(%s) needs to be defined before using a Jacobians for this state",v);
	      trans_syntax_error_report_fn(buf);
	    } else {
	      sAppend(&sb, "__PDStateVar__[%d*(__NROWPD__)+",tb.id);
	    }
	  }
	  }
	  
        } else {
          // New statement
	  aType(TJAC);
          sb.o = 0; sbDt.o = 0;
          sbt.o = 0;
	  sAppend(&sbDt,"__PDStateVar_%s_SeP_",v);
	  sAppend(&sbt,"df(%s)/dy(",v);
	  if (new_de(v)){
	    updateSyntaxCol();
	    sprintf(buf,"d/dt(%s) needs to be defined before using a Jacobians for this state",v);
            trans_syntax_error_report_fn(buf);
	  } else {
	    sAppend(&sb,"__PDStateVar__[%d*(__NROWPD__)+",tb.id);
	  }
	  new_or_ith(v);
	  tb.cdf = tb.ix;
        }
        Free(v);
        continue;
      }
      if ((nodeHas(dfdy) || nodeHas(dfdy_rhs)) && i == 4){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	ii = 0;
	if (strstr(v,"THETA[") != NULL){
	  good_jac=0;
	  sprintf(buf,"_THETA_%.*s_",(int)(strlen(v))-7,v+6);
	  sAppend(&sbt, "%s)",v);
	  sAppendN(&sb, "0]", 2);
	  sAppend(&sbDt, "%s__",buf);
	  ii = 1;
	} else if (strstr(v,"ETA[") != NULL) {
	  good_jac=0;
	  sprintf(buf,"_ETA_%.*s_",(int)(strlen(v))-5,v+4);
          sAppend(&sbt, "%s)",v);
          sAppendN(&sb, "0]",2);
	  sAppend(&sbDt, "%s__",buf);
          ii = 1;
        } else {
	  sAppend(&sbDt, "%s__",v);
          sAppend(&sbt, "%s)",v);
	  new_or_ith(v);
	  if (tb.lh[tb.ix] == 9){
	    new_de(v);
	    sAppend(&sb, "%d]",tb.id);
	  } else {
	    sAppendN(&sb, "0]",2);
	    good_jac = 0;
	  }
        }
        if (strcmp("dfdy",name) == 0){
          aAppendN(" = ", 3);
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
      if (nodeHas(selection_statement) && i==1) {
	sb.o = 0; sbDt.o = 0; sbt.o = 0;
        aAppendN("if (", 4);
        sAppendN(&sbt,"if (", 4);
        continue;
      }
      if (nodeHas(selection_statement) && i==3) {
        sAppend(&sb, " {", 2);
	sAppend(&sbDt, " {", 2);
        sAppend(&sbt, "{", 1);
	aType(TLOGIC);
	addLine(&sbPm, "%s\n", sb.s);
	addLine(&sbPmDt, "%s\n", sbDt.s);
	sAppend(&sbNrm, "%s\n", sbt.s);
        continue;
      }
      if (nodeHas(selection_statement__8) && i==0) {
	sb.o = 0; sbDt.o = 0; sbt.o = 0;
	aType(TLOGIC);
	aAppendN("}\nelse {", 8);
	sAppendN(&sbt,"}\nelse {", 8);
	addLine(&sbPm, "%s\n", sb.s);
	addLine(&sbPmDt, "%s\n", sbDt.s);
	sAppend(&sbNrm, "%s\n", sbt.s);
        continue;
      }

      if (nodeHas(power_expression) && i==0) {
        aAppendN(",", 1);
        sAppendN(&sbt, "^", 1);
      }
      if (!rx_syntax_star_pow && i == 1 &&nodeHas(power_expression)){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        if (!strcmp("**",v)){
	  updateSyntaxCol();
          trans_syntax_error_report_fn(NEEDPOW);
        }
        Free(v);
      }
      if (nodeHas(transit2) && i == 0){
        aAppendN("_transit3P(t, _cSub, ", 21);
        sAppendN(&sbt,"transit(", 8);
        rx_podo = 1;
      }
      if (nodeHas(transit3) && i == 0){
        aAppendN("_transit4P(t, _cSub, ", 21);
        sAppendN(&sbt,"transit(", 8);
        rx_podo = 1;
      }
      if ((nodeHas(fbio) || nodeHas(alag) || 
	   nodeHas(dur) || nodeHas(rate) ||
	   nodeHas(cmt_statement)) && i==2) {
        /* sprintf(sb.s, "__DDtStateVar__[%d] = InfusionRate(%d) +", tb.nd, tb.nd); */
        /* sb.o = strlen(sb.s); */
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	int hasLhs=0;
	if (nodeHas(cmt_statement)){
	  new_or_ith(v);
	  if (tb.lh[tb.ix] || tb.ini[tb.ix]){
	    hasLhs=1;
	    tb.ini[tb.ix]=2;
	  }
	  if (!strcmp("depot", v)){
	    tb.hasDepotCmt = 1;
	  } else if (!strcmp("central", v)){
	    tb.hasCentralCmt = 1;
	  }
	}
        sprintf(tb.ddt, "%s",v);
        if (new_de(v)){
	  if (rx_syntax_require_ode_first){
	    if (nodeHas(cmt_statement)){
	    } else if (!strcmp("depot", v)){
	      tb.hasDepot = 1;
	    } else if (!strcmp("central", v)){
	      tb.hasCentral = 1;
	    } else {
	      updateSyntaxCol();
	      sprintf(buf,ODEFIRST,v);
	      trans_syntax_error_report_fn(buf);
	    }
	  }
	  tb.statei++;
	  if (nodeHas(fbio)){
	    sb.o=0;sbDt.o=0; sbt.o=0;
	    sAppend(&sb, "_f[%d] = ", tb.nd);
	    sAppend(&sbDt, "_f[%d] = ", tb.nd);
	    sAppend(&sbt, "f(%s)=", v);
	    if (foundF == 0) needSort+=1;// & 1 when F
	    foundF=1;
	    aType(FBIO);
	  } else if (nodeHas(alag)){
	    sb.o=0; sbDt.o=0; sbt.o=0;
	    sAppend(&sb, "_alag[%d] = ", tb.nd);
	    sAppend(&sbDt, "_alag[%d] = ", tb.nd);
	    sAppend(&sbt, "alag(%s)=", v);
	    if (foundLag == 0) needSort+=2; // & 2 when alag
	    foundLag=1;
	    aType(ALAG); 
	  } else if (nodeHas(dur)){
	    sb.o=0;sbDt.o=0; sbt.o=0;
	    sAppend(&sb, "_dur[%d] = ", tb.nd);
	    sAppend(&sbDt, "_dur[%d] = ", tb.nd);
	    sAppend(&sbt, "dur(%s)=", v);
	    if (foundDur == 0) needSort+=4;// & 4 when dur
	    foundDur=1;
	    aType(DUR);
          } else if (nodeHas(rate)){
	    sb.o=0;sbDt.o=0; sbt.o=0;
	    sAppend(&sb, "_rate[%d] = ", tb.nd);
	    sAppend(&sbDt, "_rate[%d] = ", tb.nd);
	    sAppend(&sbt, "rate(%s)=", v);
	    if (foundRate == 0) needSort+=8;// & 8 when rate
	    foundRate=1;
	    aType(RATE);
          } else if (nodeHas(cmt_statement)){
	    sb.o=0;sbDt.o=0; sbt.o=0;
	    sAppend(&sbt, "cmt(%s)", v);
	    sAppend(&sbNrm, "%s;\n", sbt.s);
	  }
          new_or_ith(v);
	  aProp(tb.nd);
          /* Rprintf("%s; tb.ini = %d; tb.ini0 = %d; tb.lh = %d\n",v,tb.ini[tb.ix],tb.ini0[tb.ix],tb.lh[tb.ix]); */
          tb.lh[tb.ix] = 9;
	  if (hasLhs){
	    tb.lh[tb.ix] = 19;
	  }	  
          tb.di[tb.nd] = tb.ix;
          sprintf(tb.de+tb.pos_de, "%s,", v);
          tb.pos_de += strlen(v)+1;
          tb.deo[++tb.nd] = tb.pos_de;
        } else {
          new_or_ith(v);
	  aProp(tb.ix);
          /* printf("de[%d]->%s[%d]\n",tb.id,v,tb.ix); */
          if (nodeHas(fbio)){
	    sb.o=0;sbDt.o=0; sbt.o=0;
	    sAppend(&sb, "_f[%d] = ", tb.id);
	    sAppend(&sbDt, "_f[%d] = ", tb.id);
	    sAppend(&sbt, "f(%s)=", v);
	    if (foundF == 0) needSort+=1;// & 1 when F
	    foundF=1;
	    aType(FBIO);
          } else if (nodeHas(alag)){
	    sb.o=0;sbDt.o=0; sbt.o=0;
	    sAppend(&sb, "_alag[%d] = ", tb.id);
	    sAppend(&sbDt, "_alag[%d] = ", tb.id);
	    sAppend(&sbt, "alag(%s)=", v);
	    if (foundLag == 0) needSort+=2; // & 2 when alag
	    foundLag=1;
	    aType(ALAG);
          } else if (nodeHas(dur)){
	    sb.o=0;sbDt.o=0; sbt.o=0;
	    sAppend(&sb, "_dur[%d] = ", tb.id);
	    sAppend(&sbDt, "_dur[%d] = ", tb.id);
	    sAppend(&sbt, "dur(%s)=", v);
	    if (foundDur == 0) needSort+=4;// & 4 when dur
	    foundDur=1;
	    aType(DUR);
          } else if (nodeHas(rate)){
	    sb.o=0;sbDt.o=0; sbt.o=0;
	    sAppend(&sb, "_rate[%d] = ", tb.id);
	    sAppend(&sbDt, "_rate[%d] = ", tb.id);
	    sAppend(&sbt, "rate(%s)=", v);
	    if (foundRate == 0) needSort+=8;// & 8 when rate
	    foundRate=1;
	    aType(RATE);
          }
        }
        Free(v);
        continue;
      }
      if (nodeHas(derivative) && i==5) {
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	if (!strcmp("+", v) || 
	    !strcmp("-", v)){
          // = + is output  or = InfusionRate + is outupt.
        } else {
	  // = + is output  or = InfusionRate + is outupt.
          aAppendN("+ ", 2);
        }
	Free(v);
	continue;
      }
      if (nodeHas(derivative) && i==2) {
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
        sprintf(tb.ddt, "%s",v);
        if (new_de(v)){
	  tb.statei++;
	  if (strncmp(v, "rx__sens_", 3) == 0){
	    tb.sensi++;
	  }
	  if (rx_syntax_allow_dots == 0 && strstr(v, ".")){
	    updateSyntaxCol();
	    trans_syntax_error_report_fn(NODOT);
	  }
	  sb.o =0; sbDt.o =0;
	  aType(TDDT);
	  aProp(tb.nd);
          sAppend(&sb, "__DDtStateVar__[%d] = ((double)(_ON[%d]))*(_IR[%d] ", tb.nd, tb.nd, tb.nd);
	  sAppend(&sbDt, "__DDtStateVar_%d__ = ((double)(_ON[%d]))*(_IR[%d] ", tb.nd, tb.nd, tb.nd);
	  sbt.o=0;
          sAppend(&sbt, "d/dt(%s)", v);
	  new_or_ith(v);
          /* Rprintf("%s; tb.ini = %d; tb.ini0 = %d; tb.lh = %d\n",v,tb.ini[tb.ix],tb.ini0[tb.ix],tb.lh[tb.ix]); */
          if (!rx_syntax_allow_assign_state &&
	      ((tb.ini[tb.ix] == 1 && tb.ini0[tb.ix] == 0) ||
	       tb.lh[tb.ix] == 1)){
	    updateSyntaxCol();
            sprintf(buf,"Cannot assign state variable %s; For initial condition assignment use '%s(0) = #'.\n  Changing states can break sensitivity analysis (for nlmixr glmm/focei).\n  To override this behavior set 'options(RxODE.syntax.assign.state = TRUE)'",v,v);
            trans_syntax_error_report_fn0(buf);
          }
	  tb.lh[tb.ix] = 9;
          tb.di[tb.nd] = tb.ix;
          sprintf(tb.de+tb.pos_de, "%s,", v);
          tb.pos_de += (int)strlen(v)+1;
	  Free(v);
	  xpn = d_get_child(pn,4);
          v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	  tb.idu[tb.nd] = 1;
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
	  sb.o =0; sbDt.o =0;
	  if (tb.idu[tb.id] == 0){
	    sAppend(&sb, "__DDtStateVar__[%d] = ((double)(_ON[%d]))*(_IR[%d] ", tb.id, tb.id, tb.id);
	    sAppend(&sbDt, "__DDtStateVar_%d__ = ((double)(_ON[%d]))*(_IR[%d] ", tb.id, tb.id, tb.id);
	  } else {
	    sAppend(&sb, "__DDtStateVar__[%d] = ((double)(_ON[%d]))*(", tb.id, tb.id);
	    sAppend(&sbDt, "__DDtStateVar_%d__ = ((double)(_ON[%d]))*(", tb.id, tb.id);
	  }
	  tb.idu[tb.id]=1;
	  aType(TDDT);
	  aProp(tb.id);
	  sbt.o=0;
	  sAppend(&sbt, "d/dt(%s)", v);
	  Free(v);
          xpn = d_get_child(pn,4);
          v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
          if (!strcmp("~",v)){
            tb.idi[tb.id] = 1;
	    sAppendN(&sbt, "~", 1);
          } else {
	    // Don't switch idi back to 0; Once the state is ignored,
	    // keep it ignored.
	    sAppendN(&sbt, "=", 1);
	  }
        }
        Free(v);
        continue;
      }
      if (nodeHas(der_rhs)) {
	switch(sbPm.lType[sbPm.n]){
	case TMTIME:
	  updateSyntaxCol();
	  trans_syntax_error_report_fn("Modeling times cannot depend on state values");
	  break;
	case FBIO:
	  updateSyntaxCol();
	  trans_syntax_error_report_fn("Bioavailability cannot depend on state values");
	  break;
	case ALAG:
	  updateSyntaxCol();
	  trans_syntax_error_report_fn("Absorption Lag-time cannot depend on state values");
	  break;
	case RATE:
	  updateSyntaxCol();
	  trans_syntax_error_report_fn("Model-based rate cannot depend on state values");
	  break;
	case DUR:
	  updateSyntaxCol();
	  trans_syntax_error_report_fn("Model-based duration cannot depend on state values");
	  break;
	default:
	  {
	    updateSyntaxCol();
	    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	    if (new_de(v)){
	      sprintf(buf2,"d/dt(%s)",v);
	      updateSyntaxCol();
	      sprintf(buf,"Tried to use d/dt(%s) before it was defined",v);
	      updateSyntaxCol();
	      trans_syntax_error_report_fn(buf);
	    } else {
	      if (sbPm.lType[sbPm.n] == TJAC){
		sAppend(&sb,   "__DDtStateVar_%d__", tb.id);
		sAppend(&sbDt, "__DDtStateVar_%d__", tb.id);
	      } else {
		sAppend(&sb,   "__DDtStateVar__[%d]", tb.id);
		sAppend(&sbDt, "__DDtStateVar_%d__", tb.id);
		aType(TDDT);
	      }
	      aProp(tb.id);
	      sAppend(&sbt, "d/dt(%s)", v);
	    }
	    Free(v);
	  } 
	}
        
        continue;
      }

      if (nodeHas(ini0f) && rx_syntax_allow_ini && i == 0){
	foundF0=1;
	aType(TF0);
	sb.o =0; sbDt.o=0; sbt.o = 0;
	char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	sAppend(&sb,  "%s",v);
	sAppend(&sbDt,"%s",v);
	sAppend(&sbt, "%s(0)",v);
	Free(v);
      }

      if ((i==0 && (nodeHas(assignment) || nodeHas(ini) || nodeHas(ini0))) ||
	  (i == 2 && nodeHas(mtime))){
        char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	tb.ddt[0]='\0';
        if ((rx_syntax_allow_ini && nodeHas(ini)) || nodeHas(ini0)){
	  sb.o =0; sbDt.o =0;
          /* aAppendN("(__0__)", 7); */
	  aType(TINI);
          for (k = 0; k < (int)strlen(v); k++){
            if (v[k] == '.'){
                aAppendN("_DoT_", 5);
		if (rx_syntax_allow_dots == 0){
	          updateSyntaxCol();
		  trans_syntax_error_report_fn(NODOT);
		}
            } else {
              sPut(&sb, v[k]);
	      sPut(&sbDt, v[k]);
            }
          }
          if (nodeHas(ini) && !new_de(v)){
	    if (tb.idu[tb.id] == 0){
	      new_or_ith(v);
	      tb.lh[tb.ix] = 19;
	    } else {
	      updateSyntaxCol();
	      sprintf(buf,"Cannot assign state variable %s; For initial condition assigment use '%s(0) ='",v,v);
	      trans_syntax_error_report_fn(buf);
	    }
          }
          if (!rx_syntax_allow_ini0 && nodeHas(ini0)){
            sprintf(buf,NOINI0,v);
	    updateSyntaxCol();
            trans_syntax_error_report_fn(buf);
          }
        } else {
          sb.o = 0; sbDt.o = 0;
          for (k = 0; k < (int)strlen(v); k++){
            if (v[k] == '.'){
	      aAppendN("_DoT_", 5);
	      if (rx_syntax_allow_dots == 0){
		updateSyntaxCol();
		trans_syntax_error_report_fn(NODOT);
	      }
            } else {
              sPut(&sb, v[k]);
	      sPut(&sbDt, v[k]);
            }
          }
          if (!new_de(v)){
	    if (tb.idu[tb.id] == 0){
	      // Change to 19 for LHS w/stateExtra
	      new_or_ith(v);
	      tb.lh[tb.ix] = 19;
	    } else {
	      sprintf(buf,"Cannot assign state variable %s; For initial condition assigment use '%s(0) ='",v,v);
	      updateSyntaxCol();
	      trans_syntax_error_report_fn(buf);
	      
	    }
            
          }
	  aType(TASSIGN);
        }
	if (nodeHas(ini0)){
	  sbt.o=0;
	  sAppend(&sbt,"%s(0)",v);
	} else if (nodeHas(mtime)){
	  sbt.o=0;
	  sAppend(&sbt, "mtime(%s)", v);
	  needSort=1;
	  aType(TMTIME);
	  nmtime++;
	} else {
	  sbt.o=0;
	  sAppend(&sbt, "%s", v);
	}
	new_or_ith(v);
	aProp(tb.ix);
	if (nodeHas(mtime)){
	  tb.lh[tb.ix] = 1;
	  tb.mtime[tb.ix] = 1;
	} else if (nodeHas(assignment)  || (!rx_syntax_allow_ini && nodeHas(ini))){
          tb.lh[tb.ix] = 1;
        } else if (nodeHas(ini) || nodeHas(ini0)){
          if (tb.ini[tb.ix] == 0){
            // If there is only one initialzation call, then assume
            // this is a parameter with an initial value.
            tb.ini[tb.ix] = 1;
            if (nodeHas(ini0)){
	      tb.ini0[tb.ix] = 1;
	      Free(v);
	      xpn = d_get_child(pn, 3);
	      v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	      sscanf(v, "%lf", &d);
	      tb.iniv[tb.ix] = d;
	      tb.ini_i++;
            } else {
	      tb.ini0[tb.ix] = 0;
              if (strncmp(v,"rx_",3)==0){
                tb.lh[tb.ix] = 1;
              } else {
		Free(v);
		xpn = d_get_child(pn, 2);
		v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
		sscanf(v, "%lf", &d);
		tb.iniv[tb.ix] = d;
		tb.ini_i++;
	      }
            }
          } else {
            // There is more than one call to this variable, it is a
            // conditional variable
	    /* Rprintf("Duplicate %s; %d %d\n", v, tb.lh[tb.ix], tb.ini0[tb.ix]); */
	    if (tb.lh[tb.ix] != 1){
	      tb.lh[tb.ix] = 1;
	      if (nodeHas(ini0) && tb.ini0[tb.ix] == 1){
		sprintf(buf,"Cannot have conditional initial conditions for %s",v);
		updateSyntaxCol();
		trans_syntax_error_report_fn(buf);
	      } else if (tb.ini0[tb.ix] == 1){
		tb.iniv[tb.ix] = NA_REAL;
		tb.ini_i--;
	      } else if (tb.ini[tb.ix] == 1){
		tb.iniv[tb.ix] = NA_REAL;
		tb.ini_i--;
	      }
	    }
	    tb.ini0[tb.ix] = 0;	      
          }
        }
        Free(v);
      }
    }

    if (nodeHas(assignment) || nodeHas(ini) || nodeHas(dfdy) ||
        nodeHas(ini0) || nodeHas(ini0f) || nodeHas(fbio) || nodeHas(alag) || nodeHas(rate) || 
	nodeHas(dur) || nodeHas(mtime)){
      addLine(&sbPm,     "%s;\n", sb.s);
      addLine(&sbPmDt,   "%s;\n", sbDt.s);
      sAppend(&sbNrm, "%s;\n", sbt.s);
    } else if (nodeHas(derivative)){
      addLine(&sbPm,     "%s);\n", sb.s);
      addLine(&sbPmDt,   "%s);\n", sbDt.s);
      sAppend(&sbNrm, "%s;\n", sbt.s);
    }

    if (!rx_syntax_assign && (nodeHas(assignment) || nodeHas(ini) || nodeHas(ini0) || nodeHas(ini0f) || nodeHas(mtime))){
      if (nodeHas(mtime)){
	i = 4;
      } else if (nodeHas(ini0)){
        i = 2;
      } else {
        i = 1;
      }
      D_ParseNode *xpn = d_get_child(pn,i);
      char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
      if (!strcmp("<-",v)){
	updateSyntaxCol();
        trans_syntax_error_report_fn(NOASSIGN);
      }
      Free(v);
    }
    
    if (nodeHas(selection_statement)){
      sb.o = 0; sbDt.o = 0; sbt.o = 0;
      aAppendN("}", 1);
      sAppendN(&sbt,"}", 1);
      aType(TLOGIC);
      addLine(&sbPm,   "%s\n", sb.s);
      addLine(&sbPmDt, "%s\n", sbDt.s);
      sAppend(&sbNrm,  "%s\n", sbt.s);
    }
    
    if (nodeHas(power_expression)) {
      aAppendN(")", 1);
    }

  }
}

void retieve_var(int i, char *buf) {
  int len;

  len = tb.vo[i+1] - tb.vo[i] - 1;
  strncpy(buf, tb.ss+tb.vo[i], len);
  buf[len] = 0;
}

void err_msgP(int chk, const char *msg, int code, D_Parser *p)
{
  if(!chk) {
    free_D_Parser(p);
    error("%s",msg);
  }
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
  buf[0]='\0';
  buf1[0]='\0';
  buf2[0]='\0';
  sAppend(&sbOut, "%s", pre_str);
  if (scenario == 0 || scenario == 2){
    // show_ode = 1 dydt
    // show_ode = 2 Jacobian
    // show_ode = 3 Ini statement
    // show_ode = 0 LHS
    // show_ode = 5 functional bioavailibility
    // show_ode == 6 functional lag
    // show_ode == 7 functional rate
    // show_ode == 8 functional duration
    // show_ode == 9 functional mtimes
    if (show_ode == 2 || show_ode == 0){
      //__DDtStateVar_#__
      for (i = 0; i < tb.nd; i++){
	if (scenario == 0){
	  sAppend(&sbOut,"  __DDtStateVar_%d__,\n",i);
	} else {
	  sAppend(&sbOut,"  (void)__DDtStateVar_%d__;\n",i);
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
	    sAppend(&sbOut,"  __PDStateVar_%s_SeP_%s__,\n",buf1,buf2);
          } else {
	    sAppend(&sbOut,"  (void)__PDStateVar_%s_SeP_%s__;\n",buf1,buf2);
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
      sAppendN(&sbOut,"  ", 2);
      for (k = 0; k < (int)strlen(buf); k++){
        if (buf[k] == '.'){
          sAppend(&sbOut,"_DoT_");
          if (rx_syntax_allow_dots == 0){
	    updateSyntaxCol();
            trans_syntax_error_report_fn(NODOT);
          }
        } else {
          sPut(&sbOut,buf[k]);
        }
      }
      if (!strcmp("rx_lambda_", buf) || !strcmp("rx_yj_", buf)){
	sAppendN(&sbOut, "__", 2);
      }
      if (i <tb.nv-1)
        sAppendN(&sbOut, ",\n", 2);
      else
        sAppendN(&sbOut, ";\n", 2);
      break;
    case 2: // Case 2 is for suppressing all the warnings for the variables by using (void)var;
      // See https://stackoverflow.com/questions/1486904/how-do-i-best-silence-a-warning-about-unused-variables
      sAppend(&sbOut,"  ");
      sAppend(&sbOut,"(void)");
      for (k = 0; k < (int)strlen(buf); k++){
        if (buf[k] == '.'){
          sAppendN(&sbOut,"_DoT_", 5);
          if (rx_syntax_allow_dots == 0){
	    updateSyntaxCol();
            trans_syntax_error_report_fn(NODOT);
          }
        } else {
          sPut(&sbOut,buf[k]);
        }
      }
      if (!strcmp("rx_lambda_", buf) || !strcmp("rx_yj_", buf)){
        sAppendN(&sbOut, "__", 2);
      }
      sAppendN(&sbOut, ";\n", 2);
      break;
    case 1:
      // Case 1 is for declaring the par_ptr.
      sAppendN(&sbOut,"  ", 2);
      for (k = 0; k < (int)strlen(buf); k++){
        if (buf[k] == '.'){
          sAppendN(&sbOut,"_DoT_", 5);
          if (rx_syntax_allow_dots == 0){
	    updateSyntaxCol();
            trans_syntax_error_report_fn(NODOT);
          }
        } else {
          sPut(&sbOut, buf[k]);
        }
      }
      sAppend(&sbOut, " = _PP[%d];\n", j++);
      break;
    default: break;
    }
  }
  sAppend(&sbOut, "%s", post_str);
}

void print_aux_info(char *model, const char *prefix, const char *libname, const char *pMd5, const char *timeId,
		    const char *libname2){
  int i, j, islhs,pi = 0,li = 0, o=0, statei = 0, sensi=0, normi=0,fdi=0,
    in_str=0;
  char buf[512], buf2[512];
  buf[0]='\0';buf2[0]='\0';
  for (i=0; i<tb.nv; i++) {
    islhs = tb.lh[i];
    if (islhs>1 && islhs != 19) continue;      /* is a state var */
    retieve_var(i, buf);
    if (islhs == 1 || islhs == 19){
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
  int nExtra=0;
  for (i=0; i<tb.nd; i++) {                     /* name state vars */
    retieve_var(tb.di[i], buf);
    if (tb.idu[i] == 1){
      if (strncmp(buf, "rx__sens_", 9) == 0){
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
    } else {
      sprintf(s_aux_info+o, "    SET_STRING_ELT(extraState, %d, mkChar(\"%s\"));\n", nExtra++, buf);
      o = (int)strlen(s_aux_info);
    }
    
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
  sAppend(&sbOut, "extern SEXP %smodel_vars(){\n  int pro=0;\n", prefix);
  sAppend(&sbOut, "  SEXP _mv = PROTECT(_rxGetModelLib(\"%smodel_vars\"));pro++;\n", prefix);
  sAppendN(&sbOut, "  if (!_rxIsCurrentC(_mv)){\n", 28);
  sAppendN(&sbOut, "    SEXP lst      = PROTECT(allocVector(VECSXP, 20));pro++;\n", 60);
  sAppendN(&sbOut, "    SEXP names    = PROTECT(allocVector(STRSXP, 20));pro++;\n", 60);
  sAppendN(&sbOut, "    SEXP sNeedSort = PROTECT(allocVector(INTSXP,1));pro++;\n", 59);
  sAppendN(&sbOut, "    int *iNeedSort  = INTEGER(sNeedSort);\n", 42);
  sAppend(&sbOut, "    iNeedSort[0] = %d;\n", needSort);
  sAppendN(&sbOut, "    SEXP sMtime = PROTECT(allocVector(INTSXP,1));pro++;\n", 56);
  sAppendN(&sbOut, "    int *iMtime  = INTEGER(sMtime);\n", 36);
  sAppend(&sbOut,  "    iMtime[0] = %d;\n", nmtime);
  sAppendN(&sbOut, "    SEXP sExtraCmt = PROTECT(allocVector(INTSXP,1));pro++;\n", 59);
  sAppendN(&sbOut, "    int *iExtraCmt  = INTEGER(sExtraCmt);\n", 42);
  sAppend(&sbOut,  "    iExtraCmt[0] = %d;\n", extraCmt);
  sAppend(&sbOut, "    SEXP params   = PROTECT(allocVector(STRSXP, %d));pro++;\n",pi);
  sAppend(&sbOut, "    SEXP lhs      = PROTECT(allocVector(STRSXP, %d));pro++;\n",li);
  sAppend(&sbOut, "    SEXP state    = PROTECT(allocVector(STRSXP, %d));pro++;\n",statei);
  sAppend(&sbOut, "  SEXP extraState = PROTECT(allocVector(STRSXP, %d));pro++;\n",nExtra);
  sAppend(&sbOut, "    SEXP stateRmS = PROTECT(allocVector(INTSXP, %d));pro++;\n",statei);
  sAppendN(&sbOut, "    SEXP timeInt = PROTECT(allocVector(INTSXP, 1));pro++;\n", 58);
  sAppend(&sbOut, "    INTEGER(timeInt)[0] = %s;\n", timeId);
  sAppend(&sbOut, "    SEXP sens     = PROTECT(allocVector(STRSXP, %d));pro++;\n",sensi);
  sAppend(&sbOut, "    SEXP normState= PROTECT(allocVector(STRSXP, %d));pro++;\n",statei-sensi);
  sAppend(&sbOut, "    SEXP fn_ini   = PROTECT(allocVector(STRSXP, %d));pro++;\n",fdi);
  sAppend(&sbOut, "    SEXP dfdy     = PROTECT(allocVector(STRSXP, %d));pro++;\n",tb.ndfdy);
  sAppendN(&sbOut, "    SEXP tran     = PROTECT(allocVector(STRSXP, 20));pro++;\n", 60);
  sAppendN(&sbOut, "    SEXP trann    = PROTECT(allocVector(STRSXP, 20));pro++;\n", 60);
  sAppendN(&sbOut, "    SEXP mmd5     = PROTECT(allocVector(STRSXP, 2));pro++;\n", 59);
  sAppendN(&sbOut, "    SEXP mmd5n    = PROTECT(allocVector(STRSXP, 2));pro++;\n", 59);
  sAppendN(&sbOut, "    SEXP model    = PROTECT(allocVector(STRSXP, 1));pro++;\n", 59);
  sAppendN(&sbOut, "    SEXP modeln   = PROTECT(allocVector(STRSXP, 1));pro++;\n", 59);
  sAppendN(&sbOut, "    SEXP version    = PROTECT(allocVector(STRSXP, 3));pro++;\n", 61);
  sAppendN(&sbOut, "    SEXP versionn   = PROTECT(allocVector(STRSXP, 3));pro++;\n", 61);
  
  sAppend(&sbOut,  __VER_0__);
  sAppend(&sbOut,  __VER_1__);
  sAppend(&sbOut,  __VER_2__);

  sAppendN(&sbOut, "    SET_STRING_ELT(versionn,0,mkChar(\"version\"));\n", 50);
  sAppendN(&sbOut, "    SET_STRING_ELT(versionn,1,mkChar(\"repo\"));\n", 47);
  sAppendN(&sbOut, "    SET_STRING_ELT(versionn,2,mkChar(\"md5\"));\n", 46);

  sAppend(&sbOut, "%s",s_aux_info);
  // Save for outputting in trans
  tb.fdn = fdi;
  tb.pi = pi;
  tb.li = li;
  tb.sensi  = sensi;
  sAppendN(&sbOut, "    SET_STRING_ELT(modeln,0,mkChar(\"normModel\"));\n", 50);
  sAppendN(&sbOut, "    SET_STRING_ELT(model,0,mkChar(\"", 35);
  in_str=0;
  for (i = 0; i < sbNrm.o; i++){
    if (sbNrm.s[i] == '"'){
      if (in_str==1){
	in_str=0;
      } else {
	in_str=1;
      }
      sAppendN(&sbOut, "\\\"", 2);
    } else if (sbNrm.s[i] == '\''){
      if (in_str==1){
	in_str=0;
      } else {
	in_str=1;
      }
      sAppendN(&sbOut, "'", 1);
    } else if (sbNrm.s[i] == ' '){
      if (in_str==1){
	sAppendN(&sbOut, " ", 1);
      }
    } else if (sbNrm.s[i] == '\n'){
      sAppendN(&sbOut, "\\n", 2);
    } else if (sbNrm.s[i] == '\t'){
      sAppendN(&sbOut, "\\t", 2);
    } else if (sbNrm.s[i] == '\\'){
      sAppendN(&sbOut, "\\\\", 2);
    } else if (sbNrm.s[i] >= 33  && sbNrm.s[i] <= 126){ // ASCII only
      sPut(&sbOut, sbNrm.s[i]);
    }
  }
  sAppendN(&sbOut, "\"));\n", 5);
  
  s_aux_info[0] = '\0';
  o    = 0;
  tb.ini_i = getInits(s_aux_info, &o); 

  /* tb.ini_i = length(ini); */
  /* for (i = 0; i < tb.ini_i; i++){ */
  /*   sprintf(s_aux_info+o,"    SET_STRING_ELT(inin,%d,mkChar(\"%s\"));\n",i, CHAR(STRING_ELT(inin, i))); */
  /*   o = (int)strlen(s_aux_info); */
  /*   sprintf(s_aux_info+o,"    REAL(ini)[%d] = %.16f;\n",i, REAL(ini)[i]); */
  /*   o = (int)strlen(s_aux_info); */
  /* } */
  
  sAppend(&sbOut, "    SEXP ini    = PROTECT(allocVector(REALSXP,%d));pro++;\n",tb.ini_i);
  sAppend(&sbOut, "    SEXP inin   = PROTECT(allocVector(STRSXP, %d));pro++;\n",tb.ini_i);
  sAppend(&sbOut, "%s",s_aux_info);
  // Vector Names
  sAppendN(&sbOut, "    SET_STRING_ELT(names,0,mkChar(\"params\"));\n", 46);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  0,params);\n", 36);

  sAppendN(&sbOut, "    SET_STRING_ELT(names,1,mkChar(\"lhs\"));\n", 43);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  1,lhs);\n", 33);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(names,2,mkChar(\"state\"));\n", 45);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  2,state);\n", 35);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(names,3,mkChar(\"trans\"));\n", 45);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  3,tran);\n", 34);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(names,4,mkChar(\"model\"));\n", 45);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  4,model);\n", 35);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(names,5,mkChar(\"ini\"));\n", 43);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  5,ini);\n", 33);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(names,6,mkChar(\"podo\"));\n", 44);
  sAppend(&sbOut, "    SET_VECTOR_ELT(lst,   6,ScalarLogical(%d));\n",rx_podo);

  sAppendN(&sbOut, "    SET_STRING_ELT(names,7,mkChar(\"dfdy\"));\n", 44);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  7,dfdy);\n", 34);

  sAppendN(&sbOut, "    SET_STRING_ELT(names,8,mkChar(\"sens\"));\n", 44);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  8,sens);\n", 34);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(names,9,mkChar(\"fn.ini\"));\n", 46);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  9,fn_ini);\n", 36);

  sAppendN(&sbOut, "    SET_STRING_ELT(names,10,mkChar(\"state.ignore\"));\n", 53);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  10,stateRmS);\n", 39);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(names,11,mkChar(\"version\"));\n", 48);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  11,version);\n", 38);

  sAppendN(&sbOut, "    SET_STRING_ELT(names,12,mkChar(\"normal.state\"));\n", 53);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  12,normState);\n", 40);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(names,13,mkChar(\"needSort\"));\n", 49);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  13,sNeedSort);\n", 40);

  sAppendN(&sbOut, "    SET_STRING_ELT(names,14,mkChar(\"nMtime\"));\n", 47);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  14,sMtime);\n", 37);

  sAppendN(&sbOut, "    SET_STRING_ELT(names,15,mkChar(\"extraCmt\"));\n", 49);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  15,sExtraCmt);\n", 40);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(names, 16, mkChar(\"stateExtra\"));\n", 53);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  16, extraState);\n", 42);

  sAppendN(&sbOut, "    SET_STRING_ELT(names, 17, mkChar(\"dvid\"));\n", 47);
  sAppend(&sbOut,   "    SEXP sDvid = PROTECT(allocVector(INTSXP,%d));pro++;\n", tb.dvidn);
  
  for (int di = 0; di < tb.dvidn; di++){
    sAppend(&sbOut, "    INTEGER(sDvid)[%d] = %d;\n",di, tb.dvid[di]);
  }
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst, 17, sDvid);\n", 36);


  sAppendN(&sbOut, "    SET_STRING_ELT(names,18,mkChar(\"timeId\"));\n", 47);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  18,timeInt);\n", 38);

  sAppendN(&sbOut, "    SET_STRING_ELT(names,19,mkChar(\"md5\"));\n", 43);
  sAppendN(&sbOut, "    SET_VECTOR_ELT(lst,  19,mmd5);\n", 34);

  // const char *rxVersion(const char *what)
  
  // md5 values
  sAppendN(&sbOut, "    SET_STRING_ELT(mmd5n,0,mkChar(\"file_md5\"));\n", 48);
  sAppend(&sbOut, "    SET_STRING_ELT(mmd5,0,mkChar(\"%s\"));\n",md5);
  sAppendN(&sbOut, "    SET_STRING_ELT(mmd5n,1,mkChar(\"parsed_md5\"));\n", 50);
  sAppend(&sbOut, "    SET_STRING_ELT(mmd5,1,mkChar(\"%s\"));\n", pMd5);
  
  // now trans output
  sAppendN(&sbOut, "    SET_STRING_ELT(trann,0,mkChar(\"lib.name\"));\n", 48);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 0,mkChar(\"%s\"));\n", libname);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,1,mkChar(\"jac\"));\n", 43);
  if (found_jac == 1 && good_jac == 1){
    sAppendN(&sbOut, "    SET_STRING_ELT(tran,1,mkChar(\"fulluser\"));\n", 47); // Full User Matrix
  } else {
    sAppendN(&sbOut, "    SET_STRING_ELT(tran,1,mkChar(\"fullint\"));\n", 46); // Full Internal Matrix
  }
  sAppendN(&sbOut, "    SET_STRING_ELT(trann,2,mkChar(\"prefix\"));\n", 46);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 2,mkChar(\"%s\"));\n", prefix);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,3,mkChar(\"dydt\"));\n", 44);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 3,mkChar(\"%sdydt\"));\n", prefix);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,4,mkChar(\"calc_jac\"));\n", 48);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 4,mkChar(\"%scalc_jac\"));\n", prefix);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,5,mkChar(\"calc_lhs\"));\n", 48);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 5,mkChar(\"%scalc_lhs\"));\n", prefix);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,6,mkChar(\"model_vars\"));\n", 50);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 6,mkChar(\"%smodel_vars\"));\n", prefix);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(trann,7,mkChar(\"theta\"));\n", 45);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 7,mkChar(\"%stheta\"));\n", prefix);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(trann,8,mkChar(\"inis\"));\n", 44);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 8,mkChar(\"%sinis\"));\n", prefix);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(trann,  9,mkChar(\"dydt_lsoda\"));\n", 52);
  sAppend(&sbOut, "    SET_STRING_ELT(tran,   9,mkChar(\"%sdydt_lsoda\"));\n", prefix);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(trann,10,mkChar(\"calc_jac_lsoda\"));\n", 55);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 10,mkChar(\"%scalc_jac_lsoda\"));\n", prefix);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(trann,11,mkChar(\"ode_solver_solvedata\"));\n", 61);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 11,mkChar(\"%sode_solver_solvedata\"));\n", prefix);
  
  sAppendN(&sbOut, "    SET_STRING_ELT(trann,12,mkChar(\"ode_solver_get_solvedata\"));\n", 65);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 12,mkChar(\"%sode_solver_get_solvedata\"));\n", prefix);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,13,mkChar(\"dydt_liblsoda\"));\n", 54);
  sAppend(&sbOut, "    SET_STRING_ELT(tran, 13,mkChar(\"%sdydt_liblsoda\"));\n", prefix);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,14,mkChar(\"F\"));\n", 42);
  sAppend(&sbOut,  "    SET_STRING_ELT(tran, 14,mkChar(\"%sF\"));\n", prefix);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,15,mkChar(\"Lag\"));\n", 44);
  sAppend(&sbOut,  "    SET_STRING_ELT(tran, 15,mkChar(\"%sLag\"));\n", prefix);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,16,mkChar(\"Rate\"));\n", 45);
  sAppend(&sbOut,  "    SET_STRING_ELT(tran, 16,mkChar(\"%sRate\"));\n", prefix);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,17,mkChar(\"Dur\"));\n", 44);
  sAppend(&sbOut,  "    SET_STRING_ELT(tran, 17,mkChar(\"%sDur\"));\n", prefix);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,18,mkChar(\"mtime\"));\n", 46);
  sAppend(&sbOut,  "    SET_STRING_ELT(tran, 18,mkChar(\"%smtime\"));\n", prefix);

  sAppendN(&sbOut, "    SET_STRING_ELT(trann,19,mkChar(\"assignFuns\"));\n", 51);
  sAppend(&sbOut,  "    SET_STRING_ELT(tran, 19,mkChar(\"%sassignFuns\"));\n", prefix);

  sAppendN(&sbOut, "    setAttrib(tran, R_NamesSymbol, trann);\n", 43);
  sAppendN(&sbOut, "    setAttrib(mmd5, R_NamesSymbol, mmd5n);\n", 43);
  sAppendN(&sbOut, "    setAttrib(model, R_NamesSymbol, modeln);\n", 45);
  sAppendN(&sbOut, "    setAttrib(ini, R_NamesSymbol, inin);\n", 41);
  sAppendN(&sbOut, "    setAttrib(version, R_NamesSymbol, versionn);\n", 49);
  sAppendN(&sbOut, "    setAttrib(lst, R_NamesSymbol, names);\n", 42);
  sAppendN(&sbOut, "    SEXP cls = PROTECT(allocVector(STRSXP, 1));pro++;\n", 54);
  sAppendN(&sbOut, "    SET_STRING_ELT(cls, 0, mkChar(\"rxModelVars\"));\n", 51);
  sAppendN(&sbOut, "    classgets(lst, cls);\n", 25);
  sAppendN(&sbOut, "    _assign_ptr(lst);\n", 22);
  sAppendN(&sbOut, "    UNPROTECT(pro);\n", 20);
  
  sAppendN(&sbOut, "    return lst;\n", 16);
  sAppendN(&sbOut, "  } else {\n", 11);
  sAppendN(&sbOut, "    UNPROTECT(pro);\n", 20);
  sAppendN(&sbOut, "    return _mv;\n", 16);
  sAppendN(&sbOut, "  }\n", 4);
  sAppendN(&sbOut, "}\n", 2);

  sAppend(&sbOut,"extern void %sdydt_lsoda(int *neq, double *t, double *A, double *DADT)\n{\n  %sdydt(neq, *t, A, DADT);\n}\n", prefix, prefix);
  sAppend(&sbOut, "extern int %sdydt_liblsoda(double t, double *y, double *ydot, void *data)\n{\n  int *neq = (int*)(data);\n  %sdydt(neq, t, y, ydot);\n  return(0);\n}\n",
	  prefix,prefix);
  sAppend(&sbOut,"extern void %scalc_jac_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){\n  // Update all covariate parameters\n  %scalc_jac(neq, *t, A, JAC, *nrowpd);\n}\n",
	  prefix, prefix);
  sAppend(&sbOut,"\n//Create function to call from R's main thread that assigns the required functions. Sometimes they don't get assigned.\nextern void %sassignFuns(){\n  _assignFuns();\n}\n", prefix);
  sAppend(&sbOut,"\n//Initialize the dll to match RxODE's calls\nvoid R_init0_%s(){\n  // Get C callables on load; Otherwise it isn't thread safe\n", libname2);
  sAppendN(&sbOut, "  _assignFuns();\n", 17);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sassignFuns\", (DL_FUNC) %sassignFuns);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%stheta\", (DL_FUNC) %stheta);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sinis\",(DL_FUNC) %sinis);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sdydt\",(DL_FUNC) %sdydt);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%scalc_lhs\",(DL_FUNC) %scalc_lhs);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%scalc_jac\",(DL_FUNC) %scalc_jac);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sdydt_lsoda\", (DL_FUNC) %sdydt_lsoda);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%scalc_jac_lsoda\", (DL_FUNC) %scalc_jac_lsoda);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sode_solver_solvedata\", (DL_FUNC) %sode_solver_solvedata);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sode_solver_get_solvedata\", (DL_FUNC) %sode_solver_get_solvedata);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sF\", (DL_FUNC) %sF);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sLag\", (DL_FUNC) %sLag);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sRate\", (DL_FUNC) %sRate);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sDur\", (DL_FUNC) %sDur);\n", libname, prefix, prefix);
sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%smtime\", (DL_FUNC) %smtime);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sdydt_liblsoda\", (DL_FUNC) %sdydt_liblsoda);\n", libname, prefix, prefix);
  sAppend(&sbOut,"}\n//Initialize the dll to match RxODE's calls\nvoid R_init_%s(DllInfo *info){\n  // Get C callables on load; Otherwise it isn't thread safe\n  R_init0_%s();", libname2, libname2);
  sAppend(&sbOut, "\n  static const R_CallMethodDef callMethods[]  = {\n    {\"%smodel_vars\", (DL_FUNC) &%smodel_vars, 0},\n    {NULL, NULL, 0}\n  };\n",
  	  prefix, prefix);
  sAppendN(&sbOut, "\n  R_registerRoutines(info, NULL, callMethods, NULL, NULL);\n  R_useDynamicSymbols(info,FALSE);\n}\n", 97);
  sAppend(&sbOut, "\nvoid R_unload_%s (DllInfo *info){\n  // Free resources required for single subject solve.\n  SEXP _mv = PROTECT(_rxGetModelLib(\"%smodel_vars\"));\n",
	  libname2, prefix);
  sAppend(&sbOut, "  if (!isNull(_mv)){\n    _rxRmModelLib(\"%smodel_vars\");\n  }\n  UNPROTECT(1);\n}\n", prefix);
}


void codegen(char *model, int show_ode, const char *prefix, const char *libname, const char *pMd5, const char *timeId, const char *fixInis0, const char *fixInis1, const char *libname2) {
  if (show_ode == 4) {
    print_aux_info(model, prefix, libname, pMd5, timeId, libname2);
  } else {
    int i, j, k;
    char buf[64];
    buf[0]='\0';
    if (show_ode == 1){
      sAppendN(&sbOut,"#include <RxODE_model_shared.h>\n",32);
      int mx = maxSumProdN;
      if (SumProdLD > mx) mx = SumProdLD;
      sAppend(&sbOut,"#define __MAX_PROD__ %d\n", mx);
      int baseSize = tb.statei-tb.nExtra+extraCmt - tb.sensi;
      if (tb.sensi > 0){
	// This converts CMT to user CMT in model
	// Hence CMT = 4 could translate in data to 44 with sensi=10
	// Then cmt=44 translates back to cmt-10 or 4.
	// This makes the sensitivity equations insensitive to CMT changes that occur in FOCEi
	sAppend(&sbOut,"#define _CMT ((abs(CMT)<=%d) ? CMT : ((CMT<0) ? CMT+%d: CMT-%d))\n",
	      baseSize, tb.sensi, tb.sensi);
      } else {
	sAppendN(&sbOut,"#define _CMT CMT\n", 17);
      }
      sAppend(&sbOut, "extern void  %sode_solver_solvedata (rx_solve *solve){\n  _solveData = solve;\n}\n",prefix);
      sAppend(&sbOut, "extern rx_solve *%sode_solver_get_solvedata(){\n  return _solveData;\n}\n", prefix);
      sAppend(&sbOut, "SEXP %smodel_vars();\n", prefix);
      sAppend(&sbOut, "%s\nextern double* %stheta(double *theta){\n  %s\n  return _theta;\n}\n", fixInis0, prefix, fixInis1);
      sAppendN(&sbOut,"\n", 1);
      sAppendN(&sbOut, "\n// prj-specific differential eqns\nvoid ", 40);
      sAppend(&sbOut, "%sdydt(int *_neq, double t, double *__zzStateVar__, double *__DDtStateVar__)\n{\n  int _cSub = _neq[1];\n", prefix);
    } else if (show_ode == 2){
      sAppend(&sbOut, "// Jacobian derived vars\nvoid %scalc_jac(int *_neq, double t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {\n  int _cSub=_neq[1];\n", prefix);
    } else if (show_ode == 3){
      sAppend(&sbOut,  "// Functional based initial conditions.\nvoid %sinis(int _cSub, double *__zzStateVar__){\n", prefix);
      if (foundF0){
	sAppendN(&sbOut, "  double t=0;\n", 14);
      }
    } else if (show_ode == 5){
      if (foundF){
	sAppend(&sbOut,  "// Functional based bioavailability (returns amount)\ndouble %sF(int _cSub,  int _cmt, double _amt, double t){\n  double _f[%d]={1};\n  (void)_f;\n",
		prefix, tb.nd);
      } else {
	sAppend(&sbOut,  "// Functional based bioavailability\ndouble %sF(int _cSub,  int _cmt, double _amt, double t){\n return _amt;\n",
		prefix);
      }
    } else if (show_ode == 6){
      if (foundLag){
	sAppend(&sbOut,  "// Functional based absorption lag\ndouble %sLag(int _cSub,  int _cmt, double t){\n  double _alag[%d]={0};\n  (void)_alag;\n",
		prefix, tb.nd);
      } else {
	sAppend(&sbOut,  "// Functional based absorption lag\ndouble %sLag(int _cSub,  int _cmt, double t){\n return t;\n",
		prefix);
      }
    } else if (show_ode == 7){
      if (foundRate){
	sAppend(&sbOut,  "// Modeled zero-order rate\ndouble %sRate(int _cSub,  int _cmt, double _amt, double t){\n  double _rate[%d]={0};\n  (void)_rate;\n",
		prefix, tb.nd);
      } else {
	sAppend(&sbOut,  "// Modeled zero-order rate\ndouble %sRate(int _cSub,  int _cmt, double _amt, double t){\n return 0.0;\n",
		prefix);
      }
    } else if (show_ode == 8){
      if (foundDur){
	sAppend(&sbOut,  "// Modeled zero-order duration\ndouble %sDur(int _cSub,  int _cmt, double _amt, double t){\n  double _dur[%d]={0};\n  (void)_dur;\n",
		prefix, tb.nd);
      } else {
	sAppend(&sbOut,  "// Modeled zero-order duration\ndouble %sDur(int _cSub,  int _cmt, double _amt, double t){\n return 0.0;\n",
		prefix);
      }
    } else if (show_ode == 9){
      if (nmtime){
	sAppend(&sbOut,  "// Model Times\nvoid %smtime(int _cSub, double *_mtime){\n  double t = 0;\n  ",
		prefix);
      } else {
	sAppend(&sbOut,  "// Model Times\nvoid %smtime(int _cSub, double *_mtime){\n",
		prefix);
      }
      
    } else {
      sAppend(&sbOut,  "// prj-specific derived vars\nvoid %scalc_lhs(int _cSub, double t, double *__zzStateVar__, double *_lhs) {\n", prefix);
    }
    if ((show_ode == 2 && found_jac == 1 && good_jac == 1) ||
	(show_ode != 2 && show_ode != 3 && show_ode != 5  && show_ode != 8 &&
	 show_ode !=0 && show_ode != 9) ||
	(show_ode == 8 && foundDur) ||
	(show_ode == 7 && foundRate) ||
	(show_ode == 6 && foundLag) ||
	(show_ode == 5 && foundF) ||
	(show_ode == 3 && foundF0) || 
	(show_ode == 0 && tb.li) ||
	(show_ode == 9 && nmtime)){
      prnt_vars(0, 0, "  double ", "\n",show_ode);     /* declare all used vars */
      if (maxSumProdN > 0 || SumProdLD > 0){
	int mx = maxSumProdN;
	if (SumProdLD > mx) mx = SumProdLD;
	sAppend(&sbOut,  "  double _p[%d], _input[%d];\n", mx, mx);
	sAppend(&sbOut,  "  double _pld[%d];\n", mx);
      }
      else prnt_vars(2, 0, "  (void)t;\n", "\n",show_ode);     /* declare all used vars */
      if (maxSumProdN){
	sAppendN(&sbOut,  "  (void)_p;\n  (void)_input;\n", 28);
	if (SumProdLD){
	  sAppendN(&sbOut,  "  (void)_pld;\n", 14);
	}
      }
      if (show_ode == 3){
	sAppendN(&sbOut, "  _update_par_ptr(0.0, _cSub, _solveData, _idx);\n", 49);
      } else if (show_ode == 6 || show_ode == 7 || show_ode == 8 || show_ode == 9){
	sAppendN(&sbOut, "  _update_par_ptr(NA_REAL, _cSub, _solveData, _idx);\n", 53);
      } else {
	sAppendN(&sbOut, "  _update_par_ptr(t, _cSub, _solveData, _idx);\n", 47);
      }
      prnt_vars(1, 1, "", "\n",show_ode);                   /* pass system pars */
      if (show_ode != 7 && show_ode != 5 &&
	  show_ode != 6 && show_ode != 8 && show_ode != 9){
	for (i=0; i<tb.nd; i++) {                   /* name state vars */
	  retieve_var(tb.di[i], buf);
	  sAppendN(&sbOut, "  ", 2);
	  for (k = 0; k < (int)strlen(buf); k++){
	    if (buf[k] == '.'){
	      sAppendN(&sbOut, "_DoT_", 5);
	      if (rx_syntax_allow_dots == 0){
		updateSyntaxCol();
		trans_syntax_error_report_fn(NODOT);
	      }
	    } else {
	      sPut(&sbOut, buf[k]);
	    }
	  }
	  sAppend(&sbOut, " = __zzStateVar__[%d]*((double)(_ON[%d]));\n", i, i);
	}
	sAppendN(&sbOut, "\n", 1);
      }
    }
    if ((foundDur && show_ode == 8) ||
	(foundRate && show_ode == 7) ||
	(foundLag && show_ode == 6) ||
	(foundF && show_ode == 5) ||
	(foundF0 && show_ode == 3) ||
	(show_ode == 0 && tb.li) ||
	(show_ode == 9 && nmtime) ||
	(show_ode == 2 && found_jac == 1 && good_jac == 1) ||
	(show_ode != 9 && show_ode != 0 && show_ode != 2 && show_ode != 3 && show_ode != 5 && show_ode != 6  && show_ode != 7 && show_ode != 8)){
      for (i = 0; i < sbPm.n; i++){
	switch(sbPm.lType[i]){
	case TMTIME:
	case TASSIGN:
	  sAppend(&sbOut,"  %s",show_ode == 1 ? sbPm.line[i] : sbPmDt.line[i]);
	  break;
	case TINI:
	  // See if this is an ini or a reclaimed expression.
	  if (sbPm.lProp[i] >= 0 ){
	    tb.ix = sbPm.lProp[i];
	    if (tb.lh[tb.ix] == 1){
	      sAppend(&sbOut,"  %s",show_ode == 1 ? sbPm.line[i] : sbPmDt.line[i]);
	    }
	  }	  
	  break;
	case TF0:
	  // functional ini
	  if (show_ode == 3) sAppend(&sbOut,"  %s",sbPmDt.line[i]);
	  break;	  
	case FBIO:
	  if (show_ode == 5) sAppend(&sbOut,"  %s", sbPmDt.line[i]);
	  break;
	case ALAG:
	  if (show_ode == 6) sAppend(&sbOut, "  %s", sbPmDt.line[i]);
	  break;
	case RATE:
	  if (show_ode == 7) sAppend(&sbOut, "  %s", sbPmDt.line[i]);
	  break;
	case DUR:
	  if (show_ode == 8) sAppend(&sbOut,"  %s", sbPmDt.line[i]);
	  break;
	case TJAC:
	  if (show_ode == 0) sAppend(&sbOut, "  %s", sbPmDt.line[i]);
	  else if (show_ode == 2)  sAppend(&sbOut, "  %s", sbPm.line[i]);
	  break;
	case TDDT:
	  // d/dt()
	  if (show_ode != 3 && show_ode != 5 && show_ode != 6 &&
	      show_ode != 7 && show_ode != 8 && show_ode != 9){
	    sAppend(&sbOut, "  %s", show_ode == 1 ? sbPm.line[i] : sbPmDt.line[i]);
	  }
	  break;
	case PPRN:
	  // Rprintf
	  if (show_ode == 1){
	    sAppend(&sbOut, "  %s", show_ode == 1 ? sbPm.line[i] : sbPmDt.line[i]);
	  }
	  break;
	case TLOGIC:
	  sAppend(&sbOut,"  %s",show_ode == 1 ? sbPm.line[i] : sbPmDt.line[i]);
	  break;
	default:
	  Rprintf("Ignored line: \"%s\"\n\tlProp: %d\tlType: %d\n", sbPm.line[i],
		  sbPm.lProp[i], sbPm.lType[i]);
	}
      }
      // End statements
      switch (show_ode){
      case 8:
	// RATE
	sAppendN(&sbOut, "\n  return _dur[_cmt];\n", 22);
	break;
      case 7:
	// DUR
	sAppendN(&sbOut, "\n  return _rate[_cmt];\n", 23);
	break;
      case 6:
	// Alag
	sAppendN(&sbOut, "\n  return t + _alag[_cmt];\n", 27);
	break;
      case 5:
	sAppendN(&sbOut, "\n  return _f[_cmt]*_amt;\n", 25);
	break;
      }
    }
    if (show_ode == 1){
      sAppendN(&sbOut,  "  (&_solveData->subjects[_cSub])->dadt_counter[0]++;\n}\n\n", 56);
    } else if (show_ode == 2){
      //sAppendN(&sbOut, "  free(__ld_DDtStateVar__);\n");
      sAppendN(&sbOut,  "  (&_solveData->subjects[_cSub])->jac_counter[0]++;\n", 52);
      sAppendN(&sbOut,  "}\n", 2);
    } else if (show_ode == 3){
      if (foundF0){
	for (i = 0; i < tb.nd; i++){
	  retieve_var(tb.di[i], buf);
	  sAppend(&sbOut, "  __zzStateVar__[%d]=((double)(_ON[%d]))*(",i,i);
	  for (k = 0; k < (int)strlen(buf); k++){
	    if (buf[k] == '.'){
	      sAppendN(&sbOut, "_DoT_", 5);
	      if (rx_syntax_allow_dots == 0){
		updateSyntaxCol();
		trans_syntax_error_report_fn(NODOT);
	      }
	    } else {
	      sPut(&sbOut, buf[k]);
	    }
	  }
	  sAppendN(&sbOut,  ");\n", 3);
	}
      }
      sAppendN(&sbOut,  "}\n", 2);
    } else if (show_ode == 5 || show_ode == 6 || show_ode == 7 || show_ode == 8){
      sAppendN(&sbOut,  "}\n", 2);
    } else if (show_ode == 0 && tb.li){
      sAppendN(&sbOut,  "\n", 1);
      for (i=0, j=0; i<tb.nv; i++) {
	if (tb.lh[i] != 1 && tb.lh[i] != 19) continue;
	retieve_var(i, buf);
	sAppend(&sbOut,  "  _lhs[%d]=", j);
	for (k = 0; k < (int)strlen(buf); k++){
	  if (buf[k] == '.'){
	    sAppendN(&sbOut, "_DoT_", 5);
	    if (rx_syntax_allow_dots == 0){
	      updateSyntaxCol();
	      trans_syntax_error_report_fn(NODOT);
	    }
	  } else {
	    sPut(&sbOut, buf[k]);
	  }
	}
	sAppendN(&sbOut,  ";\n", 2);
	j++;
      }
      sAppendN(&sbOut,  "}\n", 2);
    } else if (show_ode == 9 && nmtime){
      sAppendN(&sbOut,  "\n", 1);
      for (i=0, j=0; i<tb.nv; i++) {
	if (tb.mtime[i] != 1) continue;
	retieve_var(i, buf);
	sAppend(&sbOut,  "  _mtime[%d]=", j);
	for (k = 0; k < (int)strlen(buf); k++){
	  if (buf[k] == '.'){
	    sAppendN(&sbOut, "_DoT_", 5);
	    if (rx_syntax_allow_dots == 0){
	      updateSyntaxCol();
	      trans_syntax_error_report_fn(NODOT);
	    }
	  } else {
	    sPut(&sbOut, buf[k]);
	  }
	}
	sAppendN(&sbOut,  ";\n", 2);
	j++;
      }
      sAppendN(&sbOut,  "}\n", 2);
    } else {
      sAppendN(&sbOut,  "}\n", 2);
    }
  }
}
  
void reset (){
  // Reset sb/sbt string buffers
  sIni(&sb);
  sIni(&sbDt);
  sIni(&sbt);
  sIni(&sbNrm);
  lineIni(&sbPm);
  lineIni(&sbPmDt);

  // Reset Arrays
  memset(tb.ss,		0, 64*MXSYM*sizeof(char));
  memset(tb.de,		0, 64*MXSYM*sizeof(char));
  memset(tb.deo,	0, MXSYM*sizeof(int));
  memset(tb.vo,		0, MXSYM*sizeof(int));
  memset(tb.lh,		0, MXSYM*sizeof(int));
  memset(tb.ini,	0, MXSYM*sizeof(int));
  memset(tb.mtime,	0, MXSYM*sizeof(int));
  memset(tb.di,		0, MXDER*sizeof(int));
  memset(tb.fdi,        0, MXDER*sizeof(int));
  memset(tb.dy,		0, MXSYM*sizeof(int));
  memset(tb.sdfdy,	0, MXSYM*sizeof(int));
  memset(tb.idu,        0, MXDER*sizeof(int));
  memset(tb.idi,        0, MXDER*sizeof(int));
  memset(tb.dvid,       0, MXDER*sizeof(int));
  // Reset integers
  tb.dvidn      = 0;
  tb.nv		= 0;
  tb.ix		= 0;
  tb.id		= 0;
  tb.fn		= 0;
  tb.nd		= 0;
  tb.pos	= 0;
  tb.pos_de	= 0;
  tb.ini_i	= 0;
  tb.nExtra     = 0;
  tb.statei	= 0;
  tb.sensi	= 0;
  tb.li		= 0;
  tb.pi		= 0;
  tb.cdf	= 0;
  tb.ndfdy	= 0;
  tb.maxtheta   = 0;
  tb.maxeta     = 0;
  tb.fdn        = 0;
  tb.linCmt     = 0;
  tb.isPi       = 0;
  tb.ini_i      = 0;
  tb.hasDepot   = 0;
  tb.hasCentral = 0;
  tb.hasKa      = 0;
  tb.hasDepotCmt = 0;
  tb.hasCentralCmt = 0;
  // reset globals
  good_jac = 1;
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

  foundLag=0;
  foundRate=0;
  foundDur=0;
  foundF=0;
  foundF0=0;
  nmtime=0;
  Free(md5);
  syntaxErrorExtra=0;
  lastSyntaxErrorLine=0;
  gBufLast=0;
  lastStrLoc=0;
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
static void rxSyntaxError(struct D_Parser *ap);

void trans_internal(char* parse_file, int isStr){
  char buf1[512], buf2[512], bufe[2048];
  int i,j,found,islhs;
  D_ParseNode *pn;
  /* any number greater than sizeof(D_ParseNode_User) will do;
     below 1024 is used */
  D_Parser *p = new_D_Parser(&parser_tables_RxODE, 1024);
  curP = p;
  p->save_parse_tree = 1;
  p->syntax_error_fn = rxSyntaxError;
  if (isStr){
    gBuf = parse_file;
  } else {
      gBuf = r_sbuf_read(parse_file);
      err_msgP((intptr_t) gBuf, "error: empty buf for FILE_to_parse\n", -2, p);
  }
  sIni(&sbNrm);
  lineIni(&sbPm);
  lineIni(&sbPmDt);
  
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
	trans_syntax_error_report_fn0(bufe);
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
          if (islhs>1 && tb.lh[i] != 19) continue; /* is a state var */
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
        trans_syntax_error_report_fn0(bufe);
      }
    }
  } else {
    rx_syntax_error = 1;
  }
  free_D_Parser(p);
}

SEXP _RxODE_trans(SEXP parse_file, SEXP extra_c, SEXP prefix, SEXP model_md5, SEXP parseStr,
		  SEXP isEscIn){
  char *in;
  char buf[1024], buf2[512], df[128], dy[128];
  int i, j, islhs, pi=0, li=0, ini_i = 0,k=0, l=0, m=0, p=0;
  // Make sure buffers are initialized.
  buf2[0]='\0';
  buf[0]='\0';
  df[0]='\0';
  dy[0]='\0';

  isEsc=INTEGER(isEscIn)[0];

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
  rx_syntax_require_ode_first = R_get_option("RxODE.syntax.require.ode.first",1);
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
    Free(md5);
    md5 = rc_dup_str(CHAR(STRING_ELT(model_md5,0)),0);
    if (strlen(md5)!= 32){
      md5 = Calloc(1,char);
      md5[0] = '\0';
    }
  } else {
    Free(md5);
    md5 = Calloc(1,char);
    md5[0] = '\0';
  }
  
  trans_internal(in, isStr);
  extraCmt = 0;
  if (tb.linCmt){
    if (tb.hasKa){
      extraCmt=2;
    } else {
      extraCmt=1;
    }
    if (tb.hasDepotCmt){
      trans_syntax_error_report_fn0("cmt(depot) does not work with linCmt()");
    }
    if (tb.hasCentralCmt) {
      trans_syntax_error_report_fn0("cmt(central) does not work with linCmt()");
    }
  } else {
    if (tb.hasDepot && rx_syntax_require_ode_first){
      sprintf(buf,ODEFIRST,"depot");
      trans_syntax_error_report_fn0(buf);
    } else if (tb.hasCentral && rx_syntax_require_ode_first){
      sprintf(buf,ODEFIRST,"central");
      trans_syntax_error_report_fn0(buf);
    }
  }
  for (i=0; i<tb.nv; i++) {
    islhs = tb.lh[i];
    if (islhs>1 && islhs != 19) continue;      /* is a state var */
    if (islhs == 1 || islhs == 19){
      li++;
    } else {
      pi++;
    }
  }
  tb.pi=pi;
  tb.li=li;
  
  int pro = 0;
  SEXP lst   = PROTECT(allocVector(VECSXP, 18));pro++;
  SEXP names = PROTECT(allocVector(STRSXP, 18));pro++;

  SEXP sNeedSort = PROTECT(allocVector(INTSXP,1));pro++;
  int *iNeedSort  = INTEGER(sNeedSort);
  iNeedSort[0] = needSort;
  
  SEXP sMtime = PROTECT(allocVector(INTSXP,1));pro++;
  int *iMtime  = INTEGER(sMtime);
  iMtime[0] = (int)nmtime;
  
  SEXP tran  = PROTECT(allocVector(STRSXP, 20));pro++;
  SEXP trann = PROTECT(allocVector(STRSXP, 20));pro++;

  int offCmt=0,nExtra = 0;
  for (int i = 0; i < tb.statei; i++){
    if (offCmt == 0 && tb.idu[i] == 0){
      offCmt = 1;
      nExtra++;
      retieve_var(tb.di[i], buf);
    } else if (offCmt == 1 && tb.idu[i] == 1){
      // There is an compartment that doesn't have a derivative
      if (tb.linCmt == 0){
	UNPROTECT(pro);
	error("Compartment '%s' needs differential equations defined", buf);
      } else if (!strcmp("depot", buf) || !strcmp("central", buf)) {
      } else {
	UNPROTECT(pro);
	error("Compartment '%s' needs differential equations defined", buf);
      }
    } else if (offCmt == 1 && tb.idu[i] == 0){
      nExtra++;
    }
  }
  tb.nExtra=nExtra;

  SEXP state      = PROTECT(allocVector(STRSXP,tb.statei-tb.nExtra));pro++;
  SEXP stateRmS   = PROTECT(allocVector(INTSXP,tb.statei-tb.nExtra));pro++;
  int *stateRm    = INTEGER(stateRmS);
  SEXP extraState = PROTECT(allocVector(STRSXP,nExtra));pro++;
  
  SEXP sens     = PROTECT(allocVector(STRSXP,tb.sensi));pro++;
  SEXP normState= PROTECT(allocVector(STRSXP,tb.statei-tb.sensi-nExtra));pro++;
  
  SEXP fn_ini   = PROTECT(allocVector(STRSXP, tb.fdn));pro++;

  SEXP dfdy = PROTECT(allocVector(STRSXP,tb.ndfdy));pro++;
  
  SEXP params = PROTECT(allocVector(STRSXP, tb.pi));pro++;
  SEXP lhs    = PROTECT(allocVector(STRSXP, tb.li));pro++;

  SEXP inin  = PROTECT(allocVector(STRSXP, tb.isPi + tb.ini_i));pro++;
  SEXP ini   = PROTECT(allocVector(REALSXP, tb.isPi + tb.ini_i));pro++;

  SEXP version  = PROTECT(allocVector(STRSXP, 3));pro++;
  SEXP versionn = PROTECT(allocVector(STRSXP, 3)); pro++;
  
  SET_STRING_ELT(versionn,0,mkChar("version"));
  SET_STRING_ELT(versionn,1,mkChar("repo"));
  SET_STRING_ELT(versionn,2,mkChar("md5"));

  SET_STRING_ELT(version,0,mkChar(__VER_ver__));
  SET_STRING_ELT(version,1,mkChar(__VER_repo__));
  SET_STRING_ELT(version,2,mkChar(__VER_md5__));
  setAttrib(version,   R_NamesSymbol, versionn);

  ini_i=0;
  int redo = 0;
  for (i = 0; i < tb.nv; i++){
    retieve_var(i, buf);
    if (tb.ini[i] == 1 && tb.lh[i] != 1){
      if (tb.isPi && !strcmp("pi", buf)) {
	redo=1;
	tb.isPi=0;
	break;
      }
      SET_STRING_ELT(inin,ini_i,mkChar(buf));
      REAL(ini)[ini_i++] = tb.iniv[i];
    }
  }
  if (tb.isPi){
    SET_STRING_ELT(inin,ini_i,mkChar("pi"));
    REAL(ini)[ini_i++] = M_PI;
  } else if (redo){
    inin  = PROTECT(allocVector(STRSXP, tb.ini_i));pro++;
    ini   = PROTECT(allocVector(REALSXP, tb.ini_i));pro++;
    ini_i=0;
    for (i = 0; i < tb.nv; i++){
      retieve_var(i, buf);
      if (tb.ini[i] == 1 && tb.lh[i] != 1){
	if (tb.isPi && !strcmp("pi", buf)) {
	  redo=1;
	  tb.isPi=0;
	  break;
	}
	SET_STRING_ELT(inin,ini_i,mkChar(buf));
	REAL(ini)[ini_i++] = tb.iniv[i];
      }
    }
  }
  tb.ini_i = ini_i;

  setAttrib(ini,   R_NamesSymbol, inin);  
  
  SEXP model  = PROTECT(allocVector(STRSXP,1));pro++;
  SEXP modeln = PROTECT(allocVector(STRSXP,1));pro++;
  k=0;j=0;l=0;m=0,p=0;
  for (i=0; i<tb.nd; i++) {                     /* name state vars */
    retieve_var(tb.di[i], buf);
    if (tb.idu[i] == 1){
      if (strncmp(buf,"rx__sens_", 9) == 0){
	SET_STRING_ELT(sens,j++,mkChar(buf));
	SET_STRING_ELT(state,k++,mkChar(buf));
	stateRm[k-1]=tb.idi[i];
      } else {
	SET_STRING_ELT(normState,m++,mkChar(buf));
	SET_STRING_ELT(state,k++,mkChar(buf));
	stateRm[k-1]=tb.idi[i];
      }
      if (tb.fdi[i]){
	SET_STRING_ELT(fn_ini,l++,mkChar(buf));
      }
    } else {
      if (tb.fdi[i]){
	UNPROTECT(pro);
	error("Initialization of non-ODE compartment '%s' makes no sense", buf);
      }
      SET_STRING_ELT(extraState, p++, mkChar(buf));
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
    if (islhs>1 && islhs != 19) continue;      /* is a state var */
    retieve_var(i, buf);
    if (islhs == 1 || islhs == 19){
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
  setInits(ini);

  SET_STRING_ELT(names,0,mkChar("params"));
  SET_VECTOR_ELT(lst,  0,params);
  
  SET_STRING_ELT(names,1,mkChar("lhs"));
  SET_VECTOR_ELT(lst,  1,lhs);
  
  SET_STRING_ELT(names,2,mkChar("state"));
  SET_VECTOR_ELT(lst,  2,state);

  SET_STRING_ELT(names,3,mkChar("trans"));
  SET_VECTOR_ELT(lst,  3,tran);

  SET_STRING_ELT(names,4,mkChar("model"));
  SET_VECTOR_ELT(lst,  4,model);

  SET_STRING_ELT(names,5,mkChar("ini"));
  SET_VECTOR_ELT(lst,  5,ini);

  SET_STRING_ELT(names,6,mkChar("podo"));
  SET_VECTOR_ELT(lst,  6,ScalarLogical(rx_podo));

  SET_STRING_ELT(names,7,mkChar("dfdy"));
  SET_VECTOR_ELT(lst,  7,dfdy);

  SET_STRING_ELT(names,8,mkChar("sens"));
  SET_VECTOR_ELT(lst,  8,sens);
  
  SET_STRING_ELT(names,9,mkChar("fn.ini"));
  SET_VECTOR_ELT(lst,  9,fn_ini);

  SET_STRING_ELT(names,10,mkChar("state.ignore"));
  SET_VECTOR_ELT(lst,  10,stateRmS);

  SET_STRING_ELT(names,11,mkChar("version"));
  SET_VECTOR_ELT(lst,  11,version);

  SET_STRING_ELT(names,12,mkChar("normal.state"));
  SET_VECTOR_ELT(lst,  12,normState);
  
  SET_STRING_ELT(names,13,mkChar("needSort"));
  SET_VECTOR_ELT(lst,  13,sNeedSort);

  SET_STRING_ELT(names,14,mkChar("nMtime"));
  SET_VECTOR_ELT(lst,  14,sMtime);

  SET_STRING_ELT(names, 15, mkChar("extraCmt"));
  SEXP sExtraCmt = PROTECT(allocVector(INTSXP,1));pro++;
  INTEGER(sExtraCmt)[0] = extraCmt;
  SET_VECTOR_ELT(lst, 15, sExtraCmt);

  SET_STRING_ELT(names, 16, mkChar("stateExtra"));
  SET_VECTOR_ELT(lst,  16, extraState);

  SET_STRING_ELT(names, 17, mkChar("dvid"));
  SEXP sDvid = PROTECT(allocVector(INTSXP,tb.dvidn));pro++;
  for (i = 0; i < tb.dvidn; i++) INTEGER(sDvid)[i]=tb.dvid[i];
  SET_VECTOR_ELT(lst,  17, sDvid);

  sprintf(buf,"%.*s", (int)strlen(model_prefix)-1, model_prefix);
  SET_STRING_ELT(trann,0,mkChar("lib.name"));
  SET_STRING_ELT(tran,0,mkChar(buf));
  
  SET_STRING_ELT(trann,1,mkChar("jac"));
  if (found_jac == 1 && good_jac == 1){
    SET_STRING_ELT(tran,1,mkChar("fulluser")); // Full User Matrix
  } else {
    SET_STRING_ELT(tran,1,mkChar("fullint")); // Full Internal Matrix
  }
  
  SET_STRING_ELT(trann,2,mkChar("prefix"));
  SET_STRING_ELT(tran,2,mkChar(buf));

  sprintf(buf,"%sdydt",model_prefix);
  SET_STRING_ELT(trann,3,mkChar("dydt"));
  SET_STRING_ELT(tran,3,mkChar(buf)) ;

  sprintf(buf,"%scalc_jac",model_prefix);
  SET_STRING_ELT(trann,4,mkChar("calc_jac"));
  SET_STRING_ELT(tran, 4,mkChar(buf));

  sprintf(buf,"%scalc_lhs",model_prefix);
  SET_STRING_ELT(trann,5,mkChar("calc_lhs"));
  SET_STRING_ELT(tran, 5,mkChar(buf));

  sprintf(buf,"%smodel_vars",model_prefix);
  SET_STRING_ELT(trann,6,mkChar("model_vars"));
  SET_STRING_ELT(tran, 6,mkChar(buf));

  sprintf(buf,"%stheta",model_prefix);
  SET_STRING_ELT(trann,7,mkChar("theta"));
  SET_STRING_ELT(tran, 7,mkChar(buf));

  sprintf(buf,"%sinis",model_prefix);
  SET_STRING_ELT(trann,8,mkChar("inis"));
  SET_STRING_ELT(tran, 8,mkChar(buf));

  sprintf(buf,"%sdydt_lsoda",model_prefix);
  SET_STRING_ELT(trann,9,mkChar("dydt_lsoda"));
  SET_STRING_ELT(tran, 9,mkChar(buf));

  sprintf(buf,"%scalc_jac_lsoda",model_prefix);
  SET_STRING_ELT(trann,10,mkChar("calc_jac_lsoda"));
  SET_STRING_ELT(tran, 10,mkChar(buf));

  sprintf(buf,"%sode_solver_solvedata",model_prefix);
  SET_STRING_ELT(trann,11,mkChar("ode_solver_solvedata"));
  SET_STRING_ELT(tran, 11,mkChar(buf));
  
  sprintf(buf,"%sode_solver_get_solvedata",model_prefix);
  SET_STRING_ELT(trann,12,mkChar("ode_solver_get_solvedata"));
  SET_STRING_ELT(tran, 12,mkChar(buf));

  sprintf(buf,"%sdydt_liblsoda",model_prefix);
  SET_STRING_ELT(trann,13,mkChar("dydt_liblsoda"));
  SET_STRING_ELT(tran, 13,mkChar(buf));

  sprintf(buf,"%sF",model_prefix);
  SET_STRING_ELT(trann,14,mkChar("F"));
  SET_STRING_ELT(tran, 14,mkChar(buf));

  sprintf(buf,"%sLag",model_prefix);
  SET_STRING_ELT(trann,15,mkChar("Lag"));
  SET_STRING_ELT(tran, 15,mkChar(buf));

  sprintf(buf,"%sRate",model_prefix);
  SET_STRING_ELT(trann,16,mkChar("Rate"));
  SET_STRING_ELT(tran, 16,mkChar(buf));

  sprintf(buf,"%sDur",model_prefix);
  SET_STRING_ELT(trann,17,mkChar("Dur"));
  SET_STRING_ELT(tran, 17,mkChar(buf));

  sprintf(buf,"%smtime",model_prefix);
  SET_STRING_ELT(trann,18,mkChar("mtime"));
  SET_STRING_ELT(tran, 18,mkChar(buf));

  sprintf(buf,"%sassignFuns",model_prefix);
  SET_STRING_ELT(trann,19,mkChar("assignFuns"));
  SET_STRING_ELT(tran, 19,mkChar(buf));

  SET_STRING_ELT(modeln,0,mkChar("normModel"));
  SET_STRING_ELT(model,0,mkChar(sbNrm.s));
  
  setAttrib(tran,  R_NamesSymbol, trann);
  setAttrib(lst,   R_NamesSymbol, names);
  setAttrib(model, R_NamesSymbol, modeln);
  SEXP cls = PROTECT(allocVector(STRSXP, 1));pro++;
  SET_STRING_ELT(cls, 0, mkChar("rxModelVars"));
  classgets(lst, cls);
  
  UNPROTECT(pro);
  if (rx_syntax_error){
    if(!rx_suppress_syntax_info){
      if (gBuf[gBufLast] != '\0'){
	gBufLast++;
	REprintf("\n:%03d: ", lastSyntaxErrorLine);
	for (; gBuf[gBufLast] != '\0'; gBufLast++){
	  if (gBuf[gBufLast] == '\n'){
	    REprintf("\n:%03d: ", ++lastSyntaxErrorLine);
	  } else{
	    REprintf("%c", gBuf[gBufLast]);
	  }
	}
      }
      if (isEsc)REprintf("\n\033[1m================================================================================\033[0m\n");
      else REprintf("\n================================================================================\n");
    }
    error("Syntax Errors (see above)");
  }
  /* Free(sbPm); Free(sbNrm); */
  return lst;
}

SEXP _RxODE_parseModel(SEXP type){
  if (!sbPm.o){
    error("Model no longer loaded in memory.");
  }
  int iT = INTEGER(type)[0];
  SEXP pm;
  switch (iT){
  case 1:
    pm = PROTECT(allocVector(STRSXP, sbPmDt.n));
    for (int i = 0; i < sbPmDt.n; i++){
      SET_STRING_ELT(pm, i, mkChar(sbPmDt.line[i]));
    }
    break;
  default:
    pm = PROTECT(allocVector(STRSXP, sbPm.n));
    for (int i = 0; i < sbPm.n; i++){
      SET_STRING_ELT(pm, i, mkChar(sbPm.line[i]));
    }
    break;
  }
  UNPROTECT(1);
  return pm;
}

SEXP _RxODE_codeLoaded(){
  SEXP pm = PROTECT(allocVector(INTSXP, 1));
  if (!sbPm.o || !sbNrm.o){
    INTEGER(pm)[0]=0;
  } else {
    INTEGER(pm)[0]=1;
  }
  UNPROTECT(1);
  return pm;
}

SEXP _RxODE_isLinCmt(){
  SEXP ret = PROTECT(allocVector(INTSXP, 1));
  INTEGER(ret)[0]=tb.linCmt;
  UNPROTECT(1);
  return ret;
}

SEXP _RxODE_codegen(SEXP c_file, SEXP prefix, SEXP libname,
		    SEXP pMd5, SEXP timeId, SEXP fixInis){
  if (!sbPm.o || !sbNrm.o){
    error("Nothing in output queue to write");
  }
  if (!isString(c_file) || length(c_file) != 1){
    error("c_file should only be 1 file");
  }
  if (length(libname) != 2){
    error("libname needs 2 elements");
  }
  fpIO = fopen(CHAR(STRING_ELT(c_file,0)), "wb");
  err_msg((intptr_t) fpIO, "error opening output c file\n", -2);
  sIniTo(&sbOut, (int)((sbPm.sN)*5.3));
  // show_ode = 1 dydt
  // show_ode = 2 Jacobian
  // show_ode = 3 Ini statement
  // show_ode = 0 LHS
  // show_ode = 5 functional bioavailibility
  // show_ode = 6 functional rate
  
  gCode(1); // d/dt()
  gCode(2); // jac
  gCode(3); // ini()
  gCode(0); //
  gCode(5);
  gCode(6);
  gCode(7);
  gCode(8);
  gCode(9); // mtime
  gCode(4); // Registration
  fclose(fpIO);
  return R_NilValue;
}

char *getLine (char *src, int line, int *lloc)
{
  int cur = 1, col=0, i;
  for(i = 0; src[i] != '\0' && cur != line; i++){
    if(src[i] == '\n') cur++;
  }
  for(col = 0; src[i + col] != '\n' && src[i + col] != '\0'; col++);
  *lloc=i+col;
  char *buf = Calloc(col + 1, char);
  memcpy(buf, src + i, col);
  buf[col] = '\0';
  return buf;
}

static void rxSyntaxError(struct D_Parser *ap) {
  if (!rx_suppress_syntax_info){
    if (lastSyntaxErrorLine == 0){
      if (isEsc)REprintf("\033[1mRxODE Model Syntax Error:\n================================================================================\033[0m");
      else REprintf("RxODE Model Syntax Error:\n================================================================================");
      lastSyntaxErrorLine=1;
    }
    char *buf;
    Parser *p = (Parser *)ap;
    for (; lastSyntaxErrorLine < p->user.loc.line; lastSyntaxErrorLine++){
      buf = getLine(gBuf, lastSyntaxErrorLine, &gBufLast);
      REprintf("\n:%03d: %s", lastSyntaxErrorLine, buf);
    }
    char *after = 0;
    ZNode *z = p->snode_hash.last_all ? p->snode_hash.last_all->zns.v[0] : 0;
    while (z && z->pn->parse_node.start_loc.s == z->pn->parse_node.end)
      z = (z->sns.v && z->sns.v[0]->zns.v) ? z->sns.v[0]->zns.v[0] : 0;
    if (z && z->pn->parse_node.start_loc.s != z->pn->parse_node.end)
      after = rc_dup_str(z->pn->parse_node.start_loc.s, z->pn->parse_node.end);
    if (after){
      if (isEsc) REprintf("\n\n\033[1mRxODE syntax error after\033[0m '\033[35m\033[1m%s\033[0m':\n",  after);
      else REprintf("\n\nRxODE syntax error after '%s':\n",  after);
    }
    else{
      if (isEsc) REprintf("\n\n\033[1mRxODE syntax error\033[0m:\n");
      else REprintf("\n\nRxODE syntax error:\n");
    }

    buf = getLine(gBuf, p->user.loc.line, &gBufLast);
    if (lastSyntaxErrorLine < p->user.loc.line) lastSyntaxErrorLine++;
    if (isEsc) REprintf("\033[1m:%03d:\033[0m ", p->user.loc.line);
    else REprintf(":%03d: ", p->user.loc.line);
    /* REprintf("      "); */
    /* REprintf("%s\n", buf); */
    int col = 0, len= strlen(buf), lenv, i;
    for (i = 0; i < p->user.loc.col; i++){
      REprintf("%c", buf[i]);
      if (i == len-2) { i++; break;}
    }
    if (isEsc) REprintf("\033[35m\033[1m%c\033[0m", buf[i++]);
    else REprintf("%c", buf[i++]);
    for (; i < len; i++){
      REprintf("%c", buf[i]);
    }
    REprintf("\n      ");
    
    if (after){
      lenv = strlen(after);
      while (col != len && strncmp(buf + col, after, lenv) != 0) col++;
      if (col == len) col = 0;
      if (col){
	for (int i = 0; i < col; i++){
	  REprintf(" ");
	  if (i == len-2) { i++; break;}
	}
	len = p->user.loc.col - col;
	if (len > 0 && len < 40){
	  for (int i = len; i--;) REprintf("~");
	}
	if (isEsc) REprintf("\033[35m\033[1m^\033[0m");
	else REprintf("^");
      } else {
	for (int i = 0; i < p->user.loc.col; i++){
	  REprintf(" ");
	  if (i == len-2) { i++; break;}
	}
	if (isEsc) REprintf("\033[35m\033[1m^\033[0m");
	else REprintf("^");
      }
    } else {
      for (int i = 0; i < p->user.loc.col; i++){
	REprintf(" ");
	if (i == len-2) { i++; break;}
      }
      if (isEsc) REprintf("\033[35m\033[1m^\033[0m");
      else REprintf("^");
    }
    Free(buf);
    if (after) Free(after);
  }
  rx_syntax_error = 1;
}

static void trans_syntax_error_report_fn0(char *err){
  if (!rx_suppress_syntax_info){
    if (lastSyntaxErrorLine == 0){
      if (isEsc) REprintf("\033[1mRxODE Model Syntax Error:\n================================================================================\033[0m");
      else REprintf("RxODE Model Syntax Error:\n================================================================================");
      lastSyntaxErrorLine=1;
    }
    if (isEsc) REprintf("\n\033[1m:ERR:\033[0m %s:\n",  err);
    else REprintf("\n:ERR: %s:\n", err);
  }
  rx_syntax_error = 1;
}

static void trans_syntax_error_report_fn(char *err) {
  if (!rx_suppress_syntax_info){
    if (lastSyntaxErrorLine == 0){
      if (isEsc) REprintf("\033[1mRxODE Model Syntax Error:\n================================================================================\033[0m");
      else REprintf("RxODE Model Syntax Error:\n================================================================================");
      lastSyntaxErrorLine=1;
    }
    Parser *p = (Parser *)curP;
    char *buf;
    for (; lastSyntaxErrorLine < p->user.loc.line; lastSyntaxErrorLine++){
      buf = getLine(gBuf, lastSyntaxErrorLine, &gBufLast);
      REprintf("\n:%03d: %s", lastSyntaxErrorLine, buf);
    }
    if (lastSyntaxErrorLine < p->user.loc.line){
      REprintf("\n");
      lastSyntaxErrorLine++;
    }
    if (isEsc) REprintf("\n\033[1m:%03d:\033[0m %s:\n", p->user.loc.line, err);
    else REprintf("\n:%03d: %s:\n", p->user.loc.line, err);
    buf = getLine(gBuf, p->user.loc.line, &gBufLast);
    REprintf("      ");
    int i, len = strlen(buf);
    for (i = 0; i < p->user.loc.col; i++){
      REprintf("%c", buf[i]);
      if (i == len-2) { i++; break;}
    }
    if (isEsc) REprintf("\033[35m\033[1m%c\033[0m", buf[i++]);
    else REprintf("%c", buf[i++]);
    for (; i < len; i++){
      REprintf("%c", buf[i]);
    }
    REprintf("\n      ");
    Free(buf);
    for (int i = 0; i < p->user.loc.col; i++){
      REprintf(" ");
      if (i == len-2) { i++; break;}
    }
    if (isEsc) REprintf("\033[35m\033[1m^\033[0m");
    else REprintf("^");
    if (syntaxErrorExtra > 0 && syntaxErrorExtra < 40){
      for (int i = syntaxErrorExtra; i--;) REprintf("~");
    }
    syntaxErrorExtra=0;
  }
  rx_syntax_error = 1;
}


void updateSyntaxCol(){
  int i = lastStrLoc, lineNum=1, colNum=0;
  for(i = 0; gBuf[i] != '\0' && lastStr != gBuf + i; i++){
    if(gBuf[i] == '\n'){
      lineNum++;
      colNum=0;
    } else {
      colNum++;
    }
  }
  lastStrLoc=i;
  Parser *p = (Parser *)curP;
  p->user.loc.line=lineNum;
  p->user.loc.col=colNum;
}
