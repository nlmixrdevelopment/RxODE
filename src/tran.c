#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>   /* dj: import intptr_t */
#include "ode.h"
#include "sbuf.h"
#include "getOption.h"
#include "parseLinCmt.h"
#include "tran.h"
#include <errno.h>
#include <dparser.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif
#include "tran.g.d_parser.c"
#include "../inst/include/RxODE.h"
#define max(a,b) (a)>(b) ? (a):(b)
#define min(a,b) (a)<(b) ? (a):(b)
#define MXSYM 50000
#define MXDER 5000
#define MXLEN 12000
/* #define SBUF_MXBUF 5 */
#define MXLINE 100
/* #define MXLINE 5 */
#define SBPTR sb.s+sb.o
#define SBTPTR sbt.s+sbt.o
#define NV tb.ss.n

#define STRINGIFY(...) STRINGIFY_AUX(__VA_ARGS__)
#define STRINGIFY_AUX(...) #__VA_ARGS__

#define ENDLINE tb.ixL=-1; tb.didEq=0; tb.NEnd=NV;

#define gCode(i) (&sbOut)->s[0]='\0';		\
  (&sbOut)->o=0;				\
  codegen(gBuf, i, CHAR(STRING_ELT(prefix,0)),	\
	  CHAR(STRING_ELT(libname, 0)),		\
	  CHAR(STRING_ELT(pMd5,0)),		\
	  CHAR(STRING_ELT(timeId, 0)),		\
	  CHAR(STRING_ELT(libname, 1)));					\
  writeSb(&sbOut, fpIO);

#define aAppendN(str, len) sAppendN(&sb, str, len); sAppendN(&sbDt, str, len);
#define aProp(prop) curLineProp(&sbPm, prop); curLineProp(&sbPmDt, prop); curLineProp(&sbNrmL, prop);
#define aType(type) curLineType(&sbPm, type); curLineType(&sbPmDt, type); curLineType(&sbNrmL, type);

#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#if (__STDC_VERSION__ >= 199901L)
#include <stdint.h>
#endif

void RSprintf(const char *format, ...);

// from mkdparse_tree.h
typedef void (print_node_fn_t)(int depth, char *token_name, char *token_value, void *client_data);

int syntaxErrorExtra = 0;
int isEsc=0;
const char *lastStr;
int lastStrLoc=0;

SEXP _goodFuns;
char * rc_dup_str(const char *s, const char *e);
vLines _dupStrs;

int rx_syntax_error = 0, rx_suppress_syntax_info=0, rx_podo = 0, rx_syntax_require_ode_first = 1;

extern D_ParserTables parser_tables_RxODE;

unsigned int found_jac = 0, nmtime=0;
int rx_syntax_assign = 0, rx_syntax_star_pow = 0,
  rx_syntax_require_semicolon = 0, rx_syntax_allow_dots = 0,
  rx_syntax_allow_ini0 = 1, rx_syntax_allow_ini = 1, rx_syntax_allow_assign_state = 0,
  maxSumProdN = 0, SumProdLD = 0, good_jac=1, extraCmt=0, gnini=0;

sbuf s_aux_info;
sbuf s_inits;

#define notLHS 0
#define isLHS 1
#define isState 9
#define isSuppressedLHS 10
#define isSuppressedParam 11
#define isLhsStateExtra 19
#define isLHSparam 70

symtab tb;

static inline void addSymbolStr(char *value) {
  addLine(&(tb.ss),"%s",value);
  if (tb.depotN == -1 && !strcmp("depot", value)) {
    tb.depotN = NV-1;
  } else if (tb.centralN && !strcmp("central", value)){
    tb.centralN = NV-1;
  }
}


sbuf sb, sbDt; /* buffer w/ current parsed & translated line */
sbuf sbt;

sbuf firstErr;
int firstErrD=0;

vLines sbPm, sbPmDt, sbNrmL;
sbuf sbNrm;
vLines depotLines, centralLines;

const char *model_prefix = NULL;
const char *me_code = NULL;
const char *md5 = NULL;
int badMd5 = 0;
int foundF=0,foundLag=0, foundRate=0, foundDur=0, foundF0=0, needSort=0;

sbuf sbOut;

static FILE *fpIO;

int lastSyntaxErrorLine=0;
void trans_syntax_error_report_fn(char *err);
static void trans_syntax_error_report_fn0(char *err);
char *getLine (char *src, int line, int *lloc);
void updateSyntaxCol();

#include "parseVars.h"

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
  int function_name;
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
  int selection_statement__9;
  int break_statement;
  int sign;
  int sum;
  int theta0;
  int theta0_noout;
  int theta;
  int cmt_statement;
  int param_statement;
  int dvid_statementI;
  int ifelse;
  int ifelse_statement;
  int mat0;
  int matF;
  int equality_str1;
  int equality_str2;
  int simfun_statement;
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
  ni->function_name=-1;
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
  ni->selection_statement__9 = -1;
  ni->break_statement = -1;
  ni->sign = -1;
  ni->sum = -1;
  ni->theta = -1;
  ni->theta0 = -1;
  ni->theta0_noout = -1;
  ni->cmt_statement = -1;
  ni->param_statement = -1;
  ni->dvid_statementI = -1;
  ni->ifelse = -1;
  ni->ifelse_statement=-1;
  ni->mat0 = -1;
  ni->matF = -1;
  ni->equality_str1 = -1;
  ni->equality_str2 = -1;
  ni->simfun_statement = -1;
}

static inline int new_de(const char *s){
  int i;
  if (!strcmp("cmt", s)) Rf_errorcall(R_NilValue, _("'cmt' cannot be a state or lhs expression"));
  if (!strcmp("dvid", s)) Rf_errorcall(R_NilValue, _("'dvid' cannot be a state or lhs expression"));
  if (!strcmp("addl", s)) Rf_errorcall(R_NilValue, _("'addl' cannot be a state or lhs expression"));
  if (!strcmp("ii", s)) Rf_errorcall(R_NilValue, _("'ii' cannot be a state or lhs expression"));
  if (!strcmp("ss", s)) Rf_errorcall(R_NilValue, _("'ss' cannot be a state or lhs expression"));
  if (!strcmp("amt", s)) Rf_errorcall(R_NilValue, _("'amt' cannot be a state or lhs expression"));
  if (!strcmp("dur", s)) Rf_errorcall(R_NilValue, _("'dur' cannot be a state or lhs expression"));
  if (!strcmp("rate", s)) Rf_errorcall(R_NilValue, _("'rate' cannot be a state or lhs expression"));
  if (!strcmp("Rprintf", s)) Rf_errorcall(R_NilValue, _("'Rprintf' cannot be a state"));
  if (!strcmp("printf", s)) Rf_errorcall(R_NilValue, _("'printf' cannot be a state"));
  if (!strcmp("print", s)) Rf_errorcall(R_NilValue, _("'print' cannot be a state"));
  for (i=0; i<tb.de.n; i++) {
    if (!strcmp(tb.de.line[i], s)) {
      tb.id = i;
      return 0;
    }
  }
  if (tb.de.n + 1 > tb.allocD){
    tb.allocD+=MXDER;
    tb.di=Realloc(tb.di, tb.allocD, int);
    tb.idi=Realloc(tb.idi, tb.allocD, int);
    tb.idu=Realloc(tb.idu, tb.allocD, int);
    tb.dvid=Realloc(tb.dvid, tb.allocD, int);
  }
  return 1;
}

static inline int nodeTime(char *value) {
  if (!strcmp("time",value)){
    aAppendN("t", 1);
    sAppendN(&sbt, "t", 1);
    return 1;
  }
  return 0;
}

static inline int nodePodo(char *value) {
  if (!strcmp("podo",value)){
    aAppendN("_solveData->subjects[_cSub].podo", 32);
    sAppendN(&sbt, "podo", 4);
    rx_podo = 1;
    return 1;
  }
  return 0;
}

static inline int nodeCmt(char *value) {
  if (!strcmp("CMT",value)){
    aAppendN("_CMT", 4);
    sAppendN(&sbt, "CMT", 3);
    return 1;
  }
  return 0;
}

static inline int nodeTlast(char *value) {
  if (!strcmp("tlast",value)){
    aAppendN("_solveData->subjects[_cSub].tlast", 33);
    sAppendN(&sbt, "tlast", 5);
    return 1;
  }
  return 0;
}

static inline int nodePtr(char *value) {
  if (!strcmp("rx__PTR__",value)){
    aAppendN("_solveData, _cSub", 17);
    sAppendN(&sbt, "rx__PTR__", 9);
    return 1;
  }
  return 0;
}

static inline int nodeNaN(char *value){
  if (!strcmp("NaN",value)){
    aAppendN("NAN", 3);
    sAppendN(&sbt,"NaN", 3);
    return 1;
  }
  return 0;
}

static inline int nodeNA(char *value) {
  if (!strcmp("NA",value)){
    aAppendN("NA_REAL", 7);
    sAppendN(&sbt,"NA", 2);
    return 1;
  }
  return 0;
}

static inline int nodeInf(char *value) {
  if (!strcmp("Inf",value)){
    if (sbt.o > 0 && sbt.s[sbt.o-1] == '-'){
      sb.o--; sbDt.o--;
      aAppendN("R_NegInf", 8);
    } else {
      aAppendN("R_PosInf", 8);
    }
    sAppendN(&sbt,"Inf", 3);
    return 1;
  }
  return 0;
}

static inline int nodeFunGamma(char *value) {
  if (!strcmp("gamma",value)){
    aAppendN("lgammafn", 8);
    sAppendN(&sbt, "lgammafn", 8);
    return 1;
  }
  return 0;
}

static inline int nodeFunLfactorial(char *value) {
  if (!strcmp("lfactorial",value)){
    aAppendN("lgamma1p", 8);
    sAppendN(&sbt, "lgamma1p", 8);
    return 1;
  }
  return 0;
}

static inline int nodeFunLog(char *value) {
  if (!strcmp("log",value)){
    aAppendN("_safe_log", 9);
    sAppendN(&sbt, "log", 3);
    return 1;
  }
  return 0;
}

static inline int nodeFunAbs(char *value) {
  if (!strcmp("abs",value)){
    aAppendN("fabs", 4);
    sAppendN(&sbt,"abs", 3);
    return 1;
  }
  return 0;
}

static inline int nodeFunLinCmt(char *value) {
  if (!strcmp("linCmt",value)) {
    if (tb.linCmt == 0){
      aAppendN("linCmt", 6);
      aProp(-100);
      sAppendN(&sbt,"linCmt", 6);
      tb.linCmt=1;
    } else {
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("only one 'linCmt()' per model"));
    }
    return 1;
  }
  return 0;
}

static inline int nodeFunLinCmtA(char *value) {
  if (!strcmp("linCmtA",value)){
    aAppendN("linCmtA", 7);
    sAppendN(&sbt,"linCmtA", 7);
    tb.linCmt=2;
    return 1;
  }
  return 0;
}

static inline int nodeFunLinCmtB(char *value) {
  if (!strcmp("linCmtB",value)){
    aAppendN("linCmtB", 7);
    sAppendN(&sbt,"linCmtB", 7);
    tb.linCmt=2;
    return 1;
  }
  return 0;
}

static inline int nodeFunLinCmtC(char *value){
  if (!strcmp("linCmtC",value)){
    aAppendN("linCmtC", 7);
    sAppendN(&sbt,"linCmtC", 7);
    tb.linCmt=2;
    return 1;
  }
  return 0;
}

void wprint_node(int depth, char *name, char *value, void *client_data) {
  int i;
  nodeInfo ni;
  niReset(&ni);
  int tmp = nodeTime(value) ||
    nodePodo(value) ||
    nodeCmt(value) ||
    nodeTlast(value) ||
    nodePtr(value) ||
    nodeNaN(value) ||
    nodeNA(value) ||
    nodeInf(value);
  if (!tmp && nodeHas(identifier)) {
    tmp = nodeFunGamma(value) ||
      nodeFunLfactorial(value) ||
      nodeFunLog(value) ||
      nodeFunAbs(value) ||
      nodeFunLinCmt(value) ||
      nodeFunLinCmtA(value) ||
      nodeFunLinCmtB(value) ||
      nodeFunLinCmtC(value);
  }
  if (!tmp) {
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
int gBufFree=0;
int gBufLast;
D_Parser *curP=NULL;
D_ParseNode *_pn = 0;

void freeP(){
  if (_pn){
    free_D_ParseTreeBelow(curP,_pn);
    free_D_ParseNode(curP,_pn);
  }
  _pn=0;
  if (curP != NULL){
    free_D_Parser(curP);
  }
  curP = NULL;
}

int toInt(char *v2){
  errno = 0;
  char *v3 = v2;
  char *endptr = NULL;
  long lagNoL = strtol(v3, &endptr, 10);
  int lagNo;
  if (errno == 0 && v3 && !*endptr){
    lagNo = (int)(lagNoL);
  } else {
    lagNo = NA_INTEGER;
  }
  errno=0;
  return lagNo;
}

int skipDouble=0;

int allSpaces(char *v2) {
  int iii=0;
  int allSpace=1;
  while(v2[iii] != '\0'){
    if (!isspace(v2[iii])){
      allSpace=0;
      break;
    }
  }
  return allSpace;
}

int depotAttr=0, centralAttr=0;

sbuf _gbuf, _mv;

#include "parseFuns.h"

////////////////////////////////////////////////////////////////////////////////
// parsing properties (logical expressions)
static inline int isIdentifier(nodeInfo ni, const char *name) {
  return nodeHas(identifier) || nodeHas(identifier_r) ||
    nodeHas(identifier_r_no_output)  ||
    nodeHas(theta0_noout) ||
    nodeHas(theta0);
}

static inline int isTbsVar(const char *value) {
  return !strcmp("rx_lambda_", value) || !strcmp("rx_yj_", value) ||
    !strcmp("rx_low_", value) || !strcmp("rx_hi_", value);
}

static inline int isDefiningParameterRecursively(const char *value) {
  return tb.ix == tb.ixL && tb.didEq==1 &&
    !strcmp(value, tb.ss.line[tb.ix]);
}

static inline int isOperatorOrPrintingIdentifier(nodeInfo ni, const char *name){
  return nodeHas(identifier) ||
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
    !strcmp(">", name);
}

static inline int isSkipChild(nodeInfo ni, const char *name, int i) {
  return ((i == 3 || i == 4 || i < 2) &&
	  (nodeHas(derivative) ||nodeHas(fbio) || nodeHas(alag) ||
	   nodeHas(rate) || nodeHas(dur))) ||
    ((i == 3 || i < 2) && nodeHas(der_rhs)) ||
    (nodeHas(dfdy)     && i< 2)  ||
    (nodeHas(dfdy_rhs) && i< 2) ||
    (nodeHas(dfdy)     && i == 3) ||
    (nodeHas(dfdy_rhs) && i == 3) ||
    (nodeHas(dfdy)     && i == 5) ||
    (nodeHas(dfdy_rhs) && i == 5) ||
    (nodeHas(dfdy)     && i == 6) ||
    (nodeHas(ini0)     && i == 1) ||
    (nodeHas(dvid_statementI) && i != 0) ||
    ((nodeHas(theta) || nodeHas(eta)) && i != 2) ||
    (nodeHas(mtime) && (i == 0 || i == 1 || i == 3)) ||
    (nodeHas(cmt_statement) && (i == 0 || i == 1 || i == 3)) ||
    (i != 2 && (nodeHas(mat0) || nodeHas(matF)));
}

////////////////////////////////////////////////////////////////////////////////
// Parsing pieces
static inline void handleIdentifier(nodeInfo ni, char *name, char *value) {
  // Handles identifiers, add it as a symbol if needed.
  if (isIdentifier(ni, name)) {
    if (new_or_ith(value)){
      // If it is new, add it
      addSymbolStr(value);
      // Ignored variables
      if (isTbsVar(value)){
	// If it is Transform both sides, suppress printouts
	tb.lh[NV-1] = isSuppressedParam; // Suppress param printout.
      }
    } else if (isDefiningParameterRecursively(value)){
      // This is x = x*exp(matt)
      // lhs defined in terms of a parameter
      if (tb.lh[tb.ix] == isSuppressedLHS){
	tb.lh[tb.ix] = notLHS;
      } else {
	tb.lh[tb.ix] = isLHSparam;
      }
    }
  }
}

static inline void handleOperatorsOrPrintingIdentifiers(int depth, print_node_fn_t fn, void *client_data,
							nodeInfo ni, char *name, char *value) {
  if (isOperatorOrPrintingIdentifier(ni, name))
    fn(depth, name, value, client_data);
  if (!strcmp("=", name)){
    tb.didEq=1;
    fn(depth, name, value, client_data);
  }

  // Operator synonyms
  if (!strcmp("<-",name)){
    aAppendN(" =", 2);
    sAppendN(&sbt, "=", 1);
    tb.didEq=1;
  } else if (!strcmp("~",name)){
    // Suppress LHS calculation with ~
    aAppendN(" =", 2);
    sAppendN(&sbt, "~", 1);
    tb.lh[tb.ix] = isSuppressedLHS; // Suppress LHS printout.
    tb.didEq=1;
  } else if (!strcmp("=", name)){
    tb.didEq=1;
  } else if (!strcmp("|",name)){
    aAppendN(" ||", 3);
    sAppendN(&sbt, "||", 2);
  } else if (!strcmp("&",name)){
    aAppendN(" &&", 3);
    sAppendN(&sbt, "&&", 2);
  }
}

static inline void handleIndLinMat0(nodeInfo ni, char *name) {
  if (nodeHas(mat0)){
    aType(TMAT0);
    sb.o =0; sbDt.o =0;
    sAppend(&sb, "_mat[%d] = ", tb.matn++);
  }
}

static inline void handleIndLinMatf(nodeInfo ni, char *name) {
  if (nodeHas(matF)){
    aType(TMATF);
    sb.o =0; sbDt.o =0;
    sAppend(&sb, "_matf[%d] = _IR[%d] + ", tb.matnf, tb.matnf);
    tb.matnf++;
  }
}

static inline int handleIfElse(nodeInfo ni, char *name, int i) {
  if (nodeHas(ifelse)){
    if (i == 0){
      return 1;
    } else if (i == 1){
      aAppendN("((", 2);
      sAppendN(&sbt,"ifelse(", 7);
      return 1;
    } else if (i == 3){
      aAppendN(") ? (", 5);
      sAppendN(&sbt,",", 1);
      return 1;
    } else if (i == 5){
      aAppendN(") : (", 5);
      sAppendN(&sbt,",", 1);
      return 1;
    } else if (i == 7){
      aAppendN("))", 2);
      sAppendN(&sbt,")", 1);
      return 1;
    }
  }
  if (nodeHas(ifelse_statement)){
    if (i == 0){
      return 1;
    } else if (i == 1){
      aAppendN("if (", 4);
      sAppendN(&sbt, "if (", 4);
      return 1;
    } else if (i == 3){
      aType(TLOGIC);
      aAppendN(") {", 3);
      sAppendN(&sbt,") {", 3);
      addLine(&sbPm, "%s\n", sb.s);
      addLine(&sbPmDt, "%s\n", sbDt.s);
      sAppend(&sbNrm, "%s\n", sbt.s);
      addLine(&sbNrmL, "%s\n", sbt.s);
      sb.o=0;sbDt.o=0; sbt.o=0;
      return 1;
    } else if (i == 5){
      sb.o=0;sbDt.o=0; sbt.o=0;
      aType(TLOGIC);
      aAppendN("}\nelse {", 8);
      sAppendN(&sbt,"}\nelse {", 1);
      addLine(&sbPm, "%s\n", sb.s);
      addLine(&sbPmDt, "%s\n", sbDt.s);
      sAppend(&sbNrm, "%s\n", sbt.s);
      addLine(&sbNrmL, "%s\n", sbt.s);
      return 1;
    } else if (i == 7){
      sb.o=0;sbDt.o=0; sbt.o=0;
      aType(TLOGIC);
      aAppendN("}", 1);
      sAppendN(&sbt,"}", 1);
      addLine(&sbPm, "%s\n", sb.s);
      addLine(&sbPmDt, "%s\n", sbDt.s);
      sAppend(&sbNrm, "%s\n", sbt.s);
      addLine(&sbNrmL, "%s\n", sbt.s);
      return 1;
    } else if (i == 8){
      return 1;
    }
  }
  return 0;
}

static inline int handleStringEqualRhs(nodeInfo ni, char *name, int i, D_ParseNode *xpn) {
  if (nodeHas(equality_str1)){
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    switch(i) {
    case 0:
      // string
      aAppendN("_cmp1(", 6);
      sAppend(&sb, "%s, ", v);
      sAppend(&sbDt, "%s, ", v);
      sAppend(&sbt, "%s", v);
      /* Free(v); */
      return 1;
    case 1:
      if (!strcmp(v, "==")) {
	aAppendN("1, ", 3);
      } else {
	aAppendN("0, ", 3);
      }
      sAppend(&sbt, "%s", v);
      /* Free(v); */
      return 1;
    case 2:
      // identifier_r
      // val, valstr
      if (!strcmp(v, "id") || !strcmp(v, "ID") || !strcmp(v, "Id")){
	aAppendN("(&_solveData->subjects[_cSub])->idReal, \"ID\")", 45);
	sAppendN(&sbt, "ID", 2);
      } else {
	if (new_or_ith(v)) addSymbolStr(v);
	sAppend(&sb, "%s, \"%s\")", v, v);
	sAppend(&sbDt, "%s, \"%s\")", v, v);
	sAppend(&sbt, "%s", v);
      }
      return 1;
    }
  }
  return 0;
}

static inline int handleStringEqualLhs(nodeInfo ni, char *name, int i, D_ParseNode *xpn) {
  if (nodeHas(equality_str2)){
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    switch(i) {
    case 0:
      aAppendN("_cmp2(", 6);
      if (!strcmp(v, "id") || !strcmp(v, "ID") || !strcmp(v, "Id")){
	aAppendN("(&_solveData->subjects[_cSub])->idReal, \"ID\", ", 46);
	sAppendN(&sbt, "ID", 2);
      } else {
	if (new_or_ith(v)) addSymbolStr(v);
	sAppend(&sb, "%s, \"%s\", ", v, v);
	sAppend(&sbDt, "%s, \"%s\", ", v, v);
	sAppend(&sbt, "%s", v);
      }
      return 1;
    case 1:
      if (!strcmp(v, "==")) {
	aAppendN("1, ", 3);
      } else {
	aAppendN("0, ", 3);
      }
      sAppend(&sbt, "%s", v);
      return 1;
    case 2:
      sAppend(&sb, "%s)", v);
      sAppend(&sbDt, "%s)", v);
      sAppend(&sbt, "%s", v);
      return 1;
    }
  }
  return 0;
}

static inline int handleStringEqualityStatements(nodeInfo ni, char *name, int i, D_ParseNode *xpn) {
  return handleStringEqualRhs(ni, name, i, xpn) ||
    handleStringEqualLhs(ni, name, i, xpn);
}

static inline int handleDvidStatement(nodeInfo ni, char *name, D_ParseNode *xpn, D_ParseNode *pn) {
  if (nodeHas(dvid_statementI)){
    if (tb.dvidn == 0){
      // dvid->cmt translation
      sb.o=0;sbDt.o=0; sbt.o=0;
      xpn = d_get_child(pn,2);
      char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
      tb.dvid[0]=atoi(v);
      /* Free(v); */
      if (tb.dvid[0] == 0){
	updateSyntaxCol();
	trans_syntax_error_report_fn(ZERODVID);
      }
      sAppend(&sbt, "dvid(%d", tb.dvid[0]);
      xpn = d_get_child(pn,3);
      tb.dvidn = d_get_number_of_children(xpn)+1;
      D_ParseNode *xpn2;
      for (int i = 0; i < tb.dvidn-1; i++){
	xpn2 = d_get_child(xpn, i);
	v = (char*)rc_dup_str(xpn2->start_loc.s, xpn2->end);
	tb.dvid[i+1]=atoi(v+1);
	if (tb.dvid[i+1] == 0){
	  /* Free(v); */
	  updateSyntaxCol();
	  trans_syntax_error_report_fn(ZERODVID);
	}
	sAppend(&sbt, ",%d", tb.dvid[i+1]);
	/* Free(v); */
      }
      sAppend(&sbNrm, "%s);\n", sbt.s);
      addLine(&sbNrmL, "%s);\n", sbt.s);
      /* Free(v); */
      return 1;
    } else {
      updateSyntaxCol();
      trans_syntax_error_report_fn(ZERODVID);
    }
    return 1;
  }
  return 0;
}

static inline int handleTheta(nodeInfo ni, char *name, D_ParseNode *xpn) {
  if (nodeHas(theta)){
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    sPrint(&_gbuf,"_THETA_%s_",v);
    int ii = strtoimax(v,NULL,10);
    if (ii > tb.maxtheta){
      tb.maxtheta =ii;
    }
    if (new_or_ith(_gbuf.s)){
      addSymbolStr(_gbuf.s);
    }
    sAppend(&sb,"_THETA_%s_",v);
    sAppend(&sbDt,"_THETA_%s_",v);
    sAppend(&sbt,"THETA[%s]",v);
    /* Free(v); */
    return 1;
  }
  return 0;
}

static inline int handleEta(nodeInfo ni, char *name, D_ParseNode *xpn) {
  if (nodeHas(eta)){
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    int ii = strtoimax(v,NULL,10);
    if (ii > tb.maxeta){
      tb.maxeta =ii;
    }
    sPrint(&_gbuf,"_ETA_%s_",v);
    if (new_or_ith(_gbuf.s)){
      addSymbolStr(_gbuf.s);
    }
    sAppend(&sb, "_ETA_%s_",v);
    sAppend(&sbDt, "_ETA_%s_",v);
    sAppend(&sbt,"ETA[%s]",v);
    /* Free(v); */
    return 1;
  }
  return 0;
}

static inline void handleSafeZero(nodeInfo ni, char *name, int i, int *safe_zero, D_ParseNode *xpn) {
  if (nodeHas(mult_part)){
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    if (i == 0){
      if (!strcmp("/",v)){
	aAppendN("safe_zero(", 10);
	*safe_zero = 1;
      } else {
	*safe_zero = 0;
      }
    }
    if (i == 1){
      if (*safe_zero){
	aAppendN(")", 1);
      }
      *safe_zero = 0;
    }
    /* Free(v); */
  }
}
#include "parseDfdy.h"

static inline int assertLogicalNoWhileElse(nodeInfo ni, char *name, int i, D_ParseNode *pn, D_ParseNode *xpn, int *isWhile) {
  if (nodeHas(selection_statement) && i== 0 ) {
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    *isWhile = !strcmp("while", v);
    /* Free(v); */
    if (*isWhile) {
      D_ParseNode *xpn2 = d_get_child(pn, 5);
      v = (char*)rc_dup_str(xpn2->start_loc.s, xpn2->end);
      if (v[0] == 0) {
      } else {
	updateSyntaxCol();
	trans_syntax_error_report_fn(_("'while' cannot be followed by 'else' (did you mean 'if'/'else')"));
      }
    }
    return 1;
  }
  return 0;
}

static inline int handleLogicalIfOrWhile(nodeInfo ni, char *name, int i, D_ParseNode *pn, D_ParseNode *xpn, int *isWhile) {
  if (nodeHas(selection_statement) && i==1) {
    sb.o = 0; sbDt.o = 0; sbt.o = 0;
    if (*isWhile) {
      sAppendN(&sb, "_itwhile=0;\nwhile (", 19);
      sAppendN(&sbDt, "_itwhile=0;\nwhile (", 19);
      sAppendN(&sbt,"while (", 7);
      tb.nwhile++;
    } else {
      sAppendN(&sb, "if (", 4);
      sAppendN(&sbDt, "if (", 4);
      sAppendN(&sbt,"if (", 4);
    }
    return 1;
  }
  return 0;
}

static inline int handleLogicalBreak(nodeInfo ni, char *name, int i, D_ParseNode *pn, D_ParseNode *xpn, int *isWhile) {
  if (nodeHas(break_statement) && i == 0) {
    if (tb.nwhile > 0) {
      aType(TLOGIC);
      sb.o = 0; sbDt.o = 0; sbt.o = 0;
      /* aType(100); */
      aAppendN("break;", 6);
      sAppendN(&sbt, "break;", 6);
      addLine(&sbPm, "%s\n", sb.s);
      addLine(&sbPmDt, "%s\n", sbDt.s);
      sAppend(&sbNrm, "%s\n", sbt.s);
      addLine(&sbNrmL, "%s\n", sbt.s);
      ENDLINE;
    } else {
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'break' can only be used in  'while' statement"));
    }
    return 1;
  }
  return 0;
}

static inline int handleLogicalBeginParen(nodeInfo ni, char *name, int i, D_ParseNode *pn, D_ParseNode *xpn, int *isWhile) {
  if (nodeHas(selection_statement) && i==3) {
    aType(TLOGIC);
    /* aType(100); */
    aAppendN("{", 1);
    sAppendN(&sbt, "{", 1);
    addLine(&sbPm, "%s\n", sb.s);
    addLine(&sbPmDt, "%s\n", sbDt.s);
    sAppend(&sbNrm, "%s\n", sbt.s);
    addLine(&sbNrmL, "%s\n", sbt.s);
    ENDLINE;
    return 1;
  }
  return 0;
}

static inline int handleLogicalElse(nodeInfo ni, char *name, int i, D_ParseNode *pn, D_ParseNode *xpn, int *isWhile) {
  if (nodeHas(selection_statement__9) && i==0) {
    sb.o = 0; sbDt.o = 0; sbt.o = 0;
    aType(TLOGIC);
    aAppendN("}\nelse {", 8);
    sAppendN(&sbt,"}\nelse {", 8);
    addLine(&sbPm, "%s\n", sb.s);
    addLine(&sbPmDt, "%s\n", sbDt.s);
    sAppend(&sbNrm, "%s\n", sbt.s);
    addLine(&sbNrmL, "%s\n", sbt.s);
    ENDLINE;
    return 1;
  }
  return 0;
}

static inline int handleLogicalExpr(nodeInfo ni, char *name, int i, D_ParseNode *pn, D_ParseNode *xpn, int *isWhile) {
  int tmp = assertLogicalNoWhileElse(ni, name, i, pn, xpn, isWhile) ||
    handleLogicalIfOrWhile(ni, name, i, pn, xpn, isWhile) ||
    handleLogicalBreak(ni, name, i, pn, xpn, isWhile) ||
    handleLogicalBeginParen(ni, name, i, pn, xpn, isWhile) ||
    handleLogicalElse(ni, name, i, pn, xpn, isWhile);
  (void)tmp;
  return 0;
}

static inline int handleCmtPropertyFbio(nodeInfo ni, char *name, char *v) {
  if (nodeHas(fbio)){
    sb.o=0;sbDt.o=0; sbt.o=0;
    sAppend(&sb, "_f[%d] = ", tb.id);
    sAppend(&sbDt, "_f[%d] = ", tb.id);
    sAppend(&sbt, "f(%s)=", v);
    tb.curPropN=tb.id;
    if (foundF == 0) needSort+=1;// & 1 when F
    foundF=1;
    aType(FBIO);
    return 1;
  }
  return 0;
}

static inline int handleCmtPropertyAlag(nodeInfo ni, char *name, char *v) {
  if (nodeHas(alag)){
    sb.o=0;sbDt.o=0; sbt.o=0;
    sAppend(&sb, "_alag[%d] = ", tb.id);
    sAppend(&sbDt, "_alag[%d] = ", tb.id);
    sAppend(&sbt, "alag(%s)=", v);
    tb.curPropN=tb.id;
    if (foundLag == 0) needSort+=2; // & 2 when alag
    foundLag=1;
    aType(ALAG);
    return 1;
  }
  return 0;
}

static inline int handleCmtPropertyDur(nodeInfo ni, char *name, char *v) {
  if (nodeHas(dur)) {
    sb.o=0;sbDt.o=0; sbt.o=0;
    sAppend(&sb, "_dur[%d] = ", tb.id);
    sAppend(&sbDt, "_dur[%d] = ", tb.id);
    sAppend(&sbt, "dur(%s)=", v);
    tb.curPropN=tb.id;
    if (foundDur == 0) needSort+=4;// & 4 when dur
    foundDur=1;
    aType(DUR);
    return 1;
  }
  return 0;
}

static inline int handleCmtPropertyRate(nodeInfo ni, char *name, char *v) {
  if (nodeHas(rate)){
    sb.o=0;sbDt.o=0; sbt.o=0;
    sAppend(&sb, "_rate[%d] = ", tb.id);
    sAppend(&sbDt, "_rate[%d] = ", tb.id);
    sAppend(&sbt, "rate(%s)=", v);
    tb.curPropN=tb.id;
    if (foundRate == 0) needSort+=8;// & 8 when rate
    foundRate=1;
    aType(RATE);
    return 1;
  }
  return 0;
}

static inline int handleCmtPropertyCmtOrder(nodeInfo ni, char *name, char *v) {
  if (nodeHas(cmt_statement)){
    sb.o=0;sbDt.o=0; sbt.o=0;
    sAppend(&sbt, "cmt(%s)", v);
    sAppend(&sbNrm, "%s;\n", sbt.s);
    addLine(&sbNrmL, "%s;\n", sbt.s);
    return 1;
  }
  return 0;
}

static inline int isCmtLhsStatement(nodeInfo ni, char *name, char *v) {
  int hasLhs = 0;
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
  return hasLhs;
}

static inline int add_deCmtProp(nodeInfo ni, char *name, char *v, int hasLhs, int fromWhere) {
  if (hasLhs == fromCMTprop) { // 1 only
    if (tb.lh[tb.ix] == isSuppressedLHS || tb.lh[tb.ix] == 29) {
      tb.lh[tb.ix] = 29;
    } else {
      tb.lh[tb.ix] = isLhsStateExtra;
    }
    new_or_ith(v);
    return 1;
  }
  return 0;
}

static inline int add_deState(nodeInfo ni, char *name, char *v, int hasLhs, int fromWhere) {
  if (fromWhere == fromDDT && strncmp(v, "rx__sens_", 3) == 0) {
    tb.sensi++;
  }
  if (rx_syntax_allow_dots == 0 && strstr(v, ".")) {
    updateSyntaxCol();
    trans_syntax_error_report_fn(NODOT);
  }
  new_or_ith(v);
  if (!rx_syntax_allow_assign_state &&
      ((tb.ini[tb.ix] == 1 && tb.ini0[tb.ix] == 0) ||
       (tb.lh[tb.ix] == isLHS || tb.lh[tb.ix] == isLHSparam))){
    updateSyntaxCol();
    sPrint(&_gbuf,_("cannot assign state variable %s; For initial condition assignment use '%s(0) = #'.\n  Changing states can break sensitivity analysis (for nlmixr glmm/focei).\n  To override this behavior set 'options(RxODE.syntax.assign.state = TRUE)'"),v,v);
    trans_syntax_error_report_fn0(_gbuf.s);
  }
  tb.lh[tb.ix] = isState;
  return 1;
}

static inline void add_de(nodeInfo ni, char *name, char *v, int hasLhs, int fromWhere) {
  tb.statei++;
  tb.id=tb.de.n;
  if (fromWhere == fromCMTprop && !nodeHas(cmt_statement)) {
    if (rx_syntax_require_ode_first) {
      if (!strcmp("depot", v)) {
	tb.hasDepot = 1;
      } else if (!strcmp("central", v)) {
	tb.hasCentral = 1;
      } else {
	updateSyntaxCol();
	sPrint(&_gbuf,ODEFIRST,v);
	trans_syntax_error_report_fn(_gbuf.s);
      }
    }
  }
  int tmp = add_deCmtProp(ni, name, v, hasLhs, fromWhere) ||
    add_deState(ni, name, v, hasLhs, fromWhere);
  (void) tmp;
  tb.di[tb.de.n] = tb.ix;
  addLine(&(tb.de),"%s",v);
}

static inline int handleCmtProperty(nodeInfo ni, char *name, int i, D_ParseNode *xpn) {
  if ((nodeHas(fbio) || nodeHas(alag) ||
       nodeHas(dur) || nodeHas(rate) ||
       nodeHas(cmt_statement)) && i==2) {
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    int hasLhs=isCmtLhsStatement(ni, name, v);
    if (new_de(v)){
      add_de(ni, name, v, hasLhs, fromCMTprop);
      aProp(tb.de.n);
      handleCmtPropertyCmtOrder(ni, name, v);
    } else {
      new_or_ith(v);
      aProp(tb.ix);
      /* printf("de[%d]->%s[%d]\n",tb.id,v,tb.ix); */
    }
    int tmp = handleCmtPropertyFbio(ni, name, v) ||
      handleCmtPropertyAlag(ni, name, v) ||
      handleCmtPropertyDur(ni, name, v) ||
      handleCmtPropertyRate(ni, name, v);
    (void) tmp;
    return 1;
  }
  return 0;
}

static inline int handleDdtAssign(nodeInfo ni, char *name, int i, D_ParseNode *pn, D_ParseNode *xpn) {
  if (nodeHas(derivative) && i==2) {
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    if (new_de(v)) {
      add_de(ni, name, v, 0, fromDDT);
    }
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
    /* Free(v); */
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
    return 1;
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
    /* Free(v); */
    return 1;
  }
  return 0;
}

static inline int handleRemainingAssignmentsIniProp(nodeInfo ni, char *name, int i, D_ParseNode *pn, D_ParseNode *xpn, char *v) {
  if ((rx_syntax_allow_ini && nodeHas(ini)) || nodeHas(ini0)) {
    sb.o =0; sbDt.o =0;
    /* aAppendN("(__0__)", 7); */
    aType(TINI);
    doDot2(&sb, &sbDt, v);
    if (nodeHas(ini) && !new_de(v)){
      if (tb.idu[tb.id] == 0){
	new_or_ith(v);
	if (tb.lh[tb.ix] == isSuppressedLHS || tb.lh[tb.ix] == 29){
	  tb.lh[tb.ix] = 29;
	} else {
	  tb.lh[tb.ix] = isLhsStateExtra;
	}
      } else {
	updateSyntaxCol();
	sPrint(&_gbuf,"Cannot assign state variable %s; For initial condition assigment use '%s(0) ='",v,v);
	trans_syntax_error_report_fn(_gbuf.s);
      }
    }
    if (!rx_syntax_allow_ini0 && nodeHas(ini0)){
      sPrint(&_gbuf,NOINI0,v);
      updateSyntaxCol();
      trans_syntax_error_report_fn(_gbuf.s);
    }
    return 1;
  }
  return 0;
}

static inline void handleRemainingAssignmentsRestProp(nodeInfo ni, char *name, int i, D_ParseNode *pn, D_ParseNode *xpn, char *v) {
  sb.o = 0; sbDt.o = 0;
  doDot2(&sb, &sbDt, v);
  if (!new_de(v)){
    if (tb.idu[tb.id] == 0){
      // Change to 19 for LHS w/stateExtra
      new_or_ith(v);
      if (tb.lh[tb.ix] == isSuppressedLHS || tb.lh[tb.ix] == 29){
	tb.lh[tb.ix] = 29;
      } else {
	tb.lh[tb.ix] = isLhsStateExtra;
      }
    } else {
      sPrint(&_gbuf,"Cannot assign state variable %s; For initial condition assigment use '%s(0) ='",v,v);
      updateSyntaxCol();
      trans_syntax_error_report_fn(_gbuf.s);
    }
  }
  aType(TASSIGN);
}

static inline int handleRemainingAssignmentsCalcPropMtime(nodeInfo ni, char *name){
  if (nodeHas(mtime)){
    tb.lh[tb.ix] = isLHS;
    tb.mtime[tb.ix] = 1;
    return 1;
  }
  return 0;
}
static inline int handleRemainingAssignmentsCalcPropComplexAssign(nodeInfo ni, char *name, char *v) {
  if (nodeHas(assignment)  || (!rx_syntax_allow_ini && nodeHas(ini))) {
    if (tb.ix+1 == NV && tb.NEnd != NV){
      // New assignment
      tb.ixL = tb.ix;
      tb.lh[tb.ix] = isLHS;
    } else if (tb.ix < 0){
      sPrint(&_gbuf,"cannot assign protected variable '%s'",v);
      updateSyntaxCol();
      trans_syntax_error_report_fn(_gbuf.s);
    } else {
      if (tb.lh[tb.ix] == notLHS){
	tb.lh[tb.ix] = isLHSparam;
      } else {
	tb.lh[tb.ix] = isLHS;
      }
      tb.ixL=-1;
    }
    return 1;
  }
  return 0;
}
static inline int handleRemainingAssignmentsCalcPropIni(nodeInfo ni, char *name, D_ParseNode *pn, char *v) {
  if (nodeHas(ini) || nodeHas(ini0)) {
    D_ParseNode *xpn;
    double d;
    if (tb.ini[tb.ix] == 0){
      // If there is only one initialzation call, then assume
      // this is a parameter with an initial value.
      tb.ini[tb.ix] = 1;
      if (nodeHas(ini0)){
	tb.ini0[tb.ix] = 1;
	xpn = d_get_child(pn, 3);
	/* Free(v); */
	v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	sscanf(v, "%lf", &d);
	tb.iniv[tb.ix] = d;
	tb.ini_i++;
      } else {
	tb.ini0[tb.ix] = 0;
	if (strncmp(v,"rx_",3)==0){
	  tb.lh[tb.ix] = isLHS;
	} else {
	  xpn = d_get_child(pn, 2);
	  /* Free(v); */
	  v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	  sscanf(v, "%lf", &d);
	  tb.iniv[tb.ix] = d;
	  tb.ini_i++;
	}
      }
      /* Free(v); */
      return 1;
    } else {
      // There is more than one call to this variable, it is a
      // conditional variable
      /* Rprintf("Duplicate %s; %d %d\n", v, tb.lh[tb.ix], tb.ini0[tb.ix]); */
      if (tb.lh[tb.ix] != isLHS){
	tb.lh[tb.ix] = isLHS;
	if (nodeHas(ini0) && tb.ini0[tb.ix] == 1){
	  sPrint(&_gbuf,"cannot have conditional initial conditions for '%s'",v);
	  updateSyntaxCol();
	  trans_syntax_error_report_fn(_gbuf.s);
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
  return 0;
}
static inline int handleRemainingAssignmentsCalcProps(nodeInfo ni, char *name, int i, D_ParseNode *pn, D_ParseNode *xpn, char *v) {
  if (!handleRemainingAssignmentsIniProp(ni, name, i, pn, xpn, v)){
    handleRemainingAssignmentsRestProp(ni, name, i, pn, xpn, v);
  }
  new_or_ith(v);
  aProp(tb.ix);
  if (!(handleRemainingAssignmentsCalcPropMtime(ni, name) ||
	handleRemainingAssignmentsCalcPropComplexAssign(ni, name, v))) {
    return handleRemainingAssignmentsCalcPropIni(ni, name, pn, v);
  }
  return 0;
}


static inline int handleRemainingAssignments(nodeInfo ni, char *name, int i, D_ParseNode *pn, D_ParseNode *xpn) {
  if (nodeHas(ini0f) && rx_syntax_allow_ini && i == 0){
    foundF0=1;
    aType(TF0);
    sb.o =0; sbDt.o=0; sbt.o = 0;
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    sAppend(&sb,  "%s",v);
    sAppend(&sbDt,"%s",v);
    sAppend(&sbt, "%s(0)",v);
  }

  if ((i==0 && (nodeHas(assignment) || nodeHas(ini) || nodeHas(ini0))) ||
      (i == 2 && nodeHas(mtime))){
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    int ret1 = handleRemainingAssignmentsCalcProps(ni, name, i, pn, xpn, v);
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
    if (ret1) return 1;
  }
  return 0;
}

static inline int handleDdtRhs(nodeInfo ni, char *name, D_ParseNode *xpn) {
  if (nodeHas(der_rhs)) {
    switch(sbPm.lType[sbPm.n]){
    case TMTIME:
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("modeling times cannot depend on state values"));
      break;
    case FBIO:
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("bioavailability cannot depend on state values"));
      break;
    case ALAG:
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("absorption lag-time cannot depend on state values"));
      break;
    case RATE:
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("model-based rate cannot depend on state values"));
      break;
    case DUR:
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("model-based duration cannot depend on state values"));
      break;
    case TMAT0:
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("model-based matricies cannot depend on state values"));
    default:
      {
	updateSyntaxCol();
	char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
	if (new_de(v)){
	  /* sPrint(&buf2,"d/dt(%s)",v); */
	  updateSyntaxCol();
	  sPrint(&_gbuf,"Tried to use d/dt(%s) before it was defined",v);
	  updateSyntaxCol();
	  trans_syntax_error_report_fn(_gbuf.s);
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
      }
    }
    return 1;
  }
  return 0;
}

static inline int isLineAssignmentStatement(nodeInfo ni, char *name) {
  return nodeHas(assignment) || nodeHas(ini) || nodeHas(dfdy) ||
      nodeHas(ini0) || nodeHas(ini0f) || nodeHas(fbio) || nodeHas(alag) || nodeHas(rate) ||
    nodeHas(dur) || nodeHas(mtime);
}

static inline char * getLineAfterAssign(char *c) {
  while ((*c != '=') && (*c != '~')) {
    c++;
  }
  while ((*c == '=') || (*c == '~') || (*c == ' ')){
    c++;
  }
  return c;
}

static inline int isLineAssigmentProperty(nodeInfo ni, char *name, int *isDepot) {
  return (nodeHas(rate) || nodeHas(alag) || nodeHas(fbio) || nodeHas(dur)) &&
    ((*isDepot = (tb.depotN == tb.di[tb.curPropN])) ||
     (tb.centralN == tb.di[tb.curPropN]));
}

static inline int finalizeLineAssign(nodeInfo ni, char *name, D_ParseNode *pn) {
  if (isLineAssignmentStatement(ni, name)) {
    int isDepot;
    if (isLineAssigmentProperty(ni, name, &isDepot)) {
      char *c = getLineAfterAssign(sb.s);
      if (isDepot){
	curLineType(&depotLines, sbPm.lType[sbPm.n]);
	addLine(&depotLines, "%s", c);
      } else {
	curLineType(&centralLines, sbPm.lType[sbPm.n]);
	addLine(&centralLines, "%s", c);
      }
      /* RSprintf("c: %s, lType: %d\n", c, sbPm.lType[sbPm.n], isDepot); */
    } else {
      addLine(&sbNrmL, "%s;\n", sbt.s);
    }
    addLine(&sbPm,     "%s;\n", sb.s);
    addLine(&sbPmDt,   "%s;\n", sbDt.s);
    sAppend(&sbNrm, "%s;\n", sbt.s);
    ENDLINE;
    return 1;
  }
  return 0;
}

static inline int finalizeLineMat(nodeInfo ni, char *name) {
  if ((nodeHas(mat0) || nodeHas(matF))){
    addLine(&sbPm,     "%s;\n", sb.s);
    addLine(&sbPmDt,   "%s;\n", sbDt.s);
    ENDLINE;
    return 1;
  }
  return 0;
}

static inline int finalizeLineDdt(nodeInfo ni, char *name) {
  if (nodeHas(derivative)){
    addLine(&sbPm,     "%s);\n", sb.s);
    addLine(&sbPmDt,   "%s);\n", sbDt.s);
    sAppend(&sbNrm, "%s;\n", sbt.s);
    addLine(&sbNrmL, "%s;\n", sbt.s);
    ENDLINE;
    return 1;
  }
  return 0;
}

static inline int finalizeLineParam(nodeInfo ni, char *name) {
  if (nodeHas(param_statement)) {
    sbDt.o = 0; sbt.o = 0;
    sAppend(&sbNrm, "%s;\n", sbt.s);
    addLine(&sbNrmL, "%s;\n", sbt.s);
    ENDLINE;
    return 1;
  }
  return 0;
}

static inline int finalizeLineSelectionStatement(nodeInfo ni, char *name, int isWhile) {
  if (nodeHas(selection_statement)){
    sb.o = 0; sbDt.o = 0; sbt.o = 0;
    aType(TLOGIC);
    /* aType(300); */
    if (isWhile) {
      sAppendN(&sb,   "if (_itwhile > _solveData->maxwhile) {_solveData->whileexit=1;break;}\n}\n", 72);
      sAppendN(&sbDt, "if (_itwhile > _solveData->maxwhile) {_solveData->whileexit=1;break;}\n}\n", 72);
      sAppendN(&sbt, "}", 1);
    } else {
      sAppendN(&sb, "}", 1);
      sAppendN(&sbDt, "}", 1);
      sAppendN(&sbt, "}", 1);
    }
    addLine(&sbPm,   "%s\n", sb.s);
    addLine(&sbPmDt, "%s\n", sbDt.s);
    sAppend(&sbNrm,  "%s\n", sbt.s);
    addLine(&sbNrmL, "%s\n", sbt.s);
    ENDLINE;
    return 1;
  }
  return 0;
}

static inline void assertLineEquals(nodeInfo ni, char *name, D_ParseNode *pn){
  if (!rx_syntax_assign && (nodeHas(assignment) || nodeHas(ini) || nodeHas(ini0) || nodeHas(ini0f) || nodeHas(mtime))){
    int i;
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
    /* Free(v); */
  }
}

static inline int finalizeLinePower(nodeInfo ni, char *name) {
  if (nodeHas(power_expression)) {
    aAppendN(")", 1);
    return 1;
  }
  return 0;
}

static inline void finalizeLine(nodeInfo ni, char *name, D_ParseNode *pn, int isWhile, int i) {
  if (isWhile) {
    tb.nwhile--;
  }
  int tmp = finalizeLineAssign(ni, name, pn) ||
    finalizeLineMat(ni, name) ||
    finalizeLineDdt(ni, name) ||
    finalizeLineParam(ni, name) ||
    finalizeLineSelectionStatement(ni, name, isWhile) ||
    finalizeLinePower(ni, name);
  (void) tmp;
  assertLineEquals(ni, name, pn);
}

////////////////////////////////////////////////////////////////////////////////
// assertions
static inline int assertNoRAssign(nodeInfo ni, char *name, D_ParseNode *pn, int i){
  if (!rx_syntax_assign  &&
      ((i == 4 && nodeHas(derivative)) ||
       (i == 6 && nodeHas(dfdy)))) {
    D_ParseNode *xpn = d_get_child(pn,i);
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    if (!strcmp("<-",v)){
      updateSyntaxCol();
      trans_syntax_error_report_fn(NOASSIGN);
    }
    /* Free(v); */
    /* continue; */
    return 1;
  }
  return 0;
}

static inline void assertEndSemicolon(nodeInfo ni, char *name, int i, D_ParseNode *xpn) {
  if (rx_syntax_require_semicolon && nodeHas(end_statement) && i == 0){
    if (xpn->start_loc.s ==  xpn->end){
      updateSyntaxCol();
      trans_syntax_error_report_fn(NEEDSEMI);
    }
  }
}

static inline int parseNodePossiblySkipRecursion(nodeInfo ni, char *name, D_ParseNode *pn, D_ParseNode *xpn,
						 int *i, int nch, int *depth) {
  if (assertNoRAssign(ni, name, pn, *i) ||
      isSkipChild(ni, name, *i))  return 1;

  // Inductive linearization matrices
  handleIndLinMat0(ni, name);
  handleIndLinMatf(ni, name);
  // Determine if this is a function and change depth flag if needed
  setFunctionFlag(ni, name, *i, depth);
  if (handleIfElse(ni, name, *i) ||
      // simeta()/simeps()
      handleSimFunctions(ni, name, i, nch, pn) ||
      handleStringEqualityStatements(ni, name, *i, xpn) ||
      handleDvidStatement(ni, name, xpn, pn) ||
      handleFunctions(ni, name, i, depth, nch, xpn, pn) ||
      handleTheta(ni, name, xpn) ||
      handleEta(ni, name, xpn)) return 1;

  if (nodeHas(param_statement) && i == 0) {
    sAppendN(&sbt,"param", 5);
    sbDt.o = 0;
  }
  return 0;
}

static inline int parseNodeAfterRecursion(nodeInfo ni, char *name, D_ParseNode *pn, D_ParseNode *xpn,
					  int *i, int nch, int *depth, int *safe_zero,
					  int *ii, int *found, int *isWhile) {
  assertEndSemicolon(ni, name, *i, xpn);
  handleSafeZero(ni, name, *i, safe_zero, xpn);  // protect against divide by zeros
  if (handlePrintf(ni, name, *i, xpn) ||
      handleJac(ni, name, *i, xpn, ii, found) ||
      handleLogicalExpr(ni, name, *i, pn, xpn, isWhile) ||
      handleCmtProperty(ni, name, *i, xpn) ||
      handleDdtAssign(ni, name, *i, pn, xpn) ||
      handleDdtRhs(ni, name, xpn)) return 1;
  if (*i==0 && nodeHas(power_expression)) {
    aAppendN("),", 2);
    sAppendN(&sbt, "^", 1);
  }
  if (!rx_syntax_star_pow && *i == 1 && nodeHas(power_expression)){
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    if (!strcmp("**",v)){
      updateSyntaxCol();
      trans_syntax_error_report_fn(NEEDPOW);
    }
  }
  handleRemainingAssignments(ni, name, *i, pn, xpn);
  return 0;
}

void wprint_parsetree(D_ParserTables pt, D_ParseNode *pn, int depth, print_node_fn_t fn, void *client_data) {
  char *name = (char*)pt.symbols[pn->symbol].name;
  nodeInfo ni;
  niReset(&ni);
  int nch = d_get_number_of_children(pn), i, ii, found, safe_zero = 0;
  char *value = (char*)rc_dup_str(pn->start_loc.s, pn->end);
  // Add symbol, check/flag if recursive
  handleIdentifier(ni, name, value);
  // Add (double) in front of function arguments
  handleFunctionArguments(name, depth);
  // print/change identifier/operator and change operator information (if needed)
  handleOperatorsOrPrintingIdentifiers(depth, fn, client_data, ni, name, value);

  if (nch != 0) {
    int isWhile=0;
    if (nodeHas(power_expression)) {
      aAppendN("Rx_pow(_as_dbleps(", 18);
    }
    for (i = 0; i < nch; i++) {
      D_ParseNode *xpn = d_get_child(pn, i);

      if (parseNodePossiblySkipRecursion(ni, name, pn, xpn, &i, nch, &depth)) continue;

      // Recursively parse tree
      wprint_parsetree(pt, xpn, depth, fn, client_data);

      parseNodeAfterRecursion(ni, name, pn, xpn, &i, nch, &depth, &safe_zero,
			      &ii, &found, &isWhile);
    }
    finalizeLine(ni, name, pn, isWhile, i);
  }
}

void err_msg(int chk, const char *msg, int code)
{
  if(!chk) {
    Rf_errorcall(R_NilValue, "%s",msg);
  }
}
// show_ode = 1 dydt
// show_ode = 2 Jacobian
// show_ode = 3 Ini statement
// show_ode = 0 LHS
// show_ode = 5 functional bioavailibility
// show_ode == 6 functional lag
// show_ode == 7 functional rate
// show_ode == 8 functional duration
// show_ode == 9 functional mtimes
// show_ode == 10 ME matrix
// show_ode == 11 Inductive vector
// show_ode == 12 initialize lhs to last value
// show_ode == 13 #define lags for lhs values
// show_ode == 14 #define lags for params/covs
// show_ode == 15 #define sync lhs for simeps
// show_ode == 16 #define sync lhs for simeps

#define ode_lhs 0
#define ode_dydt 1
#define ode_jac  2
#define ode_ini 3
#define ode_simeps 15
#define ode_simeta 16

// Scenarios
#define print_double 0
#define print_populateParameters 1
#define print_void 2
#define print_lastLhsValue  3
#define print_lhsLags 4
#define print_paramLags 5
#define print_simeps 15
#define print_simeta 16


static inline void printDdtDefine(int show_ode, int scenario) {
  if (show_ode == ode_jac || show_ode == ode_lhs){
    //__DDtStateVar_#__
    // These will be defined and used in Jacobian or LHS functions
    for (int i = 0; i < tb.de.n; i++){
      if (scenario == print_double){
	sAppend(&sbOut,"  double  __DDtStateVar_%d__;\n",i);
      } else {
	sAppend(&sbOut,"  (void)__DDtStateVar_%d__;\n",i);
      }
    }
  }
}

static inline void printPDStateVar(int show_ode, int scenario) {
  // Now get Jacobain information  __PDStateVar_df_dy__ if needed
  char *buf1, *buf2;
  if (show_ode != ode_ini && show_ode != ode_simeps){
    for (int i = 0; i < tb.ndfdy; i++){
      buf1 = tb.ss.line[tb.df[i]];
      buf2 = tb.ss.line[tb.dy[i]];
      // This is for dydt/ LHS/ or jacobian for df(state)/dy(parameter)
      if (show_ode == ode_dydt || show_ode == ode_lhs || tb.sdfdy[i] == 1){
	if (scenario == print_double){
	  sAppend(&sbOut,"  double __PDStateVar_%s_SeP_%s__;\n",buf1,buf2);
	} else {
	  sAppend(&sbOut,"  (void)__PDStateVar_%s_SeP_%s__;\n",buf1,buf2);
	}
      }
    }
  }
}

static inline int isStateLhsI(int i) {
  if (tb.lh[i] == isState){
    int doCont=0;
    for (int j = 0; j < tb.de.n; j++) {
      if (tb.di[j] == i) {
	if (!tb.idu[j]) doCont = 1;
	break;
      }
    }
    if (doCont) return 1;
  }
  return 0;
}

static inline int shouldSkipPrintLhsI(int scenario, int lhs, int i) {
  if (isStateLhsI(i)) return 1;
  switch(scenario){
  case print_paramLags:
    return (tb.lag[i] == notLHS || tb.lh[i] == isState);
  case print_lhsLags:
    return (tb.lag[i] == 0 ||tb.lh[i] != isLHS);
  case print_lastLhsValue:
    return !(tb.lh[i] == isLHS || tb.lh[i] == isLhsStateExtra || tb.lh[i] == isLHSparam);
  }
  return (lhs && tb.lh[i]>0 && tb.lh[i] != isLHSparam);
}

static inline void printParamLags(char *buf, int *j) {
  sAppendN(&sbOut, "#define diff_", 13);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "1(x) (x - _getParCov(_cSub, _solveData, %d, (&_solveData->subjects[_cSub])->idx - 1))\n", *j);
  sAppendN(&sbOut, "#define diff_", 13);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "(x,y) (x - _getParCov(_cSub, _solveData, %d, (&_solveData->subjects[_cSub])->idx - (y)))\n", *j);
  sAppendN(&sbOut, "#define first_", 14);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "1(x) _getParCov(_cSub, _solveData, %d, NA_INTEGER)\n", *j);
  sAppendN(&sbOut, "#define last_", 13);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "1(x) _getParCov(_cSub, _solveData, %d, (&_solveData->subjects[_cSub])->n_all_times - 1)\n", *j);
  sAppendN(&sbOut, "#define lead_", 13);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "1(x) _getParCov(_cSub, _solveData, %d, (&_solveData->subjects[_cSub])->idx + 1)\n", *j);
  sAppendN(&sbOut, "#define lead_", 13);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "(x, y) _getParCov(_cSub, _solveData, %d, (&_solveData->subjects[_cSub])->idx + (y))\n", *j);
  sAppendN(&sbOut, "#define lag_", 12);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "1(x) _getParCov(_cSub, _solveData, %d, (&_solveData->subjects[_cSub])->idx - 1)\n", *j);
  sAppendN(&sbOut, "#define lag_", 12);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "(x,y) _getParCov(_cSub, _solveData, %d, (&_solveData->subjects[_cSub])->idx - (y))\n", *j);
  j[0]=j[0]+1;
}

static inline void printLhsLag(char *buf, int *j) {
  sAppendN(&sbOut, "#define lead_", 13);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "1(x) _solveData->subjects[_cSub].lhs[%d]\n", *j);
  sAppendN(&sbOut, "#define lead_", 13);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "(x,y) _solveData->subjects[_cSub].lhs[%d]\n", *j);
  sAppendN(&sbOut, "#define diff_", 13);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "1(x) _solveData->subjects[_cSub].lhs[%d]\n", *j);
  sAppendN(&sbOut, "#define diff_", 13);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "(x,y) _solveData->subjects[_cSub].lhs[%d]\n", *j);
  sAppendN(&sbOut, "#define lag_", 12);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "1(x) _solveData->subjects[_cSub].lhs[%d]\n", *j);
  sAppendN(&sbOut, "#define lag_", 12);
  doDot(&sbOut, buf);
  sAppend(&sbOut, "(x, y) _solveData->subjects[_cSub].lhs[%d]\n", *j);
  j[0] = j[0]+1;
}

static inline void printLastLhsValue(char *buf, int *j) {
  sAppendN(&sbOut, "  ", 2);
  doDot(&sbOut, buf);
  sAppend(&sbOut, " = _PL[%d];\n", *j);
  j[0] = j[0]+1;
}

static inline void printDoubleDeclaration(char *buf) {
  sAppendN(&sbOut,"  double ", 9);
  doDot(&sbOut, buf);
  if (!strcmp("rx_lambda_", buf) || !strcmp("rx_yj_", buf) ||
      !strcmp("rx_hi_", buf) || !strcmp("rx_low_", buf)){
    sAppendN(&sbOut, "__", 2);
  }
  sAppendN(&sbOut, ";\n", 2);
}

static inline void printVoidDeclaration(char *buf) {
  sAppend(&sbOut,"  ");
  sAppend(&sbOut,"(void)");
  doDot(&sbOut, buf);
  if (!strcmp("rx_lambda_", buf) || !strcmp("rx_yj_", buf) ||
      !strcmp("rx_low_", buf) || !strcmp("rx_hi_", buf)){
    sAppendN(&sbOut, "__", 2);
  }
  sAppendN(&sbOut, ";\n", 2);
}

static inline void printPopulateParameters(char *buf, int *j) {
  sAppendN(&sbOut,"  ", 2);
  doDot(&sbOut, buf);
  sAppend(&sbOut, " = _PP[%d];\n", *j);
  j[0] = j[0]+1;
}

static inline void printSimEps(char *buf, int *j) {
  sAppend(&sbOut,"  if (_solveData->svar[_svari] == %d) {", *j);
  doDot(&sbOut, buf);
  sAppend(&sbOut, " = _PP[%d];}; ", *j);
  j[0] = j[0]+1;
}

static inline void printSimEta(char *buf, int *j) {
  sAppend(&sbOut,"  if (_solveData->ovar[_ovari] == %d) {", *j);
  doDot(&sbOut, buf);
  sAppend(&sbOut, " = _PP[%d];}; ", *j);
  j[0] = j[0]+1;
}


/* when prnt_vars() is called, user defines the behavior in "case" */
void prnt_vars(int scenario, int lhs, const char *pre_str, const char *post_str, int show_ode) {
  int i, j;
  char *buf;
  sAppend(&sbOut, "%s", pre_str);
  if (scenario == print_double || scenario == print_void){
    printDdtDefine(show_ode, scenario);
    printPDStateVar(show_ode, scenario);
  }
  for (i=0, j=0; i<NV; i++) {
    if (shouldSkipPrintLhsI(scenario, lhs, i)) continue;
    buf = tb.ss.line[i];
    switch(scenario) {
    case print_paramLags: // Case 5 is for using #define lag_var(x)
      printParamLags(buf, &j);
      break;
    case print_lhsLags: // Case 4 is for using #define lag_var(x)
      printLhsLag(buf, &j);
      break;
    case print_lastLhsValue: // Case 3 is for using the last lhs value
      printLastLhsValue(buf, &j);
      break;
    case print_double:   // Case 0 is for declaring the variables
      printDoubleDeclaration(buf);
      break;
    case print_void: // Case 2 is for suppressing all the warnings for the variables by using (void)var;
      // See https://stackoverflow.com/questions/1486904/how-do-i-best-silence-a-warning-about-unused-variables
      printVoidDeclaration(buf);
      break;
    case print_populateParameters:
      // Case 1 is for declaring the par_ptr.
      printPopulateParameters(buf, &j);
      break;
    case print_simeps:
      // Case 15 is for declaring eps the sync parameters
      printSimEps(buf, &j);
      break;
    case print_simeta:
      // Case 16 is for declaring eta the sync parameters
      printSimEta(buf, &j);
      break;
    default: break;
    }
  }
  sAppend(&sbOut, "%s", post_str);
}

static inline void printCModelVars(const char *prefix) {
  sAppend(&sbOut, "extern SEXP %smodel_vars(){\n  int pro=0;\n", prefix);
  sAppend(&sbOut, "  SEXP _mv = PROTECT(_rxGetModelLib(\"%smodel_vars\"));pro++;\n", prefix);
  sAppendN(&sbOut, "  if (!_rxIsCurrentC(_mv)){\n", 28);
  sAppendN(&sbOut, "    SEXP hash    = PROTECT(allocVector(STRSXP, 1));pro++;\n", 58);
  sAppend(&sbOut, "#define __doBuf__  sprintf(buf, \"", _mv.o+1);
  int off=0;
  int off2 = 0;
  for (int i = 0; i < _mv.o; i++){
    if (off != 0 && off % 4095 == 0) {
      sAppend(&sbOut, "\"); \\\n sprintf(buf+%d, \"", off2);
    }
    off++;
    off2++;
    if (_mv.s[i] == '%'){
      sAppendN(&sbOut, "%%", 2);
      off++;
    } else if (_mv.s[i] == '?') {
      // Avoid digrahps/trigraphs
      sAppendN(&sbOut, "\\?", 2);
    } else if (_mv.s[i] == '"'){
      sAppendN(&sbOut, "\\\"", 2);
    } else if (_mv.s[i] == '\''){
      sAppendN(&sbOut, "'", 1);
    } else if (_mv.s[i] == ' '){
      sAppendN(&sbOut, " ", 1);
    } else if (_mv.s[i] == '\n'){
      sAppendN(&sbOut, "\\n", 2);
    } else if (_mv.s[i] == '\t'){
      sAppendN(&sbOut, "\\t", 2);
    } else if (_mv.s[i] == '\\'){
      sAppendN(&sbOut, "\\\\", 2);
    } else if (_mv.s[i] >= 33  && _mv.s[i] <= 126){ // ASCII only
      sPut(&sbOut, _mv.s[i]);
    }
  }
  sAppendN(&sbOut, "\");\n", 4);
  sAppend(&sbOut,"    char buf[%d];\n    __doBuf__\n#undef __doBuf__\n", off+1);
  sAppendN(&sbOut,"    SET_STRING_ELT(hash, 0, mkChar(buf));\n", 42);
  sAppendN(&sbOut, "    SEXP lst      = PROTECT(_rxQr(hash));pro++;\n", 48);
  sAppendN(&sbOut, "    _assign_ptr(lst);\n", 22);
  sAppendN(&sbOut, "    UNPROTECT(pro);\n", 20);

  sAppendN(&sbOut, "    return lst;\n", 16);
  sAppendN(&sbOut, "  } else {\n", 11);
  sAppendN(&sbOut, "    UNPROTECT(pro);\n", 20);
  sAppendN(&sbOut, "    return _mv;\n", 16);
  sAppendN(&sbOut, "  }\n", 4);
  sAppendN(&sbOut, "}\n", 2);
}

static inline void printRInit(const char *libname, const char *libname2, const char *prefix) {
  sAppend(&sbOut,"\n//Create function to call from R's main thread that assigns the required functions. Sometimes they don't get assigned.\nextern void %sassignFuns(){\n  _assignFuns();\n}\n", prefix);
  sAppend(&sbOut,"\n//Initialize the dll to match RxODE's calls\nvoid R_init0_%s(){\n  // Get C callables on load; Otherwise it isn't thread safe\n", libname2);
  sAppendN(&sbOut, "  _assignFuns();\n", 17);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sassignFuns\", (DL_FUNC) %sassignFuns);\n", libname, prefix, prefix);
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
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sME\", (DL_FUNC) %sME);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sIndF\", (DL_FUNC) %sIndF);\n", libname, prefix, prefix);
  sAppend(&sbOut, "  R_RegisterCCallable(\"%s\",\"%sdydt_liblsoda\", (DL_FUNC) %sdydt_liblsoda);\n", libname, prefix, prefix);
  sAppend(&sbOut,"}\n//Initialize the dll to match RxODE's calls\nvoid R_init_%s(DllInfo *info){\n  // Get C callables on load; Otherwise it isn't thread safe\n  R_init0_%s();", libname2, libname2);
  sAppend(&sbOut, "\n  static const R_CallMethodDef callMethods[]  = {\n    {\"%smodel_vars\", (DL_FUNC) &%smodel_vars, 0},\n    {NULL, NULL, 0}\n  };\n",
  	  prefix, prefix);
  sAppendN(&sbOut, "\n  R_registerRoutines(info, NULL, callMethods, NULL, NULL);\n  R_useDynamicSymbols(info,FALSE);\n}\n", 97);
  sAppend(&sbOut, "\nvoid R_unload_%s (DllInfo *info){\n  // Free resources required for single subject solve.\n  SEXP _mv = PROTECT(_rxGetModelLib(\"%smodel_vars\"));\n",
	  libname2, prefix);
  sAppend(&sbOut, "  if (!isNull(_mv)){\n    _rxRmModelLib(\"%smodel_vars\");\n  }\n  UNPROTECT(1);\n}\n", prefix);
}

void print_aux_info(char *model, const char *prefix, const char *libname, const char *pMd5, const char *timeId,
		    const char *libname2){
  sbuf bufw;
  sIniTo(&bufw, 1024);
  sClear(&s_aux_info);
  /* char bufw[1024]; */
  printCModelVars(prefix);

  sAppend(&sbOut,"extern void %sdydt_lsoda(int *neq, double *t, double *A, double *DADT)\n{\n  %sdydt(neq, *t, A, DADT);\n}\n", prefix, prefix);
  sAppend(&sbOut, "extern int %sdydt_liblsoda(double __t, double *y, double *ydot, void *data)\n{\n  int *neq = (int*)(data);\n  %sdydt(neq, __t, y, ydot);\n  return(0);\n}\n",
	  prefix,prefix);
  sAppend(&sbOut,"extern void %scalc_jac_lsoda(int *neq, double *t, double *A,int *ml, int *mu, double *JAC, int *nrowpd){\n  // Update all covariate parameters\n  %scalc_jac(neq, *t, A, JAC, *nrowpd);\n}\n",
	  prefix, prefix);

  printRInit(libname, libname2, prefix);

  sFree(&bufw);
}


void codegen(char *model, int show_ode, const char *prefix, const char *libname, const char *pMd5, const char *timeId, const char *libname2) {
  if (show_ode == 4) {
    print_aux_info(model, prefix, libname, pMd5, timeId, libname2);
  } else {
    int i, j;
    char *buf;
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
      // Now define lhs lags
      prnt_vars(print_lhsLags, 1, "", "", 13);
      // And covariate/parameter lags
      prnt_vars(print_paramLags, 1, "", "", 15);
      // Add sync PP define
      prnt_vars(print_simeps, 1, "#define _SYNC_simeps_ for (int _svari=_solveData->neps; _svari--;){", "}\n", 15);
      prnt_vars(print_simeta, 1, "#define _SYNC_simeta_ for (int _ovari=_solveData->neta; _ovari--;){", "}\n", 16);
      sAppendN(&sbOut,"#include \"extraC.h\"\n", 20);
      sAppendN(&sbOut,"#include <RxODE_model_shared.c>\n", 32);
      sAppend(&sbOut, "extern void  %sode_solver_solvedata (rx_solve *solve){\n  _solveData = solve;\n}\n",prefix);
      sAppend(&sbOut, "extern rx_solve *%sode_solver_get_solvedata(){\n  return _solveData;\n}\n", prefix);
      sAppend(&sbOut, "SEXP %smodel_vars();\n", prefix);
      sAppendN(&sbOut,"\n", 1);
      sAppendN(&sbOut, "\n// prj-specific differential eqns\nvoid ", 40);
      sAppend(&sbOut, "%sdydt(int *_neq, double __t, double *__zzStateVar__, double *__DDtStateVar__)\n{\n  int _itwhile = 0;\n  (void)_itwhile;\n  int _cSub = _neq[1];\n  double t = __t + _solveData->subjects[_neq[1]].curShift;\n  (void)t;\n  ", prefix);
    } else if (show_ode == 2){
      sAppend(&sbOut, "// Jacobian derived vars\nvoid %scalc_jac(int *_neq, double __t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {\n  int _itwhile = 0;\n  (void)_itwhile;\n    int _cSub=_neq[1];\n  double t = __t + _solveData->subjects[_neq[1]].curShift;\n  (void)t;\n  ", prefix);
    } else if (show_ode == 3){
      sAppend(&sbOut,  "// Functional based initial conditions.\nvoid %sinis(int _cSub, double *__zzStateVar__){\n  int _itwhile = 0;\n  (void)_itwhile;\n  \n", prefix);
      if (foundF0){
	sAppendN(&sbOut, "  double t=0;\n", 14);
      }
    } else if (show_ode == 5){
      if (foundF){
	int nnn = tb.de.n;
	if (tb.linCmt){
	  if (tb.hasKa){
	    nnn+=2;
	  } else {
	    nnn+=1;
	  }
	}
	sAppend(&sbOut,  "// Functional based bioavailability (returns amount)\ndouble %sF(int _cSub,  int _cmt, double _amt, double __t, double *__zzStateVar__){\n  int _itwhile = 0;\n  (void)_itwhile;\n  double *_f=_solveData->subjects[_cSub].cF;\n  (void)_f;\n  double t = __t + _solveData->subjects[_cSub].curShift;\n  (void)t;\n  ",
		prefix, nnn);
	for (int jjj = nnn; jjj--;){
	  sAppend(&sbOut, "  _f[%d]=1.0;\n",jjj);
	}
      } else {
	sAppend(&sbOut,  "// Functional based bioavailability\ndouble %sF(int _cSub,  int _cmt, double _amt, double __t, double *__zzStateVar__){\n return _amt;\n",
		prefix);
      }
    } else if (show_ode == 6){
      if (foundLag){
	int nnn = tb.de.n;
	if (tb.linCmt){
	  if (tb.hasKa){
	    nnn+=2;
	  } else {
	    nnn+=1;
	  }
	}
	sAppend(&sbOut,  "// Functional based absorption lag\ndouble %sLag(int _cSub,  int _cmt, double __t){\n  int _itwhile = 0;\n  (void)_itwhile;\n  double *restrict _alag = _solveData->subjects[_cSub].alag;\n  (void)_alag; \n  double t = __t + _solveData->subjects[_cSub].curShift;\n  (void)t;\n  ",
		prefix, nnn);
	for (int jjj = nnn; jjj--;){
	  sAppend(&sbOut, "  _alag[%d]=0.0;\n",jjj);
	}
      } else {
	sAppend(&sbOut,  "// Functional based absorption lag\ndouble %sLag(int _cSub,  int _cmt, double __t, double *__zzStateVar__){\n return __t;\n",
		prefix);
      }
    } else if (show_ode == 7){
      if (foundRate){
	int nnn = tb.de.n;
	if (tb.linCmt){
	  if (tb.hasKa){
	    nnn+=2;
	  } else {
	    nnn+=1;
	  }
	}
	sAppend(&sbOut,  "// Modeled zero-order rate\ndouble %sRate(int _cSub,  int _cmt, double _amt, double __t){\n    int _itwhile = 0;\n  (void)_itwhile;\n  double *restrict _rate= _solveData->subjects[_cSub].cRate;\n  (void)_rate;\n   double t = __t + _solveData->subjects[_cSub].curShift;\n  (void)t;\n  ",
		prefix, nnn);
	for (int jjj = nnn; jjj--;){
	  sAppend(&sbOut, "  _rate[%d]=0.0;\n",jjj);
	}
      } else {
	sAppend(&sbOut,  "// Modeled zero-order rate\ndouble %sRate(int _cSub,  int _cmt, double _amt, double __t, double *__zzStateVar__){\n return 0.0;\n",
		prefix);
      }
    } else if (show_ode == 8){
      if (foundDur){
	int nnn = tb.de.n;
	if (tb.linCmt){
	  if (tb.hasKa){
	    nnn+=2;
	  } else {
	    nnn+=1;
	  }
	}
	sAppend(&sbOut,  "// Modeled zero-order duration\ndouble %sDur(int _cSub,  int _cmt, double _amt, double __t){\n  int _itwhile = 0;\n  (void)_itwhile;\n double *restrict _dur = _solveData->subjects[_cSub].cDur;\n  (void)_dur;\n    double t = __t + _solveData->subjects[_cSub].curShift;\n  (void)t;\n  ",
		prefix, nnn);
	for (int jjj = nnn; jjj--;){
	  sAppend(&sbOut, "  _dur[%d]=0.0;\n",jjj);
	}
      } else {
	sAppend(&sbOut,  "// Modeled zero-order duration\ndouble %sDur(int _cSub,  int _cmt, double _amt, double __t){\n return 0.0;\n",
		prefix);
      }
    } else if (show_ode == 9){
      if (nmtime){
	sAppend(&sbOut,  "// Model Times\nvoid %smtime(int _cSub, double *_mtime){\n  int _itwhile = 0;\n  (void)_itwhile;\n  double t = 0;\n  ",
		prefix);
      } else {
	sAppend(&sbOut,  "// Model Times\nvoid %smtime(int _cSub, double *_mtime){\n",
		prefix);
      }
    } else if (show_ode == 10){
      sAppend(&sbOut, "// Matrix Exponential (%d)\nvoid %sME(int _cSub, double _t, double __t, double *_mat, const double *__zzStateVar__){\n  int _itwhile = 0;\n  (void)_itwhile;\n  double t = __t + _solveData->subjects[_cSub].curShift;\n  (void)t;\n  ",
	      tb.matn, prefix);
    } else if (show_ode == 11){
      sAppend(&sbOut, "// Inductive linearization Matf\nvoid %sIndF(int _cSub, double _t, double __t, double *_matf){\n int _itwhile = 0;\n  (void)_itwhile;\n  double t = __t + _solveData->subjects[_cSub].curShift;\n  (void)t;\n  ", prefix);
    } else {
      sAppend(&sbOut,  "// prj-specific derived vars\nvoid %scalc_lhs(int _cSub, double __t, double *__zzStateVar__, double *_lhs) {\n    int _itwhile = 0;\n  (void)_itwhile;\n  double t = __t + _solveData->subjects[_cSub].curShift;\n  (void)t;\n  ", prefix);
    }
    if ((show_ode == 2 && found_jac == 1 && good_jac == 1) ||
	(show_ode != 2 && show_ode != 3 && show_ode != 5  && show_ode != 8 &&
	 show_ode != 7 && show_ode != 6 &&
	 show_ode !=0 && show_ode != 9 && show_ode != 10 && show_ode != 11) ||
	(show_ode == 8 && foundDur) ||
	(show_ode == 7 && foundRate) ||
	(show_ode == 6 && foundLag) ||
	(show_ode == 5 && foundF) ||
	(show_ode == 3 && foundF0) ||
	(show_ode == 0 && tb.li) ||
	(show_ode == 9 && nmtime) ||
	(show_ode == 10 && tb.matn) ||
	(show_ode == 11 && tb.matnf)){
      prnt_vars(print_double, 0, "", "\n",show_ode);     /* declare all used vars */
      if (maxSumProdN > 0 || SumProdLD > 0){
	int mx = maxSumProdN;
	if (SumProdLD > mx) mx = SumProdLD;
	sAppend(&sbOut,  "  double _p[%d], _input[%d];\n", mx, mx);
	sAppend(&sbOut,  "  double _pld[%d];\n", mx);
	sAppend(&sbOut,  "  for (int ddd=%d; ddd--;){_p[ddd]=_input[ddd]=_pld[ddd]=0.0;}", mx);

      }
      else prnt_vars(print_void, 0, "  (void)t;\n", "\n",show_ode);     /* declare all used vars */
      if (maxSumProdN){
	sAppendN(&sbOut,  "  (void)_p;\n  (void)_input;\n", 28);
	if (SumProdLD){
	  sAppendN(&sbOut,  "  (void)_pld;\n", 14);
	}
      }
      prnt_vars(print_lastLhsValue, 0,"","\n", 12);
      if (show_ode == 3){
	sAppendN(&sbOut, "  _update_par_ptr(0.0, _cSub, _solveData, _idx);\n", 49);
      } else if (show_ode == 6 || show_ode == 7 || show_ode == 8 || show_ode == 9){
	// functional lag, rate, duration, mtime
	sAppendN(&sbOut, "  _update_par_ptr(NA_REAL, _cSub, _solveData, _idx);\n", 53);
      } else if (show_ode == 11 || show_ode == 10){
	sAppendN(&sbOut, "  _update_par_ptr(_t, _cSub, _solveData, _idx);\n", 48);
      } else {
	sAppendN(&sbOut, "  _update_par_ptr(__t, _cSub, _solveData, _idx);\n", 49);
      }
      prnt_vars(print_populateParameters, 1, "", "\n",show_ode);                   /* pass system pars */
      if (show_ode != 9 && show_ode != 11){
	for (i=0; i<tb.de.n; i++) {                   /* name state vars */
	  buf = tb.ss.line[tb.di[i]];
	  if(tb.idu[i] != 0){
	    if (show_ode == 6 || show_ode == 8 || show_ode == 7){
	      sAppendN(&sbOut, "  ", 2);
	      doDot(&sbOut, buf);
	      sAppend(&sbOut, " = NA_REAL;\n", i, i);
	    } else {
	      // stateExtra
	      sAppendN(&sbOut, "  ", 2);
	      doDot(&sbOut, buf);
	      sAppend(&sbOut, " = __zzStateVar__[%d]*((double)(_ON[%d]));\n", i, i);
	    }
	  } else {
	    break;
	  }
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
	case TLIN:
	  if (show_ode != 10 && show_ode != 11 &&
	      show_ode != 5 && show_ode != 6 &&
	      show_ode != 7 && show_ode !=8){
	    sAppend(&sbOut,"  %s",show_ode == 1 ? sbPm.line[i] : sbPmDt.line[i]);
	  }
	  break;
	case TMTIME:
	case TASSIGN:
	  if (show_ode != 10 && show_ode != 11){
	    sAppend(&sbOut,"  %s",show_ode == 1 ? sbPm.line[i] : sbPmDt.line[i]);
	  }
	  break;
	case TINI:
	  // See if this is an ini or a reclaimed expression.
	  if (show_ode != 10 && show_ode != 11){
	    if (sbPm.lProp[i] >= 0 ){
	      tb.ix = sbPm.lProp[i];
	      if (tb.lh[tb.ix] == isLHS || tb.lh[tb.ix] == isLHSparam){
		sAppend(&sbOut,"  %s",show_ode == 1 ? sbPm.line[i] : sbPmDt.line[i]);
	      }
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
	      show_ode != 7 && show_ode != 8 && show_ode != 9 &&
	      show_ode !=10 && show_ode != 11){
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
	  if (show_ode != 10 && show_ode != 11){
	    sAppend(&sbOut,"  %s",show_ode == 1 ? sbPm.line[i] : sbPmDt.line[i]);
	  }
	  break;
	case TMAT0:
	  if (show_ode == 10){
	    sAppend(&sbOut,"  %s", sbPm.line[i]);
	  }
	  break;
	case TMATF:
	  if (show_ode == 11){
	    sAppend(&sbOut,"  %s", sbPm.line[i]);
	  }
	  break;
	default:
	  RSprintf("line Number: %d\n", i);
	  RSprintf("type: %d\n", sbPm.lType[i]);
	  RSprintf("line: %s\n", sbPm.line[i]);
	  RSprintf("PmDt Line: %s\n", sbPmDt.line[i]);
	  RSprintf("Prop: %d\n", sbPm.lProp[i]);
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
	sAppendN(&sbOut, "\n  return t + _alag[_cmt] - _solveData->subjects[_cSub].curShift;\n", 66);
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
	for (i = 0; i < tb.de.n; i++) {
	  if (tb.idu[i]) {
	    buf=tb.ss.line[tb.di[i]];
	    sAppend(&sbOut, "  __zzStateVar__[%d]=((double)(_ON[%d]))*(",i,i);
	    doDot(&sbOut, buf);
	    sAppendN(&sbOut,  ");\n", 3);
	  }
	}
      }
      sAppendN(&sbOut,  "}\n", 2);
    } else if (show_ode == 5 || show_ode == 6 || show_ode == 7 || show_ode == 8){
      sAppendN(&sbOut,  "}\n", 2);
    } else if (show_ode == 0 && tb.li){
      sAppendN(&sbOut,  "\n", 1);
      for (i=0, j=0; i<NV; i++) {
	if (tb.lh[i] != isLHS && tb.lh[i] != isLhsStateExtra && tb.lh[i] != isLHSparam) continue;
	buf = tb.ss.line[i];
	sAppend(&sbOut,  "  _lhs[%d]=", j);
	doDot(&sbOut, buf);
	sAppendN(&sbOut,  ";\n", 2);
	j++;
      }
      sAppendN(&sbOut,  "}\n", 2);
    } else if (show_ode == 9 && nmtime){
      sAppendN(&sbOut,  "\n", 1);
      for (i=0, j=0; i<NV; i++) {
	if (tb.mtime[i] != 1) continue;
	buf = tb.ss.line[i];
	sAppend(&sbOut,  "  _mtime[%d]=", j);
	doDot(&sbOut, buf);
	sAppendN(&sbOut,  ";\n", 2);
	j++;
      }
      sAppendN(&sbOut,  "}\n", 2);
    } else {
      sAppendN(&sbOut,  "}\n", 2);
    }
  }
}

sbuf _bufw, _bufw2;

void parseFree(int last){
  sFree(&sb);
  sFree(&sbDt);
  sFree(&sbt);
  sFree(&sbNrm);
  sFree(&s_aux_info);
  sFree(&s_inits);
  sFree(&_bufw);
  sFree(&_bufw2);
  sFree(&firstErr);
  sFree(&_gbuf);
  sFree(&_mv);
  lineFree(&sbPm);
  lineFree(&sbPmDt);
  lineFree(&sbNrmL);
  lineFree(&(tb.ss));
  lineFree(&(tb.de));
  lineFree(&depotLines);
  lineFree(&centralLines);
  lineFree(&_dupStrs);
  Free(tb.lh);
  Free(tb.lag);
  Free(tb.ini);
  Free(tb.mtime);
  Free(tb.iniv);
  Free(tb.ini0);
  Free(tb.di);
  Free(tb.idi);
  Free(tb.idu);
  Free(tb.dvid);
  Free(tb.df);
  Free(tb.dy);
  Free(tb.sdfdy);
  freeP();
  if (last){
    Free(gBuf);
    sFree(&sbOut);
    freeP();
    sFree(&_bufw);
    sFree(&_bufw2);
  }
}
void reset (){
  // Reset sb/sbt string buffers
  parseFree(0);
  sIniTo(&_bufw, 1024);
  sIniTo(&_bufw2, 2100);
  sIniTo(&sb, MXSYM);
  sIniTo(&sbDt, MXDER);
  sIniTo(&sbt, SBUF_MXBUF);
  sIniTo(&sbNrm, SBUF_MXBUF);
  sIniTo(&s_aux_info, 64*MXSYM);
  sIniTo(&_gbuf, 1024);
  sIni(&_mv);
  sIniTo(&firstErr, SBUF_MXBUF);
  firstErrD=0;

  sIniTo(&s_inits, MXSYM);

  lineIni(&sbPm);
  lineIni(&sbPmDt);
  lineIni(&sbNrmL);
  lineIni(&depotLines);
  lineIni(&centralLines);
  lineIni(&_dupStrs);

  lineIni(&(tb.ss));
  lineIni(&(tb.de));

  tb.lh		= Calloc(MXSYM, int);
  tb.ini	= Calloc(MXSYM, int);
  tb.mtime	= Calloc(MXSYM, int);
  tb.iniv	= Calloc(MXSYM, double);
  tb.ini0	= Calloc(MXSYM, int);
  tb.di		= Calloc(MXDER, int);
  tb.idi	= Calloc(MXDER, int);
  tb.idu	= Calloc(MXDER, int);
  tb.lag	= Calloc(MXSYM, int);
  tb.dvid	= Calloc(MXDER, int);
  tb.thread     = 1; // Thread safe flag
  tb.dvidn      = 0;
  tb.ix		= 0;
  tb.id         = 0;
  tb.fn		= 0;
  tb.ixL        = -1;
  tb.didEq      = 0;
  tb.NEnd       = -1;
  tb.pos_de	= 0;
  tb.ini_i	= 0;
  tb.statei	= 0;
  tb.nExtra     = 0;
  tb.sensi	= 0;
  tb.li		= 0;
  tb.sli	= 0;
  tb.pi		= 0;
  tb.isPi       = 0;
  tb.isNA       = 0;
  tb.linCmt     = 0;
  tb.linCmtN    = -100;
  tb.linCmtFlg  = 0;
  tb.df		= Calloc(MXSYM, int);
  tb.dy		= Calloc(MXSYM, int);
  tb.sdfdy	= Calloc(MXSYM, int);
  tb.cdf	= 0;
  tb.ndfdy	= 0;
  tb.maxtheta   = 0;
  tb.hasCmt     = 0;
  tb.maxeta     = 0;
  tb.hasDepot   = 0;
  tb.hasCentral = 0;
  tb.hasDepotCmt = 0;
  tb.hasCentralCmt = 0;
  tb.hasKa      = 0;
  tb.allocS	= MXSYM;
  tb.allocD	= MXDER;
  tb.matn	= 0;
  tb.matnf	= 0;
  tb.ncmt	= 0;
  tb.linB	= 0;
  tb.curPropN	= 0;
  tb.depotN	= -1;
  tb.centralN	= -1;
  tb.linExtra   = false;
  tb.nwhile     = 0;
  tb.nInd       = 0;
  tb.simflg     = 0;
  // Reset Arrays
  // Reset integers
  NV		= 0;

  // reset globals
  good_jac = 1;
  found_jac = 0;
  rx_syntax_error = 0;
  rx_suppress_syntax_info=0;
  rx_podo = 0;
  rx_syntax_assign = 0;
  rx_syntax_star_pow = 0;
  rx_syntax_require_semicolon = 0;
  rx_syntax_allow_dots = 0;
  rx_syntax_allow_ini0 = 1;
  rx_syntax_allow_ini = 1;

  maxSumProdN = 0;
  SumProdLD = 0;

  foundDur=0;
  foundF0=0;
  nmtime=0;
  syntaxErrorExtra=0;
  lastSyntaxErrorLine=0;
  foundF=0;
  foundLag=0;
  foundRate=0;
  gBufLast=0;
  lastStrLoc=0;
  lastSyntaxErrorLine=0;
  needSort=0;
  nmtime=0;
  syntaxErrorExtra=0;
  extraCmt=0;
  gnini = 0;
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
      Rf_errorcall(R_NilValue, _("IO error writing parsed C file"));
    } else{
      totalWritten += written; // add the written bytes
    }
  }
  if (totalWritten != sbb->o) {
    fclose(fp);
    Rf_errorcall(R_NilValue, _("IO error writing parsed C file"));
  }
}
static void rxSyntaxError(struct D_Parser *ap);

void trans_internal(const char* parse_file, int isStr){
  char *buf1, *buf2, bufe[2048];
  int i,j,found,islhs;
  freeP();
  curP = new_D_Parser(&parser_tables_RxODE, sizeof(D_ParseNode_User));
  curP->save_parse_tree = 1;
  curP->error_recovery = 1;
  curP->initial_scope = NULL;
  curP->syntax_error_fn = rxSyntaxError;
  if (isStr){
    if (gBufFree) Free(gBuf);
    // Should be able to use gBuf directly, but I believe it cause
    // problems with R's garbage collection, so duplicate the string.
    gBuf = (char*)(parse_file);
    gBufFree=0;
  } else {
    if (gBufFree) Free(gBuf);
    gBuf = rc_sbuf_read(parse_file);
    gBufFree=1;
    err_msg((intptr_t) gBuf, "error: empty buf for FILE_to_parse\n", -2);
  }
  sFree(&sbNrm);
  sIniTo(&sbNrm, SBUF_MXBUF);
  lineIni(&sbPm);
  lineIni(&sbPmDt);
  lineIni(&sbNrmL);
  // do not free these, they remain until next parse for quick parsing of linCmt() models
  lineIni(&depotLines);
  lineIni(&centralLines);

  _pn= dparse(curP, gBuf, (int)strlen(gBuf));
  if (!_pn || curP->syntax_errors) {
    rx_syntax_error = 1;
  } else {
    wprint_parsetree(parser_tables_RxODE, _pn, 0, wprint_node, NULL);
    // Determine Jacobian vs df/dvar
    for (i=0; i<tb.ndfdy; i++) {                     /* name state vars */
      buf1=tb.ss.line[tb.df[i]];
      found=0;
      for (j=0; j<tb.de.n; j++) {                     /* name state vars */
        buf2=tb.ss.line[tb.di[j]];
	if (!strcmp(buf1, buf2)){
	  found=1;
          break;
	}
      }
      if (!found){
	buf2=tb.ss.line[tb.dy[i]];
	snprintf(bufe, 2048, NOSTATE,buf1,buf2,buf1);
	trans_syntax_error_report_fn0(bufe);
      }
      // Now the dy()
      buf1=tb.ss.line[tb.dy[i]];
      found=0;
      for (j=0; j<tb.de.n; j++) {                     /* name state vars */
        buf2=tb.ss.line[tb.di[j]];
        if (!strcmp(buf1, buf2)){
          found=1;
          break;
        }
      }
      if (!found){
	for (j=0; j<NV; j++) {
          islhs = tb.lh[j];
	  buf2=tb.ss.line[j];
          if (islhs>1 && tb.lh[i] != isLhsStateExtra) continue; /* is a state var */
          buf2=tb.ss.line[j];
          if ((islhs != 1 || tb.ini[j] == 1) &&!strcmp(buf1, buf2)){
	    found=1;
	    // This is a df(State)/dy(Parameter)
	    tb.sdfdy[i] = 1;
	    break;
	  }
        }
      }
      if (!found){
        buf2=tb.ss.line[tb.df[i]];
      	buf2=tb.ss.line[tb.dy[i]];
      	snprintf(bufe,2048,NOSTATEVAR,buf1,buf2,buf2);
        trans_syntax_error_report_fn0(bufe);
      }
    }
  }
}

SEXP _RxODE_rxQs(SEXP);
SEXP _RxODE_rxQr(SEXP);

SEXP _RxODE_trans(SEXP parse_file, SEXP prefix, SEXP model_md5, SEXP parseStr,
		  SEXP isEscIn, SEXP inME, SEXP goodFuns){
  _goodFuns = goodFuns;
  const char *in = NULL;
  char *buf, *df, *dy;

  int i, j, islhs, pi=0, li=0, sli = 0, ini_i = 0,k=0, m=0, p=0;
  // Make sure buffers are initialized.
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

  if (isString(prefix) && length(prefix) == 1){
    model_prefix = CHAR(STRING_ELT(prefix,0));
  } else {
    Rf_errorcall(R_NilValue, _("model prefix must be specified"));
  }

  if (isString(inME) && length(inME) == 1){
    me_code = CHAR(STRING_ELT(inME,0));
  } else {
    freeP();
    Rf_errorcall(R_NilValue, _("extra ME code must be specified"));
  }

  if (isString(model_md5) && length(model_md5) == 1){
    md5 = CHAR(STRING_ELT(model_md5,0));
    badMd5 = 0;
    if (strlen(md5)!= 32){
      badMd5=1;
    }
  } else {
    badMd5=1;
  }

  in = CHAR(STRING_ELT(parse_file,0));
  trans_internal(in, isStr);
  extraCmt = 0;
  if (tb.linCmt){
    if (tb.hasKa){
      extraCmt=2;
    } else {
      extraCmt=1;
    }
    if (tb.hasDepotCmt){
      trans_syntax_error_report_fn0(_("'cmt(depot)' does not work with 'linCmt()'"));
    }
    if (tb.hasCentralCmt) {
      trans_syntax_error_report_fn0("'cmt(central)' does not work with 'linCmt()'");
    }
  } else {
    if (tb.hasDepot && rx_syntax_require_ode_first){
      sPrint(&_bufw2, ODEFIRST, "depot");
      trans_syntax_error_report_fn0(_bufw2.s);
    } else if (tb.hasCentral && rx_syntax_require_ode_first){
      sPrint(&_bufw2, ODEFIRST, "depot");
      trans_syntax_error_report_fn0(_bufw2.s);
    }
  }
  for (i=0; i<NV; i++) {
    islhs = tb.lh[i];
    if (islhs>1 && islhs != isLhsStateExtra && islhs != isLHSparam && islhs != isSuppressedLHS) continue;      /* is a state var */
    if (islhs == isSuppressedLHS){
      sli++;
    } else if (islhs == isLHS || islhs == isLhsStateExtra || islhs == isLHSparam){
      li++;
      if (islhs == isLHSparam) pi++;
    } else {
      pi++;
    }
  }
  tb.pi=pi;
  tb.li=li;
  tb.sli=sli;

  int pro = 0;
  SEXP lst   = PROTECT(allocVector(VECSXP, 20));pro++;
  SEXP names = PROTECT(allocVector(STRSXP, 20));pro++;

  SEXP sNeedSort = PROTECT(allocVector(INTSXP,1));pro++;
  int *iNeedSort  = INTEGER(sNeedSort);
  iNeedSort[0] = needSort;

  SEXP sLinCmt = PROTECT(allocVector(INTSXP,11));pro++;
  INTEGER(sLinCmt)[0] = tb.ncmt;
  INTEGER(sLinCmt)[1] = tb.hasKa;
  INTEGER(sLinCmt)[2] = tb.linB;
  INTEGER(sLinCmt)[3] = tb.maxeta;
  INTEGER(sLinCmt)[4] = tb.maxtheta;
  INTEGER(sLinCmt)[6] = tb.linCmtN;
  INTEGER(sLinCmt)[7] = tb.linCmtFlg;
  INTEGER(sLinCmt)[8] = tb.nInd;
  INTEGER(sLinCmt)[9] = tb.simflg;
  INTEGER(sLinCmt)[10]= tb.thread;

  SEXP sLinCmtN = PROTECT(allocVector(STRSXP, 11));pro++;
  SET_STRING_ELT(sLinCmtN, 0, mkChar("ncmt"));
  SET_STRING_ELT(sLinCmtN, 1, mkChar("ka"));
  SET_STRING_ELT(sLinCmtN, 2, mkChar("linB"));
  SET_STRING_ELT(sLinCmtN, 3, mkChar("maxeta"));
  SET_STRING_ELT(sLinCmtN, 4, mkChar("maxtheta"));
  SET_STRING_ELT(sLinCmtN, 5, mkChar("hasCmt"));
  SET_STRING_ELT(sLinCmtN, 6, mkChar("linCmt"));
  SET_STRING_ELT(sLinCmtN, 7, mkChar("linCmtFlg"));
  SET_STRING_ELT(sLinCmtN, 8, mkChar("nIndSim"));
  SET_STRING_ELT(sLinCmtN, 9, mkChar("simflg"));
  SET_STRING_ELT(sLinCmtN, 10, mkChar("thread"));
  setAttrib(sLinCmt,   R_NamesSymbol, sLinCmtN);

  SEXP sMtime = PROTECT(allocVector(INTSXP,1));pro++;
  int *iMtime  = INTEGER(sMtime);
  iMtime[0] = (int)nmtime;

  SEXP tran  = PROTECT(allocVector(STRSXP, 22));pro++;
  SEXP trann = PROTECT(allocVector(STRSXP, 22));pro++;

  int offCmt=0,nExtra = 0;
  for (int i = 0; i < tb.statei; i++){
    if (offCmt == 0 && tb.idu[i] == 0){
      offCmt = 1;
      nExtra++;
      buf=tb.ss.line[tb.di[i]];
    } else if (offCmt == 1 && tb.idu[i] == 1){
      // There is an compartment that doesn't have a derivative
      if (tb.linCmt == 0){
	UNPROTECT(pro);
	char *v = rc_dup_str(buf, 0);
	sprintf(buf, "compartment '%s' needs differential equations defined", v);
	updateSyntaxCol();
	trans_syntax_error_report_fn(buf);
      } else if (!strcmp("depot", buf) || !strcmp("central", buf)) {
      } else {
	UNPROTECT(pro);
	char *v = rc_dup_str(buf, 0);
	sprintf(buf, _("compartment '%s' needs differential equations defined"), v);
	updateSyntaxCol();
	trans_syntax_error_report_fn(buf);
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

  SEXP dfdy = PROTECT(allocVector(STRSXP,tb.ndfdy));pro++;

  SEXP params = PROTECT(allocVector(STRSXP, tb.pi));pro++;
  SEXP lhs    = PROTECT(allocVector(STRSXP, tb.li));pro++;
  SEXP slhs   = PROTECT(allocVector(STRSXP, tb.sli));pro++;

  SEXP inin  = PROTECT(allocVector(STRSXP, tb.isPi + tb.ini_i));pro++;
  SEXP ini   = PROTECT(allocVector(REALSXP, tb.isPi + tb.ini_i));pro++;
  for (int i=tb.isPi + tb.ini_i;i--;) REAL(ini)[i] = NA_REAL;

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
  for (i = 0; i < NV; i++){
    buf=tb.ss.line[i];
    if (tb.ini[i] == 1 && tb.lh[i] != isLHS){
      if (tb.isPi && !strcmp("pi", buf)) {
	redo=1;
	tb.isPi=0;
	break;
      }
      sAppend(&s_aux_info,"    SET_STRING_ELT(inin,%d,mkChar(\"%s\"));\n",ini_i,
	      buf);
      SET_STRING_ELT(inin,ini_i,mkChar(buf));
      if (ISNA(tb.iniv[i])){
	sAppend(&s_aux_info,"    REAL(ini)[%d] = NA_REAL;\n",ini_i);
      } else if (ISNAN(tb.iniv[i])){
	sAppend(&s_aux_info,"    REAL(ini)[%d] = R_NaN;\n",ini_i);
      } else if (R_FINITE(tb.iniv[i])){
	sAppend(&s_aux_info,"    REAL(ini)[%d] = %.16f;\n",ini_i, tb.iniv[i]);
      } else if (tb.iniv[i] > 0){
	sAppend(&s_aux_info,"    REAL(ini)[%d] = R_PosInf;\n",ini_i);
      } else {
	sAppend(&s_aux_info,"    REAL(ini)[%d] = R_NegInf;\n",ini_i);
      }
      REAL(ini)[ini_i++] = tb.iniv[i];
    }
  }
  if (tb.isPi){
    sAppend(&s_aux_info,"    SET_STRING_ELT(inin,%d,mkChar(\"pi\"));\n",ini_i);
    sAppend(&s_aux_info,"    REAL(ini)[%d] = M_PI;\n",ini_i);
    SET_STRING_ELT(inin,ini_i,mkChar("pi"));
    REAL(ini)[ini_i++] = M_PI;
  } else if (redo){
    inin  = PROTECT(allocVector(STRSXP, tb.ini_i));pro++;
    ini   = PROTECT(allocVector(REALSXP, tb.ini_i));pro++;
    for (int i = tb.ini_i; i--;) REAL(ini)[i] = NA_REAL;
    ini_i=0;
    for (i = 0; i < NV; i++){
      buf=tb.ss.line[i];
      if (tb.ini[i] == 1 && tb.lh[i] != isLHS){
	if (tb.isPi && !strcmp("pi", buf)) {
	  redo=1;
	  tb.isPi=0;
	  break;
	}
	sAppend(&s_aux_info,"    SET_STRING_ELT(inin,%d,mkChar(\"%s\"));\n",ini_i,
		buf);
	SET_STRING_ELT(inin,ini_i,mkChar(buf));
	if (ISNA(tb.iniv[i])){
	  sAppend(&s_aux_info,"    REAL(ini)[%d] = NA_REAL;\n",ini_i);
	} else if (ISNAN(tb.iniv[i])){
	  sAppend(&s_aux_info,"    REAL(ini)[%d] = R_NaN;\n",ini_i);
	} else if (R_FINITE(tb.iniv[i])){
	  sAppend(&s_aux_info,"    REAL(ini)[%d] = %.16f;\n",ini_i, tb.iniv[i]);
	} else if (tb.iniv[i] > 0){
	  sAppend(&s_aux_info,"    REAL(ini)[%d] = R_PosInf;\n",ini_i);
	} else {
	  sAppend(&s_aux_info,"    REAL(ini)[%d] = R_NegInf;\n",ini_i);
	}
	REAL(ini)[ini_i++] = tb.iniv[i];
      }
    }
  }
  tb.ini_i = ini_i;

  setAttrib(ini,   R_NamesSymbol, inin);

  SEXP model  = PROTECT(allocVector(STRSXP,2));pro++;
  SEXP modeln = PROTECT(allocVector(STRSXP,2));pro++;
  k=0;j=0;m=0,p=0;
  for (i=0; i<tb.de.n; i++) {                     /* name state vars */
    buf=tb.ss.line[tb.di[i]];
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
    } else {
      SET_STRING_ELT(extraState, p++, mkChar(buf));
    }
  }
  for (i=0; i<tb.ndfdy; i++) {                     /* name state vars */
    df=tb.ss.line[tb.df[i]];
    dy=tb.ss.line[tb.dy[i]];
    int foundIt=0;
    for (j = 1; j <= tb.maxtheta;j++){
      sPrint(&_bufw,"_THETA_%d_",j);
      if (!strcmp(dy,_bufw.s)){
        sPrint(&_bufw,"THETA[%d]",j);
	foundIt=1;
	break;
      }
    }
    if (!foundIt){
      for (j = 1; j <= tb.maxeta;j++){
	sPrint(&_bufw,"_ETA_%d_",j);
	if (!strcmp(dy,_bufw.s)){
	  sPrint(&_bufw,"ETA[%d]",j);
	}
      }
    }
    if (!foundIt){
      sClear(&_bufw);
      sPrint(&_bufw,"%s",dy);
    }
    sPrint(&_bufw2,"df(%s)/dy(%s)",df,_bufw.s);
    SET_STRING_ELT(dfdy,i,mkChar(_bufw2.s));
  }
  li=0, pi=0, sli = 0;
  for (i=0; i<NV; i++) {
    islhs = tb.lh[i];
    if (islhs == isSuppressedLHS){
      SET_STRING_ELT(slhs, sli++, mkChar(tb.ss.line[i]));
    }
    if (islhs>1 && islhs != isLhsStateExtra && islhs != isLHSparam) {
      if (tb.lag[i] != 0){
	buf=tb.ss.line[i];
	if (islhs == isState){
	  sPrint(&_bufw, _("state '%s': 'lag', 'lead', 'first', 'last', 'diff' not legal"), buf);
	  trans_syntax_error_report_fn0(_bufw.s);
	} else if (islhs == 10 || islhs == 11){
	  sPrint(&_bufw, _("suppress '%s': 'lag', 'lead', 'first', 'last', 'diff' not legal"), buf);
	  trans_syntax_error_report_fn0(_bufw.s);
	}
      }
      continue;
    }      /* is a state var */
    buf=tb.ss.line[i];
    if (tb.lag[i] != 0){
      if (islhs == 70){
	sPrint(&_bufw, _("redefined '%s': 'lag', 'lead', 'first', 'last', 'diff' not legal"), buf);
	trans_syntax_error_report_fn0(_bufw.s);
      } else if (islhs == 1 && tb.lag[i] != 1){
	sPrint(&_bufw, _("lhs '%s': only 'lag(%s,1)' and 'diff(%s,1)' supported"), buf, buf, buf);
	trans_syntax_error_report_fn0(_bufw.s);
      }
    }
    if (islhs == 1 || islhs == 19 || islhs == 70){
      SET_STRING_ELT(lhs, li++, mkChar(buf));
      if (islhs == 70) {
	if (!strcmp("CMT", buf)) {
	  tb.hasCmt = 1;
	}
	SET_STRING_ELT(params, pi++, mkChar(buf));
      }
    } else {
      int foundIt=0;
      for (j = 1; j <= tb.maxtheta;j++){
	sPrint(&_bufw,"_THETA_%d_",j);
	if (!strcmp(buf, _bufw.s)){
	  sPrint(&_bufw,"THETA[%d]",j);
	  foundIt=1;
	  break;
	}
      }
      if (!foundIt){
	for (j = 1; j <= tb.maxeta;j++){
	  sPrint(&_bufw,"_ETA_%d_",j);
	  if (!strcmp(buf, _bufw.s)){
	    sPrint(&_bufw,"ETA[%d]",j);
	    foundIt=1;
	    break;
	  }
	}
      }
      if (!foundIt){
	sPrint(&_bufw, "%s", buf);
      }
      if (!strcmp("CMT", _bufw.s)) {
	tb.hasCmt = 1;
      }
      SET_STRING_ELT(params, pi++, mkChar(_bufw.s));
    }
  }
  INTEGER(sLinCmt)[5] = tb.hasCmt;
  tb.ini_i = length(ini);
  sPrint(&s_inits,"%s", s_aux_info.s);
  gnini = length(ini);

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

  SET_STRING_ELT(names,9,mkChar("state.ignore"));
  SET_VECTOR_ELT(lst,  9,stateRmS);

  SET_STRING_ELT(names,10,mkChar("version"));
  SET_VECTOR_ELT(lst,  10,version);

  SET_STRING_ELT(names,11,mkChar("normal.state"));
  SET_VECTOR_ELT(lst,  11,normState);

  SET_STRING_ELT(names,12,mkChar("needSort"));
  SET_VECTOR_ELT(lst,  12,sNeedSort);

  SET_STRING_ELT(names,13,mkChar("nMtime"));
  SET_VECTOR_ELT(lst,  13,sMtime);

  SET_STRING_ELT(names, 14, mkChar("extraCmt"));
  SEXP sExtraCmt = PROTECT(allocVector(INTSXP,1));pro++;
  INTEGER(sExtraCmt)[0] = extraCmt;
  SET_VECTOR_ELT(lst, 14, sExtraCmt);

  SET_STRING_ELT(names, 15, mkChar("stateExtra"));
  SET_VECTOR_ELT(lst,   15, extraState);

  SET_STRING_ELT(names, 16, mkChar("dvid"));
  SEXP sDvid = PROTECT(allocVector(INTSXP,tb.dvidn));pro++;
  for (i = 0; i < tb.dvidn; i++) INTEGER(sDvid)[i]=tb.dvid[i];
  SET_VECTOR_ELT(lst,  16, sDvid);

  SET_STRING_ELT(names, 17, mkChar("indLin"));
  SEXP matLst = PROTECT(allocVector(VECSXP, 0));pro++;
  SET_VECTOR_ELT(lst,  17, matLst);

  SET_STRING_ELT(names, 18, mkChar("flags"));
  SET_VECTOR_ELT(lst,   18, sLinCmt);

  SET_STRING_ELT(names, 19, mkChar("slhs"));
  SET_VECTOR_ELT(lst,   19, slhs);

  sPrint(&_bufw,"%.*s", (int)strlen(model_prefix)-1, model_prefix);

  SET_STRING_ELT(trann,0,mkChar("lib.name"));
  SET_STRING_ELT(tran,0,mkChar(_bufw.s));

  SET_STRING_ELT(trann,1,mkChar("jac"));
  if (found_jac == 1 && good_jac == 1){
    SET_STRING_ELT(tran,1,mkChar("fulluser")); // Full User Matrix
  } else {
    SET_STRING_ELT(tran,1,mkChar("fullint")); // Full Internal Matrix
  }

  SET_STRING_ELT(trann,2,mkChar("prefix"));
  SET_STRING_ELT(tran,2,mkChar(_bufw.s));

  sPrint(&_bufw,"%sdydt",model_prefix);
  SET_STRING_ELT(trann,3,mkChar("dydt"));
  SET_STRING_ELT(tran,3,mkChar(_bufw.s)) ;

  sPrint(&_bufw,"%scalc_jac",model_prefix);
  SET_STRING_ELT(trann,4,mkChar("calc_jac"));
  SET_STRING_ELT(tran, 4,mkChar(_bufw.s));

  sPrint(&_bufw,"%scalc_lhs",model_prefix);
  SET_STRING_ELT(trann,5,mkChar("calc_lhs"));
  SET_STRING_ELT(tran, 5,mkChar(_bufw.s));

  sPrint(&_bufw,"%smodel_vars",model_prefix);
  SET_STRING_ELT(trann,6,mkChar("model_vars"));
  SET_STRING_ELT(tran, 6,mkChar(_bufw.s));

  sPrint(&_bufw,"%stheta",model_prefix);
  SET_STRING_ELT(trann,7,mkChar("theta"));
  SET_STRING_ELT(tran, 7,mkChar(_bufw.s));

  sPrint(&_bufw,"%sinis",model_prefix);
  SET_STRING_ELT(trann,8,mkChar("inis"));
  SET_STRING_ELT(tran, 8,mkChar(_bufw.s));

  sPrint(&_bufw,"%sdydt_lsoda",model_prefix);
  SET_STRING_ELT(trann,9,mkChar("dydt_lsoda"));
  SET_STRING_ELT(tran, 9,mkChar(_bufw.s));

  sPrint(&_bufw,"%scalc_jac_lsoda",model_prefix);
  SET_STRING_ELT(trann,10,mkChar("calc_jac_lsoda"));
  SET_STRING_ELT(tran, 10,mkChar(_bufw.s));

  sPrint(&_bufw,"%sode_solver_solvedata",model_prefix);
  SET_STRING_ELT(trann,11,mkChar("ode_solver_solvedata"));
  SET_STRING_ELT(tran, 11,mkChar(_bufw.s));

  sPrint(&_bufw,"%sode_solver_get_solvedata",model_prefix);
  SET_STRING_ELT(trann,12,mkChar("ode_solver_get_solvedata"));
  SET_STRING_ELT(tran, 12,mkChar(_bufw.s));

  sPrint(&_bufw,"%sdydt_liblsoda",model_prefix);
  SET_STRING_ELT(trann,13,mkChar("dydt_liblsoda"));
  SET_STRING_ELT(tran, 13,mkChar(_bufw.s));

  sPrint(&_bufw,"%sF",model_prefix);
  SET_STRING_ELT(trann,14,mkChar("F"));
  SET_STRING_ELT(tran, 14,mkChar(_bufw.s));

  sPrint(&_bufw,"%sLag",model_prefix);
  SET_STRING_ELT(trann,15,mkChar("Lag"));
  SET_STRING_ELT(tran, 15,mkChar(_bufw.s));

  sPrint(&_bufw,"%sRate",model_prefix);
  SET_STRING_ELT(trann,16,mkChar("Rate"));
  SET_STRING_ELT(tran, 16,mkChar(_bufw.s));

  sPrint(&_bufw,"%sDur",model_prefix);
  SET_STRING_ELT(trann,17,mkChar("Dur"));
  SET_STRING_ELT(tran, 17,mkChar(_bufw.s));

  sPrint(&_bufw,"%smtime",model_prefix);
  SET_STRING_ELT(trann,18,mkChar("mtime"));
  SET_STRING_ELT(tran, 18,mkChar(_bufw.s));

  sPrint(&_bufw,"%sassignFuns",model_prefix);
  SET_STRING_ELT(trann,19,mkChar("assignFuns"));
  SET_STRING_ELT(tran, 19,mkChar(_bufw.s));

  sPrint(&_bufw,"%sME",model_prefix);
  SET_STRING_ELT(trann,20,mkChar("ME"));
  SET_STRING_ELT(tran, 20,mkChar(_bufw.s));

  sPrint(&_bufw,"%sIndF",model_prefix);
  SET_STRING_ELT(trann,21,mkChar("IndF"));
  SET_STRING_ELT(tran, 21,mkChar(_bufw.s));

  SET_STRING_ELT(modeln,0,mkChar("normModel"));
  SET_STRING_ELT(model,0,mkChar(sbNrm.s));

  SET_STRING_ELT(modeln,1,mkChar("indLin"));
  SET_STRING_ELT(model,1,mkChar(me_code));

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
	RSprintf("\n:%03d: ", lastSyntaxErrorLine);
	for (; gBuf[gBufLast] != '\0'; gBufLast++){
	  if (gBuf[gBufLast] == '\n'){
	    RSprintf("\n:%03d: ", ++lastSyntaxErrorLine);
	  } else{
	    RSprintf("%c", gBuf[gBufLast]);
	  }
	}
      }
      if (isEsc){
	RSprintf("\n\033[1m================================================================================\033[0m\n");
      }
      else {
	RSprintf("\n================================================================================\n");
      }
    }
    if (firstErrD == 1) {
      firstErrD=0;
      Rf_errorcall(R_NilValue, firstErr.s);
    } else {
      Rf_errorcall(R_NilValue, _("syntax errors (see above)"));
    }
  }
  return lst;
}

SEXP _RxODE_parseModel(SEXP type){
  if (!sbPm.o){
    Rf_errorcall(R_NilValue, _("model no longer loaded in memory"));
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
		    SEXP pMd5, SEXP timeId, SEXP mvLast){
  if (!sbPm.o || !sbNrm.o){
    Rf_errorcall(R_NilValue, _("nothing in output queue to write"));
  }
  if (!isString(c_file) || length(c_file) != 1){
    Rf_errorcall(R_NilValue, _("c_file should only be 1 file"));
  }
  if (length(libname) != 2){
    Rf_errorcall(R_NilValue, _("libname needs 2 elements"));
  }
  fpIO = fopen(CHAR(STRING_ELT(c_file,0)), "wb");
  err_msg((intptr_t) fpIO, "error opening output c file\n", -2);
  if (badMd5){
    SET_STRING_ELT(VECTOR_ELT(mvLast, RxMv_md5), 0, mkChar(""));
  } else {
    SET_STRING_ELT(VECTOR_ELT(mvLast, RxMv_md5), 0, mkChar(md5));
  }
  SET_STRING_ELT(VECTOR_ELT(mvLast, RxMv_model), 1, mkChar(me_code));
  int pro = 0;
  SEXP trans = PROTECT(VECTOR_ELT(mvLast, RxMv_trans)); pro++;
  sbuf buf;
  sIni(&buf);
  if (strcmp(CHAR(STRING_ELT(trans, 0)), CHAR(STRING_ELT(libname, 0)))) {
    SET_STRING_ELT(trans, 0, STRING_ELT(libname, 0)); // libname
    SET_STRING_ELT(trans, 2, STRING_ELT(prefix, 0)); // prefix
    const char *curPrefix = CHAR(STRING_ELT(prefix,0));
    sPrint(&buf, "%sdydt", curPrefix);
    SET_STRING_ELT(trans, 3, mkChar(buf.s)); // dydt
    sPrint(&buf, "%scalc_jac", curPrefix);
    SET_STRING_ELT(trans, 4, mkChar(buf.s)); // calc_jac
    sPrint(&buf, "%scalc_lhs", curPrefix);
    SET_STRING_ELT(trans, 5, mkChar(buf.s)); // calc_lhs
    sPrint(&buf, "%smodel_vars", curPrefix);
    SET_STRING_ELT(trans, 6, mkChar(buf.s)); // model_vars
    sPrint(&buf, "%stheta", curPrefix);
    SET_STRING_ELT(trans, 7, mkChar(buf.s)); // theta
    sPrint(&buf, "%sinis", curPrefix);
    SET_STRING_ELT(trans, 8, mkChar(buf.s)); // inis
    sPrint(&buf, "%sdydt_lsoda", curPrefix);
    SET_STRING_ELT(trans, 9, mkChar(buf.s)); // dydt_lsoda
    sPrint(&buf, "%scalc_jac_lsoda", curPrefix);
    SET_STRING_ELT(trans, 10, mkChar(buf.s)); // calc_jac_lsoda
    sPrint(&buf, "%sode_solver_solvedata", curPrefix);
    SET_STRING_ELT(trans, 11, mkChar(buf.s)); // ode_solver_solvedata
    sPrint(&buf, "%sode_solver_get_solvedata", curPrefix);
    SET_STRING_ELT(trans, 12, mkChar(buf.s)); // ode_solver_get_solvedata
    sPrint(&buf, "%sdydt_liblsoda", curPrefix);
    SET_STRING_ELT(trans, 13, mkChar(buf.s)); // dydt_liblsoda
    sPrint(&buf, "%sF", curPrefix);
    SET_STRING_ELT(trans, 14, mkChar(buf.s)); // F
    sPrint(&buf, "%sLag", curPrefix);
    SET_STRING_ELT(trans, 15, mkChar(buf.s)); // Lag
    sPrint(&buf, "%sRate", curPrefix);
    SET_STRING_ELT(trans, 16, mkChar(buf.s)); // Rate
    sPrint(&buf, "%sDur", curPrefix);
    SET_STRING_ELT(trans, 17, mkChar(buf.s)); // Dur
    sPrint(&buf, "%smtime", curPrefix);
    SET_STRING_ELT(trans, 18, mkChar(buf.s)); // mtime
    sPrint(&buf, "%sassignFuns", curPrefix);
    SET_STRING_ELT(trans, 19, mkChar(buf.s)); // assignFuns
    sPrint(&buf, "%sME", curPrefix);
    SET_STRING_ELT(trans, 20, mkChar(buf.s)); // ME
    sPrint(&buf, "%sIndF", curPrefix);
    SET_STRING_ELT(trans, 21, mkChar(buf.s)); // IndF
  }
  sPrint(&_mv, "%s", CHAR(STRING_ELT(PROTECT(_RxODE_rxQs(mvLast)), 0))); pro++;
  UNPROTECT(pro);
  sFree(&buf);
  //SET_STRING_ELT(tran, 0, mkChar());
  sFree(&sbOut);
  sIniTo(&sbOut, (int)((sbPm.sN)*5.3));
  // show_ode = 1 dydt
  // show_ode = 2 Jacobian
  // show_ode = 3 Ini statement
  // show_ode = 0 LHS
  // show_ode = 5 functional bioavailibility
  // show_ode = 6 functional rate
  if (tb.linCmt != 0) {
    char *buf;
    int badCentral=false, badDepot=false;
    for (int i=tb.de.n; i--;) {                     /* name state vars */
      buf=tb.ss.line[tb.di[i]];
      if (tb.hasKa == 1 && !strcmp(buf,"depot")){
	badDepot=true;
      } else if (!strcmp(buf, "central")) {
	badCentral=true;
      }
    }
    if (badCentral && badDepot){
      fclose(fpIO);
      reset();
      Rf_errorcall(R_NilValue, _("linCmt() and ode have 'central' and 'depot' compartments, rename ODE 'central'/'depot'"));
    } else if (badCentral) {
      fclose(fpIO);
      reset();
      Rf_errorcall(R_NilValue, _("linCmt() and ode has a 'central' compartment, rename ODE 'central'"));
    } else if (badDepot) {
      fclose(fpIO);
      reset();
      Rf_errorcall(R_NilValue, _("linCmt() and ode has a 'depot' compartment, rename ODE 'depot'"));
    }
    (&sbOut)->s[0]='\0';
    if (tb.hasKa == 1) {
      sAppend(&sbOut, "#define _DEPOT_ %d\n", tb.statei);
      sAppend(&sbOut, "#define _CENTRAL_ %d\n", tb.statei+1);
    } else if (tb.hasCentral == 1) {
      if (tb.hasDepot){
	fclose(fpIO);
	reset();
	Rf_errorcall(R_NilValue, _("linCmt() does not have 'depot' compartment without a 'ka'"));
	return R_NilValue;
      }
      sAppend(&sbOut, "#define _CENTRAL_ %d\n", tb.statei);
    }
    writeSb(&sbOut, fpIO);
  }
  gCode(1); // d/dt()
  gCode(2); // jac
  gCode(3); // ini()
  gCode(0); //
  gCode(5);
  gCode(6);
  gCode(7);
  gCode(8);
  gCode(9); // mtime
  gCode(10); //mat
  gCode(11); //matF
  gCode(4); // Registration
  fclose(fpIO);
  reset();
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
      if (isEsc){
	RSprintf(_("\033[1mRxODE model syntax error:\n================================================================================\033[0m"));
      }
      else {
	RSprintf(_("RxODE Model Syntax Error:\n================================================================================"));
      }
      lastSyntaxErrorLine=1;
    }
    char *buf;
    Parser *p = (Parser *)ap;
    for (; lastSyntaxErrorLine < p->user.loc.line; lastSyntaxErrorLine++){
      buf = getLine(gBuf, lastSyntaxErrorLine, &gBufLast);
      RSprintf("\n:%03d: %s", lastSyntaxErrorLine, buf);
      Free(buf);

    }
    char *after = 0;
    ZNode *z = p->snode_hash.last_all ? p->snode_hash.last_all->zns.v[0] : 0;
    while (z && z->pn->parse_node.start_loc.s == z->pn->parse_node.end)
      z = (z->sns.v && z->sns.v[0]->zns.v) ? z->sns.v[0]->zns.v[0] : 0;
    if (z && z->pn->parse_node.start_loc.s != z->pn->parse_node.end)
      after = rc_dup_str(z->pn->parse_node.start_loc.s, z->pn->parse_node.end);
    if (after){
      if (isEsc){
	RSprintf(_("\n\n\033[1mRxODE syntax error after\033[0m '\033[35m\033[1m%s\033[0m':\n"),  after);
      }
      else {
	RSprintf(_("\n\nRxODE syntax error after '%s'\n"),  after);
      }
      if (firstErrD == 0) {
	sAppend(&firstErr, _("RxODE syntax error after '%s':\n"), after);
      }
    }
    else{
      if (isEsc){
	RSprintf(_("\n\n\033[1mRxODE syntax error\033[0m:\n"));
      }
      else{
	RSprintf(_("\n\nRxODE syntax error:\n"));
      }
      if (firstErrD == 0) {
	sAppendN(&firstErr, "RxODE syntax error:\n", 20);
      }
    }
    buf = getLine(gBuf, p->user.loc.line, &gBufLast);
    if (lastSyntaxErrorLine < p->user.loc.line) lastSyntaxErrorLine++;
    if (isEsc) {
      RSprintf("\033[1m:%03d:\033[0m ", p->user.loc.line);
    }
    else {
      RSprintf(":%03d: ", p->user.loc.line);
    }
    if (firstErrD == 0) {
      sAppend(&firstErr, ":%03d: ", p->user.loc.line);
    }
    int col = 0, len= strlen(buf), lenv, i;
    for (i = 0; i < p->user.loc.col; i++){
      RSprintf("%c", buf[i]);
      if (firstErrD == 0) {
	sAppend(&firstErr, "%c", buf[i]);
      }
      if (i == len-2) { i++; break;}
    }
    if (isEsc) {
      RSprintf("\033[35m\033[1m%c\033[0m", buf[i++]);
    }
    else {
      RSprintf("%c", buf[i++]);
    }
    if (firstErrD == 0) {
      sAppend(&firstErr, "%c", buf[i-1]);
    }
    for (; i < len; i++){
      RSprintf("%c", buf[i]);
      if (firstErrD == 0) {
	sAppend(&firstErr, "%c", buf[i]);
      }
    }
    RSprintf("\n      ");
    if (firstErrD == 0) {
      sAppendN(&firstErr, "\n      ", 7);
    }
    if (after){
      lenv = strlen(after);
      while (col != len && strncmp(buf + col, after, lenv) != 0) col++;
      if (col == len) col = 0;
      if (col){
	for (int i = 0; i < col; i++){
	  RSprintf(" ");
	  if (firstErrD == 0) {
	    sAppendN(&firstErr, " ", 1);
	  }
	  if (i == len-2) { i++; break;}
	}
	len = p->user.loc.col - col;
	if (len > 0 && len < 40){
	  for (int i = len; i--;) {
	    RSprintf("~");
	    if (firstErrD == 0) {
	      sAppendN(&firstErr, "~", 1);
	    }
	  }
	}
	if (isEsc) {
	  RSprintf("\033[35m\033[1m^\033[0m");
	}
	else {
	  RSprintf("^");
	}
	if (firstErrD == 0) {
	  sAppendN(&firstErr, "^", 1);
	}
      } else {
	for (int i = 0; i < p->user.loc.col; i++){
	  RSprintf(" ");
	  if (firstErrD == 0) {
	    sAppendN(&firstErr, " ", 1);
	  }
	  if (i == len-2) { i++; break;}
	}
	if (isEsc) {
	  RSprintf("\033[35m\033[1m^\033[0m");
	}
	else {
	  RSprintf("^");
	}
	if (firstErrD == 0) {
	  sAppendN(&firstErr, "^", 1);
	}

      }
    } else {
      for (int i = 0; i < p->user.loc.col; i++){
	RSprintf(" ");
	if (firstErrD == 0) {
	  sAppendN(&firstErr, " ", 1);
	}
	if (i == len-2) { i++; break;}
      }
      if (isEsc) {
	RSprintf("\033[35m\033[1m^\033[0m");
      }
      else {
	RSprintf("^");
      }
      if (firstErrD == 0) {
	sAppendN(&firstErr, "^", 1);
      }
    }
    Free(buf);
    if (firstErrD == 0) {
      firstErrD = 1;
      sAppendN(&firstErr, "\nmore errors could be listed above", 34);
    }
  }
  rx_syntax_error = 1;
}

static void trans_syntax_error_report_fn0(char *err){
  if (!rx_suppress_syntax_info){
    if (lastSyntaxErrorLine == 0){
      if (isEsc) {
	RSprintf(_("\033[1mRxODE model syntax error:\n================================================================================\033[0m"));
      }
      else {
	RSprintf(_("RxODE model syntax error:\n================================================================================"));
      }
      lastSyntaxErrorLine=1;
    }
    if (isEsc) {
      RSprintf("\n\033[1m:ERR:\033[0m %s:\n",  err);
    }
    else {
      RSprintf("\n:ERR: %s:\n", err);
    }
  }
  rx_syntax_error = 1;
}

void trans_syntax_error_report_fn(char *err) {
  if (!rx_suppress_syntax_info){
    if (lastSyntaxErrorLine == 0){
      if (isEsc) {
	RSprintf(_("\033[1mRxODE model syntax error:\n================================================================================\033[0m"));
      }
      else {
	RSprintf(_("RxODE model syntax error:\n================================================================================"));
      }
      lastSyntaxErrorLine=1;
    }
    Parser *p = (Parser *)curP;
    char *buf;
    for (; lastSyntaxErrorLine < p->user.loc.line; lastSyntaxErrorLine++){
      buf = getLine(gBuf, lastSyntaxErrorLine, &gBufLast);
      RSprintf("\n:%03d: %s", lastSyntaxErrorLine, buf);
      Free(buf);
    }
    if (lastSyntaxErrorLine < p->user.loc.line){
      RSprintf("\n");
      lastSyntaxErrorLine++;
    }
    if (isEsc) {
      RSprintf("\n\033[1m:%03d:\033[0m %s:\n", p->user.loc.line, err);
    }
    else {
      RSprintf("\n:%03d: %s:\n", p->user.loc.line, err);
    }
    buf = getLine(gBuf, p->user.loc.line, &gBufLast);
    RSprintf("      ");
    int i, len = strlen(buf);
    for (i = 0; i < p->user.loc.col; i++){
      RSprintf("%c", buf[i]);
      if (i == len-2) { i++; break;}
    }
    if (isEsc) {
      RSprintf("\033[35m\033[1m%c\033[0m", buf[i++]);
    }
    else {
      RSprintf("%c", buf[i++]);
    }
    for (; i < len; i++){
      RSprintf("%c", buf[i]);
    }
    RSprintf("\n      ");
    Free(buf);
    for (int i = 0; i < p->user.loc.col; i++){
      RSprintf(" ");
      if (i == len-2) { i++; break;}
    }
    if (isEsc) {
      RSprintf("\033[35m\033[1m^\033[0m");
    }
    else {
      RSprintf("^");
    }
    if (syntaxErrorExtra > 0 && syntaxErrorExtra < 40){
      for (int i = syntaxErrorExtra; i--;) {
	RSprintf("~");
      }
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

////////////////////////////////////////////////////////////////////////////////
// linCmtParse

// Taken from dparser and changed to use Calloc
char * rc_dup_str(const char *s, const char *e) {
  lastStr=s;
  int l = e ? e-s : (int)strlen(s);
  syntaxErrorExtra=min(l-1, 40);
  addLine(&_dupStrs, "%.*s", l, s);
  return _dupStrs.line[_dupStrs.n-1];
}
