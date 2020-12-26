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
#include "genModelVars.h"
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

#define MXSYM 50000
#define MXDER 5000
#define MXLEN 12000
/* #define SBUF_MXBUF 5 */
#define MXLINE 100
/* #define MXLINE 5 */


#define ENDLINE tb.ixL=-1; tb.didEq=0; tb.NEnd=NV;

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
vLines _dupStrs;

int rx_syntax_error = 0, rx_suppress_syntax_info=0, rx_podo = 0, rx_syntax_require_ode_first = 1;

extern D_ParserTables parser_tables_RxODE;

unsigned int found_jac = 0, nmtime=0;
int rx_syntax_assign = 0, rx_syntax_star_pow = 0,
  rx_syntax_require_semicolon = 0, rx_syntax_allow_dots = 0,
  rx_syntax_allow_ini0 = 1, rx_syntax_allow_ini = 1, rx_syntax_allow_assign_state = 0,
  maxSumProdN = 0, SumProdLD = 0, good_jac=1, extraCmt=0;

sbuf s_inits;

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

int lastSyntaxErrorLine=0;
char *getLine (char *src, int line, int *lloc);
void updateSyntaxCol();

#include "parseVars.h"

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

sbuf _bufw, _bufw2;

void parseFree(int last){
  sFree(&sb);
  sFree(&sbDt);
  sFree(&sbt);
  sFree(&sbNrm);
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
}

static void rxSyntaxError(struct D_Parser *ap);

static inline void assertCorrectDfDy() {
  char *buf1, *buf2, bufe[2048];
  int i, j, found, islhs;
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

void trans_internal(const char* parse_file, int isStr){
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
    assertCorrectDfDy();
  }
}



SEXP _RxODE_trans(SEXP parse_file, SEXP prefix, SEXP model_md5, SEXP parseStr,
		  SEXP isEscIn, SEXP inME, SEXP goodFuns){
  _goodFuns = goodFuns;
  const char *in = NULL;
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
  SEXP lst = PROTECT(generateModelVars());
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
  UNPROTECT(1);
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

void trans_syntax_error_report_fn0(char *err){
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

static inline void printSyntaxErrorHeader() {
  if (lastSyntaxErrorLine == 0){
    if (isEsc) {
      RSprintf(_("\033[1mRxODE model syntax error:\n================================================================================\033[0m"));
    }
    else {
      RSprintf(_("RxODE model syntax error:\n================================================================================"));
    }
    lastSyntaxErrorLine=1;
  }
}

static inline void printPriorLines(Parser *p) {
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
}

static inline void printErrorInfo(Parser *p, char *err, char *after, int printLine) {
  if (printLine) {
    if (isEsc) {
      RSprintf("\n\033[1m:%03d:\033[0m %s:\n", p->user.loc.line, err);
    }
    else {
      RSprintf("\n:%03d: %s:\n", p->user.loc.line, err);
    }
  } else {
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
  }
}

static inline void printErrorLineHighlightPoint(Parser *p) {
  char *buf = getLine(gBuf, p->user.loc.line, &gBufLast);
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

void trans_syntax_error_report_fn(char *err) {
  if (!rx_suppress_syntax_info){
    printSyntaxErrorHeader();
    Parser *p = (Parser *)curP;
    printPriorLines(p);
    printErrorInfo(p, err, 0, 1);
    printErrorLineHighlightPoint(p);
  }
  rx_syntax_error = 1;
}

static inline void printLineNumberAlone(Parser *p) {
  if (isEsc) {
    RSprintf("\033[1m:%03d:\033[0m ", p->user.loc.line);
  }
  else {
    RSprintf(":%03d: ", p->user.loc.line);
  }
  if (firstErrD == 0) {
    sAppend(&firstErr, ":%03d: ", p->user.loc.line);
  }
}

static inline void printErrorLineHighlight1(Parser *p, char *buf, char *after, int len) {
  int i;
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
}

static inline int printErrorLineHighligt2afterCol(Parser *p, char *buf, char *after, int len, int col) {
  if (!col || col == len) return 0;
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
  return 1;
}

static inline void printErrorLineHighligt2after(Parser *p, char *buf, char *after, int len) {
  int col = 0, lenv = strlen(after);
  while (col != len && strncmp(buf + col, after, lenv) != 0) col++;
  if (col == len) col = 0;
  if (!printErrorLineHighligt2afterCol(p, buf, after, len, col)) {
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
}

static inline void printErrorLineHighlight2(Parser *p, char *buf, char *after, int len) {
  if (after){
    printErrorLineHighligt2after(p, buf, after, len);
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
}

static inline void printErrorLineHiglightRegion(Parser *p, char *after) {
  char *buf = getLine(gBuf, p->user.loc.line, &gBufLast);
  if (lastSyntaxErrorLine < p->user.loc.line) lastSyntaxErrorLine++;
  printLineNumberAlone(p);
  int len= strlen(buf);
  printErrorLineHighlight1(p, buf, after, len);
  printErrorLineHighlight2(p, buf, after, len);
  Free(buf);
}


static void rxSyntaxError(struct D_Parser *ap) {
  if (!rx_suppress_syntax_info){
    printSyntaxErrorHeader();
    Parser *p = (Parser *)ap;
    printPriorLines(p);
    char *after = 0;
    ZNode *z = p->snode_hash.last_all ? p->snode_hash.last_all->zns.v[0] : 0;
    while (z && z->pn->parse_node.start_loc.s == z->pn->parse_node.end)
      z = (z->sns.v && z->sns.v[0]->zns.v) ? z->sns.v[0]->zns.v[0] : 0;
    if (z && z->pn->parse_node.start_loc.s != z->pn->parse_node.end)
      after = rc_dup_str(z->pn->parse_node.start_loc.s, z->pn->parse_node.end);
    printErrorInfo(p, 0, after, 0);
    printErrorLineHiglightRegion(p, after);

    if (firstErrD == 0) {
      firstErrD = 1;
      sAppendN(&firstErr, "\nmore errors could be listed above", 34);
    }
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
