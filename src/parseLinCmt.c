#include "parseLinCmt.h"

char errLin[errLinLen];
int errOff = 0;

int _linCmtParsePro=0;

SEXP _linCmtParse(SEXP vars0, SEXP inStr, SEXP verboseSXP) {
  const char *first = "linCmtB(rx__PTR__, t, ";
  const char *mid0 = "0, ";
  const char *end1 = "rx_tlag, rx_F, rx_rate, rx_dur,";
  const char *end2 = ", yrx_tlag2, rx_F2, rx_rate2, rx_dur2)";
  int verbose = 0;
  if (TYPEOF(verboseSXP) == LGLSXP) {
    verbose = INTEGER(verboseSXP)[0];
  }
  int type = TYPEOF(inStr);
  if (type == STRSXP) {
    int len = Rf_length(inStr);
    if (len > 0) {
      first = CHAR(STRING_ELT(inStr, 0));
    }
    if (len > 1) {
      mid0 = CHAR(STRING_ELT(inStr, 1));
    }
    if (len > 2) {
      end1 = CHAR(STRING_ELT(inStr, 2));
    }
    if (len > 3) {
      end2 = CHAR(STRING_ELT(inStr, 3));
    }
  }
  linCmtStruct lin;
  linCmtIni(&lin);
  lin.mid = mid0;
  lin.vars = vars0;
  for (int i = Rf_length(lin.vars); i--;){
    linCmtStr(&lin, CHAR(STRING_ELT(lin.vars, i)), &i);
  }
  linCmtAdjustPars(&lin);
  lin.trans =-1;
  lin.ncmt = -1;
  sIni(&(lin.ret0));
  sIni(&(lin.ret));
  if (lin.cl != -1) {
    lin.trans = 1;
    if (lin.vss != -1) {
      lin.ncmt = 2;
      lin.trans = 3;
      if (lin.vStyle != -1) {
	errLin[0] = '\0';
	errOff=0;
	snprintf(errLin + errOff, errLinLen-errOff, "cannot mix 'Vss' and '");
	errOff+=22;
	linCmtVStr(lin.vStyle);
	snprintf(errLin + errOff, errLinLen-errOff, "' volumes");
	errOff+=9;
	parseFree(0);
	sFree(&(lin.ret0));
	sFree(&(lin.ret));
	Rf_errorcall(R_NilValue, errLin);
      }
      if (lin.v == -1) {
	parseFree(0);
	sFree(&(lin.ret0));
	sFree(&(lin.ret));
	Rf_errorcall(R_NilValue, _("cannot figure out a central volume"));
      }
      if (lin.cl2 == -1) {
	parseFree(0);
	sFree(&(lin.ret0));
	sFree(&(lin.ret));
	Rf_errorcall(R_NilValue, _("cannot figure out distributional clearance"));
      }
      sAppend(&(lin.ret0), "%d, %s", lin.trans, lin.mid);
      sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.cl)));
      sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.v)));
      sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.cl2)));
      sAppend(&(lin.ret0), "%s, 0.0, 0.0, ", CHAR(STRING_ELT(lin.vars, lin.vss)));
    } else {
      if (lin.v == -1) {
	parseFree(0);
	sFree(&(lin.ret0));
	sFree(&(lin.ret));
	Rf_errorcall(R_NilValue, _("cannot figure out a central volume"));
      }
      lin.ncmt = 1;
      lin.trans = 1;
      sAppend(&(lin.ret0), "%d, %s", lin.trans, lin.mid);
      sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.cl)));
      sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.v)));
      if (lin.v2 != -1 || lin.cl2 != -1) {
	lin.ncmt = 2;
	if (lin.cl2 == -1) {
	  parseFree(0);
	  sFree(&(lin.ret0));
	  sFree(&(lin.ret));
	  Rf_errorcall(R_NilValue, _("cannot figure out distributional clearance"));
	}
	if (lin.v2 == -1) {
	  parseFree(0);
	  sFree(&(lin.ret0));
	  sFree(&(lin.ret));
	  Rf_errorcall(R_NilValue, _("cannot figure out distributional volume"));
	}
	sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.cl2)));
	sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.v2)));
	if (lin.v3 != -1 || lin.cl3 != -1) {
	  lin.ncmt = 3;
	  if (lin.cl3 == -1) {
	    parseFree(0);
	    sFree(&(lin.ret0));
	    sFree(&(lin.ret));
	    Rf_errorcall(R_NilValue, _("cannot figure out 2nd distributional clearance"));
	  }
	  if (lin.v3 == -1) {
	    parseFree(0);
	    sFree(&(lin.ret0));
	    sFree(&(lin.ret));
	    Rf_errorcall(R_NilValue, _("cannot figure out 2nd distributional volume"));
	  }
	  sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.cl3)));
	  sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.v3)));
	} else {
	  sAppendN(&(lin.ret0), "0.0, 0.0, ", 10);
	}
      } else {
	sAppendN(&(lin.ret0), "0.0, 0.0, 0.0, 0.0, ", 20);
      }
    }
    if (verbose) RSprintf(_("Detected %d-compartment model in terms of clearance"), lin.ncmt);
  } else if (lin.kel != -1) {
    if (lin.v == -1) {
      parseFree(0);
      sFree(&(lin.ret0));
      sFree(&(lin.ret));
      Rf_errorcall(R_NilValue, _("cannot figure out a central volume"));
    }
    lin.ncmt = 1;
    lin.trans = 2;
    sAppend(&(lin.ret0), "%d, %s", lin.trans, lin.mid);
    sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.kel)));
    sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.v)));
    if (lin.k12 != -1 || lin.k21 != -1) {
      if (lin.k12 == -1) {
	if (lin.cmtc == 1){
	  parseFree(0);
	  sFree(&(lin.ret0));
	  sFree(&(lin.ret));
	  Rf_errorcall(R_NilValue, _("'k12' not found when 'k21' present"));
	} else {
	  parseFree(0);
	  sFree(&(lin.ret0));
	  sFree(&(lin.ret));
	  Rf_errorcall(R_NilValue, _("'k23' not found when 'k32' present"));
	}
      }
      if (lin.k21 == -1) {
	if (lin.cmtc == 1){
	  parseFree(0);
	  sFree(&(lin.ret0));
	  sFree(&(lin.ret));
	  Rf_errorcall(R_NilValue, _("'k21' not found when 'k12' present"));
	} else {
	  parseFree(0);
	  sFree(&(lin.ret0));
	  sFree(&(lin.ret));
	  Rf_errorcall(R_NilValue, _("'k32' not found when 'k23' present"));
	}
      }
      lin.ncmt = 2;
      sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.k12)));
      sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.k21)));
      if (lin.k13 != -1 || lin.k31 != -1) {
	if (lin.k13 == -1) {
	  if (lin.cmtc == 1){
	    sFree(&(lin.ret0));
	    sFree(&(lin.ret));
	    parseFree(0);
	    Rf_errorcall(R_NilValue, _("'k13' not found when 'k31' present"));
	  } else {
	    sFree(&(lin.ret0));
	    sFree(&(lin.ret));
	    parseFree(0);
	    Rf_errorcall(R_NilValue, _("'k24' not found when 'k42' present"));
	  }
	}
	if (lin.k31 == -1) {
	  if (lin.cmtc == 1){
	    sFree(&(lin.ret0));
	    sFree(&(lin.ret));
	    parseFree(0);
	    Rf_errorcall(R_NilValue, _("'k31' not found when 'k13' present"));
	  } else {
	    sFree(&(lin.ret0));
	    sFree(&(lin.ret));
	    parseFree(0);
	    Rf_errorcall(R_NilValue, _("'k42' not found when 'k24' present"));
	  }
	}
	lin.ncmt = 3;
	sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.k13)));
	sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.k31)));
      } else {
	sAppendN(&(lin.ret0), "0.0, 0.0, ", 10);
      }
    } else if (lin.k31 != -1 || lin.k13 != -1){
      if (lin.cmtc == 1){
	sFree(&(lin.ret0));
	sFree(&(lin.ret));
	parseFree(0);
	Rf_errorcall(R_NilValue, _("'k13' or 'k31' present when 'k12' and 'k21' not present"));
      } else {
	sFree(&(lin.ret0));
	sFree(&(lin.ret));
	parseFree(0);
	Rf_errorcall(R_NilValue, _("'k24' or 'k42' present when 'k23' and 'k32' not present"));
      }
    } else {
      sAppendN(&(lin.ret0), "0.0, 0.0, 0.0, 0.0, ", 20);
    }
    if (verbose) RSprintf(_("detected %d-compartment model in terms of micro-constants"), lin.ncmt);
  } else if (lin.aob != -1) {
    lin.ncmt = 2;
    lin.trans = 5;
    if (lin.v == -1) {
      parseFree(0);
      sFree(&(lin.ret0));
      sFree(&(lin.ret));
      Rf_errorcall(R_NilValue, _("cannot figure out a central volume"));
    }
    if (lin.alpha == -1) {
      parseFree(0);
      sFree(&(lin.ret0));
      sFree(&(lin.ret));
      Rf_errorcall(R_NilValue, _("need an 'alpha' with 'aob'"));
    }
    if (lin.beta == -1) {
      parseFree(0);
      sFree(&(lin.ret0));
      sFree(&(lin.ret));
      Rf_errorcall(R_NilValue, _("need a 'beta' with 'aob'"));
    }
    sAppend(&(lin.ret0), "%d, %s", lin.trans, lin.mid);
    sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.alpha)));
    sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.v)));
    sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.beta)));
    sAppend(&(lin.ret0), "%s, 0.0, 0.0, ", CHAR(STRING_ELT(lin.vars, lin.aob)));
    if (verbose) RSprintf(_("detected %d-compartment model in terms of 'alpha' and 'aob'"), lin.ncmt);
  } else if (lin.k21 != -1) {
    lin.ncmt = 2;
    lin.trans = 4;
    if (lin.v == -1) {
      parseFree(0);
      sFree(&(lin.ret0));
      sFree(&(lin.ret));
      Rf_errorcall(R_NilValue, _("cannot figure out a central volume"));
    }
    if (lin.alpha == -1) {
      parseFree(0);
      sFree(&(lin.ret0));
      sFree(&(lin.ret));
      Rf_errorcall(R_NilValue, _("need an 'alpha'"));
    }
    if (lin.beta == -1) {
      parseFree(0);
      sFree(&(lin.ret0));
      sFree(&(lin.ret));
      Rf_errorcall(R_NilValue, _("need a 'beta'"));
    }
    sAppend(&(lin.ret0), "%d, %s", lin.trans, lin.mid);
    sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.alpha)));
    sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.v)));
    sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.beta)));
    sAppend(&(lin.ret0), "%s, 0.0, 0.0, ", CHAR(STRING_ELT(lin.vars, lin.k21)));
    if (verbose) {
      if (lin.cmtc == 1) {
	RSprintf(_("detected %d-compartment model in terms of 'alpha' or 'k21'"), lin.ncmt);
      } else {
	RSprintf(_("detected %d-compartment model in terms of 'alpha' or 'k32'"), lin.ncmt);
      }
    }
  } else if (lin.alpha != -1) {
    lin.ncmt = 1;
    if (lin.a != -1){
      lin.trans = 10;
    } else {
      lin.trans = 11;
    }
    sAppend(&(lin.ret0), "%d, %s", lin.trans, lin.mid);
    sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.alpha)));
    if (lin.a != -1) {
      sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.a)));
    } else {
      if (lin.v == -1) {
	parseFree(0);
	sFree(&(lin.ret0));
	sFree(&(lin.ret));
	Rf_errorcall(R_NilValue, _("cannot figure out a central volume"));
      }
      sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.v)));
    }
    if (lin.beta != -1 || lin.b != -1) {
      lin.ncmt =2;
      if (lin.beta == -1) {
	parseFree(0);
	sFree(&(lin.ret0));
	sFree(&(lin.ret));
	Rf_errorcall(R_NilValue, _("need a 'beta'"));
      }
      if (lin.b == -1) {
	parseFree(0);
	sFree(&(lin.ret0));
	sFree(&(lin.ret));
	Rf_errorcall(R_NilValue, _("need a 'b'"));
      }
      sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.beta)));
      sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.b)));
      if (lin.gamma != -1 || lin.c != -1) {
	lin.ncmt = 3;
	if (lin.gamma == -1) {
	  parseFree(0);
	  sFree(&(lin.ret0));
	  sFree(&(lin.ret));
	  Rf_errorcall(R_NilValue, _("need a 'gamma'"));
	}
	if (lin.c == -1) {
	  parseFree(0);
	  sFree(&(lin.ret0));
	  sFree(&(lin.ret));
	  Rf_errorcall(R_NilValue, _("need a 'c'"));
	}
	sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.gamma)));
	sAppend(&(lin.ret0), "%s, ", CHAR(STRING_ELT(lin.vars, lin.c)));
      } else {
	sAppendN(&(lin.ret0), "0.0, 0.0, ", 10);
      }
    } else if (lin.gamma != -1 || lin.c != -1) {
      parseFree(0);
      sFree(&(lin.ret0));
      sFree(&(lin.ret));
      Rf_errorcall(R_NilValue, _("a 'gamma' or 'c' specified without 'b' or 'beta'"));
    } else {
      sAppendN(&(lin.ret0), "0.0, 0.0, 0.0, 0.0, ", 20);
    }
    if (verbose) {
      if (lin.a != -1){
	RSprintf(_("detected %d-compartment model in terms of 'alpha' and central volume"), lin.ncmt);
      } else {
	RSprintf(_("detected %d-compartment model in terms of 'alpha' and 'a'"), lin.ncmt);
      }
    }
  }
  sAppend(&(lin.ret), "%s", first);
  sAppend(&(lin.ret), "%d, %s", lin.ncmt, lin.ret0.s);
  sAppend(&(lin.ret), "%s", end1);
  if (lin.ka == -1) {
    sAppendN(&(lin.ret), "0.0", 3);
    if (verbose) RSprintf("\n");
  } else {
    sAppend(&(lin.ret), "%s", CHAR(STRING_ELT(lin.vars, lin.ka)));
    if (verbose) RSprintf(_(" with first order absorption\n"));
  }
  sAppend(&(lin.ret), "%s", end2);
  int pro = 0;
  SEXP strV = PROTECT(allocVector(STRSXP, 1)); pro++;
  SEXP lst = PROTECT(allocVector(VECSXP, 3)); pro++;
  SEXP lstN = PROTECT(allocVector(STRSXP, 3)); pro++;

  SEXP transSXP = PROTECT(allocVector(INTSXP, 1)); pro++;
  INTEGER(transSXP)[0] = lin.trans;

  SEXP ncmtSXP = PROTECT(allocVector(INTSXP, 1)); pro++;
  INTEGER(ncmtSXP)[0] = lin.ncmt;

  SET_STRING_ELT(strV, 0, mkChar(lin.ret.s));
  SET_VECTOR_ELT(lst,  0, strV);
  SET_STRING_ELT(lstN, 0, mkChar("str"));

  SET_STRING_ELT(lstN, 1, mkChar("ncmt"));
  SET_VECTOR_ELT(lst,  1, ncmtSXP);

  SET_STRING_ELT(lstN, 2, mkChar("trans"));
  SET_VECTOR_ELT(lst,  2, transSXP);

  setAttrib(lst, R_NamesSymbol, lstN);
  sFree(&(lin.ret0));
  sFree(&(lin.ret));
  UNPROTECT(pro);
  if (lin.trans == -1) {
    UNPROTECT(_linCmtParsePro);
    _linCmtParsePro=0;
    parseFree(0);
    Rf_errorcall(R_NilValue, _("could not figure out linCmt() model"));
  }
  _linCmtParsePro=0;
  return lst;
}

SEXP _RxODE_linCmtGen(SEXP linCmt, SEXP vars, SEXP linCmtSens, SEXP verbose) {
  /* SEXP ret = PROTECT(allocVector(STRSXP, 1)); */
  int i = 0;
  linCmtGenStruct linG;
  linCmtGenIni(&linG);
  if (tb.hasKa){
    // depot, central
    for (i = 0; i < depotLines.n; i++){
      switch(depotLines.lType[i]){
      case FBIO:
	sClear(&(linG.d_F));
	sAppend(&(linG.d_F), "%s, ", depotLines.line[i]);
	break;
      case ALAG:
	sClear(&(linG.d_tlag));
	sAppend(&(linG.d_tlag), "%s, ", depotLines.line[i]);
	break;
      case RATE:
	sClear(&(linG.d_rate1));
	sAppend(&(linG.d_rate1), "%s, ", depotLines.line[i]);
	break;
      case DUR:
	sClear(&(linG.d_dur1));
	sAppend(&(linG.d_dur1), "%s, ", depotLines.line[i]);
	break;
      default:
	RSprintf("unknown depot line(%d): %s \n", depotLines.lType[i], depotLines.line[i]);
      }
    }
    for (i = 0; i < centralLines.n; i++){
      switch(centralLines.lType[i]){
      case FBIO:
	sClear(&(linG.d_F2));
	sAppend(&(linG.d_F2), "%s, ", centralLines.line[i]);
	break;
      case ALAG:
	sClear(&(linG.d_tlag2));
	sAppend(&(linG.d_tlag2), ", %s, ", centralLines.line[i]);
	break;
      case RATE:
	sClear(&(linG.d_rate2));
	sAppend(&(linG.d_rate2), "%s, ", centralLines.line[i]);
	break;
      case DUR:
	sClear(&(linG.d_dur2));
	sAppend(&(linG.d_dur2), "%s)", centralLines.line[i]);
	break;
      }
    }
  } else {
    for (i = 0; i < depotLines.n; i++){
      switch(depotLines.lType[i]){
      case FBIO:
	sAppendN(&(linG.last), "'f(depot)' ", 11);
	break;
      case ALAG:
	sAppendN(&(linG.last), "'alag(depot)' ", 14);
	break;
      case RATE:
	sAppend(&(linG.last), "'rate(depot)' ", 14);
	break;
      case DUR:
	sAppend(&(linG.last), "'dur(depot)' ", 13);
	break;
      default:
	RSprintf("unknown depot line(%d): %s \n", depotLines.lType[i], depotLines.line[i]);
      }
    }
    if (linG.last.o) {
      errLin[0] = '\0';
      errOff=0;
      snprintf(errLin, errLinLen, "%s does not exist without a 'depot' compartment, specify a 'ka' parameter", linG.last.s);
      errOff=strlen(errLin);
      sFree(&(linG.d_tlag));
      sFree(&(linG.d_tlag2));
      sFree(&(linG.d_F));
      sFree(&(linG.d_F2));
      sFree(&(linG.d_rate1));
      sFree(&(linG.d_dur1));
      sFree(&(linG.d_rate2));
      sFree(&(linG.d_dur2));
      sFree(&(linG.last));
      parseFree(0);
      Rf_errorcall(R_NilValue, errLin);
    }
    // central only
    for (i = 0; i < centralLines.n; i++){
      switch(centralLines.lType[i]){
      case FBIO:
	sClear(&(linG.d_F));
	sAppend(&(linG.d_F), "%s, ", centralLines.line[i]);
	break;
      case ALAG:
	sClear(&(linG.d_tlag));
	sAppend(&(linG.d_tlag), "%s, ", centralLines.line[i]);
	break;
      case RATE:
	sClear(&(linG.d_rate1));
	sAppend(&(linG.d_rate1), "%s, ", centralLines.line[i]);
	break;
      case DUR:
	sClear(&(linG.d_dur1));
	sAppend(&(linG.d_dur1), "%s, ", centralLines.line[i]);
	break;
      }
    }
  }
  int pro=0;
  SEXP inStr = PROTECT(Rf_allocVector(STRSXP, 4)); pro++;
  int doSens = 0;
  if (TYPEOF(linCmtSens) == INTSXP){
    doSens = INTEGER(linCmtSens)[0];
  }
  sAppend(&(linG.last), "%s%s%s%s", linG.d_tlag.s, linG.d_F.s, linG.d_rate1.s, linG.d_dur1.s);
  SET_STRING_ELT(inStr, 2, mkChar(linG.last.s));
  sClear(&(linG.last));
  sAppend(&(linG.last), "%s%s%s%s",linG.d_tlag2.s, linG.d_F2.s,  linG.d_rate2.s, linG.d_dur2.s);
  SET_STRING_ELT(inStr, 3, mkChar(linG.last.s));
  sClear(&(linG.last));
  if (doSens == 2){
    sAppend(&(linG.last), "linCmtB(rx__PTR__, t, %d, ", INTEGER(linCmt)[0]);
    SET_STRING_ELT(inStr, 0, mkChar(linG.last.s));
    SET_STRING_ELT(inStr, 1, mkChar("0, "));
  } else {
    if (doSens == 1){
      sAppend(&(linG.last), "linCmtA(rx__PTR__, t, %d, ", INTEGER(linCmt)[0]);
    } else if (doSens == 3) {
      sAppend(&(linG.last), "linCmtC(rx__PTR__, t, %d, ", INTEGER(linCmt)[0]);
    }
    SET_STRING_ELT(inStr, 0, mkChar(linG.last.s));
    SET_STRING_ELT(inStr, 1, mkChar(""));
  }
  _linCmtParsePro=pro;
  SEXP linCmtP = PROTECT(_linCmtParse(vars, inStr, verbose)); pro++;
  for (i = 0; i < sbNrmL.n; i++){
    if (sbNrmL.lProp[i]== -100){
      char *line = sbNrmL.line[i];
      if (line[0] != '\0') {
	while (strncmp(line, "linCmt(", 7)){
	  if (line[0] == '\0') {
	    linCmtGenFree(&linG);
	    UNPROTECT(pro);
	    parseFree(0);
	    Rf_errorcall(R_NilValue, _("linCmt() bad parse"));
	    return R_NilValue;
	  }
	  else sPut(&(linG.last2), line[0]);
	  line++;
	}
      }
      if (strlen(line) > 7) line +=7;
      else {
	linCmtGenFree(&linG);
	UNPROTECT(pro);
	parseFree(0);
	Rf_errorcall(R_NilValue, _("linCmt() bad parse"));
	return R_NilValue;
      }
      sAppend(&(linG.last2), "%s", CHAR(STRING_ELT(VECTOR_ELT(linCmtP, 0), 0)));
      while (line[0] != ')'){
	if (line[0] == '\0') {
	  linCmtGenFree(&linG);
	  UNPROTECT(pro);
	  parseFree(0);
	  Rf_errorcall(R_NilValue, _("linCmt() bad parse"));
	  return R_NilValue;
	}
	if (line[0] == '('){
	  linCmtGenFree(&linG);
	  UNPROTECT(pro);
	  parseFree(0);
	  Rf_errorcall(R_NilValue, _("linCmt() cannot have any extra parentheses in it"));
	  return R_NilValue;
	}
	line++;
      }
      if (line[0] != '\0') sAppend(&(linG.last2), "%s", ++line);
    } else {
      sAppend(&(linG.last2), "%s", sbNrmL.line[i]);
    }
  }
  SEXP ret = PROTECT(Rf_allocVector(STRSXP,1)); pro++;
  SET_STRING_ELT(ret, 0, mkChar(linG.last2.s));
  linCmtGenFree(&linG);
  UNPROTECT(pro);
  return ret;
}
