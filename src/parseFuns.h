////////////////////////////////////////////////////////////////////////////////
// RxODE parsing function routines

static inline int isAtFunctionArg(const char *name) {
  return !strcmp("(", name) ||
    !strcmp(")", name) ||
    !strcmp(",", name);
}

static inline void handleFunctionArguments(char *name, int depth) {
  if (isAtFunctionArg(name)) {
    sPut(&sb, name[0]);
    sPut(&sbDt, name[0]);
    if (!skipDouble && !(strcmp(",", name)) && depth == 1){
      aAppendN("(double)", 8);
      skipDouble=0;
    }
    sPut(&sbt, name[0]);
  }
}


static inline void setFunctionFlag(nodeInfo ni, char *name, int i, int *depth) {
  tb.fn = (i==0 && (nodeHas(function)) ? 1 : 0);
  if (tb.fn == 0) tb.fn = (i==0 && (nodeHas(function_name)) ? 2 : 0);
  if (tb.fn == 1) *depth = 0;
}


static inline int handleSimFunctions(nodeInfo ni, char *name, int *i, int nch,
				     D_ParseNode *pn){
  if (nodeHas(simfun_statement) && *i == 0) {
    *i = nch; // done
    if (tb.thread != 1) tb.thread = 2;
    sb.o=0;sbDt.o=0; sbt.o=0;
    D_ParseNode *xpn = d_get_child(pn, 0);
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    aType(TLOGIC);
    if (!strcmp("simeta", v)) {
      foundF0=1;
      if ((tb.simflg & 1) == 0) tb.simflg += 1;
    } else {
      if ((tb.simflg & 2) == 0) tb.simflg += 2;
    }
    sAppend(&sb, "%s(_cSub);\n  _SYNC_%s_;", v, v);
    sAppend(&sbDt, "%s(_cSub);\n  _SYNC_%s_;", v, v);
    sAppend(&sbt, "%s();", v);
    addLine(&sbPm, "%s\n", sb.s);
    addLine(&sbPmDt, "%s\n", sbDt.s);
    sAppend(&sbNrm, "%s\n", sbt.s);
    addLine(&sbNrmL, "%s\n", sbt.s);
    /* Free(v); */
    ENDLINE;
    return 1;
  }
  return 0;
}

typedef struct transFunctions {
  int isNorm;
  int isExp;
  int isF;
  int isGamma;
  int isBeta;
  int isPois;
  int isT;
  int isUnif;
  int isWeibull;
  int isNormV;
  int isCauchy;
  int isLead;
  int isFirst;
  int isLast;
  int isDiff;
  int isLinB;
  int isPnorm;
  int isTad;
  int isTafd;
  int isTlast;
  int isTfirst;
  int isInd;
  nodeInfo ni;
  char *name;
  int *i;
  int *depth;
  int nch;
  D_ParseNode *xpn;
  D_ParseNode *pn;
  char *v;
} transFunctions;

static inline void transFunctionsIni(transFunctions *tf) {
  tf->isNorm=0;
  tf->isExp=0;
  tf->isF=0;
  tf->isGamma=0;
  tf->isBeta=0;
  tf->isPois=0;
  tf->isT=0;
  tf->isUnif=0;
  tf->isWeibull=0;
  tf->isNormV=0;
  tf->isCauchy=0;
  tf->isLead=0;
  tf->isFirst=0;
  tf->isLast=0;
  tf->isDiff=0;
  tf->isLinB=0;
  tf->isPnorm=0;
  tf->isTad=0;
  tf->isTafd=0;
  tf->isTlast = 0;
  tf->isTfirst = 0;
  tf->isInd=0;
}

transFunctions _tf;

static inline int handleFunctionDosenum(transFunctions *tf) {
  if (!strcmp("dosenum", tf->v)) {
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    if (ii == 1){
      D_ParseNode *xpn = d_get_child(tf->pn, 2);
      char *v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
      if (allSpaces(v2)){
	aAppendN("(double)(_solveData->subjects[_cSub].dosenum)", 45);
	sAppendN(&sbt, "dosenum()", 9);
      } else {
	updateSyntaxCol();
	trans_syntax_error_report_fn(_("'dosenum' does not currently take arguments 'dosenum()'"));
      }
    } else {
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'dosenum' does not currently take arguments 'dosenum()'"));
    }
    tf->i[0] = tf->nch;
    return 1;
  }
  return 0;
}

static inline int handleFunctionTad(transFunctions *tf) {
  if ((tf->isTad = !strcmp("tad", tf->v)) || (tf->isTafd = !strcmp("tafd", tf->v)) ||
      (tf->isTlast = !strcmp("tlast", tf->v)) || (tf->isTfirst = !strcmp("tfirst", tf->v))) {
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    if (ii == 1){
      D_ParseNode *xpn = d_get_child(tf->pn, 2);
      char *v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
      if (allSpaces(v2)){
	// tad overall
	sAppend(&sb, "_%s0()", tf->v);
	sAppend(&sbDt, "_%s0()", tf->v);
      } else {
	sAppend(&sb, "_%s1(", tf->v);
	sAppend(&sbDt, "_%s1(", tf->v);
	if (new_de(v2)){
	  if (!strcmp("depot", v2)){
	    tb.hasDepot = 1;
	    aAppendN("_DEPOT_)", 8);
	  } else if (!strcmp("central", v2)){
	    tb.hasCentral = 1;
	    aAppendN("_CENTRAL_)", 10);
	  } else if (rx_syntax_require_ode_first){
	    updateSyntaxCol();
	    sPrint(&_gbuf,ODEFIRST,v2);
	    trans_syntax_error_report_fn(_gbuf.s);
	    /* Free(v2); */
	    /* Free(tf->v); */
	    return 1;
	  } else {
	    tb.statei++;
	    sAppend(&sb, "%d)", tb.de.n);
	    sAppend(&sbDt, "%d)", tb.de.n);
	  }
	} else {
	  new_or_ith(v2);
	  sAppend(&sb, "%d)", tb.id);
	  sAppend(&sbDt, "%d)", tb.id);
	}
	// tad(cmt)
      }
      sAppend(&sbt, "%s(%s)", tf->v, v2);
      tf->i[0] = tf->nch;
      return 1;
    }
  }
  return 0;
}

static inline int handleFunctionLogit(transFunctions *tf) {
  if (!strcmp("logit", tf->v) || !strcmp("expit", tf->v) ||
      !strcmp("invLogit", tf->v) || !strcmp("logitInv", tf->v)){
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    if (ii == 1){
      D_ParseNode *xpn = d_get_child(tf->pn, 2);
      char *v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
      if (allSpaces(v2)){
	updateSyntaxCol();
	sPrint(&_gbuf, _("'%s' takes 1-3 arguments '%s(x,low,high)'"),
	       tf->v, tf->v);
	/* Free(v2); */
	trans_syntax_error_report_fn(_gbuf.s);
      }
      /* Free(v2); */
      sAppend(&sb, "_%s1(", tf->v);
      sAppend(&sbDt,"_%s1(", tf->v);
      sAppend(&sbt, "%s(", tf->v);
    } else if (ii == 2) {
      sAppend(&sb, "_%s2(", tf->v);
      sAppend(&sbDt,"_%s2(", tf->v);
      sAppend(&sbt, "%s(", tf->v);
    } else if (ii == 3) {
      sAppend(&sb, "%s(", tf->v);
      sAppend(&sbDt,"%s(", tf->v);
      sAppend(&sbt, "%s(", tf->v);
    } else {
      updateSyntaxCol();
      sPrint(&_gbuf, _("'%s' takes 1-3 arguments '%s(x,low,high)'"),
	     tf->v, tf->v);
      trans_syntax_error_report_fn(_gbuf.s);
    }
    tf->i[0] = 1;// Parse next arguments
    tf->depth[0]=1;
    return 1;
  }
  return 0;
}

static inline int handleFunctionSum(transFunctions *tf) {
  if (!strcmp("prod",tf->v) || !strcmp("sum",tf->v) || !strcmp("sign",tf->v) ||
      !strcmp("max",tf->v) || !strcmp("min",tf->v)){
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    if (!strcmp("prod", tf->v)){
      sAppend(&sb, "_prod(_p, _input, _solveData->prodType, %d, (double) ", ii);
      sAppend(&sbDt, "_prod(_p, _input, _solveData->prodType, %d, (double) ", ii);
      if (maxSumProdN < ii){
	maxSumProdN = ii;
      }
    } else if (!strcmp("sum", tf->v)){
      sAppend(&sb, "_sum(_p, _pld, -__MAX_PROD__, _solveData->sumType, %d, (double) ", ii);
      sAppend(&sbDt, "_sum(_p, _pld, -__MAX_PROD__, _solveData->sumType, %d, (double) ", ii);
      if (SumProdLD < ii){
	SumProdLD = ii;
      }
    } else {
      sAppend(&sb, "_%s(%d, (double) ", tf->v, ii);
      sAppend(&sbDt, "_%s(%d, (double) ", tf->v, ii);
    }
    sAppend(&sbt, "%s(", tf->v);
    /* Free(tf->v); */
    tf->i[0] = 1;// Parse next arguments
    tf->depth[0]=1;
    return 1;
  }
  return 0;
}

static inline int isDiffFunction(transFunctions *tf) {
  return !strcmp("lag", tf->v) || (tf->isLead = !strcmp("lead", tf->v)) ||
    (tf->isDiff = !strcmp("diff", tf->v)) ||
    (tf->isFirst = !strcmp("first", tf->v)) ||
    (tf->isLast = !strcmp("last", tf->v));
}
static inline int getFunctionNargs(transFunctions *tf, int node){
  int ii = d_get_number_of_children(d_get_child(tf->pn,node))+1;
  if (ii == 1) {
    D_ParseNode *xpn = d_get_child(tf->pn, node-1);
    char *v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    if (allSpaces(v2)) {
      return 0;
    }
  }
  return ii;
}

static inline int assertCorrectDiffArgs(transFunctions *tf, int nargs, int *lagNo) {
  if (nargs != 1 && (tf->isFirst || tf->isLast)) {
    updateSyntaxCol();
    sPrint(&_gbuf, _("'%s' takes 1 argument %s(parameter), you have %d"),
	   tf->v, tf->v, nargs);
    trans_syntax_error_report_fn(_gbuf.s);
    return 1;
  }
  if (!(tf->isFirst || tf->isLast)) {
    if (!(nargs == 1 || nargs == 2)) {
      updateSyntaxCol();
      sPrint(&_gbuf, _("'%s' takes 1-2 arguments %s(parameter, k), you have %d"),
	     tf->v, tf->v, nargs);
      trans_syntax_error_report_fn(_gbuf.s);
      return 1;
    }
    // Check lag
    if (nargs == 2) {
      D_ParseNode *xpn = d_get_child(tf->pn, 3);
      char *v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
      if (strlen(v2) > 2){
	*lagNo = toInt(v2+1);
	if (tf->isLead && *lagNo != NA_INTEGER) *lagNo = -*lagNo;
      }
      if (*lagNo == NA_INTEGER){
	updateSyntaxCol();
	sPrint(&_gbuf, _("'%s(parameter, k)' requires k to be an integer"), tf->v);
	trans_syntax_error_report_fn(_gbuf.s);
	return 1;
      } else if (tf->isDiff && *lagNo <= 0){
	updateSyntaxCol();
	sPrint(&_gbuf, _("'%s(parameter, k)' requires k to be an integer >= 1"), tf->v);
	trans_syntax_error_report_fn(_gbuf.s);
	return 1;
      }
    }
  }
  return 0;
}

static inline int handleFunctionDiff(transFunctions *tf) {
  if (isDiffFunction(tf)) {
    int nargs = getFunctionNargs(tf, 3);
    // lag(par, 1) => lag_par(1)
    // lag(par) => lag_par(1)
    // Header what lag_par means.
    char *v2;
    int lagNo =0;
    D_ParseNode *xpn;
    if (assertCorrectDiffArgs(tf, nargs, &lagNo)) return 1;
    switch(nargs) {
    case 1:
      xpn = d_get_child(tf->pn, 2);
      v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
      lagNo=0;
      tb.fn=0;
      lagNo = 1;
      if (tf->isLead) lagNo=-1;
      if (tf->isFirst || tf->isLast) lagNo=NA_INTEGER;
      if (new_or_ith(v2)){
	addSymbolStr(v2);
	tb.lag[NV-1] = lagNo;
      } else {
	tb.lag[tb.ix] = lagNo;
      }
      tb.fn=1;
      sAppend(&sb,"%s_", tf->v);
      sAppend(&sbDt,"%s_", tf->v);
      doDot2(&sb, &sbDt, v2);
      sAppendN(&sb, "1(", 2);
      sAppendN(&sbDt, "1(", 2);
      sAppend(&sbt, "%s(", tf->v);
      break;
    case 2:
      // Check lag(x, 1);  Its OK with lhs, but nothing else is...
      xpn = d_get_child(tf->pn, 2);
      v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
      tb.fn=0;
      if (new_or_ith(v2)){
	addSymbolStr(v2);
	tb.lag[NV-1] = lagNo;
      } else {
	tb.lag[tb.ix] = lagNo;
      }
      tb.fn=1;
      if (lagNo == 0){
	doDot2(&sb, &sbDt, v2);
	/* Free(v2); */
	/* Free(tf->v); */
	tf->i[0] = 4;// skip next arguments
	tf->depth[0]=1;
	return 1;
      } else {
	skipDouble=1;
	sAppend(&sb,   "%s_", tf->v);
	sAppend(&sbDt,   "%s_", tf->v);
	doDot2(&sb, &sbDt, v2);
	sAppendN(&sb,   "(", 1);
	sAppendN(&sbDt,   "(", 1);
	sAppend(&sbt,  "%s(", tf->v);
      }
      break;
    default:

      break;
    }
    tf->i[0]     = 1;// Parse next arguments
    tf->depth[0] =1;
    return 1;
  }
  return 0;
}

static inline int handleFunctionPnorm(transFunctions *tf) {
  if ((tf->isPnorm = !strcmp("pnorm", tf->v)) ||
      !strcmp("qnorm", tf->v)){
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    if (ii == 1) {
      D_ParseNode *xpn = d_get_child(tf->pn, 2);
      char *v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
      int allSpace=allSpaces(v2);
      /* Free(v2); */
      if (allSpace){
	updateSyntaxCol();
	if (tf->isPnorm){
	  trans_syntax_error_report_fn(_("'pnorm' in RxODE takes 1-3 arguments pnorm(q, mean, sd)"));
	} else {
	  trans_syntax_error_report_fn(_("'qnorm' in RxODE takes 1-3 arguments pnorm(p, mean, sd)"));
	}
      } else {
	sAppend(&sb, "_%s1(", tf->v);
	sAppend(&sbDt,"_%s1(", tf->v);
	sAppend(&sbt, "%s(", tf->v);
      }
    } else if (ii == 2) {
      sAppend(&sb,"_%s2(", tf->v);
      sAppend(&sbDt,"_%s2(", tf->v);
      sAppend(&sbt, "%s(", tf->v);
    } else if (ii == 3) {
      sAppend(&sb,"_%s3(", tf->v);
      sAppend(&sbDt,"_%s3(", tf->v);
      sAppend(&sbt, "%s(", tf->v);
    } else {
      updateSyntaxCol();
      if (tf->isPnorm){
	trans_syntax_error_report_fn(_("'pnorm' in RxODE takes 1-3 arguments pnorm(q, mean, sd)"));
      } else {
	trans_syntax_error_report_fn(_("'qnorm' in RxODE takes 1-3 arguments pnorm(p, mean, sd)"));
      }
    }
    tf->i[0] = 1;// Parse next arguments
    tf->depth[0]=1;
    return 1;
  }
  return 0;
}

static inline int handleFunctionTransit(transFunctions *tf) {
  if (!strcmp("transit", tf->v)) {
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    if (ii == 2){
      aAppendN("_transit3P(t, _cSub, ", 21);
      sAppendN(&sbt,"transit(", 8);
      rx_podo=1;
    } else if (ii == 3){
      aAppendN("_transit4P(t, _cSub, ", 21);
      sAppendN(&sbt,"transit(", 8);
      rx_podo = 1;
    } else {
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'transit' takes 2-3 arguments transit(n, mtt, bio)"));
    }
    tf->i[0] = 1;// Parse next arguments
    tf->depth[0] = 1;
    return 1;
  }
  return 0;
}

static inline int isRxnormOrRelatedNode(transFunctions *tf) {
  return (tf->isNorm = !strcmp("rnorm", tf->v) ||
	  !strcmp("rxnorm", tf->v) || (tf->isInd = !strcmp("rinorm", tf->v))) ||
    (tf->isNormV = !strcmp("rnormV", tf->v) ||
     !strcmp("rxnormV", tf->v) || (tf->isInd = !strcmp("rinormV", tf->v))) ||
    (tf->isCauchy = !strcmp("rxcauchy", tf->v) || (tf->isInd = !strcmp("ricauchy", tf->v)) ||
     !strcmp("rcauchy", tf->v)) ||
    (tf->isF = !strcmp("rxf", tf->v) ||
     !strcmp("rf", tf->v) || (tf->isInd = !strcmp("rif", tf->v))) ||
    (tf->isGamma = !strcmp("rxgamma", tf->v) ||
     !strcmp("rgamma", tf->v) || (tf->isInd = !strcmp("rigamma", tf->v))) ||
    (tf->isBeta = !strcmp("rxbeta", tf->v) ||
     !strcmp("rbeta", tf->v) || (tf->isInd = !strcmp("ribeta", tf->v))) ||
    (tf->isUnif = !strcmp("rxunif", tf->v) ||
     !strcmp("runif", tf->v) || (tf->isInd = !strcmp("riunif", tf->v))) ||
    (tf->isWeibull = !strcmp("rxweibull", tf->v) ||
     !strcmp("rweibull", tf->v) || (tf->isInd = !strcmp("riweibull", tf->v)));
}

static inline int assertCorrectRxnormArgs2(transFunctions *tf, int nargs) {
  if (nargs != 2) {
    if (tf->isF) {
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'rif'/'rxf'/'rf' takes 2 arguments 'rxf(df1, df2)'"));
      return 1;
    }
    if (tf->isBeta) {
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'ribeta'/'rxbeta'/'rbeta' takes 2 arguments 'rxbeta(shape1, shape2)'"));
      return 1;
    }
  }
  return 0;
}

static inline int assertCorrectRxnormArgs12(transFunctions *tf, int nargs) {
  if (!(nargs == 1 || nargs == 2)) {
    if (tf->isGamma){
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'rigamma'/'rxgamma'/'rgamma' takes 1-2 arguments 'rxgamma(shape, rate)'"));
      return 1;
    }
    if (tf->isWeibull){
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'riweibull'/'rxweibull'/'rweibull' takes 1-2 arguments 'rxweibull(shape, scale)'"));
      return 1;
    }
  }
  return 0;
}

static inline int assertCorrectRxnormArgs02(transFunctions *tf, int nargs) {
  if (nargs > 2) {
    updateSyntaxCol();
    if (tf->isNormV) {
      trans_syntax_error_report_fn(_("'rinormV'/'rxnormV'/'rnormV' takes 0-2 arguments 'rxnormV(mean, sd)'"));
      return 1;
    }
    if (tf->isNorm){
      trans_syntax_error_report_fn(_("'rinorm'/'rxnorm'/'rnorm' takes 0-2 arguments 'rxnorm(mean, sd)'"));
      return 1;
    }
    if (tf->isUnif){
      trans_syntax_error_report_fn(_("'riunif'/'rxunif'/'runif' takes 0-2 arguments 'rxunif(min, max)'"));
      return 1;
    }
    if (tf->isCauchy) {
      trans_syntax_error_report_fn(_("'ricauchy'/'rxcauchy'/'rcauchy' takes 0-2 arguments 'rxcauchy(location, scale)'"));
      return 1;
    }
  }
  return 0;
}

static inline int assertCorrectRxnormArgs(transFunctions *tf, int nargs) {
  return assertCorrectRxnormArgs2(tf, nargs) ||
    assertCorrectRxnormArgs12(tf, nargs) ||
    assertCorrectRxnormArgs02(tf, nargs);
}

static inline int handleFunctionRxnorm(transFunctions *tf) {
  if (isRxnormOrRelatedNode(tf)) {
    if (tb.thread != 0) tb.thread = 2;
    int nargs = getFunctionNargs(tf, 3);
    if (assertCorrectRxnormArgs(tf, nargs)) return 1;
    switch (nargs) {
    case 0:
      if (tf->isInd) {
	// rxnorm()
	sAppend(&sb,"%s(&_solveData->subjects[_cSub], %d, 0.0, 1.0",  tf->v, tb.nInd);
	sAppend(&sbDt,"%s(&_solveData->subjects[_cSub], %d, 0.0, 1.0", tf->v, tb.nInd++);
	foundF0=1;
      } else {
	sAppend(&sb,"%s(&_solveData->subjects[_cSub], 0.0, 1.0", tf->v);
	sAppend(&sbDt,"%s(&_solveData->subjects[_cSub], 0.0, 1.0", tf->v);
      }
      sAppend(&sbt, "%s(", tf->v);
      break;
    case 1:
      if (tf->isInd) {
	sAppend(&sb,"%s1(%d,", tf->v, tb.nInd);
	sAppend(&sbDt,"%s1(%d,", tf->v, tb.nInd++);
	foundF0=1;
      } else {
	sAppend(&sb,"%s1(", tf->v);
	sAppend(&sbDt,"%s1(", tf->v);
      }
      sAppend(&sbt, "%s(", tf->v);
      break;
    case 2:
      if (tf->isInd){
	foundF0=1;
	sAppend(&sb,"%s(&_solveData->subjects[_cSub], %d, ", tf->v, tb.nInd);
	sAppend(&sbDt,"%s(&_solveData->subjects[_cSub], %d, ", tf->v, tb.nInd++);
      } else {
	sAppend(&sb,"%s(&_solveData->subjects[_cSub], ", tf->v);
	sAppend(&sbDt,"%s(&_solveData->subjects[_cSub], ", tf->v);
      }
      sAppend(&sbt, "%s(", tf->v);
    }
    tf->i[0] = 1;// Parse next arguments
    tf->depth[0]=1;
    return 1;
  }
  return 0;
}

static inline int handleFunctionRchisq(transFunctions *tf) {
  if (!strcmp("rchisq", tf->v) ||
      !strcmp("rxchisq", tf->v) ||
      (tf->isInd = !strcmp("richisq", tf->v)) ||
      (tf->isExp = !strcmp("rxexp", tf->v) ||
       !strcmp("rexp", tf->v) ||
       (tf->isInd = !strcmp("riexp", tf->v))) ||
      (tf->isT = !strcmp("rxt", tf->v) ||
       !strcmp("rt", tf->v) ||
       (tf->isInd = !strcmp("rit", tf->v)))) {
    if (tb.thread != 0) tb.thread = 2;
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    if (ii != 1){
      sPrint(&_gbuf, _("'%s' takes 1 arguments"), tf->v);
      updateSyntaxCol();
      trans_syntax_error_report_fn(_gbuf.s);
    } else {
      D_ParseNode *xpn = d_get_child(tf->pn, 2);
      char *v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
      int allSpace=allSpaces(v2);
      /* Free(v2); */
      if (allSpace){
	if (tf->isExp){
	  if (tf->isInd) {
	    sAppend(&sb,"%s(&_solveData->subjects[_cSub], %d, 1.0", tf->v, tb.nInd);
	    sAppend(&sbDt,"%s(&_solveData->subjects[_cSub], %d, 1.0", tf->v, tb.nInd++);
	    foundF0=1;
	  } else {
	    sAppend(&sb,"%s(&_solveData->subjects[_cSub], 1.0", tf->v);
	    sAppend(&sbDt,"%s(&_solveData->subjects[_cSub], 1.0", tf->v);
	  }
	  sAppend(&sbt, "%s(", tf->v);
	} else {
	  sPrint(&_gbuf, _("'%s' takes 1 argument"), tf->v);
	  updateSyntaxCol();
	  trans_syntax_error_report_fn(_gbuf.s);
	}
      } else if (tf->isT){
	if (tf->isInd) {
	  sAppend(&sb,"rit_(&_solveData->subjects[_cSub], %d, ", tb.nInd);
	  sAppend(&sbDt,"rit_(&_solveData->subjects[_cSub], %d, ", tb.nInd++);
	  foundF0=1;
	  sAppendN(&sbt, "rit(", 4);
	} else {
	  sAppendN(&sb,"rxt_(&_solveData->subjects[_cSub], ", 35);
	  sAppendN(&sbDt,"rxt_(&_solveData->subjects[_cSub], ", 35);
	  sAppendN(&sbt, "rxt(", 4);
	}
      } else {
	if (tf->isInd) {
	  sAppend(&sb,"%s(&_solveData->subjects[_cSub], %d, ", tf->v, tb.nInd);
	  sAppend(&sbDt,"%s(&_solveData->subjects[_cSub], %d, ", tf->v, tb.nInd++);
	  foundF0=1;
	} else {
	  sAppend(&sb,"%s(&_solveData->subjects[_cSub], ", tf->v);
	  sAppend(&sbDt,"%s(&_solveData->subjects[_cSub], ", tf->v);
	}
	sAppend(&sbt, "%s(", tf->v);
      }
    }
    tf->i[0] = 1;// Parse next arguments
    tf->depth[0]=1;
    return 1;
  }
  return 0;
}

static inline int isRgeomOrRelated(transFunctions *tf) {
  return (!strcmp("rxgeom", tf->v) ||
	  !strcmp("rgeom", tf->v) ||
	  (tf->isInd = !strcmp("rigeom", tf->v)) ||
	  (tf->isPois= !strcmp("rxpois", tf->v) ||
	   !strcmp("rpois", tf->v) ||
	   (tf->isInd = !strcmp("ripois", tf->v))));
}

static inline int assertCorrectGeomArgs(transFunctions *tf, int nargs) {
  if (nargs != 1) {
    updateSyntaxCol();
    if (tf->isPois){
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'ripois'/'rxpois'/'rpois' takes 1 argument 'rxpois(lambda)'"));
    } else {
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'rigeom'/'rxgeom'/'rgeom' takes 1 argument 'rxgeom(prob)'"));
    }
    return 1;
  }
  return 0;
}

static inline int handleFunctionRgeom(transFunctions *tf) {
  if (isRgeomOrRelated(tf)) {
    if (tb.thread != 0) tb.thread = 2;
    int nargs = getFunctionNargs(tf, 3);
    if (assertCorrectGeomArgs(tf, nargs)) return 1;
    if (tf->isInd) {
      sAppend(&sb, "(double)%s(&_solveData->subjects[_cSub], %d, ", tf->v, tb.nInd);
      sAppend(&sbDt, "(double)%s(&_solveData->subjects[_cSub], %d, ", tf->v, tb.nInd++);
      foundF0=1;
    } else {
      sAppend(&sb, "(double)%s(&_solveData->subjects[_cSub], ", tf->v);
      sAppend(&sbDt, "(double)%s(&_solveData->subjects[_cSub], ", tf->v);
    }
    sAppend(&sbt, "%s(", tf->v);
    tf->i[0] = 1;// Parse next arguments
    tf->depth[0]=1;
    return 1;
  }
  return 0;
}

static inline int handleFunctionRbinom(transFunctions *tf){
  if (!strcmp("rbinom", tf->v) ||
      !strcmp("rxbinom", tf->v) ||
      (tf->isInd = !strcmp("ribinom", tf->v))) {
    if (tb.thread != 0) tb.thread = 2;
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    if (ii != 2){
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'ribinom'/'rbinom'/'rxbinom' takes 2 arguments 'rxbinom(size, prob)'"));
    } else {
      if (tf->isInd){
	sAppend(&sb,   "(double)ribinom(&_solveData->subjects[_cSub], %d, (int)" , tb.nInd);
	sAppend(&sbDt, "(double)ribinom(&_solveData->subjects[_cSub], %d, (int)", tb.nInd++);
	sAppendN(&sbt, "ribinom(", 8);
      } else {
	aAppendN("(double)rxbinom(&_solveData->subjects[_cSub], (int)", 51);
	sAppendN(&sbt, "rxbinom(", 8);
      }
    }
    tf->i[0]     = 1;// Parse next arguments
    tf->depth[0] =1;
    return 1;
  }
  return 0;
}

static inline int handleFunctionIsNan(transFunctions *tf) {
  if (!strcmp("is.nan", tf->v)) {
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    D_ParseNode *xpn = d_get_child(tf->pn, 2);
    char *v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    int allSpace=allSpaces(v2);
    /* Free(v2); */
    if (ii != 1 || (ii == 1 && allSpace)) {
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'is.nan' takes 1 argument"));
    }
    sAppendN(&sb, "isnan", 5);
    sAppendN(&sbDt, "isnan", 5);
    sAppendN(&sbt, "is.nan", 6);
    return 1;
  }
  return 0;
}

static inline int handleFunctionIsNa(transFunctions *tf) {
  if (!strcmp("is.na", tf->v)) {
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    D_ParseNode *xpn = d_get_child(tf->pn, 2);
    char *v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    int allSpace=allSpaces(v2);
    /* Free(v2); */
    if (ii != 1 || (ii == 1 && allSpace)) {
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'is.na' takes 1 argument"));
    }
    sAppendN(&sb, "ISNA", 4);
    sAppendN(&sbDt, "ISNA", 4);
    sAppendN(&sbt, "is.na", 5);
    return 1;
  }
  return 0;
}

static inline int handleFunctionIsFinite(transFunctions *tf) {
  if (!strcmp("is.finite", tf->v)) {
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    D_ParseNode *xpn = d_get_child(tf->pn, 2);
    char *v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    int allSpace=allSpaces(v2);
    /* Free(v2); */
    if (ii != 1 || (ii == 1 && allSpace)) {
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'is.finite' takes 1 argument"));
    }
    sAppendN(&sb, "R_FINITE", 8);
    sAppendN(&sbDt, "R_FINITE", 8);
    sAppendN(&sbt, "is.finite", 9);
    return 1;
  }
  return 0;
}

static inline int handleFunctionIsInfinite(transFunctions *tf) {
  if (!strcmp("is.infinite", tf->v)) {
    int ii = d_get_number_of_children(d_get_child(tf->pn,3))+1;
    D_ParseNode *xpn = d_get_child(tf->pn, 2);
    char *v2 = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    int allSpace=allSpaces(v2);
    /* Free(v2); */
    if (ii != 1 || (ii == 1 && allSpace)) {
      updateSyntaxCol();
      trans_syntax_error_report_fn(_("'is.infinite' takes 1 argument"));
    }
    if (sbt.o > 0 && sbt.s[sbt.o-1] == '!'){
      sb.o--;sbDt.o--;
      sAppendN(&sb, "R_FINITE", 8);
      sAppendN(&sbDt, "R_FINITE", 8);
    } else {
      sAppendN(&sb, "!R_FINITE", 9);
      sAppendN(&sbDt, "!R_FINITE", 9);
    }
    sAppendN(&sbt, "is.infinite", 11);
    return 1;
  }
  return 0;
}

static inline void handleFunctionLinCmtAlag(transFunctions *tf, D_ParseNode *xpn1, D_ParseNode *xpn2) {
  // 10 tlag
  xpn2 = d_get_child(xpn1, 10+tf->isLinB);
  char* v2 = (char*)rc_dup_str(xpn2->start_loc.s+2, xpn2->end);
  if (!((!strcmp(v2, "0") || !strcmp(v2, "0.0") ||
	 !strcmp(v2, "0.") || !strcmp(v2, "")))) {
    // has interesting tlag
    int ixL = tb.ixL;
    int didEq = tb.didEq;
    if (foundLag == 0) needSort+=2; // & 2 when alag
    foundLag=1;
    aType(ALAG);
    addLine(&sbPm, "_alag[(&_solveData->subjects[_cSub])->linCmt] = %s;\n", v2);
    addLine(&sbPmDt, "_alag[(&_solveData->subjects[_cSub])->linCmt] = %s;\n", v2);
    addLine(&sbNrmL, "");
    ENDLINE;
    tb.ixL= ixL; tb.didEq=didEq;
  }
}

static inline void handleFunctionLinCmtF1(transFunctions *tf, D_ParseNode *xpn1, D_ParseNode *xpn2) {
  // 11 f1
  xpn2 = d_get_child(xpn1, 11+tf->isLinB);
  char* v2 = (char*)rc_dup_str(xpn2->start_loc.s+2, xpn2->end);
  if (!((!strcmp(v2, "1") || !strcmp(v2, "1.0") ||
	 !strcmp(v2, "1.") || !strcmp(v2, "")))) {
    // has interesting f1
    int ixL = tb.ixL;
    int didEq = tb.didEq;
    if (foundF == 0) needSort+=1;// & 1 when F
    foundF=1;
    aType(FBIO);
    addLine(&sbPm, "_f[(&_solveData->subjects[_cSub])->linCmt] = %s;\n", v2);
    addLine(&sbPmDt, "_f[(&_solveData->subjects[_cSub])->linCmt] = %s;\n", v2);
    addLine(&sbNrmL, "");
    /* sAppend(&sbNrm, "%s;\n", sbt.s); */
    ENDLINE;
    tb.ixL= ixL; tb.didEq=didEq;
  }
}

static inline void handleFunctionLinCmtDur1(transFunctions *tf, D_ParseNode *xpn1, D_ParseNode *xpn2) {
  // 13 dur1
  xpn2 = d_get_child(xpn1, 13+tf->isLinB);
  char* v2 = (char*)rc_dup_str(xpn2->start_loc.s + 2, xpn2->end);
  if (!((!strcmp(v2, "0") || !strcmp(v2, "0.0") ||
	 !strcmp(v2, "0.")) || !strcmp(v2, ""))) {
    // has interesting rate
    int ixL = tb.ixL;
    int didEq = tb.didEq;
    if (foundDur == 0) needSort+=4;// & 4 when dur
    foundDur=1;
    aType(DUR);
    addLine(&sbPm, "_dur[(&_solveData->subjects[_cSub])->linCmt] = %s;\n", v2);
    addLine(&sbPmDt, "_dur[(&_solveData->subjects[_cSub])->linCmt] = %s;\n", v2);
    addLine(&sbNrmL, "");
    /* sAppend(&sbNrm, "%s;\n", sbt.s); */
    ENDLINE;
    tb.ixL= ixL; tb.didEq=didEq;
  }
}

static inline void handleFunctionLinCmtRate1(transFunctions *tf, D_ParseNode *xpn1, D_ParseNode *xpn2) {
  // 12 rate1
  xpn2 = d_get_child(xpn1, 12+tf->isLinB);
  char* v2 = (char*)rc_dup_str(xpn2->start_loc.s+2, xpn2->end);
  if (!((!strcmp(v2, "0") || !strcmp(v2, "0.0") ||
	 !strcmp(v2, "0.") || !strcmp(v2, "")))) {
    // has interesting rate
    int ixL = tb.ixL;
    int didEq = tb.didEq;
    if (foundRate == 0) needSort+=8;// & 8 when rate
    foundRate=1;
    aType(RATE);
    addLine(&sbPm, "_rate[(&_solveData->subjects[_cSub])->linCmt] = %s;\n", v2);
    addLine(&sbPmDt, "_rate[(&_solveData->subjects[_cSub])->linCmt] = %s;\n", v2);
    addLine(&sbNrmL, "");
    /* sAppend(&sbNrm, "%s;\n", sbt.s); */
    ENDLINE;
    tb.ixL= ixL; tb.didEq=didEq;
  }
}

static inline void handleFunctionLinCmtKa(transFunctions *tf, D_ParseNode *xpn1, D_ParseNode *xpn2) {
  // 14 -- ka
  xpn2 = d_get_child(xpn1, 14+tf->isLinB);
  char* v2 = (char*)rc_dup_str(xpn2->start_loc.s+2, xpn2->end);
  if (!((!strcmp(v2, "0") || !strcmp(v2, "0.0") ||
	 !strcmp(v2, "0.")))) {
    tb.hasKa=1;
  }
  /* Free(v2); */
  // lag2 = 15
  xpn2 = d_get_child(xpn1, 15+tf->isLinB);
  v2 = (char*)rc_dup_str(xpn2->start_loc.s+2, xpn2->end);
  if (!((!strcmp(v2, "0") || !strcmp(v2, "0.0") ||
	 !strcmp(v2, "0.") || !strcmp(v2, "")))) {
    // has interesting tlag
    int ixL = tb.ixL;
    int didEq = tb.didEq;
    if (foundLag == 0) needSort+=2; // & 2 when alag
    foundLag=1;
    aType(ALAG);
    addLine(&sbPm, "_alag[(&_solveData->subjects[_cSub])->linCmt+1] = %s;\n", v2);
    addLine(&sbPmDt, "_alag[(&_solveData->subjects[_cSub])->linCmt+1] = %s;\n", v2);
    addLine(&sbNrmL, "");
    /* sAppend(&sbNrm, "%s;\n", sbt.s); */
    ENDLINE;
    tb.ixL= ixL; tb.didEq=didEq;
  }
}

static inline void handleFunctionLinCmtF2(transFunctions *tf, D_ParseNode *xpn1, D_ParseNode *xpn2) {
  // f2 = 16 ; This is 1 instead of zero
  xpn2 = d_get_child(xpn1, 16+tf->isLinB);
  char* v2 = (char*)rc_dup_str(xpn2->start_loc.s+2, xpn2->end);
  if (!((!strcmp(v2, "1") || !strcmp(v2, "1.0") ||
	 !strcmp(v2, "1.") || !strcmp(v2, "")))) {
    // has interesting f1
    int ixL = tb.ixL;
    int didEq = tb.didEq;
    if (foundF == 0) needSort+=1;// & 1 when F
    foundF=1;
    aType(FBIO);
    addLine(&sbPm, "_f[(&_solveData->subjects[_cSub])->linCmt+1] = %s;\n", v2);
    addLine(&sbPmDt, "_f[(&_solveData->subjects[_cSub])->linCmt+1] = %s;\n", v2);
    addLine(&sbNrmL, "");
    /* sAppend(&sbNrm, "%s;\n", sbt.s); */
    ENDLINE;
    tb.ixL= ixL; tb.didEq=didEq;
  }
}

static inline void handleFunctionLinCmtRate2(transFunctions *tf, D_ParseNode *xpn1, D_ParseNode *xpn2) {
  // rate2 = 17
  xpn2 = d_get_child(xpn1, 17+tf->isLinB);
  char* v2 = (char*)rc_dup_str(xpn2->start_loc.s+2, xpn2->end);
  if (!((!strcmp(v2, "0") || !strcmp(v2, "0.0") ||
	 !strcmp(v2, "0.") || !strcmp(v2, "")))) {
    // has interesting rate
    int ixL = tb.ixL;
    int didEq = tb.didEq;
    if (foundRate == 0) needSort+=8;// & 8 when rate
    foundRate=1;
    aType(RATE);
    addLine(&sbPm, "_rate[(&_solveData->subjects[_cSub])->linCmt+1] = %s;\n", v2);
    addLine(&sbPmDt, "_rate[(&_solveData->subjects[_cSub])->linCmt+1] = %s;\n", v2);
    addLine(&sbNrmL, "");
    /* sAppend(&sbNrm, "%s;\n", sbt.s); */
    ENDLINE;
    tb.ixL= ixL; tb.didEq=didEq;
  }
}

static inline void handleFunctionLinCmtDur2(transFunctions *tf, D_ParseNode *xpn1, D_ParseNode *xpn2) {
  xpn2 = d_get_child(xpn1, 18+tf->isLinB);
  char* v2 = (char*)rc_dup_str(xpn2->start_loc.s+2, xpn2->end);
  if (!((!strcmp(v2, "0") || !strcmp(v2, "0.0") ||
	 !strcmp(v2, "0.") || !strcmp(v2, "")))) {
    // has interesting rate
    int ixL = tb.ixL;
    int didEq = tb.didEq;
    if (foundDur == 0) needSort+=4;// & 4 when dur
    foundDur=1;
    aType(DUR);
    addLine(&sbPm, "_dur[(&_solveData->subjects[_cSub])->linCmt+1] = %s;\n", v2);
    addLine(&sbPmDt, "_dur[(&_solveData->subjects[_cSub])->linCmt+1] = %s;\n", v2);
    addLine(&sbNrmL, "");
    /* sAppend(&sbNrm, "%s;\n", sbt.s); */
    ENDLINE;
    tb.ixL= ixL; tb.didEq=didEq;
  }
}

static inline int handleFunctionLinCmt(transFunctions *tf) {
  if (!strcmp("linCmtA", tf->v) || !strcmp("linCmtC", tf->v) ||
      (tf->isLinB=!strcmp("linCmtB", tf->v))) {
    D_ParseNode *xpn1 = d_get_child(tf->pn, 3);
    D_ParseNode *xpn2 = d_get_child(xpn1, 1);
    char *v2 = (char*)rc_dup_str(xpn2->start_loc.s, xpn2->end);
    tb.linCmtN = toInt(v2+1);
    xpn2 = d_get_child(xpn1, 2);
    v2 = (char*)rc_dup_str(xpn2->start_loc.s, xpn2->end);
    tb.ncmt = toInt(v2+1);
    if (tf->isLinB) tf->isLinB=1;
    tb.linB = tf->isLinB;
    if (!tb.linExtra) {
      handleFunctionLinCmtAlag(tf, xpn1, xpn2);
      handleFunctionLinCmtF1(tf, xpn1, xpn2);
      handleFunctionLinCmtDur1(tf, xpn1, xpn2);
      handleFunctionLinCmtRate1(tf, xpn1, xpn2);
      handleFunctionLinCmtKa(tf, xpn1, xpn2);
      handleFunctionLinCmtF2(tf, xpn1, xpn2);
      handleFunctionLinCmtRate2(tf, xpn1, xpn2);
      handleFunctionLinCmtDur2(tf, xpn1, xpn2);
      tb.linExtra=true; // Only first call
    }
    aType(TLIN);
    if (tb.linB){
      xpn2 = d_get_child(xpn1, 4);
      v2 = (char*)rc_dup_str(xpn2->start_loc.s+2, xpn2->end);
      int tmp = toInt(v2);
      if (tmp > 0) {
	tmp--;
	tmp = 1 << tmp;
	if ((tb.linCmtFlg & tmp) == 0){
	  tb.linCmtFlg += tmp;
	}
      }
    }
    return 1;
  }
  return 0;
}

static inline int handleFunctionsExceptLinCmt(transFunctions *tf) {
  return handleFunctionDosenum(tf) ||
    handleFunctionTad(tf) ||
    handleFunctionSum(tf) ||
    handleFunctionLogit(tf) ||
    handleFunctionDiff(tf) ||
    handleFunctionPnorm(tf) ||
    handleFunctionTransit(tf) ||
    handleFunctionRxnorm(tf) ||
    handleFunctionRchisq(tf) ||
    handleFunctionRgeom(tf) ||
    handleFunctionRbinom(tf) ||
    handleFunctionIsNan(tf) ||
    handleFunctionIsNa(tf) ||
    handleFunctionIsFinite(tf) ||
    handleFunctionIsInfinite(tf);
}

static inline void handleBadFunctions(transFunctions *tf) {
  // Split out to handle anticipated automatic conversion of R
  // functions to C
  int foundFun = 0;
  for (int j = length(_goodFuns); j--;){
    if (!strcmp(CHAR(STRING_ELT(_goodFuns, j)),tf->v)){
      foundFun = 1;
      j=0;
      break;
    }
  }
  if (foundFun == 0){
    sPrint(&_gbuf, _("function '%s' is not supported in RxODE"), tf->v);
    updateSyntaxCol();
    trans_syntax_error_report_fn(_gbuf.s);
  }
}

static inline int handleFunctions(nodeInfo ni, char *name, int *i, int *depth, int nch, D_ParseNode *xpn, D_ParseNode *pn) {
  if (tb.fn == 1) {
    transFunctions *tf = &_tf;
    transFunctionsIni(tf);
    tf->ni = ni;
    tf->name = name;
    tf->i = i;
    tf->depth = depth;
    tf->nch = nch;
    tf->xpn = xpn;
    tf->pn = pn;
    tf->v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    if (handleFunctionsExceptLinCmt(tf)) {
      return 1;
    } else if (handleFunctionLinCmt(tf)){
      return 0;
    } else {
      handleBadFunctions(tf);
    }
  }
  return 0;
}

static inline int handlePrintf(nodeInfo ni, char *name, int i, D_ParseNode *xpn) {
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
      addLine(&sbNrmL, "%s;\n", sbt.s);
      ENDLINE
        }
    /* Free(v); */
    return 1;
  }
  return 0;
}
