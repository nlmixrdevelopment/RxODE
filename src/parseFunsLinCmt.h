#include "needSortDefines.h"

static inline void handleFunctionLinCmtAlag(transFunctions *tf, D_ParseNode *xpn1, D_ParseNode *xpn2) {
  // 10 tlag
  xpn2 = d_get_child(xpn1, 10+tf->isLinB);
  char* v2 = (char*)rc_dup_str(xpn2->start_loc.s+2, xpn2->end);
  if (!((!strcmp(v2, "0") || !strcmp(v2, "0.0") ||
	 !strcmp(v2, "0.") || !strcmp(v2, "")))) {
    // has interesting tlag
    int ixL = tb.ixL;
    int didEq = tb.didEq;
    if (foundLag == 0) needSort+=needSortAlag; // & 2 when alag
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
    if (foundF == 0) needSort+=needSortF;// & 1 when F
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
    if (foundDur == 0) needSort+= needSortDur;// & 4 when dur
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
    if (foundRate == 0) needSort+= needSortRate;// & 8 when rate
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
    if (foundLag == 0) needSort+= needSortAlag; // & 2 when alag
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
    if (foundF == 0) needSort+= needSortF;// & 1 when F
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
    if (foundRate == 0) needSort+= needSortRate;// & 8 when rate
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
    if (foundDur == 0) needSort+= needSortDur;// & 4 when dur
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
