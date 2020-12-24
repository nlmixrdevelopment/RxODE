static inline void linCmtK(linCmtStruct *lin, const char *in, int *index) {
  // ke, kel, k10 or k20
  //
  if (in[1] == '\0') {
    lin->kel = *index;
    return;
  }
  if ((in[1] == 'e' || in[1] == 'E')) {
    if (in[2] == '\0') {
      if (lin->kel != -1) {
	parseFree(0);
	Rf_errorcall(R_NilValue, _("Ambiguous 'kel'"));
      }
      lin->kel = *index;
      return;
    }
    if ((in[2] == 'l' || in[2] == 'L') && in[3] == '\0') {
      if (lin->kel != -1) {
	parseFree(0);
	Rf_errorcall(R_NilValue, _("Ambiguous 'kel'"));
      }
      lin->kel = *index;
      return;
    }
  }
  if ((in[1] == 'a' || in[1] == 'A') &&
      in[2] == '\0') {
    lin->ka = *index;
    return;
  }
  // Support: k12, k21, k13, k31 when central compartment is 1
  // Also support:  k23, k32, k24, k42 when central compartment is 2
  if (in[1] == '1') {
    if (in[2] == '0' && in[3] == '\0') {
      linCmtCmt(lin, 1);
      if (lin->kel != -1) {
	parseFree(0);
	Rf_errorcall(R_NilValue, _("Ambiguous 'kel'"));
      }
      lin->kel = *index;
      return;
    }
    if (in[2] == '2' && in[3] == '\0') {
      linCmtCmt(lin, 1);
      lin->k12 = *index;
      return;
    }
    if (in[2] == '3' && in[3] == '\0') {
      linCmtCmt(lin, 1);
      lin->k13 = *index;
      return;
    }
  }
  if (in[1] == '2') {
    if (in[2] == '0' && in[3] == '\0') {
      linCmtCmt(lin, 2);
      if (lin->kel != -1) {
	parseFree(0);
	Rf_errorcall(R_NilValue, _("Ambiguous 'kel'"));
      }
      lin->kel = *index;
      return;
    }
    if (in[2] == '1' && in[3] == '\0') {
      linCmtCmt(lin, 1);
      lin->k21 = *index;
      return;
    }
    if (in[2] == '3' && in[3] == '\0') {
      linCmtCmt(lin, 2);
      // k23->k12
      lin->k12 = *index;
      return;
    }
    if (in[2] == '4' && in[3] == '\0') {
      linCmtCmt(lin, 2);
      // k24->k13
      lin->k13 = *index;
      return;
    }
  }
  if (in[1] == '3')  {
    if (in[2] == '1' && in[3] == '\0') {
      linCmtCmt(lin, 1);
      lin->k31 = *index;
      return;
    }
    if (in[2] == '2' && in[3] == '\0') {
      linCmtCmt(lin, 2);
      // k32->k21
      lin->k21 = *index;
      return;
    }
  }
  if (in[1] == '4' && in[2] == '2' && in[3] == '\0') {
    linCmtCmt(lin, 2);
    // k42 -> k31
    lin->k31 = *index;
  }
}

#define linCmtCl1style 1
#define linCmtCld1style 2
#define linCmtQstyle 3

static inline void linCmtClStr(const int style){
  switch(style){
  case linCmtCl1style:
    snprintf(errLin + errOff, errLinLen - errOff, "Cl#");
    errOff+=3;
    break;
  case linCmtCld1style:
    snprintf(errLin + errOff, errLinLen - errOff, "Cld#");
    errOff+=4;
    break;
  case linCmtQstyle:
    snprintf(errLin + errOff, errLinLen - errOff, "Q");
    errOff++;
    break;
  }
}

static inline void linCmtClStyle(linCmtStruct *lin, const int style) {
  if (lin->clStyle == -1) {
    lin->clStyle = style;
  }
  if (lin->clStyle != style) {
    errLin[0] = '\0';
    errOff=0;
    snprintf(errLin + errOff, errLinLen - errOff, "cannot mix '");
    errOff+=12;
    linCmtClStr(lin->clStyle);
    snprintf(errLin + errOff, errLinLen - errOff, "' and '");
    errOff+=7;
    linCmtClStr(style);
    snprintf(errLin + errOff, errLinLen - errOff, "' clearance styles");
    errOff+=18;
    parseFree(0);
    Rf_errorcall(R_NilValue, errLin);
  }
}

static inline void linCmtQ(linCmtStruct *lin, const char *in, int *index) {
  if (in[1] == '\0') {
    // Q
    linCmtClStyle(lin, linCmtQstyle);
    lin->cl1 = *index;
    return;
  }
  if (in[1] == '1' && in[2] == '\0') {
    // Q1
    linCmtClStyle(lin, linCmtQstyle);
    lin->cl2 = *index;
    return;
  }
  if (in[1] == '2' && in[2] == '\0') {
    // Q2
    linCmtClStyle(lin, linCmtQstyle);
    lin->cl3 = *index;
    return;
  }
  if (in[1] == '3' && in[2] == '\0') {
    // Q3
    linCmtClStyle(lin, linCmtQstyle);
    lin->cl4 = *index;
    return;
  }
}


static inline void linCmtC(linCmtStruct *lin, const char *in, int *index) {
  // CL, CL1, CL2, CL3, CL4
  // CLD, CLD1, CLD2, CLD3, CL4
  if (in[1] == '\0'){
    // C
    lin->c = *index;
    return;
  }
  if (in[1] == 'l' || in[1] == 'L') {
    if (in[2] == '\0') {
      lin->cl = *index;
      return;
    }
    if (in[2] == '1' && in[3] == '\0') {
      linCmtClStyle(lin, linCmtCl1style);
      lin->cl1 = *index;
      return;
    }
    if (in[2] == '2' && in[3] == '\0') {
      linCmtClStyle(lin, linCmtCl1style);
      lin->cl2 = *index;
      return;
    }
    if (in[2] == '3' && in[3] == '\0') {
      linCmtClStyle(lin, linCmtCl1style);
      lin->cl3 = *index;
      return;
    }
    if (in[2] == '4' && in[3] == '\0') {
      linCmtClStyle(lin, linCmtCl1style);
      lin->cl4 = *index;
      return;
    }
    if (in[2] == 'd' || in[2] == 'D') {
      if (in[3] == '\0') {
	linCmtClStyle(lin, linCmtCld1style);
	lin->cl1 = *index;
	return;
      }
      if (in[3] == '1' && in[4] == '\0') {
	linCmtClStyle(lin, linCmtCld1style);
	  lin->cl2 = *index;
	return;
      }
      if (in[3] == '2' && in[4] == '\0') {
	linCmtClStyle(lin, linCmtCld1style);
	lin->cl3 = *index;
	return;
      }
      if (in[3] == '3' && in[4] == '\0') {
	linCmtClStyle(lin, linCmtCld1style);
	lin->cl4 = *index;
	return;
      }
    }
  }
}
