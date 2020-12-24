#ifndef __PARSELINCMT_H__
#define __PARSELINCMT_H__
#pragma once
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>   /* dj: import intptr_t */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <unistd.h>
#include <errno.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif
#include "../inst/include/RxODE.h"
#include "tran.h"
#include "sbuf.h"

typedef struct linCmtStruct {
  int ka;

  int k12;
  int k21;

  int k13;
  int k31;

  int kel;

  int a;
  int b;
  int c;

  int aob;

  int alpha;
  int beta;
  int gamma;

  int cl;

  int cl1;
  int cl2;
  int cl3;
  int cl4;

  int v;

  int v1;
  int v2;
  int v3;
  int v4;

  int vp;
  int vp1;
  int vp2;
  int vp3;

  int vss;

  int cmtc;

  int clStyle;
  int vStyle;

  int trans;
  int ncmt;

  sbuf ret0;
  sbuf ret;
  const char *mid;
  SEXP vars;
} linCmtStruct;


static inline void linCmtIni(linCmtStruct *lin){
  lin->ka   = -1;

  lin->k12  = -1;
  lin->k21  = -1;

  lin->k13  = -1;
  lin->k31  = -1;

  lin->cmtc = -1;
  lin->kel  = -1;

  lin->cl = -1;

  lin->cl1 = -1;
  lin->cl2 = -1;
  lin->cl3 = -1;
  lin->cl4 = -1;

  lin->v = -1;

  lin->v1 = -1;
  lin->v2 = -1;
  lin->v3 = -1;
  lin->v4 = -1;

  lin->vp = -1;
  lin->vp1 = -1;
  lin->vp2 = -1;
  lin->vp3 = -1;

  lin->vss = -1;

  lin->a = -1;
  lin->b = -1;
  lin->c = -1;

  lin->aob = -1;

  lin->alpha = -1;
  lin->beta  = -1;
  lin->gamma = -1;

  lin->clStyle=-1;
  lin->vStyle = -1;
}

static inline void linCmtCmt(linCmtStruct *lin, const int cmt){
  if (lin->cmtc == -1) {
    lin->cmtc = cmt;
  }
  if (lin->cmtc != cmt){
    parseFree(0);
    Rf_errorcall(R_NilValue, _("inconsistent central compartment numbers, not sure if central compartment no is '1' or '2'"));
  }
}

#define errLinLen 150
extern char errLin[errLinLen];
extern int errOff;


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
#define linCmtVdStyle 1
#define linCmtVtStyle 2
#define linCmtVpStyle 3
#define linCmtVnStyle 4

static inline void linCmtVStr(const int style){
  switch(style){
  case linCmtVdStyle:
    snprintf(errLin + errOff, errLinLen-errOff, "Vd");
    errOff+=2;
    break;
  case linCmtVtStyle:
    snprintf(errLin + errOff, errLinLen-errOff, "Vt");
    errOff+=2;
    break;
  case linCmtVpStyle:
    snprintf(errLin + errOff, errLinLen-errOff, "Vp");
    errOff+=2;
    break;
  case linCmtVnStyle:
    snprintf(errLin + errOff, errLinLen-errOff, "V#");
    errOff+=2;
    break;
  }
}

static inline void linCmtVStyle(linCmtStruct *lin, int style) {
  if (lin->vStyle == -1) {
    lin->vStyle = style;
  }
  if (lin->vStyle != style) {
    errLin[0] = '\0';
    errOff = 0;
    snprintf(errLin + errOff, errLinLen-errOff, "cannot mix '");
    errOff += 12;
    linCmtVStr(lin->vStyle);
    snprintf(errLin + errOff, errLinLen-errOff, "' and '");
    errOff += 7;
    linCmtVStr(style);
    snprintf(errLin + errOff, errLinLen-errOff, "' volume styles");
    errOff += 15;
    parseFree(0);
    Rf_errorcall(R_NilValue, errLin);
  }
}

static inline void linCmtV(linCmtStruct *lin, const char *in, int *index) {
  if ((in[1] == 's' || in[1] == 'S') &&
      (in[2] == 's' || in[2] == 'S') &&
      in[3] == '\0') {
    lin->vss = *index;
    return;
  }
  // v1, v2, v3, v4
  // vt1, vt2, vt3, vt4
  // vp1, vp2, vp3, vp4
  if (in[1] == '\0') {
    lin->v = *index;
    return;
  }
  if (in[1] == '1' && in[2] == '\0') {
    linCmtVStyle(lin, 4);
    lin->v1 = *index;
    return;
  }
  if (in[1] == '2' && in[2] == '\0') {
    linCmtVStyle(lin, 4);
    lin->v2 = *index;
    return;
  }
  if (in[1] == '3' && in[2] == '\0') {
    linCmtVStyle(lin, 4);
    lin->v3 = *index;
    return;
  }
  if (in[1] == '4' && in[2] == '\0') {
    linCmtVStyle(lin, 4);
    lin->v4 = *index;
    return;
  }
  if ((in[1] == 'c' || in[1] == 'C') &&
      in[2] == '\0') {
    lin->v = *index;
    return;
  }
  int vType = (in[1] == 'd' || in[1] == 'D')*1 +
    (in[1] == 't' || in[1] == 'T')*2 +
    (in[1] == 'p' || in[1] == 'P')*3;
  if (vType) {
    linCmtVStyle(lin, vType);
    // vp, vp1
    // vp, vp2
    // vp, vp3
    // vp1, vp2
    // vp2, vp3
    if (in[2] == '\0') {
      lin->vp = *index;
      return;
    }
    if (in[2] == '1' && in[3] == '\0') {
      lin->vp1 = *index;
      return;
    }
    if (in[2] == '2' && in[3] == '\0') {
      lin->vp2 = *index;
      return;
    }
    if (in[2] == '3' && in[3] == '\0') {
      lin->vp3 = *index;
      return;
    }
  }
}

static inline void linCmtB(linCmtStruct *lin, const char *in, int *index) {
  if (in[1] == '\0') {
    lin->b = *index;
    return;
  }
  if ((in[1] == 'E' || in[1] == 'e') &&
      (in[2] == 'T' || in[2] == 't') &&
      (in[3] == 'A' || in[3] == 'a') &&
      in[4] == '\0') {
    lin->beta = *index;
    return;
  }
}

static inline void linCmtA(linCmtStruct *lin, const char *in, int *index) {
  if (in[1] == '\0') {
    lin->a = *index;
    return;
  }
  if ((in[1] == 'O' || in[1] == 'o') &&
      (in[2] == 'B' || in[2] == 'b') &&
      in[3] == '\0') {
    lin->aob = *index;
    return;
  }
  if ((in[1] == 'L' || in[1] == 'l') &&
      (in[2] == 'P' || in[2] == 'p') &&
      (in[3] == 'H' || in[3] == 'h') &&
      (in[4] == 'A' || in[4] == 'a') &&
      in[5] == '\0') {
    lin->alpha = *index;
    return;
  }
}

static inline void linCmtStr(linCmtStruct *lin, const char *in, int *index) {
  if (in[0] == 'v' || in[0] == 'V') {
    linCmtV(lin, in, index);
    return;
  }
  if (in[0] == 'c' || in[0] == 'C') {
    linCmtC(lin, in, index);
    return;
  }
  if (in[0] == 'k' || in[0] == 'K') {
    linCmtK(lin, in, index);
    return;
  }
  if (in[0] == 'Q' || in[0] == 'q') {
    linCmtQ(lin, in, index);
    return;
  }
  if (in[0] == 'A' || in[0] == 'a') {
    linCmtA(lin, in, index);
  }
  if (in[0] == 'B' || in[0] == 'b') {
    linCmtB(lin, in, index);
  }
  if ((in[0] == 'G' || in[0] == 'g') &&
      (in[1] == 'A' || in[1] == 'a') &&
      (in[2] == 'M' || in[2] == 'm') &&
      (in[3] == 'M' || in[3] == 'm') &&
      (in[4] == 'A' || in[4] == 'a') &&
      in[5] == '\0') {
    lin->gamma = *index;
    return;
  }
}

static inline void linCmtAdjustPars(linCmtStruct *lin) {
  if (lin->clStyle == linCmtQstyle || lin->clStyle == linCmtCld1style){
    // cl,
    if (lin->cl == -1){
      if (lin->clStyle == linCmtCld1style){
	parseFree(0);
	Rf_errorcall(R_NilValue, _("'Cld' parameterization needs 'Cl'"));
      } else {
	parseFree(0);
	Rf_errorcall(R_NilValue, _("'Q' parameterization needs 'Cl'"));
      }
    }
    if (lin->cl1 != -1) {
      if (lin->cl2  != -1) {
	// Cl, Q, Q1
	if (lin->clStyle == linCmtQstyle){
	  parseFree(0);
	  Rf_errorcall(R_NilValue, _("cannot mix 'Q' and 'Q1'"));
	} else {
	  parseFree(0);
	  Rf_errorcall(R_NilValue, _("cannot mix 'Cld' and 'Cld1'"));
	}
      } else if (lin->cl3 != -1) {
	// Cl, Q (cl1->cl2), Q2 (cl3->cl3)
	lin->cl2 = lin->cl1;
	lin->cl1 = -1;
      } else if (lin->cl4 != -1){
	// Cl, Q, Q3
	if (lin->clStyle == linCmtQstyle){
	  parseFree(0);
	  Rf_errorcall(R_NilValue, _("cannot mix 'Q' and 'Q3'"));
	} else {
	  parseFree(0);
	  Rf_errorcall(R_NilValue, _("cannot mix 'Cld' and 'Cld3'"));
	}
      } else {
	// Cl, Q (cl1->cl2), Q2 (cl3->cl3)
	lin->cl2 = lin->cl1;
	lin->cl1 = -1;
      }
    } else if (lin->cl2  != -1) {
      // Cl, Q1
      if (lin->cl4 != -1) {
	if (lin->clStyle == linCmtQstyle){
	  parseFree(0);
	  Rf_errorcall(R_NilValue, _("cannot mix 'Q1' and 'Q3'"));
	} else {
	  parseFree(0);
	  Rf_errorcall(R_NilValue, _("cannot mix 'Cld1' and 'Cld3'"));
	}
      }
    } else if (lin->cl3 != -1){
      lin->cl2 = lin->cl3;
      lin->cl3 = lin->cl4;
    }
  } else {
    if (lin->cl1 != -1){
      // Cl1, Cl2, Cl3
      // -> cl, cl2, cl3
      if (lin->cl != -1) {
	// cl, cl1,
	if (lin->cl2 == -1){
	  if (lin->cl4 != -1){
	    parseFree(0);
	    Rf_errorcall(R_NilValue, _("error parsing higher 'cl'"));
	  }
	  lin->cl4 = lin->cl3;
	  lin->cl3 = lin->cl2;
	  lin->cl2 = lin->cl1;
	  lin->cl1 = -1;
	} else {
	  parseFree(0);
	  Rf_errorcall(R_NilValue, _("cannot mix 'Cl' and 'Cl1'"));
	}
      } else {
	linCmtCmt(lin, 1);
	lin->cl = lin->cl1;
	lin->cl1 = -1;
	if (lin->cl4 != -1){
	  parseFree(0);
	  Rf_errorcall(R_NilValue, _("specified clearance for 4th compartment, which does not make sense in this context"));
	}
      }
    } else if (lin->cl2 != -1){
      if (lin->cl == -1){
	//  Cl2, Cl3, Cl4
	// -> Cl, cl2, cl3
	linCmtCmt(lin, 2);
      } else if (lin->cl4 != -1) {
	// Cl, Cl2, Cl3 keeps the same;  Cl4 doesn't make sense
	parseFree(0);
	Rf_errorcall(R_NilValue, _("specified clearance for 4th compartment, which does not make sense in this context"));
      }
    } else if (lin->cl != -1){
      if (lin->cl3 != -1){
	// Cl, Cl3, Cl4
	//-> Cl, Cl2, cl3
	lin->cl2 = lin->cl3;
	lin->cl3 = lin->cl4;
	lin->cl4 = -1;
      }
    }
  }
  if (lin->v != -1) {
    if (lin->v1 != -1){
      parseFree(0);
      Rf_errorcall(R_NilValue, _("Cannot specify 'v1' and 'vc'"));
    }
    if (lin->v4 != -1){
      parseFree(0);
      Rf_errorcall(R_NilValue, _("Cannot specify 'v4' and 'vc'"));
    }
    if (lin->v2 != -1) {
      // v, v2, v3; Central Compartment is 1
      linCmtCmt(lin, 1);
      linCmtVStyle(lin, 4); // V#
    } else if (lin->v3 != -1) {
      // v, v3, v4; Central compartment is 2
      linCmtCmt(lin, 2);
      linCmtVStyle(lin, 4); // V#
      lin->v2 = lin->v3;
      lin->v3 = lin->v4;
    } else if (lin->vp != -1){
      lin->v2 = lin->vp;
      if (lin->vp1 != -1){
	// v, vp, vp, vp1
	lin->v3 = lin->vp1;
      } else if (lin->vp2 != -1) {
	// v, vp, vp, vp2
	lin->v3 = lin->vp2;
      } else if (lin->vp3 != -1) {
	// v, vp, vp, vp3
	linCmtCmt(lin, 1);
	linCmtCmt(lin, 2);
      }
    } else if (lin->vp1 != -1) {
	// v, vp1, vp2
	lin->v2 = lin->vp1;
	lin->v3 = lin->vp2;
    } else if (lin->vp2 != -1) {
	lin->v2 = lin->vp2;
	lin->v3 = lin->vp3;
    }
  } else if (lin->v1 != -1) {
    linCmtCmt(lin, 1);
    lin->v = lin->v1;
    if (lin->v2 != -1) {
      // v1, v2, v3; Central Compartment is 1
      linCmtCmt(lin, 1);
      linCmtVStyle(lin, 4); // V#
    } else if (lin->v3 != -1) {
      // v1, v3, v4; Central compartment is 2
      linCmtCmt(lin, 2);
    } else if (lin->vp != -1){
      // v1, vp,
      lin->v2 = lin->vp;
      if (lin->vp1 != -1){
	// v1, vp, vp1
	lin->v3 = lin->vp1;
      } else if (lin->vp2 != -1) {
	// v, vp, vp2
	lin->v3 = lin->vp2;
      } else if (lin->vp3 != -1) {
	linCmtCmt(lin, 2);
      }
    }
  } else if (lin->v2 != -1){
    linCmtCmt(lin, 2);
    lin->v = lin->v2;
    lin->v2 = -1;
    if (lin->v3 != -1) {
      // v2, v3, v4; Central compartment is 2
      lin->v2 = lin->v3;
      lin->v3 = lin->v4;
    } else if (lin->vp != -1){
      // v2, vp,
      lin->v2 = lin->vp;
      if (lin->vp1 != -1){
	// v2, vp, vp1
	lin->v3 = lin->vp1;
      } else if (lin->vp2 != -1) {
	// v2, vp, vp2
	lin->v3 = lin->vp2;
      } else if (lin->vp3 != -1) {
	linCmtCmt(lin, 2);
      }
    }
  }
  if (lin->cl != -1 && lin->v != -1) {
    if (lin->cl2 != -1) {
      if (lin->v2 == -1 && lin->vss == -1) {
	parseFree(0);
	Rf_errorcall(R_NilValue, _("can find distributional clearance but not peripheral volume"));
      }
    }
    if (lin->v2 != -1) {
      if (lin->cl2 == -1) {
	parseFree(0);
	Rf_errorcall(R_NilValue, _("can find peripheral volume but not distributlin->v2 ional clearance"));
      }
    }
    if (lin->cl3 != -1) {
      if (lin->v3 == -1) {
	parseFree(0);
	Rf_errorcall(R_NilValue, _("can find 2nd distributional clearance but not 2nd peripheral volume"));
      }
    }
    if (lin->v3 != -1) {
      if (lin->cl3 == -1) {
	parseFree(0);
	Rf_errorcall(R_NilValue, _("can find 2nd peripheral volume but not 2nd distributional clearance"));
      }
    }
  }
  if (lin->v != -1 && lin->v2 != -1) {
    if (lin->v == lin->v2) {
      parseFree(0);
      Rf_errorcall(R_NilValue, _("cannot distinguish between central and peripheral volumes"));
    }
  }
  if (lin->v2 != -1 && lin->v3 != -1) {
    if (lin->v2 == lin->v3) {
      parseFree(0);
      Rf_errorcall(R_NilValue, _("cannot distinguish between 1st and 2nd peripheral volumes"));
    }
  }

  if (lin->cl != -1 && lin->cl2 != -1) {
    if (lin->cl == lin->cl2) {
      parseFree(0);
      Rf_errorcall(R_NilValue, _("cannot distinguish between central and peripheral clearances"));
    }
  }
  if (lin->cl2 != -1 && lin->cl3 != -1) {
    if (lin->cl2 == lin->cl3) {
      parseFree(0);
      Rf_errorcall(R_NilValue, _("cannot distinguish between 1st and 2nd distributional clearances"));
    }
  }
}

SEXP _linCmtParse(SEXP vars, SEXP inStr, SEXP verboseSXP);
SEXP _RxODE_linCmtGen(SEXP linCmt, SEXP vars, SEXP linCmtSens, SEXP verbose);

typedef struct linCmtGenStruct {
  sbuf last;
  sbuf d_tlag;
  sbuf d_tlag2;
  sbuf d_F;
  sbuf d_F2;
  sbuf d_rate1;
  sbuf d_dur1;
  sbuf d_rate2;
  sbuf d_dur2;
  sbuf last2;
} linCmtGenStruct;

static inline void linCmtGenIni(linCmtGenStruct *linG) {
  sIni(&(linG->last));
  sIni(&(linG->d_tlag));
  sIni(&(linG->d_tlag2));
  sIni(&(linG->d_F));
  sIni(&(linG->d_F2));
  sIni(&(linG->d_rate1));
  sIni(&(linG->d_dur1));
  sIni(&(linG->d_rate2));
  sIni(&(linG->d_dur2));
  sIni(&(linG->last2));

  sAppendN(&(linG->d_tlag),"0.0, ", 5);
  sAppendN(&(linG->d_tlag2), ", 0.0, ", 7);
  sAppendN(&(linG->d_F), "1.0, ", 5);
  sAppendN(&(linG->d_F2), "1.0, ", 5);
  sAppendN(&(linG->d_rate1), "0.0, ", 5);
  sAppendN(&(linG->d_dur1), "0.0, ", 5);
  sAppendN(&(linG->d_rate2), "0.0, ", 5);
  sAppendN(&(linG->d_dur2), "0.0)", 4);
}

static inline void linCmtGenFree(linCmtGenStruct *linG) {
  sFree(&(linG->last));
  sFree(&(linG->d_tlag));
  sFree(&(linG->d_tlag2));
  sFree(&(linG->d_F));
  sFree(&(linG->d_F2));
  sFree(&(linG->d_rate1));
  sFree(&(linG->d_dur1));
  sFree(&(linG->d_rate2));
  sFree(&(linG->d_dur2));
  sFree(&(linG->last2));
}

#endif // __PARSELINCMT_H__
