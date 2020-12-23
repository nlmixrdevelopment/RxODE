#ifndef __TRAN_H__
#define __TRAN_H__
extern int rx_syntax_allow_dots;
void updateSyntaxCol();
void trans_syntax_error_report_fn(char *err);
void parseFree(int last);
void RSprintf(const char *format, ...);

/* char s_aux_info[64*MXSYM*4]; */
typedef struct symtab {
  vLines ss; // Symbol string or symbol lines
  /* char ss[64*MXSYM]; */                     /* symbol string: all vars*/
  vLines de;             /* symbol string: all Des*/
  int *lh;        /*
lhs symbols?
=0 not LHS
=1 LHS
=9 if a state var;
=10 if suppressed lhs;
=11 suppress parameter printout;
=19 is LHS with stateExtra
=29
=70 LHS + param
*/
  int *ini;        /* initial variable assignment =2 if there are two assignments */
  int *mtime;
  double *iniv;        /* Initial values */
  int *ini0;        /* state initial variable assignment =2 if there are two assignments */
  int *di;        /* ith of state vars */
  int *idi;       /* should ith state variable be ignored 0/1 */
  int *idu;       /* Has the ith state been used in a derivative expression? */
  int *lag;  // Lag number (if present)
  int *dvid;
  int dvidn;
  int ix;                       /* ith of curr symbol */
  int id;                       /* ith of curr symbol */
  int fn;                       /* curr symbol a fn?*/
  int ixL;// New assignment index
  int didEq;
  int NEnd;
  int pos_de;
  int ini_i; // #ini
  int statei; // # states
  int nExtra;
  int sensi;
  int li; // # lhs
  int sli; // # suppressed lhs
  int pi; // # param
  int isPi; // # pi?
  int isNA; // # pi?
  int linCmt; // Unparsed linear compartment
  int linCmtN; // Unparsed linear compartment
  int linCmtFlg; // Linear compartment flag
  // Save Jacobian information
  int *df;
  int *dy;
  int *sdfdy;
  int cdf;
  int ndfdy;
  int maxtheta;
  int hasCmt;
  int maxeta;
  int hasDepot;
  int hasCentral;
  int hasDepotCmt;
  int hasCentralCmt;
  int hasKa;
  int allocS;
  int allocD;
  int matn;
  int matnf;
  int ncmt;
  int linB;
  // curPropN
  int curPropN;
  int depotN;
  int centralN;
  // linCmt extras
  bool linExtra;
  int nwhile;
  int nInd;
  int simflg;
  int thread;
} symtab;

extern symtab tb;

extern vLines depotLines;
extern vLines centralLines;
extern vLines sbNrmL;


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
#define TMAT0 19
#define TMATF 20
#define TLIN 21


// new de type
#define fromDDT 2
#define fromCMTprop 1


#define NOASSIGN _("'<-' not supported, use '=' instead or set 'options(RxODE.syntax.assign = TRUE)'")
#define NEEDSEMI _("lines need to end with ';'\n     to match R's handling of line endings set 'options(RxODE.syntax.require.semicolon = FALSE)'")
#define NEEDPOW _("'**' not supported, use '^' instead or set 'options(RxODE.syntax.star.pow = TRUE)'")
#define NOINI0 _("'%s(0)' for initialization not allowed\n to allow set 'options(RxODE.syntax.allow.ini0 = TRUE)'")
#define NOSTATE _("defined 'df(%s)/dy(%s)', but '%s' is not a state")
#define NOSTATEVAR _("defined 'df(%s)/dy(%s)', but '%s' is not a state or variable")
#define ODEFIRST _("ODEs compartment 'd/dt(%s)' must be defined before changing/accessing its properties (f/alag/rate/dur/tad/tafd)\nIf you want to change this set 'options(RxODE.syntax.require.ode.first = FALSE).\nBe warned this may number compartments based on first occurance of property or ODE")
#define ZERODVID _("'dvid()' cannot have zeros in it")
#define ONEDVID _("RxODE only supports one dvid() statement per model")

#endif // __TRAN_H__
