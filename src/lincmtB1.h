
#ifndef linCmtB1_header
#define linCmtB1_header
#define safe_zero(a) ((a) == 0 ? DOUBLE_EPS : (a))
#define A1 A[0]
#define A2 A[1]
#define A1last Alast[0]
#define A2last Alast[1]

static inline void oneCmtRateSSr1D(double *A, double *r1, double *k10) {
#define A1k10 A[1]
A1 = (*r1)/safe_zero((*k10));
A1k10 = -(*r1)/safe_zero((((*k10)) * ((*k10))));
#undef A1k10
}

static inline void oneCmtRateSSD(double *A, double *tinf, double *tau, double *r1, double *k10) {
#define A1k10 A[1]
double rx0= (*tau) - (*tinf);
double rx1= exp(-(*tau) * (*k10));
double rx2= exp(-(*k10) * (*tinf));
double rx3= 1 - rx1;
double rx4= 1 - rx2;
double rx5= exp(-(*k10) * (rx0));
double rx6= (*k10) * (rx3);
double rx7= rx5 * (*r1);
A1 = rx7 * (rx4)/safe_zero((rx6));
A1k10 = -rx5 * (*r1) * (rx4)/safe_zero(((((*k10)) * ((*k10))) * (rx3))) + exp(-(*k10) * (*tinf) - (*k10) * (rx0)) * (*r1) * (*tinf)/safe_zero((rx6)) - rx7 * (rx0) * (rx4)/safe_zero((rx6)) - exp(-(*k10) * (rx0) - (*tau) * (*k10)) * (*r1) * (*tau) * (rx4)/safe_zero(((*k10) * ((rx3) * (rx3))));
#undef A1k10
}

static inline void oneCmtRateD(double *A, double *Alast, double *t, double *b1, double *r1, double *k10) {
#define A1k10 A[1]
#define A1lastk10 Alast[1]
double rx0= exp(-(*t) * (*k10));
double rx1= 1 - rx0;
double rx3= (*r1) * (rx1);
A1 = (*b1) + rx0 * A1last + rx3/safe_zero((*k10));
double rx2= (*t) * rx0;
A1k10 = rx0 * A1lastk10 - rx3/safe_zero((((*k10)) * ((*k10)))) - rx2 * A1last + rx2 * (*r1)/safe_zero((*k10));
#undef A1k10
#undef A1lastk10
}

static inline void oneCmtBolusSSD(double *A, double *tau, double *b1, double *k10) {
#define A1k10 A[1]
double rx0= exp(-(*tau) * (*k10));
double rx1= 1 - rx0;
A1 = (*b1)/safe_zero((rx1));
A1k10 = -rx0 * (*b1) * (*tau)/safe_zero(((rx1) * (rx1)));
#undef A1k10
}

static inline void oneCmtBolusD(double *A, double *Alast, double *t, double *b1, double *k10) {
#define A1k10 A[1]
#define A1lastk10 Alast[1]
double rx0= exp(-(*t) * (*k10));
A1 = (*b1) + rx0 * A1last;
A1k10 = rx0 * A1lastk10 - (*t) * rx0 * A1last;
#undef A1k10
#undef A1lastk10
}

static inline void oneCmtKaRateSSr1D(double *A, double *r1, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
A1 = (*r1)/safe_zero((*ka));
A1ka = -(*r1)/safe_zero((((*ka)) * ((*ka))));
A2 = (*r1)/safe_zero((*k20));
A2ka = 0;
A2k20 = -(*r1)/safe_zero((((*k20)) * ((*k20))));
#undef A1ka
#undef A2ka
#undef A2k20
}

static inline void oneCmtKaRateSSr2D(double *A, double *r2, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
A1 = 0;
A1ka = 0;
A2 = (*r2)/safe_zero((*k20));
A2ka = 0;
A2k20 = -(*r2)/safe_zero((((*k20)) * ((*k20))));
#undef A1ka
#undef A2ka
#undef A2k20
}

static inline void oneCmtKaRateSStr1D(double *A, double *tinf, double *tau, double *r1, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
double rx0= (*r1)/safe_zero((*ka));
double rx4= (*tau) - (*tinf);
double rx7= exp(-(*ka) * (*tau));
double rx8= exp(-(*ka) * (*tinf));
double rx12= 1 - rx7;
double rx15= rx8 * (*r1);
double rx25= rx0 - rx15/safe_zero((*ka));
A1 = exp(-(*ka) * (rx4)) * (rx25)/safe_zero((rx12));
double rx5= ((*ka)) * ((*ka));
double rx9= (*ka) * (rx4);
double rx22= rx15 * (*tinf);
double rx24= rx22/safe_zero((*ka));
double rx26= exp(-(*ka) * (*tau) - rx9);
double rx28= rx15/safe_zero((rx5));
double rx31= rx26 * (*tau);
A1ka = exp(-(*ka) * (rx4)) * (-(*r1)/(rx5) + rx28 + rx24)/safe_zero((rx12)) - exp(-(*ka) * (rx4)) * (rx4) * (rx25)/safe_zero((rx12)) - rx31 * (rx25)/safe_zero(((rx12) * (rx12)));
double rx1= (*r1)/safe_zero((*k20));
double rx2= (*ka) * (*k20);
double rx13= exp(-(*k20) * (*tinf));
double rx16= rx13 * (*r1);
double rx19= exp(-(*k20) * (rx4));
double rx20= rx16 * (*ka);
double rx27= rx15/safe_zero((-(*k20) + (*ka)));
double rx33= rx1 + rx27;
A2 = rx19 * (rx33 - rx20/(rx2 - (((*k20)) * ((*k20)))))/safe_zero((1 - exp(-(*tau) * (*k20)))) + (*ka) * (rx19/(1 - exp(-(*tau) * (*k20))) - exp(-(*ka) * (rx4))/(rx12)) * (rx25)/safe_zero((-(*k20) + (*ka)));
double rx6= ((*k20)) * ((*k20));
double rx23= (-(*k20) + (*ka)) * (-(*k20) + (*ka));
double rx38= rx15/safe_zero((rx23));
A2ka = (rx19/(1 - exp(-(*tau) * (*k20))) - exp(-(*ka) * (rx4))/(rx12)) * (rx25)/safe_zero((-(*k20) + (*ka))) + rx19 * (-rx13 * (*r1)/(rx2 - (rx6)) - rx38 - rx22/(-(*k20) + (*ka)) + rx20 * (*k20)/((rx2 - (((*k20)) * ((*k20)))) * (rx2 - (((*k20)) * ((*k20))))))/safe_zero((1 - exp(-(*tau) * (*k20)))) - (*ka) * (rx19/(1 - exp(-(*tau) * (*k20))) - exp(-(*ka) * (rx4))/(rx12)) * (rx25)/safe_zero((rx23)) + (*ka) * (rx19/(1 - exp(-(*tau) * (*k20))) - exp(-(*ka) *      (rx4))/(rx12)) * (-(*r1)/(rx5) + rx28 + rx24)/safe_zero((-(*k20) + (*ka))) + (*ka) * (exp(-(*ka) * (rx4)) * (rx4)/(rx12) + rx31/((rx12) * (rx12))) * (rx25)/safe_zero((-(*k20) + (*ka)));
double rx3= (*tau) * (*k20);
double rx29= exp(-(*k20) * (rx4) - rx3);
double rx32= rx29 * (*tau);
A2k20 = rx19 * (-(*r1)/(rx6) + rx38 + rx20 * (*tinf)/(rx2 - (rx6)) + rx20 * (-2 * (*k20) + (*ka))/((rx2 - (((*k20)) * ((*k20)))) * (rx2 - (((*k20)) * ((*k20))))))/safe_zero((1 - exp(-(*tau) * (*k20)))) + (*ka) * (rx19/(1 - exp(-(*tau) * (*k20))) - exp(-(*ka) * (rx4))/(rx12)) * (rx25)/safe_zero((rx23)) + (*ka) * (-rx19 * (rx4)/(1 - exp(-(*tau) * (*k20))) - rx32/((1 - exp(-(*tau) * (*k20))) * (1 - exp(-(*tau) * (*k20))))) * (rx25)/safe_zero((-(*k20) + (*ka))) -      rx19 * (rx4) * (rx33 - rx20/(rx2 - (rx6)))/safe_zero((1 - exp(-(*tau) * (*k20)))) - rx32 * (rx33 - rx20/(rx2 - (rx6)))/safe_zero(((1 - exp(-(*tau) * (*k20))) * (1 - exp(-(*tau) * (*k20)))));
#undef A1ka
#undef A2ka
#undef A2k20
}

static inline void oneCmtKaRateSStr2D(double *A, double *tinf, double *tau, double *r2, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
A1 = 0;
A1ka = 0;
double rx0= (*r2)/safe_zero((*k20));
double rx1= (*ka) * (*k20);
double rx2= (*tau) - (*tinf);
double rx5= exp(-(*tau) * (*k20));
double rx6= exp(-(*k20) * (*tinf));
double rx7= 1 - rx5;
double rx8= rx6 * (*r2);
double rx9= exp(-(*k20) * (rx2));
double rx11= rx8 * (-(*k20) + (*ka));
A2 = rx9 * (rx0 - rx11/(rx1 - (((*k20)) * ((*k20)))))/safe_zero((rx7));
double rx3= ((*k20)) * ((*k20));
A2ka = rx9 * (-rx6 * (*r2)/(rx1 - (rx3)) + rx8 * (*k20) * (-(*k20) + (*ka))/((rx1 - (((*k20)) * ((*k20)))) * (rx1 - (((*k20)) * ((*k20))))))/safe_zero((rx7));
A2k20 = rx9 * (-(*r2)/(rx3) + rx8/(rx1 - (rx3)) + rx8 * (*tinf) * (-(*k20) + (*ka))/(rx1 - (rx3)) + rx8 * (-2 * (*k20) + (*ka)) * (-(*k20) + (*ka))/((rx1 - (((*k20)) * ((*k20)))) * (rx1 - (((*k20)) * ((*k20))))))/safe_zero((rx7)) - rx9 * (rx2) * (rx0 - rx11/(rx1 - (rx3)))/safe_zero((rx7)) - exp(-(*k20) * (rx2) - (*tau) * (*k20)) * (*tau) * (rx0 - rx11/(rx1 - (rx3)))/safe_zero(((rx7) * (rx7)));
#undef A1ka
#undef A2ka
#undef A2k20
}

static inline void oneCmtKaRateD(double *A, double *Alast, double *t, double *b1, double *b2, double *r1, double *r2, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
#define A1lastka Alast[2]
#define A2lastka Alast[3]
#define A2lastk20 Alast[4]
double rx3= (*ka) * A1last;
double rx5= exp(-(*t) * (*ka));
double rx8= (*r1) - rx3;
double rx17= rx5 * (rx8);
A1 = (*b1) + (*r1)/safe_zero((*ka)) - rx17/safe_zero((*ka));
double rx4= ((*ka)) * ((*ka));
double rx9= (*t) * rx5;
double rx19= rx9 * (rx8);
double rx20= rx5 * (-(*ka) * A1lastka - A1last);
A1ka = -(*r1)/safe_zero((rx4)) + rx17/safe_zero((rx4)) - rx20/safe_zero((*ka)) + rx19/safe_zero((*ka));
double rx0= (*r1) + (*r2);
double rx1= (*r1) * (*ka);
double rx2= (*ka) * (*k20);
double rx6= exp(-(*t) * (*k20));
double rx10= (*r2) * (-(*k20) + (*ka));
double rx12= (-A1last - A2last) * (*ka);
A2 = (*b2) + (rx0)/safe_zero((*k20)) - rx6 * ((((*k20)) * ((*k20))) * A2last + rx1 + rx10 + rx12 * (*k20))/safe_zero((rx2 - (((*k20)) * ((*k20))))) + rx17/safe_zero((-(*k20) + (*ka)));
double rx7= ((*k20)) * ((*k20));
double rx16= (-(*k20) + (*ka)) * (-(*k20) + (*ka));
double rx23= rx17/safe_zero((rx16));
A2ka = -rx6 * (rx0 + (rx7) * A2lastka + (-A1last - A2last) * (*k20) + (-A1lastka - A2lastka) * (*ka) * (*k20))/safe_zero((rx2 - (rx7))) - rx23 + rx20/safe_zero((-(*k20) + (*ka))) - rx19/safe_zero((-(*k20) + (*ka))) + rx6 * (*k20) * ((rx7) * A2last + rx1 + rx10 + rx12 * (*k20))/safe_zero(((rx2 - (((*k20)) * ((*k20)))) * (rx2 - (((*k20)) * ((*k20))))));
A2k20 = -(rx0)/safe_zero((rx7)) - rx6 * (-(*r2) + 2 * (*k20) * A2last + (rx7) * A2lastk20 + rx12 + (-0 - A2lastk20) * (*ka) * (*k20))/safe_zero((rx2 - (rx7))) + rx23 + (*t) * rx6 * ((rx7) * A2last + rx1 + rx10 + rx12 * (*k20))/safe_zero((rx2 - (rx7))) + rx6 * (-2 * (*k20) + (*ka)) * ((rx7) * A2last + rx1 + rx10 + rx12 * (*k20))/safe_zero(((rx2 - (((*k20)) * ((*k20)))) * (rx2 - (((*k20)) * ((*k20)))))) - 0/safe_zero((-(*k20) +      (*ka)));
#undef A1ka
#undef A2ka
#undef A2k20
#undef A1lastka
#undef A2lastka
#undef A2lastk20
}

static inline void oneCmtKaSSb1D(double *A, double *tau, double *b1, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
double rx1= exp(-(*ka) * (*tau));
double rx3= 1 - rx1;
A1 = (*b1)/safe_zero((rx3));
A1ka = -rx1 * (*b1) * (*tau)/safe_zero(((rx3) * (rx3)));
double rx0= (*ka) * (*b1);
double rx2= exp(-(*tau) * (*k20));
double rx4= 1 - rx2;
A2 = rx0 * (1.0/(-(rx3)) + (1.0/(rx4)))/safe_zero((-(*k20) + (*ka)));
double rx5= 1/safe_zero((rx3));
double rx6= 1/safe_zero((rx4));
double rx7= (-(*k20) + (*ka)) * (-(*k20) + (*ka));
double rx8= (rx6);
double rx10= rx0 * (-(rx5) + rx8);
double rx11= rx10/safe_zero((rx7));
A2ka = (*b1) * (-(rx5) + rx8)/safe_zero((-(*k20) + (*ka))) - rx11 + rx1 * (*ka) * (*b1) * (*tau)/safe_zero(((-(*k20) + (*ka)) * ((rx3) * (rx3))));
A2k20 = rx11 - rx2 * (*ka) * (*b1) * (*tau)/safe_zero(((-(*k20) + (*ka)) * ((rx4) * (rx4))));
#undef A1ka
#undef A2ka
#undef A2k20
}

static inline void oneCmtKaSSb2D(double *A, double *tau, double *b2, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
A1 = 0;
A1ka = 0;
double rx0= exp(-(*tau) * (*k20));
double rx1= 1 - rx0;
A2 = (*b2)/safe_zero((rx1));
A2ka = 0;
A2k20 = -rx0 * (*b2) * (*tau)/safe_zero(((rx1) * (rx1)));
#undef A1ka
#undef A2ka
#undef A2k20
}

static inline void oneCmtKaD(double *A, double *Alast, double *t, double *b1, double *b2, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
#define A1lastka Alast[2]
#define A2lastka Alast[3]
#define A2lastk20 Alast[4]
double rx0= exp(-(*t) * (*ka));
A1 = (*b1) + rx0 * A1last;
double rx2= (*t) * rx0;
A1ka = rx0 * A1lastka - rx2 * A1last;
double rx1= exp(-(*t) * (*k20));
double rx4= rx1 - rx0;
double rx6= (*ka) * (rx4);
double rx7= rx6 * A1last;
A2 = (*b2) + rx1 * A2last + rx7/safe_zero((-(*k20) + (*ka)));
double rx5= (-(*k20) + (*ka)) * (-(*k20) + (*ka));
double rx8= rx7/safe_zero((rx5));
A2ka = rx1 * A2lastka + (rx4) * A1last/safe_zero((-(*k20) + (*ka))) - rx8 + rx6 * A1lastka/safe_zero((-(*k20) + (*ka))) + rx2 * (*ka) * A1last/safe_zero((-(*k20) + (*ka)));
double rx3= (*t) * rx1;
A2k20 = rx1 * A2lastk20 - rx3 * A2last + rx8 + 0/safe_zero((-(*k20) + (*ka))) - rx3 * (*ka) * A1last/safe_zero((-(*k20) + (*ka)));
#undef A1ka
#undef A2ka
#undef A2k20
#undef A1lastka
#undef A2lastka
#undef A2lastk20
}
 

#undef A1
#undef A2
#undef A1last
#undef A2last
#undef safe_zero
#endif
