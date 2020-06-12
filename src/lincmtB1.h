
#ifndef linCmtB1_header
#define linCmtB1_header
#define A1 A[0]
#define A2 A[1]
#define A1last Alast[0]
#define A2last Alast[1]

static inline void oneCmtRateSSr1D(double *A, double *r1, double *k10) {
#define A1k10 A[1]
A1=(*r1)/(*k10);
A1k10=-(*r1)/((*k10)*(*k10));
#undef A1k10
}

static inline void oneCmtRateSSD(double *A, double *tinf, double *tau, double *r1, double *k10) {
#define A1k10 A[1]
double rx0=(*tau)-(*tinf);
double rx1=exp(-(*tau)*(*k10));
double rx2=exp(-(*k10)*(*tinf));
double rx3=1-rx1;
double rx4=1-rx2;
double rx5=exp(-(*k10)*(rx0));
double rx6=(*k10)*(rx3);
double rx7=rx5*(*r1);
A1=rx7*(rx4)/(rx6);
A1k10=-rx5*(*r1)*(rx4)/((((*k10))*((*k10)))*(rx3))+exp(-(*k10)*(*tinf)-(*k10)*(rx0))*(*r1)*(*tinf)/(rx6)-rx7*(rx0)*(rx4)/(rx6)-exp(-(*k10)*(rx0)-(*tau)*(*k10))*(*r1)*(*tau)*(rx4)/((*k10)*(((rx3))*((rx3))));
#undef A1k10
}

static inline void oneCmtRateD(double *A, double *Alast, double *t, double *b1, double *r1, double *k10) {
#define A1k10 A[1]
#define A1lastk10 Alast[1]
double rx0=exp(-(*t)*(*k10));
double rx1=1-rx0;
double rx3=(*r1)*(rx1);
A1=(*b1)+rx0*A1last+rx3/(*k10);
double rx2=(*t)*rx0;
A1k10=rx0*A1lastk10-rx3/(((*k10))*((*k10)))-rx2*A1last+rx2*(*r1)/(*k10);
#undef A1k10
#undef A1lastk10
}

static inline void oneCmtBolusSSD(double *A, double *tau, double *b1, double *k10) {
#define A1k10 A[1]
double rx0=exp(-(*tau)*(*k10));
double rx1=1-rx0;
A1=1*(*b1)/(rx1);
A1k10=-1*rx0*(*b1)*(*tau)/(((rx1))*((rx1)));
#undef A1k10
}

static inline void oneCmtBolusD(double *A, double *Alast, double *t, double *b1, double *k10) {
#define A1k10 A[1]
#define A1lastk10 Alast[1]
double rx0=exp(-(*t)*(*k10));
A1=(*b1)+rx0*A1last;
A1k10=rx0*A1lastk10-(*t)*rx0*A1last;
#undef A1k10
#undef A1lastk10
}

static inline void oneCmtKaRateSSr1D(double *A, double *r1, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
A1=(*r1)/(*ka);
A1ka=-(*r1)/((*ka)*(*ka));
A2=(*r1)/(*k20);
A2ka=0;
A2k20=-(*r1)/((*k20)*(*k20));
#undef A1ka
#undef A2ka
#undef A2k20
}

static inline void oneCmtKaRateSSr2D(double *A, double *r2, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
A1=0;
A1ka=0;
A2=(*r2)/(*k20);
A2ka=0;
A2k20=-(*r2)/((*k20)*(*k20));
#undef A1ka
#undef A2ka
#undef A2k20
}

static inline void oneCmtKaRateSStr1D(double *A, double *tinf, double *tau, double *r1, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
double rx0=(*r1)/(*ka);
double rx4=(*tau)-(*tinf);
double rx6=exp(-(*ka)*(*tau));
double rx7=exp(-(*ka)*(*tinf));
double rx11=1-rx6;
double rx14=rx7*(*r1);
double rx23=rx0-rx14/(*ka);
A1=exp(-(*ka)*(rx4))*(rx23)/(rx11);
double rx5=(((*ka))*((*ka)));
double rx8=(*ka)*(rx4);
double rx21=rx14*(*tinf);
double rx22=rx21/(*ka);
double rx25=exp(-(*ka)*(*tau)-rx8);
double rx27=rx14/rx5;
double rx30=rx25*(*tau);
A1ka=exp(-(*ka)*(rx4))*(-(*r1)/rx5+rx27+rx22)/(rx11)-exp(-(*ka)*(rx4))*(rx4)*(rx23)/(rx11)-rx30*(rx23)/(((rx11))*((rx11)));
double rx1=(*r1)/(*k20);
double rx2=(*ka)*(*k20);
double rx9=(((*k20))*((*k20)));
double rx12=exp(-(*k20)*(*tinf));
double rx15=rx12*(*r1);
double rx18=exp(-(*k20)*(rx4));
double rx19=rx15*(*ka);
double rx20=rx2-rx9;
double rx26=rx14/(-(*k20)+(*ka));
double rx32=rx1+rx26;
double rx38=rx19/(rx20);
double rx41=rx32-rx38;
A2=rx18*(rx41)/(1-exp(-(*tau)*(*k20)))+(*ka)*(rx18/(1-exp(-(*tau)*(*k20)))-exp(-(*ka)*(rx4))/(rx11))*(rx23)/(-(*k20)+(*ka));
double rx24=(((-(*k20)+(*ka)))*((-(*k20)+(*ka))));
double rx37=rx14/rx24;
A2ka=(rx18/(1-exp(-(*tau)*(*k20)))-exp(-(*ka)*(rx4))/(rx11))*(rx23)/(-(*k20)+(*ka))+rx18*(-rx12*(*r1)/(rx20)-rx37-rx21/(-(*k20)+(*ka))+rx19*(*k20)/(((rx20))*((rx20))))/(1-exp(-(*tau)*(*k20)))-(*ka)*(rx18/(1-exp(-(*tau)*(*k20)))-exp(-(*ka)*(rx4))/(rx11))*(rx23)/rx24+(*ka)*(rx18/(1-exp(-(*tau)*(*k20)))-exp(-(*ka)*(rx4))/(rx11))*(-(*r1)/rx5+rx27+rx22)/(-(*k20)+(*ka))+(*ka)*(exp(-(*ka)*(rx4))*(rx4)/(rx11)+rx30/(((rx11))*((rx11))))*(rx23)/(-(*k20)+(*ka));
double rx3=(*tau)*(*k20);
double rx28=exp(-(*k20)*(rx4)-rx3);
double rx31=rx28*(*tau);
A2k20=rx18*(-(*r1)/rx9+rx37+rx19*(*tinf)/(rx20)+rx19*(-2*(*k20)+(*ka))/(((rx20))*((rx20))))/(1-exp(-(*tau)*(*k20)))+(*ka)*(rx18/(1-exp(-(*tau)*(*k20)))-exp(-(*ka)*(rx4))/(rx11))*(rx23)/rx24+(*ka)*(-rx18*(rx4)/(1-exp(-(*tau)*(*k20)))-rx31/(((1-exp(-(*tau)*(*k20))))*((1-exp(-(*tau)*(*k20))))))*(rx23)/(-(*k20)+(*ka))-rx18*(rx4)*(rx41)/(1-exp(-(*tau)*(*k20)))-rx31*(rx41)/(((1-exp(-(*tau)*(*k20))))*((1-exp(-(*tau)*(*k20)))));
#undef A1ka
#undef A2ka
#undef A2k20
}

static inline void oneCmtKaRateSStr2D(double *A, double *tinf, double *tau, double *r2, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
A1=0;
A1ka=0;
double rx0=(*r2)/(*k20);
double rx1=(*ka)*(*k20);
double rx2=(*tau)-(*tinf);
double rx3=(((*k20))*((*k20)));
double rx4=exp(-(*tau)*(*k20));
double rx5=exp(-(*k20)*(*tinf));
double rx6=1-rx4;
double rx7=rx5*(*r2);
double rx8=exp(-(*k20)*(rx2));
double rx9=rx1-rx3;
double rx10=rx7*(-(*k20)+(*ka));
double rx11=rx10/(rx9);
double rx13=rx0-rx11;
A2=rx8*(rx13)/(rx6);
A2ka=rx8*(-rx5*(*r2)/(rx9)+rx7*(*k20)*(-(*k20)+(*ka))/(((rx9))*((rx9))))/(rx6);
A2k20=rx8*(-(*r2)/rx3+rx7/(rx9)+rx7*(*tinf)*(-(*k20)+(*ka))/(rx9)+rx7*(-2*(*k20)+(*ka))*(-(*k20)+(*ka))/(((rx9))*((rx9))))/(rx6)-rx8*(rx2)*(rx13)/(rx6)-exp(-(*k20)*(rx2)-(*tau)*(*k20))*(*tau)*(rx13)/(((rx6))*((rx6)));
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
double rx3=(*ka)*A1last;
double rx4=exp(-(*t)*(*ka));
double rx7=(*r1)-rx3;
double rx15=rx4*(rx7);
A1=(*b1)+(*r1)/(*ka)-rx15/(*ka);
double rx5=(((*ka))*((*ka)));
double rx8=(*t)*rx4;
double rx18=rx8*(rx7);
double rx19=rx4*(-(*ka)*A1lastka-A1last);
A1ka=-(*r1)/rx5+rx15/rx5-rx19/(*ka)+rx18/(*ka);
double rx0=(*r1)+(*r2);
double rx1=(*r1)*(*ka);
double rx2=(*ka)*(*k20);
double rx6=exp(-(*t)*(*k20));
double rx9=(*r2)*(-(*k20)+(*ka));
double rx10=(((*k20))*((*k20)));
double rx11=(-A1last-A2last)*(*ka);
double rx12=rx10*A2last;
double rx13=rx2-rx10;
double rx17=rx12+rx1;
double rx20=rx17+rx9;
A2=(*b2)+(rx0)/(*k20)-rx6*(rx20+rx11*(*k20))/(rx13)+rx15/(-(*k20)+(*ka));
double rx16=(((-(*k20)+(*ka)))*((-(*k20)+(*ka))));
double rx21=rx15/rx16;
A2ka=-rx6*(rx0+rx10*A2lastka+(-A1last-A2last)*(*k20)+(-A1lastka-A2lastka)*(*ka)*(*k20))/(rx13)-rx21+rx19/(-(*k20)+(*ka))-rx18/(-(*k20)+(*ka))+rx6*(*k20)*(rx20+rx11*(*k20))/(((rx13))*((rx13)));
A2k20=-(rx0)/rx10-rx6*(-(*r2)+2*(*k20)*A2last+rx10*A2lastk20+rx11+(-0-A2lastk20)*(*ka)*(*k20))/(rx13)+rx21+(*t)*rx6*(rx20+rx11*(*k20))/(rx13)+rx6*(-2*(*k20)+(*ka))*(rx20+rx11*(*k20))/(((rx13))*((rx13)))-rx4*(*ka)*0/(-(*k20)+(*ka));
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
double rx1=exp(-(*ka)*(*tau));
double rx3=1-rx1;
A1=1*(*b1)/(rx3);
A1ka=-1*rx1*(*b1)*(*tau)/(((rx3))*((rx3)));
double rx0=(*ka)*(*b1);
double rx2=exp(-(*tau)*(*k20));
double rx4=1-rx2;
double rx5=(1.0/(rx3));
double rx6=(1.0/(rx4));
double rx7=1*rx6;
double rx10=rx0*(-1*rx5+rx7);
A2=rx10/(-(*k20)+(*ka));
double rx8=(((-(*k20)+(*ka)))*((-(*k20)+(*ka))));
double rx11=rx10/rx8;
A2ka=(*b1)*(-1*rx5+rx7)/(-(*k20)+(*ka))-rx11+1*rx1*(*ka)*(*b1)*(*tau)/((-(*k20)+(*ka))*(((rx3))*((rx3))));
A2k20=rx11-1*rx2*(*ka)*(*b1)*(*tau)/((-(*k20)+(*ka))*(((rx4))*((rx4))));
#undef A1ka
#undef A2ka
#undef A2k20
}

static inline void oneCmtKaSSb2D(double *A, double *tau, double *b2, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
A1=0;
A1ka=0;
double rx0=exp(-(*tau)*(*k20));
double rx1=1-rx0;
A2=1*(*b2)/(rx1);
A2ka=0;
A2k20=-1*rx0*(*b2)*(*tau)/(((rx1))*((rx1)));
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
double rx0=exp(-(*t)*(*ka));
A1=(*b1)+rx0*A1last;
double rx2=(*t)*rx0;
A1ka=rx0*A1lastka-rx2*A1last;
double rx1=exp(-(*t)*(*k20));
double rx4=rx1-rx0;
double rx6=(*ka)*(rx4);
double rx7=rx6*A1last;
A2=(*b2)+rx1*A2last+rx7/(-(*k20)+(*ka));
double rx5=(((-(*k20)+(*ka)))*((-(*k20)+(*ka))));
double rx8=rx7/rx5;
A2ka=rx1*A2lastka+(rx4)*A1last/(-(*k20)+(*ka))-rx8+rx6*A1lastka/(-(*k20)+(*ka))+rx2*(*ka)*A1last/(-(*k20)+(*ka));
double rx3=(*t)*rx1;
A2k20=rx1*A2lastk20-rx3*A2last+rx8+rx6*0/(-(*k20)+(*ka))-rx3*(*ka)*A1last/(-(*k20)+(*ka));
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
#endif
