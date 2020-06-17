
#ifndef linCmtB1_header
#define linCmtB1_header
#define A1 A[0]
#define A2 A[1]
#define A1last Alast[0]
#define A2last Alast[1]

static inline void oneCmtRateSSr1D(double *A, double *r1, double *k10) {
#define A1k10 A[1]
A1=(*r1)/(*k10);
A1k10=-(*r1)/(((*k10))*((*k10)));
#undef A1k10
}

static inline void oneCmtRateSSD(double *A, double *tinf, double *tau, double *r1, double *k10) {
#define A1k10 A[1]
rx0=(*tau)-(*tinf);
rx1=exp(-(*tau)*(*k10));
rx2=exp(-(*k10)*(*tinf));
rx3=1-rx1;
rx4=1-rx2;
rx5=exp(-(*k10)*(rx0));
rx6=(*k10)*(rx3);
rx7=rx5*(*r1);
A1=rx7*(rx4)/(rx6);
A1k10=-rx5*(*r1)*(rx4)/((((*k10))*((*k10)))*(rx3))+exp(-(*k10)*(*tinf)-(*k10)*(rx0))*(*r1)*(*tinf)/(rx6)-rx7*(rx0)*(rx4)/(rx6)-exp(-(*k10)*(rx0)-(*tau)*(*k10))*(*r1)*(*tau)*(rx4)/((*k10)*((rx3)*(rx3)));
#undef A1k10
}

static inline void oneCmtRateD(double *A, double *Alast, double *t, double *b1, double *r1, double *k10) {
#define A1k10 A[1]
#define A1lastk10 Alast[1]
rx0=exp(-(*t)*(*k10));
rx1=1-rx0;
rx3=(*r1)*(rx1);
A1=(*b1)+rx0*A1last+rx3/(*k10);
rx2=(*t)*rx0;
A1k10=rx0*A1lastk10-rx3/(((*k10))*((*k10)))-rx2*A1last+rx2*(*r1)/(*k10);
#undef A1k10
#undef A1lastk10
}

static inline void oneCmtBolusSSD(double *A, double *tau, double *b1, double *k10) {
#define A1k10 A[1]
rx0=exp(-(*tau)*(*k10));
rx1=1-rx0;
A1=(*b1)/(rx1);
A1k10=-rx0*(*b1)*(*tau)/((rx1)*(rx1));
#undef A1k10
}

static inline void oneCmtBolusD(double *A, double *Alast, double *t, double *b1, double *k10) {
#define A1k10 A[1]
#define A1lastk10 Alast[1]
rx0=exp(-(*t)*(*k10));
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
A1ka=-(*r1)/(((*ka)) * ((*ka)));
A2=(*r1)/(*k20);
A2ka=0;
A2k20=-(*r1)/(((*k20)) * ((*k20)));
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
A2k20=-(*r2)/(((*k20))*((*k20)));
#undef A1ka
#undef A2ka
#undef A2k20
}

static inline void oneCmtKaRateSStr1D(double *A, double *tinf, double *tau, double *r1, double *ka, double *k20) {
#define A1ka A[2]
#define A2ka A[3]
#define A2k20 A[4]
rx0=(*r1)/(*ka);
rx4=(*tau)-(*tinf);
rx7=exp(-(*ka)*(*tau));
rx8=exp(-(*ka)*(*tinf));
rx12=1-rx7;
rx15=rx8*(*r1);
rx25=rx0-rx15/(*ka);
A1=exp(-(*ka)*(rx4))*(rx25)/(rx12);
rx5=((*ka))*((*ka));
rx9=(*ka)*(rx4);
rx22=rx15*(*tinf);
rx24=rx22/(*ka);
rx26=exp(-(*ka)*(*tau)-rx9);
rx28=rx15/(rx5);
rx31=rx26*(*tau);
A1ka=exp(-(*ka)*(rx4))*(-(*r1)/(rx5)+rx28+rx24)/(rx12)-exp(-(*ka)*(rx4))*(rx4)*(rx25)/(rx12)-rx31*(rx25)/((rx12)*(rx12));
rx1=(*r1)/(*k20);
rx2=(*ka)*(*k20);
rx13=exp(-(*k20)*(*tinf));
rx16=rx13*(*r1);
rx19=exp(-(*k20)*(rx4));
rx20=rx16*(*ka);
rx27=rx15/(-(*k20)+(*ka));
rx33=rx1+rx27;
A2=rx19*(rx33-rx20/(rx2-(((*k20))*((*k20)))))/(1-exp(-(*tau)*(*k20)))+(*ka)*(rx19/(1-exp(-(*tau)*(*k20)))-exp(-(*ka)*(rx4))/(rx12))*(rx25)/(-(*k20)+(*ka));
rx6=((*k20))*((*k20));
rx23=(-(*k20)+(*ka))*(-(*k20)+(*ka));
rx38=rx15/(rx23);
A2ka=(rx19/(1-exp(-(*tau)*(*k20)))-exp(-(*ka)*(rx4))/(rx12))*(rx25)/(-(*k20)+(*ka))+rx19*(-rx13*(*r1)/(rx2-(rx6))-rx38-rx22/(-(*k20)+(*ka))+rx20*(*k20)/((rx2-(((*k20))*((*k20))))*(rx2-(((*k20))*((*k20))))))/(1-exp(-(*tau)*(*k20)))-(*ka)*(rx19/(1-exp(-(*tau)*(*k20)))-exp(-(*ka)*(rx4))/(rx12))*(rx25)/(rx23)+(*ka)*(rx19/(1-exp(-(*tau)*(*k20)))-exp(-(*ka)*(rx4))/(rx12))*(-(*r1)/(rx5)+rx28+rx24)/(-(*k20)+(*ka))+(*ka)*(exp(-(*ka)*(rx4))*(rx4)/(rx12)+rx31/((rx12)*(rx12)))*(rx25)/(-(*k20)+(*ka));
rx3=(*tau)*(*k20);
rx29=exp(-(*k20)*(rx4)-rx3);
rx32=rx29*(*tau);
A2k20=rx19*(-(*r1)/(rx6)+rx38+rx20*(*tinf)/(rx2-(rx6))+rx20*(-2*(*k20)+(*ka))/((rx2-(((*k20))*((*k20))))*(rx2-(((*k20))*((*k20))))))/(1-exp(-(*tau)*(*k20)))+(*ka)*(rx19/(1-exp(-(*tau)*(*k20)))-exp(-(*ka)*(rx4))/(rx12))*(rx25)/(rx23)+(*ka)*(-rx19*(rx4)/(1-exp(-(*tau)*(*k20)))-rx32/((1-exp(-(*tau)*(*k20)))*(1-exp(-(*tau)*(*k20)))))*(rx25)/(-(*k20)+(*ka))-rx19*(rx4)*(rx33-rx20/(rx2-(rx6)))/(1-exp(-(*tau)*(*k20)))-rx32*(rx33-rx20/(rx2-(rx6)))/((1-exp(-(*tau)*(*k20)))*(1-exp(-(*tau)*(*k20))));
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
rx0=(*r2)/(*k20);
rx1=(*ka)*(*k20);
rx2=(*tau)-(*tinf);
rx5=exp(-(*tau)*(*k20));
rx6=exp(-(*k20)*(*tinf));
rx7=1-rx5;
rx8=rx6*(*r2);
rx9=exp(-(*k20)*(rx2));
rx11=rx8*(-(*k20)+(*ka));
A2=rx9*(rx0-rx11/(rx1-(((*k20))*((*k20)))))/(rx7);
rx3=((*k20))*((*k20));
A2ka=rx9*(-rx6*(*r2)/(rx1-(rx3))+rx8*(*k20)*(-(*k20)+(*ka))/((rx1-(((*k20))*((*k20))))*(rx1-(((*k20))*((*k20))))))/(rx7);
A2k20=rx9*(-(*r2)/(rx3)+rx8/(rx1-(rx3))+rx8*(*tinf)*(-(*k20)+(*ka))/(rx1-(rx3))+rx8*(-2*(*k20)+(*ka))*(-(*k20)+(*ka))/((rx1-(((*k20))*((*k20))))*(rx1-(((*k20))*((*k20))))))/(rx7)-rx9*(rx2)*(rx0-rx11/(rx1-(rx3)))/(rx7)-exp(-(*k20)*(rx2)-(*tau)*(*k20))*(*tau)*(rx0-rx11/(rx1-(rx3)))/((rx7)*(rx7));
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
rx3=(*ka)*A1last;
rx5=exp(-(*t)*(*ka));
rx8=(*r1)-rx3;
rx17=rx5*(rx8);
A1=(*b1)+(*r1)/(*ka)-rx17/(*ka);
rx4=((*ka))*((*ka));
rx9=(*t)*rx5;
rx19=rx9*(rx8);
rx20=rx5*(-(*ka)*A1lastka-A1last);
A1ka=-(*r1)/(rx4)+rx17/(rx4)-rx20/(*ka)+rx19/(*ka);
rx0=(*r1)+(*r2);
rx1=(*r1)*(*ka);
rx2=(*ka)*(*k20);
rx6=exp(-(*t)*(*k20));
rx10=(*r2)*(-(*k20)+(*ka));
rx12=(-A1last-A2last)*(*ka);
A2=(*b2)+(rx0)/(*k20)-rx6*((((*k20))*((*k20)))*A2last+rx1+rx10+rx12*(*k20))/(rx2-(((*k20))*((*k20))))+rx17/(-(*k20)+(*ka));
rx7=((*k20))*((*k20));
rx16=(-(*k20)+(*ka))*(-(*k20)+(*ka));
rx23=rx17/(rx16);
A2ka=-rx6*(rx0+(rx7)*A2lastka+(-A1last-A2last)*(*k20)+(-A1lastka-A2lastka)*(*ka)*(*k20))/(rx2-(rx7))-rx23+rx20/(-(*k20)+(*ka))-rx19/(-(*k20)+(*ka))+rx6*(*k20)*((rx7)*A2last+rx1+rx10+rx12*(*k20))/((rx2-(((*k20))*((*k20))))*(rx2-(((*k20))*((*k20)))));
A2k20=-(rx0)/(rx7)-rx6*(-(*r2)+2*(*k20)*A2last+(rx7)*A2lastk20+rx12+(-0-A2lastk20)*(*ka)*(*k20))/(rx2-(rx7))+rx23+(*t)*rx6*((rx7)*A2last+rx1+rx10+rx12*(*k20))/(rx2-(rx7))+rx6*(-2*(*k20)+(*ka))*((rx7)*A2last+rx1+rx10+rx12*(*k20))/((rx2-(((*k20))*((*k20))))*(rx2-(((*k20))*((*k20)))))-0/(-(*k20)+(*ka));
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
rx1=exp(-(*ka)*(*tau));
rx3=1-rx1;
A1=(*b1)/(rx3);
A1ka=-rx1*(*b1)*(*tau)/((rx3)*(rx3));
rx0=(*ka)*(*b1);
rx2=exp(-(*tau)*(*k20));
rx4=1-rx2;
A2=rx0*(1.0/(-(rx3))+(1.0/(rx4)))/(-(*k20)+(*ka));
rx5=1/(rx3);
rx6=1/(rx4);
rx7=(-(*k20)+(*ka))*(-(*k20)+(*ka));
rx8=(rx6);
rx10=rx0*(-(rx5)+rx8);
rx11=rx10/(rx7);
A2ka=(*b1)*(-(rx5)+rx8)/(-(*k20)+(*ka))-rx11+rx1*(*ka)*(*b1)*(*tau)/((-(*k20)+(*ka))*((rx3)*(rx3)));
A2k20=rx11-rx2*(*ka)*(*b1)*(*tau)/((-(*k20)+(*ka))*((rx4)*(rx4)));
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
rx0=exp(-(*tau)*(*k20));
rx1=1-rx0;
A2=(*b2)/(rx1);
A2ka=0;
A2k20=-rx0*(*b2)*(*tau)/((rx1)*(rx1));
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
rx0=exp(-(*t)*(*ka));
A1=(*b1)+rx0*A1last;
rx2=(*t)*rx0;
A1ka=rx0*A1lastka-rx2*A1last;
rx1=exp(-(*t)*(*k20));
rx4=rx1-rx0;
rx6=(*ka)*(rx4);
rx7=rx6*A1last;
A2=(*b2)+rx1*A2last+rx7/(-(*k20)+(*ka));
rx5=(-(*k20)+(*ka))*(-(*k20)+(*ka));
rx8=rx7/(rx5);
A2ka=rx1*A2lastka+(rx4)*A1last/(-(*k20)+(*ka))-rx8+rx6*A1lastka/(-(*k20)+(*ka))+rx2*(*ka)*A1last/(-(*k20)+(*ka));
rx3=(*t)*rx1;
A2k20=rx1*A2lastk20-rx3*A2last+rx8+0/(-(*k20)+(*ka))-rx3*(*ka)*A1last/(-(*k20)+(*ka));
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
