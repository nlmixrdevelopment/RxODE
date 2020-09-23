
#ifndef linCmtB1g_header
#define linCmtB1g_header
#include "dual.h"

static inline dualN oneCmtRateSSr1G(double *A, double *r1, dualN *k10) {
dualN rx_dn0;
dualN A1;
assignD(&A1,divd2((*r1),k10,&rx_dn0));
A[0] = A1.f;
A[1] = A1.grad[dP1];
A[2] = A1.grad[dV1];
return A1;
}
static inline dualN oneCmtRateSSG(double *A, double *tinf, double *tau, double *r1, dualN *k10) {
dualN rx_dn0,rx_dn1;
dualN eiK;
assignD(&eiK,expDX(prod2dX(negD(k10,&rx_dn0),(*tinf))));
dualN eK;
assignD(&eK,div2X(expDX(prodd2X(((*tau)-(*tinf)),negD(k10,&rx_dn0))),(subtrd2X(1,expDX(prod2dX(negD(k10,&rx_dn1),(*tau)))))));
dualN A1;
assignD(&A1,div2X(prod2X(prodd2X((*r1),(subtrd2(1,&eiK,&rx_dn0))),&eK),(k10)));
A[0] = A1.f;
A[1] = A1.grad[dP1];
A[2] = A1.grad[dV1];
return A1;
}
static inline dualN oneCmtRateG(double *A, double *Alast, double *t, double *b1, double *r1, dualN *k10) {
dualN rx_dn0,rx_dn1,rx_dn2;
dualN A1last;
A1last.f = Alast[0];
A1last.grad[dKa] = Alast[1];
dualN eT;
assignD(&eT,expDX(prod2dX(negD(k10,&rx_dn0),(*t))));
dualN A1;
assignD(&A1,add2dX(add2X(prod2X(divd2((*r1),k10,&rx_dn0),(subtrd2(1,&eT,&rx_dn1))),prod2(&eT,&A1last,&rx_dn2)),(*b1)));
A[0] = A1.f;
A[1] = A1.grad[dP1];
A[2] = A1.grad[dV1];
return A1;
}
static inline dualN oneCmtBolusSSG(double *A, double *tau, double *b1, dualN *k10) {
dualN rx_dn0;
dualN eT;
assignD(&eT,divd2X(1,(subtrd2X(1,expDX(prod2dX(negD(k10,&rx_dn0),(*tau)))))));
dualN A1;
assignD(&A1,prodd2((*b1),&eT,&rx_dn0));
A[0] = A1.f;
A[1] = A1.grad[dP1];
A[2] = A1.grad[dV1];
return A1;
}
static inline dualN oneCmtKaRateSSr1G(double *A, double *r1, dualN *ka, dualN *k20) {
dualN rx_dn0;
// has ka
dualN A1;
assignD(&A1,divd2((*r1),ka,&rx_dn0));
dualN A2;
assignD(&A2,divd2((*r1),k20,&rx_dn0));
A[0] = A1.f;
A[1] = A2.f;
A[2] = A1.grad[dKa];
A[3] = A1.grad[dP1];
A[4] = A1.grad[dV1];
A[5] = A2.grad[dKa];
A[6] = A2.grad[dP1];
A[7] = A2.grad[dV1];
return A2;
}
static inline dualN oneCmtKaRateSSr2G(double *A, double *r2, dualN *ka, dualN *k20) {
dualN rx_dn0;
(void)(ka);
dualN A1;
iniD(0, -1, &A1);
dualN A2;
assignD(&A2,divd2((*r2),k20,&rx_dn0));
A[0] = A1.f;
A[1] = A2.f;
A[2] = A1.grad[dKa];
A[3] = A1.grad[dP1];
A[4] = A1.grad[dV1];
A[5] = A2.grad[dKa];
A[6] = A2.grad[dP1];
A[7] = A2.grad[dV1];
return A2;
}
static inline dualN oneCmtKaRateSStr1G(double *A, double *tinf, double *tau, double *r1, dualN *ka, dualN *k20) {
dualN rx_dn0,rx_dn1,rx_dn2,rx_dn3,rx_dn4,rx_dn5,rx_dn6,rx_dn7,rx_dn8,rx_dn9;
// has ka
dualN eKa;
assignD(&eKa,div2X(expDX(prodd2X(((*tau)-(*tinf)),negD(ka,&rx_dn0))),(subtrd2X(1,expDX(prodd2(-(*tau),ka,&rx_dn1))))));
dualN eiKa;
assignD(&eiKa,expDX(prod2dX(negD(ka,&rx_dn0),(*tinf))));
dualN eiK;
assignD(&eiK,expDX(prod2dX(negD(k20,&rx_dn0),(*tinf))));
dualN eK;
assignD(&eK,div2X(expDX(prodd2X(((*tau)-(*tinf)),negD(k20,&rx_dn0))),(subtrd2X(1,expDX(prodd2(-(*tau),k20,&rx_dn1))))));
dualN A1;
assignD(&A1,prod2X2(&eKa,(subtr2X(divd2((*r1),ka,&rx_dn0),div2X(prod2d(&eiKa,(*r1),&rx_dn1),ka)))));
dualN A2;
assignD(&A2,add2X(prod2X2(&eK,(subtr2X(add2X(divd2((*r1),k20,&rx_dn0),div2X(prod2d(&eiKa,(*r1),&rx_dn1),(add2X(negD(k20,&rx_dn2),ka)))),div2X(prod2X(prod2d(&eiK,(*r1),&rx_dn3),ka),(subtr2X(prod2(k20,ka,&rx_dn4),prod2(k20,k20,&rx_dn5))))))),div2X(prod2X(prod2X2(ka,(subtr2(&eK,&eKa,&rx_dn6))),(subtr2X(divd2((*r1),ka,&rx_dn7),div2X(prod2d(&eiKa,(*r1),&rx_dn8),ka)))),(add2X(negD(k20,&rx_dn9),ka)))));
A[0] = A1.f;
A[1] = A2.f;
A[2] = A1.grad[dKa];
A[3] = A1.grad[dP1];
A[4] = A1.grad[dV1];
A[5] = A2.grad[dKa];
A[6] = A2.grad[dP1];
A[7] = A2.grad[dV1];
return A2;
}
static inline dualN oneCmtKaRateSStr2G(double *A, double *tinf, double *tau, double *r2, dualN *ka, dualN *k20) {
dualN rx_dn0,rx_dn1,rx_dn2,rx_dn3,rx_dn4;
// has ka
dualN eiK;
assignD(&eiK,expDX(prod2dX(negD(k20,&rx_dn0),(*tinf))));
dualN eK;
assignD(&eK,div2X(expDX(prodd2X(((*tau)-(*tinf)),negD(k20,&rx_dn0))),(subtrd2X(1,expDX(prod2dX(negD(k20,&rx_dn1),(*tau)))))));
dualN A1;
iniD(0, -1, &A1);
dualN A2;
assignD(&A2,prod2X2(&eK,(subtr2X(divd2((*r2),k20,&rx_dn0),div2X(prod2X(prod2d(&eiK,(*r2),&rx_dn1),(add2X(negD(k20,&rx_dn2),ka))),(subtr2X(prod2(k20,ka,&rx_dn3),prod2(k20,k20,&rx_dn4))))))));
A[0] = A1.f;
A[1] = A2.f;
A[2] = A1.grad[dKa];
A[3] = A1.grad[dP1];
A[4] = A1.grad[dV1];
A[5] = A2.grad[dKa];
A[6] = A2.grad[dP1];
A[7] = A2.grad[dV1];
return A2;
}
static inline dualN oneCmtKaRateG(double *A, double *Alast, double *t, double *b1, double *b2, double *r1, double *r2, dualN *ka, dualN *k20) {
dualN rx_dn0,rx_dn1,rx_dn2,rx_dn3,rx_dn4,rx_dn5,rx_dn6,rx_dn7,rx_dn8;
// has ka
dualN A1last;
A1last.f = Alast[0];
A1last.grad[dKa] = Alast[2];
A1last.grad[dP1] = Alast[3];
A1last.grad[dV1] = Alast[4];
dualN A2last;
A2last.f = Alast[1];
A2last.grad[dKa] = Alast[5];
A2last.grad[dP1] = Alast[6];
A2last.grad[dV1] = Alast[7];
dualN eKa;
assignD(&eKa,expDX(prod2dX(negD(ka,&rx_dn0),(*t))));
dualN e20;
assignD(&e20,expDX(prod2dX(negD(k20,&rx_dn0),(*t))));
dualN A1;
assignD(&A1,add2dX(subtr2X(divd2((*r1),ka,&rx_dn0),div2X((prod2X((subtrd2X((*r1),prod2(ka,&A1last,&rx_dn1))),&eKa)),ka)),(*b1)));
dualN A2;
assignD(&A2,add2dX(add2X(subtr2X(div2X((prod2X((subtrd2X((*r1),prod2(ka,&A1last,&rx_dn0))),&eKa)),(subtr2(ka,k20,&rx_dn1))),div2X((prod2X((add2X(add2X(add2X(prod2dX((subtr2(ka,k20,&rx_dn2)),(*r2)),prod2d(ka,(*r1),&rx_dn3)),prod2X(prod2X((subtr2X(negD(&A2last,&rx_dn4),&A1last)),k20),ka)),prod2X(prod2(k20,&A2last,&rx_dn5),k20))),&e20)),(subtr2X(prod2(ka,k20,&rx_dn6),prod2(k20,k20,&rx_dn7))))),divd2(((*r2)+(*r1)),k20,&rx_dn8)),(*b2)));
A[0] = A1.f;
A[1] = A2.f;
A[2] = A1.grad[dKa];
A[3] = A1.grad[dP1];
A[4] = A1.grad[dV1];
A[5] = A2.grad[dKa];
A[6] = A2.grad[dP1];
A[7] = A2.grad[dV1];
return A2;
}
static inline dualN oneCmtKaSSb1G(double *A, double *tau, double *b1, dualN *ka, dualN *k20) {
dualN rx_dn0,rx_dn1,rx_dn2;
// has ka
dualN eKa;
assignD(&eKa,divd2X(1,(subtrd2X(1,expDX(prodd2(-(*tau),ka,&rx_dn0))))));
dualN eK;
assignD(&eK,divd2X(1,(subtrd2X(1,expDX(prodd2(-(*tau),k20,&rx_dn0))))));
dualN A1;
assignD(&A1,prod2d(&eKa,(*b1),&rx_dn0));
dualN A2;
assignD(&A2,div2X(prod2X(prod2d(ka,(*b1),&rx_dn0),(subtr2(&eK,&eKa,&rx_dn1))),(add2X(negD(k20,&rx_dn2),ka))));
A[0] = A1.f;
A[1] = A2.f;
A[2] = A1.grad[dKa];
A[3] = A1.grad[dP1];
A[4] = A1.grad[dV1];
A[5] = A2.grad[dKa];
A[6] = A2.grad[dP1];
A[7] = A2.grad[dV1];
return A2;
}
static inline dualN oneCmtKaSSb2G(double *A, double *tau, double *b2, dualN *ka, dualN *k20) {
dualN rx_dn0;
(void)(ka);
dualN eK;
assignD(&eK,divd2X(1,(subtrd2X(1,expDX(prodd2(-(*tau),k20,&rx_dn0))))));
dualN A1;
iniD(0, -1, &A1);
dualN A2;
assignD(&A2,prod2d(&eK,(*b2),&rx_dn0));
A[0] = A1.f;
A[1] = A2.f;
A[2] = A1.grad[dKa];
A[3] = A1.grad[dP1];
A[4] = A1.grad[dV1];
A[5] = A2.grad[dKa];
A[6] = A2.grad[dP1];
A[7] = A2.grad[dV1];
return A2;
}
#endif
