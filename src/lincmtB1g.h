
#ifndef linCmtB1g_header
#define linCmtB1g_header
#include "dual.h"

static inline dualN oneCmtRateSSr1G(double *A, double *r1, dualN k10) {
dualN A1=divd2((*r1),k10);
A[0] = A1.f;
A[1] = A1.grad[dP1];
A[2] = A1.grad[dV1];
return A1;
}
static inline dualN oneCmtRateSSG(double *A, double *tinf, double *tau, double *r1, dualN k10) {
dualN eiK=expD(prod2d(negD(k10),(*tinf)));
dualN eK=div2(expD(prodd2(((*tau)-(*tinf)),negD(k10))),(subtrd2(1,expD(prod2d(negD(k10),(*tau))))));
dualN A1=div2(prod2(prodd2((*r1),(subtrd2(1,eiK))),eK),(k10));
A[0] = A1.f;
A[1] = A1.grad[dP1];
A[2] = A1.grad[dV1];
return A1;
}
static inline dualN oneCmtRateG(double *A, double *Alast, double *t, double *b1, double *r1, dualN k10) {
dualN A1last = iniD(Alast[0],-1);
A1last.grad[dKa] = Alast[1];
dualN eT=expD(prod2d(negD(k10),(*t)));
dualN A1=add2d(add2(prod2(divd2((*r1),k10),(subtrd2(1,eT))),prod2(A1last,eT)),(*b1));
A[0] = A1.f;
A[1] = A1.grad[dP1];
A[2] = A1.grad[dV1];
return A1;
}
static inline dualN oneCmtBolusSSG(double *A, double *tau, double *b1, dualN k10) {
dualN eT=divd2(1,(subtrd2(1,expD(prod2d(negD(k10),(*tau))))));
dualN A1=prodd2((*b1),eT);
A[0] = A1.f;
A[1] = A1.grad[dP1];
A[2] = A1.grad[dV1];
return A1;
}
static inline dualN oneCmtKaRateSSr1G(double *A, double *r1, dualN ka, dualN k20) {
// has ka
dualN A1=divd2((*r1),ka);
dualN A2=divd2((*r1),k20);
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
static inline dualN oneCmtKaRateSSr2G(double *A, double *r2, dualN ka, dualN k20) {
(void)(ka);
dualN A1=iniD(0, -1);
dualN A2=divd2((*r2),k20);
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
static inline dualN oneCmtKaRateSStr1G(double *A, double *tinf, double *tau, double *r1, dualN ka, dualN k20) {
// has ka
dualN eKa=div2(expD(prodd2(((*tau)-(*tinf)),negD(ka))),(subtrd2(1,expD(prodd2(-(*tau),ka)))));
dualN eiKa=expD(prod2d(negD(ka),(*tinf)));
dualN eiK=expD(prod2d(negD(k20),(*tinf)));
dualN eK=div2(expD(prodd2(((*tau)-(*tinf)),negD(k20))),(subtrd2(1,expD(prodd2(-(*tau),k20)))));
dualN A1=prod2(eKa,(subtr2(divd2((*r1),ka),div2(prod2d(eiKa,(*r1)),ka))));
dualN A2=add2(prod2(eK,(subtr2(add2(divd2((*r1),k20),div2(prod2d(eiKa,(*r1)),(add2(negD(k20),ka)))),div2(prod2(prod2d(eiK,(*r1)),ka),(subtr2(prod2(ka,k20),prod2(k20,k20))))))),div2(prod2(prod2(ka,(subtr2(eK,eKa))),(subtr2(divd2((*r1),ka),div2(prod2d(eiKa,(*r1)),ka)))),(add2(negD(k20),ka))));
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
static inline dualN oneCmtKaRateSStr2G(double *A, double *tinf, double *tau, double *r2, dualN ka, dualN k20) {
// has ka
dualN eiK=expD(prod2d(negD(k20),(*tinf)));
dualN eK=div2(expD(prodd2(((*tau)-(*tinf)),negD(k20))),(subtrd2(1,expD(prod2d(negD(k20),(*tau))))));
dualN A1=iniD(0, -1);
dualN A2=prod2(eK,(subtr2(divd2((*r2),k20),div2(prod2(prod2d(eiK,(*r2)),(add2(negD(k20),ka))),(subtr2(prod2(ka,k20),prod2(k20,k20)))))));
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
static inline dualN oneCmtKaRateG(double *A, double *Alast, double *t, double *b1, double *b2, double *r1, double *r2, dualN ka, dualN k20) {
// has ka
dualN A1last = iniD(Alast[0],-1);
A1last.grad[dKa] = Alast[2];
A1last.grad[dP1] = Alast[3];
A1last.grad[dV1] = Alast[4];
dualN A2last = iniD(Alast[1],-1);
A2last.grad[dKa] = Alast[5];
A2last.grad[dP1] = Alast[6];
A2last.grad[dV1] = Alast[7];
dualN eKa=expD(prod2d(negD(ka),(*t)));
dualN e20=expD(prod2d(negD(k20),(*t)));
dualN A1=add2d(subtr2(divd2((*r1),ka),div2((prod2((subtrd2((*r1),prod2(A1last,ka))),eKa)),ka)),(*b1));
dualN A2=add2d(add2(subtr2(div2((prod2((subtrd2((*r1),prod2(A1last,ka))),eKa)),(subtr2(ka,k20))),div2((prod2((add2(add2(add2(prod2d((subtr2(ka,k20)),(*r2)),prod2d(ka,(*r1))),prod2(prod2((subtr2(negD(A2last),A1last)),k20),ka)),prod2(prod2(A2last,k20),k20))),e20)),(subtr2(prod2(k20,ka),prod2(k20,k20))))),divd2(((*r2)+(*r1)),k20)),(*b2));
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
static inline dualN oneCmtKaSSb1G(double *A, double *tau, double *b1, dualN ka, dualN k20) {
// has ka
dualN eKa=divd2(1,(subtrd2(1,expD(prodd2(-(*tau),ka)))));
dualN eK=divd2(1,(subtrd2(1,expD(prodd2(-(*tau),k20)))));
dualN A1=prod2d(eKa,(*b1));
dualN A2=div2(prod2(prod2d(ka,(*b1)),(subtr2(eK,eKa))),(add2(negD(k20),ka)));
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
static inline dualN oneCmtKaSSb2G(double *A, double *tau, double *b2, dualN ka, dualN k20) {
(void)(ka);
dualN eK=divd2(1,(subtrd2(1,expD(prodd2(-(*tau),k20)))));
dualN A1=iniD(0, -1);
dualN A2=prod2d(eK,(*b2));
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
