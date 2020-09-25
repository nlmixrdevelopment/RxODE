
#ifndef linCmtB2g_header
#define linCmtB2g_header
#include "dual.h"

static inline dualN twoCmtRateSSr1G(double *A, double *r1, dualN k10, dualN k12, dualN k21) {
dualN E1=add2(k10,k12);
dualN s=add2(E1,k21);
dualN sqr=sqrtD(subtr2(prod2(s,s),prodd2(4,(subtr2(prod2(E1,k21),prod2(k12,k21))))));
dualN lambda1=prodd2(0.5,(add2(s,sqr)));
dualN lambda2=prodd2(0.5,(subtr2(s,sqr)));
dualN l12=divd2(1,(prod2(lambda1,lambda2)));
dualN A1=prod2(prodd2((*r1),k21),l12);
dualN A2=prod2(prodd2((*r1),k12),l12);
double *_cur=A+2;
A[0]=A1.f;
A[1]=A2.f;
_cur += saveD(_cur, 0, 2, &A1);
saveD(_cur, 0, 2, &A2);
return A1;
}
static inline dualN twoCmtRateSSG(double *A, double *tinf, double *tau, double *r1, dualN k10, dualN k12, dualN k21) {
dualN E1=add2(k10,k12);
dualN E2=k21;
dualN s=add2(E1,E2);
dualN sqr=sqrtD(subtr2(prod2(s,s),prodd2(4,(subtr2(prod2(E1,E2),prod2(k12,k21))))));
dualN lambda1=prodd2(0.5,(add2(s,sqr)));
dualN lambda2=prodd2(0.5,(subtr2(s,sqr)));
dualN eTi1=expD(prodd2(-(*tinf),lambda1));
dualN eTi2=expD(prodd2(-(*tinf),lambda2));
dualN eT1=div2(expD(prodd2(((*tau)-(*tinf)),negD(lambda1))),(subtrd2(1,expD(prodd2(-(*tau),lambda1)))));
dualN eT2=div2(expD(prodd2(((*tau)-(*tinf)),negD(lambda2))),(subtrd2(1,expD(prodd2(-(*tau),lambda2)))));
dualN A1=div2((subtr2(prod2(eT1,(add2(subtr2(prod2(E2,(add2(div2((subtr2(prod2d(eTi1,(*r1)),prod2d(eTi2,(*r1)))),(add2(negD(lambda1),lambda2))),prod2(prodd2((*r1),E2),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2))))))))),prod2(lambda1,(add2(div2((subtr2(prod2d(eTi1,(*r1)),prod2d(eTi2,(*r1)))),(add2(negD(lambda1),lambda2))),prod2(prodd2((*r1),E2),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2)))))))))),prod2(prod2(prodd2((*r1),k12),k21),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2))))))))),prod2(eT2,(add2(subtr2(prod2(E2,(add2(div2((subtr2(prod2d(eTi1,(*r1)),prod2d(eTi2,(*r1)))),(add2(negD(lambda1),lambda2))),prod2(prodd2((*r1),E2),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2))))))))),prod2(lambda2,(add2(div2((subtr2(prod2d(eTi1,(*r1)),prod2d(eTi2,(*r1)))),(add2(negD(lambda1),lambda2))),prod2(prodd2((*r1),E2),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2)))))))))),prod2(prod2(prodd2((*r1),k12),k21),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2))))))))))),(add2(negD(lambda1),lambda2)));
dualN A2=div2((subtr2(prod2(eT1,(subtr2(add2(prod2(k12,(add2(div2((subtr2(prod2d(eTi1,(*r1)),prod2d(eTi2,(*r1)))),(add2(negD(lambda1),lambda2))),prod2(prodd2((*r1),E2),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2))))))))),prod2(prod2(prodd2((*r1),E1),k12),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2))))))),prod2(prod2(prodd2((*r1),k12),lambda1),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2))))))))),prod2(eT2,(subtr2(add2(prod2(k12,(add2(div2((subtr2(prod2d(eTi1,(*r1)),prod2d(eTi2,(*r1)))),(add2(negD(lambda1),lambda2))),prod2(prodd2((*r1),E2),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2))))))))),prod2(prod2(prodd2((*r1),E1),k12),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2))))))),prod2(prod2(prodd2((*r1),k12),lambda2),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eTi1,(prod2((subtr2(lambda1,lambda2)),lambda1)))),div2(eTi2,(prod2((subtr2(lambda1,lambda2)),lambda2))))))))))),(add2(negD(lambda1),lambda2)));
double *_cur=A+2;
A[0]=A1.f;
A[1]=A2.f;
_cur += saveD(_cur, 0, 2, &A1);
saveD(_cur, 0, 2, &A2);
return A1;
}
static inline dualN twoCmtRateG(double *A, double *Alast, double *t, double *b1, double *r1, dualN k10, dualN k12, dualN k21) {
dualN A1last;
dualN A2last;
double *_cur = Alast+2;
_cur += restoreD(Alast[0], _cur, 0, 2, &A1last);
restoreD(Alast[1], _cur, 0, 2, &A2last);
dualN E1=add2(k10,k12);
dualN E2=k21;
dualN s=add2(E1,E2);
dualN sqr=sqrtD(subtr2(prod2(s,s),prodd2(4,(subtr2(prod2(E1,E2),prod2(k12,k21))))));
dualN lambda1=prodd2(0.5,(add2(s,sqr)));
dualN lambda2=prodd2(0.5,(subtr2(s,sqr)));
dualN eT1=expD(prodd2(-(*t),lambda1));
dualN eT2=expD(prodd2(-(*t),lambda2));
dualN l12=(subtr2(lambda1,lambda2));
dualN l21=(subtr2(lambda2,lambda1));
dualN c10=(add2(add2d(prod2(A1last,E2),(*r1)),prod2(A2last,k21)));
dualN c11=div2((subtr2(c10,prod2(A1last,lambda1))),l21);
dualN c12=div2((subtr2(c10,prod2(A1last,lambda2))),l21);
dualN A1term1=subtr2(prod2(c11,eT1),prod2(c12,eT2));
dualN A1term2=prod2(prodd2((*r1),E2),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eT1,(prod2(lambda1,l12)))),div2(eT2,(prod2(lambda2,l12))))));
dualN A1=add2d(add2(A1term1,A1term2),(*b1));
dualN c20=(add2(prod2(A2last,E1),prod2(A1last,k12)));
dualN c21=div2((subtr2(c20,prod2(A2last,lambda1))),l21);
dualN c22=div2((subtr2(c20,prod2(A2last,lambda2))),l21);
dualN A2term1=subtr2(prod2(c21,eT1),prod2(c22,eT2));
dualN A2term2=prod2(prodd2((*r1),k12),(subtr2(add2(divd2(1,(prod2(lambda1,lambda2))),div2(eT1,(prod2(lambda1,l12)))),div2(eT2,(prod2(lambda2,(subtr2(lambda1,lambda2))))))));
dualN A2=add2(A2term1,A2term2);
_cur=A+2;
A[0]=A1.f;
A[1]=A2.f;
_cur += saveD(_cur, 0, 2, &A1);
saveD(_cur, 0, 2, &A2);
return A1;
}
static inline dualN twoCmtBolusSSG(double *A, double *tau, double *b1, dualN k10, dualN k12, dualN k21) {
dualN E2=k21;
dualN s=add2(add2(k12,k21),k10);
dualN sqr=sqrtD(subtr2(prod2(s,s),prod2(prodd2(4,k21),k10)));
dualN lambda1=prodd2(0.5,(add2(s,sqr)));
dualN lambda2=prodd2(0.5,(subtr2(s,sqr)));
dualN eL1=divd2(1,(subtrd2(1,expD(prodd2(-(*tau),lambda1)))));
dualN eL2=divd2(1,(subtrd2(1,expD(prodd2(-(*tau),lambda2)))));
dualN A1=div2((subtr2(prod2(eL1,(subtr2(prodd2((*b1),E2),prodd2((*b1),lambda1)))),prod2(eL2,(subtr2(prodd2((*b1),E2),prodd2((*b1),lambda2)))))),(add2(negD(lambda1),lambda2)));
dualN A2=div2((subtr2(prod2(prod2d(eL1,(*b1)),k12),prod2(prod2d(eL2,(*b1)),k12))),(add2(negD(lambda1),lambda2)));
double *_cur=A+2;
A[0]=A1.f;
A[1]=A2.f;
_cur += saveD(_cur, 0, 2, &A1);
saveD(_cur, 0, 2, &A2);
return A1;
}
static inline dualN twoCmtKaRateSSr1G(double *A, double *r1, dualN ka, dualN k20, dualN k23, dualN k32) {
//has ka
dualN s=add2(add2(k23,k32),k20);
dualN beta=prodd2(0.5,(subtr2(s,sqrtD(subtr2(prod2(s,s),prod2(prodd2(4,k32),k20))))));
dualN alpha=div2(prod2(k32,k20),beta);
dualN A1=divd2((*r1),ka);
dualN A2=div2(prodd2((*r1),k32),(prod2(beta,alpha)));
dualN A3=div2(prodd2((*r1),k23),(prod2(beta,alpha)));
A[0]=A1.f;
A[1]=A2.f;
A[2]=A3.f;
double *_cur = A+3;
_cur += saveD(_cur, 1, 2, &A1);
_cur += saveD(_cur, 1, 2, &A2);
saveD(_cur, 1, 2, &A3);
return A2;
}
static inline dualN twoCmtKaRateSSr2G(double *A, double *r2, dualN ka, dualN k20, dualN k23, dualN k32) {
(void)(ka);
dualN s=add2(add2(k23,k32),k20);
dualN beta=prodd2(0.5,(subtr2(s,sqrtD(subtr2(prod2(s,s),prod2(prodd2(4,k32),k20))))));
dualN alpha=div2(prod2(k32,k20),beta);
dualN A1;
iniD(0, -1, &A1);
dualN A2=div2(prodd2((*r2),k32),(prod2(beta,alpha)));
dualN A3=div2(prodd2((*r2),k23),(prod2(beta,alpha)));
A[0]=A1.f;
A[1]=A2.f;
A[2]=A3.f;
double *_cur = A+3;
_cur += saveD(_cur, 1, 2, &A1);
_cur += saveD(_cur, 1, 2, &A2);
saveD(_cur, 1, 2, &A3);
return A2;
}
static inline dualN twoCmtKaRateSStr1G(double *A, double *tinf, double *tau, double *r1, dualN ka, dualN k20, dualN k23, dualN k32) {
//has ka
dualN s=add2(add2(k23,k32),k20);
dualN beta=prodd2(0.5,(subtr2(s,sqrtD(subtr2(prod2(s,s),prod2(prodd2(4,k32),k20))))));
dualN alpha=div2(prod2(k32,k20),beta);
dualN eA=div2(expD(prodd2(((*tau)-(*tinf)),negD(alpha))),(subtrd2(1,expD(prod2d(negD(alpha),(*tau))))));
dualN eB=div2(expD(prodd2(((*tau)-(*tinf)),negD(beta))),(subtrd2(1,expD(prod2d(negD(beta),(*tau))))));
dualN eiA=expD(prod2d(negD(alpha),(*tinf)));
dualN eiB=expD(prod2d(negD(beta),(*tinf)));
dualN alpha2=prod2(alpha,alpha);
dualN alpha3=prod2(alpha2,alpha);
dualN beta2=prod2(beta,beta);
dualN beta3=prod2(beta2,beta);
dualN ka2=prod2(ka,ka);
dualN eKa=div2(expD(prodd2(((*tau)-(*tinf)),negD(ka))),(subtrd2(1,expD(prod2d(negD(ka),(*tau))))));
dualN eiKa=expD(prod2d(negD(ka),(*tinf)));
dualN A1=prod2(eKa,(subtr2(divd2((*r1),ka),div2(prod2d(eiKa,(*r1)),ka))));
dualN A2=add2(div2((subtr2(prod2(eA,(add2(add2(prod2(negD(alpha),(add2(subtr2(add2(div2(prodd2((*r1),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),(add2(negD(k32),ka))),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),(add2(negD(alpha),k32))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),(add2(negD(beta),k32))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3)))))),prod2(k32,(add2(subtr2(subtr2(div2(prodd2((*r1),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),k23),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),k23),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),k23),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))),prod2(k32,(add2(subtr2(add2(div2(prodd2((*r1),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),(add2(negD(k32),ka))),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),(add2(negD(alpha),k32))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),(add2(negD(beta),k32))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))))),prod2(eB,(add2(add2(prod2(negD(beta),(add2(subtr2(add2(div2(prodd2((*r1),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),(add2(negD(k32),ka))),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),(add2(negD(alpha),k32))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),(add2(negD(beta),k32))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3)))))),prod2(k32,(add2(subtr2(subtr2(div2(prodd2((*r1),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),k23),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),k23),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),k23),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))),prod2(k32,(add2(subtr2(add2(div2(prodd2((*r1),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),(add2(negD(k32),ka))),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),(add2(negD(alpha),k32))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),(add2(negD(beta),k32))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))))))),(add2(negD(alpha),beta))),prod2(prod2(ka,(add2(add2(div2(prod2(eA,(add2(negD(alpha),k32))),(prod2((add2(negD(alpha),beta)),(add2(negD(alpha),ka))))),div2(prod2(eB,(add2(negD(beta),k32))),(prod2((add2(negD(beta),ka)),(subtr2(alpha,beta)))))),div2(prod2(eKa,(subtr2(k32,ka))),(prod2((subtr2(beta,ka)),(subtr2(alpha,ka)))))))),(subtr2(divd2((*r1),ka),div2(prod2d(eiKa,(*r1)),ka)))));
dualN A3=add2(div2((subtr2(prod2(eA,(add2(add2(prod2(negD(alpha),(add2(subtr2(subtr2(div2(prodd2((*r1),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),k23),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),k23),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),k23),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3)))))),prod2(k23,(add2(subtr2(add2(div2(prodd2((*r1),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),(add2(negD(k32),ka))),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),(add2(negD(alpha),k32))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),(add2(negD(beta),k32))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))),prod2((add2(k20,k23)),(add2(subtr2(subtr2(div2(prodd2((*r1),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),k23),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),k23),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),k23),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))))),prod2(eB,(add2(add2(prod2(negD(beta),(add2(subtr2(subtr2(div2(prodd2((*r1),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),k23),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),k23),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),k23),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3)))))),prod2(k23,(add2(subtr2(add2(div2(prodd2((*r1),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),(add2(negD(k32),ka))),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),(add2(negD(alpha),k32))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),(add2(negD(beta),k32))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))),prod2((add2(k20,k23)),(add2(subtr2(subtr2(div2(prodd2((*r1),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiKa,(*r1)),k23),(add2(add2(prod2(beta,alpha),prod2(ka,(subtr2(negD(alpha),beta)))),ka2)))),div2(prod2(prod2(prod2d(eiA,(*r1)),ka),k23),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2(prod2d(eiB,(*r1)),ka),k23),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))))))),(add2(negD(alpha),beta))),prod2(prod2(prod2(ka,k23),(add2(add2(div2(eA,(prod2((add2(negD(alpha),beta)),(add2(negD(alpha),ka))))),div2(eB,(prod2((add2(negD(beta),ka)),(subtr2(alpha,beta)))))),div2(eKa,(prod2((subtr2(beta,ka)),(subtr2(alpha,ka)))))))),(subtr2(divd2((*r1),ka),div2(prod2d(eiKa,(*r1)),ka)))));
A[0]=A1.f;
A[1]=A2.f;
A[2]=A3.f;
double *_cur = A+3;
_cur += saveD(_cur, 1, 2, &A1);
_cur += saveD(_cur, 1, 2, &A2);
saveD(_cur, 1, 2, &A3);
return A2;
}
static inline dualN twoCmtKaRateSStr2G(double *A, double *tinf, double *tau, double *r2, dualN ka, dualN k20, dualN k23, dualN k32) {
//has ka
dualN E2=add2(k20,k23);
dualN E3=k32;
dualN s=add2(add2(k23,k32),k20);
dualN beta=prodd2(0.5,(subtr2(s,sqrtD(subtr2(prod2(s,s),prod2(prodd2(4,k32),k20))))));
dualN alpha=div2(prod2(k32,k20),beta);
dualN eA=div2(expD(prodd2(((*tau)-(*tinf)),negD(alpha))),(subtrd2(1,expD(prod2d(negD(alpha),(*tau))))));
dualN eB=div2(expD(prodd2(((*tau)-(*tinf)),negD(beta))),(subtrd2(1,expD(prod2d(negD(beta),(*tau))))));
dualN eiA=expD(prod2d(negD(alpha),(*tinf)));
dualN eiB=expD(prod2d(negD(beta),(*tinf)));
dualN alpha2=prod2(alpha,alpha);
dualN alpha3=prod2(alpha2,alpha);
dualN beta2=prod2(beta,beta);
dualN beta3=prod2(beta2,beta);
dualN A1;
iniD(0, -1, &A1);
dualN A2=div2((subtr2(prod2(eA,(add2(subtr2(prod2(E3,(add2(subtr2(div2(prodd2((*r2),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(add2(prod2(negD(k32),alpha),prod2(ka,(add2(negD(alpha),k32)))),alpha2))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(add2(prod2(negD(k32),beta),prod2(ka,(add2(negD(beta),k32)))),beta2))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3)))))),prod2(alpha,(add2(subtr2(div2(prodd2((*r2),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(add2(prod2(negD(k32),alpha),prod2(ka,(add2(negD(alpha),k32)))),alpha2))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(add2(prod2(negD(k32),beta),prod2(ka,(add2(negD(beta),k32)))),beta2))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))),prod2(k32,(add2(subtr2(div2(prodd2((*r2),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(prod2(negD(k23),alpha),prod2(ka,k23)))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(prod2(negD(k23),beta),prod2(ka,k23)))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))))),prod2(eB,(add2(subtr2(prod2(E3,(add2(subtr2(div2(prodd2((*r2),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(add2(prod2(negD(k32),alpha),prod2(ka,(add2(negD(alpha),k32)))),alpha2))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(add2(prod2(negD(k32),beta),prod2(ka,(add2(negD(beta),k32)))),beta2))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3)))))),prod2(beta,(add2(subtr2(div2(prodd2((*r2),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(add2(prod2(negD(k32),alpha),prod2(ka,(add2(negD(alpha),k32)))),alpha2))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(add2(prod2(negD(k32),beta),prod2(ka,(add2(negD(beta),k32)))),beta2))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))),prod2(k32,(add2(subtr2(div2(prodd2((*r2),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(prod2(negD(k23),alpha),prod2(ka,k23)))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(prod2(negD(k23),beta),prod2(ka,k23)))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))))))),(add2(negD(alpha),beta)));
dualN A3=div2((subtr2(prod2(eA,(add2(subtr2(prod2(E2,(add2(subtr2(div2(prodd2((*r2),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(prod2(negD(k23),alpha),prod2(ka,k23)))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(prod2(negD(k23),beta),prod2(ka,k23)))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3)))))),prod2(alpha,(add2(subtr2(div2(prodd2((*r2),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(prod2(negD(k23),alpha),prod2(ka,k23)))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(prod2(negD(k23),beta),prod2(ka,k23)))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))),prod2(k23,(add2(subtr2(div2(prodd2((*r2),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(add2(prod2(negD(k32),alpha),prod2(ka,(add2(negD(alpha),k32)))),alpha2))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(add2(prod2(negD(k32),beta),prod2(ka,(add2(negD(beta),k32)))),beta2))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))))),prod2(eB,(add2(subtr2(prod2(E2,(add2(subtr2(div2(prodd2((*r2),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(prod2(negD(k23),alpha),prod2(ka,k23)))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(prod2(negD(k23),beta),prod2(ka,k23)))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3)))))),prod2(beta,(add2(subtr2(div2(prodd2((*r2),k23),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(prod2(negD(k23),alpha),prod2(ka,k23)))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(prod2(negD(k23),beta),prod2(ka,k23)))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))),prod2(k23,(add2(subtr2(div2(prodd2((*r2),k32),(prod2(beta,alpha))),div2(prod2(prod2d(eiA,(*r2)),(add2(add2(prod2(negD(k32),alpha),prod2(ka,(add2(negD(alpha),k32)))),alpha2))),(add2(add2(prod2(negD(beta),alpha2),prod2(ka,(subtr2(prod2(beta,alpha),alpha2)))),alpha3)))),div2(prod2(prod2d(eiB,(*r2)),(add2(add2(prod2(negD(k32),beta),prod2(ka,(add2(negD(beta),k32)))),beta2))),(subtr2(add2(prod2(beta2,alpha),prod2(ka,(add2(prod2(negD(beta),alpha),beta2)))),beta3))))))))))),(add2(negD(alpha),beta)));
A[0]=A1.f;
A[1]=A2.f;
A[2]=A3.f;
double *_cur = A+3;
_cur += saveD(_cur, 1, 2, &A1);
_cur += saveD(_cur, 1, 2, &A2);
saveD(_cur, 1, 2, &A3);
return A2;
}
static inline dualN twoCmtKaRateG(double *A, double *Alast, double *t, double *b1, double *b2, double *r1, double *r2, dualN ka, dualN k20, dualN k23, dualN k32) {
// has ka
dualN A1last;
dualN A2last;
dualN A3last;
double *_cur = Alast+3;
_cur += restoreD(Alast[0], _cur, 1, 2, &A1last);
_cur += restoreD(Alast[1], _cur, 1, 2, &A2last);
restoreD(Alast[2], _cur, 1, 2, &A3last);
dualN E2=add2(k20,k23);
dualN s=add2(add2(k23,k32),k20);
dualN beta=prodd2(0.5,(subtr2(s,sqrtD(subtr2(prod2(s,s),prod2(prodd2(4,k32),k20))))));
dualN alpha=div2(prod2(k32,k20),beta);
dualN eKa=expD(prod2d(negD(ka),(*t)));
dualN eA=expD(prod2d(negD(alpha),(*t)));
dualN eB=expD(prod2d(negD(beta),(*t)));
dualN ka2=prod2(ka,ka);
dualN alpha2=prod2(alpha,alpha);
dualN alpha3=prod2(alpha2,alpha);
dualN beta2=prod2(beta,beta);
dualN beta3=prod2(beta2,beta);
dualN A1=subtr2(addd2((*b1),divd2((*r1),ka)),div2((prod2((subtrd2((*r1),prod2(A1last,ka))),eKa)),ka));
dualN A2=add2(subtr2(add2(addd2((*b2),div2((prod2((add2(subtr2(prod2d((subtr2(ka,k32)),(*r1)),prod2(A1last,ka2)),prod2(prod2(A1last,k32),ka))),eKa)),(add2(add2(ka2,prod2((subtr2(negD(beta),alpha)),ka)),prod2(alpha,beta))))),div2((prod2((subtr2(add2(add2(add2(prod2d((add2(subtr2(prod2((subtr2(k32,beta)),ka),prod2(beta,k32)),beta2)),(*r2)),prod2d(prod2((subtr2(k32,beta)),ka),(*r1))),prod2((add2(prod2(prod2((subtr2(subtr2(negD(A3last),A2last),A1last)),beta),k32),prod2((add2(A2last,A1last)),beta2))),ka)),prod2(prod2((add2(A3last,A2last)),beta2),k32)),prod2(A2last,beta3))),eB)),(add2(subtr2(prod2((subtr2(beta2,prod2(alpha,beta))),ka),beta3),prod2(alpha,beta2))))),div2((prod2((subtr2(add2(add2(add2(prod2d((add2(subtr2(prod2((subtr2(k32,alpha)),ka),prod2(alpha,k32)),alpha2)),(*r2)),prod2d(prod2((subtr2(k32,alpha)),ka),(*r1))),prod2((add2(prod2(prod2((subtr2(subtr2(negD(A3last),A2last),A1last)),alpha),k32),prod2((add2(A2last,A1last)),alpha2))),ka)),prod2(prod2((add2(A3last,A2last)),alpha2),k32)),prod2(A2last,alpha3))),eA)),(add2(subtr2(prod2((subtr2(prod2(alpha,beta),alpha2)),ka),prod2(alpha2,beta)),alpha3)))),div2((add2(prod2d(k32,(*r2)),prod2d(k32,(*r1)))),(prod2(alpha,beta))));
dualN A3=add2(subtr2(add2(div2(negD((prod2((subtr2(prod2d(k23,(*r1)),prod2(prod2(A1last,k23),ka))),eKa))),(add2(add2(ka2,prod2((subtr2(negD(beta),alpha)),ka)),prod2(alpha,beta)))),div2((prod2((add2(subtr2(add2(add2(add2(prod2d((subtr2(prod2(k23,ka),prod2(beta,k23))),(*r2)),prod2d(prod2(k23,ka),(*r1))),prod2((subtr2(add2(prod2(prod2((subtr2(negD(A2last),A1last)),beta),k23),prod2(A3last,beta2)),prod2(prod2(A3last,E2),beta))),ka)),prod2(prod2(A2last,beta2),k23)),prod2(A3last,beta3)),prod2(prod2(A3last,E2),beta2))),eB)),(add2(subtr2(prod2((subtr2(beta2,prod2(alpha,beta))),ka),beta3),prod2(alpha,beta2))))),div2((prod2((add2(subtr2(add2(add2(add2(prod2d((subtr2(prod2(k23,ka),prod2(alpha,k23))),(*r2)),prod2d(prod2(k23,ka),(*r1))),prod2((subtr2(add2(prod2(prod2((subtr2(negD(A2last),A1last)),alpha),k23),prod2(A3last,alpha2)),prod2(prod2(A3last,E2),alpha))),ka)),prod2(prod2(A2last,alpha2),k23)),prod2(A3last,alpha3)),prod2(prod2(A3last,E2),alpha2))),eA)),(add2(subtr2(prod2((subtr2(prod2(alpha,beta),alpha2)),ka),prod2(alpha2,beta)),alpha3)))),div2((add2(prod2d(k23,(*r2)),prod2d(k23,(*r1)))),(prod2(alpha,beta))));
A[0]=A1.f;
A[1]=A2.f;
A[2]=A3.f;
_cur = A+3;
_cur += saveD(_cur, 1, 2, &A1);
_cur += saveD(_cur, 1, 2, &A2);
saveD(_cur, 1, 2, &A3);
return A2;
}
static inline dualN twoCmtKaSSb1G(double *A, double *tau, double *b1, dualN ka, dualN k20, dualN k23, dualN k32) {
//has ka
dualN E2=add2(k20,k23);
dualN E3=k32;
dualN e2e3=add2(E2,E3);
dualN s=sqrtD(subtr2(prod2(e2e3,e2e3),prodd2(4,(subtr2(prod2(E2,E3),prod2(k23,k32))))));
dualN lambda1=prodd2(0.5,(add2(e2e3,s)));
dualN lambda2=prodd2(0.5,(subtr2(e2e3,s)));
dualN eKa=divd2(1,(subtrd2(1,expD(prodd2(-(*tau),ka)))));
dualN eL1=divd2(1,(subtrd2(1,expD(prodd2(-(*tau),lambda1)))));
dualN eL2=divd2(1,(subtrd2(1,expD(prodd2(-(*tau),lambda2)))));
dualN A1=prod2d(eKa,(*b1));
dualN A2=prod2(prod2d(ka,(*b1)),(add2(add2(div2(prod2(eL1,(subtr2(E3,lambda1))),(prod2((add2(negD(lambda1),lambda2)),(subtr2(ka,lambda1))))),div2(prod2(eL2,(subtr2(E3,lambda2))),(prod2((subtr2(lambda1,lambda2)),(subtr2(ka,lambda2)))))),div2(prod2(eKa,(subtr2(E3,ka))),(prod2((add2(negD(ka),lambda2)),(add2(negD(ka),lambda1))))))));
dualN A3=prod2(prod2(prod2d(ka,(*b1)),k23),(add2(add2(div2(eL1,(prod2((add2(negD(lambda1),lambda2)),(subtr2(ka,lambda1))))),div2(eL2,(prod2((subtr2(lambda1,lambda2)),(subtr2(ka,lambda2)))))),div2(eKa,(prod2((add2(negD(ka),lambda2)),(add2(negD(ka),lambda1))))))));
A[0]=A1.f;
A[1]=A2.f;
A[2]=A3.f;
double *_cur = A+3;
_cur += saveD(_cur, 1, 2, &A1);
_cur += saveD(_cur, 1, 2, &A2);
saveD(_cur, 1, 2, &A3);
return A2;
}
static inline dualN twoCmtKaSSb2G(double *A, double *tau, double *b2, dualN ka, dualN k20, dualN k23, dualN k32) {
(void)(ka);
dualN E2=add2(k20,k23);
dualN E3=k32;
dualN e2e3=add2(E2,E3);
dualN s=sqrtD(subtr2(prod2(e2e3,e2e3),prodd2(4,(subtr2(prod2(E2,E3),prod2(k23,k32))))));
dualN lambda1=prodd2(0.5,(add2(e2e3,s)));
dualN lambda2=prodd2(0.5,(subtr2(e2e3,s)));
dualN eL1=divd2(1,(subtrd2(1,expD(prodd2(-(*tau),lambda1)))));
dualN eL2=divd2(1,(subtrd2(1,expD(prodd2(-(*tau),lambda2)))));
dualN A1;
iniD(0, -1, &A1);
dualN A2=div2((subtr2(prod2(eL1,(subtr2(prodd2((*b2),E3),prodd2((*b2),lambda1)))),prod2(eL2,(subtr2(prodd2((*b2),E3),prodd2((*b2),lambda2)))))),(add2(negD(lambda1),lambda2)));
dualN A3=div2((subtr2(prod2(prod2d(eL1,(*b2)),k23),prod2(prod2d(eL2,(*b2)),k23))),(add2(negD(lambda1),lambda2)));
A[0]=A1.f;
A[1]=A2.f;
A[2]=A3.f;
double *_cur = A+3;
_cur += saveD(_cur, 1, 2, &A1);
_cur += saveD(_cur, 1, 2, &A2);
saveD(_cur, 1, 2, &A3);
return A2;
}
#endif
