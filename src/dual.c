#include <R.h>
#include <Rversion.h>
#include <Rinternals.h>
#include "../inst/include/RxODE.h"

#include "dual.h"

#undef beta
void DtwoCmtKaRate(double *A, double *Alast,
		   double *t, double *b1, double *b2,
		   double *r1, double *r2,
		   double *kad,  double *k20d, 
		   double *k23d, double *k32d) {
  dualN ka  = iniD(*kad,  0, 4);
  dualN k20 = iniD(*k20d, 2, 4);
  dualN k23 = iniD(*k23d, 3, 4);
  dualN k32 = iniD(*k32d, 4, 4);
  dualN A1last = iniD(Alast[0],-1,4);
  A1last.grad[0] = Alast[3];
  A1last.grad[1] = 0.0;
  A1last.grad[2] = Alast[4];
  A1last.grad[3] = Alast[5];
  dualN A2last = iniD(Alast[1],-1,4);
  A2last.grad[0] = Alast[6];
  A2last.grad[1] = Alast[7];
  A2last.grad[2] = Alast[8];
  A1last.grad[3] = Alast[9];
  dualN A3last = iniD(Alast[2],-1,4);
  A3last.grad[0] = Alast[10];
  A3last.grad[1] = Alast[11];
  A3last.grad[2] = Alast[12];
  A3last.grad[3] = Alast[13];
  dualN E2=add2(k20,k23);
  dualN s=add2(add2(k23,k32),k20);
  dualN beta=prodd2(0.5,(subtr2(s,sqrtD(subtr2(prod2(s,s),prod2(prodd2(4,k32),k20))))));
  dualN alpha=div2(prod2(k32,k20),beta);
  dualN eKa=expD(prodd2((*t),negD(ka)));
  dualN eA=expD(prodd2((*t),negD(alpha)));
  dualN eB=expD(prodd2((*t),negD(beta)));
  dualN ka2=prod2(ka,ka);
  dualN alpha2=prod2(alpha,alpha);
  dualN alpha3=prod2(alpha2,alpha);
  dualN beta2=prod2(beta,beta);
  dualN beta3=prod2(beta2,beta);
  dualN A1=subtr2(addd2((*b1),divd2((*r1),ka)),div2((prod2((subtrd2((*r1),prod2(A1last,ka))),eKa)),ka));
  dualN A2=add2(subtr2(add2(addd2((*b2),div2((prod2((add2(subtr2(prodd2((*r1),(subtr2(ka,k32))),prod2(A1last,ka2)),prod2(prod2(A1last,k32),ka))),eKa)),(add2(add2(ka2,prod2((subtr2(negD(beta),alpha)),ka)),prod2(alpha,beta))))),div2((prod2((subtr2(add2(add2(add2(prodd2((*r2),(add2(subtr2(prod2((subtr2(k32,beta)),ka),prod2(beta,k32)),beta2))),prodd2((*r1),prod2((subtr2(k32,beta)),ka))),prod2((add2(prod2(prod2((subtr2(subtr2(negD(A3last),A2last),A1last)),beta),k32),prod2((add2(A2last,A1last)),beta2))),ka)),prod2(prod2((add2(A3last,A2last)),beta2),k32)),prod2(A2last,beta3))),eB)),(add2(subtr2(prod2((subtr2(beta2,prod2(alpha,beta))),ka),beta3),prod2(alpha,beta2))))),div2((prod2((subtr2(add2(add2(add2(prodd2((*r2),(add2(subtr2(prod2((subtr2(k32,alpha)),ka),prod2(alpha,k32)),alpha2))),prodd2((*r1),prod2((subtr2(k32,alpha)),ka))),prod2((add2(prod2(prod2((subtr2(subtr2(negD(A3last),A2last),A1last)),alpha),k32),prod2((add2(A2last,A1last)),alpha2))),ka)),prod2(prod2((add2(A3last,A2last)),alpha2),k32)),prod2(A2last,alpha3))),eA)),(add2(subtr2(prod2((subtr2(prod2(alpha,beta),alpha2)),ka),prod2(alpha2,beta)),alpha3)))),div2((add2(prodd2((*r2),k32),prodd2((*r1),k32))),(prod2(alpha,beta))));
  dualN A3=add2(subtr2(add2(div2(negD((prod2((subtr2(prodd2((*r1),k23),prod2(prod2(A1last,k23),ka))),eKa))),(add2(add2(ka2,prod2((subtr2(negD(beta),alpha)),ka)),prod2(alpha,beta)))),div2((prod2((add2(subtr2(add2(add2(add2(prodd2((*r2),(subtr2(prod2(k23,ka),prod2(beta,k23)))),prodd2((*r1),prod2(k23,ka))),prod2((subtr2(add2(prod2(prod2((subtr2(negD(A2last),A1last)),beta),k23),prod2(A3last,beta2)),prod2(prod2(A3last,E2),beta))),ka)),prod2(prod2(A2last,beta2),k23)),prod2(A3last,beta3)),prod2(prod2(A3last,E2),beta2))),eB)),(add2(subtr2(prod2((subtr2(beta2,prod2(alpha,beta))),ka),beta3),prod2(alpha,beta2))))),div2((prod2((add2(subtr2(add2(add2(add2(prodd2((*r2),(subtr2(prod2(k23,ka),prod2(alpha,k23)))),prodd2((*r1),prod2(k23,ka))),prod2((subtr2(add2(prod2(prod2((subtr2(negD(A2last),A1last)),alpha),k23),prod2(A3last,alpha2)),prod2(prod2(A3last,E2),alpha))),ka)),prod2(prod2(A2last,alpha2),k23)),prod2(A3last,alpha3)),prod2(prod2(A3last,E2),alpha2))),eA)),(add2(subtr2(prod2((subtr2(prod2(alpha,beta),alpha2)),ka),prod2(alpha2,beta)),alpha3)))),div2((add2(prodd2((*r2),k23),prodd2((*r1),k23))),(prod2(alpha,beta))));
  A[0] = A1.f;
  A[1] = A2.f;
  A[2] = A3.f;
  A[3] = A1.grad[0];
  A[4] = A1.grad[2];
  A[5] = A1.grad[3];
  A[6] = A2.grad[0];
  A[7] = A2.grad[1];
  A[8] = A2.grad[2];
  A[9] = A2.grad[3];
  A[10] = A3.grad[0];
  A[11] = A3.grad[1];
  A[12] = A3.grad[2];
  A[13] = A3.grad[3];
}
