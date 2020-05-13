library(tidyverse)

itdepends::dep_usage_pkg("RxODE") %>%
    count(pkg, sort = TRUE) %>%
    as.data.frame


if (FALSE){
  ## Creating optimized models for linCmt Advan
  ## rate 1-3 cmt
  m <- RxODE({
    A1=r1/ka-((r1+(-b1-A1last)*ka)*exp(-ka*t))/ka
    A2=(((ka-k32)*r1+(-b1-A1last)*ka^2+(b1+A1last)*k32*ka)*exp(-ka*t))/(ka^2+(-beta-alpha)*ka+alpha*beta)+((((k32-beta)*ka-beta*k32+beta^2)*r2+(k32-beta)*ka*r1+((-b2-b1-A3last-A2last-A1last)*beta*k32+(b2+b1+A2last+A1last)*beta^2)*ka+(b2+A3last+A2last)*beta^2*k32+(-b2-A2last)*beta^3)*exp(-beta*t))/((beta^2-alpha*beta)*ka-beta^3+alpha*beta^2)-((((k32-alpha)*ka-alpha*k32+alpha^2)*r2+(k32-alpha)*ka*r1+((-alpha*b2-alpha*b1+(-A3last-A2last-A1last)*alpha)*k32+alpha^2*b2+alpha^2*b1+(A2last+A1last)*alpha^2)*ka+(alpha^2*b2+(A3last+A2last)*alpha^2)*k32-alpha^3*b2-A2last*alpha^3)*exp(-alpha*t))/((alpha*beta-alpha^2)*ka-alpha^2*beta+alpha^3)+(k32*r2+k32*r1)/(alpha*beta)
    A3=-((k23*r1+(-b1-A1last)*k23*ka)*exp(-ka*t))/(ka^2+(-beta-alpha)*ka+alpha*beta)+(((k23*ka-beta*k23)*r2+k23*ka*r1+((-b2-b1-A2last-A1last)*beta*k23+A3last*beta^2-A3last*E2*beta)*ka+(b2+A2last)*beta^2*k23-A3last*beta^3+A3last*E2*beta^2)*exp(-beta*t))/((beta^2-alpha*beta)*ka-beta^3+alpha*beta^2)-(((k23*ka-alpha*k23)*r2+k23*ka*r1+((-alpha*b2-alpha*b1+(-A2last-A1last)*alpha)*k23+A3last*alpha^2-A3last*E2*alpha)*ka+(alpha^2*b2+A2last*alpha^2)*k23-A3last*alpha^3+A3last*E2*alpha^2)*exp(-alpha*t))/((alpha*beta-alpha^2)*ka-alpha^2*beta+alpha^3)+(k23*r2+k23*r1)/(alpha*beta)
  })

  message(rxOptExpr(rxNorm(m)))

  m <- RxODE({
    A1=r1/ka-((r1+(-b1-A1last)*ka)*exp(-ka*t))/ka
    A2=-(((lam3^3+(-ka-k42-k32)*lam3^2+((k42+k32)*ka+k32*k42)*lam3-k32*k42*ka)*r2+(-ka*lam3^2+(k42+k32)*ka*lam3-k32*k42*ka)*r1+(-b2-A2last)*lam3^4+((b2+b1+A2last+A1last)*ka+(b2+A4last+A2last)*k42+(b2+A3last+A2last)*k32)*lam3^3+(((-b2-b1-A4last-A2last-A1last)*k42+(-b2-b1-A3last-A2last-A1last)*k32)*ka+(-b2-A4last-A3last-A2last)*k32*k42)*lam3^2+(b2+b1+A4last+A3last+A2last+A1last)*k32*k42*ka*lam3)*exp(-lam3*t))/(lam3^4+(-lam2-lam1-ka)*lam3^3+((lam1+ka)*lam2+ka*lam1)*lam3^2-ka*lam1*lam2*lam3)+(((lam2^3+(-ka-k42-k32)*lam2^2+((k42+k32)*ka+k32*k42)*lam2-k32*k42*ka)*r2+(-ka*lam2^2+(k42+k32)*ka*lam2-k32*k42*ka)*r1+(-b2-A2last)*lam2^4+((b2+b1+A2last+A1last)*ka+(b2+A4last+A2last)*k42+(b2+A3last+A2last)*k32)*lam2^3+(((-b2-b1-A4last-A2last-A1last)*k42+(-b2-b1-A3last-A2last-A1last)*k32)*ka+(-b2-A4last-A3last-A2last)*k32*k42)*lam2^2+(b2+b1+A4last+A3last+A2last+A1last)*k32*k42*ka*lam2)*exp(-lam2*t))/((lam2^3+(-lam1-ka)*lam2^2+ka*lam1*lam2)*lam3-lam2^4+(lam1+ka)*lam2^3-ka*lam1*lam2^2)-(((lam1^3+(-ka-k42-k32)*lam1^2+((k42+k32)*ka+k32*k42)*lam1-k32*k42*ka)*r2+(-ka*lam1^2+(k42+k32)*ka*lam1-k32*k42*ka)*r1+(-b2-A2last)*lam1^4+((b2+b1+A2last+A1last)*ka+(b2+A4last+A2last)*k42+(b2+A3last+A2last)*k32)*lam1^3+(((-b2-b1-A4last-A2last-A1last)*k42+(-b2-b1-A3last-A2last-A1last)*k32)*ka+(-b2-A4last-A3last-A2last)*k32*k42)*lam1^2+(b2+b1+A4last+A3last+A2last+A1last)*k32*k42*ka*lam1)*exp(-lam1*t))/(((lam1^2-ka*lam1)*lam2-lam1^3+ka*lam1^2)*lam3+(ka*lam1^2-lam1^3)*lam2+lam1^4-ka*lam1^3)-(((ka^2+(-k42-k32)*ka+k32*k42)*r1+(-b1-A1last)*ka^3+((b1+A1last)*k42+(b1+A1last)*k32)*ka^2+(-b1-A1last)*k32*k42*ka)*exp(-ka*t))/(((lam1-ka)*lam2-ka*lam1+ka^2)*lam3+(ka^2-ka*lam1)*lam2+ka^2*lam1-ka^3)+(k32*k42*r2+k32*k42*r1)/(lam1*lam2*lam3)
    A3=(((k23*lam3^2+(-k23*ka-k23*k42)*lam3+k23*k42*ka)*r2+(k23*k42*ka-k23*ka*lam3)*r1+A3last*lam3^4+(-A3last*ka-A3last*k42+(-b2-A2last)*k23-A3last*E2)*lam3^3+((A3last*k42+(b2+b1+A2last+A1last)*k23+A3last*E2)*ka+(-A3last*k24+(b2+A4last+A2last)*k23+A3last*E2)*k42)*lam3^2+(A3last*k24+(-b2-b1-A4last-A2last-A1last)*k23-A3last*E2)*k42*ka*lam3)*exp(-lam3*t))/(lam3^4+(-lam2-lam1-ka)*lam3^3+((lam1+ka)*lam2+ka*lam1)*lam3^2-ka*lam1*lam2*lam3)-(((k23*lam2^2+(-k23*ka-k23*k42)*lam2+k23*k42*ka)*r2+(k23*k42*ka-k23*ka*lam2)*r1+A3last*lam2^4+(-A3last*ka-A3last*k42+(-b2-A2last)*k23-A3last*E2)*lam2^3+((A3last*k42+(b2+b1+A2last+A1last)*k23+A3last*E2)*ka+(-A3last*k24+(b2+A4last+A2last)*k23+A3last*E2)*k42)*lam2^2+(A3last*k24+(-b2-b1-A4last-A2last-A1last)*k23-A3last*E2)*k42*ka*lam2)*exp(-lam2*t))/((lam2^3+(-lam1-ka)*lam2^2+ka*lam1*lam2)*lam3-lam2^4+(lam1+ka)*lam2^3-ka*lam1*lam2^2)+(((k23*lam1^2+(-k23*ka-k23*k42)*lam1+k23*k42*ka)*r2+(k23*k42*ka-k23*ka*lam1)*r1+A3last*lam1^4+(-A3last*ka-A3last*k42+(-b2-A2last)*k23-A3last*E2)*lam1^3+((A3last*k42+(b2+b1+A2last+A1last)*k23+A3last*E2)*ka+(-A3last*k24+(b2+A4last+A2last)*k23+A3last*E2)*k42)*lam1^2+(A3last*k24+(-b2-b1-A4last-A2last-A1last)*k23-A3last*E2)*k42*ka*lam1)*exp(-lam1*t))/(((lam1^2-ka*lam1)*lam2-lam1^3+ka*lam1^2)*lam3+(ka*lam1^2-lam1^3)*lam2+lam1^4-ka*lam1^3)+(((k23*ka-k23*k42)*r1+(-b1-A1last)*k23*ka^2+(b1+A1last)*k23*k42*ka)*exp(-ka*t))/(((lam1-ka)*lam2-ka*lam1+ka^2)*lam3+(ka^2-ka*lam1)*lam2+ka^2*lam1-ka^3)+(k23*k42*r2+k23*k42*r1)/(lam1*lam2*lam3)
    A4=(((k24*lam3^2+(-k24*ka-k24*k32)*lam3+k24*k32*ka)*r2+(k24*k32*ka-k24*ka*lam3)*r1+A4last*lam3^4+(-A4last*ka-A4last*k32+(-b2-A2last)*k24-A4last*E2)*lam3^3+((A4last*k32+(b2+b1+A2last+A1last)*k24+A4last*E2)*ka+((b2+A3last+A2last)*k24-A4last*k23+A4last*E2)*k32)*lam3^2+((-b2-b1-A3last-A2last-A1last)*k24+A4last*k23-A4last*E2)*k32*ka*lam3)*exp(-lam3*t))/(lam3^4+(-lam2-lam1-ka)*lam3^3+((lam1+ka)*lam2+ka*lam1)*lam3^2-ka*lam1*lam2*lam3)-(((k24*lam2^2+(-k24*ka-k24*k32)*lam2+k24*k32*ka)*r2+(k24*k32*ka-k24*ka*lam2)*r1+A4last*lam2^4+(-A4last*ka-A4last*k32+(-b2-A2last)*k24-A4last*E2)*lam2^3+((A4last*k32+(b2+b1+A2last+A1last)*k24+A4last*E2)*ka+((b2+A3last+A2last)*k24-A4last*k23+A4last*E2)*k32)*lam2^2+((-b2-b1-A3last-A2last-A1last)*k24+A4last*k23-A4last*E2)*k32*ka*lam2)*exp(-lam2*t))/((lam2^3+(-lam1-ka)*lam2^2+ka*lam1*lam2)*lam3-lam2^4+(lam1+ka)*lam2^3-ka*lam1*lam2^2)+(((k24*lam1^2+(-k24*ka-k24*k32)*lam1+k24*k32*ka)*r2+(k24*k32*ka-k24*ka*lam1)*r1+A4last*lam1^4+(-A4last*ka-A4last*k32+(-b2-A2last)*k24-A4last*E2)*lam1^3+((A4last*k32+(b2+b1+A2last+A1last)*k24+A4last*E2)*ka+((b2+A3last+A2last)*k24-A4last*k23+A4last*E2)*k32)*lam1^2+((-b2-b1-A3last-A2last-A1last)*k24+A4last*k23-A4last*E2)*k32*ka*lam1)*exp(-lam1*t))/(((lam1^2-ka*lam1)*lam2-lam1^3+ka*lam1^2)*lam3+(ka*lam1^2-lam1^3)*lam2+lam1^4-ka*lam1^3)+(((k24*ka-k24*k32)*r1+(-b1-A1last)*k24*ka^2+(b1+A1last)*k24*k32*ka)*exp(-ka*t))/(((lam1-ka)*lam2-ka*lam1+ka^2)*lam3+(ka^2-ka*lam1)*lam2+ka^2*lam1-ka^3)+(k24*k32*r2+k24*k32*r1)/(lam1*lam2*lam3)
  })

  message(rxOptExpr(rxNorm(m)))

  ## Now 1 compartment oral
  m <-RxODE({
    A1=(b1+A1last)*exp(-ka*t)
    A2=(((b2+b1+A2last+A1last)*ka+(-b2-A2last)*k20)*exp(-k20*t))/(ka-k20)-((b1+A1last)*ka*exp(-ka*t))/(ka-k20)
  })

  message(rxOptExpr(rxNorm(m)))

  ## 2 compartment oral

  m <- RxODE({
    A1=(b1+A1last)*exp(-ka*t)
    A2=-(((b1+A1last)*ka^2+(-b1-A1last)*k32*ka)*exp(-ka*t))/(ka^2+(-beta-alpha)*ka+alpha*beta)-((((b2+b1+A3last+A2last+A1last)*k32+(-b2-b1-A2last-A1last)*beta)*ka+(-b2-A3last-A2last)*beta*k32+(b2+A2last)*beta^2)*exp(-beta*t))/((beta-alpha)*ka-beta^2+alpha*beta)+((((b2+b1+A3last+A2last+A1last)*k32-alpha*b2-alpha*b1+(-A2last-A1last)*alpha)*ka+((-A3last-A2last)*alpha-alpha*b2)*k32+alpha^2*b2+A2last*alpha^2)*exp(-alpha*t))/((beta-alpha)*ka-alpha*beta+alpha^2);
    A3=((b1+A1last)*k23*ka*exp(-ka*t))/(ka^2+(-beta-alpha)*ka+alpha*beta)-((((b2+b1+A2last+A1last)*k23-A3last*beta+A3last*E2)*ka+(-b2-A2last)*beta*k23+A3last*beta^2-A3last*E2*beta)*exp(-beta*t))/((beta-alpha)*ka-beta^2+alpha*beta)+((((b2+b1+A2last+A1last)*k23-A3last*alpha+A3last*E2)*ka+(-alpha*b2-A2last*alpha)*k23+A3last*alpha^2-A3last*E2*alpha)*exp(-alpha*t))/((beta-alpha)*ka-alpha*beta+alpha^2)
  })

  message(rxOptExpr(rxNorm(m)))

  ## 3 compartment oral
  m <- RxODE({
    A1=(b1+A1last)*exp(-ka*t);
    A2=(((b2+A2last)*lam3^3+((-b2-b1-A2last-A1last)*ka+(-b2-A4last-A2last)*k42+(-b2-A3last-A2last)*k32)*lam3^2+(((b2+b1+A4last+A2last+A1last)*k42+(b2+b1+A3last+A2last+A1last)*k32)*ka+(b2+A4last+A3last+A2last)*k32*k42)*lam3+(-b2-b1-A4last-A3last-A2last-A1last)*k32*k42*ka)*exp(-lam3*t))/(lam3^3+(-lam2-lam1-ka)*lam3^2+((lam1+ka)*lam2+ka*lam1)*lam3-ka*lam1*lam2)-(((b2+A2last)*lam2^3+((-b2-b1-A2last-A1last)*ka+(-b2-A4last-A2last)*k42+(-b2-A3last-A2last)*k32)*lam2^2+(((b2+b1+A4last+A2last+A1last)*k42+(b2+b1+A3last+A2last+A1last)*k32)*ka+(b2+A4last+A3last+A2last)*k32*k42)*lam2+(-b2-b1-A4last-A3last-A2last-A1last)*k32*k42*ka)*exp(-lam2*t))/((lam2^2+(-lam1-ka)*lam2+ka*lam1)*lam3-lam2^3+(lam1+ka)*lam2^2-ka*lam1*lam2)+(((b2+A2last)*lam1^3+((-b2-b1-A2last-A1last)*ka+(-b2-A4last-A2last)*k42+(-b2-A3last-A2last)*k32)*lam1^2+(((b2+b1+A4last+A2last+A1last)*k42+(b2+b1+A3last+A2last+A1last)*k32)*ka+(b2+A4last+A3last+A2last)*k32*k42)*lam1+(-b2-b1-A4last-A3last-A2last-A1last)*k32*k42*ka)*exp(-lam1*t))/(((lam1-ka)*lam2-lam1^2+ka*lam1)*lam3+(ka*lam1-lam1^2)*lam2+lam1^3-ka*lam1^2)+(((b1+A1last)*ka^3+((-b1-A1last)*k42+(-b1-A1last)*k32)*ka^2+(b1+A1last)*k32*k42*ka)*exp(-ka*t))/(((lam1-ka)*lam2-ka*lam1+ka^2)*lam3+(ka^2-ka*lam1)*lam2+ka^2*lam1-ka^3)
    A3=((A3last*lam3^3+(-A3last*ka-A3last*k42+(-b2-A2last)*k23-A3last*E2)*lam3^2+((A3last*k42+(b2+b1+A2last+A1last)*k23+A3last*E2)*ka+(-A3last*k24+(b2+A4last+A2last)*k23+A3last*E2)*k42)*lam3+(A3last*k24+(-b2-b1-A4last-A2last-A1last)*k23-A3last*E2)*k42*ka)*exp(-lam3*t))/(lam3^3+(-lam2-lam1-ka)*lam3^2+((lam1+ka)*lam2+ka*lam1)*lam3-ka*lam1*lam2)-((A3last*lam2^3+(-A3last*ka-A3last*k42+(-b2-A2last)*k23-A3last*E2)*lam2^2+((A3last*k42+(b2+b1+A2last+A1last)*k23+A3last*E2)*ka+(-A3last*k24+(b2+A4last+A2last)*k23+A3last*E2)*k42)*lam2+(A3last*k24+(-b2-b1-A4last-A2last-A1last)*k23-A3last*E2)*k42*ka)*exp(-lam2*t))/((lam2^2+(-lam1-ka)*lam2+ka*lam1)*lam3-lam2^3+(lam1+ka)*lam2^2-ka*lam1*lam2)+((A3last*lam1^3+(-A3last*ka-A3last*k42+(-b2-A2last)*k23-A3last*E2)*lam1^2+((A3last*k42+(b2+b1+A2last+A1last)*k23+A3last*E2)*ka+(-A3last*k24+(b2+A4last+A2last)*k23+A3last*E2)*k42)*lam1+(A3last*k24+(-b2-b1-A4last-A2last-A1last)*k23-A3last*E2)*k42*ka)*exp(-lam1*t))/(((lam1-ka)*lam2-lam1^2+ka*lam1)*lam3+(ka*lam1-lam1^2)*lam2+lam1^3-ka*lam1^2)-(((b1+A1last)*k23*ka^2+(-b1-A1last)*k23*k42*ka)*exp(-ka*t))/(((lam1-ka)*lam2-ka*lam1+ka^2)*lam3+(ka^2-ka*lam1)*lam2+ka^2*lam1-ka^3)
    A4=((A4last*lam3^3+(-A4last*ka-A4last*k32+(-b2-A2last)*k24-A4last*E2)*lam3^2+((A4last*k32+(b2+b1+A2last+A1last)*k24+A4last*E2)*ka+((b2+A3last+A2last)*k24-A4last*k23+A4last*E2)*k32)*lam3+((-b2-b1-A3last-A2last-A1last)*k24+A4last*k23-A4last*E2)*k32*ka)*exp(-lam3*t))/(lam3^3+(-lam2-lam1-ka)*lam3^2+((lam1+ka)*lam2+ka*lam1)*lam3-ka*lam1*lam2)-((A4last*lam2^3+(-A4last*ka-A4last*k32+(-b2-A2last)*k24-A4last*E2)*lam2^2+((A4last*k32+(b2+b1+A2last+A1last)*k24+A4last*E2)*ka+((b2+A3last+A2last)*k24-A4last*k23+A4last*E2)*k32)*lam2+((-b2-b1-A3last-A2last-A1last)*k24+A4last*k23-A4last*E2)*k32*ka)*exp(-lam2*t))/((lam2^2+(-lam1-ka)*lam2+ka*lam1)*lam3-lam2^3+(lam1+ka)*lam2^2-ka*lam1*lam2)+((A4last*lam1^3+(-A4last*ka-A4last*k32+(-b2-A2last)*k24-A4last*E2)*lam1^2+((A4last*k32+(b2+b1+A2last+A1last)*k24+A4last*E2)*ka+((b2+A3last+A2last)*k24-A4last*k23+A4last*E2)*k32)*lam1+((-b2-b1-A3last-A2last-A1last)*k24+A4last*k23-A4last*E2)*k32*ka)*exp(-lam1*t))/(((lam1-ka)*lam2-lam1^2+ka*lam1)*lam3+(ka*lam1-lam1^2)*lam2+lam1^3-ka*lam1^2)-(((b1+A1last)*k24*ka^2+(-b1-A1last)*k24*k32*ka)*exp(-ka*t))/(((lam1-ka)*lam2-ka*lam1+ka^2)*lam3+(ka^2-ka*lam1)*lam2+ka^2*lam1-ka^3)
  })
  
  message(rxOptExpr(rxNorm(m)))

  ## Now 1-3 compartment bolus

  ## 2 compartment infusion
  m <- RxODE({
    A1=(((A1last*E2+r1+A2last*k21)-A1last*lambda1)*exp(-t*lambda1)-((A1last*E2+r1+A2last*k21)-A1last*lambda2)*exp(-t*lambda2))/(lambda2-lambda1) + r1*E2*(1/(lambda1*lambda2)+exp(-t*lambda1)/(lambda1*(lambda1-lambda2))-exp(-t*lambda2)/(lambda2*(lambda1-lambda2)))
    A2=(((A2last*E1+A1last*k12)-A2last*lambda1)*exp(-t*lambda1)-((A2last*E1+A1last*k12)-A2last*lambda2)*exp(-t*lambda2))/(lambda2-lambda1)+r1*k12*(1/(lambda1*lambda2)+exp(-t*lambda1)/(lambda1*(lambda1-lambda2))-exp(-t*lambda2)/(lambda2*(lambda1-lambda2)))
  })

  message(rxOptExpr(rxNorm(m)))

  m <- RxODE({
    A1=A1last*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3))) + exp(-t*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2)) + r1*((E2*E3)/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    A2 = A2last*(exp(-t*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3))) + exp(-t*lambda1)*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2)) + r1*k12*(E3/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E3-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E3-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E3-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
    A3=A3last*(exp(-t*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))+exp(-t*lambda1)*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2))+r1*k13*(E2/(lambda1*lambda2*lambda3)-exp(-t*lambda1)*(E2-lambda1)/(lambda1*(lambda2-lambda1)*(lambda3-lambda1))-exp(-t*lambda2)*(E2-lambda2)/(lambda2*(lambda1-lambda2)*(lambda3-lambda2))-exp(-t*lambda3)*(E2-lambda3)/(lambda3*(lambda1-lambda3)*(lambda2-lambda3)))
  })

  message(rxOptExpr(rxNorm(m)))

  m <- RxODE({
    A1 = (((A1last*E2+A2last*k21)-A1last*lambda1)*exp(-t*lambda1)-((A1last*E2+A2last*k21)-A1last*lambda2)*exp(-t*lambda2))/(lambda2-lambda1);
    A2 = (((A2last*E1+A1last*k12)-A2last*lambda1)*exp(-t*lambda1)-((A2last*E1+A1last*k12)-A2last*lambda2)*exp(-t*lambda2))/(lambda2-lambda1)
  })

  message(rxOptExpr(rxNorm(m)))

  m <- RxODE({
    A1= A1last*(exp(-t*lambda1)*(E2-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E2-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E2-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3))) + exp(-t*lambda1)*(C-B*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(B*lambda2-C)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(B*lambda3-C)/((lambda1-lambda3)*(lambda3-lambda2))
    A2=A2last*(exp(-t*lambda1)*(E1-lambda1)*(E3-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E3-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E3-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))+exp(-t*lambda1)*(I-A1last*k12*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k12*lambda2-I)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k12*lambda3-I)/((lambda1-lambda3)*(lambda3-lambda2))
    A3=A3last*(exp(-t*lambda1)*(E1-lambda1)*(E2-lambda1)/((lambda2-lambda1)*(lambda3-lambda1))+exp(-t*lambda2)*(E1-lambda2)*(E2-lambda2)/((lambda1-lambda2)*(lambda3-lambda2))+exp(-t*lambda3)*(E1-lambda3)*(E2-lambda3)/((lambda1-lambda3)*(lambda2-lambda3)))+exp(-t*lambda1)*(J-A1last*k13*lambda1)/((lambda1-lambda2)*(lambda1-lambda3))+exp(-t*lambda2)*(A1last*k13*lambda2-J)/((lambda1-lambda2)*(lambda2-lambda3))+exp(-t*lambda3)*(A1last*k13*lambda3-J)/((lambda1-lambda3)*(lambda3-lambda2))
  })

  message(rxOptExpr(rxNorm(m)))
  
}
