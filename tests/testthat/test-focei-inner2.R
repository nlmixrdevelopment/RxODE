## library(RxODE)
library(testthat)

#rxPermissive({
devtools::load_all()

context("Wang 2007 -- prop -- Inner Test")

## Another example (pred  Wang2007)
Wang2007 <- structure(list(ID = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L, 10L, 10L), TIME = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1), DV = c(10.68, 3.6837, 10.402, 6.454, 9.8814, 5.8565, 9.3408, 5.6209, 10.082, 6.7583, 9.8938, 6.5049, 9.8908, 6.9557, 10.234, 6.4488, 9.9882, 6.7112, 9.6736, 6.6402)), class = "data.frame", row.names = c(NA, -20L))
Wang2007$EVID <- 0
Wang2007$AMT <- NA

dat2 <- Wang2007[Wang2007$TIME == 0, ]
dat2$EVID <- 101
dat2$AMT <- 10;
dat2$DV <- NA;

dat2 <- rbind(dat2, Wang2007)
dat2 <- dat2[(order(dat2$ID, -dat2$EVID, dat2$TIME)), ]

write.table(dat2 %>% select(-EVID), "Wang2007amt.txt", row.names=FALSE, na=".")


df <- structure(list(ID = 1:10, ETA1 = c(0.0715450734486005, 0.00570452879884164, 0.02494718484644, 0.0319502378033361, -0.00478736927193241, 0.00397892241288246,  -0.0118035085000642, 0.00588012198219381, -0.00313671906679696,  -0.00066641434061578),
                     OBJI = c(5.20178380742227, 3.77604874005119, 3.72057773742475, 3.79581329127695, 3.84803907890983, 3.7725678919448, 3.92866450763499, 3.76425408681763, 3.83096377569697, 3.81883716853077)),
                class = "data.frame", row.names = c(NA, -10L))

ofv <- 39.4575500857102

inner.pred.prop <- RxODE("rx_pred_=10*exp(-THETA[1]*t*exp(ETA[1]));\nrx__sens_rx_pred__BY_ETA_1___=-10*THETA[1]*t*exp(ETA[1])*exp(-THETA[1]*t*exp(ETA[1]));\nrx_r_=100*Rx_pow_di(THETA[2],2)*exp(-2*THETA[1]*t*exp(ETA[1]));\nrx__sens_rx_r__BY_ETA_1___=-200*THETA[1]*Rx_pow_di(THETA[2],2)*t*exp(ETA[1])*exp(-2*THETA[1]*t*exp(ETA[1]));\n")

rxInv <- rxSymInvCholCreate(matrix(c(.04)))

Wang2007a <- Wang2007[Wang2007$ID == 1, ];

ini <- RxODE:::.foceiSetup(inner.pred.prop, Wang2007a, c(0.5, sqrt(.1)), etaMat=matrix(0.07154500), rxInv = rxInv,
                           odeOpt=foceiControl(maxOuterIterations=0, maxInnerIterations=0));

expect_equal(round(RxODE:::foceiOfv(ini), 3), round(5.20178378955669, 3))

rxSolveFree();

## Use NONMEM's selected ETAs.
etaMat <- matrix(c(0.07154500, 0.00570450, 0.02494700, 0.03195000, -0.00478740, 0.00397890, -0.01180400, 0.00588010,
                   -0.00313670,
                   -0.00066641), ncol=1)

ini <- RxODE:::.foceiSetup(inner.pred.prop, Wang2007, c(0.5, sqrt(.1)), etaMat=etaMat, rxInv = rxInv,
                           odeOpt=foceiControl(maxOuterIterations=0, maxInnerIterations=0));

df <- structure(list(ID = 1:10,
                     ETA1 = c(0.071545, 0.0057045, 0.024947, 0.03195, -0.0047874, 0.0039789, -0.011804, 0.0058801, -0.0031367,  -0.00066641),
                     OBJI = c(5.20178378955669, 3.77604873381688, 3.72057769601626, 3.79581323733943, 3.84803907238135, 3.77256788710782, 3.92866440452002, 3.76425408205748, 3.83096377975989, 3.8188371694598)),
                class = "data.frame", row.names = c(NA,  -10L));

expect_equal(round(RxODE:::foceiOfv(ini), 3), round(ofv, 3));
expect_equal(round(RxODE:::foceiEtas(), 3), round(df, 3));

rxSolveFree();


ini <- RxODE:::.foceiSetup(inner.pred.prop, Wang2007, c(0.5, sqrt(.1)), rxInv = rxInv);

expect_equal(round(RxODE:::foceiOfv(ini), 3), round(ofv, 3));


expect_equal(round(RxODE:::foceiEtas(), 3), round(df, 3));

rxSolveFree();


## Now Check Additive
## if (FALSE){

##     mypar1 = function ()
##     {
##         ke = theta[1] * exp(eta[1]);
##     }

##     mypar2 = function ()
##     {
##         k = theta[1] * exp(eta[1]);
##         v = 1
##     }

##     mod <- RxODE({
##         ipre = 10 * exp(-ke * t)
##     })

##     err3 <- function(f){
##         return(add(0.1));
##     }

##     pred <- function() ipre

##     m2b <- rxSymPySetupPred(mod, pred, mypar1, err3)

##     ## From inner.add <- m2b$inner

##     err4 <- function(f){
##         return(add(.1) + prop(0.2));
##     }

##     m2b <- rxSymPySetupPred(mod, pred, mypar1, err4)

##     err4 <- function(f){
##         return(add(.1) + prop(0.2));
##     }

##     err5 <- function(f){
##         return(2 * (THETA[2] + THETA[3] * ipre) ^ 2)
##     }

##     m2b <- rxSymPySetupPred(mod, pred, mypar1, err5)

##     ## add+ prop  = a*e +b*f*e
##     ## linearized fo on e = a+b*f
## }

context("Wang 2007 -- add -- Inner Test")

inner.add <- RxODE("rx_pred_=10*exp(-THETA[1]*t*exp(ETA[1]));\nrx__sens_rx_pred__BY_ETA_1___=-10*THETA[1]*t*exp(ETA[1])*exp(-THETA[1]*t*exp(ETA[1]));\nrx_r_=Rx_pow_di(THETA[2],2);\nrx__sens_rx_r__BY_ETA_1___=0;\n")


ini <- RxODE:::.foceiSetup(inner.add, Wang2007, c(0.5, sqrt(.1)), rxInv = rxInv);

ofv <- -2.05879522891223
df <- structure(list(ID = 1:10, ETA1 = c(0.583231076063623, -0.10164742615365, 0.0538695583777434, 0.113993027748088, -0.182309811413365, -0.11507939655624,  -0.235037475634943, -0.100276674118276, -0.169771301855972, -0.150905282267601),
                     OBJI = c(11.9606738350415, -1.19254345514792, -2.78878510918067, 1.77277064630383, -2.04088896880191, -2.61087183257036, -1.27355910316162,  -2.2689961915846, -2.24616212916238, -1.37043292064807)),
                class = "data.frame", row.names = c(NA, -10L));

expect_equal(round(RxODE:::foceiOfv(ini), 3), round(ofv, 3));

expect_equal(round(RxODE:::foceiEtas(), 3), round(df, 3));

rxSolveFree();

context("Wang 2007 -- add+prop -- Inner Test")

inner.add.prop <- RxODE("rx_pred_=10*exp(-THETA[1]*t*exp(ETA[1]));\nrx__sens_rx_pred__BY_ETA_1___=-10*THETA[1]*t*exp(ETA[1])*exp(-THETA[1]*t*exp(ETA[1]));\nrx_r_=Rx_pow_di(THETA[2],2)+100*Rx_pow_di(THETA[3],2)*exp(-2*THETA[1]*t*exp(ETA[1]));\nrx__sens_rx_r__BY_ETA_1___=-100*THETA[1]*Rx_pow_di(THETA[3],2)*t*exp(ETA[1])*exp(-2*THETA[1]*t*exp(ETA[1]));\n")

ini <- RxODE:::.foceiSetup(inner.add.prop, Wang2007, c(0.5, sqrt(.1), sqrt(.2)), rxInv = rxInv);

expect_equal(round(RxODE:::foceiOfv(ini), 3), round(51.7412183645535, 3));

expect_equal(round(RxODE:::foceiEtas(), 3), structure(list(ID = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                                           ETA1 = c(0.045, 0.01, 0.022, 0.026, 0, 0.01, 0, 0.01, 0, 0.008),
                                                           OBJI = c(5.802, 5.088, 5.059, 5.096, 5.127, 5.087, 5.17, 5.082, 5.119, 5.111)),
                                                      row.names = c(NA, -10L), class = "data.frame"))
rxSolveFree();



## Lognormal -- log scale
Wang2007l <- Wang2007
Wang2007l$DV <- log(Wang2007l$DV)

ini <- RxODE:::.foceiSetup(inner.add, Wang2007l, c(0.5, sqrt(.1)), rxInv=rxInv)

expect_equal(round(RxODE:::foceiOfv(ini), 3), round(6205.61885097915, 3))
expect_equal(round(RxODE:::foceiEtas(), 3), structure(list(ID = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                                           ETA1 = c(1.122, 0.997, 1.019, 1.028, 0.987, 0.996, 0.981, 0.998, 0.989, 0.991),
                                                           OBJI = c(618.183, 613.574, 622.858, 632.166, 617.722, 621.16, 620.271, 616.082, 619.258, 624.345)),
                                                      row.names = c(NA, -10L), class = "data.frame"))

rxSolveFree();


context("Focei Inner");

d <- as.data.frame(list(time = c(72, 95.99, 96, 119.99, 120, 143.99, 144,
                                     144.25, 144.5, 144.75, 145, 145.5, 146, 146.5, 147, 148, 150,
                                     152, 156, 160, 164, 167.99, 168, 191.99, 192, 215.99, 216, 216.25,
                                     216.5, 216.75, 217, 217.5, 218, 218.5, 219, 220, 222, 224, 228,
                                     232, 236, 240, 252, 264, 276, 288),
                            evid = c(101, 0, 101, 0,
                                     101, 0, 101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 101,
                                     0, 101, 0, 101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0),
                            amt = c(60000, NA, 60000, NA, 60000, NA, 60000,
                                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 60000,
                                    NA, 60000, NA, 60000, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                    NA, NA, NA, NA, NA, NA, NA, NA, NA),
                            dv = c(NA, 263.6, NA, 164.7,
                                   NA, 287.3, NA, 1248.7, 1211.5, 1017.7, 1690.1, 1029.8, 890.7,
                                   598.4, 1009.3, 1159.8, 742.2, 724.6, 728.2, 509.7, 243.1, 259.9,
                                   NA, 242.2, NA, 281.4, NA, 1500.1, 1281.4, 1200.2, 1378.8, 1373.2,
                                   582.9, 960.2, 720.3, 852.6, 950.3, 654.7, 402.5, 456, 346.5,
                                   268.2, 134.2, 42.6, 25.9, 14.6)))
    d$id <- 1;


inner <- RxODE("d/dt(centr)~-centr*exp(ETA[1]+THETA[1])*exp(-ETA[2]-THETA[2]);\nd/dt(rx__sens_centr_BY_ETA_1___)~-centr*exp(ETA[1]+THETA[1])*exp(-ETA[2]-THETA[2])-rx__sens_centr_BY_ETA_1___*exp(ETA[1]+THETA[1])*exp(-ETA[2]-THETA[2]);\nd/dt(rx__sens_centr_BY_ETA_2___)~centr*exp(ETA[1]+THETA[1])*exp(-ETA[2]-THETA[2])-rx__sens_centr_BY_ETA_2___*exp(ETA[1]+THETA[1])*exp(-ETA[2]-THETA[2]);\nrx_pred_=centr*exp(-ETA[2]-THETA[2]);\nrx__sens_rx_pred__BY_ETA_1___=rx__sens_centr_BY_ETA_1___*exp(-ETA[2]-THETA[2]);\nrx__sens_rx_pred__BY_ETA_2___=-centr*exp(-ETA[2]-THETA[2])+rx__sens_centr_BY_ETA_2___*exp(-ETA[2]-THETA[2]);\nrx_r_=Rx_pow_di(THETA[3],2)*Rx_pow_di(centr,2)*exp(-2*ETA[2]-2*THETA[2]);\nrx__sens_rx_r__BY_ETA_1___=2*Rx_pow_di(THETA[3],2)*centr*rx__sens_centr_BY_ETA_1___*exp(-2*ETA[2]-2*THETA[2]);\nrx__sens_rx_r__BY_ETA_2___=-2*Rx_pow_di(THETA[3],2)*Rx_pow_di(centr,2)*exp(-2*ETA[2]-2*THETA[2])+2*Rx_pow_di(THETA[3],2)*centr*rx__sens_centr_BY_ETA_2___*exp(-2*ETA[2]-2*THETA[2]);\n");

    theta <- structure(c(1.6, 4.5, sqrt(0.1)), .Names=sprintf("THETA[%d]", 1:3));
    eta <- structure(c(-0.147736086922763, -0.294637022436797), .Names=sprintf("ETA[%d]", 1:2))

    rxInv <- rxSymInvCholCreate(matrix(c(0.1, 0, 0, 0.1), nrow=2))

    etaMat <- t(matrix(eta))

    ## rxSolve setupOnly needs to be called first.
    ini <- RxODE:::.foceiSetup(inner, d, theta, rxInv = rxInv, etaMat = etaMat,
                               odeOpt=foceiControl(maxOuterIterations=0, maxInnerIterations=0));

    testthat::expect_equal(RxODE:::foceiLik(ini), -209.467627968204)


testthat::expect_equal(RxODE:::foceiInnerLp(rep(0, 2)), c(31.4177868473072, 56.8293635831368))

    RxODE:::likInner(rep(0, 2))

testthat::expect_equal(RxODE:::foceiInnerLp(eta), c(-0.000367366009438985, -2.46513261989989e-05))

    RxODE:::likInner(c(-0.000367366009438985, -2.46513261989989e-05))

    RxODE:::likInner(c(-0.1477353, -0.2946360))

    RxODE:::foceiInnerLp(c(-0.1477353, -0.2946360))

    ## print(numDeriv::grad(RxODE:::foceiLik, ini))

    ## print(RxODE:::foceiNumericGrad(ini))

rxSolveFree();

theta <- structure(c(1.6, 4.5, sqrt(0.1)), .Names=sprintf("THETA[%d]", 1:3));
eta <- structure(rep(0, 2), .Names=sprintf("ETA[%d]", 1:2))
etaMat <- t(matrix(eta))

ini <- RxODE:::.foceiSetup(inner, d, theta, rxInv = rxInv, odeOpt=foceiControl(printInner=50L, epsilon=1e-8))

print(RxODE:::foceiLik(ini))

RxODE:::foceiEtas()


## print(numDeriv::grad(RxODE:::foceiLik, rep(1, 5)))

## print(numDeriv::grad(RxODE:::foceiLik, rep(1, 5), method="simple"))

## print(numDeriv::grad(RxODE:::foceiLik, rep(1, 5)))

## print(RxODE:::foceiNumericGrad(rep(1, 5)))
rxSolveFree();

## inner.pred.prop <- RxODE("d/dt(ipre)~-ipre*THETA[1];\nd/dt(rx__sens_ipre_BY_ETA_1___)~-THETA[1]*rx__sens_ipre_BY_ETA_1___;\nrx_pred_=ipre;\nrx__sens_rx_pred__BY_ETA_1___=rx__sens_ipre_BY_ETA_1___;\nrx_r_=Rx_pow_di(THETA[2],2)*Rx_pow_di(ipre,2);\nrx__sens_rx_r__BY_ETA_1___=2*Rx_pow_di(THETA[2],2)*ipre*rx__sens_ipre_BY_ETA_1___;\n")

## dat2 <- Wang2007[Wang2007$TIME == 0, ]
## dat2$EVID <- 101
## dat2$AMT <- 10;
## dat2 <- rbind(dat2, Wang2007)
## dat2 <- dat2[(order(dat2$ID, -dat2$EVID, dat2$TIME)), ]


## ini <- RxODE:::.foceiSetup(inner.pred.prop, dat2, c(0.5, sqrt(.1)), rxInv = rxInv)

## print(round(RxODE:::foceiLik(ini), 3))
## ## expect_equal(round(RxODE:::foceiLik(ini), 3), 39.458)

## print(RxODE:::foceiNumericGrad(ini))

## rxSolveFree();


#},  silent=TRUE)
