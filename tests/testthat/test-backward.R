context("Test Backward Compatability")

rxPermissive({
    demo <- structure(list(ID = c(1012, 1012, 1012, 1012, 1012, 1012, 1012, 1012), TIME = c(588.5, 600.5, 612.5, 624.5, 636.5, 648.5, 660.5, 672.5), DOSE = c(2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000), AMT = c(1000, 1000, 1000, 1000, 1000, 1000, 1000, 0), TAD = c(0, 0, 0, 0, 0, 0, 0, 12), CL = c(5.851496056, 5.851496056, 5.851496056, 5.851496056, 5.851496056, 5.851496056, 5.851496056, 5.851496056), V = c(49.3930186, 49.3930186, 49.3930186, 49.3930186, 49.3930186, 49.3930186, 49.3930186, 49.3930186), KA = c(3.320205555, 3.320205555, 3.320205555, 3.320205555, 3.320205555, 3.320205555, 3.320205555, 3.320205555), TAD4 = c(0, 12, 24, 36, 48, 60, 72, 84)), .Names = c("ID", "TIME", "DOSE", "AMT", "TAD", "CL", "V", "KA", "TAD4"), row.names = c(NA, -8L), class = c("data.frame"), sorted = c("ID", "TIME"));

    ode1KA <- "
d/dt(abs)    = -KA*abs;
d/dt(centr)  =  KA*abs-(CL/V)*centr;
C1=centr/V;
"

    StepSize=1
    Extension=6

    mod1KA <- RxODE(model=ode1KA)

    params <- demo[1,c("CL","V","KA")]

    ev<-eventTable()
    DOSi<-as.data.frame(demo[demo$AMT>0,])
    for (j in 1:length(DOSi$AMT)){
        dos<-DOSi[j,]
        ev$add.dosing(dose=as.numeric(dos$AMT),nbr.doses=1,dosing.to=1,rate=NULL,start.time=as.numeric(dos$TAD4))
    }

    timei<-demo$TIME
    minimum<-min(timei)
    maximum<-max(timei)+Extension
    times<-sort(unique(c(timei,seq(minimum,maximum,StepSize))))
    ev$add.sampling(times)

    x <- as.data.frame(mod1KA$run(params, ev))

    test_that("run works with a data frame param", {
        expect_equal(class(x), "data.frame")
    })


    ## test old solving.


    event.table <- ev$get.EventTable()
    modelVars <- mod1KA$get.modelVars()
    state_vars <- modelVars$state;
    neq <- length(state_vars);
    lhs_vars <- modelVars$lhs;
    nlhs <- length(lhs_vars);

    ntime <- dim(event.table)[1];
    ret <- rep(0.0, ntime * neq);
    lhs <- rep(0.0, ntime * nlhs);
    rc <- as.integer(0)

    inits <- rep(0.0, neq)

    cmpMgr <- mod1KA$cmpMgr;

    cmpMgr$dynLoad()

    atol <- 1e-8
    rtol <- 1e-6
    transit_abs <- FALSE
    stiff <- TRUE

    xx <- .C(cmpMgr$ode_solver,
             as.integer(neq),
             as.double(params),
             as.double(event.table$time),
             as.integer(event.table$evid),
             length(event.table$time),
             as.double(inits),
             as.double(event.table$amt[event.table$evid > 0]),
             as.double(ret),
             as.double(atol),
             as.double(rtol),
             as.integer(stiff),
             as.integer(transit_abs),
             as.integer(nlhs),
             as.double(lhs),
             rc);


    x2 <- cbind(matrix(xx[[8]], ncol=neq, byrow=T),
                matrix(xx[[14]], ncol = nlhs, byrow = TRUE))
    colnames(x2) <- c(state_vars, lhs_vars);

    x2 <- cbind(time=event.table$time, x2)[ev$get.obs.rec(), ];
    x2 <- as.data.frame(x2)

    test_that("Old routine works correctly", {
        expect_equal(x2, x);
    })
}, silent=TRUE)
