rxPermissive({

    op <- options()
    ## To keep basename the same
    options("RxODE.basename.print"="basename",
            "RxODE.dll.print"="dll",
            "RxODE.c.print"="cfile")

    ## test_path("test-print.txt")
    verify_output(tempfile(), {

        mod <- RxODE({
            a = 6
            b = 0.6 + a / 100
            kel = b * 0.01
            V = 10
            d/dt(intestine) = -a*intestine
            d/dt(blood)     = a*intestine - b*blood
            l2 <- linCmt()
        })

        print(mod)

        summary(mod)
        str(mod) ## fragile

        print(mod$cmpMgr$rxDll())

        coef(mod)

        print(rxModelVars(mod))

        print(rxC(mod))
        summary(rxC(mod)) # too fragile

        print(mod$.rxDll)

        et <- eventTable(time.units="days")
        et$add.sampling(seq(0,10,by=1/24))
        et$add.dosing(dose=2/24,rate=2,start.time=0,
                      nbr.doses=10,dosing.interval=1)

        print(et)
        summary(et)
        str(et)
        print(attr(class(et), ".RxODE.lst"))
        str(attr(class(et), ".RxODE.lst"))

        print(format(structure(0:7,class="rxEvid")))
        print(structure(0:7,class="rxEvid"))

        print(format(structure(c(-2, -1, 0, 1, 2),class="rxRateDur")))
        print(structure(c(-2, -1, 0, 1, 2),class="rxRateDur"))

        pk <- solve(mod,et);
        print(pk)

        print(pk, width=40)

        print(pk, bound="k")

        print(pk, n=10)

        options(RxODE.display.tbl=FALSE)
        print(pk)

        options(RxODE.display.tbl=TRUE)

        summary(pk)
        str(pk)

        print(summary(pk, bound="k"))




        ## Now "destroy" the RxODE solved object and change the printout
        tmp <- class(pk)
        pk$matt <- 4
        class(tmp) <- tmp
        print(pk)


        rxUnload(mod)
        print(mod)
        summary(mod)
        str(mod)

        rxDelete(mod)
        print(mod)
        summary(mod)
        str(mod)

        mod <- RxODE(mod, indLin=TRUE)

        print(mod)
        summary(mod)
        ## str(mod)

        rxUnload(mod)
        print(mod)
        summary(mod)
        str(mod)

        rxDelete(mod)
        print(mod)
        summary(mod)
        str(mod)

        set.seed(42)
        tmp <- matrix(rnorm(5^2), 5, 5)
        mcov <- tcrossprod(tmp, tmp)
        v <- rxSymInvCholCreate(mcov, "sqrt")

        print(v)
        str(v)

        class(v) <- NULL
        v$env$theta <- NULL
        class(v) <- "rxSymInvCholEnv"
        print(v)

        mod <- RxODE("
a = 6
b = 0.6
d/dt(intestine) = -a*intestine
d/dt(blood)     = a*intestine - b*blood
")


    et <- eventTable()
    et$add.dosing(dose=2/24,rate=2,start.time=0,
                  nbr.doses=10,dosing.interval=1)
    et <- et %>% et(0.05,evid=2) %>%
        et(amt=3,time=0.5,cmt=out) %>%
        et(amt=3,time=0.1,cmt=intestine,ss=1,ii=3) %>%
        et(amt=3,time=0.3,cmt=intestine,ss=2,ii=3) %>%
        et(time=0.2,cmt="-intestine") %>%
        as.data.frame

    ett1 <- RxODE:::etTrans(et, mod, keepDosingOnly=TRUE)

    print(ett1)


})

    options(op)

}, test="print")
