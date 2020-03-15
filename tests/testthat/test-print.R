rxPermissive({

    op <- options()
    ## To keep basename the same
    options("RxODE.basename.print"="basename",
            "RxODE.dll.print"="dll")

    ## test_path("test-print.txt")
    verify_output(tempfile(), {

        mod <- RxODE({
            a = 6
            b = 0.6
            d/dt(intestine) = -a*intestine
            d/dt(blood)     = a*intestine - b*blood
        })

        print(mod)
        summary(mod)
        str(mod)

        print(rxModelVars(mod))

        print(rxC(mod))
        summary(rxC(mod))

        print(rxDll(mod))

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
        summary(pk)
        str(pk)

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
        str(mod)

        rxUnload(mod)
        print(mod)
        summary(mod)
        str(mod)

        rxDelete(mod)
        print(mod)
        summary(mod)
        str(mod)
    })

    options(op)

}, test="print")
