rxPermissive({
    context("Compartment order & extra CMTs with cmt()");
    test_that("cmt() syntax makes sense", {
        mod <- RxODE({
            a  <-  6
            b  <-  0.6
            cmt(blood); # cmt = 1 now
            d/dt(intestine) <- -a*intestine
            d/dt(blood) <- a*intestine - b*blood
        })

        expect_equal(c("blood", "intestine"), rxState(mod));

        mod <- RxODE({
            a  <-  6
            b  <-  0.6
            d/dt(intestine) <- -a*intestine
            d/dt(blood) <- a*intestine - b*blood
        })

        expect_equal(c("intestine", "blood"), rxState(mod));

        expect_error(RxODE({
            a  <-  6
            b  <-  0.6
            cmt(matt); # cmt = 1 now
            d/dt(intestine) <- -a*intestine
            d/dt(blood) <- a*intestine - b*blood
        }))

        tmp <- RxODE({
            a  <-  6
            b  <-  0.6
            cmt(blood); # cmt = 1 now
            cmt(intestine); # cmt = 2 now
            cmt(matt)
            d/dt(intestine) <- -a*intestine
            d/dt(blood) <- a*intestine - b*blood
        })
        expect_equal(tmp$stateExtra,"matt")

        context("Compartment melding with dvid");

        w  <- RxODE({
            ktr <- exp(tktr + eta.ktr)
            ka <- exp(tka + eta.ka)
            cl <- exp(tcl + eta.cl)
            ##
            logit=exp(poplogit+eta.emax)
            emax = logit/(1+logit)
            ec50 =  exp(tec50 + eta.ec50)
            kout = exp(tkout + eta.kout)
            e0 = exp(te0 + eta.e0)
            ##
            DCP = center/v
            PD=1-emax*DCP/(ec50+DCP)
            ##
            effect(0) = e0
            kin = e0*kout
            ##
            d/dt(depot) = -ktr * depot
            d/dt(gut) =  ktr * depot -ka * gut
            d/dt(center) =  ka * gut - cl / v * center
            d/dt(effect) = kin*PD -kout*effect
            ##
            cp = center / v
        })

        tmp  <- RxODE:::etTrans(warfarin, w,addCmt=TRUE)

        w2  <- warfarin
        w2$dvid  <- as.integer(w2$dvid)

        w2  <- warfarin
        w2$dvid  <- paste(w2$dvid)
        tmp  <- RxODE:::etTrans(w2, w,addCmt=TRUE)

        expect_equal(tmp$CMT[tmp$ID==1],
                     c(1L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 5L, 6L, 5L, 6L, 5L,  6L, 6L, 6L, 6L))

        tmp  <- RxODE:::etTrans(warfarin, w,addCmt=TRUE)

        expect_equal(tmp$CMT[tmp$ID==1],
                     c(1L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 5L, 6L, 5L, 6L, 5L,  6L, 6L, 6L, 6L))

        tmp  <- RxODE({
            d/dt(depot) = -ka * depot
            d/dt(center) = ka * depot - cl/v * center
            cp = center/v
            cp2 = cp + nlmixr_extra_par
            cmt(central);
            cmt(c20);
        })

        expect_equal(class(tmp), "RxODE")

        context("Check lhs allowed stateExtra while preserving lhs properties.")

        tmp  <- RxODE({
            d/dt(depot) = -ka * depot
            d/dt(center) = ka * depot - cl/v * center
            cp = center/v;
            cmt(cp);
        })

        expect_equal(class(tmp), "RxODE")
        expect_equal(tmp$lhs, "cp")
        expect_equal(tmp$stateExtra, "cp")

        tmp  <- RxODE({
            d/dt(depot) = -ka * depot
            d/dt(center) = ka * depot - cl/v * center
            cmt(cp);
            cp = center/v;
        })

        expect_equal(class(tmp), "RxODE")
        expect_equal(tmp$lhs, "cp")
        expect_equal(tmp$stateExtra, "cp")

        tmp  <- RxODE({
            d/dt(depot) = -ka * depot
            d/dt(center) = ka * depot - cl/v * center
            cmt(cp);
            cp = 3;
        })

        expect_equal(class(tmp), "RxODE")
        expect_equal(tmp$lhs, "cp")
        expect_equal(tmp$stateExtra, "cp")

        tmp  <- RxODE({
            d/dt(depot) = -ka * depot
            d/dt(center) = ka * depot - cl/v * center
            cp = 3;
            cmt(cp);
        })

        expect_equal(class(tmp), "RxODE")
        expect_equal(tmp$lhs, "cp")
        expect_equal(tmp$stateExtra, "cp")

    })
},cran=TRUE)




