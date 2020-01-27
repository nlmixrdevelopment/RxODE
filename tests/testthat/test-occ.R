rxPermissive({
    context("Occasion tests")

    mod <- RxODE({
        eff(0) = 1
        C2 = centr/V2*(1+prop.err);
        C3 = peri/V3;
        CL =  TCl*exp(eta.Cl + iov.Cl)
        KA = TKA * exp(eta.Ka + iov.Ka)
        d/dt(depot) =-KA*depot;
        d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  =                    Q*C2 - Q*C3;
        d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
    })

    mod.eta <- RxODE({
        eff(0) = 1
        C2 = centr/V2*(1+prop.err);
        C3 = peri/V3;
        CL =  TCl*exp(ETA[1] + iov.Cl)
        KA = TKA * exp(ETA[2] + iov.Ka)
        d/dt(depot) =-KA*depot;
        d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  =                    Q*C2 - Q*C3;
        d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
    })

    test_that("Expanding IOV", {
        expect_error(rxExpandOcc(mod,1,"iov.Cl"))
        expect_equal(rxExpandOcc(mod,2,c("iov.Cl", "iov.Ka")),
                     "iov.Ka=(rxOCC==2)*ETA[4]+(rxOCC==1)*ETA[2];\niov.Cl=(rxOCC==2)*ETA[3]+(rxOCC==1)*ETA[1];\neff(0)=1;\nC2=centr/V2*(1+prop.err);\nC3=peri/V3;\nCL=TCl*exp(eta.Cl+iov.Cl);\nKA=TKA*exp(eta.Ka+iov.Ka);\nd/dt(depot)=-KA*depot;\nd/dt(centr)=KA*depot-CL*C2-Q*C2+Q*C3;\nd/dt(peri)=Q*C2-Q*C3;\nd/dt(eff)=Kin-Kout*(1-C2/(EC50+C2))*eff;\n")

        expect_equal(rxExpandOcc(mod,2,c("iov.Ka", "iov.Cl")),
                     "iov.Cl=(rxOCC==2)*ETA[4]+(rxOCC==1)*ETA[2];\niov.Ka=(rxOCC==2)*ETA[3]+(rxOCC==1)*ETA[1];\neff(0)=1;\nC2=centr/V2*(1+prop.err);\nC3=peri/V3;\nCL=TCl*exp(eta.Cl+iov.Cl);\nKA=TKA*exp(eta.Ka+iov.Ka);\nd/dt(depot)=-KA*depot;\nd/dt(centr)=KA*depot-CL*C2-Q*C2+Q*C3;\nd/dt(peri)=Q*C2-Q*C3;\nd/dt(eff)=Kin-Kout*(1-C2/(EC50+C2))*eff;\n")

        expect_equal(rxExpandOcc(mod.eta,3,c("iov.Cl", "iov.Ka")),
                     "iov.Ka=(rxOCC==3)*ETA[8]+(rxOCC==2)*ETA[6]+(rxOCC==1)*ETA[4];\niov.Cl=(rxOCC==3)*ETA[7]+(rxOCC==2)*ETA[5]+(rxOCC==1)*ETA[3];\neff(0)=1;\nC2=centr/V2*(1+prop.err);\nC3=peri/V3;\nCL=TCl*exp(ETA[1]+iov.Cl);\nKA=TKA*exp(ETA[2]+iov.Ka);\nd/dt(depot)=-KA*depot;\nd/dt(centr)=KA*depot-CL*C2-Q*C2+Q*C3;\nd/dt(peri)=Q*C2-Q*C3;\nd/dt(eff)=Kin-Kout*(1-C2/(EC50+C2))*eff;\n")

    })

})
