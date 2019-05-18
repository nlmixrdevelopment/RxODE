## test ODE parsing for syntax errors
library("RxODE")
context("Test Parsing of models")
rxPermissive({

    badParse <- function(desc,code){
        test_that(desc,{
            tmp <- normalizePath(tempfile(), mustWork = FALSE);
            zz <- file(tmp, open = "wt")
            sink(zz)
            sink(zz, type = "message")
            on.exit({sink(); sink(type = "message");close(zz);unlink(tmp);});
            expect_error(RxODE(code));
        })
    }

    goodParse <- function(desc,code){
        test_that(desc,{
            tmp <- normalizePath(tempfile(), mustWork = FALSE);
            zz <- file(tmp, open = "wt")
            sink(zz)
            sink(zz, type = "message")
            on.exit({sink(); sink(type = "message");close(zz);unlink(tmp);});
            rx <- RxODE(code);
            expect_equal(class(rx),"RxODE");
            rxDelete(rx);
        })
    }

    equivSyntax <- function(desc,code1,code2){
        test_that(desc,{
            tmp <- normalizePath(tempfile(), mustWork = FALSE);
            zz <- file(tmp, open = "wt")
            sink(zz)
            sink(zz, type = "message")
            on.exit({sink(); sink(type = "message");close(zz);unlink(tmp);});
            rx1 <- RxODE(code1);
            rx2 <- RxODE(code2);
            expect_equal(rxMd5(rx1)["parsed_md5"],rxMd5(rx2)["parsed_md5"])
            rxDelete(rx1);
            rxDelete(rx2);
        })
    }

    badParse('incorrect d/dt operator','d/dt(y = 1);')

    ## Statements don't require ; now.
    options(RxODE.syntax.require.semicolon=FALSE);rxSyncOptions();
    goodParse('comments must be outside statements #1','d/dt(y) = 1   # bad comment;')
    goodParse('missing end of statement ";" dosen\'t cause errors',
              paste(sep = "\n",
                    'd/dt(depot) = -ka * depot',
                    'd/dt(centr) =  ka * depot - kout * centr;'))

    options(RxODE.syntax.require.semicolon=TRUE);rxSyncOptions();
    badParse('comments must be outside statements #2','d/dt(y) = 1   # bad comment;')
    badParse('missing end of statement ";"',
             paste(sep = "\n",
                   'd/dt(depot) = -ka * depot',
                   'd/dt(centr) =  ka * depot - kout * centr;'))

    options(RxODE.syntax.require.semicolon=FALSE);rxSyncOptions();
    badParse('arithmetic syntax error',
             paste(sep = "\n",
                   '# comment, just to show error in line 3',
                   'd/dt(y) = -ka;',
                   'C1 = /y;'))
    ## added ** operator
    options(RxODE.syntax.star.pow=TRUE);rxSyncOptions();
    goodParse('existing operator **',
              code = paste(sep = "\n",
                           'd/dt(y) = -ka;',
                           'C1 = ka *  y**2;'))
    options(RxODE.syntax.star.pow=FALSE);rxSyncOptions();
    badParse('existing operator **',
             code = paste(sep = "\n",
                          'd/dt(y) = -ka;',
                          'C1 = ka *  y**2;'))


    badParse('unexistent operator %',
             code = paste(sep = "\n",
                          'remainder = 4 % 3;',
                          'd/dt(y) = -ka;',
                          'C1 = ka * y;')
             );

    badParse(desc = 'incorrect "if" statement',
             code = paste(sep = "\n",
                          'if(comed==0){',
                          '   F = 1.0;',
                          'else {',         # missing "}"'
                          '   F = 0.75;',
                          '};',
                          'd/dt(y) = F * y;')
             )

    badParse(desc = 'illegal variable name (starting w. a digit)',
             code = paste(sep = "\n",
                          'F = 0.75;',
                          '12foo_bar = 1.0/2.0;',
                          'd/dt(y) = F * y;')
             )

    options(RxODE.syntax.allow.dots=TRUE);rxSyncOptions();

    goodParse(desc = 'dot in variable name (ini0)',
              code = paste(sep = "\n",
                           'F = 0.75;',
                           'foo.bar = 1.0/2.0;',
                           'd/dt(y) = F * y;')
              )

    goodParse(desc = 'dot in variable name in d/dt()',
              code = paste(sep = "\n",
                           'd/dt(y_1) = F * y;',   # okay
                           'd/dt(y.1) = F * y;')   # not okay
              )

    goodParse(desc = 'leading dot in variable name',
              code = paste(sep = "\n",
                           'F = 0.75;',
                           '.foo.bar = 0.5;',
                           'd/dt(y) = F * y;')
              )

    goodParse(desc = 'leading dot in variable name (ini0)',
              code = paste(sep = "\n",
                           'F = 0.75;',
                           '.foo.bar = 1.0/2.0;',
                           'd/dt(y) = F * y;')
              )

    goodParse(desc = 'leading dot in variable name in d/dt()',
              code = paste(sep = "\n",
                           'd/dt(y_1) = F * y;',   # okay
                           'd/dt(.y.1) = F * y;')   # not okay
              )

    goodParse(desc = 'leading dot in variable name',
              code = paste(sep = "\n",
                           'F = 0.75;',
                           '.foo.bar = 0.5;',
                           'd/dt(y) = F * y;')
              )

    options(RxODE.syntax.allow.dots=FALSE); rxSyncOptions();
    badParse(desc = 'dot in variable name (ini0)',
             code = paste(sep = "\n",
                          'F = 0.75;',
                          'foo.bar = 1.0/2.0;',
                          'd/dt(y) = F * y;')
             )

    badParse(desc = 'dot in variable name in d/dt()',
             code = paste(sep = "\n",
                          'd/dt(y_1) = F * y;',   # okay
                          'd/dt(y.1) = F * y;')   # not okay
             )

    badParse(desc = 'leading dot in variable name',
             code = paste(sep = "\n",
                          'F = 0.75;',
                          '.foo.bar = 0.5;',
                          'd/dt(y) = F * y;')
             )

    badParse(desc = 'leading dot in variable name (ini0)',
             code = paste(sep = "\n",
                          'F = 0.75;',
                          '.foo.bar = 1.0/2.0;',
                          'd/dt(y) = F * y;')
             )

    badParse(desc = 'leading dot in variable name in d/dt()',
             code = paste(sep = "\n",
                          'd/dt(y_1) = F * y;',   # okay
                          'd/dt(.y.1) = F * y;')   # not okay
             )

    badParse(desc = 'leading dot in variable name',
             code = paste(sep = "\n",
                          'F = 0.75;',
                          '.foo.bar = 0.5;',
                          'd/dt(y) = F * y;')
             )

    badParse(desc = 'Assignment with <<- not supported',
             'd/dt(y_1) <<- F*y')

    options(RxODE.syntax.assign=TRUE); rxSyncOptions();
    goodParse(desc = 'Assignment with <- supported #1',
              'd/dt(y_1) <- F*y')

    goodParse(desc = 'Assignment with <- supported #2',
              'y_1(0) <- 1;d/dt(y_1) = F*y_1')

    goodParse(desc = 'Assignment with <- supported #3',
              'y_2 <- 1;d/dt(y_1) = F*y')

    goodParse(desc = 'Assignment with <- supported #4',
              'y_2 <- 1+7;d/dt(y_1) = F*y')

    goodParse(desc = 'Assignment with <- supported #7',
              'd/dt(y_1) = F*y; df(y_1)/dy(y_1) <- 0')

    options(RxODE.syntax.assign=FALSE); rxSyncOptions();
    badParse(desc = 'Assignment with <- not supported #1',
             'd/dt(y_1) <- F*y')

    badParse(desc = 'Assignment with <- not supported #2',
             'y_1(0) <- 1;d/dt(y_1) = F*y')

    badParse(desc = 'Assignment with <- not supported #3',
             'y_2 <- 1;d/dt(y_1) = F*y')

    badParse(desc = 'Assignment with <- not supported #4',
             'y_2 <- 1+7;d/dt(y_1) = F*y')

    badParse(desc = 'Assignment with <- not supported #7',
             'd/dt(y_1) = F*y; df(y_1)/dy(F) <- y')
    equivSyntax(desc = "time and t are equivalent",
                "d/dt(depot) = time^2","d/dt(depot) = t^2")

    badParse(desc = "Assignment of state varaible #1","
a       = 1.0E4+0
x       = 1
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
")

    badParse(desc = "Assignment of state varaible #2","
a       = 1.0E4+0
x       = 1+2
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
")

    goodParse(desc = "Assignment of state varaible #3","
a       = 1.0E4+0
x(0)    = 1+2
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
")

    badParse(desc = "Assignment of state varaible #4","
a       = 1.0E4+0
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
x       = 1
")


    badParse(desc = "Assignment of state varaible #5","
a       = 1.0E4+0
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
x       = 1+2
")

    goodParse(desc = "Initial Condition Assignment #1","
a       = 1.0E4+0
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
x(0)     = 1
")

    goodParse(desc = "Initial Condition Assignment #2","
a       = 1.0E4+0
x(0)    = 1
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
")

    badParse(desc = "Initial conditions based on if/then statements","
if (a > 1){
x(0)    = 1
} else {
x(0) = 2
}
a       = 1.0E4+0
x(0)    = 1
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
")


    badParse(desc = "RHS d/dt(x) before defined","
d/dt
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
")
    options(RxODE.syntax.allow.ini0=FALSE); rxSyncOptions();
    badParse(desc = 'y_1(0) unsupported when RxODE.syntax.allow.ini0=FALSE',
              'y_1(0) = 1;d/dt(y_1) = F*y_1')
    options(RxODE.syntax.allow.ini0=TRUE); rxSyncOptions();

    badParse(desc="Defining df(var1)/dy(var2) where var1 is not a state variable.",
             "
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Jacobian
df(y)/dy(dy)  = 1
df(dy)/dy(y)  = -2*dy*mu*y - 1
df(dy)/dy(dy) = mu*(1-y^2)
df(mu)/dy(y)=0;
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
")

    badParse(desc="Defining df(var1)/dy(var2) where var1 is not a state variable.",
             "
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Jacobian
df(y)/dy(dy)  = 1
df(dy)/dy(y)  = -2*dy*mu*y - 1
df(dy)/dy(dy) = mu*(1-y^2)
df(mu)/dy(y)=0;
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
")

    goodParse(desc="Defining df(var1)/dy(var2) where var1 is a state variable.",
          "
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Jacobian
df(y)/dy(dy)  = 1
df(dy)/dy(y)  = -2*dy*mu*y - 1
df(dy)/dy(dy) = mu*(1-y^2)
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
")

    goodParse(desc="Defining df(var1)/dy(var2) where var2 is a variable.",
              "
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Jacobian
df(y)/dy(dy)  = 1
df(dy)/dy(y)  = -2*dy*mu*y - 1
df(dy)/dy(dy) = mu*(1-y^2)
df(dy)/dy(mu) = (1-y^2)*dy
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
")

    badParse(desc="Defining df(var1)/dy(var2) where var2 is a calculated value.",
              "
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Jacobian
df(y)/dy(dy)  = 1
df(dy)/dy(y)  = -2*dy*mu*y - 1
df(dy)/dy(dy) = mu*(1-y^2)
df(dy)/dy(mu) = (1-y^2)*dy
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1+bad ## nonstiff; 10 moderately stiff; 1000 stiff
")

    goodParse(desc="a*b/c^2",
              "d/dt(x)=a*b/c^2*x")

    goodParse(desc="a*b/c^2/d",
              "d/dt(x)=a*b/c^2/d*x")

    goodParse(desc="Transit as a compartment",
              "d/dt(transit) = -(1/mtt) * transit
        d/dt(depot)   =  (1/mtt) * transit - ka * depot
        d/dt(center)  = ka * depot - (cl/v1) * center - (q/v1) * center + (q/v2) * periph
        d/dt(periph)  = (q / v1) * center - (q / v2) * periph
        cp = center / v1")

    goodParse(desc="=+ parsing",
              "C2 = +centr/V2;
              C3 = peri/V3;
              d/dt(depot) =-KA*depot;
              d/dt(centr) = +KA*depot - CL*C2 - Q*C2 + Q*C3;
              d/dt(peri)  =                    Q*C2 - Q*C3;
              d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;")

    for (v in c("rate", "amt", "ss", "addl", "evid", "ii")){
        badParse(desc=sprintf("bad variables: %s",v),
                 sprintf("var=%s",v))
    }
    for (v in c("tlast")){
        goodParse(desc=sprintf("good variables: %s",v),
                 sprintf("var=%s",v))
    }
    for (v in c("abs")){
        goodParse(desc=sprintf("good functions: %s",v),
                 sprintf("var=%s(x)",v))
    }

    badParse(desc="No duplicate dvid()","a=b;dvid(1,2,3);dvid(3,4,5)");

    badParse(desc="Bad dvid(0)","a=b;dvid(0);");
    badParse(desc="Bad dvid(0, 1)","a=b;dvid(1,0);");

    goodParse(desc="THETA/ETA parsing","a=THETA[1]+ETA[2];");

    goodParse(desc="Duplicate d/dt(x)", "d/dt(depot) = -depot*ka;\nd/dt(depot) = d/dt(depot)+0")

    goodParse(desc="Duplicate d/dt(x)", "d/dt(depot) ~ -depot*ka;\nd/dt(depot) = d/dt(depot)+0")
    goodParse(desc="Duplicate d/dt(x)", "d/dt(depot) = -depot*ka;\nd/dt(depot) ~ d/dt(depot)+0")

    goodParse(desc="pi Parse", "a = pi+e2")

    for (v in c("f", "F", "alag", "lag", "rate", "dur")){
        badParse(sprintf("%s cannot depend on d/dt(state)", v),
                 sprintf("d/dt(depot)=-depot*ka;\nd/dt(central)=ka*depot-kel*central\n%s(depot)=d/dt(central)+3", v))
    }

    badParse(desc="Using d/dt(x) before defined", "y=d/dt(x)+3")

    goodParse(desc="Jacobain with theta and eta",
              "d/dt(x)=(THETA[1]+ETA[1])*x\ndf(x)/dy(THETA[1]) = 1\ndf(x)/dy(ETA[1]) = 1\n")


    for (v in c("f", "F", "alag", "lag", "rate", "dur")){
        badParse(sprintf("%s cannot depend on jacobain info", v),
                 sprintf("d/dt(x)=(THETA[1]+ETA[1])*x\ndf(x)/dy(THETA[1]) = 1\ndf(x)/dy(ETA[1]) = 1\n%s(x)=df(x)/dt(ETA[1])+3", v))
    }

    badParse("mtime cannot depend on df(x)/dt(ETA[1])",
             "d/dt(x)=(THETA[1]+ETA[1])*x\ndf(x)/dy(THETA[1]) = 1\ndf(x)/dy(ETA[1]) = 1\nmtime(z)=df(x)/dt(ETA[1])+3")
    badParse("mtime cannot depend on d/dt(x)",
             "d/dt(x)=(THETA[1]+ETA[1])*x\ndf(x)/dy(THETA[1]) = 1\ndf(x)/dy(ETA[1]) = 1\nmtime(z)=d/dt(x)+3")

    goodParse("functional initialization ok",
              "x(0) = y + 3\nd/dt(x) =-3*z")

    ## test_that("large model parse",{
    ##     skip_on_cran()
    ##     if (file.exists("model_default.txt")){
    ##         env  <- new.env();
    ##         m  <- microbenchmark::microbenchmark(assign("model", RxODE("model_default.txt"),env),times=1);
    ##         tmp  <- summary(m, "s")
    ##         context(sprintf("Large parse works (#61; %s s)", tmp$median))
    ##         expect_true(inherits(env$model, "RxODE"))
    ##     }
    ## })

    ## Bad parse:
    for (var in c('rate', 'dur', 'dvid')){
        badParse(sprintf("Cannot use the reserved %s variable", var),
                  sprintf("x = %s + 3"))
    }

    ## 'rate' and 'dur' can be data items, so they cannot be variables in an RxODE model
    for (var in c("lag", "alag", "f", "F")){
        goodParse(sprintf("Parsing of %s as a variable and function work.", var),
                  sprintf("d/dt(x) = -k*x;%s(x) = %s;", var, var))
    }

    options(RxODE.syntax.assign = TRUE)

    goodParse("mix lincmt with lags etc", "popCl <- 1
    popV <- 20
    popKa <- 1
    popVp <- 10
    popQ <- 2
    bsvCl <-0
    bsvV <- 0
    bsvKa <-0
    bsvVp <- 0
    bsvQ <-0
    popKeo <- 1.4
    bsvKeo <- 0
    popE0 <- 0
    popEmax <- 1
    popEC50 <- 5
    popGamma <- 1
    bsvE0 <- 0
    bsvEmax <- 0
    bsvEC50 <- 0
    ##
    cl ~ popCl * exp(bsvCl)
    v ~ popV * exp(bsvV)
    ka ~ popKa * exp(bsvKa)
    q ~ popQ * exp(bsvQ)
    vp ~ popVp * exp(bsvVp)
    keo ~ popKeo * exp(bsvKeo)
    popLagDepot <- 0
    popLagCentral <- 0
    popRateCentral <- 0
    popDurCentral <- 0
    bsvLagDepot <- 0
    bsvLagCentral <- 0
    bsvRateCentral <- 0
    bsvDurCentral <- 0
    lag(depot) <- popLagDepot * exp(bsvLagDepot)
    lag(central) <- popLagCentral * exp(bsvLagCentral)
    rate(central) <- popRateCentral * exp(bsvRateCentral)
    dur(central) <- popDurCentral * exp(bsvDurCentral)
    cp <- linCmt()
    d/dt(ce) = keo*(cp-ce)
    effect = E0 - Emax*(Ce^gamma)/((Ce^gamma)+(Ec50^gamma));")

    badParse("Still cannot take undefined compartments",
    "popCl <- 1
    popV <- 20
    popKa <- 1
    popVp <- 10
    popQ <- 2
    bsvCl <-0
    bsvV <- 0
    bsvKa <-0
    bsvVp <- 0
    bsvQ <-0
    popKeo <- 1.4
    bsvKeo <- 0
    popE0 <- 0
    popEmax <- 1
    popEC50 <- 5
    popGamma <- 1
    bsvE0 <- 0
    bsvEmax <- 0
    bsvEC50 <- 0
    ##
    cl ~ popCl * exp(bsvCl)
    v ~ popV * exp(bsvV)
    ka ~ popKa * exp(bsvKa)
    q ~ popQ * exp(bsvQ)
    vp ~ popVp * exp(bsvVp)
    keo ~ popKeo * exp(bsvKeo)
    popLagDepot <- 0
    popLagCentral <- 0
    popRateCentral <- 0
    popDurCentral <- 0
    bsvLagDepot <- 0
    bsvLagCentral <- 0
    bsvRateCentral <- 0
    bsvDurCentral <- 0
    lag(depot) <- popLagDepot * exp(bsvLagDepot)
    lag(central) <- popLagCentral * exp(bsvLagCentral)
    rate(central) <- popRateCentral * exp(bsvRateCentral)
    dur(central) <- popDurCentral * exp(bsvDurCentral)
    dur(matt) <- 3
    cp <- linCmt()
    d/dt(ce) = keo*(cp-ce)
    effect = E0 - Emax*(Ce^gamma)/((Ce^gamma)+(Ec50^gamma));")

    badParse("cmt(depot) doesn't work with linCmt()", "popCl <- 1
    cmt(depot)
    popV <- 20
    popKa <- 1
    popVp <- 10
    popQ <- 2
    bsvCl <-0
    bsvV <- 0
    bsvKa <-0
    bsvVp <- 0
    bsvQ <-0
    popKeo <- 1.4
    bsvKeo <- 0
    popE0 <- 0
    popEmax <- 1
    popEC50 <- 5
    popGamma <- 1
    bsvE0 <- 0
    bsvEmax <- 0
    bsvEC50 <- 0
    ##
    cl ~ popCl * exp(bsvCl)
    v ~ popV * exp(bsvV)
    ka ~ popKa * exp(bsvKa)
    q ~ popQ * exp(bsvQ)
    vp ~ popVp * exp(bsvVp)
    keo ~ popKeo * exp(bsvKeo)
    popLagDepot <- 0
    popLagCentral <- 0
    popRateCentral <- 0
    popDurCentral <- 0
    bsvLagDepot <- 0
    bsvLagCentral <- 0
    bsvRateCentral <- 0
    bsvDurCentral <- 0
    lag(depot) <- popLagDepot * exp(bsvLagDepot)
    lag(central) <- popLagCentral * exp(bsvLagCentral)
    rate(central) <- popRateCentral * exp(bsvRateCentral)
    dur(central) <- popDurCentral * exp(bsvDurCentral)
    cp <- linCmt()
    d/dt(ce) = keo*(cp-ce)
    effect = E0 - Emax*(Ce^gamma)/((Ce^gamma)+(Ec50^gamma));")

    badParse("cmt(central) doesn't work with linCmt()", "popCl <- 1
    cmt(central)
    popV <- 20
    popKa <- 1
    popVp <- 10
    popQ <- 2
    bsvCl <-0
    bsvV <- 0
    bsvKa <-0
    bsvVp <- 0
    bsvQ <-0
    popKeo <- 1.4
    bsvKeo <- 0
    popE0 <- 0
    popEmax <- 1
    popEC50 <- 5
    popGamma <- 1
    bsvE0 <- 0
    bsvEmax <- 0
    bsvEC50 <- 0
    ##
    cl ~ popCl * exp(bsvCl)
    v ~ popV * exp(bsvV)
    ka ~ popKa * exp(bsvKa)
    q ~ popQ * exp(bsvQ)
    vp ~ popVp * exp(bsvVp)
    keo ~ popKeo * exp(bsvKeo)
    popLagDepot <- 0
    popLagCentral <- 0
    popRateCentral <- 0
    popDurCentral <- 0
    bsvLagDepot <- 0
    bsvLagCentral <- 0
    bsvRateCentral <- 0
    bsvDurCentral <- 0
    lag(depot) <- popLagDepot * exp(bsvLagDepot)
    lag(central) <- popLagCentral * exp(bsvLagCentral)
    rate(central) <- popRateCentral * exp(bsvRateCentral)
    dur(central) <- popDurCentral * exp(bsvDurCentral)
    cp <- linCmt()
    d/dt(ce) = keo*(cp-ce)
    effect = E0 - Emax*(Ce^gamma)/((Ce^gamma)+(Ec50^gamma));")

}, silent=TRUE);
