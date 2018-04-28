## test ODE parsing for syntax errors
library("RxODE")
context("Test Parsing of models")
rxPermissive({

    badParse <- function(desc,code){
        test_that(desc,{
            tmp <- tempfile();
            sink(tmp)
            on.exit({sink(); unlink(tmp);});
            expect_error(RxODE(code));
        })
    }

    goodParse <- function(desc,code){
        test_that(desc,{
            tmp <- tempfile();
            sink(tmp)
            on.exit({sink(); unlink(tmp);});
            rx <- RxODE(code);
            expect_equal(class(rx),"RxODE");
            rxDelete(rx);
        })
    }

    equivSyntax <- function(desc,code1,code2){
        test_that(desc,{
            tmp <- tempfile();
            sink(tmp);
            on.exit({sink(); unlink(tmp);});
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
    goodParse('comments must be outside statements','d/dt(y) = 1   # bad comment;')
    goodParse('missing end of statement ";" dosen\'t cause errors',
              paste(sep = "\n",
                    'd/dt(depot) = -ka * depot',
                    'd/dt(centr) =  ka * depot - kout * centr;'))

    options(RxODE.syntax.require.semicolon=TRUE);rxSyncOptions();
    badParse('comments must be outside statements','d/dt(y) = 1   # bad comment;')
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

    goodParse(desc = 'Assignment with <- supported #5',
              'd/dt(y_1) = F*y; jac(y_1,y_1)<-0')

    goodParse(desc = 'Assignment with <- supported #6',
              'd/dt(y_1) = F*y; df/dy(y_1,y_1)<-0')

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

    badParse(desc = 'Assignment with <- not supported #5',
             'd/dt(y_1) = F*y; jac(y_1,y_1)<-y')

    badParse(desc = 'Assignment with <- not supported #6',
             'd/dt(y_1) = F*y; df/dy(y_1,y_1)<-y')

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
jac(mu,y)=0;
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

}, silent=TRUE);
