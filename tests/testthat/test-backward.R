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

    ode_solve <- cmpMgr$ode_solver;
    xx <- mod1KA$dll$.c(ode_solve,
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

    ## Dynmodel routines
    ode <- "
   dose=200;
   pi = 3.1415926535897931;
   if (t<=0) {
      fI = 0;
   } else {
      fI = F*dose*sqrt(MIT/(2.0*pi*CVI2*t^3))*exp(-(t-MIT)^2/(2.0*CVI2*MIT*t));
   }
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(centr) = fI - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =              Q*C2 - Q*C3;
"

    sys1 <- RxODE(model = ode)

    ev <- eventTable()
    ev$add.sampling(c(0, c(15, 30, 60, 90, 120, 150, 210, 270, 330, 360, 390, 420, 450, 480)))
    theta <- structure(c(190, 0.65, 0.92, 0.0793, 0.64, 0.292, 9.63), .Names = c("MIT",  "CVI2", "F", "CL", "V2", "Q", "V3"))

    val1 <- sys1$solve(theta, ev, atol=1e-6, rtol=1e-6)

    ## Prior RxODE solving...
    val2 <- structure(c(0, 15, 30, 60, 90, 120, 150, 210, 270, 330, 360, 390, 420, 450, 480, 0, 0.00468853094648131, 0.373067597127217, 2.06008521781962, 3.11868121404469, 3.61478640740494, 3.77791523493719, 3.57583622912201, 3.07240077291895, 2.50781365737478, 2.23759663817273, 1.98421582221308, 1.75065831378882, 1.53813326790626, 1.34666851432174, 0, 0.00275724533602532, 0.779135520764223, 12.5294880193912, 29.5561651547891, 42.805783957353, 50.756288454544, 54.6956530842046, 50.0020491493438, 42.249424648696, 38.1457972477322, 34.1486661506102, 30.3627294500035, 26.8473328340894, 23.6307416646032, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 3.14159265358979, 0, 0.00555514066938253, 0.241306787665833, 0.863322464865012, 0.937389549382366, 0.809244619745276, 0.654268731228509, 0.409232361102111, 0.256995004576163, 0.164606161831775, 0.132752793640541, 0.107575895938534, 0.0875618379140851, 0.0715643817688998, 0.0587112500960777, 0, 0.00732582960387705, 0.582918120511276, 3.21888315284315, 4.87293939694483, 5.64810376157022, 5.90299255458936, 5.58724410800314, 4.80062620768586, 3.9184588396481, 3.49624474714489, 3.10033722220793, 2.73540361529503, 2.40333323110352, 2.10416955362771, 0, 0.00028631831111374, 0.0809071153441561, 1.30108909858683, 3.06917602853469, 4.44504506306884, 5.27064262248639, 5.6797147543307, 5.19232078394016, 4.38727151076802, 3.96114197795765, 3.54607125136139, 3.15293140706163, 2.78788502950045, 2.45386725489129), .Dim = c(15L, 8L), .Dimnames = list(NULL, c("time", "centr", "peri", "dose", "pi", "fI", "C2", "C3")));

    val2 <- val2[, dimnames(val1)[[2]]]

    test_that("pi values will not corrupt solver", {
        expect_equal(round(val1, 3), round(val2, 3))
    })

    ## Now test without pi

    ode <- "
   dose=200;
   if (t<=0) {
      fI = 0;
   } else {
      fI = F*dose*sqrt(MIT/(2.0*pi*CVI2*t^3))*exp(-(t-MIT)^2/(2.0*CVI2*MIT*t));
   }
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(centr) = fI - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =              Q*C2 - Q*C3;
"
    sys1 <- RxODE(model = ode)

    val3 <- sys1$solve(theta, ev, atol=1e-6, rtol=1e-6)

    test_that("pi is not needed...", {
        expect_equal(round(val1, 3), round(val3, 3))
    })

}, silent=TRUE)
