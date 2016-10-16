library(RxODE)
library(dplyr)
library(digest)
context("Make sure solved rxSolve work with dplyr objects")

test.dir <- tempfile("Rxmult-")

## RxODE instance 1
m1 <- 
    RxODE(
        model = '
         C2 = centr/V2;
         C3 = peri/V3;
         d/dt(depot) =-KA*depot;
         d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
         d/dt(peri)  =                    Q*C2 - Q*C3;
         d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;',
        modName = "inst1",
        wd = test.dir
    )

test_that("RxODE instance 1 is created",{
    expect_equal(class(m1),"RxODE");
})

et1 <- eventTable(amount.units="ug", time.units = "hours")
et1$add.dosing(dose=10000, nbr.doses=5, dosing.interval = 24)
et1$add.sampling(0:24)
et1$add.sampling(24);
et1$add.sampling(seq(from = 24+8, to = 5*24, by = 8))


test_that("RxODE event table 1 was created",{
    expect_equal(class(et1), "EventTable")
    expect_equal(et1$get.nobs(),38);
    expect_equal(length(et1$get.dosing()[,1]), 5);
})

o1.first <- rxSolve(m1, params = c(KA=.291, CL=18.6, V2=40.2, Q=10.5, V3=297.0,
                                   Kin=1.0, Kout=1.0, EC50=200.0),
                    events = et1, 
                    inits = c(0, 0, 0, 1))

test_that("filter works",{
    expect_equal((o1.first %>% filter(time <= 5))$time,0:5);
})


test_that("distinct works",{
    expect_equal(sum((o1.first %>% distinct())$time == 24),1);
});


test_that("slice works",{
    expect_equal((o1.first %>% slice(2:4))$time,1:3);
});

test_that("top n works",{
    expect_equal((o1.first %>% top_n(3,depot))$time,c(48,72,96));
})

sn <- sample_n(o1.first,3);
sf <- sample_frac(o1.first,0.5);

test_that("sample n works",{
    expect_equal(length(sn$time),3);
})

test_that("sample frac works",{
    expect_equal(length(sf$time),19);
})

test_that("top n works",{
    expect_equal(digest(top_n(o1.first,3,depot),"sha512"),
                 "489db21b97dbccf5847401bee73f862f99bc08da84a2e7d040818a6dbb3cd14c85ac846cbb22421867ec8d3694ccecaf49189583f67054c3f8e002b0e3eca62d")
})

test_that("select works",{
    expect_equal(names(o1.first %>% select(time,depot)),c("time","depot"))
})

test_that("summerise works",{
    expect_equal((o1.first %>% summarize(avg=mean(depot)))$avg, mean(o1.first$depot));
});

test_that("summerise each works",{
    expect_equal(sum(as.data.frame(o1.first %>% group_by(time) %>% summarize_each(funs(first)))$time == 24),1);
})

test_that("summarize count works",{
    expect_equal((o1.first %>% filter(time==24) %>% count())$n,2)
})

test_that("mutate works",{
    expect_equal((o1.first %>% mutate(time=time+1))$time,et1$get.sampling()$time+1);
});

test_that("mutate each",{
    expect_equal((o1.first %>% filter(time==24 | time == 23) %>% mutate_each(funs(max)))$time,rep(24,3))
});

test_that("transmute works",{
    expect_equal(names(o1.first %>% transmute(C1=C2+C3,time=time)),c("C1","time"));
})

test_that("rename works",{
    expect_equal(names(o1.first %>% rename(Cdepot=C2)),c("time","depot", "centr", "peri","eff","Cdepot","C3"))
})

test_that("arrange works",{
    expect_equal(digest(o1.first %>% arrange(C2),"sha512"),
                 "ce3cbb41991b6cefd97b21fb72e0e09e01efe4238bd80315bd08bf52689a439127ec923a56d46aa6701896841af1845e1e73574f124f2d1c5d82b2510d71d5f2");
})
