
library(dplyr)

et <- eventTable(time.units="days")
    et$add.sampling(seq(0,10,by=1/24))
    et$add.dosing(dose=2/24,rate=2,start.time=0,
                  nbr.doses=10,dosing.interval=1)

et1 <- et$get.EventTable() %>% filter(is.na(amt) | amt > 0) %>%
    mutate(evid=ifelse(evid==0, 0, 90101)) %>%
    mutate(amt=ifelse(is.na(amt),NA_real_,0.4))

et2 <- et1 %>% filter(evid !=0) %>% mutate(evid=70101)

et3 <- rbind(et1, et2) %>% arrange(time,-evid)

mod.rate <- RxODE({
        a = 6
        b = 0.6
        f = 1
        ri = 10
        li = 0
        d/dt(intestine) = -a*intestine
        f(intestine) = f
        rate(intestine) = ri
        lag(intestine) = li
        d/dt(blood)     = a*intestine - b*blood
    })

library(ggplot2)

f <- rxSolve(mod.rate, et3, c(ri=2,f=2,li=0.3)); ggplot(f,aes(time,intestine)) + geom_line()


et <- eventTable(time.units="days")
et$add.sampling(seq(0,10,by=1/24))
et$add.dosing(dose=2/24,rate=2,start.time=0,
              nbr.doses=1,dosing.interval=1)

et1 <- et$get.EventTable() %>% filter(is.na(amt) | amt > 0) %>%
    mutate(evid=ifelse(evid==0, 0, 90101)) %>%
    mutate(amt=ifelse(is.na(amt),NA_real_,0.4))

et2 <- et1 %>% filter(evid !=0) %>% mutate(evid=70101)

et3 <- rbind(et1, et2, data.frame(time=1,evid=101,amt=0.2)) %>% arrange(time,-evid)

f <- rxSolve(mod.rate, et3, c(ri=2,f=2,li=0.3)); ggplot(f,aes(time,intestine)) + geom_line()

f <- rxSolve(mod.rate, et3, c(ri=0.5,f=2,li=0.3)); ggplot(f,aes(time,intestine)) + geom_line()
