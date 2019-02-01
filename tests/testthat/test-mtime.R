
library(dplyr)

mod <- RxODE({
    ka = 0.04
    kel = 0.6
    mt1 = 5
    mt2 = 7
    mtime(t1) = mt1
    mtime(t2) = mt2
    if (t >=t1){
        ka = 0.4
    }
    if (t >=t2){
        ka = 4
    }
    d/dt(intestine) = -ka*intestine
    d/dt(blood)     = ka*intestine - kel*blood
})


et <- eventTable() %>% add.dosing(dose=3, nbr.doses=1) %>%
    add.sampling(seq(0, 48, length.out=200))

et <- rbind(et$get.EventTable(),data.frame(time=0,evid=c(10L, 11L), amt=NA_real_)) %>%
    arrange(time,-evid)


s <- solve(mod,et)

library(ggplot2);
ggplot(s,aes(time,blood)) + geom_line()
