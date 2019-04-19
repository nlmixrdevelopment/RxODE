## ---- echo=FALSE---------------------------------------------------------
options(cli.unicode=FALSE, crayon.enabled=FALSE);
options(knitr.table.format = "html")
htmltools::img(src = knitr::image_uri("logo.png"), 
               alt = 'RxODE', 
               style = 'position:absolute; top:0; right:0; padding:10px; border: 0;')

## ----results="asis"------------------------------------------------------
library(RxODE)

Vtpol2 <- RxODE({
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
})

et <- eventTable();
et$add.sampling(seq(0, 10, length.out=200));
et$add.dosing(20, start.time=0);

s1 <- Vtpol2 %>%  solve(et, method="lsoda")
rxHtml(s1)

## ----results="asis"------------------------------------------------------
s2 <- Vtpol2 %>%  solve(c(mu=1000), et)
rxHtml(s2)

