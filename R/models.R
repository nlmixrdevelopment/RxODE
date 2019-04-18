##' RxODE one compartment model (solved)
##' @examples
##' pk1cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot # Oral
##' pk1cmt %>% solve(et(amt=100,cmt=central,timeUnits="hr")) %>% plot # IV
"pk1cmt"

##' RxODE two compartment model (solved)
##' @examples
##' pk2cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot # Oral
##' pk2cmt %>% solve(et(amt=100,cmt=central,timeUnits="hr")) %>% plot # IV
"pk2cmt"

##' RxODE three compartment model (solved)
##' @examples
##' pk3cmt %>% solve(et(timeUnits="hr") %>% et(amt=100)) %>% plot # Oral
##' pk3cmt %>% solve(et(amt=100,cmt=central,timeUnits="hr")) %>% plot # IV
"pk3cmt"

