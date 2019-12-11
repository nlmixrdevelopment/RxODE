library(tidyverse)

itdepends::dep_usage_pkg("RxODE") %>%
    count(pkg, sort = TRUE) %>%
    as.data.frame
