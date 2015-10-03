library("RxODE")

# create dosing and observation (sampling) events
# QD 50mg dosing, 5 days followed by 25mg 5 days

qd <- eventTable(amount.units = "mg", time.units = "days")

qd$add.dosing(dose=50, nbr.doses=5, dosing.interval = 1, do.sampling=FALSE)

# sample the system's drug amounts hourly the first day, then every 12 hours
# for the next 4 days
qd$add.sampling(seq(from = 0, to = 1, by = 1/24))
qd$add.sampling(seq(from = 1, to = 5, by = 12/24))

print(qd$get.dosing())     # table of dosing records
print(qd$get.nobs())   # number of observation (not dosing) records

# BID dosing, 5 days

bid <- eventTable("mg", "days")  # only dosing
bid$add.dosing(dose=10000, nbr.doses=2*5, 
   dosing.interval = 12, do.sampling=FALSE)

# Use the copy() method to create a copy (clone) of an existing
# event table (simple assignments just create a new reference to
# the same event table object (closure)).

bid.ext <- bid$copy()      # three-day extension for a 2nd cohort
bid.ext$add.dosing(dose = 5000, nbr.doses = 2*3,
   dosing.interval = 12, do.sampling = TRUE)


