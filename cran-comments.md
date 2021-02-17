# CRAN Comments

We see:

```
> ## This example uses `rxbinom` directly in the model
>
> rx <- RxODE({
+   a = rxbinom(1,0.5)
+ })
=================================================================
==3422613==ERROR: AddressSanitizer: odr-violation (0x7fceb25427c0):
   [1] size=8 '_compareFactorVal'
/data/gannet/ripley/R/packages/incoming/RxODE.Rcheck/RxODE/include/RxODE_model_shared.c:70:27
   [2] size=8 '_compareFactorVal'
/data/gannet/ripley/R/packages/incoming/RxODE.Rcheck/RxODE/include/RxODE_model_shared.c:70:27
These globals were registered at these points:
   [1]:
     #0 0x7fcecd4abcc8  (/lib64/libasan.so.6+0x37cc8)
     #1 0x7fcecdebf801 in call_init.part.0
(/lib64/ld-linux-x86-64.so.2+0x11801)

   [2]:
     #0 0x7fcecd4abcc8  (/lib64/libasan.so.6+0x37cc8)
     #1 0x7fcecdebf801 in call_init.part.0
(/lib64/ld-linux-x86-64.so.2+0x11801)

```

We believe we fixed this by defining new variables for every ODE model
compiled. Our tests show that this new approach does not have any ODR
violations.  We also checked the binary symbols in 2 different
compiled objects and made sure they were different.

Please check on your machines if this is correct.

On our end, we cannot reproduce your results; On our machine the
ODR-check for the last version was clean.

Our guess is the extra `-march=native` means that our chip doesn't
match the testing machine. By your back-trace (which we appreciate) we
believe that the error also points to the same memory location, same
type of object and same place in the code. If I understand ODR
correctly, this should be valid by the ODR rules.

Our guess it the prior error is due to the over-optimization of
`-march=native`, causing either buggy code or false positive odr
violations.

