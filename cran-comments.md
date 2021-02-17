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
