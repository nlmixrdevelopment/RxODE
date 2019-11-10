# RxODE xxxx
* Added ability to prune branching with `rxPrune`. This converts
  `if`/`else` or `ifelse` to single line statements without any
  `if`/`then` branching within them.
* Added ability to take more complex conditional expressions, including:
  - `ifelse(expr, yes, no)`
  - `x = (x==1)*1 + (!(x==1))*2`
  - `if (logic){ expr} else if (logic) {expr} else {}`.  The preferred
    syntax is still only `if`/`else` and the corresponding parsed code
    reflects this preference.
* Switched to `symengine` instead of using `sympy`
  - Remove dependence on python.
  - Since symengine is C-based and doesn't require the python
    interface it is much faster than `sympy`, though some functions in
    `sympy` are no longer accessible.
	
* Added new ODE solving method "indLin", or inductive linearization.
  When the full model is a linear ODE system this becomes simply the
  matrix exponential solution.  Currently this requires a different
  setup.

* Added arbitrary function definition to RxODE using `rxFun`
  - Requires function, arguments and corresponding C-code
  - Derivatives (if required) can be added to the derivative table `rxD`

# RxODE v0.9.1-8
* Added the ability to have an input parameter to be assigned to a new
  value (Issue #135)
* Added LINPACK authors as contributors
* Added a `NEWS.md` file to track changes to the package
