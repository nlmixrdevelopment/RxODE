# RxODE 0.9.0-9
* Added ability to prune branching with `rxPrune`. This converts
  `if`/`else` or `ifelse` to single line statements without any
  `if`/`then` branching within them.
* Added ability to take more complex conditional expressions, including:
  - `ifelse(expr, yes, no)`
  - `x = (x==1)*1 + (!(x==1))*2`
  - `if (logic){ expr} else if (logic) {expr} else {}`.  The preferred
    syntax is still only `if`/`else` and the corresponding parsed code
    reflects this preference.
* Added a `NEWS.md` file to track changes to the package.
