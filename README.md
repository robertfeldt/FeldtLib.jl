FeldtLib.jl
===========

This is where I collect useful Julia code for which I have not yet created a package, or where existing Julia packages (by others) do not yet (to my knowledge) include it.

# Installation

Just install from the github by calling:

    Pkg.clone("https://github.com/robertfeldt/FeldtLib.jl")

from a Julia repl.

# Contents

## Comparing two set of samples (empirical distributions)

* Baumgartner-Weis-Schindler statistic and test: 
  - `bws_statistic(x, y)`
  - `bws_test_sampled(x, y, numsamples = 10000)`

* There is also an add-on layer to these tests to try to give them the same interface as in HypothesisTests.jl:
  - `BWSTest(x, y)`
  - `pvalue(t::BWSTest)`

* Utilities used in these comparisons:
  - `ranks(x, y)` to calculate the combined ranks of all values in x and y
  - `factorial`
  - `combinations`