FeldtLib.jl
===========

This is where I collect useful Julia code for which I have not yet created a package, or where existing Julia packages (by others) do not yet (to my knowledge) include it.

# Installation

Just install from the github by calling:

    Pkg.clone("https://github.com/robertfeldt/FeldtLib.jl")

from a Julia repl.

# Contents

## Comparing two empirical distributions

* Baumgartner-Weis-Schindler test: 
  - `baumgartner_weis_schindler_statistic(x, y)`
  - `bws_test_sampled(x, y, numsamples = 10000)`