# Specialist_generalist
## Rate of specialist and generalist in a time fluctuating environment with variable rate of fluctuation

*Update of 15 Mars 2021*

To this point, there are working R code (`R/mod_burger.R` and `R/sim_burger.R`) and c++ code (`C++/mod_burger_codes/ModBurger.cpp`) for the model of Burger (see the article of [Burger(2002)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6AB0F0FA0C6991AAC624A26549D110EC/S0016672302005682a.pdf/div-class-title-fluctuating-environments-and-the-role-of-mutation-in-maintaining-quantitative-genetic-variation-div.pdf)). From this code a markdown explain a bit more in details how it works (`R/mod_burger.rmd`). Also an haploid mode has been added where we basically removed the recombination in `R/haploid_mod_burger.R`.

Can't get the stats to work as I want, in particular, the ratio gives values above 1, so I probably misunderstood the genetic variance value (since no formula in the paper).

Build a new Rmd that explain how to use on R the results produce on C++.

The code idea by Justin is in creation, but the basic idea works, need to add the outputs and should be ready to plot the first results

**Now working on the code to get results representation and comparison with the article**.
- the results are outputed in a file that can be read on the R file `plotcpp.R`
- Want to check the stats before the simulation (quite long even in C++)

**Also should do following**:
- The genetics variance formula
- Add **halpoid** dor the burger model on C++
- Continue on the code of Justin
- comment correctly the Burger code




