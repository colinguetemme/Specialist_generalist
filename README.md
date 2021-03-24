# Specialist_generalist
## Rate of specialist and generalist in a time fluctuating environment with variable rate of fluctuation

*Update of 24 Mars 2021*

### The models

To this point, there are working R code (`R/mod_burger.R` and `R/sim_burger.R`) and c++ code (`C++/mod_burger_codes/ModBurger.cpp`) for the model of Burger (see the article of [Burger(2002)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6AB0F0FA0C6991AAC624A26549D110EC/S0016672302005682a.pdf/div-class-title-fluctuating-environments-and-the-role-of-mutation-in-maintaining-quantitative-genetic-variation-div.pdf)). From this code a markdown explain a bit more in details how it works (`R/mod_burger.rmd`). Also an haploid mode has been added where we basically removed the recombination in `R/haploid_mod_burger.R`.

Can't get the stats to work as I want, in particular, the ratio gives values above 1, so I probably misunderstood the genetic variance value (since no formula in the paper).

### TO DO

**Now working on the code to get results representation and comparison with the article**.
- the results are outputed in a file that can be read on the R file `plotcpp.R`
- Want to check the stats before the simulation (quite long even in C++)

**Also should do following**:
- The genetics variance formula
- Add **halpoid** dor the burger model on C++ (add but to test)
- Get a common mesure of the specialisation of the species for both model
- Problem for Justin's model for high number of simulation, since the mutation accumulates, soon become very high number of mutation
- Get a proper output R code
