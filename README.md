# Specialist_generalist
## Rate of specialist and generalist in a time fluctuating environment with variable rate of fluctuation

*update of 3 february 2021*

To this point, there is a working R code (R/mod_burger.R and R/sim_burger.R) for the model of Burger (see the article of [Burger(2002)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6AB0F0FA0C6991AAC624A26549D110EC/S0016672302005682a.pdf/div-class-title-fluctuating-environments-and-the-role-of-mutation-in-maintaining-quantitative-genetic-variation-div.pdf)). Fron this code a markdown explain a bit more in details how it works (R/mod_burger.rmd). Also an haploid mode has been added where we basically removed the recombination in R/haploid_mod_burger.R.

The translation in C++ is almost finish with the initialisation of parameters, the environment and the recombinations working (C++/ModBurger.cpp).
