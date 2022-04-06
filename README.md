SimpleBounce
=============

## Overview
The SimpleBounce package calcualtes the Euclidean action of the bounce solution which contributes to the false vacuum decay.
The algorithm is based on the flow equation which is proposed in [arXiv:1907.02417](https://arxiv.org/abs/1907.02417).
For more details of the code, please read [arXiv:1908.10868](https://arxiv.org/abs/1908.10868).

## How to use
The following sample codes are complied by ``make``.

#### ``sample1.x``
This calculate the bounce action for a single scalar field model with V(phi) = (phi^2)/2 - (phi^3)/3.
The source file is ``sample1.cc``

#### ``benchmark/compare_with_cosmotransitions/``
``run_simplebounce.sh`` calculate the bounce solution and the bounce action for several sample models,
which are taken from
Table 1 in [arXiv:1901.03714](https://arxiv.org/abs/1901.03714)
and Eq. (40-43) in [arXiv:1906.10829](https://arxiv.org/abs/1906.10829).
You can compare the results and performance with [CosmoTransitions](http://clwainwright.github.io/CosmoTransitions/)
 by using ``benchmark/compare_with_cosmotransitions/run_cosmotransitions.py``.

#### ``benchmark/change_n_tau/``
``change_n_tau.x`` calculate the bounce action and measure the runtime for different parameter settings.
``run.sh`` plots the results.

## References
* R. Sato, "A Simple Gradient Flow Equation for the Bounce Solution," [arXiv:1907.02417 \[hep-ph\]](https://arxiv.org/abs/1907.02417), [Phys.Rev.D 101 (2020) 1, 016012](https://doi.org/10.1103/PhysRevD.101.016012)
* R. Sato, "SimpleBounce : a simple package for the false vacuum decay." [arXiv:1908.10868 \[hep-ph\]](https://arxiv.org/abs/1908.10868), [Comput.Phys.Commun. 258 (2021) 107566](https://doi.org/10.1016/j.cpc.2020.107566)
