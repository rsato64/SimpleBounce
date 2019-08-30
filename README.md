SimpleBounce
=============

## Overview
The SimpleBounce package calcualtes the Euclidean action for the bounce solution which contribute to the false vacuum decay.
The algorithm is based on the flow equation which is proposed in [arXiv:1907.02417](https://arxiv.org/abs/1907.02417).
For more details, please read [the documentation](https://www.dropbox.com/s/9ydtrjlsq480823/manual.pdf?dl=0).

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
* R. Sato, "A Simple Gradient Flow Equation for the Bounce Solution," [arXiv:1907.02417 \[hep-ph\]](https://arxiv.org/abs/1907.02417)
* R. Sato, "SimpleBounce : a simple package for the false vacuum decay." [arXiv:1908.10868 \[hep-ph\]](https://arxiv.org/abs/1908.10868)
