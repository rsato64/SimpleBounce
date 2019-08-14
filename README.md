SimpleBounce
=============

## Overview
The SimpleBounce package calcualtes the Euclidean action for the bounce solution which contribute to the false vacuum decay.
The algorithm is based on the flow equation which is proposed in [arXiv:1907.02417](https://arxiv.org/abs/1907.02417).
For more details, please read [the documentation](https://www.dropbox.com/s/9ydtrjlsq480823/manual.pdf?dl=0).

## How to use
The sample codes are complied by ``make``. The executable files are ``sample*.x``.
They calculate the bounce solution and the bounce action for several sample models,
which are taken from
Table 1 in [arXiv:1901.03714](https://arxiv.org/abs/1901.03714)
and Eq. (40-43) in [arXiv:1906.10829](https://arxiv.org/abs/1906.10829).
You can compare the results and performance with [CosmoTransitions](http://clwainwright.github.io/CosmoTransitions/)
 by using ``comparison/run_sample.sh`` and ``comparison/run_ct.py``.
