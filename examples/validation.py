#!/usr/bin/env python3

import numpy as np
import scipy as sp

from pypssalib import SSA
from pypssalib import Testcase
from pypssalib import homoreaction_pdf
from pypssalib import pSSAlib

omega = 1.0
k1 = 0.016
k2 = 10.0
A0 = 25
samples = 10000

pssa = pSSAlib(SSA.SPDM)
finalpops = pssa.sample_testcase_populations(Testcase.homoreaction, [k1, k2, omega, A0], 1000, samples)

bins, counts = np.unique(finalpops.squeeze(), axis=0, return_counts=True)

empirical_pdf = counts / len(finalpops)
analytical_pdf = homoreaction_pdf(bins, k1, k2, omega)

result = np.stack((bins, empirical_pdf, analytical_pdf), axis=1, dtype=float)
KL_div = sp.special.rel_entr(analytical_pdf, empirical_pdf)

print(f"Kullback-Leibler divergence over {samples} samples is {sum(KL_div)}")
