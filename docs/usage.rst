=====
Usage
=====

A smsmaple usage illustrating how to simulate a trajectory for a single birth-death process::

   import pypssalib as m

   pssa = m.pSSAlib(m.SSA.SPDM)
   pssa.sample_testcase_trajectory(m.Testcase.sbd, [1, 1, 10, 0], 2100, time_start=2000.0, time_step=1)

Another example that computes the Kullback-Liebler divergence to quantify similarity between empirical and analytical PDFs for the Homoreaction test case [1]_:

.. literalinclude:: ../examples/validation.py
   :lines: 2-

.. [1] Erban R, Chapman SJ. Stochastic modelling of reaction—diffusion processes: algorithms for bimolecular reactions. Physical Biology. 2009;6(4):046001 doi: `10.1088/1478-3975/6/4/046001 <https://doi.org/10.1088/1478-3975/6/4/046001>`_