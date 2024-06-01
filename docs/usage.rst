=====
Usage
=====

A smsmaple usage illustrating how to simulate a trajectory for a single birth-death process::

   import pypssalib

   pssa = pypssalib.pSSAlib()
   pssa.run_testcase(m.Testcase.sbd, [1, 1, 10, 0], 2000, 2100.0, 1)
