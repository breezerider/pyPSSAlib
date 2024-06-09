
Changelog
=========

Unreleased
----------

* Update usage example
* Fix deadlock in sampling methods when an atomic flag was not reset

0.1.0-dev0 (2024-06-09)
-----------------------

* Python bindings for pSSAlib in pSSAlib class
* Set simulation algorithm (one of DirectMethod, PartialPropensityDirectMethod, SortingPartialPropensityDirectMethod, CompositionRejectionSampling) via method property
* Collect individual trajectories via sample_testcase_trajectory method
* Collect populations at a certain time point via sample_testcase_trajectory method
* Provide following test cases: CyclicLinearChain, ColloidalAggregation, Homoreaction, SingleBirthDeath, TwoComponentSystem, EnzymaticDegration
* Implement interface to analytic PDF for Homoreaction model

