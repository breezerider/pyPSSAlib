========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|

    * - build
      - |github-actions|

    * - package
      - | |license| |version| |wheel| |supported-versions|
        | |commits-since|

.. |docs| image:: https://readthedocs.org/projects/pypssalib/badge/?style=flat
    :target: https://pypssalib.readthedocs.io/
    :alt: Documentation Status

.. |github-actions| image:: https://github.com/breezerider/pypssalib/actions/workflows/github-actions.yml/badge.svg
    :alt: GitHub Actions Build Status
    :target: https://github.com/breezerider/pypssalib/actions

.. |license| image:: https://img.shields.io/badge/license-BSD-green?style=flat
    :alt: PyPI Package license
    :target: https://test.pypi.org/project/pypssalib

.. |version| image:: https://img.shields.io/badge/test.pypi-v0.0.0-informational?style=flat
    :alt: PyPI Package latest release
    :target: https://test.pypi.org/project/pypssalib

.. |wheel| image:: https://img.shields.io/badge/wheel-yes-success?style=flat
    :alt: PyPI Wheel
    :target: https://test.pypi.org/project/pypssalib

.. |supported-versions| image:: https://img.shields.io/badge/python-3.8_|_3.9_|_3.10|_3.11-informational?style=flat
    :alt: Supported Python versions
    :target: https://test.pypi.org/project/pypssalib

.. |commits-since| image:: https://img.shields.io/github/commits-since/breezerider/pypssalib/v0.0.0.svg
    :alt: Commits since latest release
    :target: https://github.com/breezerider/pypssalib/compare/v0.0.0...main

.. end-badges

Python bindings for `pSSAlib <https://github.com/breezerider/pSSAlib>`_ supporting following stoachastic simulation algorithms (SSAs):

* Gillespie’s direct method (DM) [1]_ as a reference;
* partial-propensity direct method (PDM) [2]_;
* sorting partial-propensity direct method (SPDM) [2]_;
* partial-propensity SSA with Composition-Rejection Sampling (PSSA-CR) [3]_;

Package provides an interface to sample individual trajectories as well as steady-stat populations. A number of build-in models are provided:

* Cyclic Linear Chain (CLC) model [2]_
* Colloidal Aggregation (CA) model [2]_
* Homoreaction model [4]_
* Single-Species Birth-Death (SBD) model
* The bacterial Two-Component System [5]_
* The Enzymatic Degradation process [6]_


Installation
============

Get latest released version from `Test PyPI <https://test.pypi.org/>`_::

    pip install -i https://test.pypi.org/simple/ pypssalib

You can also install the in-development version with::

    pip install https://github.com/breezerider/pypssalib/archive/main.zip


Documentation
=============


https://pypssalib.readthedocs.io/


License
=======

- Source code: `BSD-3-Clause <https://choosealicense.com/licenses/bsd-3-clause/>`_ license unless noted otherwise in individual files/directories
- Documentation: `Creative Commons Attribution-ShareAlike 4.0 <https://creativecommons.org/licenses/by-sa/4.0/>`_ license


Development
===========

To run all the tests issue this command in a terminal::

    tox

.. [1] Gillespie DT. Exact stochastic simulation of coupled chemical reactions. The Journal of Physical Chemistry. 1977;81(25):2340–2361. doi: `10.1021/j100540a008 <https://doi.org/10.1021/j100540a008>`_
.. [2] Ramaswamy R, Gonzalez-Segredo N, Sbalzarini IF. A new class of highly efficient exact stochastic simulation algorithms for chemical reaction networks. J Chem Phys. 2009;130(24):244104 doi: `10.1063/1.3154624 <https://doi.org/10.1063/1.3154624>`_
.. [3] Ramaswamy R, Sbalzarini IF. A partial-propensity variant of the composition-rejection stochastic simulation algorithm for chemical reaction networks. The Journal of Chemical Physics. 2010;132(4):044102 doi: `10.1063/1.3297948 <https://doi.org/10.1063/1.3297948>`_
.. [4] Erban R, Chapman SJ. Stochastic modelling of reaction—diffusion processes: algorithms for bimolecular reactions. Physical Biology. 2009;6(4):046001 doi: `10.1088/1478-3975/6/4/046001 <https://doi.org/10.1088/1478-3975/6/4/046001>`_
.. [5] Kim, J.-R. & Cho, K.-H. The multi-step phosphorelay mechanism of unorthodox two-component systems in e. coli realizes ultrasensitivity to stimuli while maintaining robustness to noises. Comput. Biol. Chem. 2006, doi: `10.1016/j.compbiolchem.2006.09.004 <https://doi.org/10.1016/j.compbiolchem.2006.09.004>`_
.. [6] Fröhlich, F. et al. Inference for stochastic chemical kinetics using moment equations and system size expansion. PLOS Comput. Biol. 2016, doi: `10.1371/journal.pcbi.1005030 <https://doi.org/10.1371/journal.pcbi.1005030>`_
