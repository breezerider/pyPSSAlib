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

Python bindings for `pSSAlib <https://github.com/breezerider/pSSAlib>`_.

Installation
============

Get latest released version from `PyPI <https://pypi.org/>`_::

    pip install pypssalib

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
