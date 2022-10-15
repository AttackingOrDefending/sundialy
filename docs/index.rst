sundialy: A sundial library for python
======================================


.. image:: https://badge.fury.io/py/sundialy.svg
    :target: https://badge.fury.io/py/sundialy
    :alt: PyPI version
.. image:: https://github.com/AttackingOrDefending/sundialy/actions/workflows/tests.yml/badge.svg
    :target: https://github.com/AttackingOrDefending/sundialy/actions/workflows/tests.yml
    :alt: Tests
.. image:: https://github.com/AttackingOrDefending/sundialy/actions/workflows/build.yml/badge.svg
    :target: https://github.com/AttackingOrDefending/sundialy/actions/workflows/build.yml
    :alt: Build
.. image:: https://github.com/AttackingOrDefending/sundialy/actions/workflows/mypy.yml/badge.svg
    :target: https://github.com/AttackingOrDefending/sundialy/actions/workflows/mypy.yml
    :alt: Mypy
.. image:: https://codecov.io/gh/AttackingOrDefending/sundialy/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/AttackingOrDefending/sundialy
    :alt: codecov
.. image:: https://readthedocs.org/projects/sundialy/badge/?version=latest
    :target: https://sundialy.readthedocs.io/en/latest/?badge=latest
    :alt: Docs

Installing
----------

Download and install the latest release:

>>> pip install sundialy

Features
--------

* Includes mypy typings

* Import the library

>>> from sundialy import AnalemmaticHorizontal

* Create a sundial
>>> sundial = AnalemmaticHorizontal(latitude=34, longitude=-118)  # An analemmatic sundial for Los Angeles.

* View how it will look

>>> sundial.create_sundial("sundial.jpg")

.. image:: ../images/sundial.jpg
|br|

* Use the tools

>>> from sundialy.tools import SPA
>>> spa_results = SPA(2020, 12, 31, 23, 59, 59, 0, 0, 0, 0, pressure=1000, temperature=10, omega=0, gamma=0)

Contents
--------

.. toctree::
    :maxdepth: 1

    analemmatic
    tools
.. |br| raw:: html

   <br />
