# sundialy
[![PyPI version](https://badge.fury.io/py/sundialy.svg)](https://badge.fury.io/py/sundialy)
[![Tests](https://github.com/AttackingOrDefending/sundialy/actions/workflows/tests.yml/badge.svg)](https://github.com/AttackingOrDefending/sundialy/actions/workflows/tests.yml)
[![Build](https://github.com/AttackingOrDefending/sundialy/actions/workflows/build.yml/badge.svg)](https://github.com/AttackingOrDefending/sundialy/actions/workflows/build.yml)
[![Mypy](https://github.com/AttackingOrDefending/sundialy/actions/workflows/mypy.yml/badge.svg)](https://github.com/AttackingOrDefending/sundialy/actions/workflows/mypy.yml)
[![codecov](https://codecov.io/gh/AttackingOrDefending/sundialy/branch/main/graph/badge.svg)](https://codecov.io/gh/AttackingOrDefending/sundialy)
[![Docs](https://readthedocs.org/projects/sundialy/badge/?version=latest)](https://sundialy.readthedocs.io/en/latest/?badge=latest)

A sundial library for python.

## Installing

Download and install the latest release:

    pip install sundialy

## Features

* Includes mypy typings

* Create a sundial
```python
from sundialy import AnalemmaticHorizontal

# An analemmatic sundial for Los Angeles.
sundial = AnalemmaticHorizontal(latitude=34, longitude=-118)
```

* View how it will look

```python
sundial.create_sundial("sundial.jpg")
```
![Sundial](images/sundial.jpg)

* Use the tools

```python
from sundialy.tools import SPA

spa_results = SPA(2020, 12, 31, 23, 59, 59, 0, 0, 0, 0, pressure=1000, temperature=10, omega=0, gamma=0)
```


## License
sundialy is licensed under the MIT License. Check out LICENSE for the full text.
