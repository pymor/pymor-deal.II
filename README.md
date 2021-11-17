Requirements
============

`deal.II` with development headers

Usage
=====

```
git clone --recurse-submodules https://github.com/pymor/pymor-deal.II
cd pymor-deal.II
python -m pip install pymor
python -m pip install -e .
python test/demo.py
```


Development notes
=================

We have a `pre-commit` config that is also run/checked in CI.
Python source is formatted with `black`.
C++ source is formatted with `clang-format`.
