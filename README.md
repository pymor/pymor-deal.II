# pyMOR bindings for deal.II

This repository contains basic pyMOR binding for the [deal.II](https://dealii.org) Finite Element library.

## Installation
A working deal.II installation is required. On Debian/Ubuntu, deal.II can be installed as follows:
```
sudo apt install libdeal.ii-dev
```

Install and run demo script:
```
python -m pip install -e .
python ./src/pymor_dealii/pymor/demo.py
```
