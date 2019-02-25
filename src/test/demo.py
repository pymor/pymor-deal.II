import pytest
import pathlib
import pickle
from pathlib import Path
import numpy as np

def test_demo_results():
    fn = Path(__file__).resolve().parent / 'demo_result.pickle'
    good_result = pickle.load(open(str(fn), 'rb'))

    from pymor_dealii.pymor.demo import run
    result = run()

    compare = ['errors', 'basis_sizes', 'rel_errors', 'estimates']
    for key in compare:
        assert np.allclose(result[key], good_result[key])
