import pytest
import pathlib
import pickle
from pathlib import Path

def test_demo_results():
    fn = Path(__file__).resolve().parent / 'demo_result.pickle'
    good_result = pickle.load(open(str(fn), 'rb'))

    from pydealii.pymor.demo import run
    result = run()

    assert result == good_result
