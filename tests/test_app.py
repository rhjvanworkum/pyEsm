"""
Main Test module for the py_esm module
"""

from py_esm.app import PyEsm

def test_py_esm():
    """
    Tests the dummy function in the py_esm module
    """
    assert PyEsm.run() == "Hello World"
