"""Pytest configuration: ensure src/ is on sys.path for development imports."""
import sys
import os

_src = os.path.join(os.path.dirname(__file__), "..", "src")
if _src not in sys.path:
    sys.path.insert(0, os.path.abspath(_src))
