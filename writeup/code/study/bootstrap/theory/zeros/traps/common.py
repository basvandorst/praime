"""Shared helpers (read-only use of the repo's certs folder: core.py, kbessel.py)."""
import sys, os
import os as _os
CERTS = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "convex", "certs")
if CERTS not in sys.path:
    sys.path.insert(0, CERTS)
HERE = os.path.dirname(os.path.abspath(__file__))
