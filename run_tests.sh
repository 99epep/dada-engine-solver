#!/bin/bash
set -e
cd /tmp/test_dada_solver
PYTHONPATH=src python3 -m pytest tests/test_doty.py -v
