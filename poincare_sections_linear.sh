#!/usr/bin/env bash

source activate.sh
venv/bin/python Source/main_poincare_sections_linear.py
venv/bin/python Source/plot_poincare_sections.py
