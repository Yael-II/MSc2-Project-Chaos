# Order and Chaos in a 2D potential

Moussouni, Yaël (MSc student) & Bhat, Junaid Ramzan (MSc student)

Université de Strasbourg, CNRS, Observatoire astronomique de Strasbourg, 
UMR 7550, F-67000, Strasbourg, France

## Requirements

The project requires `python` (tested with version 3.13.1), with the `venv` module and `pip`, a `bash` interpreter (`/usr/bin/env bash` by default), and at least 450 kiB of available space.

## Installation

Place all the content of the archive in a directory. You can use:
```bash
git clone https://github.com/Yael-II/MSc2-Project-Chaos
```
Then, in this directory, to create the virtual environment, run:
```bash
python3 -m venv ./venv
```
Then, to install the requirements, run:
```bsah
source activate.sh && pip install -r requirements.txt && deactivate
```
Finally, you will have to authorize execution of the shell scripts with
```bash
chmod u+x *.sh
```

## Usage

##

## Main Poincaré Sections (Linear and Parallel)

...

The result of the simulation is saved in the Output directory, with the prefix `poincare_section_linear_` followed by 1/E (e.g. 12 for E = 1/12). The ASCII file contains all the Poincare section points (y on the first line, v on the second line) separated by blank spaces.
