# Analyse DLS with CONTIN (for Python)

[![Testing](https://github.com/BAMresearch/analyse_dls_with_contin/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/BAMresearch/analyse_dls_with_contin/actions/workflows/tests.yml)
[![Documentation](https://github.com/BAMresearch/analyse_dls_with_contin/actions/workflows/doc.yml/badge.svg?branch=main)](https://github.com/BAMresearch/analyse_dls_with_contin/actions/workflows/doc.yml)

This repository contains Python code and a Jupyter Notebook
running the original [CONTIN program by S. Provencher](http://dx.doi.org/10.1016/0010-4655(82\)90174-6)
on every DLS measurement (dynamic light scattering, aka. photon correlation spectroscopy, PCS)
read from `*.ASC` files at the specified angle found in the given subfolders.
The expected input file format is `ALV-7004 CGS-8F Data` which is found at the first line of each file.

## Requirements

For packages that are used here and which might be missing the a clean default python installation or *venv*,
please see the [requirements](requirements.txt) file. Install the listed packages by running:

    pip install -r requirements.txt

## Problems & Help

If there are any errors with the code, please look for existing similar reports and possible solutions
or open [a **new issue** on the project page](https://github.com/BAMresearch/analyse_dls_with_contin/issues).

### Online Documentation

https://bamresearch.github.io/analyse_dls_with_contin

#### Documentation as PDF

https://bamresearch.github.io/analyse_dls_with_contin/analyse-dls-with-contin.pdf
