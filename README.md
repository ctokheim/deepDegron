# deepDegron

While a few UPS substrate mutations can be implicated in cancer based on known degrons, systematic investigation requires better degron annotation. To address this challenge, we developed a protein sequence-based model, deepDegron, that leverages data from the recently published high throughput global protein stability (GPS, (Koren et al., 2018; Timms et al., 2019)) assay of the n-terminal and c-terminal proteome to predict degrons. GPS measures the protein stability impact of peptides when attached to proteins (Yen et al., 2008), as measured by FACS-sorting of cells based on a fluorescent reporter protein (GFP, green) compared to a control reporter with no peptide attached (dsRed, red). Because the peptides consisted of known sequences and could contain degrons, deepDegron can learn sequence-rules of degron impact on protein stability.  

[![Build Status](https://travis-ci.org/ctokheim/deepDegron.svg?branch=master)](https://travis-ci.org/ctokheim/deepDegron)

## Documentation

For more documentation, please visit our [website documentation](https://deepdegron.readthedocs.io/en/latest/index.html).

## Installation 

We recommend that you use python 3.7 to run deepDegron. 

### pip

The easiest way to install deepDegron is through [pip](https://pip.pypa.io/en/stable/).

```bash
$ pip install deepDegron
```

### From source

As a first step, please change to the top-level directory in the deepDegron source code.

You can install the dependencies of deepDegron using [conda](https://docs.conda.io/en/latest/). Once you have conda installed, create an environment to run deep degron using the following commands:

```bash
$ conda env create -f environment.yml  # install dependencies
$ source activate deepDegron  # activate environment
$ pyensembl install --release 75 --species human  # download human reference data
$ python setup.py install  # install deepDegron
```

An alternative way to install the python dependencies is to use pip.

```bash
$ python -m pip install --upgrade pip
$ pip install -r requirements.txt  # install required packages
$ pyensembl install --release 75 --species human  # download human reference data
$ python setup.py install  # install deepDegron
```
