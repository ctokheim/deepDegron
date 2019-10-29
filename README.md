# deepDegron

While a few UPS substrate mutations can be implicated in cancer based on known degrons, systematic investigation requires better degron annotation. To address this challenge, we developed a protein sequence-based model, deepDegron, that leverages data from the recently published high throughput global protein stability (GPS, (Koren et al., 2018; Timms et al., 2019)) assay of the n-terminal and c-terminal proteome to predict degrons. GPS measures the protein stability impact of peptides when attached to proteins (Yen et al., 2008), as measured by FACS-sorting of cells based on a fluorescent reporter protein (GFP, green) compared to a control reporter with no peptide attached (dsRed, red). Because the peptides consisted of known sequences and could contain degrons, deepDegron can learn sequence-rules of degron impact on protein stability.  

![](https://github.com/ctokheim/deepDegron/workflows/Python%20application/badge.svg)

## Installation

deepDegron is a python package, so please install the package dependencies via pip.

```bash
$ python -m pip install --upgrade pip
$ pip install -r requirements.txt  # install required packages
$ pyensembl install --release 75 --species human  # download human reference data
```
