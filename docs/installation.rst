Installation
------------

.. image:: https://travis-ci.org/ctokheim/deepDegron.svg?branch=master
    :target: https://travis-ci.org/ctokheim/deepDegron

deepDegron has only been tested on linux operating systems. We recommend that you use **python 3.7** to run deepDegron.


Installing from source
~~~~~~~~~~~~~~~~~~~~~~

Releases
++++++++

First download the deepDegron source code on on `github <https://github.com/ctokheim/deepDegron>`_.

Once downloaded, please change to the top-level directory in the deepDegron source code.

Conda Environment
+++++++++++++++++

We recommend using `conda <https://conda.io/docs/>`_ to install the deepDegron dependencies.

.. code-block:: bash

   $ conda env create -f environment.yml  # create environment for deepDegron
   $ source activate deepDegron # activate environment for deepDegron
   $ pyensembl install --release 75 --species human  # download human reference data
   $ python setup.py install  # install deepDegron

Make sure the deepDegron environment is activated when you want to run deepDegron.

Pip
+++

An alternative way to install the python dependencies is to use pip.

.. code-block:: bash

   $ python -m pip install --upgrade pip
   $ pip install -r requirements.txt  # install required packages
   $ pyensembl install --release 75 --species human  # download human reference data
   $ python setup.py install  # install deepDegron
