FAQ
===

**Who should I contact if I encounter a problem?**

If you believe your problem may be encountered by other users,
please post the question on `biostars <https://www.biostars.org/>`_.
Check to make sure your question has not been already answered 
by looking at posts with the tag `deepDegron <https://www.biostars.org/t/deepDegron>`_.
Otherwise, create a new post with the deepDegron tag. We will be checking
biostars for questions. You may also contact me directly at
ctokheim AT ds DOT dfci DOT harvard DOT edu.

**Can I run deepDegron using mutations annotated on hg19 and/or hg38?**

Yes, deepDegron supports both hg19 and hg38. You will, however, need to download the relevant reference data for pyensembl first. We recommend ensembl release 75 for hg19:

.. code-block:: bash
    
    $ pyensembl install --release 75 --species human  # download hg19 human reference data 

For hg38, we recommend ensembl release 95:

.. code-block:: bash
    
    $ pyensembl install --release 95 --species human  # download hg38 human reference data 

To correctly specify which reference genome you are using, please supply the relevant ensembl release number to --ensembl-release flag in deepDegron.

.. code-block:: bash
    
    $ deepDegron_test [options] --ensembl-release 75 # for hg19
    $ deepDegron_test [options] --ensembl-release 95 # for hg38

**Where can I obtain the training data for deepDegron?**

You can obtain the set of mutations used for training from github for `c-terminal degrons <https://raw.githubusercontent.com/ctokheim/deepDegron/master/train_data/gps_cterminal_degron_screen.txt>`_ and `n-terminal degrons <https://raw.githubusercontent.com/ctokheim/deepDegron/master/train_data/gps_nterminal_degron_screen.txt>`_.
