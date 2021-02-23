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

**Where can I obtain the training data for CHASMplus?**

You can obtain the set of mutations used for training from `here <http://karchinlab.org/data/CHASMplus/formatted_training_list.txt.gz>`_.

**I want to compare my method to CHASMplus. How should I do it?**

I recommend using the precomputed scores available through OpenCRAVAT [see :ref:`quickstart-ref`]. Scores in the precompute were generated using gene-hold out cross-validation, so there is no issue when evaluating performance about training set overlap leading to overfitting. However, the scores do reflect training based on data from The Cancer Genome Atlas (TCGA). If a new method is trained using more data than is available from the TCGA, then it is recommended to create a new CHASMplus model based on the larger data set by using the CHASMplus source code.
