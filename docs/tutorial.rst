.. _tut-ref:

Tutorial
========

In this tutorial we will be investigating somatic mutations found in GATA3 in breast cancer samples from The Cancer Genome Atlas (TCGA). Note, this analysis could equally apply to other types of variants, such as germline or de novo variants, as well.

Scoring mutations for impacting degrons
---------------------------------------

deepDegron computes a degron potential score to represent the liklihood a protein sequence contains a degron. Mutations may lead to a change in degron potential. The difference between the degron potential of the mutant compared to wildtype sequence is what we call "delta degron potential". The more negative this score is, the more deepDegron predicts a degron has likely been disrupted by a mutation.

The first step to score mutations is to download the trained c-terminal degron models, as well as the necessary data file of mutations.
Here, we are using mutations found in the GATA3 gene in breast cancer. 

.. code-block:: bash
    
    $ wget https://github.com/ctokheim/deepDegron/raw/master/models/cterm/neural_network_pos_specific.pickle
    $ wget https://github.com/ctokheim/deepDegron/raw/master/models/cterm/neural_network_bag_of_amino_acids.pickle
    $ wget https://raw.githubusercontent.com/ctokheim/deepDegron/master/tests/data/gata3_mutations.txt

To obtain delta degron potential scores for mutations, you can use the deepDegron_score command.

.. code-block:: bash
    
    $ deepDegron_score -i gata3_mutations.txt -c models/cterm/neural_network_pos_specific.pickle,models/cterm/neural_network_bag_of_amino_acids.pickle -o GATA3_delta_degron_potential.txt

You will notice the output file will contain all mutations that impact the c-terminal protein sequence of GATA3, as well as the delta degron potential scores for each mutation. Notice for GATA3 that many are frameshift indels that have a very negative delta degron potential, indicating likely degron loss. You results should match the results seen `here <https://raw.githubusercontent.com/ctokheim/deepDegron/master/docs/GATA3_delta_degron_potential.txt>`_.


Statistical test
----------------

DeepDegron also can test whether there is a significant enrichment for mutations that likely lead to degron loss. This helps to avoid non-significant cases where a degron loss mutation may have happened by chance, and may not play a role in a given phenotype/condition.

To run the deepDegron statistical test on GATA3, use the deepDegron_test command.  

.. code-block:: bash
    
    $ deepDegron_test -i gata3_mutations.txt -ns 100 -c models/cterm/neural_network_pos_specific.pickle,models/cterm/neural_network_bag_of_amino_acids.pickle -o GATA3_result.txt

Because in this example we are analyzing c-terminal degrons, we supplied the trained deepDegron models using the -c flag. However, for analyzing n-terminal degrons, the -n flag should be used. Additionally, for this toy example, we used only 100 simulations (-ns parameter), but in practical applications this should be much larger (e.g. 10,000). Note, increasing the number of simulations increases precision but has a longer run time.

Your result should show a delta degron potential of -23 and a p-value of 0.0 (beyond resolution of the 100 simulations) for the GATA3 data. It should match the results available `here <https://raw.githubusercontent.com/ctokheim/deepDegron/master/docs/GATA3_result.txt>`_.
