.. deepDegron documentation master file, created by
   sphinx-quickstart on Sun Dec  1 16:01:06 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

deepDegron: de novo prediction of degrons in protein sequence
=============================================================

:Author: Collin Tokheim, Shirley Liu
:Contact: ctokheim # ds DOT dfci DOT harvard DOT edu
:Lab: `Liu Lab <https://liulab-dfci.github.io/>`_ 
:Source code: `GitHub <https://github.com/ctokheim/deepDegron>`_
:Q&A: `Biostars (tag: deepDegron) <https://www.biostars.org/t/deepDegron/>`_ 

The Ubiquitin-Proteasome System (UPS) is the primary means for selective protein degradation in cells. While the UPS may contribute upwards of 19% of mutated driver genes in cancer, a systems-level understanding of the UPS is lacking. The regulatory specificty of the UPS is thought to be governed by E3 ligases recognizing short amino acid sequence motifs, known as degrons, on substrate proteins. However, only a handful of E3 ligases has known degron motifs, hampering our capability to understand UPS regulation in normal physiology and disease.

deepDegron is a machine learning method to systematically predict the potential for a protein sequence to contain a degron. Leveraging this capability, deepDegron also allows the user to predict whether a mutation likely disrupts a degron, which may lead to increased protein stability. Furthermore, it also includes a statistical test to examine for enrichment of mutations leading to degron loss in a gene. Currently, deepDegron supports predictions for degrons at the c-terminus and n-terminus of proteins. Future updates may expand to the full proteome.

Contents:

.. toctree::
   :maxdepth: 3

   quickstart_opencravat
   download
   installation
   faq

Citation
--------

Please cite our paper:

Tokheim, C., Wang, X., Timms, R.T., Zhang, B., Mena, E.L., Wang, B., Chen, C., Ge, J., Chu, J., Zhang, W., et al. (2021). Systematic characterization of mutations altering protein degradation in human cancers. Mol Cell. `link <https://www.cell.com/molecular-cell/fulltext/S1097-2765(21)00040-X>`_
