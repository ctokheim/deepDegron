.. _format-ref:

File formats
============

Mutations
+++++++++

Mutations are provided in a Mutation Annotation Format (MAF) file (specification `here <https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/>`_). 
Columns can be in any order, and only a few columns in the MAF file
are needed. The following is a list of the required columns.

* Hugo_Symbol 
* Chromosome
* Start_Position
* End_Position
* Reference_Allele
* Tumor_Seq_Allele2 
* Tumor_Sample_Barcode
* Variant_Classification

The remaining columns in the MAF specification can be 
left empty or not included. 

Only coding variants found in the Variant_Classification column will be used, which includes the following: 'Missense_Mutation', 'Silent', 'Nonsense_Mutation', 'Splice_Site', 'Nonstop_Mutation', 'Translation_Start_Site', 'Frame_Shift_Ins', 'Frame_Shift_Del', 'In_Frame_Ins', or 'In_Frame_Del'. 
