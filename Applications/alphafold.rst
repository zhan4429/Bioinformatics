Alphafold
==============================

Introduction
~~~~~~~~
``Alphafold`` is a protein structure prediction tool developed by DeepMind (Google). It uses a novel machine learning approach to predict 3D protein structures from primary sequences alone. The source code is available on `Github`_. It has been deployed in all RCAC clusters, supporting both CPU and GPU.   

It also relies on a huge database. The full database (~2.2TB) has been downloaded and setup for users.  

Protein struction prediction by alphafold is performed in the following steps:

* Search the amino acid sequence in uniref90 database by jackhmmer (using CPU)
* Search the amino acid sequence in  mgnify database by jackhmmer (using CPU)
* Search the amino acid sequence in pdb70 database (for monomers) or pdb_seqres database (for multimers) by hhsearch (using CPU)
* Search the amino acid sequence in bfd database and uniclust30 (updated to uniref30 since v2.3.0) database by hhblits (using CPU)
* Search structure templates in pdb_mmcif database (using CPU)
* Search the amino acid sequence in uniprot database (for multimers) by jackhmmer (using CPU)
* Predict 3D structure by machine learning (using CPU or GPU)
* Structure optimisation with OpenMM (using CPU or GPU)

