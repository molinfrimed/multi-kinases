# KUALA: Kinase drUgs mAchine Learning frAmework

Multi-target priority scores are contained in the results folder.
Each row contains:
* Kinase Uniprot ID
* ChEMBL ID of the ligand predicted as active
* MTPS score for kinase-ligand pair
* Status of the repurposing choice

## Machine learning models 
Trained models to predict activity for our set of kinases can be downloaded from Zenodo https://doi.org/10.5281/zenodo.6554043

## Results
A comprehensive overview of repurposable drugs for each kinase is reported on figure below and is freely accessible on [shinyapps.io](https://molinfrimed.shinyapps.io/kuala-demo/).

![Repurposable drugs distribution](https://github.com/molinfrimed/multi-kinases/blob/main/results/kinase_repurposing_distr.png?raw=true)

## Tutorial
KUALA models have been trained by using molecular descriptors computed with [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) software. Each molecule is represented as canonical SMILES.

The following steps are useful to correctly predict activity for our set of kinases:
1. Download KUALA models from Zenodo: https://doi.org/10.5281/zenodo.6867485 (PaDEL-Descriptors is also available as an extension of KNIME)
2. Compute all available molecular descritors (a comprehensive list of mandatory descriptors is reported in mandatory-list.txt)
3. Set path variables properly within kuala-demo.R script and execute it
4. Collect your results, default contained in kuala_predicted_ligands_activity.txt file in your current working directory
