# **HOPPred**

A computational method for the prediction of peptide hormones peptide sequences using an ensemble of machine learning and similarity-based methods

# Introduction
Peptide hormones serve as genome-encoded signal transduction molecules that play essential roles in multicellular organisms, and their dysregulation can lead to various health problems. In this study, we propose a method for predicting hormonal peptides with high accuracy. The dataset used for training, testing, and evaluating our models consisted of 1174 hormonal and 1174 non-hormonal peptide sequences. Initially, we developed similarity-based methods utilizing BLAST and MERCI software. Although these similarity-based methods provided a high probability of correct prediction, they had limitations, such as no hits or prediction of limited sequences. To overcome these limitations, we further developed machine and deep learning-based models. Our logistic regression-based model achieved a maximum AUROC of 0.93 with an accuracy of 86% on an independent/validation dataset. To harness the power of similarity-based and machine learning-based models, we developed an ensemble method that achieved an AUROC of 0.96 with an accuracy of 89.79% and a Matthews correlation coefficient (MCC) of 0.8 on the validation set. To facilitate researchers in predicting and designing hormone peptides, we developed a web-based server called HOPPred. This server offers a unique feature that allows the identification of hormone-associated motifs within hormone peptides. The server can be accessed at: https://webs.iiitd.edu.in/raghava/hoppred/.Please read/cite the content about the HOPPred for complete information including algorithm behind the approach.

# Reference
Kaur, D., Arora, A., Vigneshwar, P., & Raghava, G. P. (2023). Prediction of peptide
hormones using an ensemble of machine learning and similarity-based methods. bioRxiv, (),
2023.05.15.540764. Accessed April 25, 2024. https://doi.org/10.1101/2023.05.15.540764.

# Standalone
The standalone version is available to download at https://webs.iiitd.edu.in/raghava/hoppred/hoppred.zip
The Standalone version of hoppred is written in python3 and following libraries are necessary for the successful run:

scikit-learn
Pandas
Numpy
For the successful run of the standalone, please download and install the latest version for BLAST from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/.

# Minimum Usage
To know about the available option for the stanadlone, type the following command:
```
python3 hoppred.py -h
```

To run the example, type the following command:
```
python3 hoppred.py -i peptide.fa
```
This will predict if the submitted sequences are hormonal peptides or non-hormonal peptides. It will use other parameters by default. It will save the output in "outfile.csv" in CSV (comma seperated variables).
```
usage: hoppred.py [-h] -i INPUT [-o OUTPUT] [-t THRESHOLD] [-m {1,2}] [-d {1,2}]
```
```
Please provide following arguments

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence in FASTA format or
                        single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.5
  -m {1,2}, -- model Model
                        Model: 1: AAC based RF, 2: Hybrid, by default 1

```
**Input File:** It allows users to provide input in two formats; i) FASTA format (standard) (e.g. peptide.fa) and ii) Simple Format. In case of simple format, file should have one one peptide sequence in a single line in single letter code (eg. peptide.seq). 

**Output File:** Program will save result in CSV format, in case user do not provide output file name, it will be stored in outfile.csv.

**Threshold:** User should provide threshold between 0 and 1, please note score is proportional to hormonal potential of peptide.

**Models:** In this program, two models have beeen incorporated:
**i) Model 1** for predicting given input peptide/protein sequence as Hormonal and Non-hormonal peptide/proteins using Logistic Regression based on top 50 features of the peptide/proteins; 

**ii) Model 2** for predicting given input peptide/protein sequence as Hormonal and non-hormonal peptide/proteins using Hybrid approach, which is the ensemble of Logistic Regression+ BLAST+ MERCI. It combines the scores generated from machine learning (LR), MERCI, and BLAST as Hybrid Score, and the prediction is based on Hybrid Score.

HOPPred Package Files
=======================
It contain following files, brief descript of these files given below

INSTALLATION  	: Installation instructions

LICENSE       	: License information

envfile : This file provide the path information for BLAST and MERCI commands,and data required to run BLAST and MERCI

Database: This folder contains the BLAST database and IgE motif files

progs : This folder contains the program to run MERCI

README.md     	: This file provide information about this package

hoppred.py 	: Main python program 

rf_model        : Model file required for running Machine-learning model

MERCI_motif_locator.pl                    : Perl script for locating motifs using MERCI

peptide.fa	: Example file contain peptide sequences in FASTA format

peptide.seq	: Example file contain peptide sequences in simple format

protein.fa	: Example file contain protein sequences in FASTA format 




