# Coalescent_Bf_institute_project
The repository of my 2nd semester project in Bioinfirmatics Institute.
Goal:
- To simulate soome trees to obtain a number of statictics and allele-frequency spectra for the coalescent with different amounts of selection.

Objectives:
- To find an apporpriate simulator (beta-tree)
- To modify it in a way to fit our goal (done)
- To simulate a bunch of trees and to calculate said statistics and AFSs with different amounts of selection (see pictures)

System requirements:
python 3 with Biopython, numpy and matplotlib

Instructions for running code:
Just run the script Maincode.py and eventualy you'll get all the same pictures I store in this repository.

The repository description:
It contains my code used for obtaining distributions of different statistics and allele-frequency spectra - Maincode.py (may be not the best name, but it's something I came up with at the time). It also contains the folder with betatree - a simulator written in Neher laboratory, which I downloaded from his repository -
https://github.com/neherlab/betatree
And modified it so that it now provides a possibility of switching between selection and neutrality periods. There really wasn't much to do, actually, so it's mostly original.
There are also pictures I got. The name of a file corresponds to a statistics it represent. Not surprisingly. On each picture green line represents full selection mode (Bolthausen-Sznitman coalescent), red line - full neutrality mode (Kingsman coalescent) and the rest - intermediate states with step 0.1, the darker, the more neutral.
