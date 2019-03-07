Readme for datasets and analyses for the Harmonious Combinations paper

The data for three experiments and associated layout information are in three directories. The datasets for producing analyses and figures are located in each of the folders, and duplicated separately in the SupplementaryDataTables section. There is a readme in the SupplementaryDataTables directory explaining the headers and what they mean.

As of February 2019, this code worked with the current version of R and packages I have installed on my machine. Have a look at the sessionInfo.txt file for the versions of different packages I used for these analyses. 

Especially problematic might be executing flowCore-related packages. This package bricks any code implementing older versions of itself via radical updates every 6 months, so at some point I just got fed up and stopped updating my version of R and flowCore. As a result, my version of R is a bit old, so could lead to problems with your replication of some plots or analyses.

Dates of experiments are in YYMMDD format. 

You might be confused by the random pWXYZ names of the plots. At some point I become overwhelmed with copies of all basically-the-same-but-slightly-different plots, so to keep track of the versions I just make a 4-letter random code for the figure.


