Author of this document: Xiao Zang
Author of the analyses: Rui Zhong and Xiao Zang

The following workflow for the primary screens assumes the raw data are stored in the folder /data and the information about the siRNA library is in GeneID.RData.

I. Run a separate model for each replicate plate
1. Preliminary analysis
Code to use: 
preliminary_analysis.R

File to use: 
GeneID.RData
/data


Explanation: 
The first step is to convert xlsx in /data to csv and store them in /csv. The second step is to analyze each library plate and output the analysis results to /preliminary_analysis. The third step is to combine the analysis result tables into preliminary_summary.csv.

Note that preliminary_analysis.R run Poisson regression without random effects, but the result we eventually used for screening is produced by the mixed effects models in first_analysis.R and second_analysis.R. However, the plate layout information and cell numbers summarized in preliminary_summary.csv is used in second_analysis.R for convenience.


2. The first analysis for the primary screening
Code to use:
first_analysis.R

Files to use:
GeneID.RData
/csv


Explanation:
Analysis is done using a mixed effects model in first_analysis.R. Output for each library plate is first written in /first_analysis and then combined to a first_analysis.csv.

II. The second analysis for the primary screening: Integrated model with interaction between viral infection and gene knockdown
Code to use:
second_analysis.R

Files to use:
preliminary_summary.csv
GeneID.RData
/csv

Explanation:
Analysis is done using a mixed effects model with interaction term. The result is output to second_analysis.csv

