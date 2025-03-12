## README for: A Quantitative Survey of Cogradient and Countergradient Variation in Nature
Authors: Albecker, MA, Bittar, T, Trussell G., Lotterhos, KE. 

This document provides information for the data and code used in the manuscript, "A Quantitative Survey of Cogradient and Countergradient Variation in Nature" by Dr. Molly Albecker, Dr. Thais Bittar, Dr. Geoff Trussell and Dr. Katie Lotterhos. Author information is included in the manuscript. 

# Study Information
The purpose of this study was to apply a quantitative measure to phenotypic data that were collected by other studies to test hypotheses about natural prevalence and patterns of Covariance between genotypic and environmental effects on phenotypes (CovGE). 


# List of files and folders: 
1. Data file - "Extraction_Initialization.csv" - Meta data on studies selected for quantitative review
2. Data file - "Meta_data_checked.csv" - Compiled dataset with phenotypic, genotypic, and environmental data for each study included in quantitative review. 
3. Data file - "Meta_output_Nov27_2024.csv" - Result dataframe with estimates for CovGE and Delta_GxE for each phenotype/study
4. Code file - "CovGE_Meta_Analysis.R" - R code including all analyses, data formatting, and plotting. 
5. Code file - "CovarianceDataFunctions.R" - R companion code storing all essential functions for analyses.

# Description of Contents: 
There were three major phases to this study which are pertinent for data storage and collection. 

1. The first phase was a quantitative review in which the literature was systematically surveyed, assessed for inclusion, and relevant data extracted. This process resulted in the dataset "Extraction_Initialize.csv", which contains metadata from all studies initially considered for data collection. In some cases, studies were later deemed inappropriate and excluded from subsequent phases, but their metadata remains in this file. Once data collection was completed for the 102 studies that met the inclusion criteria, all phenotypic data—used in the second phase—were compiled into a single .csv file, "Meta_data_checked.csv". 

2. The second phase involved data analysis, where the phenotypic data and metadata were imported into R, formatted, and processed using a for-loop to estimate CovGE and ΔG×E for each phenotype independently. This loop fully implements CovGE estimation following the methods outlined in Albecker et al. (2022) and also generates confidence intervals for CovGE and p-values for ΔG×E. This analysis was executed in lines 1–220 of the "CovGE_Meta_Analysis.R" script, with essential functions housed in the companion script "CovarianceDataFunctions.R". Results of this analysis were saved as, "Meta_output_Nov27_2024.csv" 

3. The third phase involved testing the questions/hypotheses outlined in the manuscript (e.g., Is CnGV or CoGV more prevalent in nature?). These analyses were conducted in the "CovGE_Meta_Analysis.R" script, following the conclusion of the primary analysis on Line 220. This phase also included all plotting and formatting code. Any modifications made to plots in PowerPoint were documented with notes. All code is fully annotated for clarity and reproducibility.

# Variables Defined: 

## Extraction_Initialize: 
1. Study_ID - unique numeric identifier for each study
2. First_Author - name of the first author on the paper
3. Data_file_name - compiled column with study id, first name, and phenotype 4. G_match_E - Inclusion criteria was that studies must show that phenotypes from a particular genotype (G) were exposed to at least one environmental treatment that matched their native environment (E). G_match_E indicates the study did so. 
5. natural_env_type - Were environments presented in categorical or continuous terms
6. GxE_sig - did authors find significant GxE for this particular phenotype? NA is not reported, True is yes, False is no. 
7. phylum_division - taxonomic classification 
8. genus - taxonomic classification
9. species - taxonomic classification
10. gen_number - How many genotypes tested in the study
11. gen_number_knownenvs - if they studies more genotypes than they provided information for, how many genotypes have defined native environments? 
12. phenotype - what phenotype is being reported
13. phenotype_unit - what unit of measurement is phenotype
14. life_stage - what life stage is being studied
15. exp_env_type - are environmental treatments continous or categorical
16. experimental_comparison - briefly, what is being compared in study?
17. Experimental_environment - what are treatments manipulatign in study?
18. Experimental_env_unit - what units of measurement were experimental treatments
19. natural_env - What natural environment measure is being associated with experimental treatment 
20. natural_env_units - what units of measurement correspond to the natural environment
21. Study_ID_phenotype - identifier combining study id and phenotype
22. Raw_data - did the study provide raw data, means, box plots, what? 
23. Proportion - is the phenotypic data presented as a proportion? 
24. Mention_CGV - do authors mention CovGE in the paper? 
25. Data_source - where was data collected from in the paper? 
26. figure_number - if data grabbed from figure, which one?
27. table_Number - if data grabbed from table, which one? 
28. Comments - extra notes if necessary
29. Checked? - noted if someone QA/QC'd data for that phenotype

## : Meta_data_checked
1. Study_ID_phenotype - identifier combining study id and phenotype
2. First_Author - first name of author of study
3. gen_factor - Assigned name for genotype with "G_" and numeric
4. nat_env_factor - The native environment of each genotype (started with N_ and switched to E_ later - in code, all are switched to E_).
5. exp_env_factor - Assigned categorical label for each environmental treatment - one must correspond with a native environment. These all begin with "E_" and end with numeric. 
6. phen_n - How many samples were collected in each genotype/enviornment
7. phen_data - the data for the phenotype from the genotype in that environment. 
8. phen_sd - if means and standard deviation were collected, enter SD in this column. 
9. phen_mean_SE - if means and standard error were collected, enter SE in this column. 
10. phen_mean_lowCI_095 - if means and confidence intervals were collected, enter lower CI value in this column. 
11. phen_mean_highCI_095 - if means and confidence intervals were collected, enter upper CI value in this column. 
12. source - lists original, unique xlsx file 

## Meta_output_Nov27_2024
1. Artifact - numbered column which is artifact of creating dataframe in R
2. Index - Number for each phenotype output
3. Time - Amount of seconds to run analysis
4. Data_file_name - Unique identifier for each study/author/phenotype. Used to match across studies to pull in relevant metadata
5. First.Author - name of first author (mostly used to make sure no oddities in data formatting)
6. Data.type - were data presented as "raw" or "means"? 
7. Phylum - taxonomic classification
8. Genus - taxonomic classification
9. Species - taxonomic classification
10. Order - taxonomic classification
11. Thermy - thermal physiology - ecto/endotherm or plant
12. Total_Sample_Size - total number of samples across genotypes and environments
13. Phenotype - phenotype measured
14. Phenotype_units - units of phenotype measurements
15. Morphological.Trait - is the phenotype a morphological trait? T - yes, F = no
16. Growth.Rate - is the phenotype a growth-related trait? T - yes, F = no
17. Developmental.time - is the phenotype a development trait? T - yes, F = no
18. Metabolic.Physiology - is the phenotype related to metabolic trait? T - yes, F = no
19. Trait.class - derelict column, all NAs
20. Experimental_comparison - brief description of environmental comparison in study
21. Environment.type - into what of 9 environmental classifications does this comparison fit? 
22. Study_Continent - on what continent did this study occur? 
23. Study_Country - in what country did this study occur? Multiple countries marked as "multiple". 
24. Study_location - as narrow location as possible 
25. Mention.CGV - do authors mention CovGE in original paper? 
26. Sig.GxE - do authors find significant GxE in this phenotype in original paper? 
27. Covariance - CovGE estimate
28. Covariance_LCI - CovGE lower confidence interval
29. Covariance_UCI - CovGE upper confidence interval
30. Covariance_Pvalue - Signifiance (pvalue) of CovGE estimate
31. GxE_Estimate - Delta GxE estimate
32. GxE_LCI - Delta GxE lower confidence interval - not used because unreliable
33. GxE_UCI - Delta GxE upper confidence interval - not used because unreliable. 
34. GxE_Pvalue - Pvalue for significance of Delta GxE
35. Omega2 - Omega squared estimate (Variance partioning method to generate idea of how much variation is explained by GxE) - not used here.
36. Omega2_LCI - Omega Squared lower confidence interval - not used here
37. Omega2_UCI - Omega Squared upper confidence interval - not used here
38. Omega2_Pvalue - Omega Squared Pvalue - not used here
39. G_diff - average difference in phenotypes across genotypes - not used here

