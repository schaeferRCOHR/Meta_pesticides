# Repo for the submitted preprint: How pesticide exposure and effects match with the intention of European pesticide regulation – a mini review

This repo contains the R code and data to reproduce the analysis

Written by Ralf B. Schäfer

## Content overview: 
  - Authorisation_information.csv: Authorisation status of the pesticides covered in Beaumelle et al. 2023
  
  - Biodiv_overview.csv: Table with the synthetic studies and cases on pesticide effects
  
  - Exposure_overview.csv: Table with larger spatial scale studies on comparisons of pesticide exposure predictions with measurements
  
  - Script.R: Code for the analysis
  
## Metadata for variables in the files
| Variable Name | Explanation Text |
| :--- | :--- |
| Authors | Authors of publication |
| Title	| Title of publication |
| Reason for exclusion	| NA = Study not excluded, otherwise reason given |
| Case_ID	| Assigned unique ID |
| Biome	| Aquatic or terrestrial |
| Organism group	| Type of organism group (e.g. animal, plant) |
| Response	| Type of response (e.g. mortality, reproduction) |
| Response level	| Population or community response level |
| Study_ID	| Unique ID per study |
| Number of pesticides	| Number of pesticides considered in study |
| Pesticide_identity	| Individual pesticide or type (e.g. herbicide) |
| Compartment |	Study compartment: Soil or water |
| Soil depth sampling (cm) |	Sampling depth |
| Number of independent sites | Number of independent sites in study |
| Region/Spatial unit |	Study region or spatial unit |
| Number of studies	| Number of studies evaluated in synthesis or meta-analysis |
| Number of observations |	Number of observations in study, synthesis or meta-analysis |
| Observation unit	| Analytical level at which exceedances are evaluated (e.g. level of sites or pesticides) |
| Regulatory threshold	| Type of threshold (e.g. PEC, RAC) |
| Response measure	| Type of response measure  |
| Maximum	| Maximum value of response measure |
| Qualifier Minimum	| Qualifier for minimum |
| Minimum	| Minimum value of response measure |
| Mean	| Mean value of response measure |
| Median	| Median value of response measure |
| Effect direction	| Direction (positive, neutral or negative) of effect |
| Unit	| Unit of response measure |
| Values_larger_threshold_abs 	| Number of values larger than threshold |
| Threshold	| Value of threshold |
| Values_larger_threshold_rel	| Percentage of values larger than threshold |
| Unit_rel	| Unit of values larger than threshold |
| Variance measure	| Type of variance measure (e.g. standard deviation) |
| Lower value	| Lower value of variance measure |
| Upper value	| Upper value of variance measure |
| Stat_Sign_Based_on_CI	| Confidence interval intersects with 0? |
| Only_authorised	| Only authorised pesticides in 10/2025 considered in cases |
