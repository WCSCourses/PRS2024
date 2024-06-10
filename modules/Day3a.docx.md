# Advanced Polygenic Risk Score Analyses

## Day 3 - Polygenic Risk Score Analyses Workshop 2024

## Table of Contents

  1. [Key Learning Outcomes](#Key-learning-outcomes)
  2. [Base and Target datasets](#Base-and-target-datasets)
  3. [Downloading Datasets](#Downloading-Datasets)
  4. [Method for Calculating PRS](#Method-for-calculating-PRS)
  5. [Exercise 1 Estimating R<sup>2</sup>](#exercise-1-estimating-r2)
  6. [Exercise 2 Visualising and comparing R<sup>2</sup>](#exercise-2-visualising-r2)
       
## Day 3a practical
## Key Learning Outcomes
After completing this practical, you should be able to:
  1. Compute and analyse ancestry-matched and unmatched PRS using PRSice-2.
  2. Understand and interpret the results from PRSise-2 derived scores. 
  3. Understand and identify the impact of ancestry on the predictive utility of PRS.
  4. Understand and identify the impact of sample size on the predictive utility of PRS.
  5. Understand the challenges and limitations of applying PRS in populations with diverse genetic backgrounds.


## Base and Target datasets 
In this practical, we will compute a PRS for systolic blood pressure (SBP) and assess it performance across European and African ancestry datasets to clearly illustrate the portability problem. 
We will assess the predictive utility of 3 scores:

**ANCESTRY-MATCHED**
1. EUR base - EUR target: Utilise European summary statistics as the training data and individual-level genotyped data from Europeans as the target dataset. 
2. AFR base - AFR target: Utilise African summary statistics as the training data and individual-level genotyped data from Africans as the target dataset.

**ANCESTRY-UNMATCHED**
1. EUR base - AFR target: Utilise European summary statistics as the training data and individual-level genotyped data from Africans as the target dataset. 

Please note that the sample sizes of the individual-level target data are as follows: 
Europeans (n = 500) and Africans (n = 650). 

* Note: This is simulated data with no real-life biological meaning or implication. 

|**Dataset**|**Source**|**Description**|
|:---:|:---:|:---:|
|Base dataset (EUR, AFR)|simulated |GWAS summary stats of SBP|
|Target dataset (EUR, AFR)|simulated |EUR (n = 500), AFR (n = 650)| 


## Downloading Datasets

The required datasets for this part of the practical can be downloaded at:       

[Download Backup WeTransfer](https://we.tl/t-YO8YbDjDK7)        

The data will be downloaded into your "Downloads" folder. You will need to move it to right directory, using the following commands.


```
cd data
mv ~/Downloads/Day_3.zip .
unzip Day_3.zip
cd Day_3
```

All required material for this practical is found in the **data/Day_3a** directory. Relevant materials that you should see there at the start of the practical are as follows:

 📂: Base_Data
  - GIANT_Height.txt,
  - cad.add.txt,
  - cad.add.readme.

 📂: Target_Data
  - TAR.fam
  - TAR.bim
  - TAR.bed
  - TAR.height
  - TAR.cad 
  - TAR.covariate
    
  📁: Reference files
   - Homo_sapiens.GRCh38.86.gt
   - Sets.gmt
     
 🛠️: Software
  - PRSice.R 
  - PRSice_linux
  - nagelkerke.R
  - Quantile.R

## Method for calculating PRS
For this practical we will use PRSice-2. PRSice-2 is one of the dedicated PRS calculation and analysis programs that makes use of a sequence of PLINK functions. The tools utilises the standard clumping and thresholding (C+T) approach.

## Exercise 1 Estimating R<sup>2</sup> 

### Key code parameters

The parameters listed in this table remain consistent across various scenarios, but the specific values may change based on the dataset and analysis scenario. 

For illustrative purposes, this table uses the first scenario, EUR base and EUR target population:

<details>
<summary>Click to view the parameters table</summary>

<table>
  <tr>
    <th>Parameter</th>
    <th>Value</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>prsice</td>
    <td>PRSice_xxx</td>
    <td>Informs PRSice.R that the location of the PRSice binary, xxx is the operating system (mac or linux)</td>
  </tr>
  <tr>
    <td>base</td>
    <td>EUR-SBP-simulated.sumstats.prscsx</td>
    <td>Specifies the GWAS summary statistics file for input</td>
  </tr>
  <tr>
    <td>A1</td>
    <td>A1</td>
    <td>Column name for the effect allele in the GWAS summary statistics</td>
  </tr>
  <tr>
    <td>p value</td>
    <td>P</td>
    <td>Column name for the p-values of SNPs in the GWAS summary statistics</td>
  </tr>
  <tr>
    <td>no clump</td>
    <td>-</td>
    <td>Instructs PRSice to skip the clumping process, which is used to remove SNPs in linkage disequilibrium</td>
  </tr>
  <tr>
    <td>beta</td>
    <td>-</td>
    <td>Indicates that the effect sizes are given in beta coefficients (linear regression coefficients)</td>
  </tr>
  <tr>
    <td>snp</td>
    <td>SNP</td>
    <td>Column name for SNP identifiers in the GWAS summary statistics</td>
  </tr>
  <tr>
    <td>score</td>
    <td>sum</td>
    <td>Specifies that the score calculation should sum the product of SNP effect sizes and their genotype counts</td>
  </tr>
  <tr>
    <td>target</td>
    <td>EUR_1kg.hm3.only.csx</td>
    <td>Specifies the genotype data file for the target sample</td>
  </tr>
  <tr>
    <td>binary-target</td>
    <td>F</td>
    <td>Indicates that the phenotype of interest is not binary (e.g., quantitative trait), F for no</td>
  </tr>
  <tr>
    <td>pheno</td>
    <td>sbp_eur_1kg.sim_pheno</td>
    <td>Specifies the file containing phenotype data</td>
  </tr>
  <tr>
    <td>pheno-col</td>
    <td>pheno100</td>
    <td>Column name in the phenotype file that contains the phenotype data to be used for PRS calculation</td>
  </tr>
  <tr>
    <td>cov*</td>
    <td>EUR.covariate</td>
    <td>Specifies the file containing covariate data for the analysis</td>
  </tr>
  <tr>
    <td>base-maf*</td>
    <td>MAF:0.01</td>
    <td>Filter out SNPs with MAF < 0.01 in the GWAS summary statistics, using information in the MAF column</td>
  </tr>
  <tr>
    <td>base-info*</td>
    <td>INFO:0.8</td>
    <td>Filter out SNPs with INFO < 0.8 in the GWAS summary statistics, using information in the INFO column</td>
  </tr>
  <tr>
    <td>stat*</td>
    <td>OR</td>
    <td>Column name for the odds ratio (effect size) in the GWAS summary statistics</td>
  </tr>
  <tr>
    <td>or*</td>
    <td>-</td>
    <td>Inform PRSice that the effect size is an Odd Ratio</td>
  </tr>
  <tr>
    <td>thread</td>
    <td>8</td>
    <td>Specifies the number of computing threads to use for the analysis</td>
  </tr>
  <tr>
    <td>out</td>
    <td>SBP_trial.eur.eur</td>
    <td>Specifies the name for the output files generated by PRSice</td>
  </tr>
</table>

*Note: These parameters are not used within this exercise but will likely be included when conducting your own analyses.

</details>

### Code

#### Scenario 1: Predicting from EUR training to EUR target data
First create the output directory
```
mkdir /home/manager/data/Data_Day4/data/out
```


```sh
Rscript /home/manager/PRSice_linux/PRSice.R \
--prsice /home/manager/PRSice_linux/PRSice \
--base /home/manager/data/Data_Day4/data/EUR-SBP-simulated.sumstats.prscsx \
--A1 A1 \
--pvalue P \
--no-clump \
--beta \
--snp SNP \
--score sum \
--target /home/manager/data/Data_Day4/data/EUR_1kg.hm3.only.csx \
--binary-target F \
--pheno /home/manager/data/Data_Day4/data/sbp_eur_1kg.sim_pheno \
--pheno-col pheno100 \
--thread 8 \
--out /home/manager/data/Data_Day4/data/out/SBP_trial.eur.eur  
```

<details>
  <summary>How many files have been generated from this code and what does each file show you?</summary>
  This code generates xxx?? files. 
  Each files serves a different purpose in the analysis and interpretation of derived PRS.
  These are outlined below: 
  
  1. Summary Statistics File: Provides a high-level overview of the PRS analysis results, allowing for a quick assessment of the PRS model's performance.
  2. Log File: Useful for debugging and detailed tracking of the computational steps undertaken during the PRS calculation.
  3. High-resolution Plots: Assist in visually assessing the sensitivity of PRS performance to the choice of P-value cutoffs and in selecting the optimal model.
  4. Best-fit PRS File: Essential for understanding which genetic variants contribute most significantly to the trait under study in the specific dataset used.
  5. Score Files: Allow for the practical application of the PRS in predicting or studying the trait within the study population; these scores are what you would correlate with actual phenotypic outcomes to validate or   use the PRS.
</details>


Examine the bar plot indicating the R<sup>2</sup> at each p-value threshold.

**Figure 1.1: BARPLOT generated by PRSice**  
```sh
![Figure 1.1](/images/Day2.docxfolder/images-018.png)
```
<details>
  <summary>Which P-value threshold generates the "best-fit" PRS?</summary>
  P-value threshold of xxx.
</details>

<details> 
  <summary>What does "best-fit" mean?</summary>
"Best-fit" refers to the P-value threshold at which the PRS accounts for the highest proportion of variance in the phenotype compared to other thresholds tested, suggesting an optimal balance between including informative SNPs and excluding noise from less relevant variants.
</details>

<details> 
<summary>Why does this matter?</summary>
Choosing the optimal P-value threshold is crucial because it affects the sensitivity and specificity of the PRS. A threshold that is too lenient (high P-value) might include too many SNPs, adding noise and possibly diluting the predictive power of the score. Conversely, a threshold that is too stringent (low P-value) might exclude potentially informative SNPs, reducing the ability of the PRS to capture the genetic architecture of the trait.
</details>

<details>
<summary>How much phenotypic variation does the "best-fit" PRS explain? What does this mean in very simple terms?</summary>
R<sup>2</sup> = 0.078 (7.8%).
This R<sup>2</sup> value means that out of the total variability observed in the trait across the population (under study), 7.8% can be attributed to the genetic variants included in this PRS.
</details>

<details>
<summary>What is the significance of the association (p-value) between the "best-fit" PRS and trait?</summary>
The p-value is xxxx. 
A p-value below 0.05 indicates statistically significant evidence that the PRS at this threshold significantly explains phenotypic variance and captures genuine genetic associations with the phenotype.
</details>

<details>
<summary>Which file provides a summary of the "best-fit" PRS?</summary>
SBP_trial.eur.eur.summary
</details>

View this summary output file: 
```sh
cat /Volumes/Chris-1/PRS24/Data_Day4/out/SBP_trial.eur.eur.summary
```

<details>
<summary>How many SNPs are included in the "best-fit" PRS? </summary>
Number of SNPs = xxx.
</details>  

<details>
<summary>Why is the number of SNPs important in PRS calcualtion? </summary>
The number of SNPs is crucial in PRS calculations because it directly influences the accuracy and predictive power of the score. 
Including a larger number of informative SNPs generally increases the likelihood of capturing significant genetic effects associated with the trait. However, it's important to balance quantity with quality; too many SNPs might include noise, which can dilute the true genetic signal and decrease the overall effectiveness of the PRS in distinguishing phenotypic outcomes.
</details>  

#### Scenario 2: Predicting from AFR training to AFR target data
```sh
Rscript /home/manager/PRSice_linux/PRSice.R \
--prsice /home/manager/PRSice_linux/PRSice \
--base /home/manager/data/Data_Day4/data/AFR-SBP-simulated.sumstats.prscsx \
--A1 A1 \
--pvalue P \
--no-clump \
--beta \
--snp SNP \
--score sum \
--target /home/manager/data/Data_Day4/data/AFR_1kg.hm3.only.csx \
--binary-target F \
--pheno /home/manager/data/Data_Day4/data/sbp_afr_1kg.sim_pheno \
--pheno-col pheno50 \
--thread 8 \
--out /home/manager/data/Data_Day4/out/SBP_trial.afr.afr
```
View the output file:
```sh
cat /home/manager/data/Data_Day4/out/SBP_trial.afr.afr.summary
```

<details>
  <summary>Which P-value threshold generates the "best-fit" PRS?</summary>
  ANSWER.
</details>

<details>
  <summary>How many SNPs are included in the "best-fit" PRS explain?</summary>
  Number of SNPs = xxx.
</details>

<details>
  <summary>How much phenotypic variation does the "best-fit" PRS explain?</summary>
  R<sup>2</sup> = 0.02 (2.0%).
</details>

#### Scenario 3: Predicting from EUR training to AFR target data
```sh
Rscript /home/manager/PRSice_linux/PRSice.R \
--prsice /home/manager/PRSice_linux/PRSice \
--base /home/manager/data/Data_Day4/data/EUR-SBP-simulated.sumstats.prscsx \
--A1 A1 \
--pvalue P \
--no-clump \
--beta \
--snp SNP \
--score sum \
--target /home/manager/data/Data_Day4/data/AFR_1kg.hm3.only.csx \
--binary-target F \
--pheno /home/manager/data/Data_Day4/data/sbp_afr_1kg.sim_pheno \
--pheno-col pheno50 \
--thread 8 \
--out /home/manager/data/Data_Day4/out/SBP_trial.afr.by.eur
```

View the output file: 
```sh
cat /home/manager/data/Data_Day4/out/SBP_trial.afr.eur.summary 
```

<details>
  <summary>Which P-value threshold generates the "best-fit" PRS?</summary>
  ANSWER.
</details>

<details>
  <summary>How many SNPs are included in the "best-fit" PRS explain?</summary>
  Number of SNPs = xxx.
</details>

<details>
  <summary>How much phenotypic variation does the "best-fit" PRS explain?</summary>
  R<sup>2</sup> = 0.078 (7.8%).
</details>


## Exercise 2 Visualising and comparing R<sup>2</sup>

In this exercise, we will analyse and compare the phenotypic variance explained (R<sup>2</sup>) by PRS across different combinations of base and target ancestries. 

Combine the summary files and visualise the performance of each PRS 
```sh
Rscript /home/manager/data/Data_Day4/software/R --- to check \

# Load necessary libraries
library(ggplot2)
library(RColorBrewer)

# Create a function to read the files and add ancestry information
read_and_label <- function(file, ancestry) {
  data <- read.table(file, header = TRUE, sep = "\t")
  data$Ancestry <- ancestry
  return(data)
}

# Read each file with the corresponding ancestry information
EUR_EUR <- read_and_label("SBP_trial.eur.eur.summary", "EUR_EUR")
AFR_AFR <- read_and_label("SBP_trial.afr.afr.summary", "AFR_AFR")
EUR_AFR <- read_and_label("SBP_trial.eur.by.afr.summary", "EUR_AFR")


# Combine all data into one dataframe
all_data <- rbind(EUR_EUR, AFR_AFR, EUR_AFR)

# Create a bar graph with different colors for each ancestry
ggplot(all_data, aes(x = Ancestry, y = PRS.R2, fill = Ancestry)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "R2 Values by Ancestry", x = "Ancestry", y = "R2 Value") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.position = "none")
 ```

Examine the bar plot indicating the R<sup>2</sup> for each base:target ancestry pair.

**Figure 2.1: BARPLOT - Phenotypic variance explained across ancestries **  
```sh
![Figure 2.1](/images/Day2.docxfolder/images-xx.png)
```

<details>
  <summary>Which base:target pair has the highest phenotypic variance explained?</summary>
  ANSWER.
</details>

<details>
  <summary>Which base:target pair has the lowest phenotypic variance explained?</summary>
  ANSWER.
</details>

<details>
  <summary>Explain the results?</summary>
   xxx.
</details>

<details>
  <summary>Are all results as expected? Explain why this may be the case.</summary>
  R<sup>2</sup> = 0.025 (2.5%).
</details>





> 
> ‼️ Note that all target phenotype data in this worshop are **simulated**. They have no specific biological meaning and are for demonstration purposes only. 
> 
---
<a href="#top">[Back to Top](#table-of-contents)</a>

