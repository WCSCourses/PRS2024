## Day 3b practical

7. [Introduction to cross-ancestry PRS computation](#Cross-ancestry-PRS-computation)
  8. [Cross-ancestry PRS analysis using PRS-CSx](#Cross-ancestry-PRS-computation)


## Introduction to Cross-Ancestry PRS computation
Before starting the practical the following commands will need to be run from within your virtual machine:

(1) conda create -n "PRScsx" python=3.7

(2) conda activate PRScsx

(3) pip install scipy

(4) pip install h5py

 The goal of this practical is to provide you with basic understanding and experience of running the PRS-CSx software. After completing this practical, you should:
* Be able to perform cross-population descriptives.
* Be familiar with running cross-ancestry PRS analyses using PRS-CSx.
* Understand how to evaluate linear models using Akaike’s Information Criterion

#### 1. The 1000 Genomes datasetS
The data we will be working with are coming from the 1000 Genomes Project reference panel. The data relates to individuals from 26 different source populations around the world. For simplicity, the populations have been collapsed into 5 broader continental super-populations: East Asian, European, South Asian, Amerindian, African ((EAS, EUR, SAS, EUR and AFR)). The scripts used to download and process the 1000Genomes data for the purposes of this course will be provided in the course appendix at the end of this week. 

#### 2. Cross-population allele frequency
Genetic variation is conveyed using allelic frequencies. Allele frequency is shaped by evolutionary forces and drift.  Here we compare profiles of allele frequency across the five ancestral populations. Global differences in allelic frequency has important implications for the portability of PRS across populations. Using plink it is possible to generate allele frequency statistics for each SNP, across populations, using the annotations provided in the file pop_info.pheno. In _/home/manager/data/Data_Day4_:
```sh
./software/plink_linux --bfile ./data/chr1-22 --freq --within ./data/pop_info.pheno
```
Population-stratified allele frequencies are reported in the output file plink.frq.strat.
For each population, print the numbers of total SNPs to screen, as follows:
**AFR**
```sh
grep -F 'AFR' plink.frq.strat | wc -l
```
From there we can print the number of SNPs with minor allele frequencies greater than 0 (and are hence potentially available for genomic analyes).
```sh
grep -F 'AFR' plink.frq.strat | awk '$6 >0' | wc -l
```

**EUR**
```sh
grep -F 'EUR' plink.frq.strat | wc -l
```
Number of SNPs with MAF > 0 in EUR.
```sh
grep -F 'EUR' plink.frq.strat | awk '$6 >0' | wc -l
```

**EAS**
```sh
grep -F 'EAS' plink.frq.strat | wc -l
```
Number of SNPs with MAF > 0 in EAS.
```sh
grep -F 'EAS' plink.frq.strat | awk '$6 >0' | wc -l
```
**SAS**
```sh
grep -F 'SAS' plink.frq.strat | wc -l
```
Number of SNPs with MAF > 0 in SAS.
```sh
grep -F 'SAS' plink.frq.strat | awk '$6 >0' | wc -l
```
**AFR**
```sh
grep -F 'AFR' plink.frq.strat | wc -l
```
Number of SNPs with MAF > 0 in EAS.
```sh
grep -F 'AFR' plink.frq.strat | awk '$6 >0' | wc -l
```

#### **Questions**
##### (i) Which population contains the most SNPs?
##### (ii) What  is the significance of the observed population order?  
&nbsp;
#### 3. Distribution of allele frequencies
```sh
R
library(dplyr)
library(ggplot2)
freq <-read.table("plink.frq.strat", header =T)
plotDat <- freq %>%
  mutate(AlleleFrequency = cut(MAF, seq(0, 1, 0.25))) %>%
  group_by(AlleleFrequency, CLST) %>%
  summarise(FractionOfSNPs = n()/nrow(freq) * 100)

ggplot(na.omit(plotDat),
       aes(AlleleFrequency, FractionOfSNPs, group = CLST, col = CLST)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 12)) +
  ggtitle("Distribution of allele frequency across genome")
```
#### **Questions**
##### (i) Which population has the most SNPs?
##### (ii) What  is the significance of the observed population ordering?
##### (iii) What is the reason behind these two features?
&nbsp;
#### 4. Calculation of Fst
Fst is a formal metric which is used to convey the level of genetic divergence between populations (on a scale between 0 and 1), using information derived from a set of genome-wide and mutually independent SNPs. Fst between parwise populations is estimated efficiently in Plink-2. So first we need to download Plink-2 using the command below:
```sh
sudo apt install plink2
```
A higher Fst corresponds to greater divergence between populations. Use the following command to calculate Fst:
```sh
plink2 --bfile ./data/chr1-22 --fst POP --pheno ./data/pop_info.pheno
```
Check the output file plink2.fst.summary
#### **Questions**
##### (i) Which population pairs have the highest Fst ?
##### (ii) For which populations is Fst smallest??
&nbsp;
## Introduction to PRS-CSx
#### 5. Background to PRS-CSX
PRS-CSx is a Python based command line tool that integrates GWAS summary statistics and external LD reference panels from multiple populations to improve cross-population polygenic prediction. We will be using simulated trait data pertaininng to systolic blood pressure (SBP) to explore PRS performance using 2 target populations that consist of 650 individuals of African ancestry and 500 individuals of European ancestry. Please note when running PRSice that the object of the flag "--prsice" will change according to whether plink is being called within the linux-like environment of the virtual machine (PRSice_linux) or a mac (PRSice_mac). Both executables can be found in the _/home/manager/data/Data_Day4_ directory. 

#### 6. Extraction of SNPs
PRS-CSx uses only HAPMAP3 SNPs therefore we produce a set of plink files containing this SNP subset.
```sh
./software/plink_linux --bfile ./data/1kg.eur.dbSNP153.hg38 --extract ./data/csxsnp --make-bed --out ./data/EUR_1kg.hm3.only.csx

./software/plink_linux --bfile ./data/1kg.afr.dbSNP153.hg38 --extract ./data/csxsnp --make-bed --out ./data/AFR_1kg.hm3.only.csx
```

#### 7. Running PRS-CSx
To model the coupling of effect sizes at individual SNPs across ancestries PRS-CSx uses an MCMC (Bayesian) sampling algorithm to determine values of the global shrinkage parameter ("phi") by Maximum likelihood. For samples of mixed or admixed genetic ancestry (which ours are not) the optimal value of the shrinkage parameter is estimated autonomously from the data. Here we use the value of phi (1e-4), which is suggested by the software authors, given that our trait is simulated to have a relatively small number (N=110) causal variants, distributed genome-wide.
To save time, we will be running the analyses across chromosome 15, rather than the entire genome. The commands needed to run PRS-CSx are contained in two script files, located in /home/manager/data/Data_Day4/scripts. The file run_prscsx_afr-target.sh is used to estimate optimal SNP weights for predicting into the African target population. 

NB - move the snp file to the reference folder

```
cp /home/manager/PRScsx/snpinfo_mult_1kg_hm3 /home/manager/data/Data_Day4/reference/csx    
```

NB - the ld block files must be moved to /home/manager/data/Data_Day4/reference/csx to /home/manager/PRScsx/

```
mv /home/manager/PRScsx/*.tar.gz /home/manager/PRScsx/
```
then extracted with tar

```
tar -xvfz ld...tar.gz
```
be careful of space, your vm has 100GB cap. remove the .gz files once extracted


**Script contents**:
```sh
python /home/manager/data/Data_Day4/software/PRScsx.py \
--ref_dir=/home/manager/data/Data_Day4/reference/csx \
--bim_prefix=/home/manager/data/Data_Day4/data/AFR_1kg.hm3.only.csx \
--sst_file=/home/manager/data/Data_Day4/data/EUR-SBP-simulated.sumstats.prscsx,/home/manager/data/Data_Day4/data/AFR-SBP-simulated.sumstats.prscsx \
--n_gwas=25732,4855 \
--pop=EUR,AFR \
--chrom=15 \
--phi=1e-4 \
--out_dir=/home/manager/data/Data_Day4/out/csx \
--out_name=afr.target.csx
```
The file run_prscsx_eur-target.sh is used to estimate optimal SNP weights for predicting into the European target population:
**Script contents**:
```sh
python /home/manager/data/Data_Day4/software/PRScsx.py \
--ref_dir=/home/manager/data/Data_Day4/reference/csx \
--bim_prefix=/home/manager/data/Data_Day4/data/EUR_1kg.hm3.only.csx \
--sst_file=/home/manager/data/Data_Day4/data/EUR-SBP-simulated.sumstats.prscsx,/home/manager/data/Data_Day4/data/AFR-SBP-simulated.sumstats.prscsx \
--n_gwas=25732,4855 \
--pop=EUR,AFR \
--chrom=15 \
--phi=1e-4 \
--out_dir=/home/manager/data/Data_Day4/out/csx \
--out_name=eur.target.csx
```
Prior to running each script you will need to update the first line with a reference (i.e the pathname) to your home directory. From _/home/manager/data/Data_Day4_ the scripts can then be run as follows:
```sh
./scripts/run_prscsx_afr-target.sh
./scripts/run_prscsx_eur-target.sh
```
#### **Questions**
##### (i) How many results files do you see in the output directory?
##### (ii) What does each file correspond to?
##### (ii) In which column are the adjusted SNP weights contained?      
&nbsp;
#### 8. Processing
The next step would be to collate adjusted SNP weight information from multiple chromosomes and relabel the consolidated files for clarity. The following code works regardless of whether the preceding PRS-CSx analysis was performed for single or multiple chromosomes. In _/home/manager/data/Data_Day4/out/csx_
```sh
for file in afr.target.csx_AFR_*; do
cat $file >> posteriors.afr.by.afr
done

for file in afr.target.csx_EUR_*; do
cat $file >> posteriors.afr.by.eur
done

for file in eur.target.csx_AFR_*; do
cat $file >> posteriors.eur.by.afr
done

for file in eur.target.csx_EUR_*; do
cat $file >> posteriors.eur.by.eur
done
```

#### 9. Data processing in R
The next stage of data processing, done in R will create a new set of summary statistics files, containing the adjusting SNP weights from PRS-CSx. Staying in _/home/manager/data/Data_Day4/out/csx_ :
```sh
R
# Load posteriors effect sizes from PRS-CSx output

afr.afr<-read.table("posteriors.afr.by.afr", col.names=c("chr","SNP","bp","a1","a2","post.afr.afr"))
afr.by.eur<-read.table("posteriors.afr.by.eur", col.names=c("chr","SNP","bp","a1","a2","post.afr.by.eur"))
eur.eur<-read.table("posteriors.eur.by.eur", col.names=c("chr","SNP","bp","a1","a2","post.eur.eur"))
eur.by.afr<-read.table("posteriors.eur.by.afr", col.names=c("chr","SNP","bp","a1","a2","post.eur.by.afr"))

# Load original summary statistics
afr.SBP<-read.table("/home/manager/data/Data_Day4/data/AFR-SBP-simulated.sumstats.prscsx", header=T)
eur.SBP<-read.table("/home/manager/data/Data_Day4/data/EUR-SBP-simulated.sumstats.prscsx", header=T)

# Combine the posterior derived weights and summary statistics
eur.SBP.merge1<-merge(x=eur.SBP, y=eur.eur[c(2,6)], by="SNP", all.y=T)
eur.SBP.merge<-merge(x=eur.SBP.merge1, y=afr.by.eur[c(2,6)], by="SNP", all.y=T)
afr.SBP.merge1<-merge(x=afr.SBP, y=afr.afr[c(2,6)], by="SNP", all.y=T)
afr.SBP.merge<-merge(x=afr.SBP.merge1, y=eur.by.afr[c(2,6)], by="SNP", all.y=T)

# Save files
write.table(afr.SBP.merge, "/home/manager/data/Data_Day4/afr.SBP.posterior.sumstats", quote=F, row.names=F)
write.table(eur.SBP.merge, "/home/manager/data/Data_Day4/eur.SBP.posterior.sumstats", quote=F, row.names=F)

q()
```

#### 10. Use PRSice to apply the adjusted SNP effects to target phenotypes
We first need to create a list of SNPs to be used as input, based on the previous PRS-CSx analysis. This is done in location _/home/manager/data/Data_Day4_.
```sh
awk 'NR>1 {print $1}' afr.SBP.posterior.sumstats > snps.afr.posterior
awk 'NR>1 {print $1}' eur.SBP.posterior.sumstats > snps.eur.posterior
```
The next series of commands implement in-ancestry and cross-ancestry prediction of the systolic blood pressure phenotype using the PRS-CSx optimised SNP weights. **Do not forget** to exchange _PRSice_linux_ for _PRSice_mac_ if running the commands in a Mac environment. The phenotypic files sbp_afr_1kg.sim_pheno and sbp_eur_1kg.sim_pheno contain data on systolic blood pressure with simulated heritabilities that vary from 10%, 20%, 33%, 50% and 100%. To compensate for the fact that the adjusted weights generated are based on chromosome 15, rather than the entire genome, we will be using the **_100%_** trait version (pheno100) for this analysis.


#### 15. Create summary files
From location _/home/manager/data/Data_Day4/out/prscsx_prsice_
```sh
awk '{print $1,"\t",$4,"\t",$8,"\t",$10}' SBP.afr.afr.summary > quicksum.afr-afr
awk '{print $1,"\t",$4,"\t",$8,"\t",$10}' SBP.afr.by.eur.summary > quicksum.afr-by-eur
awk '{print $1,"\t",$4,"\t",$8,"\t",$10}' SBP.eur.by.afr.summary > quicksum.eur-by-afr
awk '{print $1,"\t",$4,"\t",$8,"\t",$10}' SBP.eur.eur.summary > quicksum.eur-eur
```
#### **Questions**
##### **(i) What information is being summarized in the output files??** 
##### **(ii) For each prediction model what is the R2 and  value???** 
&nbsp;

#### 16. Complete the remaining PRS analysis in R
##### Step 1: Load and prepare data
```sh 
R
# Load libraries
sudo apt install cmake
install.packages("AICcmodavg")
install.packages("fmsb")
library(AICcmodavg)
library(fmsb)

afr.afr<-read.table("/home/manager/data/Data_Day4/out/prscsx_prsice/SBP.afr.afr.best", header=T)
afr.eur<-read.table("/home/manager/data/Data_Day4/out/prscsx_prsice/SBP.afr.by.eur.best", header=T)
eur.eur<-read.table("/home/manager/data/Data_Day4/out/prscsx_prsice/SBP.eur.eur.best", header=T)
eur.afr<-read.table("/home/manager/data/Data_Day4/out/prscsx_prsice/SBP.eur.by.afr.best", header=T)

colnames(afr.afr)[4]<-"afr.afr"   # african prediction using african PRS
colnames(afr.eur)[4]<-"afr.eur"   # african prediction using european PRS
colnames(eur.eur)[4]<-"eur.eur" # european prediction using european PRS
colnames(eur.afr)[4]<-"eur.afr" # european prediction using african PRS

# Merge PRSs according to target ancestry and source population
combined.afr<-merge(x=afr.afr, y=afr.eur[c(2,4)], by="IID", all.x=T)
combined.eur<-merge(x=eur.eur, y=eur.afr[c(2,4)], by="IID", all.x=T)

# Rescale scores to have mean = ‘0’ and standard deviation = ‘1’
combined.afr[4:5] <- as.data.frame(scale(combined.afr[4:5]))
combined.eur[4:5] <- as.data.frame(scale(combined.eur[4:5]))

# Load phenotype data
sbp.eur<-read.table("/home/manager/data/Data_Day4/data/SBP_eur_1kg.sim_pheno", header=T)
sbp.afr<-read.table("/home/manager/data/Data_Day4/data/SBP_afr_1kg.sim_pheno", header=T)

# Merge phenotypes and PRS
combined.eur.pheno<-merge(x=combined.eur[-c(2,3)], y=sbp.eur[,c(2,3)], by="IID", all.x=T)
combined.afr.pheno<-merge(x=combined.afr[-c(2,3)], y=sbp.afr[,c(2,3)], by="IID", all.x=T)
```
##### Step 2: Model selection using Akaike’s Information Criterion
In the final step we want to determine whether the use of the multi-ancestry PRS formulation generated by PRS-CSx performs better at predicting systolic blood pressure, compared to the single-ancestry PRS formulation. To do this we use (AIC Akaike’s Information Criterion). AIC is used to assess the performance of a competing set of regression models. We use it to compare the performance of models that contain different numbers of predictors, given that R2 is inflated by the inclusion of redundant terms in a model. 

##### 16b. Model evaluation
```sh 
# EUROPEAN
model1.eur <- glm(pheno100 ~ eur.eur, data = combined.eur.pheno, family=gaussian)
model2.eur <- glm(pheno100 ~ eur.afr, data = combined.eur.pheno, family=gaussian)
model3.eur <- glm(pheno100 ~ eur.eur + eur.afr, data = combined.eur.pheno, family=gaussian)

models.eur <- list(model1.eur, model2.eur, model3.eur) # define model set
mod.names.eur <- c('eur.eur', 'eur.afr', 'eur.combined') # add model names
aictab(cand.set = models.eur, modnames = mod.names.eur) # calculate model AICs model

# AFRICAN
model1.afr <- glm(pheno100 ~ afr.afr, data = combined.afr.pheno, family=gaussian)
model2.afr <- glm(pheno100 ~ afr.eur, data = combined.afr.pheno, family=gaussian)
model3.afr <- glm(pheno100 ~ afr.afr + afr.eur, data = combined.afr.pheno, family=gaussian)

models.afr <- list(model1.afr, model2.afr, model3.afr) # define model set
mod.names.afr <- c('afr.afr', 'afr.eur', 'afr.combined') # add model names
aictab(cand.set = models.afr, modnames = mod.names.afr) # calculate model AICs

# NagelkerkeR2 calculations in R
# African target
NagelkerkeR2(model1.afr)
NagelkerkeR2(model2.afr)
NagelkerkeR2(model3.afr)

# European target
NagelkerkeR2(model1.eur)
NagelkerkeR2(model2.eur)
NagelkerkeR2(model3.eur)

For each ancestry which predictive model performs best and worst?
Does combining the two sets of adjusted scores perform consistently better than modelling each one separately?
```
#### **Questions**
##### **(i) For each ancestry which model performs best and which performs worst?**   
##### **(ii) Does linearly combining adjusted scores (European plus African) always result in better performance compared to single-ancestry prediction?**   


**New Day 3b**

**Set up PRS-CSX environment (Not needed)**
---------------
Before starting the practical the following fixes will need to be applied from within your virtual machine
```
# conda create -n "PRScsx" python=3.7
# conda activate PRScsx
# pip install scipy
# pip install h5py
```

**Step 1: Set up environment**
------------------------------
First change to the working directory with the data for this practical 
```sh
cd /home/manager/data/Data_Day4/data
```
Make a directory called **hm3_by_ancestry** within the data folder, and move a folder back out of the data folder

```sh
mkdir hm3_by_ancestry
cd ..
```

**AFR**
```
for chr in {1..22}; do \
/home/manager/data/Data_Day4/software/plink_linux \
	--bfile /home/manager/data/Data_Day4/data/AFR_1kg.hm3.only.csx \
	--chr $chr \
	--make-bed \
	--out /home/manager/data/Data_Day4/data/hm3_by_ancestry/AFR_1kg.hm3.chr${chr}_only.csx;
done
```

**EUR**
```
for chr in {1..22}; do \
/home/manager/data/Data_Day4/software/plink_linux \
	--bfile /home/manager/data/Data_Day4/data/EUR_1kg.hm3.only.csx \
	--chr $chr \
	--make-bed \
	--out /home/manager/data/Data_Day4/data/hm3_by_ancestry/EUR_1kg.hm3.chr${chr}_only.csx;
done
```

**Set up the necessary environment variables for threading and verify they are set correctly.**
```
export N_THREADS=10
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS
```
**Verify the variables are set**
```
echo $N_THREADS
echo $MKL_NUM_THREADS
echo $NUMEXPR_NUM_THREADS
echo $OMP_NUM_THREADS
```
<br>

Step 2: Run CSX. Derive new SNPs weights trained on European and African summary stats
--------------------------------------------------------------------------------------

**Generate job file containing the threaded PRScsx commands.**
```
for chr in {21..22}; do
  echo "python3 ./software/PRScsx.py \
	--ref_dir=./reference/csx \
	--bim_prefix=./data/hm3_by_ancestry/AFR_1kg.hm3.chr${chr}_only.csx \
	--sst_file=./data/sumstats_by_chr/EUR-SBP-simulated.sumstats.chr${chr},./data/sumstats_by_chr/AFR-SBP-simulated.sumstats.chr${chr} \
	--n_gwas=25732,4855 \
	--chrom=$chr \
	--n_iter=1000 \
	--n_burnin=500 \
	--thin=5 \
	--pop=EUR,AFR \
	--phi=1e-4 \
	--out_dir=./out/csx \
	--out_name=afr.target_chr${chr}.csx" >> multithread_job.sh
done
```
**Check the job file contains the correct commands.**
```
head multithread_job.sh
```
**Run the Job File with GNU Parallel:**
```
parallel --verbose --jobs $N_THREADS < multithread_job.sh
```
**If the above command does not work you can create an alternative version of the same script by running the following:**
```
./script/create_multithread.sh
```
<br>

Step 3: Combine CSX-derived SNP weights across chromosomes (Currently Excludes Chromosome 3)
--------------------------------------------------------------------------------------------

**Load necessary library**
```
library(dplyr)
```
**Define the path to the directory containing the PRS-CSX output files**
```
path <- "/Users/iyegbc01/Dropbox/StatGen_OReilly_group/PRS_Summer_School/PRS_Summer_school_Wellcome_Trust_2023/Day4_WCS_2024_draft/out/csx/"
```
**Define the ancestry you want to combine ("EUR" or "AFR")**
```
ancestry <- "EUR"
```
**Initialize an empty data frame to store the combined data**
```
combined_data <- data.frame()
```
**Loop through chromosomes 1 to 22, (currently excluding chromosome 3)**
```
for (chr in setdiff(1:22, 3)) {
  # Construct the file name
  file_name <- paste0("afr.target_chr", chr, ".csx_", ancestry, "_pst_eff_a1_b0.5_phi1e-04_chr", chr, ".txt")
  file_path <- file.path(path, file_name)
  
  # Check if file exists before reading
  if (file.exists(file_path)) {
	# Read the data from the file
	data <- read.table(file_path, header = FALSE, sep = "\t", col.names = c("CHR", "rsid", "pos", "ref", "alt", "beta"))
	
	# Combine the data
	combined_data <- rbind(combined_data, data)
  } else {
	warning(paste("File not found:", file_path))
  }
}

```
**Write the combined data to a new file**
```
output_file <- file.path(path, paste0("combined_", ancestry, "_pst_eff.txt"))
write.table(combined_data, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
```
### Task: Replace 'ancestry <- "EUR" ' with 'ancestry <- "AFR" ' and repeat the subsequent steps shown above
<br>

# Step 4: Merge genotype-phenotype data
-------------------------------------

**Prepare data**

```
# Load libraries
library(data.table)
library(ggplot2)
library(snpStats)

# Define the path to the directory containing the PLINK files and phenotypic data
plink_path <- "/Users/iyegbc01/Dropbox/StatGen_OReilly_group/PRS_Summer_School/PRS_Summer_school_Wellcome_Trust_2023/Day4_WCS_2024_draft/data/"

# Read PLINK files and phenotype data into R
bim_file <- file.path(plink_path, "AFR_1kg.hm3.only.csx.bim")
fam_file <- file.path(plink_path, "AFR_1kg.hm3.only.csx.fam")
bed_file <- file.path(plink_path, "AFR_1kg.hm3.only.csx.bed")
pheno_file <- file.path(plink_path, "sbp_afr_1kg.sim_pheno")

# Read the genotype data using snpStats
geno_data <- read.plink(bed_file, bim_file, fam_file)

# Extract SNP IDs from geno_data$map$snp.name
snp_ids <- geno_data$map$snp.name
if (is.null(snp_ids) || !is.character(snp_ids)) {
  stop("SNP IDs are missing or not in the correct format.")
}

# Convert the genotype data to a matrix and then to a data.table
geno_matrix <- as(geno_data$genotypes, "matrix")
geno_df <- data.table(geno_matrix)
setnames(geno_df, snp_ids)

# Add IID column from the fam file
geno_df[, IID := geno_data$fam$member]

# Read the phenotypic data**
pheno_data <- fread(pheno_file, sep = " ", header = TRUE)
```

**Merge phenotype and genotype data**
```
# Do merge
combined_data <- merge(pheno_data, geno_df, by = "IID")

# Keep only genotype columns
geno <- combined_data[, !names(combined_data) %in% c("FID", "IID", "pheno100", "pheno20", "pheno33", "pheno10"), with = FALSE]
phen <- combined_data$pheno100

# Convert geno to numeric matrix**
geno <- as.matrix(geno)
mode(geno) <- "numeric"

# Convert phen to vector
phen <- as.vector(phen)
```
<br>

Step 5: Split data into validation and test sets 
------------------------------------------------
**Specify Proportion**
```
# Here we specify that 40% of all IDs will be used to construct the validation group
set.seed(154)
vali_proportion <- 0.4
vali_size <- round(nrow(geno) * vali_proportion)
vali_indices <- sample(1:nrow(geno), vali_size, replace = FALSE)
test_indices <- setdiff(1:nrow(geno), vali_indices)
```
**Subsetting of individuals**
```
X_vali <- geno[vali_indices, , drop=FALSE]
y_vali <- phen[vali_indices]
X_test <- geno[test_indices, , drop=FALSE]
y_test <- phen[test_indices]
```
<br>

Step 6: Prepare the regression model input using the CSX-derived AFR and EUR weights
-------------------------------------------------------------------------------------
```
# Read the merged CSX output files
AFR_betas <- fread(file.path(plink_path, "../out/csx/combined_AFR_pst_eff.txt"), sep = "\t", header = TRUE)
EUR_betas <- fread(file.path(plink_path, "../out/csx/combined_EUR_pst_eff.txt"), sep = "\t", header = TRUE)

# Assuming the beta files have columns: "CHR", "rsid", "pos", "ref", "alt", "beta"
overlap_prs <- merge(AFR_betas, EUR_betas, by = "rsid", suffixes = c("_afr", "_eur"))

# Filter overlap_prs to include only SNPs present in X_vali
overlap_prs <- overlap_prs[rsid %in% colnames(X_vali)]

# Ensure that X_vali and X_test only contain SNPs present in W_afr and W_eur
common_snps <- intersect(colnames(X_vali), overlap_prs$rsid)
X_vali <- X_vali[, common_snps, drop=FALSE]
X_test <- X_test[, common_snps, drop=FALSE]

# Reorder the columns of X_vali and X_test to match the order of SNPs in overlap_prs
X_vali <- X_vali[, match(overlap_prs$rsid, colnames(X_vali)), drop=FALSE]
X_test <- X_test[, match(overlap_prs$rsid, colnames(X_test)), drop=FALSE]

# Extract the overlapping rsid and their corresponding betas
W_afr <- overlap_prs$beta_afr
W_eur <- overlap_prs$beta_eur

# Convert W_afr and W_eur to numeric vectors
W_afr <- as.numeric(W_afr)
W_eur <- as.numeric(W_eur)
```
<br>

Step 7: Prepare the variant weights matrices as vectors
-------------------------------------------------------
```
# Pre-check the alignment between the different objects
if (ncol(X_vali) != length(W_afr) || ncol(X_vali) != length(W_eur)) {
  stop("Dimensions of X_vali and W_afr/W_eur do not match.")
}

# In the validation sample: 

# (i) Compute XWafr_vali
XWafr_vali <- X_vali %*% W_afr
# (ii) Convert XWafr_vali to have zero mean and unit variance
XWafr_vali_z <- scale(XWafr_vali)

# (iii) Compute XWeur_vali
XWeur_vali <- X_vali %*% W_eur
# (iv) Convert XWeur_vali to have zero mean and unit variance
XWeur_vali_z <- scale(XWeur_vali)

# (v) Combine the normalized matrices
XW_vali <- cbind(XWafr_vali_z, XWeur_vali_z)

# Fit the model
model <- lm(scale(y_vali) ~ XWafr_vali_z + XWeur_vali_z - 1) # '- 1' removes the intercept

# Obtain the regression parameters
a_hat <- coef(model)[1]
b_hat <- coef(model)[2]
print(paste("a_hat =", a_hat))
print(paste("b_hat =", b_hat))
```
<br>

Step 8: Predict phenotype on validation and test dataset
--------------------------------------------------------
**Generate a linear combination of AFR and EUR PRSs for each individual**
```
# Each ancestry component is weighted by the regression coefficient of that ancestry, in the preceding step
y_hat_vali <- a_hat * XWafr_vali_z + b_hat * XWeur_vali_z

# In the test sample: 

# Compute XWafr_test and XWeur_test
XWafr_test <- X_test %*% W_afr
XWafr_test_z <- scale(XWafr_test)

XWeur_test <- X_test %*% W_eur
XWeur_test_z <- scale(XWeur_test)

# y_hat in the test sample
y_hat <- a_hat * XWafr_test_z + b_hat * XWeur_test_z
```
<br>

Step 9: Plot phenotype distributions of validation and test data:
---------------------------------------------------------------------------
**Check that both distributions are approximately normal**
```
library(ggplot2)

# Create data frames for validation and test sets
vali_data <- data.frame(trait = y_vali, dataset = "Validation")
test_data <- data.frame(trait = y_test, dataset = "Test")

# Combine both data frames
combined_data <- rbind(vali_data, test_data)

# Plot the distributions
ggplot(combined_data, aes(x = trait, fill = dataset)) +
  geom_density(alpha = 0.5) +
  labs(title = "Trait Distributions for Validation and Test Sets", x = "Trait Value", y = "Density") +
  theme_minimal()
```
<br>

Step 10: Plot true values against predicted values
--------------------------------------------------
**The next steps use standard normal phenotype data to reduce scale differences between PRS and trait values**
```
min_true <- min(min(y_vali), min(y_test))
max_true <- max(max(y_vali), max(y_test))
min_pred <- min(min(y_hat_vali), min(y_hat))
max_pred <- max(max(y_hat_vali), max(y_hat))

pdf("true_against_pred.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

plot(scale(y_vali), y_hat_vali, pch = 19, col = rgb(0, 0, 0, 0.5), xlab = 'True Values', ylab = 'Predicted Values', main = 'Validation Dataset')
abline(0, 1, col = 'red', lty = 2)

plot(scale(y_test), y_hat, pch = 19, col = rgb(0, 0, 0, 0.5), xlab = 'True Values', ylab = 'Predicted Values', main = 'Test Dataset')
abline(0, 1, col = 'red', lty = 2)

dev.off()

```
Step 11: Calculate deviance-based _R<sup>2</sup>_
-------------------------------------
```
# Calculate the deviance (SS_res)
deviance <- sum((scale(y_test) - y_hat) ^ 2)
print(paste("deviance =", deviance))

# Calculate the mean of the scaled y_test
y_test_mean <- mean(scale(y_test))

# Calculate the null deviance (SS_tot)
deviance_null <- sum((scale(y_test) - y_test_mean)^ 2)
print(paste("deviance_null =", deviance_null))

# Calculate R2
R2 <- 1 - (deviance / deviance_null)
print(paste("R2 =", R2))
```














