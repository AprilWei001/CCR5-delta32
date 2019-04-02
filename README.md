# CCR5-delta32

Fig1.m reads three data files, HazardDelta32Bootstrap.txt, HazardDelta32Observed.txt, MatchMAFSNPsForHWE.txt, which are also provided under the same depository.

HazardDelta32Bootstrap.txt contains the hazard rate per year for each year (age 41 to 77) estimated from 1000 bootstrap, each bootstrap occupies 37 consecuetive rows, and four columns. The columns are hazard rates of individuals who are in the heterozygous state (delta32/+),homozygous state (delta32/delta32 and +/+), delta32/delta32, and +/+.

HazardDelta32Observed.txt contains the observed hazard rate per year for each year (age 41 to 77) estimated from the data. There are 37 rows total, and four columns in each row. The genotypic order here is different from the one in the bootstrap file. The first column is for delta32/delta32, followed by homozygous state (delta32/delta32 and +/+) ,heterozygous state (delta32/+),and +/+.

MatchMAFSNPsForHWE.txt This file contains the SNPs whose MAF matches with delta32. Each row is one SNP. The first column is the SNP-ID, and followed by the numbers of individuals who are in the homozygous state for the minor allele, in the heterozygous state, and in the homozygous state for the major allele. 

CCR5FigS1.m reads files, AgeAndAncestryInformation.txt, CCRgenesVCF.txt, HazardDelta32Observed.txt, HazardBritishnonBritishData.txt, deathRateUK.txt, Delta32SurvivalProbCorrected.txt. Unfortunately, AgeAndAncestryInformation.txt, contains individual birth/death/ancestry information extract from the phenotype data from the UK Biobank, and CCRgenesVCF.txt contains the genotype information for each individual for SNPs in the CCR5 gene. These two files can be shared with the right permission from the UK Biobank. The other four files are provided here.

HazardBritishnonBritishData.txt contains the estimated hazard rate per year for each year (age 41 to 77) estimated from all British ancestry individuals (column one) and from all non-British ancestry individuals (column two).

deathRateUK.txt contains the average hazard rate per year in the general population in UK from 39 to 78.

Delta32SurvivalProbCorrected.txt contains the hazard rate per year for each year (age 41 to 77) after the correction. The ages are in the first row, followed by the hazards for delta32/delta32, delta32/+, and +/+ genotypes.

BootstrapHazard.m shows how we bootstrap the individuals, and how we estimate the hazard rate pear year for each year. This script reads in AgeAndAncestryInformation_version2.txt, AgeAndAncestryInformation.txt, CCRgenesVCF.txt, which can be shared with the right permission from the UK Biobank.
