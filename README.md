# AD-genomicLLM
Associating genotype to imaging and clinical phenotypes of Alzheimer’s disease by leveraging genomic large language model

 ![model](https://github.com/SUwonglab/AD-genomicLLM/blob/main/workflow.png)

In this work, we propose a novel computational framework that leverages genomic large language models (LLMs) to enhance the association analysis between genetic variants and Alzheimer's disease (AD)-related phenotypes, including imaging and clinical features.

 # Enviornment
- Python==3.9.0
- java==1.8.0
- TensorFlow==2.12.0
- TensorFlow-hub==0.12.0
- pyfasta==0.5.2
- scikit-learn==1.0.2

Apart from the above softwares/packages, please also make sure the following softwares is installed properly: 1) [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) for processing sMRI images. 2) [vcftools](https://vcftools.github.io/examples.html) for processing WGS .vcf file.

# Instructions
We provide detailed step-by-step instructions for running our pipeline.

## Processing genotype data

In our study, we downloaded whole genome sequencing (WGS) data from [ADNI database](https://adni.loni.usc.edu/), which provides .vcf file (gzip compressed) for each chromosome. Here, we take the chr19 and use APOE as a demonstration case study. 

**Step 1: remove indels**

```shell
vcftools --gzvcf ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr19.vcf.gz  --remove-indels --recode --recode-INFO-all --out SNPs_ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr19
```
**Step 2: genotype to haplotype**





## Associating genotype to disease phenotype

**Step 1: Data Preparing**


