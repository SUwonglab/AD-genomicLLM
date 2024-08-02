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

Apart from the above softwares/packages, please also make sure the following softwares is installed properly: 1) [vcftools](https://vcftools.github.io/examples.html) for processing WGS .vcf file. 2) [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html) for genotype phasing. 3) [vcf2diploid](http://alleleseq.gersteinlab.org/tools.html) for constructing personal genome. 4) [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) for processing sMRI images. 

# Instructions
We provide detailed step-by-step instructions for running our pipeline.

## Processing genotype data

In our study, we downloaded whole genome sequencing (WGS) data from [ADNI database](https://adni.loni.usc.edu/), which provides .vcf file (gzip compressed) for each chromosome. Here, we take the chr19 and use gene APOE as a demonstration case study. 

**Step 1: remove indels**

```shell
vcftools --gzvcf ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr19.vcf.gz  --remove-indels --recode --recode-INFO-all --out SNPs_ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr19
[gzvcf] - input compressed vcf file
[out] - output file name (prefix)
```
**Step 2: genotype to haplotype**

```shell
java -jar preprocess/beagle.22Jul22.46e.jar gt=SNPs_ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr19.recode.vcf out=SNPs_ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr19.recode_hap map=plink.GRCh37.map/plink.chr19.GRCh37.map
[jar] - path the the beagle java program
[gt] - input vcf file from Step 1
[out] - output file name (prefix)
```
Note that the beagle .jar file is from [here](https://faculty.washington.edu/browning/beagle/beagle.html) and the plink full map files are from [here](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/).

**Step 3: personal genome construction**

```shell
java -jar vcf2diploid_v0.2.6a/vcf2diploid.jar -outDir fasta/chr19  -id  003_S_1057 -chr hg19/chr19.fa -vcf SNPs_ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr19.recode_hap.vcf.gz
[outDir] - output directory
[id] - personal ID
[chr] - reference genome for a chromosome
[vcf] - input vcf file from Step 2
```

Note that the above command can only construct the personal genome (both maternal and paternal) *per individual per chromosome*. For the construction of multiple individuals, the above command should be iterated over all individuals. The reference genome for a chromosome can be downloaded from [here](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/).

Above the above three steps, one should get `chr[CID]_[PID]_maternal.fa` and `chr[CID]_[PID]_paternal.fa` in the `fasta/chr[CID]` folder where `CID`,`PID` denotes chromosome ID and personal ID, respectively. 

In the above example case, it is `chr19_003_S_1057_maternal.fa` and `chr19_003_S_1057_paternal.fa`.



## Associating genotype to disease phenotype

**Step 1: Data Preparing**


