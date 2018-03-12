# TCGA_MAF_GeneList-vs-Random
Test Frequency of SNVs in Gene Lists versus Random Gene Sets for Increased Mutation

## How-to
Download the open-access TCGA MAFs for your preferred variant caller (MuTect2) from [GDC portal](https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.access%22%2C%22value%22%3A%5B%22open%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22WXS%22%5D%7D%7D%5D%7D) and save into 'data' directory

The data for multiple projects (disease types) should be in a single file, and the test data included should be removed. 

Run from base project directory:
```
sh scripts/parse_TCGA_MAF-VCF-genes.sh \
   <gatk-4.0.2.1/gatk> \
   <picard-tools-2.5.0/picard.jar> \
   <R/x86_64-pc-linux-gnu-library/3.5> \
   <samtools-1.5/samtools> \
   create_R_output
```
N.B. that the arguments to the bash script must be full paths, these are indicative, YMMV with different versions

N.B.B. the final argument can be any string. Leaving blank means that the data is not processed by the R script, so only the VCFs are made, which is handuy if that's what you want to do!