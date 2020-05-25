#!/usr/bin/env nextflow

params.help = ""

if (params.help) {
  log.info ''
  log.info '-------------------------------------------------------'
  log.info 'NEXTFLOW 19.10 IMPLEMENT TCGA_MAF_Genelist-vs-Bootstrap'
  log.info '-------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run main.nf'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    -profile    Configuration profile (required: standard,singularity)'
  log.info '    --fastaHttp       STRING      link for download of preferred genome'
  log.info '    --genelists      STRING      path to set of .txt genelists; TSV format, two columns headed: Entrez_Gene_ID,	Gene_Name'
  log.info ''
  exit 1
}


##inputs
BASEDIR=$(pwd $0)
FASTA=$1
GATK4=$2
PICARDTOOLS=$3
RLIBPATH=$4
SAMTOOLS=$5
DRYRUN=$6

process vcf2maf_fa{

  output:
  file('fasta.fa') into fasta

  script:
  """
  wget -O fasta.fa ${params.fastaLink}
  wget -O maf2vcf.pl https://raw.githubusercontent.com/mskcc/vcf2maf/master/maf2vcf.pl
  """
}

process prep_fa {

  input:
  file(fa) from fasta

  output:
  tuple file('chr1-22-MT-X-Y.fa'), file('chr1-22-MT-X-Y.fa.fai'), file('chr1-22-MT-X-Y.fa.dict') into faidict
  
  script:
  """
  ##Making REF file (chr1-22, chrMT, chrX, chrY)"
  samtools faidx $fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrMT chrX chrY > chr1-22-MT-X-Y.fa
  samtools faidx chr1-22-MT-X-Y.fa
  samtools dict chr1-22-MT-X-Y.fa > chr1-22-MT-X-Y.fa.dict
  """
}
GENELISTS=$(ls data/geneLists/*.txt | perl -ane 'chomp;print "$_ "' | perl -ane '$sc=scalar(@F);for($i=0;$i<@F;$i++){$j=$i+1;if($j==$sc){print "$F[0]\n"}else{print $F[0] . ",";}}')


REFD=$REF".dict"
GDC=$(ls data/gdc*gz | sed 's/\.tar.gz//')
if [[ -d $GDC ]];then rm -rf $GDC;fi

##open all data into dirs, renaming by TCGA_XYZ disease types
echo "Untar GDC data"
mkdir $GDC
tar -C $GDC -xf $GDC".tar.gz"
rm $GDC"/MANIFEST.txt"

for DIR in $(ls $GDC); do

  FILEGZ=$GDC/$DIR/*gz;
  FILEMAF=$(echo $FILEGZ | sed 's/\.gz//');
  NAME=$(ls $FILEGZ | perl -ane '@s=split(/\//);@s1=split(/\./,$s[-1]);print "$s1[0].$s1[1].$s1[2]\n";');

  if [[ $NAME =~ "TCGA" ]]; then
    if [[ -d data/$NAME ]]; then
      echo "data/$NAME extant, please remove to rerun"
      exit 127
    fi
    echo "Working on: "$NAME;
    mkdir data/$NAME;
    gunzip $FILEGZ
    mv $FILEMAF data/$NAME/$NAME".somatic.maf";
    perl scripts/maf2vcf.pl \
      --input-maf data/$NAME/$NAME".somatic.maf" \
      --output-dir data/$NAME \
      --ref-fasta $REF

    java -jar $PICARDTOOLS SortVcf \
      I=data/$NAME/$NAME".somatic.vcf" \
      O=data/$NAME/$NAME".somatic.sort.vcf" \
      SD=$REFD

    tail -n +2 data/$NAME/$NAME".somatic.pairs.tsv" | \
      cut -f 1 > data/$NAME/$NAME".somatic.tumour.args"

    $GATK4 SelectVariants \
      -R $REF \
      -V data/$NAME/$NAME".somatic.sort.vcf" \
      -O data/$NAME/$NAME".somatic.tumour.vcf" \
      -sn data/$NAME/$NAME".somatic.tumour.args"

    if [[ $DRYRUN != "" ]]; then
      #R run
      Rscript --vanilla scripts/GRanges_geneLists_test.R \
        $RLIBPATH \
        $BASEDIR \
        data/$NAME \
        results/$NAME \
        data/geneLists \
        run
      Rscript --vanilla scripts/protein_domain_mutations.R \
        results/$NAME/${NAME}.somatic.tumour.vcf.glVcfGrList.RData

    fi
  fi
done
