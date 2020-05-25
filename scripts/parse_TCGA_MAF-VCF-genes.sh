#! /bin/bash

##inputs
BASEDIR=$(pwd $0)
FASTA=$1
GATK4=$2
PICARDTOOLS=$3
RLIBPATH=$4
SAMTOOLS=$5
DRYRUN=$6

##get vcf2maf
curl "https://raw.githubusercontent.com/mskcc/vcf2maf/master/maf2vcf.pl" > "scripts/maf2vcf.pl"

GENELISTS=$(ls data/geneLists/*.txt | perl -ane 'chomp;print "$_ "' | perl -ane '$sc=scalar(@F);for($i=0;$i<@F;$i++){$j=$i+1;if($j==$sc){print "$F[0]\n"}else{print $F[0] . ",";}}')

if [[ ! -e $FASTA ]];then
  echo $FASTA" not found, exiting"
  exit 127
fi
REF=$(echo $FASTA | sed 's/{.fa,.fasta}/chr1-22-MT-X-Y.fa/')
if [[ ! -e $REF ]];then
  echo "Making REF file (chr1-22, chrMT, chrX, chrY)"
  $SAMTOOLS faidx $FASTA chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrMT chrX chrY > $REF
  $SAMTOOLS faidx $REF
  $SAMTOOLS dict $REF > $REF.dict
  DICT=$(echo $REF | sed 's/fa/dict/')
  ln -s $REF.dict $DICT
fi

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
      Rscript --vanilla scripts/annotate_protein_domains.R \
        results/$NAME/${NAME}.somatic.tumour.vcf.grPIDList.RData
    fi
  fi
done
