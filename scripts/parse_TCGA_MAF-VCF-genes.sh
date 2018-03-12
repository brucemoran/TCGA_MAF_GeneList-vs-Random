#! /bin/bash

##inputs
GATK4=$1
PICARDTOOLS=$2
RLIBPATH=$3
SAMTOOLS=$4
DRYRUN=$5

##get vcf2maf
curl "https://raw.githubusercontent.com/mskcc/vcf2maf/master/maf2vcf.pl" > "scripts/maf2vcf.pl"

GENELISTS=$(ls data/geneLists/*.txt | perl -ane 'chomp;print "$_ "' | perl -ane '$sc=scalar(@F);for($i=0;$i<@F;$i++){$j=$i+1;if($j==$sc){print "$F[0]\n"}else{print $F[0] . ",";}}')

REFBASE=${HOME}"/.vep/homo_sapiens/91_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.chr.1.fa"
if [[ ! -e $REFBASE ]];then
    echo $REFBASE" not found, exiting"
    exit 127
fi
REF=$(echo $REFBASE | sed 's/fa/chr1-22-MT-X-Y.fa/')
if [[ ! -e $REF ]];then
    echo "Making REF file (chr1-22, chrMT, chrX, chrY)"
    $SAMTOOLS faidx $REFBASE chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
	      chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrMT chrX chrY > \
	      $REF
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
		    data/$NAME \
		    $GENELISTS
	fi

	rm -rf $GDC

    fi
done
