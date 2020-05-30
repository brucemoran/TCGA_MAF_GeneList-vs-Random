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
  log.info '    -profile    Configuration profile (required: standard,conda)'
  log.info '    --faPath      STRING      path to bgzip\'d fasta, hg38 from VEP seems to work fine'
  log.info '    --geneLists      STRING      path to set of .txt genelists; TSV format, two columns headed: Entrez_Gene_ID,	Gene_Name'
  log.info ''
  exit 1
}

fasta = Channel.fromPath("$params.faPath")

process prep_fa {
  label 'low_mem'

  input:
  file(fa) from fasta

  output:
  tuple file('chr1-22-X-Y.fa.gz'), file('chr1-22-X-Y.fa.gz.fai'), file('chr1-22-X-Y.dict') into faidict

  script:
  """
  ##Making REF file (chr1-22, chrMT, chrX, chrY)
  CHRTEST=\$(gunzip -c $fa | head -n1)
  if [[ ! \$CHRTEST =~ "^>chr" ]];then
    samtools faidx $fa 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y | sed 's/>/>chr/g' > chr1-22-X-Y.fa
  else
    samtools faidx $fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > chr1-22-X-Y.fa
  fi

  bgzip chr1-22-X-Y.fa
  samtools faidx chr1-22-X-Y.fa.gz
  samtools dict chr1-22-X-Y.fa.gz > chr1-22-X-Y.dict
  """
}

//Java task memory allocation via task.memory
javaTaskmem = { it.replace(" GB", "g") }

process dl_maf {
  label 'low_mem'

  output:
  file("GDCdata/**/*maf.gz") into mafs

  script:
  """
  Rscript --vanilla ${workflow.projectDir}/bin/TCGAbiolinks.maf_dl.R
  """
}

mafs
    .flatten()
    .set { flat_mafs }

process prep_maf {
  label 'med_mem'
  publishDir "analysis/$tcga/maf", mode: "copy"

  input:
  each file(maf) from flat_mafs
  tuple file(fa), file(fai), file(dict) from faidict

  output:
  tuple val(tcga), file("${tcga}.mutect2.somatic.vcf") into genelistR

  script:
  taskmem = javaTaskmem("${task.memory}")
  tcga = "${maf}".split("\\.mutect")[0].replaceAll("\\.", "-")
  """
  gunzip -c $maf > $tcga".mutect2.maf"
  maf2vcf.pl \
    --input-maf $tcga".mutect2.maf" \
    --output-dir $tcga \
    --ref-fasta $fa

  picard SortVcf -Xmx$taskmem \
    I=./$tcga/$tcga".mutect2.vcf" \
    O=$tcga".mutect2.sort.vcf" \
    SD=$dict

  tail -n +2 $tcga/$tcga".mutect2.pairs.tsv" | \
    cut -f 1 > $tcga".mutect2.somatic.args"

  gatk SelectVariants \
    -R $fa \
    -V $tcga".mutect2.sort.vcf" \
    -O $tcga".mutect2.somatic.vcf" \
    -sn $tcga".mutect2.somatic.args"
  """
}

process run_Rs {
  label 'low_mem'
  publishDir "analysis/$tcga/", mode: "copy"
  maxForks 1
  //issue with biomart and multiple access

  input:
  tuple val(tcga), file(vcf) from genelistR

  output:
  file('*') into completed
  tuple val(tcga), file("results"), file("bootstraps") into results

  script:
  """
  Rscript --vanilla ${workflow.projectDir}/bin/GRanges_geneLists_test.call.R \
    ${workflow.projectDir}/bin/GRanges_geneLists_test.func.R \
    $tcga \
    $vcf \
    ${params.geneLists}

  Rscript --vanilla ${workflow.projectDir}/bin/protein_domain_mutations.call.R \
    ${workflow.projectDir}/bin/protein_domain_mutations.func.R \
    $tcga".glVcfGrList.RData" \
    ${params.geneLists}

  chmod a+x ./results/*
  chmod a+x ./bootstraps/*
  """
}

process pack_out {
  label 'low_mem'
  publishDir "analysis/$tcga/", mode: "copy"

  input:
  tuple val(tcga), file(results), file(bootstraps) from results

  output:
  file('*') into resulted

  script:
  """
  wget -O textToExcelXLSX.pl https://raw.githubusercontent.com/brucemoran/perl/master/XLSX/textToExcelXLSX.pl

  ##write XLSX fro SNV data
  RESTXT=\$(ls results/*txt)
  perl ./textToExcelXLSX.pl \
    \$RESTXT \
    $tcga".mutect2.genelists.SNV"

  ##write genelist data
  GLTYPE=\$(ls bootstraps | grep -v p_ | grep -v xlsx | \
    perl -ane 'chomp;@s=split(/\\mutect\\./); \$so=\$s[-1]; \$so=~s/\\_[0-9]*.txt\$//; print "\$so\\n";' | sort | uniq)

  for GL in \$GLTYPE; do
    GLTYPEALL=\$(ls bootstraps/\$GL*.txt)
    perl ./textToExcelXLSX.pl \
      \$GLTYPEALL \
      $tcga"."\$GL".bootstraps.SNV"
  done

  zip -r $tcga".mutect2.results.zip" *xlsx *padjBH.nm.pdf
  """
}
