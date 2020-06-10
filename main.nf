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

process anno_R {
  label 'low_mem'

  output:
  tuple file("biomaRt_all.extensent.RData"), file("EnsDb.Hsapiens.v86.txs.PFAM.db.RData") into anno_vcf

  script:
  """
  Rscript --vanilla ${workflow.projectDir}/bin/annotations.call.R \
    ${workflow.projectDir}/bin/annotations.func.R
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
  label 'high_mem'

  input:
  each file(maf) from flat_mafs
  tuple file(fa), file(fai), file(dict) from faidict

  output:
  tuple val(tcga), file("${tcga}.mutect2.somatic.vcf") into ( genelistR1, genelistR2 )

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

process prep_Rs {
  label 'high_mem'
  publishDir "analysis/$tcga/", mode: "copy", pattern: "!EnsDb.Hsapiens.v86.txs.PFAM.db.RData, !biomaRt_all.extensent.RData"

  input:
  tuple val(tcga), file(vcf) from genelistR1
  file(vcf_rdata) from anno_vcf

  output:
  tuple val(tcga), file("${tcga}.glVcfGrList.RData"),  file("EnsDb.Hsapiens.v86.txs.PFAM.db.RData"), file('./results'), file("bootstraps") into prepd
  tuple val(tcga), file("${tcga}.genelists.geom_density-padjBH.nm.pdf") into results_prep

  script:
  """
  Rscript --vanilla ${workflow.projectDir}/bin/GRanges_geneLists_test.call.R \
    ${workflow.projectDir}/bin/GRanges_geneLists_test.func.R \
    $tcga \
    $vcf \
    ${params.geneLists}

  chmod a+x ./results/*
  chmod a+x ./bootstraps/*
  """
}

prepd
    .join(genelistR2)
    .groupTuple()
    .map { it -> tuple(it[0], it[1..-1].flatten()).flatten() }
    .set { prepdtup }

process run_Rs {
  label 'high_mem'
  publishDir "analysis/$tcga/", mode: "copy", pattern: "!EnsDb.Hsapiens.v86.txs.PFAM.db.RData"

  input:
  tuple val(tcga), file(rdatagrv), file(rdataens), file(results), file(bootstraps), file(vcf) from prepdtup

  output:
  file('*') into completed
  tuple val(tcga), file("results"), file("bootstraps") into results_run

  script:
  """
  Rscript --vanilla ${workflow.projectDir}/bin/protein_domain_mutations.call.R \
    ${workflow.projectDir}/bin/protein_domain_mutations.func.R \
    $tcga".glVcfGrList.RData" \
    ${params.geneLists}

  chmod a+x ./results/*
  chmod a+x ./bootstraps/*
  """
}

results_run
    .join(results_prep)
    .groupTuple()
    .map { it -> tuple(it[0], it[1..-1].flatten()).flatten() }
    .set { results_pr }

process pack_out {
  label 'low_mem'
  publishDir "analysis/$tcga/", mode: "copy", pattern: "*.zip"

  input:
  tuple val(tcga), file(results), file(bootstraps), file(pdf) from results_pr

  output:
  file('*') into resulted

  script:
  """
  wget -O textToExcelXLSX.pl https://raw.githubusercontent.com/brucemoran/perl/master/XLSX/textToExcelXLSX.pl

  ##write XLSX fro SNV data
  RESTXT=\$(ls results/*txt)
  perl ./textToExcelXLSX.pl \
    \$RESTXT \
    $tcga".mutect2.genelists"

  ##write genelist data
  GLTYPE=\$(ls bootstraps | grep -v p_ | grep -v xlsx | \
    perl -ane 'chomp;@s=split(/\\mutect\\./); \$so=\$s[-1]; \$so=~s/\\_[0-9]*.txt\$//; print "\$so\\n";' | sort | uniq)

  for GL in \$GLTYPE; do
    GLTYPEALL=\$(ls bootstraps/\$GL*.txt)
    perl ./textToExcelXLSX.pl \
      \$GLTYPEALL \
      $tcga"."\$GL".bootstraps"
  done

  zip -r $tcga".mutect2.results.zip" *xlsx *padjBH.nm.pdf *.geom_density_ridges.pdf
  """
}
