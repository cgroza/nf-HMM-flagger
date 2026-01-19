process bam_to_fastq {
  container "docker://mobinasri/flagger:v1.2.0"
  input:
  tuple val(sample_name), path(sample_reads)
  output:
  tuple val(sample_name), path("${sample_reads.baseName}.fq.gz")

  script:
  """
  samtools fastq -@ ${task.cpus} ${sample_reads} | gzip > ${sample_reads.baseName}.fq.gz
  """
}

process make_dip_asm {
  container "docker://mobinasri/flagger:v1.2.0"
  input:
  tuple val(sample_name), path(hap1), path(hap2)

  output:
  tuple val(sample_name), path("${sample_name}_dip.fa.gz")

  script:
  """
  zcat ${hap1} ${hap2} | bgzip > "${sample_name}_dip.fa.gz"
  """
}
process map_asm {
  container 'library://cgroza/collection/graffite:latest'
  input:
  tuple val(sample_name), path(dip_asm), path(reads)

  output:
  tuple val(sample_name), path("${sample_name}.bam")

  script:
  """
  meryl count k=15 output merylDB ${dip_asm}
  meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

  winnowmap -t ${task.cpus} -W repetitive_k15.txt -ax map-pb -Y -L --eqx --cs -I8g ${dip_asm} ${reads} | \
    samtools sort -@ ${task.cpus} -OBAM -o ${sample_name}.bam -
  """
}

process deepvariant {
  container "docker://google/deepvariant:1.3.0"
  input:
  tuple val(sample_name), path(dip_asm), path(bam)
  output:
  tuple val(sample_name), path("${sample_name}_dip.vcf")

  script:
  """
  samtools faidx ${dip_asm}
  samtools index ${bam}
  /opt/deepvariant/bin/run_deepvariant \
    --model_type="PACBIO" \
    --ref="${dip_asm}" \
    --reads="${bam}" \
    --output_vcf="${sample_name}_dip_all.vcf" \
    --make_examples_extra_args="keep_supplementary_alignments=true,min_mapping_quality=0" \
    // --call_variants_extra_args="use_openvino=true" \
    --num_shards=${task.cpus/2} \
    --dry_run=false

  bcftools view -Ov -f PASS -m2 -M2 -v snps -e 'FORMAT/VAF<~{vafCutoff} | FORMAT/GQ<~{qCutoff}' ${sample_name}_dip_all.vcf > ${sample_name}_dip.vcf
  """
}

process filter_alt_reads {
  container "docker://mobinasri/flagger:v1.2.0"

  input:
  tuple val(sample_name), path(bam), path(snps_vcf)

  output:
  tuple val(sample_name), path("${sample_name}_filtered.bam")

  script:
  """
  filter_alt_reads \
    -i ${bam}" \
    -o "${sample_name}_filtered.bam" \
    -f "removed.bam" \
    -v "${snps_vcf}" \
    -t ${task.cpus} \
    -m 1000 \
    -r 0.4
  """
}

process run_flagger {
container "docker://mobinasri/flagger:v1.2.0"
  input:
  tuple val(sample_name), path(dip_asm), path(bam)
  output:
  tuple val(sample_name), path("${sample_name}_hmm_flagger_outputs")

"""
  samtools depth -aa -Q 0 ${INPUT_DIR}/read_alignment.bam > ${INPUT_DIR}/read_alignment.depth

  bam2cov --bam ${bam} \
    --output coverage_file.cov.gz \
    --threads ${task.cpus}

  hmm_flagger \
    --input coverage_file.cov.gz \
    --outputDir ${sample_name}_hmm_flagger_outputs  \
    --alphaTsv /home/programs/config/alpha_optimum_trunc_exp_gaussian_w_4000_n_50.tsv \
    --labelNames Err,Dup,Hap,Col \
    --threads ${task.cpus}
"""
}

workflow {
  Channel.fromPath(params.samples).splitCsv(header:true).map{row ->
    [row.sample, file(row.hap1, checkIfExists:true), file(row.hap2, checkIfExists:true), file(row.reads, checkIfExists:true)]}.set{input_ch}
  bam_to_fastq(input_ch.map{[it[0], it[3]]}).set{fq_ch}
  make_dip_asm(input_ch.map{[it[0], it[1], it[2]]}).set{dip_ch};
  map_asm(dip_ch.combine(fq_ch, by: 0)).set{bam_ch}
  deepvariant(dip_ch.combine(bam_ch, by: 0)).set{vcf_ch}
  filter_alt_reads(bam_ch.combine(vcf_ch, by: 0)).set{bam_filter_ch}
  run_flagger(dip_ch.combine(bam_filter_ch, by: 0))
}
