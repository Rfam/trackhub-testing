process fetch_models {
  output:
  path("Rfam.cm")

  """
  wget "http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
  gzip -d Rfam.cm
  """
}

process randomize {
  input:
  path(all)

  output:
  path("models/*.cm")

  """
  mkdir model-list models
  cmstat $all | grep -v '^#' | awk '{ print \$3 }' | shuf > all.list
  split --lines 100 --additional-suffix='.list' all.list model-list/
  find model-list -name '*.list' | parallel cmfetch -o models/{/}.cm -f $all {}
  """
}

process fetch_gca {
  tag { "$ucsc" }

  input:
  tuple val(species), val(ucsc), val(genome_accession)

  output:
  tuple val(ucsc), path("${ucsc}.fasta"), emit: genome
  tuple val(ucsc), path('stats'), emit: stats

  """
  datasets download genome accession \
    --exclude-genomic-cds \
    --exclude-gff3 \
    --exclude-protein \
    --exclude-rna \
    --assembly-level chromosome \
    --filename genome.zip \
    $genome_accession

  unzip genome.zip
  find ncbi_dataset/data/$genome_accession/ -name '$genome_accession*_genomic.fna' > genomes
  lines="\$(wc -l genomes | awk '{ print \$1 }')"
  if [[ "\$lines" -ne 1 ]]; then
    echo "Too many genomes"
    exit 1
  fi
  head -1 genomes | xargs -I {} mv {} '${ucsc}.fasta'
  esl-seqstat '${ucsc}.fasta' > stats
  """
}

process fetch_gcf {
  tag { "$ucsc" }

  input:
  tuple val(species), val(ucsc), val(genome_accession)

  output:
  tuple val(ucsc), path("${ucsc}.fasta"), emit: genome
  tuple val(ucsc), path('stats'), emit: stats

  """
  datasets download genome accession \
    --exclude-genomic-cds \
    --exclude-gff3 \
    --exclude-protein \
    --exclude-rna \
    --assembly-level chromosome \
    --filename genome.zip \
    $genome_accession
  unzip genome.zip
  find ncbi_dataset/data -name 'chr*.fna' | grep -v '.scaf.fna' | xargs cat > '${ucsc}.fasta'
  esl-seqstat '${ucsc}.fasta' > stats
  """
}

process split_genome {
  tag { "$ucsc" }

  input:
  tuple val(ucsc), path(genome)

  output:
  tuple val(ucsc), path('parts/*.fasta')

  """
  split-sequences --max-nucleotides 2e5 '$genome' parts/
  """
}

process compute_z_value {
  tag { "$ucsc" }

  input:
  tuple val(ucsc), path('stats')

  output:
  path('z_score')

  """
  echo -n "$ucsc," > z_score
  grep 'Total # residues: ' stats | cut -d: -f2 | awk '{ print "scale = 6; (2.0000 * " \$1 ") / 1000000" }' | bc >> z_score
  """
}

process cmsearch {
  tag { "$ucsc:$genome:$cms" }
  /* cpus 8 */
  /* memory 20.GB */

  input:
  tuple val(ucsc), path(genome), path(cms), val(z_value)

  output:
  tuple val(ucsc), path("${ucsc}.txt"), emit: output
  tuple val(ucsc), path("${ucsc}.tblout"), emit: tblout

  """
  cmsearch \
    -o ${ucsc}.txt \
    -Z $z_value \
    --cpu 8 \
    --tblout ${ucsc}.tblout \
    --cut_ga \
    --rfam \
    --nohmmonly \
    $cms \
    $genome
  """
}

process save_tblout {
  publishDir "$baseDir/results"

  input:
  tuple val(ucsc), path('raw*.tblout')

  output:
  tuple path("${ucsc}.tblout.gz"), path("${ucsc}.bed")

  """
  cat raw*.tblout > ${ucsc}-merged.tblout
  cmsearch-deoverlap.pl ${ucsc}-merged.tblout
  mv ${ucsc}-merged.tblout.deoverlapped ${ucsc}.tblout

  tblout2bigBed.pl ${ucsc}.tblout $ucsc
  mv ${ucsc}.tblout.bed ${ucsc}.bed
  gzip ${ucsc}.tblout
  """
}

process save_output {
  publishDir "$baseDir/results"

  input:
  tuple val(ucsc), path(output)

  output:
  path("${ucsc}.output.gz")

  """
  cat $output > ${ucsc}.output
  gzip ${ucsc}.output
  """
}

workflow build_track_hub {
  fetch_models | randomize | flatten | set { models }

  Channel.fromPath("genomes.txt") \
  | splitCsv \
  | branch {
    gca: it[2].startsWith("GCA")
    gcf: it[2].startsWith("GCF")
    other: true
  }
  | set { specs }

  specs.other \
  | collectFile(name: 'unknown-genomes.txt', newLine: true) \
  | subscribe {
    println "Unknown genomes are saved to file: $it"
    exit 1
  }

  specs.gca | fetch_gca
  specs.gcf | fetch_gcf

  fetch_gca.out.genome.mix(fetch_gcf.out.genome) \
  | split_genome \
  | transpose \
  | set { genome_chunks }

  fetch_gca.out.stats.mix(fetch_gcf.out.stats) \
  | compute_z_value \
  | splitCsv \
  | set { z_values }

  genome_chunks \
  | combine(models) \
  | join(z_values) \
  | cmsearch

  cmsearch.out.output | groupTuple | save_output
  cmsearch.out.tblout | groupTuple | save_tblout
}

workflow {
  build_track_hub()
}
