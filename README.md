# fastq2bam_snake

Simple snakemake pipeline to align high-throughput sequencing reads (in fastq) to a reference genome (final alignments in bam format).
Uses conda environments and snakemake wrappers, so the only dependencies are:

* conda
* snakemake

## Steps:

1. trim reads with Trimmomatic (default settings)
2. map reads to reference genome with BWA
3. filter unmapped or otherwise failing reads (with Samtools)
4. merge all mapped reads for a sample in a single bam file (with Samtools)
5. remove PCR and optical duplicated reads with Picard

Currently not functional:
6. run QC on raw and trimmed fastq, and final bams, and summarize results with MultiQC

## To make it work

* put necessary paths in config.yaml
* fill in sample information in the samples file (samples.tsv by default)

## To run

* locally: `snakemake --use-conda --cores N`, where N is the number of CPU cores to use
* on SLURM cluster (on the main node): `snakemake --use-conda -j N --printshellcmds --rerun-incomplete --latency-wait 30 --cluster-config cluster_config.yaml --cluster "sbatch --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem={cluster.mem} --partition={cluster.partition} --hint={cluster.hint}"`, where N is the number of concurrent jobs submitted and cluster_config.yaml contains resource info for jobs (ntasks, cpus-per-task, mem, partition, etc)
