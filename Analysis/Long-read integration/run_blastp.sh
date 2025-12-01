#!/bin/sh
#SBATCH --partition=multinode
#SBATCH --constraint=[x2680|x2695]
#SBATCH --nodes=2
#SBATCH --ntasks=36
#SBATCH --ntasks-per-core=1
#SBATCH --time=144:00:00
#SBATCH --mem=108g

##sbatch run_talon.sh

module load samtools
module load TransDecoder

pwd

blastp -query transcripts.fasta.transdecoder_dir/longest_orfs.pep \
	-db ../uniref100/uniprot_sprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
	-num_threads 36 > blastp.outfmt6
