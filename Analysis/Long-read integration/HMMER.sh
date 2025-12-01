#!/bin/bash
set -e
module load hmmer
hmmsearch --cpu $SLURM_CPUS_PER_TASK -E 1e-10 --domtblout pfam.domtblout /data/Choi_lung/scLongreads/TALON_workspace/pfam/Pfam-A.hmm \
	/data/Choi_lung/scLongreads/TALON_workspace/test10s/transcripts.fasta.transdecoder_dir/longest_orfs.pep
