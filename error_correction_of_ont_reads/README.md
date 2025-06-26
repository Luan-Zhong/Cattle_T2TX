The HERRO correction was done using: (from this github page: `https://github.com/lbcb-sci/herro`)

1. Preprocess the reads

        preprocess.sh all_TxW_ONT_R10_q10_l10k.fastq preprocessed_all_TxW_ONT_R10_q10_l10k 32 16

2. minimap2 alignment and batching

        scripts/create_batched_alignments.sh preprocessed_all_TxW_ONT_R10_q10_l10k read_ids.txt 16 dir_bath_aln

3. Run HERRO correction

        singularity run --nv herro/herro.sif inference --read-alns alignment_run1 -t 16 -m herro/model/model_v0.1.pt -b 16 ONT_ppr.fastq.gz ONT_herro.fasta

Refer to:

- Error correction of ONT reads in the Supplementary Information (Supplementary Methods)
