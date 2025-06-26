The repeat library can be found here: `RepeatMaskerLib.h5.bostaurus-rm.withsats.fa`

The satellite repeats can be found here: `final_sats.fa`

Run repeatmasker:

    infile="wagyu_chr.list"
    lib="RepeatMasker/custom_lib/RepeatMaskerLib.h5.bostaurus-rm.withsats.fa"
    
    mapfile -t files < "$infile"
    file=${files[$SLURM_ARRAY_TASK_ID]}
    
    base=$(basename "$file" .fasta)
    echo -e "$file"
    /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/tools/repeatmasker/RepeatMasker/RepeatMasker -pa 24 -gff -lib "$lib" -dir rm_custom_lib "$file" -e ncbi

Clean the delimeters:

    tail -n +4 $1 | sed 's/[[:space:]]\+/\t/g' | sed 's/^\t\+//' > $(dirname $1)/$(basename $1).tsv

Remove overlaps:

    #!/bin/bash
    for i in rm_custom_lib/*fa.out; do
      RepeatMasker/util/RM2Bed.py --ovlp_resolution 'higher_score' --max_divergence 40 "$i" "$(dirname $i)"/"$(basename $i)".bed
    done

Futher analyse result using `final_repeatmasker_2.R`
