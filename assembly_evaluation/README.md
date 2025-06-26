## Homologous chromosomes

To get the homologous chromosomes from a draft assembly using ARS-UCD2.0 as karyotyping reference

    #!/bin/bash
    
    export QC="scripts/assembly_initialqc"
    
    ref="REFERENCES/ARS-UCD_BLRC/ARS_UCD_v2.0.fa"
    map="assembly.scfmap"
    path="assembly.paths.tsv"
    
    echo "FOR HAPLOTYPE 1"
    qry1="assembly.haplotype1.fasta"
    dir1="hap1"
    "$QC"/assembly_initialqc.sh -r "$ref" -q "$qry1" -o "$dir1" -m "$map" -p "$path" -t "32"
    
    echo "DONE :)"

## Merqury QV

For quality and completeness

    bash merqury .sh

## QUAST

For contig N50 etc.

    quast.py UOA_Wagyu_1.fa -o quast_wagyu

## compleasm

For BUSCO completeness score

    bash compleasm.sh


Refer to:

- Assembly evaluation in the Methods of the main text
- Assembly evaluation in the Supplementary Information (Supplementary Methods)
