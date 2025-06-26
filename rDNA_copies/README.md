The rDNA sequence in NCBI DQ222453.1 were searched to the Wagyu genome using:

    #!/bin/bash
    ref="$1"
    qry="$(realpath "$2")"
    #makeblastdb -in "$ref" -dbtype nucl -parse_seqids
    blastn -db "$ref" -query "$qry" -outfmt 7 -evalue 1e-10 -perc_identity 85 -out out/"$(basename "$ref" .fa)"_vs_"$(basename "$qry" .fasta)".85.blast.tsv

Further analysis were done with `rdna_blast_2.R`
