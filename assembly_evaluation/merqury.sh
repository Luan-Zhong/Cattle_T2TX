#!/bin/bash

export MERQURY="/hpcfs/groups/phoenix-hpc-avsci/Lloyd_Low/Tuli_x_Wagyu_data/merqury/merqury_meryl_filt_issue/merqury"
export PATH=$PATH:/hpcfs/groups/phoenix-hpc-avsci/Lloyd_Low/Tuli_x_Wagyu_data/merqury/merqury_meryl_filt_issue/meryl-1.4.1/bin

hap1=$1
hap2=$2
out=$3

child="merqury/meryl_new/child.meryl"
mat="merqury/meryl_new/maternal.inherited.gt5.meryl"
pat="merqury/meryl_new/paternal.inherited.gt7.meryl"

#./merqury.sh <read-db.meryl> [<mat.meryl> <pat.meryl>] <asm1.fasta> [asm2.fasta] <out>

echo "Starting merqury..."
echo -e "bash $MERQURY/merqury.sh $child $mat $pat $hap1 $hap2 $out"
bash $MERQURY/merqury.sh $child $mat $pat $hap1 $hap2 $out
echo "Done."
