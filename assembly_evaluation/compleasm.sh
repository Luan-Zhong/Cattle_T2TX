#!/bin/bash

file="UOA_Wagyu_1.fa"
dir=.
odb="compleasm/singularity/mb_downloads"
echo "$file"
out=$(basename "$file" .fa)

# singularity build compleasm_v${VERSION}.sif docker://huangnengcsu/compleasm:v${VERSION}
# singularity exec docker://huangnengcsu/compleasm:v${VERSION} compleasm download cetartiodactyla
VERSION=0.2.6
mkdir -p "$dir"/compleasm_"$out"
singularity exec compleasm_v0.2.6.sif compleasm run -t 8 -l cetartiodactyla -L "$odb" -a "$file" -o "$dir"/compleasm_"$out"
