## GENOME ASSEMBLY

1-2. To assemble the genomes for run1 and run2 with 58x HiFi and 107x ONT (filtered to >40kb)

      verkko -d Verkko2C2 \
        --hifi TulixWagyu_HiFi/tencells.fastq.gz \
        --nano TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_107x.fastq.gz \
        --hap-kmers meryl/maternal_compress.k21.hapmer.meryl meryl/paternal_compress.k21.hapmer.meryl trio \
        --snakeopts "--cores 200"

      verkko -d Verkko2C0 \
        --hifi TulixWagyu_HiFi/tencells.fastq.gz \
        --nano TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_107x.fastq.gz \
        --hap-kmers meryl/maternal_compress.k21.hapmer.meryl meryl/paternal_compress.k21.hapmer.meryl trio \
        --snakeopts "--cores 200"

3-4. To assemble the genomes for run3 and run4 with 58x HiFi and 121x ONT (filtered to >40kb)

      verkko -d Verkko2C2 \
        --hifi TulixWagyu_HiFi/tencells.fastq.gz \
        --nano TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_107x.fastq.gz TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_14x.fastq.gz \
        --hap-kmers meryl/maternal_compress.k21.hapmer.meryl meryl/paternal_compress.k21.hapmer.meryl trio \
        --snakeopts "--cores 200"

      verkko -d Verkko2C0 \
        --hifi TulixWagyu_HiFi/tencells.fastq.gz \
        --nano TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_107x.fastq.gz TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_14x.fastq.gz \
        --hap-kmers meryl/maternal_compress.k21.hapmer.meryl meryl/paternal_compress.k21.hapmer.meryl trio \
        --snakeopts "--cores 200"

5-6. To assemble the genomes for run5 and run6 with 58x HiFi, 57x HERRO corrected reads and 121x ONT (filtered to >40kb)

      verkko -d Verkko2C2 \
        --hifi TulixWagyu_HiFi/tencells.fastq.gz TulixWagyu_HERRO/herro_57x_q10_l10k.fastq.gz \
        --nano TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_107x.fastq.gz TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_14x.fastq.gz \
        --hap-kmers meryl/maternal_compress.k21.hapmer.meryl meryl/paternal_compress.k21.hapmer.meryl trio \
        --snakeopts "--cores 200"

      verkko -d Verkko2C0 \
        --hifi TulixWagyu_HiFi/tencells.fastq.gz TulixWagyu_HERRO/herro_57x_q10_l10k.fastq.gz \
        --nano TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_107x.fastq.gz TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_14x.fastq.gz \
        --hap-kmers meryl/maternal_compress.k21.hapmer.meryl meryl/paternal_compress.k21.hapmer.meryl trio \
        --snakeopts "--cores 200"

7-8. To assemble the genomes for run7 and run8 with 58x HiFi, 99x HERRO corrected reads and 121x ONT (filtered to >40kb)

      verkko -d Verkko2C2 \
        --hifi TulixWagyu_HiFi/tencells.fastq.gz TulixWagyu_HERRO/herro_99x_q10_l10k.fastq.gz \
        --nano TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_107x.fastq.gz TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_14x.fastq.gz \
        --hap-kmers meryl/maternal_compress.k21.hapmer.meryl meryl/paternal_compress.k21.hapmer.meryl trio \
        --snakeopts "--cores 200"

      verkko -d Verkko2C0 \
        --hifi TulixWagyu_HiFi/tencells.fastq.gz TulixWagyu_HERRO/herro_99x_q10_l10k.fastq.gz \
        --nano TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_107x.fastq.gz TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_14x.fastq.gz \
        --hap-kmers meryl/maternal_compress.k21.hapmer.meryl meryl/paternal_compress.k21.hapmer.meryl trio \
        --snakeopts "--cores 200"

9. To assemble the genomes using hifiasm with 58x HiFi, 57x HERRO corrected reads and 121x ONT (filtered to >40kb)

            hifiasm -o hifiasm_run1 -t 32 \
                  -1 /hpcfs/users/a1223107/Tuli_x_Wagyu_data/hifiasm_based_assemblies/yak_file/pat.yak \
                  -2 /hpcfs/users/a1223107/Tuli_x_Wagyu_data/hifiasm_based_assemblies/yak_file/mat.yak \
                  --ul TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_107x.fastq.gz,TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_14x.fastq.gz \
                  TulixWagyu_HERRO/herro_57x_q10_l10k.fastq.gz tencells.fastq.gz

9. To assemble the genomes using hifiasm with 58x HiFi and 121x ONT (filtered to >40kb)

            hifiasm -o hifiasm_run1 -t 32 \
                  -1 /hpcfs/users/a1223107/Tuli_x_Wagyu_data/hifiasm_based_assemblies/yak_file/pat.yak \
                  -2 /hpcfs/users/a1223107/Tuli_x_Wagyu_data/hifiasm_based_assemblies/yak_file/mat.yak \
                  --ul TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_107x.fastq.gz,TulixWagyu_ONT_Dorado/ONT_TxW_Filt40k_14x.fastq.gz \
                  tencells.fastq.gz

## POLISHING with DEEPVARIANT

1. First is to align the PacBio HiFi reads with the draft assembly

            ref="Wagyu_haplotype1_v1.unpolished.fa"
            hifi="F1_HiFi_longreads/tencells.fastq.gz"
            baseref=$(basename "$ref" .fa*)
            
            minimap2 -t 46 -a -asm20 "$ref" "$hifi" | samtools sort -o mmp_"$baseref"_aln_hifi.sorted.bam -@46
            samtools index mmp_"$baseref"_aln_hifi.sorted.bam -@46
            samtools faidx "$ref"

2. Next is to run deepvariant

            export TMPDIR=.
            ulimit -u 10000
            
            singularity run --bind /usr/lib/locale/ \
              deepvariant_1.6.1.sif \
                /opt/deepvariant/bin/run_deepvariant \
                --model_type PACBIO \
                --ref Wagyu_haplotype1_v1.unpolished.fa \
                --reads mmp_Wagyu_haplotype1_v1.unpolished.fa_aln_hifi.sorted.bam \
                --output_vcf Wagyu_haplotype1_v1.dv.vcf.gz \
                --num_shards $(nproc) \
                --intermediate_results_dir intermediate_results_dir

3. Then make consensus

            wagyu="Wagyu_haplotype1_v1.unpolished.fa"
            vcfwagyu="/Wagyu_haplotype1_v1.dv.vcf.gz"
            
            wagyu_base=$(basename "$wagyu" .unpolished.fa)
            vcfwagyu_base=$(basename "$vcfwagyu" .vcf.gz)
            
            bcftools view -f PASS "$vcfwagyu" | bcftools view -i 'GT="1/1"' -o "$vcfwagyu_base".PASS.homoalt.vcf.gz
            bcftools index "$vcfwagyu_base".PASS.homoalt.vcf.gz
            bcftools consensus -f "$wagyu" "$vcfwagyu_base".PASS.homoalt.vcf.gz > "$wagyu_base".polished.fa
