The raw sequencing reads generated were:

| Sequencing   data                        |  Coverage       | N50 read length (bp) |  Mean read quality              | Number of reads      |
|------------------------------------------|-----------------|----------------------|---------------------------------|----------------------|
| PacBio HiFi reads                        |           58.1  | 20277                |                           22.9  |                      |
| Proximo Hi-C short reads   (paired-end)  |           42.2  | 150                  |                           36.0  |         421,348,414  |
| ONT reads R9.4.1                         |           72.8  | 53267                |                           13.3  |                      |
| ONT reads R10.4.1                        |         156.0   | 48867                |                           17.4  |                      |
| ONT ultra long > 100Kb                   |       18.3      | 111671               |                           13.8  |                      |
| Illumina paired-end short   reads (F1)   |           81.8  | 150                  |                           34.0  |         815,812,681  |
| Illumina paired-end short   reads (Sire) |           40.0  | 150                  |                           36.0  |         400,803,174  |
| Illumina paired-end short   reads (Dam)  |           38.0  | 150                  |                           36.0  |         381,262,284  |

To assemble the genomes for run1 and run2

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

To assemble the genomes for run1 and run2

- Sequencing of the trio (Tuli x Wagyu) in the main text    
- Sequencing of the trio (Tuli x Wagyu) in the Supplementary Information (Supplementary Methods)
