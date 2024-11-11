STAR --runThreadN 40 --runMode genomeGenerate \
    --genomeDir human_reference_genome/index_file_STAR \
    --genomeFastaFiles human_reference_genome/files_from_genecode/GRCh38.p14.genome.fa \
    --sjdbGTFfile human_reference_genome/files_from_genecode/gencode.v45.chr_patch_hapl_scaff.annotation.gtf \
    human_reference_genome/files_from_genecode/gencode.v45.long_noncoding_RNAs.gtf \
    human_reference_genome/files_from_genecode/gencode.v45.tRNAs.gtf \
    human_reference_genome/files_from_genecode/gencode.v45.polyAs.gtf