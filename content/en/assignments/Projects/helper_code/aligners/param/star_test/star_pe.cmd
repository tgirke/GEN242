STAR --runThreadN 4 \ 
     --genomeDir ./data \
     --readFilesIn $fq1 $fq2 \ 
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --readFilesCommand zcat --outFileNamePrefix $base"_"
