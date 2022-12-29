library(tidyverse)
library(rjson)
library(Biostrings)

indir <- '../data'

Vdir <- dir(indir)

# pangolin
Vf <- str_glue('{indir}/{Vdir}/{Vdir}/{Vdir}_pangolin_lineage.csv')
TB <- tibble()
for (f in Vf){
    if (file.exists(f)){
        TB <- bind_rows(TB, read_csv(f))
    }
}
write_excel_csv(TB, 'merged_pangolin.csv')

# consensus fa
Vf <- str_glue('{indir}/{Vdir}/{Vdir}/{Vdir}_consensus.fa')

Vdna <- DNAStringSet()
for (f in Vf){
    if (file.exists(f)){
        dna <- readDNAStringSet(f, format='fasta')
        Vdna <- c(Vdna, dna)
    }
}

writeXStringSet(Vdna, 'merged_consensus.fasta', format='fasta')


# log
Vf <- str_glue('{indir}/{Vdir}/HAVoC_run.log')

Llog <- list()
for (f in Vf){
    if (file.exists(f)){

        # name
        name <- f %>% str_split('/') %>% map_chr(~.x[2])
        # parse log
        Vlog <- file(f,open="r") %>% readLines
        # reads
        
        Vtotal_reads <- Vlog %>% grep('total reads', ., value=TRUE) %>% 
            str_split(': ') %>% map_chr(~.x[2])
        names(Vtotal_reads) <- c('R1_reads_raw', 'R2_reads_raw', 'R1_reads_filtered', 'R2_reads_filtered')
        # q30
        Vq30 <- Vlog %>% grep('Q30 bases', ., value=TRUE) %>% str_split(': ') %>% map_chr(~.x[2])
        names(Vq30) <- c('R1_q30_raw', 'R2_q30_raw', 'R1_q30_filtered', 'R2_q30_filtered')
        # dup
        dup_rate <- c('dup_rate' = Vlog %>% grep('Duplication rate', ., value=TRUE) %>% 
            str_split(': ') %>% map_chr(~.x[2]))
        # align rate
        align_reads <- c('align_reads' = Vlog %>% grep('Found', ., value=TRUE) %>% 
                grep(' mapped reads', ., value=TRUE) %>% 
                str_split(' ') %>% map_chr(~.x[2]))

        Llog[[name]] <- c(Vtotal_reads, Vq30, dup_rate, align_reads)
    }
}

F <- file('merge_log.json', "w" )
Llog %>% toJSON %>% writeLines(F)
close(F)

