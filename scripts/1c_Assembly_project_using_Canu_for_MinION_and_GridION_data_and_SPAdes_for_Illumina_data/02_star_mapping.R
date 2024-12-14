
library(purrr)
library(ShortRead)


#unlink('data_processed/STAR_mapping/',recursive = T,force = T)
dir.create('data_processed/STAR_mapping/')


fastqlist <- list.files('data_processed/fastq/',full.names = T)
forward <- fastqlist[grepl('_1',fastqlist)]
reverse <- fastqlist[grepl('_2',fastqlist)]

samplename <- map_chr(.x = basename(forward),.f = function(x) strsplit(x,split='_1')[[1]][1])

data.frame(samplename,forward,reverse)


Nfiles <- length(samplename)

for(i in 1:Nfiles)
{
  myarg <- paste0('run --rm -v ./:/data --cpus=6 zavolab/star:2.7.1a sh -c "STAR --genomeDir /indexdir/mm_grcm39_index/  --readFilesIn /data/',forward[i],' /data/',reverse[i],' --readFilesCommand gunzip -c --outFileNamePrefix /data/data_processed/STAR_mapping/',samplename[i],' --outSAMtype BAM   SortedByCoordinate --outSAMattributes NM --runThreadN 8','"')   
  system2(command='docker',args=myarg)  
  print(i)
}


