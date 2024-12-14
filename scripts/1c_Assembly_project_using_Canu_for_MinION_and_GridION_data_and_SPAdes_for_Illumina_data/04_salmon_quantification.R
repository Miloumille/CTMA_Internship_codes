library(purrr)
library(ShortRead)

fastqlist <- list.files('data_processed//fastq/',full.names = T)
forward <- fastqlist[grepl('_1',fastqlist)]
reverse <- fastqlist[grepl('_2',fastqlist)]

samplename <- map_chr(.x = basename(forward),.f = function(x) strsplit(x,split='_1')[[1]][1])

data.frame(samplename,forward,reverse)

################ quantification against mouse

c('/media/ctma/Disk2/01_fichier/preprocessing/00_database/rnaseq/salmon_index/GRCm39_index/')

for(i in 1:12)
{
  myarg <- paste0('run --rm -v ./:/data --rm -v /media/ctma/Disk2/01_fichier/preprocessing/00_database/rnaseq/salmon_index:/indexdir --cpus=12 combinelab/salmon:1.10.3 sh -c "salmon quant -i /indexdir/GRCm39_index -l A -1 /data/',forward[i],' -2 /data/',reverse[i],' --validateMappings -p 8 -o /data/data_processed/salmon_quant/',samplename[i],'"')   
  system2(command='docker',args=myarg)  
  print(i)
}




