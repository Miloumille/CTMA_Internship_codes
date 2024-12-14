samplelist <- read.csv('data/sample_list.csv')

runs <- samplelist$Run
N_runs <- length(runs)

for(i in 1:N_runs)
{
  myarg <- paste0('run --rm -v ./:/data1 --cpus=6 staphb/sratoolkit:3.0.7 sh -c " fastq-dump -X 100 --split-files --gzip -O /data1/data_processed/fastq ',runs[i],'"')  
  print(myarg)
  system2(command='docker',args=myarg)
}
