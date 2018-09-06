#Mouse project script: going from raw seqs to an OTU table

#make new directory and go there
mkdir mouse
cd mouse/

#load most recent version of R
module swap GNU GNU/4.9
module load R/3.3.2
R

# download list of accessions
download.file('https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB26446&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,sra_ftp,sra_galaxy,cram_index_ftp,cram_index_galaxy&download=txt', 'Reese2018_accessions.txt')

#read and simplify into list of URLS.
urls <- read.table('Reese2018_accessions.txt', header = TRUE)
urls <- unlist(strsplit(as.character(urls$submitted_ftp), ';'))
urls <- paste0('http://', urls) 

#save table, unquoted
write.table(urls, 'Reese2018_accessions.txt', col.names = F, row.names = F, quote = F)

#quit R
quit(save = 'no')

# download the files from ENA, uncompress
wget -i Reese2018_accessions.txt
gunzip ERR*.fastq.gz

#save list of headers with sample names
head -n 1 ERR*.fastq > headers.txt

#R for data munging
R

#load list of sample headers (necessary to identify samples)
heds <- readChar('headers.txt', file.info('headers.txt')$size)
heds <- unlist(strsplit(heds, '==> '))
heds <- heds[c(2:length(heds))]
sample_names <- gsub(' .*', '', heds)
ids <- gsub('.*\\.1 ', '', heds)
ids <- gsub('.*\\.1 ', '', ids)
ids <- gsub('_.*', '', ids)
unwanted_fastqs <- paste(sample_names[grepl('D|swab', ids)], collapse = ' ')
sample_dict <- data.frame(sample_name = sample_names, sample_id = ids)

#write
write(unwanted_fastqs, 'unwanted.fastqs.txt')
write.csv(sample_dict, 'sample_dict.csv')

#back to bash
quit(save = 'no')

#remove unwanted fastqs
xargs rm < unwanted.fastqs.txt

# rename the files to have an 'R1' and 'R2' in them
rename '_' '_R' *.fastq

#Start Usearch analysis
#merge paired ends
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout merged.fastq -tabbedout merged_report.txt -alnout merged_aln.txt

#filter low quality or short reads
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter merged.fastq -fastq_maxee 1 -fastq_trunclen 250 -fastaout merged_filtered.fa

#get unique samples
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques merged_filtered.fa -fastaout merged_filtered_uniques.fa -sizeout

#cluster into otus
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus merged_filtered_uniques.fa -otus otus.fa -uparseout otus_uparse.txt -relabel OTU

#identify zotus
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -unoise3 merged_filtered_uniques.fa -zotus zotus.fa -tabbedout zotus_report.txt

#fix lower case zotu problem
sed -i 's/Zotu/ZOTU/g' zotus.fa

#create otu table
#note that the query file is pre-filtering. Sometimes even low quality reads can be successfully mapped.
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -otutab merged.fastq -otus otus.fa -uc otu_map.uc -otutabout otu_table.txt -biomout otu_jsn.biom -notmatchedfq otu_unmapped.fq

#create zotu table
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -otutab merged.fastq -zotus zotus.fa -uc zotu_map.uc -otutabout zotu_table.txt -biomout zotu_jsn.biom -notmatchedfq zotu_unmapped.fq

#download databases for classifications
#I'm using SILVA on the mouse data, even though it is "not recommended" by USEARCH website, because 1) the LTP routh has too many 'unclassified' genera, and 2) even if the taxonomy has some errors, I'm just looking for broad trends, and just need a list of genera to use when searching for traits.
wget "http://drive5.com/sintax/silva_16s_v123.fa.gz"
gunzip silva_16s_v123.fa.gz

#classifying otus
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sintax otus.fa -db ltp_16s_v123.fa -tabbedout otus_taxonomy.sintax -strand both

#classifying zotus
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sintax zotus.fa -db ltp_16s_v123.fa -tabbedout merged_zotus_taxonomy.sintax -strand both

