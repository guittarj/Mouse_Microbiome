wd <- 'C:\\Users\\John\\Google Drive\\Mouse_Microbiome_Shared_Files\\Data\\'
setwd(wd)

#load otu data
x <- fread('mouse_combined_merged_OTU_table.txt') %>%
  rename(otu = `#OTU ID`) %>%
  gather(sample, abun, -otu)

#load metadata 
x1 <- fread('Reese_2018_04_29_20151231_meta.txt') %>% as.data.frame()
names(x1) <- x1[3,]
x1 <- x1[c(6:nrow(x1)), ]

x2 <- fread('Reese_2018_04_29_20150930_meta.txt') %>% as.data.frame()
names(x2) <- x2[3,]
x2 <- filter(x2, sample_alias != '')
x2 <- x2[c(6:nrow(x2)), ]

#merge metadata

cutsub <- function(strings, pattern, cut) {
  # function to extract a substring in a vector of strings...
  tmp <- str_extract(strings, pattern)
  # then trim that extraction...
  tmp <- gsub(cut, '', tmp)
  #then convert to numeric...
  tmp <- as.numeric(tmp)
  return(tmp)
}

meta <- rbind(x1, x2) %>%
  rename(x = sample_description) %>%
  transmute(sample = sample_alias, 
            mouse = cutsub(x, "mouse \\d+", "mouse "),
            cohoused = grepl('cohoused=y', x),
            partner = cutsub(x, "cohoused=y with mouse \\d+", 'cohoused=y with mouse '),
            pair = as.numeric(factor(pmax(mouse, partner))),
            hour = ifelse(grepl('day', x), cutsub(x, "day \\d+", 'day ') * 24,
                     ifelse(grepl('hour', x), cutsub(x, "hour \\d+", 'hour '), 0)),
            cells_per_gram = as.numeric(`organism count`) / as.numeric(`amount or size of sample collected`), 
            Abx = `miscellaneous parameter` == 'antibiotic')


#remove samples that are in metadata but did not sequence correctly
meta <- filter(meta, !sample %in% as.character(c(1.02, 1.097, 10.087, 11.087, 13.087, 14.087, 15.073, 16.073, 3.02, 5.087, 8.008)))

#fix another mapping error
meta$hour[meta$mouse < 20 & meta$hour == 144] <- 120
meta$hour[meta$mouse < 20 & meta$hour == 168] <- 144
meta$hour[meta$mouse < 20 & meta$hour == 192] <- 168
meta$hour[meta$mouse < 20 & meta$hour == 216] <- 192
meta$hour[meta$mouse < 20 & meta$hour == 23 ] <- 216

#load list of sample headers
heds <- fread('C:\\Users\\John\\Desktop\\headers.csv')
heds$yo2 <- grepl('@', heds$yo)
heds <- data.frame(filename = heds$yo[!heds$yo2], header = heds$yo[heds$yo2])
heds$header <- as.character(heds$header)
heds$sample <- unlist(lapply(strsplit(heds$header, ' '), function(x) x[2]))
heds$sample <- unlist(lapply(strsplit(heds$sample, '_'), function(x) x[1]))
heds$sample <- ifelse(grepl('9454', heds$sample), substr(heds$sample, 1, nchar(heds$sample) - 2), heds$sample)
heds$experiment <- heds$sample %in% meta$sample
heds$filename_short <- substr(heds$filename, 1, 10)

#add filename to meta table
meta <- meta %>% mutate(filename = heds$filename_short[match(sample, heds$sample)])

#fix sample names in otu table
x <- x %>% 
  mutate(sample = meta$sample[match(x$sample, meta$filename)]) %>%
  select(sample, otu, abun)

#metadata entries with missing samples
meta$sample[!meta$sample %in% heds$sample] %>% sort()

#samples missing from metadata
heds$sample[!heds$experiment] %>% sort()

#remove filename from meta table and pare down to only rows with extant samples
meta$filename <- NULL
meta <- meta %>% mutate(sample_exists = sample %in% x$sample)

#filenames I need to move out of the folder
#paste(as.character(heds$filename)[!heds$experiment], collapse = ' ')

###Taxonomy


saveRDS(x, 'otu_table_usearch_aug19.RDS')
saveRDS(meta, 'metadata_aug19.RDS')
