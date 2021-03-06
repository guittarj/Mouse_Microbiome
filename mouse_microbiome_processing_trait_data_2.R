#Mouse gut microbiome traits -- data gathering and cleaning
#John Guittar

##Overview
#1. load raw taxonomic data from usearch output
#2. load taxonomy from LTP (Living Tree Project)
#3. remove any ambiguities/special characters
#4. gather and combine trait data from various sources
#5. filter trait data to only relevant taxa (either from usearch or LTP)
#6. load tree of both usearch and LTP taxa, created with usearch -calc_distmx and -cluster_aggd
#6.5 Do not infer unkonwn trait values yet! that's for another script
#7. save clean taxonomic, trait, and phylogenetic data

######################################################################
######################################################################
######################################################################

#load taxonomic data from usearch output
#remove any Genus/Species names with numbers or special characters
#set working directory, load packages
wd <- 'C:\\Users\\John\\Google Drive\\Mouse_Microbiome_Shared_Files\\Data\\'
setwd(wd)
source('C:\\Users\\John\\Documents\\msu\\Mouse_Microbiome\\custom_functions.R')
loadpax(c('data.table','tidyverse','castor','phytools'))

#load raw taxonomy data from usearch
#meta <- readRDS('metadata_aug19.RDS')
tax <- read.table('combined_merged_otus_taxonomy.sintax', fill = TRUE) %>%
  select(otu = V1, tax = V4) %>%
  separate(tax, c('Kingdom','Phylum','Class','Order','Family','Genus','Species'), sep = ',', fill = 'right') %>%
  mutate_all(funs(gsub('\\w:', '', .))) %>%
  mutate_all(funs(ifelse(is.na(.), 'unclassified', .)))

# Clean up taxonomy names -- no ambiguous genera names (important when querying bacdive)
library(stringr)
tax$Genus <- unlist(lapply(strsplit(tax$Genus, '_'), function(x) if(x[1] %in% tax$Family) x[2] else x[1]), tax$Genus)
tax$Genus[grepl('\\d', tax$Genus)] <- 'unclassified'
tax$Genus <- gsub('\\[|\\]', '', tax$Genus)
tax$Species <- unlist(lapply(strsplit(tax$Species, '_'), function(x) x[length(x)]))
tax$Species[grepl('\\d|\\.', tax$Species)] <- 'unclassified'
tax[is.na(tax)] <- 'unclassified'

saveRDS(tax, 'Reese_Tax_Aug2018.RDS')
tax <- readRDS('Reese_Tax_Aug2018.RDS')

#load our OTU data; 
#determine our taxonomy data (to know which species to focus on while webcrawling, etc)
otus <- readRDS('otu_table_usearch_aug19.RDS')
otus_wide <- spread(otus, otu, abun)

### Clean
####First: Edits drawn from from Barberan et al 2016 script
####Then: my own edits

#read table
ijsem<-read.delim(paste0(wd, 'IJSEM_pheno_db_v1.0.txt'), sep="\t", header=T, check.names=F, fill=T,
                  na.strings=c("NA", "", "Not indicated", " Not indicated","not indicated", "Not Indicated", "n/a", "N/A", "Na", "Not given", "not given","Not given for yeasts", "not indicated, available in the online version", "Not indicated for yeasts", "Not Stated", "Not described for yeasts", "Not determined", "Not determined for yeasts"))

#simplify column names
colnames(ijsem)<-c("Habitat", "Year", "DOI", "rRNA16S", "GC", "Oxygen",
                   "Length", "Width", "Motility", "Spore", "MetabAssays", "Genus", "Species", "Strain", "pH_optimum", "pH_range", "Temp_optimum", "Temp_range", "Salt_optimum", "Salt_range", "Pigment", "Shape", "Aggregation", "FirstPage", "CultureCollection", "CarbonSubstrate", "Genome", "Gram", "Subhabitat", "Biolog")

#clean Habitat column
levels(ijsem$Habitat)[levels(ijsem$Habitat)=="freshwater (river, lake, pond)"]<-"freshwater"
levels(ijsem$Habitat)[levels(ijsem$Habitat)=="freshwater sediment (river, lake, pond"]<-"freshwater sediment"

#clean Oxygen column
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="aerobic"]<-"obligate aerobe"
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="anerobic"]<-"obligate anerobe"
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="microerophile"]<-"microaerophile"

#clean pH_optimum column
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
ijsem$pH_optimum<-as.character(ijsem$pH_optimum)
ijsem$pH_optimum<-sapply(ijsem$pH_optimum, simplify=T, function(x){mean(swan(unlist(strsplit(x, split="-", fixed=T))), na.rm = T)})

#remove pH values <0 and >10
ijsem$pH_optimum[ijsem$pH_optimum<0 | ijsem$pH_optimum>10]<-NA

#clean Temp_optimum column
ijsem$Temp_optimum<-as.character(ijsem$Temp_optimum)

#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
ijsem$Temp_optimum<-sapply(ijsem$Temp_optimum, simplify=T, function(x){mean(swan(unlist(strsplit(x, split="-", fixed=T))))})

#clean Salt_optimum column
ijsem$Salt_optimum<-as.character(ijsem$Salt_optimum)
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs

ijsem$Salt_optimum<-sapply(ijsem$Salt_optimum, simplify=T, function(x){mean(swan(unlist(strsplit(x, split="-", fixed=T))))})
#there are some formatting issues that should be solved


############## Now, my additional edits

#Assign oxygen score
ijsem$Oxygen_score <- c(5,4,3,2,1)[match(ijsem$Oxygen, 
                                         c('obligate aerobe','microaerophile','facultative anaerobe','facultative aerobe','anaerobic'))]

#trun Aggregation into a binary
ijsem$Aggregation <- c(0,1,1)[match(ijsem$Aggregation, c('none','chain','clump'))]

#clean Length and Width
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
x <- as.character(ijsem$Length)
x <- gsub(' in diameter| \\(diameter\\)|mm|Indicated|diameter: | ', '', x)
x <- gsub(',', '.', x, fixed = TRUE)
x[grepl(">", x)] <- NA
x[x == ''] <- NA
x[x == c('0.61.6')] <- '0.6-1.6'
x[x == c('0.4 0.9 ')] <- '0.4-0.9'
x <- gsub('��� *', '-', x)

x <- ifelse(grepl("-[a-zA-Z]+", x) & !is.na(x), 
            as.character(format(as.Date(x, "%d-%b"), "%d.%m")), 
            x)
x <- ifelse(grepl("[a-zA-Z]+-", x) & !is.na(x), 
            as.character(format(as.Date(paste(1, x), "%d %b-%y"), "%m.%y")), 
            x)

#if it has a dash, assume it was a range and take the mean
filt <- grepl("-", x) & !is.na(x)
x[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[2]]))))

#if it has a backslash, assume it was a decimal point
filt <- grepl("\\/", x) & !is.na(x)
x[filt] <- paste(sapply(strsplit(x[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(x[filt], '\\/'), function(x) x[[2]]), sep = '.')

ijsem$Length <- swan(x)


######################

x <- as.character(ijsem$Width)
x <- gsub(' |\\(|\\)|not.+|indiameter|Filament.+', '', x)
x <- gsub(',', '.', x, fixed = TRUE)
x[grepl("<", x)] <- NA
x[x == ''] <- NA
x <- gsub('��� *', '-', x)
x <- ifelse(grepl("-[a-zA-Z]+", x) & !is.na(x), 
            as.character(format(as.Date(x, "%d-%b"), "%d.%m")), 
            x)
x <- ifelse(grepl("[a-zA-Z]+-", x) & !is.na(x), 
            as.character(format(as.Date(paste(1, x), "%d %b-%y"), "%m.%y")), 
            x)

#if it has a dash, assume it was a range and take the mean
filt <- grepl("-", x) & !is.na(x)
x[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[2]]))))

#if it has a backsclash, assume it was a decimal point
filt <- grepl("\\/", x) & !is.na(x)
x[filt] <- paste(sapply(strsplit(x[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(x[filt], '\\/'), function(x) x[[2]]), sep = '.')

ijsem$Width <- as.numeric(x)

#############################################

#Turn motility into a binary
ijsem$Motility <- c(0,1,1,1,1)[match(ijsem$Motility, c('non-motile', 'flagella',
                                                       'motile, but unspecified structure', 'gliding', 'axial filament'))]

#turn spore into binary
ijsem$Spore <- c(0,1)[match(ijsem$Spore, c('no','yes'))]

#turn gram status into a binary
ijsem$Gram <- c(0, 0.5, 1)[match(ijsem$Gram, c('negative','variable','positive'))]

# fix pH
# NOTE: THE WAY THAT THEY CALCULATED PH WHEN THERE IS A RANGE - BY TAKING THE MEAN -- IS INCORRECT
# MUST UNLOG PH BEFORE CALCULATING MEAN, AND THEN RE-LOG.
# I HAVEN'T DONT THAT YET HERE... 

x <- as.character(ijsem$pH_range)
x <- gsub(' |\\(|\\)|at 6\\.0 but not 5\\.5|>|<|not inidcated', '', x)
x <- gsub('MinimumpHof5|t|4\\.0\\+|5\\.0\\+;notat4\\.5|14\\.08', '', x)
x[x %in% c('13.05','4.3','7.1','7.2','5.5','7.5','30.04','5.1')] <- NA
x[x == '5.7and6.8'] <- '5.7-6.8'
x[x == '4.0and8.5'] <- '4-8.5'
x[x == '5,6,10,12,'] <- '5-12'
x[x == '5and9.5'] <- '5-9.5'
x[x == '5.59.5'] <- '5.5-9.5'

# remove anything with an inequality
x[grepl("<|>", x)] <- NA
x[x == ''] <- NA
x <- gsub('��� *', '-', x)
x <- ifelse(grepl("-[a-zA-Z]+", x) & !is.na(x), 
            as.character(format(as.Date(x, "%d-%b"), "%d.%m")), 
            x)

filt <- grepl("\\/", x) & !is.na(x)
x[filt] <- paste(sapply(strsplit(x[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(x[filt], '\\/'), function(x) x[[2]]), sep = '.')

x[!grepl("-", x) & !is.na(x) & !grepl("\\.", x)] <- NA
filt <- !grepl("-", x) & !is.na(x) & grepl("\\.", x)
x[filt] <- gsub("\\.", "-", x[filt])
lows <- swan(sapply(strsplit(x, '-'), function(x) x[[1]]))
highs <- swan(sapply(strsplit(x, '-'), function(x) if (length(x) > 1) x[[2]] else NA))

ijsem$pH_lows <- lows
ijsem$pH_highs <- highs

####GC
x <- as.character(ijsem$GC)
x <- gsub('��� *', '-', x)
x <- gsub(',', '\\.', x)
x[x %in% c('DQ987877','marine organism','BAOL01000001','40.50%','GU323338')] <- NA
x <- gsub('O�.+|\\+/-.+|%|o', '', x)
filt <- grepl("-", x) & !is.na(x)
x[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(x[filt], '-'), function(x) x[[2]]))))

ijsem$GC <- swan(x)

# at the moment, i only want these traits...
ijsem <- select(ijsem, Genus, Species, GC_content = GC, 
                Oxygen_tolerance = Oxygen_score, Length, Width, 
                Motility = Motility, Spore, 
                pH_optimum, Temp_optimum, Salt_optimum, Aggregation_score = Aggregation, 
                Gram_positive = Gram) %>%
  gather(trait, val, -Genus, -Species) %>%
  filter(!is.na(val)) %>%
  filter(!Genus %in% c("10.1099/ijs.0.006320-0","10.1099/ijs.0.02424-0")) %>%
  mutate_if(is.factor, as.character)

###########################################################
# now, sporulation data from Browne et al 2016

spo <- read.csv('Browne2016_sporulationTable.csv', header = T, skip = 1)

# fix long colnames
colnames(spo) <- substr(colnames(spo), 1, 2)

spo <- spo %>%
  transmute(
    Genus = unlist(lapply(strsplit(as.character(Sp), " |_"), function(x) x[[1]])),
    Species = unlist(lapply(strsplit(as.character(Sp), " |_"), function(x) x[[2]])),
    trait = 'Spore_score',
    val = si)

####################################### Exploring BacDrive

if (FALSE) {
  #IFF i want just the genera in the mouse experiment:
  library(RCurl)
  library(rjson)
  
  #needed <- sort(as.character(unique(filter(tax, is.na(Gram_positive))$Genus)))
  needed <- sort(unique(tax$Genus_clean))
  needed <- needed[!needed %in% c("unclassified", needed[grep('\\d', needed)])]
  
  bacdat1 <- bac_search(needed[c(199:224)]) #re-calculate bacdata
  bacdat2 <- bac_search(needed[c(100:199)]) #re-calculate bacdata
  bacdat3 <- bac_search(needed[c(1:99)]) #re-calculate bacdata
  bacdat <- c(bacdat1, bacdat2, bacdat3)
  saveRDS(bacdat, file = paste0(wd, 'bacdive_genera_', Sys.Date(), '.RDS'))
  
} 

if (FALSE) {
  
  #OR, if I want all bac IDS: (Note this kept failing halfway through so I no longer recommend it)
  library(RCurl)
  library(rjson)
  
  needed <- read.csv(paste0(wd, 'export_bacdive_ids.csv'), header = FALSE)$V1
  needed <- needed[c(9515:length(needed))]
  
  #tmp
  bacdat <- readRDS('C:\\Users\\John\\Google Drive\\Mouse_Microbiome_Shared_Files\\Data\\bacdat9515.RDS')
  
  bacdat <- bac_search_ID(needed, bacdat = bacdat) #re-calculate bacdata
  saveRDS(bacdat, file = paste0(wd, 'bacdive_bacIDs_', Sys.Date(), '.RDS'))
  
} 

bacdat <- readRDS('bacdive_genera_2018-08-22.RDS') #load bacdata
          
#taxonomic info
tax_bacdat <- sapply(bacdat, function(x)
  sapply(x, function(y)
    sapply(y$taxonomy_name$strains_tax_PNU, function(sp) sp$species_epithet)))

na.pad <- function(x) ifelse(is.null(x), NA, x)
  
names(bacdat) <- paste0('bacID', 1:length(bacdat))
tax_bacdat <- lapply(bacdat, function(z)
  data.frame(
    Phylum = na.pad(z[[1]]$taxonomy_name$strains_tax_PNU[[1]]$phylum), 
    Class = na.pad(z[[1]]$taxonomy_name$strains_tax_PNU[[1]]$class), 
    Family = na.pad(z[[1]]$taxonomy_name$strains_tax_PNU[[1]]$family),
    Genus = na.pad(z[[1]]$taxonomy_name$strains_tax_PNU[[1]]$genus),
    Species = na.pad(z[[1]]$taxonomy_name$strains_tax_PNU[[1]]$species_epithet),
    stringsAsFactors = FALSE))

tax_bacdat <- bind_rows(tax_bacdat)
tax_bacdat$id <- paste0('bacID', c(1:nrow(tax_bacdat)))

#oxygen tolerance
o2 <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$oxygen_tolerance, function(val) val$oxygen_tol)))
o2 <- bac_extract(o2, 'Oxygen_tolerance')
o2$val <- c(5,5,4,3,2,1,1)[match(o2$val, 
                             c('obligate aerobe','aerobe','facultative aerobe','microaerophile',
                               'facultative anaerobe','anaerobe','obligate anaerobe'))]

#gram stain
gs <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$cell_morphology, function(val) val$gram_stain)))
gs <- bac_extract(gs, 'Gram_positive')
gs$val <- c(1, 0.5, 0)[match(gs$val, c('positive','variable','negative'))]

#spore forming
spofo <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$spore_formation, function(val) val$ability)))
spofo <- bac_extract(spofo, "Spore")
spofo$val <- ifelse(spofo$val == 'TRUE', 1, 0)

#motility
motility <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$cell_morphology, function(val) val$motility)))
motility <- bac_extract(motility, 'Motility')
motility$val <- ifelse(motility$val == TRUE, 1, 0)

#cell length
cell_length <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$cell_morphology, function(val) val$cell_len)))
cell_length <- lapply(cell_length, unlist)
cell_length <- lapply(cell_length, function(x) ifelse(grepl('<|>', x), NA, x))
cell_length <- lapply(cell_length, function(x) ifelse(grepl('-', x), strsplit(x, '-'), x))
cell_length <- bac_extract(cell_length, 'Length')
cell_length$val <- as.numeric(cell_length$val)

#cell width
cell_width <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$cell_morphology, function(val) val$cell_width)))
cell_width <- lapply(cell_width, unlist)
cell_width <- lapply(cell_width, function(x) ifelse(grepl('<|>', x), NA, x))
cell_width <- lapply(cell_width, function(x) ifelse(grepl('-', x), strsplit(x, '-'), x))
cell_width <- bac_extract(cell_width, 'Length')
cell_width$val <- as.numeric(cell_width$val)

#aggregation (\any data at all? Am I doing it right?)
aggreg <- sapply(bacdat, function(genus) 
  sapply(genus, function(entry) 
    sapply(entry$morphology_physiology$multicellular_morphology, function(val) val$ability)))
aggreg <- bac_extract(aggreg, 'Aggregation_score')
aggreg$val <- ifelse(aggreg$val == TRUE, 1, 0)

bd <- rbind(o2, gs, spofo, motility, cell_length, cell_length, aggreg)

###############################

## Genome size, GC content, and Gene Number from NCBI...
genos <- fread('NCBI_prokaryotes.txt') %>% 
  filter(Status == 'Complete Genome') %>%
  transmute(
    sp = gsub("'|\\[|\\]", "", `Organism/Name`),
    Genus = sapply(sapply(sp, strsplit, ' '), function(x) x[[1]]),
    Species = sapply(sapply(sp, strsplit, ' '), function(x) x[[2]]),
    Genome_Mb = as.numeric(ifelse(`Size (Mb)` == '-', NA, `Size (Mb)`)),
    GC_content = `GC%`, 
    Gene_number = as.numeric(ifelse(Genes == '-', NA, Genes))) %>%
  gather(trait, val, Genome_Mb, GC_content, Gene_number) %>%
  select(-sp) %>%
  filter(!is.na(val))

####################
# 16s gene copy numbers from rrnDB

rrnDB <- read.csv('rrnDB-5.3.csv') %>%
  mutate(
    sp = gsub("'|\\[|\\]", "", NCBI.scientific.name),
    Genus = sapply(sapply(sp, strsplit, ' '), function(x) x[[1]]),
    Species = sapply(sapply(sp, strsplit, ' '), function(x) ifelse(length(x) > 1, x[[2]], NA)),
    trait = 'Copies_16S') %>%
  transmute(Genus, Species, trait, val = X16S.gene.count) %>%
  filter(!is.na(val)) %>%
  mutate(Species = ifelse(is.na(Species), 'unclassified', Species))

################ IgA
iga <- read.csv('Palm2014_IGA.csv', stringsAsFactors = FALSE)
iga <- t(iga)
colnames(iga) <- iga[1, ]
iga <- as.data.frame(iga[c(2:nrow(iga)), ])
iga[, c(4:ncol(iga))] <- sapply(iga[, c(4:ncol(iga))], function(x) as.numeric(as.character(x)))

iga <- iga %>%
  gather(tax, val, -var, -status, -subj) %>%
  spread(var, val) %>%
  mutate(tax = gsub("\\s|.__|\\[|\\]|Other", "", tax)) %>%
  separate(tax, sep = ';', c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  filter(!(Genus == 'unclassified' | Genus == '')) %>%
  mutate(Species = ifelse(Species == '', 'unclassified', Species)) %>%
  group_by(Genus, Species, trait = 'IgA') %>%
  summarise(val = log(mean(ici) + 1)) %>%
  as.data.frame()

#################
#b-vitamins
bvit <- read.csv('Bvitamins_Magnusdottir2015.csv', stringsAsFactors = FALSE) %>%
  mutate(Genus = sapply(sapply(tax, strsplit, ' '), function(x) x[[1]]),
         Species = sapply(sapply(tax, strsplit, ' '), function(x) x[[2]])) %>%
  mutate(Species = ifelse(Species == 'sp.', 'unclassified', Species)) %>%
  select(-NCBI.Taxonomy.ID, -tax, -Body.Site) %>%
  gather(Bvit, val, -Genus, -Species) %>%
  group_by(Genus, Species, trait = 'B_vitamins') %>%
  summarise(val = length(unique(Bvit[val > 0]))) %>%
  as.data.frame()

################################

#load traits from JGI
# currently I drop otu-level taxonomy, but I could probably figure it out with
# gg to genbank accession id: 
# ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5_accessions.txt.gz
# NCBI ID to genbank accession (a huge file):
# ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/ # not sure which one
jgi <- fread("GOLD-328_Organism_Metadata_02232018.txt", quote = "")

#remove duplicate oxygen column
jgi <- jgi[, -9]

#convert null to NA
jgi[jgi == '(null)'] <- NA

#clean columns
jgi <- jgi %>%
  transmute(
    Species = SPECIES, 
    Width = CELL_DIAMETER, 
    Length = CELL_LENGTH, 
    Gram_positive = GRAM_STAIN,
    Motility = MOTILITY,
    Oxygen_tolerance = OXYGEN_REQUIREMENT, 
    pH_optimum = PH,
    Spore = SPORULATION)

#Fix taxonomy
jgi <- separate(jgi, Species, c("Genus","Species"), extra = 'merge', fill = "right") %>%
  mutate(Species = gsub("sp. ", "", Species))

# deal with width - take averages of ranges when necessary.
jgi <- jgi %>%
  mutate(
    Width = ifelse(Width == '803nm', 0.803, Width),
    Width = gsub(" |[[:alpha:]]|�|�|\\?|\\&#956;|`", "", Width)
  ) %>%
  separate(Width, c("Width0","Width1"), fill = 'right', sep = '-') %>%
  mutate(Width = ifelse(is.na(Width0), Width0, 
                        ifelse(is.na(Width1), as.numeric(Width0),
                               (as.numeric(Width1) + as.numeric(Width0)) / 2))) %>%
  mutate(Width = as.numeric(Width)) %>%
  select(-Width0, -Width1)

# same with length
jgi <- jgi %>%
  mutate(
    Length = ifelse(Length == '803nm', 0.803, Length),
    Length = gsub(" |[[:alpha:]]|�|�|\\?|\\&#956;|`|\\&#8197;", "", Length),
    Length = gsub(".5.0", "5.0", Length, fixed = T),
    Length = gsub("1.55.0", "1.55", Length, fixed = T),
    Length = gsub("0.81.2", "0.8-1.2", Length, fixed = T)
    
  ) %>%
  separate(Length, c("Length0","Length1"), fill = 'right', sep = '-') %>%
  mutate(Length = ifelse(is.na(Length0), Length0, 
                         ifelse(is.na(Length1), as.numeric(Length0),
                                (as.numeric(Length1) + as.numeric(Length0)) / 2))) %>%
  mutate(Length = as.numeric(Length)) %>%
  select(-Length0, -Length1)

# gram positive, motility, oxygen tolerance, sporulation
jgi <- jgi %>%
  mutate(
    Gram_positive = c(1,0)[match(Gram_positive, c('Gram+','Gram-'))],
    Motility = c(1,1,0,0)[match(Motility, c('Motile','Chemotactic','Non-motile','Nonmotile'))],
    Oxygen_tolerance = c(5,5,4,3,2,1,1)[match(jgi$Oxygen_tolerance, 
                                              c('Obligate aerobe','Aerobe','Microaerophilic','Facultative','Facultative anaerobe','Anaerobe','Obligate anaerobe'))],
    Spore = c(1,1,0)[match(Spore, c('Non-sporulating','Nonsporulating','Sporulating'))])

#ph optimum
jgi <- jgi %>%
  mutate(
    pH_optimum = gsub(" |~", "", pH_optimum),
    pH_optimum = ifelse(pH_optimum %in% c('acido-sensible','Notknown'), NA, pH_optimum)) %>%
  separate(pH_optimum, c('pH0','pH1'), fill = 'right', sep = '-') %>%
  mutate(pH_optimum = ifelse(is.na(pH0), pH0, 
                             ifelse(is.na(pH1), as.numeric(pH0),
                                    (as.numeric(pH1) + as.numeric(pH0)) / 2))) %>%
  mutate(pH_optimum = as.numeric(pH_optimum)) %>%
  select(-pH0, -pH1)

jgi <- jgi %>%
  gather(trait, val, -Genus, -Species) %>%
  filter(!is.na(val))


###################################################

# put it all together
x <- bind_rows(
  mutate(ijsem, source = 'IJSEM'),
  mutate(spo, source = 'Browne2016'),
  mutate(bd, source = 'BacDrive'),
  mutate(genos, source = 'NCBI'),
  mutate(rrnDB, source = 'rrnDB'),
  mutate(iga, source = 'Palm2014'),
  mutate(bvit, source = 'Mag2015'),
  mutate(jgi, source = 'JGI')
)

uncleaned <- nrow(x)
### I looked at all the points plotted using 
#
#ggplot(x, aes(x = val)) + geom_density() + facet_wrap(~trait, scales = 'free')
# identifying and removing outliers
# remove inordinately large values...

##At some point, I would like/may want to use this slightly more stringent outlier filter
#x <- x %>%
#  filter(!(trait == 'Copies_16S' & val > 20)) %>%
#  filter(!(trait == 'GC_content' & (val > 80 | val < 20))) %>%
#  filter(!(trait == 'Gene_number' & val > 11000)) %>%
#  filter(!(trait == 'Genome_Mb' & val > 14)) %>%
#  filter(!(trait == 'Length' & (val > 12 | val <= 0))) %>%
#  filter(!(trait == 'pH_optimum' & val < 4)) %>%
#  filter(!(trait == 'Salt_optimum' & val > 25)) %>%
#  filter(!(trait == 'Temp_optimum' & val > 80)) %>%
#  filter(!(trait == 'Width' & (val > 4 | val <=0)))

x <- x %>%
  filter(!(trait == 'Copies_16S' & val > 20)) %>%
  filter(!(trait == 'GC_content' & (val > 80 | val < 20))) %>%
  filter(!(trait == 'Gene_number' & val > 11000)) %>%
  filter(!(trait == 'Genome_Mb' & val > 14)) %>%
  filter(!(trait == 'Length' & val > 30)) %>%
  filter(!(trait == 'pH_optimum' & val < 2.5)) %>%
  filter(!(trait == 'Salt_optimum' & val > 25)) %>%
  filter(!(trait == 'Temp_optimum' & val > 80)) %>%
  filter(!(trait == 'Width' & val > 8))

print(paste(uncleaned - nrow(x), "outliers of", nrow(x), "data points were removed"))

#for plotting source-wise comparisons
x1 <- x

#calculate means among species. Log length/width data
x <- x %>%
  mutate(Species = ifelse(Species == 'sp.', 'unclassified', Species)) %>%
  group_by(Genus, Species, trait) %>%
  summarise(val = mean(val, na.rm = T)) %>%
  mutate(val = ifelse(trait %in% c('Length','Width'), log(val), val))


####################################################################
# filter out all rows without species-level data that can be mapped onto greengenes OR
# doesn't exist in our data
# note, I lose a lot of data here; eventually, would like to use sequences to connect trait data to reference tree
mouse_tax <- readRDS('Reese_Tax_Aug2018.RDS')

x <- x %>%
  group_by(Genus, trait) %>%
  mutate(genus_mean = mean(val))


for(i in 1) {
  for(i in unique(x$trait)) {
    mouse_tax[[i]] <- ifelse(mouse_tax$Genus %in% x$Genus[x$trait == i],
      ifelse(paste(mouse_tax$Genus, mouse_tax$Species) %in% paste(x$Genus, x$Species)[x$trait == i],
             x$val[x$Genus == mouse_tax$Genus & x$Species == mouse_tax$Species & x$trait == i],
             x$genus_mean[x$Genus == mouse_tax$Genus & x$Species == mouse_tax$Species & x$trait == i]),
      NA)
  }
}



x <- x %>%
  spread(trait, val) %>%
  inner_join(ggtax[, c('Genus','Species','otu')], by = c('Genus', 'Species')) %>%
  mutate(otu = gsub('otu', '', otu)) %>%
  filter(Species != 'unclassified' | otu %in% tax$otu)


# Add a last minute manual entry
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3346390/
x$Aggregation_score[x$Genus == 'Bifidobacterium' & x$Species != 'unclassified'] <- 1 

# Merge spore scores (if no score and ijsem says 0, then 0; 
# otherwise the median of spore score when we know for ijsem spore == 1
x <- x %>% 
  mutate(
    Sporulation = ifelse(is.na(Spore_score), 
                         ifelse(Spore < 0.5, 0, median(Spore_score[tax$Spore > 0], na.rm = T)), Spore_score)) %>%
  ungroup() %>%
  select(-Spore, -Spore_score)

#add bugbase trait preditions
bugdat <- fread('default_traits_precalculated_bugbase.txt') %>%
  mutate(otu = as.character(V1), V1 = NULL)
  

x <- full_join(x, bugdat, by = 'otu')

#drop unwanted/redundant traits
x <- x %>%
  mutate(
    Aerobic = NULL,
    Aggregation_score = NULL,
    Contains_Mobile_Elements = NULL,
    Facultatively_Anaerobic = NULL,
    Genome_Mb = NULL,
    Gram_positive = rowMeans(cbind(Gram_positive, Gram_Positive), na.rm=TRUE),
    Gram_Positive = NULL,
    Gram_Negative = NULL,
    Oxygen_tolerance = NULL,
    pH_optimum = NULL,
    Potentially_Pathogenic = NULL,
    Salt_optimum = NULL,
    Stress_Tolerant = NULL,
    Width = NULL) %>%
  rename(
    Obligate_anaerobe = Anaerobic,
    Forms_biofilms = Forms_Biofilms)

saveRDS(x, file = 'traits_sparse.RDS')

#rm(list = ls())

#########PROTRAIT DATA SANDBOX
#prot <- fread('data\\ProTraits_binaryIntegratedPr0.90.txt')
#prot[prot == '?'] <- NA
#prot[, c(2:ncol(prot))] <- apply(prot[,c(2:ncol(prot))], 2, as.numeric)
#prot <- as.data.frame(prot)
#prot <- gather(prot, substrate, use, -Organism_name, -Tax_ID)
#prot <- filter(prot, !is.na(use))
#OTU ID doesn't work: table(tax$otu %in% paste0('otu',prot$Tax_ID))
# Genus species is a little better:
#table(paste(tax$Genus, tax$Species) %in% prot$Organism_name)



###########working on mouse tree stuff

library(castor)
library(phytools)
library(stringr)
library(ggtree)

tre <- read_tree(write.tree(read.tree('C:/Users/John/Desktop/LTP_mouse_clean.tree')))
tre$tip.label <- unlist(lapply(strsplit(tre$tip.label, '_'), function(x) x[[1]]))
tax_LTP <- read.table('C:/Users/John/Desktop/LTPs132_SSU.csv', sep = '\t')
tax_LTP <- tax_LTP[!duplicated(tax_LTP$V1), ]
tax_LTP <- transmute(tax_LTP, id = V1, 
                     Genus = unlist(lapply(strsplit(as.character(V5), ' ', fixed = TRUE), function(x) x[1])), 
                     Species = word(gsub("subsp.*", "", V5), 2),
                     binom = V5, tax = V10, notes = V9)
tax_me <- readRDS('C:\\Users\\John\\Google Drive\\Mouse_Microbiome_Shared_Files\\Data\\Reese_Tax_Aug2018.RDS')

iwant <- unique(tax_LTP$Genus)

tax <- bind_rows(transmute(tax_LTP, otu = id, Genus, Species, Binomial = binom) %>% mutate_all(as.character),
                 transmute(tax_me, otu, Genus, Species, Binomial = paste(Genus, Species)) %>% mutate_all(as.character)) %>%
    mutate(source = ifelse(grepl('OTU', otu), 'AR', 'LTP'))

genus_freqs <- table(tax$Genus) %>% as.data.frame() %>% arrange(desc(Freq))
iwant <- as.character(tax$otu[tax$Genus %in% genus_freqs$Var1[c(2:5)]])
iwant <- c(tax$otu[tax$source == 'AR'], sample(tax$otu[tax$source != 'AR'], sum(tax$source == 'AR')))





tmp <- drop.tip(tre, tre$tip.label[!tre$tip.label %in% iwant])
tax_tmp <- tax[match(tmp$tip.label, tax$otu), ] %>%
  select()

#color by data source
tmp <- groupOTU(tmp, tax_tmp)
tmp <- groupOTU(tmp, grep('OTU', tmp$tip.label))
attr(tmp, 'group') <- factor(ifelse(grep('OTU', tmp$tip.label), 'Mouse', 'LTP'))
#tmp <- groupOTU(tmp, match(tax$otu[tax$Genus == 'Streptomyces'], tmp$tip.label))

#tmp <- groupOTU(tmp, grupe)

#color by Genus
#tmp_genera <- as.character(tax$Genus[match(tmp$tip.label, tax$otu)])
#tmp <- groupOTU(tmp, sample(c(1:3), length(tmp$tip.label), replace = TRUE))
#ggtree(tmp, aes(color=group))
ggtree(tmp, aes(color=group), layout='circular') +
  theme(legend.position = c(0.5,0.4), 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys
