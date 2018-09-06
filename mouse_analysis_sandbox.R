###############################################################
################# PLAYING WITH TREES #####################
###############################################################

#load taxonomies for labeling tree and exploring results.

tax$source <- ifelse(paste(tax$Genus, tax$Species) %in% paste(tax_mouse$Genus, tax_mouse$Species), 'AR', 'LTP')
iwant <- c(tax$otu[tax$source == 'AR'], sample(tax$otu[tax$source != 'AR'], 3*sum(tax$source == 'AR')))
tmp <- drop.tip(tre, tre$tip.label[!tre$tip.label %in% iwant])

#color by data source
tmp <- groupOTU(tmp, grep('OTU', tmp$tip.label))

#color by Genus
ggtree(tmp, aes(color=group))
ggtree(tmp, aes(color=group), layout='circular') +
  theme(legend.position = c(0.5,0.4), 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys
