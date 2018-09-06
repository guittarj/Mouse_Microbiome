#Mouse gut microbiome traits -- inferring unknown traits using phylogeny
#John Guittar

#load packages and set working directory
wd <- 'C:\\Users\\John\\Google Drive\\Mouse_Microbiome_Shared_Files\\Data\\'
setwd(wd)
source('C:\\Users\\John\\Documents\\msu\\Mouse_Microbiome\\custom_functions.R')
loadpax(c('tidyverse','castor','phytools','stringr','ggtree','phylosignal','kableExtra'))

#load tree with LTP and mouse tips, calculated with usearch
tre <- read_tree(write.tree(read.tree('LTP_mouse.tree')))

#load traits
traits <- as.data.frame(readRDS('traits_sparse.RDS')) %>%
  mutate(binomial = paste(Genus, Species))

trait_names <- c(
  "Aggregation_score" = 'Aggregation score',
  "B_vitamins" = "B vitamins",       
  "Copies_16S" = "16S gene copies",
  "GC_content" = "GC content",       
  "Gene_number" = "Genes",
  "Genome_Mb" = "Genome size",
  "Gram_positive" = "Gram-positive",     
  "IgA" = 'IgA binding affinity',             
  "Length" = "Length",
  "Motility" = "Motility",
  "Oxygen_tolerance" = "Oxygen tolerance",
  "pH_optimum" = "pH optimum",
  "Salt_optimum" = "Salinity optimum",      
  "Sporulation" = "Sporulation score",      
  "Temp_optimum" = "Temperature optimum",
  "Width" = "Width"
)

#load taxonomy file (used to link otus and traits)
#files were created in mouse_processing_trait_data.R
tax_LTP <- readRDS('LTP_tax.RDS')
tax_mouse <- readRDS('mouse_tax_SILVA.RDS') %>%
  select(otu, Genus, Species)
tax <- bind_rows(tax_LTP, tax_mouse) %>%
  mutate(binomial = paste(Genus, Species))

#load OTU table for plotting abundance-weighted coverage
otus <- readRDS('otu_table_usearch_aug19.RDS')

#join traits and OTUs using curated latin binomials
tax_traits <- tax %>% left_join(select(traits, -Genus, -Species), by = c("binomial"))

##Estimating unknown trait values, and evaluating confidence in those estimations
#loop through each trait...
if (TRUE) {
  for(i in trait_names) {
    
    #identify trait values of each tip
    tip_states <- tax_traits[[i]][match(tre$tip.label, tax_traits$otu)]
    
    #estimate unknown trait values
    HSP <- hsp_independent_contrasts(tre, tip_states)
    SBT <- hsp_subtree_averaging(tre, tip_states)
    SCP <- hsp_squared_change_parsimony(tre, tip_states)
    
    #save temporary data for this trait 
    dat_trait <- data.frame(
      otu = tre$tip.label, 
      trait = i,
      HSP = HSP$states[1:length(tip_states)],
      SBT = SBT$states[1:length(tip_states)], 
      SCP = SCP$states[1:length(tip_states)], 
      stringsAsFactors = FALSE)
    
    #identify nearest otu with directly observed trait data
    nt <- find_nearest_tips(tre, target_tips = tre$tip.label[!is.na(tip_states)])
    dat_trait$nearest_measured_otu <- tre$tip.label[nt$nearest_tip_per_tip]
    dat_trait$dist <- nt$nearest_distance_per_tip
    
    ##manually calculate confidence of estimates
    #create tree with only measured tips
    tre_tmp <- drop.tip(tre, tre$tip.label[is.na(tip_states)])
    
    #reduce the total number of potential tips to 2000 because getting all possible pairs is overkill,
    #labs <- tre_tmp$tip.label
    #if (length(labs)>2000) tre_tmp <- drop.tip(tre_tmp, labs[sample(seq_along(labs), length(labs) - 2000)])
    
    #create list of all possible combinations of measured otus, and their distances
    obs_trait <- as.data.frame(t(combn(tre_tmp$tip.label, 2)), stringsAsFactors = FALSE) %>%
      rename(otu1 = V1, otu2 = V2)
    j <- get_pairwise_distances(tre_tmp, obs_trait$otu1, obs_trait$otu2, 
                                as_edge_counts=FALSE, check_input=FALSE)
    
    #attach trait data
    #filter OTU pairs to only those with the shortest distances > 0
    #calculate difference in trait values
    obs_trait <- data.frame(obs_trait, dist = j) %>%
      filter(dist > 0) %>%
      select(otu1, otu2, dist) %>%
      mutate(
        trait = i,
        dist_bin = dist %/% 0.001 * 0.001,
        delta = abs(tax_traits[[i]][match(otu1, tax_traits$otu)] -
                      tax_traits[[i]][match(otu2, tax_traits$otu)]))
    
    #calculate null distances (i.e. expected delta by randomly selecting two observed OTUs)
    null_trait <- obs_trait %>%
      group_by(trait, dist_bin) %>%
      summarise(mean = mean(delta))
    
    #remove pairwise distances larger than dist <= 0.15, because we already know those are 'random')
    # drop data beyond 1000 randomly selected pairwise differences for each 0.001 increment in dist
    obs_trait <- obs_trait %>%
      filter(dist <= 0.2) %>%
      group_by(dist_bin) %>%
      filter(seq_along(dist) %in% sample(seq_along(dist), min(length(dist), 1e4))) %>%
      ungroup()
    
    if (i == trait_names[1]) {
      nulls <- null_trait
      dat <- dat_trait 
      obs <- obs_trait
    } else {
      nulls <- bind_rows(nulls, null_trait)
      dat <- bind_rows(dat, dat_trait)
      obs <- bind_rows(obs, obs_trait)
    }
    
    print(paste('done with', i)); flush.console()
    
    if (i == tail(trait_names, 1)) {
      saveRDS(nulls, 'delta_nulls.RDS')
      saveRDS(dat, 'trait_predictions.RDS')
      saveRDS(obs, 'phylogenetic_and_trait_distances.RDS')
    }
  }
} else {
  nulls <- readRDS('delta_nulls.RDS')
  dat <- readRDS('trait_predictions.RDS')
  obs <- readRDS('phylogenetic_and_trait_distances.RDS')
}

#Pearson correlations among exploring differences among different methods of hidden state character predictions
kable(cor(dat[, c('HSP','SBT','SCP')]), caption = 'Pearson correlations among methods of hidden character predictions')

dat %>% 
  gather(method, val, HSP, SBT, SCP) %>%
  ggplot(aes(x = val, fill = method, alpha = 0.5)) + 
  geom_density() + 
  facet_wrap(~trait, scales = 'free') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank())

# before fitting models, I revert to a condensed version of Gram-positive data, in which I take the average delta for each 0.001 increment of dist (i.e., each 'dist_bin'). I do this because the logistic model (obviously the right model) wasn't converging.
mods <- obs %>%
  filter(trait != 'Gram_positive') %>%
  bind_rows(obs %>% 
              filter(trait == 'Gram_positive') %>% 
              group_by(trait, dist = dist_bin) %>% 
              mutate(delta = mean(delta)) %>% 
              distinct(dist, delta, .keep_all = TRUE))

#fit various models to delta~dist for each trait
#try_default ensures that do() finishes even if the model does not converge, or errors
mods <- mods %>%
  group_by(trait) %>%
  do(
    null = lm(delta ~ 1, data = .),
    linear = lm(delta ~ dist, data = .),
    linear_asym = plyr::try_default(nls(delta~SSasympOrig(dist, Asym, lrc), data = .),
                                    default = NA, quiet = TRUE),
    log = lm(delta~log(dist), data = .),
    logistic = plyr::try_default(nls(delta ~ SSlogis(dist, Asym, xmid, scal), data = .), default = NA, quiet = TRUE))

#model names for plotting later
mod_names <- c(null = 'Null', linear = 'Linear regression', linear_asym = 'Asymptotic regression',
               log = 'Logarithmic regression', logistic = 'Logistic regression')


#remove failed models 
mods <- mods %>%
  gather(type, mod, -trait) %>%
  rowwise() %>%
  filter(class(mod) != 'logical')

#calculate AIC and coefficients to see if the model was so bad that it predicted a negative slope
#if coefficient is negative, only keep null model
#normalize AIC values to lowest value for each set of models
mods <- mods %>%
  mutate(
    AIC = AIC(mod),
    est = ifelse(class(mod) == 'nls', summary(mod)$coef[[3]], summary(mod)$coef[[2]])) %>%
  ungroup() %>%
  mutate(type = mod_names[match(type, names(mod_names))]) %>%
  group_by(trait) %>%
  filter(if (min(est) < 0) type == 'Null' else est > 0) %>%
  mutate(AIC = AIC - min(AIC)) %>%
  arrange(trait, AIC)

#create dataframe/table with AIC scores for all models
mod_aic <- mods %>%
  select(Trait = trait, type, AIC) %>%
  spread(type, AIC) %>%
  select(Trait, Null, `Linear regression`, `Asymptotic regression`, `Logarithmic regression`, `Logistic regression`)

#create list of best mods
best_mods <- mods %>%
  filter(AIC == 0) %>% 
  mutate(best = paste(trait, type)) %>% 
  pull(sort(best))

#predict fits for models
#identify whether each model is the best fitting
#drop any foolish deltas less than 0
preds <- mods %>%
  rowwise() %>%
  do(data.frame(stringsAsFactors = FALSE,
                trait = .$trait,
                type = .$type,
                dist = seq(0, 0.2, length.out = 1000),
                delta = predict(.$mod, data.frame(dist = seq(0, 0.2, length.out = 1000))))) %>%
  filter(delta >= 0) %>%
  mutate(best = paste(trait, type) %in% best_mods) %>%
  ungroup()

#create a dataframe of null expectations of delta for each trait
#To remove phylogenetic effects, null delta is calculated as the mean pairwise trait-based distance for all randoly selected pairs of otus with greater than 10% difference in their 16S V4 region -- except for Salt optimum and ph optimum, which have no observable phylogenetic effect
nulls_tmp <- nulls %>%
  group_by(trait) %>%
  filter(dist_bin > 0.1 | trait %in% c('Salt_optimum','pH_optimum')) %>%
  summarise(null = mean(mean))

#append null expectations for each trait
#create a delta threshold at the point when the model reaches the delta null expectation
#identitify the maximum phylogenetic distance that falls below that threshold (I included a zero to avoid -Inf results)
preds_bin <- preds %>%
  group_by(trait, type) %>%
  mutate(
    dist_bin = dist %/% 0.001 * 0.001,
    threshold = nulls_tmp$null[match(trait, nulls_tmp$trait)]*.9,
    max_dist = max(dist[delta < threshold], 0),
    max_dist = ifelse(max_dist == max(dist), 0, max_dist)) %>%
  ungroup() %>%
  mutate(trait = trait_names[match(trait, names(trait_names))])

#rather than plot the 8M+ pairwise deltas, which crashes my computer, and looks messy...
#calculate means and sd. Add a row for null deltas. Clean up trait names.
obs_tmp <- obs %>% 
  mutate(dist_bin = dist %/% 0.001 * 0.001) %>% 
  group_by(trait, dist_bin) %>%
  summarise(
    mean = mean(delta), 
    sd = sd(delta), 
    null = nulls_tmp$null[match(trait[1], nulls_tmp$trait)]) %>%
  ungroup() %>%
  mutate(trait = trait_names[match(trait, names(trait_names))])

#plot without standard deviations
ggplot(obs_tmp, aes(x = dist_bin, y = mean)) +
  geom_point() +
  geom_line(aes(y = delta, color = type), data = filter(preds_bin, best & type != 'Null'), lwd = 1.5) +
  geom_hline(aes(yintercept = null), linetype = 3) +
  geom_vline(aes(xintercept = max_dist), data = filter(preds_bin, best), lty = 2) +
  scale_color_discrete(name = 'Best model') +
  facet_wrap(~trait, scales = 'free') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank()) +
  labs(x = 'Phylogenetic distance between tips', y = 'Mean difference in traits')

#plot with standard deviations
ggplot(obs_tmp, aes(x = dist_bin, y = mean)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0, color = 'grey') +
  geom_point() +
  geom_line(aes(y = delta, color = type), data = filter(preds_bin, best & type != 'Null'), lwd = 1.5) +
  geom_hline(aes(yintercept = null), linetype = 3) +
  geom_vline(aes(xintercept = max_dist), data = filter(preds_bin, best), lty = 2) +
  scale_color_discrete(name = 'Best model') +
  facet_wrap(~trait, scales = 'free') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank()) +
  labs(x = 'Phylogenetic distance between tips', y = 'Mean difference in traits')

#data frame of maximum distances used to infer traits, for each trait
max_dists <- preds_bin %>% 
  filter(best) %>% 
  select(trait, max_dist) %>% 
  distinct()

kable(max_dists, caption = "Phylogenetic distances at which inference becomes meaningless")

#determine if nearest_measured_otu is close enough to phylogenetically infer trait values
#append overall OTU abundances
tmp <- dat %>%
  filter(grepl("OTU", otu)) %>%
  mutate(
    trait = trait_names[match(trait, names(trait_names))],
    max_dist = max_dists$max_dist[match(trait, max_dists$trait)],
    Type = ifelse(dist > max_dist, "Insufficient data",
                  ifelse(dist > 0, "Phylogenetically inferred", 'Directly observed')),
    Type = factor(Type, levels = c('Insufficient data','Phylogenetically inferred','Directly observed'))) %>%
  left_join(otus %>% group_by(otu) %>% summarise(abun = sum(abun)), by = 'otu') %>%
  group_by(trait) %>%
  mutate(relabun = abun / sum(abun)) %>%
  arrange(trait, Type, relabun)

#summary barplot (number of OTUS by trait/data type)
ggplot(tmp, aes(x = trait, fill = Type)) + 
  geom_bar(color = 'black') +
  coord_flip() +
  scale_fill_manual(values = c('lightgrey','skyblue','blue')) +
  labs(x = "Number of OTUs", y = '') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE))

#summary barplot (abundances of OTUS by trait/data type)
ggplot(tmp, aes(x = trait, y = relabun, fill = Type, color = Type)) + 
  geom_bar(stat = 'identity') +
  coord_flip() +
  scale_fill_manual(values = c('lightgrey','#6666ff','blue')) +
  scale_color_manual(values = c('#d9d9d9','#b3b3ff','#6666ff'), guide = FALSE) +
  labs(x = "Percent Relative abundance across all samples", y = '') +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE))
