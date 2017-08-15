##workflows:
#https://f1000research.com/articles/5-1492/v1
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html

#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

setwd("C:\\Users\\arisely\\Dropbox\\PhD\\Microbiome\\CHAPTER 3 - PHYSIOLOGY\\ANALYSIS\\R\\Phyloseq analysis")

library(phyloseq)

##other packages you may need, althoguh I think a few are included with the phyloseq package so you may not need to load them

library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(ape)
library(gridExtra)
library(ade4)
library(plyr)
library(tidyr)
library(data.table)
library(stringr)





source("miseqR.R")

theme_set(theme_bw())


################# ONLY OTUS WITH OVER 10 SEQUENCES (AND WITH PHYLO TREE) ##################################################################

setwd("C:\\Users\\arisely\\Dropbox\\PhD\\Microbiome\\CHAPTER 3 - PHYSIOLOGY\\ANALYSIS\\R\\Phyloseq analysis")

# Assign variables for imported data

sharedfile = "shorebird.microbiome.shared"
taxfile = "shorebird.microbiome.taxonomy"

##import mothur shared and tax files

mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

##import metadata
map <- read.csv("shorebird.metadata.csv", header=T, row.names=1)

##make meta data into phyloseq format

map <- sample_data(map)
str(map)

head(map)
tail(map)



##merge the metadata into the phyloseq object

moth_merge <- merge_phyloseq(mothur_data, map)

##import phylogenetic tree

tree<-read.tree("shorebird.microbiome.tree")

##make tree into phyloseq format

tree<-phy_tree(tree)

##merge tree into current phyloseq object so all object now contained in moth_merge

moth_merge <-merge_phyloseq(moth_merge, tree)

moth_merge
sample_data(moth_merge)
##everything looks fie and no errors!

##rename columns to relevant taxonomic group

colnames(tax_table(moth_merge))

colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", 
                                     "Order", "Family", "Genus")

data_all <- moth_merge %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "mitochondria" &
      Class   != "Chloroplast"
  )

data_all




################################################################look at read library


sample_sum_df <- data.frame(sum = sample_sums(data_all))
plot(sample_sum_df$sum)


ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  ylab("Frequency")+
  theme(axis.title.y = element_blank())

##Broome samples didn't work :()

##############################################################################################################################

##get rid of 9185 as lowest read count
##get rid of 8688 as mostly unclassified taxa... not sure why

data_all<-data_all %>%
  subset_samples(Group != "Env_Broome_A")

data_all<-data_all %>%
  subset_samples(Group != "Env_Broome_E")


data_all<-prune_taxa(taxa_sums(data_all)>0, data_all)

##remove OTUs that are common in the neg controls
neg_control<-data_all %>% subset_samples(type == "neg")

neg_control <- prune_taxa(taxa_sums(neg_control) > 5, neg_control)
neg_control

otu_table(neg_control)


##97 OTUs which have over 1 reads

##we should get rid of these

badtaxa<-taxa_names(neg_control)
alltaxa<-taxa_names(data_all)
alltaxa1 <- alltaxa[!(alltaxa %in% badtaxa)]

data_all = prune_taxa(alltaxa1, data_all)


##get rid of negative control altogether


data_all<-data_all%>%subset_samples(species!="Neg")

#delete replicates

data_all<-data_all%>%subset_samples(remove_replicate=="No")
sample_data(data_all)

##just birds

birds<-data_all%>%subset_samples(type=="bird")
birds <- prune_taxa(taxa_sums(birds) > 0, birds)

smin <- min(sample_sums(birds))

##rarefy
birds <-rarefy_even_depth(birds, sample.size = min(sample_sums(birds)),
                         rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)




##for migration analysis
migration<-birds%>%subset_samples(species=="Red-necked stint"|species=="Curlew sandpiper")
migration<-migration%>%subset_samples(month=="August"|month=="September")
migration <- prune_taxa(taxa_sums(migration) > 0, migration)
migration

sample_sum_df <- data.frame(sum = sample_sums(migration))
plot(sample_sum_df$sum)
sum(sample_sum_df)


##QUESTION  How does migration effect the microbiome? Red necked stint and Curlew Sandpiper
curlew.b<-birds%>%subset_samples(species=="Curlew sandpiper"& site=="Broome")
curlew.b<- prune_taxa(taxa_sums(curlew.b) > 0, curlew.b)
curlew.b.a<-curlew.b%>%subset_samples(age=="Adult")
curlew.b.a<- prune_taxa(taxa_sums(curlew.b.a) > 0, curlew.b.a)
curlew.b.sy<-curlew.b%>%subset_samples(age=="Second year")
curlew.b.sy<- prune_taxa(taxa_sums(curlew.b.sy) > 0, curlew.b.sy)

stint.b<-birds%>%subset_samples(site=="Broome" & species=="Red-necked stint")
stint.b<- prune_taxa(taxa_sums(stint.b) > 0, stint.b)
stint.b.a<-stint.b%>%subset_samples(age=="Adult")
stint.b.a<- prune_taxa(taxa_sums(stint.b.a) > 0, stint.b.a)
stint.b.sy<-stint.b%>%subset_samples(age=="Second year")
stint.b.sy<- prune_taxa(taxa_sums(stint.b.sy) > 0, stint.b.sy)
stint.f<-birds%>%subset_samples(month=="September")
stint.f<-stint.f%>%subset_samples(species=="Red-necked stint")
stint.f<- prune_taxa(taxa_sums(stint.f) > 0, stint.f)
stint.f.sy<-stint.f%>%subset_samples(age=="Second year")

stint.f.sy<- prune_taxa(taxa_sums(stint.f.sy) > 0, stint.f.sy)
stint.f.a<-stint.f%>%subset_samples(age=="Adult")
stint.f.a<- prune_taxa(taxa_sums(stint.f.a) > 0, stint.f.a)

##########################################  Prevalence  ##############################################
##prevance whole group

prev0 = apply(X = otu_table(migration),
              MARGIN = ifelse(taxa_are_rows(migration), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(migration),
                    tax_table(migration))
keepPhyla = table(prevdf$Phylum)[(table(prevdf$Phylum) > 5)]
prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))

prevdf$Prevalence<-(prevdf$Prevalence/77)*100


head(prevdf)
str(prevdf)
summary(prevdf$Prevalence)
###only keep those phyla which occur in at least 5 samples

prevdf1


# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(migration)
prevalenceThreshold

migration1 = prune_taxa((prev0 > prevalenceThreshold), migration)
migration1

ggplot(prevdf1, aes(TotalAbundance, Prevalence)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum)+theme_bw() +
  theme(  axis.text.y = element_text(size=12), axis.title=element_text(size=14))

write.csv(prevdf, "prevalence_all_samples.csv")

##for supplementary table on core microbiomes for each group

##########################################################################

##for each group 

prev0 = apply(X = otu_table(stint.f.a),
              MARGIN = ifelse(taxa_are_rows(stint.f.a), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(stint.f.a),
                    tax_table(stint.f.a))


prevdf<-prevdf[order(-prevdf$TotalAbundance),]
head(prevdf)

write.csv(prevdf, "stint.f.a.csv")

##############################

prev0 = apply(X = otu_table(stint.f.sy),
              MARGIN = ifelse(taxa_are_rows(stint.f.sy), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(stint.f.sy),
                    tax_table(stint.f.sy))


prevdf<-prevdf[order(-prevdf$TotalAbundance),]
head(prevdf)

write.csv(prevdf, "stint.f.sy.csv")

########################

prev0 = apply(X = otu_table(stint.b.sy),
              MARGIN = ifelse(taxa_are_rows(stint.b.sy), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(stint.b.sy),
                    tax_table(stint.b.sy))


prevdf<-prevdf[order(-prevdf$TotalAbundance),]
head(prevdf)

write.csv(prevdf, "stint.b.sy.csv")

########################

prev0 = apply(X = otu_table(stint.b.a),
              MARGIN = ifelse(taxa_are_rows(stint.b.a), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(stint.b.a),
                    tax_table(stint.b.a))


prevdf<-prevdf[order(-prevdf$TotalAbundance),]
head(prevdf)

write.csv(prevdf, "stint.b.a.csv")

##########################


prev0 = apply(X = otu_table(curlew.b.a),
              MARGIN = ifelse(taxa_are_rows(curlew.b.a), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(curlew.b.a),
                    tax_table(curlew.b.a))


prevdf<-prevdf[order(-prevdf$TotalAbundance),]
head(prevdf)

write.csv(prevdf, "curlew.b.a.csv")

#######################


prev0 = apply(X = otu_table(curlew.b.sy),
              MARGIN = ifelse(taxa_are_rows(curlew.b.sy), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(curlew.b.sy),
                    tax_table(curlew.b.sy))


prevdf<-prevdf[order(-prevdf$TotalAbundance),]
head(prevdf)

write.csv(prevdf, "curlew.b.sy.csv")

########################create lefse files########################################################################

##broome stint
stint.b.lefse <- stint.b %>% transform_sample_counts(function(x) {x/sum(x)} )

stint.b.lefse1<-data.frame(otu_table(stint.b.lefse))
write.csv(stint.b.lefse1, "stint.b.lefse.csv")

stint.b.tax<-data.frame(tax_table(stint.b.lefse))
write.csv(stint.b.tax, "stint.b.tax.csv")

##flinders sep  stint

stint.f.sep<-stint%>%subset_samples(month=="September")
stint.f.sep<- prune_taxa(taxa_sums(stint.f.sep) > 0, stint.f.sep)

stint.f.lefse <- stint.f.sep %>% transform_sample_counts(function(x) {x/sum(x)} )

stint.f.lefse1<-data.frame(otu_table(stint.f.lefse))
write.csv(stint.f.lefse1, "stint.f.lefse.csv")

stint.f.tax<-data.frame(tax_table(stint.f.lefse))
write.csv(stint.f.tax, "stint.f.tax.csv")

#curlew broome

curlew.b<-curlew%>%subset_samples(site=="Broome")
curlew.b<- prune_taxa(taxa_sums(curlew.b) > 0, curlew.b)
curlew.b.lefse <- curlew.b %>% transform_sample_counts(function(x) {x/sum(x)} )
curlew.b.lefse1<-data.frame(otu_table(curlew.b.lefse))
write.csv(curlew.b.lefse1, "curlew.b.lefse.csv")

curlew.b.tax<-data.frame(tax_table(curlew.b.lefse))
write.csv(curlew.b.tax, "curlew.b.tax.csv")



#############################################stacked barplots #################################################################


#################### migration

migration_phylum <- migration %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

#filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
head(migration_phylum)


###############################################################

###########plot

migration_phylum$Phylum<-factor(migration_phylum$Phylum)
#broome_phylum$Order<-factor(broome_phylum$Order)
unique(migration_phylum$Phylum )
migration_phylum$Phylum<-factor(migration_phylum$Phylum)



migration_phylum$Sample<-factor(migration_phylum$Sample, levels = migration_phylum$Sample[order(migration_phylum$species, migration_phylum$site )])

unique(migration_phylum$Phylum)
migration_phylum$Phylum<-factor(migration_phylum$Phylum, level = c("Actinobacteria", 
                                                                   "Firmicutes",
                                                                   "Proteobacteria",
                                                                   "Bacteroidetes",
                                                                   "Deferribacteres",
                                                                   "Fusobacteria",
                                                                   "Tenericutes",        
                                                                   "Spirochaetae" ,
                                                                   "Deinococcus-Thermus",
                                                                   "Acidobacteria",
                                                                   "Armatimonadetes",
                                                                   "Bacteria_unclassified",
                                                                   "Chloroflexi",
                                                                   "Cyanobacteria",
                                                                   "FBP",
                                                                   "Fibrobacteres",
                                                                   "Gemmatimonadetes",
                                                                   "Gracilibacteria",
                                                                   "Ignavibacteriae",
                                                                   "Lentisphaerae",
                                                                   "Microgenomates",
                                                                   "Parcubacteria",
                                                                   "Peregrinibacteria",       
                                                                   "Planctomycetes",
                                                                   "Saccharibacteria",                   
                                                                   "SR1_(Absconditabacteria)",
                                                                   "TM6_(Dependentiae)",
                                                                   "Verrucomicrobia")) 
                                                                  
                                                                   
                                                                


unique(migration_phylum$Sample)
colors <- c("navyblue","green4","orchid","royalblue","black","red","brown","grey","bisque",
            "white","white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white")


names(colors) <- levels(migration_phylum$Phylum)
colScale <- scale_fill_manual(name = "Phylum",values = colors)

migration_phylum$Sample<-factor(migration_phylum$Sample, level = c("8682", "8629", "8615", "8685", "8610", "8662", "8560", "8607", "8671", "8589", "8657", "8573", "8625", "8627", "8639", "8549", "8645", "8632", "8660", "8649", "8601",
                                                                   "8658", "8638", "8646", "8634", "8606", "8612", "8581", "8618", "8609", "8628", "8672", "8605", "8689", "8659", "8616", "8626", "8604", "8569", "8677", "8588", "8675",
                                                                   "8551", "8613", "8546", "8683", "8590", "8557", "8558", "8674", "8559", "8680", "8630", "8563", "8586", "8556", "8561", "8644", "8667", "8566", "8665", "8633", "8564",
                                                                   "8583", "8575", "8571", "8565", "8578", "8570", "8684", "8688", "8655", "8552", "8548", "8679", "8663", "8686"))


ggplot(migration_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=12), axis.title=element_text(size=14))+
  facet_wrap(species~site~age, scales="free", nrow=3, ncol=2)+colScale


############repeat but for Family

migration_family <- migration %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa 
  arrange(Family)                                      # Sort data frame alphabetically by family
#                          
head(migration_family)

migration_family$Family<-factor(migration_family$Family)
#broome_phylum$Order<-factor(broome_phylum$Order)
unique(migration_family$Family)

#write.csv(migration_family,"family.csv")

##all common Corynebacteriales_unclassified are genus Corynebacterium (we ran sequences through ARB):

migration_family$Family[migration_family$Family=="Corynebacteriales_unclassified"]<-"Corynebacteriaceae"


###########plot
migration_family$Family<-factor(migration_family$Family, level = c("Corynebacteriaceae",  
                                                                   "Acidaminococcaceae","Clostridiaceae_1","Clostridiales_unclassified" , "Enterococcaceae" ,"Erysipelotrichaceae"  ,"Firmicutes_unclassified" , "Lachnospiraceae"  ,"Peptostreptococcaceae" , "Ruminococcaceae" , "Staphylococcaceae" ,
                                                                   "Comamonadaceae", "Desulfovibrionaceae","Enterobacteriaceae","Gammaproteobacteria_unclassified" ,"Helicobacteraceae"  ,  "Methylobacteriaceae" ,"Oxalobacteraceae" ,"Rhizobiales_unclassified" , "Rhodobacteraceae",   "Sphingomonadaceae" , "Succinivibrionaceae",
                                                                   "Bacteroidaceae" , "Bacteroidales_unclassified" ,"Chitinophagaceae" ,"Flavobacteriaceae","Porphyromonadaceae","Rikenellaceae" ,
                                                                   "Deferribacteraceae" ,
                                                                   "Fusobacteriaceae",
                                                                   "Anaeroplasmataceae", "Mycoplasmataceae",
                                                                   "Brachyspiraceae" , 
                                                                   "Deinococcaceae",
                                                                   
                                                                   "0319-6G20",	
                                                                   "Acetobacteraceae",	
                                                                   
                                                                   "Acidimicrobiaceae",	
                                                                   "Acidimicrobiales_unclassified",	
                                                                   "Acidobacteriaceae_(Subgroup_1)",	
                                                                   "Actinobacteria_unclassified",	
                                                                   "Actinomycetaceae",	
                                                                   "Aeromonadaceae",	
                                                                   "Aeromonadales_unclassified",	
                                                                   "Alcaligenaceae",	
                                                                   "Alphaproteobacteria_unclassified",	
                                                                   "Alteromonadaceae",	
                                                                   "Anaerolineaceae",	
                                                                   
                                                                   "Anaplasmataceae",	
                                                                   "Ardenticatenales_fa",	
                                                                   "Arenicellaceae",	
                                                                   "Armatimonadales_unclassified",	
                                                                   "Aurantimonadaceae",	
                                                                   "Bacillaceae",	
                                                                   "Bacilli_unclassified",	
                                                                   "Bacteria_unclassified",	
                                                                   "Bacteriovoracaceae",	
                                                                   
                                                                   "Bacteroidales_S24-7_group",	
                                                                   
                                                                   "Bacteroidetes_unclassified",	
                                                                   "Bartonellaceae",	
                                                                   "BD2-11_terrestrial_group_fa",	
                                                                   "BD7-8_marine_group_fa",	
                                                                   "Bdellovibrionaceae",	
                                                                   "Beijerinckiaceae",	
                                                                   "Bifidobacteriaceae",	
                                                                   "Blastocatellaceae_(Subgroup_4)",	
                                                                   
                                                                   "Bradymonadaceae",	
                                                                   "Bradymonadales_fa",	
                                                                   "Bradyrhizobiaceae",	
                                                                   "Brevibacteriaceae",	
                                                                   "Brucellaceae",	
                                                                   "Burkholderiaceae",	
                                                                   "Burkholderiales_unclassified",	
                                                                   "Caldilineaceae",	
                                                                   "Campylobacteraceae",	
                                                                   "Campylobacterales_unclassified",	
                                                                   "Cardiobacteriaceae",	
                                                                   "Carnobacteriaceae",	
                                                                   "Caulobacteraceae",	
                                                                   "Cellulomonadaceae",	
                                                                   "Cellvibrionaceae",	
                                                                   "Cellvibrionales_unclassified",	
                                                                   
                                                                   "Chloroflexaceae",	
                                                                   "Chloroflexi_unclassified",	
                                                                   "Chromatiaceae",	
                                                                   "Clostridia_unclassified",	
                                                                   
                                                                   "Clostridiaceae_4",	
                                                                   
                                                                   "Clostridiales_vadinBB60_group",	
                                                                   
                                                                   "Coriobacteriaceae",	
                                                                   
                                                                   "Coxiellaceae",	
                                                                   "Cryomorphaceae",	
                                                                   "CS-B046_fa",	
                                                                   "Cyanobacteria_unclassified",	
                                                                   "Cyclobacteriaceae",	
                                                                   "Cytophagaceae",	
                                                                   "Cytophagales_unclassified",	
                                                                   "d142_fa",	
                                                                   "DA111",	
                                                                   "DBS1",	
                                                                   
                                                                   
                                                                   "Deltaproteobacteria_unclassified",	
                                                                   "Dermabacteraceae",	
                                                                   "Dermacoccaceae",	
                                                                   "Dermatophilaceae",	
                                                                   "Desulfobacteraceae",	
                                                                   "Desulfobulbaceae",	
                                                                   
                                                                   "Desulfovibrionales_unclassified",	
                                                                   "Desulfuromonadaceae",	
                                                                   "Desulfuromonadales_unclassified",	
                                                                   "DEV007",	
                                                                   "Dietziaceae",	
                                                                   "Eel-36e1D6",	
                                                                   "Elev-16S-1332",	
                                                                   
                                                                   
                                                                   "Entomoplasmatales_Incertae_Sedis",	
                                                                   "Epsilonproteobacteria_unclassified",	
                                                                   
                                                                   "Erythrobacteraceae",	
                                                                   "Eubacteriaceae",	
                                                                   "Euzebyaceae",	
                                                                   "Family_XI",	
                                                                   "Family_XII",	
                                                                   "Family_XIII",	
                                                                   "FamilyI",	
                                                                   "FamilyII",	
                                                                   "FBP_fa",	
                                                                   "FD035",	
                                                                   
                                                                   "Flammeovirgaceae",	
                                                                   
                                                                   "Flavobacteriales_unclassified",	
                                                                   
                                                                   "Fusobacteriales_unclassified",	
                                                                   
                                                                   "Gastranaerophilales_fa",	
                                                                   "Gemmatimonadaceae",	
                                                                   "Geodermatophilaceae",	
                                                                   "Gitt-GS-136_fa",	
                                                                   "Gracilibacteria_fa",	
                                                                   "Granulosicoccaceae",	
                                                                   "Hahellaceae",	
                                                                   "Halanaerobiales_unclassified",	
                                                                   "Haliangiaceae",	
                                                                   "Halieaceae",	
                                                                   "Halomonadaceae",	
                                                                   
                                                                   "Herpetosiphonaceae",	
                                                                   "HOC36_fa",	
                                                                   "Holosporaceae",	
                                                                   "Hydrogenophilaceae",	
                                                                   "Hyphomicrobiaceae",	
                                                                   "Hyphomonadaceae",	
                                                                   "Iamiaceae",	
                                                                   "Intrasporangiaceae",	
                                                                   "JG30-KF-CM45_fa",	
                                                                   "JG34-KF-361",	
                                                                   "JTB255_marine_benthic_group",	
                                                                   "KD4-96_fa",	
                                                                   "KF-JG30-B3",	
                                                                   "Kineosporiaceae",	
                                                                   
                                                                   "Lactobacillaceae",	
                                                                   "Lactobacillales_unclassified",	
                                                                   "LD29",	
                                                                   "Legionellaceae",	
                                                                   "Lentisphaeraceae",	
                                                                   "Leptospiraceae",	
                                                                   "Leptotrichiaceae",	
                                                                   "Leuconostocaceae",	
                                                                   "Listeriaceae",	
                                                                   "Longimicrobiaceae",	
                                                                   
                                                                   "Methylophilaceae",	
                                                                   "MgMjR-022",	
                                                                   "Microbacteriaceae",	
                                                                   "Micrococcaceae",	
                                                                   "Micrococcales_unclassified",	
                                                                   "Microgenomates_unclassified",	
                                                                   "Milano-WF1B-44_fa",	
                                                                   "Mollicutes_unclassified",	
                                                                   "Moraxellaceae",	
                                                                   "MSB-1E8",	
                                                                   "Mycobacteriaceae",	
                                                                   
                                                                   "NB1-n_fa",	
                                                                   "Neisseriaceae",	
                                                                   "Nitriliruptoraceae",	
                                                                   "Nitrosomonadaceae",	
                                                                   "Nocardiaceae",	
                                                                   "Nocardioidaceae",	
                                                                   "NS9_marine_group",	
                                                                   "Obscuribacterales_fa",	
                                                                   "Oceanospirillaceae",	
                                                                   "Oceanospirillales_unclassified",	
                                                                   "Oligoflexaceae",	
                                                                   "OM1_clade",	
                                                                   "OM182_clade",	
                                                                   "OM190_fa",	
                                                                   
                                                                   "Paenibacillaceae",	
                                                                   "Parcubacteria_fa",	
                                                                   "Parvularculaceae",	
                                                                   "Pasteurellaceae",	
                                                                   "PeM15_fa",	
                                                                   
                                                                   "Peregrinibacteria_fa",	
                                                                   "Phycisphaeraceae",	
                                                                   "Phyllobacteriaceae",	
                                                                   "Pla3_lineage_fa",	
                                                                   "Planctomycetaceae",	
                                                                   "Planococcaceae",	
                                                                   
                                                                   "Prevotellaceae",	
                                                                   "Propionibacteriaceae",	
                                                                   "Proteobacteria_unclassified",	
                                                                   "Pseudoalteromonadaceae",	
                                                                   "Pseudomonadaceae",	
                                                                   "Pseudonocardiaceae",	
                                                                   "Psychromonadaceae",	
                                                                   "Puniceicoccaceae",	
                                                                   "PYR10d3_fa",	
                                                                   "Rhizobiaceae",	
                                                                   "Rhizobiales_Incertae_Sedis",	
                                                                   
                                                                   
                                                                   "Rhodobiaceae",	
                                                                   "Rhodocyclaceae",	
                                                                   "Rhodospirillaceae",	
                                                                   "Rhodospirillales_Incertae_Sedis",	
                                                                   "Rhodospirillales_unclassified",	
                                                                   "Rhodothermaceae",	
                                                                   "Rickettsiaceae",	
                                                                   "Rickettsiales_Incertae_Sedis",	
                                                                   "Rickettsiales_unclassified",	
                                                                   
                                                                   "Rs-D42",	
                                                                   "Rubritaleaceae",	
                                                                   
                                                                   "S0134_terrestrial_group_fa",	
                                                                   "S085_fa",	
                                                                   "Saccharibacteria_fa",	
                                                                   "Sandaracinaceae",	
                                                                   "Saprospiraceae",	
                                                                   "SAR116_clade",	
                                                                   "SAR324_clade(Marine_group_B)_fa",	
                                                                   "Shewanellaceae",	
                                                                   "Solibacteraceae_(Subgroup_3)",	
                                                                   "Solirubrobacteraceae",	
                                                                   "Sphingobacteriaceae",	
                                                                   "Sphingobacteriales_unclassified",	
                                                                   
                                                                   "Sphingomonadales_unclassified",	
                                                                   "Spirochaetaceae",	
                                                                   "Spirochaetales_unclassified",	
                                                                   "Spongiibacteraceae",	
                                                                   "SR-FBR-L83",	
                                                                   "SR1_(Absconditabacteria)_fa",	
                                                                   
                                                                   "Streptococcaceae",	
                                                                   "Streptomycetaceae",	
                                                                   "Subgroup_21_fa",	
                                                                   "Subgroup_22_fa",	
                                                                   "Subgroup_23_fa",	
                                                                   "Subgroup_6_fa",	
                                                                   "Subgroup_7_fa",	
                                                                   
                                                                   "Surface_1",	
                                                                   "Sva0071_fa",	
                                                                   "Sva0725",	
                                                                   "Sva0996_marine_group",	
                                                                   "Sva1033",	
                                                                   "Syntrophaceae",	
                                                                   "T9d",	
                                                                   "Tepidisphaeraceae",	
                                                                   "Thermoactinomycetaceae",	
                                                                   "Thiotrichaceae",	
                                                                   "TM146",	
                                                                   "TM6_(Dependentiae)_fa",	
                                                                   "TRA3-20_fa",	
                                                                   "Trueperaceae",	
                                                                   "uncultured",	
                                                                   "uncultured_fa",	
                                                                   "Unknown_Family",	
                                                                   "Veillonellaceae",	
                                                                   "Verrucomicrobia_fa",	
                                                                   "Verrucomicrobia_unclassified",	
                                                                   "Verrucomicrobiaceae",	
                                                                   "Verrucomicrobiales_unclassified",	
                                                                   "Vibrionaceae",	
                                                                   "X35_fa",	
                                                                   "Xanthobacteraceae",	
                                                                   "Xanthomonadaceae",	
                                                                   "Xanthomonadales_Incertae_Sedis",	
                                                                   "Xanthomonadales_unclassified"))

                                                                  
                                                                 
                                                                  
                                                                      

##only families that make up more than 5% total abundance are coloured. The rest are coloured white

colors <- c("navy",
            "olivedrab4","olivedrab2","forestgreen","darkgreen","yellow3","yellow1","yellowgreen","lawngreen","springgreen3","aquamarine",
            "mistyrose","pink","palevioletred1","palevioletred4", "orchid1","orchid4", "maroon1","magenta","deeppink3","red4","purple",
            
            "lightskyblue", "cyan","skyblue4", "skyblue3","darkslategray4","blue",
            "black",
            "red",
            "chocolate1","chocolate",
            "grey",
            "bisque",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white","white","white")


names(colors) <- levels(migration_family$Family)
colScale <- scale_fill_manual(name = "Family",values = colors)


migration_family$Sample<-factor(migration_family$Sample, level = c("8682", "8629", "8615", "8685", "8610", "8662", "8560", "8607", "8671", "8589", "8657", "8573", "8625", "8627", "8639", "8549", "8645", "8632", "8660", "8649", "8601",
                                                                   "8658", "8638", "8646", "8634", "8606", "8612", "8581", "8618", "8609", "8628", "8672", "8605", "8689", "8659", "8616", "8626", "8604", "8569", "8677", "8588", "8675",
                                                                   "8551", "8613", "8546", "8683", "8590", "8557", "8558", "8674", "8559", "8680", "8630", "8563", "8586", "8556", "8561", "8644", "8667", "8566", "8665", "8633", "8564",
                                                                   "8583", "8575", "8571", "8565", "8578", "8570", "8684", "8688", "8655", "8552", "8548", "8679", "8663", "8686"))




head(migration_family)
ggplot(migration_family, aes(x = Sample, y = Abundance, fill = Family)) + 
  
  geom_bar(stat = "identity", position="fill") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=12), axis.title=element_text(size=14))+
  ggtitle("Family abundance per sample")+facet_wrap(species~site~age, scales="free", nrow=3, ncol=2)+colScale + theme(legend.position="none")






################################################ ORDINATION ANLYSIS#############################################

set.seed(1)

migration_pcoa <- ordinate(
  physeq = migration, 
  method = "NMDS", 
  distance = "unifrac"
)


# Plot 
plot_ordination(
  physeq = migration,
  ordination = migration_pcoa) +
  
  geom_point(aes(colour = migration, shape = migration), alpha = 1, size = 3, stroke = 2) +
  theme_bw()+scale_color_manual(values=c("royalblue","turquoise1", 
                                         "forestgreen","lawngreen","violetred2","plum1"))+
  scale_shape_manual(values=c(17,2,15,0,16,1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16))

##stats
migration_unifrac <- phyloseq::distance(migration, method = "unifrac")


# make a data frame from the sample_data
sampledf <- data.frame(sample_data(migration))

# Adonis test
adonis(migration_unifrac ~ site+species+age, data = sampledf)

###############################################

##quantify variation##############################



#curlew sand#########################################

curlew_unifrac <- phyloseq::distance(curlew.b, method = "unifrac")

curlew_metadata<-data.frame(sample_data(curlew.b))
curlew_age<-as.factor(curlew_metadata$age)

curlew_betadisper<-betadisper(curlew_unifrac,curlew_age, type = c("median","centroid"), bias.adjust = FALSE,
                              sqrt.dist = FALSE, add = FALSE)

anova(curlew_betadisper)

####################stint broome

stint.b_unifrac <- phyloseq::distance(stint.b, method = "unifrac")
stint.b_metadata<-data.frame(sample_data(stint.b))
stint.b_age<-as.factor(stint.b_metadata$age)

stint.b_betadisper<-betadisper(stint.b_unifrac,stint.b_age, type = c("median","centroid"), bias.adjust = FALSE,
                               sqrt.dist = FALSE, add = FALSE)

anova(stint.b_betadisper)

############## stint flinders

stint.f_unifrac <- phyloseq::distance(stint.f, method = "unifrac")
#stint.f_bray <- phyloseq::distance(stint.f, method = "unifrac")
stint.f_metadata<-data.frame(sample_data(stint.f))
stint.f_age<-as.factor(stint.f_metadata$age)

stint.f_betadisper<-betadisper(stint.f_bray,stint.f_age)

anova(stint.f_betadisper)



################################bray curtis

set.seed(1)

migration_pcoa <- ordinate(
  physeq = migration, 
  method = "NMDS", 
  distance = "bray"
)



# Plot 
plot_ordination(
  physeq = migration,
  ordination = migration_pcoa) +
  
  geom_point(aes(colour = migration, shape = migration), alpha = 1, size = 3, stroke = 2) +
  theme_bw()+scale_color_manual(values=c("royalblue","turquoise1", 
                                         "forestgreen","lawngreen","violetred2","plum1"))+
  scale_shape_manual(values=c(17,2,15,0,16,1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16))


##stats

migration_bray <- phyloseq::distance(migration, method = "bray")


# make a data frame from the sample_data
sampledf <- data.frame(sample_data(migration))

# Adonis test
adonis(migration_bray ~ site+species+age, data = sampledf)

##quantify variation##############################

#overall

metadata<-data.frame(sample_data(migration))
age<-as.factor(metadata$age)
species<-as.factor(metadata$species)

betadisper<-betadisper(migration_bray,age, type = c("median","centroid"), bias.adjust = FALSE,
           sqrt.dist = FALSE, add = FALSE)

anova(betadisper)

betadisper<-betadisper(migration_bray,species, type = c("median","centroid"), bias.adjust = FALSE,
                       sqrt.dist = FALSE, add = FALSE)

anova(betadisper)

#curlew sand#########################################

curlew_bray <- phyloseq::distance(curlew.b, method = "bray")
#curlew_bray <- phyloseq::distance(curlew.b, method = "unifrac")
curlew_metadata<-data.frame(sample_data(curlew.b))
curlew_age<-as.factor(curlew_metadata$age)

curlew_betadisper<-betadisper(curlew_bray,curlew_age, type = c("median","centroid"), bias.adjust = FALSE,
                       sqrt.dist = FALSE, add = FALSE)

anova(curlew_betadisper)

####################stint broome

stint.b_bray <- phyloseq::distance(stint.b, method = "bray")
#stint.b_bray <- phyloseq::distance(stint.b, method = "unifrac")
stint.b_metadata<-data.frame(sample_data(stint.b))
stint.b_age<-as.factor(stint.b_metadata$age)

stint.b_betadisper<-betadisper(stint.b_bray,stint.b_age, type = c("median","centroid"), bias.adjust = FALSE,
                              sqrt.dist = FALSE, add = FALSE)

anova(stint.b_betadisper)

############## stint flinders

stint.f_bray <- phyloseq::distance(stint.f, method = "bray")
#stint.f_bray <- phyloseq::distance(stint.f, method = "unifrac")
stint.f_metadata<-data.frame(sample_data(stint.f))
stint.f_age<-as.factor(stint.f_metadata$age)

stint.f_betadisper<-betadisper(stint.f_bray,stint.f_age)

anova(stint.f_betadisper)


#########################################PCA

set.seed(1)

migration_pcoa <- ordinate(
  physeq = migration, 
  method = "MDS", 
  distance = "bray"
)



plot_ordination(
  physeq = migration,
  ordination = migration_pcoa) +
  
  geom_point(aes(colour = migration, shape = migration), alpha = 1, size = 3, stroke = 2) +
  theme_bw()+scale_color_manual(values=c("royalblue","turquoise1", 
                                         "forestgreen","lawngreen","violetred2","plum1"))+
  scale_shape_manual(values=c(17,2,15,0,16,1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16))

#####################################################

set.seed(1)

migration_pcoa <- ordinate(
  physeq = migration, 
  method = "MDS", 
  distance = "unifrac"
)


plot_ordination(
  physeq = migration,
  ordination = migration_pcoa) +
  
  geom_point(aes(colour = migration, shape = migration), alpha = 1, size = 3, stroke = 2) +
  theme_bw()+scale_color_manual(values=c("royalblue","turquoise1", 
                                         "forestgreen","lawngreen","violetred2","plum1"))+
  scale_shape_manual(values=c(17,2,15,0,16,1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16))

#############################################################


##remove corynebacteria to see if any other differences

migration1<-subset_taxa(migration, Order != "Corynebacteriales")
migration1<- prune_taxa(taxa_sums(migration1)>0, migration1)

set.seed(1)

migration_pcoa <- ordinate(
  physeq = migration1, 
  method = "NMDS", 
  distance = "bray"
)



# Plot 
plot_ordination(
  physeq = migration1,
  ordination = migration_pcoa) +
  
  geom_point(aes(colour = migration, shape = migration), alpha = 1, size = 3, stroke = 2) +
  theme_bw()+scale_color_manual(values=c("royalblue","turquoise1", 
                                         "forestgreen","lawngreen","violetred2","plum1"))+
  scale_shape_manual(values=c(17,2,15,0,16,1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16))

##stats

migration_bray <- phyloseq::distance(migration1, method = "bray")
migration_bray <- phyloseq::distance(migration1, method = "unifrac")


# make a data frame from the sample_data
sampledf <- data.frame(sample_data(migration1))

# Adonis test
adonis(migration_bray ~ site+species+age, data = sampledf)

##########################################################

migration_phylum <- migration1 %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

#filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
head(migration_phylum)


###############################################################

###########plot

migration_phylum$Phylum<-factor(migration_phylum$Phylum)
#broome_phylum$Order<-factor(broome_phylum$Order)
unique(migration_phylum$Phylum )
migration_phylum$Phylum<-factor(migration_phylum$Phylum)



migration_phylum$Sample<-factor(migration_phylum$Sample, levels = migration_phylum$Sample[order(migration_phylum$species, migration_phylum$site )])

unique(migration_phylum$Phylum)
migration_phylum$Phylum<-factor(migration_phylum$Phylum, level = c("Actinobacteria", 
                                                                   "Firmicutes",
                                                                   "Proteobacteria",
                                                                   "Bacteroidetes",
                                                                   "Deferribacteres",
                                                                   "Fusobacteria",
                                                                   "Tenericutes",        
                                                                   "Spirochaetae" ,
                                                                   "Deinococcus-Thermus",
                                                                   "Acidobacteria",
                                                                   "Armatimonadetes",
                                                                   "Bacteria_unclassified",
                                                                   "Chloroflexi",
                                                                   "Cyanobacteria",
                                                                   "FBP",
                                                                   "Fibrobacteres",
                                                                   "Gemmatimonadetes",
                                                                   "Gracilibacteria",
                                                                   "Ignavibacteriae",
                                                                   "Lentisphaerae",
                                                                   "Microgenomates",
                                                                   "Parcubacteria",
                                                                   "Peregrinibacteria",       
                                                                   "Planctomycetes",
                                                                   "Saccharibacteria",                   
                                                                   "SR1_(Absconditabacteria)",
                                                                   "TM6_(Dependentiae)",
                                                                   "Verrucomicrobia")) 





unique(migration_phylum$Sample)
colors <- c("navyblue","green4","orchid","royalblue","black","red","brown","grey","bisque",
            "white","white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white","white")


names(colors) <- levels(migration_phylum$Phylum)
colScale <- scale_fill_manual(name = "Phylum",values = colors)

migration_phylum$Sample<-factor(migration_phylum$Sample, level = c("8682", "8629", "8615", "8685", "8610", "8662", "8560", "8607", "8671", "8589", "8657", "8573", "8625", "8627", "8639", "8549", "8645", "8632", "8660", "8649", "8601",
                                                                   "8658", "8638", "8646", "8634", "8606", "8612", "8581", "8618", "8609", "8628", "8672", "8605", "8689", "8659", "8616", "8626", "8604", "8569", "8677", "8588", "8675",
                                                                   "8551", "8613", "8546", "8683", "8590", "8557", "8558", "8674", "8559", "8680", "8630", "8563", "8586", "8556", "8561", "8644", "8667", "8566", "8665", "8633", "8564",
                                                                   "8583", "8575", "8571", "8565", "8578", "8570", "8684", "8688", "8655", "8552", "8548", "8679", "8663", "8686"))


ggplot(migration_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=12), axis.title=element_text(size=14))+
  facet_wrap(species~site~age, scales="free", nrow=3, ncol=2)+colScale

############################################################################################################

########################################################################################################

####alpha diversity

##use unscaled data


plot_richness(migration, color = "migration")



min_lib <- min(sample_sums(migration))

# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(migration)
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(migration)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(migration)


###

set.seed(3)

for (i in 1:100) {
  # Subsample
  r <-  rarefy_even_depth(migration, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed"))) 
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
}

SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

s <- data.frame(sample_data(migration))

alphadiv_rich_observed <- merge(rich_stats, s, by = "row.names") 
alphadiv_even_inverse <- merge(even_stats, s, by = "row.names") 

#summary stats

tapply(alphadiv_rich_observed$mean, alphadiv_rich_observed$migration, summary)

unique(alphadiv_rich_observed$migration)
group1a<-subset(alphadiv_rich_observed, migration=="Broome-curlewsand-adult")
sd(group1a$mean)
group1b<-subset(alphadiv_rich_observed, migration=="Broome-curlewsand-SY")
sd(group1b$mean)
group2a<-subset(alphadiv_rich_observed, migration=="Broome-stint-adult")
sd(group2a$mean)
group2b<-subset(alphadiv_rich_observed, migration=="Broome-stint-SY")
sd(group2b$mean)
group3a<-subset(alphadiv_rich_observed, migration=="Flinders-stint-adult")
sd(group3a$mean)
group3b<-subset(alphadiv_rich_observed, migration=="Flinders-stint-SY")
sd(group3b$mean)



ggplot(alphadiv_rich_observed, aes(x = migration, y = mean))+geom_boxplot(width=0.5)+
  geom_jitter(width=0.1, alpha = 0.3, size=3)+ylab("Shannon index")+
  theme(axis.text=element_text(size=12, angle=45, hjust=1), axis.title=element_text(size=14,face="bold"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##test diversity differences between migration status

##RNS Broome

alphadiv_RNS_Broome<-subset(alphadiv_rich_observed, site=="Broome" & species=="Red-necked stint")

t.test(alphadiv_RNS_Broome$mean~ alphadiv_RNS_Broome$age)

#RNS flinders

alphadiv_RNS_flinders<-subset(alphadiv_rich_observed, site=="Flinders" & species=="Red-necked stint")

t.test(alphadiv_RNS_flinders$mean~ alphadiv_RNS_flinders$age)

#CS Broome

alphadiv_CS_Broome<-subset(alphadiv_rich_observed, site=="Broome"& species=="Curlew sandpiper")

t.test(alphadiv_CS_Broome$mean~ alphadiv_CS_Broome$age)

########################################################################


#############################FUNCTION/PICRUST 97 % similarity

##level2 collapsed functions

function_all <- read.csv("C:/Users/arisely/Dropbox/PhD/Microbiome/CHAPTER 3 - PHYSIOLOGY/ANALYSIS/PiCrust/Collapsed function paths/function_all_relabundance.csv")

picrust<-t(function_all)

colnames(picrust) = picrust[1, ] # the first row will be the header
picrust = picrust[-1, ] 

picrust<-data.frame(picrust)



picrust <- setDT(picrust, keep.rownames = TRUE)[]

colnames(picrust)[1]<-"Sample"

picrust1<- gather(picrust, functions, abundance, Amino.Acid.Metabolism:Xenobiotics.Biodegradation.and.Metabolism, factor_key=TRUE)

str(picrust1)
picrust1$Sample<-as.factor(picrust1$Sample)
picrust1$abundance<-as.numeric(picrust1$abundance)

unique(picrust1$functions)

picrust1$functions1[picrust1$functions=="Amino.Acid.Metabolism"]<-"Metabolism - amino acids"
picrust1$functions1[picrust1$functions=="Biosynthesis.of.Other.Secondary.Metabolites"]<-"Metabolism - secondary metabolites"
picrust1$functions1[picrust1$functions=="Cancers"]<-"Human diseases - cancer"
picrust1$functions1[picrust1$functions=="Carbohydrate.Metabolism"]<-"Metabolism - carbohydrates"
picrust1$functions1[picrust1$functions=="Cardiovascular.Diseases"]<-"Human diseases - cardiovascular"
picrust1$functions1[picrust1$functions=="Cell.Communication"]<-"Cell signaling"
picrust1$functions1[picrust1$functions=="Cell.Growth.and.Death"]<-"Cellular processes"
picrust1$functions1[picrust1$functions=="Cell.Motility"]<-"Cellular processes"
picrust1$functions1[picrust1$functions=="Cellular.Processes.and.Signaling"]<-"Cellular processes"
picrust1$functions1[picrust1$functions=="Circulatory.System"]<-"Organismal systems - circulatory"
picrust1$functions1[picrust1$functions=="Digestive.System"]<-"Organismal systems - digestive"
picrust1$functions1[picrust1$functions=="Endocrine.System"]<-"Organismal systems - endocrine"
picrust1$functions1[picrust1$functions=="Energy.Metabolism"]<- "Metabolism - energy"
picrust1$functions1[picrust1$functions=="Environmental.Adaptation"]<-"Organismal systems - environmental adaptation"
picrust1$functions1[picrust1$functions=="Enzyme.Families"]<-"Enzyme families"
picrust1$functions1[picrust1$functions=="Excretory.System"]<-"Organismal systems - excretory"
picrust1$functions1[picrust1$functions=="Folding..Sorting.and.Degradation"]<-"Genetic processes"
picrust1$functions1[picrust1$functions=="Genetic.Information.Processing"]<-"Genetic processes"
picrust1$functions1[picrust1$functions=="Glycan.Biosynthesis.and.Metabolism"]<-"Metabolism - glycan"
picrust1$functions1[picrust1$functions=="Immune.System"]<-"Organismal systems - immune function"
picrust1$functions1[picrust1$functions=="Immune.System.Diseases"]<-"Human diseases - immune system"
picrust1$functions1[picrust1$functions=="Infectious.Diseases"]<-"Human diseases - infectious"
picrust1$functions1[picrust1$functions=="Lipid.Metabolism"]<-"Metabolism - lipids"
picrust1$functions1[picrust1$functions=="Membrane.Transport"]<-"Cell signaling"
picrust1$functions1[picrust1$functions=="Metabolic.Diseases"]<-"Human diseases - metabolic"
picrust1$functions1[picrust1$functions=="Metabolism"]<- "Metabolism - general"
picrust1$functions1[picrust1$functions=="Metabolism.of.Cofactors.and.Vitamins"]<-"Metablism - vitamins"
picrust1$functions1[picrust1$functions=="Metabolism.of.Other.Amino.Acids"]<-"Metabolism - amino acids"
picrust1$functions1[picrust1$functions=="Metabolism.of.Terpenoids.and.Polyketides"]<-"Metabolism - terpenoids"
picrust1$functions1[picrust1$functions=="Nervous.System"]<-"Organismal systems - nervous"
picrust1$functions1[picrust1$functions=="Neurodegenerative.Diseases"]<-"Human diseases - neurodegerative"
picrust1$functions1[picrust1$functions=="Nucleotide.Metabolism"]<-"Metabolism - nucleotides"
picrust1$functions1[picrust1$functions=="Poorly.Characterized"]<-"Poorly characterized"
picrust1$functions1[picrust1$functions=="Replication.and.Repair"]<-"Genetic processes"
picrust1$functions1[picrust1$functions=="Sensory.System"]<-"Organismal systems - sensory"
picrust1$functions1[picrust1$functions=="Signal.Transduction"]<-"Cell signaling"
picrust1$functions1[picrust1$functions=="Signaling.Molecules.and.Interaction"]<-"Cell signaling"
picrust1$functions1[picrust1$functions=="Transcription"]<-"Genetic processes"
picrust1$functions1[picrust1$functions=="Translation"]<-"Genetic processes"
picrust1$functions1[picrust1$functions=="Transport.and.Catabolism"]<-"Cellular processes"
picrust1$functions1[picrust1$functions=="Xenobiotics.Biodegradation.and.Metabolism"]<- "Metabolism - xenobiotics"

unique(picrust1$functions1)

picrust1$functions1<-factor(picrust1$functions1, level = c("Metablism - vitamins"                         ,
                                                           "Metabolism - amino acids"                     ,
                                                           "Metabolism - carbohydrates"                   ,
                                                           "Metabolism - energy"                          ,
                                                           "Metabolism - general"                         ,
                                                           "Metabolism - glycan"                          ,
                                                           "Metabolism - lipids"                          ,
                                                           "Metabolism - nucleotides"                     ,
                                                           "Metabolism - secondary metabolites"           ,
                                                           "Metabolism - terpenoids"                      ,
                                                           "Metabolism - xenobiotics"                     ,
                                                           "Cell signaling"                               ,
                                                           "Cellular processes"                           ,
                                                           "Genetic processes"                            ,
                                                           "Enzyme families"                              ,
                                                           "Human diseases - cancer"                      ,
                                                           "Human diseases - cardiovascular"              ,
                                                           "Human diseases - immune system"               ,
                                                           "Human diseases - infectious"                  ,
                                                           "Human diseases - metabolic"                   ,
                                                           "Human diseases - neurodegerative"             ,
                                                           "Organismal systems - circulatory"             ,
                                                           "Organismal systems - digestive"               ,
                                                           "Organismal systems - endocrine"               ,
                                                           "Organismal systems - environmental adaptation",
                                                           "Organismal systems - excretory"               ,
                                                           "Organismal systems - immune function"         ,
                                                           "Organismal systems - nervous"                 ,
                                                           "Organismal systems - sensory"                 ,
                                                           "Poorly characterized"))

                                                           
colors <- c(
            "mistyrose","pink","palevioletred1","palevioletred4", "orchid1","orchid4", "deeppink","magenta","deeppink3","red4","purple",
            "forestgreen","yellowgreen",
            "royalblue",
            "bisque",
            "red", "red","red", "red","red","red",
            "yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow",
            "grey")
            
            
names(colors) <- levels(picrust1$function1)
colScale <- scale_fill_manual(name = "functions1",values = colors)

picrust1$Sample<-factor(picrust1$Sample, levels = picrust1$Sample[order(picrust1$Species, picrust1$Site, picrust1$Age)])



ggplot(picrust1, aes(x = Sample, y = abundance, fill = functions1)) + 
  
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=12), axis.title=element_text(size=14))+
 colScale + facet_wrap(Species~Site~Age, scales="free", nrow=3, ncol=2)

##########################

picrust<- read.csv("picrust_unrarefied.csv")
str(picrust)

picrust1<- gather(picrust, functions, abundance, Trichloro.ethane.degradation:p53.signaling.pathway, factor_key=TRUE)
picrust1$Sample<-as.factor(picrust1$Sample)
picrust1$abundance<-as.numeric(picrust1$abundance)

picrust2<-subset(picrust1, group=="two")

insulin<-picrust2%>%filter(str_detect(functions,'Insulin')==T)

ggplot(insulin, aes(x = age, y = abundance))+geom_boxplot()+geom_point()
t.test(insulin$abundance~ insulin$age)

picrust3 <- ddply(picrust2, c("functions", "age"), summarise, mean.abund=mean(abundance))

picrust3 <- spread(picrust3, age, mean.abund)

ggplot(picrust2, aes(x = functions, y = abundance, fill=age))+geom_boxplot()

picrust4<-picrust2%>%filter(str_detect(functions,'etabolism')==T)
ggplot(picrust4, aes(x = functions, y = abundance, fill=age))+geom_boxplot()+theme(axis.text.x = element_text(angle=60, hjust=1))+ylim(0,500000)

picrust5<-picrust2%>%filter(str_detect(functions,'synthesis')==T)
ggplot(picrust5, aes(x = functions, y = abundance, fill=age))+geom_boxplot()+theme(axis.text.x = element_text(angle=60, hjust=1))

picrust6<-picrust2%>%filter(str_detect(functions,'Fat')==T)
ggplot(picrust6, aes(x = functions, y = abundance, fill=age))+geom_boxplot()+theme(axis.text.x = element_text(angle=60, hjust=1))

ggplot(picrust, aes(x = age, y = stick.rate))+geom_boxplot()

##################################################################################

##rarefied to 1000

picrust<- read.csv("picrust_rarefied1000.csv")
picrust1<- gather(picrust, functions, abundance, Trichloro.ethane.degradation:p53.signaling.pathway, factor_key=TRUE)
picrust1$Sample<-as.factor(picrust1$Sample)
picrust1$abundance<-as.numeric(picrust1$abundance)

picrust2<-subset(picrust1, group=="two")

insulin<-picrust2%>%filter(str_detect(functions,'Insulin')==T)

ggplot(insulin, aes(x = age, y = abundance))+geom_boxplot()+geom_point()
t.test(insulin$abundance~ insulin$age)

picrust3 <- ddply(picrust2, c("functions", "age"), summarise, mean.abund=mean(abundance))

picrust3 <- spread(picrust3, age, mean.abund)

ggplot(picrust2, aes(x = functions, y = abundance, fill=age))+geom_boxplot()

picrust4<-picrust2%>%filter(str_detect(functions,'etabolism')==T)
ggplot(picrust4, aes(x = functions, y = abundance, fill=age))+geom_boxplot()+theme(axis.text.x = element_text(angle=60, hjust=1))+ylim(0,75000)

picrust5<-picrust2%>%filter(str_detect(functions,'synthesis')==T)
ggplot(picrust5, aes(x = functions, y = abundance, fill=age))+geom_boxplot()+theme(axis.text.x = element_text(angle=60, hjust=1))

picrust6<-picrust2%>%filter(str_detect(functions,'Fat')==T)
ggplot(picrust6, aes(x = functions, y = abundance, fill=age))+geom_boxplot()+theme(axis.text.x = element_text(angle=60, hjust=1))

ggplot(picrust, aes(x = age, y = stick.rate))+geom_boxplot()



#########################################################################

##what OTUs are even included in this analysis? Does it capture Corynebacterium?

#################PICRUST 97% #########ONLY INCLUDE OTUS WITH ASSIGNED REFERENCE OTUS


# Assign variables for imported data

#for 97% similarity of closed ref OTUS

sharedfile = "function97.shared"
taxfile = "function97.taxonomy"

#for 95% similarity

sharedfile = "function95.shared"
taxfile = "function95.taxonomy"

##import mothur shared and tax files

mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

##import metadata
map <- read.csv("shorebird.metadata.csv", header=T, row.names=1)

##make meta data into phyloseq format

map <- sample_data(map)
str(map)

head(map)
tail(map)



##merge the metadata into the phyloseq object

moth_merge <- merge_phyloseq(mothur_data, map)

colnames(tax_table(moth_merge))

colnames(tax_table(moth_merge)) <- c("Kingdom", "Phylum", "Class", 
                                     "Order", "Family", "Genus")

data_all <- moth_merge 
data_all

tax_table(moth_merge)


################################################################look at read library


sample_sum_df <- data.frame(sum = sample_sums(data_all))
plot(sample_sum_df$sum)


ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  ylab("Frequency")+
  theme(axis.title.y = element_blank())

##Broome samples didn't work :()

##############################################################################################################################

##get rid of 9185 as lowest read count
##get rid of 8688 as mostly unclassified taxa... not sure why

data_all<-data_all %>%
  subset_samples(Group != "Env_Broome_A")

data_all<-data_all %>%
  subset_samples(Group != "Env_Broome_E")


data_all<-prune_taxa(taxa_sums(data_all)>0, data_all)


##get rid of negative control altogether


data_all<-data_all%>%subset_samples(species!="Neg")

#delete replicates

data_all<-data_all%>%subset_samples(remove_replicate=="No")
sample_data(data_all)

##just birds

birds<-data_all%>%subset_samples(type=="bird")
birds <- prune_taxa(taxa_sums(birds) > 0, birds)

smin <- min(sample_sums(birds))




##for migration analysis
migration<-birds%>%subset_samples(species=="Red-necked stint"|species=="Curlew sandpiper")
migration<-migration%>%subset_samples(month=="August"|month=="September")
migration <- prune_taxa(taxa_sums(migration) > 0, migration)
migration



#############stacked barplots

migration_phylum <- migration %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

#filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
head(migration_phylum)



###############################################################

###########plot

migration_phylum$Phylum<-factor(migration_phylum$Phylum)
#broome_phylum$Order<-factor(broome_phylum$Order)
unique(migration_phylum$Phylum )
migration_phylum$Phylum<-factor(migration_phylum$Phylum)



migration_phylum$Sample<-factor(migration_phylum$Sample, levels = migration_phylum$Sample[order(migration_phylum$species, migration_phylum$site )])

unique(migration_phylum$Phylum)
migration_phylum$Phylum<-factor(migration_phylum$Phylum, level = c("p__Actinobacteria", 
                                                                   "p__Firmicutes", 
                                                                   "p__Proteobacteria",
                                                                   "p__Bacteroidetes", 
                                                                   "p__Fusobacteria", 
                                                                   "p__Tenericutes", 
                                                                   "p__Spirochaetes",
                                                                   
                                                                    "p__[Thermi]",        
                                                                   "p__Acidobacteria",    
                                                                   "p__Armatimonadetes",  
                                                                   "p__Chlorobi",         
                                                                   "p__Chloroflexi",      
                                                                   "p__Cyanobacteria",    
                                                                   "p__FBP",              
                                                                   "p__Gemmatimonadetes", 
                                                                   "p__Lentisphaerae",    
                                                                   "p__Nitrospirae",      
                                                                   "p__OD1",              
                                                                   "p__Planctomycetes",   
                                                                   "p__SAR406",           
                                                                   "p__SR1",              
                                                                   "p__Synergistetes",   
                                                                   "p__TM6",              
                                                                   "p__TM7",              
                                                                   "p__Verrucomicrobia")) 





unique(migration_phylum$Sample)
colors <- c("navyblue","green4","orchid","royalblue","black","red","brown",
            "white","white","white","white","white","white","white","white","white","white","white",
            "white","white","white","white","white","white","white")


names(colors) <- levels(migration_phylum$Phylum)
colScale <- scale_fill_manual(name = "Phylum",values = colors)

migration_phylum$Sample<-factor(migration_phylum$Sample, level = c("8682", "8629", "8615", "8685", "8610", "8662", "8560", "8607", "8671", "8589", "8657", "8573", "8625", "8627", "8639", "8549", "8645", "8632", "8660", "8649", "8601",
                                                                   "8658", "8638", "8646", "8634", "8606", "8612", "8581", "8618", "8609", "8628", "8672", "8605", "8689", "8659", "8616", "8626", "8604", "8569", "8677", "8588", "8675",
                                                                   "8551", "8613", "8546", "8683", "8590", "8557", "8558", "8674", "8559", "8680", "8630", "8563", "8586", "8556", "8561", "8644", "8667", "8566", "8665", "8633", "8564",
                                                                   "8583", "8575", "8571", "8565", "8578", "8570", "8684", "8688", "8655", "8552", "8548", "8679", "8663", "8686"))


ggplot(migration_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=12), axis.title=element_text(size=14))+
  facet_wrap(species~site~age, scales="free", nrow=3, ncol=2)+colScale
