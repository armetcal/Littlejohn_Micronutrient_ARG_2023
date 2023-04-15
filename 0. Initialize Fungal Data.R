library(phyloseq)
library(tidyverse)
library(here)

metadata = read.csv('Data Files/metadata.csv') %>% `rownames<-`(.$Sample)

species <- read.delim("Data Files/all_taxonomic_abundances.tsv",check.names=F) %>% as.data.frame() %>% 
  separate(`Consensus Lineage`,sep = '; ', 
           into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),remove = F) %>% 
  select(-`#OTU ID`) %>% 
  mutate(`Consensus Lineage` = `Consensus Lineage` %>% str_replace_all('; ','.')) %>% 
  column_to_rownames('Consensus Lineage')

names(species)[1:40] = sapply(names(species)[1:40], function(x) {
  w = which(metadata$MDB_ID == str_sub(x,end=10))
  return(metadata$Sample[w])
})

ps = phyloseq(sample_data(metadata),
                    otu_table(species %>% select(contains('LM'),contains("CON")) %>% as.matrix,taxa_are_rows = T),
                    tax_table(species %>% select(-contains('LM'),-contains("CON")) %>% as.matrix))

# Fungi only
ps = ps %>% subset_taxa(Phylum %in% c('p__Ascomycota','p__Basidiomycota','p__Microsporidia'))
ps.rel = ps %>% microbiome::transform('compositional')

# Format taxonomy table
Dttable = tax_table(ps)
Dttable = sub(Dttable, pattern = "^[urkpcofgs]__", replacement = "")

Dttable[Dttable == ""]  <-  NA

Dttable = as.data.frame(Dttable@.Data)

Dttable$Species = ifelse(!is.na(Dttable$Genus) & !is.na(Dttable$Species),
                          paste(as.character(Dttable$Genus), as.character(Dttable$Species)),
                          Dttable$Species)
tax_table(ps)  <- as.matrix(Dttable)

rm(Dttable,metadata,species)
