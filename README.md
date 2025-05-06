# ERV-Filtering-Code

Making new files with EVEs from South American primate families:
```R
setwd("/Users/zoegeorge/masterscode/eveDigest/")

# Getting interesting clusters for South American primate families

new_world_fams <- "Callitrichidae|Cebidae|Aotidae|Pitheciidae|Atelidae"
interesting_eves <- all_eves[which(grepl(new_world_fams, all_eves$host_family)),]

# Finding the clusters of those EVEs
my_clus <- unique(interesting_eves$cluster.representative)
all_clus_eves <- all_eves[which(all_eves$cluster.representative %in% my_clus),]
my_clus_info <- final[which(final$cluster.representative %in% my_clus),]

write.csv(all_clus_eves, "/home/biology/biol0273/Projects/eveDigest/mining_results/all_primates/zoe_cluster_candidates/data_frames/target_eves.csv")
write.csv(my_clus_info, "/home/biology/biol0273/Projects/eveDigest/mining_results/all_primates/zoe_cluster_candidates/data_frames/target_clusters_info.csv")

all_clus_eves <- read.csv("data/target_eves.csv", sep = ",", header = T, row.names = NULL)
my_clus_info <- read.csv("data/target_clusters_info.csv", sep = ",", header = T, row.names = NULL)
```

Removing clusters with labels that look like host proteins:
```R
reg_exp <- "akt|p120|fusion"
my_clus_info <- my_clus_info[which(!grepl(reg_exp, my_clus_info$rvdb_label_common)),]

reg_exp <- "actin|zinc|tyrosine-protein"
my_clus_info <- my_clus_info[which(!grepl(reg_exp, my_clus_info$nr_label_common)),]
```

Extract clusters with gags:
```R
gag_clusters <- my_clus_info[which(grepl("gag", my_clus_info$rvdb_label_common, ignore.case = T)),]
```

Filter for 0 framshifts, stop codons or introns:
```R
intact <- eves_target_cluster[which(eves_target_cluster$stop_codons==0 & 
                                      eves_target_cluster$frameshifts == 0 &
                                      eves_target_cluster$introns==0),]
```

Remaining filtering steps: 
```R
# Filter for between 9 and 100 intact sequences and a peptide length of over 200
gag_clusters <- my_clus_info[which(grepl("gag",
                                         my_clus_info$rvdb_label_common, ignore.case = T) 
                                          & my_clus_info$n.intact > 0
                                   &my_clus_info$mean.peptide.length > 200
                                   & my_clus_info$n.intact > 9 & 
                                     my_clus_info$n.intact < 100),] # More than three intact sequences
```

Looking at clusters of interest:
```R
# Identifying clusters of interest
eves_in_target_cluster <- all_clus_eves[which(all_clus_eves$cluster.representative=="CAUYPB010090649.1:69547-71010|Alouatta"),]
table(eves_in_target_cluster$host_family)

# Looking at values in target clusters

# Frequency of each unique value in 'introns' column in new world monkeys
introns_nw <- table(clus_nw_fams$introns)
print(introns_nw)

# Frequency of each unique value in 'introns' column in other monkeys
introns_other <- table(clus_other_fams$introns)
print(introns_other)

# Frequency of each unique value in 'frameshifts' column in new world monkeys
frameshifts_nw <- table(clus_nw_fams$frameshifts)
print(frameshifts_nw)

# Frequency of each unique value in 'frameshifts' column in other monkeys
frameshifts_other <- table(clus_other_fams$frameshifts)
print(frameshifts_other)

# Frequency of each unique value in 'stop_codons' column in new world monkeys
stop_codons_nw <- table(clus_nw_fams$stop_codons)
print(stop_codons_nw)

# Frequency of each unique value in 'stop_codons' column in other monkeys
stop_codons_other <- table(clus_other_fams$stop_codons)
print(stop_codons_other)


dim(clus)

table(clus$host_species)
num_species <- length(unique(clus$host_species))
num_species

#Finding the number of host genus in cluster
table(clus$host_genus)
num_genus <- length(unique(clus$host_genus))
num_genus

genus_counts <- table(clus$host_genus)

# Create a bar plot
par(mfrow=c(2,3))
barplot(genus_counts, 
        main = "Frequency of Host genus",
        ylab = "Frequency",
        las = 2,          # Makes species names on the x-axis vertical for readability
        col = "skyblue",  # Sets color of the bars
        border = "black") # Sets color of the bar borders

# Inspecting host families
table(clus$host_family)
num_family <- length(unique(clus$host_family))
num_family

family_counts <- table(clus$host_family)


# Create a bar plot
barplot(family_counts, 
        main = "Frequency of Host family",
        ylab = "Frequency",
        las = 2,          # Makes species names on the x-axis vertical for readability
        col = "skyblue",  # Sets color of the bars
        border = "black") # Sets color of the bar borders

# creating a histogram to look at peptide lengths
hist(clus$peptide_len)

# Making a frequency plot to look at the rvdb labels

rvdb_label <- table(clus$prot_label_rvdb)

# Create a bar plot
barplot(rvdb_label, 
        main = "Frequency of rvdb label",
        ylab = "Frequency",
        las = 2,          # Makes species names on the x-axis vertical for readability
        col = "skyblue",  # Sets color of the bars
        border = "black") # Sets color of the bar borders

#Making a frequcny plot to look at the nr labels 

nr_label <- table(clus$prot_label_nr)

# Create a bar plot
barplot(nr_label, 
        main = "Frequency of nr label",
        ylab = "Frequency",
        las = 2,          # Makes species names on the x-axis vertical for readability
        col = "skyblue",  # Sets color of the bars
        border = "black") # Sets color of the bar borders

```
 
