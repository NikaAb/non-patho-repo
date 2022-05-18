# Comparer 6 répertoires de IGH Homo sapiens avec 3 types des données: COVID +++ (S24,S26), COVID+ (M5, M6) et COVID - (H3, H4).
# Les datas ont été téléchargés dans le ireceptor-public-archives(chercher des information sur les données).
# Les séquences ont été analysés par IMGT/HighV-quest et la sortie a été récuperé en format AIRR (vquest_airr.tsv).

library(dplyr)
library(immuneREF)
library(stringr)
library(tidyr)




getwd()
rep <- c("H3")
#, "H4", "M5", "M6", "S24", "S26")
# Adresse du dossier qui contienne le fichier AIRR pour le répertoire H3
setwd("/Users/equipe.lefranc/Documents/Lorena Project/COVID/DONES")
list.files()
for (i in 1:length(rep)){
  setwd("C:/Users/equipe.lefranc/Documents")
        a <- "C:/Users/equipe.lefranc/Documents/Lorena Project/COVID/DONES/"
        b <-rep[i]
        rep_directory <- paste(a,b, sep = "")
    # On va stocker le fichier AIRR dans un tableau et on remplace les données vides
    #par NA et on va separer les données par tabulation
      file_name <- paste0(rep_directory,"/vquest_airr.tsv")

      
local_var <- read.table(file = file_name , header = TRUE, fill = TRUE,
                        sep = "\t", na.strings=c(""," ","NA"))


#head(_donnees_brut)

#names(_donnees_brut)

#remplacer le contennu de la colonne séquence des acides aminés par
#séquence_aligment_aa et aussi le contennu de la séquence par 
# sequence_alignment (sans gaps) car le format demandé par immuneREF est:
#“sequence_aa”: Full amino acid VDJ sequence
#“sequence”: Full nucleotide VDJ sequence
local_var$sequence_aa <- gsub("\\.","", 
                              as.character(local_var$sequence_alignment_aa))
local_var$sequence <- gsub("\\.","", 
                           as.character(local_var$sequence_alignment))

# On recupere les colonnes nécessaires pour lancer immuneREF
local_var <- select(local_var, sequence_id, sequence_aa, sequence, 
                    junction_aa, junction, v_call, d_call, j_call)
 


#Ajouter la colonne de fréquence et regroupe les séquences identiques

local_var <- local_var %>% 
  
  group_by(sequence_aa, sequence, junction_aa, junction, v_call, d_call, 
           j_call) %>%
  
  summarise(
    
    n = n(),
    
  ) %>%
  filter(
    
    n >= 0
  )


# Calculer la fréquence et vérifier si la sume des fréquences est égale à 1
freqs <- local_var$n/sum(local_var$n)
sum(freqs)


#Creer un data frame
local_var <- data.frame(local_var)


#Mis en format les noms des genes de V_call. J'enleve les nom de spece er ce qui
#est apres l'étoile (allele) [HACERLO POR PARTES Y UNIRLO AL FINAL!!!!!!!!!!]

dfv <- data.frame(Name = c(local_var$v_call)) 

dfv[c('specie', 'v_call')] <- str_split_fixed(dfv$Name," ", 2)

dfv[c('v_call', 'allele v_call')] <- str_split_fixed(dfv$v_call,"\\*", 2)

dfv <- dfv[c('v_call', 'allele v_call')]

dfv$v_call <- NULL
dfv$Name <- NULL
dfv$v_call <- NULL


#Separer les données de D

dfd <- data.frame(Name = c(local_var$d_call)) 

dfd[c('specie', 'd_call')] <- str_split_fixed(dfd$Name," ", 2)

dfd[c('d_call', 'allele d_call')] <- str_split_fixed(dfd$d_call,"\\*", 2)

dfd <- dfd[c('d_call', 'allele d_call')]


dfd$d_call <- NULL
dfd$Name <- NULL
dfd$d_call <- NULL
dfd$`allele d_call` <- NULL


#Separer les données de J

dfj <- data.frame(Name = c(local_var$j_call)) 

dfj[c('specie', 'j_call')] <- str_split_fixed(dfj$Name," ", 2)

dfj[c('j_call', 'allele j_call')] <- str_split_fixed(dfj$j_call,"\\*", 2)

dfj <- dfj[c('j_call', 'allele j_call')]

dfj$j_call <- NULL
dfj$Name <- NULL
dfj$j_call <- NULL
dfj$`allele j_call` <- NULL


repertoire <- cbind(local_var, dfv, dfd,  dfj, freqs)
repertoire$n <- NULL


#Convertir le format des colonnes de charactere à factor du répertoire 

repertoire <- data.frame(repertoire)
repertoire$sequence_aa=as.factor(repertoire$sequence_aa)
repertoire$sequence=as.factor(repertoire$sequence)                                 
repertoire$junction_aa=as.factor(repertoire$junction_aa)                                  
repertoire$junction=as.factor(repertoire$junction)
repertoire$v_call=as.factor(repertoire$v_call)
repertoire$d_call=as.factor(repertoire$d_call)
repertoire$j_call=as.factor(repertoire$j_call)
dfv$`allele v_call` <- NULL

#Le format de nom de répertoire est le suivant:
#nomdespece_locus_chaine_nomrépertoire
hs_ig_h_ <- repertoire

# Si il y a un probleme avec compatibilitycheck, on peut Unifier les genes pour
#voir si ils sont dans le germline list c'est pas obligatoire de le faire.
hs_ig_h_v <-hs_ig_h_ %>% group_by(v_call) %>% filter (! duplicated(v_call))
hs_ig_h_j <-hs_ig_h_ %>% group_by(j_call) %>% filter (! duplicated(j_call))

list_germline_genes$hs$ig$h$V$gene
hs_ig_h_v$v_call
hs_ig_h_v[65,]$sequence
list_germline_genes$hs$ig$h$V$sequence
hs_ig_h_j$j_call

setdiff(list_germline_genes$hs$ig$h$V$gene, hs_ig_h_v$v_call)
setdiff(list_germline_genes$hs$ig$h$J$gene, hs_ig_h_j$j_call)


#Si il manque des genes dans le germline_list, on peut les ajouter
df_species_VDJ<-list()

df_species_VDJ[["V1"]] <- data.frame(gene=c("IGHV2-70D"),allele="14",sequence=c("caggtcaccttgaaggagtctggtcctgcgctggtgaaacccacacagaccctcacactgacctgcaccttctctgggttctcactcagcactagtggaatgcgtgtgagctggatccgtcagcccccaggtaaggccctggagtggcttgcacgcattgattgggatgatgataaattctacagcacatctctgaagaccaggctcaccatctccaaggacacctccaaaaaccaggtggtccttacaatgaccaacatggaccctgtggacacagccacgtattactgtgcacggatac"),species="hs",frequency_uniform=c(0.01), frequency=c(4.147725e-06))
df_species_VDJ[["V2"]] <- data.frame(gene=c("IGHV3-41"),allele="02",sequence=c("gaggtgcagctggtggagtctgggggaggcttggtccagcctggggggtccctgagactctcctgtgcagcctcaggattctcctttagtagctatggcatgagctgggtccgccaggctccagggaaggggctggactgagtggcacatatctggaatgatggaagtcagaaatactatgcagactctgtgaagggccgattcacaatctccagagacaattctaagagcatgctctatctgcaaatggacagtctgaaagctaaggacacggccatgtattactgtaccaga"),species="hs",frequency_uniform=c(0.01), frequency=c(4.147725e-06))
df_species_VDJ[["V3"]] <- data.frame(gene=c("IGHV3-43D"),allele="04",sequence=c("gaagtgcagctggtggagtctgggggagtcgtggtacagcctggggggtccctgagactctcctgtgcagcctctggattcacctttgatgattatgccatgcactgggtccgtcaagctccggggaagggtctggagtgggtctctcttattagttgggatggtggtagcacatactatgcagactctgtgaagggtcgattcaccatctccagagacaacagcaaaaactccctgtatctgcaaatgaacagtctgagagctgaggacaccgccttgtattactgtgcaaaagata"),species="hs",frequency_uniform=c(0.01), frequency=c(4.147725e-06))

list_germline_genes$hs$ig$h$V[nrow(list_germline_genes$hs$ig$h$V) + 1,] <- df_species_VDJ[["V1"]]
list_germline_genes$hs$ig$h$V[nrow(list_germline_genes$hs$ig$h$V) + 1,] <- df_species_VDJ[["V2"]]
list_germline_genes$hs$ig$h$V[nrow(list_germline_genes$hs$ig$h$V) + 1,] <- df_species_VDJ[["V3"]]
list_germline_genes$hs$ig$h$V$gene


#suprimer les séquences qui ont des valeurs NA (annotation incomplete)
hs_ig_h_ <- na.omit(hs_ig_h_)
assign(paste("hs_ig_h_", rep[i], sep = ""), hs_ig_h_)

rm(hs_ig_h_)
rm(local_var)
rm(df_species_VDJ)
rm(dfd)
rm(dfj)
rm(dfv)
rm(freqs)
rm(file_name)
rm(repertoire)
}

ls()

repertoires <- list(hs_ig_h_H3 = hs_ig_h_H3, hs_ig_h_H4= hs_ig_h_H4, hs_ig_h_M5= hs_ig_h_M5, hs_ig_h_M6= hs_ig_h_M6, hs_ig_h_S24= hs_ig_h_S24, hs_ig_h_S26= hs_ig_h_S26)

compatibility_check(repertoire = repertoires[[1]] , species = "hs", receptor = "igh") #pour lancer ça on doit activer toujours avant la library(immuneREF)

#repertoires[[1]]

#ANALYSIS OF SINGLE FEATURES
#Setting up of reference dataframe

library(immuneREF)

#Ajouter bien les donees le meme qu'avant
setwd("/Users/equipe.lefranc/Documents/Lorena Project/COVID/DONES")


list_repertoires <- repertoires

repertoire_names <- names(list_repertoires)
repertoire_lengths <- sapply(list_repertoires,nrow)
repertoire_species <- sapply(strsplit(repertoire_names, "\\_"),function(x) x[[1]])
repertoire_receptor <- sapply(strsplit(repertoire_names, "\\_"),function(x) x[[2]])
repertoire_chain <- sapply(strsplit(repertoire_names, "\\_"),function(x) x[[3]])

input_data_ref<-data.frame(sample_id = repertoire_names,
                           nb_sequences = repertoire_lengths,
                           species = repertoire_species,
                           receptor = repertoire_receptor,
                           chain = repertoire_chain,row.names=c(1:length(repertoire_names)))

#Calculate repertoire overlap

overlap_layer<-repertoire_overlap(list_repertoires,basis="CDR3_aa")

# Subsampling of repertoires

subsample_size<-10000

list_repertoires<-subset_input_repertoires(list_repertoires=list_repertoires,
                                           subset_size=subsample_size,
                                           random=FALSE)

# Calculate all features for each repertoire

repertoires_analyzed<-list()

for(i in 1:length(list_repertoires)){
  repertoires_analyzed[[repertoire_names[i]]]<-calc_characteristics(
    repertoire_df=list_repertoires[[i]],
    species=repertoire_species[i],
    receptor=repertoire_receptor[i],
    chain=repertoire_chain[i],
    identifier_rep=repertoire_names[i])
}

#Parallelization

library(foreach)
library(doParallel)
registerDoParallel(10) 

repertoires_analyzed<-foreach(i=1:length(list_repertoires)) %dopar% {#
  repertoires_analyzed_loop<-calc_characteristics(
    repertoire_df=list_repertoires[[i]],
    species=repertoire_species[i],
    receptor=repertoire_receptor[i],
    chain=repertoire_chain[i],
    identifier_rep=repertoire_names[i])
  
  return(repertoires_analyzed_loop)
}
names(repertoires_analyzed)<-repertoire_names

#SIMILARITY SCORE CALCULATION

list_single_layers<-calculate_similarities(repertoires_analyzed=repertoires_analyzed,overlap_layer=overlap_layer)

list_single_layers<-calculate_similarities_parallel(repertoires_analyzed,overlap_layer)

#CONDENSING LAYERS INTO MULTI-LAYERS NETWORK

cormat <- condense_layers(list_single_layers,
                          weights = c(1,1,1,1,1,1),
                          method = "standard")

#VISUALIZATION OF RESULTS
# Draw heatmap of immuneREF layers

dir.create("figures3")

list_all_layers <- list_single_layers
list_all_layers[["Condensed"]] <- cormat

annotation_list<-list()
annotation_list[["categories"]]<-data.frame(Species=input_data_ref$species,
                                            Receptor = input_data_ref$receptor)

annotation_list[["colors"]]<-list(Species=c(mm='#ffffbf',hs='#fc8d59'),
                                  Receptor=c(ig='#91bfdb'))

print_heatmap_sims(list_similarity_matrices=list_all_layers,
                   annotation_list=annotation_list,
                   path_figure="figures3")


#Calculate network features of condensed immuneREF layer
network_features <- analyze_similarity_network(cormat)

# Draw Global Similarity Plots of immuneREF layers

categories_list<-list()
categories_list[["categories"]]<-input_data_ref
categories_list[["color"]]<-c("white",'#91bfdb','#ffffbf')
categories_list[["subset"]]<-"species"

print_global_similarity(list_similarity_matrices=list_all_layers,
                        categories_list = categories_list,
                        path_figure="figures3")

# Plot local similarity per category and identify max and min locally similar repertoires

max_min_reps<-print_local_similarity(list_similarity_matrices=list_all_layers,
                                     categories_list = categories_list,
                                     path_figure="figures3")


# Radar plot to visualize similarity across all 6 layers

radar_list<-list()

radar_list[["hs_ig_h_H3"]]<-repertoires_analyzed[["hs_ig_h_H3"]]
radar_list[["hs_ig_h_H4"]]<-repertoires_analyzed[["hs_ig_h_H4"]]
radar_list[["hs_ig_h_M5"]]<-repertoires_analyzed[["hs_ig_h_M5"]]
radar_list[["hs_ig_h_M6"]]<-repertoires_analyzed[["hs_ig_h_M6"]]
radar_list[["hs_ig_h_S24"]]<-repertoires_analyzed[["hs_ig_h_S24"]]
radar_list[["hs_ig_h_S26"]]<-repertoires_analyzed[["hs_ig_h_S26"]]

comparison_list<-list(roi=names(radar_list),
                      roi_names=c(
                        "HEALTHY 3",
                        "HEALTHY 4",
                        "COVID 5",
                        "COVID 6",
                        "SEVERE COVID 24",
                        "SEVERE COVID 26"),
                      ref="hs_ig_h_H3",
                      plot_names=c("HEALTHY 3","HEALTHY 4", "COVID 5", "COVID 6", "SEVERE COVID 24", "SEVERE COVID 26"),
                      colors=c("grey","blue", "red", "yellow", "pink", "green"))

comparison_list2<-list(roi=names(radar_list),
                       roi_names=c(
                         "HEALTHY 3",
                         "HEALTHY 4",
                         "COVID 5",
                         "COVID 6",
                         "SEVERE COVID 24",
                         "SEVERE COVID 26"),
                       ref="hs_ig_h_H4",
                       plot_names=c("HEALTHY 3","HEALTHY 4", "COVID 5", "COVID 6", "SEVERE COVID 24", "SEVERE COVID 26"),
                       colors=c("grey","blue", "red", "yellow", "pink", "green"))


print_repertoire_radar(list_similarity_matrices=list_single_layers,
                       to_compare=comparison_list,
                       path_figure="figures3",
                       name_plot="tutorial")
print_repertoire_radar(list_similarity_matrices=list_single_layers,
                       to_compare=comparison_list2,
                       path_figure="figures3",
                       name_plot="tutorial")

# Classical repertoire analysis of maximally and minimally similar repertoires per category

#mm_igh<-list()

#mm_igh[["mm_ig_h_2_0__0_0_0_A"]]<-repertoires_analyzed[["mm_ig_h_2_0__0_0_0_A"]]

#mm_igh[["mm_ig_h_4_0__0_0_0_A"]]<-repertoires_analyzed[["mm_ig_h_4_0__0_0_0_A"]]

#print_repertoire_comparison(list_repertoires=mm_igh,name_plots="mm_igh",aa_freq_length=14,path_figure="figures3")

hs_igh<-list()

hs_igh[["hs_ig_h_H3"]]<-repertoires_analyzed[["hs_ig_h_H3"]]

hs_igh[["hs_ig_h_H4"]]<-repertoires_analyzed[["hs_ig_h_H4"]]

hs_igh[["hs_ig_h_M5"]]<-repertoires_analyzed[["hs_ig_h_M5"]]

hs_igh[["hs_ig_h_M6"]]<-repertoires_analyzed[["hs_ig_h_M6"]]

hs_igh[["hs_ig_h_S24"]]<-repertoires_analyzed[["hs_ig_h_S24"]]

hs_igh[["hs_ig_h_S26"]]<-repertoires_analyzed[["hs_ig_h_S26"]]

print_repertoire_comparison(list_repertoires=hs_igh,
                            name_plots="hs_igh",
                            aa_freq_length=17,
                            path_figure="figures3")








