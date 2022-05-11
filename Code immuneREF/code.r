#QUICKSTART

library(immuneREF)

similarity_networks <- immuneREF_quickstart(repertoire_list = tutorial_repertoires)

network_features <- analyze_similarity_network(similarity_networks[["Condensed"]])

pheatmap::pheatmap(similarity_networks[["Condensed"]],scale='row')

#INPUT FORMAT
getwd()

#Ajouter bien les donees on bessoin les resultats de donees tsv
setwd("/Users/equipe.lefranc/Documents/Lorena Project/Pratique 1/DATA/Repertoire2/R2")

data <- read.table(file = "vquest_airr.tsv", header = TRUE, fill = TRUE)

head(data)

names(data)

library(dplyr)

select(data, sequence_aa, sequence, junction_aa, junction, v_call, d_call, j_call)
library(immuneREF)

freq <- data [, 2]

df <- data.frame(freq) 

table(df) 

install.packages("tidyverse")

library(tidyverse)

df %>%
  
  group_by(freq) %>%
  
  tally()
library(immuneREF)
compatibility_check(repertoire = tutorial_repertoires[[1]], species = "mm", receptor = "igh") #pour lancer ça on doit activer toujours avant la library(immuneREF)

#ANALYSIS OF SINGLE FEATURES
#Setting up of reference dataframe

library(immuneREF)

#Ajouter bien les donees le meme qu'avant
setwd("/Users/equipe.lefranc/Documents/Lorena Project/Pratique 1/DATA/Repertoire2/R2")

list_simulated_repertoires <- tutorial_repertoires

repertoire_names <- names(list_simulated_repertoires)
repertoire_lengths <- sapply(list_simulated_repertoires,nrow)
repertoire_species <- sapply(strsplit(repertoire_names, "\\_"),function(x) x[[1]])
repertoire_receptor <- sapply(strsplit(repertoire_names, "\\_"),function(x) x[[2]])
repertoire_chain <- sapply(strsplit(repertoire_names, "\\_"),function(x) x[[3]])

input_data_ref<-data.frame(sample_id = repertoire_names,
                           nb_sequences = repertoire_lengths,
                           species = repertoire_species,
                           receptor = repertoire_receptor,
                           chain = repertoire_chain,row.names=c(1:length(repertoire_names)))

#Calculate repertoire overlap

overlap_layer<-repertoire_overlap(list_simulated_repertoires,basis="CDR3_aa")

# Subsampling of repertoires

subsample_size<-10000

list_simulated_repertoires<-subset_input_repertoires(list_repertoires=list_simulated_repertoires,
                                                     subset_size=subsample_size,
                                                     random=FALSE)

#Calculation of the remaining 5 features

repertoires_analyzed<-list()
for(i in 1:length(list_simulated_repertoires)){
  repertoires_analyzed[[repertoire_names[i]]]<-calc_characteristics(
    repertoire_df=list_simulated_repertoires[[i]],
    species=repertoire_species[i],
    receptor=repertoire_receptor[i],
    chain=repertoire_chain[i],
    identifier_rep=repertoire_names[i])
}

#Parallelization

library(foreach)
library(doParallel)
registerDoParallel(10) 

repertoires_analyzed<-foreach(i=1:length(list_simulated_repertoires)) %dopar% {#
  repertoires_analyzed_loop<-calc_characteristics(
    repertoire_df=list_simulated_repertoires[[i]],
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

dir.create("figures")

list_all_layers <- list_single_layers
list_all_layers[["Condensed"]] <- cormat

annotation_list<-list()
annotation_list[["categories"]]<-data.frame(Species=input_data_ref$species,
                                            Receptor = input_data_ref$receptor)

annotation_list[["colors"]]<-list(Species=c(mm='#ffffbf',hs='#fc8d59'),
                                  Receptor=c(ig='#91bfdb'))

print_heatmap_sims(list_similarity_matrices=list_all_layers,
                   annotation_list=annotation_list,
                   path_figure="figures")


#Calculate network features of condensed immuneREF layer
network_features <- analyze_similarity_network(cormat)

# Draw Global Similarity Plots of immuneREF layers

categories_list<-list()
categories_list[["categories"]]<-input_data_ref
categories_list[["color"]]<-c("white",'#91bfdb','#ffffbf')
categories_list[["subset"]]<-"species"

print_global_similarity(list_similarity_matrices=list_all_layers,
                        categories_list = categories_list,
                        path_figure="figures")

# Plot local similarity per category and identify max and min locally similar repertoires

max_min_reps<-print_local_similarity(list_similarity_matrices=list_all_layers,
                                     categories_list = categories_list,
                                     path_figure="figures")


# Radar plot to visualize similarity across all 6 layers

radar_list<-list()
radar_list[["mm_ig_h_2_0__0_0_0_A"]]<-repertoires_analyzed[["mm_ig_h_2_0__0_0_0_A"]]
radar_list[["mm_ig_h_4_0__0_0_0_A"]]<-repertoires_analyzed[["mm_ig_h_4_0__0_0_0_A"]]
radar_list[["hs_ig_h_2_0__0_0_0_A"]]<-repertoires_analyzed[["hs_ig_h_2_0__0_0_0_A"]]
radar_list[["hs_ig_h_4_0__0_0_0_A"]]<-repertoires_analyzed[["hs_ig_h_4_0__0_0_0_A"]]

comparison_list<-list(roi=names(radar_list),
                      roi_names=c(
                        "Murine A",
                        "Murine B",
                        "Human A",
                        "Human B"),
                      ref="mm_ig_h_2_0__0_0_0_A",
                      plot_names=c("Murine A", "Murine B","Human A","Human B"),
                      colors=c("grey","blue",'red',"green"))


print_repertoire_radar(list_similarity_matrices=list_single_layers,
                       to_compare=comparison_list,
                       path_figure="figures",
                       name_plot="tutorial")

# Classical repertoire analysis of maximally and minimally similar repertoires per category

mm_igh<-list()

mm_igh[["mm_ig_h_2_0__0_0_0_A"]]<-repertoires_analyzed[["mm_ig_h_2_0__0_0_0_A"]]

mm_igh[["mm_ig_h_4_0__0_0_0_A"]]<-repertoires_analyzed[["mm_ig_h_4_0__0_0_0_A"]]

print_repertoire_comparison(list_repertoires=mm_igh,name_plots="mm_igh",aa_freq_length=14,path_figure="figures")

hs_igh<-list()

hs_igh[["hs_ig_h_2_0__0_0_0_A"]]<-repertoires_analyzed[["hs_ig_h_2_0__0_0_0_A"]]

hs_igh[["hs_ig_h_4_0__0_0_0_A"]]<-repertoires_analyzed[["hs_ig_h_4_0__0_0_0_A"]]

print_repertoire_comparison(list_repertoires=hs_igh,
                            name_plots="hs_igh",
                            aa_freq_length=17,
                            path_figure="figures")

