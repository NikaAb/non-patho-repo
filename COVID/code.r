# Comparer 6 répertoires de IGH Homo sapiens avec 3 types des données: COVID +++ (S24,S26), COVID+ (M5, M6) et COVID - (H3, H4).
# Les datas ont été téléchargés dans le ireceptor-public-archives(chercher des information sur les données).
# Les séquences ont été analysés par IMGT/HighV-quest et la sortie a été récuperé en format AIRR (vquest_airr.tsv).

library(dplyr)
library(immuneREF)
library(stringr)
library(tidyr)


getwd()

# Adresse du dossier qui contienne le fichier AIRR pour le répertoire H3
setwd("/Users/equipe.lefranc/Documents/Lorena Project/COVID/DONES/H3")

# On va stocker le fichier AIRR dans un tableau et on remplace les données vides
#par NA et on va separer les données par tabulation
H3_donnees_brut <- read.table(file = "vquest_airr.tsv", header = TRUE,
                              fill = TRUE, sep = "\t",
                              na.strings=c(""," ","NA"))

#head(H3_donnees_brut)

#names(H3_donnees_brut)

#remplacer le contennu de la colonne séquence des acides aminés par
#séquence_aligment_aa et aussi le contennu de la séquence par sequence_alignment 
#(sans gaps) car le format demandé par immuneREF est:
#“sequence_aa”: Full amino acid VDJ sequence
#“sequence”: Full nucleotide VDJ sequence

H3_donnees_brut$sequence_aa <- gsub("\\.","",
                            as.character(H3_donnees_brut$sequence_alignment_aa))
H3_donnees_brut$sequence <- gsub("\\.","",
                            as.character(H3_donnees_brut$sequence_alignment))

# On recupere les colonnes nécessaires pour lancer immuneREF
H3_format_immuneREF <- select(H3_donnees_brut, sequence_id, sequence_aa,
                              sequence, junction_aa, junction, v_call, d_call,
                              j_call)

#Ajouter la colonne de fréquence et regroupe les séquences identiques

H3_format_immuneREF <- H3_format_immuneREF %>% 
  
  group_by(sequence_aa, sequence, junction_aa, junction, v_call, d_call, 
           j_call) %>%
  
  summarise(
    
    n = n(),
    
  ) %>%
  filter(
    
    n >= 0
  )

# Calculer la fréquence et vérifier si la sume des fréquences est égale à 1
freqs <- H3_format_immuneREF$n/sum(H3_format_immuneREF$n)
sum(freqs)

#Creer un data frame
H3_df <- data.frame(H3_format_immuneREF)

#Mis en format les noms des genes de V_call. J'enleve les nom de spece er ce qui
#est apres l'étoile (allele) [HACERLO POR PARTES Y UNIRLO AL FINAL!!!!!!!!!!]

dfvH3 <- data.frame(Name = c(H3_df$v_call)) 

dfvH3[c('specie', 'v_callH3')] <- str_split_fixed(dfvH3$Name," ", 2)

dfvH3[c('v_call', 'allele v_call')] <- str_split_fixed(dfvH3$v_callH3,"\\*", 2)

dfvH3 <- dfvH3[c('v_call', 'allele v_call')]

repH3.v <- cbind(H3_df, dfvH3)

repH3.v$v_call <- NULL
repH3.v$Name <- NULL
repH3.v$v_callH3 <- NULL

#Separer les données de D

dfdH3 <- data.frame(Name = c(H3_df$d_call)) 



dfdH3[c('specie', 'd_callH3')] <- str_split_fixed(dfdH3$Name," ", 2)

dfdH3[c('d_call', 'allele d_call')] <- str_split_fixed(dfdH3$d_callH3,"\\*", 2)

dfdH3 <- dfdH3[c('d_call', 'allele d_call')]

repH3.d <- cbind(repH3.v, dfdH3)

repH3.d$d_call <- NULL
repH3.d$Name <- NULL
repH3.d$d_callH3 <- NULL
repH3.d$`allele d_call` <- NULL

#Separer les données de J

dfjH3 <- data.frame(Name = c(H3_df$j_call)) 

library(stringr)

dfjH3[c('specie', 'j_callH3')] <- str_split_fixed(dfjH3$Name," ", 2)

dfjH3[c('j_call', 'allele j_call')] <- str_split_fixed(dfjH3$j_callH3,"\\*", 2)

dfjH3 <- dfjH3[c('j_call', 'allele j_call')]

repH3.j <- cbind(repH3.d, dfjH3, freqs)

repH3.j$j_call <- NULL
repH3.j$Name <- NULL
repH3.j$j_callH3 <- NULL
repH3.j$`allele j_call` <- NULL
repH3.j$n <- NULL

#Convertir le format des colonnes de charactere à factor du répertoire H3

repertoireH3 <- data.frame(repH3.j)
repertoireH3$sequence_aa=as.factor(repertoireH3$sequence_aa)
repertoireH3$sequence=as.factor(repertoireH3$sequence)                                 
repertoireH3$junction_aa=as.factor(repertoireH3$junction_aa)                                  
repertoireH3$junction=as.factor(repertoireH3$junction)
repertoireH3$v_call=as.factor(repertoireH3$v_call)
repertoireH3$d_call=as.factor(repertoireH3$d_call)
repertoireH3$j_call=as.factor(repertoireH3$j_call)
repH3.v$`allele v_call` <- NULL

#Le format de nom de répertoire est le suivant:
#nomdespece_locus_chaine_nomrépertoire
hs_ig_h_H3 <- repertoireH3

# Si il y a un probleme avec compatibilitycheck, on peut Unifier les genes pour
#voir si ils sont dans le germline list c'est pas obligatoire de le faire.
hs_ig_h_H3v <-hs_ig_h_H3 %>% group_by(v_call) %>% filter (! duplicated(v_call))
hs_ig_h_H3j <-hs_ig_h_H3 %>% group_by(j_call) %>% filter (! duplicated(j_call))

list_germline_genes$hs$ig$h$V$gene
hs_ig_h_H3v$v_call
hs_ig_h_H3v[65,]$sequence
list_germline_genes$hs$ig$h$V$sequence
hs_ig_h_H3j$j_call

setdiff(list_germline_genes$hs$ig$h$V$gene, hs_ig_h_H3v$v_call)
setdiff(list_germline_genes$hs$ig$h$J$gene, hs_ig_h_H3j$j_call)

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
hs_ig_h_H3 <- na.omit(hs_ig_h_H3)


#INPUT FORMAT
getwd()

#Ajouter bien les donees on bessoin les resultats de donees tsv POUR H4
setwd("/Users/equipe.lefranc/Documents/Lorena Project/COVID/DONES/H4")


dataH4 <- read.table(file = "vquest_airr.tsv", header = TRUE, fill = TRUE, sep = "\t", na.strings=c(""," ","NA"))

filtered_dataH4<-dataH4[!is.na(dataH4$j_call),]

#head(dataH4)

#names(dataH4)

library(dplyr)

tabH4 <- select(filtered_dataH4, sequence_id, sequence_aa, sequence, junction_aa, junction, v_call, d_call, j_call)
library(immuneREF)

#Calculer la freq

tabFH4 <- tabH4 %>% 
  
  group_by(sequence_id, sequence_aa, sequence, junction_aa, junction, v_call, d_call, j_call) %>%
  
  summarise(
    
    n = n(),
    
  ) %>%
  filter(
    
    n >= 0
  )

freqs <- tabFH4$n/sum(tabFH4$n)
sum(freqs)

RH4 <- data.frame(tabFH4)

# Put the sequence_alignment_aa sequence, without gaps, in the sequence_aa coloumn  

library(stringr)


# Version two : without loop

#gsub is for removing gaps from the sequences in sequence_alignment_aa column
RH4$sequence_aa <- gsub("\\.","",as.character(filtered_dataH4$sequence_alignment_aa))


#Separer les données de V

dfvH4 <- data.frame(Name = c(RH4$v_call)) 

library(stringr)

dfvH4[c('specie', 'v_callH4')] <- str_split_fixed(dfvH4$Name," ", 2)

dfvH4[c('v_call', 'allele v_call')] <- str_split_fixed(dfvH4$v_callH4,"\\*", 2)

dfvH4 <- dfvH4[c('v_call', 'allele v_call')]

repH4.v <- cbind(RH4, dfvH4)

repH4.v$v_call <- NULL
repH4.v$Name <- NULL
repH4.v$v_callH4 <- NULL

#Separer les données de D

dfdH4 <- data.frame(Name = c(RH4$d_call)) 

library(stringr)

dfdH4[c('specie', 'd_callH4')] <- str_split_fixed(dfdH4$Name," ", 2)

dfdH4[c('d_call', 'allele d_call')] <- str_split_fixed(dfdH4$d_callH4,"\\*", 2)

dfdH4 <- dfdH4[c('d_call', 'allele d_call')]

repH4.d <- cbind(repH4.v, dfdH4)

repH4.d$d_call <- NULL
repH4.d$Name <- NULL
repH4.d$d_callH4 <- NULL
repH4.d$`allele d_call` <- NULL

#Separer les données de J

dfjH4 <- data.frame(Name = c(RH4$j_call)) 

library(stringr)

dfjH4[c('specie', 'j_callH4')] <- str_split_fixed(dfjH4$Name," ", 2)

dfjH4[c('j_call', 'allele j_call')] <- str_split_fixed(dfjH4$j_callH4,"\\*", 2)

dfjH4 <- dfjH4[c('j_call', 'allele j_call')]

repH4.j <- cbind(repH4.d, dfjH4, freqs)

repH4.j$j_call <- NULL
repH4.j$Name <- NULL
repH4.j$j_callH4 <- NULL
repH4.j$`allele j_call` <- NULL
repH4.j$n <- NULL

#Changer le format de répertoires

repertoireH4 <- data.frame(repH4.j)
repertoireH4$sequence_aa=as.factor(repertoireH4$sequence_aa)
repertoireH4$sequence=as.factor(repertoireH4$sequence)                                 
repertoireH4$junction_aa=as.factor(repertoireH4$junction_aa)                                  
repertoireH4$junction=as.factor(repertoireH4$junction)
repertoireH4$v_call=as.character(repertoireH4$v_call)
repertoireH4$d_call=as.factor(repertoireH4$d_call)
repertoireH4$j_call=as.factor(repertoireH4$j_call)
repH4.v$`allele v_call` <- NULL
hs_ig_h_H4 <- repertoireH4

#Unifier les genes pour voir si ils sont dans le germline list
hs_ig_h_H4v <-hs_ig_h_H4 %>% group_by(v_call) %>% filter (! duplicated(v_call))
hs_ig_h_H4j <-hs_ig_h_H4 %>% group_by(j_call) %>% filter (! duplicated(j_call))

list_germline_genes$hs$ig$h$V$gene
hs_ig_h_H4v$v_call
hs_ig_h_H4v[65,]$sequence
list_germline_genes$hs$ig$h$V$sequence
hs_ig_h_H4j$j_call

setdiff(list_germline_genes$hs$ig$h$V$gene, hs_ig_h_H4v$v_call)
setdiff(list_germline_genes$hs$ig$h$J$gene, hs_ig_h_H4j$j_call)

#replace also the sequence with the v domain sequence

# with gaps
hs_ig_h_H4$sequence_aa <- filtered_dataH4$sequence_alignment_aa   
hs_ig_h_H4$sequence <- filtered_dataH4$sequence_alignment


# without gaps
hs_ig_h_H4$sequence_aa <- gsub("\\.","",as.character(filtered_dataH4$sequence_alignment_aa))
hs_ig_h_H4$sequence <- gsub("\\.","",as.character(filtered_dataH4$sequence_alignment))

#suprimer NA
library(tidyr)
hs_ig_h_H4 <- na.omit(hs_ig_h_H4)


#INPUT FORMAT
getwd()

#Ajouter bien les donees on bessoin les resultats de donees tsv POUR M5
setwd("/Users/equipe.lefranc/Documents/Lorena Project/COVID/DONES/M5")


dataM5 <- read.table(file = "vquest_airr.tsv", header = TRUE, fill = TRUE, sep = "\t", na.strings=c(""," ","NA"))

filtered_dataM5<-dataM5[!is.na(dataM5$j_call),]

#head(dataM5)

#names(dataM5)

library(dplyr)

tabM5 <- select(filtered_dataM5, sequence_id, sequence_aa, sequence, junction_aa, junction, v_call, d_call, j_call)
library(immuneREF)

#Calculer la freq

tabFM5 <- tabM5 %>% 
  
  group_by(sequence_id, sequence_aa, sequence, junction_aa, junction, v_call, d_call, j_call) %>%
  
  summarise(
    
    n = n(),
    
  ) %>%
  filter(
    
    n >= 0
  )

freqs <- tabFM5$n/sum(tabFM5$n)
sum(freqs)

RM5 <- data.frame(tabFM5)

# Put the sequence_alignment_aa sequence, without gaps, in the sequence_aa coloumn  

library(stringr)


# Version two : without loop

#gsub is for removing gaps from the sequences in sequence_alignment_aa column
RM5$sequence_aa <- gsub("\\.","",as.character(filtered_dataM5$sequence_alignment_aa))


#Separer les données de V

dfvM5 <- data.frame(Name = c(RM5$v_call)) 

library(stringr)

dfvM5[c('specie', 'v_callM5')] <- str_split_fixed(dfvM5$Name," ", 2)

dfvM5[c('v_call', 'allele v_call')] <- str_split_fixed(dfvM5$v_callM5,"\\*", 2)

dfvM5 <- dfvM5[c('v_call', 'allele v_call')]

repM5.v <- cbind(RM5, dfvM5)

repM5.v$v_call <- NULL
repM5.v$Name <- NULL
repM5.v$v_callM5 <- NULL

#Separer les données de D

dfdM5 <- data.frame(Name = c(RM5$d_call)) 

library(stringr)

dfdM5[c('specie', 'd_callM5')] <- str_split_fixed(dfdM5$Name," ", 2)

dfdM5[c('d_call', 'allele d_call')] <- str_split_fixed(dfdM5$d_callM5,"\\*", 2)

dfdM5 <- dfdM5[c('d_call', 'allele d_call')]

repM5.d <- cbind(repM5.v, dfdM5)

repM5.d$d_call <- NULL
repM5.d$Name <- NULL
repM5.d$d_callM5 <- NULL
repM5.d$`allele d_call` <- NULL

#Separer les données de J

dfjM5 <- data.frame(Name = c(RM5$j_call)) 

library(stringr)

dfjM5[c('specie', 'j_callM5')] <- str_split_fixed(dfjM5$Name," ", 2)

dfjM5[c('j_call', 'allele j_call')] <- str_split_fixed(dfjM5$j_callM5,"\\*", 2)

dfjM5 <- dfjM5[c('j_call', 'allele j_call')]

repM5.j <- cbind(repM5.d, dfjM5, freqs)

repM5.j$j_call <- NULL
repM5.j$Name <- NULL
repM5.j$j_callM5 <- NULL
repM5.j$`allele j_call` <- NULL
repM5.j$n <- NULL

#Changer le format de répertoires

repertoireM5 <- data.frame(repM5.j)
repertoireM5$sequence_aa=as.factor(repertoireM5$sequence_aa)
repertoireM5$sequence=as.factor(repertoireM5$sequence)                                 
repertoireM5$junction_aa=as.factor(repertoireM5$junction_aa)                                  
repertoireM5$junction=as.factor(repertoireM5$junction)
repertoireM5$v_call=as.character(repertoireM5$v_call)
repertoireM5$d_call=as.factor(repertoireM5$d_call)
repertoireM5$j_call=as.factor(repertoireM5$j_call)
repM5.v$`allele v_call` <- NULL
hs_ig_h_M5 <- repertoireM5

#Unifier les genes pour voir si ils sont dans le germline list
hs_ig_h_M5v <-hs_ig_h_M5 %>% group_by(v_call) %>% filter (! duplicated(v_call))
hs_ig_h_M5j <-hs_ig_h_M5 %>% group_by(j_call) %>% filter (! duplicated(j_call))

list_germline_genes$hs$ig$h$V$gene
hs_ig_h_M5v$v_call
hs_ig_h_M5v[65,]$sequence
list_germline_genes$hs$ig$h$V$sequence
hs_ig_h_M5j$j_call

setdiff(list_germline_genes$hs$ig$h$V$gene, hs_ig_h_M5v$v_call)
setdiff(list_germline_genes$hs$ig$h$J$gene, hs_ig_h_M5j$j_call)

#replace also the sequence with the v domain sequence

# with gaps
hs_ig_h_M5$sequence_aa <- filtered_dataM5$sequence_alignment_aa   
hs_ig_h_M5$sequence <- filtered_dataM5$sequence_alignment


# without gaps
hs_ig_h_M5$sequence_aa <- gsub("\\.","",as.character(filtered_dataM5$sequence_alignment_aa))
hs_ig_h_M5$sequence <- gsub("\\.","",as.character(filtered_dataM5$sequence_alignment))

#suprimer NA
library(tidyr)
hs_ig_h_M5 <- na.omit(hs_ig_h_M5)


#INPUT FORMAT
getwd()

#Ajouter bien les donees on bessoin les resultats de donees tsv POUR M6
setwd("/Users/equipe.lefranc/Documents/Lorena Project/COVID/DONES/M6")


dataM6 <- read.table(file = "vquest_airr.tsv", header = TRUE, fill = TRUE, sep = "\t", na.strings=c(""," ","NA"))

filtered_dataM6<-dataM6[!is.na(dataM6$j_call),]

#head(dataM6)

#names(dataM6)

library(dplyr)

tabM6 <- select(filtered_dataM6, sequence_id, sequence_aa, sequence, junction_aa, junction, v_call, d_call, j_call)
library(immuneREF)

#Calculer la freq

tabFM6 <- tabM6 %>% 
  
  group_by(sequence_id, sequence_aa, sequence, junction_aa, junction, v_call, d_call, j_call) %>%
  
  summarise(
    
    n = n(),
    
  ) %>%
  filter(
    
    n >= 0
  )

freqs <- tabFM6$n/sum(tabFM6$n)
sum(freqs)

RM6 <- data.frame(tabFM6)

# Put the sequence_alignment_aa sequence, without gaps, in the sequence_aa coloumn  

library(stringr)


# Version two : without loop

#gsub is for removing gaps from the sequences in sequence_alignment_aa column
RM6$sequence_aa <- gsub("\\.","",as.character(filtered_dataM6$sequence_alignment_aa))


#Separer les données de V

dfvM6 <- data.frame(Name = c(RM6$v_call)) 

library(stringr)

dfvM6[c('specie', 'v_callM6')] <- str_split_fixed(dfvM6$Name," ", 2)

dfvM6[c('v_call', 'allele v_call')] <- str_split_fixed(dfvM6$v_callM6,"\\*", 2)

dfvM6 <- dfvM6[c('v_call', 'allele v_call')]

repM6.v <- cbind(RM6, dfvM6)

repM6.v$v_call <- NULL
repM6.v$Name <- NULL
repM6.v$v_callM6 <- NULL

#Separer les données de D

dfdM6 <- data.frame(Name = c(RM6$d_call)) 

library(stringr)

dfdM6[c('specie', 'd_callM6')] <- str_split_fixed(dfdM6$Name," ", 2)

dfdM6[c('d_call', 'allele d_call')] <- str_split_fixed(dfdM6$d_callM6,"\\*", 2)

dfdM6 <- dfdM6[c('d_call', 'allele d_call')]

repM6.d <- cbind(repM6.v, dfdM6)

repM6.d$d_call <- NULL
repM6.d$Name <- NULL
repM6.d$d_callM6 <- NULL
repM6.d$`allele d_call` <- NULL

#Separer les données de J

dfjM6 <- data.frame(Name = c(RM6$j_call)) 

library(stringr)

dfjM6[c('specie', 'j_callM6')] <- str_split_fixed(dfjM6$Name," ", 2)

dfjM6[c('j_call', 'allele j_call')] <- str_split_fixed(dfjM6$j_callM6,"\\*", 2)

dfjM6 <- dfjM6[c('j_call', 'allele j_call')]

repM6.j <- cbind(repM6.d, dfjM6, freqs)

repM6.j$j_call <- NULL
repM6.j$Name <- NULL
repM6.j$j_callM6 <- NULL
repM6.j$`allele j_call` <- NULL
repM6.j$n <- NULL

#Changer le format de répertoires

repertoireM6 <- data.frame(repM6.j)
repertoireM6$sequence_aa=as.factor(repertoireM6$sequence_aa)
repertoireM6$sequence=as.factor(repertoireM6$sequence)                                 
repertoireM6$junction_aa=as.factor(repertoireM6$junction_aa)                                  
repertoireM6$junction=as.factor(repertoireM6$junction)
repertoireM6$v_call=as.character(repertoireM6$v_call)
repertoireM6$d_call=as.factor(repertoireM6$d_call)
repertoireM6$j_call=as.factor(repertoireM6$j_call)
repM6.v$`allele v_call` <- NULL
hs_ig_h_M6 <- repertoireM6

#Unifier les genes pour voir si ils sont dans le germline list
hs_ig_h_M6v <-hs_ig_h_M6 %>% group_by(v_call) %>% filter (! duplicated(v_call))
hs_ig_h_M6j <-hs_ig_h_M6 %>% group_by(j_call) %>% filter (! duplicated(j_call))

list_germline_genes$hs$ig$h$V$gene
hs_ig_h_M6v$v_call
hs_ig_h_M6v[65,]$sequence
list_germline_genes$hs$ig$h$V$sequence
hs_ig_h_M6j$j_call

setdiff(list_germline_genes$hs$ig$h$V$gene, hs_ig_h_M6v$v_call)
setdiff(list_germline_genes$hs$ig$h$J$gene, hs_ig_h_M6j$j_call)

#replace also the sequence with the v domain sequence

# with gaps
hs_ig_h_M6$sequence_aa <- filtered_dataM6$sequence_alignment_aa   
hs_ig_h_M6$sequence <- filtered_dataM6$sequence_alignment


# without gaps
hs_ig_h_M6$sequence_aa <- gsub("\\.","",as.character(filtered_dataM6$sequence_alignment_aa))
hs_ig_h_M6$sequence <- gsub("\\.","",as.character(filtered_dataM6$sequence_alignment))

#suprimer NA
library(tidyr)
hs_ig_h_M6 <- na.omit(hs_ig_h_M6)


#INPUT FORMAT
getwd()

#Ajouter bien les donees on bessoin les resultats de donees tsv POUR S24
setwd("/Users/equipe.lefranc/Documents/Lorena Project/COVID/DONES/S24")


dataS24 <- read.table(file = "vquest_airr.tsv", header = TRUE, fill = TRUE, sep = "\t", na.strings=c(""," ","NA"))

filtered_dataS24<-dataS24[!is.na(dataS24$j_call),]

#head(dataS24)

#names(dataS24)

library(dplyr)

tabS24 <- select(filtered_dataS24, sequence_id, sequence_aa, sequence, junction_aa, junction, v_call, d_call, j_call)
library(immuneREF)

#Calculer la freq

tabFS24 <- tabS24 %>% 
  
  group_by(sequence_id, sequence_aa, sequence, junction_aa, junction, v_call, d_call, j_call) %>%
  
  summarise(
    
    n = n(),
    
  ) %>%
  filter(
    
    n >= 0
  )

freqs <- tabFS24$n/sum(tabFS24$n)
sum(freqs)

RS24 <- data.frame(tabFS24)

# Put the sequence_alignment_aa sequence, without gaps, in the sequence_aa coloumn  

library(stringr)


# Version two : without loop

#gsub is for removing gaps from the sequences in sequence_alignment_aa column
RS24$sequence_aa <- gsub("\\.","",as.character(filtered_dataS24$sequence_alignment_aa))


#Separer les données de V

dfvS24 <- data.frame(Name = c(RS24$v_call)) 

library(stringr)

dfvS24[c('specie', 'v_callS24')] <- str_split_fixed(dfvS24$Name," ", 2)

dfvS24[c('v_call', 'allele v_call')] <- str_split_fixed(dfvS24$v_callS24,"\\*", 2)

dfvS24 <- dfvS24[c('v_call', 'allele v_call')]

repS24.v <- cbind(RS24, dfvS24)

repS24.v$v_call <- NULL
repS24.v$Name <- NULL
repS24.v$v_callS24 <- NULL

#Separer les données de D

dfdS24 <- data.frame(Name = c(RS24$d_call)) 

library(stringr)

dfdS24[c('specie', 'd_callS24')] <- str_split_fixed(dfdS24$Name," ", 2)

dfdS24[c('d_call', 'allele d_call')] <- str_split_fixed(dfdS24$d_callS24,"\\*", 2)

dfdS24 <- dfdS24[c('d_call', 'allele d_call')]

repS24.d <- cbind(repS24.v, dfdS24)

repS24.d$d_call <- NULL
repS24.d$Name <- NULL
repS24.d$d_callS24 <- NULL
repS24.d$`allele d_call` <- NULL

#Separer les données de J

dfjS24 <- data.frame(Name = c(RS24$j_call)) 

library(stringr)

dfjS24[c('specie', 'j_callS24')] <- str_split_fixed(dfjS24$Name," ", 2)

dfjS24[c('j_call', 'allele j_call')] <- str_split_fixed(dfjS24$j_callS24,"\\*", 2)

dfjS24 <- dfjS24[c('j_call', 'allele j_call')]

repS24.j <- cbind(repS24.d, dfjS24, freqs)

repS24.j$j_call <- NULL
repS24.j$Name <- NULL
repS24.j$j_callS24 <- NULL
repS24.j$`allele j_call` <- NULL
repS24.j$n <- NULL

#Changer le format de répertoires

repertoireS24 <- data.frame(repS24.j)
repertoireS24$sequence_aa=as.factor(repertoireS24$sequence_aa)
repertoireS24$sequence=as.factor(repertoireS24$sequence)                                 
repertoireS24$junction_aa=as.factor(repertoireS24$junction_aa)                                  
repertoireS24$junction=as.factor(repertoireS24$junction)
repertoireS24$v_call=as.character(repertoireS24$v_call)
repertoireS24$d_call=as.factor(repertoireS24$d_call)
repertoireS24$j_call=as.factor(repertoireS24$j_call)
repS24.v$`allele v_call` <- NULL
hs_ig_h_S24 <- repertoireS24

#Unifier les genes pour voir si ils sont dans le germline list
hs_ig_h_S24v <-hs_ig_h_S24 %>% group_by(v_call) %>% filter (! duplicated(v_call))
hs_ig_h_S24j <-hs_ig_h_S24 %>% group_by(j_call) %>% filter (! duplicated(j_call))

list_germline_genes$hs$ig$h$V$gene
hs_ig_h_S24v$v_call
hs_ig_h_S24v[65,]$sequence
list_germline_genes$hs$ig$h$V$sequence
hs_ig_h_S24j$j_call

setdiff(list_germline_genes$hs$ig$h$V$gene, hs_ig_h_S24v$v_call)
setdiff(list_germline_genes$hs$ig$h$J$gene, hs_ig_h_S24j$j_call)

#replace also the sequence with the v domain sequence

# with gaps
hs_ig_h_S24$sequence_aa <- filtered_dataS24$sequence_alignment_aa   
hs_ig_h_S24$sequence <- filtered_dataS24$sequence_alignment


# without gaps
hs_ig_h_S24$sequence_aa <- gsub("\\.","",as.character(filtered_dataS24$sequence_alignment_aa))
hs_ig_h_S24$sequence <- gsub("\\.","",as.character(filtered_dataS24$sequence_alignment))

#suprimer NA
library(tidyr)
hs_ig_h_S24 <- na.omit(hs_ig_h_S24)


#INPUT FORMAT
getwd()

#Ajouter bien les donees on bessoin les resultats de donees tsv POUR S26
setwd("/Users/equipe.lefranc/Documents/Lorena Project/COVID/DONES/S26")


dataS26 <- read.table(file = "vquest_airr.tsv", header = TRUE, fill = TRUE, sep = "\t", na.strings=c(""," ","NA"))

filtered_dataS26<-dataS26[!is.na(dataS26$j_call),]

#head(dataS26)

#names(dataS26)

library(dplyr)

tabS26 <- select(filtered_dataS26, sequence_id, sequence_aa, sequence, junction_aa, junction, v_call, d_call, j_call)
library(immuneREF)

#Calculer la freq

tabFS26 <- tabS26 %>% 
  
  group_by(sequence_id, sequence_aa, sequence, junction_aa, junction, v_call, d_call, j_call) %>%
  
  summarise(
    
    n = n(),
    
  ) %>%
  filter(
    
    n >= 0
  )

freqs <- tabFS26$n/sum(tabFS26$n)
sum(freqs)

RS26 <- data.frame(tabFS26)

# Put the sequence_alignment_aa sequence, without gaps, in the sequence_aa coloumn  

library(stringr)


# Version two : without loop

#gsub is for removing gaps from the sequences in sequence_alignment_aa column
RS26$sequence_aa <- gsub("\\.","",as.character(filtered_dataS26$sequence_alignment_aa))


#Separer les données de V

dfvS26 <- data.frame(Name = c(RS26$v_call)) 

library(stringr)

dfvS26[c('specie', 'v_callS26')] <- str_split_fixed(dfvS26$Name," ", 2)

dfvS26[c('v_call', 'allele v_call')] <- str_split_fixed(dfvS26$v_callS26,"\\*", 2)

dfvS26 <- dfvS26[c('v_call', 'allele v_call')]

repS26.v <- cbind(RS26, dfvS26)

repS26.v$v_call <- NULL
repS26.v$Name <- NULL
repS26.v$v_callS26 <- NULL

#Separer les données de D

dfdS26 <- data.frame(Name = c(RS26$d_call)) 

library(stringr)

dfdS26[c('specie', 'd_callS26')] <- str_split_fixed(dfdS26$Name," ", 2)

dfdS26[c('d_call', 'allele d_call')] <- str_split_fixed(dfdS26$d_callS26,"\\*", 2)

dfdS26 <- dfdS26[c('d_call', 'allele d_call')]

repS26.d <- cbind(repS26.v, dfdS26)

repS26.d$d_call <- NULL
repS26.d$Name <- NULL
repS26.d$d_callS26 <- NULL
repS26.d$`allele d_call` <- NULL

#Separer les données de J

dfjS26 <- data.frame(Name = c(RS26$j_call)) 

library(stringr)

dfjS26[c('specie', 'j_callS26')] <- str_split_fixed(dfjS26$Name," ", 2)

dfjS26[c('j_call', 'allele j_call')] <- str_split_fixed(dfjS26$j_callS26,"\\*", 2)

dfjS26 <- dfjS26[c('j_call', 'allele j_call')]

repS26.j <- cbind(repS26.d, dfjS26, freqs)

repS26.j$j_call <- NULL
repS26.j$Name <- NULL
repS26.j$j_callS26 <- NULL
repS26.j$`allele j_call` <- NULL
repS26.j$n <- NULL

#Changer le format de répertoires

repertoireS26 <- data.frame(repS26.j)
repertoireS26$sequence_aa=as.factor(repertoireS26$sequence_aa)
repertoireS26$sequence=as.factor(repertoireS26$sequence)                                 
repertoireS26$junction_aa=as.factor(repertoireS26$junction_aa)                                  
repertoireS26$junction=as.factor(repertoireS26$junction)
repertoireS26$v_call=as.character(repertoireS26$v_call)
repertoireS26$d_call=as.factor(repertoireS26$d_call)
repertoireS26$j_call=as.factor(repertoireS26$j_call)
repS26.v$`allele v_call` <- NULL
hs_ig_h_S26 <- repertoireS26

#Unifier les genes pour voir si ils sont dans le germline list
hs_ig_h_S26v <-hs_ig_h_S26 %>% group_by(v_call) %>% filter (! duplicated(v_call))
hs_ig_h_S26j <-hs_ig_h_S26 %>% group_by(j_call) %>% filter (! duplicated(j_call))

list_germline_genes$hs$ig$h$V$gene
hs_ig_h_S26v$v_call
hs_ig_h_S26v[65,]$sequence
list_germline_genes$hs$ig$h$V$sequence
hs_ig_h_S26j$j_call

setdiff(list_germline_genes$hs$ig$h$V$gene, hs_ig_h_S26v$v_call)
setdiff(list_germline_genes$hs$ig$h$J$gene, hs_ig_h_S26j$j_call)

#replace also the sequence with the v domain sequence

# with gaps
hs_ig_h_S26$sequence_aa <- filtered_dataS26$sequence_alignment_aa   
hs_ig_h_S26$sequence <- filtered_dataS26$sequence_alignment


# without gaps
hs_ig_h_S26$sequence_aa <- gsub("\\.","",as.character(filtered_dataS26$sequence_alignment_aa))
hs_ig_h_S26$sequence <- gsub("\\.","",as.character(filtered_dataS26$sequence_alignment))

#suprimer NA
library(tidyr)
hs_ig_h_S26 <- na.omit(hs_ig_h_S26)

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
rlang::last_error()
rlang::last_trace()
