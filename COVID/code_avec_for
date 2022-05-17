# Comparer 6 répertoires de IGH Homo sapiens avec 3 types des données: COVID +++ (S24,S26), COVID+ (M5, M6) et COVID - (H3, H4).
# Les datas ont été téléchargés dans le ireceptor-public-archives(chercher des information sur les données).
# Les séquences ont été analysés par IMGT/HighV-quest et la sortie a été récuperé en format AIRR (vquest_airr.tsv).

library(dplyr)
library(immuneREF)
library(stringr)
library(tidyr)
setwd("C:/Users/equipe.lefranc/Documents")
a<-"C:/Users/equipe.lefranc/Documents/Lorena Project/COVID/DONES +"/H3
setwd(a)

getwd()
rep <- c("H3", "H4", "M5", "M6", "S24", "S26")
# Adresse du dossier qui contienne le fichier AIRR pour le répertoire H3
setwd("/Users/equipe.lefranc/Documents/Lorena Project/COVID/DONES")
list.files()

for (i in 1:length(rep)){
  print(i)
  print(rep[i])
}
  rep[1] <- setwd("H3")
  print(# On va stocker le fichier AIRR dans un tableau et on remplace les données vides
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
}
   
