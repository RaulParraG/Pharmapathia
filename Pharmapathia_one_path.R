### Title: "Pharmapathia - Genérico repositorio"
#Author: "Raul Parra Garces"
#date: "05/11/2022"
# Resumen: Pharmapathia pretende aplicar una capa adicional al algoritmo ya desarrollado de HiPathia, 
#de modo que pueda integrar la información procedente de la acción farmacogenética y farmacocinética de fármacos, 
#desde un punto de vista mecanístico. Para ello se implementarán las vías de metabolización disponibles en PharmGKB".

#Instalacion de dependencias y librerias
if(!require(plyr)){
  install.packages("plyr")
  library(plyr)
}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(tidyverse)){
  install.packages("tidyverse")
  library(tidyverse)
}
if(!require(readr)){
  install.packages("readr")
  library(readr)
}
if(!require(data.table)){
  install.packages("data.table")
  library(data.table)
}
if(!require(purrr)){
  install.packages("purrr")
  library(purrr)
}
if(!require(rebus)){
  install.packages("rebus")
  library(rebus)
}
if(!require(conflicted)){
  install.packages("conflicted")
  library(conflicted)
}
if(!require(do)){
  install.packages("do")
  library(do)
}
if(!require(org.Hs.eg.db)){
  install.packages("org.Hs.eg.db")
  library(org.Hs.eg.db)
}
if(!require(AnnotationDbi)){
  install.packages("AnnotationDbi")
  library(AnnotationDbi)
}


##################### FORMATO PARA RUTA UNICA
# Farmaco. Aqui se carga el archivo para parsear

pathfiles<-list.files(path="./dir_que_contiene_los_files_de_pathways", full.names=T)
codigo por pathfiles[i] y todo dentro de un  
for(i in 1:length(pathfiles)){
Farmaco <- read_tsv(pathfiles[i],
                 col_names = TRUE, na = c("","NA")) %>%
  dplyr::select(c("From", "Reaction Type", "To", "Controller")) %>%
  separate(Controller, into = c("Gene1", "Gene2", "Gene3"), sep = ", ") ## Verificar si ", " o ","

######### Unificar las palabras de tipo de reacción
Farmaco[Farmaco == "Biochemical Reaction" | Farmaco == "Transport" | Farmaco == "Leads To"] <- "activation"
Farmaco[Farmaco == "Inhibition"] <- "inhibition"

##### Separo los elementos únicos de las columnas "From" y "To" para crear el diccionario de moléculas
Farmaco_mol <- unique(c(Farmaco$From, Farmaco$To)) %>%
  as.data.frame()
colnames(Farmaco_mol)[1] <- "molecules"
Farmaco_mol <- dplyr::mutate(Farmaco_mol, coded_mol = paste("N-PAxxx", str_pad(row_number(), width = 2, pad = "0"), sep = "-"))  
Farmaco_mol <- add_column(Farmaco_mol, "NA" = NA, .after="coded_mol")
}

## Verificar que "N-PAxxx" coincide con la ruta seleccionada

### Diccionario moleculas 
# Columna "From"
Farmaco <- add_column(Farmaco, map_mol = NA, .after="From")
Farmaco[,2] <- plyr::mapvalues(Farmaco$From, from = Farmaco_mol[,1], to= Farmaco_mol[,2])

# Columna "To"
Farmaco <- add_column(Farmaco, map_mol2 = NA, .after="To")
Farmaco[,5] <- plyr::mapvalues(Farmaco$To, from = Farmaco_mol[,1], to= Farmaco_mol[,2])

##### Diccionario de genes 
Farmaco_gen <- unique_no.NA(c(Farmaco$Gene1, Farmaco$Gene2, Farmaco$Gene3)) %>%
  as.data.frame()
colnames(Farmaco_gen)[1] <- "genes"


#### Anotacion de EntrezID
Farmaco_gen <- AnnotationDbi::select(org.Hs.eg.db, keys = Farmaco_gen$genes,
                                  column = "ENTREZID", keytype = "SYMBOL")
colnames(Farmaco_gen)[2] <- "EntrezID"

Farmaco_gen <- dplyr::mutate(Farmaco_gen, coded_gen = paste("N-PAxxx", str_pad(row_number(), width = 3, pad = "0"), sep = "-")) 
## Verificar que "N-PAxxx" coincide con la ruta seleccionada

Farmaco_dic <- as.data.frame(rbind(c(Farmaco_mol, Farmaco_gen)))

Farmaco <- add_column(Farmaco, "Entrez1"=NA, "coded_gene1"=NA, .after="Gene1")
Farmaco <- add_column(Farmaco, "Entrez2"=NA, "coded_gene2"=NA, .after="Gene2")
Farmaco <- add_column(Farmaco, "Entrez3"=NA, "coded_gene3"=NA, .after="Gene3")
Farmaco$Entrez1 <- mapvalues(Farmaco$Gene1, from = Farmaco_gen[,1], to= Farmaco_gen[,2])
Farmaco$Entrez2 <- mapvalues(Farmaco$Gene2, from = Farmaco_gen[,1], to= Farmaco_gen[,2])
Farmaco$Entrez3 <- mapvalues(Farmaco$Gene3, from = Farmaco_gen[,1], to= Farmaco_gen[,2])

Farmaco$coded_gene1 <- mapvalues(Farmaco$Gene1, from = Farmaco_gen[,1], to= Farmaco_gen[,3])
Farmaco$coded_gene2 <- mapvalues(Farmaco$Gene2, from = Farmaco_gen[,1], to= Farmaco_gen[,3])
Farmaco$coded_gene3 <- mapvalues(Farmaco$Gene3, from = Farmaco_gen[,1], to= Farmaco_gen[,3])

pattern_gv <- one_or_more(DGT) %R% END
Gene_extracted2 <- as.data.frame(sapply(stringr::str_extract_all(Farmaco$coded_gene2, pattern = pattern_gv),paste0, collapse = ''))
Farmaco$coded_gene2 <- mapvalues(Farmaco$coded_gene2, from = Farmaco$coded_gene2, to= Gene_extracted2[,1])
Gene_extracted3<- as.data.frame(sapply(stringr::str_extract_all(Farmaco$coded_gene3, pattern = pattern_gv),paste0, collapse = ''))
Farmaco$coded_gene3 <- mapvalues(Farmaco$coded_gene3, from = Farmaco$coded_gene3, to= Gene_extracted3[,1])


#### Unifico nuevamente la Farmaco
Farmaco_coded <- unite(Farmaco, "coded_gene", c(coded_gene1, coded_gene2, coded_gene3), sep = " ", remove = TRUE, na.rm = TRUE) %>%
  dplyr::select(c("map_mol", "Reaction Type", "map_mol2", "coded_gene"))
### Completo la ruta con los pasos intermedios 
#1
No_NA <- which(!Farmaco_coded[,4]=="")
Farmaco_expanded <- Farmaco_coded[-No_NA,]

for (i in 1:length(No_NA)) {
  temp1 <- cbind(Farmaco_coded[No_NA[i],1], Farmaco_coded[No_NA[i],2], Farmaco_coded[No_NA[i],4], Farmaco_coded[No_NA[i],4])
  temp2 <- cbind(Farmaco_coded[No_NA[i],4], Farmaco_coded[No_NA[i],2], Farmaco_coded[No_NA[i],3], Farmaco_coded[No_NA[i],4])
  colnames(temp1)=colnames(temp2) <- colnames(Farmaco_expanded)
  Farmaco_expanded <- rbind(Farmaco_expanded, temp1, temp2)}
Farmaco_expanded <- dplyr::select(Farmaco_expanded, c("map_mol", "Reaction Type", "map_mol2")) 

############# Guardo las Farmacos en formato .txt
## Farmaco referencia
write.table(Farmaco, file = "Farmaco_ref.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)

## Farmaco interacciones (SIF)
write.table(Farmaco_expanded, file = "Farmaco_pk.sif", sep = "\t",
            row.names = FALSE, col.names = FALSE)

## Farmaco att
Farmaco_att <- unique(c(Farmaco_expanded$map_mol, Farmaco_expanded$map_mol2)) %>%
  as.data.frame()
colnames(Farmaco_att)[1] <- "ID"
#Farmaco_att <- as.data.frame(setdiff(Farmaco_expanded[,1], Farmaco_expanded[,3]))
Farmaco_att <- add_column(Farmaco_att, "label"=NA, "X"=NA, "Y"=NA, "color"= NA, "shape"= NA, "type"= NA, "label.cex"=NA, "label.color"=NA, "width"=NA, "height"=NA, "genesList"=NA,.after="ID") 
Farmaco_att$color <- "white"
Farmaco_att$shape

#Comprobación de archivos
setdiff(Farmaco_att$ID,c(Farmaco_expanded$map_mol, Farmaco_expanded$map_mol2))
setdiff(c(Farmaco_expanded$map_mol, Farmaco_expanded$map_mol2), Farmaco_att$ID)

write.table(Farmaco_att, file = "Farmaco1.tsv", sep = "\t",
            row.names = FALSE, col.names = FALSE)

## Farmaco Gene List
write.table(Farmaco_gen[,2], file = "FarmacogeneList", sep = "\t",
            row.names = FALSE, col.names = FALSE)
  

sessionInfo()
