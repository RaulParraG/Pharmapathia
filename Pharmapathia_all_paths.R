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
##################################################################################################################################################################

## Ruta de las tablas PharmGKB
table_pgkb<- "./dir_que_contiene_los_files_de_pathways"
## Ruta guardado
save_path <- "./dir_donde_guardar_los_pathways"

##################################################################################################################################################################
## Obtener tabla asociación Pathway - Drug
names <- read.table(text = list.files(table_pgkb,
                             include.dirs= FALSE, 
                             full.names=TRUE), skip = 2,
                             sep = "/")

Names_PGKB <- data.frame(names, skip = 2) %>%
  select("V8") %>%
  separate(V8, into = c("Pathway", "Drug"), sep = "-") %>%
  separate(Drug, into = c("Drug", "else", sep = "_"))
colnames(Names_PGKB)[4] <- "pkpd"
Names_PGKB<- dplyr::select(Names_PGKB, c("Pathway", "Drug", "pkpd"))

write.table(Names_PGKB, file = "Pathways_ref.tsv", sep = "\t", row.names = FALSE)


############# PARA TODOS LOS ARCHIVOS
##### Crear un diccionario para moléculas y genes a partir de todos las rutas en PharmGKB

files <- list.files(table_pgkb, pattern="PA.*\\.tsv$", full.name=TRUE)
names(files) <- str_remove(basename(files), "\\.tsv$")

#### Crear una tabla con todos las Moleculas
Molecules <- files |>
  map(read_tsv) |>
  map(mutate, across(everything(), as.character)) |>
  dplyr::bind_rows(.id="filename") |>
  pivot_longer(c(From, To), values_to="molecules") |>
  select(filename, molecules) |>
  distinct() |>
  group_by(filename) |>
  mutate(coded_mol = str_c(filename, str_pad(seq(1, n(), 1), side="left", width=2, pad=0), sep="-")) |>
  ungroup()

Molecules <- 
  separate(Molecules, coded_mol, into = c("code", "else"), sep = "-") %>%
  select(-"else") %>%
  dplyr::mutate(Molecules, coded_mol = paste("N", code, str_pad(row_number(), width = 2, pad = "0"), sep = "-"))

#### Crear una tabla con todos los genes
Genes <- files |>
  map(read_tsv) |>
  map(mutate, across(everything(), as.character)) |>
  dplyr::bind_rows(.id="filename") |>
  pivot_longer(Controller, values_to="genes") |>
  select(filename, genes) |>
  distinct() |>
  group_by(filename) |>
  ungroup()

Genes <- separate(Genes, genes, into = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6", "Gene7", "Gene8", "Gene9"), sep = ",") 
Genes_u <- unique_no.NA(c(Genes$Gene1, Genes$Gene2, Genes$Gene3, Genes$Gene4, Genes$Gene5, Genes$Gene6, Genes$Gene7, Genes$Gene8, Genes$Gene9)) %>%
         as.data.frame()

No_esp <- which(!Genes_u[,1]=="")
Genes_u <- Genes_u[No_esp,] %>%
  as.data.frame()
colnames(Genes_u)[1] <- "genes"

pattern_gv <- one_or_more(UPPER) %R% one_or_more(DGT) %R% one_or_more(UPPER) %R% one_or_more(DGT)

genes_v <- sapply(stringr::str_extract_all(Genes_u$genes, pattern = pattern_gv),paste0, collapse = ' ')
Genes_t <- as.data.frame(genes_v)
Genes_t <- unique_no.NA(Genes_t$genes_v) %>% as.data.frame()
No_esp_gent <- which(!Genes_t[,1]=="")
Genes_t <- Genes_t[No_esp_gent,] %>%
  as.data.frame()
colnames(Genes_t)[1] <- "genes"
Genes_t <- AnnotationDbi::select(org.Hs.eg.db, keys = Genes_t$genes,
                                  column = "ENTREZID", keytype = "SYMBOL")
colnames(Genes_u)[3] <- "EntrezID"

No_entrez <- which(!is.na(Genes_t[,2]))
Genes_t <- Genes_t[No_entrez,] %>%
  as.data.frame()

