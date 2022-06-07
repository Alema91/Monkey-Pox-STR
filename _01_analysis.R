#!/usr/bin/env Rscript

# library -------

library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(readxl, quietly = TRUE, warn.conflicts = FALSE)
library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(readODS, quietly = TRUE, warn.conflicts = FALSE)
library(readr, quietly = TRUE, warn.conflicts = FALSE)
library(janitor, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(viridis, quietly = TRUE, warn.conflicts = FALSE)
library(rgr, quietly = TRUE, warn.conflicts = FALSE)
library(ggfortify, quietly = TRUE, warn.conflicts = FALSE)


# data --------

### Direccion ----

dir_excel <- list.files(path = "Data", pattern = "xlsx", full.names = TRUE, recursive = TRUE, include.dirs = FALSE)

### long encoded ----
long_full <- read_excel(dir_excel[1], sheet = 1)
long_encoded <- read_excel(dir_excel[2], sheet = 1)
datos <- data.frame(long_encoded[, c(1:3)], long_full[, 2], long_encoded[, c(4)])
datos$sequence_2 <- as.numeric(as.factor(datos$Sequence_OG))
datos$sequence_3 <- factor(paste(datos$Sample_name_num, datos$STR_mark, sep = "_"))

# datos$sequence_og_3<- factor(paste(datos$Sample_name_num, datos$STR_mark, datos$sequence_og_2, sep="_"))
# datos$sequence_og_3<- factor(datos$Sample_name_num)

# filter

# filter_datos<- subset(datos, AlleleFrequency >= 0.5 & STR_mark == "STR9")
# filter_datos<- subset(datos[,c(2,6,7)], AlleleFrequency >= 0.5)
filter_datos <- subset(datos[, c(2, 6, 7)])
head(filter_datos)
# filter_datos<- subset(datos[,c(2,6,7)])


# histograma de las frecuencia alelica

ggplot(filter_datos, aes(AlleleFrequency)) +
    geom_histogram()

# wide table

pivot_table <- pivot_wider(filter_datos, names_from = "sequence_3", values_from = "AlleleFrequency", id_cols = "sequence_2")


pivot_table[is.na(pivot_table)] <- 0.0001
pivot_table$Sequence <- factor(pivot_table$sequence_2)

# Intento de PCA

matrix_table <- as.matrix(pivot_table[, -1])
rownames(matrix_table) <- NULL

modelo_pca1 <- prcomp(t(matrix_table), scale = FALSE)
plot(modelo_pca1)
summary(modelo_pca1)



############ 33

# Pruebas

df <- iris[1:4]
data_df <- iris
str(data_df)