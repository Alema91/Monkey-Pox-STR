#!/usr/bin/env Rscript

# library -------

library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(readxl, quietly = TRUE, warn.conflicts = FALSE)
library(janitor, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(forcats, quietly = TRUE, warn.conflicts = FALSE)
library(gridExtra, quietly = TRUE, warn.conflicts = FALSE)
library(patchwork, quietly = TRUE, warn.conflicts = FALSE)
library(reshape2, quietly = TRUE, warn.conflicts = FALSE)
library(ggfortify, quietly = TRUE, warn.conflicts = FALSE)

# data --------

### Direccion ----

dir_excel <- list.files(path = "Data", pattern = "xlsx", full.names = TRUE, recursive = TRUE, include.dirs = FALSE)

### long encoded ----
long_full <- read_excel(dir_excel[1], sheet = 1)
long_encoded <- read_excel(dir_excel[2], sheet = 1)

nombres_df <- as.data.frame(str_split(long_full$Sample_name, "_", simplify = TRUE))
datos$sequence_nueva <- factor(paste(nombres_df$V1, nombres_df$V2, sep = "_"))

datos <- data.frame(long_encoded[, c(1:3)], long_full[, 2], long_encoded[, c(4)])
datos$sequence_2 <- as.numeric(as.factor(datos$Sequence_OG))
datos$sequence_3 <- factor(paste(datos$Sample_name_num, datos$STR_mark, sep = "_"))

# datos$sequence_og_3<- factor(paste(datos$Sample_name_num, datos$STR_mark, datos$sequence_og_2, sep="_"))
# datos$sequence_og_3<- factor(datos$Sample_name_num)

# filter

# filter_datos<- subset(datos, AlleleFrequency >= 0.5 & STR_mark == "STR9")
# filter_datos<- subset(datos[,c(2,6,7)], AlleleFrequency >= 0.5)
filter_datos <- subset(datos[, c(2, 6, 7)])

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

##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

# matrix con los valores de STR y muestras

datos <- read.csv2("bin/abundant_alelles.tsv", sep = "\t")
colnames(datos)

# tengo que cambiar la tabla a long y despues otra vez a wide

df_melt <- melt(datos, id.var = c("X"))
df_melt$value <- as.numeric(as.factor(df_melt$value))

df_wide <- pivot_wider(df_melt, names_from = "variable", values_from = "value", id_cols = "X")

# Intento de PCA

df_wide <- data.frame(df_wide)

matrix_table <- as.matrix(df_wide[, -1])
row.names(matrix_table) <- as.character(df_wide[, 1])

modelo_pca1 <- prcomp(matrix_table, scale = FALSE)
plot(modelo_pca1)
summary(modelo_pca1)

# Plot PCA

autoplot(modelo_pca1, data = df_wide, loadings.label = TRUE, loadings.label.size = 3, loadings = TRUE, shape = T, label = T, label.size = 2.5)
ggsave("plots/pca_str_samples.png", width = 18, height = 18, dpi = 300, units = c("cm"))

# heatmap

heatmap(matrix_table)
ggsave("plots/heatmap_str_samples.png", width = 18, height = 18, dpi = 300, units = c("cm"))

# intento de hierathical clustering

# dataset

matrix_table <- as.matrix(df_wide[, -1])
row.names(matrix_table) <- as.character(df_wide[, 1])

# Finding distance matrix
distance_mat <- dist(matrix_table, method = "euclidean")

# Fitting Hierarchical clustering Model
# to training dataset
set.seed(240)
hier_cl <- hclust(distance_mat, method = "average")
hier_cl

# Plotting dendrogram
plot(hier_cl)

# plot heatmap
heatmap(matrix_table)