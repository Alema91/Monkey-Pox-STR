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
library(pheatmap, quietly = TRUE, warn.conflicts = FALSE)
library(compositions, quietly = TRUE, warn.conflicts = FALSE)
library(ggfortify, quietly = TRUE, warn.conflicts = FALSE)

# matrix con los valores de muestras vs STR_alelos

datos <- read.csv2("data/df_prueba_clr_2.csv", sep = "\t", header = T)
muestras<- as.character(datos[,1])

numeric_datos<- datos[,-c(1)]
indx <- sapply(numeric_datos, is.factor)
numeric_datos[indx] <- lapply(numeric_datos[indx], function(x) as.numeric(as.character(x)))

ndatos<- data.frame (
    samples = muestras,
    numeric_datos
)

# normalizacion por CLR

matrix_datos<- as.matrix(numeric_datos)
clr_datos<- clr(matrix_datos)

df_clr<- data.frame (
    samples = muestras,
    clr_datos)

# guardamos los datos transformados

write.table(df_clr, "data/df_clr_str.csv", sep = "\t", col.names = T, quote = F)

# plots con CLR

# heatmap

modelo_pca1 <- prcomp(clr_datos, scale = FALSE)
plot(modelo_pca1)
summary(modelo_pca1)

# Plot PCA

autoplot(modelo_pca1, data = df_wide, loadings.label = TRUE, loadings.label.size = 3, loadings = TRUE, shape = T, label = T, label.size = 2.5)
ggsave("plots/pca_str_samples.png", width = 18, height = 18, dpi = 300, units = c("cm"))

# dataset

row.names(clr_datos) <- df_clr$samples

# Finding distance matrix
distance_mat <- dist(clr_datos, method = "euclidean")

# Fitting Hierarchical clustering Model
# to training dataset
set.seed(240)
hier_cl <- hclust(distance_mat, method = "average")

# Plotting dendrogram
png(file = "C:/Users/alber/Desktop/Desarrollo/Monkey-Pox-STR/plots/hierclust.png", width = 600, height = 350)
plot(hier_cl)
dev.off()

########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################


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

autoplot(modelo_pca1, data = ndatos, loadings.label = TRUE, loadings.label.size = 3, loadings = TRUE, shape = T, label = T, label.size = 2.5)
ggsave("plots/pca_str_samples.png", width = 18, height = 18, dpi = 300, units = c("cm"))

# heatmap

heatmap(matrix_table)

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
png(file = "/home/alberto.lema/Documents/Desarrollo/Monkey-Pox-analysis/plots/eeeheatmap.png", width = 600, height = 350)
pheatmap(matrix_table, cluster_rows = T, cluster_cols = T)
dev.off()