#!/usr/bin/env Rscript

# library -------

library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(readxl, quietly = TRUE, warn.conflicts = FALSE)
library(tibble, quietly = TRUE, warn.conflicts = FALSE)
library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(forcats, quietly = TRUE, warn.conflicts = FALSE)
library(gridExtra, quietly = TRUE, warn.conflicts = FALSE)
library(patchwork, quietly = TRUE, warn.conflicts = FALSE)
library(reshape2, quietly = TRUE, warn.conflicts = FALSE)
library(ggfortify, quietly = TRUE, warn.conflicts = FALSE)
library(pheatmap, quietly = TRUE, warn.conflicts = FALSE)
library(viridis, quietly = TRUE, warn.conflicts = FALSE)

colorRampPalette(c("navy", "white", "firebrick3"))(50)
# utilidades
# indx <- sapply(numeric_datos, is.factor)
# numeric_datos[indx] <- lapply(numeric_datos[indx], function(x) as.numeric(as.character(x)))

# data 01
mfalleles <- read.csv2("bin/01-most_frequent_alleles_seqs.tsv", sep = "\t", header = T)
long_mfalleles <- melt(mfalleles, id.var = c("X"))
colnames(long_mfalleles) <- c("samples", "STR", "Alelles")

# Presence ausence data
paalelles <- read.csv2("bin/03-all_alleles_presence.tsv", sep = "\t", header = T)

# matrix con los valores de muestras vs STR_alelos
muestras <- as.character(paalelles[, 1])
matrix_pa <- as.matrix(paalelles[, -1])

# heatmap presencia vs ausencia
pheatmap(matrix_pa)

# CLR data
supportclr <- read.csv2("bin/06-all_alleles_supporting_reads_CLR.tsv", sep = "\t", header = T)
muestras <- as.character(supportclr[, 1])
matrix_sp <- supportclr[, -1]

indx <- sapply(matrix_sp, is.character)
matrix_sp[indx] <- lapply(matrix_sp[indx], function(x) as.numeric(x))
row.names(matrix_sp) <- muestras

# realizamos un heatmap por cada STR
long_sp <- melt(supportclr, id.var = c("X"))
colnames(long_sp) <- c("samples", "STR", "CLR")
data_sep_sp <- long_sp %>% separate(STR, c("STR", "alelle", "id"), sep = "_")
df_sep_sp <- data.frame(
    sample = data_sep_sp[, 1],
    STR_alelles = long_sp$STR,
    STR = data_sep_sp$STR,
    CLR = data_sep_sp$CLR
)

id_str <- unique(df_sep_sp$STR)
for (i in 1:length(id_str)) {
    tmp_df <- df_sep_sp[df_sep_sp$STR == id_str[i], ]
    alelos_tmp <- as.character(unique(tmp_df$STR_alelles))
    color_tmp <- colorRampPalette(c("navy", "white", "firebrick3"))(length(alelos_tmp))
    tmp_df_wide <- data.frame(pivot_wider(tmp_df, names_from = "STR_alelles", values_from = "CLR", id_cols = "sample"))
    tmp_matrix <- tmp_df_wide[, -1]
    rownames(tmp_matrix) <- muestras
    indx <- sapply(tmp_matrix, is.character) # transform to numeric
    tmp_matrix[indx] <- lapply(tmp_matrix[indx], function(x) as.numeric(x))
    nombre_str <- as.character(id_str[i])
    ruta <- paste0("/home/alberto.lema/Documents/Desarrollo/Monkey-Pox-analysis/plots/", nombre_str, ".png")
    # heatmap by STR
    png(file = ruta, width = 600, height = 350)
    pheatmap(tmp_matrix, show_rownames = T, show_colnames = F, color = color_tmp)
    dev.off()
}

# prueba con STR11

df_str11 <- df_sep_sp[df_sep_sp$STR == "STR11", ]
alelos_str11 <- as.character(unique(df_str11$STR_alelles))
color_str11 <- colorRampPalette(c("navy", "white", "firebrick3"))(length(alelos_str11))

png(file = "/home/alberto.lema/Documents/Desarrollo/Monkey-Pox-analysis/plots/STR11.png", width = 600, height = 350)
pheatmap(tmp_matrix, show_rownames = T, show_colnames = F, color = color_str11)
dev.off()

# PCA
modelo_pca1 <- prcomp(matrix_sp, scale = FALSE)
plot(modelo_pca1)
summary(modelo_pca1)

# Plot PCA

autoplot(modelo_pca1, loadings.label = TRUE, loadings.label.size = 3, loadings = TRUE, shape = T, label = T, label.size = 2.5)
ggsave("plots/pca_str_samples.png", width = 18, height = 18, dpi = 300, units = c("cm"))

# Fitting Hierarchical clustering Model
# to training dataset
set.seed(240)
hier_cl <- hclust(matrix_sp, method = "average")

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
