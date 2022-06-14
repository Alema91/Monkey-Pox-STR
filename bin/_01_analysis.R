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
library(ggdendro, quietly = TRUE, warn.conflicts = FALSE)

# utilidades
# indx <- sapply(numeric_datos, is.factor)
# numeric_datos[indx] <- lapply(numeric_datos[indx], function(x) as.numeric(as.character(x)))

# data 01 - alelos + frecuentes para cada muestra y STR
mfalleles <- read.csv2("bin/01-most_frequent_alleles_seqs.tsv", sep = "\t", header = T)
long_mfalleles <- melt(mfalleles, id.var = c("X"))
colnames(long_mfalleles) <- c("samples", "STR", "Alelles")

# calculo del 80%
table_percNone <- function(in_str, in_sample) {
    id_str <- unique(as.character(in_str))
    matrix_str_n <- matrix(0, nrow = length(id_str), ncol = 2)
    for (i in 1:length(id_str)) {
        df_tmp <- long_mfalleles[long_mfalleles$STR == id_str[i] & long_mfalleles$Alelles == "None", ]
        if (identical(df_tmp$sample, character(0))) {
            matrix_str_n[i, 1] <- id_str[i]
            matrix_str_n[i, 2] <- 0
        } else {
            matrix_str_n[i, 1] <- id_str[i]
            matrix_str_n[i, 2] <- round(as.numeric(dim(df_tmp)[1]) / 34, 2)
        }
    }
    id_sample <- unique(as.character(in_sample))
    matrix_sample_n <- matrix(0, nrow = length(id_sample), ncol = 2)
    for (i in 1:length(id_sample)) {
        df_tmp <- as.matrix(mfalleles[mfalleles$X == id_sample[i], ])
        colnames(aa) <- NULL
        matrix_sample_n[i, 1] <- id_sample[i]
        matrix_sample_n[i, 2] <- round(as.numeric(table(df_tmp[1, ] == "None")[2]) / 10, 2)
    }

    df_str_none <- data.frame(
        STR = as.character(matrix_str_n[, 1]),
        PerNone = as.numeric(matrix_str_n[, 2])
    )
    print(df_str_none)

    df_sample_none <- data.frame(
        samples = as.character(matrix_sample_n[, 1]),
        PerNone = as.numeric(matrix_sample_n[, 2])
    )
    print(df_sample_none)
}

table_percNone(long_mfalleles$STR, long_mfalleles$samples)

# Eliminamos el STR3, STR1, STR4
# Eliminamos las muestras 390, 350

mf_filter <- mfalleles[-c(1, 30), -c(3, 5, 6)]
long_mf_filter <- melt(mf_filter, id.var = c("X"))
colnames(long_mf_filter) <- c("samples", "STR", "alelles")
long_mf_filter$alelles <- as.numeric(as.factor(long_mf_filter$alelles))

wide_mf_filter <- as.data.frame(pivot_wider(long_mf_filter, names_from = "STR", values_from = "alelles", id_cols = "samples"))
matrix_mf_filter <- as.matrix(wide_mf_filter[, -1])
nombres <- c(
    "395",          "399",          "353_R_miseq", "403",
    "407",          "351_R",          "347_R", "411",
    "415",          "416",          "417", "418",
    "419",          "420",          "353_novaseq", "431",
    "435",          "438",          "353_R_novaseq", "441",
    "352",          "447",          "453", "349_R",
    "457",          "345",          "347", "349",
    "351",          "352",          "353_R_nanopore", "350_R"
)
row.names(matrix_mf_filter) <- nombres

# PCA
modelo_pca1 <- prcomp(matrix_mf_filter, scale = FALSE)
plot(modelo_pca1)
summary(modelo_pca1)

autoplot(modelo_pca1, loadings.label = TRUE, loadings.label.size = 3, loadings = TRUE, shape = T, label = T, label.size = 2.5)
ggsave("plots/pca_str_samples.png", width = 18, height = 18, dpi = 300, units = c("cm"))

# Hclustering

## Matriz de distancias
dist_mf_filter <- dist(matrix_mf_filter, method = "manhattan")
## Clúster jerárquico
hc_mf_filter <- hclust(dist_mf_filter)

## Dendrograma
png(file = "/home/alberto.lema/Documents/Desarrollo/Monkey-Pox-analysis/plots/hc_mf_filter.png", height = 400, width = 500)
ggdendrogram(hc_mf_filter, rotate = T, type = "rectangle")
dev.off()

# data 02 - alelos + frecuentes para cada muestra y STR
afalleles <- read.csv2("bin/02-most_frequent_alleles_freqs.tsv", sep = ",", header = T)
af_filter <- afalleles[-c(1, 30), -c(3, 5, 6)]
matrix_af <- af_filter[, -1]

indx <- sapply(matrix_af, is.character)
matrix_af[indx] <- lapply(matrix_af[indx], function(x) as.numeric(as.character(x)))

# heatmap
pheatmap(matrix_af, show_rownames = T, show_colnames = T, color = viridis(50))

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

# negative
`%notin%` <- Negate(`%in%`)

# Filtrar por los STR mayores del 80% de None
str_mayor80 <- c("STR3", "STR1", "STR4")
df_sep_sp_filter <- df_sep_sp[df_sep_sp$STR %notin% str_mayor80, ]

# heatmap

id_str <- unique(df_sep_sp_filter$STR)
for (i in 1:length(id_str)) {
    tmp_df <- df_sep_sp[df_sep_sp_filter$STR == id_str[i], ]
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

png(file = "/home/alberto.lema/Documents/Desarrollo/Monkey-Pox-analysis/plots/STR11_prueba.png", width = 600, height = 350)
pheatmap(tmp_matrix, show_rownames = T, show_colnames = F, color = color_str11)
dev.off()

df_wide_full <- data.frame(pivot_wider(df_sep_sp, names_from = "STR_alelles", values_from = "CLR", id_cols = "sample"))
full_matrix <- df_wide_full[, -1]
rownames(full_matrix) <- muestras
indx <- sapply(full_matrix, is.character) # transform to numeric
full_matrix[indx] <- lapply(full_matrix[indx], function(x) as.numeric(x))

modelo_pca1 <- prcomp(full_matrix, scale = FALSE)
plot(modelo_pca1)
summary(modelo_pca1)

autoplot(modelo_pca1, loadings.label = TRUE, loadings.label.size = 3, loadings = TRUE, shape = T, label = T, label.size = 2.5)
ggsave("plots/pca_str_samples.png", width = 18, height = 18, dpi = 300, units = c("cm"))



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
