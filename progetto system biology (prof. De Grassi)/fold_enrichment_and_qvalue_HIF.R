setwd("C:/Users/carme/Desktop/Bioinf/3° semestre/Systems Biology/progetto esame/")

# RISULTATI PER HIF1A

# df <- read.csv("Cleaned_Annotated_ChIP-seq_Data.csv")
df_err<- read.table("merged_HIF1A_A549.tabular", header = TRUE, sep = "\t")
df <- df_err[!grepl("\\bUTR\\b", df_err$Annotation), ]


# Assicuriamoci che le colonne numeriche siano effettivamente numeriche
df$Fold_Enrichment <- as.numeric(df$Fold_Enrichment)
df$X.log10.qvalue. <- as.numeric(df$X.log10.qvalue.)

attach(df)
annotations <- unique(Annotation)

# Test di Shapiro-Wilk per ciascun gruppo di annotation
shapiro_results_fold <- data.frame(Annotation = character(), p_value = numeric(), stringsAsFactors = FALSE)
shapiro_results_qvalue <- data.frame(Annotation = character(), p_value = numeric(), stringsAsFactors = FALSE)

for (cat in annotations) {
  subset_data <- df[Annotation == cat, ]  # Filtra i dati per ogni gruppo di annotation
  
  if (nrow(subset_data) > 3) {  # Shapiro-Wilk richiede almeno 3 dati
    # Test per Fold Enrichment
    p_fold <- shapiro.test(subset_data$Fold_Enrichment)$p.value
    shapiro_results_fold <- rbind(shapiro_results_fold, data.frame(Annotation = cat, p_value = p_fold))
    
    # Test per -log10(q-value)
    p_qvalue <- shapiro.test(subset_data$X.log10.qvalue.)$p.value
    shapiro_results_qvalue <- rbind(shapiro_results_qvalue, data.frame(Annotation = cat, p_value = p_qvalue))
  }
}

print(shapiro_results_fold)
print(shapiro_results_qvalue)

kruskal.test(Fold_Enrichment~Annotation)
kruskal.test(X.log10.qvalue.~Annotation)

# dato che 5'-UTR e 3'-UTR avevano solo 1 riga associata, il pairwise.wilcox.test dava errorre. 
# Per questo il df è stato filtrato rimuovendo quelle righe

library(dplyr)
df_filtered <- df %>%
  group_by(Annotation) %>%
  filter(n() > 1)

pairwise.wilcox.test(df_filtered$Fold_Enrichment, df_filtered$Annotation, p.adjust.method = "BH")
pairwise.wilcox.test(df_filtered$X.log10.qvalue., df_filtered$Annotation, p.adjust.method = "BH")

# pairwise.wilcox.test(Fold_Enrichment, Annotation, p.adjust.method = "BH")
# pairwise.wilcox.test(X.log10.qvalue., Annotation, p.adjust.method = "BH")


detach(df)




# RISULTATI PER HIF2A

# Caricare il dataset
# df2 <- read.csv("Final_Annotated_ChIP-seq_Data_HIF2A.csv")
df2 <- read.table("merged_HIF2A_A549.tabular", header = TRUE, sep = "\t")
df2 <- df_err2[!grepl("\\bUTR\\b", df_err2$Annotation), ]

# Assicuriamoci che le colonne numeriche siano effettivamente numeriche
df2$X_Fold_Enrichment <- as.numeric(df2$X_Fold_Enrichment)
df2$X_.log10.qvalue. <- as.numeric(df2$X_.log10.qvalue.)

# Attaccare il dataset per facilitare l'accesso alle variabili
attach(df2)

# Lista unica delle categorie di annotazione
categories <- unique(X_Annotation)

# Inizializzazione dei dataframe per i risultati
shapiro_results_enrichment <- data.frame(X_Annotation = character(), p_value = numeric(), stringsAsFactors = FALSE)
shapiro_results_qvalue <- data.frame(X_Annotation = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Test di Shapiro-Wilk per ciascun gruppo di annotazione
for (group in categories) {
  subset_data <- df2[X_Annotation == group, ]  # Filtra i dati per ogni gruppo di annotation
  
  if (nrow(subset_data) > 3) {  # Shapiro-Wilk richiede almeno 3 dati
    # Test per X_fold_enrichment
    p_enrichment <- shapiro.test(subset_data$X_Fold_Enrichment)$p.value
    shapiro_results_enrichment <- rbind(shapiro_results_enrichment, data.frame(X_Annotation = group, p_value = p_enrichment))
    
    # Test per X_.log10.qvalue.
    p_qvalue <- shapiro.test(subset_data$X_.log10.qvalue.)$p.value
    shapiro_results_qvalue <- rbind(shapiro_results_qvalue, data.frame(X_Annotation = group, p_value = p_qvalue))
  }
}

# Stampare i risultati dei test di Shapiro-Wilk
print(shapiro_results_enrichment)
print(shapiro_results_qvalue)

# Test di Kruskal-Wallis
kruskal.test(X_Fold_Enrichment ~ X_Annotation)
kruskal.test(X_.log10.qvalue. ~ X_Annotation)

# Test di Wilcoxon con correzione per confronti multipli (Benjamini-Hochberg)
pairwise.wilcox.test(X_Fold_Enrichment, X_Annotation, p.adjust.method = "BH")
pairwise.wilcox.test(X_.log10.qvalue., X_Annotation, p.adjust.method = "BH")

library(dplyr)
df2_filtered <- df2 %>%
  group_by(X_Annotation) %>%
  filter(n() > 1)

pairwise.wilcox.test(df2_filtered$X_Fold_Enrichment, df2_filtered$X_Annotation, p.adjust.method = "BH")
pairwise.wilcox.test(df2_filtered$X_.log10.qvalue., df2_filtered$X_Annotation, p.adjust.method = "BH")


detach(df2)


# RAPPRESENTAZIONE GRAFICA

library(viridis)

# Creiamo una palette di colori basata sulle annotazioni (basta una volta!)
annotation_levels <- names(which(table(df2$X_Annotation) >= 3))  # Identiche tra df e df2
color_palette <- setNames(viridis(length(annotation_levels)), annotation_levels)  # Associa colori alle categorie

# Convertiamo annotation in fattore con lo stesso ordine nei due dataset
df$Annotation <- factor(df$Annotation, levels = annotation_levels)
df2$X_Annotation <- factor(df2$X_Annotation, levels = annotation_levels)

# Calcoliamo la numerosità per ogni livello di Annotation
n_HIF1A <- table(df$Annotation)
n_HIF2A <- table(df2$X_Annotation)

#Calcolo le mediane
mediana_HIF1A_FE <- median(df$Fold_Enrichment[df$Annotation == "Distal Intergenic"], na.rm=TRUE)
mediana_HIF2A_FE <- median(df2$X_Fold_Enrichment[df2$X_Annotation == "Distal Intergenic"], na.rm=TRUE)

# Imposta il layout per due boxplot affiancati per il Fold Enrichment
par(mfrow = c(1, 2))  
par(mar = c(10, 5, 4, 2))  # Aumenta il margine inferiore per evitare il taglio delle etichette

# Boxplot per Fold Enrichment - HIF2A
boxplot_obj1 <- boxplot(Fold_Enrichment ~ Annotation, data = df, outline = FALSE, 
                        col = color_palette[levels(df$Annotation)], las = 2, cex.axis = 0.8,
                        main = "Fold Enrichment - HIF2A", ylab = "Fold Enrichment", xlab = "")

# Posizionare la numerosità sulla cima del whisker superiore
text(x = 1:length(n_HIF1A), y = boxplot_obj1$stats[5, ] + 0.9, 
     labels = paste0("n=", n_HIF1A), col = "black", cex = 0.8, font = 2)
# Linea tratteggiata alla mediana di Distal Intergenic
abline(h = mediana_HIF1A_FE, col = "red", lty = 2, lwd = 2)

# Boxplot per Fold Enrichment - HIF1A
boxplot_obj2 <- boxplot(X_Fold_Enrichment ~ X_Annotation, data = df2, outline = FALSE, 
                        col = color_palette[levels(df2$X_Annotation)], las = 2, cex.axis = 0.8,
                        main = "Fold Enrichment - HIF1A", ylab = "Fold Enrichment", xlab = "")

# Posizionare la numerosità sulla cima del whisker superiore
text(x = 1:length(n_HIF2A), y = boxplot_obj2$stats[5, ] + 1.1, 
     labels = paste0("n=", n_HIF2A), col = "black", cex = 0.8, font = 2)
# Linea tratteggiata alla mediana di Distal Intergenic
abline(h = mediana_HIF2A_FE, col = "red", lty = 2, lwd = 2)

# Aggiungiamo la legenda unica (con colori corretti)
# legend("topright", legend = annotation_levels, fill = color_palette, bty = "n", cex = 0.8)





# Imposta una nuova finestra grafica per i -log10(q-value)
par(mfrow = c(1, 2))  
par(mar = c(10, 5, 4, 2))  # Mantiene i margini per evitare il taglio delle etichette

# Boxplot per -log10(q-value) - HIF1A
boxplot(X.log10.qvalue. ~ annotation, data = df, outline = FALSE, 
        col = color_palette[levels(df$annotation)], las = 2, cex.axis = 0.8,
        main = "-log10(q-value) - HIF1A", ylab = "-log10(q-value)", xlab = "")
abline(h = mediana_HIF1A_Q, col = "red", lty = 2, lwd = 2)  # Linea tratteggiata alla mediana di Distal Intergenic

# Boxplot per -log10(q-value) - HIF2A
boxplot(X_.log10.qvalue. ~ X_annotation, data = df2, outline = FALSE, 
        col = color_palette[levels(df2$X_annotation)], las = 2, cex.axis = 0.8,
        main = "-log10(q-value) - HIF2A", ylab = "-log10(q-value)", xlab = "")
abline(h = mediana_HIF2A_Q, col = "red", lty = 2, lwd = 2)  # Linea tratteggiata alla mediana di Distal Intergenic

# Aggiungiamo di nuovo la stessa legenda per i -log10(q-value)
legend("topright", legend = annotation_levels, fill = color_palette, bty = "n", cex = 0.8)

