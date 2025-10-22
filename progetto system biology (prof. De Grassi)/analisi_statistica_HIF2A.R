setwd("C:/Users/carme/Desktop/Bioinf/3° semestre/Systems Biology/progetto esame/")
annotation_mod <- read.table("Modified_Galaxy116_Annotated_Peaks.tabular", header = TRUE, sep = "\t", quote = "", fill = TRUE, comment.char = "")

summary(annotation_mod)
attach(annotation_mod)
hist(distanceToTSS)
hist(distanceToTSS[distanceToTSS < 1000 & distanceToTSS > -1000])
length(distanceToTSS[distanceToTSS < 1000 & distanceToTSS > -1000])

summary(abs(distanceToTSS))

score_close_TSS = score[abs(distanceToTSS) < 1108]
length(score_close_TSS)

score_far_TSS = score[abs(distanceToTSS) > 66074]
length(score_far_TSS)

score_middle_TSS = score[abs(distanceToTSS)>=1108 & abs(distanceToTSS)<=66074]
length(score_middle_TSS)

summary(score_close_TSS)
summary(score_far_TSS)

# TEST DI NORMALITA'

shapiro.test(score_close_TSS) # p-value significativo mi permette di rigettare H0, quindi la distribuzione NON è normale
shapiro.test(score_far_TSS) # p-value significativo mi permette di rigettare H0, quindi la distribuzione NON è normale

ks.test(unique(score_close_TSS), "ppois", lambda = mean(score_close_TSS)) # p-value significativo mi permette di rigettare H0, quindi i dati non seguono distribuzione di Poisson

# ANOVA (invoco il TEOREMA DEL LIMITE CENTRALE, perchè distribuzione non è normale)

score_col = c(score_close_TSS, score_middle_TSS, score_far_TSS) # unisco gli score in un unico vettore
length(score_col)
dist_col = rep(c('close', 'mid', 'far'), times=c(length(score_close_TSS), length(score_middle_TSS), length(score_far_TSS))) #creo un vettore con le acategorie associate alle distanze
length(dist_col)

score_aov = aov(score_col~dist_col)
summary(score_aov) # p-value significativo, quindi c'è differenza nei confronti a coppie, ma devo capire quale specifico confronto incotra una differenza

TukeyHSD(score_aov) # test "post-hoc"
# in questo caso, il p-value è significativo nei confronti a coppie in cui compare "far". Di conseguenza questo fattore lega maggiormente porzioni distanti dal TSS
par(mfrow=c(2,2))
plot(score_aov)

# TEST DI KRUSKAL WALLIS

kruskal.test(score~annotation) #p-value significativo

annotation <- ifelse(grepl("Intron", annotation), "Intron", annotation)
annotation <- ifelse(grepl("Exon", annotation), "Exon", annotation)

pairwise.wilcox.test(score, annotation, p.adjust.method = "BH") # i p-value che identifico sono indicativi delle coppie i cui score sono utili per discriminare l'appertenenza a uno dei due gruppi

# provo l'ANOVA usando tutte le categorie di annotation e gli score del dataset di partenza
score_aov_2 = aov(score~annotation)
summary(score_aov_2)
TukeyHSD(score_aov_2)

mediana_distal <- median(annotation_mod$score[annotation_mod$annotation == "Distal Intergenic"], na.rm = TRUE)

library(viridisLite)
par(mar=c(10,5,4,2))
boxplot(score ~ annotation, outline = FALSE, 
        col = viridis(length(unique(annotation))), las = 2, cex.axis = 0.8, xlab = "", main = "score - HIF2A")
abline(h = mediana_distal, col = "red", lty = 2, lwd = 2)

detach(annotation_mod)
