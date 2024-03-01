data <-
    read.csv("genomic_data_visualization_project/xenium_breast_cancer1.csv.gz",
             row.names = 1)

data[1:10, 1:10]

pos <- data[, 1:2]
area <- data[, 3]
gexp <- data[, 4:ncol(data)]

hist(area)
hist(colSums(gexp))
hist(rowSums(gexp))

plot(area, rowSums(gexp))

hist(log10(colSums(gexp) + 1))
hist(log10(rowSums(gexp) + 1))

# gexpnorm <- log10(gexp/area * mean(area)+1)

gexpnorm <- log10(gexp/rowSums(gexp) * mean(rowSums(gexp))+1)

# gexpnorm2 <- gexp/area * mean(area)
# gexpnorm <- log10(gexpnorm2/rowSums(gexpnorm2) * mean(rowSums(gexpnorm2)) + 1)

hist(rowSums(gexpnorm))

library(ggplot2)
ggplot(data.frame(pos, area)) + 
    geom_point(aes(x = x_centroid, y = y_centroid, col = area), size = 1) +
    scale_color_gradient(high = "darkred", low = "gray")

ggplot(data.frame(pos, gexpnorm)) + 
    geom_point(aes(x = x_centroid, y = y_centroid, col = CD4)) +
    scale_color_gradient(high = "darkred", low = "gray")

pcs <- prcomp(gexpnorm)
plot(pcs$sdev[1:50], type="o")

library(Rtsne)
?Rtsne
set.seed(42)
emb <- Rtsne(gexpnorm, perplexity = 20)$Y

set.seed(17)
tw <- sapply(1:12, function(i) {
    print(i)
    kmeans(gexpnorm, centers=i, iter.max = 1000)$tot.withinss
})
plot(tw, type='o')

?kmeans
set.seed(17)
com <- as.factor(kmeans(gexpnorm, centers=10, iter.max = 1000)$cluster)

ggplot(data.frame(pos, gexpnorm, com)) + 
    geom_point(aes(x = x_centroid, y = y_centroid, col = com), size = 1)

ggplot(data.frame(pcs$x, gexpnorm, com)) + 
    geom_point(aes(x = PC1, y = PC2, col = com), size = 1)

ggplot(data.frame(emb, gexpnorm, com)) + 
    geom_point(aes(x = X1, y = X2, col = com), size = 1)

