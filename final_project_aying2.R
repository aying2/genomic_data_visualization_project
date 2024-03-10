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

gexpnorm <- log10(gexp/area * mean(area)+1)
# gexpnorm <- gexp/area * mean(area)

# gexpnorm <- log10(gexp/rowSums(gexp) * mean(rowSums(gexp))+1)

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

ggplot(data.frame(pos, gexpnorm)) + 
    geom_point(aes(x = x_centroid, y = y_centroid, col = CD4)) +
    scale_color_gradient(high = "darkred", low = "gray")

pcs <- prcomp(gexpnorm)
plot(pcs$sdev[1:50], type="o")

library(Rtsne)
?Rtsne
set.seed(39)
emb <- Rtsne(gexpnorm)$Y

set.seed(39)
tw <- sapply(1:12, function(i) {
    print(i)
    kmeans(gexpnorm, centers=i, iter.max = 1000)$tot.withinss
})
plot(tw, type='o')

tw_df <- data.frame(tw)
tw_df$k = as.numeric(rownames(tw_df))

ggplot(tw_df, aes(x = k, y = tw, grouping = 1)) + 
    geom_point(color="blue2", size= 2) + geom_line() + 
    labs(
        title =
            'total withinness vs. k',
    ) +
    ylab("total withinness")+
    theme_bw()

?kmeans
set.seed(39)
nclusters = 8
com <- as.factor(kmeans(gexpnorm, centers=nclusters, iter.max = 1000)$cluster)

ggplot(data.frame(pos, gexpnorm, com)) + 
    geom_point(aes(x = x_centroid, y = y_centroid, col = com), size = 0.8) + 
    labs(
        title =
            'cluster vs. y_centroid vs. x_centroid (physical space)',
    ) +
    labs(col = "cluster")+
    theme_bw()

clustersofinterest <- c(2, 3, 5, 6)

df <- data.frame(pos, emb, com)
df2 <- df
df2[df["com"]!=clustersofinterest[1], "com"] = 8
df2[df["com"]==clustersofinterest[2], "com"] = clustersofinterest[2]
df2[df["com"]==clustersofinterest[3], "com"] = clustersofinterest[3]
df2[df["com"]==clustersofinterest[4], "com"] = clustersofinterest[4]

def_col = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67",
            "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")
ggplot(df2) + 
    geom_point(aes(x = x_centroid, y = y_centroid, col = com), size=0.8) +
    scale_color_manual(values = c(def_col[clustersofinterest], "gray"), 
                       labels = c("1, 3, 4, 5, 7", "2", "6"))+
    labs(col='cluster') +
    labs(title = 'cluster vs. y_centroid vs. x_centroid (physical space)') +
    
    theme_bw()

ggplot(data.frame(emb, gexpnorm, com)) + 
    geom_point(aes(x = X1, y = X2, col = com), size = 0.8) +
    labs(
        title =
            'cluster vs. X2 vs. X1 (tSNE space)',
    ) +
    labs(col = "cluster")+
    theme_bw()

com_idx = 1
df3 <- data.frame(pos, emb, com)
df3[df["com"]!=clustersofinterest[com_idx], "com"] = 8
ggplot(df3) + 
    geom_point(aes(x = x_centroid, y = y_centroid, col = com), size=0.8) +
    scale_color_manual(values = c(def_col[clustersofinterest[com_idx]], "gray"), 
                       labels = c("other", "6"))+
    labs(col='cluster') +
    labs(title = 'cluster vs. y_centroid vs. x_centroid (physical space)') +
    
    theme_bw()

volcano_plts <- list()
for (x in 1:nclusters) {
    print(x)
    clusterofinterest <- x
    
    pv <- sapply(colnames(gexpnorm), function(i) {
        print(i) ## print out gene name
        wilcox.test(gexpnorm[com == clusterofinterest, i],
                    gexpnorm[com != clusterofinterest, i])$p.val
    })
    logfc <- sapply(colnames(gexpnorm), function(i) {
        print(i) ## print out gene name
        log2(mean(gexpnorm[com == clusterofinterest, i])
             / mean(gexpnorm[com != clusterofinterest, i]))
    })
    
    df4 <- data.frame(pv = pv, logpv = -log10(pv + 1e-300), logfc)
    df4["gene"] <- rownames(df4)
    df4["diffexp"] = "Not Significant"
    df4[df4["pv"] < 0.01 &
            df4["logfc"] > 1, "diffexp"] = "Upregulated"
    df4[df4["pv"] < 0.01 &
            df4["logfc"] < -1, "diffexp"] = "Downregulated"
    df4$diffexp = as.factor(df4$diffexp)
    
    library(ggrepel)
    volcano_plts[[x]] <-
        ggplot(df4) + geom_point(aes(x = logfc, y = logpv, color = diffexp)) +
        scale_color_manual(values = c("cornflowerblue", "grey", "firebrick")) +
        geom_text_repel(
            aes(
                label = ifelse(diffexp == "chicken", as.character(gene), "")
                ,
                x = logfc,
                y = logpv
            ),
            size = 3,
            max.overlaps = Inf,
            box.padding = .5
        ) +
        labs(
            title = sprintf(
                'Cluster %d Differential Expression vs. -log10(pv + 1e-300) vs. log2(fc)',
                clusterofinterest
            )
        ) +
        ylab("-log10(pv)") +
        xlab("log2(fc)") +
        theme_bw() +
        theme(plot.title = element_text(size = 10))
}
volcano_plts[[3]]
