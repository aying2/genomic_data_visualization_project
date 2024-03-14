data <-
    read.csv("genomic_data_visualization_project/xenium_breast_cancer1.csv.gz",
             row.names = 1)

data[1:10, 1:10]

pos <- data[, 1:2]
area <- data[, 3]
gexp <- data[, 4:ncol(data)]

radius <- sqrt(area / pi)

volume <- 4 / 3 * radius ^ 3

hist(area)
hist(volume)
hist(colSums(gexp))
hist(rowSums(gexp))

plot(area, rowSums(gexp))
plot(volume, rowSums(gexp))

hist(log10(colSums(gexp) + 1))
hist(log10(rowSums(gexp) + 1))

gexpnorm <- log10(gexp/area * mean(area)+1)
# gexpnorm <- log10(gexp / volume * mean(volume) + 1)
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
plot(pcs$sdev[1:50], type = "o")

library(Rtsne)
? Rtsne
set.seed(39)
emb <- Rtsne(gexpnorm)$Y

tw <- sapply(1:12, function(i) {
    print(i)
    set.seed(39)
    kmeans(gexpnorm, centers = i, iter.max = 1000)$tot.withinss
})
plot(tw, type = 'o')

tw_df <- data.frame(tw)
tw_df$k = as.numeric(rownames(tw_df))

ggplot(tw_df, aes(x = k, y = tw, grouping = 1)) +
    geom_line() + geom_point(color = "cornflowerblue", size = 2) +
    labs(title =
             'total withinness vs. k', ) +
    ylab("total withinness") +
    scale_x_continuous(breaks = seq(0, 15, 2)) +
    theme_bw()

library(cluster)
? silhouette
sis <- list()
for (i in 2:15) {
    print(i)
    set.seed(39)
    sis[[i]] <-
        silhouette(kmeans(gexpnorm, centers = i)$cluster, dist(gexpnorm))
}
plot(2:15,mean_sis, type = 'o')

? kmeans
set.seed(39)
nclusters = 8
com <-
    as.factor(kmeans(gexpnorm, centers = nclusters, iter.max = 1000)$cluster)

ggplot(data.frame(pos, gexpnorm, com)) +
    geom_point(aes(x = x_centroid, y = y_centroid, col = com), size = 0.8) +
    labs(title =
             'cluster vs. y_centroid vs. x_centroid (physical space)', ) +
    labs(col = "cluster") +
    theme_bw()

ggplot(data.frame(emb, gexpnorm, com)) +
    geom_point(aes(x = X1, y = X2, col = com), size = 0.8) +
    labs(title =
             'cluster vs. X2 vs. X1 (tSNE space)', ) +
    labs(col = "cluster") +
    theme_bw()


plot_clustersofinterest <-
    function(df, clustersofinterest, nclusters) {
        stopifnot(length
                  (clustersofinterest) < nclusters &&
                      length(clustersofinterest) >= 1)
        df2 <- df
        levels(df2$com) <- c(levels(df2$com), nclusters + 1)
        df2[df["com"] != clustersofinterest[1], "com"] = nclusters + 1
        
        for (i in 2:length(clustersofinterest)) {
            df2[df["com"] == clustersofinterest[i], "com"] = clustersofinterest[i]
        }
        
        library(scales)
        def_col = hue_pal()(nclusters)
        
        other_lab <- ""
        for (i in 1:nclusters) {
            if (!(i %in% clustersofinterest)) {
                sep <- ", "
                if (nchar(other_lab) == 0) {
                    sep = ""
                }
                other_lab <- paste(other_lab, i, sep = sep)
            }
        }
        
        
        p1 <- ggplot(df2) +
            geom_point(aes(x = x_centroid, y = y_centroid, col = com), size =
                           0.8) +
            scale_color_manual(
                values = c(def_col[clustersofinterest], "gray"),
                labels = c(clustersofinterest, other_lab)
            ) +
            labs(col = 'cluster') +
            labs(title = 'cluster vs. y_centroid vs. x_centroid (physical space)') +
            
            theme_bw()
        
        p2 <- ggplot(df2) +
            geom_point(aes(x = X1, y = X2, col = com), size =
                           0.8) +
            scale_color_manual(
                values = c(def_col[clustersofinterest], "gray"),
                labels = c(clustersofinterest, other_lab)
            ) +
            labs(col = 'cluster') +
            labs(title = 'cluster vs. X2 vs. X1 (tSNE space)') +
            
            theme_bw()
        
        list(p1, p2)
    }

clustersofinterest <- c(1, 3, 5, 7)
df <- data.frame(pos, emb, com)

plot_clustersofinterest(df, clustersofinterest, nclusters)[2]

plot_volcano <- function(clusterofinterest) {
    pv <- sapply(colnames(gexpnorm), function(i) {
        # print(i) ## print out gene name
        wilcox.test(gexpnorm[com == clusterofinterest, i],
                    gexpnorm[com != clusterofinterest, i])$p.val
    })
    logfc <- sapply(colnames(gexpnorm), function(i) {
        # print(i) ## print out gene name
        log2(mean(gexpnorm[com == clusterofinterest, i])
             / mean(gexpnorm[com != clusterofinterest, i]))
    })
    
    df4 <- data.frame(pv = pv, logpv = -log10(pv + 1e-300), logfc)
    df4["gene"] <- rownames(df4)
    df4["diffexp"] = "Not Significant"
    df4[df4["pv"] < 0.01 &
            df4["logfc"] > .58, "diffexp"] = "Upregulated"
    df4[df4["pv"] < 0.01 &
            df4["logfc"] < -.58, "diffexp"] = "Downregulated"
    df4$diffexp = as.factor(df4$diffexp)
    
    # print(levels(df4$diffexp))
    
    diff_col <- c()
    for (lev in levels(df4$diffexp)) {
        if (lev == "Downregulated") {
            diff_col <- c(diff_col, "cornflowerblue")
        }
        else if (lev == "Not Significant") {
            diff_col <- c(diff_col, "grey")
        }
        else if (lev == "Upregulated") {
            diff_col <- c(diff_col, "firebrick")
        }
    }
    
    # print(diff_col)
    
    # pvup <- df4[df4["diffexp"] == "Upregulated", "logpv", drop=FALSE]
    # pvup_top <- pvup[order(pvup$logpv, decreasing = TRUE),"logpv",drop=FALSE]
    # pvup_thresh <- pvup_top[5,]
    # pvup_names <- rownames(pvup_top[1:5,"logpv",drop=FALSE])
    # fcup_thresh <- sort(df4[pvup_names, "logfc"], decreasing = TRUE)[5]
    # 
    # print(pvup_names)
    # print(df4[pvup_names, "logfc"])
    # print(pvup_thresh)
    # print(fcup_thresh)
    # 
    # pvdown <- df4[df4["diffexp"] == "Downregulated", "logpv", drop=FALSE]
    # pvdown_top <- pvdown[order(pvdown$logpv, decreasing = TRUE),"logpv",drop=FALSE]
    # pvdown_thresh <- pvdown_top[5,]
    # pvdown_names <- rownames(pvdown_top[1:5,"logpv",drop=FALSE])
    # fcdown_thresh <- sort(df4[pvdown_names, "logfc"], decreasing = TRUE)[5]
    # 
    # print(pvdown_thresh)
    # print(fcdown_thresh)
    
    library(ggrepel)
    ggplot(df4) + geom_point(aes(x = logfc, y = logpv, color = diffexp)) +
        scale_color_manual(values = diff_col) +
        geom_text_repel(
            aes(
                label = ifelse(FALSE
                               , as.character(gene), "")
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
        theme_bw()
    # + theme(plot.title = element_text(size = 10))
}

plot_volcano(5)

volcano_plts <- list()
for (i in clustersofinterest) {
    print(i)
    volcano_plts[[i]] <- plot_volcano(i)
}

volcano_plts[[5]]
