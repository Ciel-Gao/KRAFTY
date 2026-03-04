# This document provides analysis for world trade data for 2010 and 2023.
# Use the data of raw chicken meat product.
# The data is sourced from the Food and Agriculture Organization (FAO) of the United Nations.

# Load necessary libraries
library(tidyverse)  
library(igraph)    
library(ggmap)     
library(data.table) 
library(Matrix)    
library(gridExtra) 
library(maps)      # Geographic mapping
library(geosphere) # Geographic calculations
library(irlba)     
library(mclust)    
library(MGLM)      
library(GGally)    
library(pracma)    
library(MASS)      
library(ggpubr)    
library(rnaturalearth)  # For world map plots
library(rnaturalearthdata)  # For world map plots
library(sf)
library(ggalluvial)
library(RColorBrewer)
library(openxlsx)


# Helper functions
class2mat <- function(tauvec){
    tauvec<- as.numeric(tauvec)
    k <- length(unique(tauvec))
    taumat <- matrix(rep(0,length(tauvec)*k),nrow=length(tauvec))
    for(i in 1:length(tauvec)){
        taumat[i,tauvec[i]]<- 1
    }
    return(taumat)
}

A_asym <- function(df, df_trade_total, product=NULL){
    df <- copy(df) 
    setDT(df)
    
    if (is.null(product)) {
        df_filtered <- df
    }else{
        df_filtered <- df %>% 
            filter(Item == product, Value > 0) 
    }
    
    countries <- unique(c(df_trade_total$`Export Countries`, df_trade_total$`Import Countries`))
    country_idx <- setNames(seq_along(countries), countries)
    
    
    A_matrix <- sparseMatrix(
        i = country_idx[df_filtered$`Export Countries`],
        j = country_idx[df_filtered$`Import Countries`],
        x = df_filtered$Value,  
        dims = c(length(countries), length(countries)),
        dimnames = list(countries, countries)
    )
    
    A <- as.matrix(A_matrix)
    
    # compute exports (row sums) and imports (col sums)
    exports <- rowSums(A, na.rm = TRUE)
    imports <- colSums(A, na.rm = TRUE)
    total_trade <- exports + imports  # per country total
    
    # first quartile threshold
    # q1 <- quantile(total_trade, 0.5, na.rm = TRUE)
    
    # keep countries at or above Q1
    keep <- names(total_trade)[total_trade > 0]
    # keep <- names(total_trade)[total_trade > q1]
    
    # drop countries below Q1 from both axes
    A_kept <- A[keep, keep, drop = FALSE]
    return(A_kept)
}

normalize_A <- function(X) {
    row_norms <- sqrt(rowSums(X^2))  
    row_norms[row_norms == 0] <- 1  
    X_normalized <- X / row_norms    
    return(X_normalized)
}

extract_latent_positions <- function(trade_matrix, k, role = "export") {
    
    #countries <- rownames(trade_matrix)
    
    trade_svd <- irlba(trade_matrix, nu=k, nv = k, maxit = 1000)
    
    u <- trade_svd$u
    s <- diag(sqrt(trade_svd$d[1:k]))
    v <- trade_svd$v
    
    if (role == "export") {
        latent_positions <- u %*% s
    } else if (role == "import") {
        latent_positions <- v %*% s
    } else {
        stop("Invalid role. Choose either 'export' or 'import'.")
    }
    
    
    return(latent_positions) # this is a matrix format, not data frame format
}

export_and_import <- function(A, k){
    US <- extract_latent_positions(A, k, role = "export")
    VS <- extract_latent_positions(A, k, role = "import")
    US_norm <- normalize_A(US)
    VS_norm <- normalize_A(VS)
    return(list(U = US_norm, V = VS_norm))
}


getElbows <- function(dat, n = 3, threshold = FALSE, plot = FALSE, main="") {
    ## Given a decreasingly sorted vector, return the given number of elbows
    ##
    ## Args:
    ##   dat: a input vector (e.g. a vector of standard deviations), or a input feature matrix.
    ##   n: the number of returned elbows.
    ##   threshold: either FALSE or a number. If threshold is a number, then all
    ##   the elements in d that are not larger than the threshold will be ignored.
    ##   plot: logical. When T, it depicts a scree plot with highlighted elbows.
    ##
    ## Return:
    ##   q: a vector of length n.
    ##
    ## Reference:
    ##   Zhu, Mu and Ghodsi, Ali (2006), "Automatic dimensionality selection from
    ##   the scree plot via the use of profile likelihood", Computational
    ##   Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006. 
    
    #  if (is.unsorted(-d))
    
    
    if (is.matrix(dat)) {
        d <- sort(apply(dat,2,sd), decreasing=TRUE)
    } else {
        d <- sort(dat,decreasing=TRUE)
    }
    
    if (!is.logical(threshold))
        d <- d[d > threshold]
    
    p <- length(d)
    if (p == 0)
        stop(paste("d must have elements that are larger than the threshold ",
                   threshold), "!", sep="")
    
    lq <- rep(0.0, p)                     # log likelihood, function of q
    for (q in 1:p) {
        mu1 <- mean(d[1:q])
        mu2 <- mean(d[-(1:q)])              # = NaN when q = p
        sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
            (p - 1 - (q < p))
        lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
            sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
    }
    
    q <- which.max(lq)
    if(is.na(q)){
        q <- p
    }
    if (n > 1 && q < (p-1)) {
        q <- c(q, q + getElbows(d[(q+1):p], n-1, plot=FALSE))
    }
    
    if (plot==TRUE) {
        if (is.matrix(dat)) {
            sdv <- d # apply(dat,2,sd)
            plot(sdv,type="b",xlab="dim",ylab="stdev",main=main)
            points(q,sdv[q],col=2,pch=19)
        } else {
            plot(dat, type="b",xlab="Index", ylab = "Singular Values", main=main)
            p + points(q,dat[q],col=2,pch=19)
        }
    }
    p <- recordPlot()
    
    return(list(q=q,p=p))
}


procrustes2 <- function(X, Y) {
    tmp <- t(X) %*% Y
    tmp.svd <- svd(tmp)
    W <- tmp.svd$u %*% t(tmp.svd$v)
    newX <- X %*% W
    return(list(newX = newX, error = norm(newX-Y, type="F"), W = W))
}

find_country_differences <- function(matrix1, matrix2) {
    countries1 <- rownames(matrix1)
    countries2 <- rownames(matrix2)
    
    only_in_matrix1 <- setdiff(countries1, countries2)
    only_in_matrix2 <- setdiff(countries2, countries1)
    
    list(
        only_in_matrix1 = only_in_matrix1,
        only_in_matrix2 = only_in_matrix2
    )
}

filter_common_countries <- function(matrix1, matrix2) {
    common_countries <- intersect(rownames(matrix1), rownames(matrix2))
    
    matrix1_filtered <- matrix1[common_countries, common_countries, drop = FALSE]
    matrix2_filtered <- matrix2[common_countries, common_countries, drop = FALSE]
    
    list(
        matrix1_filtered = matrix1_filtered,
        matrix2_filtered = matrix2_filtered
    )
}

hierachical_cluster <- function(embeddings,num_clusts){
    eH <- hc(embeddings)
    CR <- c(hclass(eH,num_clusts))
    return(CR)
}

Krafty <- function(M1, M2){
    ## M: input matrices, U or Z
    
    set.seed(123)
    Mkr <- kr(M1, M2, byrow = TRUE)
    Mkr.svd <- svd(Mkr)
    elbMkr_list <-getElbows(Mkr.svd$d, n=1, plot=TRUE)
    elbMkr <- elbMkr_list$q
    plot_elbMkr <- elbMkr_list$p
    K <- elbMkr[length(elbMkr)]
    u <- Mkr.svd$u[,1:K]
    cluster_result <- hierachical_cluster(u, K)
    # cluster_result <- kmeans(u, centers=K)$cluster
    return(list(Mkr = Mkr, Mkr.svd = Mkr.svd, K = K, cluster_result = cluster_result, figure=plot_elbMkr))
}

MASE <- function(M1, M2){
    ## M: input matrices, U or Z
    
    set.seed(123)
    Mmase <- cbind(M1, M2)
    Mmase.svd <- svd(Mmase)
    elbMmase_list <- getElbows(Mmase.svd$d, n=1, plot=TRUE)
    elbMmase <- elbMmase_list$q
    plot_elbMmase <- elbMmase_list$p
    K <- elbMmase[length(elbMmase)]
    u <- Mmase.svd$u[,1:K]
    cluster_result <- hierachical_cluster(u, K)
    # cluster_result <- kmeans(u, centers=K)$cluster
    return(list(Mmase = Mmase, Mmase.svd = Mmase.svd, K = K, cluster_result = cluster_result, figure=plot_elbMmase))
}

# This plot is for comparing joint clusters and clusters from single source.
alluvial_plot <- function(
        data,
        kr_col,
        exporter_col,
        importer_col,
        joint_label = "Joint Cluster",
        pair_label  = "Exporter x Importer",
        type = "1yr"
) {
    
    if (type == "1yr") {
        df2 <- data %>%
            dplyr::select({{ kr_col }}, {{ exporter_col }}, {{ importer_col }}) %>%
            mutate(pair = paste0("E", {{ exporter_col }}, " x I", {{ importer_col }})) %>%
            count({{ kr_col }}, pair, name = "n")
    } else if (type == "exp") {
        df2 <- data %>%
            dplyr::select({{ kr_col }}, {{ exporter_col }}, {{ importer_col }}) %>%
            mutate(pair = paste0("E", {{ exporter_col }}, " x E", {{ importer_col }})) %>%
            count({{ kr_col }}, pair, name = "n")
    } else if (type == "imp") {
        df2 <- data %>%
            dplyr::select({{ kr_col }}, {{ exporter_col }}, {{ importer_col }}) %>%
            mutate(pair = paste0("I", {{ exporter_col }}, " x I", {{ importer_col }})) %>%
            count({{ kr_col }}, pair, name = "n")
    } else {
        stop("Invalid type. Choose either '1yr' or '2yr'.")
    }
    my_cols <- c("#D97D55", "#E6D8C3", "#B8C4A9", "#D9A299", "#8CB9BD")
    
    p <- ggplot(
        df2,
        aes(y = n, axis1 = {{ kr_col }}, axis2 = pair)
    ) +
        geom_alluvium(aes(fill = {{ kr_col }}),
                      width = 1/12, alpha = 0.8, show.legend = FALSE,
                      discern = TRUE) +
        geom_stratum(width = 1/12, fill = "grey80", color = "grey40",
                     show.legend = FALSE, discern = TRUE) +
        geom_text(stat = "stratum",
                  aes(label = after_stat(stratum)),
                  size = 5, discern = TRUE) +
        # scale_fill_brewer(palette = "Blues") +
        scale_fill_manual(values = my_cols) +
        scale_x_discrete(limits = c(joint_label, pair_label), expand = c(.1, .1)) +
        labs(y = "Number of Countries", x = NULL) +
        theme_minimal() +
        theme(
            text = element_text(size = 18),
            axis.text.x  = element_text(size = 16),
            axis.text.y  = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            legend.position = "none"
        )
    
    plot(p)
    return(p)
}

# Confusion heatmap
construct_conf_mat <- function(df1, df, target_var, predicted_var, fill_limits = c(0, 60), comments) {
    
    df[[target_var]] <- factor(df[[target_var]])
    df[[predicted_var]] <- factor(df[[predicted_var]])
    
    target_levels <- levels(df[[target_var]])
    predicted_levels <- levels(df[[predicted_var]])
    
    conf_matrix <- table(df[[target_var]], df[[predicted_var]])
    cm_df <- as.data.frame(as.table(conf_matrix))
    colnames(cm_df) <- c("Target", "Predicted", "Freq")
    
    full_grid <- expand.grid(
        Target = factor(target_levels, levels = target_levels),
        Predicted = factor(predicted_levels, levels = predicted_levels)
    )
    
    cm_df <- left_join(full_grid, cm_df, by = c("Target", "Predicted")) %>%
        mutate(Freq = ifelse(is.na(Freq), 0, Freq))
    
    plot <- ggplot(cm_df, aes(x = Target, y = Predicted, fill = Freq)) +
        geom_tile(color = "white") +
        scale_fill_gradient(low = "white", high = "red", limits = fill_limits) +
        labs(x = target_var, y = predicted_var, title = paste0("Confusion Matrix (", comments, ")")) +
        theme_minimal() +
        theme(axis.text.x = element_text(hjust = 1))
    
    print(plot)
    return(conf_matrix)
}

my_cols <- c("#D97D55", "#E6D8C3", "#B8C4A9", "#D9A299", "#8CB9BD")

# --------------------------Main Analysis--------------------------

## -------------- Part1: Data preprocessing --------------
# Load trade data for 2010 and 2023
trade_2010 <- read_csv("trade_2010_money.csv", show_col_types = FALSE)
trade_2023 <- read_csv("trade_2023_money.csv", show_col_types = FALSE)


# Choose product
chicken_2010 <- trade_2010 %>% filter(Item %in% c("Meat of chickens, fresh or chilled"))
chicken_2023 <- trade_2023 %>% filter(Item %in% c("Meat of chickens, fresh or chilled"))

# Construct asymmetric trade matrices for 2010 and 2023
chicken_asym_2010 <- A_asym(chicken_2010, trade_2010)
chicken_asym_2023 <- A_asym(chicken_2023, trade_2023)

# 2010 data has 196 rows and 2023 has 198 rows, try to find the differences
differences <- find_country_differences(chicken_asym_2010, chicken_asym_2023)
print(differences)


# Filter matrices to keep only common countries
filtered_matrices <- filter_common_countries(chicken_asym_2010, chicken_asym_2023)

chicken_2010_filtered <- filtered_matrices$matrix1_filtered
chicken_2023_filtered <- filtered_matrices$matrix2_filtered

elb_2010_list <- getElbows(svd(chicken_2010_filtered)$d[1:20], n=1, plot=TRUE)
elb_2010 <- elb_2010_list$q
plot_2010 <- elb_2010_list$p

# Save the scree plot for 2010
png("screePlot_2010.png", width = 6, height = 4, units = "in", res = 300)
replayPlot(plot_2010)
dev.off()

elb_2023_list <- getElbows(svd(chicken_2023_filtered)$d[1:20], n=1, plot=TRUE)
elb_2023 <- elb_2023_list$q
plot_2023 <- elb_2023_list$p

# Save the scree plot for 2023
png("screePlot_2023.png", width = 6, height = 4, units = "in", res = 300)
replayPlot(plot_2023)
dev.off()

K_2010 <- elb_2010[length(elb_2010)]
K_2023 <- elb_2023[length(elb_2023)]

# 2010 exporter and importer
extract_2010 <- export_and_import(chicken_2010_filtered, K_2010)
US_2010 <- extract_2010$U
VS_2010 <- extract_2010$V

# 2023 exporter and importer
extract_2023 <- export_and_import(chicken_2023_filtered, K_2023)
US_2023 <- extract_2023$U
VS_2023 <- extract_2023$V

# Z matrices for both datasets
set.seed(123)
# 2010
# Z10_exp_clusters <- hierachical_cluster(US_2010,K_2010)
Z10_exp_clusters <- kmeans(US_2010, centers=K_2010)$cluster
Z10_exp <- class2mat(Z10_exp_clusters)

# Z10_imp_clusters <- hierachical_cluster(VS_2010,K_2010)
Z10_imp_clusters <- kmeans(VS_2010, centers=K_2010)$cluster
Z10_imp <- class2mat(Z10_imp_clusters)

# 2023
# Z23_exp_clusters <- hierachical_cluster(US_2023,K_2023)
Z23_exp_clusters <- kmeans(US_2023, centers=K_2023)$cluster
Z23_exp <- class2mat(Z23_exp_clusters)

Z23_imp_clusters <- kmeans(VS_2023, centers=K_2023)$cluster
# Z23_imp_clusters <- hierachical_cluster(VS_2023,K_2023)
Z23_imp <- class2mat(Z23_imp_clusters)


## -------------- Part 2: Single year analysis -- 2010 --------------
Zkr_1yr_lists <- Krafty(Z10_exp, Z10_imp)
Zkr_1yr <- Zkr_1yr_lists$Mkr
est_K_kr <- Zkr_1yr_lists$K
Zkr_1yr_clusters <- Zkr_1yr_lists$cluster_result
Zkr_1yr_figure <- Zkr_1yr_lists$figure

png("screePlot_kr_1yr.png", width = 6, height = 4, units = "in", res = 300)
replayPlot(Zkr_1yr_figure)
dev.off()


Zmase_1yr_lists <- MASE(Z10_exp, Z10_imp)
Zmase_1yr <- Zmase_1yr_lists$Mmase
est_K_mase <- Zmase_1yr_lists$K
Zmase_1yr_clusters <- Zmase_1yr_lists$cluster_result
Zmase_1yr_figure <- Zmase_1yr_lists$figure

png("screePlot_mase_1yr.png", width = 6, height = 4, units = "in", res = 300)
replayPlot(Zmase_1yr_figure)
dev.off()

ari_1yr <- adjustedRandIndex(Zkr_1yr_clusters, Zmase_1yr_clusters)

table(Z10_exp_clusters, Z10_imp_clusters)
table(Zkr_1yr_clusters, Zmase_1yr_clusters)

df_1yr_clusters <- data.frame(name = rownames(chicken_2010_filtered),
                          kr_cluster = factor(Zkr_1yr_clusters),
                          mase_cluster = factor(Zmase_1yr_clusters),
                          exporter_cluster = factor(Z10_exp_clusters), 
                          importer_cluster = factor(Z10_imp_clusters))


# Alluvial plot for KR clusters vs Exporter & Importer clusters
alluvial_1yr <- alluvial_plot(
    data = df_1yr_clusters,
    kr_col = kr_cluster,
    exporter_col = exporter_cluster,
    importer_col = importer_cluster,
    joint_label = "KRAFTY",
    pair_label  = "Exporter x Importer",
    type = "1yr"
)

ggsave("alluvial_kr_exp_imp_1yr.png", plot = alluvial_1yr, width = 12, height = 6, units = "in", dpi = 300)

# Alluvial plot for MASE vs Exp & Imp Clusters
alluvial_1yr_mase <- alluvial_plot(
    data = df_1yr_clusters,
    kr_col = mase_cluster,
    exporter_col = exporter_cluster,
    importer_col = importer_cluster,
    joint_label = "MASE",
    pair_label  = "Exporter x Importer",
    type = "1yr"
)
ggsave("alluvial_mase_exp_imp_1yr.png", plot = alluvial_1yr, width = 12, height = 6, units = "in", dpi = 300)


# Alluvial plot for MASE clusters vs KR clusters
df2 <- df_1yr_clusters %>%
    dplyr::select(kr_cluster, mase_cluster) %>%
    count(kr_cluster, mase_cluster, name = "n")

alluvial_kr_mase_1yr <- ggplot(df2, aes(y = n, axis1 = kr_cluster, axis2 = mase_cluster)) +
    geom_alluvium(aes(fill = kr_cluster), width = 1/12, alpha = 0.8,
                  show.legend = FALSE, discern = TRUE) +
    geom_stratum(width = 1/12, fill = "grey80", color = "grey40",
                 show.legend = FALSE, discern = TRUE) +
    geom_text(stat = "stratum",
              aes(label = after_stat(sub("\\.\\d+$", "", stratum))),
              size = 5, discern = TRUE) +
    scale_fill_manual(values = my_cols) +
    scale_x_discrete(limits = c("KRAFTY", "MASE"), expand = c(.1, .1)) +
    labs(y = "Number of Countries", x = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.position = "none")
alluvial_kr_mase_1yr
ggsave("alluvial_kr_mase_1yr.png", plot = alluvial_kr_mase_1yr, width = 12, height = 6, units = "in", dpi = 300)


# Load world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Convert names to lowercase for comparison
df_1yr_clusters$name <- tolower(df_1yr_clusters$name)
world$name_long <- tolower(world$name_long)

# Check which country names in df_kr_cluster_labels are NOT in the world dataset
unmatched_names_kr <- setdiff(df_1yr_clusters$name, world$name_long)
unmatched_names_wd <- setdiff(world$name_long, df_1yr_clusters$name)

df_1yr_clusters$name <- recode(df_1yr_clusters$name,
                           "bolivia (plurinational state of)" = "bolivia",
                           "iran (islamic republic of)" = "iran",
                           "china, mainland" = "china",
                           "czechia" = "czech republic",
                           "china, hong kong sar" = "hong kong",
                           "netherlands (kingdom of the)" = "netherlands",
                           "gambia" = "the gambia",
                           "syrian arab republic" = "syria",
                           "congo" = "republic of the congo",
                           "united republic of tanzania" = "tanzania",
                           "türkiye" = "turkey",
                           "united states of america" = "united states",
                           "united kingdom of great britain and northern ireland" = "united kingdom",
                           "china, taiwan province of" = "taiwan",
                           "republic of moldova" = "moldova",
                           "venezuela (bolivarian republic of)" = "venezuela",
                           "cabo verde" = "republic of cabo verde",
                           "democratic people's republic of korea" = "dem. rep. korea",
                           "eswatini" = "kingdom of eswatini",
                           "viet nam" = "vietnam",
                           "china, macao sar" = "macao",
                           "lao people's democratic republic" = "lao pdr",
                           "faroe islands" = "faeroe islands",
                           "sudan (former)" = "sudan",
                           "sao tome and principe" = "são tomé and principe",
                           "micronesia (federated states of)" = "federated states of micronesia"
)

setdiff(df_1yr_clusters$name, world$name_long)

write.xlsx(df_1yr_clusters, file = "clusters_1yr.xlsx")


# world_clustered_kr <- left_join(world, df_1yr_clusters, by = c("name_long" = "name"))
# 
# p1 <- ggplot(world_clustered_kr) +
#     geom_sf(aes(fill = factor(kr_cluster)), color = "black", size = 0.2) +  
#     scale_fill_brewer(palette = "RdYlBu", na.value = "gray80") +  
#     theme_minimal() +
#     ggtitle("KR: Asymmetric") +
#     theme(legend.title = element_text(size = 12),
#           legend.text = element_text(size = 10))
# 
# 
# p1
# 
# p3 <- ggplot(world_clustered_kr) +
#     geom_sf(aes(fill = factor(exporter_cluster)), color = "black", size = 0.2) +  
#     scale_fill_brewer(palette = "RdYlBu", na.value = "gray80") +  
#     theme_minimal() +
#     ggtitle("Exp: Asymmetric") +
#     theme(legend.title = element_text(size = 12),
#           legend.text = element_text(size = 10))
# 
# 
# p3
# 
# p4 <- ggplot(world_clustered_kr) +
#     geom_sf(aes(fill = factor(importer_cluster)), color = "black", size = 0.2) +  
#     scale_fill_brewer(palette = "RdYlBu", na.value = "gray80") +  
#     theme_minimal() +
#     ggtitle("Imp: Asymmetric") +
#     theme(legend.title = element_text(size = 12),
#           legend.text = element_text(size = 10))
# 
# 
# p4
# 
# p2 <- ggplot(world_clustered_kr) +
#     geom_sf(aes(fill = factor(mase_cluster)), color = "black", size = 0.2) +  
#     scale_fill_brewer(palette = "RdYlBu", na.value = "gray80") +  
#     theme_minimal() +
#     ggtitle("MASE: Asymmetric") +
#     theme(legend.title = element_text(size = 12),
#           legend.text = element_text(size = 10))
# 
# 
# p2




df_check <- df_1yr_clusters %>% mutate(name = rownames(chicken_2010_filtered),
                                       exporter_value = rowSums(chicken_2010_filtered),
                                       importer_value = colSums(chicken_2010_filtered),
                                       add_value = exporter_value + importer_value)

plot_trade_heatmap <- function(A, df_clusters = NULL, log_scale = TRUE, order_by_clusters = FALSE) {
    # deps: dplyr, tidyr, tibble, ggplot2
    A <- as.matrix(A)
    A <- A + t(A)  # symmetrize
    if (is.null(rownames(A)) || is.null(colnames(A))) {
        stop("A must have rownames and colnames that are country names.")
    }
    
    # Keep these for boundaries if ordering is requested
    row_bound_pos <- col_bound_pos <- numeric(0)
    
    if (order_by_clusters) {
        if (is.null(df_clusters) || !all(c("name", "kr_cluster") %in% names(df_clusters))) {
            stop("order_by_clusters=TRUE requires df_clusters with columns: name, kr_cluster.")
        }
        
        # Factor
        levs <- sort(unique(df_clusters$kr_cluster))
        df_clusters <- df_clusters |>
            dplyr::mutate(kr_cluster = factor(kr_cluster, levels = levs))
        
        # Align sets (keep only names present everywhere)
        keep <- Reduce(intersect, list(rownames(A), colnames(A), df_clusters$name))
        if (length(keep) < 2L) stop("After aligning names, fewer than 2 countries remain. Check name matching.")
        A <- A[keep, keep, drop = FALSE]
        
        # Cluster labels in current matrix order
        kmap <- stats::setNames(df_clusters$kr_cluster, df_clusters$name)
        exp_lab0 <- kmap[rownames(A)]
        imp_lab0 <- kmap[colnames(A)]
        
        if (any(is.na(exp_lab0)) || any(is.na(imp_lab0))) {
            stop("Some countries in A are missing cluster labels in df_clusters.")
        }
        
        # Reorder rows/cols by (cluster, then name) and keep the order indices
        row_ord <- order(exp_lab0, rownames(A))
        col_ord <- order(imp_lab0, colnames(A))
        A <- A[row_ord, col_ord, drop = FALSE]
        
        exp_lab <- exp_lab0[row_ord]
        imp_lab <- imp_lab0[col_ord]
        
        # ---- Compute boundary positions between cluster blocks ----
        # positions are at i + 0.5 where cluster changes between i and i+1
        if (length(exp_lab) > 1) {
            row_changes <- which(diff(as.integer(exp_lab)) != 0)
            row_bound_pos <- row_changes + 0.5
        }
        if (length(imp_lab) > 1) {
            col_changes <- which(diff(as.integer(imp_lab)) != 0)
            col_bound_pos <- col_changes + 0.5
        }
    } else {
        # If not ordering, still fix the axis levels later based on current row/col names
        keep <- intersect(rownames(A), colnames(A))
        A <- A[keep, keep, drop = FALSE]
    }
    
    # Long format for ggplot; lock axis order to current row/col order
    df_heat <- A |>
        as.data.frame(check.names = FALSE) |>
        tibble::rownames_to_column(var = "exporter") |>
        tidyr::pivot_longer(cols = -exporter, names_to = "importer", values_to = "value") |>
        dplyr::mutate(
            value_plot = if (log_scale) log10(1 + value) else value,
            exporter   = factor(exporter, levels = rownames(A)),
            importer   = factor(importer, levels = colnames(A))
        )
    
    # Base heatmap
    p <- ggplot2::ggplot(df_heat, ggplot2::aes(x = importer, y = exporter, fill = value_plot)) +
        ggplot2::geom_tile(width = 1.05, height = 1.05) +  # <— bigger cells
        ggplot2::coord_fixed() +
        ggplot2::labs(
            title = "Trade heatmap",
            x = "Countries (cols)", y = "Countries (rows)",
            fill = if (log_scale) "log10(1 + value)" else "Volume"
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
            panel.grid  = ggplot2::element_blank()
        ) +
        ggplot2::scale_fill_gradient(low = "white", high = "black", na.value = "grey80")
    
    if (length(row_bound_pos)) { p <- p + ggplot2::geom_hline(yintercept = row_bound_pos, color = "red", linewidth = 0.4) } 
    if (length(col_bound_pos)) { p <- p + ggplot2::geom_vline(xintercept = col_bound_pos, color = "red", linewidth = 0.4) }
    
    
    print(p)
    invisible(df_heat)
    return(p)
}

hp <- plot_trade_heatmap(chicken_2010_filtered, df_clusters = df_check, log_scale = TRUE, order_by_clusters = TRUE)
ggsave("trade_heatmap_1yr.png", plot = hp, width = 20, height = 20, units = "in", dpi = 300)

# write.xlsx(df_check, file = "clusters_1yr_with_values.xlsx")

## -------------- Part 3: Multiple year analysis -- 2010 vs 2023 --------------
## -------------- kr and mase: exporters --------------
set.seed(123)
Zkr_2yr_exp_lists <- Krafty(Z10_exp, Z23_exp)
Zkr_2yr_exp <- Zkr_2yr_exp_lists$Mkr
K_kr_2yr_exp <- Zkr_2yr_exp_lists$K
Zkr_2yr_exp_clusters <- Zkr_2yr_exp_lists$cluster_result
Zkr_2yr_exp_figure <- Zkr_2yr_exp_lists$figure

# save the scree plot
png("screePlot_kr_2yr_exp.png", width = 6, height = 4, units = "in", res = 300)
replayPlot(Zkr_2yr_exp_figure)
dev.off()

Zmase_2yr_exp_lists <- MASE(Z10_exp, Z23_exp)
Zmase_2yr_exp <- Zmase_2yr_exp_lists$Mmase
K_mase_2yr_exp <- Zmase_2yr_exp_lists$K
Zmase_2yr_exp_clusters <- Zmase_2yr_exp_lists$cluster_result
Zmase_2yr_exp_figure <- Zmase_2yr_exp_lists$figure

# save the scree plot
png("screePlot_mase_2yr_exp.png", width = 6, height = 4, units = "in", res = 300)
replayPlot(Zmase_2yr_exp_figure)
dev.off()


ari_exp <- adjustedRandIndex(Zkr_2yr_exp_clusters, Zmase_2yr_exp_clusters)

df_2yr_exp_clusters <- data.frame(name = df_1yr_clusters$name,
                          kr_cluster = factor(Zkr_2yr_exp_clusters),
                          mase_cluster = factor(Zmase_2yr_exp_clusters),
                          exporter_10 = factor(Z10_exp_clusters), 
                          exporter_23 = factor(Z23_exp_clusters))

write.xlsx(df_2yr_exp_clusters, file = "clusters_2yr_exp.xlsx")

# Alluvial plot for KR clusters vs Exporter clusters in 2010 & 2023
alluvial_2yr_exp <- alluvial_plot(
    data = df_2yr_exp_clusters,
    kr_col = kr_cluster,
    exporter_col = exporter_10,
    importer_col = exporter_23,
    joint_label = "KRAFTY",
    pair_label  = "Exporter 2010 x Exporter 2023",
    type = "exp"
)

ggsave("alluvial_kr_exp_2yr.png", plot = alluvial_2yr_exp, width = 12, height = 6, units = "in", dpi = 300)

df3 <- df_2yr_exp_clusters %>%
    dplyr::select(kr_cluster, mase_cluster) %>%
    count(kr_cluster, mase_cluster, name = "n")

alluvial_kr_mase_2yr_exp <- ggplot(df3, aes(y = n, axis1 = kr_cluster, axis2 = mase_cluster)) +
    geom_alluvium(aes(fill = kr_cluster), width = 1/12, alpha = 0.8,
                  show.legend = FALSE, discern = TRUE) +
    geom_stratum(width = 1/12, fill = "grey80", color = "grey40",
                 show.legend = FALSE, discern = TRUE) +
    geom_text(stat = "stratum",
              aes(label = after_stat(sub("\\.\\d+$", "", stratum))),
              size = 5, discern = TRUE) +
    scale_fill_manual(values = my_cols) +
    scale_x_discrete(limits = c("KRAFTY", "MASE"), expand = c(.1, .1)) +
    labs(y = "Number of Countries", x = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = "none")
alluvial_kr_mase_2yr_exp
ggsave("alluvial_kr_mase_exp_2yr.png", plot = alluvial_kr_mase_2yr_exp, width = 12, height = 6, units = "in", dpi = 300)

## -------------- kr and mase: importers --------------
set.seed(123)
Zkr_2yr_imp_lists <- Krafty(Z10_imp, Z23_imp)
Zkr_2yr_imp <- Zkr_2yr_imp_lists$Mkr
K_kr_2yr_imp <- Zkr_2yr_imp_lists$K
Zkr_2yr_imp_clusters <- Zkr_2yr_imp_lists$cluster_result
Zkr_2yr_imp_figure <- Zkr_2yr_imp_lists$figure

# save the scree plot
png("screePlot_kr_2yr_imp.png", width = 6, height = 4, units = "in", res = 300)
replayPlot(Zkr_2yr_imp_figure)
dev.off()

Zmase_2yr_imp_lists <- MASE(Z10_imp, Z23_imp)
Zmase_2yr_imp <- Zmase_2yr_imp_lists$Mmase
K_mase_2yr_imp <- Zmase_2yr_imp_lists$K
Zmase_2yr_imp_clusters <- Zmase_2yr_imp_lists$cluster_result
Zmase_2yr_imp_figure <- Zmase_2yr_imp_lists$figure

# save the scree plot
png("screePlot_mase_2yr_imp.png", width = 6, height = 4, units = "in", res = 300)
replayPlot(Zmase_2yr_imp_figure)
dev.off()


ari_imp <- adjustedRandIndex(Zkr_2yr_imp_clusters, Zmase_2yr_imp_clusters)

df_2yr_imp_clusters <- data.frame(name = df_1yr_clusters$name,
                                  kr_cluster = factor(Zkr_2yr_imp_clusters),
                                  mase_cluster = factor(Zmase_2yr_imp_clusters),
                                  importer_10 = factor(Z10_imp_clusters), 
                                  importer_23 = factor(Z23_imp_clusters))

write.xlsx(df_2yr_imp_clusters, file = "clusters_2yr_imp.xlsx")

# Alluvial plot for KR clusters vs Exporter clusters in 2010 & 2023
alluvial_2yr_imp <- alluvial_plot(
    data = df_2yr_imp_clusters,
    kr_col = kr_cluster,
    exporter_col = importer_10,
    importer_col = importer_23,
    joint_label = "KRAFTY",
    pair_label  = "Importer 2010 x Importer 2023",
    type = "imp"
)

ggsave("alluvial_kr_imp_2yr.png", plot = alluvial_2yr_imp, width = 12, height = 6, units = "in", dpi = 300)

df4 <- df_2yr_imp_clusters %>%
    dplyr::select(kr_cluster, mase_cluster) %>%
    count(kr_cluster, mase_cluster, name = "n")

alluvial_kr_mase_2yr_imp <- ggplot(df4, aes(y = n, axis1 = kr_cluster, axis2 = mase_cluster)) +
    geom_alluvium(aes(fill = kr_cluster), width = 1/12, alpha = 0.8,
                  show.legend = FALSE, discern = TRUE) +
    geom_stratum(width = 1/12, fill = "grey80", color = "grey40",
                 show.legend = FALSE, discern = TRUE) +
    geom_text(stat = "stratum",
              aes(label = after_stat(sub("\\.\\d+$", "", stratum))),
              size = 5, discern = TRUE) +
    scale_fill_manual(values = my_cols) +
    scale_x_discrete(limits = c("KRAFTY", "MASE"), expand = c(.1, .1)) +
    labs(y = "Number of Countries", x = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.position = "none")
alluvial_kr_mase_2yr_imp 
ggsave("alluvial_kr_mase_imp_2yr.png", plot = alluvial_kr_mase_2yr_imp, width = 12, height = 6, units = "in", dpi = 300)


