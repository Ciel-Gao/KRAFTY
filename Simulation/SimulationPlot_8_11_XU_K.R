## Simulation results plot 6.24
## Fixed sigma=0.1, true K varying. Panel (a) with X matrix, panel (b) with U matrix [True K known and unknown]
## True K = 4, 9, 15, sigma varying. Top 3 panels with X matrix, bottom 3 panels with U matrix. [True K known and unknown]
## Combined U/X sigma =0.01 and 1, true K varying. [True K unknown]
## Combined U/X true K = 4, 9, 15, sigma varying. Definitely make a comment about ARI increasing with sigma! [True K unknown]
## Combined U/X, sigma=0.1, true K varying. [True K known]
## 4-heatmap of sigma vs. TrueK with U vs. X and KR vs. MASE fixed
## try looking at the heatmap of sigma vs. TrueK where we have error for X - error for U, with KR vs. MASE fixed; and the heatmaps of sigma vs. TrueK where we have error for MASE - error for KR, with U vs. X fixed.

library(tidyverse)
library(gridExtra)



num_K <- length(trueKvec)
num_sigma <- length(sigmavec)

# function used to organize results
construct_df <- function(results, kvec, sigmavec){
    methods <- c("KR", "MASE")
    
    df <- data.frame(
        method = character(),
        K1 = numeric(),
        K2 = numeric(),
        m = numeric(),
        trueK = numeric(),
        sigma = numeric(),
        mean = numeric(),
        lower = numeric(),
        upper = numeric(),
        stringsAsFactors = TRUE
    )
    
    k1 <- 4
    k2 <- 4
    m <- 20
    
    for (method in methods) {
        if (method == "KR") {
            meanidx <- 2
            loweridx <- 1
            upperidx <- 3
        }else{
            meanidx <- 5
            loweridx <- 4
            upperidx <- 6
        }
        for (k_idx in 1:length(kvec)) {
            k <- kvec[k_idx]
            for (sigma_idx in 1:length(sigmavec)){
                sigma <- sigmavec[sigma_idx]
                df <- rbind(df, data.frame(
                    method = method,
                    K1 = k1,
                    K2 = k2,
                    m = m,
                    trueK = k,
                    sigma = sigma,
                    mean = results[k_idx, sigma_idx, meanidx],
                    lower = results[k_idx, sigma_idx, loweridx],
                    upper = results[k_idx, sigma_idx, upperidx]
                ))
            }
        }
    }
    
    return(df)
}


Xplot <- construct_df(Xresults.known, trueKvec, sigmavec) # X results, K known
Xkplot <- construct_df(Xkestresults, trueKvec, sigmavec) # X results for estimating true K
Unknown_Xplot <- construct_df(Xresults.unknown, trueKvec, sigmavec) # X results, K unknown
Uplot <- construct_df(Uresults.known, trueKvec, sigmavec) # U results, K known
Unknown_Uplot <- construct_df(Uresults.unknown, trueKvec, sigmavec) # U results, K unknown
Ukplot <- construct_df(Ukestresults, trueKvec, sigmavec) # U results for estimating K

# X with K un/known, U with K un/known, both X and U with K unknown
Xplot <- Xplot %>% mutate(type = "knownK", matrix = "X")
Uplot <- Uplot %>% mutate(type = "knownK", matrix = "U")
Unknown_Xplot <- Unknown_Xplot %>% mutate(type = "unknownK", matrix = "X")
Unknown_Uplot <- Unknown_Uplot %>% mutate(type = "unknownK", matrix = "U")

# include all U and X results
full_df <- rbind(Xplot, Unknown_Xplot, Uplot, Unknown_Uplot) %>% 
    mutate(sigma = as.factor(sigma), trueK = as.factor(trueK),
           matrix = factor(matrix, levels = c("X", "U")),
           type = factor(type, levels = c("knownK", "unknownK")))



# function to create plots
draw_plots <- function(df, changing_variable, name, metric, combine = FALSE){
    lab_trueK <- function(x) as.list(parse(text = paste0("K = ", x))) # to control the facet labels
    lab_sigma <- function(x) as.list(parse(text = paste0("sigma^2 = ", x)))
    
    color_vals <- c("#c31e23", "#ff5a5e", "#0d7d87", "#99c6cc")
    shape_vals <- c(16, 17, 15, 18) 
    color_mapping <- setNames(color_vals, levels(df$group_id))
    shape_mapping <- setNames(shape_vals, levels(df$group_id))
    
    # metric is ARI or estimated value of true K
    p <- ggplot(df, aes(x = .data[[changing_variable]], group = group_id)) +
        geom_point(aes(y = mean, color = group_id, shape = group_id), 
                   position = position_dodge(width = 0.1), size = 1.5) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = group_id), alpha = 0.6, color = NA)+
        theme(
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            plot.title = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 13),
            legend.text = element_text(size = 12)
        ) +
        scale_color_manual(values = color_mapping) +
        scale_fill_manual(values = color_mapping)+
        scale_shape_manual(values = shape_mapping)
    
    
    if (changing_variable == "trueK" & combine == FALSE) {
        p <- p + labs(
            # title = paste0(name, ", \n", "K1=4, K2=4, m=20, sigma=0.1"),
            x = "K",
            y = metric) +
            geom_vline(xintercept = 5, linetype = "dashed", color = "red")+
            facet_wrap(~matrix, ncol = 3, labeller = label_both)+
            theme(
                axis.text.x = element_text(size = 12, face = "bold"),
                axis.text.y = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 14, face = "bold"),
                axis.title.y = element_text(size = 14, face = "bold"),
                plot.title = element_text(size = 16, face = "bold"),
                legend.title = element_text(size = 13, face = "bold"),
                legend.text = element_text(size = 12, face = "bold"),
                strip.text = element_text(size = 14, face = "bold")
            )+ 
            labs(color = "Method, Type", fill = "Method, Type", shape = "Method, Type")
        
    }else if(changing_variable == "trueK" & combine == TRUE){
        p <- p + labs(
            # title = paste0(name, ", K1=4, K2=4, m=20"),
            x = "K",
            y = metric) +
            geom_vline(xintercept = 5, linetype = "dashed", color = "red")+
            facet_wrap(~sigma, labeller = labeller(sigma=lab_sigma))+
            theme(
                axis.text.x = element_text(size = 12, face = "bold"),
                axis.text.y = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 14, face = "bold"),
                axis.title.y = element_text(size = 14, face = "bold"),
                plot.title = element_text(size = 16, face = "bold"),
                legend.title = element_text(size = 13, face = "bold"),
                legend.text = element_text(size = 12, face = "bold"),
                strip.text = element_text(size = 14, face = "bold")
            )+
            labs(color = "Method, Matrix", fill = "Method, Matrix", shape = "Method, Matrix")
        
    }else if(changing_variable == "sigma" & combine == FALSE){
        p <- p + labs(
            # title = paste0(name, ", \n", "K1=4, K2=4, m=20"),
            x = expression(sigma^2),
            y = metric) +
            facet_grid(matrix ~ trueK, labeller = labeller(matrix=label_both,
                                                           trueK=lab_trueK))+
            theme(
                axis.text.x = element_text(size = 12, face = "bold"),
                axis.text.y = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 14, face = "bold"),
                axis.title.y = element_text(size = 14, face = "bold"),
                plot.title = element_text(size = 16, face = "bold"),
                legend.title = element_text(size = 13, face = "bold"),
                legend.text = element_text(size = 12, face = "bold"),
                strip.text = element_text(size = 14, face = "bold")
            )+
            labs(color = "Method, Type", fill = "Method, Type", shape = "Method, Type")
    }else if(changing_variable == "sigma" & combine == TRUE){
        p <- p + labs(
            # title = paste0(name, ", K1=4, K2=4, m=20"),
            x = expression(sigma^2),
            y = metric) +
            facet_wrap(~trueK, ncol=3, labeller = labeller(trueK=lab_trueK))+
            theme(
                axis.text.x = element_text(size = 12, face = "bold"),
                axis.text.y = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 14, face = "bold"),
                axis.title.y = element_text(size = 14, face = "bold"),
                plot.title = element_text(size = 16, face = "bold"),
                legend.title = element_text(size = 13, face = "bold"),
                legend.text = element_text(size = 12, face = "bold"),
                strip.text = element_text(size = 14, face = "bold")
            )+
            labs(color = "Method, Matrix", fill = "Method, Matrix", shape = "Method, Matrix")
    }
    
    return(p)
}

# plot 1: Fixed sigma=0.1, true K varying. 
# Panel (a) with X matrix, panel (b) with U matrix [True K known and unknown] 

# keep records only when sigma == 0.1
df_sigma01 <- full_df %>% filter(sigma == 0.1) %>% 
    mutate(group_id = factor(paste0(method, ".", type),
                             levels = c("KR.knownK", "KR.unknownK", "MASE.knownK", "MASE.unknownK"),
                             labels = c("KR, K known", "KR, K unknown", 
                                        "MASE, K known", "MASE, K unknown")))

p1 <- draw_plots(
    df = df_sigma01,
    changing_variable = "trueK",
    name = "Joint clusters estimation: varying K with sigma fixed, X and U results",
    metric = "ARI") + ylim(0.2, 1)
plot(p1)
ggsave("New_Figures/XU_K/K_XU.png", plot = p1, width = 10, height = 5, dpi = 300)

# plot 2: True K = 4, 9, 15, sigma varying. 
# Top 3 panels with X matrix, bottom 3 panels with U matrix. [True K known and unknown]

# keep records only when K = 4, 9, 15
df_K4915 <- full_df %>% filter(trueK %in% c(4, 9, 15)) %>%
    mutate(group_id = factor(paste0(method, ".", type),
                             levels = c("KR.knownK", "KR.unknownK", "MASE.knownK", "MASE.unknownK"),
                             labels = c("KR, K known", "KR, K unknown", 
                                        "MASE, K known", "MASE, K unknown")))


p2 <- draw_plots(
    df = df_K4915,
    changing_variable = "sigma",
    name = "Joint clusters estimation: varying sigma with K fixed, X and U results",
    metric = "ARI") + ylim(0.2, 1)
plot(p2)
ggsave("New_Figures/XU_K/K_XU_sigma.png", plot = p2, width = 12, height = 5, dpi = 300)
# plot 3: Combined U/X sigma =0.01 and 1, true K varying. [True K unknown]
# keep records when k is unknown, sigma == 0.01 or 1
df_compareXU_K <- full_df %>% filter(sigma %in% c(0.01, 1), type == "unknownK") %>%
    mutate(group_id = paste0(method, ", ", matrix)) %>%
    mutate(group_id = factor(group_id, levels = c("KR, X", "KR, U", "MASE, X", "MASE, U")))


p3 <- draw_plots(
    df = df_compareXU_K,
    changing_variable = "trueK",
    name = "Joint clusters estimation: varying K, comparing X and U results",
    metric = "ARI",
    combine = TRUE) + ylim(0.2, 1)
plot(p3)
ggsave("New_Figures/XU_K/K_compareXU_K.png", plot = p3, width = 10, height = 5, dpi = 300)
# plot 4: Combined U/X true K = 4, 9, 15, sigma varying. [True K unknown]
# keep records when k is unknown
df_compareXU_sigma <- full_df %>% filter(trueK %in% c(4, 9, 15), type == "unknownK") %>%
    mutate(group_id = paste0(method, ", ", matrix)) %>%
    mutate(group_id = factor(group_id, levels = c("KR, X", "KR, U", "MASE, X", "MASE, U")))

p4 <- draw_plots(
    df = df_compareXU_sigma,
    changing_variable = "sigma",
    name = "Joint clusters estimation: varying sigma, comparing X and U results",
    metric = "ARI",
    combine = TRUE) + ylim(0.2, 1)
plot(p4)
ggsave("New_Figures/XU_K/K_compareXU_sigma.png", plot = p4, width = 10, height = 5, dpi = 300)

# plot 5: Combined U/X, sigma=0.1, true K varying. [True K known]
# keep records when k is known and sigma = 0.1
df_compareXU_knownK <- full_df %>% filter(sigma == 0.1, type == "knownK") %>%
    mutate(group_id = paste0(method, ", ", matrix)) %>%
    mutate(group_id = factor(group_id, levels = c("KR, X", "KR, U", "MASE, X", "MASE, U")))

p5 <- draw_plots(
    df = df_compareXU_knownK,
    changing_variable = "trueK",
    name = "Joint clusters estimation: varying true K, comparing X and U results, K known, \nsigma=0.1",
    metric = "ARI",
    combine = TRUE) + ylim(0.7, 1)
plot(p5)

ggsave("New_Figures/XU_K/K_compareXU_knownK.png", plot = p5, width = 6, height = 5, dpi = 300)

# plot 6: heatmap
lab_trueK <- function(x) as.list(parse(text = paste0("K = ", x))) # to control the facet labels
lab_sigma <- function(x) as.list(parse(text = paste0("sigma^2 = ", x)))

# 4-heatmap of sigma vs. TrueK with U vs. X and KR vs. MASE fixed
Xkplot <- Xkplot %>% mutate(matrix = "X")
Ukplot <- Ukplot %>% mutate(matrix = "U")
Kresults <- rbind(Xkplot, Ukplot)

KR_MASE_diff <- Kresults %>%
    select(method, matrix, trueK, sigma, mean) %>%
    pivot_wider(names_from = method, values_from = mean) %>%
    mutate(diff = MASE-KR, comparison_type = "MASE-KR", group = paste0(matrix, ":", comparison_type)) %>%
    select(group, diff, comparison_type, trueK, sigma)

XU_diff <- Kresults %>%
    select(method, matrix, trueK, sigma, mean) %>%
    pivot_wider(names_from = matrix, values_from = mean) %>%
    mutate(diff = X - U, comparison_type = "X-U", group = paste0(method, ":", comparison_type)) %>%
    select(group, diff, comparison_type, trueK, sigma)

# K_comparison <- rbind(KR_MASE_diff, XU_diff)

# my_palette <- c("#008080","#70a494","#b4c8a8","#f6edbd","#edbb8a","#de8a5a","#ca562c")
# my_palette <- c("#70284a", "#9c3f5d", "#c8586c", "#ee8a82", "#fbe6c5")

max_val <- max(abs(XU_diff$diff), na.rm = TRUE) # to control the color scale

p6 <- ggplot(XU_diff, aes(x = factor(sigma), y = factor(trueK), fill = diff)) +
    geom_tile() +
    facet_wrap(~group, ncol = 2) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, 
                         n.breaks = 15,
                         limits = c(-max_val, max_val))+
    labs(
        # title = "Comparative K Estimation Error: Cool Colors\u2192 X Better, Warm Colors \u2192 U Better",
        x = expression(sigma^2),
        y = "K",
        fill = "MAE (X) - MAE (U)"
    ) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.key.height = unit(1.2, "cm")
    )
p6
ggsave("New_Figures/XU_K/K_compareXU_estK.png", plot = p6, width = 9, height = 5, dpi = 300)

max_val <- max(abs(KR_MASE_diff$diff), na.rm = TRUE) # to control the color scale

p7 <- ggplot(KR_MASE_diff, aes(x = factor(sigma), y = factor(trueK), fill = diff)) +
    geom_tile() +
    facet_wrap(~group, ncol=2) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, 
                         limits = c(-max_val, max_val),
                         n.breaks = 15)+  
    labs(
        # title = "Comparative K Estimation Error: Warm Colors \u2192 KR Better, Cool Colors\u2192 MASE Better",
        x = expression(sigma^2),
        y = "K",
        fill = "MAE (MASE) - MAE (KR)"
    ) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.key.height = unit(1.2, "cm")
    )
p7
ggsave("New_Figures/XU_K/K_compareMethod_estK.png", plot = p7, width = 9, height = 5, dpi = 300)



