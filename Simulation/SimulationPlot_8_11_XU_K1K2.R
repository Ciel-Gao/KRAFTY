## Plots for varying K1 and K2: Use sigma=0.1 for everything
## Heatmaps of ARI for k1 and k2 with the different True K regimes, and KR-MASE for X and U fixed [True K unknown]
## Heatmaps of ARI for k1 and k2 with the different True K regimes, and KR, MASE, X, U all fixed [True K unknown]
## Heatmaps of ARI for k1 and k2 with the different True K regimes, and KR-MASE for X and U fixed [True K known]
## Heatmaps of ARI for k1 and k2 with the different True K regimes, and KR, MASE, X, U all fixed [True K known]
## Heatmaps of K estimation for k1 and k2 with the different True K regimes, and KR-MASE for X and U fixed
##Heatmaps of K estimation for k1 and k2 with the different True K regimes, and KR, MASE, X, U all fixed

library(tidyverse)
library(gridExtra)

k_category <- c("max(K1,K2)", "K1+K2", "floor(K1*K2/2)", "K1*K2")

construct_df <- function(results, kvec, sigmavec, k_category){
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
    
    m <- 101
    
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
        for (k1_idx in 1:length(kvec)) {
            k1 <- kvec[k1_idx]
            for (k2_idx in 1:length(kvec)) {
                k2 <- kvec[k2_idx]
                for (k_idx in 1: 4) {
                    k_cat <- k_category[k_idx]
                    for (sigma_idx in 1:length(sigmavec)){
                        sigma <- sigmavec[sigma_idx]
                        df <- rbind(df, data.frame(
                            method = method,
                            K1 = k1,
                            K2 = k2,
                            m = m,
                            trueK = k_cat,
                            sigma = sigma,
                            mean = results[k1_idx, k2_idx, k_idx, sigma_idx, meanidx],
                            lower = results[k1_idx, k2_idx, k_idx, sigma_idx, loweridx],
                            upper = results[k1_idx, k2_idx, k_idx, sigma_idx, upperidx]
                        ))
                    }
                }
            }
        }
    }
    
    
    
    return(df)
}


Xplot <- construct_df(Xresults.known, kvec, sigmavec, k_category) # X results, K known
Xkplot <- construct_df(Xkestresults, kvec, sigmavec, k_category) # X results for estimating true K
Unknown_Xplot <- construct_df(Xresults.unknown, kvec, sigmavec, k_category) # X results, K unknown
Uplot <- construct_df(Uresults.known, kvec, sigmavec, k_category) # U results, K known
Unknown_Uplot <- construct_df(Uresults.unknown, kvec, sigmavec, k_category) # U results, K unknown
Ukplot <- construct_df(Ukestresults, kvec, sigmavec, k_category) # # U results for estimating K


Xplot <- Xplot %>% mutate(type = "knownK", matrix = "X")
Uplot <- Uplot %>% mutate(type = "knownK", matrix = "U")
Unknown_Xplot <- Unknown_Xplot %>% mutate(type = "unknownK", matrix = "X")
Unknown_Uplot <- Unknown_Uplot %>% mutate(type = "unknownK", matrix = "U")

# include all U and X results
full_df <- rbind(Xplot, Unknown_Xplot, Uplot, Unknown_Uplot) %>% 
    filter(sigma == 0.1)%>%
    mutate(trueK = factor(trueK, levels = c("max(K1,K2)", "K1+K2", "floor(K1*K2/2)", "K1*K2")),
           K1 = as.factor(K1),
           K2 = as.factor(K2),
           matrix = factor(matrix, levels = c("X", "U")),
           type = factor(type, levels = c("knownK", "unknownK")))

# facet labels
lab_sigma <- function(x) as.list(parse(text = paste0("sigma^2 = ", x)))
lab_trueK <- function(x) as.list(parse(text = paste0("K = ", x)))

# ARI heatmap
KR_MASE_diff_ARI <- full_df %>%
    select(method, matrix, type, trueK, K1, K2, mean) %>%
    pivot_wider(names_from = method, values_from = mean) %>%
    mutate(diff = KR-MASE, comparison_type = "KR-MASE", group = paste0(matrix, ":", comparison_type)) %>%
    select(group, diff, type, comparison_type, trueK, K1, K2)

max_val <- max(abs(KR_MASE_diff_ARI$diff), na.rm = TRUE) # to control the color scale

p1 <- ggplot(KR_MASE_diff_ARI, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(rows = vars(trueK), cols = vars(type, group), labeller = labeller(.rows = lab_trueK)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, n.breaks = 10,
                         limits = c(-max_val, max_val)) +
    labs(
        x = "K1",
        y = "K2",
        fill = "ARI (KR) - ARI (MASE)"
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

p1
# ggsave("New_Figures/k1k2_XU_ARI_method.png", plot = p1, width = 12, height = 10, dpi = 300)

KR_MASE_diff_ARI_knownK <- full_df %>%
    filter(type == "knownK")%>%
    select(method, matrix, type, trueK, K1, K2, mean) %>%
    pivot_wider(names_from = method, values_from = mean) %>%
    mutate(diff = KR-MASE, comparison_type = "KR-MASE", group = paste0(matrix, ":", comparison_type)) %>%
    select(group, diff, type, comparison_type, trueK, K1, K2)

max_val <- max(abs(KR_MASE_diff_ARI_knownK$diff), na.rm = TRUE) # to control the color scale
p12 <- ggplot(KR_MASE_diff_ARI_knownK, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(rows = vars(trueK), cols = vars(type, group), labeller = labeller(.rows = label_both)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, n.breaks = 10,
                         limits = c(-max_val, max_val))+
    labs(
        title = "Comparative ARI: Cool Colors\u2192 MASE Better, Warm Colors \u2192 KR Better\n sigma=0.1, m=101",
        x = "K1",
        y = "K2",
        fill = "ARI (KR) - ARI (MASE)"
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

p12

# ggsave("k1k2_ARI_method.png", plot = p1, width = 12, height = 10, dpi = 300)

XU_diff_ARI <- full_df %>%
    select(method, matrix, type, trueK, K1, K2, mean) %>%
    pivot_wider(names_from = matrix, values_from = mean) %>%
    mutate(diff = U-X, comparison_type = "U-X", group = paste0(method, ":", comparison_type)) %>%
    select(group, diff, type, comparison_type, trueK, K1, K2)

max_val <- max(abs(XU_diff_ARI$diff), na.rm = TRUE) # to control the color scale

p2 <- ggplot(XU_diff_ARI, aes(x = K1, y = K2, fill = diff)) +
    geom_tile() +
    facet_grid(rows = vars(trueK), cols = vars(type, group), labeller = labeller(.rows = lab_trueK)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, n.breaks = 10,
                         limits = c(-max_val, max_val))+
    labs(
        # title = "Comparative ARI: Cool Colors\u2192 X Better, Warm Colors \u2192 U Better\nsigma=0.1, m=101",
        x = "K1",
        y = "K2",
        fill = "ARI (U) - ARI (X)"
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

p2
ggsave("New_Figures/XU_K1K2/k1k2_XU_ARI_matrix.png", plot = p2, width = 12, height = 10, dpi = 300)
# K estimation absolute error heatmap
Xkplot <- Xkplot %>% mutate(matrix = "X")
Ukplot <- Ukplot %>% mutate(matrix = "U")
Kresults <- rbind(Xkplot, Ukplot) %>% 
    filter(sigma==0.1)%>%
    mutate(trueK = factor(trueK, levels = c("max(K1,K2)", "K1+K2", "floor(K1*K2/2)", "K1*K2")))

KR_MASE_diff <- Kresults %>%
    select(method, matrix, trueK, K1, K2, mean) %>%
    pivot_wider(names_from = method, values_from = mean) %>%
    mutate(diff = MASE-KR, comparison_type = "MASE-KR", group = paste0(matrix, ":", comparison_type)) %>%
    select(group, diff, comparison_type, trueK, K1, K2)

XU_diff <- Kresults %>%
    select(method, matrix, trueK, K1, K2, mean) %>%
    pivot_wider(names_from = matrix, values_from = mean) %>%
    mutate(diff = X - U, comparison_type = "X-U", group = paste0(method, ":", comparison_type)) %>%
    select(group, diff, comparison_type, trueK, K1, K2)


p3 <- ggplot(KR_MASE_diff, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(group~trueK, labeller = labeller(.cols = label_both, .rows = label_value)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, n.breaks = 10)+
    labs(
        # title = "Comparative K Estimation Error: Cool Colors\u2192 MASE Better, Warm Colors \u2192 KR Better\nsigma=0.1, m=101",
        x = "K1",
        y = "K2",
        fill = "MAE (MASE) - MAE (KR)"
    ) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.key.height = unit(1.5, "cm"),       # Increase height of each legend key
        legend.key.width = unit(1, "cm"),          # Increase width
        legend.box.margin = margin(10, 10, 10, 10) # Add margin around the legend
    )
p3

KR_MASE_diff_smallK <- KR_MASE_diff %>% filter(trueK %in% c("max(k1,k2)", "k1+k2"))

max_val <- max(abs(KR_MASE_diff_smallK$diff), na.rm = TRUE) # to control the color scale
p4 <- ggplot(KR_MASE_diff_smallK, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(trueK~group, labeller = labeller(.rows = label_both, .cols = label_value)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, n.breaks = 10,
                         limits = c(-max_val, max_val))+
    labs(
        # title = "Comparative K Estimation Error: Cool Colors\u2192 MASE Better, Warm Colors \u2192 KR Better\nsmall K cases, sigma=0.1, m=101",
        x = "K1",
        y = "K2",
        fill = "MAE (MASE) - MAE (KR)"
    ) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.key.height = unit(1, "cm"),       # Increase height of each legend key
        legend.key.width = unit(1, "cm"),          # Increase width
        legend.box.margin = margin(10, 10, 10, 10) # Add margin around the legend
    )
p4

# ggsave("New_Figures/k1k2_smallK_method.png", plot = p4, width = 8, height = 5, dpi = 300)

KR_MASE_diff_largeK <- KR_MASE_diff %>% filter(trueK %in% c("floor(k1*k2/2)", "k1*k2"))
max_val <- max(abs(KR_MASE_diff_largeK$diff), na.rm = TRUE) # to control the color scale
p5 <- ggplot(KR_MASE_diff_largeK, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(trueK~group, labeller = labeller(.rows = label_both, .cols = label_value)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, n.breaks = 10, 
                         limits = c(-max_val, max_val))+
    labs(
        # title = "Comparative K Estimation Error: Cool Colors\u2192 MASE Better, Warm Colors \u2192 KR Better\nlarge K cases, sigma=0.1, m=101",
        x = "K1",
        y = "K2",
        fill = "MAE (MASE) - MAE (KR)"
    ) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.key.height = unit(1, "cm"),       # Increase height of each legend key
        legend.key.width = unit(1, "cm"),          # Increase width
        legend.box.margin = margin(10, 10, 10, 10) # Add margin around the legend
    )
p5
# ggsave("New_Figures/k1k2_largeK_method.png", plot = p5, width = 8, height = 5, dpi = 300)


max_val <- max(abs(XU_diff$diff), na.rm = TRUE) # to control the color scale
p6 <- ggplot(XU_diff, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(group~trueK, labeller = labeller(.cols = lab_trueK, .rows = label_value)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, n.breaks = 10,
                         limits = c(-max_val, max_val))+
    labs(
        # title = "Comparative K Estimation Error: Cool Colors\u2192 X Better, Warm Colors \u2192 U Better\nsigma=0.1, m=101",
        x = "K1",
        y = "K2",
        fill = "MAE (X) - MAE (U)"
    ) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        legend.key.height = unit(1, "cm"),       # Increase height of each legend key
        legend.key.width = unit(1, "cm"),          # Increase width
        legend.box.margin = margin(10, 10, 10, 10) # Add margin around the legend
    )
p6
ggsave("New_Figures/XU_K1K2/k1k2_XU_matrix.png", plot = p6, width = 12, height = 5, dpi = 300)
