## Plots for varying K1 and K2: Use sigma=0.1 for everything
## Heatmaps of ARI for k1 and k2 with the different True K regimes, and KR-MASE for Z and U fixed [True K unknown]
## Heatmaps of ARI for k1 and k2 with the different True K regimes, and KR, MASE, Z, U all fixed [True K unknown]
## Heatmaps of ARI for k1 and k2 with the different True K regimes, and KR-MASE for Z and U fixed [True K known]
## Heatmaps of ARI for k1 and k2 with the different True K regimes, and KR, MASE, Z, U all fixed [True K known]
## Heatmaps of K estimation for k1 and k2 with the different True K regimes, and KR-MASE for Z and U fixed
##Heatmaps of K estimation for k1 and k2 with the different True K regimes, and KR, MASE, Z, U all fixed

library(tidyverse)
library(gridExtra)
load("ExperimentResults7-22-25_varyingK1K2.Rdata")

k_category <- c("max(K1,K2)", "K1+K2", "floor(K1*K2/2)", "K1*K2")
# k_category <- c("$\\max (K_1,K_2)$", "$\\ K_1+K_2$", "$\\ \text{floor} (K_1*K_2/2)$", "$\\K_1*K_2$")

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


Zplot <- construct_df(Zresults.known, kvec, sigmavec, k_category) # Z results, K known
Zkplot <- construct_df(Zkestresults, kvec, sigmavec, k_category) # Z results for estimating true K
Unknown_Zplot <- construct_df(Zresults.unknown, kvec, sigmavec, k_category) # Z results, K unknown
Uplot <- construct_df(Uresults.known, kvec, sigmavec, k_category) # U results, K known
Unknown_Uplot <- construct_df(Uresults.unknown, kvec, sigmavec, k_category) # U results, K unknown
Ukplot <- construct_df(Ukestresults, kvec, sigmavec, k_category) # # U results for estimating K


Zplot <- Zplot %>% mutate(type = "knownK", matrix = "Z")
Uplot <- Uplot %>% mutate(type = "knownK", matrix = "U")
Unknown_Zplot <- Unknown_Zplot %>% mutate(type = "unknownK", matrix = "Z")
Unknown_Uplot <- Unknown_Uplot %>% mutate(type = "unknownK", matrix = "U")

# include all U and Z results
full_df <- rbind(Zplot, Unknown_Zplot, Uplot, Unknown_Uplot) %>% 
    filter(sigma == 0.1)%>%
    mutate(trueK = factor(trueK, levels = c("max(K1,K2)", "K1+K2", "floor(K1*K2/2)", "K1*K2")),
           K1 = as.factor(K1),
           K2 = as.factor(K2),
           matrix = factor(matrix, levels = c("Z", "U")),
           type = factor(type, levels = c("knownK", "unknownK")))


lab_trueK <- function(x) as.list(parse(text = paste0("K = ", x))) # to control the facet labels
# ARI heatmap
KR_MASE_diff_ARI <- full_df %>%
    select(method, matrix, type, trueK, K1, K2, mean) %>%
    pivot_wider(names_from = method, values_from = mean) %>%
    mutate(diff = KR-MASE, comparison_type = "KR-MASE", group = paste0(matrix, ":", comparison_type)) %>%
    select(group, diff, type, comparison_type, trueK, K1, K2)

max_val <- max(abs(KR_MASE_diff_ARI$diff), na.rm = TRUE) # to control the color scale

p1 <- ggplot(KR_MASE_diff_ARI, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(rows = vars(trueK), cols = vars(type, group), labeller = labeller(trueK=lab_trueK)) +
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
ggsave("New_Figures/varyingK1K2/k1k2_ARI_method.png", plot = p1, width = 12, height = 10, dpi = 300)

KR_MASE_diff_ARI_knownK <- full_df %>%
    filter(type == "knownK")%>%
    select(method, matrix, type, trueK, K1, K2, mean) %>%
    pivot_wider(names_from = method, values_from = mean) %>%
    mutate(diff = KR-MASE, comparison_type = "KR-MASE", group = paste0(matrix, ":", comparison_type)) %>%
    select(group, diff, type, comparison_type, trueK, K1, K2)

max_val <- max(abs(KR_MASE_diff_ARI_knownK$diff), na.rm = TRUE) # to control the color scale
p12 <- ggplot(KR_MASE_diff_ARI_knownK, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(rows = vars(trueK), cols = vars(type, group), labeller = labeller(trueK=lab_trueK)) +
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

ZU_diff_ARI <- full_df %>%
    select(method, matrix, type, trueK, K1, K2, mean) %>%
    pivot_wider(names_from = matrix, values_from = mean) %>%
    mutate(diff = U-Z, comparison_type = "U-Z", group = paste0(method, ":", comparison_type)) %>%
    select(group, diff, type, comparison_type, trueK, K1, K2)

max_val <- max(abs(ZU_diff_ARI$diff), na.rm = TRUE) # to control the color scale

p2 <- ggplot(ZU_diff_ARI, aes(x = K1, y = K2, fill = diff)) +
    geom_tile() +
    facet_grid(rows = vars(trueK), cols = vars(type, group), labeller = labeller(trueK=lab_trueK)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, n.breaks = 10,
                         limits = c(-max_val, max_val))+
    labs(
        # title = "Comparative ARI: Cool Colors\u2192 Z Better, Warm Colors \u2192 U Better\nsigma=0.1, m=101",
        x = "K1",
        y = "K2",
        fill = "ARI (U) - ARI (Z)"
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
ggsave("New_Figures/varyingK1K2/k1k2_ARI_matrix.png", plot = p2, width = 12, height = 10, dpi = 300)
# K estimation absolute error heatmap
Zkplot <- Zkplot %>% mutate(matrix = "Z")
Ukplot <- Ukplot %>% mutate(matrix = "U")
Kresults <- rbind(Zkplot, Ukplot) %>% 
    filter(sigma==0.1)%>%
    mutate(trueK = factor(trueK, levels = c("max(K1,K2)", "K1+K2", "floor(K1*K2/2)", "K1*K2")))

KR_MASE_diff <- Kresults %>%
    select(method, matrix, trueK, K1, K2, mean) %>%
    pivot_wider(names_from = method, values_from = mean) %>%
    mutate(diff = MASE-KR, comparison_type = "MASE-KR", group = paste0(matrix, ":", comparison_type)) %>%
    select(group, diff, comparison_type, trueK, K1, K2)

ZU_diff <- Kresults %>%
    select(method, matrix, trueK, K1, K2, mean) %>%
    pivot_wider(names_from = matrix, values_from = mean) %>%
    mutate(diff = Z - U, comparison_type = "Z-U", group = paste0(method, ":", comparison_type)) %>%
    select(group, diff, comparison_type, trueK, K1, K2)


p3 <- ggplot(KR_MASE_diff, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(group~trueK, labeller = labeller(trueK=lab_trueK, .rows = label_value)) +
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
ggsave("New_Figures/varyingK1K2/k1k2_K_method.png", plot = p3, width = 8, height = 5, dpi = 300)


KR_MASE_diff_smallK <- KR_MASE_diff %>% filter(trueK %in% c("max(K1,K2)", "K1+K2"))

max_val <- max(abs(KR_MASE_diff_smallK$diff), na.rm = TRUE) # to control the color scale
p4 <- ggplot(KR_MASE_diff_smallK, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(trueK~group, labeller = labeller(trueK=lab_trueK, .cols = label_value)) +
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

ggsave("New_Figures/varyingK1K2/k1k2_smallK_method.png", plot = p4, width = 8, height = 5, dpi = 300)

KR_MASE_diff_largeK <- KR_MASE_diff %>% filter(trueK %in% c("floor(K1*K2/2)", "K1*K2"))
max_val <- max(abs(KR_MASE_diff_largeK$diff), na.rm = TRUE) # to control the color scale
p5 <- ggplot(KR_MASE_diff_largeK, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(trueK~group, labeller = labeller(trueK=lab_trueK, .cols = label_value)) +
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
ggsave("New_Figures/varyingK1K2/k1k2_largeK_method.png", plot = p5, width = 8, height = 5, dpi = 300)


max_val <- max(abs(ZU_diff$diff), na.rm = TRUE) # to control the color scale
p6 <- ggplot(ZU_diff, aes(x = factor(K1), y = factor(K2), fill = diff)) +
    geom_tile() +
    facet_grid(group~trueK, labeller = labeller(trueK=lab_trueK, .rows = label_value)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, n.breaks = 10,
                         limits = c(-max_val, max_val))+
    labs(
        # title = "Comparative K Estimation Error: Cool Colors\u2192 Z Better, Warm Colors \u2192 U Better\nsigma=0.1, m=101",
        x = "K1",
        y = "K2",
        fill = "MAE (Z) - MAE (U)"
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
ggsave("New_Figures/varyingK1K2/k1k2_matrix.png", plot = p6, width = 12, height = 5, dpi = 300)
