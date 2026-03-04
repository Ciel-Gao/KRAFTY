## Simulation results plot 8.11
## Varying m and sigma with K fixed, X results [True K known and unknown], keep the whole thing.
## Varying m and sigma with K fixed, U results [True K known and unknown], keep the whole thing.
## Varying m and sigma, X and U [True K unknown], keep the whole thing.
## Varying m and sigma, X and U [True K known], keep the whole thing.
## Heatmaps for K estimation with trueK =6,12 and KR vs. MASE fixed, 2x2 with X
## Heatmaps for K estimation with trueK =6,12 and KR vs. MASE fixed, 2x2 with U
## Heatmaps for K estimation with trueK =6,12 and KR vs. MASE fixed, 2x2 with X-U
## Heatmaps for ARI with trueK =6,12 and KR vs. MASE fixed, 2x2 with X
## Heatmaps for ARI with trueK =6,12 and KR vs. MASE fixed, 2x2 with U
## Heatmaps for ARI with trueK =6,12 and KR vs. MASE fixed, 2x2 with X-U


library(tidyverse)
library(gridExtra)

construct_df <- function(results, kvec, sigmavec, mvec){
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
        for (m_idx in 1: length(mvec)) {
            m <- mvec[m_idx]
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
                        mean = results[m_idx, k_idx, sigma_idx, meanidx],
                        lower = results[m_idx, k_idx, sigma_idx, loweridx],
                        upper = results[m_idx, k_idx, sigma_idx, upperidx]
                    ))
                }
            }
        }
        
    }
    
    return(df)
}


Xplot <- construct_df(Xresults.known, trueKvec, sigmavec, mvec) # X results, K known
Xkplot <- construct_df(Xkestresults, trueKvec, sigmavec, mvec) # X results for estimating true K
Unknown_Xplot <- construct_df(Xresults.unknown, trueKvec, sigmavec, mvec) # X results, K unknown
Uplot <- construct_df(Uresults.known, trueKvec, sigmavec, mvec) # U results, K known
Unknown_Uplot <- construct_df(Uresults.unknown, trueKvec, sigmavec, mvec) # U results, K unknown
Ukplot <- construct_df(Ukestresults, trueKvec, sigmavec, mvec) # # U results for estimating K


# X with K un/known, U with K un/known, both X and U with K unknown
Xplot <- Xplot %>% mutate(type = "knownK", matrix = "X")
Uplot <- Uplot %>% mutate(type = "knownK", matrix = "U")
Unknown_Xplot <- Unknown_Xplot %>% mutate(type = "unknownK", matrix = "X")
Unknown_Uplot <- Unknown_Uplot %>% mutate(type = "unknownK", matrix = "U")

# include all U and X results
full_df <- rbind(Xplot, Unknown_Xplot, Uplot, Unknown_Uplot) %>% 
    mutate(sigma = as.factor(sigma), 
           trueK = as.factor(trueK),
           m = as.factor(m),
           matrix = factor(matrix, levels = c("X", "U")),
           type = factor(type, levels = c("knownK", "unknownK")))


# function to create plots
draw_plots <- function(df, changing_variable, name, metric){
    
    color_vals <- c("#c31e23", "#ff5a5e", "#0d7d87", "#99c6cc")
    shape_vals <- c(16, 17, 15, 18) 
    color_mapping <- setNames(color_vals, levels(df$group_id))
    shape_mapping <- setNames(shape_vals, levels(df$group_id))
    
    # metric is ARI or estimated value of true K
    p <- ggplot(df, aes(x = .data[[changing_variable]], group = group_id)) +
        geom_point(aes(y = mean, color = group_id, shape = group_id), position = position_dodge(width = 0.1), size = 1.5) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = group_id), alpha = 0.6, color = NA)+
        theme(
            axis.text.x = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 12, face = "bold"),
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            plot.title = element_text(size = 16, face = "bold"),
            legend.title = element_text(size = 13, face = "bold"),
            legend.text = element_text(size = 12, face = "bold")
        )+
        labs(color = "Method, Matrix", fill = "Method, Matrix", shape = "Method, Matrix") +
        scale_color_manual(values = color_mapping) +
        scale_fill_manual(values = color_mapping)+
        scale_shape_manual(values = shape_mapping)
    
    
    if (changing_variable == "m") {
        p <- p + labs(x = "p",
                      y = metric) +
            facet_grid(trueK~sigma, labeller = label_both)+
            theme(
                panel.background = element_rect(fill = "white", color = NA),
                plot.background = element_rect(fill = "white", color = NA),
                strip.background = element_rect(fill = "gray80", color = NA),
                strip.text = element_text(size = 12, face = "bold"),
                panel.grid.major = element_line(color = "gray90"),  
                panel.grid.minor = element_blank(),
                panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)
            )
    }else if(changing_variable == "sigma"){
        p <- p + labs(x = "sigma",
                      y = metric) +
            facet_wrap(~trueK, ncol = 3, labeller = label_both)+
            theme(
                panel.background = element_rect(fill = "white", color = NA),
                plot.background = element_rect(fill = "white", color = NA),
                strip.background = element_rect(fill = "gray80", color = NA),
                strip.text = element_text(size = 12, face = "bold"),
                panel.grid.major = element_line(color = "gray90"),  
                panel.grid.minor = element_blank(),
                panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)
            )
        
    }
    
    return(p)
}

# plot 1: Varying m and sigma with K fixed, X results [True K known and unknown]
X_combine <- full_df %>% filter(matrix == "X") %>% 
    mutate(group_id = factor(paste0(method, ".", type),
                             levels = c("KR.knownK", "KR.unknownK", "MASE.knownK", "MASE.unknownK"),
                             labels = c("KR, X with K known", "KR, X with K unknown", "MASE, X with K known", "MASE, X with K unknown")))
p1 <- draw_plots(
    df = X_combine,
    changing_variable = "m",
    name = "Joint clusters estimation: varying m with sigma and K fixed, X results",
    metric = "ARI") + ylim(0, 1)
plot(p1)
ggsave("New_Figures/m_X.png", plot = p1, width = 12, height = 5, dpi = 300)



# plot 2: Varying m and sigma with K fixed, U results [True K known and unknown].
U_combine <- full_df %>% filter(matrix == "U") %>% 
    mutate(group_id = factor(paste0(method, ".", type),
                             levels = c("KR.knownK", "KR.unknownK", "MASE.knownK", "MASE.unknownK"),
                             labels = c("KR, U with K known", "KR, U with K unknown", "MASE, U with K known", "MASE, U with K unknown")))


p2 <- draw_plots(
    df = U_combine,
    changing_variable = "m",
    name = "Joint clusters estimation: varying m with sigma and K fixed, U results",
    metric = "ARI") + ylim(0, 1)
plot(p2)
ggsave("New_Figures/m_U.png", plot = p2, width = 12, height = 5, dpi = 300)

# plot 3: Varying m and sigma, X and U [True K unknown].
XU_unknownK <- full_df %>% filter(type == "unknownK") %>% 
    mutate(group_id = factor(paste0(method, ".", matrix),
                             levels = c("KR.X", "KR.U", "MASE.X", "MASE.U"),
                             labels = c("KR, X with K unknown", "KR, U with K unknown", "MASE, X with K unknown", "MASE, U with K unknown")))

p3 <- draw_plots(
    df = XU_unknownK ,
    changing_variable = "m",
    name = "Joint clusters estimation: varying m with K unknown, X and U results",
    metric = "ARI") + ylim(0, 1)
plot(p3)
ggsave("New_Figures/m_XU_unknown.png", plot = p3, width = 12, height = 5, dpi = 300)

# plot 4: Varying m and sigma, X and U [True K known].
XU_knownK <- full_df %>% filter(type == "knownK") %>% 
    mutate(group_id = factor(paste0(method, ".", matrix),
                             levels = c("KR.X", "KR.U", "MASE.X", "MASE.U"),
                             labels = c("KR, X with K known", "KR, U with K known", "MASE, X with K known", "MASE, U with K known")))

p4 <- draw_plots(
    df = XU_knownK ,
    changing_variable = "m",
    name = "Joint clusters estimation: varying m with K known, X and U results",
    metric = "ARI") + ylim(0, 1)
plot(p4)
ggsave("New_Figures/m_XU_known.png", plot = p4, width = 12, height = 5, dpi = 300)

# plot 5: Heatmaps for K estimation with trueK =6,12 and KR vs. MASE fixed, 2x2 with X
Xkplot <- Xkplot %>% mutate(matrix = "X")
Ukplot <- Ukplot %>% mutate(matrix = "U")
Kresults <- rbind(Xkplot, Ukplot)

KR_MASE_diff <- Kresults %>%
    dplyr::select(method, matrix, trueK, sigma, m, mean) %>%
    pivot_wider(names_from = method, values_from = mean) %>%
    mutate(diff = MASE-KR, comparison_type = "MASE-KR", group = paste0(matrix, ":", comparison_type)) %>%
    dplyr::select(group, diff, comparison_type, trueK, sigma, m)

XU_diff <- Kresults %>%
    dplyr::select(method, matrix, trueK, sigma, m, mean) %>%
    pivot_wider(names_from = matrix, values_from = mean) %>%
    mutate(diff = X - U, comparison_type = "X-U", group = paste0(method, ":", comparison_type)) %>%
    dplyr::select(group, diff, comparison_type, trueK, sigma, m)

K_comparison <- rbind(KR_MASE_diff, XU_diff)


max_val <- max(abs(K_comparison$diff), na.rm = TRUE) # to control the color scale
p5 <- ggplot(K_comparison, aes(x = factor(sigma), y = factor(m), fill = diff)) +
    geom_tile() +
    facet_grid(trueK~group, labeller = labeller(.rows = label_both, .cols = label_value)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, 
                         limits = c(-max_val, max_val), n.breaks = 10)+
    labs(
        # title = "Comparative K Estimation Error: Cool Colors\u2192 X or MASE Better, Warm Colors \u2192 U or KR Better",
        x = "Sigma Value",
        y = "p",
        fill = "MAE (MASE - KR or X - U)"
    ) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.key.height = unit(1, "cm"),       # Increase height of each legend key
        legend.key.width = unit(1, "cm"),          # Increase width
        legend.box.margin = margin(10, 10, 10, 10) # Add margin around the legend
    )
p5
ggsave("New_Figures/m_K_XU_method.png", plot = p5, width = 12, height = 5, dpi = 300)

KR_MASE_diff_ARI <- full_df %>%
    dplyr::select(method, matrix, type, trueK, sigma, m, mean) %>%
    pivot_wider(names_from = method, values_from = mean) %>%
    mutate(diff = KR-MASE, comparison_type = "KR-MASE", group = paste0(matrix, ":", comparison_type)) %>%
    dplyr::select(group, diff, type, comparison_type, trueK, sigma, m)

XU_diff_ARI <- full_df %>%
    dplyr::select(method, matrix, type, trueK, sigma, m, mean) %>%
    pivot_wider(names_from = matrix, values_from = mean) %>%
    mutate(diff = U-X, comparison_type = "U-X", group = paste0(method, ":", comparison_type)) %>%
    dplyr::select(group, diff, type, comparison_type, trueK, sigma, m)


max_val <- max(abs(KR_MASE_diff_ARI$diff), na.rm = TRUE) # to control the color scale
p6 <- ggplot(KR_MASE_diff_ARI, aes(x = factor(sigma), y = factor(m), fill = diff)) +
    geom_tile() +
    facet_grid(rows = vars(trueK), cols = vars(type, group), labeller = labeller(.rows = label_both)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, 
                         limits = c(-max_val, max_val), n.breaks = 10)+
    labs(
        # title = "Comparative ARI: Cool Colors\u2192 MASE Better, Warm Colors \u2192 KR Better",
        x = "Sigma Value",
        y = "p",
        fill = "ARI (KR) - ARI (MASE)"
    ) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.key.height = unit(1, "cm"),       # Increase height of each legend key
        legend.key.width = unit(1, "cm"),          # Increase width
        legend.box.margin = margin(10, 10, 10, 10) # Add margin around the legend
    )
p6
ggsave("New_Figures/m_diffARI_XU_method.png", plot = p6, width = 12, height = 5, dpi = 300)

max_val <- max(abs(XU_diff_ARI$diff), na.rm = TRUE) # to control the color scale
p7 <- ggplot(XU_diff_ARI, aes(x = factor(sigma), y = factor(m), fill = diff)) +
    geom_tile() +
    facet_grid(rows = vars(trueK), cols = vars(type, group), labeller = labeller(.rows = label_both)) +
    scale_fill_gradient2(low = "#00004d", mid = "white", high = "#c31e23", midpoint = 0, 
                         limits = c(-max_val, max_val), n.breaks = 10)+
    labs(
        # title = "Comparative ARI: Cool Colors\u2192 X Better, Warm Colors \u2192 U Better",
        x = "Sigma Value",
        y = "p",
        fill = "ARI (U) - ARI (X)"
    ) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.key.height = unit(1, "cm"),       # Increase height of each legend key
        legend.key.width = unit(1, "cm"),          # Increase width
        legend.box.margin = margin(10, 10, 10, 10) # Add margin around the legend
    )
p7
ggsave("New_Figures/m_diffARI_XU.png", plot = p7, width = 12, height = 5, dpi = 300)

my_palette <- c("#fbe6c5","#ee8a82","#c8586c","#9c3f5d", "#70284a")

ARI_unknownK <- full_df %>% filter(type=="unknownK")
ARI_knownK <- full_df %>% filter(type=="knownK")
p8 <- ggplot(ARI_unknownK , aes(x = factor(sigma), y = factor(m), fill = mean)) +
    geom_tile() +
    facet_grid(rows = vars(method), cols = vars(matrix, trueK), labeller = labeller(.cols = label_both)) +
    scale_fill_gradientn(colors = my_palette) +  
    labs(
        # title = "ARI, K unknown",
        x = "Sigma Value",
        y = "p",
        fill = "Mean ARI"
    ) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.key.height = unit(1, "cm"),       # Increase height of each legend key
        legend.key.width = unit(1, "cm"),          # Increase width
        legend.box.margin = margin(10, 10, 10, 10) # Add margin around the legend
    )

p8
ggsave("New_Figures/m_ARI_unknown.png", plot = p8, width = 12, height = 5, dpi = 300)

p9 <- ggplot(ARI_knownK , aes(x = factor(sigma), y = factor(m), fill = mean)) +
    geom_tile() +
    facet_grid(rows = vars(method), cols = vars(matrix, trueK), labeller = labeller(.cols = label_both)) +
    scale_fill_gradientn(colors = my_palette) +  
    labs(
        # title = "ARI, K known",
        x = "Sigma Value",
        y = "p",
        fill = "Mean ARI"
    ) +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.1),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.key.height = unit(1, "cm"),       # Increase height of each legend key
        legend.key.width = unit(1, "cm"),          # Increase width
        legend.box.margin = margin(10, 10, 10, 10) # Add margin around the legend
    )
p9
ggsave("New_Figures/m_ARI_known.png", plot = p9, width = 12, height = 5, dpi = 300)
