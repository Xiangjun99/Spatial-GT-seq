
library(copykit)
library(ggplot2)
library(dplyr)

load("1708_copykit_addRNAcelltype_adt_V3_260122.rdata")

copykit

subclone_pal = c(
  'Diploid' = "#5085C4",
  'C1' = "#EB545C",
  'C2' = "#F6EC1B",
  'C3' = '#87C55F',
  'C4' = '#CC66FF')

p_spatial <-ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=ident)) + 
  scale_color_manual(values = subclone_pal)  + 
  #scale_color_gradientn(colours = c("black", "green")) + 
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  ggtitle("UMAP") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(shape = 16, size = 4.0)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(0.013, 0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(0.013, 0.013))) +
  coord_equal(xlim=c(0,73),ylim=c(73,1)) +
  theme(plot.title = element_text(hjust = 0.8, size = 25, face = "bold"),
        #axis.text=element_text(size=20),
        #axis.title=element_text(size=20,face="bold"),
        #axis.ticks =element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #axis.text.x = element_blank(), axis.text.y = element_blank(),
        #axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +   
  theme(plot.background = element_rect(fill = "white"))

p_spatial

###calculate proportion
library(dplyr)

row_prop <- test %>%
  dplyr::mutate(
    A = as.numeric(A),
    B = as.numeric(B)
  ) %>%
  dplyr::filter(ident %in% c("C1", "C2", "C3")) %>%
  dplyr::group_by(B, ident) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  dplyr::ungroup()


row_prop_wide <- row_prop %>%
  dplyr::select(B, ident, proportion) %>%
  tidyr::pivot_wider(
    names_from = ident,
    values_from = proportion,
    values_fill = 0
  )

heat_df <- row_prop_wide %>%
  dplyr::arrange(B) %>%
  tidyr::pivot_longer(
    cols = c("C1", "C2", "C3"),
    names_to = "clone",
    values_to = "proportion"
  ) %>%
  dplyr::mutate(
    B = factor(B, levels = rev(sort(unique(B)))) 
  )

#total
p_heat <- ggplot(heat_df, aes(x = clone, y = B, fill = proportion)) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("brewer_spectra")
) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  )
p_heat

p_heat <- p_heat +
  coord_fixed(ratio = 0.25) +
  scale_y_discrete(expand = c(0, 0))

library(patchwork)

p_spatial 

p_spatial + p_heat +
  plot_layout(widths = c(4, 1))

###split heatmap
heat_C1 <- heat_df %>% filter(clone == "C1")
heat_C2 <- heat_df %>% filter(clone == "C2")
heat_C3 <- heat_df %>% filter(clone == "C3")

p1 <- ggplot(heat_C1, aes(x = clone, y = B, fill = proportion)) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("brewer_red")) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_blank()) +
  scale_y_discrete(drop=FALSE) +
  coord_fixed(ratio = 0.2) 
p1

p2 <- ggplot(heat_C2, aes(x = clone, y = B, fill = proportion)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = jdb_palette("white_orange"),
    limits = c(0, 1),    
    oob = scales::squish 
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_blank()
  ) +
  scale_y_discrete(drop = FALSE) +
  coord_fixed(ratio = 0.2)  

p2

p3 <- ggplot(heat_C3, aes(x = clone, y = B, fill = proportion)) +
  geom_tile() +
  scale_fill_gradientn(colors = jdb_palette("brewer_green"),
                       limits = c(0, 1),    
                       oob = scales::squish
                       ) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_blank()) +
  scale_y_discrete(drop=FALSE)+
  coord_fixed(ratio = 0.2)
p3

p_heat <- p1 | p2 | p3 +
  #coord_fixed(ratio = 0.5) +    
  scale_y_discrete(expand = c(0, 0))


line_df <- row_prop_wide %>%
  pivot_longer(
    cols = c("C1", "C2", "C3"),
    names_to = "clone",
    values_to = "proportion",
    values_drop_na = FALSE
  ) %>%
  complete(B = 1:72, clone = c("C1","C2","C3"), fill = list(proportion = 0))
head(line_df)

line_df <- line_df %>%
  mutate(
    B = 73 - B
  )
line_df <- line_df %>%
  dplyr::filter(B >= 1 & B <= 65)

library(ggplot2)

p5 <- ggplot(line_df, aes(x = B, y = proportion, color = clone)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.5, size = 1.5) +
  scale_x_continuous(name = "Spatial Row (B)", limits = c(1,72), expand = c(0,0)) +
  scale_y_continuous(name = "Proportion of Clone", limits = c(0,1), expand = c(0,0)) +
  scale_color_manual(values = subclone_pal) +
  theme_minimal(base_size = 14) +  
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid.major = element_line(color = "grey30"), 
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "white", size = 12),
    axis.title = element_text(color = "white", size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(color="white", size=14)
  )
p5
ggsave("1708_clone_frequency_line.pdf", plot = p5, width = 8, height = 6, dpi = 300)
