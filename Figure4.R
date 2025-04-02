library("ggplot2")
library('dplyr')
library('tidyr')
library('ggpattern')
library('ggtext')
library('stringr')

# Load and prep data
dat <- read.csv("mgCST.csv") %>%
  drop_na() %>%
  mutate(mgCST = as.numeric(gsub("mgCST ", "", mgCST)))

dat$site <- factor(dat$site, levels = c("Vagina", "Endocervix"))
dat$Treatment <- factor(dat$Treatment, levels = c("Ct Persistence", "Ct Clearance", "No Treatment"))

# Group by for plotting
total <- dat %>%
  group_by(site, Treatment, Timepoint) %>%
  summarise(total_count = n(), .groups = 'drop')

grouped <- dat %>%
  left_join(total, by = c('site', 'Treatment', 'Timepoint')) %>%
  group_by(mgCST, site, Treatment, Timepoint) %>%
  summarise(
    mean_score = mean(mgCST.score),
    max = max(mgCST.score),
    count = n(),
    total_count = first(total_count),
    .groups = 'drop'
  ) %>%
  mutate(
    percentage = round((count / total_count) * 100, 2),
    plt = paste0(percentage, "% (", count, "/", total_count, ")")
  )

# Optional: replace mgCST numbers with scientific names
legend_labs <- c("Lactobacillus crispatus 1", "Lactobacillus crispatus 4", "Lactobacillus crispatus 6", 
                 "Lactobacillus gasseri 2", "Lactobacillus iners 1", "Lactobacillus iners 2", 
                 "Lactobacillus iners 3", "Lactobacillus iners 5", "Lactobacillus jensenii 1", 
                 "Lactobacillus jensenii 2", "Gardnerella vaginalis 1", "Gardnerella vaginalis 2", 
                 "Gardnerella vaginalis 3", "Gardnerella vaginalis 4", "Gardnerella vaginalis 5", 
                 "Bifidobacterium breve", "Bifidobacterium dentium")

i_legend_labs <- lapply(legend_labs, function(x) bquote(italic(.(x))))

# Styling maps
pattern_map <- c("Baseline" = "none", "Follow-up" = "stripe")
angle_map <- c("Baseline" = 0, "Follow-up" = 30)
custom_palette <- c("#E6194B", "#3CB44B", "#FFE119", "#0082C8", "#F58231", 
                    "#911EB4", "#46F0F0", "#F032E6", "#D2F53C", "#FABEBE", 
                    "#008080", "#E6BEFF", "#AA6E28", "#800000", "#808000", 
                    "#BFEF45", "#FFD8B1") 

# Build the plot
fig2 <- ggplot(dat) +
  geom_bar_pattern(
    data = grouped,
    aes(x = factor(mgCST), y = mean_score, fill = factor(mgCST),
        pattern = Timepoint, pattern_angle = Timepoint),
    stat = 'identity',
    position = position_dodge2(preserve = 'single', width = 0.95),
    pattern_density = 0.1, pattern_fill = "black", color = "black"
  ) +
  scale_fill_manual(
    values = custom_palette,
    labels = i_legend_labs,
    guide = guide_legend(override.aes = list(pattern = "none"))
  ) +
  geom_boxplot(
    aes(x = factor(mgCST), y = mgCST.score, group = interaction(mgCST, Timepoint)),
    outlier.shape = NA,
    position = position_dodge2(preserve = 'single', width = 0.95)
  ) +
  geom_dotplot(
    aes(x = factor(mgCST), y = mgCST.score, group = interaction(mgCST, Timepoint)),
    binaxis = 'y', stackdir = 'center', dotsize = 0.4, stackratio = 1, alpha = 0.8,
    position = position_dodge2(preserve = 'single', width = 0.95),
    show.legend = FALSE
  ) +
  # Add labels above bars
  geom_text(
    data = grouped,
    aes(x = factor(mgCST), y = -0.04, label = plt, group = Timepoint),
    size = 3,
    angle = 30,
    hjust = 1,
    vjust = 1,
    position = position_dodge(width = 0.95),
    color = "black"
  ) +
  coord_cartesian(ylim = c(-0.2, 1.0), clip = "off") +
  scale_pattern_manual(values = pattern_map,
                       guide = guide_legend(override.aes = list(fill = "white", color = "black", pattern_density = 0.05))) +
  scale_pattern_angle_manual(values = angle_map) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  facet_grid(Treatment ~ site, labeller = labeller(Treatment = c(
    "Ct Persistence" = "*Ct* Persistence",
    "Ct Clearance" = "*Ct* Clearance",
    "No Treatment" = "No Treatment Control"
  ))) +
  labs(x = 'mgCST', y = 'mgCST Score', fill = "mgSs") +
  scale_x_discrete(
    position = 'top',
    expand = expansion(mult = c(0.07, 0.05))
    ) +
  theme_bw() +
  theme(
    legend.position = 'none',
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    strip.text.y = element_markdown(size = 11),
    strip.text.x = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "pt")
  )

# Display and save
print(fig2)
ggsave("mgCST2.png", fig2, width = 20, height = 10, dpi = 300)
