library('dplyr')
library('ggplot2')
library('tidyr')
library('ggtext')
library('ggpattern')

summary_df <- read.csv("suppl_data_FINAL.csv")

summary_df <- summary_df %>%
  mutate(followup = ifelse(grepl('F', Sample_ID), 'Follow-up', 'Baseline'),
         site = case_when(grepl('V', Sample_ID) ~ "Vagina",
                          grepl('C', Sample_ID) ~ "Endocervix",
                          grepl('R', Sample_ID) ~ "Rectum"),
         sample_type = case_when(grepl('Pos to pos', sample_type) ~ "pos-pos",
                                 grepl('Pos to neg', sample_type) ~ "pos-neg",
                                 grepl('Neg to neg', sample_type) ~ "neg-neg")) %>%
  filter(!is.na(sample_type)) %>%
  select(-c("Paired.endocervical..vaginal.and.rectal.datasets"))

#fig 10a
df_long <- summary_df %>%
  pivot_longer(
    cols = c("ermX", "ermB", "rplV", "tetM"), 
    names_to = "Gene",
    values_to = "Presence"
  ) %>%
  mutate(
    Gene = case_when(
      Gene == "ermX" ~ "*Gv erm*X",
      Gene == "ermB" ~ "*Ln erm*B",
      Gene == "tetM" ~ "*Ln tet*(M)", 
      Gene == "rplV" ~ "*Ct rpl*V",
    )
  ) %>%
  group_by(site, sample_type, followup,  Gene) %>%
  summarise(percentage = mean(Presence)*100, total = sum(Presence), .groups = "drop")

pattern_map <- c("Baseline" = "none",
                 "Follow-up" = "stripe")

angle_map <- c("Baseline" = 0, 
               "Follow-up" = 30)

col_pal <- c("*Ct rpl*V" = "#F10C81", 
                 "*Ln erm*B" = "#468d26",
                  "*Ln tet*(M)" = "#0a589f",
                 "*Gv erm*X" = "#b7a2f3")

facet_labs <- c("*Ct* Persistence", "*Ct* Clearance", "No Treatment Control")

df_long <- df_long %>%
  filter(!total == 0)

fig_10a <- ggplot(df_long, aes(x = Gene, y = percentage, fill = Gene, pattern = followup, pattern_angle = followup)) +
  geom_bar_pattern(stat = "identity", position = position_dodge2(preserve = "single"), pattern_density = 0.1, pattern_fill = "black", color = "black") +
  facet_grid(factor(site, levels = c("Vagina", "Endocervix", "Rectum")) ~ factor(sample_type, levels = c("pos-pos", "pos-neg", "neg-neg"), labels = facet_labs)) +
  scale_pattern_manual(values = pattern_map, guide = guide_legend(override.aes = list(fill = "white", color = "black", pattern_density = 0.05))) +
  scale_pattern_angle_manual(values = angle_map) +
  scale_fill_manual(values = col_pal, guide = guide_legend(override.aes = list(pattern = "none"))) +
  labs(x = "", y = "Percentage of Samples", fill = "AMR Gene") +
  scale_x_discrete(labels=c("b" = "Baseline", "f" = "Follow-up")) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_markdown(size = 10),
        legend.title = element_markdown(size = 10),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 11, color = "black"),
        strip.text.x = element_markdown(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text.x = element_markdown(size = 10, color = "black", angle = -35),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("fig_10a.png", fig_10a, width = 8, height = 8)


#10b

df_long2 <- summary_df %>%
  pivot_longer(
    cols = c("ermX", "ermB", "rplV", "tetM"), 
    names_to = "Gene",
    values_to = "Presence"
  ) %>%
  mutate(
    Gene = case_when(
      Gene == "ermX" ~ "*Gv erm*X",
      Gene == "ermB" ~ "*Ln erm*B",
      Gene == "tetM" ~ "*Ln tet*(M)", 
      Gene == "rplV" ~ "*Ct rpl*V",
    )
  ) %>%
  group_by(site, sample_type, followup, Gene) %>%
  summarise(percentage = mean(Presence)*100, count = sum(Presence), total_count = n(), .groups = "drop")

baseline_df <- df_long2 %>%
  filter(followup == 'Baseline') %>%
  rename(baseline_count = count, baseline_total_count = total_count)

followup_df <- df_long2 %>%
  filter(followup == 'Follow-up') %>%
  rename(followup_count = count, followup_total_count = total_count)

joined_df <- baseline_df %>%
  full_join(followup_df, by = c("site", "sample_type", "Gene")) %>%
  mutate(total_count = coalesce(baseline_total_count, followup_total_count)) %>%
  select(site, sample_type, Gene, baseline_count, followup_count, total_count)

test <- joined_df %>%
  mutate(baseline_present = baseline_count,
         followup_present = followup_count,
         baseline_absent = total_count - baseline_present,
         followup_absent = total_count - followup_present) %>%
  select(site, sample_type, Gene, baseline_present, followup_present, baseline_absent, followup_absent)


results <- test %>%
  rowwise() %>%
  mutate(
    contingency_matrix = list(matrix(c(baseline_present, baseline_absent,followup_present, followup_absent), nrow = 2, byrow = TRUE)),
    test_result = list(tryCatch(mcnemar.test(contingency_matrix), error = function(e) NA)), 
    p_value = test_result$p.value) %>%
  ungroup() %>%
  select(site, sample_type, Gene, test_result, contingency_matrix, p_value)


filtered_results <- results %>%
  filter(p_value < 0.05) %>%
  select(site, sample_type, Gene, p_value)


x_filtered <- df_long2 %>%
  inner_join(filtered_results, by = c("site", "sample_type", "Gene")) %>%
  filter(!count == 0)


fig_10b <- ggplot(x_filtered, aes(x = Gene, y = percentage, fill = Gene, pattern = followup, pattern_angle = followup)) +
  geom_bar_pattern(stat = "identity", position = position_dodge2(preserve = "single"), width = 0.4, pattern_density = 0.1, pattern_fill = "black", color = "black") +
  facet_grid(factor(site, levels = c("Vagina", "Endocervix", "Rectum")) ~ factor(sample_type, levels = c("pos-pos", "pos-neg", "neg-neg"), labels = facet_labs)) +
  scale_pattern_manual(values = pattern_map, guide = guide_legend(override.aes = list(fill = "white", color = "black", pattern_density = 0.05))) +
  scale_pattern_angle_manual(values = angle_map) +
  scale_fill_manual(values = col_pal, guide = guide_legend(override.aes = list(pattern = "none"))) +
  labs(x = "", y = "Percentage of Samples", fill = "AMR Gene") +
  scale_x_discrete(labels=c("b" = "Baseline", "f" = "Follow-up")) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_markdown(size = 10),
        legend.title = element_markdown(size = 10),
        axis.title = element_text(size = 12),
        axis.text.y = element_text(size = 11, color = "black"),
        strip.text.x = element_markdown(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text.x = element_markdown(size = 10, color = "black", angle = -35),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("fig_10b.png", fig_10b, width = 8, height = 8)
