library(MCMCglmm)
library(readxl)
library(tidyverse)

setwd("Z:/Fly_infections")
data <- read_excel("immune_genes_summary.tsv.xlsx")

data <- data[,c("Hcam", "Hconf", "Sdef", "HOG", "FB_ID", "Receptor", "Signalling", "Effector")]

# convert to long format
data_long <- data %>%
  pivot_longer(cols = c(Hcam, Hconf, Sdef), 
               names_to = "species", 
               values_to = "present") %>%
  mutate(
    present = ifelse(present == "Yes", 1, 0),  # convert to binary
    geneID = HOG
  )

data_long <- data_long %>%
  mutate(
    category = case_when(
      Receptor == "Yes" ~ "Receptor",
      Signalling == "Yes" ~ "Signalling",
      Effector == "Yes" ~ "Effector",
      TRUE ~ "Unknown"
    )
  )

data_long <- as.data.frame(data_long)
data_long$geneID <- as.factor(data_long$geneID)
data_long$species <- as.factor(data_long$species)
data_long$category <- as.factor(data_long$category)

prior <- list(R = list(V = 1, fix = 1),
              G = list(G1 = list(V = 1, nu = 0.002)))

# Fit model
model <- MCMCglmm(
  present ~ species * category,
  random = ~ geneID,
  data = data_long,
  family = "categorical",
  prior = prior,
  nitt = 10000000, burnin = 1000000, thin = 1000
)

summary(model)
plot(model$Sol)

library(tidyverse)

species_levels <- c("Hcam", "Hconf", "Sdef")
category_levels <- c("Effector", "Receptor", "Signalling", "Unknown")

pred_grid <- expand.grid(species = species_levels,
                         category = category_levels,
                         stringsAsFactors = FALSE)
pred_grid$species <- factor(pred_grid$species, levels = species_levels)
pred_grid$category <- factor(pred_grid$category, levels = category_levels)

X <- model.matrix(~ species * category, data = pred_grid)
posterior_samples <- model$Sol 


# Multiply posterior_samples %*% t(X) to get log-odds predictions
log_odds <- as.matrix(posterior_samples) %*% t(X)

# log odds to predicted probability
probs <- plogis(log_odds)

probs_df <- as.data.frame(t(probs))

# Add metadata
probs_df <- cbind(pred_grid, probs_df)


summary_df <- probs_df %>%
  pivot_longer(cols = starts_with("V"), names_to = "draw", values_to = "prob") %>%
  group_by(species, category) %>%
  summarise(
    mean_prob = mean(prob),
    lower_95 = quantile(prob, 0.025),
    upper_95 = quantile(prob, 0.975),
    .groups = "drop"
  )


species_summary <- probs_df %>%
  pivot_longer(cols = starts_with("V"), names_to = "draw", values_to = "prob") %>%
  group_by(species, draw) %>%
  summarise(prob = mean(prob), .groups = "drop") %>%
  group_by(species) %>%
  summarise(
    mean_prob = mean(prob),
    lower_95 = quantile(prob, 0.025),
    upper_95 = quantile(prob, 0.975)
  )

category_summary <- probs_df %>%
  pivot_longer(cols = starts_with("V"), names_to = "draw", values_to = "prob") %>%
  group_by(category, draw) %>%
  summarise(prob = mean(prob), .groups = "drop") %>%
  group_by(category) %>%
  summarise(
    mean_prob = mean(prob),
    lower_95 = quantile(prob, 0.025),
    upper_95 = quantile(prob, 0.975)
  )


fileConn <- file("Immune_repetoire_model_results.txt", open = "w")

writeLines(capture.output(summary(model)), fileConn)

writeLines(capture.output(summary_df), fileConn)

writeLines(capture.output(species_summary), fileConn)

writeLines(capture.output(category_summary), fileConn)

close(fileConn)

# Pivot to long format
probs_long <- probs_df %>%
  pivot_longer(
    cols = starts_with("V"), 
    names_to = "draw", 
    values_to = "prob"
  )


box_plot <- ggplot(probs_long, aes(x = category, y = prob, fill = species)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6,
               outlier.size = 0.5, lwd = 0.3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(ylim = c(0.97, 1.00)) +
  ylab("Gene recovery probability") +
  xlab("Immune gene category") +
  scale_fill_brewer(palette = "Set2", name = "Species") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1))

ggsave("Immune_gene_recov_prob.svg", plot = box_plot, width = 8.27/1.5, height = 11.69/3, units = "in")

