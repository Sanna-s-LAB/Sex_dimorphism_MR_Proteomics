library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(ggVennDiagram)
library(patchwork)
library(extrafont)

loadfonts(device = "pdf")

# Load data

lifelines <- read_excel("~/Supplementary_Tables.xlsx", sheet = 7)
lifelines <- subset(lifelines, select=c("exposure","outcome",
                                        "Consistency_GLGC_Lifelines_sexspecific",
                                        "Consistency_QCochran_Lifelines_sexspecific"))
lifelines <- lifelines %>% distinct()

df <- read_excel("~/Supplementary_Tables.xlsx", sheet = 5)

df <- merge(df, lifelines, all.x = TRUE)

# Create Cochran column
df$Cochran[df$Q_Cochran_pvalue < 0.05] <- "yes"
df$Cochran[is.na(df$Cochran)] <- "no"

df$Consistency_GLGC_Lifelines_sexspecific[is.na(df$Consistency_GLGC_Lifelines_sexspecific)] <- "no"
df$Consistency_QCochran_Lifelines_sexspecific[is.na(df$Consistency_QCochran_Lifelines_sexspecific)] <- "no"

###########################################################
########   FOREST PLOT OF SEX-SPECIFIC RESULTS  ###########
###########################################################

df1 <- subset(df, SexSpecific == "women-only" | SexSpecific == "men-only")

df_long <- df1 %>%
  pivot_longer(
    cols = c(b.men, b.women, b.all,
             se.men, se.women, se.all,
             fdr.men, fdr.women, fdr.all),
    names_to = c(".value", "Sex"),
    names_pattern = "(.*)\\.(men|women|all)"
  ) %>%
  mutate(
    Sex = recode(Sex,
                 men = "Men",
                 women = "Women",
                 all = "All"),
    
    ci_lower = b - 1.96 * se,
    ci_upper = b + 1.96 * se,
    
    sig = case_when(
      Sex == "Men" & SexSpecific == "men-only" ~ "Significant",
      Sex == "Women" & SexSpecific == "women-only" ~ "Significant",
      is.na(fdr) ~ "NA",
      TRUE ~ "Not significant"
    ),
    
    ExposureOutcome = paste(Assay),
    
    highlight_GLGC = if_else(Consistency_GLGC_Lifelines_sexspecific == "yes", TRUE, FALSE),
    
    symbol = case_when(
      Sex == "Men" & Consistency_GLGC_Lifelines_sexspecific == "yes" &
        SexSpecific == "men-only" ~ "*",
      Sex == "Women" & Consistency_GLGC_Lifelines_sexspecific == "yes" &
        SexSpecific == "women-only" ~ "*",
      TRUE ~ ""
    ),
    
    highlight_Cochran = if_else(Cochran == "yes", TRUE, FALSE),
    
    color_group = case_when(
      Sex == "Women" & SexSpecific == "women-only" ~ "Women_significant",
      Sex == "Men" & SexSpecific == "men-only" ~ "Men_significant",
      Sex == "All" ~ "All",
      TRUE ~ "Not_significant"
    ),
    
    color_group = factor(color_group,
                         levels = c("Women_significant",
                                    "Men_significant",
                                    "All",
                                    "Not_significant"))
  ) %>%
  group_by(outcome) %>%
  mutate(
    ExposureOutcome = factor(
      ExposureOutcome,
      levels = unique(ExposureOutcome[order(b)])
    )
  ) %>%
  ungroup()

###########################################################
# PLOTS 
###########################################################

plots <- df_long %>%
  group_split(outcome) %>%
  lapply(function(significant_sub) {
    
    significant_sub <- significant_sub %>%
      dplyr::mutate(
        ExposureOutcome = as.character(ExposureOutcome),
        ExposureOutcome = factor(
          ExposureOutcome,
          levels = rev(sort(unique(ExposureOutcome)))
        )
      )
    
    min_x <- min(significant_sub$ci_lower, na.rm = TRUE) -
      0.05 * diff(range(significant_sub$ci_lower, na.rm = TRUE))
    
    ggplot(significant_sub,
           aes(x = b, y = ExposureOutcome, color = color_group)) +
      
      geom_text(aes(x = min_x, label = symbol),
                color = "black", size = 7, vjust = 0.72) +
      
      geom_point(position = position_dodge(width = 0.6), size = 2.5) +
      
      geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                     position = position_dodge(width = 0.6),
                     height = 0.5) +
      
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      
      scale_color_manual(
        values = c(
          "Women_significant" = "#e41a1c",
          "Men_significant" = "#377eb8",
          "All" = "black",
          "Not_significant" = "gray70"
        ),
        labels = c(
          "Women_significant" = "Women (Significant)",
          "Men_significant" = "Men (Significant)",
          "All" = "All",
          "Not_significant" = "Not significant"
        )
      ) +
      
      labs(
        title = unique(significant_sub$outcome),
        x = "Effect (IVW Beta ± 95% CI)",
        y = "Exposure",
        color = ""
      ) +
      
      theme_minimal(base_size = 13) +
      theme(
        axis.text.y = element_text(size = 8),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_markdown(
          hjust = 0.5,
          size = 13
        )
      ) +
      
      guides(color = guide_legend(override.aes = list(size = 3)))
  })

###########################################################
# SAVE
###########################################################

names <- c("a","b","c","d","e")

for (i in 1:5) {
  ggsave(
    filename = paste0("~/Figure_2", names[i], ".pdf"),
    plot = plots[[i]],
    device = cairo_pdf,
    width = 7,
    height = 13,
    units = "in",
    dpi = 300
  )
}






###########################################################
#########   FOREST PLOT OF COCHRAN'S Q RESULTS  ###########
###########################################################

library(dplyr)
library(tidyr)
library(ggplot2)

df2 <- subset(df, Q_Cochran_pvalue < 0.05)

###########################################################
# LONG FORMAT (INCLUDING MEN + WOMEN + ALL)
###########################################################

df_long2 <- df2 %>%
  pivot_longer(
    cols = c(b.men, b.women, b.all,
             se.men, se.women, se.all),
    
    names_to = c(".value", "Sex"),
    names_pattern = "(.*)\\.(men|women|all)"
  ) %>%
  
  mutate(
    Sex = recode(Sex,
                 men = "Men",
                 women = "Women",
                 all = "All"),
    
    ci_lower = b - 1.96 * se,
    ci_upper = b + 1.96 * se,
    
    ExposureOutcome = paste(Assay),
    
    symbol = case_when(
      Consistency_QCochran_Lifelines_sexspecific == "yes" ~ "*",
      TRUE ~ ""
    )
  ) %>%
  
  group_by(outcome) %>%
  mutate(
    ExposureOutcome = factor(
      ExposureOutcome,
      levels = unique(ExposureOutcome[order(b)])
    )
  ) %>%
  ungroup()

###########################################################
# PLOTs
###########################################################

plots <- df_long2 %>%
  group_split(outcome) %>%
  lapply(function(df_sub) {
    
    min_x <- min(df_sub$ci_lower, na.rm = TRUE) -
      0.05 * diff(range(df_sub$ci_lower, na.rm = TRUE))
    
    ggplot(df_sub,
           aes(x = b,
               y = ExposureOutcome,
               group = Sex,
               color = Sex)) +
      
      geom_text(aes(x = min_x, label = symbol),
                color = "black",
                size = 7,
                vjust = 0.72) +
      
      geom_point(position = position_dodge(width = 0.6),
                 size = 2.5) +
      
      geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                     position = position_dodge(width = 0.6),
                     height = 0.5) +
      
      geom_vline(xintercept = 0,
                 linetype = "dashed",
                 color = "gray40") +
      
      scale_color_manual(
        values = c(
          "Men" = "#377eb8",    # blue
          "Women" = "#e41a1c",  # red
          "All" = "black"       # black reference
        )
      ) +
      
      labs(
        title = unique(df_sub$outcome),
        x = "Effect (IVW Beta ± 95% CI)",
        y = "Exposure",
        color = ""
      ) +
      
      theme_minimal(base_size = 13) +
      theme(
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
      ) +
      
      guides(color = guide_legend(override.aes = list(size = 3)))
  })

###########################################################
# SAVE 
###########################################################

names <- c("a","b","c","d","e")

for (i in 1:5) {
  ggsave(
    filename = paste0(
      "~/nas/MR_onlycis_UKBB/MR_res_proxies/Figure_3",
      names[i],
      ".pdf"
    ),
    plot = plots[[i]],
    device = cairo_pdf,
    width = 6,
    height = 8,
    units = "in",
    dpi = 300
  )
}






######################################################################
########   VENN DIAGRAM OF MALES AND FEMALES MAIN RESULTS  ###########
######################################################################

library(readxl)
library(dplyr)
library(tidyr)
library(ggVennDiagram)
library(patchwork)
library(ggplot2)

names(df)
df_women <- df %>%
  filter(SexSpecific %in% c("women-only") ) %>%
  dplyr::select(exposure, outcome) %>%
  distinct()

df_men <- df %>%
  filter(SexSpecific %in% c("men-only")) %>%
  dplyr::select(exposure, outcome) %>%
  distinct()

# Create lists of proteins splitted by outcome
sets_women <- split(df_women$exposure, df_women$outcome)
sets_men <- split(df_men$exposure, df_men$outcome)


venn_women <- ggVennDiagram(sets_women, label_alpha = 0, set_size = 4, label = "count",label_size = 5) +
  scale_fill_gradient(name = "protein\ncount", low = "#FEECEC", high = "#D73027") +
  ggtitle("Women") +
  theme(text = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5, size = 19))  # Increase title size

venn_men <- ggVennDiagram(sets_men, label_alpha = 0, set_size = 4, label = "count",label_size = 5) +
  scale_fill_gradient(name = "protein\ncount", low = "#F4FAFE", high = "#4981BF") +
  ggtitle("Men") +
  theme(text = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5, size = 19))  # Increase title size


# Space
venn_plot <- venn_women + plot_spacer() + venn_men + plot_layout(widths = c(1, 0.12, 1))
venn_plot
# Save
ggsave("~/Venn_diagram.pdf", plot = venn_plot, width = 30, height = 15, units = "cm", dpi = 300)



######################################################################
########   BAR PLOT OF MALES AND FEMALES MAIN RESULTS  ###########
######################################################################

Proteins_selected  <- df %>%
  dplyr::filter(SexSpecific %in% c("men-only", "women-only")) %>%
  dplyr::select("Assay",  "outcome", "SexSpecific") 

Proteins_plot_unique <- unique(Proteins_selected )


HDL_w <- filter(Proteins_plot_unique, outcome == "HDL" & SexSpecific == "women-only")
HDL_m <- filter(Proteins_plot_unique, outcome == "HDL" & SexSpecific == "men-only")

nonHDL_w <- filter(Proteins_plot_unique, outcome == "nonHDL" & SexSpecific == "women-only")
nonHDL_m <- filter(Proteins_plot_unique, outcome == "nonHDL" & SexSpecific == "men-only")

LDL_w <- filter(Proteins_plot_unique, outcome == "LDL" & SexSpecific == "women-only")
LDL_m <- filter(Proteins_plot_unique, outcome == "LDL" & SexSpecific == "men-only")

TC_w <- filter(Proteins_plot_unique, outcome == "TC" & SexSpecific == "women-only")
TC_m <- filter(Proteins_plot_unique, outcome == "TC" & SexSpecific == "men-only")

TG_w <- filter(Proteins_plot_unique, outcome == "TG" & SexSpecific == "women-only")
TG_m <- filter(Proteins_plot_unique, outcome == "TG" & SexSpecific == "men-only")



# Define custom outcome order
custom_order <- c("HDL", "nonHDL", "LDL", "TC", "TG")



# Prepare data
protein_counts <- Proteins_plot_unique %>%
  dplyr::filter(SexSpecific %in% c("women-only", "men-only")) %>%
  dplyr::mutate(Sex = ifelse(SexSpecific == "women-only", "Women", "Men")) %>%
  dplyr::group_by(outcome, Sex) %>%
  dplyr::summarise(Count = n(), .groups = 'drop') %>%
  mutate(outcome = factor(outcome, levels = custom_order)) 

# Create plot object
p <- ggplot(protein_counts, aes(x = outcome, y = Count, fill = Sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.7),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = c("Women" = "#D73027", "Men" = "#4981BF")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Lipid Outcome",
    y = "Protein Count",
    fill = "Sex"
  ) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Export to PDF
pdf("~/BarPlot.pdf", width = 8, height = 6)
print(p)
dev.off()




######################################################################
########   VENN DIAGRAM OF MALES AND FEMALES MAIN RESULTS ###########
########                    after MVMR                   ###########
######################################################################


df <- read_excel("~/Supplementary_Tables.xlsx", sheet = 6)

library(readxl)
library(dplyr)
library(tidyr)
library(ggVennDiagram)
library(patchwork)
library(ggplot2)

names(df)
df_women <- df %>%
  filter(SexSpecific %in% c("women-only") & (Significant==1) ) %>%
  dplyr::select(Exposure, Outcome) %>%
  distinct()

df_men <- df %>%
  filter(SexSpecific %in% c("men_only") & (Significant==1) ) %>%
  dplyr::select(Exposure, Outcome) %>%
  distinct()

# Create lists of proteins splitted by outcome
sets_women <- split(df_women$Exposure, df_women$Outcome)
sets_men <- split(df_men$Exposure, df_men$Outcome)


venn_women <- ggVennDiagram(sets_women, label_alpha = 0, set_size = 4, label = "count",label_size = 5) +
  scale_fill_gradient(name = "protein\ncount", low = "#FEECEC", high = "#D73027") +
  ggtitle("Women") +
  theme(text = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5, size = 19))  # Increase title size

venn_men <- ggVennDiagram(sets_men, label_alpha = 0, set_size = 4, label = "count",label_size = 5) +
  scale_fill_gradient(name = "protein\ncount", low = "#F4FAFE", high = "#4981BF") +
  ggtitle("Men") +
  theme(text = element_text(size = 12), 
        plot.title = element_text(hjust = 0.5, size = 19))  # Increase title size


# Space
venn_plot <- venn_women + plot_spacer() + venn_men + plot_layout(widths = c(1, 0.12, 1))
venn_plot
# Save
ggsave("~/Venn_diagram_postMVMR.pdf", plot = venn_plot, width = 30, height = 15, units = "cm", dpi = 300)



######################################################################
########   BAR PLOT OF MALES AND FEMALES MAIN RESULTS  ###########
######################################################################
names(df)
Proteins_selected  <- df %>%
  dplyr::filter(SexSpecific %in% c("men_only", "women-only") & (Significant==1)) %>%
  dplyr::select("Assay",  "Outcome", "SexSpecific") 

Proteins_plot_unique <- unique(Proteins_selected )


HDL_w <- filter(Proteins_plot_unique, Outcome == "HDL" & SexSpecific == "women-only")
HDL_m <- filter(Proteins_plot_unique, Outcome == "HDL" & SexSpecific == "men_only")

nonHDL_w <- filter(Proteins_plot_unique, Outcome == "nonHDL" & SexSpecific == "women-only")
nonHDL_m <- filter(Proteins_plot_unique, Outcome == "nonHDL" & SexSpecific == "men_only")

LDL_w <- filter(Proteins_plot_unique, Outcome == "LDL" & SexSpecific == "women-only")
LDL_m <- filter(Proteins_plot_unique, Outcome == "LDL" & SexSpecific == "men_only")

TC_w <- filter(Proteins_plot_unique, Outcome == "TC" & SexSpecific == "women-only")
TC_m <- filter(Proteins_plot_unique, Outcome == "TC" & SexSpecific == "men_only")

TG_w <- filter(Proteins_plot_unique, Outcome == "TG" & SexSpecific == "women-only")
TG_m <- filter(Proteins_plot_unique, Outcome == "TG" & SexSpecific == "men_only")



# Define custom outcome order
custom_order <- c("HDL", "nonHDL", "LDL", "TC", "TG")



# Prepare data
protein_counts <- Proteins_plot_unique %>%
  dplyr::filter(SexSpecific %in% c("women-only", "men_only")) %>%
  dplyr::mutate(Sex = ifelse(SexSpecific == "women-only", "Women", "Men")) %>%
  dplyr::group_by(Outcome, Sex) %>%
  dplyr::summarise(Count = n(), .groups = 'drop') %>%
  mutate(outcome = factor(Outcome, levels = custom_order)) 

# Create plot object
p <- ggplot(protein_counts, aes(x = Outcome, y = Count, fill = Sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.7),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = c("Women" = "#D73027", "Men" = "#4981BF")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Lipid Outcome",
    y = "Protein Count",
    fill = "Sex"
  ) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Export to PDF
pdf("~/BarPlot_afterMVMR.pdf", width = 8, height = 6)
print(p)
dev.off()
