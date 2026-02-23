# Boiler Plate


library(cowplot)
library(ggplot2)
library(tidyverse)
library(readxl)
library(dplyr)
library(ggpubr)
library(tidymodels)
library(tidytext)
library(embed)
library(ggfortify)
library(caret)
library(ggdendro)
library(tidyclust)
library(scales)
library(ggforce)
library(pheatmap)
library(ggplotify)
library(rstatix)
library(viridis)
library(ggpp)
library(patchwork)


symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

path <- getwd() 
# %>%
#   paste0("/Data Copy")

here::i_am("Src/Boilerplate.R")

# 
hc_spec <- hier_clust(
  # num_clusters = 4,
  linkage_method = "average"
)


df_merged_all <- read.csv(
  here::here("Data Copy", "All_data_merged_redone.csv"),
  check.names = FALSE
  )



cytokine_names <- colnames(df_merged_all)[3:25]
celltype_names <- colnames(df_merged_all)[28:43]
colnames_for_log <-  c(cytokine_names,celltype_names)

labels <- c("Severity", "Day_0_BW", "BAL.protein", "BAL.Cell.counts", "viral_burden_grouped","Viral.Burden","Hysteresis", "Deviation", "weights_grouped", "Sex", "Day.x"
)
severity_features <- c("Viral.Burden", "BAL.protein", "BAL.Cell.counts", "Hysteresis")


inclination <- c("Complex", "Pro", "Pro", "Pro", "Anti", "Complex", "Pro", "Anti", "Complex", "Pro", "Complex", "Pro", "Anti", "Pro", "Pro", "Pro", "Pro", "Pro", "Pro", "Pro", "Pro", "Pro", "Pro", "Complex", "Complex", "Complex", "Pro", "Anti", "Pro", "Anti", "Pro", "Pro", "Pro", "Pro", "Anti", "Pro", "Complex", "Anti", "Complex" )

lymph_myl <- c("Myeloid", "Myeloid", "Myeloid", "Lymphoid", "Both", "Myeloid", "Myeloid", "Lymphoid", "Lymphoid", "Myeloid", "Lymphoid", "Lymphoid", "Both", "Both", "Both", "Lymphoid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Both", "Both", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Lymphoid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid")

df_inclination <- data.frame(`Inflammatory nature` = factor(inclination),
                             `Lymph_Myl` = factor(lymph_myl) %>% as.numeric(),
                             `varname` = colnames_for_log)
rownames(df_inclination) <- colnames_for_log


temp_arr <- df_inclination %>% arrange(Lymph_Myl)

temp_ord <- temp_arr$varname

df_inclination_2 <- df_inclination %>%
  # filter(varname %in% colnames_for_log_2) %>% 
  select(Inflammatory.nature, Lymph_Myl)


df_male_Master_ID <- read_excel(
  here::here("Data Copy", "Exp_16_Master_ID_Male_Copy.xlsx"),
  sheet = "Data_Reordered_BW",
  # sheet = "Body weights",
  # range = "B2:H77"
)

df_female_Master_ID <- read_excel(
  here::here("Data Copy", "Exp_22_Master_ID_Female_Copy.xlsx"),
  sheet = "Data_Reordered_BW",
)


t1 <- read.csv(
  here::here("Data Copy", "Male_merged_dataset.csv")
)
t2 <- read.csv(
  here::here("Data Copy", "Female_merged_dataset.csv")
)

df_t <- bind_rows(
  t1 %>% select(Mouse_ID, Viral.Burden, BAL.Cell.counts, BAL.protein, Hysteresis) %>% 
    mutate(Sex = "Male"),
  t2 %>% select(Mouse_ID, Viral.Burden, BAL.Cell.counts, BAL.protein, Hysteresis) %>% 
    mutate(Sex = "Female",
           Mouse_ID = Mouse_ID +75),
) %>% 
  mutate(Viral.Burden = if_else(Mouse_ID %in% 76:90, 0, Viral.Burden)
  )

df_merged_all_2 <- dplyr::left_join(
  df_merged_all,
  df_t %>% select(-Sex),
  by = "Mouse_ID"
)

# replacing severity label column NAs with zeroes

df_female_bws <- df_female_Master_ID %>% 
  mutate(`Day 0` = 100) %>% 
  pivot_longer(cols = starts_with("Day "),
               names_to = "DPI",
               values_to = "BWs") %>% 
  mutate(DPI = case_when( DPI == "Day 0" ~ 0,
                          DPI == "Day 3" ~ 3,
                          DPI == "Day 4" ~ 4,
                          DPI == "Day 5" ~ 5,
                          DPI == "Day 7" ~ 7,
                          DPI == "Day 8" ~ 8,)) %>% 
  filter(BWs >0) %>% 
  mutate(BWs = BWs/100,
         Group = Days) %>% 
  select(-Days)



df_male_bws <- df_male_Master_ID %>% 
  pivot_longer(cols = starts_with("Day "),
               names_to = "DPI",
               values_to = "BWs") %>%
  mutate(DPI = case_when( DPI == "Day 0" ~ 0,
                          DPI == "Day 3" ~ 3,
                          DPI == "Day 4" ~ 4,
                          DPI == "Day 5" ~ 5,
                          DPI == "Day 7" ~ 7,
                          DPI == "Day 8" ~ 8,)) %>% 
  filter(BWs >0) %>% 
  mutate(Mouse_ID = ID) %>% 
  select(-ID)

df_bws_merged <- bind_rows( df_male_bws %>% mutate(Sex = "Male"), 
                            df_female_bws %>% mutate(Mouse_ID = Mouse_ID+75,
                                                     Sex = "Female"),
) %>% 
  mutate(control = case_when(Group == "Day 0" ~ "Control",
                             Group != "Day 0" ~ Sex)) %>% 
  filter(Mouse_ID > 15)





df_merged_all_2 <- df_merged_all_2 %>% 
  mutate(across(c("Viral.Burden", "BAL.protein", "BAL.Cell.counts", "Hysteresis"), ~ifelse(is.na(.),0,.))) %>% 
  select(-Severity)


df_merged_severity_2 <- df_merged_all_2 %>%
  # mutate(across(c("Viral.Burden", "BAL.protein", "BAL.Cell.counts"), ~ log10(1+.))) %>%
  group_by(Sex, Day.x) %>%
  filter(Day.x>0) %>% 
  # mutate(across(severity_features, ~ . - mean(.))) %>% 
  
  
  mutate(
    bool_viral = (Viral.Burden > median(Viral.Burden)),
    bool_balprotein =   (BAL.protein > median(BAL.protein)), 
    bool_balcellcount =  (BAL.Cell.counts > median(BAL.Cell.counts)),
    bool_hysteresis =  (Hysteresis > median(Hysteresis)),
    
    severe_criteria = bool_viral + bool_balprotein + bool_balcellcount + bool_hysteresis,
    ) %>% 
  ungroup() %>% 
  group_by(Sex, Day.x) %>%
  mutate(
    Severity_2 = case_when(
      severe_criteria >= 2 ~ "severe",
      severe_criteria < 2 ~ "mild",
      )
    ) %>%
  select(-severe_criteria, -starts_with("bool")) %>%
  ungroup()

df_merged_all_2 <-
  df_merged_all_2 %>% 
  left_join(
    df_merged_severity_2 %>% select(Mouse_ID, Severity_2),
    by = join_by(Mouse_ID)
  ) %>% 
  mutate(status = ifelse(Day.x == 0, "control", Severity_2)) 



# Filter and select the data for males and females
df_diff <- df_merged_all_2 %>%
  select(colnames_for_log, Sex, Mouse_ID, Severity_2, Day.x, Day_0_BW) %>%
  filter(Sex %in% c("Male", "Female"))
# mutate(across(all_of(colnames_for_log), ~ log10(0.00003 + .) ))

# Calculate the mean cytokine levels at Day.x == 0 for each sex
mean_cols <- df_diff %>%
  filter(Day.x == 0) %>%
  group_by(Sex) %>%
  summarise(across(all_of(colnames_for_log), mean, na.rm = TRUE)) 

# Normalize the cytokine levels by the mean values
df_diff <- df_diff %>%
  left_join(mean_cols, by = "Sex", suffix = c("", "_mean")) %>%
  mutate(across(all_of(colnames_for_log), ~ . / get(paste0(cur_column(), "_mean")))) %>%
  select(-ends_with("_mean")) %>% 
  mutate(across(all_of(colnames_for_log), ~ log10(0.003 + .) )) %>% 
  ungroup()

df_diff_2 <- df_diff %>% 
  filter(Day.x>0)%>% 
  left_join(
    df_merged_severity_2 %>% select(Mouse_ID,
                                    Viral.Burden,
                                    BAL.protein,
                                    BAL.Cell.counts,
                                    Hysteresis,
                                    ),
    by = "Mouse_ID"
  )

df_bws_merged_2 <- df_bws_merged %>% left_join(df_merged_severity_2 %>% select(Mouse_ID, Severity_2 ),
                                               by = "Mouse_ID")

#df_merged_all_2 should be the default dataframe for use everywhere