library("tidyverse")


###Format Long

PanelB_long <- PanelB_measurements_all_new_areas %>% 
  pivot_longer(names_to = "IHC2",
               values_to = "Result",
               -(1:6))

View(PanelB_long)

write_csv(PanelB_long, "data_output/PanelB_long.CSV")

PanelA_measurements_all_new_areas_redo_csv %>% paste(`% AF555+: AF750+: AF488+`)


PanelA_long <- PanelA_measurements_all_new_areas_redo_csv %>% 
  pivot_longer(names_to = "IHC2",
               values_to = "Result",
               -(1:7))

View(PanelA_long)

write_csv(PanelA_long, "data_output/PanelA_long.CSV")


PanelB_long_redo <- PanelB_redo_threshold %>% 
  pivot_longer(names_to = "IHC2",
               values_to = "Result",
               -(1:3))

View(PanelB_long_redo)

write_csv(PanelB_long_redo, "data_output/PanelB_long_redo.CSV")

### Graph


Multiplex_all_data_long_format %>% 
  filter(IHC %in% c("CD15", "CD20", "CD117")) %>% 
  group_by(Classification, IHC) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
            y = Result)) +
  geom_boxplot() +
  facet_wrap(~ IHC, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))


Multiplex_all_data_long_format %>% 
  filter(IHC %in% c("CD15", "CD20", "CD117")) %>% 
  group_by(Phase, IHC) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Phase,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ IHC, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))


Multiplex_all_data_long_format %>% 
  filter(IHC %in% c("CD15", "CD20", "CD117")) %>% 
  group_by(Type, IHC) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Type,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ IHC, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))



Multiplex_all_data_long_format %>% 
  group_by(Phase, Cells_group) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Phase,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ Cells_group, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))


Multiplex_all_data_long_format %>% 
  group_by(Type, Cells_group) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Type,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ Cells_group, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))



Multiplex_all_data_long_format %>% 
  filter(Cells_group != "Neutrophils")
  group_by(Type, Cells_group) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Cells_group,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ Type, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))

  
  Multiplex_all_data_long_format %>% 
    filter(Cells_group == "Neutrophils") %>% 
  group_by(Type) %>% 
    summarise(Result) %>% 
    ggplot(aes(x = Type,
               y = Result)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
    labs(title = "Neutrophils distribution",
         x = "Groups of patients",
         y = "Neutrophils (%)")
  
  
  Multiplex_all_data_long_format %>% 
    filter(Cells_group == "B_cells") %>% 
    group_by(Type) %>% 
    summarise(Result) %>% 
    ggplot(aes(x = Type,
               y = Result)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
    labs(title = "B_cells distribution",
         x = "Groups of patients",
         y = "B_cells (%)")
  
  
  Multiplex_all_data_long_format %>% 
    filter(Cells_group == "T_cells") %>% 
    group_by(Type) %>% 
    summarise(Result) %>% 
    ggplot(aes(x = Type,
               y = Result)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
    labs(title = "T_cells distribution",
         x = "Groups of patients",
         y = "T_cells (%)")
  
  Multiplex_all_data_long_format %>% 
    filter(Cells_group == "Monocytes") %>% 
    group_by(Type) %>% 
    summarise(Result) %>% 
    ggplot(aes(x = Type,
               y = Result)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
    labs(title = "Monocytes distribution",
         x = "Groups of patients",
         y = "Monocytes (%)")
  
  
  Multiplex_all_data_long_format %>% 
    filter(Cells_group == "M2_Macrophages") %>% 
    group_by(Type) %>% 
    summarise(Result) %>% 
    ggplot(aes(x = Type,
               y = Result)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
    labs(title = "M2_Macrophages distribution",
         x = "Groups of patients",
         y = "M2_Macrophages (%)")
  
  
  Multiplex_all_data_long_format %>% 
    filter(Cells_group == "M1_Macrophages") %>% 
    group_by(Type) %>% 
    summarise(Result) %>% 
    ggplot(aes(x = Type,
               y = Result)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
    labs(title = "M1_Macrophages distribution",
         x = "Groups of patients",
         y = "M1_Macrophages (%)")
  
  
  Multiplex_all_data_long_format %>% 
    filter(Cells_group == "NK_cells") %>% 
    group_by(Type) %>% 
    summarise(Result) %>% 
    ggplot(aes(x = Type,
               y = Result)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
    labs(title = "NK_cells distribution",
         x = "Groups of patients",
         y = "NK_cells (%)")
  
  
  Multiplex_all_data_long_format %>% 
    filter(Cells_group == "Dendritic_cells") %>% 
    group_by(Type) %>% 
    summarise(Result) %>% 
    ggplot(aes(x = Type,
               y = Result)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
    labs(title = "Dendritic_cells distribution",
         x = "Groups of patients",
         y = "Dendritic_cells (%)")
  
  Multiplex_all_data_long_format %>% 
    filter(Cells_group == "Helper1_T_cells") %>% 
    group_by(Type) %>% 
    summarise(Result) %>% 
    ggplot(aes(x = Type,
               y = Result)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
    labs(title = "Helper1_T_cells distribution",
         x = "Groups of patients",
         y = "Helper1_T_cells (%)")
  
  
  Multiplex_all_data_long_format %>% 
    filter(Cells_group == "Helper2_T_cells") %>% 
    group_by(Type) %>% 
    summarise(Result) %>% 
    ggplot(aes(x = Type,
               y = Result)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
    labs(title = "Helper2_T_cells distribution",
         x = "Groups of patients",
         y = "Helper2_T_cells (%)")
  
  Multiplex_all_data_long_format %>% 
    filter(Cells_group == "Helper2_T_cells") %>% 
    group_by(Type) %>% 
    summarise(Result) %>% 
    ggplot(aes(x = Type,
               y = Result)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
    labs(title = "Helper2_T_cells distribution",
         x = "Groups of patients",
         y = "Helper2_T_cells (%)")
  
  
  Multiplex_all_data_long_format %>% 
  group_by(Type, Classification, Cells_group) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ Cells_group, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))


Multiplex_all_data_long_format %>% 
  group_by(Classification, Cells_group) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Cells_group,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ Classification, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))

Multiplex_all_data_long_format %>% 
  group_by(Classification, Cells_subgroup) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Cells_subgroup,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ Classification, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))


Multiplex_all_data_long_format %>% 
  filter(!is.na(Cells)) %>% 
  group_by(Phase, Cells) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Phase,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ Cells, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))



Multiplex_all_data_long_format %>% 
  filter(Type == "Adenomyosis") %>% 
  filter(Cells_group == "M2_Macrophages") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))


Multiplex_all_data_long_format %>% 
  filter(Type == "Transition") %>% 
  filter(Cells_group == "M2_Macrophages") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))


Multiplex_all_data_long_format %>% 
  filter(Classification %in% c("Control_eutopic", "Transition_eutopic", "AD_eutopic")) %>% 
  filter(Cells_group == "M2_Macrophages") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))


Multiplex_all_data_long_format %>% 
  group_by(IHC, Type) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Type,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ IHC, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90))


Multiplex_all_data_long_format %>% 
  filter(IHC %in% c("CD117", "CD15", "CD20", "CD1a", "CD3_GATA3", "CD3_CD8", 
                    "CD3", "CD68", "CD3_Tbet", "CD68_CD86", "CD68_CD163", "NKp46")) %>% 
  group_by(IHC, Type) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Type,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ IHC, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Inflammatory cells distribution according to diagnosis",
       x = "Diagnosis",
       y = "Inflammatory cells (%)")


Multiplex_all_data_long_format %>% 
  filter(IHC %in% c("CD117", "CD15", "CD20", "CD1a", "CD3_GATA3", "CD3_CD8", 
                    "CD3", "CD68", "CD3_Tbet", "CD68_CD86", "CD68_CD163", "NKp46")) %>% 
  filter(Type == "Adenomyosis") %>% 
  group_by(IHC, Phase) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Phase,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ IHC, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of inflammatory cells in adenomyosis patients, according to the phase of the cylcle",
       x = "Phase of the cycle",
       y = "Inflammatory cells (%)")


Multiplex_all_data_long_format %>% 
  filter(IHC %in% c("CD117", "CD15", "CD20", "CD1a", "CD3_GATA3", "CD3_CD8", 
                    "CD3", "CD68", "CD3_Tbet", "CD68_CD86", "CD68_CD163", "NKp46")) %>% 
  filter(Type == "Transition") %>% 
  group_by(IHC, Phase) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Phase,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ IHC, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of inflammatory cells in transition patients, according to the phase of the cylcle",
       x = "Phase of the cycle",
       y = "Inflammatory cells (%)")

Multiplex_all_data_long_format %>% 
  filter(IHC %in% c("CD117", "CD15", "CD20", "CD1a", "CD3_GATA3", "CD3_CD8", 
                    "CD3", "CD68", "CD3_Tbet", "CD68_CD86", "CD68_CD163", "NKp46")) %>% 
  filter(Type == "Control") %>% 
  group_by(IHC, Phase) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Phase,
             y = Result)) +
  geom_boxplot() +
  facet_wrap(~ IHC, scales = "free_y") +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of inflammatory cells in controls, according to the phase of the cylcle",
       x = "Phase of the cycle",
       y = "Inflammatory cells (%)")

Multiplex_all_data_long_format %>% 
  filter(IHC == "CD3") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of T cells",
       x = "Phase of the cycle",
       y = "T cells (%)")


Multiplex_all_data_long_format %>% 
  filter(IHC == "CD68") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of macrophages",
       x = "Phase of the cycle",
       y = "Macrophages (%)")

Multiplex_all_data_long_format %>% 
  filter(IHC == "CD15") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of neutrophils",
       x = "Phase of the cycle",
       y = "neutrophils (%)")

Multiplex_all_data_long_format %>% 
  filter(IHC == "NKp46") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of NK cells",
       x = "Phase of the cycle",
       y = "NK cells (%)")



Multiplex_all_data_long_format %>% 
  filter(IHC == "CD68_CD86") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of M1 macrophages",
       x = "Phase of the cycle",
       y = "M1 macrophages (%)")
  

Multiplex_all_data_long_format %>% 
  filter(IHC == "CD3_CD8") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of cytotoxic T cells",
       x = "Phase of the cycle",
       y = "Cytotoxic T cells (%)")


Multiplex_all_data_long_format %>% 
  filter(IHC == "CD3_GATA3") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of Helper 2 T cells",
       x = "Phase of the cycle",
       y = "Helper 2 T cells (%)")

Multiplex_all_data_long_format %>% 
  filter(IHC == "CD117") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of mast cells",
       x = "Phase of the cycle",
       y = "Mast cells (%)")


Multiplex_all_data_long_format %>% 
  filter(IHC == "CD68_CD163") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of M2 macrophages",
       x = "Phase of the cycle",
       y = "M2 macrophages (%)")


Multiplex_all_data_long_format %>% 
  filter(IHC == "CD20") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of B cells",
       x = "Phase of the cycle",
       y = "B cells (%)")


Multiplex_all_data_long_format %>% 
  filter(IHC == "CD1a") %>% 
  group_by(Classification) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  labs(title = "Distribution of dendritic cells",
       x = "Phase of the cycle",
       y = "dendritic cells (%)")


Multiplex_all_data_long_format %>% 
  filter(Classification %in% c("AD_myometrium", "Control_myometrium", "Transition_myometrium")) %>% 
  filter(IHC %in% c("CD117", "CD15", "CD20", "CD1a", "CD3_GATA3", "CD3_CD8", 
                    "CD3", "CD68", "CD3_Tbet", "CD68_CD86", "CD68_CD163", "NKp46")) %>% 
  group_by(Classification, IHC) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  facet_wrap(~ IHC, scales = "free_y")


Multiplex_all_data_long_format %>% 
  filter(Classification %in% c("AD_eutopic", "Control_eutopic", "Transition_eutopic")) %>% 
  filter(IHC %in% c("CD117", "CD15", "CD20", "CD1a", "CD3_GATA3", "CD3_CD8", 
                    "CD3", "CD68", "CD3_Tbet", "CD68_CD86", "CD68_CD163", "NKp46")) %>% 
  group_by(Classification, IHC) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  facet_wrap(~ IHC, scales = "free_y")


Multiplex_all_data_long_format %>% 
  filter(Classification %in% c("AD_ectopic_transition", "AD_ectopic_transition_dilated", "Transition_ectopic", "Transition_ectopic_dilated")) %>% 
  filter(IHC %in% c("CD117", "CD15", "CD20", "CD1a", "CD3_GATA3", "CD3_CD8", 
                    "CD3", "CD68", "CD3_Tbet", "CD68_CD86", "CD68_CD163", "NKp46")) %>% 
  group_by(Classification, IHC) %>% 
  summarise(Result) %>% 
  ggplot(aes(x = Classification,
             y = Result)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size = rel(1.5), angle = 90)) +
  facet_wrap(~ IHC, scales = "free_y")







