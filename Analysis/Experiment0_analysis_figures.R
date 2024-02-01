library(tidyverse)
library(patchwork)

data <- read_csv('/home/sam/scRNAseq/Xenium/ClassVsSubclass/Merged_Experiments_FINAL.csv')  %>%
  mutate(Class = as.factor(Class),
         Cluster = as.factor(Cluster),
         NumEpochs = as.factor(NumEpochs),
         L1Lambda = as.factor(L1Lambda),
         EarlyStopping = as.factor(EarlyStopping),
         deltaRate = as.factor(deltaRate),
         ClassWeight = as.factor(ClassWeight),
         Experiment = as.factor(Experiment))

FF_df <- data %>%
  filter(Experiment=="FF") %>%
  filter(Layer1==125*3, Layer2==125*3, Layer3==125*3)

p1 <- ggplot(FF_df, aes(x=L1Lambda, y = mean_F1_class, fill = Class, color = Class)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_F1_class - sd_F1_class, ymax = mean_F1_class + sd_F1_class), 
                width = 0.1) +
  ggtitle("FF Hyperparameter Grid Search Class Performance") +
  facet_wrap(~NumEpochs*EarlyStopping, nrow=4)+ 
  ylab("Mean F1") +
  theme_minimal()

p2 <- ggplot(FF_df, aes(x=L1Lambda, y = mean_F1, fill = Class)) +
  geom_boxplot() +
  # ggtitle("FF Hyperparameter Grid Search SubClass Performance") +
  facet_wrap(~NumEpochs*EarlyStopping, nrow=4) + 
  ylab("Mean F1") +
  theme_minimal()

p3 <- ggplot(FF_df, aes(x=L1Lambda, y = sd_F1**2, fill = Class)) +
  geom_boxplot() +
  ggtitle("FF Hyperparameter Grid Search SubClass Variance") +
  scale_y_log10() +
  facet_wrap(~NumEpochs*EarlyStopping, nrow=4) + 
  ylab("F1 Variance") +
  theme_minimal()

optimal_FF <- FF_df %>%
  filter(NumEpochs==50, EarlyStopping==10, L1Lambda==0.001)

p4.1 <- ggplot(optimal_FF, aes(x=Class, y = mean_F1_class, fill = Class, color = Class)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_F1_class - sd_F1_class, ymax = mean_F1_class + sd_F1_class), 
                width = 0.1) +
  ggtitle("Optimal FF Class Performance") +
  ylab("Mean F1") +
  theme_minimal()
  

p4.2 <- ggplot(optimal_FF, aes(x=Class, y = mean_F1, fill = Class)) +
  geom_boxplot() +
  ggtitle("Optimal FF SubClass Performance") +
  ylab("Mean F1") +
  theme_minimal()

p4.3 <- ggplot(optimal_FF, aes(x=Class, y = sd_F1**2, fill = Class)) +
  geom_boxplot() +
  ggtitle("Optimal FF SubClass Variance") +
  ylab("F1 Variance") +
  theme_minimal()

p1
p2
p3

supFig15 <- p4.1+p4.2+p4.3
ggsave('/home/sam/Poleg-Polsky/ICML/Figures/SupFig14.pdf', p2, device = "pdf", width = 20, height = 20)
ggsave('/home/sam/Poleg-Polsky/ICML/Figures/SupFig15.pdf', supFig15, device = "pdf", width = 15, height = 10)

