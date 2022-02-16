library(tidyverse)
# The simulations consist of 3 estimates of the mean, c(muhat1, muhat2, muhat3), each of
# dimension 3.
# muhat1: naive imputation by the mean
# muhat2: known missingness mechanism
# muhat3: unknown missingness mechanism: needs to be estimated
load("./R/sims/sim_MCAR.Rdata")
mcar <- sim_MCAR
rm(sim_MCAR)

load("./R/sims/sim_trunc.Rdata")
mnar <- sim_trunc
rm(sim_trunc)

mcar <- mcar[,c(1,4,7)]
mcar <- as_tibble(mcar)
colnames(mcar) <- c("naive", "known", "unknown")

mnar <- mnar[,c(1,4,7)]
mnar <- as_tibble(mnar)
colnames(mnar) <- c("naive", "known", "unknown")

mcar <- mcar %>% mutate(missType="mcar")
mnar <- mnar %>% mutate(missType="mnar")

sims <- rbind(mnar, mcar)

sims <- sims %>% pivot_longer(1:3, names_to ="estimator", values_to="bias")
sims$estimator <- factor(sims$estimator, levels = c("naive", "known", "unknown"), labels=c("Mean Impute", "Known Mechanism", "Estimated Mechanism"))
sims$missType <- factor(sims$missType, levels=c("mcar", "mnar"), labels=c("MCAR", "Truncated"))

sims1 <- subset(sims, estimator=="Mean Impute")
sims2 <- subset(sims, estimator=="Mean Impute" | estimator=="Known Mechanism")
sims3 <- sims

plot_sim <- function(sim){
  ggplot(sim, aes(x=missType, y=bias)) +
  geom_boxplot() +
  xlab("Missingness Mechanism")+
  ylab("Bias")+
  geom_hline(yintercept=0, col="red") +
  facet_wrap(~ estimator) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.title.x = element_text(margin=margin(t=20))
        ) +
  theme(axis.text = element_text(size = 12))+
  theme(strip.text.x = element_text(size = 15))
}
p1 <- plot_sim(sims1)
p2 <- plot_sim(sims2)
p3 <- plot_sim(sims3)

ggsave(filename="./../Presentation/images/plot_bias_1.jpg", plot=p1, width=16,  height=12, units="cm")
ggsave(filename="./../Presentation/images/plot_bias_2.jpg", plot=p2, width=20, height=12, units="cm")
ggsave(filename="./../Presentation/images/plot_bias_3.jpg", plot=p3, width=24, height=12, units="cm")
