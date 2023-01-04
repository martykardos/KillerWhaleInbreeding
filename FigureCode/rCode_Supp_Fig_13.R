load("~/Documents/orca/analyses_31May2022/demograhpy/fig_16.Rdata")
# code for fig 16
g0 = ggplot(data=g0_input, aes(age, pred_invlogit, group=draws)) +
  geom_line(col = "dark blue",alpha=0.01,size=0.01) +
  geom_line(aes(age, est),col="dark blue") +
  #geom_line(col = "dark blue") +
  theme_bw() +
  ylab("Female survival") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

g1 = ggplot(g1_input, aes(age, pred_invlogit, group=draws)) +
  geom_line(col = "dark blue",alpha=0.01,size=0.01) +
  geom_line(aes(age, est),col="dark blue") +
  theme_bw() +
  ylab("Male survival") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

g2 = ggplot(g2_input, aes(age, pred_invlogit,group=draws)) +
  #geom_ribbon(aes(ymin=lo,ymax=hi),fill="#a6cef9",alpha = 0.2) +
  geom_line(col = "dark blue",alpha=0.01,size=0.01) +
  geom_line(aes(age, mean_est), col = "dark blue") +
  #geom_line(aes(Froh_raw, med), col = viridis(1), linetype=2) +
  xlim(1,60) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))+
  ylab("Fecundity") +
  xlab("Age") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# align figures
library(patchwork)
pdf("Figure_16.pdf")
g0/g1/g2
dev.off()