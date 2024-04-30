# =============================================== #
#             Supplementary figure 4              #
# =============================================== #

  # Packages ----
  library("ggplot2")

  # General ----
  rm(list = ls())
  load("./CRC/Processed_data/data.org.K849.Rdata")
  L <- length(data.rel)
  Study <- c()
  Group <- c()
  seqdepth <- c()
  p.values <- c()
  for(l in 1:L){
    x <- rowSums(data.rel[[l]]$Y[data.rel[[l]]$X==1,])
    y <- rowSums(data.rel[[l]]$Y[data.rel[[l]]$X==0,])
    Group <- c(Group, rep("case", length(x)), rep("control", length(y)))
    Study <- c(Study, rep(paste0("CRC",l), length(x) + length(y)))
    seqdepth <- c(seqdepth,x,y)
    w.test <- wilcox.test(x, y, alternative = "two.sided")
    p.values <- c(p.values, w.test$p.value)
  }
  df1 <- data.frame(seqdepth = seqdepth, Group = Group, Study = Study)

  ### Generate figures ----
  pdf("./figures/FigS4.pdf", width = 6, height = 5, bg = "white")

  ggplot(data = df1, aes(x=Study, y=seqdepth, color = Group)) +
    geom_boxplot(position=position_dodge(0.8), width=0.6) +
    scale_color_manual(
      breaks =c("case", "control"),
      values = c("red", "blue"))  +
    ggtitle("") + xlab("Study") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "right",  text = element_text(size=15),
          legend.key = element_rect(fill = "white"),
          legend.title = element_blank())

  dev.off()
