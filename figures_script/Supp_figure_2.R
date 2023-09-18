  ##############################################
#                                                #
#             Supplementary figure 2             #
#                                                #
  ##############################################

  library("ggplot2")
  
  rm(list = ls())
  # General
  load("./CRC_Real/CRC_all_K849/prepare_data/data.rel.all.Rdata")
  
  Study <- c()
  Group <- c()
  seqdepth <- c()
  p.values <- c()
  for(l in 1:5){
    x <- rowSums(data.rel[[l]]$Y[data.rel[[l]]$X==1,])
    y <- rowSums(data.rel[[l]]$Y[data.rel[[l]]$X==0,])
    Group <- c(Group, rep("case", length(x)), rep("control", length(y)))
    Study <- c(Study, rep(paste0("CRC",l), length(x) + length(y)))
    seqdepth <- c(seqdepth,x,y)
    w.test <- wilcox.test(x, y, alternative = "two.sided")
    p.values <- c(p.values, w.test$p.value)
  }
  df1 <- data.frame(seqdepth = seqdepth, Group = Group, Study = Study)
  
  
  # Customizing the output
  pdf("./figures/Supp_2.pdf",         # File name
      width = 6, height = 5, # Width and height in inches
      bg = "white")   
  
  ggplot(data = df1, aes(x=Study, y=seqdepth, color = Group)) +
    geom_boxplot(position=position_dodge(0.8), width=0.6) +
    scale_color_manual(
      breaks =c("case", "control"),
      values = c("red", "blue"))  +
    ggtitle("") +
    xlab("Study") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "right",  text = element_text(size=15),
          legend.key = element_rect(fill = "white"),
          legend.title = element_blank()) 
  
  # Closing the graphical device
  dev.off()