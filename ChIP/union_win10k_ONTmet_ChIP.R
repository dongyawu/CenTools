library(RColorBrewer)
library(ggplot2)
library(Cairo)
library(zoo)
library(dplyr)
library(cowplot
args<-commandArgs(T)  #SL170
prefix = args[1]

chipfile = paste(args[1],"log2ratio_10k.bdg",sep="_")
chip <- read.table(chipfile)
ONTmet_file = paste(args[1],"window_ONTmet_ChIPLOG2_10k.txt",sep="_")
ONTmet <- read.table(ONTmet_file)

names(chip) <- c("chr", "start", "end", "logchip", "type")
names(ONTmet) <- c("chr", "start", "end", "number", "type")
chip$type <- factor(chip$type, levels=c("B1","B2"))
ONTmet$type <- factor(ONTmet$type, levels=c("CG","CHG","CHH"))
df_newchip <- chip %>%mutate(avg_scale4 = rollmean(logchip, k=4, fill=NA, align = 'right'))
df_newONTmet <- ONTmet %>%mutate(avg_scale4 = rollmean(number, k=4, fill=NA, align = 'right'))

plot_all <- 
  ggplot()+
  geom_bar(data=df_newONTmet, mapping=aes(x=start, y=avg_scale4*6, fill=type), stat = 'identity', position = 'identity')+
  geom_line(data=df_newchip, mapping=aes(x=start, y=avg_scale4), group=1) + 
  scale_fill_manual(values = c("CG"="red", "CHG"="blue", "CHH"="green","B1"="black", "B2"="black"))+
  facet_wrap(~chr, scales = "free_x", ncol = 3)+
  scale_y_continuous(name = "log2(chip)", sec.axis = sec_axis(~./6, name = "DNAmethylation"))+
  theme(strip.background = element_blank(), panel.background = element_rect(fill=NA, color="black"), legend.position = "none")+  ### 去掉分面板的框和底纹
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=10))
  
chip <- read.table("CW15_log2chip_10k_cent1M.bdg")
names(chip) <- c("chr", "start", "end", "number", "type", "fillcolor")
df_chip <- chip %>%mutate(avg_scale4 = rollmean(number, k=4, fill=NA, align = 'right'))

monomer <- read.table("CW15_window_monomerNUM_10k_cent1M.txt")
names(monomer) <- c("chr", "start", "end", "number", "type", "fillcolor")

df_met2 <- read.table("CW15_ONTwindow_10k_cent1M.txt")
names(df_met2) <- c("chr", "start", "end", "number", "type")
df_CG2 <- df_met2[df_met2$type == "CG", ]
df_met_CG_slide2 <- df_CG2 %>%mutate(avg_scale4 = rollmean(number, k=4, fill=NA, align = 'right'))
df_CHG2 <- df_met2[df_met2$type == "CHG", ]
df_met_CHG_slide2 <- df_CHG2 %>%mutate(avg_scale4 = rollmean(number, k=4, fill=NA, align = 'right'))
df_CHH2 <- df_met2[df_met2$type == "CHH", ]
df_met_CHH_slide2 <- df_CHH2 %>%mutate(avg_scale4 = rollmean(number, k=4, fill=NA, align = 'right'))
df_all2 <- df_met2[df_met2$type == "all_met", ]
df_met_all_slide2 <- df_all2 %>%mutate(avg_scale4 = rollmean(number, k=4, fill=NA, align = 'right'))


met_count2 <- read.table("CW15_ONTcount_10k_cent1M.bdg")
names(met_count2) <- c("chr", "start", "end", "number", "type")

p1<-ggplot() +
  geom_line(data = df_chip[df_chip$fillcolor == "B1", ], aes(x = start, y = avg_scale4), color = "black") +
  facet_wrap(~ fillcolor + chr, scales = "free_x", ncol = 12)+
  theme(strip.background = element_blank(), panel.background = element_rect(fill=NA, color="black"), legend.position = "none")+ 
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=6), axis.text.x = element_blank(), strip.text = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  ylab("CENH3")


p2<-ggplot() +
  geom_line(data = df_chip[df_chip$fillcolor == "B2", ], aes(x = start, y = avg_scale4), color = "black") +
  facet_wrap(~ fillcolor + chr, scales = "free_x", ncol = 12)+
  theme(strip.background = element_blank(), panel.background = element_rect(fill=NA, color="black"), legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=6), axis.text.x = element_blank(), strip.text = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  ylab("CENH3")

p3 <- ggplot() +
  geom_bar(data = met_count2, aes(x = start, y = number, fill = type),stat = "identity")+
  scale_fill_manual(values=c("CG" = "#336600", "CHG" = "#663300", "CHH" = "#003366"))+
  facet_wrap(~ chr, scales = "free_x", ncol = 12)+
  theme(strip.background = element_blank(), panel.background = element_rect(fill=NA, color="black"), legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=6), axis.text.x = element_blank(), strip.text = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  ylab("Met count")
  

p4 <- ggplot() +
  geom_bar(data = df_met_all_slide2, aes(x = start, y = avg_scale4),  fill = "grey", stat = 'identity', position = 'identity', width = 0.02) +
  facet_wrap(~ type + chr, scales = "free_x", ncol = 12)+
  theme(strip.background = element_blank(), panel.background = element_rect(fill=NA, color="black"), legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=6), axis.text.x = element_blank(), strip.text = element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  ylab("All met")

p5 <- ggplot() +
  geom_bar(data = df_met_CG_slide2, aes(x = start, y = avg_scale4), fill = "#336600", stat = 'identity', position = 'identity', width = 0.02) +
  facet_wrap(~ type + chr, scales = "free_x", ncol = 12)+
  theme(strip.background = element_blank(), panel.background = element_rect(fill=NA, color="black"), legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=6), axis.text.x = element_blank(), strip.text = element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  ylab("CG")


p6 <- ggplot() +
  geom_bar(data = df_met_CHG_slide2, aes(x = start, y = avg_scale4),  fill = "#663300", stat = 'identity', position = 'identity', width = 0.02) +
  facet_wrap(~ type + chr, scales = "free_x", ncol = 12)+
  theme(strip.background = element_blank(), panel.background = element_rect(fill=NA, color="black"), legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=6), axis.text.x = element_blank(), strip.text = element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  ylab("CHG")

p7 <- ggplot() +
  geom_bar(data = df_met_CHH_slide2, aes(x = start, y = avg_scale4),  fill = "#003366", stat = 'identity', position = 'identity', width = 0.02) +
  facet_wrap(~ type + chr, scales = "free_x", ncol = 12)+
  theme(strip.background = element_blank(), panel.background = element_rect(fill=NA, color="black"), legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=6), axis.text.x = element_blank(), strip.text = element_blank())+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0, 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  ylab("CHH")



p8<-ggplot() + geom_bar(data = monomer, aes(x = start, y = number, fill=fillcolor), stat = 'identity', position = 'identity') +
  scale_fill_manual(values=c("monomer_zheng" = "#893F8B", "monomer_fu" = "#008C8C"))+
  facet_wrap(~ type + chr, scales = "free_x", ncol = 12)+
  theme(strip.background = element_blank(), panel.background = element_rect(fill=NA, color="black"), legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=6), axis.text.x = element_text(size=6), strip.text = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+
  ylab("Monomer")


plot_all <- plot_grid(p_met, p1,p2,p3, p4,p5, p6,p7,p8, align = "v", ncol = 1, rel_heights = c(3/19, 2/19,2/19,2/19,2/19,2/19,2/19,2/19,2/19))

output <- paste(args[1],"window_ONTmet_ChIPLOG2_10k.pdf",sep="_")

CairoPDF(file=output, width=10, height = 8)
plot_all
dev.off()
