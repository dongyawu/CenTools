library(RColorBrewer)
library(ggplot2)
library(Cairo)
library(cowplot)
args<-commandArgs(T)
prefix = args[1]  #Chr11

chip <- paste(prefix, "log2ratio_window2k_range500k.txt", sep="_")
anno <- paste(prefix, "AnnoTracks_window2k_range500k.txt", sep="_")


######### chip peaks
df_chip <- read.table(chip, sep="\t")
names(df_chip) <- c("chr", "start", "end", "prop", "type")
df_chip$type <- factor(df_chip$type, levels=c("B1", "B2", "box"))
df_chip$chr <- factor(df_chip$chr, levels=c("Obar NH284 Chr01: 15511100-17006245", "Ogla NH265 Chr01: 15615346-17661539", "Or4 CW03 Chr01: 17043650-18996996", "Or3 CW15 Chr01: 17853466-19388036", "GJtrp SL177 Chr01: 16781127-18554733", "GJtmp SL170 Chr01: 16306697-18014535", "GJtmp NIP Chr01: 16364227-17991793", "Or2 CW11 Chr01: 16933911-19210069", "Or1 CW06 Chr01: 16984379-19014746", "AUS SL121 Chr01: 16581032-18600791", "XI1A NJ11 Chr01: 16662389-18770468"))

######### anno tracks
df_anno <- read.table(anno, sep="\t")
names(df_anno) <- c("chr", "id", "start", "end", "y", "type")
df_anno$type <- factor(df_anno$type, levels=c("F1", "SF1","F2", "SF2","F3", "SF3","F4", "SF4","F5", "SF5","F6", "SF6","F7", "SF7","F8", "SF8","F9", "SF9","F10", "SF10","F11", "SF11","F12", "SF12","F13", "SF13","F14", "SF14","F15", "SF15", "sz22", "RIRE7", "CRM", "Interval", "box"))
df_anno$chr <- factor(df_anno$chr, levels=c("Obar NH284 Chr01: 15511100-17006245", "Ogla NH265 Chr01: 15615346-17661539", "Or4 CW03 Chr01: 17043650-18996996", "Or3 CW15 Chr01: 17853466-19388036", "GJtrp SL177 Chr01: 16781127-18554733", "GJtmp SL170 Chr01: 16306697-18014535", "GJtmp NIP Chr01: 16364227-17991793", "Or2 CW11 Chr01: 16933911-19210069", "Or1 CW06 Chr01: 16984379-19014746", "AUS SL121 Chr01: 16581032-18600791", "XI1A NJ11 Chr01: 16662389-18770468"))


plot_all <- 
  ggplot() +
  geom_rect(data=df_anno, mapping=aes(xmin=start, xmax=end, ymin=y-0.5, ymax=y+0.5, fill=type))+
  geom_line(data=df_chip, mapping=aes(x=start,y=prop, color=type), linewidth=0.1)+
  geom_rect(data=df_anno[df_anno$type == "box", ], mapping=aes(xmin=start, xmax=end, ymin=y-0.5, ymax=y+0.5), fill="transparent", color="grey", linewidth=0.05)+
  geom_rect(data=df_chip[df_chip$type == "box", ], mapping=aes(xmin=start, xmax=end, ymin=0, ymax=8), fill="transparent", color="grey", linewidth=0.05)+
  scale_fill_manual(values = c("box"="transparent", "F1"="#3f66a1", "SF1"="#3f66a1","F2"="#9D5427", "SF2"="#9D5427","F3"="#85A7CC", "SF3"="#85A7CC","F4"="#D0AF62", "SF4"="#D0AF62","F5"="#D67C1B", "SF5"="#D67C1B","F6"="#3897c5", "SF6"="#3897c5","F7"="#a874b5", "SF7"="#a874b5","F8"="#7ec0b4", "SF8"="#7ec0b4","F9"="#6DAB30", "SF9"="#6DAB30","F10"="#096858", "SF10"="#096858","F11"="#d66c54", "SF11"="#d66c54","F12"="#b13f73", "SF12"="#b13f73","F13"="#afa7d8", "SF13"="#afa7d8","F14"="#6566a9", "SF14"="#6566a9","F15"="#893f8b", "SF15"="#893f8b", "sz22"= "#5B86A8", "RIRE7"= "#AC7DB3","CRM"= "#FC8C80", "Interval"= "#F4B617"))+
  scale_color_manual(values = c("box"="grey", "B1" = "#008C8C","B2" = "#E85827", "F1"="#3f66a1", "SF1"="#3f66a1","F2"="#9D5427", "SF2"="#9D5427","F3"="#85A7CC", "SF3"="#85A7CC","F4"="#D0AF62", "SF4"="#D0AF62","F5"="#D67C1B", "SF5"="#D67C1B","F6"="#3897c5", "SF6"="#3897c5","F7"="#a874b5", "SF7"="#a874b5","F8"="#7ec0b4", "SF8"="#7ec0b4","F9"="#6DAB30", "SF9"="#6DAB30","F10"="#096858", "SF10"="#096858","F11"="#d66c54", "SF11"="#d66c54","F12"="#b13f73", "SF12"="#b13f73","F13"="#afa7d8", "SF13"="#afa7d8","F14"="#6566a9", "SF14"="#6566a9","F15"="#893f8b", "SF15"="#893f8b", "sz22"= "#5B86A8", "RIRE7"= "#AC7DB3","CRM"= "#FC8C80", "Interval"= "#F4B617"))+
  facet_wrap(~chr, ncol = 1, strip.position = "right", scales = "free_y") +
  theme(strip.text.y = element_text(angle = 0, hjust = 1), strip.placement = "outside", plot.title = element_text(hjust = 1), plot.title.position = "plot")+
  scale_y_continuous(limits = c(-6, 8), breaks = seq(-6, 8, 2))+
  scale_x_continuous(expand = c(0, 0))+
  theme(strip.background = element_blank(), panel.background = element_blank(), legend.position = "none",axis.title.x = element_blank())+
  ylab("log2(chip/input)")+
  #theme(strip.text = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(size=10))+
  ggtitle(prefix)

  

output <- paste(prefix, "monomer_LTR_INV.pdf", sep="_")

CairoPDF(file=output, width=8, height = 12)
plot_all
dev.off()

