library(tidyverse)
library(readxl)

total_mutations_9_29_23 <- read_excel("Downloads/total_mutations_9.29.23.xlsx")
muts <- total_mutations_9_29_23 %>% unite(Patient_ID:Tumor_ID, col = "Patient_ID", sep = "-") %>%
  pivot_longer(SNV:indel, names_to = "type", values_to = "total")  %>% filter(total>0)
muts$Patient_ID <- as.character(muts$Patient_ID)
muts$BRCA <- factor(muts$BRCA, levels = c("Wildtype", "BRCA1", "BRCA2", "BRCA1&2"))

write.table(muts, file = "POCROC_muts_10.4.23.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = T,col.names = T, 
            qmethod = c("escape", "double"),
            fileEncoding = "")


#to remove the samples that were suspect in the initial analysis
fewer_muts <- muts %>% filter_at(vars(Patient_ID), all_vars(!. %in% c("12645-886", 
                                                                      "17836-1356", 
                                                                      "17836-1399",
                                                                      "26862-2562_RtOv", 
                                                                      "26862-4860", 
                                                                      "38264-3228", 
                                                                      "40729-3360RtOv", 
                                                                      "40729-4336")))

fewer_muts$BRCA_vs <- factor(fewer_muts$BRCA_vs, levels = c("Wildtype", "Mutant"))

patient_id_colors <- c("18975" ="#0000FF", "21739" ="#FF0000", "22421" ="#00FF00",
                       "27481" ="#000033","27561"= "#FF00B6", "28601"= "#005300",
                       "41245" ="#FFD300","16030"= "#009FFF", "17687"="#9A4D42",
                       "18915"="#00FFBE", "16006"="#783FC1","25258"= "#1F9698",
                       "15815"="#FFACFD", "24892"="#B1CC71","30961"= "#F1085C",
                       "41323"="#FE8F42", "21020"="#DD00FF", "45097" ="#204254",
                       "32761"="#720055", "25236"="#766C95", "47573"="#02AD24",
                       "24487"="#C8FF00", "40729"="#886C00", "29764"="#FFB79F",
                       "48986"="#858567", "19155"="#A10300", "26862"="#14F9FF",
                       "37306"="#00479E", "17836" = "white", "12645" ="red", 
                       "34001" = "#8944a3", "38264" = "#a76b56", "22720" = "#7e8071", 
                       "16015" = "purple", "40814" = "grey")

grid_colors <- c('#7fcdbb', '#2c7fb8')

#Before/after dotplot by patient
fewer_muts %>% mutate(person = gsub("\\-.*", "", Patient_ID)) %>% 
  group_by(person, recurrence) %>% 
  mutate(type_total=sum(total)) %>% 
  ggplot() +
  geom_point(aes(x= recurrence, y= type_total, color=person)) +
  geom_line(aes(x = recurrence, y = type_total, group = person)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95)) +
  scale_fill_manual(values=grid_colors) +
  theme(strip.text.x = element_text(size = 10)) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0,50000)) +
  ggtitle("Mutation Counts") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(fill=guide_legend(title="mutation")) +
  scale_color_manual(values = patient_id_colors)

#primary only or split by BRCA status
fewer_muts %>% mutate(person = gsub("\\-.*", "", Patient_ID)) %>% group_by(person, recurrence) %>% mutate(type_total=sum(total)) %>% 
  #filter(recurrence == "Primary") %>%
  ggplot() +
  geom_point(aes(x= BRCA_vs, y= type_total, color=person)) +
  #geom_line(aes(x = BRCA_vs, y = type_total, group = person)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95)) +
  scale_fill_manual(values=grid_colors) +
  theme(strip.text.x = element_text(size = 10)) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0,50000)) +
  ggtitle("Mutation Counts") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(fill=guide_legend(title="mutation")) +
  scale_color_manual(values = patient_id_colors)

#primary only or split by BRCA status
fewer_muts %>% mutate(person = gsub("\\-.*", "", Patient_ID)) %>% group_by(person, recurrence) %>% mutate(type_total=sum(total)) %>% 
  #filter(recurrence == "Primary") %>%
  ggplot() +
  geom_point(aes(x= recurrence, y= type_total, color = person)) +
  #scale_fill_manual(values=grid_colors) +
  geom_line(aes(x = recurrence, y = type_total, group = person, color = BRCA_vs)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.95)) +
 
  theme(strip.text.x = element_text(size = 10)) +
  scale_y_continuous(expand = c(0.01, 0), limits = c(0,50000)) +
  ggtitle("Mutation Counts") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(fill=guide_legend(title="mutation")) +
  scale_color_manual(values = patient_id_colors)





load("~/Downloads/HNF4a-ESCA.rda")
redo <- Mutation_Table_Sorted_10.7.22
redo$person <- gsub("\\-.*", "", redo$Patient_ID)
redo$Tumor <- redo$recurrence
redo$Type = redo$BRCA

library(tidyr)
data <- fewer_muts

library(ggpubr)
Mutation_Table_Sorted_10.7.22 <- read.delim("~/Mutation_Table_Sorted_10.7.22.txt")
redo <- total_mutations_HRD_labels_3_17_23

library(beeswarm)
redo <- muts


# Legend
legend("topright", legend = c("Yes", "No"),
       col = 1:2, pch = 19)

redo3 <- redo[redo$recurrence == "Primary",]
redo4 <- redo[redo$recurrence == "Recurrent",]
BRCA.df <- redo[redo$BRCA_vs == "Mutant",]
WT.df <- redo[redo$BRCA_vs == "Wildtype",]

#reoder such that wildtype comes first
redo3$BRCA_vs <- factor(redo3$BRCA_vs,levels=c("Wildtype","Mutant"))
redo4$BRCA_vs <- factor(redo4$BRCA_vs,levels=c("Wildtype","Mutant"))
redo$BRCA_vs <- factor(redo$BRCA_vs,levels=c("Wildtype","Mutant"))

df <- c("#edf8fb", "#b3cde3")

beeswarm(redo$total ~ redo$BRCA_vs, 
         pch = 19,  
         pwcol = as.factor(redo$type),
         #col = rgb(red = 1, green = 0, blue = 1, alpha = 0.133),
         #pwcol = c("black", rep("grey15", 120)),
         #pwpch = c(23, rep(1, 25)),
         #pwbg = c("red", rep("transparent", 25)),
         #col = as.numeric(redo2$recurrence),
         #pwbg = c("red", rep("transparent", )),
         method="swarm", 
         cex=0.5,
         xlab="Tumor Type", 
         ylab="Mutations", 
         corral="wrap", 
         vertical = T)


legend("topright", legend = c("INDEL", "SNV"),
       col = 1:2, pch = 19, fill = guide_legend(reverse = T))


POCROC_muts_11.17.22 <- read.delim("~/POCROC_muts_3.17.23.txt")
new.df <- POCROC_muts_11.17.22
 
HRD.df <- new.df[new.df$HR_status == "HRD",]
HRP.df <- new.df[new.df$HR_status == "HRP",]
rec <- new.df[new.df$recurrence == "Recurrent",]
rec$HRD_status <- factor(rec$HRD_status,levels=c("HRP","HRD"))
prim <- new.df[new.df$recurrence == "Primary",]
prim$HRD_status <- factor(prim$HRD_status,levels=c("HRP","HRD"))
new.df$HRD_status <- factor(new.df$HRD_status,levels=c("HRP","HRD"))

colorPOCROC <- c("#fdae61","#abd9e9")

#a50026
#d73027
#f46d43
#fdae61
#fee090
#ffffbf
#e0f3f8
#abd9e9
#74add1
#4575b4
#313695

# this will plot using baseR ... unfortunately not very customizable

beeswarm(new.df$total ~ new.df$HR_status, 
         pch = 19,  
         pwcol = as.factor(new.df$type),
         #col = rgb(red = 1, green = 0, blue = 1, alpha = 0.133),
         #pwcol = c("black", rep("grey15", 120)),
         #pwpch = c(23, rep(1, 25)),
         #pwbg = c("red", rep("transparent", 25)),
         #col = as.numeric(redo2$recurrence),
         #pwbg = c("red", rep("transparent", )),
         method="swarm", 
         cex=0.5,
         xlab="Tumor Type", 
         ylab="Mutations", 
         corral="wrap", 
         vertical = T)

library(ggplot2)
library(ggbeeswarm)

# very basic ggplot
#ggplot(new.df, aes(y = total, x = HRD_status, color = type)) +
#  geom_beeswarm(cex = 3,
#                priority = "density")




ggplot(new.df, aes(y = Total, x = HR_status, color = recurrence)) +
  geom_quasirandom(cex = 3,
                   #method = "pseudorandom",
                ) +
  # removes white grid lines on gray box
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   #removes legend title 
   #theme(legend.title = element_blank()) + 
  # blows up axes titles, and legend text. does not blow up tick mark labels
  theme(axis.title = element_text(size = 20)) + theme(legend.text = element_text(size = 18)) +
  # makes white background with thin gridlines, blows up tick mark text
  theme_linedraw() + theme(axis.text = element_text(size = 18)) +
scale_color_manual(values = c("#fdae61", "#4575b4"))



#a50026
#d73027
#f46d43
#fdae61
#fee090
#ffffbf
#e0f3f8
#abd9e9
#74add1
#4575b4
#313695

library(beeswarm)
df <- data.frame(x = c(LETTERS), y = "1", 
                 z = c(rnorm(26, 11, 4)))
beeswarm(z ~ y, data = df,
         pwcol = c(1, rep(2, 25)), pwpch = c(1, rep(2, 25)), corral = "wrap", method = "center", 
         xlab = "", ylab = "variable", las=1
)

