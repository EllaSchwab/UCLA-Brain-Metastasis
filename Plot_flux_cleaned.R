#Author: Ella Schwab
#Written December 2022
#This script was used to read in raw bioluminescence output from an IVIS in vivo imaging system,
#from a pilot study to determine efficacy of a novel small molecule for brain metastasis in mice.
#Mice were injected with cancer cell lines through intracardiac injection and then bioluminescence was recorded each week 
#This script cleans, processes, and plots the data
#BLI (bioluminescence in the brain), mouse weights, and weekly survival plots are included


install.packages(c('dplyr', 'readxl','tidyverse', 
                   'installr', 'extrafont', 'plyr',
                   'swimplot', 'lubridate', 'survminer', 'gridExtra','patchwork',dependencies = TRUE))


#load multiple packages at once
my_packages <- c('dplyr', 'readxl','tidyverse',
                 'installr', 'extrafont', 'ggplot2', 'plyr', 
                 'swimplot', 'lubridate', 'survminer', 'gridExtra', 'patchwork')  
lapply(my_packages, require, character.only = TRUE) 

###############################################################################

filenames <- list.files("Bioluminescence")  # Identify day file names
filenames

#read in all the files as a separate data frame
for(j in 1:length(filenames)) { 
  #assigns a name (substring of filenames names) to each new file being read in
  assign(substr(filenames[j],1,nchar(filenames[j])),  #from 1 to nchar (the last char)   
         read.delim(paste0("Bioluminescence/", filenames[j]),
                    header = TRUE, sep = "\t", dec = ".", check.names = FALSE))
  
} 

#make a list of dataframes that end with .txt (all dataframes that are text files)
dflist <- mget(ls(pattern = "\\.txt$"))

#rbind and combine all the dataframes  
CTable <- do.call("rbind", dflist)

#extract the date from the original text files
CTable$Date <- substr(CTable$`Date and Time`,1,10)
#format as a date
CTable$Date <- as.Date(CTable$Date , format = "%m/%d/%Y") #ex: 09/25/2022

#keep only certain columns
CTable <- subset(CTable, 
                 select = c("ROI Label", "Total Flux [p/s]", "Date", "Image Number"))

#change the name of "ROI Label" to "ID"
colnames(CTable)[1] ="ID"

#read in file linking mouse ID to group
GroupID <- read.csv(file = 'GroupIDs.csv', header = TRUE)
CTable <- merge(CTable, GroupID, by = "ID")

#organize df by date
CTable <- CTable %>% arrange(CTable$Date)

#delete obs for XX from 11/09/22
CTable <- CTable[!(CTable$ID == "XX" & CTable$Date == "2022-11-09"),]

#change StartTrt to date 
CTable$StartTrt <-as.Date(CTable$StartTrt, format = "%m/%d/%Y")

#calculate week based on date
start = CTable$StartTrt #initialize the start date of the study

#check that 11/09/22 is correctly input as week 0
floor(as.numeric(difftime(CTable$Date[220], CTable$StartTrt[220], units = "weeks")))

#iterate over every row in the date column
#ceiling rounds to the nearest week (upper)
#difftime find the number of weeks from the start to date at i 
#as numeric converts the weird difftime output to a single number
for(i in 1:length(CTable$Date)){
  CTable$Week[i] <- floor(as.numeric(difftime(CTable$Date[i],start[i], units = "weeks"))) 
  
}

#Grab highest BLI values for each week
CTable <- CTable %>%
  arrange(ID,desc(`Total Flux [p/s]`)) %>%
  group_by(ID, Week) %>%
  slice(1) #takes the the highest x values (x=1 to take the first row i.e. the highest value)


#make a separate df for baseline values
BaseTable <- CTable[CTable$Group %in% "base",]
#calculate avereage baseline flux from BaseTable and bind it back to CTable
meanbase<-  aggregate(`Total Flux [p/s]` ~ Week, BaseTable, mean) #calc mean flux
colnames(meanbase)[2] <- "MeanBaseFlux"

CTable <- merge(CTable, meanbase, by = "Week")
#adjust the baseline values
CTable$adjflux <- CTable$`Total Flux [p/s]` - CTable$MeanBaseFlux
#remove baseline values from working dataset
CTable <- CTable[CTable$Group != "base",]


#calculate group means and sd by week 
data_msd <- CTable %>%  # Get mean & standard deviation by group
  group_by(Week, Group) %>%
  summarise_at(vars(`Total Flux [p/s]`),
               list(rawmean = mean,
                    rawsd = sd)) %>% 
  as.data.frame(data_msd)

CTable <- merge(CTable, data_msd, by = c("Week", "Group"))

###############################################################################
#plot the flux plot in a presentable way

windowsFonts(Arial=windowsFont("Arial"))
fluxplot<-ggplot(CTable, aes(x = Week, y = `Total Flux [p/s]`, color = Group)) + 
  geom_jitter( position = position_jitter(0.05), size = 1.5) +
  geom_line(data = CTable, aes(x = Week, y= rawmean), size = 1)+
  geom_errorbar(data= CTable, aes(x = Week, ymin= rawmean, ymax= rawmean), 
                position = position_dodge(0), width=0.5, size = 1.5,) + 
  scale_color_manual(labels = c("Group 1: Control", "Group 2: XXXX XXXmg/kg", "Group 3: XXXX XXXmg/kg", "Group 4: XXXXXX XXXmg/kg"), 
                     name = "Treatment Group",
                     values = c( "red2", "dodgerblue3", "darkseagreen3", "darkorange"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+ 
 
   labs(y= "Bioluminescence (photons/s)", x = "Week") +
  
  #limits to shrink distance
  scale_x_continuous(breaks= seq(0,max(CTable$Week), 1)) +
  scale_y_log10(limits = c(1e5, 1e10), 
               breaks = c(1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11) ,expand = c(0.003, 0)) +
  theme(axis.text.x=element_text(family="sans", size=9, color = "black"), 
        axis.text.y=element_text(family="sans", size=9, color = "black"),
        axis.title.x = element_text(family="sans", face="bold", size=11, color = "black"),
        axis.title.y = element_text(family="sans", face="bold", size=11, color = "black"),
        legend.key=element_blank()) 

fluxplot


####PLOT EACH MOUSE BLI SEPARATELY##########################################

cbp2 <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
          "#F39B7F", "#B254A5", "#91D1C2", "#D7A7C7")

#make separate datasets for each group of mice
BLIGroup_1<- CTable[CTable$Group == "1",]
BLIGroup_2<- CTable[CTable$Group == "2",]
BLIGroup_3<- CTable[CTable$Group == "3",]
BLIGroup_4<- CTable[CTable$Group == "4",]

#function to plot each mouse in an individual plot
BLImouseplot <- function(d, t){
  ggplot(d, aes(x = Week, y = `Total Flux [p/s]`, color = ID)) + 
    geom_point(size = 2) +   
    geom_line(data = d, aes(x = Week, y= `Total Flux [p/s]`, group = ID), size = 1)+
    geom_vline(aes(xintercept = 12), linetype = "dotted") + 
    labs(y= "Bioluminescence (photons/s)", x = "Week", color = "Mouse ID") +
    ggtitle(t) +
    scale_colour_manual(values=cbp2) + 
    #limits to shrink distance
    scale_x_continuous(breaks= seq(0,max(CTable$Week), 1))+
    scale_y_log10(limits = c(1e5, 1e10), 
                  breaks = c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11) ,expand = c(0.003, 0)) +
    theme(axis.text.x=element_text(family="sans", size=9, color = "black"), 
          axis.text.y=element_text(family="sans", size=9, color = "black"),
          axis.title.x = element_text(family="sans", face="bold", size=11, color = "black"),
          axis.title.y = element_text(family="sans", face="bold", size=11, color = "black"),
          legend.key=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) 
    
    
}

BG1m<- BLImouseplot(BLIGroup_1,"Group 1: Control")
BG2m<- BLImouseplot(BLIGroup_2, "Group 2: XXXX XXXmg/kg")
BG3m<- BLImouseplot(BLIGroup_3,  "Group 3: XXXX XXXmg/kg")
BG4m<- BLImouseplot(BLIGroup_4, "Group 4: XXXX XXXmg/kg")



grid.arrange(arrangeGrob(BG1m, BG2m, BG3m, BG4m, ncol = 2)) 
###################################################################################
#Weight plot

weights <- read.csv(file = 'weight12.csv', header = T)

weights$fGroup <- as.factor(weights$Group)
weights$xDate <-as.Date(weights$Date, format = "%m/%d/%Y")
weights$xICDate  <-as.Date(weights$ICDate, format = "%m/%d/%Y")
weights <-na.omit(weights)


#calculate days
for(i in 1:length(weights$xICDate)){
  weights$Days[i] <- as.numeric(difftime(weights$xDate[i], weights$xICDate, units = "days")) 
  
}


# weights<- weights %>%
#   arrange(xDate)
# 
 b <- unique(weights$Days)

weights_msd <- weights %>%   # Get mean & standard deviation by group
  group_by(Days, Group) %>%
  summarise_at(vars(Percent),
               list(mean = mean,
                    sd = sd)) %>% 
  as.data.frame(weights_msd)

weights <- merge(weights, weights_msd, by = c("Days", "Group"))


#make separate datasets for each group of mice
Group_1<- weights[weights$Group == "1",]
Group_2<- weights[weights$Group == "2",]
Group_3<- weights[weights$Group == "3",]
Group_4<- weights[weights$Group == "4",]


#function to make weight plot without red 90
#c = color 
wts <- function(c){
  
totalwt <<- ggplot(weights, aes(x = Days, y = Percent, color = fGroup)) + 
  #geom_point(size = 2) +   
  geom_jitter( position = position_jitter(0.25), size = 1.5) +
  geom_errorbar(data = weights, aes(x = Days, ymin=mean, ymax=mean), 
                width=0.5, size = 1.3, ) +
  #geom_errorbar(data = weights, aes(x = Date, ymin=  mean-sd, ymax=mean+sd), 
                #width=0.5, size = 1.3, ) +
  geom_line(data = weights, aes(x = Days, y= mean, group = Group), size = 1) +
  scale_color_manual(labels = c("Group 1: Control", "Group 2: XXXX XXXmg/kg", "Group 3: XXXX XXXmg/kg", "Group 4: XXXX XXXmg/kg"), 
                     name = "Treatment Group",
                     values = c( "red2", "dodgerblue3", "darkseagreen3", " dark orange"))+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major.y = element_line( size=.1, color="black" ),
        axis.line = element_line(colour = "black")) +
  
  labs(y= "Body Weight (%)", x = "Time Since IC Injection (Days)") +
  geom_hline(yintercept = 90.0, size = 0.7, color = "red") +
  
  
  #limits to shrink distance
  scale_x_continuous(breaks = seq(0,max(weights$Days) + 1, 5)) + 
  scale_y_continuous(limits = c(70, max(weights$Percent) + 5),
                     breaks = seq(70,(max(weights$Percent) + 5), 5)) + 

  theme(axis.text.x=element_text(family="sans", size=9, color = "black"), 
        axis.text.y=element_text(family="sans", size=9, color = c),
        axis.ticks.y = element_line(color = c),
        axis.title.x = element_text(family="sans", face="bold", size=11, color = "black"),
        axis.title.y = element_text(family="sans", face="bold", size=11, color = "black"),
        legend.key=element_blank())

totalwt

}

#wts("black") #test function

wts("red")


#initialize red labeling 
yLabVals <- as.numeric(ggplot_build(totalwt)$layout$panel_params[[1]]$y$get_labels())
yLabs <- ifelse(yLabVals == 90, "red", "black") 

#print plot with appropriate red labeling at 90
wts(yLabs)


#function to plot each group in an individual plot
wtplot1 <- function(d, c, t){
  ggplot(d, aes(x = Days, y = Percent)) + 
    geom_jitter( position = position_jitter(0.25), size = 1.5) +
    geom_errorbar(data = d, aes(x = Days, ymin=mean, ymax=mean,), 
                  width=1, size = 1.3, color = c ) +
    geom_errorbar(data = d, aes(x = Days, ymin= mean-sd, ymax= mean + sd), 
                  width=2, size = 1.3, color = c) +
    
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(),
          panel.background = element_blank(), 
          panel.grid.major.y = element_line( size=.1, color="black" ),
          axis.line = element_line(colour = "black")) +
          
          
          ggtitle(t)+
    
          labs(y= "Body Weight (%)", x = "Time Since IC Injection (Days)") +
          geom_hline(yintercept = 90.0, size = 0.7, color = "red") +
    
    
    #limits to shrink distance
    scale_x_continuous(breaks =  seq(0,max(weights$Days) + 3, 10), limits = c(10, max(weights$Days) + 3)) + 
    scale_y_continuous(limits = c(70, 110),
                       breaks = seq(70,110,5)) + 
          
    
    theme(axis.text.x=element_text(family="sans", size=9, color = "black"), 
          axis.text.y=element_text(family="sans", size=9, color = yLabs),
          axis.title.x = element_text(family="sans", face="bold", size=11, color = "black"),
          axis.title.y = element_text(family="sans", face="bold", size=11, color = "black"),
          axis.ticks.y = element_line(color = yLabs),
          legend.key=element_blank()) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_line( size=.1, color="black" ),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 9)) 
  
}



G1<- wtplot1(Group_1, "red2",  "Group 1: Control")
G2<- wtplot1(Group_2, "dodgerblue3", "Group 2: XXXX XXXmg/kg")
G3<- wtplot1(Group_3, "darkseagreen3", "Group 3: XXXX XXXmg/kg")
G4<- wtplot1(Group_4, "darkorange", "Group 4: XXXX XXXmg/kg")



grid.arrange(arrangeGrob(G1, G2, G3, G4, ncol = 2))



#function to plot each mouse in an individual plot
mouseplot <- function(d, t){
  ggplot(d, aes(x = Days, y = Percent, color = ID)) + 
    geom_point(size = 2) +   
    geom_line(data = d, aes(x = Days, y= Percent, group = ID), size = 1)+
    labs(y= "Body Weight (%)", x = "Time Since IC Injection (Days)", color = "Mouse ID") +
    
    #limits to shrink distance
    scale_x_continuous(breaks =   seq(0,max(weights$Days) + 3, 10), limits = c(10, max(weights$Days) + 3)) + 
    scale_y_continuous(limits = c(70, 110),
                       breaks = seq(70,110,5),
                       minor_breaks = seq(70,110,1)) +
    ggtitle(t)+
    geom_hline(yintercept = 90.0, size = 0.7, color = "red") +
    
    
    theme(axis.text.x=element_text(family="sans", size=9, color = "black"), 
          axis.text.y=element_text(family="sans", size=9, color = yLabs),
          axis.title.x = element_text(family="sans", face="bold", size=11, color = "black"),
          axis.title.y = element_text(family="sans", face="bold", size=11, color = "black"),
          axis.ticks.y = element_line(color = yLabs),
          legend.key=element_blank()) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major.y = element_line( size=.1, color="black" ),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size = 9)) 
        
}

G1m<- mouseplot(Group_1,"Group 1: Control")
G2m<- mouseplot(Group_2, "Group 2: XXXX XXXmg/kg")
G3m<- mouseplot(Group_3,  "Group 3: XXXX XXXmg/kg")
G4m<- mouseplot(Group_4, "Group 4: XXXX XXXmg/kg")



grid.arrange(arrangeGrob(G1m, G2m, G3m, G4m, ncol = 2))

###########################SURVIVAL PLOTS#######################################

#Swimmer plot

#read in data file
swimdata <- read.csv(file = 'swimdata.RE.csv', header = TRUE, na.strings = "")

#convert Edate and Date to a date var and calculate days 
swimdata$ICDate <- as.Date(swimdata$ICDate, format = "%m/%d/%Y")
swimdata$Edate <- as.Date(swimdata$Edate, format = "%m/%d/%Y")

library(data.table)
swimdata$Edate <- fifelse(is.na(swimdata$Event), as.Date(Sys.Date(),format = "%m/%d/%Y"), swimdata$Edate)


#calculate days
for(i in 1:length(swimdata$ICDate)){
  swimdata$Days[i] <- as.numeric(difftime(swimdata$Edate[i], swimdata$ICDate, units = "days")) 
  
}


swimdata.AE <- swimdata[!is.na(swimdata$Event),]

swimdata$fGroup <- as.factor(swimdata$Group)
swimdata.AE$fGroup <- as.factor(swimdata.AE$Group)


swimdata$Event <- ifelse(is.na(swimdata$Event), "Continued Response", swimdata$Event)
swimcont <- subset(swimdata, Event == "Continued Response")


swimtrt <- swimmer_plot(df= swimdata,id='ID',end='Days', name_fill='fGroup',
                        id_order='Group', col="black",alpha=0.75,width=.8 ) 

swimAE <- swimtrt + swimmer_points(df_points = swimdata.AE,id='ID',time ='Days',name_shape =
                                     'Event',size=2.5 , fill = 'white', col='black') 


swimplot <- swimAE + swimmer_arrows(df_arrows = swimcont ,id='ID',arrow_start='Days',
                                    cont = 'Event',name_col='fGroup',type =
                                      "open", show.legend = FALSE, cex= 1, arrow_positions = c(0.2,1.7))

swimplot <-  swimplot +
  
  scale_fill_manual(labels = c("Group 1: Control", "Group 2: XXXX XXXmg/kg",
                                "Group 3: XXXX XXXmg/kg", "Group 4:XXXX XXXmg/kg"), 
                     name = "Treatment Group",
                     values = c("1" = "red2", "2" = "dodgerblue3", "3" = "darkseagreen3", "4" =  "darkorange")) +

  
  scale_color_manual(labels = c("Group 1: Control", "Group 2: XXXX XXXmg/kg",
                                "Group 3: XXXX XXXmg/kg", "Group 4:XXXX XXXmg/kg"), 
                     name = "Treatment Group",
                     values = c("1" = "red2", "2" = "dodgerblue3", "3" = "darkseagreen3", "4" =  "darkorange")) +
  
  scale_shape_manual(name="Mortality Event", values = c(15, 18))+ 
  
  labs(y= "Time Since IC Injection (Days)", x = "Mouse ID") +
  
  theme(axis.text.x=element_text(family="sans", size=9, color = "black"), 
        axis.text.y=element_text(family="sans", size=9, color = "black"),
        axis.title.x = element_text(family="sans", face="bold", size=11, color = "black"),
        axis.title.y = element_text(family="sans", face="bold", size=11, color = "black"),
        legend.key=element_blank()) +
  
  
  coord_flip(clip = 'on', ylim=c(35,max(swimdata$Days) + 5)) +
  scale_y_continuous(breaks= seq(35,max(swimdata$Days) + 20,10), 
                     limits = c(0, max(swimdata$Days) + 10))

swimplot


#Kaplan Meier Curve

#read in survival dataset
survival <- read.csv(file = 'KaplanMeierP3.csv', header = TRUE)

# #keep columns only where there are no na values (colsum of NAs is equal to 0)
# survival <- survival[ , colSums(is.na(survival))==0]
#change the name of "Serial Time in Days" to "Time"
colnames(survival)[2] ="Time"


require("survival")
fit <- survfit(Surv(Time, Status) ~ Group, data = survival)

# Pairwise survdiff
res <- pairwise_survdiff(Surv(Time, Status) ~ Group,
                         data = survival, p.adjust.method = "BH")
res

multicomp <- res$p.value
multicomp <- as.data.frame(multicomp)


round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

multicomp <- round_df(multicomp, digits = 4)


p1 <- ggsurvplot(
  fit,
  data = survival,
  risk.table = TRUE, # Add risk table
  pval = TRUE,       # Add p-value
  palette =
    c("red2", "dodgerblue3", "darkseagreen3", "darkorange"),# custom color palettes
  xlim = c(0,max(survival$Time) + 1),  
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 10,     # break X axis in time intervals by 10.
  ggtheme = theme_light() ,     # Change ggplot2 theme
  legend = "top",
  legend.title = "Treatment Group",
  
  risk.table.height = 0.25, 
  
  
  
  risk.table.y.text.col = T, #colour risk table text annotations.
  risk.table.y.text = FALSE,
  legend.labs =
    c("Group 1: Control", "Group 2: XXXX XXXmg/kg", 
      "Group 3: XXXX XXXmg/kg", "Group 4: XXXX XXXmg/kg"),# Change legend labels
  # show bars instead of names in text annotations
  
) 





