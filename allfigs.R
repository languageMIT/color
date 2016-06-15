# Produce all figures and useful numbers.
# Based on GetEntropy.R by Julian,
# universal_snake.R and color_heatmap.R and focal_locations.R by Richard,
# and color_analysis.Rmd by Kyle

library(MASS)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(hexbin)

# todo analysis using only the modal term for each chip
# spearman correlatios of average surprisal across en es ts
# education correlations


# We've determined the Julian-style data contains duplicated subjects so we should
# use the Kyle-style data. (Meeting on 2016-06-02.)
LABELLING_DATA = "Kyle"  # "Julian" or "Kyle"

SPECIAL_FOCAL_COLORS = c("blue", "green", "red", "yellow")


TSIMANE_COLORS = c("black", "white", "red", "blue", "green", "yellow",
                   "grey", "purple / violet", "orange", "brown", "yellowish",
		   "yellow", "orange", "greenish", "brownish")

ENGLISH_COLORS = c("black", "white", "red", "blue", "green", "yellow",
                   "grey", "purple / violet", "orange", "brown", "pink", "celeste")

TSIMANE_SUBJECTS_UNDER_18 = c(8, 15, 21, 25, 58, 59, 1509)

get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


# Load data -------------------------------

RowColumnToCode = function(Row, Column) {
  NewRow = LETTERS[-Row]
  paste(NewRow, Column, sep="")
}

CodeToRowColumn = function(Code) {
  letter = str_sub(Code, 1, 1)
  column = str_sub(Code, 2, str_length(Code))
  row = -match(letter, LETTERS)
  str_c(row, "_", column)
}

assert = stopifnot

if (LABELLING_DATA == "Julian") {

# Tsimane data open
d1<-read.csv("s19_66_half_grid.csv") %>%
  dplyr::select(-order) %>%
  mutate(Experiment="Half Grid") %>% tbl_df
# Subject 0 is Salomon.
d2<-read.csv("s0-18_big_grid.csv") %>%
  mutate(Experiment="Full Grid") %>% tbl_df %>% filter(subject!="0")
#d3<-read.csv("7_subjs_double_half_grid.csv") %>%
#  filter(first_second==1) %>%
#  dplyr::select(-order,-first_second) %>%
#  mutate(Experiment="Double Grid") %>% tbl_df
#d3<-dplyr::rename(d3,grid_location=grid_.location)
TsimaneData<-rbind(d1,d2)
rm(d1,d2)
KeepChips<-as.data.frame(table(TsimaneData$grid_location)==max(table(TsimaneData$grid_location)))
names(KeepChips)<-"Keep"
KeepChips$Chip<-rownames(KeepChips)
row.names(KeepChips) <- NULL
KeepChips<-KeepChips[which(KeepChips$Keep),]
TsimaneData$grid_location<-as.character(TsimaneData$grid_location)
TsimaneData<-filter(TsimaneData,grid_location %in% KeepChips$Chip)
assert(dim(table(TsimaneData$grid_location)) == 80)
rm(KeepChips)
TDO<-filter(TsimaneData,!is.na(color),color!="")

# Load Tsimane data fixed
TDF<-read.csv("Tsimane_WCS.csv") %>% tbl_df

# Load English data open
EDO<-read.csv("USColorChipLabels nov6 2014.csv") %>% tbl_df

# Load English data fixed
EDF<-read.csv("English_WCS.csv") %>% tbl_df

# Load Spanish data open
SDO <- read.csv("Spanish_open.csv") %>% tbl_df

# Load Spanish data closed
SDF <- read.csv("Spanish_WCS.csv") %>% tbl_df


# Have all dataset match in columns and names
EDF<-dplyr::select(EDF,-order) %>% dplyr::rename(grid_location=location)
EDO<-dplyr::select(EDO,-order) %>% dplyr::rename(grid_location=location)
SDF<-dplyr::select(SDF,-Experiment)
SDO<-dplyr::select(SDO,-Experiment)
TDF<-dplyr::select(TDF,-Experiment)
TDO<-dplyr::select(TDO,-Experiment)
length(unique(EDF$grid_location))
length(unique(EDO$grid_location))
length(unique(SDF$grid_location))
length(unique(SDO$grid_location))
length(unique(TDF$grid_location))
length(unique(TDO$grid_location))

CommonChips<-EDF$grid_location %>% intersect(EDO$grid_location) %>%
  intersect(TDO$grid_location) %>% intersect(TDO$grid_location) %>%
  intersect(SDO$grid_location) %>% intersect(SDO$grid_location)

TDO$Experiment="Tsimane_Open"
TDF$Experiment="Tsimane_Fixed"
EDO$Experiment="English_Open"
EDF$Experiment="English_Fixed"
SDO$Experiment="Spanish_Open"
SDF$Experiment="Spanish_Fixed"

 load("sparse_chips_to_use.rda")

ColorData<-tbl_df(rbind(TDO,TDF,EDO,EDF,SDO,SDF))
ColorData<-filter(ColorData,grid_location %in% sparse_chips_to_use)

} else if (LABELLING_DATA == "Kyle") {

 load("sparse_chips_to_use.rda")

 ColorData = read.csv("color labeling tsimane english spanish oct 2015 v4.csv") %>%
   rename(grid_location=location) %>%
   select(subject, color, Language, Task, grid_location) %>%
   unite(Experiment, Language, Task) %>%
   filter(grid_location %in% sparse_chips_to_use)

}


m = read.csv("Munsell_WCS_codes.csv") %>%
  rename(grid_location=Our.code,
         MunsellCode=OtHer.code,
	 WCSCode=WCS)


# Read demographics --------------------------------------------

#demo2014 = read.csv("demographics2014.csv", stringsAsFactors=F) %>% mutate(Year=2014)
#demo2015 = read.csv("demographics2015.csv", stringsAsFactors=F) %>% mutate(Year=2015)
#demo2015 %>%
#  rename(Subject=Subject..) %>%
#  mutate(year=2015) %>%
#  select(Subject, Age, Sex, Spanish, Year) %>%
#  rbind(demo2014 %>% select(Subject, Age, Sex, Spanish, Year)) -> ts_demo
#
#ColorData %>%
##  separate(Experiment, into=c("Language", "Task"), sep="_") %>%
#  filter(Language == "Tsimane'") %>%
#  rename(Subject=subject) %>%
#  mutate(Subject=as.character(Subject)) %>%
#  inner_join(ts_demo) -> ColorDataWithDemographics
#  

# Compute conditional entropy ----------------------------------

P_ChipGivenWord <- plyr::ddply(ColorData,c("Experiment","color"),function(x){return(as.data.frame(prop.table(table(x$grid_location))))})
P_WordGivenChip <- plyr::ddply(ColorData,c("Experiment","grid_location"),function(x){return(as.data.frame(prop.table(table(x$color))))})
names(P_ChipGivenWord)<-c("Experiment","color","grid_location","Frequency")
names(P_WordGivenChip)<-c("Experiment","grid_location","color","Frequency")

# Get each chip's score
GetChipScore<-function(x){
  # Sum over all words.
  Val=0
  for (CurrCol in unique(x$color)){
    ChipWord=filter(P_ChipGivenWord,color==CurrCol,grid_location==x$grid_location[1],Experiment==x$Experiment[1])$Frequency
    WordChip=filter(P_WordGivenChip,grid_location==x$grid_location[1],color==CurrCol,Experiment==x$Experiment[1])$Frequency
    Val=Val+WordChip*log2(1.0/ChipWord)
  }
  return(Val)
}

ChipEntropies <- plyr::ddply(ColorData,c("Experiment","grid_location"),GetChipScore)
names(ChipEntropies)<-c("Experiment","grid_location","Entropy")
ChipNo<-length(unique(ChipEntropies$grid_location))
# This line assumes that each subdataset has the same number of chips
ChipEntropies <- ChipEntropies %>% mutate(Probability=1/ChipNo, ExpEnt=Entropy*Probability)

LanguageEntropies <- ChipEntropies %>%
  group_by(Experiment) %>%
  summarise(Uncertainty=sum(ExpEnt)) 

write.csv(ChipEntropies, "output/Latest_ChipEntropies.csv")

ChipEntropies <- ChipEntropies %>%
  separate(Experiment, into=c("Language", "Type"), sep="_") %>%
  mutate(RowColumn=CodeToRowColumn(grid_location)) %>%
  separate(RowColumn, into=c("Row", "Column"), sep="_") %>%
  mutate(Row=as.numeric(Row), Column=as.numeric(Column))

# Conditional entropy per subject----------------------


# Spearman correlations -------------------------------

en_es = cor.test(filter(ChipEntropies, Language == "English")$Entropy,
                 filter(ChipEntropies, Language == "Spanish")$Entropy,
          	 method='spearman')[4]$estimate

en_ts = cor.test(filter(ChipEntropies, Language == "English")$Entropy,
                 filter(ChipEntropies, Language == "Tsimane")$Entropy,
          	 method='spearman')[4]$estimate

es_ts = cor.test(filter(ChipEntropies, Language == "Spanish")$Entropy,
                 filter(ChipEntropies, Language == "Tsimane")$Entropy,
          	 method='spearman')[4]$estimate


# Informativity per language --------------------------

ChipEntropies %>%
  group_by(Language, Type) %>%
  summarise(LangEntropy=mean(Entropy)) %>%
  write.csv("output/Latest_EntropyByLanguage.csv")


# Modal informativity --------------------------

# Suppose every chip is just named according to its modal label.
ModalColorData = ColorData %>%
  group_by(Experiment, grid_location) %>%
  summarise(color=get_mode(color)) %>%
  ungroup()


Modal_P_ChipGivenWord <- plyr::ddply(ModalColorData,c("Experiment","color"),
     function(x){return(as.data.frame(prop.table(table(x$grid_location))))})
Modal_P_WordGivenChip <- plyr::ddply(ModalColorData,c("Experiment","grid_location"),
     function(x){return(as.data.frame(prop.table(table(x$color))))})
names(Modal_P_ChipGivenWord)<-c("Experiment","color","grid_location","Frequency")
names(Modal_P_WordGivenChip)<-c("Experiment","grid_location","color","Frequency")

# Get each chip's score
GetModalChipScore<-function(x){
  # Sum over all words.
  Val=0
  for (CurrCol in unique(x$color)){
    ChipWord=filter(Modal_P_ChipGivenWord,
                   color==CurrCol,
		   grid_location==x$grid_location[1],
		   Experiment==x$Experiment[1])$Frequency
    WordChip=filter(Modal_P_WordGivenChip,
                    grid_location==x$grid_location[1],
		    color==CurrCol,
		    Experiment==x$Experiment[1])$Frequency
    Val=Val+WordChip*log2(1.0/ChipWord)
  }
  return(Val)
}

ModalChipEntropies <- plyr::ddply(ModalColorData,c("Experiment","grid_location"),GetModalChipScore)
names(ModalChipEntropies)<-c("Experiment","grid_location","Entropy")
ModalChipNo<-length(unique(ModalChipEntropies$grid_location))
# This line assumes that each subdataset has the same number of chips
ModalChipEntropies <- ModalChipEntropies %>%
  mutate(Probability=1/ModalChipNo, ExpEnt=Entropy*Probability)

write.csv(ModalChipEntropies, "output/Latest_ModalChipEntropies.csv")

ModalChipEntropies <- ModalChipEntropies %>%
  separate(Experiment, into=c("Language", "Type"), sep="_") %>%
  mutate(RowColumn=CodeToRowColumn(grid_location)) %>%
  separate(RowColumn, into=c("Row", "Column"), sep="_") %>%
  mutate(Row=as.numeric(Row), Column=as.numeric(Column))

ModalChipEntropies %>%
  group_by(Language, Type) %>%
  summarise(LangEntropy=mean(Entropy)) %>%
  write.csv("output/Latest_ModalEntropyByLanguage.csv")


# By-subject average surprisal ------------------------




# Load WCS --------------------------------------------

munsell<- read.csv("chip.csv",header=T) %>% tbl_df %>% mutate(NumRow=as.numeric(GridRow))
answers<- read.csv("term.csv",header=T) %>% tbl_df
languages<- read.csv("lang.csv",header=T) %>% tbl_df %>% dplyr::select(Lang_Id,Lang_Name)
answers<-dplyr::rename(answers,Lang_Id=Language)

# Add chip data.
answers<-merge(munsell,answers,by.y="Chip",by.x="ChipNumber") %>% tbl_df
answers<-tbl_df(merge(languages,answers,by="Lang_Id"))

# Keep only what you need
answers<-dplyr::select(answers,Lang_Name,FielConcat,Term)
WCS<-dplyr::rename(answers,Language=Lang_Name,Chip=FielConcat,Color=Term)

WCS$Color<-as.character(WCS$Color)
# some colors are called NA and R thinks they're missing values
WCS[which(is.na(WCS$Color)),]$Color="NA"

# clean up
rm(answers,languages,munsell)

#Standardize with main dataset
names(WCS)<-c("Experiment","grid_location","color")
WCS$grid_location<-as.character(WCS$grid_location)

# Reduce WCS chips to those in our other analyses -------------------------------------
# Load code conversion chart
ConversionChart<-read.csv("Munsell_WCS_codes.csv") %>% tbl_df
names(ConversionChart) = c("BoliviaCode", "MainCode", "WCSCode")

ConversionChart$BoliviaCode<-as.character(ConversionChart$BoliviaCode)
ConversionChart$WCSCode<-as.character(ConversionChart$WCSCode)
ConversionChart <- filter(ConversionChart,BoliviaCode %in% unique(ColorData$grid_location)) %>%
  dplyr::rename(grid_location=WCSCode)

WCS <- filter(WCS,grid_location %in% unique(ConversionChart$grid_location))

assert(length(unique(WCS$grid_location)) == 80)

WCS <- inner_join(WCS,
                  ConversionChart,
		  by=c("grid_location")) %>%
         dplyr::select(-MainCode,-grid_location) %>%
	 dplyr::rename(grid_location=BoliviaCode)
	 
WCS %>% dplyr::group_by(Experiment) %>% dplyr::summarise(Chips=length(unique(grid_location))) %>%
  dplyr::select(Chips) %>% table # good!
  
rm(ConversionChart)

# Compute WCS mode conditional entropy ----------------------

WCS_Modes<-tbl_df(plyr::ddply(WCS,c("Experiment","grid_location"),function(x){return(data.frame(table(x$color)))}))
names(WCS_Modes)<-c("Experiment","grid_location","color","frequency")

# Get modes
WCS_Modes <- WCS_Modes %>% dplyr::group_by(Experiment,grid_location) %>% dplyr::top_n(1,frequency)
# Check total number of color words
#res<-WCS_Modes %>% dplyr::group_by(Experiment) %>% dplyr::summarise(ColorWords=length(unique(color)))
# Delete ties
WCS_Modes <- WCS_Modes %>% dplyr::group_by(Experiment,grid_location) %>% dplyr::summarise(color=color[1])

WCS_Mode_ChipGivenWord <- plyr::ddply(WCS_Modes,c("Experiment","color"),function(x){return(as.data.frame(prop.table(table(x$grid_location))))})
names(WCS_Mode_ChipGivenWord)<-c("Experiment","color","grid_location","Probability")

WCS_Modes<-full_join(WCS_Modes,WCS_Mode_ChipGivenWord)
WCS_Modes <- WCS_Modes %>% mutate(Prior=1/length(unique(WCS_Modes$grid_location)),ExpEnt=Prior*log2(1.0/Probability))

WordNumbers <- WCS_Modes %>% dplyr::group_by(Experiment) %>% dplyr::summarise(WordNumber = length(unique(color)))

WCS_Modes <- WCS_Modes %>% dplyr::group_by(Experiment) %>% dplyr::summarise(ConEnt=sum(ExpEnt))

WCS_Modes <- full_join(WordNumbers, WCS_Modes)

#WCS_Modes %>% filter(WordNumber %in% c(7,8)) %>%
#  ggplot(aes(x=WordNumber,y=ConEnt,label=Experiment))+geom_text()

# Compute WCS conditional entropy --------------------------

WCS_ChipGivenWord <- plyr::ddply(WCS,c("Experiment","color"),function(x){return(as.data.frame(prop.table(table(x$grid_location))))})
WCS_WordGivenChip <- plyr::ddply(WCS,c("Experiment","grid_location"),function(x){return(as.data.frame(prop.table(table(x$color))))})
names(WCS_ChipGivenWord)<-c("Experiment","color","grid_location","Frequency")
names(WCS_WordGivenChip)<-c("Experiment","grid_location","color","Frequency")


# Get each chip's score
WCS_GetChipScore<-function(x){
  # Sum over all words.
  Val=0
  for (CurrCol in unique(x$color)){
    ChipWord=filter(WCS_ChipGivenWord,
                    color==CurrCol,
		    grid_location==x$grid_location[1],
		    Experiment==x$Experiment[1])$Frequency
    WordChip=filter(WCS_WordGivenChip,
                    grid_location==x$grid_location[1],
		    color==CurrCol,
		    Experiment==x$Experiment[1])$Frequency
    #print(x$Experiment[1])
    #print(x$grid_location[1])
    if (WordChip != 0){
      Val=Val+WordChip*log2(1.0/ChipWord)
    }
  }
  return(Val)
}

WCS = filter(WCS, !is.na(Experiment))
WCS_ChipEntropies <- plyr::ddply(WCS,c("Experiment","grid_location"),WCS_GetChipScore)
names(WCS_ChipEntropies)<-c("Experiment","grid_location","Entropy")

ChipNo<-length(unique(WCS_ChipEntropies$grid_location))
# This line assumes that each subdataset has the same number of chips: CONFIRMED ITS CORRECT
WCS_ChipEntropies <- WCS_ChipEntropies %>% mutate(Probability=1/ChipNo, ExpEnt=Entropy*Probability)

write.csv(WCS_ChipEntropies, "output/Latest_WCS_ChipEntropies.csv")

ChipEntropies %>%
  filter(Type == "Open") %>%
  select(Language, grid_location, Entropy) %>%
  rbind(WCS_ChipEntropies %>%
        rename(Language=Experiment) %>%
	select(Language, grid_location, Entropy)) %>%
  spread(Language, Entropy) %>%
  inner_join(m) %>%
  rename(OurCode=grid_location) %>%
  write.csv("output/Latest_All_ChipEntropies.csv")


# Load focal location data ----------------------------

d5 = read.csv("Tsimane2015_2_focal_colors_object_colors.csv")
d4 = read.csv("Tsimane2014_focal_colors_object_colors.csv")

d4 %>%
  filter(Toss == "N") %>%
  gather(term, Code, 13:30) %>%
  select(Subject, term, Code) %>%
  filter(Code != "") %>%
  separate(term, into=c("tmp1", "tmp2", "term"), sep="_") %>%
  select(-tmp1, -tmp2) -> d4_choices

d5 %>%
  gather(term, Code, 14:25) %>%
  select(Subject.number, term, Code) %>%
  filter(Code != "") %>%
  separate(term, into=c("tmp1", "tmp2", "term", "mod"), sep="_") %>%
  filter(is.na(mod)) %>% # disregard _2 terms
  mutate(Subject = Subject.number) %>%
  select(-tmp1, -tmp2, -mod, -Subject.number) -> d5_choices

d = rbind(d4_choices, d5_choices) %>% filter(!is.na(term))

focal_ts = d %>% mutate(Language = "Tsimane")
focal_ts = filter(focal_ts, !str_detect(term, "_")) # drop _2 colors

get_en_es_focals = function(filename, lang) {
  d = read.csv(filename) %>% select(-X) %>% mutate(Subject=subject) %>% select(-subject)
  d = filter(d, str_detect(thing, "focal"))
  d %>%
    mutate(Language=lang) %>%
    rename(term=thing, Code=color) %>%
    filter(Code != "") %>%
    mutate(term=str_replace(term, "focal ", ""))
 }	

focal_en = get_en_es_focals("english_objects.csv", "English")
focal_es = get_en_es_focals("spanish_objects.csv", "Spanish")

focal_locations = rbind(focal_ts, focal_en, focal_es)

focal_locations %>%
  group_by(Language, term, Code) %>%
  summarise(freq=length(Code)) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(Z=sum(freq)) %>%
  mutate(p=freq/Z) %>%
  select(-freq, -Z) %>%
  write.csv("output/focal_P_ChipGivenWord.csv")

write.csv(focal_locations, "output/focal_locations.csv")

focal_locations = focal_locations %>%
  mutate(RowColumn=CodeToRowColumn(Code)) %>%
  separate(RowColumn, into=c("Row", "Column"), sep="_") %>%
  mutate(Row=as.numeric(Row), Column=as.numeric(Column))


# For each focal term, how many subjects know that term?
# first need to convert the color indices in ColorData into terms
ColorData %>%
  separate(Experiment, into=c("Language", "Task"), sep="_") %>%
  filter(Task == "Open") %>%
  select(subject, color, Language) %>%
  mutate(color=as.character(color)) %>%
  group_by(Language) %>%
  mutate(color=ifelse(Language == "Tsimane",
                      TSIMANE_COLORS[as.numeric(color)],
		      ENGLISH_COLORS[as.numeric(color)])) %>%
  ungroup() %>%		      
  filter(!is.na(color)) %>%
  filter(color %in% focal_locations$term) %>%
  group_by(Language, color) %>%
  filter(!(Language == "English" & color == "celeste")) %>%
  mutate(n=length(unique(subject))) %>%
  ungroup() %>%
  group_by(Language) %>%
  mutate(Z=length(unique(subject))) %>%
  ungroup() %>%
  select(-subject) %>%
  unique() %>%
  mutate(prop_subjects=n/Z) %>%
  write.csv("output/how_many_subjects_use_terms.csv")


# Numbers of subjects ---------------------------------

rts = read.csv("english_rts.csv") %>%
  select(-Native_speaker, -X) %>%
  mutate(Language="English") %>%
  rbind(read.csv("tsimane_rts.csv") %>% mutate(Language="Tsimane"))

ColorData %>%
  separate(Experiment, into=c("Language", "Task"), sep="_") %>%
  group_by(Language, Task) %>%
  summarise(NumSubjects=length(unique(subject))) %>%
  ungroup() ->
  labelling_num_subjects


focal_locations %>%
  group_by(Language) %>%
  summarise(NumSubjects=length(unique(Subject))) %>%
  ungroup() %>%
  mutate(Task="FocalColors") ->
  objects_num_subjects

rts %>%
  group_by(Language) %>%
  summarise(NumSubjects=length(unique(subject))) %>%
  ungroup() %>%
  mutate(Task="RT") ->
  rts_num_subjects

num_subjects = rbind(labelling_num_subjects, objects_num_subjects, rts_num_subjects) %>%
  arrange(Task, Language) 

write.csv(num_subjects, "output/num_subjects.csv")

# Chips used in each experiment -----------------------

m %>%
  mutate(InLabelling=(grid_location %in% ColorData$grid_location) | grid_location == "C3",
         InObjects=T,
	 InRT=grid_location %in% rts$color_object.label) %>%
  rename(OurCode=grid_location) %>%	 
  write.csv("output/chips_in_experiments.csv")	 
	 
  

# Plot Heatmap ----------------------------------------

LETTERS = c("A", "B", "C", "D", "E", "F", "G", "H")

load("cc.rda") # Load the chip code -> hex code conversion table



ts_color_indices = focal_locations[focal_locations$Language == "Tsimane",]$term %>% as.numeric()
ts_colors = TSIMANE_COLORS[ts_color_indices]
focal_locations[focal_locations$Language == "Tsimane",]$term = ts_colors

special_focal_locations = filter(focal_locations, term %in% SPECIAL_FOCAL_COLORS)
special_focal_density = group_by(special_focal_locations, Language, Row, Column, term) %>% summarise(value=n())
write.csv(special_focal_density, "output/special_focal_density.csv") # to be used by contours.py

# E10 is missing, so interpolate it from E8, E12, D9, D11, F9, F11
interp_sources = c("E8", "E12", "D9", "D11", "F9", "F11")
interp_rows = ChipEntropies %>%
  filter(grid_location %in% interp_sources) %>%
  group_by(Language, Type) %>%
  summarise(Entropy=mean(Entropy)) %>%
  ungroup() %>%
  mutate(grid_location="E10", Probability=0.0125)
  
  
ChipEntropiesInterp = ChipEntropies %>%
  select(Language, Type, Entropy, grid_location, Probability) %>%
  rbind(interp_rows) %>%
  mutate(RowColumn=CodeToRowColumn(grid_location)) %>%
  separate(RowColumn, into=c("Row", "Column"), sep="_") %>%
  mutate(Row=as.numeric(Row), Column=as.numeric(Column))

# plot heatmap of entropy 
(
ggplot(ChipEntropiesInterp, aes(x=Column, y=Row, fill=Entropy))
    + facet_grid(Type~Language)
    + geom_hex(stat="identity")
    + scale_y_continuous(breaks=c(-1,-2,-3,-4,-5,-6,-7,-8),
                         labels=c("A", "B", "C", "D", "E", "F", "G", "H"))
    + scale_fill_continuous(name="Average surprisal", low="black", high="white")
    + theme_bw()
    + theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
)

ggsave("output/heatmaps.pdf", height=5, width=13)

# Contours in ggplot are uncontrollable, so do them in matplotlib instead
##
##kde2d.weighted <- function (x, y, w, h, n = 25, lims = c(range(x), range(y))) {
##  nx <- length(x)
##  if (length(y) != nx)
##    stop("data vectors must be the same length")
##  if (length(w) != nx & length(w) != 1)
##    stop("weight vectors must be 1 or length of data")
##  gx <- seq(lims[1], lims[2], length = n) # gridpoints x
##  gy <- seq(lims[3], lims[4], length = n) # gridpoints y
##  if (missing(h))
##    h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
##  if (missing(w))
##    w <- numeric(nx)+1;
##  h <- h/4
##  ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
##  ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
##  z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*%
##       t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) # z is the density
##  return(list(x = gx, y = gy, z = z))
##}
##
##getLevel <- function(x, y, z, prob) {
##    kk <- kde2d.weighted(x,y,z)
##    dx <- diff(kk$x[1:2])
##    dy <- diff(kk$y[1:2])
##    sz <- sort(kk$z)
##    c1 <- cumsum(sz) * dx * dy
##    approx(c1, sz, xout = 1 - prob)$y
##}
##
##L20 <- with(special_focal_density, getLevel(Row, Column, value, .2))
##L40 <- with(special_focal_density, getLevel(Row, Column, value, .4))
##L60 <- with(special_focal_density, getLevel(Row, Column, value, .6))
##L80 <- with(special_focal_density, getLevel(Row, Column, value, .8))
##
### plot contours of focal locations
##(
##ggplot(special_focal_density, aes(x=Column, y=Row, z=value, color=term))
##    + facet_grid(~Language)
##    + scale_y_continuous(breaks=c(-1,-2,-3,-4,-5,-6,-7,-8),
##                         labels=c("A", "B", "C", "D", "E", "F", "G", "H"))
##    + theme_bw()
##    + theme(panel.grid.major = element_blank(),
##            panel.grid.minor = element_blank())
##    + geom_contour(stat="density2d", bins=4, breaks=c(L20, L40, L60, L80))
##    + scale_color_manual(name="Focal color", values=SPECIAL_FOCAL_COLORS)
##)
##
##ggsave("output/focal_contours.pdf", height=5, width=13)

# Save data of for each chip how many times it is chosen as focal for any color
focal_locations %>%
  mutate(Code=RowColumnToCode(Row, Column)) %>%
  select(-Row, -Column) %>%
  group_by(Language, Code) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  write.csv("output/overall_focal_density.csv")


### heatmap + contours
##(
##ggplot(ChipEntropiesInterp, aes(x=Column, y=Row, fill=Entropy))
##    + facet_grid(Type~Language)
##    + geom_hex(stat="identity")
##    + scale_y_continuous(breaks=c(-1,-2,-3,-4,-5,-6,-7,-8),
##                         labels=c("A", "B", "C", "D", "E", "F", "G", "H"))
##    + scale_fill_continuous(name="Average surprisal", low="black", high="white")
##    + theme_bw()
##    + theme(panel.grid.major = element_blank(),
##            panel.grid.minor = element_blank())
##    + geom_density2d(data=special_focal_density,
##                     aes(x=Column, y=Row, z=value, fill=1, color=term),
##		     h=3, bins=4)
##    + scale_color_manual(name="Focal color", values=SPECIAL_FOCAL_COLORS)	    
##)
##
##ggsave("output/heatmaps_with_focal_contours.pdf", height=5, width=13)


# Snake plot ----------------------------------------

d = ChipEntropies %>%
  filter(Type == "Fixed") %>%
  select(Language, grid_location, Entropy) 
d = WCS_ChipEntropies %>%
      rename(Language=Experiment) %>%
      select(Language, grid_location, Entropy) %>%
      rbind(d)
d = rename(d, Code=grid_location)

d_open = ChipEntropies %>% rename(Code=grid_location) %>% filter(Type == "Open") %>% select(-Type)

# For some reason group_by does the wrong thing in combination with rank(), so do it the bad old way...
d$score_rank = NA
for (lang in unique(d$Language)) {
  d[d$Language == lang,]$score_rank = rank(d[d$Language == lang,]$Entropy, ties.method=c("random"))
}

d_open$score_rank = NA
for (lang in unique(d$Language)) {
  d_open[d_open$Language == lang,]$score_rank = rank(d_open[d_open$Language == lang,]$Entropy, ties.method=c("random"))
}

d %>%
  ggplot(aes(x=score_rank,
             y=Entropy,
	     xmin=score_rank-.5,
	     xmax=score_rank+.5,
	     ymin=Entropy-.1,
	     ymax=Entropy+.1,
	     color=Code,
	     fill=Code)) +
    theme_bw(18) +
    scale_color_manual(values=cc, guide=F) +
    scale_fill_manual(values=cc, guide=F) +
    facet_wrap(~Language, ncol=4) +
    geom_rect() +
    xlab("Chips (rank ordered)") + 
    ylab("Average surprisal")

ggsave("output/WCS_snake.pdf", width=7, height=50, limitsize=F)

d %>%
  filter(Language %in% c("English", "Spanish", "Tsimane")) %>%
  ggplot(aes(x=score_rank,
             y=Entropy,
	     xmin=score_rank-.5,
	     xmax=score_rank+.5,
	     ymin=Entropy-.1,
	     ymax=Entropy+.1,
	     color=Code,
	     fill=Code)) +
    theme_bw(18) +
    scale_color_manual(values=cc, guide=F) +
    scale_fill_manual(values=cc, guide=F) +
    facet_wrap(~Language, ncol=4) +
    geom_rect() +
    xlab("Chips (rank ordered)") + 
    ylab("Average surprisal") +
    ylim(0, 6)

ggsave("output/snake.pdf", height=3, width=10)

d_open %>%
  filter(Language %in% c("English", "Spanish", "Tsimane")) %>%
  ggplot(aes(x=score_rank,
             y=Entropy,
	     xmin=score_rank-.5,
	     xmax=score_rank+.5,
	     ymin=Entropy-.1,
	     ymax=Entropy+.1,
	     color=Code,
	     fill=Code)) +
    theme_bw(18) +
    scale_color_manual(values=cc, guide=F) +
    scale_fill_manual(values=cc, guide=F) +
    facet_wrap(~Language, ncol=4) +
    geom_rect() +
    xlab("Chips (rank ordered)") + 
    ylab("Average surprisal") +
    ylim(0, 6)

ggsave("output/snake_open.pdf", height=3, width=10)


# Tapestry plot ----------------------------------------

d_to_plot = d %>%
  mutate(Language=factor(Language)) %>%
  group_by(Language) %>%
  mutate(lang_score=sum(Entropy), min_score=min(Entropy), max_score=max(Entropy)) %>%
  ungroup() %>%
  mutate(Language=reorder(factor(Language), -lang_score))

d_to_plot %>%
  ggplot(aes(x=score_rank, y=Language, color=Code, fill=Code)) +
  scale_color_manual(values=cc, guide=F) +
  scale_fill_manual(values=cc, guide=F) +
  geom_tile() +
  xlab("Chips by increasing average surprisal =>") +
  ylab("<= Languages by increasing total surprisal") +
  theme_classic() +
  theme(axis.ticks=element_blank(),
        axis.text.x=element_blank(),
	axis.text.y = element_text(size=6,
	   face=ifelse(levels(d_to_plot$Language) %in% c("Tsimane", "English", "Spanish"),
	               "bold", "plain"),
           color=ifelse(levels(d_to_plot$Language) %in% c("Tsimane", "English", "Spanish"),
	               "blue", "black")))

ggsave("output/WCS_tapestry.pdf", width=7, height=8)


# Tornado plots ----------------------------------------

HINC = 0.02
VINC = 0.4

d_to_plot %>%
  ggplot(aes(x=Entropy,
             y=Language,
	     color=Code,
	     fill=Code,
	     xmin=Entropy-HINC,
	     xmax=Entropy+HINC,
	     ymin=as.numeric(Language)-VINC,
	     ymax=as.numeric(Language)+VINC)) +
    scale_color_manual(values=cc, guide=F) +
    scale_fill_manual(values=cc, guide=F) +
    geom_rect() +
    theme_bw() +
    theme(axis.text.y = element_text(size=6,
	    face=ifelse(levels(d_to_plot$Language) %in% c("Tsimane", "English", "Spanish"),
	                "bold", "plain"),
            color=ifelse(levels(d_to_plot$Language) %in% c("Tsimane", "English", "Spanish"),
	                 "blue", "black")))

ggsave("output/WCS_tornado.pdf", height=8, width=7)

only_scores = d_to_plot %>%
  select(Language, lang_score, min_score, max_score) %>%
  unique()

only_scores %>%
  ggplot(aes(x=lang_score/80,
             y=Language,
	     xmin=min_score,
	     xmax=max_score)) +
    geom_errorbarh() +
    geom_point(color="red") +
    xlab("Average surprisal per chip") +
    theme_bw() +
    theme(axis.text.y = element_text(size=6,
	    face=ifelse(levels(d_to_plot$Language) %in% c("Tsimane", "English", "Spanish"),
	                "bold", "plain"),
            color=ifelse(levels(d_to_plot$Language) %in% c("Tsimane", "English", "Spanish"),
	                 "blue", "black")))

ggsave("output/WCS_simple_tornado.pdf", height=8.5, width=7)

d %>%
  group_by(score_rank) %>%
  summarise(Code=get_mode(Code)) %>%
  mutate(y=1) %>%
  ggplot(aes(x=score_rank, y=y, color=Code, fill=Code)) +
    scale_color_manual(values=cc, guide=F) +
    scale_fill_manual(values=cc, guide=F) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
	  axis.text.x = element_blank(),
	  axis.ticks.x = element_blank()) +
    ylab("") +
    xlab("Colors by increasing average surprisal =>")

ggsave("output/WCS_strip.pdf", height=1, width=10)


# Diamond plots ----------------------------------------

blanky = theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
	       panel.background = element_blank(),
	       axis.line = element_line(colour = "black"),
	       axis.text.x = element_text(size=18),
	       axis.text.y=element_text(size=18),
	       axis.title.y=element_text(size=18))

ccd = data.frame(cc) %>% mutate(grid_location=rownames(.)) %>% rename(Hex=cc)

emp.foc = read.csv("empirical_focal3.csv", stringsAsFactors=F)
names(emp.foc) = c("Language", "term", "grid_location", "Hex")
emp.foc$term = as.character(emp.foc$term)

# for each chip, figure out which term was most used for that chip, and show the focal or modal color of that term

    x = ColorData %>%
      separate(Experiment, into=c("Language", "Task"), sep="_") %>%
      mutate(RowColumn=CodeToRowColumn(grid_location)) %>%
      separate(RowColumn, into=c("Row", "Column"), sep="_") %>%
      mutate(Row=as.numeric(Row), Column=as.numeric(Column)) 

    x.sum <- group_by(x, grid_location, Language, Task, Row, Column) %>%
      mutate(location.n = n()) %>%
      group_by(grid_location, color, location.n, Language, Task, Row, Column) %>%
      summarise(n=n()) %>%
      ungroup() %>%
      group_by(grid_location, location.n, Language, Task, Row, Column) %>%
      mutate(n.max=max(n)) %>%
      filter(n == n.max)
      
    x.sum$modal.pct = x.sum$n / x.sum$location.n

    # If we have a focal hex for a term in emp.foc, use that, otherwise use modal chip hex
    x.sum$mode.hex = NA
    for (i in 1:nrow(x.sum)) {
      lang = x.sum[i,]$Language
      term = x.sum[i,]$color
      if (term %in% emp.foc[emp.foc$Language == lang,]$term) {
        x.sum[i,]$mode.hex = as.character(emp.foc[emp.foc$Language == lang & emp.foc$term == term,]$Hex)
      } else {
        loc = get_mode(filter(ColorData,
	                      Experiment == paste(lang, "Open", sep="_"),
			      color == term)$grid_location)
        x.sum[i,]$mode.hex = as.character(ccd[ccd$grid_location == loc,]$Hex)
      }
    }

diamond = function(x, ym=15) {
  x = rename(x, location=grid_location)
  x$color = paste(substr(x$Language, 1, 3), ": ", x$color, sep="")
  x$location = as.character(x$location)
  x$row = substr(x$location, 2, nchar(x$location))
  x$col = substr(x$location, 1, 1)
  nums = unique(x$n)
  names(nums) = nums
  x$row = as.numeric(as.character(x$row))
  key = unique(x[, c("mode.hex", "color")])
  rownames(key) = as.character(key$color)
  cols = as.character(key$mode.hex)
  names(cols) = as.character(key$color)
  x$col = as.factor(x$col)
  x$col = factor(x$col, levels=rev(levels(x$col)))
  x$Language = factor(x$Language, levels=c("English", "Spanish", "Tsimane"))
  
  ggplot(x, aes(x=row,
                y=as.numeric(col),
		fill=color,
		size=(modal.pct))) +
	geom_point(shape=23, alpha=1) +
	scale_x_continuous(breaks=seq(1, 20)) + 
	scale_fill_manual(values=cols) +
        #scale_size_continuous(range=c(0, ym), limits=c(0, 1)) +
	scale_radius(range=c(0, ym)) + 
        theme(axis.ticks = element_line(size = 1, color="black")) +
        blanky +
	theme(legend.position="none", 
              axis.text.x = element_text(color="black", size=15),
	      axis.text.y = element_text(color="black", size=15),
              axis.title.x = element_blank(),
	      axis.title.y = element_blank()) + 
        guides(fill = "none") +
	scale_y_continuous(limits=c(0, 9),
	                   breaks=seq(1, 8),
			   labels=rev(c("A", "B", "C", "D", "E", "F", "G", "H")))
  
}

diamond(filter(x.sum, Language == "English", Task == "Open")) 
ggsave("output/english_diamond.png", width=710/100 * 1.18, height=249/100 * 1.4)

diamond(filter(x.sum, Language == "Spanish", Task == "Open")) 
ggsave("output/spanish_diamond.png", width=710/100 * 1.18, height=249/100 * 1.4)

diamond(filter(x.sum, Language == "Tsimane", Task == "Open"))
ggsave("output/tsimane_diamond.png", width=710/100 * 1.18, height=249/100 * 1.4)


diamond(filter(x.sum, Language == "English", Task == "Fixed")) 
ggsave("output/english_fixed_diamond.png", width=710/100 * 1.18, height=249/100 * 1.4)

diamond(filter(x.sum, Language == "Spanish", Task == "Fixed")) 
ggsave("output/spanish_fixed_diamond.png", width=710/100 * 1.18, height=249/100 * 1.4)

diamond(filter(x.sum, Language == "Tsimane", Task == "Fixed"))
ggsave("output/tsimane_fixed_diamond.png", width=710/100 * 1.18, height=249/100 * 1.4)

