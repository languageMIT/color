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
library(Hmisc)
library(lme4)
library(tidyverse)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# We've determined the Julian-style data contains duplicated subjects so we should
# use the Kyle-style data. (Meeting on 2016-06-02.)
LABELLING_DATA = "Kyle"  # "Julian", "Kyle", or "Richard"

SPECIAL_FOCAL_COLORS = c("blue", "green", "red", "yellow")

TSIMANE_COLORS = c("black", "white", "red", "blue", "green", "yellow2",
                   "grey", "purple / violet", "orange", "brown", "yellowish",
		   "yellow", "orange", "greenish", "brownish")

ENGLISH_COLORS = c("black", "white", "red", "blue", "green", "yellow",
                   "grey", "purple / violet", "orange", "brown", "pink", "celeste")

RESTRICT_COOL = F
FORBIDDEN_COLS = c(7, 9, 11, 13, 15)
RESTRICT_FOCAL = F
PRIOR = "Uniform" # "Uniform", "Foreground", "Background", ...
RESTRICT_CIELAB = F
BLACK_WHITE = F
INTERPOLATE_ENTROPY = F

get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

is_close = function(x, y) {
  abs(x - y) < .00001
}

blanky = theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
	       panel.background = element_blank(),
	       axis.line = element_line(colour = "black"),
	       axis.text.x = element_text(size=18),
	       axis.text.y=element_text(size=18),
	       axis.title.y=element_text(size=18))


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

 ColorData = read.csv("color labeling tsimane english spanish oct 2015 v4.csv", stringsAsFactors=F) %>%
   rename(grid_location=location) %>%
   select(subject, color, Language, Task, grid_location) %>%
   unite(Experiment, Language, Task) %>%
   filter(grid_location %in% sparse_chips_to_use)

} else if (LABELLING_DATA == "Richard") {

 load("sparse_chips_to_use.rda")

 ColorData = read.csv("colorlabeling.csv") %>%
   rename(grid_location=location) %>%
   filter(!Toss) %>%
   select(subject, color, Language, Task, grid_location) %>%
   unite(Experiment, Language, Task) %>%
   filter(grid_location %in% sparse_chips_to_use) 
}

if (RESTRICT_COOL) {
   ColorData %>%
     mutate(rc=CodeToRowColumn(grid_location)) %>%
     separate(rc, into=c("Row", "Column"), sep="_") %>%
     mutate(Row=as.numeric(Row), Column=as.numeric(Column)) %>%
     filter(!(Column %in% FORBIDDEN_COLS)) %>%
     mutate(grid_location=RowColumnToCode(Row, Column)) %>%
     select(-Row, -Column) -> ColorData    
}



m = read.csv("Munsell_WCS_codes.csv") %>%
  rename(grid_location=Our.code,
         MunsellCode=OtHer.code,
	 WCSCode=WCS)


## Add black/white ---------------------------------------
## For the color naming task, everyone answered consistently for white/black.
## Add this data in two new dummy grid locations: I1 (black), I2 (white)
BLACK_LOCATION = "I1"
WHITE_LOCATION = "I2"

Black = ColorData %>%
    group_by(subject, Experiment) %>%
    summarise(color="1", grid_location=BLACK_LOCATION) %>%
    ungroup()

White = ColorData %>%
    group_by(subject, Experiment) %>%
    summarise(color="2", grid_location=WHITE_LOCATION) %>%
    ungroup()

if (BLACK_WHITE) {
    ColorData = ColorData %>% bind_rows(Black, White)
}

# Load focal location data ----------------------------

d5 = read.csv("Tsimane2015_2_focal_colors_object_colors.csv", stringsAsFactors=F)
d4 = read.csv("Tsimane2014_focal_colors_object_colors.csv", stringsAsFactors=F)

d4 %>%
  filter(Toss == "N") %>%
  gather(term, Code, 13:30) %>%
  select(Subject, term, Code) %>%
  filter(Code != "") %>%
  separate(term, into=c("tmp1", "tmp2", "term"), sep="_") %>%
  mutate(year=2014) %>%
  select(-tmp1, -tmp2) -> d4_choices

d5 %>%
  gather(term, Code, 14:25) %>%
  select(Subject.number, term, Code) %>%
  filter(Code != "") %>%
  separate(term, into=c("tmp1", "tmp2", "term", "mod"), sep="_") %>%
  filter(is.na(mod)) %>% # disregard _2 terms
  mutate(Subject = Subject.number) %>%
  mutate(year=2015) %>%
  select(-tmp1, -tmp2, -mod, -Subject.number) -> d5_choices

d = rbind(d4_choices, d5_choices) %>% filter(!is.na(term))

focal_ts = d %>% mutate(Language = "Tsimane")
focal_ts = filter(focal_ts, !str_detect(term, "_")) # drop _2 colors
assert((focal_ts %>% select(year, Subject) %>% unique() %>% nrow()) == 99)
focal_ts = focal_ts %>% select(-year)

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
    summarise(freq=n()) %>%
    ungroup() %>%
  group_by(Language, term) %>%
    mutate(Z=sum(freq)) %>%
    mutate(p=freq/Z) %>%
    select(-freq) %>%
    write.csv("output/focal_P_ChipGivenWord.csv")

write.csv(focal_locations, "output/focal_locations.csv")

focal_locations = focal_locations %>%
  mutate(RowColumn=CodeToRowColumn(Code)) %>%
  separate(RowColumn, into=c("Row", "Column"), sep="_") %>%
  mutate(Row=as.numeric(Row), Column=as.numeric(Column))

focal_chips = focal_locations %>%
  group_by(Language, term) %>%
    summarise(Code=get_mode(Code)) %>%
    ungroup() 

if (BLACK_WHITE) {
  focal_b = focal_chips %>%
    group_by(Language) %>%
    summarise(term="1", Code="I1") %>%
    ungroup()
  
  
  focal_w = focal_chips %>%
    group_by(Language) %>%
    summarise(term="2", Code="I2") %>%
    ungroup()
  
  focal_chips = focal_chips %>%
    rbind(focal_b, focal_w)
}

# For each focal term, how many subjects know that term?
# first need to convert the color indices in ColorData into terms
ColorData %>%
  separate(Experiment, into=c("Language", "Task"), sep="_") %>%
  filter(Task == "Open") %>%
  select(subject, color, Language) %>%
  mutate(color=as.character(color)) %>%
  group_by(Language) %>%
    mutate(color=ifelse(Language == "Tsimane",
                        TSIMANE_COLORS[as.numeric(as.character(color))],
    		        ENGLISH_COLORS[as.numeric(as.character(color))])) %>%
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


cielab_chips = c("A6", "A8", "B9", "C2", "C10", "D19", "E12", "E14", "E16", "E18", "F3", "F5", "F7", "F11", "F13", "F15", "F17", "F19", "G2", "G8", "G14", "G16", "H1", "H15")
assert(length(cielab_chips) == 24)


if (RESTRICT_FOCAL) {

   #allowed_chips = filter(focal_chips, Language == "Spanish")$Code %>% as.character()
   allowed_chips = c("E14", "H3", "E8", "A14", "E2", "D19", "G18", "F1", "B5")
   ColorData = filter(ColorData, as.character(grid_location) %in% allowed_chips)
}


if (RESTRICT_CIELAB) {
   ColorData = filter(ColorData, as.character(grid_location) %in% cielab_chips)
}

# Read demographics --------------------------------------------

demo2014 = read.csv("demographics2014.csv", stringsAsFactors=F) %>% mutate(Year=2014, Task="Open")
demo2015 = read.csv("demographics2015.csv", stringsAsFactors=F) %>% mutate(Year=2015, Task="Fixed")
demo2015 %>%
  rename(Subject=Subject..) %>%
  mutate(year=2015) %>%
  select(Subject, Age, Sex, Spanish, Year, Task, Education) %>%
  rbind(demo2014 %>% select(Subject, Age, Sex, Spanish, Year, Task, Education)) -> ts_demo

ColorData %>%
  separate(Experiment, into=c("Language", "Task"), sep="_") %>%
  filter(Language == "Tsimane") %>%
  rename(Subject=subject) %>%
  mutate(Subject=as.character(Subject)) %>%
  inner_join(ts_demo) -> ColorDataWithDemographics
  

# Compute conditional entropy ----------------------------------

pixels = read.csv("foreground_background_pixels.csv")
in_grid = ColorData %>% select(grid_location) %>% unique()
if (PRIOR == "Uniform") {
  prior = in_grid %>% mutate(Prior=1/n())
} else if (PRIOR == "Foreground") {
  prior = in_grid %>%
    inner_join(pixels) %>%
    mutate(PriorStar=(foreground+1)/(foreground+background+2),
           Z=sum(PriorStar),
	   Prior=PriorStar/Z) %>%
    select(-background, -foreground)
} else if (PRIOR == "Background") {
  prior = in_grid %>%
    inner_join(pixels) %>%
    mutate(PriorStar=(background+1)/(foreground+background+2),
           Z=sum(PriorStar),
	   Prior=PriorStar/Z) %>%    
    select(-background, -foreground)
} else {
  assert(PRIOR %in% c("Uniform", "Foreground", "Background"))
}

# Get the joint distribution of words and chips

ComputeConditionalEntropy_Old = function(ColorData) {
    # Old function which I keep around as a reference

P_WordGivenChip = ColorData %>%
  group_by(Experiment, grid_location, color) %>%
    summarise(Frequency=n()) %>%
    ungroup() %>%
  group_by(Experiment, grid_location) %>%
    mutate(Z=sum(Frequency)) %>%
    ungroup() %>%
  mutate(Probability=Frequency/Z) %>%
  select(Experiment, grid_location, color, Probability, Frequency) %>%
  unique()

if (PRIOR == "Uniform") {    
  in_grid = ColorData %>% select(grid_location) %>% distinct()
  prior = in_grid %>% mutate(Prior=1/n())
  }
    
  P_ChipGivenWord = P_WordGivenChip %>%  # we have p(w|c) for each w and c
    rename(Likelihood=Frequency) %>%     # use Frequency rather than Probability to increase precision
    inner_join(prior) %>%                # now we merge in p(c)
    mutate(Pstar=Likelihood*Prior) %>%   # calculate p(c)p(w|c), which is prop. to p(w, c)
    group_by(Experiment, color) %>%                  
      mutate(Z=sum(Pstar)) %>%           # now calculate p(w) = \sum_c p(w, c)
      ungroup() %>%        
    mutate(Probability=Pstar/Z) %>%      # p(c|w) = p(c)p(w|c) / p(w)
    select(Experiment, grid_location, color, Probability)


# If uniform prior, do it this way to check
if (PRIOR == "Uniform") {

  P_ChipGivenWord_Direct = ColorData %>%
    group_by(Experiment, grid_location, color) %>%
      summarise(Frequency=n()) %>%
      ungroup() %>%
    group_by(Experiment, color) %>%
      mutate(Z=sum(Frequency)) %>%
      ungroup() %>%
    mutate(Probability=Frequency/Z) %>%
    select(Experiment, grid_location, color, Probability) %>%
    distinct()

  assert(is_close(P_ChipGivenWord_Direct$Probability, P_ChipGivenWord$Probability))
}

# Get each chip's score
GetChipScore<-function(x){
  # Sum over all words. This is very slow, so I rewrote this function below using dplyr stuff.
  Val=0
  for (CurrCol in unique(x$color)){
    ChipWord=filter(P_ChipGivenWord,
                    color==CurrCol,
		    grid_location==x$grid_location[1],
		    Experiment==x$Experiment[1])$Probability
    WordChip=filter(P_WordGivenChip,
                    grid_location==x$grid_location[1],
		    color==CurrCol,
		    Experiment==x$Experiment[1])$Probability
    Val=Val+WordChip*log2(1.0/ChipWord)
  }
  return(Val)
}

ChipEntropies <- plyr::ddply(ColorData, c("Experiment","grid_location"), GetChipScore)
names(ChipEntropies)<-c("Experiment","grid_location","Entropy")
ChipNo<-length(unique(ChipEntropies$grid_location))
ChipEntropies <- ChipEntropies %>% inner_join(prior)
ChipEntropies
}

ComputeConditionalEntropy = function(ColorData) {

P_WordGivenChip = ColorData %>%
  group_by(Experiment, grid_location, color) %>%
    summarise(Frequency=n()) %>%
    ungroup() %>%
  group_by(Experiment, grid_location) %>%
    mutate(Z=sum(Frequency)) %>%
    ungroup() %>%
  mutate(Probability=Frequency/Z) %>%
  select(Experiment, grid_location, color, Probability, Frequency) %>%
  unique()

if (PRIOR == "Uniform") {    
  in_grid = ColorData %>% select(grid_location) %>% distinct()
  prior = in_grid %>% mutate(Prior=1/n())
  }
    
  P_ChipGivenWord = P_WordGivenChip %>%  # we have p(w|c) for each w and c
    rename(Likelihood=Frequency) %>%     # use Frequency rather than Probability to increase precision
    inner_join(prior) %>%                # now we merge in p(c)
    mutate(Pstar=Likelihood*Prior) %>%   # calculate p(c)p(w|c), which is prop. to p(w, c)
    group_by(Experiment, color) %>%                  
      mutate(Z=sum(Pstar)) %>%           # now calculate p(w) = \sum_c p(w, c)
      ungroup() %>%        
    mutate(Probability=Pstar/Z) %>%      # p(c|w) = p(c)p(w|c) / p(w)
    select(Experiment, grid_location, color, Probability)


# If uniform prior, do it this way to check
if (PRIOR == "Uniform") {

  P_ChipGivenWord_Direct = ColorData %>%
    group_by(Experiment, grid_location, color) %>%
      summarise(Frequency=n()) %>%
      ungroup() %>%
    group_by(Experiment, color) %>%
      mutate(Z=sum(Frequency)) %>%
      ungroup() %>%
    mutate(Probability=Frequency/Z) %>%
    select(Experiment, grid_location, color, Probability) %>%
    distinct()

  assert(is_close(P_ChipGivenWord_Direct$Probability, P_ChipGivenWord$Probability))
}

    ChipEntropies = P_WordGivenChip %>%
        select(Experiment, grid_location, color, Probability) %>%
        rename(P_WordGivenChip=Probability) %>%
        inner_join(P_ChipGivenWord) %>%
        rename(P_ChipGivenWord=Probability) %>%
        group_by(Experiment, grid_location) %>%
          summarise(Entropy=sum(P_WordGivenChip * log2(1/P_ChipGivenWord))) %>%
          ungroup()

    ChipEntropies %>% inner_join(prior)
        
}

ChipEntropiesOld = ComputeConditionalEntropy_Old(ColorData)
ChipEntropies = ComputeConditionalEntropy(ColorData)
assert(all(is_close(ChipEntropiesOld$Entropy, ChipEntropies$Entropy)))

LanguageEntropies <- ChipEntropies %>%
  group_by(Experiment) %>%
  summarise(Uncertainty=sum(Entropy*Prior))

write.csv(ChipEntropies, "output/Latest_ChipEntropies.csv")
write.csv(LanguageEntropies, "output/Latest_EntropyByLanguage.csv")

ChipEntropies %>%
  select(-Prior) %>%
  group_by(Experiment) %>%
    arrange(Entropy) %>%
    mutate(X=1:n()) %>%
    ungroup() %>%
  select(-Entropy) %>%
  spread(Experiment, grid_location) %>%
  select(-X) %>%
  select(English_Open, Spanish_Open, Tsimane_Open, English_Fixed, Spanish_Fixed, Tsimane_Fixed) %>%
  write.csv("output/Latest_ChipsRankOrdered.csv")

ChipEntropies <- ChipEntropies %>%
  separate(Experiment, into=c("Language", "Type"), sep="_") %>%
  mutate(RowColumn=CodeToRowColumn(grid_location)) %>%
  separate(RowColumn, into=c("Row", "Column"), sep="_") %>%
  mutate(Row=as.numeric(Row), Column=as.numeric(Column))

# Bootstrap conditional entropy ----------------------------------

BootstrapResampleColorData = function(ColorData) {
    ColorData %>%
        mutate(color=as.character(color)) %>%
        group_by(Experiment, grid_location) %>%
          sample_frac(replace=T) %>%
          ungroup()
}

num_bootstrap_samples = 100

ResampledChipEntropies = NA
for (i in 1:num_bootstrap_samples) {
    new_data = BootstrapResampleColorData(ColorData)
    new_entropies = ComputeConditionalEntropy(new_data) %>% mutate(sample=i)
    if (i == 1) {
        ResampledChipEntropies = new_entropies
    } else {
        ResampledChipEntropies = rbind(ResampledChipEntropies, new_entropies)
    }
    ResampledChipEntropies %>% separate(Experiment, into=c("Language", "Type"), sep="_")
}

ResampledChipEntropies %>%
  group_by(Experiment, grid_location) %>%
    summarise(lower=quantile(Entropy, .05), upper=quantile(Entropy, .95)) %>%
    ungroup() %>%
  separate(Experiment, into=c("Language", "Type")) %>%
  inner_join(ChipEntropies) -> ChipEntropiesWithCI


# Spearman correlations -------------------------------

en_es = cor.test(filter(ChipEntropies, Language == "English")$Entropy,
                 filter(ChipEntropies, Language == "Spanish")$Entropy,
          	 method='spearman')

en_ts = cor.test(filter(ChipEntropies, Language == "English")$Entropy,
                 filter(ChipEntropies, Language == "Tsimane")$Entropy,
          	 method='spearman')

es_ts = cor.test(filter(ChipEntropies, Language == "Spanish")$Entropy,
                 filter(ChipEntropies, Language == "Tsimane")$Entropy,
          	 method='spearman')


# Modal informativity --------------------------

# Suppose every chip is just named according to its modal label.
ModalColorData = ColorData %>%
  group_by(Experiment, grid_location) %>%
    summarise(color=get_mode(color)) %>%
    ungroup()

ModalChipEntropies = ComputeConditionalEntropy(ModalColorData)

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
WCS = answers %>%
  select(Lang_Name, FielConcat, Term, Subject) %>%
  rename(Language=Lang_Name, Chip=FielConcat, Color=Term)

WCS$Color<-as.character(WCS$Color)
# some colors are called NA and R thinks they're missing values
WCS[which(is.na(WCS$Color)),]$Color="NA"


#Standardize with main dataset
names(WCS)<-c("Experiment","grid_location","color", "subject")
WCS$grid_location<-as.character(WCS$grid_location)

WCS_with_subjects = WCS
WCS = select(WCS, -subject)

# Reduce WCS chips to those in our other analyses -------------------------------------
# Load code conversion chart
ConversionChart<-read.csv("Munsell_WCS_codes.csv") %>% tbl_df
names(ConversionChart) = c("BoliviaCode", "MainCode", "WCSCode")

ConversionChart$BoliviaCode<-as.character(ConversionChart$BoliviaCode)
ConversionChart$WCSCode<-as.character(ConversionChart$WCSCode)
ConversionChart <- filter(ConversionChart,BoliviaCode %in% unique(ColorData$grid_location)) %>%
  dplyr::rename(grid_location=WCSCode)

WCS_unfiltered = WCS
WCS <- filter(WCS,grid_location %in% unique(ConversionChart$grid_location))
answers <- filter(answers, FielConcat %in% unique(ConversionChart$grid_location))

print("Computing full WCS conditional entropy (may be slow)...")
FullWCSChipEntropies = WCS_unfiltered %>%
    ComputeConditionalEntropy()
print("Done")

wcs_to_srgb_d = read.csv("sRGB_Munsell320.csv") %>%
    select(R, G, B, Row, Column) %>%
    mutate(hex=rgb(R, G, B, maxColorValue=255)) %>%
    unite(Code, Row, Column, sep="")

wcs_to_srgb = wcs_to_srgb_d$hex
wcs_to_srgb = setNames(wcs_to_srgb, wcs_to_srgb_d$Code)

FullWCSChipEntropies %>%
    select(-Prior) %>%
    spread(Experiment, Entropy) %>%
    inner_join(wcs_to_srgb_d %>% rename(grid_location=Code)) %>%
    write.csv("output/Full_WCS_ChipEntropies.csv")
           

assert(length(unique(WCS$grid_location)) == 80 | RESTRICT_COOL | RESTRICT_FOCAL| RESTRICT_CIELAB | BLACK_WHITE)

WCS <- inner_join(WCS,
                  ConversionChart,
		  by=c("grid_location")) %>%
         select(-MainCode,-grid_location) %>%
	 rename(grid_location=BoliviaCode)
	 
WCS %>% dplyr::group_by(Experiment) %>% dplyr::summarise(Chips=length(unique(grid_location))) %>%
  dplyr::select(Chips) %>% table # good!
  
if (RESTRICT_COOL) {
   WCS %>%
     mutate(rc=CodeToRowColumn(grid_location)) %>%
     separate(rc, into=c("Row", "Column"), sep="_") %>%
     mutate(Row=as.numeric(Row), Column=as.numeric(Column)) %>%
     filter(!(Column %in% FORBIDDEN_COLS)) %>%
     mutate(grid_location=RowColumnToCode(Row, Column)) %>%
     select(-Row, -Column) -> WCS
}

if (RESTRICT_CIELAB) {
   WCS = filter(WCS, as.character(grid_location) %in% cielab_chips)
}

# Compute WCS mode conditional entropy ----------------------
# or don't

#WCS_Modes<-tbl_df(plyr::ddply(WCS,c("Experiment","grid_location"),function(x){return(data.frame(table(x$color)))}))
#names(WCS_Modes)<-c("Experiment","grid_location","color","frequency")
#
## Get modes
#WCS_Modes <- WCS_Modes %>% dplyr::group_by(Experiment,grid_location) %>% dplyr::top_n(1,frequency)
## Check total number of color words
##res<-WCS_Modes %>% dplyr::group_by(Experiment) %>% dplyr::summarise(ColorWords=length(unique(color)))
## Delete ties
#WCS_Modes <- WCS_Modes %>% dplyr::group_by(Experiment,grid_location) %>% dplyr::summarise(color=color[1])
#
#WCS_Mode_ChipGivenWord <- plyr::ddply(WCS_Modes,c("Experiment","color"),function(x){return(as.data.frame(prop.table(table(x$grid_location))))})
#names(WCS_Mode_ChipGivenWord)<-c("Experiment","color","grid_location","Probability")
#
#WCS_Modes<-full_join(WCS_Modes,WCS_Mode_ChipGivenWord)
#WCS_Modes <- WCS_Modes %>% mutate(Prior=1/length(unique(WCS_Modes$grid_location)))
#
#WordNumbers <- WCS_Modes %>% dplyr::group_by(Experiment) %>% dplyr::summarise(WordNumber = length(unique(color)))
#
#WCS_Modes <- WCS_Modes %>% dplyr::group_by(Experiment) %>% dplyr::summarise(ConEnt=sum(Probability*log2(1.0/Probability)))
#
#WCS_Modes <- full_join(WordNumbers, WCS_Modes)

#WCS_Modes %>% filter(WordNumber %in% c(7,8)) %>%
#  ggplot(aes(x=WordNumber,y=ConEnt,label=Experiment))+geom_text()


# Compute WCS conditional entropy --------------------------

WCS = filter(WCS, !is.na(Experiment))
print("Computing WCS chip entropies...")
WCS_ChipEntropies = ComputeConditionalEntropy(WCS)
print("Done")

write.csv(WCS_ChipEntropies, "output/Latest_WCS_ChipEntropies.csv")

All_ChipEntropies = ChipEntropies %>%
  filter(Type == "Open") %>%
  select(Language, grid_location, Entropy) %>%
  rbind(WCS_ChipEntropies %>%
        rename(Language=Experiment) %>%
	select(Language, grid_location, Entropy)) %>%
  spread(Language, Entropy) %>%
  inner_join(m) %>%
  rename(OurCode=grid_location)

All_ChipEntropies %>%
  write.csv("output/Latest_All_ChipEntropies.csv")

tap = function(x) {
  print(x)
  return(x)
}

interpolate_entropy = function(entropy, row, column) {
  d = data.frame(entropy, row, column)
  names(d) = c("entropy", "row", "column")
  d$interpolated_entropy = d$entropy
  for(i in 1:nrow(d)) {
    if(is.na(d[i,]$entropy)) {
      row = d[i,]$row
      column = d[i,]$column
      below = d$row==row+1 & d$column==column
      above = d$row==row-1 & d$column==column
      left = d$row==row & d$column==column-1
      right = d$row==row & d$column==column+1
      #print("row:")
      #print(row)
      #print("column:")
      #print(column)
      #print("interpolands:")
      d[i,]$interpolated_entropy = mean(
        c(d[above,]$entropy,
        d[below,]$entropy,
        d[left,]$entropy,
        d[right,]$entropy),
        na.rm=T
      )
    }
  }
  d$interpolated_entropy
}



if (INTERPOLATE_ENTROPY) { # optional because it's slow
  
  NUM_ROWS = 9
  NUM_COLS = 20

  rows = c()
  for(i in 1:NUM_ROWS) {
    rows = c(rows, rep(LETTERS[i], NUM_COLS))
  }
  columns = rep(1:NUM_COLS, NUM_ROWS)
  
  grid_locations = data.frame(rows, columns)
  names(grid_locations) = c("row", "column")
  grid_locations = grid_locations %>% unite(OurCode, row, column, sep="")
  
## Interpolate to the full WCS grid
All_ChipEntropies_Interpolated = grid_locations %>%
  full_join(All_ChipEntropies, by="OurCode") %>%
  filter(!(WCSCode %in% c("A0", "J0"))) %>% # white and black will not participate in the interpolation
  mutate(RowColumn=CodeToRowColumn(OurCode)) %>%
  separate(RowColumn, c("row", "column"), sep="_") %>%
  mutate(row=as.numeric(row), column=as.numeric(column)) %>%
  gather(Language, Entropy, -OurCode, -WCSCode, -row, -column) %>%
  group_by(Language) %>%
    mutate(InterpolatedEntropy=interpolate_entropy(as.numeric(Entropy), row, column)) %>%
    ungroup() %>%
  select(-Entropy, -row, -column) %>%
  spread(Language, InterpolatedEntropy) %>%
  select(-MunsellCode) %>%
  bind_rows(filter(All_ChipEntropies, WCSCode %in% c("A0", "J0")) %>% select(-MunsellCode)) %>%
  full_join(select(ConversionChart, BoliviaCode, MainCode) %>% rename(OurCode=BoliviaCode)) %>%
  rename(MunsellCode=MainCode)

write.csv(All_ChipEntropies_Interpolated, "output/All_ChipEntropies_Interpolated.csv")

}







# Numbers of subjects -----------------------------------------------

rts = read.csv("english_rts.csv") %>%
  filter(Native_speaker == "Y") %>%
  select(-Native_speaker, -X) %>%
  mutate(Language="English") %>%
  rbind(read.csv("tsimane_rts.csv") %>% mutate(Language="Tsimane"))

sedivy = read.csv("sedivy_data_fixed.csv")

sedivy %>%
  rename(Task=Presentation_type) %>%
  group_by(Language, Task) %>%
    summarise(NumSubjects=length(unique(subject))) %>%
    ungroup() ->
  sedivy_num_subjects

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

num_subjects = rbind(labelling_num_subjects,
                     objects_num_subjects,
		     rts_num_subjects,
		     sedivy_num_subjects) %>%
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

LETTERS = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")

load("cc.rda") # Load the chip code -> hex code conversion table
cc["I1"] = "#000000"
cc["I2"] = "#F5F5F5"

ts_color_indices = focal_locations[focal_locations$Language == "Tsimane",]$term %>% as.numeric()
ts_colors = TSIMANE_COLORS[ts_color_indices]
focal_locations[focal_locations$Language == "Tsimane",]$term = ts_colors

special_focal_locations = filter(focal_locations, term %in% SPECIAL_FOCAL_COLORS)
special_focal_density = group_by(special_focal_locations, Language, Row, Column, term) %>% summarise(value=n())
write.csv(special_focal_density, "output/special_focal_density.csv") # to be used by contours.py

# If E10 is missing, interpolate it from E8, E12, D9, D11, F9, F11
if (!RESTRICT_CIELAB & !("E10" %in% as.character(ChipEntropies$grid_location))) {
  interp_sources = c("E8", "E12", "D9", "D11", "F9", "F11")
  interp_rows = ChipEntropies %>%
    filter(grid_location %in% interp_sources) %>%
    group_by(Language, Type) %>%
      summarise(Entropy=mean(Entropy)) %>%
      ungroup() %>%
    inner_join(prior)
  
  ChipEntropiesInterp = ChipEntropies %>%
    select(Language, Type, Entropy, grid_location, Prior) %>%
    rbind(interp_rows) %>%
    mutate(RowColumn=CodeToRowColumn(grid_location)) %>%
    separate(RowColumn, into=c("Row", "Column"), sep="_") %>%
    mutate(Row=as.numeric(Row), Column=as.numeric(Column))
} else {
  ChipEntropiesInterp = ChipEntropies
}

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

d_open = ChipEntropiesWithCI %>%
  filter(Type == "Open") %>%
  select(Language, grid_location, Entropy) %>%
  rbind(WCS_ChipEntropies %>%
        rename(Language=Experiment) %>%
	select(Language, grid_location, Entropy)) %>%
  rename(Code=grid_location)
  
  

# For some reason group_by does the wrong thing in combination with rank(), so do it the bad old way...
d$score_rank = NA
for (lang in unique(d$Language)) {
  d[d$Language == lang,]$score_rank = rank(d[d$Language == lang,]$Entropy, ties.method=c("random"))
}

d_open$score_rank = NA
for (lang in unique(d$Language)) {
  d_open[d_open$Language == lang,]$score_rank = rank(d_open[d_open$Language == lang,]$Entropy,
                                                     ties.method=c("random"))
}

full_wcs = FullWCSChipEntropies %>%
    rename(Language=Experiment) %>%
    select(Language, grid_location, Entropy)

full_wcs$score_rank = NA
for (lang in unique(full_wcs$Language)) {
    full_wcs[full_wcs$Language == lang,]$score_rank = rank(full_wcs[full_wcs$Language == lang,]$Entropy, ties.method=c("random"))
}

d_ci = ChipEntropiesWithCI %>%
    filter(Type == "Open") %>%
    select(-Row, -Column, Prior) %>%
    rename(Code=grid_location)

d_ci$score_rank = NA
for (lang in unique(d_ci$Language)) {
  d_ci[d_ci$Language == lang,]$score_rank = rank(d_ci[d_ci$Language == lang,]$Entropy, ties.method=c("random"))
}

d_ci$is_focal = ""
for (i in 1:nrow(d_ci)) {
  lang = first(d_ci[i,]$Language)
  if (d_ci[i,]$Code %in% focal_chips[focal_chips$Language == lang,]$Code) {
    d_ci[i,]$is_focal = "*"
  }
}

d_ci %>%
  ggplot(aes(x=score_rank,
             y=Entropy,
	     xmin=score_rank-.5,
	     xmax=score_rank+.5,
	     ymin=Entropy-.1,
	     ymax=Entropy+.1,
	     color=Code,
	     fill=Code,
	     label=is_focal,
	     )) +
    geom_errorbar(aes(ymin=lower, ymax=upper), color="light grey", size=.3) +
    geom_rect() +
    geom_text(position=position_nudge(y=1/7), size=5, color="black") +
    blanky +
    scale_color_manual(values=cc, guide=F) +
    scale_fill_manual(values=cc, guide=F) +
    facet_wrap(~Language, ncol=4) +
    xlab("Chips (rank ordered)") + 
    ylab("Average surprisal") +
    ylim(0, NA)

ggsave("output/snake_open_with_ci.pdf", height=4, width=12)

d_ci %>%
  ggplot(aes(x=score_rank,
             y=Entropy,
	     xmin=score_rank-.5,
	     xmax=score_rank+.5,
	     ymin=Entropy-.1,
	     ymax=Entropy+.1,
	     color=Code,
	     fill=Code)) +
    geom_rect() +
    blanky +
    scale_color_manual(values=cc, guide=F) +
    scale_fill_manual(values=cc, guide=F) +
    facet_wrap(~Language, ncol=4) +
    xlab("Chips (rank ordered)") + 
    ylab("Average surprisal") +
    ylim(0, NA)

ggsave("output/snake_open.pdf", height=3, width=10)

d_open %>%
    select(-score_rank) %>%
    filter(Language %in% c("English", "Spanish", "Tsimane")) %>%
    spread(Language, Entropy) %>%
    ggplot(aes(x=English,
               y=Tsimane,
	       xmin=English-.05,
	       xmax=English+.05,
	       ymin=Tsimane-.05,
	       ymax=Tsimane+.05,
	       color=Code,
	       fill=Code)) +
    geom_rect() +
    theme_bw(18) +
    blanky +
    scale_color_manual(values=cc, guide=F) +
    scale_fill_manual(values=cc, guide=F) +
    xlab("Tsimane' surprisal") +
    ylab("English surprisal")

ggsave("output/tsimane_vs_english_surprisal.pdf")

d_open %>%
    select(-score_rank) %>%
    filter(Language %in% c("English", "Spanish", "Tsimane")) %>%
    spread(Language, Entropy) %>%
    ggplot(aes(x=Spanish,
               y=Tsimane,
	       xmin=Spanish-.05,
	       xmax=Spanish+.05,
	       ymin=Tsimane-.05,
	       ymax=Tsimane+.05,
	       color=Code,
	       fill=Code)) +
    theme_bw(18) +
    blanky + 
    scale_color_manual(values=cc, guide=F) +
    scale_fill_manual(values=cc, guide=F) +
    geom_rect() +
    xlab("Tsimane' surprisal") +
    ylab("Spanish surprisal")

ggsave("output/tsimane_vs_spanish_surprisal.pdf")

d_open %>%
    select(-score_rank) %>%
    filter(Language %in% c("English", "Spanish", "Tsimane")) %>%
    spread(Language, Entropy) %>%
    ggplot(aes(x=English,
               y=Tsimane,
	       xmin=English-.05,
	       xmax=English+.05,
	       ymin=Tsimane-.05,
	       ymax=Tsimane+.05,
	       color=Code,
	       fill=Code)) +
    geom_rect() +
    theme_bw(18) +
    blanky +
    scale_color_manual(values=cc, guide=F) +
    scale_fill_manual(values=cc, guide=F) +
    xlab("Tsimane' surprisal") +
    ylab("English surprisal")

ggsave("output/tsimane_vs_english_surprisal.pdf")

ChipEntropies %>%
  rename(Code=grid_location) %>%
  spread(Type, Entropy) %>%
  ggplot(aes(x=Open, y=Fixed,
             xmin=Open-.05, xmax=Open+.05,
	     ymin=Fixed-.05, ymax=Fixed+.05,
	     color=Code, fill=Code)) +
  geom_rect() +
  theme_bw(18) +
  blanky +
  scale_color_manual(values=cc, guide=F) +
  scale_fill_manual(values=cc, guide=F) +
  xlab("Open task surprisal") +
  ylab("Fixed task surprisal") +
  facet_wrap(~Language) +
  xlim(0, 6) +
  ylim(0, 6)

ggsave("output/open_vs_fixed_surprisal.pdf", height=4, width=10)

print("Open-Fixed Correlations:")
ChipEntropies %>%
  spread(Type, Entropy) %>%
  group_by(Language) %>%
    summarise(r=cor.test(Open, Fixed, method="pearson")[[4]],
              rho=cor.test(Open, Fixed, method="spearman")[[4]]) %>%
  print()	      
	    


# Tapestry plot ----------------------------------------

d_to_plot = d_open %>%
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

full_wcs_to_plot = full_wcs %>%
    mutate(Language=factor(Language)) %>%
    group_by(Language) %>%
        mutate(lang_score=sum(Entropy), min_score=min(Entropy), max_score=max(Entropy)) %>%
        ungroup() %>%
    mutate(Language=reorder(factor(Language), -lang_score))

wcs_to_srgb_d = read.csv("sRGB_Munsell320.csv") %>%
    mutate(hex=rgb(R, G, B, maxColorValue=255)) %>%
    unite(Code, Row, Column, sep="")

wcs_to_srgb = wcs_to_srgb_d$hex
wcs_to_srgb = setNames(wcs_to_srgb, wcs_to_srgb_d$Code)

full_wcs_to_plot %>%
  ggplot(aes(x=score_rank, y=Language, color=grid_location, fill=grid_location)) +
  scale_color_manual(values=wcs_to_srgb, guide=F) +
  scale_fill_manual(values=wcs_to_srgb, guide=F) +
  geom_tile() +
  xlab("Chips by increasing average surprisal =>") +
  ylab("<= Languages by increasing total surprisal") +
  theme_classic() +
  theme(axis.ticks=element_blank(),
        axis.text.x=element_blank(),
	axis.text.y = element_text(size=6))

ggsave("output/Full_WCS_tapestry.pdf", width=10, height=8)





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

d_open %>%
    filter(Language %in% c("English", "Spanish", "Tsimane")) %>%
    ggplot(aes(x=score_rank, y=Language, color=Code, fill=Code)) +
      scale_color_manual(values=cc, guide=F) +
      scale_fill_manual(values=cc, guide=F) +
      geom_tile() +
      theme_bw() +
      theme(axis.text.x = element_blank(),
	    axis.ticks.x = element_blank()) +
      ylab("") +
    xlab("Colors by increasing average surprisal =>")

ggsave("output/WCS_strip_est.pdf", width=30, height=3)
    

d_open %>%
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
      filter(Row != "I") %>% # black and white are useless for this part so drop them
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
                        
                              !(grid_location %in% c("I1", "I2")),
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
	scale_radius(range=c(0, ym), limits=c(0, 1)) + 
        theme(axis.ticks = element_line(size = 1, color="black")) +
        blanky +
	theme(axis.text.x = element_text(color="black", size=15),
	      axis.text.y = element_text(color="black", size=15),
              axis.title.x = element_blank(),
	      axis.title.y = element_blank()) + 
        guides(fill="none", size=guide_legend()) +
	scale_y_continuous(limits=c(0, 9),
	                   breaks=seq(1, 8),
			   labels=rev(c("A", "B", "C", "D", "E", "F", "G", "H")))
  
}

YM = 14
DIAMOND_WIDTH=10
DIAMOND_HEIGHT=3

diamond(filter(x.sum, Language == "English", Task == "Open"), ym=YM) 
ggsave("output/english_diamond.pdf", width=DIAMOND_WIDTH, height=DIAMOND_HEIGHT)

diamond(filter(x.sum, Language == "Spanish", Task == "Open"), ym=YM) 
ggsave("output/spanish_diamond.pdf", width=DIAMOND_WIDTH, height=DIAMOND_HEIGHT)

diamond(filter(x.sum, Language == "Tsimane", Task == "Open"), ym=YM)
ggsave("output/tsimane_diamond.pdf", width=DIAMOND_WIDTH, height=DIAMOND_HEIGHT)


diamond(filter(x.sum, Language == "English", Task == "Fixed"), ym=YM) 
ggsave("output/english_fixed_diamond.pdf", width=DIAMOND_WIDTH, height=DIAMOND_HEIGHT)

diamond(filter(x.sum, Language == "Spanish", Task == "Fixed"), ym=YM) 
ggsave("output/spanish_fixed_diamond.pdf", width=DIAMOND_WIDTH, height=DIAMOND_HEIGHT)

diamond(filter(x.sum, Language == "Tsimane", Task == "Fixed"), ym=YM)
ggsave("output/tsimane_fixed_diamond.pdf", width=DIAMOND_WIDTH, height=DIAMOND_HEIGHT)


# Individual subject surprisal analysis ------------------------------

# mostly copied in from color_analysis_kyle.Rmd

e = ColorData %>%
  separate(Experiment, into=c("Language", "Task"), sep="_") %>%
  group_by(Language, Task, grid_location) %>%
    mutate(location.n = n()) %>%
    ungroup() %>%
  group_by(Language, Task, color) %>%
    mutate(color.n = n()) %>%
    group_by(Language, Task, grid_location, color) %>% 
  mutate(n=n()) %>%
  ungroup()
  
e.ent = e %>%
  group_by(Language, Task, grid_location, color, subject) %>% 
  summarise(p.loc.given.color = unique(n/color.n), 
            p.color.given.loc = unique(n/location.n),
            inside = #p.color.given.loc * 
            log2(1/p.loc.given.color))
  
  #sum up over all words
  e.ent.loc = e.ent #summarise(group_by(ungroup(e.ent), Language, Task, subject, location), s.loc = first(inside))
  
  #average over all chips
  ent = summarise(group_by(ungroup(e.ent.loc), Language, Task, subject), e=round(mean(inside), 2))
  tsimane = filter(ent, Language == "Tsimane")

  # We lose 4 subjects in this join
  tsimane = tsimane %>%
    rename(Subject=subject) %>%
    mutate(Subject=as.character(Subject)) %>%
    inner_join(ts_demo)
    
  ggplot(tsimane, aes(y=e, x=Spanish)) +
    geom_point() +
    geom_smooth(method='lm') +
    theme_bw() +
    xlab("Spanish knowledge") +
    ylab("Uncertainty") +
    xlim(0, 11)
  ggsave("output/education_correlation.pdf", width=8 ,height=4)
  
  print(with(filter(tsimane, Task == "Open"), cor.test(e, Spanish, method='pearson')))
  print(with(filter(tsimane, Task == "Fixed"), cor.test(e, Spanish, method='pearson')))
  print(with(tsimane, cor.test(e, Spanish, method='pearson')))

  m = lm(e ~ Age + Spanish, data=tsimane)
  print(summary(m))


# Contrastive labeling -------------------------------------------------

sed = sedivy
sed$UsesColor = sed$color != "N" 
sed$UsesAdj = sed$other_adj != "N"
sed$same_noun = paste("same_noun: ", sed$same_noun)
#sed = filter(sed, same_noun == "Y")
sed.sum = group_by(sed, Language, Presentation_type, subject, same_noun) %>%
  summarise(m.color = mean(UsesColor), m.adj = mean(UsesAdj))

give.n <- function(x){
  return(c(y = -.1, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

sed.sum = group_by(sed.sum, Language, Presentation_type, subject, same_noun) %>% 
  mutate(uses.something = (m.color > 0 | m.adj > 0))

ggplot(sed.sum, aes(x=Language, y=m.color)) +
  geom_jitter(alpha=1, position = position_jitter(width = .04, height=.04)) +  
  facet_grid(same_noun ~ Presentation_type) +
  theme_bw() +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.01, colour="red") +
  stat_summary(fun.y = mean, geom = "point", colour="red") +
  ylab("proportion using color word") + 
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
  ggtitle("COLOR")
#ggsave("sedivy_colors.pdf")

ggplot(filter(sed.sum, uses.something == TRUE), aes(x=Language, y=m.color)) +
  geom_jitter(alpha=1, position = position_jitter(width = .04, height=.04)) +  
  facet_grid(same_noun ~ Presentation_type) +
  theme_bw() +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.01, colour="red") +
  stat_summary(fun.y = mean, geom = "point", colour="red") +
  ylab("proportion using color word") + 
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
  ggtitle("COLOR")
#ggsave("sedivy_colors_uses_something.pdf")

ggplot(sed.sum, aes(x=Language, y=m.adj)) +
  geom_jitter(alpha=1, position = position_jitter(width = .04, height=.04)) +  
  facet_grid(same_noun ~ Presentation_type) +
  theme_bw() +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.01, colour="red") +
  stat_summary(fun.y = mean, geom = "point", colour="red") +
  ylab("proportion using non-color adj. word") + 
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
  ggtitle("OTHER ADJ.")
#ggsave("sedivy_adj.pdf")

###make all same
sed.sum = group_by(sed, Language, Presentation_type, subject, same_noun) %>%
  summarise(m.color = mean(UsesColor), m.adj = mean(UsesAdj))

sed.sum2 = filter(ungroup(sed.sum), same_noun == "same_noun:  Y") %>%
  select(Language, subject, Presentation_type, m.color, m.adj) %>%
  gather(variable, value, -Language, -subject, -Presentation_type)
sed.sum2$V = ifelse(sed.sum2$variable == "m.adj", "Other Adjectives", "Color Adjectives")

ggplot(filter(sed.sum2, V == "Color Adjectives"), aes(x=Language, y=value)) +
  geom_jitter(alpha=.5, position = position_jitter(width = .04, height=.04)) +  
  facet_grid(Presentation_type ~ V) +
  theme_bw() +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.01, colour="#b20000", size=1,alpha=.7) +
  stat_summary(fun.y = mean, geom = "point", colour="#b20000", size=3, alpha=.7) +
  ylab("proportion using modifier") +
  xlab("") +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) 
#ggsave("sedivy_color_adj.pdf", width=4, height=8)

ggplot(filter(sed.sum2), aes(x=paste(Language, Presentation_type), y=value)) + 
  geom_jitter(alpha=.5, position = position_jitter(width = .04, height=.04)) +  
  facet_grid(. ~ V) +
  theme_bw() +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.01, colour="#b20000", size=1,alpha=.7) +
  stat_summary(fun.y = mean, geom = "point", colour="#b20000", size=3, alpha=.7) +
  ylab("proportion using modifier") +
  xlab("") +
  theme(axis.text.x = element_text(angle=90))
#ggsave("sedivy_color_adj_pres.pdf", width=4, height=8)

sed$object.to.be.labeled = gsub("green ", "", sed$object.to.be.labeled)
sed$object.to.be.labeled = gsub("red ", "", sed$object.to.be.labeled)
sed$object.to.be.labeled = gsub("pink ", "", sed$object.to.be.labeled)
sed$object.to.be.labeled = gsub("yellow ", "", sed$object.to.be.labeled)
sed$object.to.be.labeled = gsub("blue ", "", sed$object.to.be.labeled)

sed.obj = group_by(filter(sed, same_noun != "same_noun:  N"),
                   object.to.be.labeled,
		   Language,
		   Presentation_type,
		   presentation_order) %>%
          summarise(m=mean(color != "N"))
sed.obj$presentation_order = as.factor(sed.obj$presentation_order)
sed.obj$object.to.be.labeled = reorder(sed.obj$object.to.be.labeled, sed.obj$m)
sed.obj$IsNat = ifelse(grepl("pepper|apple|banana|tomato",
                             as.character(sed.obj$object.to.be.labeled)),
		       "Artificial",
		       "Natural")

ggplot(filter(sed.obj, Presentation_type == "1 at a time"),
       aes(x=object.to.be.labeled,
           y=m * 100,
	   group=paste(Language, object.to.be.labeled),
	   color=Language,
	   fill=paste(Language, presentation_order))) +
	geom_line(size=2) + 
	geom_point(size=5, shape=21) +  
	theme(axis.text.x = element_text(angle=90, hjust=1, color="black"),
	      axis.text.y = element_text(color="black")) +
        ylab("proportion using color word") +
	xlab("") +
	theme(legend.position="none") +
	scale_color_manual(values=c("gray", "black")) + 
	theme(panel.grid.major = element_blank(),
	      panel.grid.minor = element_blank(), 
 	      panel.background = element_blank(),
	      axis.line = element_line(colour = "black"),
	      axis.text.x = element_text(size=18, hjust=1, vjust=.5),
	      axis.text.y=element_text(size=18),
	      axis.title.y=element_text(size=18)) +
	scale_fill_manual(values=c("white", "lightgray", "white", "black")) +
  	scale_linetype_manual(values=c(2, 1)) +
	ylab("Use of\n color word (%)") +
	scale_y_continuous(breaks=c(0, 100)) + #seq(0, 100, 70)) +
	ylim(0, 100) + 
	facet_grid(. ~ IsNat, scales="free_x") +
	theme(strip.background = element_blank(),
	      strip.text.x = element_blank())
ggsave("output/contrastive_color_word_use.pdf", width=432/100, height=511/100)

##look at nat vs not
sed$IsNat = ifelse(grepl("pepper|apple|banana|tomato", as.character(sed$object.to.be.labeled)),
                   "Natural",
		   "Artificial")
sed.sum = group_by(sed, Language, Presentation_type, subject, same_noun, IsNat) %>%
  summarise(m.color = mean(UsesColor), m.adj = mean(UsesAdj))
sed.sum2 = filter(ungroup(sed.sum), same_noun == "same_noun:  Y", Presentation_type == "1 at a time") %>%
  select(Language, subject, m.color, m.adj, IsNat) %>%
  gather(variable, value, -Language, -subject, -IsNat)
sed.sum2$V = ifelse(sed.sum2$variable == "m.adj", "Other Adjectives", "Color Adjectives")

ggplot(sed.sum2, aes(x=Language, y=value)) +
  geom_jitter(alpha=.5, position = position_jitter(width = .04, height=.04)) +  
  facet_grid(IsNat ~ V) +
  theme_bw() +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.01, colour="#b20000", size=1,alpha=.7) +
  stat_summary(fun.y = mean, geom = "point", colour="#b20000", size=3, alpha=.7) +
  ylab("proportion using modifier") +
  xlab("") +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) 
#ggsave("natural_artificial.pdf", width=4, height=5)



sed$subject = with(sed, as.factor(paste(subject, Language)))
s = filter(sed, same_noun == "same_noun:  Y", Presentation_type == "1 at a time")
l.color = glmer(data=s, family = 'binomial',
            UsesColor ~ Language + (1|subject) + (1 + Language|object.to.be.labeled))
print("Color adjective use, 1 at a time")	    
print(summary(l.color))


s = filter(sed, same_noun == "same_noun:  Y", Presentation_type == "2 at a time" | Language == "English")
l.color = glmer(data=s, family = 'binomial',
            UsesColor ~ Language + (1|subject) + (1 + Language|object.to.be.labeled))
print("Color adjective use, 2 at a time")	    
print(summary(l.color))


s2 = filter(sed, same_noun == "same_noun:  Y", Language == "Tsimane")
l.color = glmer(data=s2, family = 'binomial',
            UsesColor ~  Presentation_type + (1|subject) + (1 + Presentation_type|object.to.be.labeled))
summary(l.color)

l_no_int = glmer(data=s2, family="binomial",
  UsesColor ~ presentation_order + IsNat + (1|subject))

print(summary(l_no_int))

l_int = glmer(data=s2, family="binomial",
  UsesColor ~ presentation_order * IsNat + (1|subject))

print(anova(l_no_int, l_int))

s2_en = filter(sed, same_noun == "same_noun:  Y", Language == "English")
l_no_int_en = glmer(data=s2_en, family="binomial",
  UsesColor ~ presentation_order + IsNat + (1|subject))
print(summary(l_no_int_en))




s = group_by(s, subject, Language) %>% mutate(uses.something = mean(UsesColor) > 0)
l.color.uses = glmer(data=filter(s, uses.something == T), family = 'binomial',
            UsesColor ~ Language + (1|subject) + (1 + Language|object.to.be.labeled))
print("Color adjective use, subjects who use at least 1 color adjective")	    
print(summary(l.color.uses))


l.adj = glmer(data=s, family = 'binomial',
              UsesAdj ~ Language + (1|subject) + (1 + Language|object.to.be.labeled))
summary(l.adj)

l.adj.uses = glmer(data=filter(s, uses.something == T), family = 'binomial',
            UsesAdj ~ Language + (1|subject) + (1 + Language|object.to.be.labeled))
summary(l.adj.uses)

s.melt = select(s, subject, Language, object.to.be.labeled, UsesColor, UsesAdj) %>%
  gather(variable, value, -subject, -Language, -object.to.be.labeled)
  
l.all = glmer(data=s.melt, family = 'binomial',
              value ~ variable * Language +
	     (1 + variable|subject) +
	     (1 + variable + Language|object.to.be.labeled))
summary(l.all)

l.class.order.ts = glmer(data=filter(sed, Language == "Tsimane"), family='binomial',
                      UsesColor ~ presentation_order * IsNat
		      + (1|subject))
print(summary(l.class.order.ts))

                  

# RT analysis ---------------------------------------------------

rts = mutate(rts, rt=(Isabel.Time + Rashida.Time)/2)

#object entropy
rtent = group_by(rts, Language, color_object.label) %>%
  mutate(label.count = n()) %>%
  group_by(Language, color_object.label, response) %>%
  summarise(response.count = n(), total.count = first(label.count), p=response.count/total.count) %>%
  ungroup() %>% group_by(Language, color_object.label) %>% summarise(ent=sum(-p*log2(p)))

rtent$IsColor = ifelse(grepl("^[A-H]", rtent$color_object.label), "color naming", "object naming")
rtent = merge(rtent, rts[, c("Language", "color_object.label", "rt")], by = c("Language", "color_object.label"))
rtent$log.rt = log(rtent$rt)
#rtent.sum$log.rt = NA

std.error = function(x) sqrt(var(x)/length(x))


rtent.sum = group_by(rtent, Language, IsColor, color_object.label) %>%
  summarise(m.ent=mean(ent),
            m.rt = mean(log.rt),
	    se = std.error(log.rt),
	    l = m.rt - 1.96 * se,
	    u = m.rt + 1.96 * se)
	    
re = arrange(rtent.sum, m.ent) %>% mutate(l=lag(m.ent), too.close= abs(l - m.ent) < .06)

re$jitter = rnorm(nrow(re), 0, .05)
rtent = merge(rtent, select(re, Language, IsColor, color_object.label, jitter),
              by = c("Language", "IsColor", "color_object.label"))
rtent$ent = rtent$ent + rtent$jitter
rtent.sum = group_by(rtent, Language, IsColor, color_object.label) %>%
  summarise(m.ent=mean(ent),
            m.rt = mean(log.rt),
	    se = std.error(log.rt),
	    l = m.rt - 1.96 * se,
	    u = m.rt + 1.96*se)

ggplot(data=rtent, aes(x=ent,y=log.rt,label=color_object.label)) + geom_point(alpha=1, shape=95, size=2) + 
   facet_grid(IsColor ~ Language, drop = T, scales="free") + theme_bw(18) +
   geom_smooth(method=lm) +
   geom_point(data=rtent.sum, aes(x=m.ent, y=m.rt), size=2, colour="red") +
   geom_errorbar(data=rtent.sum, aes(x=m.ent, ymax=u, ymin=l, y=0), colour='red')

#rtent.sum$m.ent = rtent.sum$m.ent + rnorm(nrow(rtent.sum), 0, .3)
ggplot(data=rtent, aes(x=ent,y=log.rt,label=color_object.label, group=Language, colour=Language)) +
  #geom_point(alpha=.4, shape=95, size=3) + 
  facet_grid(. ~ IsColor, drop = T, scales="free_x") +
  theme_bw(24) +
  geom_smooth(method=lm) +
  geom_point(data=rtent.sum, aes(x=m.ent, y=m.rt, colour=Language), size=3) +
  geom_errorbar(data=rtent.sum, aes(x=m.ent, ymax=u, ymin=l, y=0, colour=Language)) + 
  scale_colour_manual(values=c("red", "blue")) +
  xlab("entropy") +
  ylab("log reaction time") +
  theme(legend.position = "bottom")
#ggsave("object_rts.pdf", width=14, height=12)

rtent.sum$Language = factor(rtent.sum$Language, levels=c("Tsimane", "English"))
ggplot(data=filter(rtent.sum, IsColor == "color naming"),
  aes(x=m.ent,
      y=m.rt,
      label=color_object.label,
      group=Language,
      colour=Language)) +
  #geom_point(alpha=.4, shape=95, size=3) + 
  geom_smooth(method=lm, se=F, size=2, alpha=.8) +
  geom_errorbar(data=filter(rtent.sum, IsColor == "color naming"),
                aes(x=m.ent, ymax=u, ymin=l, y=0, colour=Language)) +
  geom_point(size=5) +
  xlab("Entropy\nColor naming") +
  ylab("Reaction time (log sec.)") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
	axis.line = element_line(colour = "black"),
	axis.text.x = element_text(size=18),
	axis.text.y=element_text(size=18),
	axis.title.y=element_text(size=18),
	axis.title.x=element_text(size=18)) + 
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4),
                     limits=c(0, 4.5)) +
  scale_y_continuous(breaks=c(0, 2), limits=c(-.3, 2)) +
  scale_color_manual(values=c("black", "lightgray"))
ggsave("output/color_naming_rts.pdf", width=432/100, height=511/100 * 1.1)

rtent.sum$Language = factor(rtent.sum$Language, levels=c("Tsimane", "English"))
ggplot(data=filter(rtent.sum, IsColor != "color naming"),
       aes(x=m.ent,y=m.rt,label=color_object.label, group=Language, colour=Language)) +
  #geom_point(alpha=.4, shape=95, size=3) + 
  geom_smooth(method=lm, se=F, size=2, alpha=.8) +
  geom_errorbar(data=filter(rtent.sum, IsColor != "color naming"),
                aes(x=m.ent, ymax=u, ymin=l, y=0, colour=Language)) +
  geom_point(size=5) +
  xlab("Entropy\nObject naming") +
  ylab("Reaction time (log sec.)") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
	panel.background = element_blank(),
	axis.line = element_line(colour = "black"),
	axis.text.x = element_text(size=18),
	axis.text.y=element_text(size=18),
	axis.title.y=element_text(size=18),
	axis.title.x=element_text(size=18)) + 
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4), limits=c(0, 4.5)) +
  scale_y_continuous(breaks=c(0, 2), limits=c(-.3, 2)) +
  scale_color_manual(values=c("black", "lightgray"))
ggsave("output/object_naming_rts.pdf", width=432/100, height=511/100 * 1.1)


rpp = unique(rtent[, c("Language", "color_object.label", "ent", "IsColor")])
rts.l = merge(rts, rpp, by = c("Language", "color_object.label"))

rts.l$Language = as.factor(rts.l$Language)
contrasts(rts.l$Language) = contr.sum(2)


l.rt.color = lmer(log(rt) ~ scale(ent)  * Language +
                  (1|subject) +
		  (1 + Language|color_object.label),
		  data=subset(rts.l, IsColor == "color naming"))
print(summary(l.rt.color))
l.rt.color.0 = lmer(log(rt) ~ scale(ent)  + Language +
                    (1|subject) +
		    (1 + Language|color_object.label),
		    data=subset(rts.l, IsColor == "color naming"))
print(anova(l.rt.color, l.rt.color.0))

l.rt.object = lmer(log(rt) ~ scale(ent)  * Language +
                   (1 + ent|subject) +
		   (1 + Language|color_object.label),
		   data=filter(rts.l, IsColor != "color naming"))
print(summary(l.rt.object))
l.rt.object.0 = lmer(log(rt) ~ scale(ent)  + Language +
                     (1 + ent|subject) +
		     (1 + Language|color_object.label),
		     data=filter(rts.l, IsColor != "color naming"))
l.rt.object.0.0 = lmer(log(rt) ~ scale(ent) + scale(ent):Language +
                       (1 + ent|subject) +
		       (1 + Language|color_object.label),
		       data=filter(rts.l, IsColor != "color naming"))
print(anova(l.rt.object, l.rt.object.0.0))
print(anova(l.rt.object, l.rt.object.0))


# Permutation test ----------------------------------------------------------

# Supporting functions for tree navigation -----------------------------
GetRoot<-function(x,Tree){
  if (Tree[x] == x){
    return(x)
  }
  else{
    return(GetRoot(Tree[x],Tree))
  }
}
Merge<-function(x,y,Tree){
  r1 <- GetRoot(x,Tree)
  r2 <- GetRoot(y,Tree)
  # Connect the trees
  Tree[r1]=r2
  return(Tree)
}
CountPools<-function(Tree){
  temp <- data.frame(Tree) %>%
    add_rownames("Id") %>%
    mutate(Root = (Id==Tree)) %>%
    dplyr::select(Root)
  return(sum(temp$Root))
}
# Create a function that tells you neighbors of the same color given an Id
GetNeighbors<-function(data,x){
  neighbors<-c()
  # Find neighbors of entry with id x in dataset
  central <- filter(data,Enumeration==x)
  crow <- central$row
  ccolumn <- central$column
  ccolor <- central$color
  # Data to the left
  ldat <- filter(data, row==crow, column<ccolumn) %>% arrange(-column)
  if (dim(ldat)[1]!=0){
    if (ldat$color[1] == ccolor){
      neighbors <- c(neighbors,ldat$Enumeration[1])
    }
  }
  # Data to the right
  rdat <- filter(data, row==crow, column>ccolumn) %>% arrange(column)
  if (dim(rdat)[1]!=0){
    if (rdat$color[1] == ccolor){
      neighbors <- c(neighbors,rdat$Enumeration[1])
    }
  }
  # Data to the top
  tdat <- filter(data, row<crow, column==ccolumn) %>% arrange(-row)
  if (dim(tdat)[1]!=0){
    if (tdat$color[1] == ccolor){
      neighbors <- c(neighbors,tdat$Enumeration[1])
    }
  }
  # Data below
  bdat <- filter(data, row>crow, column==ccolumn) %>% arrange(row)
  if (dim(bdat)[1]!=0){
    if (bdat$color[1] == ccolor){
      neighbors <- c(neighbors,bdat$Enumeration[1])
    }
  }
  # Check immediate neighbors
  t<-filter(data, row==(crow-1),column==(ccolumn-1))
  if (dim(t)[1]!=0){
    if (t$color[1] == ccolor){
      neighbors <- c(neighbors,t$Enumeration[1])
    }
  }
  t<-filter(data, row==(crow-1),column==(ccolumn+1))
  if (dim(t)[1]!=0){
    if (t$color[1] == ccolor){
      neighbors <- c(neighbors,t$Enumeration[1])
    }
  }
  t<-filter(data, row==(crow+1),column==(ccolumn-1))
  if (dim(t)[1]!=0){
    if (t$color[1] == ccolor){
      neighbors <- c(neighbors,t$Enumeration[1])
    }
  }
  t<-filter(data, row==(crow+1),column==(ccolumn+1))
  if (dim(t)[1]!=0){
    if (t$color[1] == ccolor){
      neighbors <- c(neighbors,t$Enumeration[1])
    }
  }
  return(neighbors)
}

# Count pools per participant ------------------------------------------

# CurrDat = filter(ColorData,Id==part)

GetPoolCount<-function(CurrDat){
  # Get the number of clusters in a participant's data
  # Enumerate the grid.
  CurrDat<-CurrDat %>% add_rownames("Enumeration")
  CurrDat$Enumeration <- as.numeric(CurrDat$Enumeration)
  # Create a tree list.
  Connections <- seq(1,max(CurrDat$Enumeration))
  # For each entry in the data, get the neighbors of the same color and call the merge function!
  for (item in CurrDat$Enumeration){
    neighbors <- GetNeighbors(CurrDat,item)
    for (neighbor in neighbors){
      Connections <- Merge(item,neighbor,Connections)
    }
  }
  return(CountPools(Connections))
}

ColorData = ColorData %>%
  mutate(rc=CodeToRowColumn(grid_location)) %>%
  separate(rc, into=c("row", "column"), sep="_") %>%
  mutate(row=-as.numeric(row), column=as.numeric(column)) %>%
  filter(!(Experiment == "Tsimane_Fixed" & color == 19)) # get rid of noncompliant subject

ColorData$Id = paste(ColorData$Experiment,ColorData$subject)

# Set this to something higher than 0 to do the permutation test
# The permutation test is pretty slow so it's disabled by default
samples=0
participant<-c()
cluster<-c()

for (discard in 1:samples) {
  for (part in unique(ColorData$Id)){
    print(c(discard,part))
    # Create a temporary data frame
    temp = filter(ColorData,Id==part)
    # Shuffle
    temp$color=sample(temp$color)
    ClusterCount<-GetPoolCount(temp)
    participant<-c(participant,part)
    cluster<-c(cluster,ClusterCount)
  }
}

# Save your results!
permutation_test_results = tbl_df(data.frame(participant,cluster))


# Work with objects --------------------------------------------------------

english_objects = read.csv("english_objects.csv", stringsAsFactors=F) %>%
  select(-X) %>%
  filter(!str_detect(thing, "focal"))

spanish_objects = read.csv("spanish_objects.csv", stringsAsFactors=F) %>%
  select(-X) %>%
  filter(!str_detect(thing, "focal"))

tsimane_objects = read.csv("tsimane_objects.csv", stringsAsFactors=F) %>%
  gather(thing, color, -subject, -Year, -Language) %>%
  select(-Year)

names(tsimane_objects) = c("subject", "Language", "thing", "color")

objects = rbind(tsimane_objects, english_objects, spanish_objects)

# Clean up names
objects$thing = gsub("grass.*", "grass", objects$thing)
objects$thing = gsub("yuca.*ins.*", "yuca inside", objects$thing)
objects$thing = gsub("yuca.*outs.*", "yuca outside", objects$thing)
objects$thing = gsub("e\\.b", "e b", objects$thing)

obj = group_by(objects, Language, thing) %>%
  mutate(total.labels = n()) %>%
  filter(color %in% as.character(1:11))


get_simple_ent = function(obj) {

#obj$color = ifelse(obj$color %in% as.character(1:11), obj$color, "other")

obj.sum = group_by(obj, thing, color, Language) %>%
  summarise(n= n(), total.labels=first(total.labels)) %>%
  ungroup()

  e.simple.ent = obj.sum
  e.simple.ent$p = with(e.simple.ent, n/total.labels)
  
  e.simple.ent = group_by(e.simple.ent, Language, thing) %>%
     summarise(s = sum(-p*log2(p))) %>%
     ungroup() 
  ent = ungroup(e.simple.ent) %>% group_by(Language) %>%
     summarise(m=mean(s))
  tsi = filter(e.simple.ent, Language == "Tsimane") %>%
     arrange(s) %>%
     mutate(r=rank(s))
  e.simple.ent = merge(e.simple.ent, tsi[, c("thing", "r")], by = c("thing"))

  e.simple.ent$thing = gsub(" ", "\n", e.simple.ent$thing)
    e.simple.ent$thing = as.factor(e.simple.ent$thing)
  e.simple.ent$thing = reorder(e.simple.ent$thing, e.simple.ent$r)
  e.simple.ent$Language = gsub("Tsimane", "Tsimane'", e.simple.ent$Language)
  e.simple.ent$Language = factor(e.simple.ent$Language, c("English", "Spanish", "Tsimane'"))
  e.simple.ent
}

  e.simple.ent = get_simple_ent(obj)

  num_object_bootstrap_samples = 100
  
  e.simple.ent.resampled = NA
  for (i in 1:num_object_bootstrap_samples) {
      new_data = obj %>% group_by(Language, thing) %>% sample_frac(replace=T) %>% ungroup()
      new_entropies = get_simple_ent(new_data) %>% mutate(sample=i)
      if (i == 1) {
          e.simple.ent.resampled = new_entropies
      } else {
          e.simple.ent.resampled = rbind(e.simple.ent.resampled, new_entropies)
      }
   }

   e.simple.ent.resampled %>%
     group_by(Language, thing) %>%
     summarise(upper=quantile(s, .95), lower=quantile(s, .05)) %>%
     inner_join(e.simple.ent) -> e.simple.ent
    
  filter(e.simple.ent, !str_detect(thing, "yuca")) %>%
  ggplot(aes(x=thing, y=s + .05, group=Language, colour=Language, fill=Language)) +
    geom_line(size=2) +
    ylab("entropy on object color") +
    xlab("object") +
    geom_point(size=6) + 
    #geom_bar(stat="identity", position="dodge") + 
    theme_bw(20) +
    theme(legend.position="bottom")  +
    scale_y_continuous(breaks=c(0, 1.5, 3), limits=c(0, 3.2)) +
    scale_color_manual(values=c("lightgray", "darkgray", "black")) +
    ylab("Entropy on Object Color") +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(), 
	  panel.grid.minor = element_blank(),
	  axis.line = element_line(colour = "black")) +
    theme(legend.position="none",
          axis.text.x = element_text(color="black", size=14, angle=90, hjust=1, vjust=.5),
	  axis.text.y = element_text(color="black", size=18),
	  axis.title.x = element_blank(),
	  axis.title.y=element_text(size=20)) +
    theme(axis.ticks = element_line(size = 1)) 
  ggsave("output/objects_entropy.pdf", width=(710/100) * 1.3, height=497/100)

object_entropy = group_by(e.simple.ent, Language) %>% summarise(m=mean(s))
write.csv(object_entropy, "output/object_entropy.csv")

## Calculate type-token ratios on WCS ----------------------------

Y = 20

get_ttr = function(d) {
ttrs = c()
hapaxes = c()

NUM_SAMPLES = 1000

for(i in 1:NUM_SAMPLES) {
  d %>%
    select(subject) %>%
    unique() %>%
    mutate(ok=as.numeric(as.factor(subject)) %in% sample(1:n(), Y)) %>%
    inner_join(d) %>%
    filter(ok) %>%
    group_by(color) %>%
      mutate(num_color_tokens=n()) %>%
      ungroup() -> d2
  d2 %>%
    summarise(num_tokens=length(color),
              num_types=length(unique(color)),
	      ttr=num_types/num_tokens,
	      hapax_rate=sum(num_color_tokens==1)/num_tokens) -> d_ttr
  ttrs = c(ttrs, d_ttr$ttr)
  hapaxes = c(hapaxes, d_ttr$hapax_rate)
}

d_ttr = data.frame(ttr=ttrs, hapax_rate=hapaxes, baseline=T)
d_ttr
}

max_sample = function(f, xs, num_samples, num_attempts) {
  samples = c()
  scores = c()
  for(i in 1:num_attempts) {
     the_sample = sample(xs, num_samples)
     samples = c(samples, the_sample)
     scores = c(scores, f(the_sample))
  }
  start = (which.max(scores) - 1) * num_samples + 1
  end = start + num_samples - 1
  samples[start:end]
}

assert(
  mean(replicate(100, sum(sample(1:50, 10)))) <
  mean(replicate(100, sum(max_sample(sum, 1:50, 10, 10))))
)

uniformity_criterion = function(xs) {
  mean(gap_size(xs)^2)		    
}

gap_size = function(numbers) {
  xs = numbers[order(numbers)]
  n = length(xs)
  xs[2:n] - xs[1:(n-1)]
}

mysample = function(xs, num_samples) {
  max_sample(uniformity_criterion, xs, num_samples, 10)
}

ColorDataFiltered = ColorData %>%
  group_by(Experiment) %>%
    select(subject) %>%
    unique() %>%
    mutate(ok=as.numeric(as.factor(subject)) %in% mysample(1:n(), Y)) %>%
    ungroup() %>%
  inner_join(ColorData) %>%
  filter(ok)

assert(ColorDataFiltered %>%
  group_by(Experiment) %>%
    summarise(s=length(unique(subject))) %>%
    select(s) %>%
    unique() == Y)

ts_ttr = ColorData %>% filter(Experiment == "Tsimane_Open") %>% get_ttr()
ts_ttr_fixed = ColorData %>%
  filter(Experiment == "Tsimane_Fixed", color != 19) %>% # weird erroneous one-off
  get_ttr()

answers %>%
  group_by(Lang_Name) %>%
    mutate(s=length(unique(Subject)), ok=s >= Y) %>%
    ungroup() %>%
  filter(ok) -> answers_filtered

answers_filtered %>%
  select(Lang_Name, Subject, s) %>%
  unique() %>%
  group_by(Lang_Name) %>%
    mutate(ok=Subject %in% sample(1:s, Y)) %>%
    ungroup() %>%
  filter(ok) %>%
  inner_join(answers) -> answers_Y

assert(answers_Y %>%
  group_by(Lang_Name) %>%
    summarise(s=length(unique(Subject))) %>%
    select(s) %>%
    unique() == Y)

WCS_ttr = answers_Y %>%
  group_by(Lang_Name, Term) %>%
    mutate(num_color_tokens=n()) %>%
    ungroup() %>%
  group_by(Lang_Name) %>%
    summarise(num_tokens=length(Term),
              num_types=length(unique(Term)),
	      ttr=num_types/num_tokens,
	      hapax_rate=sum(num_color_tokens==1)/num_tokens) %>%
    ungroup()

WCS_ttr %>%
  ggplot(aes(x=ttr)) +
    geom_histogram(bins=50) +
    xlab("Type-Token Ratio") +
    ylab("Languages") +
    geom_vline(color="red", aes(xintercept=mean(ts_ttr$ttr))) +
    geom_vline(color="blue", aes(xintercept=mean(ts_ttr_fixed$ttr)))

ggsave("output/wcs_ttr.pdf", width=5.6, height=3.8)

WCS_ttr %>%
  ggplot(aes(x=hapax_rate)) +
    geom_histogram(bins=50) +
    xlab("One-off rate") +
    ylab("Languages") +
    geom_vline(color="red", aes(xintercept=mean(ts_ttr$hapax_rate))) +
    geom_vline(color="blue", aes(xintercept=mean(ts_ttr_fixed$hapax_rate)))    

ggsave("output/wcs_hapax.pdf", width=5.6, height=3.8)

ts_ttr = rbind(ts_ttr, select(WCS_ttr, ttr, hapax_rate) %>% mutate(baseline=F))

# TODO analyze with mean number of terms per subject / total terms with freq 2

# Second round of lexical diversity stuff ------------------------

# 1. Lexical-overlap score:  calculate a number for each subject as follows: for the words spoken by this subject, the % of terms that were given at least once by 1/2 of the subjects (or maybe 3/4? some high threshold) average this number across subjects to get a lexical-overlap score for a language.  
# for the fixed task, we should get close to 1 here: subjects are choosing words that everyone else is also choosing. for the open task, we will get a much lower number: subjects are choosing words that few others are choosing.

# 2. Lexical-one-off (hapax) score:  calculate a number for each subject as follows: for the words spoken by this subject, the % of terms that were given by no other subjects (or maybe one other) average this number across subjects to get a one-off / hapax score for a language.  

PROP_THRESHOLD = 3/4

WCS_with_subjects_restr = answers_Y %>%
    rename(Experiment=Lang_Name,
	   color=Term,
	   grid_location=FielConcat,
	   subject=Subject) %>%
    select(Experiment, color, grid_location, subject) %>%	   
    mutate(WCS=T) %>%
    rbind(
      ColorDataFiltered %>%
	  mutate(WCS=F) %>%
	  select(Experiment, color, grid_location, subject, WCS))

WCS_with_subjects_full = answers %>%
    rename(Experiment=Lang_Name,
	   color=Term,
	   grid_location=FielConcat,
	   subject=Subject) %>%
    select(Experiment, color, grid_location, subject) %>%	   
    mutate(WCS=T) %>%
    rbind(
      ColorData %>%
	  mutate(WCS=F) %>%
	  select(Experiment, color, grid_location, subject, WCS))

NUM_TERMS_CUTOFF = 5 # need more than this many instances of a term for it to count

make_diversity = function(WCS_with_subjects) {
WCS_with_subjects %>%
    group_by(Experiment, color) %>%
        mutate(num_subjects_using=length(unique(subject))) %>%
    	ungroup() %>%
    group_by(Experiment, grid_location) %>%
        mutate(modal_term=Mode(color)) %>%
        ungroup() %>%
    group_by(Experiment) %>%
        mutate(num_modal_terms=length(unique(modal_term))) %>%
        ungroup() %>%	
    group_by(Experiment) %>%
        mutate(full_num_terms=sum(tabulate(color) > 0),
	       num_subjects=length(unique(subject)),
	       num_terms=sum(tabulate(color) > NUM_TERMS_CUTOFF)) %>%
	ungroup() %>%
    mutate(subj_hapax=num_subjects_using == 1,
           prop_subjects_using=num_subjects_using/num_subjects,
           over_threshold=prop_subjects_using > PROP_THRESHOLD) %>%
    group_by(Experiment, subject) %>%
        summarise(num_modal_terms=first(num_modal_terms),
	          num_terms=first(num_terms),
		  full_num_terms=first(full_num_terms),
		  prop_over_threshold = mean(over_threshold),
               	  prop_subj_hapax = mean(subj_hapax)) %>%
        ungroup() %>%
    group_by(Experiment) %>%
        summarise(num_modal_terms=first(num_modal_terms),
	          num_terms=first(num_terms),
	          full_num_terms=first(full_num_terms),
		  mean_prop_over_threshold=mean(prop_over_threshold),
	          mean_prop_subj_hapax=mean(prop_subj_hapax)) %>%
        ungroup() %>%
    mutate(WCS=!(Experiment %in% c("Tsimane_Open",
    			    	   "Tsimane_Fixed",
				   "English_Open",
				   "English_Fixed",
				   "Spanish_Open",
				   "Spanish_Fixed")))
				   }

WCS_diversity = make_diversity(WCS_with_subjects_restr)
WCS_diversity_full = make_diversity(WCS_with_subjects_full)


ts_diversity_open = filter(WCS_diversity, Experiment == "Tsimane_Open")
ts_diversity_fixed = filter(WCS_diversity, Experiment == "Tsimane_Fixed")

WCS_diversity[WCS_diversity$Experiment == "Tsimane_Fixed",]$num_terms = 8

gg_color_hues <- function(n) {
   hues = seq(15, 375, length=n+1)
   hcl(h=hues, l=65, c=100)[1:n]
}

hues = gg_color_hues(6)
hues = c("red", "darkred", "green", "darkgreen", "cyan", "darkcyan", "darkgrey")

WCS_diversity %>%
    mutate(WCS=ifelse(WCS, "WCS", as.character(Experiment))) %>%
    ggplot(aes(x=mean_prop_over_threshold, fill=WCS)) +
        geom_histogram() +
	xlab("Mean color-word-overlap proportion") +
    	ylab("Languages") +
	theme(legend.title=element_blank()) +
	scale_fill_manual(values=hues,
	                  labels=c("English fixed choice",
	                           "English free choice",
				   "Spanish fixed choice",
				   "Spanish free choice",
				   "Tsimane' fixed choice",
				   "Tsimane' free choice",
				   "WCS"))                            				      

ggsave("output/wcs_propcommon.pdf", width=6, height=4)

WCS_diversity %>%
    filter(WCS) %>%
    ggplot(aes(x=mean_prop_subj_hapax)) +
        geom_histogram() +
        xlab("Mean proportion of one-off color words by subject") +
        ylab("Languages") +
        geom_vline(color="red",
                   aes(xintercept=mean(ts_diversity_open$mean_prop_subj_hapax))) +
        geom_vline(color="blue",
                   aes(xintercept=mean(ts_diversity_fixed$mean_prop_subj_hapax)))

ggsave("output/wcs_proponeoff.pdf", width=5.6, height=3.8)

WCS_diversity %>%
    filter(WCS) %>%
    ggplot(aes(x=num_terms)) +
        geom_histogram() +
        xlab("Total number of terms used") +
        ylab("Languages") +
        geom_vline(color="red",
                   aes(xintercept=mean(ts_diversity_open$num_terms))) +
        geom_vline(color="blue",
                   aes(xintercept=mean(ts_diversity_fixed$num_terms)))

WCS_diversity %>%
  mutate(WCS=ifelse(WCS, "WCS", as.character(Experiment))) %>%
  group_by(WCS, num_terms) %>%
    summarise(count=n()) %>%
    ungroup() %>%
  ggplot(aes(x=num_terms, y=count, fill=WCS)) +
    geom_bar(stat="identity", position="stack") +
    xlab("Total terms") +
    ylab("Languages") +
    theme(legend.title=element_blank(),
          legend.position="bottom")
    
ggsave("output/wcs_num_terms.pdf", width=5.6, height=3.8)





# Spanish knowledge demographics --------------------------

ColorDataWithDemographics %>%
  select(Subject, Spanish) %>%
  unique() %>%
  summarise(m=mean(Spanish, na.rm=T),
            num_perfect=sum(Spanish == 11, na.rm=T),
	    tot=sum(!is.na(Spanish))) %>%
  write.csv("output/spanish_knowledge_summary.csv")

# Situating Tsimane' in the WCS ---------------------------

d_to_plot %>%
  mutate(WCS=ifelse(Language %in% c("English", "Spanish", "Tsimane"),
                    as.character(Language),
		    "WCS")) %>%
  select(Language, lang_score, WCS) %>%
  unique() %>%
  ggplot(aes(x=lang_score/80, fill=WCS)) +
    geom_histogram() + 
    xlab("Conditional entropy") +
    ylab("Number of languages") +
    theme(legend.title=element_blank())

ggsave("output/wcs_conditional_entropy.pdf", width=5.6, height=3.8)


WCS %>%
  group_by(Experiment, grid_location) %>%
    summarise(modal_term=Mode(color)) %>%
    ungroup() %>%
  group_by(Experiment) %>%
    summarise(num_modal_terms=length(unique(modal_term))) %>%
    ungroup() %>%
  mutate(Experiment=as.character(Experiment)) %>%    
  mutate(WCS="WCS") -> WCS_modal

ColorDataFiltered %>%
  select(Experiment, grid_location, color) %>%
  group_by(Experiment, grid_location) %>%
    summarise(modal_term=Mode(color)) %>%
    ungroup() %>%
  group_by(Experiment) %>%
    summarise(num_modal_terms=length(unique(modal_term))) %>%
    ungroup() %>%
  filter(Experiment %in% c("English_Open",
                           "Spanish_Open",
			   "Tsimane_Open",
			   "English_Fixed",
			   "Spanish_Fixed",
			   "Tsimane_Fixed")) %>%
  mutate(WCS=as.character(Experiment)) -> our_modal

modal = rbind(WCS_modal, our_modal)

max_modal = max(modal$num_modal_terms)

modal %>%
  group_by(WCS, num_modal_terms) %>%
    summarise(count=n()) %>%
    ungroup() %>%
  mutate(num_modal_terms=factor(num_modal_terms)) %>%    
  ggplot(aes(x=num_modal_terms, y=count, fill=WCS)) +
    geom_bar(stat="identity", position="stack") +
    xlab("Number of modal terms") +
    ylab("Number of languages") +
    theme(legend.title=element_blank()) 

ggsave("output/wcs_num_modal.pdf", width=5.6, height=3.8)

modeprop = function(x) {
  ux = unique(x)
  Z = length(x)
  top = max(tabulate(match(x, ux)))
  top / Z
}

# subject 18 in EO is weird

modeprops = ColorDataFiltered %>%
  select(Experiment, color, grid_location) %>%
  group_by(Experiment, grid_location) %>%
    summarise(mode=as.numeric(as.character(Mode(color))),
              modeprop=modeprop(color)) %>%
    ungroup() %>%
    mutate(mode=ifelse(str_detect(Experiment, "Tsimane"),
                       TSIMANE_COLORS[mode],
		       ENGLISH_COLORS[mode]))
		       
modeprops %>% write.csv("output/modeprops.csv")

ColorDataFiltered %>%
  select(Experiment, color, grid_location) %>%
  group_by(Experiment, grid_location) %>%
    summarise(mode=modeprop(color)) %>%
    ungroup() %>%
  group_by(Experiment) %>%
    summarise(m=mean(mode)) %>%
    ungroup() %>%
    print()

ChipEntropies %>%
  unite(Experiment, Language, Type, sep="_") %>%
  select(Experiment, grid_location, Entropy) %>%
  rbind(WCS_ChipEntropies %>% select(-Prior)) %>%
  inner_join(WCS_diversity) %>%
  mutate(WCS=ifelse(WCS, "WCS", as.character(Experiment))) %>%
  group_by(Experiment, WCS) %>%
    summarise(ratio=first(num_terms/(sum(Entropy)/80))) %>%
    ungroup() %>%
  ggplot(aes(x=ratio, fill=WCS)) +
    geom_histogram() +
    xlab("Number of terms / Conditional entropy") +
    ylab("Languages") +
    theme(legend.title=element_blank()) +
    scale_fill_discrete(labels=c("Tsimane' fixed choice", "Tsimane' free choice", "WCS"))

ggsave("output/wcs_term_condent_ratio.pdf", width=6, height=4)

WCS_diversity2 = WCS_diversity %>%
  select(-WCS) %>%
  inner_join(modal)

WCS_diversity2 %>%
  group_by(Experiment, WCS) %>%
    summarise(ratio=num_terms/num_modal_terms) %>%
    ungroup() %>%
  ggplot(aes(x=ratio, fill=WCS)) +
    geom_histogram() +
    xlab("Terms / Modal Terms") +
    ylab("Languages") +
    theme(legend.title=element_blank()) +
    scale_fill_manual(values=hues,
	              labels=c("English fixed choice",
	                       "English free choice",
			       "Spanish fixed choice",
			       "Spanish free choice",
			       "Tsimane' fixed choice",
			       "Tsimane' free choice",
			       "WCS"))                            		

ggsave("output/wcs_terms_over_modal_terms.pdf", width=6, height=4)

AllColorData = ColorData %>%
  select(Experiment, subject, grid_location, color) %>%
  rbind(WCS_with_subjects)

AllChipEntropies = ChipEntropies %>%
  unite(Experiment, Language, Type, sep="_") %>%
  select(Experiment, grid_location, Entropy, Prior) %>%
  rbind(WCS_ChipEntropies)

LangEntropies = AllChipEntropies %>%
  group_by(Experiment) %>%
    summarise(CondEnt=sum(Entropy * Prior)) %>%
    ungroup()

# Bootstrap the diversity measures ----------------------------

# bootstrap the Tsimane' fixed choice terms/modal_terms ratio to show that
# other languages with the same total number of terms are not from the same distribution,
# nor are other languages with the same total number of modal terms.

# resample subjects for Tsimane' fixed and free
# then do make_diversity
# then calculate ratios
# then look at WCS ratios and compare to the set of bootstrapped ratios

bootstrap_by_subject = function(d, f, num_samples, num_resamples) {
  subjects = unique(d$subject)
  if (is.na(num_resamples)) {
    num_resamples = length(subjects)
  }
  result = data.frame()
  for(i in 1:num_samples) {
    sub_subjects = sample(subjects, size=num_resamples, replace=F)
    subresult = data.frame()
    for(sub_subject in sub_subjects) {
      subresult = rbind(subresult, filter(d, subject == sub_subject))
    }
    result = rbind(result, f(subresult))
  }
  result
}

ts_fixed_boot_diversity = WCS_with_subjects_full %>%
  filter(Experiment == "Tsimane_Fixed") %>%
  bootstrap_by_subject(make_diversity, 10, NA)

ts_open_boot_diversity = WCS_with_subjects_full %>%
  filter(Experiment == "Tsimane_Open") %>%
  bootstrap_by_subject(make_diversity, 100, NA)

fixed_ratios = with(ts_fixed_boot_diversity, num_terms / num_modal_terms)
open_ratios  = with(ts_open_boot_diversity, num_terms / num_modal_terms) 

WCS_diversity_full %>%
  filter(WCS) %>%
  mutate(ratio=num_terms/num_modal_terms,
         less=ratio<min(open_ratios),
	 more=ratio>max(open_ratios),
	 significant = less | more) %>%
  select(Experiment, ratio, significant) %>%
  write.csv("output/ratio_significance.csv")

WCS_diversity_full %>%
  filter(WCS) %>%
  filter(num_modal_terms == 8) %>%
  mutate(ratio=num_terms/num_modal_terms,
         less=ratio<min(open_ratios),
	 more=ratio>max(open_ratios),
	 significant=less | more) %>%
  select(Experiment, ratio, significant) 


# WCS Uncertainty Plots -------------------------------

circle_small_filled=20
circle_filled=16
circle_empty=1
square_filled=15
square_empty=0

color_white="#FFFFFF"
color_gray0="#DEDEDE"
color_gray1="#B3B3B3"
color_gray2="#808080"
color_gray3="#363636"
color_black="#000000"
color_spanish="#A9A9A9"
color_english="#D3D3D3"

GeneratePlot<-function(data, method){
  # Shape order on data: Fixed, Open.
  # Color order on data: English, Spanish, Tsimane', WCS
  # Bevil's tones: WCS white, English light grey, Spanish dark grey, Tsimane black
  smallsize=3
  largesize=5
  if (method=="log2"){
    temp<-mutate(data,xcolors=log2(colors))
  }
  if (method=="log"){
    temp<-mutate(data,xcolors=log(colors))
  }
  if (method=="linear"){
    temp<-mutate(data,xcolors=colors)
  }
  max_x = ceiling(max(data$Uncertainty))
  max_y = ceiling(log(max(data$colors)))
  temp %>% ggplot(aes(x=xcolors,y=Uncertainty,shape=Method,color=Language,size=Type))+
    geom_point()+
    geom_point(shape = 1,size = smallsize,colour = "black")+
    theme_bw()+
    scale_shape_manual(values=c(circle_filled,square_filled))+
    scale_size_manual(values=c(largesize,smallsize))+
    scale_colour_manual(values=c(color_english,color_spanish,color_black,color_white))+
    scale_x_continuous("")+
    scale_y_continuous("",limit=c(3.5,6),breaks=c(3.5,4,4.5,5,5.5,6))+
    theme(panel.border = element_blank(),axis.line = element_line(colour="black"),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank())
}

LangEntropies %>%
  inner_join(select(WCS_diversity_full, Experiment, full_num_terms)) %>%
  separate(Experiment, into=c("Language", "Method"), sep="_") %>%
  mutate(Type=ifelse(is.na(Method), "WCS", "Non-WCS")) %>%
  rename(Uncertainty=CondEnt, colors=full_num_terms) %>%
  mutate(Language = ifelse(Language %in% c("English", "Spanish", "Tsimane"), Language, "WCS")) %>%
  GeneratePlot("log2")

ggsave("output/uncertainty.pdf", height=4, width=6)

