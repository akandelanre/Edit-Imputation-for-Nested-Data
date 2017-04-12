


###########################################################################
###########################################################################
####### MI using NDPMPM: Use rejection sampler to sample missing data #####
####### Add weighting option for sampling impossibles #####################
####### Add the Rejection sampler/sampling proposals hybrid option ########
####### Also make household head into a household level variable ##########
###########################################################################
###########################################################################

################################# START ###################################

####################### Step 1: Get True Population Data #################
rm(list = ls())
###### 1: Import Data
House <- read.csv("Data/House.csv",header=T)
Indiv <- read.csv("Data/Indiv.csv",header=T)


###### 2: Remove Households with size < 2 and > 4
House <- House[which(House$NP >= 2 & House$NP <= 4),]


###### 3: Keep only Households with TEN == 1,2,or 3 and recode 1,2 as 1 and 3 as 2
House <- House[which(House$TEN == 1 | House$TEN == 2 | House$TEN == 3),]
House$TEN[which(House$TEN == 2)] <- 1
House$TEN[which(House$TEN == 3)] <- 2


###### 4: Pick the same households in the indiv data
pick_index <- is.element(Indiv$SERIALNO,House$SERIALNO)
Indiv <- Indiv[pick_index,]


###### 5: Recode within-household variables
###### 5a: First, the relationship variable
Indiv$RELP[which(Indiv$RELP == 11 | Indiv$RELP == 12 | Indiv$RELP == 13)] <- 12 #Boarder, roommate or partner
Indiv$RELP[which(Indiv$RELP == 14 | Indiv$RELP == 15)] <- 13 #Other non-relative or foster child
Indiv$RELP[which(Indiv$RELP == 10)] <- 11 #Other relative
Indiv$RELP[which(Indiv$RELP == 9)] <- 10 #Child-in-law
Indiv$RELP[which(Indiv$RELP == 8)] <- 9 #Parent-in-law
Indiv$RELP[which(Indiv$RELP == 7)] <- 8 #Grandchild
Indiv$RELP[which(Indiv$RELP == 6)] <- 7 #Parent
Indiv$RELP[which(Indiv$RELP == 5)] <- 6 #Sibling
Indiv$RELP[which(Indiv$RELP == 4)] <- 5 #Stepchild
Indiv$RELP[which(Indiv$RELP == 3)] <- 4 #Adopted child
Indiv$RELP[which(Indiv$RELP == 2)] <- 3 #Biological child
Indiv$RELP[which(Indiv$RELP == 1)] <- 2 #Spouse
Indiv$RELP[which(Indiv$RELP == 0)] <- 1 #Household head
###### 5b: Next, the race variable
Indiv$RAC3P[which(Indiv$RAC3P == 4 | Indiv$RAC3P == 8| Indiv$RAC3P == 9 | Indiv$RAC3P == 10)] <- 6
Indiv$RAC3P[which(Indiv$RAC3P == 5)] <- 4
Indiv$RAC3P[which(Indiv$RAC3P == 7)] <- 5
Indiv$RAC3P[which(Indiv$RAC3P >= 11 & Indiv$RAC3P <= 15)] <- 7
Indiv$RAC3P[which(Indiv$RAC3P >= 16 & Indiv$RAC3P <= 59)] <- 8
Indiv$RAC3P[which(Indiv$RAC3P >= 60 & Indiv$RAC3P <= 100)] <- 9
###### 5c: Next, the hisp variable
Indiv$HISP[which(Indiv$HISP >= 5 & Indiv$HISP <= 24)] <- 5
###### 5d: Lastly, age
Indiv$AGEP <- Indiv$AGEP + 1L


###### 6: Make household head into household level data
HHhead_data <- Indiv[which(Indiv$SPORDER==1),]
Indiv_minHH <- Indiv[-which(Indiv$SPORDER==1),]


###### 7: Combine household head with household level data, rename all columns
Indiv_variables <- c("SEX","RAC3P","HISP","AGEP","RELP")
X_house <- cbind((House$NP-1L),House$TEN,HHhead_data[,Indiv_variables])
X_indiv <- Indiv_minHH[,Indiv_variables]
colnames(X_house) <- c("HHSize","Owner","HHGender","HHRace","HHHisp","HHAge","HHRelate")
colnames(X_indiv) <- c("Gender","Race","Hisp","Age","Relate")


###### 8: Save!!!
write.table(X_house, file = "Results/Population_House.txt",row.names = FALSE)
write.table(X_indiv, file = "Results/Population_Indiv.txt",row.names = FALSE)
########################## End of Step 1 ########################## 


###########################################################################
###########################################################################
#############################   BLANK SPACE   #############################
###########################################################################
###########################################################################
###########################################################################


####################### Step 2: Model Assessment ##########################
rm(list = ls())


###### 1a: Load saved results
Population_House <- read.table("Results/Population_House.txt",header=TRUE)
Population_Indiv <- read.table("Results/Population_Indiv.txt",header=TRUE)
Data_house_truth <- read.table("Results/Data_house_truth.txt",header=TRUE)
Data_indiv_truth <- read.table("Results/Data_indiv_truth.txt",header=TRUE)
dp_imput_house <- read.table("Results/dp_imput_house.txt",header=TRUE)
dp_imput_indiv <- read.table("Results/dp_imput_indiv.txt",header=TRUE)
#dp_imput_house_weighted <- 
#  read.table("Results/dp_imput_house_weighted.txt",header=TRUE)
#dp_imput_indiv_weighted <- 
#  read.table("Results/dp_imput_indiv_weighted.txt",header=TRUE)
weight_option <- FALSE

###### 1b: Define parameters
mm <- 50;
N <- nrow(Data_indiv_truth)
n <- nrow(Data_house_truth)
n_i <- as.numeric(as.character(Data_house_truth[,1]))
p <- ncol(Data_indiv_truth)
q <- ncol(Data_house_truth)
house_index <- rep(c(1:n),n_i)
H <- sort(unique(n_i))
#N_cc <- nrow(Data_indiv_cc)
#n_cc <- nrow(Data_house_cc)
#n_i_cc <- as.numeric(as.character(Data_house_cc[,1]))
#house_index_cc <- rep(c(1:n_cc),n_i_cc)
N_Pop <- nrow(Population_Indiv)
n_Pop <- nrow(Population_House)
n_i_Pop <- as.numeric(as.character(Population_House[,1]))
house_index_Pop <- rep(c(1:n_Pop),n_i_Pop)

###### 1c: Define categories for RELATE, GENDER & RACE variable
HEAD <- 1
SPOUSE <- 2
BIOLOGICALCHILD <- 3
ADOPTEDCHILD <- 4
STEPCHILD <- 5
SIBLING <- 6
PARENT <- 7
GRANDCHILD <- 8
PARENTINLAW <- 9
CHILDINLAW <- 10
WHITE <- 1; BLACK <- 2
MALE <- 1; FEMALE <- 2
NONHISP <- 1; OWN <- 1


###### 2: Define quantities/estimands of interest
Estimands <- c("Same race household, $n_i = 2$", "Same race household, $n_i = 3$", 
               "Same race household, $n_i = 4$","Spouse present","Black HH, own",
               "Spouse present, HH is White","Spouse present, HH is Black",
               "White couple", "Non-White couple, own","Same race couple", 
               "White-nonwhite couple", "Only one parent",
               "At least one biological child present","One grandchild present",
               "At least three generations present",
               "Couples with age difference less than 5",
               "HH older than spouse, White HH","White HH with Hispanic origin",
               "Only one adult female in house with at least one child less than 5", 
               "Only one adult Hispanic male in house with at least one child less than 10",
               "Only one adult Black female in house with at least one child less than 18",
               "HH over 35, no children present","Male HH, own",
               "Black HH younger than 40, own","White HH younger than 25, own",
               "Hispanic HH older than 50, own",
               "2 generations present, Black HH",
               "3 generations present, White couple",
               "At least 2 generations present, Hispanic couple",
               "At least one stepchild",
               "At least one adopted child, White couple",
               "At least one biological child, Hispanic couple",
               "At least two biological children, Black couple")


###### 3: Calculate the estimands
###### 3a: Population
Probs_Pop <- matrix(0,nrow=length(Estimands))
samp_size_Pop <- Probs_Pop
samp_size_Pop[1:length(H)] <- table(n_i_Pop)
samp_size_Pop[samp_size_Pop==0] <- n_Pop

for(kk in 1:n_Pop){
  hh_check <- Population_House[kk,(q-p+1):q]
  colnames(hh_check) <- colnames(Population_Indiv)
  hh_check <- rbind(hh_check,Population_Indiv[which(house_index_Pop==kk),])
  hh_check <- data.frame(Owner=Population_House[kk,"Owner"],hh_check)
  
  if(nrow(hh_check)==2 && hh_check[1,"Race"]==hh_check[2,"Race"]){
    Probs_Pop[1] <- Probs_Pop[1] + 1 #Same race household, n_i = 2
  }
  if(nrow(hh_check)==3 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
     hh_check[2,"Race"]==hh_check[3,"Race"]){
    Probs_Pop[2] <- Probs_Pop[2] + 1 #Same race household, n_i = 3
  }
  if(nrow(hh_check)==4 &&
     hh_check[1,"Race"]==hh_check[2,"Race"]&&
     hh_check[2,"Race"]==hh_check[3,"Race"]&&
     hh_check[3,"Race"]==hh_check[4,"Race"]){
    Probs_Pop[3] <- Probs_Pop[3] + 1 #Same race household, n_i = 4
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1){
    Probs_Pop[1 + length(H)] <- Probs_Pop[1 + length(H)] + 1 #Spouse present
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==BLACK)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs_Pop[2 + length(H)] <- Probs_Pop[2 + length(H)] + 1 #Black HH, own
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1 && 
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs_Pop[3 + length(H)] <- Probs_Pop[3 + length(H)] + 1 #Spouse present, HH is white
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
    Probs_Pop[4 + length(H)] <- Probs_Pop[4 + length(H)] + 1 #Spouse present, HH is black
  }
  if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs_Pop[5 + length(H)] <- Probs_Pop[5 + length(H)] + 1 #White Couple
  }
  if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]!=WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]!=WHITE)==1 &&
     hh_check[1,"Owner"] == OWN){
    Probs_Pop[6 + length(H)] <- Probs_Pop[6 + length(H)] + 1 #Non-White Couple, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
    Probs_Pop[7 + length(H)] <- Probs_Pop[7 + length(H)] + 1 #Same race couple
  }
  if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 &&
      sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
     (sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 &&
      sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1)){
    Probs_Pop[8 + length(H)] <- Probs_Pop[8 + length(H)] + 1 #White-nonwhite couple
  }
  if(sum(hh_check[,"Relate"]==PARENT)==1){
    Probs_Pop[9 + length(H)] <- Probs_Pop[9 + length(H)] + 1 #Only one parent   
  }
  if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1){
    Probs_Pop[10 + length(H)] <- Probs_Pop[10 + length(H)] + 1 #At least one biological child present
  }
  if(sum(hh_check[,"Relate"]==GRANDCHILD)==1){
    Probs_Pop[11 + length(H)] <- Probs_Pop[11 + length(H)] + 1 #One grandchild present
  }
  if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
     ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
             sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
     ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
    Probs_Pop[12 + length(H)] <- Probs_Pop[12 + length(H)] + 1 #At least three generations present
  }
  if(sum(abs(hh_check[hh_check[,"Relate"]==HEAD,"Age"]-hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])<5)==1){
    Probs_Pop[13 + length(H)] <- Probs_Pop[13 + length(H)] + 1 #Couples with age difference less than 5
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
     sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
    Probs_Pop[14 + length(H)] <- Probs_Pop[14 + length(H)] + 1 #HH older than spouse, white HH
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
         hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
    Probs_Pop[15 + length(H)] <- Probs_Pop[15 + length(H)] + 1 #White HH with hisp origin
  }
  if(sum(hh_check[hh_check[,"Gender"]==FEMALE,"Age"]>=21)==1 &&
     ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<5)>=1 |
             sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<5)>=1 |
             sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<5)>=1),1,0)==1){
    Probs_Pop[16 + length(H)] <- Probs_Pop[16 + length(H)] + 1
    #Only one adult female (>=21) in house with at least one child less than 5
  }
  if(sum(hh_check[hh_check[,"Gender"]==MALE & hh_check[,"Hisp"]!=NONHISP,"Age"]>=21)==1 &&
     ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<10)>=1 |
             sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<10)>=1 |
             sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<10)>=1),1,0)==1){
    Probs_Pop[17 + length(H)] <- Probs_Pop[17 + length(H)] + 1
    #Only one adult male hispanic (>=21) in house with at least one child less than 10
  }
  if(sum(hh_check[hh_check[,"Gender"]==FEMALE & hh_check[,"Race"]==BLACK,"Age"]>=21)==1 &&
     ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<18)>=1 |
             sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<18)>=1 |
             sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<18)>=1),1,0)==1){
    Probs_Pop[18 + length(H)] <- Probs_Pop[18 + length(H)] + 1
    #Only one adult female black (>=21) in house with at least one child less than 18
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>35)==1 &&
     ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)==0 &
             sum(hh_check[,"Relate"]==ADOPTEDCHILD)==0 & sum(hh_check[,"Relate"]==STEPCHILD)==0),1,0)==1){
    Probs_Pop[19 + length(H)] <- Probs_Pop[19 + length(H)] + 1 #HH over 35, no children present
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Gender"]==MALE)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs_Pop[20 + length(H)] <- Probs_Pop[20 + length(H)] + 1 #Male HH, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK,"Age"]<40)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs_Pop[21 + length(H)] <- Probs_Pop[21 + length(H)] + 1 #Black HH younger than 40, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE,"Age"]<25)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs_Pop[22 + length(H)] <- Probs_Pop[22 + length(H)] + 1 #White HH younger than 25, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP,"Age"]>50)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs_Pop[23 + length(H)] <- Probs_Pop[23 + length(H)] + 1 #Hispanic HH older than 50, own
  }
  if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
      ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
              sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
      ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)==1) &&
     sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==BLACK)==1){
    Probs_Pop[24 + length(H)] <- Probs_Pop[24 + length(H)] + 1 #2 generations present, black HH
  }
  if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
      ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
              sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
      ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)==2) &&
     sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs_Pop[25 + length(H)] <- Probs_Pop[25 + length(H)] + 1 #3 generations present, white couple
  }
  if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
      ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
              sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
      ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=1) &&
     sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Hisp"]!=NONHISP)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP)==1){
    Probs_Pop[26 + length(H)] <- Probs_Pop[26 + length(H)] + 1 #At least 2 generations present, hispanic couple
  }
  if(sum(hh_check[,"Relate"]==STEPCHILD)>=1){
    Probs_Pop[27 + length(H)] <- Probs_Pop[27 + length(H)] + 1 #At least one stepchild
  }
  if(sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 &&
     sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs_Pop[28 + length(H)] <- Probs_Pop[28 + length(H)] + 1 #At least one adopted child, white couple
  }
  if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 &&
     sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Hisp"]!=NONHISP)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP)==1){
    Probs_Pop[29 + length(H)] <- Probs_Pop[29 + length(H)] + 1 #At least one biological child, hispanic couple
  }
  if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=2 &&
     sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==BLACK)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
    Probs_Pop[30 + length(H)] <- Probs_Pop[30 + length(H)] + 1 #At least two biological children, Black couple
  }
}
Probs_Pop <- Probs_Pop/samp_size_Pop



###### 3b: Original Data
Probs <- matrix(0,nrow=length(Estimands))
samp_size <- Probs
samp_size[1:length(H)] <- table(n_i)
samp_size[samp_size==0] <- n

for(kk in 1:n){
  hh_check <- Data_house_truth[kk,(q-p+1):q]
  colnames(hh_check) <- colnames(Data_indiv_truth)
  hh_check <- rbind(hh_check,Data_indiv_truth[which(house_index==kk),])
  hh_check <- data.frame(Owner=Data_house_truth[kk,"Owner"],hh_check)
  
  if(nrow(hh_check)==2 && hh_check[1,"Race"]==hh_check[2,"Race"]){
    Probs[1] <- Probs[1] + 1 #Same race household, n_i = 2
  }
  if(nrow(hh_check)==3 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
     hh_check[2,"Race"]==hh_check[3,"Race"]){
    Probs[2] <- Probs[2] + 1 #Same race household, n_i = 3
  }
  if(nrow(hh_check)==4 &&
     hh_check[1,"Race"]==hh_check[2,"Race"]&&
     hh_check[2,"Race"]==hh_check[3,"Race"]&&
     hh_check[3,"Race"]==hh_check[4,"Race"]){
    Probs[3] <- Probs[3] + 1 #Same race household, n_i = 4
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1){
    Probs[1 + length(H)] <- Probs[1 + length(H)] + 1 #Spouse present
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==BLACK)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs[2 + length(H)] <- Probs[2 + length(H)] + 1 #Black HH, own
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1 && 
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs[3 + length(H)] <- Probs[3 + length(H)] + 1 #Spouse present, HH is white
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
    Probs[4 + length(H)] <- Probs[4 + length(H)] + 1 #Spouse present, HH is black
  }
  if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs[5 + length(H)] <- Probs[5 + length(H)] + 1 #White Couple
  }
  if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]!=WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]!=WHITE)==1 &&
     hh_check[1,"Owner"] == OWN){
    Probs[6 + length(H)] <- Probs[6 + length(H)] + 1 #Non-White Couple, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
    Probs[7 + length(H)] <- Probs[7 + length(H)] + 1 #Same race couple
  }
  if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 &&
      sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
     (sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 &&
     sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1)){
    Probs[8 + length(H)] <- Probs[8 + length(H)] + 1 #White-nonwhite couple
  }
  if(sum(hh_check[,"Relate"]==PARENT)==1){
    Probs[9 + length(H)] <- Probs[9 + length(H)] + 1 #Only one parent   
  }
  if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1){
    Probs[10 + length(H)] <- Probs[10 + length(H)] + 1 #At least one biological child present
  }
  if(sum(hh_check[,"Relate"]==GRANDCHILD)==1){
    Probs[11 + length(H)] <- Probs[11 + length(H)] + 1 #One grandchild present
  }
  if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
     ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
             sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
     ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
    Probs[12 + length(H)] <- Probs[12 + length(H)] + 1 #At least three generations present
  }
  if(sum(abs(hh_check[hh_check[,"Relate"]==HEAD,"Age"]-hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])<5)==1){
    Probs[13 + length(H)] <- Probs[13 + length(H)] + 1 #Couples with age difference less than 5
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
     sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
    Probs[14 + length(H)] <- Probs[14 + length(H)] + 1 #HH older than spouse, white HH
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
     hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
    Probs[15 + length(H)] <- Probs[15 + length(H)] + 1 #White HH with hisp origin
  }
  if(sum(hh_check[hh_check[,"Gender"]==FEMALE,"Age"]>=21)==1 &&
     ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<5)>=1 |
             sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<5)>=1 |
             sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<5)>=1),1,0)==1){
    Probs[16 + length(H)] <- Probs[16 + length(H)] + 1
    #Only one adult female (>=21) in house with at least one child less than 5
  }
  if(sum(hh_check[hh_check[,"Gender"]==MALE & hh_check[,"Hisp"]!=NONHISP,"Age"]>=21)==1 &&
     ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<10)>=1 |
             sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<10)>=1 |
             sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<10)>=1),1,0)==1){
    Probs[17 + length(H)] <- Probs[17 + length(H)] + 1
    #Only one adult male hispanic (>=21) in house with at least one child less than 10
  }
  if(sum(hh_check[hh_check[,"Gender"]==FEMALE & hh_check[,"Race"]==BLACK,"Age"]>=21)==1 &&
     ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<18)>=1 |
             sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<18)>=1 |
             sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<18)>=1),1,0)==1){
    Probs[18 + length(H)] <- Probs[18 + length(H)] + 1
    #Only one adult female black (>=21) in house with at least one child less than 18
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>35)==1 &&
     ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)==0 &
             sum(hh_check[,"Relate"]==ADOPTEDCHILD)==0 & sum(hh_check[,"Relate"]==STEPCHILD)==0),1,0)==1){
    Probs[19 + length(H)] <- Probs[19 + length(H)] + 1 #HH over 35, no children present
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Gender"]==MALE)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs[20 + length(H)] <- Probs[20 + length(H)] + 1 #Male HH, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK,"Age"]<40)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs[21 + length(H)] <- Probs[21 + length(H)] + 1 #Black HH younger than 40, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE,"Age"]<25)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs[22 + length(H)] <- Probs[22 + length(H)] + 1 #White HH younger than 25, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP,"Age"]>50)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs[23 + length(H)] <- Probs[23 + length(H)] + 1 #Hispanic HH older than 50, own
  }
  if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
     ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
             sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
     ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)==1) &&
     sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==BLACK)==1){
    Probs[24 + length(H)] <- Probs[24 + length(H)] + 1 #2 generations present, black HH
  }
  if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
      ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
              sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
      ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)==2) &&
     sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs[25 + length(H)] <- Probs[25 + length(H)] + 1 #3 generations present, white couple
  }
  if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
      ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
              sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
      ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=1) &&
     sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Hisp"]!=NONHISP)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP)==1){
    Probs[26 + length(H)] <- Probs[26 + length(H)] + 1 #At least 2 generations present, hispanic couple
  }
  if(sum(hh_check[,"Relate"]==STEPCHILD)>=1){
    Probs[27 + length(H)] <- Probs[27 + length(H)] + 1 #At least one stepchild
  }
  if(sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 &&
     sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs[28 + length(H)] <- Probs[28 + length(H)] + 1 #At least one adopted child, white couple
  }
  if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 &&
     sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Hisp"]!=NONHISP)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP)==1){
    Probs[29 + length(H)] <- Probs[29 + length(H)] + 1 #At least one biological child, hispanic couple
  }
  if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=2 &&
     sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==BLACK)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
    Probs[30 + length(H)] <- Probs[30 + length(H)] + 1 #At least two biological children, Black couple
  }
}
Probs <- Probs/samp_size
V <- (Probs*(1-Probs))/samp_size
CIntLower <- Probs + (qnorm(0.025)*sqrt(V))
CIntUpper <- Probs - (qnorm(0.025)*sqrt(V))
CInt <- cbind(CIntLower,CIntUpper)



###### 2c: Regular Model
Probs_syn <- matrix(0,nrow=length(Estimands),ncol=mm)
V_syn <- matrix(0,nrow=length(Estimands),ncol=mm)
for(k in 1:mm){
  k_imp_house <- dp_imput_house[((n*(k-1))+1):(n*k),]
  k_imp_indiv <- dp_imput_indiv[((N*(k-1))+1):(N*k),]
  new_n_i <- as.numeric(as.character(k_imp_house[,"HHSize"]))
  new_house_index <- rep(c(1:n),new_n_i)
  for(kk in 1:n){
    hh_check <- k_imp_house[kk,(q-p+1):q]
    colnames(hh_check) <- colnames(k_imp_indiv)
    hh_check <- rbind(hh_check,k_imp_indiv[which(new_house_index==kk),])
    hh_check <- data.frame(Owner=k_imp_house[kk,"Owner"],hh_check)
    
    if(nrow(hh_check)==2 && hh_check[1,"Race"]==hh_check[2,"Race"]){
      Probs_syn[1,k] <- Probs_syn[1,k] + 1 #Same race household, n_i = 2
    }
    if(nrow(hh_check)==3 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]){
      Probs_syn[2,k] <- Probs_syn[2,k] + 1 #Same race household, n_i = 3
    }
    if(nrow(hh_check)==4 &&
       hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]&&
       hh_check[3,"Race"]==hh_check[4,"Race"]){
      Probs_syn[3,k] <- Probs_syn[3,k] + 1 #Same race household, n_i = 4
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1){
      Probs_syn[1 + length(H),k] <- Probs_syn[1 + length(H),k] + 1 #Spouse present
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==BLACK)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn[2 + length(H),k] <- Probs_syn[2 + length(H),k] + 1 #Black HH, own
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && 
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn[3 + length(H),k] <- Probs_syn[3 + length(H),k] + 1 #Spouse present, HH is white
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
      Probs_syn[4 + length(H),k] <- Probs_syn[4 + length(H),k] + 1 #Spouse present, HH is black
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn[5 + length(H),k] <- Probs_syn[5 + length(H),k] + 1 #White Couple
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]!=WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]!=WHITE)==1 &&
       hh_check[1,"Owner"] == OWN){
      Probs_syn[6 + length(H),k] <- Probs_syn[6 + length(H),k] + 1 #Non-White Couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
      Probs_syn[7 + length(H),k] <- Probs_syn[7 + length(H),k] + 1 #Same race couple
    }
    if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 &&
        sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
       (sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 &&
        sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1)){
      Probs_syn[8 + length(H),k] <- Probs_syn[8 + length(H),k] + 1 #White-nonwhite couple
    }
    if(sum(hh_check[,"Relate"]==PARENT)==1){
      Probs_syn[9 + length(H),k] <- Probs_syn[9 + length(H),k] + 1 #Only one parent   
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1){
      Probs_syn[10 + length(H),k] <- Probs_syn[10 + length(H),k] + 1 #At least one biological child present
    }
    if(sum(hh_check[,"Relate"]==GRANDCHILD)==1){
      Probs_syn[11 + length(H),k] <- Probs_syn[11 + length(H),k] + 1 #One grandchild present
    }
    if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
       ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
               sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
       ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
      Probs_syn[12 + length(H),k] <- Probs_syn[12 + length(H),k] + 1 #At least three generations present
    }
    if(sum(abs(hh_check[hh_check[,"Relate"]==HEAD,"Age"]-hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])<5)==1){
      Probs_syn[13 + length(H),k] <- Probs_syn[13 + length(H),k] + 1 #Couples with age difference less than 5
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
       sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn[14 + length(H),k] <- Probs_syn[14 + length(H),k] + 1 #HH older than spouse, white HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
           hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn[15 + length(H),k] <- Probs_syn[15 + length(H),k] + 1 #White HH with hisp origin
    }
    if(sum(hh_check[hh_check[,"Gender"]==FEMALE,"Age"]>=21)==1 &&
       ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<5)>=1 |
               sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<5)>=1 |
               sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<5)>=1),1,0)==1){
      Probs_syn[16 + length(H),k] <- Probs_syn[16 + length(H),k] + 1
      #Only one adult female (>=21) in house with at least one child less than 5
    }
    if(sum(hh_check[hh_check[,"Gender"]==MALE & hh_check[,"Hisp"]!=NONHISP,"Age"]>=21)==1 &&
       ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<10)>=1 |
               sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<10)>=1 |
               sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<10)>=1),1,0)==1){
      Probs_syn[17 + length(H),k] <- Probs_syn[17 + length(H),k] + 1
      #Only one adult male hispanic (>=21) in house with at least one child less than 10
    }
    if(sum(hh_check[hh_check[,"Gender"]==FEMALE & hh_check[,"Race"]==BLACK,"Age"]>=21)==1 &&
       ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<18)>=1 |
               sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<18)>=1 |
               sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<18)>=1),1,0)==1){
      Probs_syn[18 + length(H),k] <- Probs_syn[18 + length(H),k] + 1
      #Only one adult female black (>=21) in house with at least one child less than 18
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>35)==1 &&
       ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)==0 &
               sum(hh_check[,"Relate"]==ADOPTEDCHILD)==0 & sum(hh_check[,"Relate"]==STEPCHILD)==0),1,0)==1){
      Probs_syn[19 + length(H),k] <- Probs_syn[19 + length(H),k] + 1 #HH over 35, no children present
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Gender"]==MALE)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn[20 + length(H),k] <- Probs_syn[20 + length(H),k] + 1 #Male HH, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK,"Age"]<40)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn[21 + length(H),k] <- Probs_syn[21 + length(H),k] + 1 #Black HH younger than 40, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE,"Age"]<25)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn[22 + length(H),k] <- Probs_syn[22 + length(H),k] + 1 #White HH younger than 25, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP,"Age"]>50)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn[23 + length(H),k] <- Probs_syn[23 + length(H),k] + 1 #Hispanic HH older than 50, own
    }
    if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
        ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
                sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
        ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)==1) &&
       sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==BLACK)==1){
      Probs_syn[24 + length(H),k] <- Probs_syn[24 + length(H),k] + 1 #2 generations present, black HH
    }
    if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
        ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
                sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
        ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)==2) &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn[25 + length(H),k] <- Probs_syn[25 + length(H),k] + 1 #3 generations present, white couple
    }
    if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
        ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
                sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
        ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=1) &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Hisp"]!=NONHISP)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP)==1){
      Probs_syn[26 + length(H),k] <- Probs_syn[26 + length(H),k] + 1 #At least 2 generations present, hispanic couple
    }
    if(sum(hh_check[,"Relate"]==STEPCHILD)>=1){
      Probs_syn[27 + length(H),k] <- Probs_syn[27 + length(H),k] + 1 #At least one stepchild
    }
    if(sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn[28 + length(H),k] <- Probs_syn[28 + length(H),k] + 1 #At least one adopted child, white couple
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Hisp"]!=NONHISP)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP)==1){
      Probs_syn[29 + length(H),k] <- Probs_syn[29 + length(H),k] + 1 #At least one biological child, hispanic couple
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=2 &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==BLACK)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
      Probs_syn[30 + length(H),k] <- Probs_syn[30 + length(H),k] + 1 #At least two biological children, Black couple
    }
  }
  Probs_syn[,k] <- Probs_syn[,k]/samp_size
  V_syn[,k] <- (Probs_syn[,k]*(1-Probs_syn[,k]))/samp_size
}
dp_qbar_syn <- rowMeans(Probs_syn)
dp_b_syn <- apply(Probs_syn,1,var)
dp_ubar_syn <- rowMeans(V_syn)
dp_t_syn <- dp_ubar_syn + (dp_b_syn*(mm+1)/mm) #dp_t_syn <- dp_ubar_syn + (dp_b/mm) for synthetic data
dp_r_syn <- dp_ubar_syn/dp_b_syn
dp_v_syn <- (mm-1)*((1+((mm/(mm+1))*dp_r_syn))^2) #dp_v_syn <- (mm-1)*((1+((mm)*dp_r_syn))^2) for synthetic data
CIntLower_dp_syn <- dp_qbar_syn + (qt(0.025,dp_v_syn)*sqrt(dp_t_syn))
CIntUpper_dp_syn <- dp_qbar_syn - (qt(0.025,dp_v_syn)*sqrt(dp_t_syn))
CInt_syn <- cbind(CIntLower_dp_syn,CIntUpper_dp_syn)
r_syn <- (1+(1/mm))/dp_r_syn
lambda_syn <- (r_syn+(2/(dp_v_syn + 3)))/(r_syn + 1)
RE_syn <- 1/(1 + (lambda_syn/mm))



###### 2d: Weighted Sampler
if(weight_option){
Probs_syn_weighted <- matrix(0,nrow=length(Estimands),ncol=mm)
V_syn_weighted <- matrix(0,nrow=length(Estimands),ncol=mm)
for(k in 1:mm){
  k_imp_house <- dp_imput_house_weighted[((n*(k-1))+1):(n*k),]
  k_imp_indiv <- dp_imput_indiv_weighted[((N*(k-1))+1):(N*k),]
  new_n_i <- as.numeric(as.character(k_imp_house[,"HHSize"]))
  new_house_index <- rep(c(1:n),new_n_i)
  for(kk in 1:n){
    hh_check <- k_imp_house[kk,(q-p+1):q]
    colnames(hh_check) <- colnames(k_imp_indiv)
    hh_check <- rbind(hh_check,k_imp_indiv[which(new_house_index==kk),])
    hh_check <- data.frame(Owner=k_imp_house[kk,"Owner"],hh_check)
    
    if(nrow(hh_check)==2 && hh_check[1,"Race"]==hh_check[2,"Race"]){
      Probs_syn_weighted[1,k] <- Probs_syn_weighted[1,k] + 1 #Same race household, n_i = 2
    }
    if(nrow(hh_check)==3 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]){
      Probs_syn_weighted[2,k] <- Probs_syn_weighted[2,k] + 1 #Same race household, n_i = 3
    }
    if(nrow(hh_check)==4 &&
       hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]&&
       hh_check[3,"Race"]==hh_check[4,"Race"]){
      Probs_syn_weighted[3,k] <- Probs_syn_weighted[3,k] + 1 #Same race household, n_i = 4
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1){
      Probs_syn_weighted[1 + length(H),k] <- Probs_syn_weighted[1 + length(H),k] + 1 #Spouse present
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==BLACK)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn_weighted[2 + length(H),k] <- Probs_syn_weighted[2 + length(H),k] + 1 #Black HH, own
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && 
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_weighted[3 + length(H),k] <- Probs_syn_weighted[3 + length(H),k] + 1 #Spouse present, HH is white
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
      Probs_syn_weighted[4 + length(H),k] <- Probs_syn_weighted[4 + length(H),k] + 1 #Spouse present, HH is black
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_weighted[5 + length(H),k] <- Probs_syn_weighted[5 + length(H),k] + 1 #White Couple
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]!=WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]!=WHITE)==1 &&
       hh_check[1,"Owner"] == OWN){
      Probs_syn_weighted[6 + length(H),k] <- Probs_syn_weighted[6 + length(H),k] + 1 #Non-White Couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
      Probs_syn_weighted[7 + length(H),k] <- Probs_syn_weighted[7 + length(H),k] + 1 #Same race couple
    }
    if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 &&
        sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
       (sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 &&
        sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1)){
      Probs_syn_weighted[8 + length(H),k] <- Probs_syn_weighted[8 + length(H),k] + 1 #White-nonwhite couple
    }
    if(sum(hh_check[,"Relate"]==PARENT)==1){
      Probs_syn_weighted[9 + length(H),k] <- Probs_syn_weighted[9 + length(H),k] + 1 #Only one parent   
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1){
      Probs_syn_weighted[10 + length(H),k] <- Probs_syn_weighted[10 + length(H),k] + 1 #At least one biological child present
    }
    if(sum(hh_check[,"Relate"]==GRANDCHILD)==1){
      Probs_syn_weighted[11 + length(H),k] <- Probs_syn_weighted[11 + length(H),k] + 1 #One grandchild present
    }
    if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
       ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
               sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
       ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
      Probs_syn_weighted[12 + length(H),k] <- Probs_syn_weighted[12 + length(H),k] + 1 #At least three generations present
    }
    if(sum(abs(hh_check[hh_check[,"Relate"]==HEAD,"Age"]-hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])<5)==1){
      Probs_syn_weighted[13 + length(H),k] <- Probs_syn_weighted[13 + length(H),k] + 1 #Couples with age difference less than 5
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
       sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn_weighted[14 + length(H),k] <- Probs_syn_weighted[14 + length(H),k] + 1 #HH older than spouse, white HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
           hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn_weighted[15 + length(H),k] <- Probs_syn_weighted[15 + length(H),k] + 1 #White HH with hisp origin
    }
    if(sum(hh_check[hh_check[,"Gender"]==FEMALE,"Age"]>=21)==1 &&
       ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<5)>=1 |
               sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<5)>=1 |
               sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<5)>=1),1,0)==1){
      Probs_syn_weighted[16 + length(H),k] <- Probs_syn_weighted[16 + length(H),k] + 1
      #Only one adult female (>=21) in house with at least one child less than 5
    }
    if(sum(hh_check[hh_check[,"Gender"]==MALE & hh_check[,"Hisp"]!=NONHISP,"Age"]>=21)==1 &&
       ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<10)>=1 |
               sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<10)>=1 |
               sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<10)>=1),1,0)==1){
      Probs_syn_weighted[17 + length(H),k] <- Probs_syn_weighted[17 + length(H),k] + 1
      #Only one adult male hispanic (>=21) in house with at least one child less than 10
    }
    if(sum(hh_check[hh_check[,"Gender"]==FEMALE & hh_check[,"Race"]==BLACK,"Age"]>=21)==1 &&
       ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<18)>=1 |
               sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<18)>=1 |
               sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<18)>=1),1,0)==1){
      Probs_syn_weighted[18 + length(H),k] <- Probs_syn_weighted[18 + length(H),k] + 1
      #Only one adult female black (>=21) in house with at least one child less than 18
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>35)==1 &&
       ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)==0 &
               sum(hh_check[,"Relate"]==ADOPTEDCHILD)==0 & sum(hh_check[,"Relate"]==STEPCHILD)==0),1,0)==1){
      Probs_syn_weighted[19 + length(H),k] <- Probs_syn_weighted[19 + length(H),k] + 1 #HH over 35, no children present
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Gender"]==MALE)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn_weighted[20 + length(H),k] <- Probs_syn_weighted[20 + length(H),k] + 1 #Male HH, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK,"Age"]<40)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn_weighted[21 + length(H),k] <- Probs_syn_weighted[21 + length(H),k] + 1 #Black HH younger than 40, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE,"Age"]<25)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn_weighted[22 + length(H),k] <- Probs_syn_weighted[22 + length(H),k] + 1 #White HH younger than 25, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP,"Age"]>50)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn_weighted[23 + length(H),k] <- Probs_syn_weighted[23 + length(H),k] + 1 #Hispanic HH older than 50, own
    }
    if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
        ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
                sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
        ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)==1) &&
       sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==BLACK)==1){
      Probs_syn_weighted[24 + length(H),k] <- Probs_syn_weighted[24 + length(H),k] + 1 #2 generations present, black HH
    }
    if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
        ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
                sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
        ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)==2) &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_weighted[25 + length(H),k] <- Probs_syn_weighted[25 + length(H),k] + 1 #3 generations present, white couple
    }
    if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
        ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
                sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
        ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=1) &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Hisp"]!=NONHISP)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP)==1){
      Probs_syn_weighted[26 + length(H),k] <- Probs_syn_weighted[26 + length(H),k] + 1 #At least 2 generations present, hispanic couple
    }
    if(sum(hh_check[,"Relate"]==STEPCHILD)>=1){
      Probs_syn_weighted[27 + length(H),k] <- Probs_syn_weighted[27 + length(H),k] + 1 #At least one stepchild
    }
    if(sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_weighted[28 + length(H),k] <- Probs_syn_weighted[28 + length(H),k] + 1 #At least one adopted child, white couple
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Hisp"]!=NONHISP)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP)==1){
      Probs_syn_weighted[29 + length(H),k] <- Probs_syn_weighted[29 + length(H),k] + 1 #At least one biological child, hispanic couple
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==BLACK)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
      Probs_syn_weighted[30 + length(H),k] <- Probs_syn_weighted[30 + length(H),k] + 1 #At least two biological children, Black couple
    }
  }
  Probs_syn_weighted[,k] <- Probs_syn_weighted[,k]/samp_size
  V_syn_weighted[,k] <- (Probs_syn_weighted[,k]*(1-Probs_syn_weighted[,k]))/samp_size
}
dp_qbar_syn_weighted <- rowMeans(Probs_syn_weighted)
dp_b_syn_weighted <- apply(Probs_syn_weighted,1,var)
dp_ubar_syn_weighted <- rowMeans(V_syn_weighted)
dp_t_syn_weighted <- dp_ubar_syn_weighted + (dp_b_syn_weighted*(mm+1)/mm)
#dp_t_syn_weighted <- dp_ubar_syn_weighted + (dp_b/mm) for synthetic data
dp_r_syn_weighted <- dp_ubar_syn_weighted/dp_b_syn_weighted
dp_v_syn_weighted <- (mm-1)*((1+((mm/(mm+1))*dp_r_syn_weighted))^2)
#dp_v_syn_weighted <- (mm-1)*((1+((mm)*dp_r_syn_weighted))^2) for synthetic data
CIntLower_dp_syn_weighted <- dp_qbar_syn_weighted + (qt(0.025,dp_v_syn_weighted)*sqrt(dp_t_syn_weighted))
CIntUpper_dp_syn_weighted <- dp_qbar_syn_weighted - (qt(0.025,dp_v_syn_weighted)*sqrt(dp_t_syn_weighted))
CInt_syn_weighted <- cbind(CIntLower_dp_syn_weighted,CIntUpper_dp_syn_weighted)
r_syn_weighted <- (1+(1/mm))/dp_r_syn_weighted
lambda_syn_weighted <- (r_syn_weighted+(2/(dp_v_syn_weighted + 3)))/(r_syn_weighted + 1)
RE_syn_weighted <- 1/(1 + (lambda_syn_weighted/mm))
}



###### 3: Combine and save!!!
if(weight_option){
  CompareProbs <- cbind(Probs_Pop,Probs,dp_qbar_syn,dp_qbar_syn_weighted,CInt,CInt_syn,CInt_syn_weighted,
                        lambda_syn,lambda_syn_weighted,RE_syn,RE_syn_weighted)
  colnames(CompareProbs) <- 
    c("Pop. Truth","Orig. Data Q","Model Q","Weighted Sampler Q",
      "Orig. Data L","Orig. Data U","Model L","Model U","Weighted Sampler L","Weighted Sampler U",
      "Model M.Frac","Weighted Sampler M.Frac","Model RE","Weighted Sampler RE")
} else {
  CompareProbs <- cbind(Probs_Pop,Probs,dp_qbar_syn,CInt,CInt_syn,lambda_syn,RE_syn)
  colnames(CompareProbs) <- c("Pop. Truth","Orig. Data Q","Model Q","Orig. Data L","Orig. Data U","Model L","Model U",
                              "Model M.Frac","Model RE")
}
rownames(CompareProbs) <- Estimands

write.table(CompareProbs,"Results/CompareProbs.txt",row.names = TRUE)
CompareProbs <- read.table("Results/CompareProbs.txt",header=TRUE)
round(CompareProbs,3)

#library(xtable)
#xtable(round(CompareProbs[,c(1:3)],3),digits=3)
#xtable(round(CompareProbs[,c(4:7)],3),digits=3)

########################## End of Step 1 ########################## 

###########################################################################
###########################################################################
#############################   BLANK SPACE   #############################
###########################################################################
###########################################################################
###########################################################################