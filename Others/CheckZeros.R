


rm(list = ls())
###### 1: Import Data
House <- read.csv("Data/House.csv",header=T)
Indiv <- read.csv("Data/Indiv.csv",header=T)


###### 2: Remove Households with size < 2 and > 4
House <- House[which(House$NP >= 2),]


###### 3: Keep only Households with TEN == 1,2,or 3 and recode 1,2 as 1 and 3 as 2
House <- House[which(House$TEN == 1 | House$TEN == 2 | House$TEN == 3),]
House$TEN[which(House$TEN == 2)] <- 1
House$TEN[which(House$TEN == 3)] <- 2


###### 4: Take a sample of size 2,000 Households
sample_size <- nrow(House)
samp_index <- sort(sample(1:nrow(House),sample_size,replace=F))
House <- House[samp_index,]


###### 5: Pick the same households in the indiv data
pick_index <- is.element(Indiv$SERIALNO,House$SERIALNO)
Indiv <- Indiv[pick_index,]


###### 6: Recode within-household variables
###### 6a: First, the relationship variable
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

###### 6b: Next, the race variable
Indiv$RAC3P[which(Indiv$RAC3P == 4 | Indiv$RAC3P == 8| Indiv$RAC3P == 9 | Indiv$RAC3P == 10)] <- 6
Indiv$RAC3P[which(Indiv$RAC3P == 5)] <- 4
Indiv$RAC3P[which(Indiv$RAC3P == 7)] <- 5
Indiv$RAC3P[which(Indiv$RAC3P >= 11 & Indiv$RAC3P <= 15)] <- 7
Indiv$RAC3P[which(Indiv$RAC3P >= 16 & Indiv$RAC3P <= 59)] <- 8
Indiv$RAC3P[which(Indiv$RAC3P >= 60 & Indiv$RAC3P <= 100)] <- 9

###### 6c: Next, the hisp variable
Indiv$HISP[which(Indiv$HISP >= 5 & Indiv$HISP <= 24)] <- 5

###### 6d: Lastly, age
Indiv$AGEP <- Indiv$AGEP + 1L


###### 7: Make household head and spouse into household level data
HHhead_data <- Indiv[which(Indiv$RELP==1),]
Indiv_minHH <- Indiv[-which(Indiv$RELP==1),]


###### 8: Combine Household and within-household data using the following ordering:
###### c("HHIndex","WithinHHIndex","Gender","Race","Hisp","Age","Relate","Owner")
#origdata <- data.frame(HHIndex = rep(c(1:sample_size),House$NP),WithinHHIndex = Indiv$SPORDER,
#                       Gender = Indiv$SEX,Race = Indiv$RAC3P,Hisp = Indiv$HISP,
#                       Age = Indiv$AGEP,Relate = Indiv$RELP,Owner = rep(House$TEN,House$NP))
#colnames(origdata) <- c("HHIndex","WithinHHIndex","Gender","Race","Hisp","Age","Relate","Owner")
origdata <- data.frame(HHIndex = rep(c(1:sample_size),(House$NP-1L)),
                       WithinHHIndex = Indiv_minHH$SPORDER,
                       Gender = Indiv_minHH$SEX,Race = Indiv_minHH$RAC3P,Hisp = Indiv_minHH$HISP,
                       Age = Indiv_minHH$AGEP,Relate = Indiv_minHH$RELP,
                       Owner = rep(House$TEN,(House$NP-1L)),
                       HHGender = rep(HHhead_data$SEX,(House$NP-1L)),
                       HHRace = rep(HHhead_data$RAC3P,(House$NP-1L)),
                       HHHisp = rep(HHhead_data$HISP,(House$NP-1L)),
                       HHAge = rep(HHhead_data$AGEP,(House$NP-1L)),
                       HHRelate = rep(HHhead_data$RELP,(House$NP-1L)))


#HEAD 1
#SPOUSE 2
#BIOLOGICALCHILD 3
#ADOPTEDCHILD 4
#STEPCHILD 5
#SIBLING 6
#PARENT 7
#GRANDCHILD 8
#PARENTINLAW 9
#CHILDINLAW 10

#HEAD 1
summary(origdata$HHAge) #Min head age is 16

#SPOUSE 2
summary(origdata$Age[which(origdata$Relate == 2)]) #Min head age is 16
summary(abs(origdata$HHAge[which(origdata$Relate == 2)] - origdata$Age[which(origdata$Relate == 2)])) #Max age diff between couples is 49
length(which(origdata$Gender[which(origdata$Relate == 2)] == 1 & origdata$HHGender[which(origdata$Relate == 2)] == 1))

#BIOLOGICALCHILD 3
summary((origdata$HHAge[which(origdata$Relate == 3)] - origdata$Age[which(origdata$Relate == 3)])) #Min age diff btw head and child is 7

#ADOPTEDCHILD 4
summary((origdata$HHAge[which(origdata$Relate == 4)] - origdata$Age[which(origdata$Relate == 4)])) #Min age diff btw head and adopted child is 11

#STEPCHILD 5
summary((origdata$HHAge[which(origdata$Relate == 5)] - origdata$Age[which(origdata$Relate == 5)])) #Min age diff btw head and step child is 9

#SIBLING 6
summary(abs(origdata$HHAge[which(origdata$Relate == 6)] - origdata$Age[which(origdata$Relate == 6)])) #Max age diff btw head and siblings is 37

#PARENT 7
summary(-(origdata$HHAge[which(origdata$Relate == 7)] - origdata$Age[which(origdata$Relate == 7)])) #Min age diff btw head and parent is 4

#GRANDCHILD 8
summary((origdata$HHAge[which(origdata$Relate == 8)] - origdata$Age[which(origdata$Relate == 8)])) #Min age diff btw head and grandchild is 26
pick_child_index <- origdata$HHIndex[which(origdata$Relate == 8)]
pick_spouse_index <- origdata$HHIndex[which(origdata$Relate == 2)]
keep_child_index <- pick_child_index[is.element(pick_child_index,pick_spouse_index)]
keep_spouse_index <- pick_spouse_index[is.element(pick_spouse_index,pick_child_index)]
data_child_spouse <- origdata[which(is.element(origdata$HHIndex,keep_spouse_index)),]
age_child <- data_child_spouse$Age[which(data_child_spouse$Relate==8)]
age_spouse <- data_child_spouse$Age[which(data_child_spouse$Relate==2)]
age_spouse <- rep(age_spouse,as.data.frame(table(keep_child_index))$Freq)
summary(age_spouse- age_child) #Min age diff btw spouse and head's biological child is 0
summary(age_spouse) #Min age should be 17
summary(origdata$HHAge[which(origdata$Relate == 8)]) #Min age should be 31

#PARENTINLAW 9
summary(-(origdata$HHAge[which(origdata$Relate == 9)] - origdata$Age[which(origdata$Relate == 9)])) #Min age diff btw head and parentinlaw is 4
pick_parentinlaw_index <- origdata$HHIndex[which(origdata$Relate == 9)]
pick_spouse_index <- origdata$HHIndex[which(origdata$Relate == 2)]
keep_parentinlaw_index <- pick_parentinlaw_index[is.element(pick_parentinlaw_index,pick_spouse_index)]
keep_spouse_index <- pick_spouse_index[is.element(pick_spouse_index,pick_parentinlaw_index)]
data_parentinlaw_spouse <- origdata[which(is.element(origdata$HHIndex,keep_spouse_index)),]
age_parentinlaw <- data_parentinlaw_spouse$Age[which(data_parentinlaw_spouse$Relate==9)]
age_spouse <- data_parentinlaw_spouse$Age[which(data_parentinlaw_spouse$Relate==2)]
age_spouse <- rep(age_spouse,as.data.frame(table(keep_parentinlaw_index))$Freq)
summary(age_parentinlaw - age_spouse) #Min age diff btw spouse and head's biological parentinlaw is -7

#CHILDINLAW 10
summary((origdata$HHAge[which(origdata$Relate == 10)] - origdata$Age[which(origdata$Relate == 10)])) #Min age diff btw head and adopted child is -3






