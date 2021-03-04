library(dplyr)
library(caret)
library(tidyr)


metadata_version <- "30Dec2020"
version <- "14Dec2020"
datapath <- "/home/kenneth/data/UKBB/covid19"

setwd("/home/kenneth/out/UKBB/covid19/")

# Load the input data
load(file=paste0("rerun_", metadata_version, "/df_ukbb_", version, ".RData"))

# Read the list of categorical and quantitative features with NA missing rate < 20%
selected_cat_traits <- scan(paste0(datapath, "/meta/", metadata_version, "/selected_cat_features.txt"), character(), sep="\n")
selected_cat_traits <- selected_cat_traits[selected_cat_traits != "sex"]
selected_qtl_traits <- scan(paste0(datapath, "/meta/", metadata_version, "/selected_qtl_features.txt"), character(), sep="\n")

# Read pre-processed classification of EID
inpatient.id <- scan(paste0("/home/kenneth/out/UKBB/covid19/id_", version, "/tested_case_inpatient.id"), integer(), sep="\n")
outpatient.id <- scan(paste0("/home/kenneth/out/UKBB/covid19/id_", version, "/tested_case_outpatient.id"), integer(), sep="\n")
u071.id <- scan(paste0("/home/kenneth/out/UKBB/covid19/id_", version, "/dead_U071.id"), integer(), sep="\n")

# Split dataset by Sex
df.ukbb.male <- df.ukbb.out[df.ukbb.out$sex==1,]
df.ukbb.female <- df.ukbb.out[df.ukbb.out$sex==0,]
df.ukbb.male$sex <- NULL
df.ukbb.female$sex <- NULL
rm(df.ukbb.out)


## FUNCTION

gen_cohort <- function(df.ukbb.out, selected_cat_traits, selected_qtl_traits, sex) {

df.ukbb <- df.ukbb.out
df.ukbb$eid <- rownames(df.ukbb.out)
df.test <- df.ukbb[!is.na(df.ukbb$result),]
df.test.pos <- df.test[df.test$result == 1,]
df.test.neg <- df.test[df.test$result == 0,]
df.test.pos.in <- subset(df.test.pos, eid %in% inpatient.id) 
df.test.pos.out <- subset(df.test.pos, eid %in% outpatient.id)

# Remove the resurrection subject (eid=3648154) dead on 6 May 2012 but got covid-19 on 1 May 2020
df.test.pos.out <- df.test.pos.out[df.test.pos.out$eid != 3648154,]
#resurrection_subj_idx <- which(is.na(df.test.pos.out$eid))
#df.test.pos.out <- df.test.pos.out[-resurrection_subj_idx,]

# Extract the 5 U071 subjects having -ve covid-19 test result >1 week before dead => +ve U071 case
U071.neg2pos <- df.test[df.test$eid %in% c(1986819, 2661083, 3356568, 3457484, 1819103),]

# Cohort C2 (UKBB alive subj without Covid-19)
df.nottest <- df.ukbb[is.na(df.ukbb$result),]
df.nottest.alive <- subset(df.nottest, !eid %in% u071.id)
df.test.neg.alive <- subset(df.test.neg, !eid %in% u071.id)
df.c2 <- rbind(df.nottest.alive,
	                      df.test.neg.alive)
dim(df.c2)

# cohort A1/C1 (Serious Covid-19 patients)
df.nottest.u071 <- subset(df.nottest, eid %in% u071.id)
df.test.pos.out.u071 <- subset(df.test.pos.out, eid %in% u071.id)
df.test.pos.in.u071 <- subset(df.test.pos.in, eid %in% u071.id)
df.a1 <- rbind(df.nottest.u071,
	                      df.test.pos.out.u071,
			                     df.test.pos.in,
			                     U071.neg2pos)
dim(df.a1)

# cohort A2 (Mild Covid-19 patients)
df.a2 <- subset(df.test.pos.out, !eid %in% u071.id)
dim(df.a2)

# Add label to dataset A1 and A2, and create dataset A
df.a1[, "outcome"] <- as.factor("case")
df.a2[, "outcome"] <- as.factor("ctrl")
df.a <- rbind(df.a1, df.a2)
df.a$result <- NULL
df.a$eid <- NULL
dim(df.a)

# cohort B1 (death from Covid-19 cases excluding negative covid-19 tested subject and death from Covid-19)
df.test.neg <- subset(df.test, result==0)
df.test.neg.u071 <- subset(df.test.neg, eid %in% u071.id)
df.test.pos.u071 <- subset(df.ukbb, eid %in% u071.id & !eid %in% df.test.neg.u071$eid)
df.b1 <- rbind(df.test.pos.u071, U071.neg2pos)
dim(df.b1)

# cohort B2 (Other Covid-19 tested alive cases)
df.test.pos.in.alive <- subset(df.test.pos.in, !eid %in% u071.id)
df.test.pos.out.alive <- subset(df.test.pos.out, !eid %in% u071.id)
df.test.pos.alive <- rbind(df.test.pos.in.alive, df.test.pos.out.alive)
df.b2 <- df.test.pos.alive
dim(df.b2)

# Add label to dataset B1 and B2, and create dataset B
df.b1[, "outcome"] <- "case"
df.b2[, "outcome"] <- "ctrl"
df.b <- rbind(df.b1, df.b2)
df.b$result <- NULL
df.b$eid <- NULL
dim(df.b)

# Add label to dataset C1 and C2, and create dataset C
df.c2[, "outcome"] <- "ctrl"
df.c <- rbind(df.a1, df.c2)
df.c$result <- NULL
df.c$eid <- NULL
dim(df.c)

# Cohort D
df.d <- rbind(df.b1, df.c2)
df.d$result <- NULL
df.d$eid <- NULL
dim(df.d)

# Save cohort A-D dataframe to files
save(df.a, df.b, df.c, df.d, selected_cat_traits, selected_qtl_traits, file=paste0("rerun_", metadata_version, "_", sex, "/df_cohorts_", version, ".", sex, ".RData"))
}


# Generate cohorts A-D splitted by sex
gen_cohort(df.ukbb.male, selected_cat_traits, selected_qtl_traits, "male")
gen_cohort(df.ukbb.female, selected_cat_traits, selected_qtl_traits, "female")
