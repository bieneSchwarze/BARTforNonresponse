################################################################################
################################################################################
## BART FOR PREDICTING PARTICIPATION PROBABILITIES
## 19.06.2019
## Timo Gnambs
################################################################################
################################################################################


################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################

# Clear workspace
rm(list = ls())

# Load libraries
library(haven)
library(dplyr)
library(psych)
library(MBESS)
library(lavaan)

# Load rawdata
dat1 <- read_dta("DesignData/_designDataSC3W1toW9_ID_t.dta")  # unpublished, for internal use only
dat2 <- read_sav("SC3_D_8-0-0/SC3_CohortProfile_D_8-0-0.sav")
dat3 <- read_sav("SC3_D_8-0-0/SC3_pTarget_D_8-0-0.sav")
dat4 <- read_sav("SC3_D_8-0-0/SC3_xTargetCompetencies_D_8-0-0.sav")


################################################################################
# II. RECODE DESIGN VARIABLES
################################################################################

#' Keep main sample only
table(dat1$tstud_st) 
attributes(dat1$tstud_st)$labels # A28 = main sample, A30A = refreshment sample, 
                                 # A56 = Förderschule, A63= migrant oversampling
dat <- filter(dat1, tstud_st == 1) # keep main sample

# Exclude school types: GS, SU
dat <- filter(dat, school_type %in% 3:7) # exclude GS, SU

# Keep participants at wave 1 (exclude temporary dropouts)
dat <- filter(dat, status_1 %in% "Part")

# School type
table(dat$school_type)
attributes(dat$school_type)$labels
dat$SCHOOL <- recode(unclass(dat$school_type), '3' = 'GY', '4' = 'HS', '5' = 'IG', '6' = 'MB', '7' = 'RS')
dat$SCHOOLGY <- recode(dat$SCHOOL, 'GY' = 1, .default = 0) # 1 = Gymnasium
dat$SCHOOLRS <- recode(dat$SCHOOL, 'RS' = 1, .default = 0) # 1 = Realschule
dat$SCHOOLOS <- recode(dat$SCHOOL, 'GY' = 0, 'RS' = 0, .default = 1) # 1 = other school type
table(dat$SCHOOL, useNA = "always")
table(dat$SCHOOLGY, useNA = "always")
table(dat$SCHOOLRS, useNA = "always")
table(dat$SCHOOLOS, useNA = "always")

# Sex
dat$sex <- dat$Gender - 1 # 0 = male, 1 = female
attributes(dat$Gender)$labels
table(dat$sex, useNA = "always")

# Age (in years)
dat$age <- dat$Age
describe(dat$age)

# Native language
table(dat1$NativeLanguage)
attributes(dat1$NativeLanguage)$labels
dat$lang <- recode(unclass(dat$NativeLanguage), '1' = 0, '2' = 1)
dat$lang[dat$lang == 3] <- NA
table(dat$lang, useNA = "always")

# Federal state
table(dat$State)
attributes(dat$State)$labels
dat$STATE <- recode(unclass(dat$State), '1' = "BB", '2' = 'BE', '3' = 'BW', '4' = 'BY',
                    '5' = 'HB', '6' = 'HE', '7' = 'HH', '8' = 'MV', '9' = 'NI' , '10' = 'NW',
                    '11' = 'RP', '12' = 'SH', '13' = 'SL', '14' = 'SN', '15' = 'ST', '16' = 'TH')
dat$STATEBB <- recode(dat$STATE, 'BB' = 1, .default = 0)
dat$STATEBE <- recode(dat$STATE, 'BE' = 1, .default = 0)
dat$STATEBW <- recode(dat$STATE, 'BW' = 1, .default = 0)
dat$STATEBY <- recode(dat$STATE, 'BY' = 1, .default = 0)
dat$STATEHB <- recode(dat$STATE, 'HB' = 1, .default = 0)
dat$STATEHE <- recode(dat$STATE, 'HE' = 1, .default = 0)
dat$STATEHH <- recode(dat$STATE, 'HH' = 1, .default = 0)
dat$STATEMV <- recode(dat$STATE, 'MV' = 1, .default = 0)
dat$STATENI <- recode(dat$STATE, 'NI' = 1, .default = 0)
dat$STATENW <- recode(dat$STATE, 'NW' = 1, .default = 0)
dat$STATERP <- recode(dat$STATE, 'RP' = 1, .default = 0)
dat$STATESH <- recode(dat$STATE, 'SH' = 1, .default = 0)
dat$STATESL <- recode(dat$STATE, 'SL' = 1, .default = 0)
dat$STATESN <- recode(dat$STATE, 'SN' = 1, .default = 0)
dat$STATEST <- recode(dat$STATE, 'ST' = 1, .default = 0)
dat$STATETH <- recode(dat$STATE, 'TH' = 1, .default = 0)
table(dat$STATE, useNA = "always")
apply(dat[, paste0("STATE", c("BB", "BE", "BW", "BY", "HB", "HE", "HH", "MV", "NI", 
                              "NW", "RP", "SH", "SL", "SN", "ST", "TH"))], 2, table, useNA = "always")

# Dropout indicators
table(dat$status_1) # PART = Participated, TDO = temporary dropout, FD=12 = final dropout
dat$status1 <- recode(dat$status_1, "Part" = 0, "TDO" = 1)
table(dat$status_2)
dat$status2 <- recode(dat$status_2, "Part" = 0, "TDO" = 1, .default = 2)
table(dat$status_3)
dat$status3 <- recode(dat$status_3, "Part" = 0, "TDO" = 1, .default = 2)
table(dat$status_4)
dat$status4 <- recode(dat$status_4, "Part" = 0, "TDO" = 1, .default = 2)
table(dat$status_5)
dat$status5 <- recode(dat$status_5, "Part" = 0, "TDO" = 1, .default = 2)
table(dat$status_6)
dat$status6 <- recode(dat$status_6, "Part" = 0, "TDO" = 1, .default = 2)
table(dat$status_7)
dat$status7 <- recode(dat$status_7, "Part" = 0, "TDO" = 1, .default = 2)
table(dat$status_8)
dat$status8 <- recode(dat$status_8, "Part" = 0, "TDO" = 1, "TDO_ab8" = 1, .default = 2)
table(dat$status_9)
dat$status9 <- recode(dat$status_9, "Part" = 0, "TDO" = 1, .default = 2)
# 0 = participation, 1 = temporary dropout, 2 = final dropout
apply(dat[, paste0("status", 1:9)], 2, table, useNA = "always")

# Individual field
table(dat$indNV_1) # 1 = no, 2 = yes
dat$indNV1 <- recode(unclass(dat$indNV_1), '1' = 0, '2' = 1)
table(dat$indNV_2) # 1 = no, 2 = yes
dat$indNV2 <- recode(unclass(dat$indNV_2), '1' = 0, '2' = 1)
table(dat$indNV_3) # 1 = no, 2 = yes
dat$indNV3 <- recode(unclass(dat$indNV_3), '1' = 0, '2' = 1)
table(dat$indNV_4) # 1 = no, 2 = yes
dat$indNV4 <- recode(unclass(dat$indNV_4), '1' = 0, '2' = 1)
table(dat$indNV_5) # 1 = no, 2 = yes
dat$indNV5 <- recode(unclass(dat$indNV_5), '1' = 0, '2' = 1)
table(dat$indNV_6) # 1 = no, 2 = yes
dat$indNV6 <- recode(unclass(dat$indNV_6), '1' = 0, '2' = 1)
table(dat$indNV_7) # 1 = no, 2 = yes
dat$indNV7 <- recode(unclass(dat$indNV_7), '1' = 0, '2' = 1)
table(dat$indNV_8) # 1 = no, 2 = yes
dat$indNV8 <- recode(unclass(dat$indNV_8), '1' = 0, '2' = 1)
table(dat$indNV_9) # 1 = no, 2 = yes
dat$indNV9 <- recode(unclass(dat$indNV_9), '1' = 0, '2' = 1)

# City - urban
dat$URBAN <- factor(dat$stadtland,labels = c("halbstädtisch", "ländlich", "städtisch"))
table(dat$URBAN, useNA = "always")
dat$URBANH <- ifelse(dat$URBAN == "halbstädtisch", 1, 0)
dat$URBANL <- ifelse(dat$URBAN == "ländlich", 1, 0)
table(dat$URBANH, useNA = "always")
table(dat$URBANL, useNA = "always")

# School holding
table(dat$traeger)
dat$TRAEGER <- recode(dat$traeger, "f" = 1, .default = 0) # 0 = öffentlich, 1 = privat
table(dat$TRAEGER, useNA = "always")

# School characateristics in grade 8
table(dat$anzahl_K8, useNA = "always") # Anzahl Klassen in K8
dat$CLASSES8 <- dat$anzahl_K8
table(dat$anzahl_S8, useNA = "always") # Anzahl Schüler in K8
dat$STUDENTS8 <- dat$anzahl_S8



################################################################################
# III. RECODE PERSON VARIABLES
################################################################################

# Migration background and number of books at home
dat <- left_join(select(dat,
                        ID_t, ID_i, sex, age, lang, SCHOOL, SCHOOLGY, SCHOOLRS, SCHOOLOS,
                        STATE, STATEBB, STATEBE, STATEBW, STATEBY, STATEHB, STATEHE, 
                        STATEHH, STATEMV, STATENI, STATENW, STATERP, STATESH, STATESL,
                        STATESN, STATEST, STATETH, status1, status2, status3, status4,
                        status5, status6, status7, status8, status9, indNV1, indNV2,
                        indNV3, indNV4, indNV5, indNV6, indNV7, indNV8, indNV9, URBAN,
                        URBANH, URBANL, CLASSES8, STUDENTS8, TRAEGER),
                 dat3 %>%
                   group_by(ID_t) %>%
                   summarize(mig = round(median(t400500_g1, na.rm = TRUE)),       # migration background
                             books = round(median(t34005a, na.rm = TRUE))) %>%    # number of books
                   ungroup(),
                 by = "ID_t")
dat$mig <- ifelse(dat$mig >= 1 & dat$mig <= 5, 1, 0)
table(dat$mig, useNA = "always")   # migration background
table(dat$books, useNA = "always") # number of books at home

# Competences
dat <- left_join(dat, 
                 select(dat4, ID_t, 
                        dgg5_sc3a, dgg5_sc3b,              # perceptual speed, reasoning
                        org5_sc1a, org5_sc1b, rsg5_sc3,    # orthography, reading speed
                        mag5_sc1u, reg5_sc1) %>%           # maths, reading
                 rename(speed = dgg5_sc3a, reas = dgg5_sc3b, 
                        ortho1 = org5_sc1a, ortho2 = org5_sc1b, rspeed = rsg5_sc3,
                        maths = mag5_sc1u, read = reg5_sc1),
                 by = "ID_t")
table(is.na(dat$speed))   # perceptual speed
table(is.na(dat$reas))    # reasoning
dat$ortho <- rowMeans(dat[, c("ortho1", "ortho2")])
table(is.na(dat$ortho))   # orthography
table(is.na(dat$rspeed))  # reading speed
table(is.na(dat$maths))   # mathematical competence
table(is.na(dat$read))    # reading competence

# Satisfaction, subjective health, and self-esteem
dat <- left_join(dat, 
                 filter(dat3, wave == 1) %>%
                 select(ID_t, 
                        t514001, t514002, t514003,  # life, living standards, health
                        t514004, t514005, t514006,  # family, friends, school
                        t521000, t66003a_g1) %>%    # subjective health, self-esteem
                 rename(sat1 = t514001, sat2 = t514002, sat3 = t514003, 
                        sat4 = t514004, sat5 = t514005, sat6 = t514006,
                        health = t521000, sees = t66003a_g1),
                 by = "ID_t")
table(is.na(dat$sat1))   # satisfaction with life
table(is.na(dat$sat2))   # satisfaction with current standards of living
table(is.na(dat$sat3))   # satisfaction with health
table(is.na(dat$sat4))   # satisfaction with family
table(is.na(dat$sat5))   # satisfaction with friends
table(is.na(dat$sat6))   # satisfaction with school
table(is.na(dat$health)) # subjective health
table(is.na(dat$sees))   # self-esteem

#' Family composition
dat <- left_join(dat, 
                 filter(dat3, wave == 1) %>%
                   select(ID_t, t741002) %>%       # family size
                   rename(famsize = t741002),
                 by = "ID_t")
dat$famsize[dat$famsize %in% c(0, 1, 99)] <- NA
table(dat$famsize, useNA = "always")

# School attendance
dat <- left_join(dat,
                 mutate(dat3, sickdays = t523010,
                        classrep = t725020 - 1) %>% # 0 = no, 1 = yes
                 filter(wave == 1) %>%
                   select(ID_t, sickdays, classrep),
                 by = "ID_t")
dat$sickdays[dat$sickdays > 50] <- NA
table(dat$sickdays, useNA = "always")  # number of sick days
table(dat$classrep, useNA = "always")  # class repeated (1 = yes, 0 = no)

# Self-concept and grades
dat <- left_join(dat, 
                   filter(dat3, wave == 1) %>%
                   select(ID_t, 
                          t66000a_g1, t66001a_g1, t66002a_g1,  # self-concept German/ maths / school
                          t724101, t724102) %>%                # grades German / maths
                   rename(sc1 = t66000a_g1, sc2 = t66001a_g1, sc3 = t66002a_g1, 
                          grade1 = t724101, grade2 = t724102),
                 by = "ID_t")
table(dat$sc1, useNA = "always") # self-concept German
table(dat$sc2, useNA = "always") # self-concept mathematics
table(dat$sc3, useNA = "always") # self-concept school
table(dat$grade1, useNA = "always") # grades German
table(dat$grade2, useNA = "always") # grades mathematics



################################################################################
# IV. CREATE SCHOOL VARIABLES
################################################################################

dati <- group_by(dat, ID_i) %>%
        summarize(N = n(),
                  sexM = mean(sex, na.rm = TRUE),
                  ageM = mean(age, na.rm = TRUE),
                  langM = mean(lang, na.rm = TRUE),
                  migM = mean(mig, na.rm = TRUE),
                  booksM = mean(books, na.rm = TRUE),
                  famsizeM = mean(famsize, na.rm = TRUE),
                  sickdaysM = mean(sickdays, na.rm = TRUE),
                  classrepM = mean(classrep, na.rm = TRUE),
                  grade1M = mean(grade1, na.rm= TRUE),
                  grade2M = mean(grade2, na.rm= TRUE),
                  speedM = mean(speed, na.rm = TRUE),
                  reasM = mean(reas, na.rm = TRUE),
                  orthoM = mean(ortho, na.rm = TRUE),
                  mathsM = mean(maths, na.rm = TRUE),
                  readM = mean(read, na.rm = TRUE),
                  rspeedM = mean(rspeed, na.rm = TRUE),
                  sat1M = mean(sat1, na.rm= TRUE),
                  sat2M = mean(sat2, na.rm= TRUE),
                  sat3M = mean(sat3, na.rm= TRUE),
                  sat4M = mean(sat4, na.rm= TRUE),
                  sat5M = mean(sat5, na.rm= TRUE),
                  sat6M = mean(sat6, na.rm= TRUE),
                  healthM = mean(health, na.rm= TRUE),
                  seesM = mean(sees, na.rm= TRUE),
                  sc1M = mean(sc1, na.rm= TRUE),
                  sc2M = mean(sc2, na.rm= TRUE),
                  sc3M = mean(sc3, na.rm= TRUE))
describe(dati$N)        # number of students from school in panel
describe(dati$sexM)      # percentage of female students
describe(dati$ageM)      # mean age of students
describe(dati$langM)     # percentage of students with non-German native language
describe(dati$migM)      # percentage of students with migration background
describe(dati$booksM)    # mean number of books at home
describe(dati$famsizeM)  # mean family size
describe(dati$sickdaysM) # mean number of sick days
describe(dati$classrepM) # mean number of classess repaeted
describe(dati$grade1M)   # mean grade in German
describe(dati$grade2M)   # mean grade in mathematics
describe(dati$speedM)    # mean perceptual speed
describe(dati$reasM)     # mean reasoning
describe(dati$orthoM)    # mean orthography
describe(dati$mathsM)    # mean mathematical competence
describe(dati$readM)     # mean reading competence
describe(dati$rspeedM)   # mean reading speed
describe(dati$sat1M)     # mean satisfaction with life
describe(dati$sat2M)     # mean satisfaction with current standards of living
describe(dati$sat3M)     # mean satisfaction with health
describe(dati$sat4M)     # mean satisfaction with family
describe(dati$sat5M)     # mean satisfaction with friends
describe(dati$sat6M)     # mean satisfaction with school
describe(dati$healthM)   # mean subjective health
describe(dati$seesM)     # mean self-esteem
describe(dati$sc1M)      # mean self-concept German
describe(dati$sc2M)      # mean self-concept mathematics
describe(dati$sc3M)      # mean self-concept school



################################################################################
# V. DESCRIPTIVES
################################################################################


# Sample size
nrow(dat)

# Sex
prop.table(table(dat$sex))

# Age
describe(dat$age)

# School type
prop.table(table(dat$SCHOOL))

# Descriptives

vars <- c( "sex", "age", "lang", "mig", "books", "speed", "reas", "ortho",
           "rspeed", "maths", "read", "sat1", "sat2", "sat3",
           "sat4", "sat5", "sat6", "health", "sees", "famsize", "sickdays",
           "classrep", "sc1", "sc2", "sc3", "grade1", "grade2", 
           "SCHOOLRS", "SCHOOLOS", "CLASSES8", "STUDENTS8", "TRAEGER", 
           "URBANH", "URBANL")
apply(dat[, vars], 2, function(x) {
  round(c(M = mean(x, na.rm = TRUE), 
          SD = sd(x, na.rm = TRUE), 
          Min = min(x, na.rm = TRUE), 
          Max = max(x, na.rm = TRUE)), 1)
})

vars <- c("sexM", "ageM", "migM", "booksM", "langM", "speedM", "reasM",
           "orthoM", "rspeedM", "mathsM", "readM", "sat1M", "sat2M", "sat3M",
           "sat4M", "sat5M", "sat6M", "healthM", "seesM", "famsizeM", "classrepM",
           "sc1M", "sc2M", "sc3M", "grade1M", "grade2M", "sickdaysM")
apply(dati[, vars], 2, function(x) {
  round(c(M = mean(x, na.rm = TRUE), 
          SD = sd(x, na.rm = TRUE), 
          Min = min(x, na.rm = TRUE), 
          Max = max(x, na.rm = TRUE)), 1)
})

#' **Reliability of RSES**
t <- filter(dat3, ID_t %in% datl$ID_t & dat3$wave == 1) %>%
  select(t66003a, t66003b, t66003c,
         t66003d, t66003e, t66003f,
         t66003g, t66003h, t66003i,
         t66003j)
t <- na.omit(t)
t$t66003b <- 6 - t$t66003b
t$t66003e <- 6 - t$t66003e
t$t66003f <- 6 - t$t66003f
t$t66003h <- 6 - t$t66003h
t$t66003i <- 6 - t$t66003i
ci.reliability(as.matrix(t), interval.type = "none", type="categorical")$est
rm(t)



################################################################################
# VI. HANDLE MISSING VALUES
################################################################################


vars <- c( "sex", "age", "lang", "mig", "books", "speed", "reas", "ortho",
           "rspeed", "maths", "read", "sat1", "sat2", "sat3",
           "sat4", "sat5", "sat6", "health", "sees", "famsize", "sickdays",
           "classrep", "sc1", "sc2", "sc3", "grade1", "grade2", 
           "SCHOOLRS", "SCHOOLOS", "CLASSES8", "STUDENTS8", "TRAEGER", 
           "URBANH", "URBANL")

# Percentage of missing values
apply(dat[, vars], 2, function(x) { round(prop.table(table(is.na(x))), 2) })

# Replace missing values with mean if less than 5%
dat$mig[is.na(dat$mig)] <- median(dat$mig, na.rm = TRUE)
dat$books[is.na(dat$books)] <- median(dat$books, na.rm = TRUE)
dat$speed[is.na(dat$speed)] <- median(dat$speed, na.rm = TRUE)
dat$reas[is.na(dat$reas)] <- median(dat$reas, na.rm = TRUE)
dat$maths[is.na(dat$maths)] <- mean(dat$maths, na.rm = TRUE)
dat$read[is.na(dat$read)] <- mean(dat$read, na.rm = TRUE)
dat$sat1[is.na(dat$sat1)] <- median(dat$sat1, na.rm = TRUE)
dat$sat2[is.na(dat$sat2)] <- median(dat$sat2, na.rm = TRUE)
dat$sat3[is.na(dat$sat3)] <- median(dat$sat3, na.rm = TRUE)
dat$sat4[is.na(dat$sat4)] <- median(dat$sat4, na.rm = TRUE)
dat$sat5[is.na(dat$sat5)] <- median(dat$sat5, na.rm = TRUE)
dat$sat6[is.na(dat$sat6)] <- median(dat$sat6, na.rm = TRUE)
dat$health[is.na(dat$health)] <- median(dat$health, na.rm = TRUE)
dat$famsize[is.na(dat$famsize)] <- median(dat$famsize, na.rm = TRUE)
dat$classrep[is.na(dat$classrep)] <- median(dat$classrep, na.rm = TRUE)
dat$CLASSES8[is.na(dat$CLASSES8)] <- median(dat$CLASSES8, na.rm = TRUE)
dat$STUDENTS8[is.na(dat$STUDENTS8)] <- median(dat$STUDENTS8, na.rm = TRUE)

# Create missing indicator for missing more than 5%
dat$seesNA <- is.na(dat$sees)
dat$sees[is.na(dat$sees)] <- median(dat$sees, na.rm = TRUE)
dat$grade1NA <- is.na(dat$grade1)
dat$grade1[is.na(dat$grade1)] <- median(dat$grade1, na.rm = TRUE)
dat$grade2NA <- is.na(dat$grade2)
dat$grade2[is.na(dat$grade2)] <- median(dat$grade2, na.rm = TRUE)
dat$sickdaysNA <- is.na(dat$sickdays)
dat$sickdays[is.na(dat$sickdays)] <- median(dat$sickdays, na.rm = TRUE)
dat$sc1NA <- is.na(dat$sc1)
dat$sc1[is.na(dat$sc1)] <- median(dat$sc1, na.rm = TRUE)
dat$sc2NA <- is.na(dat$sc2)
dat$sc2[is.na(dat$sc2)] <- median(dat$sc2, na.rm = TRUE)
dat$sc3NA <- is.na(dat$sc3)
dat$sc3[is.na(dat$sc3)] <- median(dat$sc3, na.rm = TRUE)

# Percentage of missing values
apply(dat[, vars], 2, function(x) { prop.table(table(is.na(x))) })
apply(dat[, vars], 2, function(x) { length(unique(x)) })



################################################################################
# VII. TRANSFORM TO LONG FORMAT
################################################################################


# Respondent-level data
dat <- left_join(dat, dati, by = "ID_i")

# Convert respondent data to long format
datl <- reshape(as.data.frame(dat), direction = "long", idvar = "ID_t", 
                varying = list(paste0("status", 1:9), paste0("indNV", 1:9)), sep = "",
                times = 1:9, timevar = "wave", v.names = c("status", "indNV"))
datl <- datl[order(datl$ID_t), ]
table(datl$status, datl$wave, useNA = "always")
table(datl$indNV, datl$wave, useNA = "always")

# Select variables for analyses
names(datl)
vars <- c("ID_t", "sex", "age", "mig", "books", "lang", "speed", "reas",
                  "ortho", "rspeed", "maths", "read", "sat1", "sat2", "sat3",
                  "sat4", "sat5", "sat6", "health", "sees", "famsize", "classrep",
                  "sc1", "sc2", "sc3", "grade1", "grade2", "sickdays",
                  "seesNA", "sc1NA", "sc2NA", "sc3NA", "grade1NA", "grade2NA",
                  "sexM", "ageM", "migM", "booksM", "langM", "speedM", "reasM",
                  "orthoM", "rspeedM", "mathsM", "readM", "sat1M", "sat2M", "sat3M",
                  "sat4M", "sat5M", "sat6M", "healthM", "seesM", "famsizeM", "classrepM",
                  "sc1M", "sc2M", "sc3M", "grade1M", "grade2M", "sickdaysM",
                  "SCHOOLRS", "SCHOOLOS", "CLASSES8", "STUDENTS8", "TRAEGER", 
                  "URBANH", "URBANL", "STATEBB", "STATEBE", "STATEBW", "STATEBY",
                  "STATEHB", "STATEHE", "STATEHH", "STATEMV", "STATENI", "STATENW",
                  "STATERP", "STATESH", "STATESL", "STATESN", "STATEST", "STATETH")
datl <- select(datl, one_of(vars), wave, status, indNV)
table(datl$wave, datl$indNV, useNA = "always")
apply(datl[, vars], 2, function(x) { prop.table(table(is.na(x))) })
apply(datl[, vars], 2, function(x) { length(unique(x)) })



################################################################################
# VIII. SAVE DATA
################################################################################


save(dat, datl, file = "data.Rdata")



################################################################################
# IX. FURTHER DESCRIPTIVES
################################################################################


# Sample size
nrow(dat)

# Sex
prop.table(table(dat$sex))

# Age
describe(dat$age)

# School type
prop.table(table(dat$SCHOOL))

# Variables
colnames(datl)

# Descriptives
apply(datl[, -c(1, 85:87)], 2, function(x) {
    round(c(M = mean(x), SD = sd(x), Min = min(x), Max = max(x)), 2)
})

# Reliability of RSES
t <- filter(dat3, ID_t %in% datl$ID_t & dat3$wave == 1) %>%
  select(t66003a, t66003b, t66003c,
         t66003d, t66003e, t66003f,
         t66003g, t66003h, t66003i,
         t66003j)
t <- na.omit(t)
t$t66003b <- 6 - t$t66003b
t$t66003e <- 6 - t$t66003e
t$t66003f <- 6 - t$t66003f
t$t66003h <- 6 - t$t66003h
t$t66003i <- 6 - t$t66003i
ci.reliability(as.matrix(t), interval.type = "none", type="categorical")$est
rm(t)

# Reliability of self-concept
t <- filter(dat3, ID_t %in% datl$ID_t & dat3$wave == 1) %>%
  select(t66000a, t66000b, t66000c,
         t66001a, t66001b, t66001c,
         t66002a, t66002b, t66002c)
t <- na.omit(t)
t$t66000a <- 6 - t$t66000a
ci.reliability(as.matrix(t)[, 1:3], interval.type = "none", type="categorical")$est # German
ci.reliability(as.matrix(t)[, 4:6], interval.type = "none", type="categorical")$est # maths
ci.reliability(as.matrix(t)[, 7:9], interval.type = "none", type="categorical")$est # school
rm(t)


# Reliability of reasoning
t <- filter(dat4, ID_t %in% datl$ID_t) %>%
  select(dgci2101_sc3g5_c, dgci2102_sc3g5_c, dgci2103_sc3g5_c, dgci2104_sc3g5_c,
         dgci2201_sc3g5_c, dgci2202_sc3g5_c, dgci2203_sc3g5_c, dgci2204_sc3g5_c,
         dgci2301_sc3g5_c, dgci2302_sc3g5_c, dgci2303_sc3g5_c, dgci2304_sc3g5_c)
t <- na.omit(t)
ci.reliability(as.matrix(t), interval.type = "none", type="categorical")$est
rm(t)

