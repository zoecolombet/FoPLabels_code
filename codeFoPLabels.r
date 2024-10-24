##doing mc monte carlo
library(data.table)
library(readxl)
library(fst)
library(gamlss)
library(demography)
#library(xlsx)
library(foreach)
library(qs) # save R objects to disk

setDTthreads(1L, restore_after_fork = FALSE)
threads_fst(1L, reset_after_fork = FALSE)

setwd("Warning labels model UK")

#cREATION OF THE POPULATION
# we only need year 2012 (init year) and ages 30-90
# Set simulation parameters
init_year <- 12L     # initial year
sim_hor <- 32L       # simulation horizon
impl_year <- 24L     # year of the implementation of the policy
min_age <- 30L       # min age for the population / need to be the age of interest (especially regarding the RR) minus lagtime. For us here, RR start a 40, we have 6y lag time, we might want to start at 34. For comparibility with UK model, lets keep 30
max_age <- 89L       # max age for the population
#should be 94  // issue with 90 as for IMD its 90+
SES <- 1:5           #SES IMD quintiles
scale_factor <- 500 # factor to downscale the country population #was 1000 before


#database needed
#population
onspop_import <- readRDS("./Population_Mortality/onspop.rds")
#mortality
onsmort_import <- readRDS("./Population_Mortality/onsmort_new.rds")
#GAMLSS bmi
tbl_bmi_withoutNRJtable <- read_fst("./Exposures/bmi_withoutNRJtable_quicker_CK.fst", as.data.table = TRUE)
#I change the GAMLSS for now as it was too long to run. REDO GAMLSS BETTER ASAP
#tbl_bmi_withoutNRJtable <- read_fst("./bmi_withoutNRJtable.fst", as.data.table = TRUE)
qread("./Exposures/bmi_withoutNRJmodelfinal_quicker.qs")$family
#TO DO
#need to change the distribution name within the code cerca line 135 (qSEP4)


#Energy potentially impacted by the policy
#how much of your daily energy is coming from product potentially labelled?
#1. percentage of energy from product from shopping in supermarkets/markets
#According to 2019 figures on household consumption expenditure, 55% of all food and beverage expenditure (including alcoholic beverages)
#was for at-home consumption compared to 45 percent spent on restaurants and other out-of-home food services.

#on the food/drink purchased, how much is eligible to label?
#https://www.theguardian.com/society/2017/apr/07/uk-eats-almost-four-times-more-packaged-food-than-fresh#comments
# around 400kcal are from fresh (not labelled) on the 2050kcal purchased per capita per day: 20% are fresh/80% are packaged
#energy_purchased_packaged <- Energykcal * 0.55 * 0.80
#We use NDNS to make estimation of energy_purchased_packaged
#and we did GAMLSS for it
#tbl_energy_purchased_packaged_table <- read_fst("./Exposures/energy_purchased_packaged_faster_table_light.fst", as.data.table = TRUE)
tbl_energy_purchased_packaged_table <- read_fst("./Exposures/energy_purchased_packaged_faster_table.fst", as.data.table = TRUE)
# tbl_energy_purchased_packaged_table_light <- tbl_energy_purchased_packaged_table[year>=init_year & year<=init_year+sim_hor & age>=min_age,]
# write_fst(tbl_energy_purchased_packaged_table_light, path = "./Exposures/energy_purchased_packaged_faster_table_light.fst", compress = 100L)
# saveRDS(tbl_energy_purchased_packaged_table_light, "./Exposures/tbl_energy_purchased_packaged_table_SERVER.rds")
qread("./Exposures/energy_purchased_packaged_modelfinal_faster.qs")$family
#TO DO
#need to change the distribution name within the code cerca line 166 (qBCTo for now)


###In this model, we have for now 13 scenarios
#Thus, we have 260 results (13 scenarios, 5 SES, 4 diseases (chd,haem,isch,obesity) = 260 possibilities) + (1 baseline, 5 SES, 4 diseases = 20)
poss_tot <- 20*5*4 + 20 # CKnote: if 13 scenarios why do we multiply by 19? #RE: there are now 19 #now 20


#######################################MONTECARLO
#Number of iterations (for testing our model, put 1:2)
montecarlo <- 100

library(bigstatsr)
mat3 <- FBM(montecarlo, (poss_tot*sim_hor))
cl <- parallel::makeCluster(4, outfile="") ##HERE you might want to change the number of cluster
#detectCores(all.tests = FALSE, logical = TRUE)
doParallel::registerDoParallel(cl)
time.begin <- Sys.time()


system.time(
  tmp <- foreach(sim = 1:montecarlo, .combine = 'rbind', .packages=c("data.table","fst","gamlss")) %dopar% {
    set.seed(sim) # CKNote: set seed for reproducibility
    setDTthreads(1L, restore_after_fork = FALSE)
    threads_fst(1L, reset_after_fork = FALSE)

    time.old <- Sys.time()

    # Import ONS population estimates
    onspop <- onspop_import

    sp <- onspop[between(as.integer(Age), min_age, max_age),
                 .(year = init_year, age = as.integer(Age), sex = factor(Sex), SES = factor(SES), popsize = `2022`)]
    head(sp)
    # We scale down the population by the scale factor. We will scale up by factor the results (for now ignoring the truncation error)
    sp[, popsize := round(popsize/scale_factor)]
    nrow(sp)
    head(sp)

    # Let's create the simulants (person-year table)
    sp <- sp[rep(1:.N, times = popsize), .(pid = .I, year, age, sex, SES)] # rep(1:3, times = 1:3)
    nrow(sp)
    head(sp)

    # Now lets project their life course for the next 20 years (simulation horizon)
    sp <- rbindlist(rep(list(sp), sim_hor), idcol = "time")
    sp[, `:=` (year = year + time - 1L,
               age = age + time - 1L,
               time = NULL)]
    sp[, table(age, year)]   # this is the close cohort
    sp <- sp[age <= max_age] # Some pruning
    setkey(sp, pid, year)

    # The structure above is very useful in dynamic microsimulation, each row is a person year

    # We will create an open cohort
    # With our previous approach, we have no one aged 20 in 2021 or 20-21 in 22
    # We need a new cohort of 20yo simulants enter the model every year
    # Let's get back to the ONS pop size estimates file
    onspop <- melt(onspop, id.vars = c("Age", "Sex","SES"),
                   variable.name = "year", value.name = "popsize")
    setDT(onspop)
    onspop[, lapply(.SD, class)]
    onspop[, year := as.integer(as.character(year)) - 2000L]

    # Let's consider the cohort of 30yo entering the model in years post 2022
    sp2 <- onspop[as.integer(Age) == min_age &
                    between(year, init_year + 1L, init_year + sim_hor - 1L),
                  .(year, age = as.integer(Age), sex = Sex, SES=factor(SES), popsize = round(popsize / scale_factor))] #was 1000 before

    # Form the person-year table
    sp2 <- sp2[rep(1:.N, times = popsize), .(pid = .I + max(sp$pid), year, age, sex, SES)]

    # Create their lifecourse as with the original cohort
    sp2 <- rbindlist(rep(list(sp2), sim_hor), idcol = "time")
    sp2[, `:=` (year = year + time - 1L,
                age = age + time - 1L,
                time = NULL)]
    setkey(sp2, pid, year)

    # some pruning is necessary
    sp2[, table(age, year)]
    sp2 <- sp2[year < init_year + sim_hor ]
    sp2[, table(age, year)]

    # And we bind the 2 cohorts
    sp <- rbind(sp, sp2)
    setkey(sp, pid, year)
    sp[, table(age, year)]

    # Now I can start adding attributes for my simulants. The order we add the parameters
    # reflects the causality network assumptions.
    # Here I just estimate the bmi I will make the rank stability assumption
    tt <- sp[, .("pid" = unique(pid))]
    #set.seed(42L)
    tt[, rn_bmi := runif(.N) * 0.995] #test here to find good number to avoid long tail
    sp[tt, on = "pid", rn_bmi := rn_bmi]
    head(sp)
    tbl <- tbl_bmi_withoutNRJtable
    head(tbl)
    sp[tbl, on = c("year","age","sex","SES"), `:=` (mu = i.mu, sigma = i.sigma, nu = i.nu, tau = i.tau)]
    head(sp)
    sp[, bmi := qSHASH(rn_bmi, mu, sigma, nu, tau)] #need to use the quantile function here! qXXX. Also CHANGE the distribution accroding to the GAMLSS fitting
    #plot(density(sp$bmi))
    #summary(sp$bmi) #14 to 50, mean at 28, seems high CHRIS
    sp[, c("rn_bmi", "mu", "sigma", "nu", "tau") := NULL]
    head(sp)
    #BMI is the best way of doing that, instead of height and weight as modelling weight with height is super complicated
    #Here I need to do a hard boundary as I have a specific range for Energy
    sp[bmi > 55, bmi := 55]
    sp[bmi < 13, bmi := 13]

    sp[, `:=` (bmi_grp="<20")]
    sp[(bmi>=20 & bmi<25), `:=` (bmi_grp="20-25")]
    sp[(bmi>=25 & bmi<30), `:=` (bmi_grp="25-30")]
    sp[(bmi>=30), `:=` (bmi_grp=">=30")]
    head(sp)
    sp[year==impl_year, prop.table(table(bmi_grp))]
    # bmi_grp
    # <20             >=30      20-25      25-30
    # 0.03162434 0.28629612 0.28222329 0.39985625

    #having energy intakes, also predicted by the bmi
    tt <- sp[, .("pid" = unique(pid))]
    tt[, rn_nrj := runif(.N) * 0.997] #test here to find good number to avoid long tail
    sp[tt, on = "pid", rn_nrj := rn_nrj]
    head(sp)
    tbl <- tbl_energy_purchased_packaged_table
    tbl[,`:=` (sex=ifelse(sex=="Men","1","2"), bmi2=bmi)]
    head(tbl)
    #sp[,`:=` (sex=factor(sex))]
    ###ASK CHRIS IF PUTTING BMI AS AN INTEGER IS THE BEST THING TO DO
    #no, we should have integer just for the matching, so create another variable with integer that we delete after the merging
    sp[,`:=` (bmi2=round(bmi))]
    sp[tbl, on = c("year","age","sex","SES","bmi2"), `:=` (mu = i.mu, sigma = i.sigma, nu = i.nu, tau = i.tau)]
    head(sp)
    sp[, energy_purchased_packaged := qBCTo(rn_nrj, mu, sigma, nu, tau)] #need to use the quantile function here! qXXX. Also CHANGE the distribution according to the GAMLSS fitting
    #plot(density(sp$energy_purchased_packaged))
    #summary(sp$energy_purchased_packaged)
    #NDNS
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # 79.48  629.46  770.24  795.85  941.87 2340.70
    sp[, c("rn_nrj", "mu", "sigma", "nu", "tau", "bmi2") := NULL]
    head(sp)

    # Now I will model mortality dependant on current bmi
    #Artilce ERFC 2011 - Separate and combined associations of body-mass index and abdominal adiposity with cardiovascular disease: collaborative analysis of 58 prospective studies
    # https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(11)60105-0/fulltext#secd19223949e1020
    # Dans le doc https://els-jbs-prod-cdn.jbs.elsevierhealth.com/cms/attachment/6d6bb978-0d17-4017-9845-2903d3b68778/gr2.gif
    # In people with BMI of 20 kg/m² or higher: HRs per 1 SD higher baseline values : 4·56 kg/m² higher BMI
    # CHD
    # 40-59 years : 1.41 (1.30-1.53)
    # 60-69 years : 1.23 (1.15-1.31)
    # 70+ years   : 1.12 (1.05-1.19)
    #
    # Ischaemic stroke
    # 40-59 years : 1.34 (1.21-1.48)
    # 60-69 years : 1.22 (1.13-1.31)
    # 70+ years   : 1.08 (0.99-1.18)

    sp[, `:=` (rr.obesity.chd_e =1, rr.obesity.chd_l =1, rr.obesity.chd_u= 1,
               rr.obesity.isch_e=1, rr.obesity.isch_l=1, rr.obesity.isch_u=1,
               rr.obesity.haem_e=1, rr.obesity.haem_l=1, rr.obesity.haem_u=1)]
    sp[(age>=40 & age<60), `:=` (rr.obesity.chd_e= 1.41, rr.obesity.chd_l= 1.30, rr.obesity.chd_u= 1.53,
                                 rr.obesity.isch_e=1.34, rr.obesity.isch_l=1.21, rr.obesity.isch_u=1.48,
                                 rr.obesity.haem_e=1,    rr.obesity.haem_l=1,    rr.obesity.haem_u=1)]
    sp[(age>=60 & age<70), `:=` (rr.obesity.chd_e= 1.23, rr.obesity.chd_l= 1.15, rr.obesity.chd_u= 1.31,
                                 rr.obesity.isch_e=1.22, rr.obesity.isch_l=1.13, rr.obesity.isch_u=1.31,
                                 rr.obesity.haem_e=1,    rr.obesity.haem_l=1,    rr.obesity.haem_u=1)]
    sp[(age>=70), `:=` (rr.obesity.chd_e= 1.12, rr.obesity.chd_l= 1.05, rr.obesity.chd_u= 1.19,
                        rr.obesity.isch_e=1.08, rr.obesity.isch_l=0.99, rr.obesity.isch_u=1.18,
                        rr.obesity.haem_e=1,    rr.obesity.haem_l=1,    rr.obesity.haem_u=1)]
    sp[, table(rr.obesity.chd_e)]
    #sp[, table(rr.obesity.chd_e,age)]
    sp[, table(rr.obesity.isch_e)]
    sp[, `:=` (rr.obesity.chd_se=((log(rr.obesity.chd_u) - log(rr.obesity.chd_e))/1.96 + (log(rr.obesity.chd_l) - log(rr.obesity.chd_e))/-1.96)/2,
               rr.obesity.isch_se=((log(rr.obesity.isch_u) - log(rr.obesity.isch_e))/1.96 + (log(rr.obesity.isch_l) - log(rr.obesity.isch_e))/-1.96)/2,
               rr.obesity.haem_se=((log(rr.obesity.haem_u) - log(rr.obesity.haem_e))/1.96 + (log(rr.obesity.haem_l) - log(rr.obesity.haem_e))/-1.96)/2)]
    head(sp)


    # CKNote: This is not ideal as the RR should be the same for all individuals of the same age and change every iteration
    # Simulation: RR
    # sp[, rn_rr := runif(.N)] # CKChange
    sp[, rn_rr := runif(1)] # CKChange, seed was set at the top as sim
    sp[, `:=` (rr.obesity.chd=qlnorm(rn_rr,meanlog = log(rr.obesity.chd_e), sdlog = rr.obesity.chd_se))]
    sp[rr.obesity.chd < 1, rr.obesity.chd := 1] # CKChange, to avoid values < 1
    sp[, rn_rr := runif(1)] # CKChange, seed was set at the top as sim
    sp[, `:=` (rr.obesity.isch=qlnorm(rn_rr,meanlog = log(rr.obesity.isch_e), sdlog = rr.obesity.isch_se))]
    sp[rr.obesity.isch < 1, rr.obesity.isch := 1] # CKChange, to avoid values < 1

    sp[, rn_rr := runif(1)] # CKChange, seed was set at the top as sim
    sp[, `:=` (rr.obesity.haem=qlnorm(rn_rr,meanlog = log(rr.obesity.haem_e), sdlog = rr.obesity.haem_se))]
    sp[rr.obesity.haem < 1, rr.obesity.haem := 1] # CKChange, to avoid values < 1

    sp[, c("rr.obesity.chd_e", "rr.obesity.chd_l", "rr.obesity.chd_u", "rr.obesity.chd_se",
           "rr.obesity.isch_e", "rr.obesity.isch_l", "rr.obesity.isch_u", "rr.obesity.isch_se",
           "rr.obesity.haem_e", "rr.obesity.haem_l","rr.obesity.haem_u", "rr.obesity.haem_se") := NULL]
    #View(sp)

    #We are assuming that we have a 6y lag time between the bmi and the RR (ERFC 2011: 5.7 years [SD 3.0–9.0]))
    #to be easier here, we are creating a bmi_lagged to obtain the rr
    #Calculating the RR according to BMI
    setkey(sp, pid, year) # NECESSARY!!!
    sp[, bmi_lagged := shift(bmi, 6), by = pid]
    #check <- sp[,.(pid,year,bmi,bmi_lagged)]

    sp[, `:=` (rr.obesity.chd.indiv=rr.obesity.chd^((bmi_lagged-20)/4.56),  #BMIss bring back to 4.56 to be comparative with rr.Obesity
               rr.obesity.isch.indiv=rr.obesity.isch^((bmi_lagged-20)/4.56),
               rr.obesity.haem.indiv=rr.obesity.haem^((bmi_lagged-20)/4.56))]
    sp[bmi_lagged<20, `:=` (rr.obesity.chd.indiv=1,
                            rr.obesity.isch.indiv=1,
                            rr.obesity.haem.indiv=1)]


    mrtlparf1 <-
      sp[, .(parf_obesity_chd = 1 - 1 / (sum(rr.obesity.chd.indiv) / .N )), keyby = .(year, age, sex, SES)] #dont need bmi grp here/ if strong trend no year to avoid double counting, for us Chris say keep year
    head(mrtlparf1)
    mrtlparf2 <-
      sp[, .(parf_obesity_isch = 1 - 1 / (sum(rr.obesity.isch.indiv) / .N )), keyby = .(year, age, sex, SES)]
    mrtlparf3 <-
      sp[, .(parf_obesity_haem = 1 - 1 / (sum(rr.obesity.haem.indiv) / .N )), keyby = .(year, age, sex, SES)]

    mrtlparf <- mrtlparf1[mrtlparf2, on = .NATURAL,]
    mrtlparf <- mrtlparf[mrtlparf3, on = .NATURAL,]
    head(mrtlparf)
    rm(mrtlparf1,mrtlparf2,mrtlparf3)

    # Get mortality rates from ons
    onsmort <- onsmort_import

    # mx = central rate of mortality for 1
    k=runif(1)
    onsmort <- onsmort[, `:=` (cvd.rate=qlnorm(k,meanlog = log(mx_total), sdlog = mx_total_se),
                                      chd.rate=qlnorm(k,meanlog = log(mx_chd), sdlog = mx_chd_se),
                                      stroke.rate=qlnorm(k,meanlog = log(mx_stroke), sdlog = mx_stroke_se),
                                      isch.rate=qlnorm(k,meanlog = log(mx_isch), sdlog = mx_isch_se),
                                      haem.rate=qlnorm(k,meanlog = log(mx_haem), sdlog = mx_haem_se))]


    onsmort <- onsmort[,c("year","Sex","Age","AgeGroups","age","SES","popsize","Country","chd.rate","isch.rate","haem.rate")]

    yeye <- onsmort[between(as.integer(age), min_age, max_age) &
                      between(year, init_year, init_year + sim_hor - 1L),
                    .(year, age = as.integer(age), sex = factor(Sex), SES = factor(SES),
                      #here you should divide by 100,000 if the rate is for 100,000. mr = mr/1e5
                      chd.rate,isch.rate,haem.rate)
    ]


    mrtlparf[yeye, on = c("year","age","sex","SES"), `:=` (chd.rate = i.chd.rate, isch.rate = i.isch.rate,  haem.rate = i.haem.rate)]


    mrtlparf[, `:=` (p0_bmi_mrtl_chd  = chd.rate  * (1 - parf_obesity_chd),# p0 is the mortality rate of the unexposed (bmi<20)
                     p0_bmi_mrtl_isch = isch.rate * (1 - parf_obesity_isch),
                     p0_bmi_mrtl_haem = haem.rate * (1 - parf_obesity_haem))]
    head(mrtlparf)
    sp[mrtlparf, on = c("year","age","sex","SES"), `:=` (
      # chd.rate = i.chd.rate,  isch.rate = i.isch.rate, haem.rate = i.haem.rate,
      p0_bmi_mrtl_chd = i.p0_bmi_mrtl_chd, p0_bmi_mrtl_isch = i.p0_bmi_mrtl_isch,
      p0_bmi_mrtl_haem = i.p0_bmi_mrtl_haem)]

    #should I rename mr_chd > mr_bmi_chd????
    sp[, `:=` (mr_chd  = p0_bmi_mrtl_chd  * rr.obesity.chd.indiv,
               mr_isch = p0_bmi_mrtl_isch * rr.obesity.isch.indiv,
               mr_haem = p0_bmi_mrtl_haem * rr.obesity.haem.indiv)]

    # CKTest: PASSED
    # sp[year >= 18, .(sum(chd.rate), sum(mr_chd))]
    # sp[year >= 18, .(sum(isch.rate), sum(mr_isch))]
    # sp[year >= 18, .(sum(haem.rate), sum(mr_haem))]
    # sp[, c("chd.rate", "isch.rate", "haem.rate") := NULL]


    #TEST <- sp[mr_chd!=p0_bmi_mrtl_chd,]

    # Let's test it
    # sp[mrtlparf, on = c("year", "age",  "sex", "SES"), `:=` (mr_chd_ons = i.chd.rate, mr_isch_ons = i.isch.rate, mr_haem_ons = i.haem.rate)]
    # sp[, .("Our deaths" = sum(mr_chd), "ONS deaths" = sum(mr_chd_ons)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_isch), "ONS deaths" = sum(mr_isch_ons)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_haem), "ONS deaths" = sum(mr_haem_ons)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_chd+mr_isch+mr_haem), "ONS deaths" = sum(mr_chd_ons+mr_isch_ons+mr_haem_ons)), keyby = year] # our expected number of deaths by year vs the ONS one
    # sp[, .("Our deaths" = sum(mr_chd+mr_isch+mr_haem), "ONS deaths" = sum(mr_chd_ons+mr_isch_ons+mr_haem_ons))] # our expected number of deaths by year vs the ONS one
    # mrtlparf[, .("Our deaths" = sum(chd.rate+isch.rate+haem.rate)), keyby = year]
    # mrtlparf[, .("Our deaths" = sum(chd.rate+isch.rate+haem.rate))]


    # Let's see who dies
    sp[, rn_mrtl := runif(.N)] # Keep this as well
    pid_toremove <- integer(0) # a place holder
    for (i in sort(unique(sp$year))) {
      #print(paste0("Year: ", i))

      # rebalance the mortality rates to account for simulants removed from the
      # population when dead
      if (length(pid_toremove) > 0L) {
        correction_factor <- sp[year == i, sum(mr_chd)]/sp[year == i & !pid %in% pid_toremove, sum(mr_chd)]
        sp[year == i, mr_chd := mr_chd * correction_factor]
        #print(paste0("correction factor: ", correction_factor))
      }

      sp[year == i & !pid %in% pid_toremove, dead_bmi_chd := rn_mrtl < mr_chd]
      pid_toremove <- sp[year <= i & (dead_bmi_chd), pid]
      #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
    }
    sp[, longdead_bmi_chd := cumsum(dead_bmi_chd), by = pid]
    sp[dead_bmi_chd == FALSE & longdead_bmi_chd == 1, dead_bmi_chd := NA] # Best for later calculations.
    # So dead == FALSE mean alive, dead == TRUE means died that year,
    # and dead == NA means died in a previous year
    #dcast(sp, year~dead_bmi_chd) # Note NA some years is decreasing because simulants get removed from the dataset because their age is > 89
    #table(sp$dead_bmi_chd)
    #table(sp$dead_bmi_chd)
    rm(correction_factor)

    #CKNote for 2nd COD do not reset pid_toremove
    for (i in sort(unique(sp$year))) {
      #print(paste0("Year: ", i))

      # rebalance the mortality rates to account for simulants removed from the
      # population when dead
      if (length(pid_toremove) > 0L) {
        correction_factor <- sp[year == i, sum(mr_isch)]/sp[year == i & !pid %in% pid_toremove, sum(mr_isch)]
        sp[year == i, mr_isch := mr_isch * correction_factor]
        #print(paste0("correction factor: ", correction_factor))
      }

      sp[year == i & !pid %in% pid_toremove, dead_bmi_isch := rn_mrtl < mr_isch]
      pid_toremove <- sp[year <= i & (dead_bmi_isch), pid]
      #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
    }
    sp[, longdead_bmi_isch := cumsum(dead_bmi_isch), by = pid]
    sp[dead_bmi_isch == FALSE & longdead_bmi_isch == 1, dead_bmi_isch := NA] # Best for later calculations.
    # So dead == FALSE mean alive, dead == TRUE means died that year,
    # and dead == NA means died in a previous year
    #dcast(sp, year~dead_bmi_isch) # Note NA some years is decreasing because simulants get removed from the dataset because their age is > 89
    rm(correction_factor)

    #CKNote for 3rd COD do not reset pid_toremove
    for (i in sort(unique(sp$year))) {
      #print(paste0("Year: ", i))

      # rebalance the mortality rates to account for simulants removed from the
      # population when dead
      if (length(pid_toremove) > 0L) {
        correction_factor <- sp[year == i, sum(mr_haem)]/sp[year == i & !pid %in% pid_toremove, sum(mr_haem)]
        sp[year == i, mr_haem := mr_haem * correction_factor]
        #print(paste0("correction factor: ", correction_factor))
      }

      sp[year == i & !pid %in% pid_toremove, dead_bmi_haem := rn_mrtl < mr_haem]
      pid_toremove <- sp[year <= i & (dead_bmi_haem), pid]
      #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
    }
    sp[, longdead_bmi_haem := cumsum(dead_bmi_haem), by = pid]
    sp[dead_bmi_haem == FALSE & longdead_bmi_haem == 1, dead_bmi_haem := NA] # Best for later calculations.
    # So dead == FALSE mean alive, dead == TRUE means died that year,
    # and dead == NA means died in a previous year
    #dcast(sp, year~dead_bmi_haem) # Note NA some years is decreasing because simulants get removed from the dataset because their age is > 24
    rm(correction_factor)

    ##############################################################################
    #                      DEFINING ALL THE MODEL SCENARIOS                      #
    ##############################################################################
    # S0. no policy
    ##############################################################################
    ##############################################################################
    ## WARNING LABELLING IN SUPERMARKET
    ##############################################################################
    ##############################################################################
    #Let's see the effect of warning labelling

    #####EFFECT ON CONSUMER
    #Impact of nutrient warning labels on energy purchase
    #https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1003765#abstract0
    #Song 2021 - meta-analysis found change in total energy purchased
    # traffic light labelling system (TLS) vs control −6.5% [-11%;-2%];
    # Nutri-Score (NS) vs control , −6% [-11%,-1%];
    # nutrient warning (NW) vs control, −12.9% [-18%; -8%]
    #seems to be the same effect in all SES, age, sex (BMI not studied)
    sp[, policy.energy.effect.TLS := rnorm(1,mean = -0.065, sd = (-0.02-(-0.11))/3.92)]
    sp[policy.energy.effect.TLS > 0, policy.energy.effect.TLS := 0]
    sp[, policy.energy.effect.NS := rnorm(1,mean = -0.06,  sd = (-0.01-(-0.11))/3.92)]
    sp[policy.energy.effect.NS > 0, policy.energy.effect.NS := 0]
    sp[, policy.energy.effect.NW := rnorm(1,mean = -0.129, sd = (-0.08-(-0.18))/3.92)]
    sp[policy.energy.effect.NW > 0, policy.energy.effect.NW := 0]

    #if we are putting NW instead of TLS, Song found that NW appeared to outperform TLS in lowering total amount of energy purchased (RMD and 95% CI: −0.064 [−0.125, −0.004]).
    sp[, policy.energy.effect.NW.vs.TLS := rnorm(1,mean = -0.064, sd = (-0.004-(-0.125))/3.92)]
    sp[policy.energy.effect.NW.vs.TLS > 0, policy.energy.effect.NW.vs.TLS := 0]
    #if we are putting NS instead of TLS, Song found no significant diff between NS and TLS in lowering total amount of energy purchased (0.01 [−0.04, 0.05]).
    #SHOULD WE DO THAT THEN?? CHRIS
    sp[, policy.energy.effect.NS.vs.TLS := rnorm(1,mean = 0.01, sd = (0.05-(-0.04))/3.92)]
    sp[policy.energy.effect.NS.vs.TLS > 0, policy.energy.effect.NS.vs.TLS := 0]


    #using Chile effect (https://www.medrxiv.org/content/10.1101/2023.11.21.23298789v1.full.pdf)
    #they observed change in energy content before and after policy (calories/capita/day)
    #PHASE 1
    #Labeled -17.1% (-18.9,-15.3)
    #Non-labeled 0.0% (-2.4, 2.3)
    #Total -8.8% (-10.5, -7.1)
    #Phase 2
    #Labeled -23.0% (-26.2,-19.8)
    #Non-labeled 7.1% (2.3, 11.9)
    #Total -8.3% (-11.6, -5.0)
    policy.energy.effect.chile.NW.p1 <- rnorm(1,mean = -0.088, sd = (-0.071-(-0.105))/3.92)
    if (policy.energy.effect.chile.NW.p1 > 0) policy.energy.effect.chile.NW.p1 <- 0
    policy.energy.effect.chile.NW.p2 <- rnorm(1,mean = -0.083, sd = (-0.050-(-0.116))/3.92)
    if (policy.energy.effect.chile.NW.p2 > 0) policy.energy.effect.chile.NW.p2 <- 0
    policy.energy.effect.chile.NW.p2.withoutcomp <- rnorm(1,mean = -0.23, sd = (-0.198-(-0.262))/3.92)
    if (policy.energy.effect.chile.NW.p2.withoutcomp > 0) policy.energy.effect.chile.NW.p2.withoutcomp <- 0
    sp[, policy.energy.effect.compensation.chile.NW := rnorm(1,mean = 0.071, sd = (0.119-(-0.023))/3.92)]
    sp[policy.energy.effect.compensation.chile.NW < 0, policy.energy.effect.compensation.chile.NW := 0]

     #####EFFECT ON REFORMULATION
    #In a paper for Chile, they found a 3.9% reduction in energy content of packaged foods post-implementation
    #of the nutrient warning label policy (Scarpelli et al., 2020, https://doi.org/10.3390/nu12082371) [also a reduction for other nutrients of concern]
    #Paper finds evidence of reformulation in response to food labels (not label specific) (Shangguan et al., 2024). This paper also supports there being some reformulation, particularly when labels are mandatory -
    #(Ganderats-Fuentes & Morgan, 2023, https://doi.org/10.3390/nu15112630). Therefore, based on lack of traffic-light specific evidence, we will
    #assume same reformulation % as NW
    sp[, policy.energy.effect.reformulation.chile.NW := rnorm(1, mean = -0.039, sd = (0.125+(0.0495))/3.92)] # CKFix: negative sd
    sp[policy.energy.effect.reformulation.chile.NW > 0, policy.energy.effect.reformulation.chile.NW := 0]

    sp[, policy.energy.effect.reformulation.NW := rnorm(1, mean = -0.039, sd = (0.125+(0.0495))/3.92)] #added by RE
    sp[policy.energy.effect.reformulation.NW > 0, policy.energy.effect.reformulation.NW := 0]
    sp[, policy.energy.effect.reformulation.TLS := rnorm(1, mean = -0.039, sd = (0.125+(0.0495))/3.92)]
    sp[policy.energy.effect.reformulation.TLS > 0, policy.energy.effect.reformulation.TLS := 0]
    sp[, policy.energy.effect.reformulation.NS := rnorm(1, mean = -0.039, sd = (0.125+(0.0495))/3.92)] #added by RE
    sp[policy.energy.effect.reformulation.NS > 0, policy.energy.effect.reformulation.NS := 0]
    sp[, policy.energy.effect.reformulation.TLSsensi := rnorm(1, mean = -0.009, sd = (0.031+(0.049))/3.92)] #added by RE
    sp[policy.energy.effect.reformulation.TLSsensi > 0, policy.energy.effect.reformulation.TLSsensi := 0]
    ###RE SENT AN EMAIL TO THE AUTHOR TO HAVE THE CONFIDENCE INTERVAL
    #RE - CIs imputed


    #############COVERAGE OF THE LABELS
    #Proportion of products that currently feature a traffic light label (voluntary) in the UK = 75%
    #https://publications.parliament.uk/pa/cm201516/cmselect/cmhealth/465/465.pdf (p.34)
    #so, 25% of the products will be labeled if the TLS become mandatory
    label_extra_coverage.TLS <- 1-0.75

    #Proportion of products that currently feature a warning label in the UK = 0%
    #Potential coverage of warning labels in the UK? no info for Uk
    #For Chile, the overall proportion of products with any “high in” warnings significantly decreased from 51% (95% confidence interval [CI] 49–52) to 44% (95% CI 42–45) after the initial implementation of the law
    #(TABLE 1 https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1003220).
    #The decrease might be higher for PHASE 2 & 3 but we don’t have the information.
    #Thus, we assume that if implemented in the UK, 51% (49–52) of the packaged products will be labelled.
    label_coverage.NW <- rnorm(1,mean = 0.51, sd = (0.52-(-0.49))/3.92)
    if (label_coverage.NW < 0) label_coverage.NW <- 0
    #for sensitivity scenarios, #we assume that 51% (95% CI 49–52) of the packaged products will be impacted by this new mandatory policy
    #(i.e. they will be above the threshold for warning) for the first 4 years and then, after reformulation of the products, only 44% (42–45) will be impacted.
    label_coverage.NW_4y <- rnorm(1,mean = 0.44, sd = (0.45-(-0.42))/3.92)
    if (label_coverage.NW_4y < 0) label_coverage.NW_4y <- 0
    #assuming that NS is not implemented in the UK, 100% of the packaged products will be targetted
    label_coverage.NS <- 1
    #assuming that NS is not implemented in the UK but that the effect will be only on the products not already targetted by TLS
    label_coverage.NS_sensi <- label_extra_coverage.TLS

    #sc0. BEFORE Traffic light labeling is implemented as a voluntary policy
    #test to estimate the energy from packaged food before the effect of TLS on 75% of the products
    sp[, sc0_label_nrj := energy_purchased_packaged/(1 + (policy.energy.effect.TLS*(0.75)))]

    #sc1. Traffic light labeling is implemented as a mandatory policy.
    sp[, sc1_label_nrj := (energy_purchased_packaged*(1-label_extra_coverage.TLS)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.TLS) * (1-0))) * label_extra_coverage.TLS)]
    sp[year<impl_year, sc1_label_nrj := energy_purchased_packaged]

    #sc1sensiR. Traffic light labelling is implemented as a mandatory policy. #added by RE
    #we assume reformulation of -3.9% energy [industry effect only]
    sp[, sc1_sensiR_label_nrj := (energy_purchased_packaged*(1-label_extra_coverage.TLS)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.reformulation.TLS) * (1-0))) * label_extra_coverage.TLS)]
    sp[year<impl_year, sc1_sensiR_label_nrj := energy_purchased_packaged]
    
    #sc1sensiR2. Traffic light labelling is implemented as a mandatory policy. #added by RE
    #we assume reformulation of -0.9% energy [industry effect only]
    sp[, sc1_sensiR2_label_nrj := (energy_purchased_packaged*(1-label_extra_coverage.TLS)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.reformulation.TLSsensi) * (1-0))) * label_extra_coverage.TLS)]
    sp[year<impl_year, sc1_sensiR2_label_nrj := energy_purchased_packaged]

    #sc1sensiboth. Traffic light labelling is implemented as a mandatory policy #added by RE
    #consumer response and reformulation combined
    sp[, sc1_sensiboth_label_nrj := (energy_purchased_packaged*(1-label_extra_coverage.TLS)) +
         ((energy_purchased_packaged +
             ((energy_purchased_packaged*policy.energy.effect.TLS) * (1-0)) +
             (energy_purchased_packaged*policy.energy.effect.reformulation.TLS)) * label_extra_coverage.TLS)]

    #sc1sensi. Traffic light labeling is implemented as a mandatory policy BUT ONLY PART OF THE POPULATION IS INFLUENCED BY THE LABEL
    #FOR TLS: 45% of the consumers that are using the traffic light label when shopping (https://www.food.gov.uk/research/behaviour-and-perception/eating-well-choosing-better-tracker-survey-wave-8-2022?print=1).
    label_extra_coverage.consumer.TLS. <- (1-0.75)*0.45
    sp[, sc1_sensi_label_nrj := (energy_purchased_packaged*(1-label_extra_coverage.consumer.TLS.)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.TLS) * (1-0))) * label_extra_coverage.consumer.TLS.)]
    sp[year<impl_year, sc1_sensi_label_nrj := energy_purchased_packaged]

    #sc2. nutrient warning labels are implemented as a mandatory policy instead of voluntary traffic light labeling
    #in this one, we assume that a new NW label will have an additional effect on 100% of the targeted products
    sp[, sc2_label_nrj := (energy_purchased_packaged*(1-label_coverage.NW)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.NW.vs.TLS) * (1-0))) * label_coverage.NW)]
    sp[year<impl_year, sc2_label_nrj := energy_purchased_packaged]

    #sc2sensiR. nutrient warning labels are implemented as a mandatory policy instead of voluntary traffic light labeling #added by RE
    #we assume reformulation of -3.9% energy [industry effect only]
    sp[, sc2_sensiR_label_nrj := (energy_purchased_packaged*(1-label_coverage.NW)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.reformulation.NW) * (1-0))) * label_coverage.NW)]
    sp[year<impl_year, sc2_sensiR_label_nrj := energy_purchased_packaged]

    #sc2sensiboth. nutrient warning labels are implemented as a mandatory policy instead of voluntary traffic light labeling #added by RE
    #consumer response and reformulation combined
    sp[, sc2_sensiboth_label_nrj := (energy_purchased_packaged*(1-label_coverage.NW)) +
         ((energy_purchased_packaged +
             ((energy_purchased_packaged*policy.energy.effect.NW.vs.TLS) * (1-0)) +
             (energy_purchased_packaged*policy.energy.effect.reformulation.NW)) * label_coverage.NW)]

    #sc2sensiT. nutrient warning labels are implemented as a mandatory policy instead of voluntary traffic light labeling
    #in this one, we assume that a new NW label will have an additional effect on 100% of the targeted products
    #we assume that 51% (95% CI 49–52) of the packaged products will be impacted by this new mandatory policy (i.e. they will be above the threshold for warning)
    #for the first 4 years and then, after reformulation of the products, only 44% (42–45) will be impacted. #now changed to 1 year
    sp[, sc2_sensiT_label_nrj := (energy_purchased_packaged*(1-label_coverage.NW_4y)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.NW.vs.TLS) * (1-0))) * label_coverage.NW_4y)]
    sp[year<impl_year, sc2_sensiT_label_nrj := energy_purchased_packaged]
    sp[year>=impl_year & year<=(impl_year+1), sc2_sensiT_label_nrj := (energy_purchased_packaged*(1-label_coverage.NW)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.NW.vs.TLS) * (1-0))) * label_coverage.NW)]

    #sc2sensiC. nutrient warning labels are implemented as a mandatory policy instead of voluntary traffic light labeling
    #in this one, we assume that a new NW label will have an additional effect on 100% of the targeted products
    #We add a compensation: increase of energy from non-labelled products by 7.1% (2.3, 11.9) as found in Tallie 2023 for Chile
    sp[, sc2_sensiC_label_nrj := ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.compensation.chile.NW))) * (1-label_coverage.NW)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.NW.vs.TLS) * (1-0))) * label_coverage.NW)]
    sp[year<impl_year, sc2_sensiC_label_nrj := energy_purchased_packaged]

    #sc2sensiTC. nutrient warning labels are implemented as a mandatory policy instead of voluntary traffic light labeling
    #in this one, we assume that a new NW label will have an additional effect on 100% of the targeted products
    #we are seeing a decrease in the number of prodcuts having the label (from 51 to 44)
    #We add a compensation: increase of energy from non-labelled products by 7.1% (2.3, 11.9) as found in Tallie 2023 for Chile
    sp[, sc2_sensiTC_label_nrj := ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.compensation.chile.NW))) * (1-label_coverage.NW_4y)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.NW.vs.TLS) * (1-0))) * label_coverage.NW_4y)]
    sp[year<impl_year, sc2_sensiTC_label_nrj := energy_purchased_packaged]
    sp[year>=impl_year & year<=(impl_year+4), sc2_sensiTC_label_nrj := ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.compensation.chile.NW))) * (1-label_coverage.NW)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.NW.vs.TLS) * (1-0))) * label_coverage.NW)]

    #sc3. NutriScore is implemented as a mandatory policy instead of voluntary traffic light labeling
    #as NS vs TLS is not significant, should I still model that??
    sp[, sc3_label_nrj := (energy_purchased_packaged*(1-label_coverage.NS)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.NS.vs.TLS) * (1-0))) * label_coverage.NS)]
    sp[year<impl_year, sc3_label_nrj := energy_purchased_packaged]

    #sc3sensi. NutriScore is implemented as a mandatory policy instead of voluntary traffic light labeling
    #but NS will have an effect only on the products non labelled by traffic light already (25% of the products) thus we are using the effect NS vs control
    sp[, sc3_sensi_label_nrj := (energy_purchased_packaged*(1-label_coverage.NS_sensi)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.NS) * (1-0))) * label_coverage.NS_sensi)]
    sp[year<impl_year, sc3_sensi_label_nrj := energy_purchased_packaged]

    #sc3sensiR. Nutri Score is implemented as a mandatory policy instead of voluntary traffic light labeling #added by RE
    #we assume reformulation of -3.9% energy [industry effect only] on 25% of products not covered by TLS
    sp[, sc3_sensiR_label_nrj := (energy_purchased_packaged*(1-label_coverage.NS_sensi)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.reformulation.NS) * (1-0))) * label_coverage.NS_sensi)]
    sp[year<impl_year, sc3_sensiR_label_nrj := energy_purchased_packaged]

    #sc3sensiboth. Nutri Score is implemented as a mandatory policy instead of voluntary traffic light labeling #added by RE
    #consumer response and reformulation combined
    sp[, sc3_sensiboth_label_nrj := (energy_purchased_packaged*(1-label_coverage.NS_sensi)) +
         ((energy_purchased_packaged +
             ((energy_purchased_packaged*policy.energy.effect.NS) * (1-0)) +
             (energy_purchased_packaged*policy.energy.effect.reformulation.NS)) * label_coverage.NS_sensi)]


    #sc4. Black octagon is implemented as a mandatory policy instead of voluntary traffic light labeling
    #in this one, we assume that a new NW label will have an additional effect on 100% of the targeted products
    #using Chile results IN TWO PHASE INCLUDING COMPENSATION
    #phase 1 - 2022-2024
    sp[, sc4_label_nrj := (energy_purchased_packaged*(1-label_coverage.NW)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.chile.NW.p1) * (1-0))) * label_coverage.NW)]
    #phase 2 - 2024 and after
    sp[year>=24, sc4_label_nrj := (energy_purchased_packaged*(1-label_coverage.NW)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.chile.NW.p2) * (1-0))) * label_coverage.NW)]
    sp[year<impl_year, sc4_label_nrj := energy_purchased_packaged]

    #sc5. Black octagon is implemented as a mandatory policy instead of voluntary traffic light labeling
    #using Chile results IN ONE PHASE ONLY INCLUDING COMPENSATION
    sp[, sc5_label_nrj := (energy_purchased_packaged*(1-label_coverage.NW)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.chile.NW.p2) * (1-0))) * label_coverage.NW)]
    sp[year<impl_year, sc5_label_nrj := energy_purchased_packaged]

    #sc5sensi. Black octagon is implemented as a mandatory policy instead of voluntary traffic light labeling
    #We do without the compensation
    sp[, sc5_sensi_label_nrj := (energy_purchased_packaged*(1-label_coverage.NW)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.chile.NW.p2.withoutcomp) * (1-0))) * label_coverage.NW)]
    sp[year<impl_year, sc5_sensi_label_nrj := energy_purchased_packaged]

    #sc6. Nutrient warning label is implemented as a mandatory policy instead of voluntary traffic light labeling
    #and this will led to a 3.9% reformulation of the energy of the products (EFFECT ON INDUSTRY NOT INDIVIDUAL)
    sp[, sc6_label_nrj := (energy_purchased_packaged*(1-label_coverage.NW)) + ((energy_purchased_packaged + ((energy_purchased_packaged*policy.energy.effect.reformulation.chile.NW) * (1-0))) * label_coverage.NW)]
    sp[year<impl_year, sc6_label_nrj := energy_purchased_packaged]


    #List all your scenarios here!
    lili_1 <- c("sc0_label","sc1_label","sc1_sensiR_label","sc1_sensiR2_label", "sc1_sensiboth_label", "sc1_sensi_label", "sc2_label", "sc2_sensiR_label", "sc2_sensiboth_label", "sc2_sensiT_label","sc2_sensiC_label","sc2_sensiTC_label",
                "sc3_label","sc3_sensi_label","sc3_sensiR_label", "sc3_sensiboth_label", "sc4_label","sc5_label","sc5_sensi_label","sc6_label")


    #change in nrj
    #sc1_label_nrjdiftest = ifelse(sc1_label_nrj > energy_purchased_packaged, 0, sc1_label_nrj - energy_purchased_packaged)

    fnrj <- function(x,y) ifelse(x > y, 0, x - y)

    lili_nrj <- paste(lili_1,"_nrj",sep="")
    lili_nrjdif <- paste(lili_1,"_nrjdif",sep="")
    sp[, (lili_nrjdif)  := lapply(.SD,fnrj,energy_purchased_packaged), .SDcols = lili_nrj]
    rm(fnrj,lili_nrj)

    #####TO DELETE
    ##we are redoing sc0 in purpose as it's a test for the death attributable to the current policy
    sp[, sc0_label_nrjdif := sc0_label_nrj - energy_purchased_packaged]



    #test <- sp[, sc1_label_nrj:sc6sensi_label_nrjdif]
    #test <- sp[,.(energy_purchased_packaged,sc1_label_nrj,sc1_label_nrjdif)]
    #test <- sp[,.(energy_purchased_packaged,sc6sensi_label_nrj,sc6sensi_label_nrjdif)]

    #linking the change in nrj and the change in weight (kg)
    #for men   = 17.7 * ((sc1_label_nrjdif*(4.2/1000))/PAL) #1 Kilocalorie = (4.2/1000) MJ
    #for women = 20.7 * ((sc1_label_nrjdif*(4.2/1000))/PAL)
    #Assuming a PAL of 1.5

    fbwss <- function(x,y) ifelse(y=="1", 17.7 * ((x*(4.2/1000))/1.5), 20.7 * ((x*(4.2/1000))/1.5)) #1 Kilocalorie = (4.2/1000) MJ
    lili_BWss <- paste(lili_1,"_BWss",sep="")

    sp[, (lili_BWss)  := lapply(.SD,fbwss,sex), .SDcols = lili_nrjdif]

    #we have an absolute change in weight and we cannot convert in a relative one, so we Assume that everyone has the same height.
    #We then calculate their weight given BMI and height (we assumed a height of 1.6m)
    sp[, wt_assumption := (bmi)*(1.6*1.6)]

    #we then apply the absolute difference to their calculated weight, and from their, calculate the relative difference
    fbwtrd <- function(x,y) (x+y)/(y)
    lili_wtrd <- paste(lili_1,"_wtrd",sep="")

    sp[, (lili_wtrd)  := lapply(.SD,fbwtrd,wt_assumption), .SDcols = lili_BWss]



    #Then we apply that relative change to the baseline BMI
    #All the calculations are independent of height so id doesn't matter what value will you use.
    #sc1_label_bmi = bmi*sc1_label_wtrd
    fbmi <- function(x,y) y*x
    lili_bmi <- paste(lili_1,"_bmi",sep="")

    sp[, (lili_bmi)  := lapply(.SD,fbmi,bmi), .SDcols = lili_wtrd]  #we are working on bmi and not bmi_lagged!! this is because we assumed that change in nrj change your bmi very quickly while your change in bmi take around 6y to change your RR
    #tt <- sp[,.(bmi,sc1_label_wtrd,sc1_label_bmi)]
    #summary(sp$bmi)
    #summary(sp$sc1_label_bmi)
    ###TEST FOR THE MODEl
    #10% DECREASE IN BMI
    #sp[, sc1_label_bmi_test := bmi*0.90)]

    sp[, c(lili_BWss) := NULL]


    #obesity prevalence
    lili_bmi_grp <- paste(lili_1,"_bmi_grp",sep="")

    fobprev <- function(x) fifelse(x<20,"<20",
                                  fifelse(x>=20 & x<25,"20-25",
                                         fifelse(x>=25 & x<30,"25-30",
                                                fifelse(x>=30,">=30", NA))))
    sp[, (lili_bmi_grp)  := lapply(.SD,fobprev), .SDcols = lili_bmi]
    #check
    # tt <- sp[,.(sc1_label_bmi,sc1_label_bmi_grp,sc2_label_bmi,sc2_label_bmi_grp)]
    # tt[,table(sc1_label_bmi_grp)]

    rm(lili_bmi_grp,lili_bmi,lili_nrjdif,lili_BWss,lili_wtrd,fbmi,fbwtrd,fobprev,fbwss)

    #New RR according to new BMI, assuming a lag time of 6 years (ERFC 2011: 5.7 years [SD 3.0–9.0]))
    #to be simpler, we will lag the bmi
    setkey(sp, pid, year) # NECESSARY!!!

    sp[, `:=` (sc0_label_bmi_lagged         = shift(sc0_label_bmi, 6),    #BMIss bring back to 4.56 to be comparative with rr.Obesity
               sc1_label_bmi_lagged         = shift(sc1_label_bmi, 6),
               sc1_sensiR_label_bmi_lagged  = shift(sc1_sensiR_label_bmi, 6),
               sc1_sensiR2_label_bmi_lagged  = shift(sc1_sensiR2_label_bmi, 6),
               sc1_sensiboth_label_bmi_lagged = shift(sc1_sensiboth_label_bmi, 6),
               sc1_sensi_label_bmi_lagged   = shift(sc1_sensi_label_bmi, 6),
               sc2_label_bmi_lagged         = shift(sc2_label_bmi, 6),
               sc2_sensiR_label_bmi_lagged  = shift(sc2_sensiR_label_bmi, 6),
               sc2_sensiboth_label_bmi_lagged = shift (sc2_sensiboth_label_bmi, 6),
               sc2_sensiT_label_bmi_lagged  = shift(sc2_sensiT_label_bmi, 6),
               sc2_sensiC_label_bmi_lagged  = shift(sc2_sensiC_label_bmi, 6),
               sc2_sensiTC_label_bmi_lagged = shift(sc2_sensiTC_label_bmi, 6),
               sc3_label_bmi_lagged         = shift(sc3_label_bmi, 6),
               sc3_sensi_label_bmi_lagged   = shift(sc3_sensi_label_bmi, 6),
               sc3_sensiR_label_bmi_lagged  = shift(sc3_sensiR_label_bmi, 6),
               sc3_sensiboth_label_bmi_lagged = shift (sc3_sensiboth_label_bmi, 6),
               sc4_label_bmi_lagged         = shift(sc4_label_bmi, 6),
               sc5_label_bmi_lagged         = shift(sc5_label_bmi, 6),
               sc5_sensi_label_bmi_lagged   = shift(sc5_sensi_label_bmi, 6),
               sc6_label_bmi_lagged         = shift(sc6_label_bmi, 6)), by = pid]

    #check <- sp[,.(pid,year,bmi,bmi_lagged,sc1mincomp_label_bmi,sc1mincomp_label_bmi_lagged)]

    lili_bmilag <- paste(lili_1,"_bmi_lagged",sep="")
    lili_lagbmichd <- paste(lili_1,"_lag.bmi.chd",sep="")
    lili_lagbmiisch <- paste(lili_1,"_lag.bmi.isch",sep="")
    lili_lagbmihaem <- paste(lili_1,"_lag.bmi.haem",sep="")

    #then we calculate the risk with the lagged bmi to take into account the 6y lag time
    frr.indiv.bmi <- function(x,y) y^((x-20)/4.56) #BMIss bring back to 4.56 to be comparative with rr.Obesity
    sp[, (lili_lagbmichd)  := lapply(.SD,frr.indiv.bmi,rr.obesity.chd), .SDcols = lili_bmilag]
    sp[, (lili_lagbmiisch) := lapply(.SD,frr.indiv.bmi,rr.obesity.isch), .SDcols = lili_bmilag]
    sp[, (lili_lagbmihaem) := lapply(.SD,frr.indiv.bmi,rr.obesity.haem), .SDcols = lili_bmilag]
    #for info, the formula for sc1_label for example then is
    #sc1_label_lag.bmi.chd = rr.obesity.chd^((sc1_label_bmi_lagged-20)/4.56)
    rm(frr.indiv.bmi)


    #the people bmi<20 have no risk, not a protective risk
    sp[sc0_label_bmi_lagged<20,      `:=` (sc0_label_lag.bmi.chd=1,sc0_label_lag.bmi.isch=1,sc0_label_lag.bmi.haem=1)]
    sp[sc1_label_bmi_lagged<20,      `:=` (sc1_label_lag.bmi.chd=1,sc1_label_lag.bmi.isch=1,sc1_label_lag.bmi.haem=1)]
    sp[sc1_sensiR_label_bmi_lagged<20,      `:=` (sc1_sensiR_label_lag.bmi.chd=1,sc1_sensiR_label_lag.bmi.isch=1,sc1_sensiR_label_lag.bmi.haem=1)]
    sp[sc1_sensiR2_label_bmi_lagged<20,      `:=` (sc1_sensiR2_label_lag.bmi.chd=1,sc1_sensiR2_label_lag.bmi.isch=1,sc1_sensiR2_label_lag.bmi.haem=1)]
    sp[sc1_sensiboth_label_bmi_lagged<20,      `:=` (sc1_sensiboth_label_lag.bmi.chd=1,sc1_sensiboth_label_lag.bmi.isch=1,sc1_sensiboth_label_lag.bmi.haem=1)]
    sp[sc1_sensi_label_bmi_lagged<20,`:=` (sc1_sensi_label_lag.bmi.chd=1,sc1_sensi_label_lag.bmi.isch=1,sc1_sensi_label_lag.bmi.haem=1)]
    sp[sc2_label_bmi_lagged<20,      `:=` (sc2_label_lag.bmi.chd=1,sc2_label_lag.bmi.isch=1,sc2_label_lag.bmi.haem=1)]
    sp[sc2_sensiR_label_bmi_lagged<20,      `:=` (sc2_sensiR_label_lag.bmi.chd=1,sc2_sensiR_label_lag.bmi.isch=1,sc2_sensiR_label_lag.bmi.haem=1)]
    sp[sc2_sensiboth_label_bmi_lagged<20,      `:=` (sc2_sensiboth_label_lag.bmi.chd=1,sc2_sensiboth_label_lag.bmi.isch=1,sc2_sensiboth_label_lag.bmi.haem=1)]
    sp[sc2_sensiT_label_bmi_lagged<20,`:=` (sc2_sensiT_label_lag.bmi.chd=1,sc2_sensiT_label_lag.bmi.isch=1,sc2_sensiT_label_lag.bmi.haem=1)]
    sp[sc2_sensiC_label_bmi_lagged<20,`:=` (sc2_sensiC_label_lag.bmi.chd=1,sc2_sensiC_label_lag.bmi.isch=1,sc2_sensiC_label_lag.bmi.haem=1)]
    sp[sc2_sensiTC_label_bmi_lagged<20,`:=` (sc2_sensiTC_label_lag.bmi.chd=1,sc2_sensiTC_label_lag.bmi.isch=1,sc2_sensiTC_label_lag.bmi.haem=1)]
    sp[sc3_label_bmi_lagged<20,      `:=` (sc3_label_lag.bmi.chd=1,sc3_label_lag.bmi.isch=1,sc3_label_lag.bmi.haem=1)]
    sp[sc3_sensi_label_bmi_lagged<20,`:=` (sc3_sensi_label_lag.bmi.chd=1,sc3_sensi_label_lag.bmi.isch=1,sc3_sensi_label_lag.bmi.haem=1)]
    sp[sc3_sensiR_label_bmi_lagged<20,      `:=` (sc3_sensiR_label_lag.bmi.chd=1,sc3_sensiR_label_lag.bmi.isch=1,sc3_sensiR_label_lag.bmi.haem=1)]
    sp[sc3_sensiboth_label_bmi_lagged<20,      `:=` (sc3_sensiboth_label_lag.bmi.chd=1,sc3_sensiboth_label_lag.bmi.isch=1,sc3_sensiboth_label_lag.bmi.haem=1)]
    sp[sc4_label_bmi_lagged<20,      `:=` (sc4_label_lag.bmi.chd=1,sc4_label_lag.bmi.isch=1,sc4_label_lag.bmi.haem=1)]
    sp[sc5_label_bmi_lagged<20,      `:=` (sc5_label_lag.bmi.chd=1,sc5_label_lag.bmi.isch=1,sc5_label_lag.bmi.haem=1)]
    sp[sc5_sensi_label_bmi_lagged<20,`:=` (sc5_sensi_label_lag.bmi.chd=1,sc5_sensi_label_lag.bmi.isch=1,sc5_sensi_label_lag.bmi.haem=1)]
    sp[sc6_label_bmi_lagged<20,      `:=` (sc6_label_lag.bmi.chd=1,sc6_label_lag.bmi.isch=1,sc6_label_lag.bmi.haem=1)]

    #cleaning var
    sp[, c(lili_bmilag) := NULL]

    ##TEST FOR A LOOP, NOT WORKING YET
    # test <- copy(sp)
    # test <- test[,.(sc1_label_bmi_lagged,sc1_label_lag.bmi.chd,sc1mincomp_label_bmi_lagged,sc1mincomp_label_lag.bmi.chd)]

    # lili_bmilagTEST <- c("sc1_label_bmi_lagged","sc1mincomp_label_bmi_lagged")
    # lili_lagchdTEST <- c("sc1_label_lag.chdTEST22","sc1mincomp_label_lag.chdTEST22")
    # lili_lagchdTEST <- c("sc1_label_lag.chd","sc1mincomp_label_lag.chd")
    #working but I can not put a list in y (instead of sc1_label_bmi_lagged in the current example...)
    #test[, (lili_lagchdTEST)  := lapply(.SD, function(x,y) ifelse(y<20 & !is.na(y),1,x), sc1_label_bmi_lagged), .SDcols = lili_lagchdTEST]
    #tt <- test[sc1_label_lag.chdTEST!=sc1_label_lag.chd,] #should be nrow=12995

    #not working at all
    #test[, mapply(function(X,Y) {(lili_lagchdTEST) := lapply(X, function(x,y) ifelse(y<20 & !is.na(y),1,x), Y)}, X=lili_lagchdTEST, Y=lili_bmilagTEST)]



    # Estimate new mortality in scenario
    lili_mr_bmi_chd <- paste(lili_1,"_mr_bmi_chd",sep="")
    lili_mr_bmi_isch <- paste(lili_1,"_mr_bmi_isch",sep="")
    lili_mr_bmi_haem <- paste(lili_1,"_mr_bmi_haem",sep="")
    mortp0 <- function(x,y) y*x
    sp[, (lili_mr_bmi_chd)  := lapply(.SD,mortp0,p0_bmi_mrtl_chd), .SDcols = lili_lagbmichd]
    sp[, (lili_mr_bmi_isch)  := lapply(.SD,mortp0,p0_bmi_mrtl_isch), .SDcols = lili_lagbmiisch]
    sp[, (lili_mr_bmi_haem)  := lapply(.SD,mortp0,p0_bmi_mrtl_haem), .SDcols = lili_lagbmihaem]
    #for info, the formula for sc1_label for example then is
    #sc1_label_mr_bmi_chd = p0_bmi_mrtl_chd * sc1_label_lag.bmi.chd
    rm(mortp0)

    #as <20  have rr of 1 already no need to do >=20 here, just p0*rr

    #tt<-sp[,.(pid,year,rr.obesity.chd,bmi,bmi_lagged,rr.obesity.chd.indiv,mr_chd,sc3_label_bmi,sc3_label_bmi_lagged,sc3_label_lag.bmi.chd,p0_bmi_mrtl_chd,sc3_label_mr_bmi_chd)]

    #cleaning var
    sp[, c(lili_lagbmichd,lili_lagbmiisch,lili_lagbmihaem) := NULL]

    lili <- c(lili_mr_bmi_chd,lili_mr_bmi_isch,lili_mr_bmi_haem)
    for (h in lili_1) {
      pid_toremove <- integer(0) # a place holder

      hchd <- paste0(h,"_mr_bmi_chd")
      hh <- sub("_mr_", "_dead_", hchd)
      hhh <- sub("_mr_", "_longdead_", hchd)
      # hh <- paste(sub("_mr.*", "", h),"_dead_bmi_chd",sep="")
      # hhh <- paste(sub("_mr.*", "", h),"_longdead_bmi_chd",sep="")
      for (i in sort(unique(sp$year))) {
        #print(paste0("Year: ", i))

        # rebalance the mortality rates to account for simulants removed from the
        # population when dead
        if (length(pid_toremove) > 0L) {
          correction_factor <- as.numeric(sp[year == i, lapply(.SD,sum), .SDcols=hchd]/sp[year == i & !pid %in% pid_toremove, lapply(.SD,sum), .SDcols=hchd])
          sp[year == i, (hchd) := lapply(.SD,"*",correction_factor), .SDcols=hchd]
          #print(paste0("correction factor: ", correction_factor))
        }

        sp[year == i & !pid %in% pid_toremove, paste(hh) := lapply(.SD,">",rn_mrtl), .SDcols=hchd]
        pid_toremove <- sp[year <= i & (sp[[hh]]), pid]
        #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
      }
      sp[, paste(hhh) := lapply(.SD,cumsum), by = pid, .SDcols=hh]
      sp[sp[[hh]] == FALSE & sp[[hhh]] == 1, paste(hh) := NA] # Best for later calculations.
      # So dead == FALSE mean alive, dead == TRUE means died that year,
      # and dead == NA means died in a previous year
      #dcast(sp, year~sp[[hh]]) # Note NA some years is decreasing because simulants get removed from the dataset because their age is > 89
      #table(sp$sc1_label_dead_bmi_chd)

      hisch <- paste0(h,"_mr_bmi_isch")
      hh <- sub("_mr_", "_dead_", hisch)
      hhh <- sub("_mr_", "_longdead_", hisch)
      # hh <- paste(sub("_mr.*", "", h),"_dead_bmi_isch",sep="")
      # hhh <- paste(sub("_mr.*", "", h),"_longdead_bmi_isch",sep="")
      for (i in sort(unique(sp$year))) {
        #print(paste0("Year: ", i))

        # rebalance the mortality rates to account for simulants removed from the
        # population when dead
        if (length(pid_toremove) > 0L) {
          correction_factor <- as.numeric(sp[year == i, lapply(.SD,sum), .SDcols=hisch]/sp[year == i & !pid %in% pid_toremove, lapply(.SD,sum), .SDcols=hisch])
          sp[year == i, (hisch) := lapply(.SD,"*",correction_factor), .SDcols=hisch]
          #print(paste0("correction factor: ", correction_factor))
        }

        sp[year == i & !pid %in% pid_toremove, paste(hh) := lapply(.SD,">",rn_mrtl), .SDcols=hisch]
        pid_toremove <- sp[year <= i & (sp[[hh]]), pid]
        #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
      }
      sp[, paste(hhh) := lapply(.SD,cumsum), by = pid, .SDcols=hh]
      sp[sp[[hh]] == FALSE & sp[[hhh]] == 1, paste(hh) := NA] # Best for later calculations.

      hhaem <- paste0(h,"_mr_bmi_haem")
      hh <- sub("_mr_", "_dead_", hhaem)
      hhh <- sub("_mr_", "_longdead_", hhaem)
      # hh <- paste(sub("_mr.*", "", h),"_dead_bmi_haem",sep="")
      # hhh <- paste(sub("_mr.*", "", h),"_longdead_bmi_haem",sep="")
      for (i in sort(unique(sp$year))) {
        #print(paste0("Year: ", i))

        # rebalance the mortality rates to account for simulants removed from the
        # population when dead
        if (length(pid_toremove) > 0L) {
          correction_factor <- as.numeric(sp[year == i, lapply(.SD,sum), .SDcols=hhaem]/sp[year == i & !pid %in% pid_toremove, lapply(.SD,sum), .SDcols=hhaem])
          sp[year == i, (hhaem) := lapply(.SD,"*",correction_factor), .SDcols=hhaem]
          #print(paste0("correction factor: ", correction_factor))
        }

        sp[year == i & !pid %in% pid_toremove, paste(hh) := lapply(.SD,">",rn_mrtl), .SDcols=hhaem]
        pid_toremove <- sp[year <= i & (sp[[hh]]), pid]
        #print(paste0("PID to remove: ", paste(pid_toremove, collapse = ",")))
      }
      sp[, paste(hhh) := lapply(.SD,cumsum), by = pid, .SDcols=hh]
      sp[sp[[hh]] == FALSE & sp[[hhh]] == 1, paste(hh) := NA] # Best for later calculations.
    }

    #cleaning var
    sp[, c(lili_mr_bmi_chd,lili_mr_bmi_isch,lili_mr_bmi_haem) := NULL]
    rm(lili_mr_bmi_chd,lili_mr_bmi_isch,lili_mr_bmi_haem, lili)


    ## Results (distributional impact)
    # CKNote: This can be done in a more efficient way, but I'm keeping it as is for now
    lili_totY <- c()
    testdata2 <- data.table()
    top_year <- (init_year+sim_hor)-1
    for (y in init_year:top_year){
      byyear <- sp[year==y,]
      testdata1 <- data.table()
      for (d in c("dead_bmi_chd","dead_bmi_isch","dead_bmi_haem")){
        for (ses in SES){
          byses <- byyear[SES==ses,]
          lili_totY <- c(lili_totY, paste("expected",d, "SES",ses, "Year", y, sep="_"))
          tt <- data.table((sum(byses[[d]]==TRUE,na.rm=T)) * scale_factor)
          names(tt) <- c(paste("expected",d, "SES",ses, "Year", y, sep="_"))
          testdata1 <- cbind(testdata1, tt)
          rm(tt)
          for (e in lili_1){
            sc <- paste(e, d,sep="_")
            lili_totY <- c(lili_totY, paste("DPPs", d, e, "SES",ses, "Year", y, sep="_"))
            tt <-  data.table((sum(byses[[d]]==TRUE,na.rm=T) - sum(byses[[sc]]==TRUE,na.rm=T)) * scale_factor)
            names(tt) <- c(paste("DPPs", d, e, "SES",ses, "Year", y, sep="_"))
            testdata1 <- cbind(testdata1, tt)
            rm(tt)
          }
          rm(byses)
        }
      }
      testdata2 <- cbind(testdata2, testdata1)
      rm(byyear,testdata1)
    }
    rm(y,d,e,ses)
    #saveRDS(lili_totY, file="lili_totY.RData")

    ##obesity prevalence
    lili_totOY <- c()
    obesity_prevalence <- data.table(V1=c("<20", ">=30", "20-25", "25-30"))
    top_year <- (init_year+sim_hor)-1
    for (y in init_year:top_year){
      byyear <- sp[year==y,]
      testdata1 <- data.table(V1=c("<20", ">=30", "20-25", "25-30"))
      for (d in c("bmi_grp")){
        for (ses in SES){
          byses <- byyear[SES==ses,]
          lili_totOY <- c(lili_totOY, paste("expected",d, "SES",ses, "Year", y, sep="_"))
          tt <- data.table(prop.table(table(byses[[d]])))
          names(tt) <- c("V1", paste("expected",d, "SES",ses, "Year", y, sep="_"))
          testdata1 <- testdata1[tt, on=c("V1")]
          rm(tt)
          for (e in lili_1){
            sc <- paste(e, d,sep="_")
            lili_totOY <- c(lili_totOY, paste("expected", d, e, "SES",ses, "Year", y, sep="_"))
            tt <- data.table(prop.table(table(byses[[sc]])))
            names(tt) <- c("V1", paste("expected", d, e, "SES",ses, "Year", y, sep="_"))
            testdata1 <- testdata1[tt, on=c("V1")]
            rm(tt)
          }
          rm(byses)
        }
      }
      obesity_prevalence <- obesity_prevalence[testdata1, on=c("V1")]
      rm(byyear,testdata1)
    }
    rm(y,d,e,ses,top_year)
    #saveRDS(lili_totOY, file="lili_totOY.RData")

    obesity_prevalence2 <- obesity_prevalence[V1==">=30",]
    obesity_prevalence2 <-obesity_prevalence2[,V1:=NULL]
    #saveRDS(obesity_prevalence2, file="obesity_prevalence2.RData")

    lili_totYtest <- c(lili_totY, lili_totOY)
    saveRDS(lili_totYtest, file="lili_totYtest.RData")

    time.new <- difftime(Sys.time(),time.old, units = "secs")
    time.end <- difftime(Sys.time(),time.begin, units= "mins")
    print(paste0("Simulation: ",sim,"   Time spent: ",round(time.new, digits=2), " seconds    Total time: ", round(time.end, digits=2)," minutes"))

    #putting the results in the mat3 for each iterations
    for (nbl in 1:ncol(mat3)){
      if (nbl <= length(lili_totY)){
        mat3[sim,nbl] <- testdata2[[nbl]]
      } else {
        t <- nbl-length(lili_totY)
        mat3[sim,nbl] <- obesity_prevalence2[[t]]
      }
    }
    rm(nbl)
  }
)
##############END OF THE monte carlo
parallel::stopCluster(cl)
#mat3[]
#View(mat3)

#####PREPARING THE RESULTS
lili_totYtest <- readRDS("lili_totYtest.RData")
resultsYY2 <- data.table()
for (nbl in 1:(length(lili_totYtest))){
  resultsYY2 <- cbind(resultsYY2,  mat3[,nbl])
}
rm(nbl)
colnames(resultsYY2) <- lili_totYtest
#View(resultsYY2)

resultsYY2$sim <- seq.int(nrow(resultsYY2))

melt_data <- melt(resultsYY2, id = c("sim"))
melt_data[,`:=` (variable=as.character(variable))]
melt_data[,`:=` (year=substring(variable,nchar(variable)-1,nchar(variable)),
                 SES=substring(variable,nchar(variable)-8,nchar(variable)-8),
                 var2=substring(variable,1,nchar(variable)-14))]
melt_data[,`:=` (year=as.numeric(as.character(year)))]
melt_data[,variable:=NULL]
dcast_data <- dcast(melt_data, sim + year + SES ~var2, mean)

resultsYY <- dcast_data

####SAVING THE RESULTS
saveRDS(resultsYY, file=paste("results.test", toString(Sys.Date()), "rds", sep = "."))

##############################################################
##############################################################
##RESULTS
##############################################################
##############################################################
rm(list=ls())
##doing mc monte carlo
library(data.table)
library(readxl)
library(fst)
library(gamlss)
library(demography)
library(xlsx)
library(foreach)
library(qs) # save R objects to disk

# Load results
resultsYY <- readRDS("results.test.2024-10-03.rds") #change based on most recent results file

#results for year 2024-2043 AS WE ARE STARTING THE POLICIES IN 2024!!!!!!!!!!!!!!!!!!!!!
resultsYYall <- resultsYY[year>=24,]
resultsYYall <- cbind(resultsYYall[,sim:DPPs_dead_bmi_isch_sc6_label],resultsYYall[,expected_dead_bmi_chd:expected_dead_bmi_isch])

#we don't want results by year or SES
setDT(resultsYYall)
resultsYYall <- resultsYYall[, c("SES","year"):=NULL]
resultsYYall <- aggregate(resultsYYall[,2:ncol(resultsYYall)],by=list(resultsYYall$sim), sum)
setnames(resultsYYall, c("Group.1"), c("sim"))
setDT(resultsYYall)

#DOING CVD, which is the sum of CHD+ISCH+HAEM
resultsYYall[, `:=` (expected_dead_cvd=expected_dead_bmi_chd+expected_dead_bmi_isch+expected_dead_bmi_haem,
                     DPPs_dead_cvd_sc0_label=(DPPs_dead_bmi_chd_sc0_label+DPPs_dead_bmi_isch_sc0_label+DPPs_dead_bmi_haem_sc0_label),
                     DPPs_dead_cvd_sc1_label=(DPPs_dead_bmi_chd_sc1_label+DPPs_dead_bmi_isch_sc1_label+DPPs_dead_bmi_haem_sc1_label),
                     DPPs_dead_cvd_sc1_sensi_label=(DPPs_dead_bmi_chd_sc1_sensi_label+DPPs_dead_bmi_isch_sc1_sensi_label+DPPs_dead_bmi_haem_sc1_sensi_label),
                     DPPs_dead_cvd_sc1_sensiR_label=(DPPs_dead_bmi_chd_sc1_sensiR_label+DPPs_dead_bmi_isch_sc1_sensiR_label+DPPs_dead_bmi_haem_sc1_sensiR_label),
                     DPPs_dead_cvd_sc1_sensiR2_label=(DPPs_dead_bmi_chd_sc1_sensiR2_label+DPPs_dead_bmi_isch_sc1_sensiR2_label+DPPs_dead_bmi_haem_sc1_sensiR2_label),
                     DPPs_dead_cvd_sc1_sensiboth_label=(DPPs_dead_bmi_chd_sc1_sensiboth_label+DPPs_dead_bmi_isch_sc1_sensiboth_label+DPPs_dead_bmi_haem_sc1_sensiboth_label),
                     DPPs_dead_cvd_sc2_label=(DPPs_dead_bmi_chd_sc2_label+DPPs_dead_bmi_isch_sc2_label+DPPs_dead_bmi_haem_sc2_label),
                     DPPs_dead_cvd_sc2_sensiT_label=(DPPs_dead_bmi_chd_sc2_sensiT_label+DPPs_dead_bmi_isch_sc2_sensiT_label+DPPs_dead_bmi_haem_sc2_sensiT_label),
                     DPPs_dead_cvd_sc2_sensiC_label=(DPPs_dead_bmi_chd_sc2_sensiC_label+DPPs_dead_bmi_isch_sc2_sensiC_label+DPPs_dead_bmi_haem_sc2_sensiC_label),
                     DPPs_dead_cvd_sc2_sensiTC_label=(DPPs_dead_bmi_chd_sc2_sensiTC_label+DPPs_dead_bmi_isch_sc2_sensiTC_label+DPPs_dead_bmi_haem_sc2_sensiTC_label),
                     DPPs_dead_cvd_sc2_sensiR_label=(DPPs_dead_bmi_chd_sc2_sensiR_label+DPPs_dead_bmi_isch_sc2_sensiR_label+DPPs_dead_bmi_haem_sc2_sensiR_label),
                     DPPs_dead_cvd_sc2_sensiboth_label=(DPPs_dead_bmi_chd_sc2_sensiboth_label+DPPs_dead_bmi_isch_sc2_sensiboth_label+DPPs_dead_bmi_haem_sc2_sensiboth_label),
                     DPPs_dead_cvd_sc3_label=(DPPs_dead_bmi_chd_sc3_label+DPPs_dead_bmi_isch_sc3_label+DPPs_dead_bmi_haem_sc3_label),
                     DPPs_dead_cvd_sc3_sensi_label=(DPPs_dead_bmi_chd_sc3_sensi_label+DPPs_dead_bmi_isch_sc3_sensi_label+DPPs_dead_bmi_haem_sc3_sensi_label),
                     DPPs_dead_cvd_sc3_sensiR_label=(DPPs_dead_bmi_chd_sc3_sensiR_label+DPPs_dead_bmi_isch_sc3_sensiR_label+DPPs_dead_bmi_haem_sc3_sensiR_label),
                     DPPs_dead_cvd_sc3_sensiboth_label=(DPPs_dead_bmi_chd_sc3_sensiboth_label+DPPs_dead_bmi_isch_sc3_sensiboth_label+DPPs_dead_bmi_haem_sc3_sensiboth_label),
                     DPPs_dead_cvd_sc4_label=(DPPs_dead_bmi_chd_sc4_label+DPPs_dead_bmi_isch_sc4_label+DPPs_dead_bmi_haem_sc4_label),
                     DPPs_dead_cvd_sc5_label=(DPPs_dead_bmi_chd_sc5_label+DPPs_dead_bmi_isch_sc5_label+DPPs_dead_bmi_haem_sc5_label),
                     DPPs_dead_cvd_sc5_sensi_label=(DPPs_dead_bmi_chd_sc5_sensi_label+DPPs_dead_bmi_isch_sc5_sensi_label+DPPs_dead_bmi_haem_sc5_sensi_label),
                     DPPs_dead_cvd_sc6_label=(DPPs_dead_bmi_chd_sc6_label+DPPs_dead_bmi_isch_sc6_label+DPPs_dead_bmi_haem_sc6_label))]


#number of expected death
resultsYYall[, as.list(quantile(expected_dead_cvd, probs = c(0.5, 0.025, 0.975)))]
# resultsYYall[, as.list(quantile(expected_dead_bmi_chd, probs = c(0.5, 0.025, 0.975)))]
# resultsYYall[, as.list(quantile(expected_dead_bmi_isch, probs = c(0.5, 0.025, 0.975)))]
# resultsYYall[, as.list(quantile(expected_dead_bmi_haem, probs = c(0.5, 0.025, 0.975)))]


############################
## LABELLING
############################
###############################################################
##      DPPs CVD by scenario
###############################################################
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc1_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc1_sensi_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc1_sensiR_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc1_sensiR2_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc1_sensiboth_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc2_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc2_sensiT_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc2_sensiC_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc2_sensiTC_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc2_sensiR_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc2_sensiboth_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc3_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc3_sensi_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc3_sensiR_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc3_sensiboth_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc4_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc5_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc5_sensi_label, probs = c(0.5, 0.025, 0.975)))]
resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc6_label, probs = c(0.5, 0.025, 0.975)))]
#resultsYYall[, as.list(quantile(DPPs_dead_cvd_sc0_label, probs = c(0.5, 0.025, 0.975)))]

###################################################
### DPPS RESULTS BY SES 
###################################################
#results for year 2024-2043 AS WE ARE STARTING THE POLICIES IN 2024!!!!!!!!!!!!!!!!!!!!!
resultsYYallSES <- resultsYY[year>=24,]
resultsYYallSES <- cbind(resultsYYallSES[,sim:DPPs_dead_bmi_isch_sc6_label],resultsYYallSES[,expected_dead_bmi_chd:expected_dead_bmi_isch])

setDT(resultsYYallSES)
resultsYYallSES <- resultsYYallSES[, c("year"):=NULL]
resultsYYallSES <- aggregate(resultsYYallSES[,3:ncol(resultsYYallSES)],by=list(resultsYYallSES$sim,resultsYYallSES$SES), sum)
setnames(resultsYYallSES, c("Group.1","Group.2"), c("sim","SES"))
setDT(resultsYYallSES)

#DOING CVD, which is the sum of CHD+ISCH+HAEM
resultsYYallSES[, `:=` (expected_dead_cvd=expected_dead_bmi_chd+expected_dead_bmi_isch+expected_dead_bmi_haem,
                        DPPs_dead_cvd_sc0_label=(DPPs_dead_bmi_chd_sc0_label+DPPs_dead_bmi_isch_sc0_label+DPPs_dead_bmi_haem_sc0_label),
                        DPPs_dead_cvd_sc1_label=(DPPs_dead_bmi_chd_sc1_label+DPPs_dead_bmi_isch_sc1_label+DPPs_dead_bmi_haem_sc1_label),
                        DPPs_dead_cvd_sc1_sensi_label=(DPPs_dead_bmi_chd_sc1_sensi_label+DPPs_dead_bmi_isch_sc1_sensi_label+DPPs_dead_bmi_haem_sc1_sensi_label),
                        DPPs_dead_cvd_sc1_sensiR_label=(DPPs_dead_bmi_chd_sc1_sensiR_label+DPPs_dead_bmi_isch_sc1_sensiR_label+DPPs_dead_bmi_haem_sc1_sensiR_label),
                        DPPs_dead_cvd_sc1_sensiR2_label=(DPPs_dead_bmi_chd_sc1_sensiR2_label+DPPs_dead_bmi_isch_sc1_sensiR2_label+DPPs_dead_bmi_haem_sc1_sensiR2_label),
                        DPPs_dead_cvd_sc1_sensiboth_label=(DPPs_dead_bmi_chd_sc1_sensiboth_label+DPPs_dead_bmi_isch_sc1_sensiboth_label+DPPs_dead_bmi_haem_sc1_sensiboth_label),
                        DPPs_dead_cvd_sc2_label=(DPPs_dead_bmi_chd_sc2_label+DPPs_dead_bmi_isch_sc2_label+DPPs_dead_bmi_haem_sc2_label),
                        DPPs_dead_cvd_sc2_sensiT_label=(DPPs_dead_bmi_chd_sc2_sensiT_label+DPPs_dead_bmi_isch_sc2_sensiT_label+DPPs_dead_bmi_haem_sc2_sensiT_label),
                        DPPs_dead_cvd_sc2_sensiC_label=(DPPs_dead_bmi_chd_sc2_sensiC_label+DPPs_dead_bmi_isch_sc2_sensiC_label+DPPs_dead_bmi_haem_sc2_sensiC_label),
                        DPPs_dead_cvd_sc2_sensiTC_label=(DPPs_dead_bmi_chd_sc2_sensiTC_label+DPPs_dead_bmi_isch_sc2_sensiTC_label+DPPs_dead_bmi_haem_sc2_sensiTC_label),
                        DPPs_dead_cvd_sc2_sensiR_label=(DPPs_dead_bmi_chd_sc2_sensiR_label+DPPs_dead_bmi_isch_sc2_sensiR_label+DPPs_dead_bmi_haem_sc2_sensiR_label),
                        DPPs_dead_cvd_sc2_sensiboth_label=(DPPs_dead_bmi_chd_sc2_sensiboth_label+DPPs_dead_bmi_isch_sc2_sensiboth_label+DPPs_dead_bmi_haem_sc2_sensiboth_label),
                        DPPs_dead_cvd_sc3_label=(DPPs_dead_bmi_chd_sc3_label+DPPs_dead_bmi_isch_sc3_label+DPPs_dead_bmi_haem_sc3_label),
                        DPPs_dead_cvd_sc3_sensi_label=(DPPs_dead_bmi_chd_sc3_sensi_label+DPPs_dead_bmi_isch_sc3_sensi_label+DPPs_dead_bmi_haem_sc3_sensi_label),
                        DPPs_dead_cvd_sc3_sensiR_label=(DPPs_dead_bmi_chd_sc3_sensiR_label+DPPs_dead_bmi_isch_sc3_sensiR_label+DPPs_dead_bmi_haem_sc3_sensiR_label),
                        DPPs_dead_cvd_sc3_sensiboth_label=(DPPs_dead_bmi_chd_sc3_sensiboth_label+DPPs_dead_bmi_isch_sc3_sensiboth_label+DPPs_dead_bmi_haem_sc3_sensiboth_label),
                        DPPs_dead_cvd_sc4_label=(DPPs_dead_bmi_chd_sc4_label+DPPs_dead_bmi_isch_sc4_label+DPPs_dead_bmi_haem_sc4_label),
                        DPPs_dead_cvd_sc5_label=(DPPs_dead_bmi_chd_sc5_label+DPPs_dead_bmi_isch_sc5_label+DPPs_dead_bmi_haem_sc5_label),
                        DPPs_dead_cvd_sc5_sensi_label=(DPPs_dead_bmi_chd_sc5_sensi_label+DPPs_dead_bmi_isch_sc5_sensi_label+DPPs_dead_bmi_haem_sc5_sensi_label),
                        DPPs_dead_cvd_sc6_label=(DPPs_dead_bmi_chd_sc6_label+DPPs_dead_bmi_isch_sc6_label+DPPs_dead_bmi_haem_sc6_label))]


#for simplicity, i only select SES 1 and 5, but feel free to delete this line if you want all SES level
resultsYYallSES <- resultsYYallSES[SES=="1" |SES=="5",]


#number of expected death for most deprived (SES==1) and least deprived (SES==5) 
resultsYYallSES[, as.list(quantile(expected_dead_cvd, probs = c(0.5, 0.025, 0.975))),by=c("SES")]

resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc1_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc1_sensi_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc1_sensiR_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc1_sensiR2_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc1_sensiboth_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc2_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc2_sensiT_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc2_sensiC_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc2_sensiTC_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc2_sensiR_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc2_sensiboth_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc3_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc3_sensi_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc3_sensiR_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc3_sensiboth_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc4_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc5_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc5_sensi_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc6_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]
#resultsYYallSES[, as.list(quantile(DPPs_dead_cvd_sc0_label, probs = c(0.5, 0.025, 0.975))),by=c("SES")]

####################################################################################
#####################RESULTS ON OBESITY
####################################################################################
#we don't want results by ses
resultsYYO <- cbind(resultsYY[,sim:SES],resultsYY[,expected_bmi_grp:expected_bmi_grp_sc6_label])
setDT(resultsYYO)
resultsYYO <- resultsYYO[, c("SES"):=NULL]
resultsYYO <- aggregate(resultsYYO[,3:ncol(resultsYYO)],by=list(resultsYYO$sim,resultsYYO$year), mean)
setnames(resultsYYO, c("Group.1","Group.2"), c("sim","year"))
setDT(resultsYYO)

#expected obesity prevalence
resultsYYO[, as.list(quantile(expected_bmi_grp*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
#difference in prevalence by scenario (-2 means a decrease in 2 percentage point in the obesity prevalence)
resultsYYO[, as.list(quantile((expected_bmi_grp_sc1_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc1_sensi_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc1_sensiR_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc1_sensiR2_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc1_sensiboth_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc2_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc2_sensiT_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc2_sensiC_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc2_sensiTC_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc2_sensiR_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc2_sensiboth_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc3_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc3_sensi_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc3_sensiR_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc3_sensiboth_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc4_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc5_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc5_sensi_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
resultsYYO[, as.list(quantile((expected_bmi_grp_sc6_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]
#resultsYYO[, as.list(quantile((expected_bmi_grp_sc0_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year")]

###BY SES
resultsYYOSES <- cbind(resultsYY[,sim:SES],resultsYY[,expected_bmi_grp:expected_bmi_grp_sc6_label])
setDT(resultsYYOSES)
resultsYYOSES <- aggregate(resultsYYOSES[,4:ncol(resultsYYOSES)],by=list(resultsYYOSES$sim,resultsYYOSES$year,resultsYYOSES$SES), mean)
setnames(resultsYYOSES, c("Group.1","Group.2","Group.3"), c("sim","year","SES"))
setDT(resultsYYOSES)

resultsYYOSES <- resultsYYOSES[SES==1 |SES==5,]

#expected obesity prevalence for most deprived (SES==1) and least deprived (SES==5)
resultsYYOSES[, as.list(quantile(expected_bmi_grp*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
#difference in prevalence by scenario (-2 means a decrease in 2 percentage point in the obesity prevalence)
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc1_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc1_sensi_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc1_sensiR_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc1_sensiR2_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc1_sensiboth_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc2_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc2_sensiT_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc2_sensiC_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc2_sensiTC_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc2_sensiR_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc2_sensiboth_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc3_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc3_sensi_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc3_sensiR_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc3_sensiboth_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc4_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc5_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc5_sensi_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc6_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]
#resultsYYOSES[, as.list(quantile((expected_bmi_grp_sc0_label-expected_bmi_grp)*100, probs = c(0.5, 0.025, 0.975))),by=c("year","SES")]



######end

