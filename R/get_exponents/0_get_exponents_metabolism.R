#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2019.05.02: Max Lindmark
#
# - This code estimates size-scaling exponents from data on the rate if it was not 
#   provided in the original paper
# 
# A. Load libraries
#
# B. Fit log-log model by species
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. LOAD LIBRARIES & READ DATA ====================================================
rm(list = ls())

# Provide package names
pkgs <- c("dplyr",
          "tidyr",
          "tidylog",
          "readxl",
          "RCurl",
          "magrittr")

# Install packages
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}

# Load all packages
lapply(pkgs, library, character.only = TRUE)

# Print package version
# script <- getURL("https://raw.githubusercontent.com/maxlindmark/scaling/master/R/functions/packageInfo.R", ssl.verifypeer = FALSE)
# eval(parse(text = script))
# pkg_info(pkgs)

# package   version
# 1   dplyr   0.8.0.1
# 2   RCurl 1.95-4.10
# 3  readxl     1.3.1
# 4 tidylog     0.1.0
# 5   tidyr     0.8.3

# Will crate a csv that one can read once data collection is finished.
# dat <- read_excel(text=GET("https://raw.githubusercontent.com/maxlindmark/scaling/master/data/metabolism_data.xlsx"))

dat <- read_excel("data/metabolism_data.xlsx")

glimpse(dat)

cols = c(1, 2, 3, 14, 15, 16, 17)
dat[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))

glimpse(dat)

unique(dat$species)


# B. FIT LOG-LOG MODEL BY SPECIES ==================================================
#** Coregonus albula ===============================================================
c_albula <- dat %>% 
  filter(species == "Coregonus albula")

sort(unique(c_albula$temp_c))

# 4C
d4 <- c_albula %>% filter(temp_c == 4)
d4
plot(log(d4$metabolic_rate) ~ log(d4$mass_g))
summary(lm(log(d4$metabolic_rate) ~ log(d4$mass_g)))

# 8C
d8 <- c_albula %>% filter(temp_c == 8)
d8
plot(log(d8$metabolic_rate) ~ log(d8$mass_g))
summary(lm(log(d8$metabolic_rate) ~ log(d8$mass_g)))

# 15C
d15 <- c_albula %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g)))


#** Coregonus fontanae =============================================================
c_fontanae <- dat %>% 
  filter(species == "Coregonus fontanae")

sort(unique(c_fontanae$temp_c))

# 4C
d4 <- c_fontanae %>% filter(temp_c == 4)
d4
plot(log(d4$metabolic_rate) ~ log(d4$mass_g))
summary(lm(log(d4$metabolic_rate) ~ log(d4$mass_g)))

# 8C
d8 <- c_fontanae %>% filter(temp_c == 8)
d8
plot(log(d8$metabolic_rate) ~ log(d8$mass_g))
summary(lm(log(d8$metabolic_rate) ~ log(d8$mass_g)))

# 15C
d15 <- c_fontanae %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g)))


#** Abramis brama ==================================================================
a_brama <- dat %>% 
  filter(species == "Abramis brama")

sort(unique(a_brama$temp_c))

# 5C
d5 <- a_brama %>% filter(temp_c == 5)
d5
plot(log(d5$metabolic_rate) ~ log(d5$mass_g))
summary(lm(log(d5$metabolic_rate) ~ log(d5$mass_g)))

# 10C
d10 <- a_brama %>% filter(temp_c == 10)
d10
plot(log(d10$metabolic_rate) ~ log(d10$mass_g))
summary(lm(log(d10$metabolic_rate) ~ log(d10$mass_g)))

# 15C
d15 <- a_brama %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g)))

# 20C
d20 <- a_brama %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g)))

# 23C
d23 <- a_brama %>% filter(temp_c == 23)
d23
plot(log(d23$metabolic_rate) ~ log(d23$mass_g))
summary(lm(log(d23$metabolic_rate) ~ log(d23$mass_g)))


#** Rutilus rutilus ================================================================
r_rutilus <- dat %>% 
  filter(species == "Rutilus rutilus")

sort(unique(r_rutilus$temp_c))

# 5C
d5 <- r_rutilus %>% filter(temp_c == 5)
d5
plot(log(d5$metabolic_rate) ~ log(d5$mass_g))
summary(lm(log(d5$metabolic_rate) ~ log(d5$mass_g)))

# 10C
d10 <- r_rutilus %>% filter(temp_c == 10)
d10
plot(log(d10$metabolic_rate) ~ log(d10$mass_g))
summary(lm(log(d10$metabolic_rate) ~ log(d10$mass_g)))

# 15C
d15 <- r_rutilus %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g)))

# 20C
d20 <- r_rutilus %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g)))

# 23C
d23 <- r_rutilus %>% filter(temp_c == 23)
d23
plot(log(d23$metabolic_rate) ~ log(d23$mass_g))
summary(lm(log(d23$metabolic_rate) ~ log(d23$mass_g)))


#** Salvelinus confluentus =========================================================
s_confluentus <- dat %>% 
  filter(species == unique(dat$species)[5]) # this is odd, can't write species name

sort(unique(s_confluentus$temp_c))

# 3C
d3 <- s_confluentus %>% filter(temp_c == 3)
d3
plot(log(d3$metabolic_rate) ~ log(d3$mass_g))
summary(lm(log(d3$metabolic_rate) ~ log(d3$mass_g))) 

# 7C
d7 <- s_confluentus %>% filter(temp_c == 7)
d7
plot(log(d7$metabolic_rate) ~ log(d7$mass_g))
summary(lm(log(d7$metabolic_rate) ~ log(d7$mass_g))) 

# 10C
d10 <- s_confluentus %>% filter(temp_c == 10)
d10
plot(log(d10$metabolic_rate) ~ log(d10$mass_g))
summary(lm(log(d10$metabolic_rate) ~ log(d10$mass_g))) # Not significant size effect

# 13C
d13 <- s_confluentus %>% filter(temp_c == 13)
d13
plot(log(d13$metabolic_rate) ~ log(d13$mass_g))
summary(lm(log(d13$metabolic_rate) ~ log(d13$mass_g))) # Not significant size effect

# 16C
d16 <- s_confluentus %>% filter(temp_c == 16)
d16
plot(log(d16$metabolic_rate) ~ log(d16$mass_g))
summary(lm(log(d16$metabolic_rate) ~ log(d16$mass_g))) # Not significant size effect

# 20C
d20 <- s_confluentus %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g))) # Not significant size effect


#** Oncorhynchus mykiss ============================================================
o_mykiss <- dat %>% 
  filter(species == "Oncorhynchus mykiss")

sort(unique(o_mykiss$temp_c)) # will round these temperatures..

o_mykiss$temp_c <- round(o_mykiss$temp_c)

sort(unique(o_mykiss$temp_c))

o_mykiss$temp_c <- ifelse(o_mykiss$temp_c == 10,
                          11,
                          o_mykiss$temp_c)

sort(unique(o_mykiss$temp_c))

# 5C
d5 <- o_mykiss %>% filter(temp_c == 5)
d5
plot(log(d5$metabolic_rate) ~ log(d5$mass_g))
summary(lm(log(d5$metabolic_rate) ~ log(d5$mass_g))) 

# 8C
d8 <- o_mykiss %>% filter(temp_c == 8)
d8
plot(log(d8$metabolic_rate) ~ log(d8$mass_g))
summary(lm(log(d8$metabolic_rate) ~ log(d8$mass_g))) 

# 11C
d11 <- o_mykiss %>% filter(temp_c == 11)
d11
plot(log(d11$metabolic_rate) ~ log(d11$mass_g))
summary(lm(log(d11$metabolic_rate) ~ log(d11$mass_g))) 

# 15C
d15 <- o_mykiss %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g))) 

# 20C
d20 <- o_mykiss %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g))) 


#** Salvelinus fontinalis ==========================================================
s_fontinalis <- dat %>% 
  filter(species == unique(dat$species)[7]) # again same issue with filtering on chr...

sort(unique(s_fontinalis$temp_c))

# 10C
d10 <- s_fontinalis %>% filter(temp_c == 10)
d10
plot(log(d10$metabolic_rate) ~ log(d10$mass_g))
summary(lm(log(d10$metabolic_rate) ~ log(d10$mass_g))) 

# 15C
d15 <- s_fontinalis %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g))) 

# 20C
d20 <- s_fontinalis %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g))) 


#** Salmo trutta ===================================================================
s_trutta <- dat %>% 
  filter(species == "Salmo trutta") 

sort(unique(s_trutta$temp_c))

# 10C
d10 <- s_trutta %>% filter(temp_c == 10)
d10
plot(log(d10$metabolic_rate) ~ log(d10$mass_g))
summary(lm(log(d10$metabolic_rate) ~ log(d10$mass_g))) 


#** Catostomus commersonii =========================================================
c_commersonii <- dat %>% 
  filter(species == "Catostomus commersonii")

sort(unique(c_commersonii$temp_c))

# 10C
d10 <- c_commersonii %>% filter(temp_c == 10)
d10
plot(log(d10$metabolic_rate) ~ log(d10$mass_g))
summary(lm(log(d10$metabolic_rate) ~ log(d10$mass_g))) 

# 15C
d15 <- c_commersonii %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g))) 

# 20C
d20 <- c_commersonii %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g))) 


#** Cyprinus carpio ================================================================
c_carpio <- dat %>% 
  filter(species == "Cyprinus carpio")

sort(unique(c_carpio$temp_c))

# 10C
d10 <- c_carpio %>% filter(temp_c == 10)
d10
plot(log(d10$metabolic_rate) ~ log(d10$mass_g))
summary(lm(log(d10$metabolic_rate) ~ log(d10$mass_g))) 

# 20C
d20 <- c_carpio %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g))) 

# 30C
d30 <- c_carpio %>% filter(temp_c == 30)
d30
plot(log(d30$metabolic_rate) ~ log(d30$mass_g))
summary(lm(log(d30$metabolic_rate) ~ log(d30$mass_g))) 

# 35C
d35 <- c_carpio %>% filter(temp_c == 35)
d35
plot(log(d35$metabolic_rate) ~ log(d35$mass_g))
summary(lm(log(d35$metabolic_rate) ~ log(d35$mass_g))) 


#** Ameiurus nebulosus =============================================================
a_nebulosus <- dat %>% 
  filter(species == "Ameiurus nebulosus")

sort(unique(a_nebulosus$temp_c))

# 10C
d10 <- a_nebulosus %>% filter(temp_c == 10)
d10
plot(log(d10$metabolic_rate) ~ log(d10$mass_g))
summary(lm(log(d10$metabolic_rate) ~ log(d10$mass_g))) 

# 20C
d20 <- a_nebulosus %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g))) 

# 30C
d30 <- a_nebulosus %>% filter(temp_c == 30)
d30
plot(log(d30$metabolic_rate) ~ log(d30$mass_g))
summary(lm(log(d30$metabolic_rate) ~ log(d30$mass_g))) 


#** Silurus meridionalis ===========================================================
s_meridionalis <- dat %>% 
  filter(species == "Silurus meridionalis")

sort(unique(s_meridionalis$temp_c))

# 10C
d10 <- s_meridionalis %>% filter(temp_c == 10)
d10
plot(log(d10$metabolic_rate) ~ log(d10$mass_g))
summary(lm(log(d10$metabolic_rate) ~ log(d10$mass_g))) 

# 15C
d15 <- s_meridionalis %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g))) 

# 20C
d20 <- s_meridionalis %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g))) 

# 25C
d25 <- s_meridionalis %>% filter(temp_c == 25)
d25
plot(log(d25$metabolic_rate) ~ log(d25$mass_g))
summary(lm(log(d25$metabolic_rate) ~ log(d25$mass_g))) 

# 30C
d30 <- s_meridionalis %>% filter(temp_c == 30)
d30
plot(log(d30$metabolic_rate) ~ log(d30$mass_g))
summary(lm(log(d30$metabolic_rate) ~ log(d30$mass_g))) 


#** Carassius auratus ==============================================================
c_auratus <- dat %>% 
  filter(species == "Carassius auratus")

sort(unique(c_auratus$temp_c))

# 10C
d10 <- c_auratus %>% filter(temp_c == 10)
d10
plot(log(d10$metabolic_rate) ~ log(d10$mass_g))
summary(lm(log(d10$metabolic_rate) ~ log(d10$mass_g))) 

# 20C
d20 <- c_auratus %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g))) 

# 30C
d30 <- c_auratus %>% filter(temp_c == 30)
d30
plot(log(d30$metabolic_rate) ~ log(d30$mass_g))
summary(lm(log(d30$metabolic_rate) ~ log(d30$mass_g))) 

# 35C
d35 <- c_auratus %>% filter(temp_c == 35)
d35
plot(log(d35$metabolic_rate) ~ log(d35$mass_g))
summary(lm(log(d35$metabolic_rate) ~ log(d35$mass_g))) 


#** Pomadasys commersonnii =========================================================
p_commersonnii <- dat %>% 
  filter(species == "Pomadasys commersonnii")

sort(unique(p_commersonnii$temp_c))

# 15C
d15 <- p_commersonnii %>% filter(temp_c == 15)
d15
plot(log(d15$metabolic_rate) ~ log(d15$mass_g))
summary(lm(log(d15$metabolic_rate) ~ log(d15$mass_g))) 

# 20C
d20 <- p_commersonnii %>% filter(temp_c == 20)
d20
plot(log(d20$metabolic_rate) ~ log(d20$mass_g))
summary(lm(log(d20$metabolic_rate) ~ log(d20$mass_g))) 

# 25C
d25 <- p_commersonnii %>% filter(temp_c == 25)
d25
plot(log(d25$metabolic_rate) ~ log(d25$mass_g))
summary(lm(log(d25$metabolic_rate) ~ log(d25$mass_g))) 


