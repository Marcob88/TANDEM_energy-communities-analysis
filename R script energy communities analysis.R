############################################################
# ECONOMETRIC ANALYSIS SCRIPT
# Enabling conditions for Energy Communities (NUTS-2)
# Assumes dataset is already cleaned and final
############################################################

# Required packages
library(dplyr) #packageVersion 1.1.4
library(lmtest) #packageVersion 0.9.40
library(sandwich)#packageVersion 3.1.1
library(car)#packageVersion 3.1.2
library(lme4)#packageVersion 1.1.37

############################################################
# DATA ASSUMPTION
############################################################

# The dataset is assumed to be provided as:
# filtered_data
# containing all variables used in the analysis

# Country identifier (NUTS prefix)

filtered_data$country_code <- substr(filtered_data$NUTS_ID, 1, 2)

############################################################
# VARIABLES USED IN THE ANALYSIS
############################################################

# Dependent variable
# EC_hab_2.1            : Energy communities per 100,000 inhabitants

# Core controls
# avg_POP_D             : Population density
# avg_EDUC              : Share of population with tertiary education
# Avg_PAT               : Patent applications per billion GDP
# avg_HH_Income         : Household net disposable income per capita

# Financial variables
# ERDF_Paym_PPS_PC      : ERDF payments per capita (PPS)
# avg_RD_total          : Total R&D expenditure (% GDP)
# Avg_RD_exp_private    : Business R&D expenditure (% GDP)

# Infrastructure variables
# Nr_facilities_Ren_cap : Renewable energy facilities per capita
# Nr_facilities_NonRen_cap : Non-renewable energy facilities per capita


# Institutional variables
# quality_of_government
# political_trust
# generalised_trust


############################################################
# INSTITUTIONAL INDEX (PCA)
############################################################

pca_vars <- filtered_data %>%
  select(quality_of_government, political_trust, generalised_trust)

pca_inst <- prcomp(pca_vars, scale. = TRUE)

# First principal component (≈80% variance explained)
filtered_data$inst_index_PC1 <- -pca_inst$x[, 1]


############################################################
# H1 – FINANCIAL CONDITIONS
############################################################

# Model 1.1 – Baseline financial effects
lm1.1 <- lm(
  log(EC_hab_2.1) ~
    log(avg_POP_D) +
    log(avg_EDUC) +
    log(Avg_PAT + 0.1) +
    log(avg_HH_Income) +
    log(ERDF_Paym_PPS_PC + 0.1) +
    log(avg_RD_total) +
    log(Avg_RD_exp_private + 0.01),
  data = filtered_data
)

# HC2 robust SE (influential points)
coeftest(lm1.1, vcov = vcovHC(lm1.1, type = "HC2"))

# Country-clustered SE (main specification)
coeftest(lm1.1, vcov = vcovCL(lm1.1, cluster = ~ country_code))

# Multicollinearity
vif(lm1.1)

# Random effects
lm1.1.RE <- lmer(
  update(formula(lm1.1), . ~ . + (1 | country_code)),
  data = filtered_data
)

# Fixed effects
lm1.1.FE <- lm(
  update(formula(lm1.1), . ~ . + factor(country_code)),
  data = filtered_data
)

############################################################
# Model 1.2 – ERDF × Macro-regions
############################################################

lm1.2 <- lm(
  log(EC_hab_2.1) ~
    log(avg_POP_D) +
    log(avg_EDUC) +
    log(Avg_PAT + 0.1) +
    log(avg_HH_Income) +
    log(ERDF_Paym_PPS_PC + 0.1) * European_region +
    log(avg_RD_total) +
    log(Avg_RD_exp_private + 0.01),
  data = filtered_data
)

coeftest(lm1.2, vcov = vcovHC(lm1.2, type = "HC2"))

vif(lm1.2)

############################################################
# Model 1.3 – ERDF × Institutional quality
############################################################

lm1.3 <- lm(
  log(EC_hab_2.1) ~
    log(avg_POP_D) +
    log(avg_EDUC) +
    log(Avg_PAT + 0.1) +
    log(avg_HH_Income) +
    log(ERDF_Paym_PPS_PC + 0.1) * inst_index_PC1 +
    log(avg_RD_total) +
    log(Avg_RD_exp_private + 0.01),
  data = filtered_data
)

coeftest(lm1.3, vcov = vcovHC(lm1.3, type = "HC2"))
coeftest(lm1.3, vcov = vcovCL(lm1.3, cluster = ~ country_code))
vif(lm1.3)

lm1.3.RE <- lmer(
  update(formula(lm1.3), . ~ . + (1 | country_code)),
  data = filtered_data
)

lm1.3.FE <- lm(
  update(formula(lm1.3), . ~ . + factor(country_code)),
  data = filtered_data
)

############################################################
# H2 – ENERGY INFRASTRUCTURE
############################################################

lm2.1 <- lm(
  log(EC_hab_2.1) ~
    log(avg_POP_D) +
    log(avg_EDUC) +
    log(Avg_PAT + 0.1) +
    log(avg_HH_Income) +
    log(Nr_facilities_NonRen_cap + 0.1) +
    log(Nr_facilities_Ren_cap + 0.1),
  data = filtered_data
)

coeftest(lm2.1, vcov = vcovHC(lm2.1, type = "HC2"))
coeftest(lm2.1, vcov = vcovCL(lm2.1, cluster = ~ country_code))
vif(lm2.1)

lm2.1.RE <- lmer(
  update(formula(lm2.1), . ~ . + (1 | country_code)),
  data = filtered_data
)

lm2.1.FE <- lm(
  update(formula(lm2.1), . ~ . + factor(country_code)),
  data = filtered_data
)

############################################################
# H3 – INSTITUTIONAL CONTEXT
############################################################

lm3.1 <- lm(
  log(EC_hab_2.1) ~
    log(avg_POP_D) +
    log(avg_EDUC) +
    log(Avg_PAT + 0.1) +
    log(avg_HH_Income) +
    inst_index_PC1,
  data = filtered_data
)

coeftest(lm3.1, vcov = vcovHC(lm3.1, type = "HC2"))
coeftest(lm3.1, vcov = vcovCL(lm3.1, cluster = ~ country_code))
vif(lm3.1)

lm3.1.RE <- lmer(
  update(formula(lm3.1), . ~ . + (1 | country_code)),
  data = filtered_data
)

lm3.1.FE <- lm(
  update(formula(lm3.1), . ~ . + factor(country_code)),
  data = filtered_data
)

############################################################
# Model 3.2 – Institutional quality × macro-regions
############################################################

lm3.2 <- lm(
  log(EC_hab_2.1) ~
    log(avg_POP_D) +
    log(avg_EDUC) +
    log(Avg_PAT + 0.1) +
    log(avg_HH_Income) +
    inst_index_PC1 * European_region,
  data = filtered_data
)

coeftest(lm3.2, vcov = vcovHC(lm3.2, type = "HC2"))
vif(lm3.2)

############################################################
# END OF SCRIPT
############################################################