# EF_Case_1_Team_100
# Installing ggplot2 for data visualization
install.packages("ggplot2")

# Reading the dataset
df <- read.csv("data_case_1_group_100.csv")

# Part 1a
# Creating the excess returns from the stock return and risk-free rate
df$EXC_RET <- df$RET - df$RF

# Rescaling Unemp such that it has unit variance (dividing by its SD to make it work)
sd_unemp <- sqrt(var(df$Unemp, na.rm = TRUE))
df$Unemp_S <- df$Unemp / sd_unemp

# Quick check, was Unemp correctly rescaled? 
var(df$Unemp_S, na.rm = TRUE)

# Are there any outliers? 
boxplot(df$EXC_RET, main = "Excess Returns")

# Treating the outliers
q <- quantile(df$EXC_RET, probs = c(0.01, 0.99), na.rm = TRUE)
df$EXC_RET_wins <- pmin(pmax(df$EXC_RET, q[1]), q[2])

#Variables needed for summary stats from Model 3
vars_needed <- c("EXC_RET", "MktRF", "SMB", "HML", "RMW", "CMA",
                 "Inflation", "Unemp_S", "VIX", "EXC_RET_wins")

# Summary stats
summary_stats_raw <- data.frame(
  Variable = vars_needed,
  Mean = sapply(df[vars_needed], function(x) mean(x, na.rm = TRUE)),
  SD   = sapply(df[vars_needed], function(x) sd(x, na.rm = TRUE)),
  Min  = sapply(df[vars_needed], function(x) min(x, na.rm = TRUE)),
  Max  = sapply(df[vars_needed], function(x) max(x, na.rm = TRUE)),
  Obs  = sapply(df[vars_needed], function(x) sum(!is.na(x)))
)

# Clean version
summary_stats_clean <- subset(summary_stats_raw, Variable != "EXC_RET")

# Printing both tables
print(summary_stats_raw)    # goes to Appendix
print(summary_stats_clean)  # Table 1

# How many observations do we have of the VIX index? 
n_vix <- sum(!is.na(df$VIX))
cat("Number of VIX observations:", n_vix, "\n")

# Part 1b
# Creating the log of MktRF (exclude negatives/zeros to avoid NaN)
df$ln_MktRF <- ifelse(df$MktRF > 0, log(df$MktRF), NA)

# Scatter plot in Journal of Finance style
library(ggplot2)

ggplot(df, aes(x = ln_MktRF, y = EXC_RET)) +
  geom_point(alpha = 0.3, size = 0.7) +
  theme_minimal(base_size = 12) +
  labs(
    x = "ln(MktRF)",
    y = "Excess Return (EXC_RET, %)",
    title = "Figure 1. Scatter plot of Excess Returns vs. ln(MktRF)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )

#Part 1c
# Model 1: Estimate the FF 3 Factor model
model1 <- lm(EXC_RET ~ MktRF + SMB + HML, data = df)

# Model 2: Estimate the FF 5 Factor model
model2 <- lm(EXC_RET ~ MktRF + SMB + HML + RMW + CMA, data = df)

# Model 3: extend Model 2 by including Inflation, Unemp_S and the VIX as independent variables
model3 <- lm(EXC_RET ~ MktRF + SMB + HML + RMW + CMA + Inflation + Unemp_S + VIX, data = df)

# Goodness of fit statistics + more
get_gof_pack3 <- function(fit) {
  s   <- summary(fit)
  n   <- nobs(fit)                
  k   <- length(coef(fit))          
  rss <- sum(residuals(fit)^2)  
  mse <- rss / n              
  aic_pack3 <- log(mse) + (2 * k / n)
  bic_pack3 <- log(mse) + (k / n) * log(n)
  data.frame(
    N        = n,
    K        = k,
    R2       = s$r.squared,
    Adj_R2   = s$adj.r.squared,
    AIC_pack3 = aic_pack3,
    BIC_pack3 = bic_pack3
  )
}

# Adding the joint F-statistic
get_fstat <- function(fit){
  s <- summary(fit)
  Fv  <- as.numeric(s$fstatistic[1])
  df1 <- as.integer(s$fstatistic[2])
  df2 <- as.integer(s$fstatistic[3])
  data.frame(
    F_stat = Fv,
    F_df   = sprintf("(%d, %d)", df1, df2),
    F_pval = pf(Fv, df1, df2, lower.tail = FALSE)
  )
}

# Example output
gof_tab <- rbind(
  cbind(Model = "FF3 (Model 1)",get_gof_pack3(model1), get_fstat(model1)),
  cbind(Model = "FF5 (Model 2)",get_gof_pack3(model2), get_fstat(model2)),
  cbind(Model = "FF5 + Macro (Model 3)", get_gof_pack3(model3), get_fstat(model3))
)

print(gof_tab)

# Table formatted with rounded and aligned figures
gof_tab_clean <- data.frame(
  Model  = gof_tab$Model,
  N      = gof_tab$N,
  K      = gof_tab$K,
  R2     = round(gof_tab$R2, 3),
  Adj_R2 = round(gof_tab$Adj_R2, 3),
  F_stat = round(gof_tab$F_stat, 2),
  F_df   = gof_tab$F_df,
  F_pval = signif(gof_tab$F_pval, 3),
  AIC    = round(gof_tab$AIC_pack3, 5),   # Pack-3
  BIC    = round(gof_tab$BIC_pack3, 5)    # Pack-3
)

print(gof_tab_clean)
      
#Part 1d
s3 <- summary(model3)$coefficients

# Grab the two coefficients
b_vix <- s3["VIX","Estimate"]
p_vix <- s3["VIX","Pr(>|t|)"]
b_un  <- s3["Unemp_S","Estimate"]
p_un  <- s3["Unemp_S","Pr(>|t|)"]

# Converting Unemp_S to effect per +1 pp in unemployment for easier interpretation
sd_unemp <- sd(df$Unemp, na.rm = TRUE)
effect_per_1pp <- b_un / sd_unemp

# Printing the coefficients and p-values
print(b_vix)
print(p_vix)
print(b_un)
print(p_un)
print(effect_per_1pp)
