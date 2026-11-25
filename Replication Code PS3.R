# =============================================================================
# Replication of Fehr & Goette (2007), "Do Workers Work More if Wages Are High?
# Evidence from a Randomized Field Experiment", AER 97(1)
# Pure replication: Tables 1, 3, 5, 6 and Figure 1A.
# =============================================================================

library(haven); library(dplyr); library(sandwich)
library(lmtest); library(ggplot2); library(knitr)

dir.create("output", showWarnings = FALSE)

t14 <- read_dta("tables1to4.dta")   # Tables 1-4: messenger x four-week block
t56 <- read_dta("tables5to6.dta")   # Tables 5-6, Figure 1: messenger x day

# -----------------------------------------------------------------------------
# Table 1: mean revenues and shifts by group and period
# Group codes inferred from sample sizes (verify against paper's Table 1):
# 1 = Flash, 2 = non-participating Veloblitz, 3 = participants (control role),
# 4 = participants being treated
# -----------------------------------------------------------------------------
table1 <- t14 %>%
  group_by(group, block) %>%
  summarise(mean_rev = mean(totrev), sd_rev = sd(totrev),
            mean_shifts = mean(shifts), sd_shifts = sd(shifts),
            N = n(), .groups = "drop")
print(as.data.frame(table1), digits = 2)

# -----------------------------------------------------------------------------
# Table 3: effect of the higher commission (paper eq. 8)
# OLS with individual dummies, SEs clustered on messenger.
# Stata: areg ..., absorb(fahrer) cluster(fahrer); in R: lm + vcovCL with
# HC1 and cadjust = TRUE to mimic Stata's small-sample correction.
# -----------------------------------------------------------------------------
reg <- function(y, rhs, data, cl) {
  m  <- lm(as.formula(paste(y, "~", rhs, "+ factor(fahrer)")), data = data)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = data[[cl]],
                                  type = "HC1", cadjust = TRUE))
  c(estimate = ct["high", 1], se = ct["high", 2],
    R2 = summary(m)$r.squared, N = nobs(m))
}
d_part  <- filter(t14, maxhigh == 1)   # participating messengers
d_vebli <- filter(t14, vebli == 1)     # all Veloblitz messengers

res <- rbind(
  "(1) revenue, participants" = reg("totrev", "high+block2+block3", d_part,  "fahrer"),
  "(2) revenue, Veloblitz"    = reg("totrev", "high+block2+block3", d_vebli, "fahrer"),
  "(4) shifts, participants"  = reg("shifts", "high+block2+block3", d_part,  "fahrer"),
  "(5) shifts, Veloblitz"     = reg("shifts", "high+block2+block3", d_vebli, "fahrer"))
paper <- data.frame(paper_estimate = c(1033.6, 1094.6, 3.99, 4.08),
                    paper_se       = c(326.9, 297.8, 1.030, 0.942))
print(cbind(round(res, 3), paper))

# Columns (3)/(6): all messengers (incl. Flash) + dummy `experiment`.
# Paper: 1035.8 (444.7) and 3.44 (1.610).
# We get about 1090 (297) and 4.06 (0.94): does NOT match; other coefficients do.
print(round(rbind(
  "(3) revenue" = reg("totrev", "high+experiment+block2+block3", t14, "fahrer"),
  "(6) shifts"  = reg("shifts", "high+experiment+block2+block3", t14, "fahrer")), 3))

# Spillover coefficient on `experiment` (paper: -54.4, insignificant)
m3 <- lm(totrev ~ high + experiment + block2 + block3 + factor(fahrer), data = t14)
print(coeftest(m3, vcov = vcovCL(m3, cluster = t14$fahrer,
                                 type = "HC1", cadjust = TRUE))[
                                   c("high", "experiment", "block2", "block3"), ])

# -----------------------------------------------------------------------------
# Elasticities (paper p. 309, 311): (effect / mean in treatment period 1) / 0.25
# Paper: 1.12 (revenue) and 1.34 (shifts), using rounded effects.
# -----------------------------------------------------------------------------
rev_t1 <- d_part %>% filter(block == 2) %>% group_by(odd) %>%
  summarise(m_rev = mean(totrev), m_sh = mean(shifts))
print(c(
  base_rev = mean(rev_t1$m_rev),
  elasticity_rev = (res["(1) revenue, participants", "estimate"] /
                      mean(rev_t1$m_rev)) / 0.25,
  base_shifts = mean(rev_t1$m_sh),
  elasticity_shifts = (res["(4) shifts, participants", "estimate"] /
                         mean(rev_t1$m_sh)) / 0.25))

# -----------------------------------------------------------------------------
# Table 5: effort per shift (log revenue per shift, fixed shifts only)
# Day effects (datum), tenure (lnten). The do-file clusters on day (datum);
# only day-clustering reproduces the published SEs.
# -----------------------------------------------------------------------------
s  <- t56 %>% filter(fixed == 1, maxhigh == 1, fahrer != 27, tag != 0, tag != 6)
fs <- paste0("f", 2:38, collapse = " + ")   # messenger dummies f2-f38
m1 <- lm(lnum ~ high + lnten + female + factor(datum), data = s)
m2 <- lm(as.formula(paste("lnum ~ high + lnten +", fs, "+ factor(datum)")), data = s)

for (m in list(m1, m2)) {
  ct <- coeftest(m, vcov = vcovCL(m, cluster = s$datum, type = "HC1", cadjust = TRUE))
  print(round(ct[c("high", "lnten"), ], 4))
  cat("R2 =", round(summary(m)$r.squared, 3), " N =", nobs(m), "\n\n")
}
# Paper: col 1 high = -0.0642 (0.030); col 2 high = -0.0601 (0.030)

# -----------------------------------------------------------------------------
# Figure 1A: mean log daily revenue by group and period
# 9 obs (messengers 7 and 10) have no group code `odd`: dropped from the figure
# only (as in Stata's `table odd block`), kept in the regressions.
# -----------------------------------------------------------------------------
fig1a <- s %>%
  filter(!is.na(odd)) %>%
  group_by(odd, block) %>%
  summarise(m = mean(lnum), se = sd(lnum) / sqrt(n()), .groups = "drop") %>%
  mutate(group = ifelse(odd == 1, "Group A", "Group B")) %>%
  ggplot(aes(factor(block), m, colour = group, group = group)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.1) +
  scale_x_discrete(labels = c("Pre-period", "Treatment 1", "Treatment 2")) +
  labs(x = NULL, y = "log(daily revenues)")
print(fig1a)
ggsave("output/figure1a.png", fig1a, width = 6.5, height = 4, dpi = 300)

# -----------------------------------------------------------------------------
# Table 6: treatment effect by loss aversion (lottery experiment)
# Same specification as Table 5 col 2 (messenger + day effects), SEs on day.
# -----------------------------------------------------------------------------
fit_cl <- function(rhs, data) {
  m <- lm(as.formula(paste("lnum ~", rhs, "+", fs, "+ factor(datum)")), data = data)
  list(m = m, V = vcovCL(m, cluster = data$datum, type = "HC1", cadjust = TRUE))
}
cb1 <- fit_cl("high_la + high_not + lnten", s)
cb2 <- fit_cl("high_la1 + high_la2 + high_not0 + lnten", s)
print(round(coeftest(cb1$m, vcov = cb1$V)[c("high_la", "high_not"), ], 4))
print(round(coeftest(cb2$m, vcov = cb2$V)[c("high_la1", "high_la2", "high_not0"), ], 4))
# Paper col 1: -0.105 (0.046), -0.0273 (0.033)
# Paper col 2: -0.0853 (0.062), -0.12 (0.053), -0.027 (0.032)
# Our col-2 SEs are slightly larger (0.064, 0.055, 0.033): small unexplained difference.
