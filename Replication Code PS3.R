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

# exp_nt: dummy for non-treated Veloblitz messengers in the experimental periods
# (participants in their control period and non-participants), as in the paper's
# regressions (3) and (6) ("Dummy for nontreated at Veloblitz"):
t14 <- t14 %>% mutate(exp_nt = as.numeric(experiment) * (1 - as.numeric(high)))

# -----------------------------------------------------------------------------
# Table 1: mean revenues and shifts by group and period
# Group codes (checked against the paper's Table 1: means and N match):
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
paper <- data.frame(paper_estimate = c(1033.6, 1094.5, 3.99, 4.08),
                    paper_se       = c(326.9, 297.8, 1.030, 0.942))
print(cbind(round(res, 3), paper))

# Columns (3)/(6): all messengers (incl. Flash) + dummy for non-treated Veloblitz.
# Paper: 1035.8 (444.7) and 3.44 (1.610).
# Column (3) matches. Column (6): treatment effect 3.44 and R2 match, SE 1.59
# (paper 1.61); the other coefficients differ slightly.
print(round(rbind(
  "(3) revenue" = reg("totrev", "high+exp_nt+block2+block3", t14, "fahrer"),
  "(6) shifts"  = reg("shifts", "high+exp_nt+block2+block3", t14, "fahrer")), 3))

# Spillover coefficient on the non-treated dummy (paper: -54.4, insignificant)
m3 <- lm(totrev ~ high + exp_nt + block2 + block3 + factor(fahrer), data = t14)
print(coeftest(m3, vcov = vcovCL(m3, cluster = t14$fahrer,
                                 type = "HC1", cadjust = TRUE))[
                                   c("high", "exp_nt", "block2", "block3"), ])

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
  labs(x = NULL, y = "log(daily revenues)", colour = NULL) +
  theme_bw()
print(fig1a)
ggsave("output/figure1a.png", fig1a, width = 6.5, height = 4, dpi = 300)
ggsave("output/figure1a.pdf", fig1a, width = 6.5, height = 4)

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

# =============================================================================
# EXPORT: save the replication results as LaTeX tables (for Overleaf)
# =============================================================================
dir.create("tables", showWarnings = FALSE)

save_tex <- function(df, file, digits = 3, align = NULL) {
  tex <- kable(df, format = "latex", booktabs = TRUE, row.names = FALSE,
               digits = digits, linesep = "", align = align)
  writeLines(tex, file.path("tables", file))
}
fmt <- function(e, s, d = 3) {      # "estimate (se)"
  paste0(formatC(e, format = "f", digits = d), " (",
         formatC(s, format = "f", digits = d), ")")
}

# Table 1 -----------------------------------------------------------------
tab1 <- table1 %>% mutate(group = as.numeric(group), block = as.numeric(block))
names(tab1) <- c("Group", "Period", "Mean revenue", "SD revenue",
                 "Mean shifts", "SD shifts", "N")
save_tex(tab1, "table1.tex", digits = 2)

# Table 3 (columns 1, 2, 4, 5 and 3, 6) ---------------------------------------
r36 <- rbind(reg("totrev", "high+exp_nt+block2+block3", t14, "fahrer"),
             reg("shifts", "high+exp_nt+block2+block3", t14, "fahrer"))
all_res <- rbind(res, r36)
all_paper <- data.frame(est = c(paper$paper_estimate, 1035.8, 3.44),
                        se  = c(paper$paper_se, 444.7, 1.610))
d3 <- c(2, 2, 3, 3, 2, 3)           # digits per row
tab3 <- data.frame(
  Column = c("(1) Revenue, participants", "(2) Revenue, Veloblitz",
             "(4) Shifts, participants", "(5) Shifts, Veloblitz",
             "(3) Revenue, all messengers", "(6) Shifts, all messengers"),
  Ours  = unname(mapply(fmt, all_res[, "estimate"], all_res[, "se"], d3)),
  Paper = c("1033.6 (326.9)", "1094.5 (297.8)", "3.99 (1.030)",    # as published
            "4.08 (0.942)", "1035.8 (444.7)", "3.44 (1.610)"),
  R2    = round(all_res[, "R2"], 3),
  N     = all_res[, "N"])
tab3 <- tab3[c(1, 2, 5, 3, 4, 6), ] # order as in the paper: (1)-(6)
save_tex(tab3, "table3.tex", align = c("l", "r", "r", "r", "r"))

# Elasticities --------------------------------------------------------------
tab_el <- data.frame(
  Outcome  = c("Revenue", "Shifts"),
  Effect   = c(res["(1) revenue, participants", "estimate"],
               res["(4) shifts, participants", "estimate"]),
  Mean_T1  = c(mean(rev_t1$m_rev), mean(rev_t1$m_sh)),
  Ours     = c((res["(1) revenue, participants", "estimate"] / mean(rev_t1$m_rev)) / 0.25,
               (res["(4) shifts, participants", "estimate"] / mean(rev_t1$m_sh)) / 0.25),
  Paper    = c(1.12, 1.34))
names(tab_el) <- c("Outcome", "Effect", "Mean, treatment period 1",
                   "Elasticity (ours)", "Elasticity (paper)")
save_tex(tab_el, "elasticities.tex", digits = 2)

# Table 5 -------------------------------------------------------------------
ct5 <- lapply(list(m1, m2), function(m)
  coeftest(m, vcov = vcovCL(m, cluster = s$datum, type = "HC1", cadjust = TRUE)))
tab5 <- data.frame(
  Specification = c("(1) day effects, tenure, gender",
                    "(2) day effects, tenure, messenger effects"),
  Ours  = c(fmt(ct5[[1]]["high", 1], ct5[[1]]["high", 2], 4),
            fmt(ct5[[2]]["high", 1], ct5[[2]]["high", 2], 4)),
  Paper = c("-0.0642 (0.030)", "-0.0601 (0.030)"),
  R2    = c(summary(m1)$r.squared, summary(m2)$r.squared),
  N     = c(nobs(m1), nobs(m2)))
save_tex(tab5, "table5.tex", align = c("l", "r", "r", "r", "r"))

# Table 6 -------------------------------------------------------------------
ct61 <- coeftest(cb1$m, vcov = cb1$V)
ct62 <- coeftest(cb2$m, vcov = cb2$V)
v6 <- c("high_la", "high_not", "high_la1", "high_la2", "high_not0")
tab6 <- data.frame(
  Column = c("(1)", "(1)", "(2)", "(2)", "(2)"),
  Variable = c("Higher wage x loss averse", "Higher wage x not loss averse",
               "Higher wage x rejects one lottery", "Higher wage x rejects both lotteries",
               "Higher wage x not loss averse"),
  Ours = c(fmt(ct61["high_la", 1],  ct61["high_la", 2],  4),
           fmt(ct61["high_not", 1], ct61["high_not", 2], 4),
           fmt(ct62["high_la1", 1], ct62["high_la1", 2], 4),
           fmt(ct62["high_la2", 1], ct62["high_la2", 2], 4),
           fmt(ct62["high_not0", 1], ct62["high_not0", 2], 4)),
  Paper = c("-0.105 (0.046)", "-0.0273 (0.033)", "-0.0853 (0.062)",
            "-0.120 (0.053)", "-0.027 (0.032)"))
save_tex(tab6, "table6.tex", align = c("l", "l", "r", "r"))

# Table 3 in the layout of the paper (all six columns, stars, fixed rows) ----

reg_full <- function(y, rhs, data) {
  m  <- lm(as.formula(paste(y, "~", rhs, "+ factor(fahrer)")), data = data)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = data$fahrer,
                                  type = "HC1", cadjust = TRUE))
  list(ct = ct, G = length(unique(data$fahrer)),
       R2 = summary(m)$r.squared, N = nobs(m))
}
cell <- function(fit, var, d_est, d_se) {
  if (!(var %in% rownames(fit$ct))) return(c("", ""))
  e <- fit$ct[var, 1]; s <- fit$ct[var, 2]
  p <- 2 * pt(-abs(e / s), df = fit$G - 1)
  stars <- ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
  c(paste0(formatC(e, format = "f", big.mark = ",", digits = d_est), stars),
    paste0("(", formatC(s, format = "f", digits = d_se), ")"))
}
fits <- list(
  reg_full("totrev", "high+block2+block3",            d_part),
  reg_full("totrev", "high+block2+block3",            d_vebli),
  reg_full("totrev", "high+exp_nt+block2+block3", t14),
  reg_full("shifts", "high+block2+block3",            d_part),
  reg_full("shifts", "high+block2+block3",            d_vebli),
  reg_full("shifts", "high+exp_nt+block2+block3", t14))
digs <- list(c(1, 1), c(1, 1), c(1, 1), c(2, 3), c(2, 3), c(2, 3))
rws  <- list(c("Treatment dummy", "high"),
             c("Dummy for nontreated at Veloblitz", "exp_nt"),
             c("Treatment period 1", "block2"),
             c("Treatment period 2", "block3"))
body <- rbind(
  c("Dependent variable", rep("Revenue", 3), rep("Shifts", 3)),
  c("Sample", rep(c("Participants", "Veloblitz", "All messengers"), 2)))
for (r in rws) {
  cc <- mapply(function(f, d) cell(f, r[2], d[1], d[2]), fits, digs)
  body <- rbind(body, c(r[1], cc[1, ]), c("", cc[2, ]))
}
body <- rbind(body,
              c("Individual fixed effects", rep("Yes", 6)),
              c("R squared", sapply(fits, function(f) formatC(f$R2, format = "f", digits = 3))),
              c("N", sapply(fits, function(f) as.character(f$N))))
body[, -1] <- gsub("^-", "$-$", body[, -1])
tab3p <- as.data.frame(body, stringsAsFactors = FALSE)
names(tab3p) <- c("", paste0("(", 1:6, ")"))
writeLines(kable(tab3p, format = "latex", booktabs = TRUE, row.names = FALSE,
                 escape = FALSE, linesep = "", align = c("l", rep("r", 6))),
           "tables/table3_paperstyle.tex")
print(tab3p)

# Software versions (for the README / reproducibility checklist) ------------
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")

cat("LaTeX tables written to the folder 'tables'.\n")

# =============================================================================
# Extension: randomization inference for Table 3, column (1)
# Reassign the A/B labels across the 42 participants (group sizes kept),
# rebuild the treatment dummy and re-estimate. Two-sided p-value.
# =============================================================================
set.seed(1)
p   <- filter(t14, maxhigh == 1)
ids <- p %>% group_by(fahrer) %>% summarise(odd = first(odd))   # A = odd 1

delta_hat <- function(odd_new) {
  g <- odd_new[match(p$fahrer, ids$fahrer)]
  # treated in block 2 if group A, in block 3 if group B
  h <- as.integer((g == 1 & p$block == 2) | (g == 0 & p$block == 3))
  coef(lm(p$totrev ~ h + p$block2 + p$block3 + factor(p$fahrer)))["h"]
}
obs  <- delta_hat(ids$odd)
sims <- replicate(2000, delta_hat(sample(ids$odd)))
print(c(observed     = unname(obs),
        null_sd      = sd(sims),
        p_value      = mean(abs(sims) >= abs(obs)),
        p_value_cons = (1 + sum(abs(sims) >= abs(obs))) / (1 + length(sims))))
# Expected: observed 1033.56, null_sd 310.82, p_value 0.0005

ri_plot <- ggplot(data.frame(sims = sims), aes(sims)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "white") +
  geom_vline(xintercept = obs, colour = "red", linewidth = 1) +
  labs(x = "Estimated treatment effect (CHF)", y = "Frequency")
print(ri_plot)
ggsave("output/figure_ri.png", ri_plot, width = 6.5, height = 4, dpi = 300)

# =============================================================================
# Robustness 1: Table 5 with SEs clustered on messenger instead of day
# =============================================================================
for (m in list(m1, m2)) {
  ct <- coeftest(m, vcov = vcovCL(m, cluster = s$fahrer, type = "HC1", cadjust = TRUE))
  print(round(ct["high", ], 4))
}
# Expected: -0.0642 (0.0271); -0.0601 (0.0301)

# =============================================================================
# Robustness 2: Table 5 without messengers 7 and 10 (no group code)
# =============================================================================
s_rob <- s %>% filter(!fahrer %in% c(7, 10))
r1 <- lm(lnum ~ high + lnten + female + factor(datum), data = s_rob)
r2 <- lm(as.formula(paste("lnum ~ high + lnten +", fs, "+ factor(datum)")), data = s_rob)
for (m in list(r1, r2)) {
  ct <- coeftest(m, vcov = vcovCL(m, cluster = s_rob$datum, type = "HC1", cadjust = TRUE))
  print(round(ct["high", ], 4))
  cat("N =", nobs(m), "\n\n")
}
# Expected: -0.0697 (0.0299); -0.0599 (0.0297); N = 1128
