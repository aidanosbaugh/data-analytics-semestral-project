# =============================================================================
# 03_analysis.R  --  Data Analytics semestral project, topic: Specialty Cafes
#
# Main analysis: research questions A1-A10 of the Specialty Cafes topic.
# Reads analysis/dataset.csv (written by 01_cleaning.R) and exports one table
# and one figure per analysis to report/tables and report/figures.
#
# Run from the REPOSITORY ROOT:  source("analysis/03_analysis.R")
# =============================================================================


# ---- Settings fixed in advance ----------------------------------------------

ALPHA <- 0.05          # significance level for every test
SEED  <- 2026          # seed for every bootstrap
B     <- 5000          # bootstrap resamples

# A1 benchmark: 4 US dollars converted to CZK.
# Exchange rate: Czech National Bank fixing of 5 November 2026,
# 1 USD = 21.646 CZK  ->  4 USD = 86.58 CZK.
USD_CZK   <- 21.646
BENCH_CZK <- 4 * USD_CZK


# ---- Hypotheses (written BEFORE looking at the data) -------------------------
#
# Direction rule: one-tailed only where the question names a direction
# ("charge more", "how much more"); two-tailed otherwise.
#
# A1  One sample. Is the average price of a filter coffee different from 4 USD?
#     H0: mu_filter  = 86.58 CZK
#     H1: mu_filter != 86.58 CZK                                  (two-tailed)
#
# A2  Two groups. Do cafes that offer filter coffee charge more for an espresso
#     than cafes that don't?  Difference = offers filter MINUS does not offer.
#     H0: mu_espresso(filter) - mu_espresso(no filter)  = 0
#     H1: mu_espresso(filter) - mu_espresso(no filter)  > 0       (one-tailed)
#
# A3  Paired. How much more does a filter coffee cost than an espresso in the
#     same cafe?  Difference = filter price MINUS espresso price, per cafe.
#     H0: mu_d  = 0
#     H1: mu_d  > 0                                               (one-tailed)
#
# A4  More groups. Does the espresso price differ between cafes with none, one
#     or both specialty features (roastery, filter coffee)?
#     H0: mu_none = mu_one = mu_both
#     H1: at least one group mean differs                         (omnibus)
#     Post hoc: all pairwise comparisons, adjusted for multiple testing.
#
# A5  Two categorical. Is a roastery more likely to offer filter coffee than
#     other cafes, and by how many percentage points?
#     H0: offering filter coffee is independent of being a roastery
#         (P(filter | roastery) - P(filter | other) = 0)
#     H1: the two are associated                                  (two-sided test)
#     The direction and size are reported as the difference in percentage
#     points, roastery MINUS other cafes.
#
# A6  Correlation. Is the filter coffee price correlated with the espresso
#     price, and in which direction and how strongly?
#     H0: rho  = 0
#     H1: rho != 0                                                (two-tailed)
#
# A7  Simple regression. By how many CZK does the filter coffee price change
#     per 10 CZK increase in the espresso price?
#     Model: filter_price = b0 + b1 * (espresso_price / 10) + e
#     H0: b1  = 0
#     H1: b1 != 0                                                 (two-tailed)
#
# A8  Multiple regression. Does the espresso price differ between cafes that
#     offer filter coffee and those that don't once you account for the
#     distance from the Orloj, and if so, how?
#     Model 1: espresso_price = b0 + b1 * offers_filter + e
#     Model 2: espresso_price = b0 + b1 * offers_filter + b2 * distance_km + e
#     H0: b1  = 0 in Model 2
#     H1: b1 != 0 in Model 2                                      (two-tailed)
#     Also: nested F-test of Model 2 vs Model 1, and the change in b1.
#
# A9  Interaction. Does the relationship between the espresso price and the
#     filter coffee price differ between independent cafes and chains?
#     Model: filter_price = b0 + b1*espresso_price + b2*chain
#                           + b3*espresso_price:chain + e
#     H0: b3  = 0
#     H1: b3 != 0                                                 (two-tailed)
#
# A10 Logistic regression. Does the probability that a cafe offers filter
#     coffee change with each additional km east of the Orloj, and by how
#     many percentage points?
#     Model: logit P(offers_filter = 1) = b0 + b1 * ew_km
#     H0: average marginal effect of ew_km  = 0
#     H1: average marginal effect of ew_km != 0                   (two-tailed)
#
# =============================================================================

# ---- A1 ----------------------------------------------------------------------
# (to be written after 02_descriptive.R has chosen the method)
