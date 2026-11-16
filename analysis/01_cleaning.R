# =============================================================================
# 01_cleaning.R  --  Data Analytics semestral project, topic: Specialty Cafes
#
# What it does: reads the raw collection sheet (analysis/raw_data.csv, never
# edited by hand), cleans it step by step, creates the variables the research
# questions need, and writes analysis/dataset.csv. It also writes a short
# cleaning log (report/tables/cleaning_table.tex) for Section 2.1 of the report.
#
# How to run: from the REPOSITORY ROOT (e.g. open the RStudio project there).
#   source("analysis/01_cleaning.R")
#
# Packages: none. Base R only, so it runs on any fresh R installation.
#
# Steps follow Appendix B.1 of the guidelines:
#   0. Read the raw file
#   1. Inspect: structure, empty rows, duplicated rows, same cafe entered twice
#   2. Categories: district, cafe type, payment, roastery
#   3. Explicit decisions: cells with several prices, NA vs blank filter price
#   4. Text -> numbers (thousands separator, "NA" text)
#   5. Create variables: offers_filter, specialty_features, is_chain, ew_km ...
#   6. Check the result and write dataset.csv
# =============================================================================


# ---- 0. Paths and reading ---------------------------------------------------

raw_path <- file.path("analysis", "raw_data.csv")
out_path <- file.path("analysis", "dataset.csv")
tab_dir  <- file.path("report", "tables")

if (!file.exists(raw_path)) {
  stop("Cannot find ", raw_path, ". Run this script from the repository root.")
}
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# colClasses = "character": read EVERYTHING as text first, so R does not guess
#   types and silently turn "1,020" or "85, 95" into NA.
# na.strings = character(0): do NOT turn the text "NA" into a missing value yet.
#   In this sheet "NA" means "the cafe does not serve this drink", which is
#   different from a blank cell ("nobody wrote anything"). We decide later.
raw <- read.csv(raw_path, colClasses = "character", na.strings = character(0),
                check.names = FALSE, fileEncoding = "UTF-8")

cat("Raw file:", nrow(raw), "rows x", ncol(raw), "columns\n")
print(names(raw))   # the original (long) column names, for the record

# Short, R-friendly names, in the order of the sheet. "_raw" marks columns we
# will convert later; the converted version gets a new name (never overwrite).
short_names <- c("collector", "visit_time", "district_raw", "district_name",
                 "cafe_name", "address", "coords", "gmaps_link",
                 "espresso_raw", "cappuccino_raw", "flatwhite_raw", "filter_raw",
                 "payment_raw", "seats_in_raw", "seats_out_raw",
                 "rating_raw", "n_ratings_raw", "type_raw", "roastery_raw",
                 "photo_link", "note", "auto_separator", "distance_raw",
                 "lat_auto", "lon_auto", "adj_rating_auto")
stopifnot(ncol(raw) == length(short_names))
names(raw) <- short_names

# Row number in the raw file (1 = first data row under the header). It never
# changes because the raw file is never edited, so we can refer to rows by it.
raw$raw_row <- seq_len(nrow(raw))

# Cleaning log: every step and how many cafes (rows) it affected.
clean_log <- data.frame(step = character(0), cafes = integer(0))
log_step <- function(step, n) {
  clean_log <<- rbind(clean_log, data.frame(step = step, cafes = as.integer(n)))
  cat(sprintf("[log] %-60s %d\n", step, n))
}
log_step("Rows in the raw file", nrow(raw))

# Helpers used several times below
is_blank   <- function(x) is.na(x) | x == ""
is_na_text <- function(x) !is.na(x) & toupper(x) == "NA"        # "NA", "na"


# ---- 1. Inspect -------------------------------------------------------------

str(raw)

# 1a. Stray spaces ("Prague 7 ", "Joseph Muchnick "). Trim them in a copy and
#     count how many rows had at least one changed cell.
d <- raw
text_cols <- setdiff(names(d), "raw_row")
d[text_cols] <- lapply(d[text_cols], trimws)
rows_trimmed <- sum(apply(d[text_cols] != raw[text_cols], 1, any))
log_step("Leading or trailing spaces removed (rows)", rows_trimmed)

# 1b. Empty rows. At the bottom of the sheet there are rows with no collector
#     and no cafe, only the default "Cash & Card" and spreadsheet "#VALUE!".
empty_row <- is_blank(d$collector) & is_blank(d$cafe_name)
cat("\nNon-blank cells per column in the empty rows (what is left in them):\n")
print(colSums(!sapply(d[empty_row, text_cols], is_blank)))
d <- d[!empty_row, ]
log_step("Empty rows dropped", sum(empty_row))

# 1c. Exactly duplicated rows (identical in every column)
exact_dup <- duplicated(d[text_cols])
cat("\nExactly duplicated rows:", sum(exact_dup), "\n")
d <- d[!exact_dup, ]
log_step("Exactly duplicated rows dropped", sum(exact_dup))

# 1d. Coordinates: one cell "50.0836433, 14.4226154" -> latitude, longitude.
#     (Needed here already, to find cafes entered twice by different students.)
coord_parts <- strsplit(d$coords, ",")
d$lat <- as.numeric(sapply(coord_parts, `[`, 1))
d$lon <- as.numeric(sapply(coord_parts, `[`, 2))
cat("\nCoordinates that failed to parse:", sum(is.na(d$lat) | is.na(d$lon)), "\n")
print(summary(d[c("lat", "lon")]))   # Prague: lat ~ 50.0-50.2, lon ~ 14.2-14.7

# 1e. The same cafe entered twice by two students. These are NOT identical rows
#     (names are spelled differently, prices sometimes differ), so 1c misses them.
#     Step 1: list every pair of cafes less than 50 m apart (candidates).
#     Step 2: check the candidates by eye (name + address) and confirm by hand.
#     Several candidates are different cafes at the same address (e.g. The Cup,
#     Atomic matcha and Imaginarni zmrzlina share one building; Cafe DECADA and
#     DECAGRAM are separate buildings according to the collector's note), so an
#     automatic rule would merge cafes that are not the same.
pairs  <- t(combn(nrow(d), 2))
metres <- 1000 * sqrt(((d$lat[pairs[, 1]] - d$lat[pairs[, 2]]) * 111)^2 +
                      ((d$lon[pairs[, 1]] - d$lon[pairs[, 2]]) * 71)^2)
close  <- !is.na(metres) & metres < 50
candidates <- data.frame(
  row_a  = d$raw_row[pairs[close, 1]], cafe_a = d$cafe_name[pairs[close, 1]],
  row_b  = d$raw_row[pairs[close, 2]], cafe_b = d$cafe_name[pairs[close, 2]],
  metres = round(metres[close]))
candidates <- candidates[order(candidates$metres), ]
cat("\nCandidate duplicates (cafes less than 50 m apart):\n")
print(candidates, row.names = FALSE)

# Confirmed by hand: same name (up to spelling) AND same address.
# Each line: raw_row of the two entries.
dup_pairs <- rbind(
  c(  2, 144),   # Starbucks Reserve, Jungmannova 750/34
  c(  9, 103),   # Cafe Louvre, Narodni 22
  c( 44,  97),   # Caffe Vescovi, Ujezd 429
  c( 47, 158),   # Cafe Savoy, Vitezna 124/5
  c( 49, 159),   # Kampa Cafe, Zborovska 620/47
  c( 57, 196),   # Caffe Milani, Jungmannovo nam. 7
  c( 94, 142),   # Duo cafe & bistro, Plzenska 182/179
  c(102, 198),   # Cafe-Cafe, Rytirska 10
  c(106, 200),   # Republica Coffee, Martinska 418/4
  c(110, 197),   # Mondieu, Rytirska 403/16
  c(111, 131)    # No26 Espresso, Rytirska 398
)
# Safety check: every confirmed pair must be one of the candidates above.
stopifnot(all(paste(dup_pairs[, 1], dup_pairs[, 2]) %in%
              paste(candidates$row_a, candidates$row_b)))

# Which entry to keep? DECISION: keep the more complete entry, i.e. the one with
# more of the 10 data fields filled in (a blank cell or "NA" counts as not filled).
# If both are equally complete, keep the one entered first in the sheet.
# Reason: a price written down (with a menu photo as evidence) is positive
# evidence that the cafe sells the drink, while "NA" may be an oversight. This
# matters for Kampa Cafe and Cafe Savoy, where one student wrote a filter price
# and the other wrote NA.
data_fields <- c("espresso_raw", "cappuccino_raw", "flatwhite_raw", "filter_raw",
                 "seats_in_raw", "seats_out_raw", "rating_raw", "n_ratings_raw",
                 "type_raw", "roastery_raw")
d$completeness <- rowSums(!sapply(d[data_fields],
                                  function(x) is_blank(x) | is_na_text(x)))

drop_rows <- apply(dup_pairs, 1, function(p) {
  comp_a <- d$completeness[d$raw_row == p[1]]
  comp_b <- d$completeness[d$raw_row == p[2]]
  if (comp_b > comp_a) p[1] else p[2]     # tie -> drop the later entry
})
dup_report <- data.frame(
  cafe      = d$cafe_name[match(dup_pairs[, 1], d$raw_row)],
  row_a     = dup_pairs[, 1], complete_a = d$completeness[match(dup_pairs[, 1], d$raw_row)],
  row_b     = dup_pairs[, 2], complete_b = d$completeness[match(dup_pairs[, 2], d$raw_row)],
  dropped   = drop_rows)
cat("\nDuplicate cafes, completeness of each entry and the row dropped:\n")
print(dup_report, row.names = FALSE)

d <- d[!d$raw_row %in% drop_rows, ]
log_step("Same cafe entered twice: second entry dropped", length(drop_rows))
log_step("Cafes after removing empty and duplicate rows", nrow(d))


# ---- 2. Categories ----------------------------------------------------------

# 2a. District: "prague 1", "Prage 7", "prauge 1" -> "Prague 1"
cat("\nDistrict labels before:\n"); print(table(d$district_raw, useNA = "ifany"))
d$district <- sub("^pr[a-z]*e\\s*", "Prague ", tolower(d$district_raw))
cat("District labels after:\n");  print(table(d$district, useNA = "ifany"))
stopifnot(all(grepl("^Prague [0-9]+$", d$district)))
log_step("District labels harmonised", sum(d$district != d$district_raw))

# 2b. Cafe type: one blank cell -> missing (NA), not a fourth category.
cat("\nCafe type before:\n"); print(table(d$type_raw, useNA = "ifany"))
d$cafe_type <- ifelse(is_blank(d$type_raw), NA, d$type_raw)
d$cafe_type <- factor(d$cafe_type,
                      levels = c("Independent", "Local Chain", "International Chain"))
stopifnot(sum(is.na(d$cafe_type)) == sum(is_blank(d$type_raw)))  # no typo lost
log_step("Blank cafe type set to missing", sum(is_blank(d$type_raw)))

# 2c. Payment: already consistent (dropdown), only checked.
cat("\nPayment:\n"); print(table(d$payment_raw, useNA = "ifany"))
d$payment <- d$payment_raw

# 2d. Roastery 0/1: must be exactly 0 or 1, never blank.
cat("\nIs roastery:\n"); print(table(d$roastery_raw, useNA = "ifany"))
stopifnot(all(d$roastery_raw %in% c("0", "1")))
d$is_roastery <- as.integer(d$roastery_raw)


# ---- 3. Explicit decisions: prices ------------------------------------------

# First look at every price cell that is not a plain number.
price_cols <- c("espresso_raw", "cappuccino_raw", "flatwhite_raw", "filter_raw")
for (col in price_cols) {
  x <- d[[col]]
  cat("\n", col, ": values that are not a plain number:\n", sep = "")
  print(table(x[!grepl("^[0-9]+(\\.[0-9]+)?$", x)], useNA = "ifany"))
}

# DECISION (default rule of the course): a cell with several prices
# ("85, 95, 105" = sizes, "55/30" = dine-in/take-away, "75/90" = filter blends)
# -> keep the FIRST listed price. The first price is the smallest size or the
# dine-in price, so it is the most comparable across cafes.
# "NA"/"na" -> the cafe does not serve the drink -> price is missing (NA).
# Blank -> nothing was recorded -> price is missing (NA).
first_price <- function(x) {
  out    <- rep(NA_real_, length(x))
  served <- !(is_blank(x) | is_na_text(x))
  out[served] <- as.numeric(sub("^([0-9]+(\\.[0-9]+)?).*$", "\\1", x[served]))
  out
}
d$espresso_price   <- first_price(d$espresso_raw)
d$cappuccino_price <- first_price(d$cappuccino_raw)
d$flatwhite_price  <- first_price(d$flatwhite_raw)
d$filter_price     <- first_price(d$filter_raw)

# Check visually: the rows with several prices, before and after.
for (col in price_cols) {
  new_col <- sub("_raw$", "_price", col)
  multi   <- grepl("[,/]", d[[col]])
  cat("\n", col, " -> ", new_col, " (cells with several prices):\n", sep = "")
  print(d[multi, c("cafe_name", col, new_col)], row.names = FALSE)
  log_step(paste0("Several prices, first kept: ", sub("_raw$", "", col)), sum(multi))
}

# The filter coffee price has THREE meanings, and our topic depends on them:
#   a price   -> the cafe offers filter coffee         (offers_filter = 1)
#   "NA"/"na" -> the cafe does not offer filter coffee (offers_filter = 0)
#   blank     -> not recorded, we do not know          (offers_filter = NA)
# DECISION: a blank is NOT treated as "does not offer", because the guidelines
# told collectors to write NA, never leave a blank; a blank means the
# information is missing.
d$offers_filter <- ifelse(is_blank(d$filter_raw), NA,
                          ifelse(is_na_text(d$filter_raw), 0L, 1L))
cat("\nFilter coffee: offered (1), not offered (0), unknown (NA):\n")
print(table(d$offers_filter, useNA = "ifany"))
log_step("Filter price blank: offers filter set to missing",
         sum(is_blank(d$filter_raw)))
log_step("Flat white marked NA (not served)", sum(is_na_text(d$flatwhite_raw)))


# ---- 4. Text -> numbers -----------------------------------------------------

# Number of ratings: some are written with a thousands separator ("1,020").
bad <- !grepl("^[0-9]*$", d$n_ratings_raw)
cat("\nNumber of ratings with characters other than digits:\n")
print(d[bad, c("cafe_name", "n_ratings_raw")], row.names = FALSE)
d$n_ratings <- as.numeric(ifelse(is_blank(d$n_ratings_raw), NA,
                                 gsub(",", "", d$n_ratings_raw)))
cat("After conversion (same rows):\n")
print(d[bad, c("cafe_name", "n_ratings_raw", "n_ratings")], row.names = FALSE)
log_step("Thousands separator removed (number of ratings)", sum(bad))

# Rating and seats: blank or "NA" -> NA (never 0).
to_number <- function(x) as.numeric(ifelse(is_blank(x) | is_na_text(x), NA, x))
d$rating    <- to_number(d$rating_raw)
d$seats_in  <- to_number(d$seats_in_raw)
d$seats_out <- to_number(d$seats_out_raw)

# Distance from the Orloj in km: already computed in the sheet, only converted.
d$distance_km <- as.numeric(d$distance_raw)


# ---- 5. Create the variables of the research questions ----------------------

# A4: number of specialty features = roastery (0/1) + offers filter (0/1),
#     so 0 = none, 1 = one, 2 = both. NA if offers_filter is unknown.
d$specialty_features <- d$is_roastery + d$offers_filter
d$specialty_group    <- factor(d$specialty_features, levels = 0:2,
                               labels = c("None", "One", "Both"))

# A9: independent cafes vs chains (local + international chains together).
d$is_chain    <- ifelse(is.na(d$cafe_type), NA,
                        as.integer(d$cafe_type != "Independent"))
d$chain_group <- factor(d$is_chain, levels = 0:1,
                        labels = c("Independent", "Chain"))

# A10 (and footnote *** of Appendix C): distances from the Orloj in km.
# Positive ew_km = east of the Orloj, positive ns_km = north of it.
ORLOJ_LAT <- 50.0871
ORLOJ_LON <- 14.4207
d$ew_km <- (d$lon - ORLOJ_LON) * 71
d$ns_km <- (d$lat - ORLOJ_LAT) * 111


# ---- 6. Check the result and write ------------------------------------------

keep <- c("raw_row", "collector", "cafe_name", "district", "cafe_type",
          "payment", "espresso_price", "cappuccino_price", "flatwhite_price",
          "filter_price", "offers_filter", "is_roastery", "specialty_features",
          "specialty_group", "is_chain", "chain_group", "seats_in", "seats_out",
          "rating", "n_ratings", "distance_km", "lat", "lon", "ew_km", "ns_km")
clean <- d[keep]

cat("\nRows: raw file", nrow(raw), "-> cleaned", nrow(clean), "\n")
cat("\nMissing values per variable:\n")
print(colSums(is.na(clean)))
cat("\nRanges of the numeric variables:\n")
print(summary(clean[sapply(clean, is.numeric)]))

# Plausibility checks: the script stops if any of these fails.
prices <- c("espresso_price", "cappuccino_price", "flatwhite_price", "filter_price")
stopifnot(
  all(clean[prices] > 0 & clean[prices] < 400, na.rm = TRUE),  # CZK, plausible
  all(clean$rating >= 1 & clean$rating <= 5,  na.rm = TRUE),
  all(clean$n_ratings >= 0, na.rm = TRUE),
  all(clean$seats_in >= 0 & clean$seats_out >= 0, na.rm = TRUE),
  all(clean$offers_filter == 1 | is.na(clean$filter_price)),   # price <=> offered
  all(is.na(clean$offers_filter) | clean$offers_filter == 0 |
        !is.na(clean$filter_price)),
  all(abs(clean$ew_km) < 20 & abs(clean$ns_km) < 20)            # inside Prague
)

log_step("Cafes in the analysis dataset", nrow(clean))

# write.csv: no row numbers, as required.
write.csv(clean, out_path, row.names = FALSE, fileEncoding = "UTF-8")
cat("\nWritten:", out_path, "\n")

# Cleaning log for Section 2.1 of the report (booktabs table, no caption:
# the caption is added in report.tex around \input{tables/cleaning_table.tex}).
tex <- c("\\begin{tabular}{lr}", "\\toprule",
         "Cleaning step & Caf\\'es \\\\", "\\midrule",
         paste0(clean_log$step, " & ", clean_log$cafes, " \\\\"),
         "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(tab_dir, "cleaning_table.tex"))
cat("Written:", file.path(tab_dir, "cleaning_table.tex"), "\n")
