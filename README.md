# Prague Café Market: Specialty Cafés

Semestral project for **Data Analytics**, <batch>, <semester/year>.
Author: Aidan Osbaugh · Team partner: <partner's name>

What do specialty cafés in Prague offer? This repository analyses data on
189 Prague cafés collected in person by the class. It looks at filter coffee
prices, espresso prices, roasteries and café location (research questions
A1–A10 of the Specialty Cafés topic).

## Repository structure

```
analysis/
  raw_data.csv        raw collection sheet (never edited by hand)
  01_cleaning.R       cleaning → dataset.csv
  02_descriptive.R    description, normality and CLT checks → D_table
  03_analysis.R       hypotheses and analyses A1–A10
  dataset.csv         cleaned data (written by 01_cleaning.R)
report/
  report.tex, report.pdf
  tables/             all tables, exported by the R scripts
  figures/            all figures, exported by the R scripts
presentation/
  slides.tex, slides.pdf
```

## How to reproduce

1. Clone the repository and open its root folder in RStudio or VS Code.
2. Run the three scripts in order, from the repository root:
   ```r
   source("analysis/01_cleaning.R")
   source("analysis/02_descriptive.R")
   source("analysis/03_analysis.R")
   ```
3. Compile `report/report.tex` and `presentation/slides.tex` with LaTeX
   (e.g. `latexmk -pdf report.tex` inside `report/`). Numbers quoted in the
   text come from `report/tables/results_macros.tex`, written by `03_analysis.R`.

Every number, table and figure in the report is produced by these scripts.

## Software

- R 4.6.1. `01_cleaning.R` and `02_descriptive.R` use base R only;
  `03_analysis.R` uses `sandwich` 3.1.3, `lmtest` 0.9.40, `car` 3.1.5 and
  `marginaleffects` 1.0.0:
  ```r
  install.packages(c("sandwich", "lmtest", "car", "marginaleffects"))
  ```
- Bootstrap seed: 2026, 5,000 resamples
- Significance level: α = .05