# Data

The PIAAC microdata are **not included** in this repository. They are freely
available from the OECD, but they are not mine to redistribute.

## How to get them

1. Download the **PIAAC Cycle 2 Public Use Files (CSV)** and the label
   dictionary from the OECD website:
   <https://www.oecd.org/en/about/programmes/piaac/piaac-data.html>
2. Use the OECD's R import script (distributed with the files) to merge the
   country files into a single data frame.
3. Save that data frame as `data/raw/PIAAC_Cycle2.rds`:

   ```r
   saveRDS(PIAAC_Cycle2, "data/raw/PIAAC_Cycle2.rds")
   ```

`R/01_data_cleaning.R` reads this file and writes the cleaned datasets
(`df_clean.rds`, `df_model.rds`, `train.rds`, `test.rds`) to `data/processed/`.

## Other sources

PPP conversion rates (2022, national currency per USD) come from the OECD
Data Explorer (National Accounts, table DSD_NAMAIN10 / DF_TABLE4) and are
hard-coded in `R/01_data_cleaning.R`, so no extra download is needed.
