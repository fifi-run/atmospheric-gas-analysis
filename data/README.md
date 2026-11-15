# Data access and provenance

## Supplied input

The uploaded coursework archive contains `MaunaLoa_miss.csv` with:
- 183 monthly observations;
- a `Date` column and five numeric variables: `CO`, `CO2`, `Methane`, `NitrousOx`, `CFC11`;
- dates from 2000-01-01 to 2015-03-01;
- 35 originally missing cells: 21 in CO and 14 in CO2.

The coursework report's variable table describes a 2000–2019 date range. This does not match the supplied CSV, which ends in March 2015. The raw file's actual dates are used in this repository's description.

## Redistribution status

The archive did not identify the exact upstream dataset URL, dataset version, original publisher or licence. For that reason, the raw CSV is not committed to this public repository.

Before redistributing the file, identify and verify its source and licence. If the source permits redistribution, add the file here as `MaunaLoa_miss.csv` and record the citation, URL, access date, version and licence in this document.

## Expected path

```text
data/MaunaLoa_miss.csv
```

The date field must be parseable as YYYY-MM-DD and all five gas variables must be numeric. Measurement units were not documented in the supplied data file; do not infer or publish units without confirming the original source metadata.
