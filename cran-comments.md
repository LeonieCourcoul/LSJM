## R CMD check results

0 errors | 0 warnings | 2 notes

## CRAN resubmission fixes

This resubmission addresses all issues raised by CRAN.

- Replaced all occurrences of `T` and `F` with `TRUE` and `FALSE` in the code and documentation.

- Added missing `\value{}` sections in the documentation of exported functions, including `plot.lsjm` and `ranef`, with a detailed description of returned objects and their structure.

- Improved all Rd examples:
  - Removed commented-out code lines.
  - Ensured that all examples are syntactically valid and reproducible.
  - Due to computational cost of the models, all examples are wrapped in `\dontrun{}` as they are not intended to be executed during R CMD check (>5 seconds).

- Replaced direct console output using `print()`/`cat()` in internal and exported functions with `message()` or `warning()`.

No changes were made to the statistical methodology or model definitions; all modifications are related to code quality, documentation, and CRAN compliance.
