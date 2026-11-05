# Verified software versions

R 4.6.0 on Apple Silicon macOS. Key versions for fitting and reporting:

| Package | Version |
|---|---|
| glmnet | 5.0 |
| mice | 3.19.0 |
| brglm2 | 1.1.0 |
| pROC | 1.19.0.1 |
| Matrix | 1.7-5 |
| pmsampsize | 1.1.3 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| purrr | 1.2.2 |
| tibble | 3.3.1 |
| readxl | 1.5.0 |
| readr | 2.2.0 |
| ggplot2 | 4.0.3 |
| ragg | 1.5.2 |
| stringr | 1.6.0 |
| jsonlite | 2.0.0 |

The Quarto notebook uses the additional dependencies in its setup chunk. The existing project lock lists glmnet 5.1 and omits brglm2 and pmsampsize. The actual refit used the versions above. For a local installation that uses separate R libraries, set `THREEPAD_ADDITIONAL_R_LIB` to the folder that contains missing packages. The isolated scripts add that library after the default R libraries.
