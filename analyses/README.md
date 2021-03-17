To run all unit tests, execute
```R
source("estimators.R")
source("coupling-functions.R")
testthat::test_dir("tests/testthat")
```
or
```R
source("<file>.R")
testthat::test_file("tests/testthat/test-<file>.R")
```
for individual tests.

Add `reporter = c("summary", "debug")` as an argument to `test_*` to pause execution at failed `expect_*` instances; for example,
```R
testthat::test_dir("tests/testthat", reporter = c("summary", "debug"))
```
