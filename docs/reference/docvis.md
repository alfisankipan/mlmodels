# U.S. Medical Expenditure Panel Survey

Cross sectional sample for 2003.

## Usage

``` r
docvis
```

## Format

A data frame with 3677 rows and 22 variables.

- offer:

  Binary: employer offers insurance

- ssiratio:

  Ratio of SSI income to total income

- age:

  Age in years

- educyr:

  Years of education

- physician:

  Number of visits to doctor

- nonphysician:

  Number of visits to health professional (not doctor)

- medicaid:

  Binary: has medicaid public insurance

- private:

  Binary: has private supplementary insurance

- female:

  Binary: female

- phylim:

  Binary: physical limitation

- actlim:

  Binary: activity limitation

- income:

  Income (in thousands)

- totchr:

  Number of chronic conditions

- insured:

  Same as `private`

- age2:

  `age` squared

- linc:

  `log(income)`

- bh:

  Binary: black or hispanic

- docvis:

  Number of doctor visits

- ldocvis:

  `log(docvis)` if `docvis > 0`

- ldocvisa:

  `log(docvis + .01)`

- one:

  Constant term (1)

- docbin:

  Binary: `docvis > 0`

## Source

<https://www.stata-press.com/data/mus2.html> `mus220mepsdocvis.dta`

## References

Cameron, A. C., and P. K. Trivedi. 2022. *Microeconometrics Using Stata,
Second Edition*. Volumes I and II. College Station, TX: Stata Press.
