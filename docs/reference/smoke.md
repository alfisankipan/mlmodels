# 1979 National Health Interview Survey

Sample of males in 1979 and early 1980 from the smoking supplement.

## Usage

``` r
smoke
```

## Format

A data frame with 807 observations and 10 variables.

- educ:

  Years of education

- cigpric:

  State's cigarette price (cents per pack)

- white:

  Binary: white

- age:

  Age in years

- income:

  Annual income in dollars

- cigs:

  Number of cigarettes smoked per day

- restaurn:

  Binary: restaurant restrictions on smoking

- lincome:

  `log(income)`

- agesq:

  `age^2`

- lcigpric:

  `log(cigprice)`

## Source

<https://cran.r-project.org/package=wooldridge>

## References

Wooldridge, J. M. (2020). *Introductory Econometrics: A Modern
Approach*. 7th Edition. Boston, MA: Cengage Learning.

Mullahy, J. (1997). "Instrumental-Variable Estimation of Count Data
Models: Applications to Models of Cigarette Smoking and Infant Health."
*Review of Economics and Statistics* 79, 586-593.
