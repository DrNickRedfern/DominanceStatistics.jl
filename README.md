# DominanceStatistics.jl

The module `DominanceStatistics.jl` includes a range of functions for calculating statistics from the difference or dominance matrix for two samples, including the Hodges-Lehmann median difference, Cliff's $d$ statistic of stochastic equality, and Vargha and Delaney's $A$ statistic.

First, we load the module:

```julia
include("./DominanceStatistics.jl")
using .DominanceStatistics
```

Next, let's create some data to play with:

```julia
sample1 = [3.1, 3.3, 3.3, 3.4, 3.8, 3.8, 3.8, 3.8, 4.1, 4.3]
sample2 = [3.0, 3.8, 4.0, 4.2, 4.7, 5.0, 5.0, 5.5, 7.2, 8.8]
```

## Dominance and difference matrices
The dominance matrix can be calculated directly if desired using `dominanceMatrix`, returning a matrix with values

- `1` where values in `sample1` are greater than values in `sample2`
- `0` where values in `sample1` are equal to values in `sample2`
- -1 where values in `sample1` are less than values in `sample2`

```julia
m = dominanceMatrix(sample1, sample2)
```
```
10×10 Matrix{Float64}:
 1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0
 1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0
 1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0
 1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0
 1.0   0.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0
 1.0   0.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0
 1.0   0.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0
 1.0   0.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0
 1.0   1.0   1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0
 1.0   1.0   1.0   1.0  -1.0  -1.0  -1.0  -1.0  -1.0  -1.0
```

The difference matrix can also be accessed directly using `differenceMatrix`, which returns a matrix containing the calculated difference between the `i`-th element of `sample1` and the `j`-th element of `sample2`.

```julia
M = differenceMatrix(sample1, sample2)
```
```
10×10 Matrix{Float64}:
 0.1  -0.7  -0.9  -1.1  -1.6  -1.9  -1.9  -2.4  -4.1  -5.7
 0.3  -0.5  -0.7  -0.9  -1.4  -1.7  -1.7  -2.2  -3.9  -5.5
 0.3  -0.5  -0.7  -0.9  -1.4  -1.7  -1.7  -2.2  -3.9  -5.5
 0.4  -0.4  -0.6  -0.8  -1.3  -1.6  -1.6  -2.1  -3.8  -5.4
 0.8   0.0  -0.2  -0.4  -0.9  -1.2  -1.2  -1.7  -3.4  -5.0
 0.8   0.0  -0.2  -0.4  -0.9  -1.2  -1.2  -1.7  -3.4  -5.0
 0.8   0.0  -0.2  -0.4  -0.9  -1.2  -1.2  -1.7  -3.4  -5.0
 0.8   0.0  -0.2  -0.4  -0.9  -1.2  -1.2  -1.7  -3.4  -5.0
 1.1   0.3   0.1  -0.1  -0.6  -0.9  -0.9  -1.4  -3.1  -4.7
 1.3   0.5   0.3   0.1  -0.4  -0.7  -0.7  -1.2  -2.9  -4.5
```

## Hodges-Lehmann median difference
To calculate the [Hodges-Lehmann median difference](https://en.wikipedia.org/wiki/Hodges–Lehmann_estimator) for two samples, we calculate the difference matrix of our two data sets and find the median value: $HL\Delta = median(X_{i}-Y_{j})$.

```julia
HLD(x, y; xname="A", yname="B", nboot=1000, conf_level=0.95, print_results=true)
```

A bootstrapped confidence interval using [Bootstrap.jl](https://github.com/juliangehring/Bootstrap.jl) is also returned, with options `Percentile`, `Basic`, `BCa`, and `Normal`. The number of bootstrap samples and the confidence level can also be set.

By setting `print_results=false` the summary of results is suppressed.

```julia
res = HLD(sample1, sample2; xname="Control", yname="Treatment", type="Percentile", nboot=200, conf_level=0.95, print_results=true);
```
```
Hodges-Lehmann median difference
Differences are calculated as Control - Treatment.
Positive differences indicate that Control > Treatment.
HLΔ = -1.15 (-1.2, -0.85)
Bias = 0.094
SE = 0.148
Type = Percentile
Confidence level = 95.0%
```
The fieldnames of the struct storing the results of the function can be accessed using:
```julia
fieldnames(typeof(res))
```
```
(:hld, :CI, :type, :bias, :se, :level)
```
Individual results can be accessed from the struct in the usual way:
```julia
res.hld
```
```
-1.15
```
The struct storing the function's outputs can be dumped in one go:
```julia
dump(res)
```
```
Main.DominanceStatistics.hldresults{Float64}
  hld: Float64 -1.15
  CI: Tuple{Float64, Float64}
    1: Float64 -1.2
    2: Float64 -0.85
  type: String "Percentile"
  bias: Float64 0.094
  se: Float64 0.148
  level: Float64 0.95
```

## Dominance effect sizes
The `dominanceEffectSizes` function returns

- Cliff's $d$ statistic of stochastic equality, $d = P(X > Y) - P(X < Y)$, along with a confidence interval, $Z$-statistic, and $p$-value.
- Vargha and Delaney's $A$ statistic, where $A = P(X > Y) + 0.5(X = Y)$, with confidence interval.

```julia
dominanceEffectSizes(x, y; xname="A", yname="B", conf_level=0.95, print_results="true")
```

This function is based on the methods described in Wilcox, R.R. (2017) *Introduction to Robust Hypothesis Testing*, 4th Edn. London: Academic Press: 192-193. 

Applying `dominanceEffectSizes` to our data:

```julia
res = dominanceEffectSizes(sample1, sample2; xname="Control", yname="Treatment", print_results=true);
```
```
Differences are calculated as Control - Treatment.
Positive differences indicate that Control > Treatment.
Cliff's d: -0.66 (-0.903, -0.097)
Z = -3.046, p = 0.0023
Vargha & Delaney's A: 0.17 (0.049, 0.451)
Confidence level = 95.0%
```

Individual elements of the struct storing the results of the function can be accessed as above.

The confidence interval for $d$ can be calculated directly from a dominance matrix and a given value of $d$. The $Z$-statistic and $p$-value are also returned.

```julia
dominanceCI(m, d=-0.66; conf_level = 0.95)
```
```
((-0.9029463650549244, -0.09737031969644462), -3.046, 0.0023)
```