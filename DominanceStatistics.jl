"""
This Julia module implements a number of functions for dominance statistics.

# Functions:
    dominanceMatrix - Compute the dominance matrix for two numeric vectors.
    differenceMatrix - Compute the difference matrix for two numeric vectors.
    HLD - Compute the Hodges-Lehmann median difference, with bootstrapped confidence intervals.
    dominanceEffectSizes - Compute a range of effect sizes.
    MannWhitneyOdds - Compute the results of the Mann-Whitney test and the MW odds ratio.

# Author:
    Nick Redfern
    nickredfernres@outlook.com
    http://computationalfilmanalysis.wordpress.com
    https://github.com/DrNickRedfern
    14/09/22
"""

module DominanceStatistics

using DataFrames, Distributions, Statistics, Bootstrap, StatsBase

export dominanceMatrix, differenceMatrix, HLD, dominanceEffectSizes, MannWhitneyOdds, dominanceCI

"""
Custom type for accesing results from HLD:
    hld: the median of the difference matrix.
    CI: bootstrapped confidence interval
    type: the tpe of bootstrapped confidence interval
    bias: bias of the bootstrap estimate
    se: standard error of the bootstrap estimate
    level: confidence level of the confidence interval
"""
struct hldresults{T<:Real}
    hld::T 
    CI::Tuple{T,T} 
    type::String 
    bias::T 
    se::T 
    level::T 
end

"""
Type for accessing results from dominanceEffectSizes:
    d: Cliff's d statistic of stochastic equality
    d_ci: confidence interval for d
    Z: Z statistic
    p_value: p-value
    A: Vargha-Delaney's A statistic
    A:_ci: confidence interval for A
    level: confidence level for the confidence intervals
"""
struct desresults{T<:Real}
    d::T 
    d_ci::Tuple{T,T} 
    Z::T 
    p_value::T 
    A::T 
    A_ci::Tuple{T,T} 
    level::T 
end

"""
   M = dominanceMatrix(x, y)

Compute the dominance matrix for two numeric vectors. Returns a named matrix M where:
    - M[i,j] = 1 if x[i] > y[j] 
    - M[i,j] = 0 if x[i] = y[j]
    - M[i,j] = -1 if x[i] < y[j]

# Arguments:
- `x::Array{Real}`: a numeric array.
- `y::Array{Real}`: a numeric array.

# Example:
```
# Load the dominanceStatistics module
include("./dominanceStatistics.jl")
using .dominanceStatistics

# Generate some data
sample1 = [3.1, 3.3, 3.3, 3.4, 3.8, 3.8, 3.8, 3.8, 4.1, 4.3]
sample2 = [3.0, 3.8, 4.0, 4.2, 4.7, 5.0, 5.0, 5.5, 7.2, 8.8]

# Calculate the matrix
dominanceMatrix(sample1, sample2)
```
"""
function dominanceMatrix(x::Vector{<:Real}, y::Vector{<:Real})

    nx = length(x)
    ny = length(y)

    dominance_matrix = Array(zeros(nx, ny))
    for i in 1:nx
        for j in 1:ny
            dominance_matrix[i, j] = sign(x[i] - y[j])
        end
    end

    return dominance_matrix

end

"""
   M = differenceMatrix(x, y)

Compute the difference matrix for two numeric vectors. Returns a named matrix M where:
    - M[i,j] = x[i] - y[j]

# Arguments:
- `x::Array{Real}`: a numeric array.
- `y::Array{Real}`: a numeric array.

# Example:
```
# Load the dominanceStatistics module
include("./dominanceStatistics.jl")
using .dominanceStatistics

# Generate some data
sample1 = [3.1, 3.3, 3.3, 3.4, 3.8, 3.8, 3.8, 3.8, 4.1, 4.3]
sample2 = [3.0, 3.8, 4.0, 4.2, 4.7, 5.0, 5.0, 5.5, 7.2, 8.8]

# Calculate the matrix
differenceMatrix(sample1, sample2)
```
"""

function differenceMatrix(x::Vector{<:Real}, y::Vector{<:Real})

    nx = length(x)
    ny = length(y)

    difference_matrix = Array(zeros(nx, ny))
    for i in 1:nx
        for j in 1:ny
            difference_matrix[i, j] = x[i] - y[j]
        end
    end

    return difference_matrix

end

"""
   Δ = HLD(x, y; xname = "A", yname = "B", nboot = 1000, conf_level = 0.95, print_results=true)

Compute the Hodges-Lehmann median diffrence from the difference matrix of two independent samples with bootstrapped confidence intervals.

# Arguments:
- `x::Array{Real}`: a numeric array.
- `y::Array{Real}`: a numeric array.
- `xname::String`: the name of the vector x.
- `yname::String`: the name of the vector y.
- `type::String`: the type of bootstrap interval. One of Percentile, Basic, BCa, Normal.
- `nboot::int64`: the number of samples for bootstrap.
- `conf_level::Real`: the confidence level for the bootstrapped intervals.
- `print_results::Bool`: should a summary of the results be printed?

# Example:
```
# Load the dominanceStatistics module
include("./dominanceStatistics.jl")
using .dominanceStatistics

# Generate some data
sample1 = [3.1, 3.3, 3.3, 3.4, 3.8, 3.8, 3.8, 3.8, 4.1, 4.3]
sample2 = [3.0, 3.8, 4.0, 4.2, 4.7, 5.0, 5.0, 5.5, 7.2, 8.8]

# Calculate the Hodges-Lehmann median difference
HLD(sample1, sample2; xname = "Treatment", yname = "Control", type = "Percentile", nboot = 1000, conf_level = 0.95, print_results = true)
```
"""
function HLD(x::Vector{<:Real}, y::Vector{<:Real}; xname="A"::String, yname="B"::String, type="Percentile"::String, nboot=1000::Int64, conf_level=0.95::Float64, print_results=true::Bool)

    v = sort(vec(differenceMatrix(x, y)))

    HLΔ = round(median(v); digits=2)
    bootΔ = bootstrap(median, v, BasicSampling(nboot))
    bs = round(Bootstrap.bias(bootΔ)[1]; digits=3)
    se = round(Bootstrap.stderror(bootΔ)[1]; digits=3)

    if type == "Percentile"
        CI = round.(confint(bootΔ, PercentileConfInt(conf_level))[1][2:3], digits=4)
    end

    if type == "Basic"
        CI = round.(confint(bootΔ, BasicConfInt(conf_level))[1][2:3]; digits=4)
    end

    if type == "BCa"
        CI = round.(confint(bootΔ, BCaConfInt(conf_level))[1][2:3]; digits=4)
    end

    if type == "Normal"
        CI = round.(confint(bootΔ, NormalConfInt(conf_level))[1][2:3]; digits=4)
    end

    if print_results == true
        println("Hodges-Lehmann median difference")
        println("Differences are calculated as $xname - $yname.")
        println("Positive differences indicate that $xname > $yname.")
        println("HLΔ = ", HLΔ, " ", CI)
        println("Bias = $bs")
        println("SE = $se")
        println("Type = $type")
        println("Confidence level = ", 100 * conf_level, "%")
    end

    hldresults(HLΔ, CI, type, bs, se, conf_level)

end

"""
   ES = dominanceEffectSizes(x, y; xname = "A", yname = "B", conf_level = 0.95, print_results = "true")

Compute the Cliff's d statistic of stochastic equality (d = P(A > B) - P(A < B)) with confidence interval and Vargha and Delaney's A statistic of stochastic superiority (A = p(A = B) + 0.5P(A = B)) with confidence interval.

# Arguments:
- `x::Array{Real}`: a numeric array.
- `y::Array{Real}`: a numeric array.
- `xname::String`: the name of the vector x.
- `yname::String`: the name of the vector y.
- `conf_level::Real`: the confidence level for the confidence intervals.
- `print_results::Bool`: should a summary of the results be printed


# Example:
```
# Load the dominanceStatistics module
include("./dominanceStatistics.jl")
using .dominanceStatistics

# Generate some data
sample1 = [3.1, 3.3, 3.3, 3.4, 3.8, 3.8, 3.8, 3.8, 4.1, 4.3]
sample2 = [3.0, 3.8, 4.0, 4.2, 4.7, 5.0, 5.0, 5.5, 7.2, 8.8]

# Calculate the effect sizes
dominanceEffectSizes(sample1, sample2; xname = "Treatment", yname = "Control", conf_level = 0.95, print_results = true)
```
"""

function dominanceEffectSizes(x::Vector{<:Real}, y::Vector{<:Real}; xname="A"::String, yname="B"::String, conf_level=0.95::Float64, print_results=true::Bool)

    dm = dominanceMatrix(x, y)

    # Cliff's d
    d = mean(dm)

    # Vargha & Delaney's A
    A = round(0.5 * (d + 1); digits=3)

    res = dominanceCI(dm, d=d; conf_level=conf_level)
    dci = (round(res[1][1]; digits=3), round(res[1][2]; digits=3))
    Aci = (round(0.5 * (res[1][1] + 1); digits=3), round(0.5 * (res[1][2] + 1); digits=3))

    if print_results == true
        println("Differences are calculated as $xname - $yname.")
        println("Positive differences indicate that $xname > $yname.")
        println("Cliff's d: ", d, " ", dci)
        println("Z = ", res[2], ", p = ", res[3])
        println("Vargha & Delaney's A: ", A, " ", Aci)
        println("Confidence level = ", 100 * conf_level, "%")
    end

    desresults(d, dci, res[2], res[3], A, Aci, conf_level)

end

"""
   ci = dominanceCI(dm, d=-0.66; conf_level = 0.95)

Compute the confidence interval for Cliff's d statistic of stochastic equality and Vargha and Delaney's A statistic of stochastic superiority. This method is based on Wilcox, R.R. (2017) Introduction to Robust Hypothesis Testing, 4th Edn. London: Academic Press: 192-193. Also returns the Z statistic and a p-value.

# Arguments:
- `dm::Matrix{Real}`: a dominance matrix.
- `d::Real`: Cliff's d statistic of stochastic equality.
- `conf_level::Real`: the confidence level for the intervals.

# Example:
```
# Load the dominanceStatistics module
include("./dominanceStatistics.jl")
using .dominanceStatistics

# Generate some data
sample1 = [3.1, 3.3, 3.3, 3.4, 3.8, 3.8, 3.8, 3.8, 4.1, 4.3]
sample2 = [3.0, 3.8, 4.0, 4.2, 4.7, 5.0, 5.0, 5.5, 7.2, 8.8]

# Calculate the confidence interval
dominanceCI(dm, d=-0.66; conf_level = 0.95)
```
"""
function dominanceCI(x::Matrix{<:Real}; d=0.0<:Real, conf_level=0.95::Float64)

    nrow = size(x, 1)
    ncol = size(x, 2)

    z = abs(quantile(Normal(), (1 - conf_level) / 2))

    colmeans = mean(x, dims=1)
    rowmeans = mean(x, dims=2)
    var_rowmeans = sum([(rowmeans[i] - d)^2 for i in 1:nrow])
    var_colmeans = sum([(colmeans[i] - d)^2 for i in 1:ncol])

    var_dominancematrix = Array(zeros(nrow, ncol))

    for i in 1:nrow
        for j in 1:ncol
            var_dominancematrix[i, j] = (x[i, j] - d)^2
        end
    end

    S = (((ncol - 1) * (var_rowmeans / (nrow - 1))) + ((nrow - 1) * (var_colmeans / (ncol - 1))) + (sum(var_dominancematrix) / ((nrow - 1) * (ncol - 1)))) / (nrow * ncol)
    Z = round(d / (sqrt(S)); digits=3)
    p = round(2 * cdf(Normal(), Z); digits=4)

    conf_int = ((d - d^3) - z * sqrt(S) * sqrt((1 - d^2)^2 + (z^2 * S))) / ((1 - d^2) + (z^2 * S)),
    ((d - d^3) + z * sqrt(S) * sqrt((1 - d^2)^2 + (z^2 * S))) / ((1 - d^2) + (z^2 * S))

    return conf_int, Z, p

end

end