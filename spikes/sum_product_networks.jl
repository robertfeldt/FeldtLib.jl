# Sum-Product Networks as described by Poon&Domingos, Gens&Domingos etc.

abstract SPNode{T} # SPN node with names of type T

immutable Sum{T} <: SPNode{T}
  children::Vector{SPNode{T}}
  weights::Vector{Float64}
end

immutable Prod{T} <: SPNode{T}
  children::Vector{SPNode{T}}
end

# Variables have names of a given type.
immutable Var{T} <: SPNode{T}
  name::T
end

immutable NVar{T} <: SPNode{T}
  name::T
  negname::T
end
negname(name::Symbol) = symbol("__N__" * string(name))
negname{I<:Integer}(name::I) = -name
NVar(name) = NVar{typeof(name)}(name, negname(name))

# A SPN with variable names of type T.
type SumProductNetwork{T}
  root::SPNode{T}
end
SPN{T}(r::SPNode{T}) = SumProductNetwork{T}(r)

# Direct calculation of probabilities
calc{T, N <: Number}(spn::SumProductNetwork{T}, vals::Dict{T,N}) = calc(spn.root, vals)
function calc{T, N <: Number}(s::Sum{T}, vals::Dict{T,N})
  sum(calcchildren(s.children, vals) .* s.weights)
end
calcchildren{T, N <: Number}(children::Vector{SPNode{T}}, vals::Dict{T,N}) = map(c->calc(c, vals), children)
calc{T, N <: Number}(s::Prod{T}, vals::Dict{T,N}) = prod(calcchildren(s.children, vals))
calc{T, N <: Number}(v::Var{T}, vals::Dict{T,N}) = haskey(vals, v.name) ? vals[v.name] : 1
function calc{T, N <: Number}(v::NVar{T}, vals::Dict{T,N})
  haskey(vals, v.negname) ? vals[v.negname] : (haskey(vals, v.name) ? (1-vals[v.name]) : 1)
end

# But we also want to be able to calculate probabilities 
# when values of some variables is not known. This is done by setting both
# the vars and negvars of the unknown variables to 1. Eas

# Simple examples
s1 = SPN(Var(1))
calc(s1, Dict{Int64,Int64}(1 => 1))

s2 = SPN(NVar(1))
calc(s2, Dict{Int64,Int64}(1 => 1))
calc(s2, Dict{Int64,Int64}(1 => 0))
calc(s2, Dict{Int64,Int64}(-1 => 1))
calc(s2, Dict{Int64,Int64}(-1 => 0))
calc(s2, Dict{Int64,Int64}())

s3 = SPN(Var(:a))
calc(s3, Dict{Symbol,Int64}(:a => 1))

# Fig 1 from cs.ubc.ca paper:
at = Sum(SPNode{Symbol}[Var(:a), NVar(:a)], Float64[0.5, 0.5])
bt = Sum(SPNode{Symbol}[Var(:b), NVar(:b)], Float64[0.5, 0.5])
ct = Sum(SPNode{Symbol}[Var(:c), NVar(:c)], Float64[0.5, 0.5])
dt = Sum(SPNode{Symbol}[Var(:d), NVar(:d)], Float64[0.5, 0.5])
left = Prod(SPNode{Symbol}[at, bt])
right = Prod(SPNode{Symbol}[ct, dt])
r = Sum(SPNode{Symbol}[left, right], Float64[0.5, 0.5])
@time calc(r, Dict{Symbol,Int64}(:a => 1, :b => 0, :c => 1, :d => 1))
calc(r, Dict{Symbol,Int64}(:a => 1, :b => 0, :c => 1))
