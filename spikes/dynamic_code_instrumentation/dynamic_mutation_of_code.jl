# A general way to do mutation testing is to rewrite the AST of a function
# so that we can then dynamically select what, if anything, to mutate.
# This is a manual test of this.

# Target function we want to mutate
function tf(x::Integer)
    a = 2*x
    res = 0
    for i in 1:a
        res += x * i
    end
    if a > 9
        res = -res
    end
    return res
end

# Abstract class for any dynamic code instrumentation scheme and the
abstract DynamicInstrumenter

# The default instrumenter has inlined no-ops as its instrumentation implementation =>
# as little overhead as possible
@inline start_execution(d::DynamicInstrumenter, funcname::Symbol, types, maxlines::Int) = nothing
@inline exec_line(d::DynamicInstrumenter, linenum::Int) = true
@inline binop_mult(d::DynamicInstrumenter, lhs, rhs) = lhs * rhs
@inline binop_plus(d::DynamicInstrumenter, lhs, rhs) = lhs + rhs
@inline binop_geq(d::DynamicInstrumenter, lhs, rhs) = lhs >= rhs
@inline binop_gt(d::DynamicInstrumenter, lhs, rhs) = lhs > rhs
@inline unop_neg(d::DynamicInstrumenter, val) = -val
@inline const_val(d::DynamicInstrumenter, val) = val
@inline const_var(d::DynamicInstrumenter, varname::Symbol, val) = val
@inline forloop_iterator(d::DynamicInstrumenter, it) = it
@inline stop_execution(d::DynamicInstrumenter, funcname::Symbol, types, retval) = retval


# The target function tf is rewritten dynamically to (this is the default instrumentation,
# to speed things up individual instrumenters might need to make leaner versions, this includes a lot 
# to allow dynamic instrumentation/mutation etc):
function tf(d::DynamicInstrumenter, x::Integer)
    start_execution(d, :tf, (Int), 7)
    if exec_line(d, 1)
        a = binop_mult(d, 2, x)
    end
    if exec_line(d, 2)
        res = const_val(d, 0)
    end
    if exec_line(d, 3)
        for i in forloop_iterator(d, 1:a)
            if exec_line(d, 4)
                res = binop_plus(d, res, binop_mult(d, x, i))
            end
        end
    end
    if exec_line(d, 5)
        if binop_gt(d, a, 9)
            if exec_line(d, 6)
                res = unop_neg(d, res)
            end
        end
    end    
    if exec_line(d, 7)
        retval = const_var(d, :res, res)
        return stop_execution(d, :tf, (Int), retval)
    end
    stop_execution(d, :tf, (Int), nothing)
end


# Instrumenter 1.
# Toy instrumenter that just measures the execution time. Not very well though since
# it is hard to separate the instrumented code from the original code...
type ExecTimingInstrumenter <: DynamicInstrumenter
    start::Float64
    elapsed::Float64
    ExecTimingInstrumenter() = new(0.0, 0.0)
end
@inline start_execution(i::ExecTimingInstrumenter, funcname::Symbol, types, maxlines::Int) = (i.start = time())
@inline function stop_execution(i::ExecTimingInstrumenter, funcname::Symbol, types, retval)
    i.elapsed = time() - i.start
    retval
end

arg1 = 5
@time r1 = tf(arg1)
eti = ExecTimingInstrumenter()
@time r2 = tf(eti, arg1)
@show (r1 == r2, eti) # Not getting the performance I expected here yet, but solvable in principle with better type info?


# Instrumenter 2.
# A more interesting instrumenter is one that counts which statements has been executed.
type StmtCoverageInstrumenter <: DynamicInstrumenter
    totallines::Int
    executions::Vector{Set{Int}}
    StmtCoverageInstrumenter() = new(-1, Set{Int}[])
end
function start_execution(s::StmtCoverageInstrumenter, funcname::Symbol, types, maxlines::Int)
    s.totallines = maxlines
    push!(s.executions, Set{Int}())
end
function exec_line(s::StmtCoverageInstrumenter, linenum::Int)
    push!(s.executions[end], linenum)
    true
end
stmt_coverages(s::StmtCoverageInstrumenter) = Float64[stmt_coverage(s, executedlines) for executedlines in s.executions]
function stmt_coverage(s::StmtCoverageInstrumenter, execlines::Set{Int})
    length(unique(execlines)) / s.totallines
end
all_covered_lines(s::StmtCoverageInstrumenter) = foldl(union, Set{Int}(), s.executions)
function stmt_coverage(s::StmtCoverageInstrumenter)
    stmt_coverage(s, all_covered_lines(s))
end

scov = StmtCoverageInstrumenter()
@time r3 = tf(scov, arg1)
@show (r1 == r3)
stmt_coverage(scov) # Should be 1.0

# But when argument is below 5 we don't cover all lines
scov2 = StmtCoverageInstrumenter()
arg2 = 4
@time r4 = tf(scov2, arg2)
@show (tf(arg2) == r4)
stmt_coverage(scov2) # should be 6/7 since only the unary negation is not executed

# Instrumenter 3.
# Mutates binary operator + to * and vice versa.
type BinopMutatingInstrumenter <: DynamicInstrumenter
end
@inline binop_mult(b::BinopMutatingInstrumenter, lhs, rhs) = lhs + rhs
@inline binop_plus(d::BinopMutatingInstrumenter, lhs, rhs) = lhs * rhs

mut = BinopMutatingInstrumenter()
@time r5 = tf(mut, arg1)
@show (r1 != r5)