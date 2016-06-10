# The basic idea is that we want to be able to instrument one and the same method multiple times 
# without them interferring with each other. Tim Holy on julia-dev proposed:
#
#   @instrument Tool f(x, y)
#
# which would get the AST for f and then create a new function
#
#   f{T <: Tool}(t::T, x, y)
#
# which when called
#
#   t = Tool(...) # construct a specific tool to collect the instrumented info
#   f(t, x, y)    # will now run the instrumented version
#
#  after which we can inspect t to get info about the execution.
#
# Below we test how this could be used assuming it is implemented in a nice way.

# Function we want to instrument. We have two different instrumentation goals: (a) to
# count the number of lines of code that is executed, (b) to count the number of methods
# that are called.
function origf(x::Int, y::Int)
	if x < 0
		return y
	elseif x > y
		res = y
		for i in 1:(x-y)
			res += y
		end
		return res
	end
end

using CodeInstrumenter

type CountExecutedLines <: CodeInstrumenter
end

@instrument CountExecutedLines origf(x::Int, y::Int)

