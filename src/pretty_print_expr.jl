IndentStringCache = Dict{Int, ASCIIString}()

indentstr(indent) = get!(IndentStringCache, indent, join([" " for i in 1:indent]))

iprint(s, indent) = print(indentstr(indent), s)

function pp(d, indent = 0)
	print(indentstr(indent))
	show(d)
end

pp(s::String, indent = 0) = print(indentstr(indent), s)

# Pretty-printing of Julia Expr's
function pp(ex::Expr, indent = 0)
	# Skip line nodes unless told to print them
	if ex.head == :line && includeLineNodes == false
		return()
	end

	pp("head: $(ex.head)\n", indent)
	pp("args:\n", indent+1)
	for arg in ex.args
		pp(arg, indent+2)
		println()
	end
end
