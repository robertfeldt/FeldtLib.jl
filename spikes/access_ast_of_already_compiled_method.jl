function f(x,y)
	x+y
end

mf1 = methods(f, (Any, Any))[end]
ast1 = Base.uncompressed_ast(mf1.func.code)
@show ast1

h(x) = x+1

function f(x::Int64, y)
	a = h(x+y)
    b = x*y
    # Insert a line in between just to see how lines numbered below...
    a/b
end

mf2 = methods(f, (Int64, Any))[end]
ast2 = Base.uncompressed_ast(mf2.func.code)
@show ast2

body = ast2.args[3]

# Now when we have the methods body Expr we can for example instrument it, here using Match...
#using Match 
#function instrument(e::Expr)
#	@match e begin
#		Expr(:block, , _)     => name
#		Expr(exprtype, _...)     => error("Cannot instrument $e ($exprtype)")
#	end
#end
#@show instrument(body)