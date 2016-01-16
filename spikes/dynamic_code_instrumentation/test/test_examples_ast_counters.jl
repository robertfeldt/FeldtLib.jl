import Base.Test

using DebugPrintln
nodebug()

type FunctionCountingAstWalker <: AstWalker
	count::Int
end
FunctionCountingAstWalker() = FunctionCountingAstWalker(0)

# We count number of visits to function signatures...
function visit_expr_function_sig(fc::FunctionCountingAstWalker, sig, funcex)
	fc.count += 1
	sig
end

astfunc1 = :(
	function f(a)
		a+1
	end
)

fc = FunctionCountingAstWalker()
visit(fc, astfunc1)
@assert fc.count == 1

astfunc2 = quote
	function f1(a)
		a+1
	end
	function f2(b)
		b+2
	end
end

fc2 = FunctionCountingAstWalker()
visit(fc2, astfunc2)
@assert fc2.count == 2


type LineCountingAstWalker <: AstWalker
	count::Int
end
LineCountingAstWalker() = LineCountingAstWalker(0)

function visit_line_node(lc::LineCountingAstWalker, l::LineNumberNode)
	lc.count += 1
	l
end

lc1 = LineCountingAstWalker()
visit(lc1, astfunc1)
@assert lc1.count == 1 # Only actual lines in the body of the funcs are counted

lc2 = LineCountingAstWalker()
visit(lc2, astfunc2)
@assert lc2.count == 4 # Only actual lines in the body of the funcs are counted
