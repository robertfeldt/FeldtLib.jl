using Match

include("debug_println.jl")

# An AstWalker walks over and visits all elements of the AST.
abstract AstWalker

visit(w::AstWalker, s::Symbol) = visit_symbol(w, s)
visit(w::AstWalker, i::Int) = visit_int(w, i)
visit(w::AstWalker, l::LineNumberNode) = visit_line_node(w, l)
visit(w::AstWalker, q::QuoteNode) = visit_quote_node(w, q)
visit(w::AstWalker, e::Expr) = visit_expr(w, e)
visit(w::AstWalker, n) = error("Cannot visit AST node $n of type $(typeof(n))")

visit_expr(w::AstWalker, e::Expr) = dispatch_on_expr_head(w, e)

function dispatch_on_expr_head(w::AstWalker, e::Expr)
	dprintln("Dispatching on $e")
	@match e begin
		Expr(:function,  [signature, body], _)  => visit_expr_function(w, signature, body, e)
		Expr(:block,  stmts, _)      			=> visit_expr_block(w, stmts, e)
		Expr(:call,  [name, args...], _)      	=> visit_expr_call(w, name, args, e)
      	_   									=> error("Cannot visit Expr $e")
  	end
end

# Default behavior is to just visit nodes and return the node itself. 
# Override for more specific behavior.
visit_int(w::AstWalker, i::Int) = i
visit_symbol(w::AstWalker, s::Symbol) = s
visit_quote_node(w::AstWalker, q::QuoteNode) = q
visit_line_node(w::AstWalker, l::LineNumberNode) = l
function visit_expr_function(w::AstWalker, sig, body, funcex)
	visit_expr_function_sig(w, sig, funcex)
	visit_expr_function_body(w, body, funcex)
	funcex
end
visit_expr_function_sig(w::AstWalker, sig, funcex) = sig
function visit_expr_function_body(w::AstWalker, body, funcex)
	visit(w, body)
	body
end
function visit_expr_block(w::AstWalker, stmts, blockex)
	dprintln("Visiting block:\n$e")
	map(s -> visit(w, s), stmts)
	blockex
end
function visit_expr_call(w::AstWalker, name, args, callex)
	map(s -> visit(w, s), args)
	callex
end