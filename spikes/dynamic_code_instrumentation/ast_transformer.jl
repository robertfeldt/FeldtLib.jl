# Top type in hierarchy of AST transformers.
abstract ASTTransformer

# ASTTransformers have methods to transform each type of Expr block. Their default
# methods do no transformation, they just leave the tree untouched.
[]

# Default transformer is identity, i.e. doesn't transform at all, but traverses into child trees.
function transform{AT <: ASTTransformer}(at::AT, ast::Expr)
	for i in 1:numsubtrees(ast)
		setchild(ast, i, transform(at, child(ast, i))
	end
	ast
end

function child(ast::Expr, i)
	if i == 1 && in(ast.head, [:call])
		ast

type DictDefinedASTTransformer <: ASTTransformer

end

