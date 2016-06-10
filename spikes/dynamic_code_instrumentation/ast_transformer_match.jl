

# An AST walker that knows about the context of the currently visited node in the AST.

type ContextAstWalker <: AstWalker
end




# The default behavior is to only copy each AST node.
transform(at::AstTransformer, n) = deepcopy(n)


# For the expressions we must traverse sub-trees.
function transform(at::AstTransformer, e::Expr)
  # println("Transforming expression: head = ", e.head, "\n ", e)

  @match e begin

      Expr(:function,  [signature, body], _)      => begin
      	# println("Function body = ", body)

      	Expr(:function, signature, ast_transform(at, body))
      end

      Expr(:block,  stmts, _)      => begin
      	# We create a new array to hold the statements since we might need to insert new ones while
      	# transforming the old ones.
      	transformed_stmts = Any[]

      	# Populate it by transforming all sub-stmts in this block
      	for i in 1:length(stmts)
      		res = transform(at, stmts[i])
      		if isa(res, Array)
      			# Add all statements if an Array was returned
      			push!(transformed_stmts, res...)
      		else
      			push!(transformed_stmts, res)
      		end
      	end

      	# And create and return the transformed block
      	Expr(:block, transformed_stmts...)
      end

      # Default case is just to copy
      _   => deepcopy(e)
  end

end


dump(a)

at = AstTransformer()
a2 = ast_transform(at, a)
dump(a2)