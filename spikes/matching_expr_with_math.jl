using Match

function rw(ex)
  @match ex begin
    Expr(:comparison, [l,op,r], _) => "Comparison, lhs = $(l)"
    _ => "No match"
  end
end

rw( :(a == 1) )