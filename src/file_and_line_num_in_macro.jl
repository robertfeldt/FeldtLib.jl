function each_line_matching(callback, filepath, regexp)
  f = open(filepath, "r")
  i = 0
  for(line in eachline(f))
    i += 1
    if ismatch(regexp, line)
      callback(i, line, line[length(match(regexp, line).match)+1:end])
    end
  end
  close(f)
end

function levenshtein(s::ASCIIString, t::ASCIIString)
  m = length(s)
  n = length(t)
  if n == 0
    return(m)
  elseif m == 0
    return(m)
  end

  d = Array(Int64, m+1, n+1)

  for(i in 1:(m+1))
    d[i, 1] = i
  end
  for(j in 1:(n+1))
    d[1, j] = j
  end
  for(j in 2:n)
    for(i in 2:m)
      if s[i-1] == t[j-1]
        d[i, j] = d[i-1, j-1]
      else
        d[i, j] = minimum([d[i-1,j]+1, d[i, j-1]+1, d[i-1, j-1]+1])
      end
    end
  end
  d[m,n]
end

# Find lines in file
function find_macro_invocations_similar_to(file, macroname, exstr)
  regex = Regex("^\s*@" * macroname)

  # For each line which invokes the macro take what comes after the macroname 
  # and skip spaces etc. Then check which line is closest to exstr.
  lines_with_distance = Any[]
  cb(linenum, line, lineaftermacro) = begin
    lam = replace(lineaftermacro, " ", "")
    push!(lines_with_distance, (line, levenshtein(lam, exstr), linenum))
  end
  each_line_matching(cb, file, regex)

  # Now hope we have the right one... It might not be if there are several same/similar.
  sort(lines_with_distance, by=(t)->t[2])[1]
end

macro mytestmacro(x)
  if !isa(x.args[1], Expr) || x.args[1].head != :line
    file = Base.source_path()
    # Ok, we need to disambiguate..
    line, dist, linenum = find_macro_invocations_similar_to(file, "mytestmacro", string(x)[3:end-1])
    println("line: $(linenum) file: $(file)")
  else
    println("line: $(x.args[1].args[1]) file: $(x.args[1].args[2])")
  end
end

@mytestmacro begin mytestfun() end
@mytestmacro 1==2
@mytestmacro 3<2