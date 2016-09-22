testfile = ARGS[1]

expandtabs(str) = replace(str, r"\t", "  ")
totab4(str) = " " ^ (div(length(expandtabs(str)), 2) * 4)

function transform_line(line::AbstractString)
    m = match(r"^(\s*)facts\(\"(.*)\"\)\s+do(\s*)\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@testset \"$(m.captures[2])\" begin\n"
    end

    m = match(r"^(\s*)context\(\"(.*)\"\)\s+do\s*\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@testset \"$(m.captures[2])\" begin\n"
    end

    m = match(r"^(\s*)@fact_throws\s+(.*)\s+([A-Za-z]+)\s*\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@test_throws $(m.captures[3]) $(m.captures[2])\n"
    end

    m = match(r"^(\s*)@fact\s+(.*)\s+-->\s+is([a-z]+)\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@test is$(m.captures[3])($(m.captures[2]))\n"
    end

    m = match(r"^(\s*)@fact\s+(.*)\s+-->\s+not\(is([a-z]+)\)\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@test !is$(m.captures[3])($(m.captures[2]))\n"
    end

    m = match(r"^(\s*)@fact\s+(.*)\s+-->\s+less_than_or_equal\((.*)\)\s*\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@test $(m.captures[2]) <= $(m.captures[3])\n"
    end

    m = match(r"^(\s*)@fact\s+(.*)\s+-->\s+greater_than_or_equal\((.*)\)\s*\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@test $(m.captures[2]) >= $(m.captures[3])\n"
    end

    m = match(r"^(\s*)@fact\s+(.*)\s+-->\s+less_than\((.*)\)\s*\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@test $(m.captures[2]) < $(m.captures[3])\n"
    end

    m = match(r"^(\s*)@fact\s+(.*)\s+-->\s+greater_than\((.*)\)\s*\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@test $(m.captures[2]) > $(m.captures[3])\n"
    end

    m = match(r"^(\s*)@fact\s+(.*)\s+-->\s+true\s*\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@test $(m.captures[2])\n"
    end

    m = match(r"^(\s*)@fact\s+(.*)\s+-->\s+(.*)\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@test $(m.captures[2]) == $(m.captures[3])\n"
    end

    m = match(r"^(\s*)(.*)\s*\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))$(m.captures[2])\n"
    end
    
    line
end

newlines = open(testfile, "r") do fh
    map(transform_line, readlines(fh))
end

# Delete so only max one empty line at end
while length(newlines) >= 2 && strip(newlines[end]) == "" && strip(newlines[end-1]) == ""
    deleteat!(newlines, length(newlines))
end

altname = tempname()
open(altname, "w") do fh
    print(fh, join(newlines, ""))
end

cmd = `diff $(testfile) $(altname)`
try
    run(cmd)
catch _e
end

# Remove file unless we should write it
if length(ARGS) < 2
    rm(altname)
elseif length(ARGS) >= 2 && in(lowercase(ARGS[2]), ["o", "-o", "--overwrite", "overwrite"])
    mv(altname, testfile; remove_destination = true)
else length(ARGS) >= 2 && in(lowercase(ARGS[2]), ["s", "-s", "--save", "write"])
    println("Altered file saved to: ", altname)
end
