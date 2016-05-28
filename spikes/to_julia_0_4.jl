testfile = ARGS[1]

totab4(str) = " " ^ (div(length(str), 2) * 4)

function transform_line(line::AbstractString)
    m = match(r"^(.+)::String\s*,(.*)\n$", line)
    if m != nothing
        return transform_line("$(m.captures[1])::AbstractString,$(m.captures[2])\n")
    end

    m = match(r"^(.+)::String\s*\)(.*)\n$", line)
    if m != nothing
        return transform_line("$(m.captures[1])::AbstractString)$(m.captures[2])\n")
    end

    m = match(r"^(.+)::String\n$", line)
    if m != nothing
        return transform_line("$(m.captures[1])::AbstractString\n")
    end

    m = match(r"^(.+)<:\s+String\s*\}(.*)\n$", line)
    if m != nothing
        return transform_line("$(m.captures[1])<: AbstractString}$(m.captures[2])\n")
    end

    m = match(r"^(.+)<:\s+String\s*,(.*)\n$", line)
    if m != nothing
        return transform_line("$(m.captures[1])<: AbstractString,$(m.captures[2])\n")
    end

    m = match(r"^(.+)<:\s+String\s+$", line)
    if m != nothing
        return transform_line("$(m.captures[1])<: AbstractString\n")
    end

    m = match(r"^(.*)Dict\{\s*String\s*,(.*)\n", line)
    if m != nothing
        return transform_line("$(m.captures[1])Dict{AbstractString,$(m.captures[2])\n")
    end

    m = match(r"^(\s*)(.+)\{\s*String\s*\}(.*)\n", line)
    if m != nothing
        return transform_line("$(m.captures[1])$(m.captures[2]){AbstractString}$(m.captures[3])\n")
    end

    m = match(r"^(.+)\s+String\[\](.*)\n$", line)
    if m != nothing
        return transform_line("$(m.captures[1]) AbstractString[]$(m.captures[2])\n")
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