testfile = ARGS[1]

totab4(str) = " " ^ (div(length(str), 2) * 4)

function transform_line(line::AbstractString)
    m = match(r"^(\s*)describe\(\"(.*)\"\)\s+do(.*)\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@testset \"$(m.captures[2])\" begin$(m.captures[3])\n"
    end

    m = match(r"^(\s*)test\(\"(.*)\"\)\s+do(.*)\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@testset \"$(m.captures[2])\" begin$(m.captures[3])\n"
    end

    m = match(r"^(\s*)@repeat\s+test\(\"(.*)\"\)\s+do(.*)\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@testset \"$(m.captures[2])\" for i in 1:NumReps $(m.captures[3])\n"
    end

    m = match(r"^(\s*)@check(\s+)(.+)\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))@test$(m.captures[2])$(m.captures[3])\n"
    end

    m = match(r"^(\s*)([^\s].*)\n", line)
    if m != nothing
        return "$(totab4(m.captures[1]))$(m.captures[2])\n"
    end
    
    line
end

newlines = open(testfile, "r") do fh
    map(transform_line, readlines(fh))
end

open(testfile * ".alt", "w") do fh
    println(fh, join(newlines, ""))
end