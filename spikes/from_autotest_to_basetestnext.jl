testfile = ARGS[1]

totab4(str) = " " ^ (div(length(str), 2) * 4)

newlines = open(testfile, "r") do fh
    map(readlines(fh)) do line
        m1 = match(r"(\s*)describe(\"(.*)\") do", line)
        m2 = match(r"(\s*)test(\"(.*)\") do", line)
        if m1 != nothing
            "$(totab4(m1.captures[1]))@testset \"$(m1.captures[2])\" begin"
        elseif m2 != nothing
            "$(totab4(m2.captures[1]))@testset \"$(m2.captures[2])\" begin"
        else
            line
        end
    end

end

open(testfile * ".alt", "w") do fh
    println(fh, join(newlines, "\n"))
end