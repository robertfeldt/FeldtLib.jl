include("helper.jl")

path_to_this_file = Base.source_path()
testdir = join(split(path_to_this_file, "/")[1:end-1], "/")
srcdir = join([testdir, "../src"], "/")

# First run the base tests
include("runbasetests.jl")

# Then run
AutoTest.run_all_tests_in_dir(testdir; regexpTestFiles = r"^test.*\.jl$")