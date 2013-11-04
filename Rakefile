MainFile = "src/FeldtLib.jl"

def runtestfile(filename)
  if File.exist?(filename)
    sh "julia -L #{MainFile} #{filename}"
  else
    # puts "No test file named #{filename} found."
  end
end

task :runtest do
  runtestfile "test/runtests.jl"
end

task :runtest do
  runtestfile "test/runtests.jl"
end

desc "Run all tests in the dir named 'test'"
task :alltests do
  sh "julia -L #{MainFile} -L test/helper.jl -e 'AutoTest.run_all_tests_in_dir(\"test\")'"
end

def filter_latest_changed_files(filenames, numLatestChangedToInclude = 1)
  filenames.sort_by{ |f| File.mtime(f) }[-numLatestChangedToInclude, numLatestChangedToInclude]
end

desc "Run only the latest changed test file"
task :t do
  latest_changed_test_file = filter_latest_changed_files Dir["test/test*.jl"]
  sh "julia -L #{MainFile} -L test/helper.jl -e 'AutoTest.run_tests_in_file(\"#{latest_changed_test_file.first}\")'"
end

task :st => :runslowtest
task :at => :alltests
task :default => :at