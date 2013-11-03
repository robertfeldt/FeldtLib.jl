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

def filter_latest_changed_files(filenames, numLatestChangedToInclude = 1)
  filenames.sort_by{ |f| File.mtime(f) }[-numLatestChangedToInclude, numLatestChangedToInclude]
end

desc "Run only the latest changed test file"
task :t do
  latest_changed_test_file = filter_latest_changed_files Dir["test/test*.jl"]
  sh "julia -L #{MainFile} -L test/helper.jl #{latest_changed_test_file.first}"
end

task :st => :runslowtest

task :default => :runtest