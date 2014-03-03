using Homebrew

# First time there is an error. Works the 2nd time though. Not sure why, but 
# seems to be a problem with GnuTLS.
try
  using GnuTLS
catch
end

using GnuTLS

using Requests

r = get("http://rubygems.org/profiles/feldt")

m = match(r"<div><strong>([0-9,]+)</strong> all time</div>", r.data)

if m == nothing

  println("ERROR! Could not identify a dl count!")
  exit(-1)

else

  count = replace(m.captures[1], ",", "")
  num_dls = int(count)
  rounded = int( floor(num_dls/1000) * 1000 )
  datestr = strftime("%B %e %Y", time())
  print("Feldt's Ruby gems has been downloaded ($(datestr)) more than $(rounded) times ($(num_dls))")

end