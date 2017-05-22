# We want to extract the trees in a BART tree created by the bartMachine R package.
# The java object is avaialble in bart_machine$java_bart_machine and we can serialize
# that with .jserialize in the rJava package. Example code:
#
# library(bartMachine)
# n = 10
# p = 1
# X = data.frame(matrix(runif(n * p), ncol = p))
# y = 10 * X[ ,1]
# bart_machine = bartMachine(X,y, num_trees = 2, serialize = TRUE)
# str = .jserialize(bart_machine$java_bart_machine)
# length(str)
# fh <- file("bart_machine_java_serialized.bin", "wb")
# writeBin(str, fh)
# close(fh)

# Now we want to parse out the tree from this serialized java object.
# We create a small java serialized object parser in Julia. It will not be
# complete but will be focused for this particular use case.
#
function parse_java_serialized_object(o::Vector{UInt8})
    p = 1
    ensure_value(o, UInt8[0xAC, 0xED], p, "STREAM_MAGIC not matched")
end

o = open(fh -> read(fh), "bart_machine_java_serialized.bin")

function ensure_value(o::Vector{UInt8}, bytes::Vector{UInt8}, pos::Int = 1, message = "")
    allmatched = all(i -> o[pos + i - 1] == bytes[i], 1:length(bytes))
    if !allmatched
        error(message * " (at position $pos)")
    else
        return allmatched
    end
end

# The basic format of the seralized object is as a series of objects that adds a few different 
# but related methods/vars and then has an array of bartMachineTreeNode objects to represent the
# root of each BART tree. Each node has left, right, split dimension and split value, plus flags 
# and intermediate values. So basically we just need to skip ahead to the first splitAttributeM
# int and then read its double splitValue and just recurse.
# Seems fairly simple in theory. :)

# Example of how to parse these objects: http://www.javaworld.com/article/2072752/the-java-serialization-algorithm-revealed.html
# Format spec: https://docs.oracle.com/javase/7/docs/platform/serialization/spec/protocol.html