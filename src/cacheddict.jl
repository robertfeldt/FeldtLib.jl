struct CachedDict
    cachefile::String
    d
    CachedDict(fn::Function, filename::String) = begin
        new(filename, fromcacheorcreate(fn, filename))
    end
end
Base.haskey(c::CachedDict, key) = haskey(c.d, key)
Base.getindex(c::CachedDict, key) = c.d[key]
Base.setindex!(c::CachedDict, value, key) = update!(c, value, key, true)

function cachedict(c::CachedDict)
    if !isdir(dirname(c.cachefile))
        mkpath(dirname(c.cachefile))
    end
    serialize(c.cachefile, c.d)
end

function update!(c::CachedDict, value, key, writetodisc = false)
    c.d[key] = value
    writetodisc && cachedict(c)
    value
end

function fromcacheorcreate(fn, filename)
    !isfile(filename) && return fn()
    try
        deserialize(filename)
    catch _err
        return fn()
    end
end
