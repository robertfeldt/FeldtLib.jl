module DebugPrintln

export dprintln, debug, nodebug

DebugFlag = true

function set_debug_flag(value)
	global DebugFlag
	DebugFlag = value
end
nodebug() = set_debug_flag(false)
debug() = set_debug_flag(true)

function dprintln(args...)
	if DebugPrintln.DebugFlag == true
		println(args...)
	end
end

end