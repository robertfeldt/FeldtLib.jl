# We need some simple way that julia packages can "call home" and indicate use.
# This is needed since our funding agencies wants proof that our software
# and solutions are of use in society.
function async_send_to_server_tcp(message, url = "127.0.0.1", port = 2001)
  @async begin
    socket = connect(url, port)
    println(socket, message)
    close(socket)
  end
end

function start_test_server(port = 2001)
  @async begin
    server = listen(port)
    while true
      socket = accept(server)
      tstr = strftime("%Y%m%d %H:%M.%S", time())
      println("$(tstr), Server received: ", readline(socket), "\n")
    end
  end
end