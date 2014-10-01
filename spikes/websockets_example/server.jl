using HttpServer
using WebSockets
using JSON

decodeMessage(msg) = bytestring(msg)

onepage = readall("./index.html")
httph = HttpHandler() do req::Request, res::Response
  Response(onepage)
end

global connections = Dict{Int,WebSocket}()

function send(msg)
  global connections
  jsonmsg = JSON.json(msg)
  if length(connections) == 0
    println("No one to send to! ", jsonmsg)
  end
  for (id, client) in connections
    println("Sending to ", client, ": ", jsonmsg)
    write(client, jsonmsg)
  end
end

# Send some random data every now and then
@async begin
  while true
    sleep(10 + 5 * rand())
    send({"data" => rand(), "matrix" => "A"})
  end
end

wsh = WebSocketHandler() do req, client
  global connections
  println("Got request from client: ", client)
  connections[client.id] = client
  while true
    msg = decodeMessage(read(client))
    println("Received from ", client.id, ": ", msg)
    try
      jsonmsg = JSON.parse(msg)
      handle_json_message(req, client, jsonmsg)
    catch error
      error("Could not parse and/or handle the msg: ", msg)
    end
  end
end
println("Now starting server!")
server = Server(httph, wsh)
run(server, 8084)
