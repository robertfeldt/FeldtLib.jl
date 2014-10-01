using HttpServer
using WebSockets
using JSON

decodeMessage(msg) = bytestring(msg)

global jsonwebsocketconnections = Dict{Int,WebSocket}()

function send(msg)
  global jsonwebsocketconnections
  jsonmsg = JSON.json(msg)
  if length(jsonwebsocketconnections) == 0
    println("No one to send to! ", jsonmsg)
  end
  for (id, client) in jsonwebsocketconnections
    println("Sending to ", id, ": ", jsonmsg)
    write(client, jsonmsg)
  end
end

function start_websocket_handler_while_processing(processingFunc; mainpage = "index.html", port = 8084)

  onepage = readall(mainpage)
  httph = HttpHandler() do req::Request, res::Response
    Response(onepage)
  end

  @async begin
    while true
      processingFunc()
    end
  end

  wsh = WebSocketHandler() do req, client
    global jsonwebsocketconnections
    println("Got request from client: ", client)
    jsonwebsocketconnections[client.id] = client
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
end

global currentrow = 1

# For now the processing is just to send some random data every now and then...
function processing_func()
  global currentrow
  println("Starting processing of row ", currentrow)
  numcols = 10
  @sync @parallel for col in 1:numcols
    sleep(5 + 10 * rand()) # Fake that processing takes time...
    send({"matrix" => "A", "action" => "set_value", 
      "value" => rand(), "row" => currentrow, "col" => col})
  end
  currentrow += 1
end

start_websocket_handler_while_processing(processing_func)