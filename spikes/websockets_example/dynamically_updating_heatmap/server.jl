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
    try
      write(client, jsonmsg)
    catch error
      println(error)
    end
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

session_timestamp() = begin
  t = time()
  strftime("%Y-%m-%dT%H:%M:%S.", t) * string(int(round(t % 1, 3) * 1000)) * "Z"
end

function rand_test_outcome()
  pr = rand()
  if pr < 0.10
    -1 # Failure
  elseif pr < 0.95
    1 # Pass
  else
    0 # Missing info
  end
end

# For now the processing is just to send some random data every now and then...
function processing_func()
  numtestcases = 10
  session = session_timestamp()
  @sync @parallel for tc in 1:numtestcases
    sleep(1 + 5 * rand()) # Fake that processing takes time...
    msgtype = (rand() < 0.80) ? "dataupdate" : "unspecified"
    send({"type" => msgtype, "name" => "tom", 
      "session" => session, "testcase" => tc,
      "value" => rand_test_outcome()})
  end
end

start_websocket_handler_while_processing(processing_func)