from fastapi import FastAPI, Request, Depends
from mcp.server.sse import SseServerTransport
from starlette.routing import Mount
import uvicorn
import dotenv
from shared_mcp_object import mcp, register_all_tools

dotenv.load_dotenv()

tools_count = register_all_tools()
print(f"Registered tools from {tools_count} modules")

app = FastAPI(docs_url=None, redoc_url=None)
sse = SseServerTransport("/messages/")
app.router.routes.append(Mount("/messages", app=sse.handle_post_message))


@app.get("/sse", tags=["MCP"])
async def handle_sse(request: Request):
    
    async with sse.connect_sse(request.scope, request.receive, request._send) as (
        read_stream,
        write_stream,
    ):
        init_options = mcp._mcp_server.create_initialization_options()

        await mcp._mcp_server.run(
            read_stream,
            write_stream,
            init_options,
        )


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)