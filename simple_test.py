#!/usr/bin/env python3
"""Simple MCP tool test to understand result format."""

import sys
import asyncio
from pathlib import Path

# Add src to path
sys.path.append('src')
from server import mcp

async def test_server_info():
    """Test server info and see result format."""
    print("Testing get_server_info...")
    try:
        result = await mcp.call_tool("get_server_info", {})
        print(f"Result type: {type(result)}")
        print(f"Result dir: {[attr for attr in dir(result) if not attr.startswith('_')]}")

        if hasattr(result, 'content'):
            print(f"Content: {result.content}")
        if hasattr(result, 'result'):
            print(f"Result: {result.result}")
        if hasattr(result, 'value'):
            print(f"Value: {result.value}")

        return True
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    asyncio.run(test_server_info())