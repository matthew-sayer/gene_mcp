from shared_mcp_object import mcp
from openai import OpenAI
from config import NVIDIA_API_KEY

@mcp.tool()
async def chat_with_model(message: str, model: str = "deepseek") -> str:
    """Chat with an AI model
    
    Args:
        message (str): The message to send to the AI.
        model (str, optional): The model to use. Options:
            - "deepseek"
            - "nemotron"
            - "palmyra" (default as this is a medical specialised model)
        
    Returns:
        str: The AI's response.
    """
    model_map = {
        "deepseek": "deepseek-ai/deepseek-r1",
        "nemotron": "nvidia/llama-3.1-nemotron-ultra-253b-v1",
        "palmyra": "writer/palmyra-med-70b-32k"
    }
    
    actual_model = model_map.get(model.lower(), model)
    
    client = OpenAI(
        base_url="https://integrate.api.nvidia.com/v1",
        api_key=NVIDIA_API_KEY
    )

    completion = client.chat.completions.create(
        model=actual_model,
        messages=[{"role":"user","content":message}],
        temperature=0.6,
        top_p=0.7,
        max_tokens=4096,
        stream=True
    )

    concatenated_response = ""
    for chunk in completion:
        if chunk.choices[0].delta.content:
            concatenated_response += chunk.choices[0].delta.content

    return concatenated_response