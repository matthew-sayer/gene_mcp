from shared_mcp_object import mcp
from config import NVIDIA_API_KEY
import requests
from colorama import Fore, Style

@mcp.tool()
async def dna_generator(sequence):
    """
    DNA generator using the EVO2-40B Model.
    Your goal is to generate a DNA sequence based on the input sequence. The input should ideally be bad genetics that require improvement.
    The model will generate a new sequence that is an improvement over the input.

    Args:
        sequence: Input DNA sequence (EXAMPLE: "ACTGACTGACTGACTG")
    
    Returns:
        Dictionary with generated DNA data.

        EXAMPLE:
        {
  "sequence": "TCGCTGATTGCTAGCACAGCAGTCTGAGATCAAACTGCAAGGCGGCAGCGAGGCTGGGGGAGGGGCGCCCACCATTGCCCAGGCTTGCTTAGGTAAACAAAG",
  "logits": null,
  "sampled_probs": [
    0.9759454131126404,
    0.991154134273529,
    0.99029141664505,
    0.9881120920181274,
    0.9856293797492981
  ],
  "elapsed_ms": 3183
}
    """
    print(f"{Fore.GREEN}╔═══════════════════════════════════════════╗{Style.RESET_ALL}")
    print(f"{Fore.GREEN}║ Starting DNA Generation with NVIDIA API   ║{Style.RESET_ALL}")
    print(f"{Fore.GREEN}║ There may be a queue, please be patient   ║{Style.RESET_ALL}")
    print(f"{Fore.GREEN}╚═══════════════════════════════════════════╝{Style.RESET_ALL}")
    
    print(f"{Fore.CYAN}Input sequence: {Fore.YELLOW}{sequence}{Style.RESET_ALL}")
    
    try:
        print(f"{Fore.BLUE}Making request to NVIDIA EVO2-40B API...{Style.RESET_ALL}")
        response = requests.post(
            url="https://health.api.nvidia.com/v1/biology/arc/evo2-40b/generate",
            headers={"Authorization": f"Bearer {NVIDIA_API_KEY}"},
            json={
                "sequence": sequence,
                "num_tokens": 16,
                "top_k": 1,
                "enable_sampled_probs": True,
            }
        )
        
        print(f"{Fore.BLUE}API response received with status code: {response.status_code}{Style.RESET_ALL}")
        
        if response.status_code == 200:
            result = response.json()
            generated_sequence = result.get("sequence", "")
            
            print(f"{Fore.GREEN}Successfully generated DNA sequence:{Style.RESET_ALL}")
            print(f"{Fore.YELLOW}Generated: {generated_sequence}{Style.RESET_ALL}")
            print(f"{Fore.BLUE}Length: {len(generated_sequence)} bases{Style.RESET_ALL}")
            
            print(f"\n{Fore.GREEN}╔═══════════════════════════════════════════╗{Style.RESET_ALL}")
            print(f"{Fore.GREEN}║ DNA Generation Complete                   ║{Style.RESET_ALL}")
            print(f"{Fore.GREEN}╚═══════════════════════════════════════════╝{Style.RESET_ALL}")

            return result
        else:
            print(f"{Fore.RED}API error: {response.status_code}{Style.RESET_ALL}")
            print(f"{Fore.RED}Details: {response.text[:200]}{Style.RESET_ALL}")
            return {"error": f"API error: {response.status_code}", "details": response.text}
    
    except Exception as e:
        print(f"{Fore.RED}Error during DNA generation: {str(e)}{Style.RESET_ALL}")
        return {"error": f"Error during DNA generation: {str(e)}"}