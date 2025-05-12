from shared_mcp_object import mcp
from config import GOOGLE_SEARCH_API_KEY, GOOGLE_CUSTOM_SEARCH_CX, NVIDIA_API_KEY
import httpx
import json
import re
import urllib.parse
import asyncio
from typing_extensions import TypedDict
from typing import List, Dict, Any, Optional
from colorama import Fore, Style
import time 

from langgraph.graph import StateGraph, END, START
from .chat import chat_with_model

TARGET_LLM_MODEL = "palmyra" #Model to be used for the LLM analysis at the end, I'm using Palmyra because this is a specialist model in the Nvidia NIM for medical usage - trained on 70b params.

class VariantInfo(TypedDict, total=False):
    rsid: str
    risk_allele: Optional[str]
    odds_ratio: Optional[float]
    beta: Optional[float]
    pvalue: Optional[float]
    reported_gene: Optional[str]
    mention_count: Optional[int]
    study_ids: Optional[List[str]]
    
class GeneticResearchState(TypedDict):
    condition: str
    rsids: List[str]
    variant_details: Dict[str, VariantInfo]
    web_research_summary: Dict[str, str]
    compiled_report: Dict[str, Any]
    final_llm_summary: Optional[str]
    error: Optional[str]
    current_step_message: Optional[str]


async def find_rsids_from_gwas_catalog(condition: str, max_results: int = 50) -> tuple[List[str], Dict[str, VariantInfo]]:
    """Query GWAS Catalog for rsIDs associated with a condition with improved reliability tracking"""
    #Get studies by disease trait from GWAS catalogue, so it's checking for study IDs to find associations within. Next step will be to get associations for each study ID.
    encoded_condition = urllib.parse.quote(condition)
    study_url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByDiseaseTrait?diseaseTrait={encoded_condition}"
    variant_details: Dict[str, VariantInfo] = {}
    
    # Track which studies mention each rsID. An RSID is basically a genetic marker so we can relate it to a DNA profile later on.
    rsid_mentions: Dict[str, List[str]] = {}
    
    print(f"{Fore.CYAN}Searching GWAS Catalog for {Fore.YELLOW}'{condition}'{Fore.CYAN} studies...{Style.RESET_ALL}")
    
    async with httpx.AsyncClient(timeout=45.0) as client:
        try:
            # Get studies to read through and find associations.
            study_response = await client.get(study_url)
            study_response.raise_for_status()
            study_data = study_response.json()
            
            if '_embedded' not in study_data or 'studies' not in study_data['_embedded']:
                print(f"{Fore.RED}No studies found for '{condition}'{Style.RESET_ALL}")
                return [], {}
            
            studies = study_data['_embedded']['studies']
            print(f"{Fore.GREEN}Found {len(studies)} studies for '{condition}'{Style.RESET_ALL}")
            
            # Grab each association for each study.
            for study in studies[:50]:  # Limit to 50 studies for performance in the demo.
                if 'accessionId' not in study:
                    continue
                    
                acc_id = study['accessionId']
                study_title = study.get('title', 'Unknown study')
                assoc_url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/{acc_id}/associations"
                
                try:
                    assoc_response = await client.get(assoc_url)
                    assoc_response.raise_for_status()
                    assoc_data = assoc_response.json()
                    
                    if '_embedded' in assoc_data and 'associations' in assoc_data['_embedded']:
                        associations = assoc_data['_embedded']['associations']
                        print(f"{Fore.BLUE}Found {len(associations)} associations for study {acc_id}{Style.RESET_ALL}")
                        
                        # Process each association so that can find what genetic markers match a certain condition as a risk factor.
                        for assoc in associations:
                            # Check for loci information which contains the risk alleles, so that will contain the RSID. DNA Files contain RSIDs so this is important for analysis later.
                            if 'loci' in assoc:
                                for locus in assoc['loci']:
                                    if 'strongestRiskAlleles' in locus:
                                        for risk_allele in locus['strongestRiskAlleles']:
                                            if 'riskAlleleName' in risk_allele:
                                                # Extract RSID from format like "rs3764147-G", using regex to match it.
                                                allele_name = risk_allele['riskAlleleName']
                                                rsid_match = re.match(r'(rs\d+)(?:-([A-Z]+))?', allele_name)
                                                
                                                if rsid_match:
                                                    rsid = rsid_match.group(1).lower()
                                                    risk_allele_letter = rsid_match.group(2) if rsid_match.group(2) else None
                                                    
                                                    # Track this rsID mention in a list of mentions. If it has not been seen before, create a new list for it so we can see all the studies it was mentioned in.
                                                    if rsid not in rsid_mentions:
                                                        rsid_mentions[rsid] = []
                                                    if acc_id not in rsid_mentions[rsid]:
                                                        rsid_mentions[rsid].append(acc_id)
                                                    
                                                    # Extract gene info if available in the loci info.
                                                    gene_name = None
                                                    if 'authorReportedGenes' in locus and locus['authorReportedGenes']:
                                                        gene_info = locus['authorReportedGenes'][0]
                                                        if gene_info.get('geneName') and gene_info.get('geneName') != "Unknown":
                                                            gene_name = gene_info.get('geneName')
                                                    
                                                    # Update or create variant info using the VariantInfo class.
                                                    # This is where we will store the rsID, risk allele, p-value, and gene name, because these determine the significance of the variant.
                                                    # We do try to be optimistic and update the variant info if we find a better p-value or more mentions in other studies.
                                                    if rsid in variant_details:
                                                        # Update existing entry with better values if available in another study. Better means more mentions, lower p-value, etc. 
                                                        current = variant_details[rsid]
                                                        
                                                        # Update mention count so we know how many studies it was mentioned in.
                                                        current["mention_count"] = len(rsid_mentions[rsid])
                                                        current["study_ids"] = rsid_mentions[rsid]
                                                        
                                                        # Keep the lowest p-value to be more reliable and optimistic.
                                                        if 'pvalue' in assoc and (current.get('pvalue') is None or float(assoc['pvalue']) < current.get('pvalue')):
                                                            current['pvalue'] = float(assoc['pvalue'])
                                                            
                                                        # Only update gene info if it was unknown before.
                                                        if not current.get('reported_gene') and gene_name:
                                                            current['reported_gene'] = gene_name
                                                            
                                                        # Only update risk allele if it was unknown before.
                                                        if not current.get('risk_allele') and risk_allele_letter:
                                                            current['risk_allele'] = risk_allele_letter
                                                    else:
                                                        # Create new variant info object if it doesn't exist.
                                                        variant_info = {
                                                            "rsid": rsid,
                                                            "risk_allele": risk_allele_letter,
                                                            "pvalue": float(assoc['pvalue']) if 'pvalue' in assoc else None,
                                                            "reported_gene": gene_name,
                                                            "mention_count": len(rsid_mentions[rsid]),
                                                            "study_ids": rsid_mentions[rsid]
                                                        }
                                                        
                                                        # Add create an odds ratio if available. An odds ratio is a measure of association between an exposure and an outcome.
                                                        if 'orPerCopyNum' in assoc:
                                                            try:
                                                                variant_info["odds_ratio"] = float(assoc['orPerCopyNum'])
                                                            except (ValueError, TypeError):
                                                                pass
                                                        
                                                        # Add beta if available. A beta is a measure of effect size in regression analysis.
                                                        if 'betaNum' in assoc:
                                                            try:
                                                                variant_info["beta"] = float(assoc['betaNum'])
                                                            except (ValueError, TypeError):
                                                                pass
                                                                
                                                        variant_details[rsid] = variant_info
                                                        print(f"{Fore.GREEN}Found rsID: {rsid} with risk allele {risk_allele_letter}{Style.RESET_ALL}")
                    
                except Exception as e:
                    print(f"{Fore.RED}Error fetching associations for study {acc_id}: {e}{Style.RESET_ALL}")
                
                # Add a slight delay to avoid rate limiting on the main GWAS API.
                await asyncio.sleep(0.3)
                    
        except Exception as e:
            print(f"{Fore.RED}Error in GWAS Catalog search: {e}{Style.RESET_ALL}")
            return [], {}
    
    # Sort rsIDs by mention count (reliability), then by p-value (significance). This will help us to find the most reliable and significant rsIDs for medical indicators.
    # Using a lambda function for efficient sorting for our reported output.
    sorted_rsids = sorted(
        variant_details.keys(),
        key=lambda x: (
            -variant_details[x].get('mention_count', 0),  # Sort by mention count (descending)
            variant_details[x].get('pvalue', 1.0) if variant_details[x].get('pvalue') is not None else 1.0  # Then by p-value (ascending)
        )
    )
    
    # Print top mentioned rsIDs for debugging and user feedback at server-level. This will NOT be available to a client with no access to the server as it's just terminal output.
    print(f"\n{Fore.YELLOW}Top mentioned rsIDs:{Style.RESET_ALL}")
    for rsid in sorted_rsids[:50]:
        mention_count = variant_details[rsid].get('mention_count', 0)
        p_value = variant_details[rsid].get('pvalue', 'Unknown')
        gene = variant_details[rsid].get('reported_gene', 'Unknown gene')
        print(f"{Fore.GREEN}rsID: {rsid} - mentioned in {mention_count} studies, p-value: {p_value}, gene: {gene}{Style.RESET_ALL}")
    
    # Return the sorted list of rsIDs and the detailed info.
    rsids = sorted_rsids[:max_results]
    
    print(f"{Fore.GREEN}Extracted {len(rsids)} unique rsIDs from GWAS Catalog{Style.RESET_ALL}")
    return rsids, variant_details

#Function to search Google for RSIDs and condition info.
async def search_web_for_rsids_or_condition( 
    search_terms: List[str], 
    is_condition_search: bool = False,
    max_retries: int = 3,
    initial_delay: float = 1.0
) -> Dict[str, str]:
    """Search the web for information about each rsID or a general condition with retry logic."""
    web_info: Dict[str, str] = {}
    if not GOOGLE_SEARCH_API_KEY or not GOOGLE_CUSTOM_SEARCH_CX:
        if is_condition_search and search_terms:
            return {search_terms[0]: "Google Search not configured for condition research."}
        return {term: "Google Search not configured" for term in search_terms}
    if not search_terms:
        return {}
        
    print(f"{Fore.CYAN}Searching web for {len(search_terms)} items...{Style.RESET_ALL}")
    
    # Use a single client session for all requests.
    async with httpx.AsyncClient(timeout=30.0) as client:
        for i, term in enumerate(search_terms):
            query = ""
            if is_condition_search:
                query = f"{term} genetic basis OR overview OR symptoms OR genetic risk factors"
                print(f"{Fore.BLUE}Searching web for condition: {term} ({i+1}/{len(search_terms)}){Style.RESET_ALL}")
            else: # rsID search.
                query = f"{term} genomics significance OR association OR function"
                print(f"{Fore.BLUE}Searching for rsID: {term} ({i+1}/{len(search_terms)}){Style.RESET_ALL}")

            params = {"key": GOOGLE_SEARCH_API_KEY, "cx": GOOGLE_CUSTOM_SEARCH_CX, "q": query}
            
            current_retry_delay = initial_delay
            for attempt in range(max_retries):
                try:
                    resp = await client.get("https://www.googleapis.com/customsearch/v1", params=params)
                    
                    if resp.status_code == 429: # Too Many Requests error
                        if attempt == max_retries - 1: # Last attempt
                            print(f"{Fore.RED}Rate limit hit (429) for {term}. Max retries reached. Giving up.{Style.RESET_ALL}")
                            resp.raise_for_status() # Raise error to be caught below and carry on.
                        
                        print(f"{Fore.YELLOW}Rate limit hit (429) for {term}. Retrying in {current_retry_delay:.1f}s... (Attempt {attempt+1}/{max_retries}){Style.RESET_ALL}")
                        await asyncio.sleep(current_retry_delay)
                        current_retry_delay *= 2 
                        continue # Retry the request.
                    
                    resp.raise_for_status() # Raise for other HTTP errors (4xx, 5xx) that may come up.
                    search_results = resp.json()
                    snippets = [item.get("snippet", "") for item in search_results.get("items", [])[:3]]
                    web_info[term] = " ".join(filter(None, snippets)).strip() if snippets else f"No significant web results found for {term}."
                    break # Success, exit retry loop
                except httpx.HTTPStatusError as e:
                    # For 429, this is hit if max_retries is reached and raise_for_status() is called.
                    web_info[term] = f"Error searching web for {term} after {attempt+1} attempts: {str(e)}"
                    print(f"{Fore.RED}HTTP Error searching for {term} (Status: {e.response.status_code}): {e}{Style.RESET_ALL}")
                    break 
                except Exception as e:
                    web_info[term] = f"General error searching web for {term}: {str(e)}"
                    print(f"{Fore.RED}General error searching for {term}: {e}{Style.RESET_ALL}")
                    break 
            
            # General delay between different search terms to prevent rate limiting.
            if i < len(search_terms) - 1:
                await asyncio.sleep(1.0) 
    
    return web_info

async def web_research_node(state: GeneticResearchState) -> Dict[str, Any]:
    print(f"\n{Fore.CYAN}=== Node: Web Research ==={Style.RESET_ALL}")
    rsids = state.get("rsids", [])
    condition = state["condition"]
    web_summary = {}

    if rsids:
        print(f"{Fore.BLUE}Performing web research for {len(rsids)} rsIDs related to {condition}.{Style.RESET_ALL}")
        web_summary = await search_web_for_rsids_or_condition(rsids, is_condition_search=False)
        message = f"Completed web research for {len(rsids)} rsIDs."
    elif condition:
        print(f"{Fore.BLUE}No rsIDs found from GWAS. Performing general web research for condition: {condition}.{Style.RESET_ALL}")
        condition_web_info = await search_web_for_rsids_or_condition([condition], is_condition_search=True)
        web_summary = condition_web_info 
        message = f"Completed general web research for condition '{condition}' as no specific rsIDs were found by GWAS."
    else:
        message = "Skipping web research, no rsIDs and no condition specified."
        print(f"{Fore.YELLOW}{message}{Style.RESET_ALL}")
        return {"web_research_summary": {}, "current_step_message": message}
        
    return {"web_research_summary": web_summary, "current_step_message": message}

def summarise_findings_core(condition: str, rsids: List[str], variant_details: Dict[str, VariantInfo], 
                          web_info: Dict[str, str]) -> Dict[str, Any]:
    """Compile all research findings into a single report"""
    return {
        "condition_researched": condition,
        "identified_rsids": rsids,
        "variant_details": variant_details,
        "web_research_summary": web_info,
    }

async def llm_analyse_report(report_json_str: str) -> Optional[str]:
    """Generate a summary of the report using the LLM"""
    if not NVIDIA_API_KEY:
        return "LLM summarisation skipped (NVIDIA_API_KEY not configured)."
    if not report_json_str:
        return "No report data to summarise."

    print(f"{Fore.CYAN}Generating LLM summary with {Fore.YELLOW}{TARGET_LLM_MODEL}{Style.RESET_ALL}")

    prompt = f"""Analyze the following genetic research report (in JSON format) about a medical condition.

Your task is to create a structured summary of the key findings. Format your response as JSON with the following structure:
{{
  "condition_summary": "Brief overview of the condition",
  "genetic_findings": {{
    "overview": "General overview of the genetic associations found",
    "variants_by_significance": [
      {{
        "rsid": "rsID",
        "risk_allele": "risk allele letter",
        "significance": "p-value",
        "effect_size": "odds ratio or beta value",
        "reported_gene": "associated gene",
        "biological_impact": "interpretation of the variant's role"
      }}
    ]
  }},
  "clinical_relevance": "How these findings might impact clinical understanding or treatment",
  "research_implications": "Suggestions for further research based on findings"
}}

Focus especially on:
- Identifying variants with strongest associations (lowest p-values)
- Noting which risk alleles increase disease risk (odds_ratio > 1 or beta > 0)
- Explaining the biological significance of mentioned genes
- Drawing connections between multiple variants if patterns exist

JSON Report:
{report_json_str}

Summary:
"""
    try:
        result = await chat_with_model(message=prompt, model=TARGET_LLM_MODEL)
        return result
    except Exception as e:
        print(f"{Fore.RED}Error during LLM summarisation: {e}{Style.RESET_ALL}")
        return f"LLM summarisation failed: {str(e)}"

#Graph node functions
async def gwas_catalog_search_node(state: GeneticResearchState) -> Dict[str, Any]:
    print(f"\n{Fore.CYAN}=== Node: GWAS Catalog Search ==={Style.RESET_ALL}")
    condition = state["condition"]
    rsids, variant_details = await find_rsids_from_gwas_catalog(condition)
    if not rsids:
        message = f"No rsIDs found in GWAS Catalog for '{condition}'. Subsequent steps might yield limited results."
        print(f"{Fore.RED}{message}{Style.RESET_ALL}")
        return {"rsids": [], "variant_details": {}, "current_step_message": message}
    print(f"{Fore.GREEN}Found {len(rsids)} rsIDs for {condition}{Style.RESET_ALL}")
    return {
        "rsids": rsids, 
        "variant_details": variant_details,
        "current_step_message": f"Found {len(rsids)} potential rsIDs from GWAS Catalog."
    }

def compile_report_node(state: GeneticResearchState) -> Dict[str, Any]:
    print(f"\n{Fore.CYAN}=== Node: Compile Report ==={Style.RESET_ALL}")
    report = summarise_findings_core(
        state["condition"],
        state.get("rsids", []),
        state.get("variant_details", {}),
        state.get("web_research_summary", {})
    )
    print(f"{Fore.GREEN}Report compiled successfully.{Style.RESET_ALL}")
    return {"compiled_report": report, "current_step_message": "Report compiled."}

async def llm_final_summary_node(state: GeneticResearchState) -> Dict[str, Any]:
    print(f"\n{Fore.CYAN}=== Node: LLM Final Summary ==={Style.RESET_ALL}")
    compiled_report = state.get("compiled_report")
    if not compiled_report:
        message = "Skipping LLM summary, no compiled report."
        print(f"{Fore.YELLOW}{message}{Style.RESET_ALL}")
        return {"final_llm_summary": message, "current_step_message": message}
    
    report_str = json.dumps(compiled_report, indent=2)
    llm_summary = await llm_analyse_report(report_str)
    print(f"{Fore.GREEN}LLM summary generated.{Style.RESET_ALL}")
    return {"final_llm_summary": llm_summary, "current_step_message": "LLM summary generated."}

#Workflow graph initialisation of the object based on the standard genetic research state that was defined above.
#This graph maintains the state of each part throughout the workflow as it runs through the nodes (steps).
workflow_graph = StateGraph(GeneticResearchState)

#Create the node definitions, which are the steps.
workflow_graph.add_node("gwas_search", gwas_catalog_search_node)
workflow_graph.add_node("web_research", web_research_node)
workflow_graph.add_node("compile_report", compile_report_node)
workflow_graph.add_node("llm_summary", llm_final_summary_node)

#Connect up the nodes, so we know what order to run it in for a reliable output. This is important as a medical research tool, as defined research methods are crucial for reliability.
workflow_graph.add_edge(START, "gwas_search")
workflow_graph.add_edge("gwas_search", "web_research")
workflow_graph.add_edge("web_research", "compile_report")
workflow_graph.add_edge("compile_report", "llm_summary")
workflow_graph.add_edge("llm_summary", END)

#Compile the graph to run the workflow.
app_graph = workflow_graph.compile()

# MCP registration
@mcp.tool()
async def deep_genetic_research(condition: str) -> dict:
    """
    Research genetic associations for a medical condition using GWAS Catalog data.
    
    Args:
        condition: The medical condition to research.
    Returns:
        A dictionary with research findings and LLM summary.
    """
    print(f"\n{Fore.GREEN}╔══════════════════════════════════════════════╗{Style.RESET_ALL}")
    print(f"{Fore.GREEN}║ Starting Genetic Research for: {Fore.YELLOW}{condition}{Fore.GREEN} ║{Style.RESET_ALL}")
    print(f"{Fore.GREEN}╚══════════════════════════════════════════════╝{Style.RESET_ALL}")
    
    initial_state: GeneticResearchState = {
        "condition": condition,
        "rsids": [],
        "variant_details": {},
        "web_research_summary": {},
        "compiled_report": {},
        "final_llm_summary": None,
        "error": None,
        "current_step_message": "Workflow initiated."
    }
    
    final_state = await app_graph.ainvoke(initial_state)
    
    print(f"\n{Fore.GREEN}╔══════════════════════════════════════════════╗{Style.RESET_ALL}")
    print(f"{Fore.GREEN}║ Research Complete: {Fore.YELLOW}{final_state.get('current_step_message')}{Fore.GREEN} ║{Style.RESET_ALL}")
    print(f"{Fore.GREEN}╚══════════════════════════════════════════════╝{Style.RESET_ALL}")
    
    return {
        "condition": final_state.get("condition"),
        "compiled_report": final_state.get("compiled_report"),
        "llm_summary": final_state.get("final_llm_summary"),
        "error": final_state.get("error")
    }