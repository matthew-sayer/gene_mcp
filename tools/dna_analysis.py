import asyncio
import json
import math
import pandas as pd
from typing import List, Dict, Any, Tuple, Optional

from colorama import Fore, Style
from shared_mcp_object import mcp

from .deep_research import VariantInfo 


async def load_patient_dna_from_file(file_path: str) -> Dict[str, Tuple[str, str]]:
    """
    Loads patient DNA data from a CSV file.
    Expected columns: rsid, allele1, allele2 (chromosome and position are optional but good to have).
    Returns a dictionary mapping rsID to (allele1, allele2).
    """
    print(f"{Fore.CYAN}Loading patient DNA data from: {Fore.YELLOW}{file_path}{Style.RESET_ALL}")
    patient_genotypes: Dict[str, Tuple[str, str]] = {}
    try:
        df = pd.read_csv(file_path, sep=',', comment='#', low_memory=False, keep_default_na=False, na_values=['', 'NA', 'N/A', 'NaN', '0'])
        
        required_cols = ['rsid', 'allele1', 'allele2'] 
        current_columns = df.columns.str.lower().tolist()
        
        # Normalize column names for checking
        df.columns = df.columns.str.lower().str.strip()

        # Check for common variations of rsid
        if 'rsid' not in df.columns:
            if 'rs id' in df.columns:
                df.rename(columns={'rs id': 'rsid'}, inplace=True)
            elif 'snp id' in df.columns: # Another common variant
                df.rename(columns={'snp id': 'rsid'}, inplace=True)
        
        if not all(col in df.columns for col in required_cols):
            missing_cols = [col for col in required_cols if col not in df.columns]
            raise ValueError(f"Patient DNA file is missing required columns (expected: {required_cols}). Found columns: {current_columns}. Missing: {missing_cols}")

        print(f"{Fore.BLUE}Processing {len(df)} rows from patient DNA file...{Style.RESET_ALL}")
        for index, row in df.iterrows():
            try:
                rsid = str(row['rsid']).strip().lower()
                allele1 = str(row['allele1']).strip().upper()
                allele2 = str(row['allele2']).strip().upper()
                
                if not rsid or rsid == 'nan': 
                    if index < 5: 
                        print(f"{Fore.YELLOW}Warning: Skipping row {index+2} due to missing rsID.{Style.RESET_ALL}")
                    continue

                patient_genotypes[rsid] = (allele1, allele2)
            except Exception as e_row:
                if index < 5: 
                    print(f"{Fore.RED}Error processing row {index+2} in patient DNA file: {e_row}. Row data: {row.to_dict()}{Style.RESET_ALL}")
                continue 
            
        print(f"{Fore.GREEN}Successfully loaded {len(patient_genotypes)} genotypes for the patient.{Style.RESET_ALL}")
        if not patient_genotypes:
            print(f"{Fore.YELLOW}Warning: No genotypes loaded. The file might be empty, in an unexpected format, or all rows had issues.{Style.RESET_ALL}")
        return patient_genotypes
    except FileNotFoundError:
        print(f"{Fore.RED}Error: Patient DNA file not found at {file_path}{Style.RESET_ALL}")
        raise
    except pd.errors.EmptyDataError:
        print(f"{Fore.RED}Error: Patient DNA file is empty at {file_path}{Style.RESET_ALL}")
        raise
    except Exception as e:
        print(f"{Fore.RED}Error loading patient DNA data: {e}{Style.RESET_ALL}")
        raise

@mcp.tool(
    name="assess_patient_genetic_risk_profile",
    description="Assesses a patient's DNA against a pre-compiled genetic risk profile for specified conditions."
)
async def assess_patient_genetic_risk_profile(
    patient_dna_file_path: str,
    condition_risk_profiles: Dict[str, Dict[str, VariantInfo]] 
) -> Dict[str, Any]:
    """
    Assesses patient's genetic data against pre-compiled risk profiles for specified conditions.

    NEVER use FAKE RSIDs, like rs-1234567, or similar unbelievable sequences that you might hallucinate. If you are not sure, DO THE DEEP RESEARCH BEFOREHAND USING OUR MCP TOOLS.

    Args:
        patient_dna_file_path: Path to the patient's DNA data file (CSV format).
                               Expected columns: rsid, allele1, allele2.
        condition_risk_profiles: A dictionary where keys are condition names (e.g., "autism")
                                 and values are dictionaries of variant details (Dict[str, VariantInfo])
                                 for that condition, typically from the deep_genetic_research tool's
                                 compiled_report.variant_details.
                                 Example: {"autism": {"rs123": {"risk_allele": "A", "odds_ratio": 1.5, ...}}}

    Returns:
        A dictionary where keys are condition names. Each value is a dictionary containing:
        - 'polygenic_risk_score': Calculated PRS for the condition.
        - 'risk_variants_found_count': Number of risk variants present in the patient.
        - 'protective_variants_found_count': Number of protective variants present in the patient.
        - 'risk_variants_details': List of specific risk variants present in the patient.
        - 'protective_variants_details': List of specific protective variants present in the patient.
        - 'summary_notes': Notes on the analysis.
    """
    print(f"{Fore.GREEN}╔═══════════════════════════════════════════════════════════════╗{Style.RESET_ALL}")
    print(f"{Fore.GREEN}║ Starting Patient Genetic Risk Assessment for conditions: {Fore.YELLOW}{', '.join(condition_risk_profiles.keys())}{Fore.GREEN} ║{Style.RESET_ALL}")
    print(f"{Fore.GREEN}╚═══════════════════════════════════════════════════════════════╝{Style.RESET_ALL}")

    analysis_results: Dict[str, Any] = {}

    try:
        patient_genotypes = await load_patient_dna_from_file(patient_dna_file_path)
        if not patient_genotypes:
            error_message = "Failed to load patient DNA data or no valid genotypes found."
            for condition_name in condition_risk_profiles.keys():
                analysis_results[condition_name] = {"error": error_message, "status": "Patient DNA Load Failed"}
            return analysis_results
    except Exception as e:
        error_message = f"Failed to load patient DNA: {str(e)}"
        for condition_name in condition_risk_profiles.keys():
            analysis_results[condition_name] = {"error": error_message, "status": "Patient DNA Load Failed"}
        return analysis_results

    for condition_name, condition_variant_details in condition_risk_profiles.items():
        print(f"\n{Fore.CYAN}--- Assessing condition: {Fore.YELLOW}{condition_name}{Fore.CYAN} ---{Style.RESET_ALL}")
        
        if not isinstance(condition_variant_details, dict):
            print(f"{Fore.RED}Error: Risk profile for '{condition_name}' is not in the expected format (should be a dictionary). Skipping.{Style.RESET_ALL}")
            analysis_results[condition_name] = {
                'error': f"Invalid risk profile format for {condition_name}.",
                'status': "Profile Error"
            }
            continue

        if not condition_variant_details:
            print(f"{Fore.YELLOW}No genetic risk factors provided in the profile for '{condition_name}'. Skipping detailed analysis.{Style.RESET_ALL}")
            analysis_results[condition_name] = {
                'polygenic_risk_score': 0.0,
                'risk_variants_found_count': 0,
                'protective_variants_found_count': 0,
                'risk_variants_details': [],
                'protective_variants_details': [],
                'summary_notes': [f"No genetic risk factors provided in the profile for {condition_name}."],
                'status': "No Profile Variants"
            }
            continue

        current_condition_score = 0.0
        patient_risk_variants_details = []
        patient_protective_variants_details = []
        condition_notes = []
        
        for rsid_profile, variant_info in condition_variant_details.items():
            rsid_lookup = rsid_profile.lower() 
            
            if rsid_lookup in patient_genotypes:
                patient_allele1, patient_allele2 = patient_genotypes[rsid_lookup]
                
                effect_allele = variant_info.get('risk_allele') 
                odds_ratio = variant_info.get('odds_ratio')

                if effect_allele and len(effect_allele) == 1 and effect_allele in "AGCT":
                    if odds_ratio is not None and odds_ratio > 0:
                        patient_has_effect_allele = (patient_allele1 == effect_allele or patient_allele2 == effect_allele)
                        
                        if patient_has_effect_allele:
                            try:
                                score_contribution = math.log(odds_ratio) 
                                current_condition_score += score_contribution

                                variant_detail_entry = {
                                    'rsid': rsid_profile,
                                    'patient_genotype': f"{patient_allele1}/{patient_allele2}",
                                    'effect_allele': effect_allele,
                                    'odds_ratio': odds_ratio,
                                    'reported_gene': variant_info.get('reported_gene'),
                                    'p_value': variant_info.get('pvalue'),
                                    'profile_mention_count': variant_info.get('mention_count', 1) 
                                }
                                if odds_ratio > 1:
                                    patient_risk_variants_details.append(variant_detail_entry)
                                elif odds_ratio < 1: 
                                    patient_protective_variants_details.append(variant_detail_entry)
                            except ValueError: 
                                condition_notes.append(f"Skipping rsID {rsid_profile} for PRS due to invalid odds_ratio for log: {odds_ratio}")
            else: 
                condition_notes.append(f"Variant {rsid_profile} from {condition_name} profile not found in patient's DNA data.")
        
        num_risk_variants_found = len(patient_risk_variants_details)
        num_protective_variants_found = len(patient_protective_variants_details)

        print(f"{Fore.BLUE}  Analysed {len(condition_variant_details)} profile variants for {condition_name}. Patient has {num_risk_variants_found} matching risk variants and {num_protective_variants_found} matching protective variants.{Style.RESET_ALL}")

        analysis_results[condition_name] = {
            'polygenic_risk_score': round(current_condition_score, 4),
            'risk_variants_found_count': num_risk_variants_found,
            'protective_variants_found_count': num_protective_variants_found,
            'risk_variants_details': patient_risk_variants_details,
            'protective_variants_details': patient_protective_variants_details,
            'summary_notes': condition_notes,
            'status': "Analysis Complete"
        }
        print(f"{Fore.CYAN}Finished assessment for {condition_name}. PRS: {current_condition_score:.4f}{Style.RESET_ALL}")
    
    print(f"\n{Fore.GREEN}╔═══════════════════════════════════════════╗{Style.RESET_ALL}")
    print(f"{Fore.GREEN}║ Patient Genetic Risk Assessment Complete  ║{Style.RESET_ALL}")
    print(f"{Fore.GREEN}╚═══════════════════════════════════════════╝{Style.RESET_ALL}")
    return analysis_results