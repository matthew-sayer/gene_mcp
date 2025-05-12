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
        df = pd.read_csv(file_path, sep=',', comment='#', low_memory=False, keep_default_na=False, na_values=['', 'NA', 'N/A', 'NaN', '0', '-'])
        
        required_cols = ['rsid', 'allele1', 'allele2'] 
        
        df.columns = df.columns.str.lower().str.strip()
        original_columns = df.columns.tolist() # For error reporting

        if 'rsid' not in df.columns:
            if 'rs id' in df.columns:
                df.rename(columns={'rs id': 'rsid'}, inplace=True)
            elif 'snp id' in df.columns: 
                df.rename(columns={'snp id': 'rsid'}, inplace=True)
            elif 'marker' in df.columns: # Another common one
                 df.rename(columns={'marker': 'rsid'}, inplace=True)
        
        if not all(col in df.columns for col in required_cols):
            missing_cols = [col for col in required_cols if col not in df.columns]
            raise ValueError(f"Patient DNA file is missing required columns (expected: {required_cols}). Found columns: {original_columns}. Missing: {missing_cols}")

        print(f"{Fore.BLUE}Processing {len(df)} rows from patient DNA file...{Style.RESET_ALL}")
        skipped_rows_count = 0
        error_rows_count = 0

        for index, row in df.iterrows():
            try:
                rsid_val = row.get('rsid')
                allele1_val = row.get('allele1')
                allele2_val = row.get('allele2')

                if pd.isna(rsid_val):
                    # Only skip if the RSID itself is missing - we need that for identification in all cases, or it becomes unreliable.
                    if skipped_rows_count < 5:
                        print(f"{Fore.YELLOW}Warning: Skipping row {index+2} due to missing rsID. Data: rsid='{rsid_val}'{Style.RESET_ALL}")
                    skipped_rows_count += 1
                    continue
                
                # Process the RSID even with missing alleles
                rsid = str(rsid_val).strip().lower()
                if not rsid or rsid == 'nan':
                    if skipped_rows_count < 5:
                        print(f"{Fore.YELLOW}Warning: Skipping row {index+2} due to invalid rsID after processing: '{rsid}'{Style.RESET_ALL}")
                    skipped_rows_count += 1
                    continue

                # Handle alleles, using placeholders for missing data
                allele1 = "?" if pd.isna(allele1_val) else str(allele1_val).strip().upper()
                allele2 = "?" if pd.isna(allele2_val) else str(allele2_val).strip().upper()

                if pd.isna(allele1_val) or pd.isna(allele2_val):
                    print(f"{Fore.BLUE}Info: Processing row {index+2} with incomplete allele data. RSID={rsid}, alleles={allele1}/{allele2}{Style.RESET_ALL}")

                valid_alleles = "AGCTDI0-?" # valid allele values
                if not all(a in valid_alleles for a in allele1) or not all(a in valid_alleles for a in allele2):
                    if error_rows_count < 5:
                        print(f"{Fore.YELLOW}Warning: Skipping row {index+2} due to invalid characters in alleles for rsID {rsid}. Alleles: ({allele1}, {allele2}){Style.RESET_ALL}")
                    error_rows_count += 1
                    continue

                patient_genotypes[rsid] = (allele1, allele2)
            except Exception as e_row:
                if error_rows_count < 5: 
                    print(f"{Fore.RED}Error processing row {index+2} in patient DNA file: {e_row}. Row data: {row.to_dict()}{Style.RESET_ALL}")
                error_rows_count += 1
                continue 
        
        if skipped_rows_count > 0:
            print(f"{Fore.YELLOW}Total skipped rows due to missing data: {skipped_rows_count}{Style.RESET_ALL}")
        if error_rows_count > 0:
            print(f"{Fore.RED}Total rows with errors during processing: {error_rows_count}{Style.RESET_ALL}")
            
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

@mcp.tool()
async def assess_patient_genetic_risk_profile(
    patient_dna_file_path: str,
    condition_risk_profiles: Dict[str, Dict[str, VariantInfo]] 
) -> Dict[str, Any]:
    """
    Assesses patient's genetic data against pre-compiled risk profiles for specified conditions.
    This tool will meticulously check every RSID provided in the condition_risk_profiles.

    NEVER use FAKE RSIDs, like rs-1234567, or similar unbelievable sequences that you might hallucinate. If you are not sure, DO THE DEEP RESEARCH BEFOREHAND USING OUR MCP TOOLS.

    Args:
        patient_dna_file_path: Path to the patient's DNA data file (CSV format).
                               Expected columns: rsid, allele1, allele2.
        condition_risk_profiles: A dictionary where keys are condition names (e.g., "autism")
                                 and values are dictionaries of variant details (Dict[str, VariantInfo])
                                 for that condition.
                                 Example: {"autism": {"rs123": {"risk_allele": "A", "odds_ratio": 1.5, ...}}}

    Returns:
        A dictionary where keys are condition names. Each value is a dictionary containing:
        - 'polygenic_risk_score': Calculated PRS for the condition.
        - 'risk_variants_found_count': Number of risk variants present in the patient.
        - 'protective_variants_found_count': Number of protective variants present in the patient.
        - 'risk_variants_details': List of specific risk variants present in the patient.
        - 'protective_variants_details': List of specific protective variants present in the patient.
        - 'found_variants_details': List of ALL variants found in patient's DNA, regardless of risk.
        - 'missing_variants_details': List of variants from the profile NOT FOUND in patient's DNA.
        - 'missing_variants_count': Count of variants from profile NOT FOUND in patient's DNA.
        - 'total_variants_in_profile': Total variants defined in the risk profile for the condition.
        - 'summary_notes': General notes on the analysis for this condition.
    """
    print(f"{Fore.GREEN}╔═══════════════════════════════════════════════════════════════╗{Style.RESET_ALL}")
    print(f"{Fore.GREEN}║ Starting Patient Genetic Risk Assessment for conditions: {Fore.YELLOW}{', '.join(condition_risk_profiles.keys())}{Fore.GREEN} ║{Style.RESET_ALL}")
    print(f"{Fore.GREEN}╚═══════════════════════════════════════════════════════════════╝{Style.RESET_ALL}")

    analysis_results: Dict[str, Any] = {}

    try:
        patient_genotypes = await load_patient_dna_from_file(patient_dna_file_path)
        if not patient_genotypes:
            error_message = "Failed to load patient DNA data or no valid genotypes found."
            print(f"{Fore.RED}{error_message}{Style.RESET_ALL}")
            for condition_name in condition_risk_profiles.keys():
                analysis_results[condition_name] = {"error": error_message, "status": "Patient DNA Load Failed"}
            return analysis_results
    except Exception as e:
        error_message = f"Failed to load patient DNA: {str(e)}"
        print(f"{Fore.RED}{error_message}{Style.RESET_ALL}")
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
                'missing_variants_details': [],
                'missing_variants_count': 0,
                'total_variants_in_profile': 0,
                'summary_notes': [f"No genetic risk factors provided in the profile for {condition_name}."],
                'status': "No Profile Variants"
            }
            continue

        current_condition_score = 0.0
        patient_risk_variants_details = []
        patient_protective_variants_details = []
        found_variants_details = []
        missing_variants_details = [] 
        condition_notes = []
        
        for rsid_profile, variant_info in condition_variant_details.items():
            rsid_lookup = rsid_profile.lower().strip() 

            if rsid_lookup in patient_genotypes:
                patient_allele1, patient_allele2 = patient_genotypes[rsid_lookup]
                
                found_variant_entry = {
                    'rsid': rsid_profile,
                    'patient_genotype': f"{patient_allele1}/{patient_allele2}",
                    'expected_risk_allele': variant_info.get('risk_allele'),
                    'odds_ratio': variant_info.get('odds_ratio'),
                    'reported_gene': variant_info.get('reported_gene'),
                    'p_value': variant_info.get('pvalue')
                }
                found_variants_details.append(found_variant_entry)
                
                effect_allele = None
                allele_type = None # 'risk' or 'protective'
                patient_has_effect_allele = False  # Initialize this variable

                if 'risk_allele' in variant_info and variant_info['risk_allele']:
                    effect_allele = str(variant_info['risk_allele']).strip().upper()
                    allele_type = 'risk'
                    
                    # Check for both direct match and complement allele match
                    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                    comp_effect_allele = complement_map.get(effect_allele)
                    
                    if comp_effect_allele:
                        # Check for both direct and complement match
                        patient_has_effect_allele = (patient_allele1 == effect_allele or patient_allele2 == effect_allele or
                                                   patient_allele1 == comp_effect_allele or patient_allele2 == comp_effect_allele)
                    else:
                        # For non-standard alleles, just do direct match
                        patient_has_effect_allele = (patient_allele1 == effect_allele or patient_allele2 == effect_allele)
                    
                elif 'protective_allele' in variant_info and variant_info['protective_allele']:
                    effect_allele = str(variant_info['protective_allele']).strip().upper()
                    allele_type = 'protective'
                    
                    # Apply the same strand-aware matching to protective alleles
                    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                    comp_effect_allele = complement_map.get(effect_allele)
                    
                    if comp_effect_allele:
                        patient_has_effect_allele = (patient_allele1 == effect_allele or patient_allele2 == effect_allele or
                                                   patient_allele1 == comp_effect_allele or patient_allele2 == comp_effect_allele)
                    else:
                        patient_has_effect_allele = (patient_allele1 == effect_allele or patient_allele2 == effect_allele)
                
                odds_ratio = variant_info.get('odds_ratio')

                if effect_allele and odds_ratio is not None and odds_ratio > 0 and patient_has_effect_allele:
                    try:
                        score_contribution = math.log(odds_ratio) 
                        current_condition_score += score_contribution

                        variant_detail_entry = {
                            'rsid': rsid_profile,
                            'patient_genotype': f"{patient_allele1}/{patient_allele2}",
                            'effect_allele': effect_allele,
                            'allele_type_in_profile': allele_type,
                            'odds_ratio': odds_ratio,
                            'reported_gene': variant_info.get('reported_gene'),
                            'p_value': variant_info.get('pvalue'),
                            'profile_mention_count': variant_info.get('mention_count', 1) 
                        }
                        
                        # Categorise based on odds_ratio primarily. Odds ratio is a measure of association between a variant and a trait.
                        if odds_ratio > 1: # Generally considered risk-increasing
                            patient_risk_variants_details.append(variant_detail_entry)
                        elif odds_ratio < 1: # Generally considered protective
                            patient_protective_variants_details.append(variant_detail_entry)
                        
                    except ValueError: 
                        condition_notes.append(f"Skipping rsID {rsid_profile} for PRS due to invalid odds_ratio for log: {odds_ratio}")
                
            else: 
                missing_variant_entry = {
                    'rsid': rsid_profile,
                    'risk_allele_in_profile': variant_info.get('risk_allele'),
                    'protective_allele_in_profile': variant_info.get('protective_allele'),
                    'odds_ratio_in_profile': variant_info.get('odds_ratio'),
                    'reported_gene_in_profile': variant_info.get('reported_gene'),
                    'p_value_in_profile': variant_info.get('pvalue')
                }
                missing_variants_details.append(missing_variant_entry)

        num_risk_variants_found = len(patient_risk_variants_details)
        num_protective_variants_found = len(patient_protective_variants_details)
        num_found_variants = len(found_variants_details)  # NEW
        num_missing_variants = len(missing_variants_details)
        total_variants_in_profile = len(condition_variant_details)

        print(f"{Fore.BLUE}  Condition '{condition_name}': Analysed {total_variants_in_profile} profile variants.")
        print(f"{Fore.BLUE}    Found in DNA: {num_found_variants} variants regardless of risk status.")
        print(f"{Fore.BLUE}    Risk assessment: {num_risk_variants_found} risk variants, {num_protective_variants_found} protective variants.")
        
        if num_missing_variants > 0:
            print(f"{Fore.YELLOW}    Not Found in Patient DNA: {num_missing_variants} variants from profile.{Style.RESET_ALL}")
        else:
            print(f"{Fore.GREEN}    All {total_variants_in_profile} profile variants were found in patient's DNA data.{Style.RESET_ALL}")


        analysis_results[condition_name] = {
            'polygenic_risk_score': round(current_condition_score, 4),
            'risk_variants_found_count': num_risk_variants_found,
            'protective_variants_found_count': num_protective_variants_found,
            'risk_variants_details': patient_risk_variants_details,
            'protective_variants_details': patient_protective_variants_details,
            'found_variants_details': found_variants_details,  # NEW: All variants found regardless of risk status
            'found_variants_count': num_found_variants,  # NEW
            'missing_variants_details': missing_variants_details,
            'missing_variants_count': num_missing_variants,
            'total_variants_in_profile': total_variants_in_profile,
            'summary_notes': condition_notes,
            'status': "Analysis Complete"
        }
        
        if num_missing_variants > 0:
             analysis_results[condition_name]['summary_notes'].append(
                f"WARNING: {num_missing_variants}/{total_variants_in_profile} variants for {condition_name} were NOT FOUND in the patient's DNA data. Results may be incomplete."
            )
            
        print(f"{Fore.CYAN}Finished assessment for {condition_name}. PRS: {current_condition_score:.4f}{Style.RESET_ALL}")
    
    print(f"\n{Fore.GREEN}╔═══════════════════════════════════════════╗{Style.RESET_ALL}")
    print(f"{Fore.GREEN}║ Patient Genetic Risk Assessment Complete  ║{Style.RESET_ALL}")
    print(f"{Fore.GREEN}╚═══════════════════════════════════════════╝{Style.RESET_ALL}")
    return analysis_results