import gseapy
import logging
import sys

# --- 0. Setup Logging ---
# Configure logging to provide clear feedback on the script's progress.
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def get_filtered_kegg_pathways(matched_gene_ids, library_name='KEGG_2021_Human', min_genes_per_pathway=3):
    """
    Fetches KEGG pathway data using gseapy and filters it based on a provided list of genes.

    Args:
        matched_gene_ids (list or set): A list or set of gene symbols (strings) that are
                                        present in both transcriptomics and proteomics datasets.
        library_name (str): The name of the gseapy library to use. Defaults to 'KEGG_2021_Human'.
                            To see available libraries, use `gseapy.get_library_name()`.
        min_genes_per_pathway (int): The minimum number of matched genes a pathway must contain
                                     to be included in the final output.

    Returns:
        dict: A dictionary where keys are pathway names and values are lists of gene symbols
              from `matched_gene_ids` that belong to that pathway. Returns an empty dict on error.
    """
    logger.info(f"--- Starting KEGG Pathway Fetch for library: {library_name} ---")

    # Ensure matched_gene_ids is a set for efficient lookup (O(1) average time complexity)
    matched_gene_set = set(matched_gene_ids)
    
    try:
        # Fetch the specified gene set library from gseapy's database
        # This returns a dictionary: {'pathway_name_1': [gene1, gene2,...], ...}
        logger.info(f"Downloading gene set library: {library_name}...")
        all_kegg_pathways = gseapy.get_library(name=library_name, organism='Human')
        logger.info(f"Successfully downloaded {len(all_kegg_pathways)} total pathways.")

    except Exception as e:
        logger.error(f"Failed to download or access the gseapy library '{library_name}'. Error: {e}")
        logger.error("Please check your internet connection and if the library name is correct.")
        logger.error("You can list available libraries with `gseapy.get_library_name()`.")
        return {}

    # --- Filtering Step ---
    logger.info("Filtering pathways based on the provided list of matched genes...")
    filtered_pathways = {}

    # Iterate through each pathway and its associated genes from the downloaded library
    for pathway_name, genes_in_pathway in all_kegg_pathways.items():
        
        # Find the intersection between genes in the current KEGG pathway and our matched genes
        # Gene symbols from gseapy are typically uppercase, so we convert our list for a robust match.
        # However, for this example we assume matching cases. For real data, normalization is key.
        genes_in_pathway_set = set(genes_in_pathway)
        common_genes = list(matched_gene_set.intersection(genes_in_pathway_set))
        
        # Check if the number of common genes meets our minimum threshold
        if len(common_genes) >= min_genes_per_pathway:
            # If it does, add the pathway and its relevant genes to our results
            # Sort the genes for consistent ordering
            filtered_pathways[pathway_name] = sorted(common_genes)
    
    logger.info(f"Filtering complete. Found {len(filtered_pathways)} pathways containing at least {min_genes_per_pathway} matched genes.")
    
    return filtered_pathways

# --- Example Usage ---
if __name__ == '__main__':
    logger.info("--- Running Example Usage ---")

    # 1. To see all available libraries in gseapy (especially useful for finding the latest KEGG version)
    try:
        logger.info("Available gseapy libraries:")
        # print(gseapy.get_library_name(organism='Human')) # Uncomment to see the full list
    except Exception as e:
        logger.warning(f"Could not retrieve library names. Might be an internet issue. {e}")


    # 2. Define a dummy list of "matched gene IDs" for demonstration.
    # In your full script, this would be the list of genes common to your omics data.
    # These are examples of genes involved in the Glycolysis pathway.
    example_matched_ids = [
        'HK1', 'HK2', 'GCK', 'HK3', 'GPI', 'PFKL', 'PFKM', 'PFKP', 'ALDOA', 
        'ALDOB', 'ALDOC', 'TPI1', 'GAPDH', 'PGK1', 'PGK2', 'PGAM1', 'PGAM2', 
        'ENO1', 'ENO2', 'ENO3', 'PKLR', 'PKM', 'LDHA', 'LDHB', 'ACSS2', 'FOXO1', 'PCK1'
        # Add a few random genes not in the pathway to simulate a real scenario
        , 'EGFR', 'TP53', 'BRCA1'
    ]
    logger.info(f"\nUsing an example list of {len(example_matched_ids)} matched gene IDs for demonstration.")

    # 3. Call the function to get the filtered pathways
    # We use KEGG_2021_Human, which is a common, stable choice.
    pathways = get_filtered_kegg_pathways(
        matched_gene_ids=example_matched_ids,
        library_name='KEGG_2021_Human',
        min_genes_per_pathway=5 # We want pathways with at least 5 of our genes
    )

    # 4. Print the results
    if pathways:
        logger.info("\n--- Results: Filtered Pathways ---")
        # Print the first 5 pathways found for brevity
        for i, (pathway_name, genes) in enumerate(pathways.items()):
            if i >= 5:
                logger.info(f"...and {len(pathways) - 5} more pathways.")
                break
            logger.info(f"\nPathway: {pathway_name}")
            logger.info(f"  Matched Genes ({len(genes)}): {genes}")
            
        # Example of accessing a specific pathway
        glycolysis_pathway = "Glycolysis / Gluconeogenesis hsa00010"
        if glycolysis_pathway in pathways:
            logger.info(f"\nSuccessfully found the '{glycolysis_pathway}' pathway.")
            logger.info(f"It contains {len(pathways[glycolysis_pathway])} of our matched genes.")

    else:
        logger.warning("No pathways matching the criteria were found.")
    logger.info("--- Example Usage Complete ---")