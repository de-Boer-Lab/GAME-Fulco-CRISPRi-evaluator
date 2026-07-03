'''Handle Loading and Validating Evaluator Input/Request Data'''

import os
import json
from collections import Counter
import functools
import pandas as pd
import random
#set seed so the enhancer shuffling is reproducible
random.seed(42)
from Bio.Seq import Seq
from config import EVALUATOR_INPUT_PATH_GENE_SEQ, EVALUATOR_INPUT_PATH_Enhancer_coordinates

class DuplicateKeysError(ValueError):
    """Raised when duplicate keys are found in a JSON object."""
    pass

# Internal helper function to detect duplicates during JSON parsing
def _detect_duplicates(pairs, duplicate_keys_state):

    """
    Detects duplicate keys during JSON parsing and counts occurrences of each key.

    This function intercepts the key-value pairs provided by `json.loads` and ensures that
    duplicate keys are flagged. It constructs the dictionary normally but counts how often
    each key appears, recording any keys that occur more than once.

    Args:
        pairs (list of tuple): A list of key-value pairs at the current level of the JSON.
        duplicate_keys_state (dict): The dictionary to update with any duplicates found.

    Returns:
        result_dict: A dictionary created from the key-value pairs.
    """

    # Use a local Counter to count occurrences of keys at this level
    local_counts = Counter()
    result_dict = {}
    for key, value in pairs:
        # Increment the count for each key
        local_counts[key] += 1
        # If the key is a duplicate, record it in the duplicate_keys dictionary
        if local_counts[key] > 1:
            duplicate_keys_state[key] = local_counts[key]
        # Add the key-value pair to the resulting dictionary
        result_dict[key] = value
    return result_dict

def _process_results(data, duplicate_keys):
    """
    Checks the duplicate_keys dictionary and prints a report.

    Args:
        data (dict): The dictionary of parsed data. 
        duplicate_keys (dict): The dictionary of duplicates.

    Returns:
        data or None: The parsed data if no duplicates. None, if duplicates are found.
    """
    # Report duplicates if any were found
    if duplicate_keys:
        print("Duplicate keys found:")
        error_messages = [f"Key: '{key}', Count: {count}" for key, count in duplicate_keys.items()]
        raise DuplicateKeysError(f"Duplicate keys found:\n" + "\n".join(error_messages))
    else:
        print("No duplicates found.")
        return data # Return the parsed data if no duplicates.


# Function to check for duplicate keys in JSON object

def check_duplicates_from_string(json_string):

    """
    Parses a JSON string to detect and report any duplicate keys at the same level in the same object.
    This function ensures that no keys are silently overwritten in dictionaries.

    The function uses a helper to track the number of times each key appears during parsing,
    leveraging the `object_pairs_hook` parameter of `json.loads()` to intercept key-value pairs
    before they are processed into a dictionary. If duplicates are detected at any level, they
    are reported with their counts. Keys reused in separate objects within arrays (e.g. lists) 
    are not considered duplicates.

    Args:
        json_string (str): The JSON content as a string to parse and check for duplicates.

    Raises:
        json.JSONDecodeError: If the string is not valid JSON.
        DuplicateKeysError: If duplicate keys are found in the JSON structure.

    Returns:
        dict: The parsed data if no errors or duplicates are found.
    """

    # Initialize a dictionary to track duplicate keys and their counts
    duplicate_keys = {}
    
    # Create a 1-argument hook callable by "freezing" the duplicate_keys dict
    # as the second argument to the helper.
    hook = functools.partial(_detect_duplicates, duplicate_keys_state=duplicate_keys)

    # Parse the JSON string using the helper to track duplicates
    data = json.loads(json_string, object_pairs_hook=hook)
    
    return _process_results(data, duplicate_keys)
    
# Function for check for duplicate keys if input file is in JSON format

def check_duplicates_from_json(json_file_path):
    """
    Parses a JSON file to detect and report any duplicate keys at the same level in the same object.
    This function ensures that no keys are silently overwritten in dictionaries.

    The function uses a helper to track the number of times each key appears during parsing,
    leveraging the `object_pairs_hook` parameter of `json.load()` to intercept key-value pairs 
    before they are processed into a dictionary. If duplicates are detected at any level, they
    are reported with their counts and paths. Keys reused in separate objects within arrays 
    (e.g. lists) are not considered duplicates.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        json.JSONDecodeError: If the file content is not valid JSON.
        DuplicateKeysError: If duplicate keys are found in the JSON structure.

    Returns:
        dict: The parsed data if no errors or duplicates are found.
    """

    # Initialize a dictionary to track duplicate keys and their counts
    duplicate_keys = {}
    
    # Create a 1-argument hook callable by "freezing" the duplicate_keys dict
    # as the second argument to the helper.
    hook = functools.partial(_detect_duplicates, duplicate_keys_state=duplicate_keys)

    # Open and parse the JSON file, using the helper to track duplicates
    with open(json_file_path, 'r') as file:
        data = json.load(file, object_pairs_hook=hook)
        
    return _process_results(data, duplicate_keys)

def shuffle_enhancer(row, gene_sequence, strand):
    #seq = row["sequence"]  # full pulled window sequence
    seq = gene_sequence
    # enhancer coords relative to the window
    enh_start = row["start"] - row["sequence_start"]
    enh_end   = row["end"] - row["sequence_start"]

    #Gene coordinates relative to the window
    gene_start_relative = row["gene_start"] - row["sequence_start"]
    gene_end_relative   = row["gene_end"] - row["sequence_start"]
    # extract enhancer substring
    enhancer_seq = list(seq[enh_start:enh_end])

    # shuffle enhancer bases
    random.shuffle(enhancer_seq)
    shuffled_enhancer = "".join(enhancer_seq)
    #print(shuffled_enhancer)
    # rebuild full sequence with shuffled enhancer in place
    shuffled_full_seq = seq[:enh_start] + shuffled_enhancer + seq[enh_end:]
    #if the gene is on the negative strand take the reverse complement of the full sequence after shuffling the enhancer
    if strand == "-":
        shuffled_full_seq = str(Seq(shuffled_full_seq).reverse_complement())

    return pd.Series({
            "shuffled_enhancer": shuffled_full_seq,
            "gene_start_relative": gene_start_relative,
            "gene_end_relative": gene_end_relative
        })

def load_enhancer_gene_pair_data(gene):
    """
    Loads and validates the input file specified in `config.py`.
    
    This function checks for file existence and output directory, handles JSON or MsgPack formats,
    and runs duplicate key validation.

    Raises:
        FileNotFoundError: If the input file specified in `EVALUATOR_INPUT_PATH` 
                           does not exist.
        ValueError: If the file type is unsupported (not .json, .msgpack, or .mpk),
                    or if the data is malformed (e.g. invalid JSON/MsgPack),
                    or if duplicate keys are found via `evaluator_utils`.
    Returns:
        data_dict (dict): The validated data dictionary if loading and validation are successful.
    """
    # Validate evaluator input file exists
    if not os.path.exists(EVALUATOR_INPUT_PATH_GENE_SEQ):
        print(f"ERROR: Evaluator input file '{EVALUATOR_INPUT_PATH_GENE_SEQ}' not found.")
        raise FileNotFoundError(f"Evaluator input file not found: {EVALUATOR_INPUT_PATH_GENE_SEQ}")

    # Validate evaluator input file exists
    if not os.path.exists(EVALUATOR_INPUT_PATH_Enhancer_coordinates):
        print(f"ERROR: Evaluator input file '{EVALUATOR_INPUT_PATH_Enhancer_coordinates}' not found.")
        raise FileNotFoundError(f"Evaluator input file not found: {EVALUATOR_INPUT_PATH_Enhancer_coordinates}")

    #Code can be altered for specific needs
    try:
        gene_sequences = pd.read_parquet(EVALUATOR_INPUT_PATH_GENE_SEQ)
        gene_sequence_current = gene_sequences.loc[gene_sequences['Gene'] == gene, 'sequence'].iloc[0]
        gene_strand_current = gene_sequences.loc[gene_sequences['Gene'] == gene, 'strand'].iloc[0]
        
        #this will pull all the e-g pairs for the current gene
        eg_coordinates = pd.read_parquet(EVALUATOR_INPUT_PATH_Enhancer_coordinates)
        eg_coordinates_current = eg_coordinates[eg_coordinates['Gene'] == gene].copy()

        eg_coordinates_current[["shuffled_enhancer", "gene_start_relative", "gene_end_relative"]] = eg_coordinates_current.apply(shuffle_enhancer, axis=1, args = (gene_sequence_current,  gene_strand_current,))
        #If the gene is on the negative strand take the reverse complement of the reference gene sequence
        if gene_strand_current == "-":
            gene_sequence_current = str(Seq(gene_sequence_current).reverse_complement())
            
        prediction_tasks_str = f"""
        [
            {{
                "name": "{gene}",
                "type": "expression",
                "cell_type": "K562",
                "scale": "linear",
                "species": "homo_sapiens"
            }}
        ]
        """
        prediction_tasks = check_duplicates_from_string(prediction_tasks_str)

        sequence_dict = {gene + '_reference_seq': gene_sequence_current}
        eg_dict = dict(zip(gene + "_" + eg_coordinates_current["Element name"], eg_coordinates_current["shuffled_enhancer"]))

        sequence_dict.update(eg_dict)

        gene_start = int(eg_coordinates_current["gene_start_relative"].unique()[0])
        gene_end = int(eg_coordinates_current["gene_end_relative"].unique()[0])
        gene_ranges_current = [gene_start, gene_end]
        prediction_ranges_dict  = {key: gene_ranges_current for key in sequence_dict.keys()}
        
        # Assemble the request payload (metadata already duplicate-checked above).
        json_evaluator = {
            "readout": "point",
            "prediction_tasks": prediction_tasks,
            "sequences": sequence_dict,
            "prediction_ranges": prediction_ranges_dict
        }

        data_dict = check_duplicates_from_string(json.dumps(json_evaluator))

        print("Input data loaded and validated successfully.")
        return data_dict
    
    except (json.JSONDecodeError, 
        DuplicateKeysError) as e:
        # Raise a general ValueError that the main script's handler
        # will catch and report cleanly
        raise ValueError(f"Input data is invalid.\nDetails: {e}") from e
