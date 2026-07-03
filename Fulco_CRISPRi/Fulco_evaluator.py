'''RESTful Test Evaluator Utilizing Flask'''

import os
import sys
import json
from requests.exceptions import RequestException, HTTPError
import glob
from data_loader import load_enhancer_gene_pair_data
from evaluator_content_handler import *
import evaluator_metrics_calculator
import pandas as pd
from config import EVALUATOR_NAME, EVALUATOR_INPUT_PATH_GENE_SEQ
from datetime import datetime, timezone

def run_evaluator(predictor_ip, predictor_port, output_dir):
    """
    Preprocesses the data, sends request, receives response,
    saves the response, and returns file path.
    """
    
    # Validate output directory exists; create if it does not
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory '{output_dir}' did not exist. Created it successfully!")
        
    # CHANGE (v1): collect per-gene Pearson r values in memory so the final average does not
    # depend on globbing/re-reading CSVs (the old glob double-counted the final summary on
    # re-runs and broke when the 'Value' column held the string "None"/"NaN")
    per_gene_r_values = []
    last_predictor_name = "UnknownPredictor"
        
    #Start a big loop here for all 30 Genes tested
    #This Evaluator will /post 30 requests due to large size
    gene_list = pd.read_parquet(EVALUATOR_INPUT_PATH_GENE_SEQ, columns=["Gene"])
    for gene_idx in range(0, gene_list.shape[0]):
        gene_current = gene_list.iloc[gene_idx]["Gene"]
        print(f"Sending sequence data for gene: {gene_current}")
        # Load and validate data, returns a JSON dictionary
        # This function will be Evaluator specific
        data_dict = load_enhancer_gene_pair_data(gene_current)
        
        #Total number of sequences being sent to the Evaluator
        #This is used to confirm the correct number of predictions will be returned from the Predictor
        total_sequences = len(data_dict["sequences"])

        # Communicate with Predictor, send request, receive predictions
        predictor_url = f"http://{predictor_ip}:{predictor_port}"
        response_payload = None
        is_success_response = False # Flag to track if we got a 200 OK
        
        try:
            #Decide the request format and response format
            req_fmt, resp_fmt = negotiate_formats(predictor_url)
            
            # This call will return a 200 OK response OR raise HTTPError (4xx/5xx)
            response = get_predictions(predictor_url, data_dict, req_fmt, resp_fmt)
            
            is_success_response = True
            print("Predictor returned 200 OK.")
            response_payload = deserialize_response(response, resp_fmt)
            
        except HTTPError as http_err:
            # We got a 4xx/5xx. Treat the payload as an error report.
            print(f"Predictor returned HTTP {http_err.response.status_code}. Processing error payload...")
            is_success_response = False
            try:
                # Deserialize the error body (always JSON)
                response_payload = deserialize_response(http_err.response, "application/json") 
            except ValueError as decode_err:
                # Handle cases where the error response itself is malformed
                print(f"Could not decode the error response body: {decode_err}", file=sys.stderr)
                # Create a fallback error payload
                response_payload = {
                    "predictor_name": "UnknownPredictor_ErrorResponse",
                    "error": [{"server_error": f"Failed to decode error response (Status {http_err.response.status_code}). Body: {http_err.response.text[:500]}..."}]
                }
                
        if response_payload is None:
            print("FATAL: No response payload received or processed.", file=sys.stderr)
            # Create a fallback payload indicating a severe issue
            response_payload = {
                "predictor_name": "UnknownPredictor_NoResponse",
                "error": [{"evaluator_error": "No response payload could be processed after request."}]
            }
        
        # Save predictions
        predictor_name = response_payload.get("predictor_name", "UnknownPredictor").replace(" ", "_")
        last_predictor_name = predictor_name  # CHANGE (v1): track for the final summary row
        # CHANGE (v1): build-aware output filename (predictor name already carries its build stamp)
        output_filename = f"{EVALUATOR_NAME}_{gene_current}_predictions_from_{predictor_name}.json"
        saved_predictions_path = os.path.join(output_dir, output_filename)
            
        # Check sequence counts before saving
        for task_idx, task in enumerate(response_payload.get("prediction_tasks", []), start=1):
            preds = task.get("predictions", {})
            if len(preds) != total_sequences:
                print(f"Warning: Task {task_idx} ('{task.get('name')}') has {len(preds)} predictions, but {total_sequences} sequences were sent to the Predictor.")
        
        try:
            with open(saved_predictions_path, 'w', encoding='utf-8') as f:
                json.dump(response_payload, f, ensure_ascii=False, indent=4)
            print(f"Raw predictions saved to {saved_predictions_path}")
        except IOError as e:
            print(f"FATAL: Could not save predictions to {saved_predictions_path}. {e}", file=sys.stderr)
            return

    # Calculate and save final metrics
        if is_success_response:
            all_lengths_match = True
            #Loop through and check all the prediction tasks
            for task_idx, task in enumerate(response_payload.get("prediction_tasks", []), start=1):
                # CHANGE (v1): predictions is a dict per spec -> default {} not [].
                preds = task.get("predictions", {})
                # CHANGE (v1): if the task reports an error, skip the length check entirely
                # (the original set the flag True but fell through to the length check, which
                # then overwrote it to False).
                if "error" in preds:
                    print(f"Task {task_idx} ('{task.get('name')}') returned an error -- skipping length check.")
                    continue
                #Otherwise length of predictions needs to == the # of sequences
                if len(preds) != total_sequences:
                    print(f"Warning: Task {task_idx} ('{task.get('name')}') has {len(preds)} predictions, but {total_sequences} sequences were sent to the Predictor.")
                    all_lengths_match = False
            if all_lengths_match:
                # CHANGE (v1): capture the returned r so we can average across genes in memory.
                gene_r = evaluator_metrics_calculator.calculate_and_save_metrics(gene_current, saved_predictions_path, output_dir)
                if gene_r is not None:
                    per_gene_r_values.append(gene_r)
            else:
                print("Skipping metric calculation because not all sequences got predictions.")
        else:
            print("Skipping metrics calculation because the Predictor did not return a 200 OK status.")

    #Finally concatenate and average the pearson r for each gene into 1 final pearson r
    print("Averaging per-gene Pearson r into one final value")
    if per_gene_r_values:
        final_pearsonr = float(pd.Series(per_gene_r_values, dtype="float64").mean())
        final_value_str = str(final_pearsonr)
    else:
        print("No per-gene correlations were computed; final average is NaN.")
        final_value_str = "NaN"
   
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")
    evaluation_output = pd.DataFrame([{
            'evaluator_name': EVALUATOR_NAME,
            'description': "Average Pearson Correlation across all genes for CRISPRi",
            'predictor_name': last_predictor_name,
            'time_stamp': timestamp,
            'metric': 'pearson_r',
            'value': final_value_str,
            'prediction_task(s)_data': "See per-gene rows in this same file"
        }])
    summary_filepath = os.path.join(output_dir, f"evaluation_summary_{EVALUATOR_NAME}.csv")
    evaluator_metrics_calculator._save_df_to_csv(evaluation_output, summary_filepath)
    

if __name__ == '__main__':
    # Check that mandatory arguments are passed to the script
    if len(sys.argv) != 4:
        print(f"Invalid arguments! Arguments must have: <container image/python script> <predictor_ip_address> <predictor_port> <mounted_output_directory>")
        sys.exit(1)
    
    # Call Evaluator here
    predictor_ip = sys.argv[1]
    predictor_port = int(sys.argv[2])
    output_dir_arg = sys.argv[3]
    
    try:
        run_evaluator(predictor_ip, int(predictor_port), output_dir_arg)
        print("Evaluation complete.")
        sys.exit(0)
        
    except (FileNotFoundError, ValueError) as e:
        print(f"FATAL ERROR (Data): {e}", file=sys.stderr)
        sys.exit(1)
    except RequestException as e:
        print(f"FATAL ERROR (Network): Could not connect to predictor at http://{predictor_ip}:{predictor_port}.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected fatal error occurred: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)