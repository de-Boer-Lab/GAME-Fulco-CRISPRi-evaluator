'''Calculate and save the final evaluation metrics.'''

# NOTE: Every evaluator will do this slightly differently depending on how the data is presented

import os
import sys
import json
import pandas as pd
import numpy as np
import itertools
from datetime import datetime, timezone
from scipy.stats import pearsonr
from config import EVALUATOR_NAME, EVALUATOR_INPUT_PATH_GENE_SEQ, EVALUATOR_INPUT_PATH_Enhancer_coordinates

def _save_df_to_csv(df, filepath):
    """
    Appends a DataFrame to a CSV file, adding a header if the file is new.
    """
    if df.empty:
        print(f"No metrics to save for {os.path.basename(filepath)}. Skipping.")
        return
    
    try:
        file_exists = os.path.isfile(filepath)
        df.to_csv(filepath, mode='a', sep='\t', header=(not file_exists), index=False)
        print(f"DEBUG: Metrics file '{filepath}' exists: {file_exists}")
        if file_exists:
            print(f"Appended metrics to {filepath}")
        else:
            print(f"Created new metrics file {filepath}")
    except IOError as e:
        print(f"\nError: Could not save metrics to {filepath}. {e}", file=sys.stderr)

def calculate_pearson_r(predictions_content, eg_measured_values, scale_prediction_actual, gene):
    """
    1- Calculates the Fraction change in gene expr.
    2- Calculated the pearsonR between measured and predicted fractional change

    Args:
        predictions_json_path (str): Path to JSON file with predictions
    Returns:
        float: The Pearson correlation coefficient (r), or None if calculation isn't possible.
    """
   
    predictions_dict = predictions_content['prediction_tasks'][0]['predictions']

    if "error" in predictions_dict:
        print("No predictions were returned for this task -> Skipping evaluation calculation")
        return None
    # Create DataFrame from Predictions
    predictions_df = pd.DataFrame(list(predictions_dict.items()), columns=['Element name', 'Predicted_Value'])
    predictions_df['Predicted_Value'] = predictions_df['Predicted_Value'].apply(
        lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x
        )
    print(predictions_df)
    #check here is there is NA is any of the prediction values
    na_rows = predictions_df[predictions_df['Predicted_Value'].isna()]
    if not na_rows.empty:
        print("NA values were found in the predictions, skipping evaluation")
        print(na_rows)
        return None

    ref_row = predictions_df[predictions_df["Element name"].str.contains("_reference_seq", case=False)]
    #First row of the predictions is the Gene reference sequence predictions
    ref_id = ref_row["Element name"].iloc[0]
    ref_pred = ref_row["Predicted_Value"].iloc[0]

    #Add column that has the reference predictions in each row
    predictions_df[ref_id] = ref_pred
    #Remove the row that has the reference sequence prediction twice
    predictions_df = predictions_df[~predictions_df["Element name"].str.contains("_reference_seq", case=False)].reset_index(drop=True)
    if scale_prediction_actual == 'log':
        #Needs to be in linear scale for this Evaluator
        print("Error: Received log scaled predictions from Predictor, this Evaluator can only handle linear")
        return None
    if scale_prediction_actual == 'linear':
        fraction_change_predicted = (predictions_df['Predicted_Value'] - predictions_df[ref_id])/predictions_df[ref_id]

    predictions_df['fraction_change_predicted'] = fraction_change_predicted
    print(predictions_df)

    eg_measured_values_current_gene = eg_measured_values[eg_measured_values['Gene'] == gene]
    
    eg_measured_values_current_gene['Element name'] = gene + "_" + eg_measured_values_current_gene['Element name']
    #Covert , to .
    eg_measured_values_current_gene['Fraction change in gene expr'] = (
        eg_measured_values_current_gene['Fraction change in gene expr']
        .str.replace(',', '.', regex=False)
        .astype(float)
        )

    merged_data = predictions_df.merge(eg_measured_values_current_gene, how='left', on='Element name')
    print(merged_data['Fraction change in gene expr'])
    print(merged_data['fraction_change_predicted'])
    try:
        r, _ = pearsonr(merged_data['fraction_change_predicted'], merged_data['Fraction change in gene expr'])
        print(f"Calculated Pearson r: {r}") 
    
        return r
    except:
        print("An error occured during the correlation calculation")

def calculate_and_save_metrics(gene, saved_predictions_path_current, output_dir):

    print(f"Using measured data from: {EVALUATOR_INPUT_PATH_Enhancer_coordinates}")

    # Define output paths
    correlation_summary_filepath = f"correlation_summary_{EVALUATOR_NAME}_{gene}.csv"
    try:
        with open(saved_predictions_path_current, 'r') as f:
            predictions_file_content = json.load(f)
    except Exception as e:
        print(f"FATAL: Could not load predictions from {saved_predictions_path_current}. {e}", file=sys.stderr)
        return

    # Now load measured values
    try:
        eg_coordinates = pd.read_parquet(EVALUATOR_INPUT_PATH_Enhancer_coordinates)
    except Exception as e:
        print(f"FATAL: Could not load measured data from {EVALUATOR_INPUT_PATH_Enhancer_coordinates}. {e}", file=sys.stderr)
        return

    try:
        # ADDITION: Construct file name after receiving predictor_name
        predictor_name_received = predictions_file_content.get("predictor_name")
        scale_prediction_actual = predictions_file_content['prediction_tasks'][0].get("scale_prediction_actual")

        predictor_name = predictor_name_received.replace(" ", "_").replace("/", "_")

        # Get UTC timestamp for predictor_nam
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")
        pearson_r = calculate_pearson_r(predictions_file_content, eg_coordinates, scale_prediction_actual, gene)
        description = f"Fulco CRISPRi ({gene})"

        prediction_task_data_onlyinfo = [{k: v for k, v in predictions_file_content["prediction_tasks"][0].items() if k != "predictions"}]

        evaluation_output = pd.DataFrame([{
            'Evaluator': EVALUATOR_NAME,
            'Description': description,
            'Predictor_name': predictor_name,
            'Time_stamp': timestamp,
            'Metric': 'pearson_r',
            'Value': str(pearson_r),
            'Prediction_task(s)_data': prediction_task_data_onlyinfo
        }])
        #print(evaluation_output)
        evaluation_output.to_csv(output_dir + '/' + correlation_summary_filepath , sep = "\t")
    except Exception as e:
        print(f"An unexpected error occurred during evaluation calculations: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()


