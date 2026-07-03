'''Calculate and save the final evaluation metrics.'''

# NOTE: Every evaluator will do this slightly differently depending on how the data is presented

import os
import sys
import json
import pandas as pd
import numpy as np
from datetime import datetime, timezone
from scipy.stats import pearsonr
from config import EVALUATOR_NAME, EVALUATOR_INPUT_PATH_Enhancer_coordinates

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
        eg_measured_values (pd.DataFrame): DataFrame containing measured values with columns ['Gene', 'Element name', 'Fraction change in gene expr']
        scale_prediction_actual (str): The scale of the predictions received from the predictor, either 'log' or 'linear'
        gene (str): The gene for which the evaluation is being calculated
    Returns:
        float: The Pearson correlation coefficient (r), or None if calculation isn't possible
        0.0: If predictions (or measured values) have zero variable -> "ran but useless"
        None: If the task is disqualified due to errors, NA values, or other issues.
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
    else:
        print(f"Warning: Unrecognized or missing scale '{scale_prediction_actual}'. Defaulting to linear math.")
        fraction_change_predicted = (predictions_df['Predicted_Value'] - predictions_df[ref_id])/predictions_df[ref_id]

    predictions_df['fraction_change_predicted'] = fraction_change_predicted
    print(predictions_df)

    eg_measured_values_current_gene = eg_measured_values[eg_measured_values['Gene'] == gene].copy()
    
    eg_measured_values_current_gene['Element name'] = gene + "_" + eg_measured_values_current_gene['Element name']
    # Convert , to .
    eg_measured_values_current_gene['Fraction change in gene expr'] = (
        eg_measured_values_current_gene['Fraction change in gene expr']
        .str.replace(',', '.', regex=False)
        .astype(float)
        )

    # CHANGE (v1): Changed to inner join to ensure only rows with both predictions and measurements are kept
    # Any missing predictions for sequences sent to Predictors will be flagged by the main Evaluator script
    # This will drop any rows that have predictions but no measurements, which is necessary for a valid correlation calculation.
    merged_data = predictions_df.merge(eg_measured_values_current_gene, how='inner', on='Element name')
    
    # CHANGE (v1): coerce to numeric and drop any rows that can't be correlated, with a count.
    merged_data['fraction_change_predicted'] = pd.to_numeric(merged_data['fraction_change_predicted'], errors='coerce')
    merged_data['Fraction change in gene expr'] = pd.to_numeric(merged_data['Fraction change in gene expr'], errors='coerce')
    before = len(merged_data)
    merged_data = merged_data.dropna(subset=['fraction_change_predicted', 'Fraction change in gene expr'])
    dropped = before - len(merged_data)
    if dropped:
        print(f"Dropped {dropped} row(s) with non-numeric/missing values before correlation.")
    
    print(merged_data['Fraction change in gene expr'])
    print(merged_data['fraction_change_predicted'])
    
    # CHANGE (v1): zero-variance handling
    # pearsonr is undefined for a constant input; a flat-line prediction "ran but is useless" -> 0.0, not NaN.
    std_predicted = merged_data['fraction_change_predicted'].std()
    std_measured = merged_data['Fraction change in gene expr'].std()
    if std_predicted == 0 or std_measured == 0:
        print("Zero variance in predicted or measured values -> assigning pearson_r = 0.0")
        return 0.0

    try:
        r, _ = pearsonr(merged_data['fraction_change_predicted'], merged_data['Fraction change in gene expr'])
        # CHANGE (v1) a NaN result maps to 0.0 ("ran but useless").
        pearson_r = 0.0 if np.isnan(r) else float(r)
        print(f"Calculated Pearson r: {pearson_r}")
        return pearson_r
    except Exception as e:
        # CHANGE (v1): no bare except; report and disqualify.
        print(f"An error occurred during the correlation calculation: {e}", file=sys.stderr)
        return None

def calculate_and_save_metrics(gene, saved_predictions_path_current, output_dir):
    """
    Computes the per-gene Pearson r, appends one row to the single evaluation summary
    CSV, and RETURNS the r value (float, 0.0, or None) so the caller can compute the
    cross-gene average without re-reading any CSVs.
    """

    print(f"Using measured data from: {EVALUATOR_INPUT_PATH_Enhancer_coordinates}")

    # Define output paths
    summary_filepath = os.path.join(output_dir, f"evaluation_summary_{EVALUATOR_NAME}.csv")
    
    try:
        with open(saved_predictions_path_current, 'r') as f:
            predictions_file_content = json.load(f)
    except Exception as e:
        print(f"FATAL: Could not load predictions from {saved_predictions_path_current}. {e}", file=sys.stderr)
        return None

    # Now load measured values
    try:
        eg_coordinates = pd.read_parquet(EVALUATOR_INPUT_PATH_Enhancer_coordinates)
    except Exception as e:
        print(f"FATAL: Could not load measured data from {EVALUATOR_INPUT_PATH_Enhancer_coordinates}. {e}", file=sys.stderr)
        return None

    try:
        # ADDITION: Construct file name after receiving predictor_name
        predictor_name_received = predictions_file_content.get("predictor_name")
        scale_prediction_actual = predictions_file_content['prediction_tasks'][0].get("scale_prediction_actual")

        predictor_name = predictor_name_received.replace(" ", "_").replace("/", "_")

        # Get UTC timestamp
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")
        pearson_r = calculate_pearson_r(predictions_file_content, eg_coordinates, scale_prediction_actual, gene)
        description = f"Fulco CRISPRi ({gene})"

        prediction_task_data_onlyinfo = [{k: v for k, v in predictions_file_content["prediction_tasks"][0].items() if k != "predictions"}]
        
        value_str = "NaN" if pearson_r is None else str(pearson_r)

        evaluation_output = pd.DataFrame([{
            'evaluator_name': EVALUATOR_NAME,
            'description': description,
            'predictor_name': predictor_name,
            'time_stamp': timestamp,
            'metric': 'pearson_r',
            'value': value_str,
            'prediction_task(s)_data': prediction_task_data_onlyinfo
        }])
        
        # CHANGE (v1): append to the single summary file (header only when new).
        _save_df_to_csv(evaluation_output, summary_filepath)

        return pearson_r

    except Exception as e:
        print(f"An unexpected error occurred during evaluation calculations: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return None



