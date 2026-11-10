# -*- coding: utf-8 -*-

import numpy as np
from tqdm import tqdm
import json
import pyBigWig
import pandas as pd
import os
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr

# ===== CONFIGURATION =====
batching = False 

# ===== FUNCTIONS =====

def compute_log_total_count_correlation(pred_sum, truth_sum):
    """Compute Pearson r on log-transformed total counts per region."""
    pred_log = np.log1p(pred_sum)
    truth_log = np.log1p(truth_sum)
    
    if np.std(pred_log) > 0 and np.std(truth_log) > 0:
        corr, _ = pearsonr(pred_log, truth_log)
        return corr
    else:
        return np.nan


def compute_jensen_shannon_distance(pred_values, truth_values, pseudocount=1e-16):
    """Compute Jensen-Shannon distance between predicted and ground truth signal distributions."""
    pred = pred_values + pseudocount
    truth = truth_values + pseudocount
    
    pred_prob = pred / np.sum(pred)
    truth_prob = truth / np.sum(truth)
    
    js_distance = jensenshannon(pred_prob, truth_prob)
    
    return js_distance if not np.isnan(js_distance) else np.nan


def compute_binned_correlation(pred_values, truth_values, bin_size):
    """Compute Pearson correlation at specified binning resolution."""
    if bin_size == 1:
        if np.std(pred_values) > 0 and np.std(truth_values) > 0:
            return np.corrcoef(pred_values, truth_values)[0, 1]
        else:
            return np.nan
    
    effective_len = (len(pred_values) // bin_size) * bin_size
    
    if effective_len < bin_size:
        return np.nan
    
    pred_binned = pred_values[:effective_len].reshape(-1, bin_size).sum(axis=1)
    truth_binned = truth_values[:effective_len].reshape(-1, bin_size).sum(axis=1)
    
    if np.std(pred_binned) > 0 and np.std(truth_binned) > 0:
        return np.corrcoef(pred_binned, truth_binned)[0, 1]
    else:
        return np.nan


def get_evaluation_window(chromosome, start, end, window_size=1000):
    """Get a window of specified size centered at the peak (midpoint of interval)."""
    center = (start + end) // 2
    window_start = center - window_size // 2
    window_end = window_start + window_size
    
    return chromosome, window_start, window_end


def get_values_from_bigwig(bw, chromosome, start, end):
    """Extract values from a bigwig file for a given region."""
    values = np.nan_to_num(bw.values(chromosome, start, end, numpy=True))
    return values


def evaluate_sample(sample_name, bw_path, pred_bw_path, bed, 
                   batch_size, bin_size, output_dir, 
                   eval_window_size=1000, jsd_pseudocount=1e-16):
    """Evaluate a single sample and return metrics."""
    
    print(f"\n{'='*60}")
    print(f"Evaluating: {sample_name}")
    print(f"{'='*60}")
    
    regions_df = pd.read_table(bed, names=['chr', 'start', 'end', 'name', 'value', 'strand', 
                                           'value1', 'value2', 'value3', 'summit'])
    truth_bw = pyBigWig.open(bw_path)
    pred_bw = pyBigWig.open(pred_bw_path)
    
    if batching and len(regions_df) > batch_size:
        batch_df = regions_df.sample(n=batch_size, random_state=42)
        batch_df.reset_index(inplace=True, drop=True)
    else:
        batch_df = regions_df
    

    all_pred_values = []  
    all_truth_values = []  
    per_sample_sum_preds = []
    per_sample_sum_targs = []
    per_sample_js_distances = []
    total_pred_sum = 0.0
    total_target_sum = 0.0
    failed_regions = 0
    
    print(f"Processing {len(batch_df)} regions...")
    print(f"Evaluation window: {eval_window_size}bp centered at peak")
    print(f"JS pseudocount: {jsd_pseudocount}")
    
    with tqdm(total=len(batch_df), desc=f"  {sample_name}") as pbar:
        for row in range(len(batch_df)):
            try:
                chrom = batch_df.at[row, 'chr']
                region_start = batch_df.at[row, 'start']
                region_end = batch_df.at[row, 'end']
                
                eval_chrom, eval_start, eval_end = get_evaluation_window(
                    chrom, region_start, region_end, eval_window_size
                )
                
                pred_values = get_values_from_bigwig(
                    pred_bw, eval_chrom, eval_start, eval_end
                )
                
                truth_values = get_values_from_bigwig(
                    truth_bw, eval_chrom, eval_start, eval_end
                )
                
                all_pred_values.extend(pred_values.tolist())
                all_truth_values.extend(truth_values.tolist())
                pred_sum = np.sum(pred_values)
                truth_sum = np.sum(truth_values)
                per_sample_sum_preds.append(pred_sum)
                per_sample_sum_targs.append(truth_sum)
                
                if np.sum(pred_values) > 0 and np.sum(truth_values) > 0:
                    jsd = compute_jensen_shannon_distance(pred_values, truth_values, jsd_pseudocount)
                    per_sample_js_distances.append(jsd)
                else:
                    per_sample_js_distances.append(np.nan)
                
                total_pred_sum += float(pred_sum)
                total_target_sum += float(truth_sum)
                
            except Exception as e:
                print(f"\n  Error processing region {row}: {e}")
                failed_regions += 1
                # Fill with NaN for failed regions (will be excluded from global correlation)
                all_pred_values.extend([np.nan] * eval_window_size)
                all_truth_values.extend([np.nan] * eval_window_size)
                per_sample_sum_preds.append(0)
                per_sample_sum_targs.append(0)
                per_sample_js_distances.append(np.nan)
            
            pbar.update(1)
    
    truth_bw.close()
    pred_bw.close()
    
    if failed_regions > 0:
        print(f"\n  Warning: {failed_regions} regions failed processing")
    
    # ===== COMPUTE GLOBAL PEARSON CORRELATION =====
    # Remove any NaN values (from failed regions)
    pred_array = np.array(all_pred_values)
    truth_array = np.array(all_truth_values)
    valid_mask = ~(np.isnan(pred_array) | np.isnan(truth_array))
    pred_clean = pred_array[valid_mask]
    truth_clean = truth_array[valid_mask]
    
    # Compute global Pearson r on all concatenated bases
    if len(pred_clean) > 1 and np.std(pred_clean) > 0 and np.std(truth_clean) > 0:
        global_pearson_r, global_pearson_pvalue = pearsonr(pred_clean, truth_clean)
    else:
        global_pearson_r = np.nan
        global_pearson_pvalue = np.nan
    
    # Also compute binned correlation if bin_size > 1
    if bin_size > 1 and len(pred_clean) >= bin_size:
        global_pearson_r_binned = compute_binned_correlation(pred_clean, truth_clean, bin_size)
    else:
        global_pearson_r_binned = global_pearson_r
    
    # Log total count correlation (per-region)
    if len(per_sample_sum_preds) > 1 and len(per_sample_sum_targs) > 1:
        log_count_corr = compute_log_total_count_correlation(
            np.array(per_sample_sum_preds),
            np.array(per_sample_sum_targs)
        )
    else:
        log_count_corr = np.nan
    
    # Compute aggregate JS distance
    valid_jsd = [j for j in per_sample_js_distances if not np.isnan(j)]
    
    results = {
        "sample_name": sample_name,
        "n_regions": len(batch_df),
        "n_failed_regions": failed_regions,
        "n_valid_regions": len(batch_df) - failed_regions,
        "n_valid_bases": len(pred_clean),
        "eval_window_size": eval_window_size,
        "jsd_pseudocount": jsd_pseudocount,
        f"global_pearson_r_{bin_size}bp": float(global_pearson_r_binned) if not np.isnan(global_pearson_r_binned) else np.nan,
        f"global_pearson_r_1bp": float(global_pearson_r) if not np.isnan(global_pearson_r) else np.nan,
        f"global_pearson_pvalue_1bp": float(global_pearson_pvalue) if not np.isnan(global_pearson_pvalue) else np.nan,
        "log_total_count_pearson_r": float(log_count_corr) if not np.isnan(log_count_corr) else np.nan,
        "mean_js_distance": float(np.mean(valid_jsd)) if valid_jsd else np.nan,
        "median_js_distance": float(np.median(valid_jsd)) if valid_jsd else np.nan,
        "total_pred_sum": float(total_pred_sum),
        "total_truth_sum": float(total_target_sum),
        "ratio_pred_truth": float(total_pred_sum / total_target_sum) if total_target_sum > 0 else np.nan,
        "total_bases_evaluated": len(pred_clean),
        "per_sample_js_distances": per_sample_js_distances,
        "per_sample_pred_sums": per_sample_sum_preds,
        "per_sample_truth_sums": per_sample_sum_targs
    }
    
    # Save individual sample results
    sample_output_dir = os.path.join(output_dir, sample_name.replace(' ', '_'))
    os.makedirs(sample_output_dir, exist_ok=True)
    
    with open(os.path.join(sample_output_dir, f"{sample_name}_results.json"), "w") as f:
        results_serializable = {}
        for k, v in results.items():
            if isinstance(v, np.ndarray):
                results_serializable[k] = v.tolist()
            elif isinstance(v, (np.float32, np.float64)):
                results_serializable[k] = float(v)
            elif isinstance(v, list):
                results_serializable[k] = [float(x) if isinstance(x, (np.floating,)) else x for x in v]
            else:
                results_serializable[k] = v
        json.dump(results_serializable, f, indent=2)
    
    # Print summary
    print(f"\n  Results for {sample_name}:")
    print(f"  Valid regions: {results['n_valid_regions']}/{results['n_regions']}")
    print(f"  Total bases evaluated: {results['total_bases_evaluated']:,}")
    print(f"  Global Pearson r (1bp): {results['global_pearson_r_1bp']:.4f}" if not np.isnan(results['global_pearson_r_1bp']) else "  Global Pearson r (1bp): N/A")
    print(f"  Global Pearson r ({bin_size}bp): {results[f'global_pearson_r_{bin_size}bp']:.4f}" if not np.isnan(results[f'global_pearson_r_{bin_size}bp']) else f"  Global Pearson r ({bin_size}bp): N/A")
    print(f"  Log total count Pearson r: {results['log_total_count_pearson_r']:.4f}" if not np.isnan(results['log_total_count_pearson_r']) else "  Log total count Pearson r: N/A")
    print(f"  Mean JS distance: {results['mean_js_distance']:.4f}" if not np.isnan(results['mean_js_distance']) else "  Mean JS distance: N/A")
    
    return results


