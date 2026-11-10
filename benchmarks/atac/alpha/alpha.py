# -*- coding: utf-8 -*-

import numpy as np
from tqdm import tqdm
import json
import pyBigWig
from alphagenome.data import genome
from alphagenome.models import dna_client
import pandas as pd
import os
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr

# ===== CONFIGURATION =====

CACHE = None
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


def get_evaluation_window(interval, window_size=1000):
    """Get a window of specified size centered at the peak (midpoint of interval)."""
    center = (interval.start + interval.end) // 2
    window_start = center - window_size // 2
    window_end = window_start + window_size
    
    return genome.Interval(
        chromosome=interval.chromosome,
        start=window_start,
        end=window_end
    )


def predict_with_reverse_complement(model, pred_interval, eval_interval, biosample, data_type):
    """
    Run prediction on both forward and reverse complement strands and average.
    
    Args:
        model: AlphaGenome model
        pred_interval: Full prediction interval (resized to model input length)
        eval_interval: Evaluation window interval
        biosample: Biosample ontology term
        data_type: Data type for prediction
    
    Returns:
        Averaged prediction values for the evaluation window
    """
    path = CACHE / f'{eval_interval.chromosome}_{eval_interval.start}_{eval_interval.end}.npz' if CACHE else None
    if path and path.exists():
        with np.load(path) as saved:
            return saved['prediction']
    # ===== FORWARD STRAND PREDICTION =====
    interval_forward = genome.Interval(
        chromosome=pred_interval.chromosome,
        start=pred_interval.start,
        end=pred_interval.end,
        strand='+'
    )
    
    pred_forward = model.predict_interval(
        interval=interval_forward,
        requested_outputs=[dna_client.OutputType[str(data_type)]],
        ontology_terms=[biosample],
    )
    
    pred_forward_eval = getattr(pred_forward, data_type.lower()).slice_by_interval(
        eval_interval, match_resolution=True
    )
    pred_values_forward = np.asarray(pred_forward_eval.values).squeeze()
    
    # ===== REVERSE STRAND PREDICTION =====
    interval_reverse = genome.Interval(
        chromosome=pred_interval.chromosome,
        start=pred_interval.start,
        end=pred_interval.end,
        strand='-'
    )
    
    pred_reverse = model.predict_interval(
        interval=interval_reverse,
        requested_outputs=[dna_client.OutputType[str(data_type)]],
        ontology_terms=[biosample],
    )
    
    pred_reverse_eval = getattr(pred_reverse, data_type.lower()).slice_by_interval(
        eval_interval, match_resolution=True
    )
    pred_values_reverse = np.asarray(pred_reverse_eval.values).squeeze()
    
    # Average forward and reverse complement predictions
    pred_values_avg = (pred_values_forward + pred_values_reverse) / 2.0
    
    if pred_values_avg.shape != (1000,) or not np.isfinite(pred_values_avg).all():
        raise ValueError('Invalid AlphaGenome prediction')
    if path:
        temporary = path.with_suffix('.tmp.npz')
        np.savez_compressed(temporary, prediction=pred_values_avg)
        temporary.replace(path)
    return pred_values_avg


def evaluate_sample(sample_name, biosample, bw_path, bed, data_type, 
                   model, seq_length, batch_size, bin_size, output_dir, 
                   eval_window_size=1000, jsd_pseudocount=1e-16):
    """Evaluate a single sample and return metrics."""
    
    print(f"\n{'='*60}")
    print(f"Evaluating: {sample_name}")
    print(f"{'='*60}")
    
    regions_df = pd.read_table(bed, names=['chr', 'start', 'end', 'name', 'value', 'strand', 
                                           'value1', 'value2', 'value3', 'summit'])
    truth_bw = pyBigWig.open(bw_path)
    
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
    print(f"Using reverse-complement averaging (forward + reverse strands)")
    print(f"JS pseudocount: {jsd_pseudocount}")
    
    with tqdm(total=len(batch_df), desc=f"  {sample_name}") as pbar:
        for row in range(len(batch_df)):
            try:
                orig_interval = genome.Interval(
                    chromosome=batch_df.at[row, 'chr'],
                    start=batch_df.at[row, 'start'],
                    end=batch_df.at[row, 'end']
                )
                
                eval_interval = get_evaluation_window(orig_interval, eval_window_size)
                pred_interval = eval_interval.resize(seq_length)
                
                pred_values = predict_with_reverse_complement(
                    model=model,
                    pred_interval=pred_interval,
                    eval_interval=eval_interval,
                    biosample=biosample,
                    data_type=data_type
                )
                
                truth_values = np.nan_to_num(truth_bw.values(
                    eval_interval.chromosome,
                    eval_interval.start,
                    eval_interval.end,
                    numpy=True
                ))
                
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
    
    if failed_regions > 0:
        print(f"\n  Warning: {failed_regions} regions failed processing")
    
    # ===== COMPUTE GLOBAL PEARSON CORRELATION =====
    # Remove any NaN values (from failed regions)
    pred_array = np.array(all_pred_values)
    truth_array = np.array(all_truth_values)
    
    valid_mask = ~(np.isnan(pred_array) | np.isnan(truth_array))
    pred_clean = pred_array[valid_mask]
    truth_clean = truth_array[valid_mask]
    
    if len(pred_clean) > 1 and np.std(pred_clean) > 0 and np.std(truth_clean) > 0:
        global_pearson_r, global_pearson_pvalue = pearsonr(pred_clean, truth_clean)
    else:
        global_pearson_r = np.nan
        global_pearson_pvalue = np.nan
    
    if bin_size > 1 and len(pred_clean) >= bin_size:
        global_pearson_r_binned = compute_binned_correlation(pred_clean, truth_clean, bin_size)
    else:
        global_pearson_r_binned = global_pearson_r
    
    if len(per_sample_sum_preds) > 1 and len(per_sample_sum_targs) > 1:
        log_count_corr = compute_log_total_count_correlation(
            np.array(per_sample_sum_preds),
            np.array(per_sample_sum_targs)
        )
    else:
        log_count_corr = np.nan
    
    valid_jsd = [j for j in per_sample_js_distances if not np.isnan(j)]
    
    results = {
        "sample_name": sample_name,
        "biosample": biosample,
        "data_type": data_type,
        "n_regions": len(batch_df),
        "n_failed_regions": failed_regions,
        "n_valid_regions": len(batch_df) - failed_regions,
        "n_valid_bases": len(pred_clean),
        "eval_window_size": eval_window_size,
        "jsd_pseudocount": jsd_pseudocount,
        "reverse_complement": True,
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
    
    # Print summary to terminal
    print(f"\n  Results for {sample_name}:")
    print(f"  Valid regions: {results['n_valid_regions']}/{results['n_regions']}")
    print(f"  Total bases evaluated: {results['total_bases_evaluated']:,}")
    print(f"  Global Pearson r (1bp): {results['global_pearson_r_1bp']:.4f}" if not np.isnan(results['global_pearson_r_1bp']) else "  Global Pearson r (1bp): N/A")
    print(f"  Global Pearson r ({bin_size}bp): {results[f'global_pearson_r_{bin_size}bp']:.4f}" if not np.isnan(results[f'global_pearson_r_{bin_size}bp']) else f"  Global Pearson r ({bin_size}bp): N/A")
    print(f"  Log total count Pearson r: {results['log_total_count_pearson_r']:.4f}" if not np.isnan(results['log_total_count_pearson_r']) else "  Log total count Pearson r: N/A")
    print(f"  Mean JS distance: {results['mean_js_distance']:.4f}" if not np.isnan(results['mean_js_distance']) else "  Mean JS distance: N/A")
    
    return results


def main():
    import argparse
    import hashlib
    from pathlib import Path
    global CACHE
    root = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description='Predict and evaluate IMR-90 ATAC with AlphaGenome.')
    parser.add_argument('--model-version', choices=['FOLD_0','FOLD_1'], required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--limit', type=int)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    CACHE = args.out / 'regions'; CACHE.mkdir(exist_ok=True)
    bed = root / 'data/atac/fold_1_test_regions.bed'
    truth = root / 'data/atac/ENCFF650PHA.bigWig'
    identity = {'model':args.model_version, 'limit':args.limit,
                'inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [Path(__file__),bed,truth]}}
    manifest = args.out / 'inputs.json'
    if manifest.exists() and json.loads(manifest.read_text()) != identity:
        raise ValueError('Inputs changed; use a new output folder')
    manifest.write_text(json.dumps(identity,indent=2)+'\n')
    if args.limit:
        selected = args.out / 'selected_regions.bed'
        pd.read_csv(bed,sep='\t',header=None).head(args.limit).to_csv(selected,sep='\t',index=False,header=False)
        bed = selected
    key = os.environ.get('ALPHAGENOME_API_KEY')
    if not key: raise ValueError('Set ALPHAGENOME_API_KEY')
    model = dna_client.create(key,model_version=dna_client.ModelVersion[args.model_version])
    result = evaluate_sample('ATAC_IMR-90','EFO:0001196',str(truth),str(bed),'ATAC',
                             model,1048576,500,1,str(args.out),1000,1e-16)
    summary = {k:v for k,v in result.items() if not k.startswith('per_sample_')}
    summary.update(model_version=args.model_version,held_out=args.model_version=='FOLD_1',
                   status='complete' if result['n_failed_regions']==0 else 'incomplete')
    (args.out/'metrics.json').write_text(json.dumps(summary,indent=2)+'\n')
    pd.DataFrame([summary]).to_csv(args.out/'summary.csv',index=False)
    if result['n_failed_regions']: raise RuntimeError('Some regions failed; rerun to retry missing predictions')


if __name__ == '__main__':
    main()
