#!/bin/bash
# Model pilot: PIE (arc-pie 1.0.0) with the replogle_wdataset recipe values (delta_p bin 1.1, evidence seed 42 / chunk 256,
# max_steps 5000, effective batch 16 x 2 GPUs x 4 = 128 kept as 16 x 1 GPU x 8); the published knowledge sources except
# ncbi_text (43 GB, not downloaded), placed in $PIE_DIR/sources_pub. Inputs from pilot_pie_prep.py <data_dir>
# $WORK/model_eval/state_D1/prep.json $PIE_DIR/cd4. Then `pie eval` (save_predictions=true, checkpoint best_auprc) and
# pilot_pie_score.py. Environment: PIE_DIR (default $WORK/model_eval/pie_data), PIE_BIN (the `pie` command).
# Usage: run_pie.sh [gpu]
PIE_DIR=${PIE_DIR:-$WORK/model_eval/pie_data}
printf 'WANDB_ENTITY=none\nWANDB_PROJECT=none\nPIE_DATA_ROOT=%s\nPIE_RUNS_ROOT=%s/pie_runs\nPIE_CACHE_DIR=%s/pie_cache\n' "$PIE_DIR" "$WORK/model_eval" "$WORK/model_eval" > $PIE_DIR/common.sh
export PIE_ENV_FILE=$PIE_DIR/common.sh CUDA_VISIBLE_DEVICES=${1:-0} NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
S=$PIE_DIR/sources_pub; P=$PIE_DIR/cd4
cd $PIE_DIR
${PIE_BIN:-pie} train experiment_name=cd4_pilot \
  "data.preprocessed_dirs=[$P/preprocessed]" data.dataset_weights=null data.split_dir=$P/splits \
  "data.source_dirs={esm2:$S/esm2,string_space:$S/string_space,depmap_gene_effect:$S/depmap_gene_effect,context_text:$S/context_text,perturbation_text:$S/perturbation_text}" \
  data.gene_text_dir=$S/gene_text data.aliases_path=$P/aliases.yaml \
  data.delta_p.bin_width_fold_change=1.1 data.evidence.seed=42 data.evidence.chunk=256 \
  trainer.devices=1 trainer.num_nodes=1 trainer.max_steps=5000 trainer.accumulate_grad_batches=8 \
  logger.enabled=false > $P/train.log 2>&1
echo "TRAIN EXIT $?" >> $P/train.log
