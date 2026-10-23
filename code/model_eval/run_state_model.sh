#!/bin/bash
# Model pilot (Fig. 7a, Supplementary Table S6): Arc models bundled with arc-state (perturb_mean, pseudobulk), same data,
# TOML and training settings as the donor-1 STATE run of run_transfer.sh (run in $WORK/model_eval/state_D1).
# Usage: run_state_model.sh <model> <gpu>; then `state tx predict --output-dir runs/cd4_d1_<model> --checkpoint final.ckpt
# --predict-only` and pilot_state_score.py <data_dir> <adata_pred.h5ad> <out_json>.
M=$1; G=$2; SD=$WORK/model_eval/state_D1; ST=${STATE_BIN:-state}
cd $SD
export CUDA_VISIBLE_DEVICES=$G NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 WANDB_MODE=disabled
$ST tx train \
  data.kwargs.toml_config_path=$SD/cd4_fewshot.toml \
  data.kwargs.embed_key=X_hvg \
  data.kwargs.output_space=gene \
  data.kwargs.pert_col=gene \
  data.kwargs.cell_type_key=cell_type \
  data.kwargs.batch_col=gem_group \
  data.kwargs.control_pert=non-targeting \
  data.kwargs.num_workers=8 \
  training.max_steps=40000 \
  training.batch_size=8 \
  use_wandb=false \
  model=$M \
  output_dir=$SD/runs \
  name=cd4_d1_$M > $SD/train_$M.log 2>&1
echo "TRAIN EXIT $?" >> $SD/train_$M.log
