#!/bin/bash
# Model pilot: STATE on SE-600M embeddings (the full State recipe): as the donor-1 STATE run (README example settings,
# model=state, same fewshot splits), with embed_key=X_state as input and output_space=gene (targets read from obsm X_hvg
# = pilot genes). Inputs from pilot_se_embed.py <data_dir> <se_dir> $WORK/model_eval/state_se (writes data/cd4_d1_se.h5ad;
# the TOML cd4_fewshot_se.toml is the donor-1 TOML pointed at that folder). Usage: run_state_se.sh [gpu]
SD=$WORK/model_eval/state_se; ST=${STATE_BIN:-state}
cd $SD
export CUDA_VISIBLE_DEVICES=${1:-0} NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 WANDB_MODE=disabled
$ST tx train \
  data.kwargs.toml_config_path=$SD/cd4_fewshot_se.toml \
  data.kwargs.embed_key=X_state \
  data.kwargs.output_space=gene \
  data.kwargs.pert_col=gene \
  data.kwargs.cell_type_key=cell_type \
  data.kwargs.batch_col=gem_group \
  data.kwargs.control_pert=non-targeting \
  data.kwargs.num_workers=8 \
  training.max_steps=40000 \
  training.batch_size=8 \
  use_wandb=false \
  model=state \
  output_dir=$SD/runs \
  name=cd4_d1_se > $SD/train.log 2>&1
echo "TRAIN EXIT $?" >> $SD/train.log
