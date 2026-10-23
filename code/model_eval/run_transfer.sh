#!/bin/bash
# Transfer-model evaluation, one run (Fig. 7, Supplementary Table S6): resting -> stimulated 48 h for one CD4 donor
# (env PILOT_DONOR), or source -> target of another pair (env PILOT_STATES, e.g. K562,RPE1 or Control,IFNg).
# With the pilot's code and settings unchanged: baselines (pilot_baselines.py); STATE Transition (pilot_state_prep.py;
# `state tx train` with the settings of the arc-state README example: model=state, 40,000 steps, batch 8, input X_hvg =
# the pilot genes); prediction; scoring (pilot_state_score.py, except in mode local, and me_score_detail.py); export of the
# predictions (me_export_preds.py) and the STATE fit check (run_fitcheck.sh). Large intermediates (STATE input h5ad,
# adata_real) are deleted at the end.
# Usage: run_transfer.sh <data_dir> <name> <gpu> [full | local | baselines]
#   data_dir  pilot-format folder (pilot_<state>.npz, pilot_design.json) from cd4_pilot_extract.py,
#             me_replogle_extract.py or me_frangieh_extract.py
#   full      CD4 donors and K562 -> RPE1; local: Frangieh folds; baselines: baselines only
# Environment: WORK (outputs under $WORK/model_eval), STATE_PYTHON (python of the arc-state environment), STATE_BIN
# (its `state` command). Originally run as run_5090_donor.sh / run_5090_pair.sh / run_5090_local.sh followed by
# run_post_5090.sh on a GPU server, with the inputs copied from the extraction server; the commands are the same.
DD=$1; N=$2; G=$3; MODE=${4:-full}
R=$WORK/model_eval; SD=$R/state_$N; O=$R/out_$N; C=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PY=${STATE_PYTHON:-python3}; ST=${STATE_BIN:-state}
export CUDA_VISIBLE_DEVICES=$G NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 WANDB_MODE=disabled PYTHONPATH=$C
mkdir -p $SD $O
$PY $C/pilot_baselines.py $DD $O/baselines.json || { echo "$(date) baselines failed"; exit 1; }
[ "$MODE" = baselines ] && { echo "$(date) ${N}_BASELINES_DONE"; exit 0; }
$PY $C/pilot_state_prep.py $DD $SD || { echo "$(date) state prep failed"; exit 1; }
cd $SD
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
  model=state \
  output_dir=$SD/runs \
  name=cd4_$N > $SD/train.log 2>&1 || { echo "$(date) train failed"; exit 1; }
echo "$(date) $N trained"
$ST tx predict --output-dir runs/cd4_$N --checkpoint final.ckpt --predict-only > $SD/predict.log 2>&1 || { echo "$(date) predict failed"; exit 1; }
P=$SD/runs/cd4_$N/eval_final.ckpt/adata_pred.h5ad
if [ "$MODE" = full ]; then
  $PY $C/pilot_state_score.py $DD $P $O/state_st_scores.json || { echo "$(date) state score failed"; exit 1; }
fi
$PY $C/me_score_detail.py $DD $O/baseline_preds.npz $P $O/detail.json || { echo "$(date) detail failed"; exit 1; }
$PY $C/me_export_preds.py $DD $O/baseline_preds.npz $P $O/preds_$N.npz || { echo "$(date) export failed"; exit 1; }
rm -f $SD/data/cd4_d1.h5ad $SD/runs/cd4_$N/eval_final.ckpt/adata_real.h5ad
bash $C/run_fitcheck.sh $SD/runs/cd4_$N $DD $R/fit_$N $O/preds_$N.npz $O/fitcheck.json $G
echo "$(date) ${N}_DONE"
