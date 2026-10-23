#!/bin/bash
# Model evaluation: STATE fit check (me_fitcheck.py). Rebuilds the run's input with pilot_state_prep.py in a work folder,
# predicts NFIT training guides (+ NCHK test guides as a reproduction check) with the trained checkpoint via a TOML
# override, scores, and deletes the rebuilt input and the predictions. The original run folder is only read.
# Usage: run_fitcheck.sh <run_dir> <data_dir> <work_dir> <preds.npz> <out_json> <gpu>; env PILOT_DONOR or PILOT_STATES;
# STATE_PYTHON and STATE_BIN as in run_transfer.sh.
RUN=$1; DD=$2; WD=$3; PP=$4; OUT=$5; G=$6
C=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); PY=${STATE_PYTHON:-python3}; ST=${STATE_BIN:-state}
export CUDA_VISIBLE_DEVICES=$G NUMPY_MADVISE_HUGEPAGE=0 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 WANDB_MODE=disabled PYTHONPATH=$C
mkdir -p $WD/fitrun
$PY $C/pilot_state_prep.py $DD $WD || { echo "$(date) prep failed"; exit 1; }
$PY $C/me_fitcheck.py toml $WD $DD || { echo "$(date) toml failed"; exit 1; }
cp $RUN/config.yaml $WD/fitrun/
$ST tx predict --output-dir $WD/fitrun --checkpoint final.ckpt --toml $WD/fit.toml --predict-only > $WD/predict_fit.log 2>&1 || { echo "$(date) predict failed"; exit 1; }
$PY $C/me_fitcheck.py score $DD $WD/fitrun/eval_final.ckpt/adata_pred.h5ad $PP $OUT || { echo "$(date) score failed"; exit 1; }
rm -rf $WD/data $WD/fitrun/eval_final.ckpt
echo "$(date) FITCHECK_DONE"
