#!/bin/bash
# Environments of the model evaluation (CUDA 12.8 wheels). The original environments were installed without version pins;
# the pins below are the versions recorded in the analysis (Supplementary Table S6 and the script docstrings).
# Prefix: $WORK/envs. Usage: mkenv_models.sh
E=$WORK/envs; mkdir -p $E; cd $E
# arc-state (STATE, perturb_mean, pseudobulk, SE-600M): python 3.11, torch 2.7.1, arc-state 0.11.1 (cell-load 0.10.4)
conda create -y -p $E/stateenv -c conda-forge --override-channels python=3.11
stateenv/bin/pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
stateenv/bin/pip install arc-state==0.11.1
# biolord 0.0.3: python 3.10
conda create -y -p $E/biolordenv -c conda-forge --override-channels python=3.10
biolordenv/bin/pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
biolordenv/bin/pip install biolord==0.0.3 scanpy
# CPA (cpa-tools 0.8.8): python 3.10
conda create -y -p $E/cpaenv -c conda-forge --override-channels python=3.10
cpaenv/bin/pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
cpaenv/bin/pip install cpa-tools==0.8.8 scanpy
# PIE (arc-pie 1.0.0): python 3.12
conda create -y -p $E/pieenv -c conda-forge --override-channels python=3.12
pieenv/bin/pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
pieenv/bin/pip install "arc-pie[all]==1.0.0"
# Stack (arc-stack) and the Stack-Large weights: python 3.10
conda create -y -p $E/stackenv -c conda-forge --override-channels python=3.10
stackenv/bin/pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
stackenv/bin/pip install arc-stack huggingface_hub
stackenv/bin/python -c "from huggingface_hub import snapshot_download as s; print(s('arcinstitute/Stack-Large', local_dir='$E/stack_weights'))"
