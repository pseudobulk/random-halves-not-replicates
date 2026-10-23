"""Model pilot: load SE-600M from its published model.safetensors (2.9 GB) instead of the 11.5 GB Lightning checkpoints
that `state emb transform` requires (the 5090 disk is nearly full). The model is built exactly as arc-state's trainer
builds it (state/emb/train/trainer.py: token_dim = embedding size, d_model = emsize, nhead, d_hid, nlayers, output_dim
from config.yaml), the weights are loaded strictly, and the remaining steps of Inference.load_model are repeated
(inference precision, protein-embedding table from protein_embeddings.pt, eval mode). Embeddings then come from arc-state's
own Inference.encode_adata. Usage: embed(se_dir, input_h5ad, batch_size) -> (n_cells, output_dim) array."""
import os
import torch
from torch import nn
from omegaconf import OmegaConf
from safetensors.torch import load_file
from state.emb.inference import Inference
from state.emb.nn.model import StateEmbeddingModel
from state.emb.utils import get_embedding_cfg, get_precision_config


def load(se_dir):
    cfg = OmegaConf.load(os.path.join(se_dir, "config.yaml"))
    model = StateEmbeddingModel(token_dim=get_embedding_cfg(cfg).size, d_model=cfg.model.emsize, nhead=cfg.model.nhead,
                                d_hid=cfg.model.d_hid, nlayers=cfg.model.nlayers, output_dim=cfg.model.output_dim,
                                dropout=0.0, compiled=False, emb_size=get_embedding_cfg(cfg).size, cfg=cfg)
    sd = load_file(os.path.join(se_dir, "model.safetensors"))
    sd = {k: v for k, v in sd.items() if not k.startswith("pe_embedding.")}   # replaced below, as in load_model
    missing, unexpected = model.load_state_dict(sd, strict=False)
    missing = [k for k in missing if not k.startswith("pe_embedding.")]
    assert not missing and not unexpected, "weights do not match the model: missing %s unexpected %s" % (missing[:5], unexpected[:5])
    device = "cuda" if torch.cuda.is_available() else "cpu"
    precision = get_precision_config(device_type=device)
    model = model.to(device).to(precision)
    pe = torch.load(os.path.join(se_dir, "protein_embeddings.pt"), map_location="cpu", weights_only=False)
    all_pe = torch.vstack(list(pe.values())) if isinstance(pe, dict) else pe
    model.pe_embedding = nn.Embedding.from_pretrained(all_pe)
    model.pe_embedding.to(model.device, dtype=precision)
    model.binary_decoder.requires_grad = False
    model.eval()
    inf = Inference(cfg=cfg, protein_embeds=pe)
    inf.model = model; inf._vci_conf = cfg
    return inf


def embed(se_dir, input_h5ad, batch_size=128):
    inf = load(se_dir)
    return inf.encode_adata(input_adata_path=input_h5ad, output_adata_path=None, emb_key="X_state", batch_size=batch_size)
