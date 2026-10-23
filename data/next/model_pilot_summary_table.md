| Model | Unexpressed A | Unexpressed B | A − B (SE) | Knockdown A | Knockdown B | A − B (SE) | Note |
|---|---|---|---|---|---|---|---|
| Context mean | 0.081 | 0.081 | 0.000 (0.012) | 0.106 | 0.106 | 0.000 (0.012) | no perturbation information |
| Copy resting effect | 0.220 | -0.001 | 0.221 (0.008) | 0.228 | 0.023 | 0.205 (0.009) |  |
| Additive | 0.226 | 0.007 | 0.220 (0.008) | 0.239 | 0.034 | 0.204 (0.009) |  |
| Linear (ridge) | 0.463 | 0.020 | 0.443 (0.020) | 0.457 | 0.042 | 0.415 (0.020) |  |
| Arc perturb_mean | 0.165 | 0.007 | 0.158 (0.009) | 0.179 | 0.026 | 0.153 (0.010) |  |
| Arc pseudobulk | 0.222 | 0.006 | 0.216 (0.012) | 0.238 | 0.030 | 0.208 (0.014) |  |
| STATE Transition (HVG) | 0.260 | 0.005 | 0.254 (0.016) | 0.260 | 0.020 | 0.241 (0.016) |  |
| State (SE-600M + ST) | 0.280 | 0.014 | 0.266 (0.016) | 0.297 | 0.038 | 0.259 (0.018) | SE loaded from safetensors |
| PIE (expression shift) | 0.228 | 0.108 | 0.120 (0.015) | 0.233 | 0.122 | 0.111 (0.018) | without ncbi_text source |
| PIE (log2 fold change) | 0.026 | 0.011 | 0.015 (0.003) | 0.029 | 0.014 | 0.015 (0.003) |  |
| Stack (zero-shot) | 0.038 | 0.003 | 0.036 (0.014) | 0.023 | -0.005 | 0.028 (0.015) | 1,672 of 2,152 genes |
| biolord (adapted) | 0.030 | 0.004 | 0.027 (0.026) | 0.020 | -0.033 | 0.053 (0.028) | fit r = 0.35 |
| biolord (tutorial) | 0.005 | 0.010 | -0.006 (0.022) | -0.035 | -0.031 | -0.005 (0.023) | fit r = 0.025; did not learn guides |
| CPA (default) | -0.001 | 0.003 | -0.003 (0.017) | -0.031 | -0.000 | -0.031 (0.020) | fit r = 0.003; did not learn guides |
| CPA (adapted) | 0.001 | -0.002 | 0.003 (0.017) | -0.021 | -0.003 | -0.018 (0.020) | fit r = 0.016; did not learn guides |

Reference: the two guides' stimulated effects agree at r = 0.041 (unexpressed) and 0.059 (knockdown).
