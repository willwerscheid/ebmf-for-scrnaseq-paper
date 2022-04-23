#!/bin/bash

#                           datfile         method      K select.genes outfile
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  3 TRUE         pijuan-sel-topics-k=5
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         3 TRUE         pijuan-sel-nmf-k=5
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    3 TRUE         pijuan-sel-snnebmf-k=5
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     3 TRUE         pijuan-sel-nnebmf-k=5
