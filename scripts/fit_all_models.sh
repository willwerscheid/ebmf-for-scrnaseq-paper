#!/bin/bash

#                           datfile         method      K  select.genes outfile
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  10 TRUE         pijuan-sel-topics-k=10
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         10 TRUE         pijuan-sel-nmf-k=10
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    10 TRUE         pijuan-sel-snnebmf-k=10
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     10 TRUE         pijuan-sel-nnebmf-k=10
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  20 TRUE         pijuan-sel-topics-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         20 TRUE         pijuan-sel-nmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    20 TRUE         pijuan-sel-snnebmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     20 TRUE         pijuan-sel-nnebmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  30 TRUE         pijuan-sel-topics-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         30 TRUE         pijuan-sel-nmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    30 TRUE         pijuan-sel-snnebmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     30 TRUE         pijuan-sel-nnebmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  30 TRUE         pijuan-sel-topics-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         30 TRUE         pijuan-sel-nmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    30 TRUE         pijuan-sel-snnebmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     30 TRUE         pijuan-sel-nnebmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  10 FALSE        pijuan-sel-topics-k=10
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         10 FALSE        pijuan-sel-nmf-k=10
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    10 FALSE        pijuan-sel-snnebmf-k=10
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     20 FALSE        pijuan-sel-nnebmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  20 FALSE        pijuan-sel-topics-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         20 FALSE        pijuan-sel-nmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    20 FALSE        pijuan-sel-snnebmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     20 FALSE        pijuan-sel-nnebmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  30 FALSE        pijuan-sel-topics-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         30 FALSE        pijuan-sel-nmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    30 FALSE        pijuan-sel-snnebmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     30 FALSE        pijuan-sel-nnebmf-k=30
