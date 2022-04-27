#!/bin/bash

#                           datfile         method      K  select.genes outfile
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  20 TRUE         pijuan-sel-topics-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         20 TRUE         pijuan-sel-nmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    20 TRUE         pijuan-sel-snnebmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     20 TRUE         pijuan-sel-nnebmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  30 TRUE         pijuan-sel-topics-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         30 TRUE         pijuan-sel-nmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    30 TRUE         pijuan-sel-snnebmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     30 TRUE         pijuan-sel-nnebmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  40 TRUE         pijuan-sel-topics-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         40 TRUE         pijuan-sel-nmf-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    40 TRUE         pijuan-sel-snnebmf-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     40 TRUE         pijuan-sel-nnebmf-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  20 FALSE        pijuan-nosel-topics-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         20 FALSE        pijuan-nosel-nmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    20 FALSE        pijuan-nosel-snnebmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     20 FALSE        pijuan-nosel-nnebmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  30 FALSE        pijuan-nosel-topics-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         30 FALSE        pijuan-nosel-nmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    30 FALSE        pijuan-nosel-snnebmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     30 FALSE        pijuan-nosel-nnebmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics  40 FALSE        pijuan-nosel-topics-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds nmf         40 FALSE        pijuan-nosel-nmf-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds snn-ebmf    40 FALSE        pijuan-nosel-snnebmf-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds nn-ebmf     40 FALSE        pijuan-nosel-nnebmf-k=40
