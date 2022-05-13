#!/bin/bash

#                           datfile  method        K  select.genes outfile
sbatch fit_one_model.sbatch deng.rds nmf-log       10 TRUE         deng-sel-nmf-log-k=10
sbatch fit_one_model.sbatch deng.rds nmf-identity  10 TRUE         deng-sel-nmf-id-k=10
sbatch fit_one_model.sbatch deng.rds fasttopics    10 TRUE         deng-sel-topics-k=10
sbatch fit_one_model.sbatch deng.rds pcmf          10 TRUE         deng-sel-pcmf-k=10
sbatch fit_one_model.sbatch deng.rds ebmf-log      10 TRUE         deng-sel-ebmf-log-k=10
sbatch fit_one_model.sbatch deng.rds ebmf-identity 10 TRUE         deng-sel-ebmf-id-k=10
