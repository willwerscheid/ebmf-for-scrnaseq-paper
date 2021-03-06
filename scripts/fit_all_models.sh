#!/bin/bash

#                           datfile         method        K  select.genes outfile
sbatch fit_one_model.sbatch pijuan-sala.rds nmf-log       20 TRUE         pijuan-sel-nmf-log-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds nmf-identity  20 TRUE         pijuan-sel-nmf-id-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics    20 TRUE         pijuan-sel-topics-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds pcmf          20 TRUE         pijuan-sel-pcmf-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds ebmf-log      20 TRUE         pijuan-sel-ebmf-log-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds ebmf-identity 20 TRUE         pijuan-sel-ebmf-id-k=20
sbatch fit_one_model.sbatch pijuan-sala.rds nmf-log       30 TRUE         pijuan-sel-nmf-log-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nmf-identity  30 TRUE         pijuan-sel-nmf-id-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics    30 TRUE         pijuan-sel-topics-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds pcmf          30 TRUE         pijuan-sel-pcmf-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds ebmf-log      30 TRUE         pijuan-sel-ebmf-log-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds ebmf-identity 30 TRUE         pijuan-sel-ebmf-id-k=30
sbatch fit_one_model.sbatch pijuan-sala.rds nmf-log       40 TRUE         pijuan-sel-nmf-log-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds nmf-identity  40 TRUE         pijuan-sel-nmf-id-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds fasttopics    40 TRUE         pijuan-sel-topics-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds pcmf          40 TRUE         pijuan-sel-pcmf-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds ebmf-log      40 TRUE         pijuan-sel-ebmf-log-k=40
sbatch fit_one_model.sbatch pijuan-sala.rds ebmf-identity 40 TRUE         pijuan-sel-ebmf-id-k=40
sbatch fit_one_model.sbatch montoro.rds     nmf-log       20 TRUE         montoro-sel-nmf-log-k=20
sbatch fit_one_model.sbatch montoro.rds     nmf-identity  20 TRUE         montoro-sel-nmf-id-k=20
sbatch fit_one_model.sbatch montoro.rds     fasttopics    20 TRUE         montoro-sel-topics-k=20
sbatch fit_one_model.sbatch montoro.rds     pcmf          20 TRUE         montoro-sel-pcmf-k=20
sbatch fit_one_model.sbatch montoro.rds     ebmf-log      20 TRUE         montoro-sel-ebmf-log-k=20
sbatch fit_one_model.sbatch montoro.rds     ebmf-identity 20 TRUE         montoro-sel-ebmf-id-k=20
sbatch fit_one_model.sbatch montoro.rds     nmf-log       30 TRUE         montoro-sel-nmf-log-k=30
sbatch fit_one_model.sbatch montoro.rds     nmf-identity  30 TRUE         montoro-sel-nmf-id-k=30
sbatch fit_one_model.sbatch montoro.rds     fasttopics    30 TRUE         montoro-sel-topics-k=30
sbatch fit_one_model.sbatch montoro.rds     pcmf          30 TRUE         montoro-sel-pcmf-k=30
sbatch fit_one_model.sbatch montoro.rds     ebmf-log      30 TRUE         montoro-sel-ebmf-log-k=30
sbatch fit_one_model.sbatch montoro.rds     ebmf-identity 30 TRUE         montoro-sel-ebmf-id-k=30
