# Behavioral Analyses

Scripts beginning with GNG_ implement behavioral measures of Go/No-Go performance:

GNG_bin… — Generate PSTHs for neuronal firing aligned to trial events

GNG_dprime — Sensitivity (d′) across conditions

GNG_lick_rt_dprime — Lick reaction time and discrimination performance

GNG_raster — Trial raster plotting

GNG_smoothed_PSTH — Smoothed peri-stimulus time histograms

GNG_trial_outcomes — Trial outcome separation

# Electrophysiological & Neural Analyses

Figure_3_ephys_kra_GLM.m — GLM-based encoding analysis of neural activity

Figure_4_ephys_PCA.m — Principal component analysis of population activity

Figure_4_ephys_SVM.m — Decoding analysis (support vector machine)

# Optogenetic & Imaging Analyses

Figure_2_optogenetics.m — Optogenetic manipulation analysis

Figure_5_MGreenLantern_fluorescence.m — Fluorescence images quantification

Figure_5_rAAV_count.m — AAV expression counting

Figure_6_MGreen_structural_imaging.m — Structural imaging metrics

# Python Utilities

analyze_axon_peaks_interactive.py — Interactive axon peak detection

analyze_bouton_density.py — Bouton density estimation

mark_axons_loop.py — Manual axon annotation

# Requirements
# Software

MATLAB (R2020a or later recommended)

Statistics and Machine Learning Toolbox

Signal Processing Toolbox (optional, for advanced PSTHs)

Python 3.7+ for the axon and bouton analysis utilities

matplotlib

numpy

scipy

tqdm (optional)

opencv-python (optional, for image inspection)

Before running Python scripts, install dependencies (example):
````
pip install matplotlib numpy scipy tqdm opencv-python
````

Analysis Workflow

Set up paths

Add repo folder to MATLAB path.

Confirm data folders are accessible (not included in the repository).

Run analyses in order
The analysis is designed to proceed figure by figure:
````
% Behavioral & recording summary
run('Figure_1_behavior_recording.m')

% Optogenetic manipulation results
run('Figure_2_optogenetics.m')

% Electrophysiology
run('Figure_3_ephys_kra_GLM.m')
run('Figure_4_ephys_PCA.m')
run('Figure_4_ephys_SVM.m')

% Imaging & structure
run('Figure_5_MGreenLantern_fluorescence.m')
run('Figure_5_rAAV_count.m')
run('Figure_6_MGreen_structural_imaging.m')
````

# Reproducibility

Random seed is set consistently for any stochastic procedures.

All analyses are script-driven with no GUI dependence.

Intermediate results can be saved by uncommenting the save() commands in each script.

