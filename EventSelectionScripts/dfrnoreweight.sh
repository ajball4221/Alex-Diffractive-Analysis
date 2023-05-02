#!/bin/sh
python CC-NuE-XSec/selection/gridSelection.py --playlist me5A_NCDiffSignal_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --mc_only --unweighted
python CC-NuE-XSec/selection/gridSelection.py --playlist me6A_NCDiffSignal_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --mc_only --unweighted
python CC-NuE-XSec/selection/gridSelection.py --playlist me6B_NCDiffSignal_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --mc_only --unweighted
python CC-NuE-XSec/selection/gridSelection.py --playlist me6C_NCDiffSignal_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --mc_only --unweighted
python CC-NuE-XSec/selection/gridSelection.py --playlist me6D_NCDiffSignal_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --mc_only --unweighted
python CC-NuE-XSec/selection/gridSelection.py --playlist me6E_NCDiffSignal_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --mc_only --unweighted
python CC-NuE-XSec/selection/gridSelection.py --playlist me6F_NCDiffSignal_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --mc_only --unweighted
python CC-NuE-XSec/selection/gridSelection.py --playlist me6G_NCDiffSignal_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --mc_only --unweighted
jobsub_q --user ajball
