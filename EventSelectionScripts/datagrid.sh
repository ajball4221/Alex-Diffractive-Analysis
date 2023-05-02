#!/bin/bash
python CC-NuE-XSec/selection/gridSelection.py --playlist me5A_Audit_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --data_only
python CC-NuE-XSec/selection/gridSelection.py --playlist me6A_Audit_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --data_only
python CC-NuE-XSec/selection/gridSelection.py --playlist me6B_Audit_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --data_only
python CC-NuE-XSec/selection/gridSelection.py --playlist me6C_Audit_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --data_only
python CC-NuE-XSec/selection/gridSelection.py --playlist me6D_Audit_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --data_only
python CC-NuE-XSec/selection/gridSelection.py --playlist me6E_Audit_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --data_only
python CC-NuE-XSec/selection/gridSelection.py --playlist me6F_Audit_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --data_only
python CC-NuE-XSec/selection/gridSelection.py --playlist me6G_Audit_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --data_only
python CC-NuE-XSec/selection/gridSelection.py --playlist me6H_Audit_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --data_only
python CC-NuE-XSec/selection/gridSelection.py --playlist me6I_Audit_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --data_only
python CC-NuE-XSec/selection/gridSelection.py --playlist me6J_Audit_nx --ntuple_tag fspline --use-sideband Pi0 Coherent Electron_Neutrino --truth --cal_POT --data_only
jobsub_q --user ajball
