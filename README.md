# EDiTS - Embayment Decompression in Two Stages
A MATLAB code (written in v. R2020a) for modeling diffusion of H2O and CO2 through rhyolitic melt embayments during two-stage, continuous magma
decompression. The example is set up to fit H2O and CO2 concentration gradients measured along a XX-um long quartz-hosted rhyolitic melt embayment that was
experimentally decompressed from 150 to XX MPa at 780 C and a rate of XX MPa/s.

Instructions for running code:

Download both emb_diffusion_ts.m (master script) and diffusion_function_ts.m (function) into a folder.

emb_diffusion_ts.m is the main code where all input parameters are defined and the misfit is calculated. On Lines 22-30, specify input parameters, including initial pressure (MPa), final pressure (MPa), embayment length (um), bubble radius, initial H2O (wt. %), melt density (kg/m^3), temperature (C), and uncertainties. On lines 37-39, enter H2O +/- CO2 concentrations (wt. % and ppm, respectively), as well as distance along the measured profile (um). All arrays go from embayment outlet to the interior. On lines 43-45, specify the range of initial decompression rates (MPa/s), final decompression rates (MPa/s), and transition pressures (MPa) to iterate through.

diffusion_function_ts.m is the function that performs the diffusion calculation. User does not need to modify.

example_output.tiff shows the best fit model profiles to the measured H2O and CO2 concentration gradients for the example embayment.
