# Margolis_foraging

The following files were used for the Margolis & Gordus manuscript:

N2_cumevents_matrix.mat: This matlab file was provided by Alejandro Lopez-Cruz, and is the quantification of reorientation events from:

López-Cruz, A., Sordillo, A., Pokala, N., Liu, Q., McGrath, P.T., and Bargmann, C.I. (2019). Parallel Multimodal Circuits Control an Innate Foraging Behavior. Neuron 102, 407-419.e8. https://doi.org/10.1016/j.neuron.2019.01.053.![image](https://github.com/user-attachments/assets/b85de0dd-aefb-4585-a6fd-fd2efd5adcfc)

The elements of N2_cumevents_matrix are the cumulative reorientions observed for 1631 worms. Rows = frames (3fps, 8100 frames= 45 min), Columns = Individuals.

Model_2025.m: This code uses the Gillespie algorithm to model the experimental data from N2_cumevents_matrix.
Fig2b.m: This file generates figure panels for figure 2b in Maroglis et al.
Fig2d.m: This file generates figure panels for figure 2d in Margolis et al.
M_values_data.m: This file iteratively generates several versions of the model with different M0 values, to generate the data used in Figure 3c of Margolis et al.
Jensen_Shannon.m: This file calculates the Jensen-Shannon divergence used in M_values_data.



