# KTB_Graded_Inh
Code for Chapter 3 of the Dissertation "Computational Models of Modulated Gamma Oscillations and Activity-Dependent Graded Inhibition in the Olfactory Bulb." 

Code is based on the model by Kersen, Tavoni and Balasubramanian (2022), (D. E. C. Kersen, G. Tavoni, and V. Balasubramanian. Connectivity and dynamics in the olfactory bulb. PLoS Comput Biol, 18(2):e1009856, 2022. ISSN 1553-7358. doi: 10.1371/journal.pcbi.1009856.). For their code, see https://github.com/dkersen/olfactory-bulb.

LFP50_*.m runs simulations of various iterations of the model

 
Fig_GPFI_Osc_Full.m calculates power spectra from the data and plots Figures 3.9 - 3.11, B.1 and B.2 from the dissertation. 

Fig_GPFI_gspiking.m calculates the average GC firing rate and spike participation and plots Figure 3.8 from the dissertation.

Fig_GPFI_mspiking.m calculates the average MC firing rate, MC spike participation, and MC spike distributions and plots Figures 3.7 and 3.12 from the dissertation.
Note: The figures for "Graded Only" used the same code, just loaded trial data from "Graded Only" trials.

Fig_SITI2_gSpiking.m plots Figure 3.3(c) from the dissertation

Fig_SITI2_mSpiking.m plots Figure 3.3(d) and (e) from the dissertation

Fig_spont_SITI2_osc_heatmaps.m plots Figure 3.3(a) from the dissertation

distance50.mat, distance_3D50.mat, distance_GC50.mat, fullNetwork50.mat, and glomeruli50.mat are data needed for running the simulations (generated from networkGenerator50.m). mitralCells50.mat contains the properties of the MCs in the network, but is not necessary for running simulations. granuleCells50 is too large to put on github, but is available at https://drive.google.com/file/d/1Yfc-Mfj-2HLsi6-qkEzOg0T0uF2mb2VF/view?usp=drive_link (and is also not necessary for running simulations).

granule.m and mitral.m define the granule cell and mitral cell classes respectively, and are as found at https://github.com/dkersen/olfactory-bulb (Kersen 2022).

networkGenerator50.m generates the network connectivity and is as found at https://github.com/dkersen/olfactory-bulb (Kersen 2022) with a few minor changes.
