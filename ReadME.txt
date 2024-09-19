Within_group
----Simulations (Simulation data using hotspot or streak patterns)
    1) Simulation data generation.R: Generate simulation data with hotpot and streak patterns, considering different parameters (Cell proportions, FC).
    2) Simulation_within_groups_hotspot_run.R: Analysis code for simulation data with hotspot pattern to identify SV genes using methods of SPARK, SPADE, MERINGUE, and DESpace_BayesSpace.
    3) Simulation_within_groups_streak_run.R: Analysis code for simulation data with streak pattern to identify SV genes using methods of SPARK, SPADE, MERINGUE, and DESpace_BayesSpace.
    4) stlearn_simu_run.py: Analysis code for simulation data to identify clusters using StLearn tool.
    5) DEspace_stlearn_simu_run.R: Analysis code for simulation data to identify SV genes using DEspace_StLearn approach, based on clusters identifed above.
    6) SpatialDE_simu_run.py: Analysis code for simulation data to identify SV genes using SpatialDE tool.
    7) SpaGCN_simu_run.py: Analysis code for simulation data to identify SV genes using SpaGCN tool.

----Simu_real_data (Simulation data using patterns from real data)
    ----SeqFISH
        1) Simulation_within_groups_SeqFISH_run.R:Analysis code for simulation data with pattern from SeqFISH data to identify SV genes using methods of SPARK, SPADE, MERINGUE, and DESpace_BayesSpace
        2) Simulation_within_groups_SeqFISH_run.swarm: Swarm file to run "Simulation_within_groups_SeqFISH_run.R"
        3) stlearn_simu_run.py: Analysis code for simulation data to identify clusters using StLearn tool.    
        4) DEspace_stlearn_simu_run.R: Analysis code for simulation data to identify SV genes using DEspace_StLearn approach, based on clusters identifed above.
        5) DEspace_stlearn_simu_run.swarm: Swarm file to run "DEspace_stlearn_simu_run.R".
        6) SpatialDE_simu_run.py: Analysis code for simulation data to identify SV genes using SpatialDE tool.
        7) SpaGCN_simu_run.py: Analysis code for simulation data to identify SV genes using SpaGCN tool.
    ----MERFISH
        1) Simulation_within_groups_MERFISH_run.R:Analysis code for simulation data with pattern from MERFISH data to identify SV genes using methods of SPARK, SPADE, MERINGUE, and DESpace_BayesSpace.
        2) Simulation_within_groups_MERFISH_run.swarm: Swarm file to run "Simulation_within_groups_MERFISH_run.R".
        3) stlearn_simu_run.py: Analysis code for simulation data to identify clusters using StLearn tool.    
        4) DEspace_stlearn_simu_run.R: Analysis code for simulation data to identify SV genes using DEspace_StLearn approach, based on clusters identifed above.
        5) DEspace_stlearn_simu_run.swarm: Swarm file to run "DEspace_stlearn_simu_run.R".
        6) SpatialDE_simu_run.py: Analysis code for simulation data to identify SV genes using SpatialDE tool.
        7) SpaGCN_simu_run.py: Analysis code for simulation data to identify SV genes using SpaGCN tool.
    ----MOB
        1) Simulation_within_groups_MOB_run.R:Analysis code for simulation data with pattern from MOB data to identify SV genes using methods of SPARK, SPADE, MERINGUE, and DESpace_BayesSpace.
        2) Simulation_within_groups_MOB_run.swarm: Swarm file to run "Simulation_within_groups_MOB_run.R".
        3) stlearn_simu_run.py: Analysis code for simulation data to identify clusters using StLearn tool.    
        4) DEspace_stlearn_simu_run.R: Analysis code for simulation data to identify SV genes using DEspace_StLearn approach, based on clusters identifed above.
        5) DEspace_stlearn_simu_run.swarm: Swarm file to run "DEspace_stlearn_simu_run.R" R code.
        6) SpatialDE_simu_run.py: Analysis code for simulation data to identify SV genes using SpatialDE tool.
        7) SpaGCN_simu_run.py: Analysis code for simulation data to identify SV genes using SpaGCN tool.

----real_data_single_group
    ----SeqFISH
        1) Real_within_groups_SeqFISH_run.R:Analysis code for SeqFISH real data to identify SV genes using methods of SPARK, SPADE, MERINGUE, and DESpace_BayesSpace.
        2) Real_within_groups_SeqFISH_runn.swarm: Swarm file to run "Real_within_groups_SeqFISH_run.R".
        3) stlearn_run.py: Analysis code for real data to identify clusters using StLearn tool.    
        4) DEspace_stlearn_run.R: Analysis code for real data  to identify SV genes using DEspace_StLearn approach, based on clusters identifed above.
        5) SpatialDE_run.py: Analysis code for real data to identify SV genes using SpatialDE tool.
        6) SpaGCN_run.py: Analysis code for real data to identify SV genes using SpaGCN tool.
    ----MOB			
        1) Real_within_groups_MOB_run.R:Analysis code for MOB real data to identify SV genes using methods of SPARK, SPADE, MERINGUE, and DESpace_BayesSpace.
        2) Real_within_groups_MOB_runn.swarm: Swarm file to run "Real_within_groups_MOB_run.R".
        3) stlearn_run.py: Analysis code for real data to identify clusters using StLearn tool.    
        4) DEspace_stlearn_run.R: Analysis code for real data  to identify SV genes using DEspace_StLearn approach, based on clusters identifed above.
        5) SpatialDE_run.py: Analysis code for real data to identify SV genes using SpatialDE tool.
        6) SpaGCN_run.py: Analysis code for real data to identify SV genes using SpaGCN tool.

Between_groups
----Simulations
    1) Simulation_hotspot_AUC.R: Analysis code to evaluate power of our method to identify SV genes between groups for hotspot pattern.
    2) Simulation_hotspot_AUC_diffCoor.R: Evaluating power of our method based on hotspot pattern with different coordinates.
    3) Simulation_hotspot_FPR_diffDirs.R: Evaluating FPR of our method based on hotspot patterns with different locatiosn and same coordinates.
    4) Simulation_streak_AUC.R: Analysis code to evaluate power of our method to identify SV genes between groups for streak pattern.
    5) Simulation_streak_AUC_diffCoor.R: Evaluating power of our method based on streak pattern with different coordinates.
    6) Simulation_streak_FPR_diffDirs.R: Evaluating FPR of our method based on streak patterns with different locatiosn and same coordinates
----Real_data_analysis
    1) Stereo-seq.R: Evaluating our method using a real spatial transcriptomic dataset with axolotl telencephalon (ARRISTA). Two different post injury stages (i.e., 2DPI, 5DPI) were included in the analysis.