# single-probe
 single-probe data processing

### filename tips
 spa_corr_spac_v2/3 -> not use axis method simulation spatial correlation using SPAC
 spa_corr_grid_spac_v2 -> add phase estimate based axis method simulation using SPAC
 spa_corr_grid_mpac_v2 -> based "spa_corr_grid_spac_v2" change to use MPAC
 spa_corr_grid_mpac_v2_1 -> add simulation h(channel parameter), and simulation the spatial correlation without signal, just using channel itself.
 
 Above files use two antennas moving from center point
 
 spa_corr_grid_mpac_v3 -> two antennas change to single ref ant
 spa_corr_grid_mpac_v3_1 -> add single-cluster model
 spa_corr_grid_mpac_v4 -> add CE influnce on signal. If it works, this h parameter would implemented on VST CE.