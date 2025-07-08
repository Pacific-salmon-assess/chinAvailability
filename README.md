# Estimates of Chinook salmon composition in SRKW prey remains and critical habitat

Data and code associated with estimates of relative abundance of Chinook salmon stocks and size classes in SRKW prey remains and critical habitat. Associated with CSAS Research Document **Southern Resident Killer Whale prey selectivity in relation to Chinook Salmon stock and size composition within Canadian critical habitat**.

Scripts:
1) functions/: utility functions
2) clean_data.R: cleaning script used to preprocess data (not runnable with publicly available data) 
2) mvtweedie_fit.R: estimates of spatiotemporal variability in size composition
3) mvtweedie_fit_size.R: estimates of spatiotemporal variability in size composition
4) mvtweedie_fit_sdmTMB.R: estimates of spatiotemporal variability in size composition using mvtweedie analogue but fit in sdmTMB (unused)
5) rkw_diet_exp.R: figures associated with SRKW prey remains
6) sep_figs.R: figures associated with pHOS and PBT marking rates
7) sampling_maps.R: figures representing sampling locations
8) size_by_stock.R: size-at-age of Chinook salmon stocks
9) esc_index.R: generates figures showing terminal abundance
10) aging_error.R: code associated with calculating aging error in supplementary material
11) selectivity_sim.R: conducts selectivity analyses
12) selectivity_sim_sensitivity.R: conducts selectivity analyses in sensitivity analysis