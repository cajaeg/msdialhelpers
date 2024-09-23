##
## TBD
##
# .mslib <-
#   readMSP(.mslib_used_file) %>%
#   as_tibble() %>%
#   mutate(PEAKS = gsub("([0-9]+) ([0-9]+)", "\\1:\\2", PEAKS)) %>%
#   rename(met_name = Name,
#          lib.ri = RI,
#          lib.spectrum = PEAKS,
#          met_class = CLASS
#   ) %>%
#   filter(!duplicated(met_name)) %>%
#   select(met_name, KEGG, lib.ri, lib.spectrum, met_class)
#
# lib <- readMSP
#
#
# .full <-
#   ## load full results from MSDIAL (from which further tables are produced below)
#   plyr::ldply(.msd_results_files, msdialhelpers::loadAlignmentResults4) %>%
#   as_tibble() %>%
#   rename(met_name = name,
#          ri = Average.RI
#   ) %>%
#   ## restrict features to identified and non-contaminant ones
#   filter(str_starts(met_name, "Unknown", negate = T)) %>%
#   mutate(neutral_mass = round(check_chemform(isotopes, form)[,"monoisotopic_mass"], 5)) %>%
#   mutate(nm_group = interaction(cutree(hclust(dist(c(-10000, neutral_mass))), h = 0.01)[-1]-1,
#                                 cutree(hclust(dist(c(-10000, rt))), h = 2)[-1]-1,
#                                 drop=T)) %>%
#   group_by(nm_group) %>% # remove duplicate annotations of same spectrum based on better dot product
#   filter(row_number(desc(total_score)) == 1) %>%
#   group_by(met_name) %>% # remove redundant spectral matches (same metabolite assigned to second and possibly more peaks)
#   filter(row_number(desc(total_score)) == 1) %>%
#   ungroup() %>%
#   left_join(.mslib, by = "met_name")
