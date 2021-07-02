### Variant calling pipeline

These scripts are included to aid in the reproducibility of the variant calls used in subsequent analyses. While all commands are included they are not currently set up to be easily run on a new system.

A brief description of each stage of the pipeline adn the accompanying R-scripts is given below.


1) Alignment was performed using the ARTIC Illumina nextflow pipeline available from https://github.com/connor-lab/ncov2019-artic-nf. The pipeline trims primer and adapter sequences prior to alignment to the reference genome using bwa (MN908947.3). The resulting bam files were then used in the subsequent steps.

2) bam2R was used to count the number of reads supporting each base type at each site for each sequencing run for each sample using the scripts
    - heron_genome_wide_pileup.R
    - heron_genome_wide_pileup_wrapper.R

3) The individual pile-up was then output into a 3D array for all samples. This produced one combined table for each replicate using the script
    - heron_combined_coverage_table.R

4) Coverage statistics were calculated and used to identifying the 521 “good” samples from which the background panel samples were selected.
    - heron_sample_amplicon_coverage_exploration.R

5) The application of the selection criteria (outlined in the [methods](https://www.biorxiv.org/content/10.1101/2020.12.23.424229v1)) was used to generate the list of samples to be used in the background panel
    - heron_good_sample_and_normal_panel_list_creation.R

6) A background paned was generated using the lust of selected samples. Sites with ≥1% VAF for non-ref bases were set as uninformative.
    - heron_background_panel_generation.R

7) The rho and mu values were calculated for the background panel.
    - heron_shearwater_normal_panel_rho_mu_calculation.R

8) The Shearwater p-values were calculated using background panel with pre-calculated rho and mu. We set a lower bound value of rho equal to 10^-3 to be conservative.
    - heron_shearwater_al28_lower_bound_precalculated_rho_farm_script.R
    - heron_shearwater_al28_lower_bound_precalculated_rho_farm_wrapper.R

9) The per sample Shearwater results were combined into a single .RDS object. 
    - heron_merge_shearwater_results.R

10) q-values were then calculated and consecutive variants were merged. We used a difference in VAF of ≤ 0.05 to merge variants.
    - heron_post_shearwater_filter_all_samples_subs_and_indels_min_depth_100x.R

11) Approaches for combining calls from two replicates. The final method used an intermediate approach requiring a qval < 0.05 in one replicate and a pval < 0.01 in the other.
    - heron_shearwater_panel_1_trinucleotide_plot_function.R

12) Calculation of read support for each position in MNVs and annotation of support from both replicates.
    - heron_shearwater_strict_union_support_from_each_read.R