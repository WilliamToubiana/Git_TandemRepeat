#############################
# Mapping 250kb Tps windows #
#############################

##Tandem repeats
setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/whole_genome_alignment")

data<-read.table("Tps-genome_250kb_to_Tcm_genome_minimap2_s40.paf", header=F, sep="\t", quote="", fill = TRUE)
data<-read.table("Tps-genome_250kb_to_Tcm_genome_minimap2_s2500.paf", header=F, sep="\t", quote="", fill = TRUE)
colnames(data) <- c("query", "query_length", "query_start", "query_end", "strand", "target", "target_length", "target_start", "target_end",
                    "number_matches", "alignment_block_length","mapping_quality","edit_distance","matching_score","alignment_score","number_ambiguous_bases",
                    "type_alignment","V18","V19","V20","V21","V22","V23","V24")

data$cov<-abs(data$query_start-data$query_end)

library(stringr)
data[c('AS', 'i', 'alignment_score2')] <- str_split_fixed(data$alignment_score, ':', 3)
data[c('NM', 'i', 'edit_distance2')] <- str_split_fixed(data$edit_distance, ':', 3)
data[c('ms', 'i', 'matching_score2')] <- str_split_fixed(data$matching_score, ':', 3)
summary(data)

data$alignment_score2<-as.numeric(data$alignment_score2)
data$edit_distance2<-as.numeric(data$edit_distance2)
data$matching_score2<-as.numeric(data$matching_score2)

data$seq_id<-(data$number_matches)/(data$number_matches+data$edit_distance2) #values of edit_distance2 include 1)nb of mismatches, 2)nb of insertions, 3)nb of deletions
summary(data)

##Optional filtering based on mapping quality

data_sub<-subset(data, mapping_quality>59)

##Summarize data per 250-kb regions (replace by data_sub if used before)
library(dplyr)
window_summary <- data %>%
  group_by(query) %>%
  summarise(
    alignment_score_avg = mean(alignment_score2),
    edit_distance_avg = mean(edit_distance2),
    matching_score_avg = mean(matching_score2),
    seq_id_avg = mean(seq_id),
  )

summary(window_summary)

window_summary[c('chr', 'start_end')] <- str_split_fixed(window_summary$query, ':', 2)
window_summary[c('window_start', 'window_end')] <- str_split_fixed(window_summary$start_end, '-', 2)
window_summary$window_start<-as.numeric(window_summary$window_start)
window_summary$window_start2<-window_summary$window_start+1
window_summary$chm_pos<-paste(window_summary$chr, window_summary$window_start2, sep="-")



#########################
### ChIP cenh3 testes ###
#########################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
cenh3<-read.table("Tps_cenh3_testes_1_GW_coverage_DR_250kb.txt",
                  header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w250000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_cenh3_R1/Tps_testes_cenh3_1_coverage_DR_100kb.txt

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/ChIP/cenh3")
input1<-read.table("Tps_input_testes_1_GW_coverage_DR_250kb.txt",
                   header=FALSE) #bedtools coverage -a genomes/Tps_chm_size_mtDNAv350_w250000.bed -b /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/mapping/Tps_testes_input_R1/Tps_testes_input1_bwa_final_DR.bam -sorted -g genomes/Tps_LRv5b_mtDNAv350.fasta.fai -mean > /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D2c/BICC/tracks/Tps_testes_input_R1/Tps_testes_input1_coverage_DR_100kb.txt

data1<-cbind(cenh3,input1[4])
colnames(data1)<-c("Scaffold_name", "start", "stop","coverage_cenh3", "coverage_input")

data1$coverage_cenh3_norm<-data1$coverage_cenh3/66295120 #normalized by the number of mapped reads (chip_cenh3_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$coverage_input_norm<-data1$coverage_input/53252668 #normalized by the number of mapped reads (chip_input_tdi_testes_1_bwa_final_DR_flagstat_out.txt)
data1$ratiocenH3_norm<-data1$coverage_cenh3_norm/data1$coverage_input_norm
data1$log2cenH3_norm<-log2(data1$ratiocenH3_norm)
data1$start2<-data1$start+1
data1$chm_pos<-paste(data1$Scaffold_name, data1$start2, sep="-")
data1$log2cenH3_norm[!is.finite(data1$log2cenH3_norm)] <- 0




#######################################
### TR Proportion along chromosomes ###
#######################################

setwd("/Users/wtoubian/Desktop/Tandem_repeats/Timema_evolution/TR_proportion_Tps")

Tps_cov_genome<-read.table("Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_GW_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_GW_sum_coverage_250kb.txt
colnames(Tps_cov_genome)<-c("chr", "start", "end", "sum_cov_GW")
Tps_cov_TR<-read.table("Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt", header=F, sep=" ", quote="") #This was generated from the command line: awk 'NR==FNR {a[$1, int($2/250000)] += $3; next} {print $1, $2, $3, a[$1, int($2/250000)]}' /nas/FAC/FBM/DEE/tschwand/asex_sinergia/D1c/wtoubian/minimap2_timema_reseq/Tps-filtered_indG_minimap2_TR_coverage.txt genomes/Tps_chm_size_mtDNAv350_w250000.bed > Tps-filtered_indG_minimap2_TR_sum_coverage_250kb.txt
colnames(Tps_cov_TR)<-c("chr", "start", "end", "sum_cov_TR")

Tps_cov_genome$start2<-Tps_cov_genome$start+1
Tps_cov_TR$start2<-Tps_cov_TR$start+1
Tps_cov_genome$chm_pos<-paste(Tps_cov_genome$chr, Tps_cov_genome$start2, sep="-")
Tps_cov_TR$chm_pos<-paste(Tps_cov_TR$chr, Tps_cov_TR$start2, sep="-")

mergeTRprop<-cbind(Tps_cov_genome, Tps_cov_TR[4])
summary(mergeTRprop)
mergeTRprop[is.na(mergeTRprop)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation
mergeTRprop$TRprop<-mergeTRprop$sum_cov_TR/mergeTRprop$sum_cov_GW
summary(mergeTRprop)
mergeTRprop[is.na(mergeTRprop)] <- 0 #NAs are retrieved when windows have no mapped reads or no TR annotation
summary(mergeTRprop)



######################################################
### Alignment centromere vs non-centromere regions ###
######################################################

##Combine TR prop and shared repeat datasets
combined<-merge(mergeTRprop, window_summary[c(2,3,4,5,11)], by="chm_pos", all.x=T)
summary(combined)

##Subdivide windows into centromere vs non-centromere based on log2 ratios
data1$centromere<-ifelse(data1$log2cenH3_norm>0.25, "centromere", "non-centromere")

##Combine TR prop, centromere and shared repeat datasets
combined2<-merge(combined, data1[c(9,11,12)], by="chm_pos", all.x=T)
summary(combined2)





##Scatterplot and density plots with TR proportion
##################################################
combined2_subset<-combined2
combined2_subset2<-subset(combined2_subset, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                  | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                  | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")
summary(combined2_subset2)

scatterPlot <- ggplot() +
  geom_point(data = combined2_subset2 %>% filter(centromere == "non-centromere"),
             aes(TRprop, seq_id_avg),
             color = "grey", alpha = 0.4) +
  geom_point(data = combined2_subset2 %>% filter(centromere == "centromere"),
             aes(TRprop, seq_id_avg),
             color = "red", alpha = 0.6) +
  geom_smooth(data = combined2_subset2,
              aes(TRprop, seq_id_avg, color = centromere),
              method = "lm", se = TRUE) +
  scale_color_manual(values = c("centromere" = "red",
                                "non-centromere" = "grey")) +
  theme_bw() +
  theme(legend.position = c(0,1),
        legend.justification = c(-5.5,-1)) +
  ylab("Sequence identity")

# Marginal density plot of x (top panel)
xdensity <- ggplot(combined2_subset2, aes(TRprop, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('red','black')) + 
  theme_bw() +
  theme(legend.position = "none")

# Marginal density plot of y (right panel)
ydensity <- ggplot(combined2_subset2, aes(seq_id_avg, fill=centromere)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('red','black')) + 
  theme_bw() +
  theme(legend.position = "none")+coord_flip()+
  xlab("Sequence identity")


blankPlot <- ggplot() + 
  geom_blank(aes(1, 1)) +
  theme_void()


library("gridExtra")
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))


#Stats
library(npANCOVA)
combined2_subset2<-na.omit(combined2_subset2)
results<-McSweeny_Porter(combined2_subset2, seq_id_avg~TRprop+centromere)
print(results$regression_equation_interaction)
print(results$regression_equation_covariate_factor)
print(results$group_effect)
print(results$interaction_effect)


#####################################################
## Distribution alignment scores along chromosomes ##
#####################################################

combined3<-subset(combined2, chr=="Tps_LRv5b_scf1" | chr=="Tps_LRv5b_scf2" | chr=="Tps_LRv5b_scf3" | chr=="Tps_LRv5b_scf4" | chr=="Tps_LRv5b_scf5"
                                | chr=="Tps_LRv5b_scf6" | chr=="Tps_LRv5b_scf7" | chr=="Tps_LRv5b_scf8" | chr=="Tps_LRv5b_scf9" | chr=="Tps_LRv5b_scf10"
                                | chr=="Tps_LRv5b_scf11" | chr=="Tps_LRv5b_scf12")

combined3$scaff_number<-as.numeric(as.character(gsub("^.*scf","", combined3$chr)))

combined3<-combined3[order(combined3$scaff_number, combined3$start2),]
combined3$cumul_start<-seq(0, by = 250000, length.out = nrow(combined3))


# Calculate scaffold boundaries (start and end positions)
library(dplyr)
scaffold_bounds <- combined3 %>%
  group_by(scaff_number) %>%
  summarise(min_pos = min(cumul_start),
            max_pos = max(cumul_start)) %>%
  mutate(fill_col = rep(c("white", "grey90"), length.out = n()))


library(ggplot2)

alignment<-ggplot(combined3, aes(x=cumul_start, y=seq_id_avg, color = centromere))+
  # Add alternating shaded backgrounds
  geom_rect(data = scaffold_bounds,
            aes(xmin = min_pos, xmax = max_pos, ymin = -Inf, ymax = Inf, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.6) +
  scale_fill_identity() +
  geom_point(size=0.5, alpha = 0.6) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "grey"))+
  ylab("Sequence identity\nmatches/(matches + mismatches and gaps)")+
  xlab("Genomic coordinates")+
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  )+ 
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))



cenh3_testes<-ggplot(combined3, 
                     aes(x = cumul_start, y = log2cenH3_norm, color = centromere)) +
  # Add alternating shaded backgrounds
  geom_rect(data = scaffold_bounds,
            aes(xmin = min_pos, xmax = max_pos, ymin = -Inf, ymax = Inf, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.6) +
  scale_fill_identity() +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "grey")) +
  geom_smooth(aes(group = scaff_number), size = 0.5, color = "white",
              method = "gam", formula = y ~ s(x, bs = "cs", k = 30)) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "black") +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  ylab("log2(CenH3/Input)") +
  xlab("Genomic coordinates") +
  ylim(-1, 3)


TR_array<-ggplot(combined3, aes(x=cumul_start, y=TRprop, color = centromere))+
  # Add alternating shaded backgrounds
  geom_rect(data = scaffold_bounds,
            aes(xmin = min_pos, xmax = max_pos, ymin = -Inf, ymax = Inf, fill = fill_col),
            inherit.aes = FALSE, alpha = 0.4) +
  scale_fill_identity() +
  geom_point(size=0.5) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "grey")) +
  ylab("TR proportion")+
  xlab("Genomic coordinates")+
  geom_hline(yintercept = 0.125, linetype = "dashed", color = "black") +
  theme(
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "none"
  )+
  geom_smooth(aes(group = scaff_number), size=0.5, colour = "white", method="gam", formula = y~s(x, bs="cs", k=30))


# Plots TRprop vs cenh3 vs avg_shared
library(cowplot)

theme_set(theme_minimal())

plot_grid(
  plot_grid(
    TR_array + theme(legend.position = "none")
    , cenh3_testes
    , alignment + theme(legend.position = "none")
    , ncol = 1
    , align = "hv")
  , plot_grid(
    get_legend(cenh3_testes)
    , ggplot()
    , get_legend(alignment)
    , ncol =1)
  , rel_widths = c(12,0)
)





###################################################################################
### Subsampling "non-centromere" windows to match distribution of TR proportion ###
###################################################################################

library(dplyr)

#set.seed(123)

# Split data
cen  <- combined2_subset2 %>% filter(centromere == "centromere")
noncen <- combined2_subset2 %>% filter(centromere == "non-centromere")

# Create TRprop bins
n_bins <- 12   # adjust for smoother matching
breaks <- quantile(combined2_subset2$TRprop, probs = seq(0, 1, length.out = n_bins + 1))

cen$bin <- cut(cen$TRprop, breaks = breaks, include.lowest = TRUE)
noncen$bin <- cut(noncen$TRprop, breaks = breaks, include.lowest = TRUE)

# Count how many centromere points in each bin
cen_bin_counts <- table(cen$bin)

# Sample non-centromere points to match counts
noncen_sub <- noncen %>%
  group_by(bin) %>%
  sample_n(size = cen_bin_counts[bin][1], replace = FALSE) %>%
  ungroup()

# Combine subsampled non-centromeres with all centromeres
matched <- bind_rows(cen, noncen_sub)


#Plots
scatterPlot <- ggplot(matched, aes(TRprop, seq_id_avg, color = centromere)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c('red','grey')) +
  geom_smooth(method = lm) +
  theme(legend.position = c(0,1), legend.justification = c(-5.5, -1)) +
  ylab("Sequence identity")

xdensity <- ggplot(matched, aes(TRprop, fill = centromere)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = c('red','black')) +
  theme(legend.position = "none")

ydensity <- ggplot(matched, aes(seq_id_avg, fill = centromere)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = c('red','black')) +
  theme(legend.position = "none") +
  coord_flip() +
  xlab("Sequence identity")

grid.arrange(
  xdensity, blankPlot,
  scatterPlot, ydensity,
  ncol = 2, nrow = 2,
  widths = c(4, 1.4),
  heights = c(1.4, 4)
)



##################################################################################################################
### Iterations and Subsampling of "non-centromere" windows to match "centromere" distribution of TR proportion ###
##################################################################################################################

library(dplyr)
library(broom)
library(ggplot2)

#set.seed(42)    # reproducible overall

# Parameters
n_iter <- 1000        # number of random subsamples (change to 1000 if you want)
n_bins  <- 12        # bins for TRprop quantile-matching (as before)

# Split data
cen    <- combined2_subset2 %>% filter(centromere == "centromere")
noncen <- combined2_subset2 %>% filter(centromere == "non-centromere")

# Precompute quantile breaks (same bins used every iteration)
breaks <- quantile(combined2_subset2$TRprop, probs = seq(0, 1, length.out = n_bins + 1),
                   na.rm = TRUE, type = 7)

# Assign bins to both sets (keeps binning constant across iterations)
cen <- cen %>% mutate(bin = cut(TRprop, breaks = breaks, include.lowest = TRUE))
noncen <- noncen %>% mutate(bin = cut(TRprop, breaks = breaks, include.lowest = TRUE))

# How many to sample per bin (target = counts in centromere)
target_per_bin <- table(cen$bin)

# quick helper to sample from noncen for a given bin counts (handles small bins)
sample_noncen_for_bin <- function(bin_name, n_target) {
  sub <- noncen %>% filter(bin == bin_name)
  n_avail <- nrow(sub)
  if (n_avail == 0) {
    # no available non-centromere in this bin: return 0 rows
    return(sub[0, ])
  }
  # if fewer available than needed, sample with replacement and warn once (but do not stop)
  replace_flag <- ifelse(n_avail < n_target, TRUE, FALSE)
  sub %>% sample_n(size = n_target, replace = replace_flag)
}

# Precompute centromere model (doesn't change each iteration)
lm_cen <- lm(seq_id_avg ~ TRprop, data = cen)
lm_cen_sum <- summary(lm_cen)
cen_stats <- list(
  slope = coef(lm_cen)["TRprop"],
  intercept = coef(lm_cen)["(Intercept)"],
  r2 = lm_cen_sum$r.squared,
  mean_seq_id = mean(cen$seq_id_avg, na.rm = TRUE),
  n = nrow(cen)
)

# Container for iteration results
results <- vector("list", n_iter)

# Main loop: iterate subsampling + fit model on sampled noncen
for (i in seq_len(n_iter)) {
  # build sampled non-centromere dataframe by sampling within each bin
  sampled_list <- mapply(
    sample_noncen_for_bin,
    names(target_per_bin),
    as.integer(target_per_bin),
    SIMPLIFY = FALSE
  )
  noncen_sub <- bind_rows(sampled_list)
  
  # In case some bins had zero available noncen, the total may be < nrow(cen).
  # You can choose to skip this iteration or proceed — here we proceed and record the actual n.
  # Fit lm on sampled noncentromeres (if enough points)
  if (nrow(noncen_sub) >= 3) {
    lm_non <- lm(seq_id_avg ~ TRprop, data = noncen_sub)
    s <- summary(lm_non)
    non_stats <- list(
      slope = coef(lm_non)["TRprop"],
      intercept = coef(lm_non)["(Intercept)"],
      r2 = s$r.squared,
      mean_seq_id = mean(noncen_sub$seq_id_avg, na.rm = TRUE),
      n = nrow(noncen_sub)
    )
  } else {
    non_stats <- list(slope = NA, intercept = NA, r2 = NA, mean_seq_id = NA, n = nrow(noncen_sub))
  }
  
  results[[i]] <- tibble(
    iter = i,
    cen_slope = cen_stats$slope,
    cen_intercept = cen_stats$intercept,
    cen_r2 = cen_stats$r2,
    cen_mean_seq_id = cen_stats$mean_seq_id,
    cen_n = cen_stats$n,
    non_slope = non_stats$slope,
    non_intercept = non_stats$intercept,
    non_r2 = non_stats$r2,
    non_mean_seq_id = non_stats$mean_seq_id,
    non_n = non_stats$n,
    slope_diff = non_stats$slope - cen_stats$slope
  )
}

res_df <- bind_rows(results)  
summary(res_df)

# Summaries
summary_stats <- res_df %>%
  summarize(
    non_slope_mean = mean(non_slope, na.rm = TRUE),
    non_slope_sd = sd(non_slope, na.rm = TRUE),
    slope_diff_mean = mean(slope_diff, na.rm = TRUE),
    slope_diff_sd = sd(slope_diff, na.rm = TRUE),
    prop_non_less_than_cen = mean(non_slope < cen_slope, na.rm = TRUE),
    median_non_n = median(non_n, na.rm = TRUE)
  )

print(summary_stats)

# Quick diagnostics plots
p1 <- ggplot(res_df, aes(x = non_slope)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = cen_stats$slope, color = "red", linetype = "dashed") +
  ggtitle("Distribution of non-centromere slopes across iterations\n(vertical dashed = centromere slope)") +
  xlab("slope (seq_id_avg ~ TRprop)")

p2 <- ggplot(res_df, aes(x = slope_diff)) +
  geom_density() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  ggtitle("Distribution of slope differences (non - cen)") +
  xlab("slope difference")

p3 <- ggplot(res_df, aes(x = non_mean_seq_id)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = cen_stats$mean_seq_id, color = "red", linetype = "dashed") +
  ggtitle("Distribution of mean(seq_id_avg) for sampled non-centromere sets") +
  xlab("mean(seq_id_avg)")

# Print plots
print(p1); print(p2); print(p3)


# Choose a representative iteration (e.g., iteration 1)
set.seed(123)
sampled_list_rep <- mapply(
  sample_noncen_for_bin,
  names(target_per_bin),
  as.integer(target_per_bin),
  SIMPLIFY = FALSE
)
noncen_rep <- bind_rows(sampled_list_rep)

# Combine cen and the representative noncen
plot_df <- bind_rows(
  cen %>% mutate(group = "centromere"),
  noncen_rep %>% mutate(group = "non-centromere")
)

#Plot of a representative iteration
xlim_range <- range(plot_df$TRprop, na.rm = TRUE)

scatterPlot_rep <- ggplot(plot_df, aes(TRprop, seq_id_avg, color = group)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("centromere" = "red", "non-centromere" = "grey")) +
  scale_x_continuous(limits = xlim_range, expand = c(0,0)) +
  geom_smooth(method = "lm", se = FALSE) +
  ylab("Sequence identity") +
  theme_bw() +
  theme(legend.position = c(0,1), legend.justification = c(-0.2,1.2))

xdensity_rep <- ggplot(plot_df, aes(TRprop, fill = group)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = c("centromere" = "red", "non-centromere" = "black")) +
  scale_x_continuous(limits = xlim_range, expand = c(0,0)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

ydensity_rep <- ggplot(plot_df, aes(seq_id_avg, fill = group)) +
  geom_density(alpha = .5) +
  scale_fill_manual(values = c("centromere" = "red", "non-centromere" = "black")) +
  coord_flip() +
  xlab("Sequence identity") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

blankPlot <- ggplot() + 
  geom_blank(aes(1, 1)) +
  theme_void()

library(gridExtra)

grid.arrange(
  xdensity_rep, blankPlot,
  scatterPlot_rep, ydensity_rep,
  ncol = 2, nrow = 2,
  widths = c(4, 1.4),
  heights = c(1.4, 4)
)



