## Figure 1

Irene Stevens

23/05/2023

# This script will generate Figure 1: Compare S. Cerevisae (BY4741) and C. albicans (sc5314) ribosome dynamics

## Load libraries

``` r
library(ggplot2)
library(gridExtra)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ lubridate 1.9.2     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.2     ✔ tidyr     1.3.0
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::combine() masks gridExtra::combine()
    ## ✖ dplyr::filter()  masks stats::filter()
    ## ✖ dplyr::lag()     masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggrepel)
```

## Figure 1A: Frame preferences in BY4741 versus SC5314.

Values represent the average of replicates from frame_stats.txt files
(see Data folder). For more details, see ‘Figure1’ R script notes in the
Figures folder.

# S. Cerevisae (BY4741) Control

Values represent the average of replicates of values from files:
y4741-ypd_ctr-r1_s1_frame_stats.txt, y4741-ypd_ctr-r2_s1_frame_stats.txt

``` r
data1 <- data.frame(frequency = c(30.75, 42.3, 27), group = c("F0", "F1", "F2"))
plot1 <- ggplot(data1, aes(x = group, y = frequency, fill = group)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = c("F0" = rgb(255, 221, 214, maxColorValue = 255), 
                                 "F1" = rgb(255, 221, 214, maxColorValue = 255), 
                                 "F2" = rgb(255, 221, 214, maxColorValue = 255))) +
    labs(x = "Control", y = "% Reads") +
    theme_classic() +
    ylim(0, 60)  # Set the same y-axis limits as plot1
plot1 <- plot1 + theme(text = element_text(size = 20), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 25))
```

## S. Cerevisae (BY4741) Cycloheximide

Values represent the average of replicates of values from files:
by4741-chx_t05-r1_s1_frame_stats.txt,
by4741-chx_t05-r2_s1_frame_stats.txt

``` r
data2 <- data.frame(frequency = c(23.84, 56.61, 19.55), group = c("F0", "F1", "F2"))
plot2 <- ggplot(data2, aes(x = group, y = frequency, fill = group)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = c("F0" = rgb(199, 236, 220, maxColorValue = 255), 
                                 "F1" = rgb(199, 236, 220, maxColorValue = 255), 
                                 "F2" = rgb(199, 236, 220, maxColorValue = 255))) +
    labs(x = "CHX", y = "% Reads", title="S. Cerevisae (BY4741)") +
    theme_classic() +
    ylim(0, 60)  # Set the same y-axis limits as plot1
plot2 <- plot2 + theme(text = element_text(size = 20), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 25))

# Combine the S. Cerevisae plots in a panel
grid.arrange(plot1, plot2, ncol = 2)
```

![](unnamed-chunk-3-1.png)<!-- -->

## Candida albicans (sc5314) Control

Values represent the average of replicates of values from files
CTX1_S21_R1_001_frame_stats.txt, CTX3_S22_R1_001_frame_stats.txt

``` r
data3 <- data.frame(frequency = c(41.6, 34.2, 24.7), group = c("F0", "F1", "F2"))
plot3 <- ggplot(data3, aes(x = group, y = frequency, fill = group)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = c("F0" = rgb(255, 135, 111, maxColorValue = 255), 
                                 "F1" = rgb(255, 135, 111, maxColorValue = 255), 
                                 "F2" = rgb(255, 135, 111, maxColorValue = 255))) +
    labs(x = "Control", y = "% Reads", title= "C.albicans (SC5314)") +
    theme_classic() + ylim(0,60)
plot3 <- plot3 + theme(text = element_text(size = 20), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 25))
```

## Candida albicans (sc5314) Cycloheximide

Values represent the average of replicates of values from files:
CHX1_S13_R1_001_frame_stats.txt, CHX3_S15_R1_001_frame_stats.txt

``` r
data4 <- data.frame(frequency = c(38.83, 33.92, 27.24), group = c("F0", "F1", "F2"))
plot4 <- ggplot(data4, aes(x = group, y = frequency, fill = group)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = c("F0" = rgb(77, 181, 188, maxColorValue = 255), 
                                 "F1" = rgb(77, 181, 188, maxColorValue = 255), 
                                 "F2" = rgb(77, 181, 188, maxColorValue = 255))) +
    labs(x = "CHX", y = "% Reads") +
    theme_classic() + ylim(0,60)
plot4 <- plot4 + theme(text = element_text(size = 20), 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 25))

#Combine the Candida plots in a panel
grid.arrange(plot3, plot4, ncol = 2)
```

![](unnamed-chunk-5-1.png)<!-- -->

## Figure 1B: Metaplots of termination site in BY4741 versus SC5314.

Counts per Million (CPM) = raw read count / total reads x 10^6 Raw read
count files: see Data folder, meta_counts_TERM.txt files Total read
count information: see Data folder, data_summary.txt files
(NumOfMapPositions) The average of 3 biological replicates is plotted

## BY4741 Control

``` r
control <- read.table("/Users/vip/Documents/CUG_reprocessed/Figure 1/A./BY4741/BY4741_Control_-60to+10.txt", header = TRUE)
ggplot(control, aes(x = Position, y = Avg_replicates)) +
  geom_line(color = rgb(255,221,214, maxColorValue = 255)) +  
  geom_vline(xintercept = -17, linetype = "dashed", color = "black") +   
  labs(x = "TERM", y = "CPM", title ="BY4741 Control") +
  theme_classic() + 
  scale_x_continuous(breaks = seq(min(control$Position), max(control$Position), by = 5))
```

![](unnamed-chunk-6-1.png)<!-- -->

## BY4741 Cycloheximide

``` r
cycloheximide <- read.table("/Users/vip/Documents/CUG_reprocessed/Figure 1/A./BY4741/BY4741_Cycloheximide_-60to+10.txt", header = TRUE)
ggplot(cycloheximide, aes(x = Position, y = average_replicates)) +
  geom_line(color = rgb(199,236,220, maxColorValue = 255)) +  
  geom_vline(xintercept = -17, linetype = "dashed", color = "black") +   
  labs(x = "TERM", y = "CPM", title ="BY4741 + Cycloheximide") +
  theme_classic() + 
  scale_x_continuous(breaks = seq(min(control$Position), max(control$Position), by = 5))
```

![](unnamed-chunk-7-1.png)<!-- -->

## SC5314 Control

``` r
Candida_control <- read.table("/Users/vip/Documents/CUG_reprocessed/Figure 1/A./Candida_Control_-60to+10.txt", header = TRUE)
ggplot(Candida_control, aes(x = Position, y = sc5314_control_CPM)) +
  geom_line(color = rgb(255,135,111, maxColorValue = 255)) +  geom_vline(xintercept = -18, linetype = "dashed", color = "black")  +   labs(x = "TERM", y = "CPM", title ="SC5314 Control") +
  theme_classic() + 
  scale_x_continuous(breaks = seq(min(control$Position), max(control$Position), by = 5))
```

![](unnamed-chunk-8-1.png)<!-- -->

## SC5314 Cycloheximide

``` r
Candida_cycloheximide <- read.table("/Users/vip/Documents/CUG_reprocessed/Figure 1/A./Candida_Cycloheximide_-60to+10.txt", header = TRUE)
ggplot(Candida_cycloheximide, aes(x = Position, y = sc5314_Cycloheximide)) +
  geom_line(color = rgb(77,181,188, maxColorValue = 255)) +  geom_vline(xintercept = -18, linetype = "dashed", color = "black")  +   labs(x = "TERM", y = "CPM", title ="SC5314 + Cycloheximide") +
  theme_classic()+ scale_x_continuous(breaks = seq(min(control$Position), max(control$Position), by = 5))
```

![](unnamed-chunk-9-1.png)<!-- -->

## Figure 1C. Metaplot at TSS in SC53144 treated with Cycloheximide

Input file was calculated as follows: Counts per million = raw counts
(meta_counts_TERM.txt)/total reads (NumOfMapPositions from
data_summary.txt) x 10^6 Plotted are Average of Cycloheximide_r1_CPM,
Cycloheximide_r2_CPM cor Treatment group and Average of Control_r1_CPM,
Control_r3_CPM for Control group

``` r
tss_start <- read.table("/Users/vip/Documents/CUG_reprocessed/Figure 1/START_metaplot_Candida.txt", header = TRUE)

ggplot(tss_start, aes(x = Position)) +
  geom_line(aes(y = Candida_Cycloheximide), color = rgb(77, 181, 188, maxColorValue = 255)) +
  geom_line(aes(y = Candida_Control), color = rgb(255, 135, 111, maxColorValue = 255)) +
  geom_vline(xintercept = -14, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = +4, linetype = "dashed", color = "grey") +
  annotate("text", x = -14, y = max(tss_start$Candida_Cycloheximide), label = "-14", vjust = -1) +
  annotate("text", x = 0, y = max(tss_start$Candida_Cycloheximide), label = "0", vjust = -1) +
  annotate("text", x = 4, y = max(tss_start$Candida_Cycloheximide), label = "+4", vjust = -1) +
  labs(x = "TSS", y = "CPM")+
  theme_classic() 
```

![](unnamed-chunk-10-1.png)<!-- -->

## Figure 1D. Fast Fourier Transform of SC5314 treated with Cycloheximide: see Candida_main.html report

\##Figure 1E. See Candida_codon_linecharts.html report

\##Figure 1F. Scatterplot of amino acid pauses in depleted versus
control conditions

Note: the values at -18 from amino_acid_pauses.txt

``` r
data <- read.table("/Users/vip/Documents/CUG_reprocessed/Figure 1/aa_pauses_ypd_-18.txt", header = TRUE, row.names = 1)

data_avg <- data %>%
  mutate(ProArgDepletion_avg = rowMeans(select(., starts_with("proline")), na.rm = TRUE),
         Control_avg = rowMeans(select(., starts_with("control")), na.rm = TRUE))

ggplot(data_avg, aes(x = Control_avg, y = ProArgDepletion_avg, label = rownames(data_avg))) +
  geom_point(size = 3, color = ifelse(rownames(data_avg) == "ARG", "red", "black")) +
  geom_text_repel(
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey',
    segment.size = 0.5,
    aes(label = rownames(data_avg)),
    color = "black"
  ) +
  labs(title = "Position -18",
       x = "Control (5Pseq Reads)",
       y = "Arg-Pro Depleted (5PSeq Reads)") +
  theme_minimal()
```

![](unnamed-chunk-11-1.png)<!-- -->

Note: the values at -12 from amino_acid_pauses.txt

``` r
data <- read.table("/Users/vip/Documents/CUG_reprocessed/Figure 1/aa_pauses_ypd_-18.txt", header = TRUE, row.names = 1)

data_avg <- data %>%
  mutate(ProArgDepletion_avg = rowMeans(select(., starts_with("proline")), na.rm = TRUE),
         Control_avg = rowMeans(select(., starts_with("control")), na.rm = TRUE))

ggplot(data_avg, aes(x = Control_avg, y = ProArgDepletion_avg, label = rownames(data_avg))) +
  geom_point(size = 3, color = ifelse(rownames(data_avg) %in% c("ARG", "PRO"), "red", "black")) +
  geom_text_repel(
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey',
    segment.size = 0.5,
    aes(label = rownames(data_avg)),
    color = ifelse(rownames(data_avg) %in% c("ARG", "PRO"), "red", "black")
  ) +
  labs(title = "Position -12",
       x = "Control (5Pseq Reads)",
       y = "Arg-Pro Depleted (5PSeq Reads)") +
  theme_minimal()
```

    ## Warning: ggrepel: 4 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](unnamed-chunk-12-1.png)<!-- -->
