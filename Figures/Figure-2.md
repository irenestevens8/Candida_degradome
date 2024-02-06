## Figure 2

Irene Stevens

23/05/2023

## This script will generate Figure 2

## Load libraries

``` r
library(ggplot2)
library(ggrepel)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ lubridate 1.9.2     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.2     ✔ tidyr     1.3.0
    ## ✔ readr     2.1.4

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ gridExtra::combine() masks dplyr::combine()
    ## ✖ dplyr::filter()      masks stats::filter()
    ## ✖ dplyr::lag()         masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

## Load Differentially Expressed Genes

See the DeSeq analysis section for how this input file was generated

``` r
ctx<-read.delim('/Users/vip/Documents/Candida/RNA-seq/DEG_sc5314_geneid.txt', sep="\t", header=T)
```

## Create volcano plot

Note: here I set log2FoldChange=0.5

``` r
ctx <- ctx %>% mutate(gene_type = case_when(log2FoldChange >= -0.5 & padj <= 0.05 ~ "up",log2FoldChange <= 0.5 & padj <= 0.05 ~ "down",TRUE ~ "non_significant"))
cols <- c("up" = "#870052", "down" = "#4DB5BC", "non_significant" = "#DDDEE0")
sizes <- c("up" = 3, "down" = 3, "non_significant" = 2)
alphas <- c("up" = 3, "down" = 3, "non_significant" = 2)

volcano<- ctx %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), fill = gene_type, size = gene_type, alpha = gene_type)) +
  geom_point(shape = 21, colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  scale_fill_manual(values = cols) +
  scale_size_manual(values = sizes) +
  scale_alpha_manual(values = alphas) +
  scale_x_continuous(breaks = c(seq(-3, 3, 2)), limits = c(-8, 8)) +
  geom_text_repel(
    data = ctx %>%
      filter(gene_id %in% c("ERG251", "ERG6", "ERG1", "ERG24", "MET13", "FRP1", "ERG11", "VPS23", "ERG2", "ERG5", "ERG4", "FMP45", "BTA1", "ERG7", "ALS7", "SOU2")),
    aes(label = gene_id),
    force = 2,
    nudge_y = 1
  ) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(breaks = c(seq(-2, 2, 1)), limits = c(-2, 2)) +
  theme_minimal() +  # Set the plot theme to minimal
  theme(plot.background = element_rect(fill = "white"))  # Change the background color to white
```

    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

``` r
volcano
```

    ## Warning: Removed 147 rows containing missing values (`geom_point()`).

![](unnamed-chunk-33-1.png)

## Plot Gene Ontology Terms downloaded from DAVID

Load DAVID enriched terms (see Suplementary Table S4)

``` r
data_up <- read.delim('/Users/vip/Documents/Candida/RNA-seq/GO_sc5314_up_genes.txt', header = TRUE)
data_up <- data_up[order(-data_up$Fold.Enrichment), ]

data_down <- read.delim('/Users/vip/Documents/Candida/RNA-seq/GO_sc5314_down_genes.txt', header = TRUE)
data_down <- data_down[order(-data_down$Fold.Enrichment), ]
```

Barplot of Enriched GO Terms in transcriptionally upregulated genes

``` r
plot_up <- ggplot(data_up, aes(x = reorder(Term, Fold.Enrichment), y = Fold.Enrichment, fill=-log10(PValue))) + 
  coord_flip() +   
  geom_bar(stat = "identity", color = "black") +  
  scale_fill_gradient(low = "white", high = "#870052") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  labs(x = "GO Terms", y = "Fold Enrichment", title = "Transcriptionally upregulated")

plot_up
```

![](unnamed-chunk-55-1.png)

Barplot of Enriched GO Terms in transcriptionally downregulated genes

``` r
plot_down <- ggplot(data_down, aes(x = reorder(Term, Fold.Enrichment), y = Fold.Enrichment, fill=-log10(PValue))) + 
  coord_flip() +   
  geom_bar(stat = "identity", color = "black") +  
  scale_fill_gradient(low = "white", high = "#4DB5BC") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +  
  labs(x = "GO Terms", y = "Fold Enrichment", title = "Transcriptionally downregulated")

plot_down
```

![](unnamed-chunk-66-1.png)

## Violin plot of amino acid frame protetion index

Load frame protection index per amino acid file

``` r
data<-read.delim('/Users/vip/Documents/CUG_reprocessed/Figure 3/sc5314_amino_acid_FPI.txt', header=T, row.names=1)
```

Calculate average of 3 replicates for fluconazole and control

``` r
data_avg <- data %>%
  mutate(Fluconazole_avg = rowMeans(select(., starts_with("fluc_rep")), na.rm = TRUE),
         Control_avg = rowMeans(select(., starts_with("control_rep")), na.rm = TRUE))
```

Reshape data and reorder levels

``` r
data_long <- data_avg %>%
  gather(key = "Condition", value = "Value", Fluconazole_avg, Control_avg)
data_long$Condition <- factor(data_long$Condition, levels = c("Fluconazole_avg", "Control_avg"))
```

Plot

``` r
ggplot(data_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "white", color = "black") +
  scale_fill_manual(values = c("#870052", "#4DB5BC")) +  
  labs(title = "Amino Acid Frame Protection Index",
       x = "Condition",
       y = "Frame Protection Index")
```

![](unnamed-chunk-100-1.png)

Perform T-test

``` r
t_test_result <- t.test(Value ~ Condition, data = data_long)

t_test_result
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  Value by Condition
    ## t = -3.8403, df = 38.367, p-value = 0.0004476
    ## alternative hypothesis: true difference in means between group Fluconazole_avg and group Control_avg is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.03457130 -0.01070931
    ## sample estimates:
    ## mean in group Fluconazole_avg     mean in group Control_avg 
    ##                    0.06968436                    0.09232467
