############################################################
## GEN242 AI Tutor for Coding Exercise Template
##
## Purpose:
##   This script provides two small exercises for demonstrating
##   coding-based AI tools such as GitHub Copilot and Claude Code.
##
## Exercises:
##   1. Simple ggplot2 example using the built-in diamonds dataset
##   2. Slightly more complex volcano plot example using a made-up
##      differential gene expression (DEG) result table
##
## Suggested use:
##   - Copy individual sections into RStudio.
##   - Use Copilot to request suggestions.
##   - Use Claude Code to explain, revise, and debug the code.
##   - Inspect and run all generated code before accepting it.
##
## RStudio Copilot shortcut:
##   Ctrl + \   request suggestion
##   Tab        accept suggestion
##   Esc        dismiss suggestion
############################################################


############################################################
## Setup
############################################################

## Load required package
## If ggplot2 is not installed, run:
## install.packages("ggplot2")

library(ggplot2)


############################################################
## Exercise 1: Simple ggplot2 example with diamonds
############################################################

## Goal:
##   Use Copilot to create a simple scatter plot showing the
##   relationship between diamond carat, price, and cut.

## Instructions:
##   1. Put the cursor after the comment below.
##   2. In RStudio, press Ctrl + \ to request a Copilot suggestion.
##   3. Press Tab to accept the suggestion.
##   4. Run the code.
##   5. Ask Copilot or Claude to explain what the code does.
##   6. Ask for one improvement to the plot.

# Create a scatter plot of carat versus price colored by cut


## Optional seed if Copilot does not suggest code automatically:
## Uncomment the next line and request a suggestion after the opening aes(.
# ggplot(diamonds, aes(


## Possible follow-up prompts for Copilot/Claude:
##
## Explain this ggplot2 code line by line.
## Improve this plot by adding transparency, labels, and a minimal theme.
## Why is alpha transparency useful in this plot?
## What limitation does this plot have when many points overlap?


############################################################
## Exercise 1B: Intentional debugging example
############################################################

## Goal:
##   Introduce a simple error and ask Copilot or Claude to fix it.

## Instructions:
##   1. Run the code below.
##   2. It contains an intentional error: the column name "carats"
##      does not exist in the diamonds dataset.
##   3. Select the code and the error message.
##   4. Ask CopilotChat or Claude:
##
##      This ggplot2 code gives an error. Explain why the error
##      happens and fix the code.

## Intentional error: "carats" should be "carat"
p_diamonds_error <- ggplot(diamonds, aes(x = carats, y = price, color = cut)) +
  geom_point(alpha = 0.3, size = 0.8) +
  labs(
    title = "Diamond price by carat and cut",
    x = "Carat",
    y = "Price",
    color = "Cut"
  ) +
  theme_minimal()

p_diamonds_error


############################################################
## Exercise 2: Made-up DEG result table
############################################################

## Goal:
##   Use a small made-up differential gene expression result table
##   to create a volcano plot function.

## Notes:
##   - This is not a real DEG analysis.
##   - The table is only for practicing plotting, debugging,
##     and AI-assisted code explanation.
##   - A real DEG table would usually come from tools such as
##     DESeq2, edgeR, or limma-voom.


## Create a small made-up DEG result table
deg_results <- data.frame(
  gene = paste0("gene", 1:30),
  log2FoldChange = c(
    3.2, -2.8, 1.5, -1.2, 0.4,
    2.4, -3.1, 0.2, 1.1, -0.7,
    4.0, -2.2, 0.8, -1.8, 2.9,
    -0.3, 1.9, -2.6, 0.5, -1.4,
    0.9, -0.2, 2.1, -3.4, 1.3,
    -1.1, 0.1, 3.6, -2.0, 0.6
  ),
  padj = c(
    0.0005, 0.002, 0.04, 0.08, 0.7,
    0.01, 0.0008, 0.9, 0.15, 0.4,
    0.0001, 0.03, 0.2, 0.04, 0.006,
    0.8, 0.049, 0.01, 0.6, 0.07,
    0.12, 0.95, 0.02, 0.0003, 0.05,
    0.09, 0.99, 0.0002, 0.025, 0.5
  )
)

## Inspect the table
head(deg_results)
str(deg_results)
summary(deg_results)


############################################################
## Exercise 2A: Ask Copilot to write a volcano plot function
############################################################

## Instructions:
##   1. Put the cursor after the comment block below.
##   2. In RStudio, press Ctrl + \ to request a Copilot suggestion.
##   3. Press Tab to accept the suggestion.
##   4. Run the function on deg_results.
##   5. Ask Copilot or Claude to explain the function line by line.

# Write a function called plot_volcano that takes a DEG result table
# with columns gene, log2FoldChange, and padj.
# The function should create a volcano plot using ggplot2.
# Plot log2FoldChange on the x-axis and -log10(padj) on the y-axis.
# Color genes as significant if padj < 0.05 and abs(log2FoldChange) > 1.
# Add vertical cutoff lines at -1 and 1 and a horizontal cutoff line at padj = 0.05.
# Return the ggplot object.


## Optional seed if Copilot does not suggest code automatically:
# plot_volcano <- function(deg_table, lfc_cutoff = 1, padj_cutoff = 0.05) {


## One possible reference solution is provided below.
## Keep it commented at first if you want students to generate their own solution.

plot_volcano_reference <- function(deg_table,
                                   lfc_cutoff = 1,
                                   padj_cutoff = 0.05) {
  ## Add a logical column indicating whether each gene passes
  ## both the adjusted p-value and fold-change cutoffs.
  deg_table$significant <- with(
    deg_table,
    padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff
  )

  ## Create the volcano plot.
  ggplot(deg_table, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant), size = 2, alpha = 0.8) +
    geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed"
    ) +
    labs(
      title = "Volcano plot of made-up DEG results",
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = "Significant"
    ) +
    theme_minimal()
}

## Run the reference function
plot_volcano_reference(deg_results)


############################################################
## Exercise 2B: Intentional debugging example
############################################################

## Goal:
##   Introduce a realistic column-name error and ask Copilot or
##   Claude to explain and fix it.

## Instructions:
##   1. Run the function below.
##   2. It contains an intentional error: it refers to "logFC",
##      but the table contains "log2FoldChange".
##   3. Select the function and the error message.
##   4. Ask CopilotChat or Claude:
##
##      This volcano plot function gives an error that object
##      'logFC' is not found. Explain why the error happens and
##      fix the function.

plot_volcano_error <- function(deg_table,
                               lfc_cutoff = 1,
                               padj_cutoff = 0.05) {
  ## Intentional error:
  ## "logFC" is not a column in deg_table.
  ## The correct column name is "log2FoldChange".
  deg_table$significant <- with(
    deg_table,
    padj < padj_cutoff & abs(logFC) > lfc_cutoff
  )

  ggplot(deg_table, aes(x = logFC, y = -log10(padj))) +
    geom_point(aes(color = significant), size = 2, alpha = 0.8) +
    geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed"
    ) +
    labs(
      title = "Volcano plot of made-up DEG results",
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = "Significant"
    ) +
    theme_minimal()
}

## Run this to produce the intentional error
plot_volcano_error(deg_results)


############################################################
## Exercise 2C: Optional more subtle data-analysis issue
############################################################

## Goal:
##   Show that code can be syntactically correct but still needs
##   careful handling of data edge cases.

## Issue:
##   In volcano plots, the y-axis often uses -log10(padj).
##   If padj is exactly zero, then -log10(0) is Inf.
##   This can create plotting problems or misleading axes.

## Create a copy of the DEG table with one adjusted p-value set to zero
deg_results_zero <- deg_results
deg_results_zero$padj[1] <- 0

## Run the reference function on the modified table
## Depending on the plotting behavior, this may generate an infinite value.
plot_volcano_reference(deg_results_zero)

## Ask Copilot or Claude:
##
## One adjusted p-value is zero, which causes an infinite -log10 value.
## How should this be handled in a volcano plot function?
## Explain the issue and suggest a safe fix.


## One possible safer version:
plot_volcano_safe <- function(deg_table,
                              lfc_cutoff = 1,
                              padj_cutoff = 0.05) {
  ## Copy the input table so the original object is not modified.
  deg_table <- as.data.frame(deg_table)

  ## Replace zero adjusted p-values only for plotting.
  ## The original padj values are kept for significance testing.
  deg_table$padj_plot <- pmax(deg_table$padj, .Machine$double.xmin)

  deg_table$significant <- with(
    deg_table,
    padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff
  )

  ggplot(deg_table, aes(x = log2FoldChange, y = -log10(padj_plot))) +
    geom_point(aes(color = significant), size = 2, alpha = 0.8) +
    geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed"
    ) +
    labs(
      title = "Volcano plot with safe handling of zero adjusted p-values",
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = "Significant"
    ) +
    theme_minimal()
}

plot_volcano_safe(deg_results_zero)


############################################################
## Suggested questions to ask for understanding
############################################################

## 1. Which parts of the generated code did you understand immediately?
## 2. Which parts required explanation?
## 3. Did the AI-generated code run without modification?
## 4. Did the AI explain the error correctly?
## 5. What assumptions does the volcano plot function make about the input table?
## 6. How would you make the function more robust for real DEG results?
## 7. Could you explain the final code without using AI?
