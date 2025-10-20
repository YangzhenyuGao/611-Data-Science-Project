# Makefile for the BIOS6611 HW5 Clustering Project
#
# Author: Yangzhenyu Gao
# Date: October 20, 2025
#
# Description:
# This Makefile automates the execution of the R analysis script.
# Its purpose is to run the script from start to finish to verify
# that the analysis completes without any errors.

# --- Variables ---
# The primary R script for the analysis.
R_SCRIPT = BIOS661HW5.R

# --- Rules ---

# Define the default goal. Typing 'make' in the terminal will execute this 'run' rule.
.DEFAULT_GOAL := run

# The primary rule to execute the analysis.
# This is a .PHONY target because 'run' is an action, not a file to be created.
.PHONY: run
run:
	@echo "--- Executing the analysis script: $(R_SCRIPT) ---"
	# The command below runs the R script from the command line.
	# The script will print its messages and generate plots in the RStudio pane if run there,
	# but no files will be saved.
	@Rscript $(R_SCRIPT)
	@echo "--- Script finished successfully. ---"