#!/bin/bash


# Create directories
mkdir -p config/docker \
		 data/raw data/processed data/metadata \
		 notebooks/exploratory notebooks/final notebooks/archive \
		 scripts/R scripts/python scripts/bash \
		 workflows/modules \
		 results/tables results/figures results/logs \
		 docs \
		 downloads \
		 tests/data_tests tests/workflow_tests \
		 .github/workflows \
		 logs

# Notify completion
echo "Folder structure created successfully."
