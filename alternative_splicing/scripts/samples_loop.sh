#!/bin/bash

# List of cell lines
cell_lines=("LINENAME_1" "LINENAME_2")

# List of condition pairs
condition_pairs=(
    "State1A State1U"
    "State1B State1U"
    "State2A State2U"
    "State2B State2U"
)

# Iterate over cell lines
echo "Analysis started"
echo ""
for cell_line in "${cell_lines[@]}"; do
    # Iterate over condition pairs
    for pair in "${condition_pairs[@]}"; do
        # Extract conditions from the pair
        cond1=$(echo "$pair" | awk '{print $1}')
        cond2=$(echo "$pair" | awk '{print $2}')

        # Run the Python command
        echo "cell line: $cell_line, cond1: $cond1, cond2: $cond2"
        python AS_analysis.py --line_name "$cell_line" --cond1 "$cond1" --cond2 "$cond2"
    done
done
