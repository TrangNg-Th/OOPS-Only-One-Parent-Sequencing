import pandas as pd
import os
import sys

# Set working directory to the script's directory
print("Current working directory ",
      os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Parse command line arguments
if len(sys.argv) != 3:
    print("Usage: python fix_PhaseSet.py <input_file> <output_file>")
    sys.exit(1)


# To run after already merged mom unphased and child phased VCFs into a single vcf
in_file = sys.argv[1]
out_file = sys.argv[2]

# Load the data
df = pd.read_csv(in_file, dtype=str)

# Add headers
df.columns = ["chrom", "pos", "PS_mom", "PS_child",
              "GT_mom", "GT_child", "DP_mom", "DP_child",
              "GQ_mom", "GQ_child",
              "ADref_mom", "ADalt_mom", "ADref_child", "ADalt_child",
              "None"]

# Drop the last empty column
df = df.drop(columns=["None"])

# Write file with tab delimiter
df.to_csv(out_file, sep="\t", index=False)