import sys
import csv
import re
import pandas as pd

#csv.field_size_limit(sys.maxsize)

INPUT_CSV = "full_forest_export.csv"
OUTPUT_CSV = "node_mutation.csv"

# Regex: right-most nXX before stop_reason=5
pattern = re.compile(
    r'(n\d+)(?=(?:(?!\bn\d+\b).)*stop_reason=5)'
)

with open(INPUT_CSV, newline="", encoding="utf-8") as fin, \
     open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as fout:

    reader = csv.DictReader(fin)
    writer = csv.writer(fout)

    # Write header
    writer.writerow(["index", "extracted_n"])

    for row in reader:
        tree = row["tree"]
        match = pattern.search(tree)

        writer.writerow([
            row["index"],
            match.group(1) if match else ""
        ])



df = pd.read_csv("node_mutation.csv")

print(df["extracted_n"])