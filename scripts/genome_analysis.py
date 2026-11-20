"""
Calculate chromosome-level statistics for the
Zymoseptoria tritici IPO323 genome assembly.

Statistics include chromosome length, GC content,
and genome assembly summary metrics.
"""

from pathlib import Path
from Bio import SeqIO
import pandas as pd

# 1. Find the genome FASTA file
data_folder = Path("../data/raw/IPO323/")
fasta_files = list(data_folder.rglob("*.fna"))

if len(fasta_files) != 1:
    raise FileNotFoundError(
        f"Expected one .fna file, but found {len(fasta_files)}. "
        "Check your downloaded genome files."
    )

fasta_file = fasta_files[0]
print(f"Reading genome file: {fasta_file}")


# 2. Read the genome sequences
records = list(SeqIO.parse(fasta_file, "fasta"))

if not records:
    raise ValueError("No sequences were found in the FASTA file.")


# 3. Calculate chromosome lengths
chromosome_lengths = [len(record.seq) for record in records]

total_length = sum(chromosome_lengths)
number_of_chromosomes = len(chromosome_lengths)

average_length = sum(chromosome_lengths) / number_of_chromosomes
sorted_lengths = sorted(chromosome_lengths)

median_length = sorted_lengths[number_of_chromosomes // 2]

if number_of_chromosomes % 2 == 0:
    median_length = (
        sorted_lengths[number_of_chromosomes // 2 - 1]
        + sorted_lengths[number_of_chromosomes // 2]
    ) / 2

shortest_chromosome = min(chromosome_lengths)
longest_chromosome = max(chromosome_lengths)


# 4. Calculate GC content
gc_count = 0
at_count = 0

for record in records:
    sequence = str(record.seq).upper()
    gc_count += sequence.count("G") + sequence.count("C")
    at_count += sequence.count("A") + sequence.count("T")

gc_content = gc_count / (gc_count + at_count) * 100


# 5. Calculate N50
def calculate_n50(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    half_total = sum(sorted_lengths) / 2
    cumulative_length = 0

    for length in sorted_lengths:
        cumulative_length += length

        if cumulative_length >= half_total:
            return length


n50 = calculate_n50(chromosome_lengths)


# 6. Print the results
print("\nGenome assembly statistics")
print("--------------------------")
print(f"Total genome size: {total_length:,} bp")
print(f"Number of chromosomes: {number_of_chromosomes}")
print(f"Average chromosome length: {average_length:,.2f} bp")
print(f"Median chromosome length: {median_length:,.2f} bp")
print(f"Shortest chromosome: {shortest_chromosome:,} bp")
print(f"Longest chromosome: {longest_chromosome:,} bp")
print(f"GC content: {gc_content:.2f}%")
print(f"N50: {n50:,} bp")


# 7. Save individual contig statistics
chromosome_data = []

for record in records:
    sequence = str(record.seq).upper()
    sequence_gc = (
        (sequence.count("G") + sequence.count("C"))
        / (sequence.count("A") + sequence.count("T")
           + sequence.count("G") + sequence.count("C"))
        * 100
    )

    chromosome_data.append({
        "chromosome_id": record.id,
        "length": len(record.seq),
        "gc_content": sequence_gc
    })

chromosome_df = pd.DataFrame(chromosome_data)

output_folder = Path("../data/processed")
output_folder.mkdir(parents=True, exist_ok=True)

chromosome_df.to_csv(
    output_folder / "chromosome_statistics.csv",
    index=False
)

print("\nChromosome statistics saved to data/processed/")

