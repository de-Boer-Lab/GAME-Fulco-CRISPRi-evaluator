import pandas as pd

# import os
# os.environ["MPLCONFIGDIR"] = "/scratch/st-cdeboer-1/iluthra/mplconfig"
import matplotlib.pyplot as plt
import seaborn as sns
import pysam


fasta_path = "/arc/project/st-cdeboer-1/Genomes/human/GRChg37/hg19.fa"
genome = pysam.FastaFile(fasta_path)

#read in the enhancer + gene pairs
eg_pairs = pd.read_csv("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Fulco_CRISPRi/evaluator_data/Supplementary Table 3a. Data for E-G pairs.csv", sep = ";", engine = "python")
print(eg_pairs)
print(eg_pairs.Gene.unique())
genes_tested = eg_pairs.Gene.unique()
print("Gene's tested")
print(genes_tested)
eg_pairs = eg_pairs.iloc[:, :6]
print(eg_pairs)

# #plot scatter plot of fraction change distribution
# x_col = "Element name"
# y_col = "Fraction change in gene expr"
# gene_col = "Gene"

# plt.figure(figsize=(10,6))
# # Use Gene to color points
# scatter = plt.scatter(
#     eg_pairs[x_col],
#     eg_pairs[y_col],
#     c=eg_pairs[gene_col].astype('category').cat.codes,   # numeric codes for colouring
#     cmap='tab20'
# )
# # Build legend manually
# handles, labels = scatter.legend_elements(prop="colors", alpha=0.9)
# plt.legend(handles, eg_pairs[gene_col].unique(), title="Gene", bbox_to_anchor=(1.05, 1), loc='upper left')

# plt.xlabel(x_col)
# plt.ylabel("Fractional change in expression")
# plt.title("Element vs Fractional Expression Change (Colored by Gene)")
# plt.tight_layout()
# plt.savefig("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/Fulco_CRISPRi/evaluator_data/eg_scatter.png", dpi=300, bbox_inches='tight')



#Read in consensus gene coordinates
gene_coordinates = pd.read_csv("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Fulco_CRISPRi/evaluator_data/Supplementary Table 5b. Human gene annotations from RefSeq collapsed to a single consensus annotation per gene.csv", sep = ";", engine = "python")
#print(gene_coordinates)

#Extract coordinates for genes that were tested
gene_coordinates_fulco = gene_coordinates[gene_coordinates["Gene"].isin(genes_tested)]
gene_coordinates.columns = ['gene_chr', 'gene_start', 'gene_end', 'Gene', 'strand']
#print(gene_coordinates_fulco)
print("30 gene info")
print(gene_coordinates_fulco)

#Merge gene locations with e-g pair data
merged_eg_pairs = pd.merge(eg_pairs, gene_coordinates, how='left', on = 'Gene')
print(merged_eg_pairs.columns)

#Calculate distance from start of "enhancer" to start of gene -> can be up or downstream
merged_eg_pairs['distance_to_gene'] = abs(merged_eg_pairs["start"] - merged_eg_pairs["gene_start"])

#Pull sequence 2MB up and downstream
merged_eg_pairs["sequence_start"] = merged_eg_pairs['gene_start'] - 2000000
merged_eg_pairs["sequence_end"] = merged_eg_pairs['gene_start'] + 2000000
print(merged_eg_pairs)

#check that enhancers are in this sequence we are pulling
mask = (merged_eg_pairs["start"] >= merged_eg_pairs["sequence_start"]) & (merged_eg_pairs["end"] <= merged_eg_pairs["sequence_end"])
merged_eg_pairs_filtered = merged_eg_pairs[mask]
#We loose 32 e-g pairs when we set sequence window to 4MB
print(merged_eg_pairs_filtered)
# print(df_invalid)

#merged_eg_pairs_filtered = merged_eg_pairs_filtered.head(5)
sequences = []
for i in range(0,len(merged_eg_pairs_filtered)):
    chrom_current = merged_eg_pairs_filtered['gene_chr'].iloc[i]
    start_current = merged_eg_pairs_filtered['sequence_start'].iloc[i]
    end_current = merged_eg_pairs_filtered['sequence_end'].iloc[i]
    seq = genome.fetch(chrom_current, start_current,  end_current)
    seq = seq.upper()
    sequences.append(seq)

merged_eg_pairs_filtered['sequence'] = sequences

print("BOB")
print(merged_eg_pairs_filtered)


gene_and_sequence = merged_eg_pairs_filtered[["Gene", "sequence", "strand"]].copy()
print(gene_and_sequence)
gene_and_sequence_unique = gene_and_sequence.drop_duplicates()
print(gene_and_sequence_unique)
print(gene_and_sequence_unique.shape)

gene_and_sequence_unique.to_parquet("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Fulco_CRISPRi/evaluator_data/fulco_reference_gene_sequence.parquet", engine='pyarrow', compression='snappy')
print(merged_eg_pairs_filtered)
merged_eg_pairs_filtered = merged_eg_pairs_filtered.drop(columns=["sequence"])
print(merged_eg_pairs_filtered)

merged_eg_pairs_filtered.to_parquet("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/new_game_dev/Evaluators/Fulco_CRISPRi/evaluator_data/fulco_evaluator_coordinates.parquet", engine='pyarrow', compression='snappy')
