import tensorflow_hub as hub
import tensorflow as tf
import numpy as np
import os
import math
from pyfasta import Fasta
import argparse

# Load the Enformer model
enformer_model = hub.load("https://tfhub.dev/deepmind/enformer/1").model

def seq_to_mat(seq):
    """Convert DNA sequence to a one-hot encoded matrix."""
    d = {'a': 0, 'A': 0, 'c': 1, 'C': 1, 'g': 2, 'G': 2, 't': 3, 'T': 3, 'N': 4, 'n': 4}
    mat = np.zeros((5, len(seq)))
    for i in range(len(seq)):
        mat[d[seq[i]], i] = 1
    mat = mat[:4, :]
    return mat

def main(args):
    gene2loc = {item.split('\t')[4]: (item.split('\t')[0], int(item.split('\t')[1])) for item in open(args.refGene_path).readlines()}

    assert args.gene_name in gene2loc, "Gene not found in refGene database"
    chr_id = gene2loc[args.gene_name][0]
    center = gene2loc[args.gene_name][1]
    start = center - 100000
    end = center + 100000

    SEQ_LENGTH = 393216
    interval = 896 * 128
    nb_regions = math.ceil((end - start - interval) / (2 * interval))

    os.makedirs(args.output_path, exist_ok=True)
    output_file = os.path.join(args.output_path, f"{args.fasta_path.split('/')[-1].split('.')[0]}.npy")

    key = "maternal" if "maternal" in args.fasta_path else "paternal"
    genome = Fasta(args.fasta_path)

    enformer_feats = []
    for coor in range(center - interval * nb_regions, center + interval * (nb_regions + 1), interval):
        seq = genome[f'{chr_id}_{key}'][(coor - SEQ_LENGTH // 2):(coor + SEQ_LENGTH // 2)]
        onehot_mat = seq_to_mat(seq).T
        onehot_mat = np.expand_dims(onehot_mat, 0)
        enformer_feats.append(enformer_model.predict_on_batch(onehot_mat)['human'])

    enformer_feats = np.squeeze(np.stack(enformer_feats))
    np.save(output_file, enformer_feats)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Enformer features for a given gene.")
    parser.add_argument('--gene_name', type=str, required=True, default='APOE', help='Name of the gene.')
    parser.add_argument('--fasta_path', type=str, required=True, default='chr19_003_S_1057_maternal.fa', help='Path to the fasta file.')
    parser.add_argument('--refGene_path', type=str, required=True, help='Path to the refGene hg19 TSS bed file.')
    parser.add_argument('--output_path', type=str, required=True, help='Directory to save the output files.')

    args = parser.parse_args()
    main(args)
