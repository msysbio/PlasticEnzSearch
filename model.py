#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 15:44:59 2023

@author: u0145079
"""

import torch
import esm
from sklearn import svm
from Bio import SeqIO

# Load ESM-1b model
model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
batch_converter = alphabet.get_batch_converter()

# Read the FASTA files and parse the sequences
def read_fasta_file(file):
    records = list(SeqIO.parse(file, "fasta"))
    return [(record.id, str(record.seq)) for record in records]

# Prepare data (sequences)
known_petase_sequences = read_fasta_file("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/Plastic-tool/Fastas_without_duplicate_titles/Fasta_improved/PET.fa")

# Convert data to model input format
batch_labels, batch_strs, batch_tokens = batch_converter(known_petase_sequences)

# Extract per-residue representations (on CPU)
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=False)

# Get the representation at the [CLS] token
known_petase_representations = [results["representations"][33][i, 0].numpy() for i, label in enumerate(batch_labels)]

# Train a one-class SVM
clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.1)
clf.fit(known_petase_representations)

# Predict for putative PETase sequences
putative_sequences = read_fasta_file("putative_sequences.fasta")
batch_labels, batch_strs, batch_tokens = batch_converter(putative_sequences)

with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=False)

putative_representations = [results["representations"][33][i, 0].numpy() for i, label in enumerate(batch_labels)]

predictions = clf.predict(putative_representations)

print(predictions)
