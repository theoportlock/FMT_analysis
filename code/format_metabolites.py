#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: Theo Portlock
For diversity analysis
'''
import pandas as pd

df = pd.read_csv(f"data/inflammeta.csv")
df = df.dropna(subset = 'Sample ID')
df = df.set_index('Sample ID')

# Filter for the metabolites
df = df[[
 'Acetic_acid',
 'Alanine',
 'Butyric_acid',
 'Ethanol',
 'Formic_acid',
 'Fumaric_acid',
 'Glutamic_acid',
 'Glycine',
 'Isoleucine',
 'Lactic_acid',
 'Leucine',
 'Methanol',
 'N_trimethylamine',
 'N_methylamine',
 'Phenylacetic_acid',
 'Propionic_acid',
 'Succinic_acid',
 'Tyrosine',
 'Uracil',
 'Valeric_acid',
 'Valine',
 'IL-17A',
 'IL-17E',
 'IL-17F',
 'IL-21',
 'IL-22',
 'IFNgamma',
 'IL-10',
 'IL-1beta',
 'IL-6',
 'TNF-alpha',
 'IL-12',
 'IL-23',
 'IL-8',
 'FABP2',
 'DLactate_faecal',
 'DLactate_plasma',
 'faecal_ammonia',
 'endotoxin',
 'INR',
 'CRP',
 'Sodium',
 'Potassium',
 'Creatinine',
 'Urea',
 'Haemoglobin',
 'MCV',
 'Platelets',
 'Neuts',
 'Lymphs',
 'neut/lymph',
 'WCC',
 'Bilirubin',
 'ALT',
 'AST',
 'ALP',
 'GGT',
 'Albumin',
 'Ammonia',
 'Lactate',
 'Ascites',
 'HE grade',
 'Diabetes',
 'TIPS',
 'Beta blocker',
 'Lactulose',
 'Rifaximin',
 'PPI',
 'Other meds',
 'Calprotectin/(µg/g)',
 'N/L ratio',
 'Calprotectin (µg/g)',
 ]]

df.index.name = 'sampleID'

df.to_csv(f'results/metabolites.tsv', sep='\t')
