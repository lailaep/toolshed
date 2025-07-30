---
description: Predicting Glycan Binding Affinities for Predicted Orbaceae Lectins
output-file: flycans_072825.html
title: Flycans 7/28/25

---


## Initialize glycowork based on original example:


```python
#| include: false
from nbdev.showdoc import *
from IPython.display import HTML
import pandas as pd
import copy
import warnings
import torch
warnings.filterwarnings("ignore")
from glycowork.glycan_data.loader import df_species, glycan_binding, df_glycan, glycomics_data_loader
from glycowork.motif.analysis import plot_embeddings, get_heatmap, characterize_monosaccharide, get_pvals_motifs, get_differential_expression, get_volcano
from glycowork.motif.processing import presence_to_matrix
from glycowork.motif.query import get_insight
from glycowork.motif.draw import GlycoDraw
from glycowork.ml.model_training import train_ml_model, analyze_ml_model, get_mismatch
from glycowork.ml.train_test_split import general_split
from glycowork.network.biosynthesis import construct_network, plot_network, evoprune_network, highlight_network, net_dic, network_alignment
from glycowork.network.evolution import distance_from_metric, dendrogram_from_distance
from glycowork.ml.models import prep_model
from glycowork.ml.inference import get_lectin_preds
#%load_ext autoreload
#%autoreload 2
```

## Setting up my data (predicted CU fimbrial lectins):


```python
#!pip install biopython
```


```python
from Bio import SeqIO

# Load protein sequences from FASTA file
sequences = []
names = []
path_to_fasta = "/Users/Snacks/software/glycowork/predicted_fim_adhesin_sequences_rn.fasta"
for record in SeqIO.parse(path_to_fasta, "fasta"):
    seq = str(record.seq)
    # Optionally truncate or pad sequences to 528 aa if needed
#    if len(seq) > 528:
#        seq = seq[:528]
#    elif len(seq) < 528:
#        seq = seq + 'A'*(528 - len(seq))  # Padding with alanines (or choose another padding strategy)

    sequences.append(seq)
    names.append(record.id)

```

Load the LectinOracle model:


```python
from glycowork.ml.models import prep_model
from glycowork.ml.inference import get_lectin_preds

# Load pretrained LectinOracle model
model = prep_model("LectinOracle_flex", 1, trained=True)
```

## Predict binding:

From OG example: "Fortunately, glycowork can be really helpful to answer such a question! The binding specificity of 1,392 lectins for a vast range of glycans is available in the glycan_binding dataset. Let's start by importing these data."


```python
#| echo: false
glycan_binding.head().style.set_properties(**{'font-size': '11pt', 'font-family': 'Helvetica','border-collapse': 'collapse','border': '1px solid black'})
```


```python
# Checking heatmap from example first
get_heatmap(glycan_binding.iloc[:,:-2], motifs = True, feature_set = ['exhaustive'], datatype = 'response', yticklabels = 1,
            xticklabels = False)
```


```python
# Check what's in glycowork.glycan_data.loader
import glycowork.glycan_data.loader as loader
print(dir(loader))

```


```python
from glycowork.glycan_data.loader import glycan_binding

# glycan_binding is a DataFrame loaded lazily from 'v12_glycan_binding.csv'
# Get glycan names (columns), excluding non-glycan metadata columns like 'Motif', 'Binding_type', etc.

glycan_columns = [col for col in glycan_binding.columns if col not in ['Motif', 'Binding_type', 'Class', 'Type', 'Family', 'Species']]
#print(glycan_columns[:50])  # see some glycan names
# Now use this as your glycan list for get_lectin_preds:
glycan_list = glycan_columns

```


```python
# Original in case:
#from glycowork.ml.inference import get_lectin_preds
#from glycowork.ml.models import prep_model
#leor = prep_model("LectinOracle_flex", 1, trained = True)
#get_lectin_preds("QWERTFVCF", ['Neu5Ac(a2-6)GalNAc', 'Gal(b1-4)Glc'], leor, flex = True)
from glycowork.ml.inference import get_lectin_preds
from glycowork.ml.models import prep_model
leor = prep_model("LectinOracle_flex", 1, trained = True)
```


```python
# Test dataset
get_lectin_preds("QWERTFVCF", ['Neu5Ac(a2-6)GalNAc', 'Gal(b1-4)Glc'], leor, flex = True)
```


```python
# Remove 'target' and 'protein' from list of glycan names

# Define known non-glycan columns you want to exclude
exclude_cols = ['target', 'protein']

# Filter glycans list
glycan_list = [col for col in glycan_columns if col not in exclude_cols]

print(glycan_list[:10])
```


```python
# Predict glycan binding for all sequences
predictions = []
for seq in sequences:
    pred = get_lectin_preds(glycans=glycan_list, 
                            prot=seq,
                            model=leor, # Could also try torch.nn.Module trained LectinOracle-type model
                            background_correction=False,
                            correction_df=None,
                            batch_size=128,
                            libr=None,
                            sort=True,
                            flex=True) # using LectinOracle-flex)
    predictions.append(pred)
```

## Combine predictions into a matrix


```python
# Convert list of prediction dataframes to a single matrix
import pandas as pd

# Assume all DataFrames have the same glycans in the same order
binding_matrix = pd.DataFrame()

#for i, pred in enumerate(predictions):
#    binding_matrix[f'seq_{i+1}'] = pred['pred'].values
for i, pred in enumerate(predictions):
    binding_matrix[names[i]] = pred['pred'].values  # use name instead of seq_i

# Set glycans as row labels
binding_matrix.index = predictions[0]['motif']
binding_matrix = binding_matrix.T  # Transpose to have glycans as columns
binding_matrix.head
```

## Plot data as heatmap (cluttered!)


```python
get_heatmap(binding_matrix,
            motifs=True,
            feature_set=['exhaustive'],
            datatype='response',
            yticklabels=1,
            xticklabels=True)
```

## Keep only glycans that have high binding in at least one lectin:


```python
import matplotlib.pyplot as plt
top_glycans = binding_matrix.columns[(binding_matrix > 1.2).any()]
reduced_matrix = binding_matrix[top_glycans]
```

## Plot data as heatmap (cleaned up?):


```python
# Replot (sequences)
plt.figure(figsize=(10, max(6, len(reduced_matrix) * 0.4)))
flycan_seq_heatmap = get_heatmap(reduced_matrix,
                motifs=False,
                feature_set=['exhaustive'],
                datatype='response',
                yticklabels=1,
                xticklabels=True)

flycan_motif_heatmap = get_heatmap(reduced_matrix,
                motifs=True,
                feature_set=['exhaustive'],
                datatype='response',
                yticklabels=1,
                xticklabels=True)
```


    <Figure size 1000x1440 with 0 Axes>



    
![png](output_25_1.png)
    



    
![png](output_25_2.png)
    

