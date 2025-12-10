'''
# To add: Conversion of viral conc to Ct value, final grouping of matrixes >3 (could maybe just do for all matrix groups), csv/excel file conversion (probably pandas)
# How does it handle s1=Ct 20, s2=Ct 23, s3=Ct 17?
    # Ct 17 and Ct 23 end up in the same group
import numpy as np
from sklearn.cluster import DBSCAN

def compute_matrix(data):
    """
    Compute pairwise distance matrix for sequences in a window using cached distances.
    """
    ids = [rec["id"] for rec in data]
    n = len(ids)
    mat = np.zeros((n, n), dtype=float)
    for i in range(n):
        ai = ids[i]
        for j in range(i + 1, n):
            aj = ids[j]
            d = abs(data[i]['ct'] - data[j]['ct'])
            mat[i, j] = mat[j, i] = d
    np.fill_diagonal(mat, 0.0)
    return mat, ids

def create_clusters(data, threshold, min_samples):
    mat, ids = compute_matrix(data)
    n = len(ids)

    db = DBSCAN(eps=threshold, min_samples=min_samples, metric='precomputed').fit(mat)
    labels = db.labels_  # -1 for noise, >=0 for clusters

    for (sid, cluster) in zip(ids, labels):
        # Divide into df based on cluster #
        # Combine with initial df
        print(sid, cluster)

# order samples within each cluster by ct value and then split into pools of 3 starting with the lowest ct

def main():
    dataset = [
        {
            "id": "sample1",
            "ct": 21
        },
                {
            "id": "sample2",
            "ct": 16
        },
                {
            "id": "sample3",
            "ct": 19
        },        
        {
            "id": "sample4",
            "ct": 32
        },    
        {
            "id": "sample5",
            "ct": 15
        },    
        {
            "id": "sample6",
            "ct": 17
        },      
        {
            "id": "sample7",
            "ct": 18
        }
    ]

    create_clusters(dataset, 1, 1)



if __name__ == "__main__":
    main()
'''

import pandas as pd

# Example dataframe
sample_column = 'WA#'
ct_coumn = 'Ct'

data = {
    f'{sample_column}': ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'],
    'Ct': [20.0, 21.5, 22.0, 24.0, 25.5, 27.0]
}
df = pd.DataFrame(data)

# Sort by Ct values
df = df.sort_values(by='Ct').reset_index(drop=True)

clusters = []
used = set()

for i in range(len(df)):
    if i in used:
        continue
    
    # Start a new cluster with the current sample
    cluster = [df.loc[i, f'{sample_column}']]
    used.add(i)
    
    # Try to add up to 2 more samples within Ct difference <= 2
    for j in range(i+1, len(df)):
        if j in used:
            continue
        if abs(df.loc[j, 'Ct'] - df.loc[i, 'Ct']) <= 2:
            cluster.append(df.loc[j, f'{sample_column}'])
            used.add(j)
        if len(cluster) == 3:  # max cluster size
            break
    
    clusters.append(cluster)

# Have each wa# listed w/ cluster number in dataframe
# Add negative control into cluster with <3 samples (lowest or highest Cts?)
# Print clusters
for idx, c in enumerate(clusters, 1):
    print(f"{idx}: {c}") 
