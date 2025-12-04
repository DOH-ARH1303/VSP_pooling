
# To add: Conversion of viral conc to Ct value, final grouping of matrixes >3 (could maybe just do for all matrix groups), csv/excel file conversion (probably pandas)
# How does it handle s1=Ct 20, s2=Ct 23, s3=Ct 17?
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
        print(sid, cluster)

# order samples within each cluster by ct value and then split into pools of 3 starting with the lowest ct

def main():
    dataset = [
        {
            "id": "sample1",
            "ct": 28
        },
                {
            "id": "sample2",
            "ct": 16
        },
                {
            "id": "sample3",
            "ct": 18
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

    create_clusters(dataset, 3, 1)



if __name__ == "__main__":
    main()
