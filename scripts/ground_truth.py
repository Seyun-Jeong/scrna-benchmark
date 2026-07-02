PBMC3K_MARKERS = {
    "Naive CD4+ T":     ["IL7R", "CCR7"],
    "CD14+ Monocyte":   ["CD14", "LYZ"],
    "B cell":           ["MS4A1"],
    "CD8+ T":           ["CD8A"],
    "NK":               ["GNLY", "NKG7"],
    "FCGR3A+ Monocyte": ["FCGR3A", "MS4A7"],
    "Dendritic":        ["FCER1A", "CST3"],
    "Platelet":         ["PPBP"],
}

def jaccard(set_a, set_b):
    a, b = set(set_a), set(set_b)
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)

def best_cell_type_match(top_genes, marker_dict=PBMC3K_MARKERS, top_n=10):
    top_genes = set(list(top_genes)[:top_n])
    scores = {ct: jaccard(top_genes, markers) for ct, markers in marker_dict.items()}
    best = max(scores, key=scores.get)
    return best, scores[best], scores
