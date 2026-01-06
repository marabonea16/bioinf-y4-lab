"""
Exercise 9.1 — Drug–Gene Bipartite Network & Drug Similarity Network

Scop:
- să construiți o rețea bipartită drug–gene plecând de la un CSV
- să proiectați layer-ul de medicamente folosind similaritatea dintre seturile de gene
- să exportați un fișier cu muchiile de similaritate între medicamente

TODO:
- încărcați datele drug–gene
- construiți dict-ul drug -> set de gene țintă
- construiți graful bipartit drug–gene (NetworkX)
- calculați similaritatea dintre medicamente (ex. Jaccard)
- construiți graful drug similarity
- exportați tabelul cu muchii: drug1, drug2, weight
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, Tuple, List

import itertools

import networkx as nx
import pandas as pd

# --------------------------
# Config — adaptați pentru handle-ul vostru
# --------------------------
HANDLE = "marabonea16"

# Input: fișier cu coloane cel puțin: drug, gene
DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")

# Output directory & files
OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_DRUG_SIMILARITY = OUT_DIR / f"drug_similarity_{HANDLE}.csv"
OUT_GRAPH_DRUG_GENE = OUT_DIR / f"network_drug_gene_{HANDLE}.gpickle"


def ensure_exists(path: Path) -> None:
    """
    TODO:
    - verificați că fișierul există
    - dacă nu, ridicați FileNotFoundError cu un mesaj clar
    """
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")


def load_drug_gene_table(path: Path) -> pd.DataFrame:
    """
    TODO:
    - citiți CSV-ul cu pandas
    - validați că există cel puțin coloanele: 'drug', 'gene'
    - returnați DataFrame-ul
    """
    df = pd.read_csv(path)
    required_cols = {'drug', 'gene'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")
    return df


def build_drug2genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    TODO:
    - construiți un dict: drug -> set de gene țintă
    - sugestie: folosiți groupby("drug") și aplicați set() pe coloana gene
    """
    drug2genes = df.groupby("drug")["gene"].apply(set).to_dict()
    return drug2genes


def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    """
    TODO:
    - construiți graful bipartit:
      - nodurile 'drug' cu atribut bipartite="drug"
      - nodurile 'gene' cu atribut bipartite="gene"
      - muchii drug-gene
    """
    G = nx.Graph()
    
    # Add drug nodes
    G.add_nodes_from(drug2genes.keys(), bipartite="drug")
    
    # Collect all genes and add gene nodes
    all_genes = set()
    for genes in drug2genes.values():
        all_genes.update(genes)
    G.add_nodes_from(all_genes, bipartite="gene")
    
    # Add edges
    edges = []
    for drug, genes in drug2genes.items():
        for gene in genes:
            edges.append((drug, gene))
    G.add_edges_from(edges)
    
    return G


def summarize_drugs(drug2genes: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    TODO:
    - construiți un DataFrame cu:
        drug, num_targets (numărul de gene țintă)
    - returnați DataFrame-ul
    """
    summary = []
    for drug, genes in drug2genes.items():
        summary.append({"drug": drug, "num_targets": len(genes)})
    return pd.DataFrame(summary)


def jaccard_similarity(s1: Set[str], s2: Set[str]) -> float:
    """
    Calculați similaritatea Jaccard între două seturi de gene:
    J(A, B) = |A ∩ B| / |A ∪ B|
    """
    if not s1 and not s2:
        return 0.0
    inter = len(s1 & s2)
    union = len(s1 | s2)
    return inter / union if union > 0 else 0.0


def compute_drug_similarity_edges(
    drug2genes: Dict[str, Set[str]],
    min_sim: float = 0.0,
) -> List[Tuple[str, str, float]]:
    """
    TODO:
    - pentru toate perechile de medicamente (combinații de câte 2),
      calculați similaritatea Jaccard între seturile de gene
    - rețineți doar muchiile cu similaritate >= min_sim
    - returnați o listă de tuple (drug1, drug2, weight)
    """
    edges = []
    drugs = list(drug2genes.keys())
    
    for drug1, drug2 in itertools.combinations(drugs, 2):
        genes1 = drug2genes[drug1]
        genes2 = drug2genes[drug2]
        sim = jaccard_similarity(genes1, genes2)
        
        if sim >= min_sim:
            edges.append((drug1, drug2, sim))
    
    return edges


def edges_to_dataframe(edges: List[Tuple[str, str, float]]) -> pd.DataFrame:
    """
    TODO:
    - transformați lista de muchii (drug1, drug2, weight) într-un DataFrame
      cu coloanele: drug1, drug2, similarity
    """
    if not edges:
        return pd.DataFrame(columns=["drug1", "drug2", "similarity"])
    
    df = pd.DataFrame(edges, columns=["drug1", "drug2", "similarity"])
    return df


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # TODO 1: verificați că fișierul de input există
    print(f"[INFO] Checking input file: {DRUG_GENE_CSV}")
    ensure_exists(DRUG_GENE_CSV)

    # TODO 2: încărcați tabelul drug-gene
    print("[INFO] Loading drug-gene table...")
    df = load_drug_gene_table(DRUG_GENE_CSV)
    print(f"[INFO] Loaded {len(df)} drug-gene interactions")

    # TODO 3: construiți mapping-ul drug -> set de gene
    print("[INFO] Building drug2genes mapping...")
    drug2genes = build_drug2genes(df)
    print(f"[INFO] Found {len(drug2genes)} unique drugs")

    # TODO 4: construiți graful bipartit și salvați-l (opțional)
    print("[INFO] Building bipartite graph...")
    G = build_bipartite_graph(drug2genes)
    print(f"[INFO] Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    
    import pickle
    with open(OUT_GRAPH_DRUG_GENE, 'wb') as f:
        pickle.dump(G, f)
    print(f"[INFO] Saved graph to {OUT_GRAPH_DRUG_GENE}")

    # TODO 5: generați și salvați sumarul pe medicamente
    print("[INFO] Generating drug summary...")
    drug_summary = summarize_drugs(drug2genes)
    drug_summary.to_csv(OUT_DRUG_SUMMARY, index=False)
    print(f"[INFO] Saved drug summary to {OUT_DRUG_SUMMARY}")
    print(drug_summary.head())

    # TODO 6: calculați similaritatea între medicamente
    print("[INFO] Computing drug similarity edges...")
    similarity_edges = compute_drug_similarity_edges(drug2genes, min_sim=0.0)
    print(f"[INFO] Found {len(similarity_edges)} similarity edges")
    
    similarity_df = edges_to_dataframe(similarity_edges)
    similarity_df = similarity_df.sort_values("similarity", ascending=False)
    similarity_df.to_csv(OUT_DRUG_SIMILARITY, index=False)
    print(f"[INFO] Saved drug similarity to {OUT_DRUG_SIMILARITY}")
    print(similarity_df.head(10))

    print("\n[SUCCESS] Exercise 9.1 completed successfully!")
