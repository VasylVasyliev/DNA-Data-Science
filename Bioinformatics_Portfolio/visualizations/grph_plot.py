import matplotlib
matplotlib.use('Agg')  # Force non-interactive mode (no window will pop up)
import networkx as nx
import matplotlib.pyplot as plt
from Bio import SeqIO
import os

def draw_overlap_graph():
    input_file = "data/grph.txt"
    output_path = "data/plots/grph_graph.png"

    # 1. Load data
    if not os.path.exists(input_file):
        print(f"File {input_file} not found!")
        return

    records = list(SeqIO.parse(input_file, "fasta"))
    k = 3
    G = nx.DiGraph()
    
    # 2. Build graph
    for s_rec in records:
        G.add_node(s_rec.id)
        for t_rec in records:
            if s_rec.id != t_rec.id and s_rec.seq.endswith(t_rec.seq[:k]):
                G.add_edge(s_rec.id, t_rec.id)
    
    # 3. Create Plot
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, k=0.5)
    nx.draw(G, pos, with_labels=True, node_color='lightgreen', 
            edge_color='gray', node_size=1000, font_size=8, arrowsize=20)
    
    # 4. Save
    plt.savefig(output_path)
    print(f"--- SUCCESS! ---")
    print(f"File created: {output_path}")

if __name__ == "__main__":
    draw_overlap_graph()