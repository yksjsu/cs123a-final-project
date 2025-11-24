from pyvis.network import Network
import random
from time import perf_counter_ns
from biorust import linear_gap_alignment, affine_gap_alignment

from Bio import Align


def calculate_id(i: int, j: int, m: int) -> int:
    """
    Calculate id, given coordinates and
    number of cols ("m").
    """

    return i * m + j


match_score = 1
mismatch_penalty = -1
gap_open_penalty = -2
gap_extend_penalty = -1
is_local = False


random.seed(None)
upper_limit = 100
seq1 = "".join(random.choices(["A", "C", "T", "G"], k=random.randint(1, upper_limit)))
seq2 = "".join(random.choices(["A", "C", "T", "G"], k=random.randint(1, upper_limit)))

seq1 = "GAT"
seq2 = "GCCC"

start = perf_counter_ns()
data = affine_gap_alignment(seq1, seq2, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty, is_local)
end = perf_counter_ns()

score = data[-1][-1][0]

print("=== Rust ===")
print("Score:", score)
print("Runtime (ns):", end - start)
print()

# get number of rows and columns for graph
n = len(data)
m = len(data[0])

# create network
net = Network(height="800px", width="100%", directed=True, heading="Alignment")
net.toggle_physics(False)
id = 0

# create nodes and edges
for i in range(n):
    for j in range(m):
        score, directions = data[i][j]

        # add node with data
        title = "{} {}".format(seq1[i - 1], seq2[j - 1]) if i != 0 and j != 0 else "N/A"
        color = "#ffffff" if i == 0 or j == 0 else "#90e0ef"
        net.add_node(id, label="{:<10}".format(score), title=title, shape="box", x = j * 100, y = i * 100, physics=False, value=id, color=color)

        edge_color = "#333333" if i == 0 or j == 0 else "#90e0ef"

        # for each of the directions, add a directed edge from source
        # to target
        if directions[0]:
            net.add_edge(id, calculate_id(i, j - 1, m), physics=False, color=edge_color)
        
        if directions[1]:
            net.add_edge(id, calculate_id(i - 1, j - 1, m), physics=False, color=edge_color)
        
        if directions[2]:
            net.add_edge(id, calculate_id(i - 1, j, m), physics=False, color=edge_color)

        id += 1

# create the sequence characters above and to the left
# of the graph
for i in range(1, n):
    net.add_node(id, label=seq1[i - 1], x=-100, y = i * 100, color="#ffffff", shape="box", physics=False)
    id += 1

for j in range(1, m):
    net.add_node(id, label=seq2[j - 1], x = j * 100, y=-100, color="#ffffff", shape="box", physics=False)
    id += 1

net.write_html("example.html")