import random
from time import perf_counter_ns

import pandas as pd
import matplotlib.pyplot as plt

import typer
from pyvis.network import Network
from Bio import Align

from biorust import linear_gap_alignment, affine_gap_alignment, reconstruct_alignment


# main app for holding the cli interface
app = typer.Typer()


def calculate_id(i: int, j: int, m: int) -> int:
    """
    Calculate id, given coordinates and
    number of cols ("m").
    """

    return i * m + j



def make_output_visual(matrix, seq1, seq2, file_output_name):
    # get number of rows and columns for graph
    n = len(matrix)
    m = len(matrix[0])

    # create network
    net = Network(height="800px", width="100%", directed=True, heading="Alignment")
    net.toggle_physics(False)
    id = 0

    # create nodes and edges
    for i in range(n):
        for j in range(m):
            score, directions = matrix[i][j]

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

    net.write_html(file_output_name)



@app.command(name="affalign")
def align_affine(seq1: str, seq2: str, file_output_name: str = "output.html", match_score: int = 1, mismatch_penalty: int = -1, gap_open_penalty: int = -2, gap_extend_penalty: int = -1, is_local: bool = True):
    """
    Command to generate a visualization for aligning two sequences, affine gap.
    """

    # measure performance of getting full scoring matrix
    # and alignment
    start = perf_counter_ns()
    score, starting_pos, matrix = affine_gap_alignment(seq1, seq2, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty, is_local)
    alignment = reconstruct_alignment(seq1, seq2, matrix, starting_pos, is_local)
    end = perf_counter_ns()

    # print basic information
    print("=== Rust ===")
    print("Score:", score)
    print("Runtime (ns):", end - start)
    print(alignment)
    print()

    make_output_visual(matrix, seq1, seq2, file_output_name)


@app.command(name="linalign")
def align_affine(seq1: str, seq2: str, file_output_name: str = "output.html", match_score: int = 1, mismatch_penalty: int = -1, linear_gap_penalty: int = -2, is_local: bool = True):
    """
    Command to generate a visualization for aligning two sequences, linear gap.
    """

    start = perf_counter_ns()
    score, starting_pos, matrix = linear_gap_alignment(seq1, seq2, match_score, mismatch_penalty, linear_gap_penalty, is_local)
    alignment = reconstruct_alignment(seq1, seq2, matrix, starting_pos, is_local)
    end = perf_counter_ns()

    print("=== Rust ===")
    print("Score:", score)
    print("Runtime (ns):", end - start)
    print(alignment)

    make_output_visual(matrix, seq1, seq2, file_output_name)


@app.command(name="testimpl")
def test_implementation():
    """
    Perform basic tests for implementation.
    """

    # set parameters
    trials = 1000
    match_score = 1
    mismatch_penalty = -1
    linear_gap_penalty = -2
    gap_open_penalty = -2
    gap_extend_penalty = -1
    is_local = True
    affine_gap = False

    # run tests
    for i in range(trials):
        upper_limit = 100

        # change whether to do local or global alignment
        is_local = not is_local

        # change whether to do affine or linear gap penalty every four times
        if (i % 4 == 0):
            affine_gap = not affine_gap
        
        # random sequence generation
        seq1 = "".join(random.choices(["A", "C", "T", "G"], k=random.randint(1, upper_limit)))
        seq2 = "".join(random.choices(["A", "C", "T", "G"], k=random.randint(1, upper_limit)))

        # get score from alignment for basic check
        score = 0
        if affine_gap:
            score, _, _ = affine_gap_alignment(seq1, seq2, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty, is_local)
        else:
            score, _, _ = linear_gap_alignment(seq1, seq2, match_score, mismatch_penalty, linear_gap_penalty, is_local)

        # construct Biopython aligner
        aligner = Align.PairwiseAligner()
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_penalty
        aligner.mode = "global" if not is_local else "local"

        if affine_gap:
            aligner.open_gap_score = gap_open_penalty
            aligner.extend_gap_score = gap_extend_penalty
        else:
            aligner.gap_score = linear_gap_penalty

        biopython_score = aligner.score(seq1, seq2)

        if biopython_score != score:
            print("Error! Incorrect implementation!")
            return

    print("All good!")


@app.command(name="compareperf")
def compare_performance(file_output_name: str = "output.csv"):
    """
    Compare performance over 100 trials, 25 each for
    global and local, linear and affine gap alignment
    over 100 and 1000-length sequences.
    """

    trials = 100
    match_score = 1
    mismatch_penalty = -1
    linear_gap_penalty = -2
    gap_open_penalty = -2
    gap_extend_penalty = -1
    is_local = True
    affine_gap = False

    # data to add later to dataframe
    perf_data = {
        "biopython": [],
        "rust": [],
        "seq1len": [],
        "seq2len": [],
        "islocal": [],
        "islinear": [],
    }

    # run tests
    for i in range(trials):
        # change whether to do local or global alignment
        is_local = not is_local

        # change whether to do affine or linear gap penalty every four times
        if ((i + 2) % 4 == 0):
            affine_gap = not affine_gap

        # random sequence generation
        seq1 = "".join(random.choices(["A", "C", "T", "G"], k=100))
        seq2 = "".join(random.choices(["A", "C", "T", "G"], k=100))


        matrix = None
        start = perf_counter_ns()
        if affine_gap:
            _, starting_pos, matrix = affine_gap_alignment(seq1, seq2, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty, is_local)
        else:
            _, starting_pos, matrix = linear_gap_alignment(seq1, seq2, match_score, mismatch_penalty, linear_gap_penalty, is_local)
        reconstruct_alignment(seq1, seq2, matrix, starting_pos, is_local)
        end = perf_counter_ns()

        rust_time = end - start

        # construct Biopython aligner
        start = perf_counter_ns()
        aligner = Align.PairwiseAligner()
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_penalty
        aligner.mode = "global" if not is_local else "local"

        if affine_gap:
            aligner.open_gap_score = gap_open_penalty
            aligner.extend_gap_score = gap_extend_penalty
        else:
            aligner.gap_score = linear_gap_penalty

        _ = aligner.align(seq1, seq2)
        end = perf_counter_ns()
        biopython_time = end - start

        perf_data["biopython"].append(biopython_time)
        perf_data["rust"].append(rust_time)
        perf_data["seq1len"].append(len(seq1))
        perf_data["seq2len"].append(len(seq2))
        perf_data["islocal"].append(is_local)
        perf_data["islinear"].append(not affine_gap)
    
    # run tests (for 1000)
    for i in range(trials):
        # change whether to do local or global alignment
        is_local = not is_local

        # change whether to do affine or linear gap penalty every four times
        if ((i + 2) % 4 == 0):
            affine_gap = not affine_gap

        # random sequence generation
        seq1 = "".join(random.choices(["A", "C", "T", "G"], k=1000))
        seq2 = "".join(random.choices(["A", "C", "T", "G"], k=1000))

        matrix = None
        start = perf_counter_ns()
        if affine_gap:
            _, starting_pos, matrix = affine_gap_alignment(seq1, seq2, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty, is_local)
        else:
            _, starting_pos, matrix = linear_gap_alignment(seq1, seq2, match_score, mismatch_penalty, linear_gap_penalty, is_local)
        reconstruct_alignment(seq1, seq2, matrix, starting_pos, is_local)
        end = perf_counter_ns()

        rust_time = end - start

        # construct Biopython aligner
        start = perf_counter_ns()
        aligner = Align.PairwiseAligner()
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_penalty
        aligner.mode = "global" if not is_local else "local"

        if affine_gap:
            aligner.open_gap_score = gap_open_penalty
            aligner.extend_gap_score = gap_extend_penalty
        else:
            aligner.gap_score = linear_gap_penalty

        _ = aligner.align(seq1, seq2)
        end = perf_counter_ns()
        biopython_time = end - start

        perf_data["biopython"].append(biopython_time)
        perf_data["rust"].append(rust_time)
        perf_data["seq1len"].append(len(seq1))
        perf_data["seq2len"].append(len(seq2))
        perf_data["islocal"].append(is_local)
        perf_data["islinear"].append(not affine_gap)

    df = pd.DataFrame(perf_data)
    df.to_csv(file_output_name)


@app.command(name="visualizecsv")
def visualize_csv(file_input_name: str = "output.csv"):
    """
    Visualize results from performance comparison.
    """

    df = pd.read_csv(file_input_name)

    # df from 100-length sequences
    hundred_df = df[df["seq1len"] == 100]
    thousand_df = df[df["seq2len"] == 1000]

    # get averages by rust vs biopython implementation
    # by type of alignment and length of sequence
    rust_avg_100_global_linear = hundred_df[(hundred_df["islocal"] == False) & (hundred_df["islinear"] == True)]["rust"].mean()
    rust_avg_100_local_linear = hundred_df[(hundred_df["islocal"] == True) & (hundred_df["islinear"] == True)]["rust"].mean()
    rust_avg_100_global_affine = hundred_df[(hundred_df["islocal"] == False) & (hundred_df["islinear"] == False)]["rust"].mean()
    rust_avg_100_local_affine = hundred_df[(hundred_df["islocal"] == True) & (hundred_df["islinear"] == False)]["rust"].mean()

    rust_avg_1000_global_linear = thousand_df[(thousand_df["islocal"] == False) & (thousand_df["islinear"] == True)]["rust"].mean()
    rust_avg_1000_local_linear = thousand_df[(thousand_df["islocal"] == True) & (thousand_df["islinear"] == True)]["rust"].mean()
    rust_avg_1000_global_affine = thousand_df[(thousand_df["islocal"] == False) & (thousand_df["islinear"] == False)]["rust"].mean()
    rust_avg_1000_local_affine = thousand_df[(thousand_df["islocal"] == True) & (thousand_df["islinear"] == False)]["rust"].mean()

    biopython_avg_100_global_linear = hundred_df[(hundred_df["islocal"] == False) & (hundred_df["islinear"] == True)]["biopython"].mean()
    biopython_avg_100_local_linear = hundred_df[(hundred_df["islocal"] == True) & (hundred_df["islinear"] == True)]["biopython"].mean()
    biopython_avg_100_global_affine = hundred_df[(hundred_df["islocal"] == False) & (hundred_df["islinear"] == False)]["biopython"].mean()
    biopython_avg_100_local_affine = hundred_df[(hundred_df["islocal"] == True) & (hundred_df["islinear"] == False)]["biopython"].mean()

    biopython_avg_1000_global_linear = thousand_df[(thousand_df["islocal"] == False) & (thousand_df["islinear"] == True)]["biopython"].mean()
    biopython_avg_1000_local_linear = thousand_df[(thousand_df["islocal"] == True) & (thousand_df["islinear"] == True)]["biopython"].mean()
    biopython_avg_1000_global_affine = thousand_df[(thousand_df["islocal"] == False) & (thousand_df["islinear"] == False)]["biopython"].mean()
    biopython_avg_1000_local_affine = thousand_df[(thousand_df["islocal"] == True) & (thousand_df["islinear"] == False)]["biopython"].mean()

    # compare average runtime from rust to biopython ratio
    categories = ["global, linear", "local, linear", "global, affine", "local, affine"]
    values = [
        rust_avg_100_global_linear,
        rust_avg_100_local_linear,
        rust_avg_100_global_affine,
        rust_avg_100_local_affine,
    ]

    values = [value / 1000000 for value in values]

    plt.bar(categories, values)
    plt.title("Average Runtime Ratio, Rust (100-seq)")
    plt.xlabel("Alignment Type")
    plt.ylabel("Rust Runtime (s)")
    plt.savefig("rust-100-sequence-comparison-rust.png")
    plt.cla()

    categories = ["global, linear", "local, linear", "global, affine", "local, affine"]
    values = [
        rust_avg_1000_global_linear,
        rust_avg_1000_local_linear,
        rust_avg_1000_global_affine,
        rust_avg_1000_local_affine,
    ]

    values = [value / 1000000 for value in values]

    plt.bar(categories, values)
    plt.title("Average Runtime Ratio, Rust (1000-seq)")
    plt.xlabel("Alignment Type")
    plt.ylabel("Rust Runtime (s)")
    plt.savefig("rust-1000-sequence-comparison.png")
    plt.cla()

    categories = ["global, linear", "local, linear", "global, affine", "local, affine"]
    values = [
        biopython_avg_100_global_linear,
        biopython_avg_100_local_linear,
        biopython_avg_100_global_affine,
        biopython_avg_100_local_affine,
    ]

    values = [value / 1000000 for value in values]

    plt.bar(categories, values)
    plt.title("Average Runtime Ratio, Biopython (100-seq)")
    plt.xlabel("Alignment Type")
    plt.ylabel("Biopython Runtime (s)")
    plt.savefig("biopython-100-sequence-comparison.png")
    plt.cla()

    categories = ["global, linear", "local, linear", "global, affine", "local, affine"]
    values = [
        biopython_avg_1000_global_linear,
        biopython_avg_1000_local_linear,
        biopython_avg_1000_global_affine,
        biopython_avg_1000_local_affine,
    ]

    values = [value / 1000000 for value in values]

    plt.bar(categories, values)
    plt.title("Average Runtime Ratio, Biopython (1000-seq)")
    plt.xlabel("Alignment Type")
    plt.ylabel("Biopython Runtime (s)")
    plt.savefig("biopython-1000-sequence-comparison.png")
    plt.cla()

app()