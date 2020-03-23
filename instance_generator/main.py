import graph
import os.path
import argparse


def main():
    vertices = None
    seed = None
    p_general = 0.5
    p_clique = 1
    directory = "."

    parser = argparse.ArgumentParser(prog="Instance_Generator", description="Generate tree_depth instances")

    parser.add_argument("--vertices", "--v", type=int, help="The amount of vertices", required=True)
    parser.add_argument("--dir", "--d", default=".", help="The directory where the exported file should go")
    parser.add_argument("--gprob", "--g", default="0.5", type=float, help="The probability that allowed edges are added")
    parser.add_argument("--cprob", "--c", default="1", type=float, help="The probability that edges on the longest path are added")
    parser.add_argument("--seed", "--s", default=None, type=int, help="The seed for the random number generator")

    result = parser.parse_args()

    vertices = result.vertices
    if vertices <= 0:
        print("Vertex amount has to be positive!")
        return

    directory = result.dir

    p_general = result.gprob
    if p_general < 0 or p_general > 1:
        print("General probability has to be between 0 and 1")
        return

    p_clique = result.cprob
    if p_clique < 0 or p_clique > 1:
        print("Clique probability has to be between 0 and 1")
        return

    seed = result.seed
    if seed is not None and (seed < 0 or seed >= 2**31):
        print("Seed has to be between 0 and 2**31-1")
        return

    g = graph.Graph(vertices, p_general=p_general, p_clique=p_clique, seed=seed)
    g.permute()

    if p_clique >= 1:
        file_name = "{}/e{}_{}_{}_{}.gr".format(directory, g.vertices, g.edge_amount, g.tree_depth, g.seed)
    else:
        file_name = "{}/r{}_{}_{}_{}.gr".format(directory, g.vertices, g.edge_amount, g.tree_depth, g.seed)

    if os.path.isfile(file_name):
        print("ERROR: File {} already exists".format(file_name))
    else:
        f = open(file_name, "w")
        f.write(g.export())
        f.close()

        f = open(file_name[:-2] + "tree", "w")
        f.write(g.export_decomposition())
        f.close()

        print("Graph has been exported to {}, with {} vertices, {} edges and treedepth {}.".format(file_name, g.vertices, g.edge_amount, g.tree_depth))


if __name__ == "__main__":
    main()