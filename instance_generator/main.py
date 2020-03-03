import graph
import os.path
import sys


def main():
    vertices = None
    seed = None
    p_general = 0.5
    p_clique = 1
    directory = "."

    if len(sys.argv) < 2:
        print("ERROR: Amount of vertices has to be specified.")
        print("Syntax: python main.py (vertices (int)) [export directory (default \".\")] [general probability (float, default 0.5)] [clique probability, default 1] [seed (int)]")
        return

    vertices = int(sys.argv[1])
    if vertices <= 0:
        print("Vertex amount has to be positive!")
        return
    if len(sys.argv) > 2:
        directory = sys.argv[2]
    if len(sys.argv) > 3:
        p_general = float(sys.argv[3])
        if p_general < 0 or p_general > 1:
            print("General probability has to be between 0 and 1")
            return
    if len(sys.argv) > 4:
        p_clique = float(sys.argv[4])
        if p_clique < 0 or p_clique > 1:
            print("Clique probability has to be between 0 and 1")
            return
    if len(sys.argv) > 5:
        seed = int(sys.argv[5])
        if seed < 0 or seed >= 2**31:
            print("Seed has to be between 0 and 2**31-1")
            return

    g = graph.Graph(vertices, p_general=p_general, p_clique=p_clique, seed=seed)

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

        print("Graph has been exported to {}, with {} vertices, {} edges and treedepth {}.".format(file_name, g.vertices, g.edge_amount, g.tree_depth))


if __name__ == "__main__":
    main()