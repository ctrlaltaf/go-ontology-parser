import csv
from pathlib import Path
import networkx as nx


def parse_obo_file(filepath):
    current = None
    in_terms_section = False
    terms = []
    # for each entry here are the possible keys. not all entries have these keys
    # {'relationship', 'is_obsolete', 'name', 'alt_id', 'synonym', 'xref'
    # 'consider', 'is_a', 'namespace', 'id', 'def', 'subset', 'comment', 'replaced_by'}

    # TODO: the relationship keys includes values such as "regulates GO:0006310 ! DNA recombination".
    # for now I am just using the the is_a key to build the hierarchy

    # parse the obo file and get only the [Term] entries
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line == "[Term]":
                in_terms_section = True
                if current:
                    terms.append(current)
                current = {}

            elif line.startswith("[") and line != "[Term]":
                in_terms_section = False
                if current:
                    terms.append(current)
                    current = None

            elif in_terms_section and current is not None:
                if ": " in line:
                    key, value = line.split(": ", 1)
                    current.setdefault(key, []).append(value)

    # store all unique keys
    keys = set()
    for term in terms:
        for key in term.keys():
            keys.add(key)
    return terms


def main():
    print("hierarchy parser")

    go_file = Path("data/go-basic.obo")
    terms = parse_obo_file(go_file)

    G = nx.DiGraph()

    for term in terms:
        go_id = term["id"][0]
        name = term["name"][0]
        namespace = term["namespace"][0]

        if 'is_obsolete' in term: # do not include obsolete terms. if there is a key, then value is always true
            continue

        # add new edges
        if go_id not in G.nodes():
            G.add_node(go_id)

        # add parent terms
        if "is_a" in term:
            for node in term["is_a"]:
                # print(node.split(" ")[0])
                parent_go_term = node.split(" ")[0]

                # if node doesnt exist then add
                if parent_go_term not in G.nodes():
                    G.add_node(parent_go_term)

                G.add_edge(
                    go_id,
                    parent_go_term,
                    name=name,
                    namespace=namespace,
                    is_direct_annotation=True,
                )

    print("number of go terms", len(G.nodes()))
    print("number of edges between go term nodes", len(G.edges()))

    output_go_network_path = Path("./output/go_hierarchy_edge_file.txt")

    with open(output_go_network_path, "w+") as f:
        headers = ["go_1", "go_2", "name", "namespace", "is_direct_annotation"]
        writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        for edge in G.edges(data=True):
            row = {
                "go_1": edge[0],
                "go_2": edge[1],
                "name": edge[2]['name'],
                "namespace": edge[2]['namespace'],
                "is_direct_annotation": edge[2]['is_direct_annotation'],
            }
            writer.writerow(row)
        f.close()


if __name__ == "__main__":
    main()
