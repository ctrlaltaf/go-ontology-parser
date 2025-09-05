from collections import defaultdict
import networkx as nx
import csv
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

def get_go_annotation_data(filepath):
    result = []

    with open(filepath, "r") as f:
        file = csv.reader(f, delimiter="\t")

        result = []

        headers = [
            "db",  # 0. Source database (e.g., UniProtKB)
            "db_object_id",  # 1. Accession ID (e.g., P02750)
            "db_object_symbol",  # . Gene/protein symbol (e.g., LRG1)
            "qualifier",  # 3. Relation qualifier (enables, located_in, etc.)
            "go_id",  # 4. GO term ID
            "reference",  # 5. Supporting reference (PubMed, Reactome, GO_REF)
            "evidence_code",  # 6. Evidence code (IEA, IDA, NAS, etc.)
            "with_from",  # 7. With/From field (additional IDs, rules)
            "aspect",  # 8. Ontology aspect: P/F/C
            "db_object_name",  # 9. Full protein name
            "db_object_synonyms",  # 10. Synonyms
            "db_object_type",  # 11. Object type (protein, gene, etc.)
            "taxon",  # 12. Taxon ID (e.g., taxon:9606)
            "date",  # 13. Annotation date (YYYYMMDD)
            "assigned_by",  # 14. Source of annotation (UniProt, Reactome, Ensembl)
            "annotation_extension",  # 15. Extra context (e.g., part_of(UBERON:...))
            "gene_product_form_id",  # 16. Specific isoform/product ID
        ]

        for row in file:
            if (
                len(row) > 1 and row[0] == "UniProtKB"
            ):  # just want UniProtKB. the rest are {'ComplexPortal', 'UniProtKB', 'RNAcentral'}
                entry = {}
                for i, value in enumerate(row):
                    entry[headers[i]] = value
                result.append(entry)

    return result


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


def rank_go_terms(gene_to_terms):
    term_to_genes = defaultdict(set)
    max_depth = 0

    # Build reverse index: GO term -> set of genes
    for gene, terms in gene_to_terms.items():
        for go_term, depth in terms:
            term_to_genes[go_term].add(gene)
            max_depth = max(max_depth, depth)

    # Score each GO term
    ranked_terms = []
    for term, genes in term_to_genes.items():
        depth = max(
            d for g, terms in gene_to_terms.items() for (t, d) in terms if t == term
        )
        weight = depth / max_depth
        frequency = len(genes)  # how many genes have this term
        score = weight * frequency
        ranked_terms.append((term, score, depth, frequency))

    # Sort by score (highest first)
    ranked_terms.sort(key=lambda x: x[1], reverse=True)
    return ranked_terms


def main():

    # Get GO terms and the hierarchy

    go_file = Path("data/go-basic.obo")
    terms = parse_obo_file(go_file)
    G = (
        nx.DiGraph()
    )  # go hierarchy network. edges are directed. if there is a directed edge from A to B, the A is more specifc and B is less specific

    for term in terms:
        go_id = term["id"][0]
        name = term["name"][0]
        namespace = term["namespace"][0]

        if (
            "is_obsolete" in term
        ):  # do not include obsolete terms. if there is a key, then value is always true
            continue

        # add new edges
        if go_id not in G.nodes():
            G.add_node(go_id, name=name)
        elif (
            "name" not in G.nodes[go_id]
        ):  # handles the case when we add a node in the edges TODO probably theres a better way to handle this
            G.nodes[go_id]["name"] = name

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
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for edge in G.edges(data=True):
            row = {
                "go_1": edge[0],
                "go_2": edge[1],
                "name": edge[2]["name"],
                "namespace": edge[2]["namespace"],
                "is_direct_annotation": edge[2]["is_direct_annotation"],
            }
            writer.writerow(row)
        f.close()

    # print("is directed acyclic?",  nx.is_directed_acyclic_graph(G))

    # get depth information for each go term in the hierarchy

    depths = {}

    # Process nodes so parents come before children
    for node in nx.topological_sort(G.reverse()):
        preds = list(G.successors(node))  # successors in original graph
        if not preds:  # no children → top/general node
            depths[node] = 0
        else:
            depths[node] = 1 + max(depths[c] for c in preds)
    nx.set_node_attributes(G, depths, "depth")
    depth_data = []
    min_depth = float("inf")
    max_depth = -float("inf")
    depth_fig_path = Path("./output/go_depth_hist.pdf")

    for n, d in G.nodes(data=True):
        curr_depth = d["depth"]
        depth_data.append(curr_depth)
        if curr_depth < min_depth:
            min_depth = curr_depth
        if curr_depth > max_depth:
            max_depth = curr_depth

    plt.figure(figsize=(10, 8))
    plt.hist(depth_data, bins=max_depth + 1, alpha=0.7)
    plt.title("Histogram of depth in GO hierarchy")
    plt.xlabel("GO depth")
    plt.ylabel("Frequency")
    plt.savefig(depth_fig_path)
    # plt.show()
    plt.close()

    # get GO annotations for proteins

    go_annotation_path = Path("data/goa_human.gaf")
    data = get_go_annotation_data(go_annotation_path)
    name_id_mapper = {}
    P = (
        nx.DiGraph()
    )  # network with nodes as GO terms. TODO edges could be TF relationships from PAthwayNet

    min_go_annotation = float("inf")
    max_go_annotation = -float("inf")

    for entry in data:
        gene_id = entry["db_object_id"]
        gene_name = entry["db_object_symbol"]
        go_qualifier = entry["qualifier"]
        go_id = entry["go_id"]
        name_id_mapper[gene_name] = gene_id

        if gene_id not in P.nodes():
            P.add_node(
                gene_id,
                gene_name=gene_name,
                go_qualifier=go_qualifier,
                go_annotations_list=[go_id],
            )
        elif go_id not in P.nodes[gene_id]["go_annotations_list"]:
            P.nodes[gene_id]["go_annotations_list"].append(go_id)

        curr_len = len(P.nodes[gene_id]["go_annotations_list"])
        if curr_len < min_go_annotation:
            min_go_annotation = curr_len
        if curr_len > max_go_annotation:
            max_go_annotation = curr_len

    print(
        "min go annotation = ",
        min_go_annotation,
        "\tmax go annotation = ",
        max_go_annotation,
    )

    # make hist of number of go annotations for all the genes

    go_hist_data = []
    hist_output_path = Path("./output/go_annotations_hist.pdf")

    for node in P.nodes(data=True):
        go_annotation_count = len(node[1]["go_annotations_list"])
        go_hist_data.append(go_annotation_count)

    plt.figure(figsize=(10, 8))
    plt.hist(go_hist_data, bins=max_go_annotation + 1, edgecolor="black", alpha=0.7)
    plt.title("Histogram of GO annotations for all genes")
    plt.xlabel("GO annotations")
    plt.ylabel("Frequency")
    plt.savefig(hist_output_path)
    # plt.show()
    plt.close()

    # get go term information from set of genes

    gene_set = [
        "P31371",
        "P07996",
        "P42830",
        "P46695",
        "O15169",
        "Q9C004",
        "P42574",
        "P06756",
        "P17676",
        "Q15646",
        "P13164",
        "P09228",
        "P20823",
        "Q92985",
        "Q7Z7D3",
        "P20591",
        "P22674",
        "Q8IVU3",
    ]  # testing with top 18 genes from 2020 paper
    gene_go_data_dict = (
        {}
    )  # gene as the key and the value are the list of GO terms and their associated information

    for gene in gene_set:
        go_annotations_list = P.nodes[gene]["go_annotations_list"]

        # get indirect annotations

        complete_annotations_set = set()
        for go_term in go_annotations_list:
            indirect_annotations = nx.descendants(G, go_term)
            complete_annotations_set.update(indirect_annotations)
            # print(len(complete_annotations_set))
        # print(complete_annotations_set)

        go_depth_list = []  # contains (go_term, depth)
        for go_term in complete_annotations_set:
            depth = G.nodes[go_term]["depth"]
            go_depth_list.append((go_term, depth))

        go_depth_list.sort(key=lambda x: x[1])

        gene_go_data_dict[gene] = go_depth_list

    ranked = rank_go_terms(gene_go_data_dict)
    go_comparison_data  = []
    go_comparison_path = Path("./output/go_comparison.pdf")

    for term, score, depth, freq in ranked:
        go_comparison_data.append({"GO_term": f"{term} {G.nodes[term]["name"]}", "score" : score})

    df = pd.DataFrame(go_comparison_data)
    df_sorted = df.sort_values("score", ascending=False)
    df_top = df_sorted.head(20)

    plt.figure(figsize=(10, len(df_top) * 0.3))
    plt.barh(df_top["GO_term"], df_top["score"], color="skyblue")
    plt.xlabel("Score")
    plt.ylabel("GO Term")
    plt.title("GO Terms Ranked by Score")
    plt.gca().invert_yaxis()  # so the highest score is at the top
    plt.tight_layout()
    plt.savefig(go_comparison_path)
    plt.show()

if __name__ == "__main__":
    main()
