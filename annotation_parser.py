import csv
from pathlib import Path


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
            if len(row) > 1 and row[0] == "UniProtKB": # just want UniProtKB. the rest are {'ComplexPortal', 'UniProtKB', 'RNAcentral'}
                entry = {}
                for i, value in enumerate(row):
                    entry[headers[i]] = value
                result.append(entry)

    return result


def main():
    print("annotation parser")

    go_annotation_path = Path("data/goa_human.gaf")

    data = get_go_annotation_data(go_annotation_path)

    # db_set = set()
    # for entry in data:
    #     db_set.add(entry['db'])
    # print(db_set)

if __name__ == "__main__":
    main()
