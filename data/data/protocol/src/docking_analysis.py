import csv


class DockingAnalysis:
    def __init__(self, ligand_file, result_file):
        self.ligands = self.load_csv(ligand_file)
        self.results = self.load_csv(result_file)

    def load_csv(self, file_path):
        records = []

        with open(file_path, "r", encoding="utf-8") as file:
            reader = csv.DictReader(file)

            for row in reader:
                records.append(row)

        return records

    def rank_ligands(self):
        ranked = sorted(
            self.results,
            key=lambda row: float(row["binding_affinity_kcal_mol"])
        )

        return ranked

    def best_candidate(self):
        return self.rank_ligands()[0]

    def generate_summary(self):
        ranked = self.rank_ligands()

        summary = "Virtual Screening Summary\n"
        summary += "=" * 40 + "\n\n"

        summary += "Ligand Ranking Based on Binding Affinity\n\n"

        for ligand in ranked:
            summary += (
                f"{ligand['ligand_name']} | "
                f"{ligand['binding_affinity_kcal_mol']} kcal/mol\n"
            )

        return summary
