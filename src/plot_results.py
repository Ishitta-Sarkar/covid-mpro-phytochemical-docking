import csv
import matplotlib.pyplot as plt


def load_docking_results(file_path):
    ligands = []
    scores = []

    with open(file_path, "r", encoding="utf-8") as file:
        reader = csv.DictReader(file)

        for row in reader:
            ligands.append(row["ligand_name"])
            scores.append(float(row["binding_affinity_kcal_mol"]))

    return ligands, scores


def plot_binding_affinity():
    ligands, scores = load_docking_results("data/docking_results.csv")

    plt.figure(figsize=(9, 5))
    plt.bar(ligands, scores)
    plt.xlabel("Ligands")
    plt.ylabel("Binding Affinity (kcal/mol)")
    plt.title("Docking Score Comparison of Selected Phytochemicals")
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig("results/binding_affinity_plot.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    plot_binding_affinity()
