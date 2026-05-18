from pathlib import Path
import matplotlib.pyplot as plt
from normalize_breakthrough import BreakthroughDataset
import os

filepath = (
        Path(__file__).resolve().parent.parent
        / "data_in"
        / "breakthrough_data.xlsx"
    )

dataset = BreakthroughDataset(filepath)

dataset.load()

for matrix in dataset.water_matrices.values():

    os.makedirs(matrix.name, exist_ok=True)

    for curve in matrix.curves.values():

        df = curve.to_dataframe()

        plt.figure()

        plt.plot(
            df["BVs"],
            df["C/Co"],
            marker="o"
        )

        plt.xlabel("Bed Volumes")
        plt.ylabel("C/Co")
        plt.title(f"{curve.name} - {matrix.name}")

        plt.savefig(f'./{matrix.name}/{curve.name}.png', bbox_inches='tight')
        plt.close()