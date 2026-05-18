from pathlib import Path
import pandas as pd
import numpy as np


class PFASCurve:

    def __init__(self, name, bed_volumes, concentrations,
                 c0_initial=None, c0_midpoint=None):

        self.name = name
        self.bed_volumes = bed_volumes
        self.concentrations = concentrations

        self.c0_initial = c0_initial
        self.c0_midpoint = c0_midpoint
        self.c0_average = np.mean([self.c0_initial, self.c0_midpoint])
        self.c0_error = np.abs(self.c0_average - self.c0_midpoint)

    def normalize(self):
        """
        Normalize concentrations by initial feed concentration.
        """
        if self.c0_initial is None:
            raise ValueError("Initial concentration is missing.")

        self.normalized_concentrations = self.concentrations / self.c0_average

    def remove_nans(self):
        """
        Remove NaN values from curve data.
        """
        mask = ~(np.isnan(self.bed_volumes) |
                 np.isnan(self.concentrations))

        self.bed_volumes = self.bed_volumes[mask]
        self.concentrations = self.concentrations[mask]

    def has_no_breakthrough(self):
        """
        Remove curves that do not exceed 20% breakthrough
        """
        if len(self.normalized_concentrations) == 0:
            return True  # No data, treat as no breakthrough and remove
        return float(np.max(self.normalized_concentrations)) < 0.2
    
    def remove_outliers_zscore(self, window=10, threshold=1.5, verbose=False):
        if len(self.normalized_concentrations) == 0:
            return

        s = pd.Series(self.normalized_concentrations)
        rolling_mean = s.rolling(window, center=True, min_periods=1).mean()
        residuals = s - rolling_mean
        std = residuals.std()

        if std == 0:
            return

        z_scores = np.abs((residuals - residuals.mean()) / std)
        mask = (z_scores < threshold).to_numpy()

        if verbose:
            removed = np.sum(~mask)
            if removed > 0:
                print(f"  [{self.name}] Removed {removed} outlier(s) at BVs: "
                    f"{self.bed_volumes[~mask].tolist()}")

        self.bed_volumes = self.bed_volumes[mask]
        self.normalized_concentrations = self.normalized_concentrations[mask]
            
    def to_dataframe(self):
        return pd.DataFrame({
            "BVs": self.bed_volumes,
            "C/Co": self.normalized_concentrations
        })
    
    def summary(self):

        print(f"Compound: {self.name}")
        print(f"Points: {len(self.bed_volumes)}")
        print(f"Average feed: {self.c0_average:.3f} ± {self.c0_error:.3f} ng/mL")


class WaterMatrix:

    def __init__(self, name):

        self.name = name

        # Dictionary of PFAS curves
        self.curves = {}

    def add_curve(self, curve):

        self.curves[curve.name] = curve

    def get_curve(self, compound_name):

        return self.curves.get(compound_name)
    
    def filter_no_breakthrough(self):
        self.curves = {
            name: curve
            for name, curve in self.curves.items()
            if not curve.has_no_breakthrough()
        }
    
    def summary(self):

        print(f"Water Matrix: {self.name}")

        print("Compounds:")

        for name in self.curves:
            print(f" - {name}")


class BreakthroughDataset:

    def __init__(self, filepath):

        self.filepath = Path(filepath)

        self.water_matrices = {}

    def load(self):

        excel_data = pd.read_excel(
            self.filepath,
            sheet_name=None
        )

        for sheet_name, df in excel_data.items():

            matrix = WaterMatrix(sheet_name)

            bed_volumes = df.iloc[2:, 0].to_numpy(dtype=float)

            compound_columns = df.columns[1:]

            for col in compound_columns:

                concentrations = df.loc[2:, col].to_numpy(dtype=float)

                curve = PFASCurve(
                    name=col,
                    bed_volumes=bed_volumes.copy(),
                    concentrations=concentrations,
                    c0_initial=df.loc[0, col],
                    c0_midpoint=df.loc[1, col]
                )

                curve.remove_nans()

                matrix.add_curve(curve)
            
            for curve in matrix.curves.values():
                curve.normalize()

            matrix.filter_no_breakthrough()

            self.add_matrix(matrix)

    def add_matrix(self, matrix):

        self.water_matrices[matrix.name] = matrix

    def get_matrix(self, matrix_name):

        return self.water_matrices.get(matrix_name)

    def get_curve(self, matrix_name, compound_name):

        matrix = self.get_matrix(matrix_name)

        if matrix is None:
            return None

        return matrix.get_curve(compound_name)

    def summary(self):

        print("Dataset Summary")

        for matrix_name in self.water_matrices:

            print(f" - {matrix_name}")


if __name__ == "__main__":

    filepath = (
        Path(__file__).resolve().parent.parent
        / "data_in"
        / "breakthrough_data.xlsx"
    )

    dataset = BreakthroughDataset(filepath)

    dataset.load()

    dataset.summary()