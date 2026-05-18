from dataclasses import dataclass
from pathlib import Path
import re
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl import Workbook
from scipy.optimize import OptimizeWarning, curve_fit
from scipy.special import expit
from scipy.stats import t

from normalize_breakthrough import BreakthroughDataset, PFASCurve, WaterMatrix


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_PATH = SCRIPT_DIR.parent / "data_in" / "breakthrough_data.xlsx"
OUT_DIR = SCRIPT_DIR / "data_out"
FIT_DATA_PATH = SCRIPT_DIR / "fit_data.xlsx"
PARAMS_PATH = SCRIPT_DIR / "thomas_params.xlsx"
DUPLICATE_PARAMS_PATH = SCRIPT_DIR / "thomas_params_duplicates.xlsx"


@dataclass(frozen=True)
class ColumnParameters:
    sorbent_mass_g: float = 45.7 / 1000
    flow_rate_l_per_bv: float = 0.119 * 10**-3


@dataclass(frozen=True)
class ResidualOutlierSettings:
    enabled: bool = True
    min_points: int = 5
    breakthrough_threshold: float = 0.2
    min_breakthrough_points: int = 5
    residual_z_threshold: float = 3.0
    max_points_to_drop: int = 0
    max_normalized_concentration: float = 1.0
    min_rsquared: float = 0.6


@dataclass
class ScipyParameter:
    value: float
    stderr: float


@dataclass
class ScipyFitResult:
    params: dict
    covariance: np.ndarray
    residual: np.ndarray
    redchi: float


@dataclass
class ThomasFitResult:
    matrix_name: str
    compound_name: str
    result: object
    c0: float
    kTh: float
    kTh_lb: float
    kTh_ub: float
    qe: float
    qe_lb: float
    qe_ub: float
    BV20: float
    BV20_lb: float
    BV20_ub: float
    rsquared: float
    p_value_a: float
    p_value_b: float

    @property
    def a_fit(self):
        return self.result.params["a"].value

    @property
    def b_fit(self):
        return self.result.params["b"].value


class ThomasModel:
    @staticmethod
    def breakthrough(x, a, b):
        return expit(a * x - b)


class ThomasFitter:
    def __init__(self, column_parameters=None, outlier_settings=None):
        self.column_parameters = column_parameters or ColumnParameters()
        self.outlier_settings = outlier_settings or ResidualOutlierSettings()

    def fit_curve(
        self,
        matrix_name: str,
        curve: PFASCurve,
        bed_volumes=None,
        normalized=None,
    ) -> ThomasFitResult:
        if bed_volumes is None or normalized is None:
            bed_volumes, normalized = self._clean_curve_values(curve)

        if len(bed_volumes) < 3:
            raise ValueError(
                f"{matrix_name} - {curve.name} has fewer than 3 valid points."
            )

        x = bed_volumes / 1000
        y = normalized

        try:
            result = self.fit_values(bed_volumes, normalized)
        except Exception as exc:
            raise ValueError(f"SciPy curve_fit failed: {exc}") from exc
        ci = self._confidence_intervals(result)

        c0 = curve.c0_average
        if pd.isna(c0) or c0 == 0:
            raise ValueError(f"{matrix_name} - {curve.name} has invalid C0: {c0}")

        a_fit = result.params["a"].value
        b_fit = result.params["b"].value

        kTh = a_fit / c0
        kTh_lb = ci["a"]["lb"] / c0
        kTh_ub = ci["a"]["ub"] / c0

        qe = self._qe(b_fit, kTh)
        qe_lb = self._qe(ci["b"]["lb"], kTh)
        qe_ub = self._qe(ci["b"]["ub"], kTh)

        bv20 = self._bv20(a_fit, b_fit)
        bv20_lb = self._bv20(ci["a"]["ub"], ci["b"]["lb"])
        bv20_ub = self._bv20(ci["a"]["lb"], ci["b"]["ub"])

        p_value_a = self._p_value(result.params["a"].value, result.params["a"].stderr, len(x))
        p_value_b = self._p_value(result.params["b"].value, result.params["b"].stderr, len(x))

        return ThomasFitResult(
            matrix_name=matrix_name,
            compound_name=curve.name,
            result=result,
            c0=c0,
            kTh=kTh,
            kTh_lb=kTh_lb,
            kTh_ub=kTh_ub,
            qe=qe,
            qe_lb=qe_lb,
            qe_ub=qe_ub,
            BV20=bv20,
            BV20_lb=bv20_lb,
            BV20_ub=bv20_ub,
            rsquared=self._rsquared(result, y),
            p_value_a=p_value_a,
            p_value_b=p_value_b,
        )

    def fit_values(self, bed_volumes, normalized):
        x = np.asarray(bed_volumes, dtype=float) / 1000
        y = np.asarray(normalized, dtype=float)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            popt, pcov = curve_fit(
                ThomasModel.breakthrough,
                x,
                y,
                p0=[0.1, 5],
                maxfev=10000,
            )

        fitted = ThomasModel.breakthrough(x, *popt)
        residual = y - fitted
        degrees_of_freedom = max(len(y) - len(popt), 1)
        redchi = float(np.sum(residual**2) / degrees_of_freedom)
        stderr = self._standard_errors(pcov, len(popt))

        return ScipyFitResult(
            params={
                "a": ScipyParameter(value=float(popt[0]), stderr=stderr[0]),
                "b": ScipyParameter(value=float(popt[1]), stderr=stderr[1]),
            },
            covariance=pcov,
            residual=residual,
            redchi=redchi,
        )

    def fitted_curve(self, fit_result: ThomasFitResult, x_values):
        return ThomasModel.breakthrough(
            x_values,
            fit_result.a_fit,
            fit_result.b_fit,
        )

    def _qe(self, b_value, kTh):
        params = self.column_parameters
        return b_value * params.flow_rate_l_per_bv / (kTh * params.sorbent_mass_g)

    def _bv20(self, a_value, b_value):
        if not np.isfinite(a_value) or not np.isfinite(b_value) or abs(a_value) < 1e-12:
            return np.nan
        return float(((b_value - np.log(4)) / a_value) * 1000)

    def _clean_curve_values(self, curve: PFASCurve):
        bed_volumes = np.asarray(curve.bed_volumes, dtype=float)
        normalized = np.asarray(curve.normalized_concentrations, dtype=float)
        mask = ~(np.isnan(bed_volumes) | np.isnan(normalized))
        return bed_volumes[mask], normalized[mask]

    def _confidence_intervals(self, result):
        return self._stderr_intervals(result)

    def _standard_errors(self, covariance, n_parameters):
        if covariance is None or covariance.shape != (n_parameters, n_parameters):
            return [np.nan] * n_parameters

        diagonal = np.diag(covariance)
        if np.any(diagonal < 0):
            return [np.nan] * n_parameters

        return [float(stderr) for stderr in np.sqrt(diagonal)]

    def _stderr_intervals(self, result):
        intervals = {}

        for param_name in ("a", "b"):
            param = result.params[param_name]
            stderr = param.stderr

            if stderr is None or pd.isna(stderr) or not np.isfinite(stderr):
                lb = param.value
                ub = param.value
            else:
                lb = param.value - 1.96 * stderr
                ub = param.value + 1.96 * stderr

            intervals[param_name] = {"lb": lb, "ub": ub}

        return intervals

    def _p_value(self, value, stderr, n_observations):
        if stderr is None or stderr == 0 or pd.isna(stderr) or not np.isfinite(stderr):
            return np.nan

        degrees_of_freedom = n_observations - 2
        if degrees_of_freedom <= 0:
            return np.nan

        t_stat = value / stderr
        return 2 * (1 - t.cdf(abs(t_stat), degrees_of_freedom))

    def _rsquared(self, result, observed):
        total_sum_squares = np.sum((observed - np.mean(observed)) ** 2)
        if total_sum_squares == 0:
            return np.nan

        residual_sum_squares = np.sum(result.residual**2)
        return 1 - residual_sum_squares / total_sum_squares


@dataclass
class DroppedPoint:
    matrix_name: str
    compound_name: str
    bed_volume: float
    normalized_concentration: float
    fitted_concentration: float
    residual: float
    residual_z: float


class ResidualOutlierFilter:
    def __init__(self, fitter: ThomasFitter, settings=None):
        self.fitter = fitter
        self.settings = settings or ResidualOutlierSettings()

    def filter_curve(self, matrix_name: str, curve: PFASCurve, bed_volumes=None, normalized=None):
        if bed_volumes is None or normalized is None:
            bed_volumes, normalized = self.fitter._clean_curve_values(curve)

        dropped_points = []

        if not self.settings.enabled:
            return bed_volumes, normalized, dropped_points

        for _ in range(self.settings.max_points_to_drop):
            if len(bed_volumes) <= self.settings.min_points:
                break

            try:
                result = self._fit_values(bed_volumes, normalized)
            except Exception:
                break

            x = bed_volumes / 1000
            fitted = ThomasModel.breakthrough(
                x,
                result.params["a"].value,
                result.params["b"].value,
            )
            residuals = normalized - fitted
            residual_z = self._robust_z_scores(residuals)
            drop_index = int(np.argmax(residual_z))

            if residuals[drop_index] <= 0:
                break

            if residual_z[drop_index] < self.settings.residual_z_threshold:
                break

            dropped_points.append(
                DroppedPoint(
                    matrix_name=matrix_name,
                    compound_name=curve.name,
                    bed_volume=bed_volumes[drop_index],
                    normalized_concentration=normalized[drop_index],
                    fitted_concentration=fitted[drop_index],
                    residual=residuals[drop_index],
                    residual_z=residual_z[drop_index],
                )
            )

            keep_mask = np.ones(len(bed_volumes), dtype=bool)
            keep_mask[drop_index] = False
            bed_volumes = bed_volumes[keep_mask]
            normalized = normalized[keep_mask]

        return bed_volumes, normalized, dropped_points

    def has_breakthrough(self, normalized):
        if len(normalized) == 0:
            return False
        breakthrough_points = normalized >= self.settings.breakthrough_threshold
        return int(np.sum(breakthrough_points)) >= self.settings.min_breakthrough_points

    def fitted_curve_is_decreasing(self, bed_volumes, normalized):
        if len(bed_volumes) < self.settings.min_points:
            return False

        try:
            result = self._fit_values(bed_volumes, normalized)
        except Exception:
            return False

        return result.params["a"].value < 0

    def cap_normalized_values(self, normalized):
        capped = np.minimum(normalized, self.settings.max_normalized_concentration)
        changed_mask = capped != normalized
        return capped, changed_mask

    def _fit_values(self, bed_volumes, normalized):
        return self.fitter.fit_values(bed_volumes, normalized)

    def _robust_z_scores(self, residuals):
        median = np.median(residuals)
        mad = np.median(np.abs(residuals - median))

        if mad == 0:
            std = np.std(residuals)
            if std == 0:
                return np.zeros(len(residuals))
            return np.abs((residuals - np.mean(residuals)) / std)

        return np.abs(0.6745 * (residuals - median) / mad)


class ThomasAnalysis:
    y_column_groups = {
        "PFCAs": ["PFBA", "PFPeA", "PFHxA", "PFHpA", "PFOA", "PFNA", "PFDA", "PFDoA", "PFUnA"],
        "PFSAs and FTS": ["PFBS", "PFHxS", "PFOS", "4:2 FTS", "6:2 FTS", "8:2 FTS"],
        "PFEAs": ["GenX", "PMPA", "PFMOAA", "NBP 1", "NBP 2-1", "NBP 2-2"],
    }

    def __init__(self, data_path=DATA_PATH, output_dir=OUT_DIR):
        self.data_path = Path(data_path)
        self.output_dir = Path(output_dir)
        self.dataset = BreakthroughDataset(self.data_path)
        self.settings = ResidualOutlierSettings()
        self.fitter = ThomasFitter(outlier_settings=self.settings)
        self.outlier_filter = ResidualOutlierFilter(self.fitter, settings=self.settings)
        self.fit_results = []
        self.filtered_curves = {}
        self.dropped_points = []

    def run(self):
        self.output_dir.mkdir(exist_ok=True)
        self.dataset.load()

        for matrix in self.dataset.water_matrices.values():
            self._fit_matrix(matrix)
            self._plot_matrix(matrix)

        self._print_summary()
        self.export_fitted_curves(FIT_DATA_PATH)
        self.export_parameters(PARAMS_PATH)
        self.export_duplicate_parameters(DUPLICATE_PARAMS_PATH)

    def _fit_matrix(self, matrix: WaterMatrix):
        for curve in matrix.curves.values():
            bed_volumes, normalized = self.fitter._clean_curve_values(curve)
            original_normalized = normalized.copy()
            normalized, capped_mask = self.outlier_filter.cap_normalized_values(normalized)

            capped_count = int(np.sum(capped_mask))
            if capped_count > 0:
                print(
                    f"Capped {matrix.name} - {curve.name}: "
                    f"{capped_count} point(s) with C/C0 > "
                    f"{self.outlier_filter.settings.max_normalized_concentration:g}; "
                    f"max original C/C0={np.max(original_normalized[capped_mask]):g} -> "
                    f"{self.outlier_filter.settings.max_normalized_concentration:g}"
                )

            bed_volumes, normalized, dropped_points = self.outlier_filter.filter_curve(
                matrix.name,
                curve,
                bed_volumes=bed_volumes,
                normalized=normalized,
            )
            self.dropped_points.extend(dropped_points)

            for dropped in dropped_points:
                print(
                    "Dropped high-residual point "
                    f"{dropped.matrix_name} - {dropped.compound_name}: "
                    f"BV={dropped.bed_volume:g}, "
                    f"C/C0={dropped.normalized_concentration:g}, "
                    f"fit={dropped.fitted_concentration:g}, "
                    f"residual={dropped.residual:g}, "
                    f"robust_z={dropped.residual_z:g}"
                )

            if not self.outlier_filter.has_breakthrough(normalized):
                print(
                    f"Skipping {matrix.name} - {curve.name}: "
                    "fewer than 3 breakthrough points after capping and residual cleanup."
                )
                continue

            try:
                fit_result = self.fitter.fit_curve(
                    matrix.name,
                    curve,
                    bed_volumes=bed_volumes,
                    normalized=normalized,
                )
            except ValueError as exc:
                print(f"Skipping {matrix.name} - {curve.name}: {exc}")
                continue

            if fit_result.rsquared < self.outlier_filter.settings.min_rsquared:
                print(
                    f"Skipping {matrix.name} - {curve.name}: "
                    f"R-squared {fit_result.rsquared:g} < "
                    f"{self.outlier_filter.settings.min_rsquared:g} after final fit."
                )
                continue

            self.filtered_curves[(matrix.name, curve.name)] = (bed_volumes, normalized)
            self.fit_results.append(fit_result)

    def _plot_matrix(self, matrix: WaterMatrix):
        for group_name, compound_names in self.y_column_groups.items():
            group_results = [
                fit
                for fit in self.fit_results
                if fit.matrix_name == matrix.name and fit.compound_name in compound_names
            ]

            if not group_results:
                continue

            fig, ax = plt.subplots(1, 1, figsize=(8, 6))
            ax.set_title(f"Thomas model fit - {matrix.name} - {group_name}", fontsize=18)

            marker_styles = ["o", "s", "^", "D", "v", "<", ">"]
            x_values = np.linspace(-50, 250, 120)

            for i, fit_result in enumerate(group_results):
                bed_volumes, normalized = self.filtered_curves[
                    (fit_result.matrix_name, fit_result.compound_name)
                ]
                x = bed_volumes / 1000
                marker = marker_styles[i % len(marker_styles)]

                scatter = ax.scatter(
                    x,
                    normalized,
                    label=fit_result.compound_name,
                    marker=marker,
                )
                ax.plot(
                    x_values,
                    self.fitter.fitted_curve(fit_result, x_values),
                    color=scatter.get_facecolor()[0],
                    linewidth=1,
                )

            ax.legend(loc="upper right", bbox_to_anchor=(1, 1), fontsize=9)
            ax.set_xlabel("Bed Volumes (x1000)", fontsize=16)
            ax.set_ylabel("C/C0", fontsize=16)
            ax.set_ylim(0, 1.5)

            fig.savefig(
                self.output_dir / f"{matrix.name}_{group_name}.png",
                bbox_inches="tight",
            )
            plt.close(fig)

    def export_fitted_curves(self, filepath):
        x_values = np.linspace(-25, 300, 1500)

        with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
            for matrix_name in self.dataset.water_matrices:
                sheet_results = [
                    fit for fit in self.fit_results if fit.matrix_name == matrix_name
                ]

                df_fit = pd.DataFrame({"x_values": x_values})
                for fit_result in sheet_results:
                    df_fit[fit_result.compound_name] = self.fitter.fitted_curve(
                        fit_result,
                        x_values,
                    )

                df_fit.to_excel(writer, sheet_name=matrix_name, index=False)

    def export_parameters(self, filepath):
        wb = Workbook()
        ws = wb.active
        ws.title = "Results"

        fieldnames = [
            "WQ",
            "PFAS",
            "kTh",
            "kTh_lb",
            "kTh_ub",
            "qe",
            "qe_lb",
            "qe_ub",
            "BV20",
            "BV20_lb",
            "BV20_ub",
            "rsquared",
            "a_pval",
            "b_pval",
            "Co",
            "n_replicates",
        ]
        ws.append(fieldnames)

        for fit_group in self._group_replicate_results():
            ws.append(
                [
                    fit_group["matrix_name"],
                    fit_group["compound_name"],
                    fit_group["kTh"],
                    fit_group["kTh_lb"],
                    fit_group["kTh_ub"],
                    fit_group["qe"],
                    fit_group["qe_lb"],
                    fit_group["qe_ub"],
                    fit_group["BV20"],
                    fit_group["BV20_lb"],
                    fit_group["BV20_ub"],
                    fit_group["rsquared"],
                    fit_group["a_pval"],
                    fit_group["b_pval"],
                    fit_group["Co"],
                    fit_group["n_replicates"],
                ]
            )

        wb.save(filepath)

    def export_duplicate_parameters(self, filepath):
        wb = Workbook()
        ws = wb.active
        ws.title = "Duplicate Results"

        fieldnames = [
            "WQ",
            "WQ_base",
            "replicate",
            "PFAS",
            "kTh",
            "kTh_lb",
            "kTh_ub",
            "qe",
            "qe_lb",
            "qe_ub",
            "BV20",
            "BV20_lb",
            "BV20_ub",
            "rsquared",
            "a_pval",
            "b_pval",
            "Co",
        ]
        ws.append(fieldnames)

        for fit in sorted(self.fit_results, key=lambda item: (item.matrix_name, item.compound_name)):
            ws.append(
                [
                    fit.matrix_name,
                    self._replicate_base_name(fit.matrix_name),
                    self._replicate_number(fit.matrix_name),
                    fit.compound_name,
                    fit.kTh,
                    fit.kTh_lb,
                    fit.kTh_ub,
                    fit.qe,
                    fit.qe_lb,
                    fit.qe_ub,
                    fit.BV20,
                    fit.BV20_lb,
                    fit.BV20_ub,
                    fit.rsquared,
                    fit.p_value_a,
                    fit.p_value_b,
                    fit.c0,
                ]
            )

        wb.save(filepath)

    def _group_replicate_results(self):
        grouped = {}

        for fit in self.fit_results:
            group_key = (
                self._replicate_base_name(fit.matrix_name),
                fit.compound_name,
            )
            grouped.setdefault(group_key, []).append(fit)

        rows = []
        for (matrix_name, compound_name), fits in grouped.items():
            rows.append(
                {
                    "matrix_name": matrix_name,
                    "compound_name": compound_name,
                    **self._group_summary(fits),
                }
            )

        return sorted(rows, key=lambda row: (row["matrix_name"], row["compound_name"]))

    def _replicate_base_name(self, matrix_name):
        return re.sub(r"-\d+$", "", matrix_name)

    def _replicate_number(self, matrix_name):
        match = re.search(r"-(\d+)$", matrix_name)
        if match is None:
            return np.nan
        return int(match.group(1))

    def _group_summary(self, fits):
        kTh_mean, kTh_lb, kTh_ub = self._positive_mean_ci_bounds(
            [fit.kTh for fit in fits],
            [(fit.kTh_lb, fit.kTh_ub) for fit in fits],
        )
        qe_mean, qe_lb, qe_ub = self._positive_mean_ci_bounds(
            [fit.qe for fit in fits],
            [(fit.qe_lb, fit.qe_ub) for fit in fits],
        )
        bv20_mean, bv20_lb, bv20_ub = self._mean_ci_bounds(
            [fit.BV20 for fit in fits],
            [(fit.BV20_lb, fit.BV20_ub) for fit in fits],
        )

        return {
            "kTh": kTh_mean,
            "kTh_lb": kTh_lb,
            "kTh_ub": kTh_ub,
            "qe": qe_mean,
            "qe_lb": qe_lb,
            "qe_ub": qe_ub,
            "BV20": bv20_mean,
            "BV20_lb": bv20_lb,
            "BV20_ub": bv20_ub,
            "rsquared": self._finite_mean([fit.rsquared for fit in fits]),
            "a_pval": self._finite_mean([fit.p_value_a for fit in fits]),
            "b_pval": self._finite_mean([fit.p_value_b for fit in fits]),
            "Co": self._finite_mean([fit.c0 for fit in fits]),
            "n_replicates": len(fits),
        }

    def _mean_ci_bounds(self, values, single_fit_bounds):
        finite_pairs = [
            (value, bounds)
            for value, bounds in zip(values, single_fit_bounds)
            if np.isfinite(value)
        ]
        finite_values = np.asarray([value for value, _ in finite_pairs], dtype=float)

        if len(finite_values) == 0:
            return np.nan, np.nan, np.nan

        mean = float(np.mean(finite_values))

        if len(finite_values) > 1:
            sem = float(np.std(finite_values, ddof=1) / np.sqrt(len(finite_values)))
            ci_half_width = float(t.ppf(0.975, len(finite_values) - 1) * sem)
            return mean, mean - ci_half_width, mean + ci_half_width

        _, (lower_bound, upper_bound) = finite_pairs[0]
        lower_bound = lower_bound if np.isfinite(lower_bound) else np.nan
        upper_bound = upper_bound if np.isfinite(upper_bound) else np.nan
        return mean, lower_bound, upper_bound

    def _positive_mean_ci_bounds(self, values, single_fit_bounds):
        finite_pairs = [
            (value, bounds)
            for value, bounds in zip(values, single_fit_bounds)
            if np.isfinite(value) and value > 0
        ]
        finite_values = np.asarray([value for value, _ in finite_pairs], dtype=float)

        if len(finite_values) == 0:
            return np.nan, np.nan, np.nan

        mean = float(np.mean(finite_values))

        if len(finite_values) > 1:
            log_values = np.log(finite_values)
            log_sem = float(np.std(log_values, ddof=1) / np.sqrt(len(log_values)))
            multiplier = float(np.exp(t.ppf(0.975, len(log_values) - 1) * log_sem))
            return mean, mean / multiplier, mean * multiplier

        _, (lower_bound, upper_bound) = finite_pairs[0]
        lower_bound = lower_bound if np.isfinite(lower_bound) else np.nan
        upper_bound = upper_bound if np.isfinite(upper_bound) else np.nan

        if np.isfinite(lower_bound):
            lower_bound = max(lower_bound, 0)
        if np.isfinite(upper_bound):
            upper_bound = max(upper_bound, mean)

        return mean, lower_bound, upper_bound

    def _finite_mean(self, values):
        finite_values = [value for value in values if np.isfinite(value)]
        if not finite_values:
            return np.nan
        return float(np.mean(finite_values))

    def _print_summary(self):
        for fit in self.fit_results:
            print(f"Results for {fit.matrix_name}_{fit.compound_name}:")
            print(f"kTh: {fit.kTh} [{fit.kTh_lb}, {fit.kTh_ub}]")
            print(f"qe: {fit.qe} [{fit.qe_lb}, {fit.qe_ub}]")
            print(f"R-squared: {fit.rsquared}")
            print()


if __name__ == "__main__":
    analysis = ThomasAnalysis()
    analysis.run()
