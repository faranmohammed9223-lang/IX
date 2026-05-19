from pathlib import Path
import pandas as pd

# 1. File paths
DATA_FOLDER = Path(".")
BREAKTHROUGH_FILE = DATA_FOLDER / "breakthrough_data.xlsx"
COLUMN_PARAMETERS_FILE = DATA_FOLDER / "column_parameters.xlsx"
PFAS_PROPERTIES_FILE = DATA_FOLDER / "PFAS_properties.xlsx"

# 2. Breakthrough matrices
def list_breakthrough_matrices():
    excel_file = pd.ExcelFile(BREAKTHROUGH_FILE)
    return excel_file.sheet_names

def read_breakthrough_matrix(matrix_name):
    breakthrough_matrix = pd.read_excel(
        BREAKTHROUGH_FILE,
        sheet_name=matrix_name,
        engine="openpyxl"
    )

    breakthrough_matrix = breakthrough_matrix.rename(
        columns={breakthrough_matrix.columns[0]: "Bed Volumes"}
    )
    return breakthrough_matrix

# 3. PFAS properties

def read_pfas_properties_table():
    pfas_properties_table = pd.read_excel(
        PFAS_PROPERTIES_FILE,
        engine="openpyxl"
    )
    return pfas_properties_table


def get_pfas_properties_for_compound(pfas_properties_table, compound_name):
    property_name_map = {
        "MW": "Molecular Weight",
        "MolarVol": "Molar Volume",
        "Density": "Density",
        "Solubility": "Solubility",
        "VaporPress": "Vapor Pressure",
        "BP": "Boiling Point",
    }

    if compound_name not in pfas_properties_table.columns:
        return {}
    compound_properties = {}
    for _, row in pfas_properties_table.iterrows():
        raw_property_name = row["compound"]
        if raw_property_name == "compound":
            continue
        full_property_name = property_name_map.get(
            raw_property_name,
            raw_property_name
        )
        compound_properties[full_property_name] = row[compound_name]
    return compound_properties

# 4. Column parameters
def read_column_parameters():
    raw_column_parameters = pd.read_excel(
        COLUMN_PARAMETERS_FILE,
        sheet_name="column",
        engine="openpyxl",
        header=None
    )

    parameter_name_map = {
        "EBCT": "Empty Bed Contact Time",
        "dp": "Particle Diameter",
        "x": "Non-Constant Diffusivity Modifier",
        "D": "Diffusivity",
        "di": "Column Diameter",
        "A": "Cross-Sectional Area",
        "p": "Water Density",
        "u": "Water Viscosity",
        "V": "Volumetric Flow Rate",
        "HLR": "Hydraulic Loading Rate",
        "Sc": "Schmidt Number",
        "Re": "Reynolds Number",
        "RexSc": "Reynolds Number Times Schmidt Number",
        "BV": "Bed Volume",
        "h ": "Bed Depth",
        "di/dp": "Column Diameter to Particle Diameter Ratio",
    }

    column_parameters = {
        "Large Column": {},
        "Small Column": {},
        "Constraints": {},
        "Media Properties": {},
    }

    current_section = None

    for _, row in raw_column_parameters.iterrows():
        parameter_code = row[0]

        if pd.isna(parameter_code):
            continue
        if parameter_code == "LC":
            current_section = "Large Column"
            continue
        if parameter_code == "SC":
            current_section = "Small Column"
            continue
        if parameter_code == "Constraints:":
            current_section = "Constraints"
            continue
        if parameter_code == "Parameter":
            current_section = "Media Properties"
            continue
        if current_section in ["Large Column", "Small Column"]:
            full_parameter_name = parameter_name_map.get(
                parameter_code,
                parameter_code
            )
            column_parameters[current_section][full_parameter_name] = {
                "code": parameter_code,
                "value": row[1],
                "unit": row[2],
                "description": row[3],
                "equation": row[5] if len(row) > 5 else None,
            }
        elif current_section == "Constraints":
            column_parameters["Constraints"][parameter_code] = {
                "value": row[1],
                "description": row[2],
            }
        elif current_section == "Media Properties":
            column_parameters["Media Properties"][parameter_code] = {
                "unit": row[1],
                "average": row[2],
                "standard_deviation": row[3],
            }

    return column_parameters

def get_column_parameter(
    column_parameters,
    section_name,
    parameter_name,
    value_name="value"
):
    return column_parameters[section_name][parameter_name][value_name]



# 5. Prepare compounds from selected breakthrough matrix
def prepare_compounds_from_breakthrough_matrix(
    breakthrough_matrix,
    pfas_properties_table
):
    initial_row = breakthrough_matrix[
        breakthrough_matrix["Bed Volumes"] == "INITIAL"
    ]

    midpoint_row = breakthrough_matrix[
        breakthrough_matrix["Bed Volumes"] == "MIDPOINT"
    ]

    bed_volume_rows = breakthrough_matrix[
        pd.to_numeric(
            breakthrough_matrix["Bed Volumes"],
            errors="coerce"
        ).notna()
    ].copy()

    compounds_data = {}

    compound_names = breakthrough_matrix.columns[1:]

    for compound_name in compound_names:
        if pd.isna(compound_name):
            continue

        compound_breakthrough_data = pd.DataFrame()

        compound_breakthrough_data["Bed Volumes"] = bed_volume_rows[
            "Bed Volumes"
        ].values

        compound_breakthrough_data["Concentration"] = bed_volume_rows[
            compound_name
        ].values

        if not initial_row.empty:
            compound_breakthrough_data["Initial Concentration"] = initial_row[
                compound_name
            ].iloc[0]

        if not midpoint_row.empty:
            compound_breakthrough_data["Midpoint Concentration"] = midpoint_row[
                compound_name
            ].iloc[0]

        compound_pfas_properties = get_pfas_properties_for_compound(
            pfas_properties_table=pfas_properties_table,
            compound_name=compound_name
        )

        compounds_data[compound_name] = {
            "breakthrough_data": compound_breakthrough_data,
            "pfas_properties": compound_pfas_properties,
        }

    return compounds_data



# 6. Main reader
def read_model_inputs(matrix_name):
    breakthrough_matrix = read_breakthrough_matrix(matrix_name)
    pfas_properties_table = read_pfas_properties_table()
    column_parameters = read_column_parameters()

    compounds_data = prepare_compounds_from_breakthrough_matrix(
        breakthrough_matrix=breakthrough_matrix,
        pfas_properties_table=pfas_properties_table
    )

    model_inputs = {
        "matrix_name": matrix_name,
        "breakthrough_matrix": breakthrough_matrix,
        "compounds_data": compounds_data,
        "column_parameters": column_parameters,
        "pfas_properties_table": pfas_properties_table,
    }

    return model_inputs



# 7. Test the reader


if __name__ == "__main__":

    print("\nAvailable breakthrough matrices:")
    print(list_breakthrough_matrices())

    selected_matrix = list_breakthrough_matrices()[0]

    data = read_model_inputs(selected_matrix)

    print("\nSelected matrix:")
    print(data["matrix_name"])
    print("\nCompounds found in selected matrix:")
    print(list(data["compounds_data"].keys()))
    print("\nColumn parameter sections:")
    print(data["column_parameters"].keys())
    print("\nLarge Column parameters:")
    print(data["column_parameters"]["Large Column"].keys())
    print("\nSmall Column parameters:")
    print(data["column_parameters"]["Small Column"].keys())
    print("\nMedia Properties:")
    print(data["column_parameters"]["Media Properties"].keys())
    first_compound = list(data["compounds_data"].keys())[0]
    print(f"\nExample compound: {first_compound}")
    print("\nBreakthrough data:")
    print(data["compounds_data"][first_compound]["breakthrough_data"].head())
    print("\nPFAS properties:")
    print(data["compounds_data"][first_compound]["pfas_properties"])