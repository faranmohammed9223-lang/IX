#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
from scipy.stats import t 
from scipy.optimize import fsolve
from openpyxl import Workbook
from pathlib import Path
import os

# Read all sheets into a dictionary of DataFrames
DATA_PATH = Path("data_cleaned.xlsx")  
all_experiments = pd.read_excel(DATA_PATH, sheet_name=None)

#create output directory
out_dir = "data_out"
os.makedirs(out_dir, exist_ok=True)

#RSSCT Parameters
x_s = 45.7/1000 #mass of sorbent (g) 
Q = 0.119*10**-3 #flow rate (L/BV)

# Access each DataFrame from the dictionary
for sheet_name, df_sheet in all_experiments.items():
    all_experiments[sheet_name] = df_sheet.rename(columns={df_sheet.columns[0]: 'BV'})
    
# Create a dictionary to store regression results
regression_results = {}

# Define the equation to be solved
def equation_bv20(x, a, b):
    return logistic_function(x, a, b) - 0.2

# Logistic function
def logistic_function(x, a, b):
    return 1 / (1 + np.exp(b - a * x))

# Columns to exclude for each sheet
columns_to_exclude = {
    'NaCl':['GenX', 'PFHxA', 'PFHpA','PFOA','PFNA','PFDA','PFDoA','PFUnA','PFBS', 'PFHxS','PFOS','6:2 FTS','8:2 FTS','NBP 1','NBP 2-1','NBP 2-2'],
    'NaCl+SO4':['PFHpA','PFOA','PFNA','PFDA','PFDoA','PFUnA','PFBS', 'PFHxS','PFOS','6:2 FTS','8:2 FTS','NBP 1','NBP 2-1','NBP 2-2'],
    'NaCl+SO4 (0.57 mM)':['PFHpA','PFOA','PFNA','PFDA','PFDoA','PFUnA','PFBS', 'PFHxS','PFOS','6:2 FTS','8:2 FTS','NBP 1','NBP 2-1','NBP 2-2'],
    'NaCl+NO3':['PFOA','PFNA','PFDA','PFDoA','PFUnA','PFBS', 'PFHxS','PFOS', '8:2 FTS','NBP 1','NBP 2-1','NBP 2-2'],
    'NaCl+NOM1':['PFOA','PFNA','PFDA','PFDoA','PFUnA','PFBS', 'PFHxS','PFOS','8:2 FTS','NBP 1','NBP 2-1','NBP 2-2'],
    'NaCl+NOM5':['PFDoA','PFUnA'],
    'NaCl+NOM1+SO4':['PFNA','PFDA','PFDoA','PFUnA','PFBS', 'PFHxS','PFOS','8:2 FTS','NBP 1','NBP 2-1','NBP 2-2'],
    'NaCl+NOM5+SO4':['PFDoA','PFUnA'],
    'NaCl+NOM1+SO4 (0.57 mM)':['PFNA','PFDA','PFDoA','PFUnA','PFBS', 'PFHxS','PFOS','8:2 FTS','NBP 1','NBP 2-1','NBP 2-2'],
    'NaCl+NOM1+NO3':['PFNA','PFDA','PFDoA','PFUnA','PFBS', 'PFHxS','PFOS','8:2 FTS','NBP 1','NBP 2-1','NBP 2-2'],
    'NaCl+SO4+NO3':['PFNA','PFDA','PFDoA','PFUnA','PFBS', 'PFHxS','PFOS','8:2 FTS','NBP 1','NBP 2-1','NBP 2-2'],
    'NaCl+NOM1+SO4+NO3':['PFDA','PFDoA','PFUnA','PFHxS','PFOS','8:2 FTS','NBP 1','NBP 2-1','NBP 2-2']}

# Define the groups of y columns you want to plot together
y_column_groups = {
    'PFCAs': ['PFBA', 'PFPeA','PFHxA','PFHpA','PFOA','PFNA','PFDA','PFDoA','PFUnA'],
    'PFSAs and FTS': ['PFBS', 'PFHxS','PFOS','4:2 FTS','6:2 FTS','8:2 FTS'],
    'PFEAs': ['GenX', 'PMPA','PFMOAA','NBP 1','NBP 2-1','NBP 2-2']
}

# Iterate over all sheets
for sheet_name in list(all_experiments.keys())[-12:]:
    df_sheet_cleaned = all_experiments[sheet_name].iloc[:-1]
    
    last_row_values = all_experiments[sheet_name].iloc[-1]
    
    # Remove columns with no breakthrough
    if sheet_name in columns_to_exclude:
        df_sheet_cleaned = df_sheet_cleaned.drop(columns=columns_to_exclude[sheet_name], errors='ignore')

    # Get the first column as x
    x_column_name = df_sheet_cleaned.columns[0]

    for group_name, y_columns in y_column_groups.items():
        # Remove columns not present in the cleaned DataFrame
        y_columns = [col for col in y_columns if col in df_sheet_cleaned.columns]

        # Skip the group if no columns are present after cleaning
        if not y_columns:
            continue

        fig, axs = plt.subplots(1, 1, figsize=(8, 6))
        axs.set_title(f'Thomas model fit - {sheet_name} - {group_name}', fontsize=18)

        for i, y_column_name in enumerate(y_columns):
            # Check if the column exists in the cleaned DataFrame
            if y_column_name not in df_sheet_cleaned.columns:
                continue

            # Remove rows with NaN values for the current y column
            df_cleaned_y = df_sheet_cleaned[[x_column_name, y_column_name]].copy()

            df_cleaned_y[x_column_name] = pd.to_numeric(df_cleaned_y[x_column_name], errors='coerce')
            df_cleaned_y[y_column_name] = pd.to_numeric(df_cleaned_y[y_column_name], errors='coerce')

            df_cleaned_y = df_cleaned_y.dropna()

            x_column = df_cleaned_y[x_column_name].values / 1000
            y_column = df_cleaned_y[y_column_name].values

            # Fit the data to the logistic function using lmfit
            model = Model(logistic_function)
            params = model.make_params(a=0.1, b=5) 
            result = model.fit(y_column, params, x=x_column)
            
            # Generate x values for plotting
            x_values = np.linspace(-50, 250, 120)

            # Calculate corresponding y values using the optimized parameters and the logistic function
            y_values = logistic_function(x_values, **result.best_values)
            
            # Calculate confidence intervals for a_fit and b_fit
            ci = result.conf_interval()
            conf_int_a = ci['a']
            conf_int_b = ci['b']

            # Extract the last row of corresponding y_column from df_sheet
            if y_column_name in last_row_values.index:
                last_row_value = last_row_values[y_column_name]

                # Divide a_fit by the last row of y_column
                kTh = result.params['a'].value/ (last_row_value)
                kTh_lb = conf_int_a[1][1]/ (last_row_value)
                kTh_ub = conf_int_a[5][1]/ (last_row_value)
                
                # Calculate qe using the formula b_fit * Q / kTh
                qe = result.params['b'].value * Q / (kTh*x_s)
                qe_lb = conf_int_b[1][1] * Q / (kTh*x_s) 
                qe_ub = conf_int_b[5][1] * Q / (kTh*x_s)                 
                                       
                #find thomas model BV estimates
                BV20 = fsolve(equation_bv20, 0, args=(result.params['a'].value, result.params['b'].value))*1000
                BV20_lb = fsolve(equation_bv20, 0, args=(conf_int_a[5][1], conf_int_b[1][1]))*1000
                BV20_ub = fsolve(equation_bv20, 0, args=(conf_int_a[1][1], conf_int_b[5][1]))*1000
                
                # Extract standard errors for a and b
                stderr_a = result.params['a'].stderr
                stderr_b = result.params['b'].stderr

                # Compute t-statistics
                t_stat_a = result.params['a'].value / stderr_a
                t_stat_b = result.params['b'].value / stderr_b

                # Degrees of freedom (df) = number of observations - number of parameters
                df = len(x_column) - 2  # 2 parameters: a and b

                # Compute p-values using the survival function (1 - CDF)
                p_value_a = 2 * (1 - t.cdf(abs(t_stat_a), df))
                p_value_b = 2 * (1 - t.cdf(abs(t_stat_b), df))
                
                # Store regression results and diagnostics
                regression_results[f'{sheet_name}_{y_column_name}'] = {
                    'a_fit': result.params['a'].value,
                    'b_fit': result.params['b'].value,
                    'stderr_a': stderr_a,
                    'stderr_b': stderr_b,
                    't_stat_a': t_stat_a,
                    't_stat_b': t_stat_b,
                    'p_value_a': p_value_a,
                    'p_value_b': p_value_b,
                    'conf_int_a': conf_int_a,
                    'conf_int_b': conf_int_b,
                    'residuals': result.residual,
                    'kTh': kTh,
                    'kTh_lb': kTh_lb,
                    'kTh_ub': kTh_ub,
                    'qe': qe,
                    'qe_lb': qe_lb,
                    'qe_ub': qe_ub,
                    'BV20': BV20[0],
                    'BV20_lb': BV20_lb[0],
                    'BV20_ub': BV20_ub[0],
                    'rsquared': 1 - result.redchi / np.var(y_column),
                    'Co': last_row_value
                    }
                
            # Scatter plot with color based on unique elements
            # Scatter plot with different marker styles
                marker_styles = ['o', 's', '^', 'D', 'v', '<', '>']  # Add more styles as needed
                marker = marker_styles[i % len(marker_styles)]
                scatter = axs.scatter(x_column, y_column, label=f'{y_column_name}', marker=marker)  # Automatic color assignment
                axs.plot(x_values, y_values, color=scatter.get_facecolor()[0], linewidth=1)

        axs.legend(loc='upper right', bbox_to_anchor=(1,1), fontsize=9)  # Add legend for the current figure
        axs.set_xlabel('Bed Volumes (x1000)', fontsize=16)  
        axs.set_ylabel('C/C\u2080', fontsize=16)
        axs.set_ylim(0,1.5)
        
        # Save the figure in the working directory
        plt.savefig(f'./{out_dir}/{sheet_name}_{group_name}.png', bbox_inches='tight')
        plt.close()
        
        #plt.show()

# Display regression results and diagnostics
for key, value in regression_results.items():
    print(f"Results for {key}:")
    print(f"kTh: {value['kTh']} [{value['kTh_lb']}, {value['kTh_ub']}]")
    print(f"qe: {value['qe']} [{value['qe_lb']}, {value['qe_ub']}]")
    #print(f"Co: {value['Co']}")
    #print(f"a: {value['a_fit']}")
    #print(f"b: {value['b_fit']}")
    #print(f"pval_a: {value['p_value_a']}")
    #print(f"pval_b: {value['p_value_b']}")
    print(f"R-squared: {value['rsquared']}")
    #print(f"BV20_lb: {value['BV20_lb']}")
    #print(f"BV20_ub: {value['BV20_ub']}")
    print("\n")

# Create an Excel writer object
with pd.ExcelWriter('fit_data.xlsx', engine='xlsxwriter') as writer:

    # Iterate over the regression_results dictionary
    for sheet_name, df_sheet_cleaned in all_experiments.items():
        # Create an empty DataFrame to store the results for the current sheet
        df_regression_results = pd.DataFrame()

        # Iterate over the regression_results dictionary for the current sheet
        for key, value in regression_results.items():
            # Extract x_values from the key
            current_sheet_name, y_column_name = key.split('_')

            # Check if the current sheet matches the iteration sheet
            if current_sheet_name == sheet_name:
                # Generate x values for plotting
                x_values = np.linspace(-25, 300, 1500)

                # Calculate corresponding y values using the optimized parameters and the logistic function
                y_values = logistic_function(x_values, value['a_fit'], value['b_fit'])

                # Create a DataFrame for the current regression results
                df_temp = pd.DataFrame({'x_values': x_values, y_column_name: y_values})

                # Append the DataFrame to df_regression_results
                if df_regression_results.empty:
                    df_regression_results = df_temp
                else:
                    df_regression_results = pd.merge(df_regression_results, df_temp, on='x_values', how='outer')

        # Write the DataFrame to the Excel file as a separate sheet
        df_regression_results.to_excel(writer, sheet_name=sheet_name, index=False)

# Path for the Excel file
excel_file_path = 'thomas_params.xlsx'

# Create a new Excel workbook and select the active worksheet
wb = Workbook()
ws = wb.active
ws.title = "Results"

# Define field names
fieldnames = ['WQ', 'PFAS', 'kTh', 'kTh_lb', 'kTh_ub', 'qe', 'qe_lb', 'qe_ub', 'BV20', 'BV20_lb', 'BV20_ub', 'rsquared', 'a_pval', 'b_pval', 'Co']

# Write the header row
ws.append(fieldnames)

# Write the data
for key, value in regression_results.items():
    experiment_column, y_column_name = key.split('_')
    ws.append([
        experiment_column,
        y_column_name,
        value['kTh'],
        value['kTh'] - value['kTh_lb'],
        value['kTh_ub'] - value['kTh'],
        value['qe'],
        value['qe'] - value['qe_lb'],
        value['qe_ub'] - value['qe'],
        value['BV20'],
        value['BV20'] - value['BV20_lb'],
        value['BV20_ub'] - value['BV20'],
        value['rsquared'],
        value['p_value_a'],
        value['p_value_b'],
        value['Co']
    ])

# Save the workbook
wb.save(excel_file_path)

# %%
