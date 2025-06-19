import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from drormd import plot


#This is how you would import it:
    # Correctly expand the path to your Desktop/dror directory

#script_dir = os.path.expanduser("~/Desktop/dror")
#sys.path.append(script_dir)

# #Now try to import

#import useful
#importlib.reload(useful)


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from drormd import plot


def plot_smooth_frequency(data, labels=None, title="Smoothed Frequency Distribution",
                          vert=None, leg=True, xaxis='Value', colors=None, xlim=None):
    """
    Plots smoothed frequency distributions (Kernel Density Estimates) for multiple datasets.

    Parameters:
        data (list of array-like): List of numerical datasets to visualize via KDE.
        labels (list of str, optional): Names for each dataset, used in the legend. Defaults to 'Dataset 1', etc.
        title (str): Title of the plot.
        vert (list of float, optional): List of x-values at which to draw vertical dashed reference lines.
        leg (bool): If True, includes a legend on the plot.
        xaxis (str): Label for the x-axis.
        colors (list of str, optional): Colors for each dataset line. Defaults to a preset palette if not given.
        xlim (list of float, optional): Specifies the x-axis limits as [xmin, xmax].

    Returns:
        None. Displays a matplotlib plot with smoothed frequency distributions.
    """
    if labels is None:
        labels = [f"Dataset {i+1}" for i in range(len(data))]

    if colors is None or len(colors) < len(data):
        default_palette = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666']
        colors = (colors or []) + default_palette[len(colors or []):len(data)]

    # Compute KDEs
    kdes = [gaussian_kde(d) for d in data]

    # Determine x range
    all_min = min(np.min(d) for d in data)
    all_max = max(np.max(d) for d in data)
    dist = all_max - all_min
    x_values = np.linspace(all_min - dist * 0.05, all_max + dist * 0.05, 1000)

    # Plot
    plt.figure(figsize=(10, 5))
    for i, kde in enumerate(kdes):
        y = kde(x_values)
        plt.plot(x_values, y, label=labels[i], linewidth=6, color=colors[i])

    if vert:
        for val in vert:
            plt.axvline(val, linestyle='dashed', color='red')

    plt.xlabel(xaxis, size=25)
    plt.ylabel("Frequency", size=30)
    plt.xticks(size=30)
    plt.yticks([])
    plt.title(title)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)

    if xlim:
        plt.xlim(xlim[0], xlim[1])

    if leg:
        plt.legend()

    plt.show()



def create_conditions_dict(user_dict):
    """
    Generates a dictionary of atom selection strings based on GPCRDB-style keys by prompting the user
    for reference residues and segment IDs.

    Parameters:
        user_dict (dict): Dictionary where keys are GPCRDB labels (e.g., '3x50', 'TM3', 'H3') and values
                          are initially empty strings. Keys may represent individual residues or whole helices.

    Returns:
        None. Prints a formatted dictionary mapping GPCRDB notation to PyMOL-style atom selection strings
        using segment ID and residue number ranges.

    Behavior:
        - Prompts the user to enter a segment ID shared by all selections.
        - For each helix inferred from keys, prompts the user for the residue number that corresponds to .50.
        - Calculates selection strings based on residue offsets and GPCRDB positions.
        - Updates the original dictionary in-place and prints it in a usable format.
    """
    # Identify the helices present
    helices_used = set()

    for key in user_dict.keys():
        if 'x' in key:
            helix = key.split('x')[0]
            if helix.isdigit():
                helices_used.add(int(helix))
        elif key.startswith('TM') and key[2:].isdigit():
            helices_used.add(int(key[2:]))

    # Prompt for segid
    segid = input("Enter the segid for the protein: ").strip()

    # Prompt for reference residue numbers (.50) for each helix
    helix_refs = {}
    for helix in sorted(helices_used):
        ref_res = int(input(f"Enter the residue number for {helix}.50: "))
        helix_refs[helix] = {'ref_res': ref_res, 'segid': segid}

    # Fill in the dictionary
    for key in user_dict.keys():
        if 'x' in key:
            helix, gpcrdb_pos = key.split('x')
            if helix.isdigit():
                helix_num = int(helix)
                gpcrdb_pos = int(gpcrdb_pos)
                offset = gpcrdb_pos - 50
                ref_res = helix_refs[helix_num]['ref_res']
                resid = ref_res + offset
                segid = helix_refs[helix_num]['segid']
                user_dict[key] = f"protein and segid {segid} and resid {resid}"
        elif key.startswith('TM') and key[2:].isdigit():
            helix_num = int(key[2:])
            ref_res = helix_refs[helix_num]['ref_res']
            segid = helix_refs[helix_num]['segid']

            # Prompt for GPCRDB range
            gpcrdb_range = input(f"Enter the GPCRDB range for TM{helix_num} (e.g., 3.44-3.56): ").strip()
            start_gpcrdb, end_gpcrdb = gpcrdb_range.split('-')
            start_gpcrdb = int(start_gpcrdb.split('.')[1])
            end_gpcrdb = int(end_gpcrdb.split('.')[1])

            # Calculate residue range
            start_res = ref_res + (start_gpcrdb - 50)
            end_res = ref_res + (end_gpcrdb - 50)

            user_dict[key] = f"protein and segid {segid} and resid {start_res} to {end_res}"
        elif key.startswith('H') and key[1:].isdigit():
            helix_num = int(key[1:])
            ref_res = helix_refs[helix_num]['ref_res']
            segid = helix_refs[helix_num]['segid']

            # Prompt for GPCRDB range
            gpcrdb_range = input(f"Enter the GPCRDB range for TM{helix_num} (e.g., 3.44-3.56): ").strip()
            start_gpcrdb, end_gpcrdb = gpcrdb_range.split('-')
            start_gpcrdb = int(start_gpcrdb.split('.')[1])
            end_gpcrdb = int(end_gpcrdb.split('.')[1])

            # Calculate residue range
            start_res = ref_res + (start_gpcrdb - 50)
            end_res = ref_res + (end_gpcrdb - 50)

            user_dict[key] = f"protein and segid {segid} and resid {start_res} to {end_res}"
        else:
            # Keep any other pre-defined entries as-is, or update manually if desired
            if user_dict[key] == '':
                user_dict[key] = f"# Please update manually for {key}"

    # --- OUTPUT ---
    print("\nGenerated dictionary:\n")
    print("generic = {")
    for k, v in user_dict.items():
        print(f"    '{k}': '{v}',")
    print("}")






def add_min_column(csv_path, column1, column2, new_column_name, output_path=None):
    """
    Given a CSV file and two column names, this function creates a new column containing
    the minimum value of the two specified columns for each row and saves the modified CSV.

    Parameters:
        csv_path (str): Path to the input CSV file.
        column1 (str): Name of the first column.
        column2 (str): Name of the second column.
        new_column_name (str): Name of the new column that will store the minimum values.
        output_path (str, optional): Path to save the modified CSV. If None, overwrites the input file.

    Returns:
        pd.DataFrame: The modified DataFrame with the new column.
    """
    # Load CSV
    df = pd.read_csv(csv_path)

    # Ensure columns exist
    if column1 not in df.columns or column2 not in df.columns:
        raise ValueError(f"Columns '{column1}' and/or '{column2}' not found in CSV.")

    # Compute minimum value row-wise
    df[new_column_name] = df[[column1, column2]].min(axis=1)

    # Save the modified CSV
    output_path = output_path or csv_path  # If no output path is given, overwrite the input file
    df.to_csv(output_path, index=False)

    return df


def make_lists(data, cond, reps, start=0, finish=1000, offset=0, skip=False, met=None, filt=False):
    """
    Extracts a combined list of metric values from multiple simulation replicates for a given condition.

    Parameters:
        data (dict): Nested dictionary containing simulation data organized as data[condition][replicate] = DataFrame.
        cond (str): Condition name (e.g., 'WT', 'mutantA') used to access data[cond].
        reps (int): Number of replicates for the condition.
        start (int): Start index for slicing each metric trace. Default is 0.
        finish (int): End index for slicing each metric trace. Default is 1000.
        offset (int): Offset added to start/finish indices. Useful for burn-in trimming.
        skip (bool): Whether to skip replicate 10 (useful for known problematic runs). Default is False.
        met (str): The metric to extract (e.g., 'chi1', 'distance', etc.). If None, nothing is extracted.
        filt (tuple or bool): Optional filtering condition as a tuple: (column, operator, threshold),
                              e.g., ('chi1', '<', 90). If False, no filtering is applied.

    Returns:
        list: Flattened list of values across all replicates for the specified metric, optionally filtered and windowed.
              Dihedral angles are wrapped to the [0, 360) or [-180, 180] range depending on content.
    """
    if met:
        met_dist = []
        for i in range(1, reps + 1):
            if i == 10 and skip == True:
                continue
            r = f'{cond}_r{i}'
            data_r = data[cond][r]
            if filt:
                if filt[1] == '<':
                    data_r_2 = data_r[data_r[filt[0]] < filt[2]].reset_index(drop=True)
                elif filt[1] == '>':
                    data_r_2 = data_r[data_r[filt[0]] > filt[2]].reset_index(drop=True)
            if 'chi' in met:
                if not filt:
                    for i in range(start + offset, finish + offset):
                        if i < len(data_r[met]):
                            val = data_r[met][i]
                            if val < 0:
                                val += 360
                            if val > 360:
                                val -= 360
                            met_dist.append(val)
                else:
                    for i in range(len(data_r_2[met])):
                        if i < len(data_r_2[met]):
                            val = data_r_2[met][i]
                            if val < 0:
                                val += 360
                            if val > 360:
                                val -= 360
                            met_dist.append(val)
            else:
                if not filt:
                    met_dist += list(data_r[met].values)[start + offset:finish + offset]
                else:
                    met_dist += list(data_r_2[met].values)
        return met_dist




def sliding_mean(data_array, window):
    new_list = []
    for i in range(data_array.shape[0]):
        if i <= window:
            indices = range(0,2*i+1)
        else:
            indices = range(max(i - window, 0),
                            min(i + window + 1, len(data_array)))
        avg = 0
        for j in indices:
            avg += data_array[j]
        avg /= float(len(indices))
        new_list.append(avg)
    return np.array(new_list)


colors_lst = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']


def plot_time_traces(conditions, data, colors=colors_lst, leg_title=''):
    """
    Plots smoothed time traces and kernel density estimates (KDEs) for a given metric across multiple conditions.

    Parameters:
        conditions (list of str): List of condition names (e.g., ['WT', 'mutantA']) corresponding to keys in `data`.
        data (dict): Nested dictionary organized as data[condition][replicate][metric] = time series array.
        colors (list of str): List of colors to use for each condition. Defaults to a predefined color palette.
        leg_title (str): Title for the legend box.

    Returns:
        None. Generates matplotlib plots for each metric across the specified conditions.

    Behavior:
        - For each metric in the first replicate of the first condition, the function:
            - Gathers all corresponding metric traces across replicates and conditions.
            - Applies special angle wrapping for metrics containing 'chi' to ensure cyclic continuity.
            - Plots time traces and KDE distributions using helper functions from `drormd.plot`.
        - Distinguishes metrics as either dihedral angles or distance-like values to label axes accordingly.
    """
    for metric in data[conditions[0]][f'{conditions[0]}_r1']:
        trace, kde = plot.setup_trace_and_kde(metric, '')
        for i, condition in enumerate(conditions):
            A = [data[condition][rep][metric]
                 for rep in data[condition]
                 if metric in data[condition][rep]]
            if not A:
                continue
            if 'chi' in metric:
                plot.dihedral_range(A, 0)
                ylabel = 'degree'
            else:
                ylabel = 'distance'
            plot.add_time_trace(trace, A, f'{condition}', colors=colors[i], lw=0, lw_smooth=1)
            plot.add_kde(kde, A, f'{condition}', flip=True, color=colors[i], burnin=0)
        trace.set_xlim(0)
        trace.set_ylabel(ylabel)
        legend = trace.legend(title=leg_title, loc='best')
        plt.show()






def get_values(data, condition, nreps, metric):
    reps = ['r'+str(i) for i in range(1,nreps+1)] 
    traces = []
    for rep in reps:
        trace = data[condition][condition+'_'+rep][metric].values
        if 'chi' in metric:
            print(type(trace))
            trace_new = []
            for val in trace:
                if val < 0:
                    val += 360
                if val > 230:
                    val -= 360
                trace_new.append(val)
            traces.append(np.array(trace_new))
        else:
            traces.append(trace)
    return traces  

