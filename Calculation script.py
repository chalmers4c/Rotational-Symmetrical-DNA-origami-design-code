# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 23:40:08 2025
Combined all the calculation into the same script
This works with even number of helix, unsure how odd number of helix in the subunit will perform
@author: chalmers
"""

import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
from tkinter import messagebox
from tkinter import ttk
import math
import pandas as pd
#%%
# Global variable adjustment function
def Global_variable_selector():
    input_window = tk.Tk()
    input_window.title("Rotational symmetry origami initial parameter calculator")
    input_window.attributes('-topmost', True)

    # Grid layout configuration (without left_frame)
    input_window.columnconfigure(0, weight=1)  # Label column
    input_window.columnconfigure(1, weight=2)  # Entry column

    # Labels and entry widgets for 6 variables
    tk.Label(input_window, text="Number of subunits: ", font = 16).grid(row=0, column=0, padx=10, pady=5, sticky="e")
    variable1_entry = tk.Entry(input_window, width = 10)
    variable1_entry.insert(0, 5)
    variable1_entry.grid(row=0, column=1, padx=10, pady=5)

    tk.Label(input_window, text="Rise per basepair (nm/bp): ", font = 16).grid(row=1, column=0, padx=10, pady=5, sticky="e")
    variable2_entry = tk.Entry(input_window, width = 10)
    variable2_entry.insert(0, 0.34)
    variable2_entry.grid(row=1, column=1, padx=10, pady=5)

    tk.Label(input_window, text="Interhelical distance between helices (nm): ", font = 16).grid(row=2, column=0, padx=10, pady=5, sticky="e")
    variable3_entry = tk.Entry(input_window, width = 10)
    variable3_entry.insert(0, 2.69)
    variable3_entry.grid(row=2, column=1, padx=10, pady=5)

    tk.Label(input_window, text="Minimum helix length of the subunit (multiples of 8) (bp):  ", font = 16).grid(row=3, column=0, padx=10, pady=5, sticky="e")
    variable4_entry = tk.Entry(input_window, width = 10)
    variable4_entry.insert(0, 48)
    variable4_entry.grid(row=3, column=1, padx=10, pady=5)

    tk.Label(input_window, text="Number of helix bundle in a subunit: ", font = 16).grid(row=4, column=0, padx=10, pady=5, sticky="e")
    variable5_entry = tk.Entry(input_window, width = 10)
    variable5_entry.insert(0, 16)
    variable5_entry.grid(row=4, column=1, padx=10, pady=5)

    # Function to retrieve values from entries and print/store them
    def store_values():
        
        variable1 = variable1_entry.get()
        variable2 = variable2_entry.get()
        variable3 = variable3_entry.get()
        variable4 = variable4_entry.get()
        variable5 = variable5_entry.get()
        
        global parameterlist
        parameterlist = [variable1, variable2, variable3, variable4, variable5]
        
        print(f"Number of subunits: {variable1}" + " subunits")
        print(f"Rise per basepair: {variable2}" + " nm/bp")
        print(f"Interhelical distance between helices: {variable3}" + " nm")
        print(f"Minimum helix length of the subunit: {variable4}" + " bp")
        print(f"Number of helix bundle in a subunit: {variable5}" + " helices")
        input_window.destroy()
        return 

    # Button to store the values from entries (placed at the lower-right corner)
    store_button = tk.Button(input_window, text="Proceed to calculation", command=store_values)
    store_button.grid(row=6, column=1, padx=10, pady=20, sticky="e")

    input_window.mainloop()
    return parameterlist

# Create the input window to select options
Global_variable_selector()
Num_subunits = float(parameterlist[0])
Rise_nm_bp = float(parameterlist[1])
Interhelical_distance = float(parameterlist[2])
Min_length_subunit = float(parameterlist[3])
Num_helix = int(parameterlist[4])
#%% Parameter_set_1 calculation
def Parameter_set_1_calculation():    
    global Subunit_angles, tan_angle, Additional_bases, length_shortest_helix
    Subunit_angles = 360/Num_subunits
    Subunit_angles_radians = Subunit_angles * (math.pi / 180)
    tan_angle = math.tan(Subunit_angles_radians/2)
    Additional_bases = Interhelical_distance/(Rise_nm_bp * tan_angle)
    length_shortest_helix = Rise_nm_bp * Min_length_subunit
    return
Parameter_set_1_calculation()
#%% Calculate the number of basepairs per helix
def Calculate_number_basepairs_per_helix():
    # Calculate the number of rows for 'helix number'
    if Num_helix % 2 == 0:
        max_helix_number = Num_helix // 2  # If Num_helix is even
    else:
        max_helix_number = (Num_helix // 2) + 1  # If Num_helix is odd
    # Generate the table with the columns helix number, helix length (bp), and round up (bp)
    data = []
    for helix_num in range(1, max_helix_number + 1):
        # Calculate the helix length (bp) using the formula
        helix_length_bp = Min_length_subunit + (helix_num - 1) * Additional_bases
        # Round up the helix length
        round_up_bp = math.ceil(helix_length_bp)
        # Append the data to the table
        data.append([helix_num, helix_length_bp, round_up_bp])
    df = pd.DataFrame(data, columns=["Helix Number", "Helix Length (bp)", "Round up (bp)"])
    return df
df_half = Calculate_number_basepairs_per_helix()
def Mirror_dataframe_to_get_full(df = df_half):
    # Check if Num_helix is even or odd
    if Num_helix % 2 == 0:
        # If Num_helix is even, mirror the entire dataframe except the first row
        mirrored_data = df[::-1].reset_index(drop=True)
    else:
        # If Num_helix is odd, mirror all rows except the last one
        mirrored_data = df[::-1].reset_index(drop=True)
        
    df_full = pd.concat([df, mirrored_data], ignore_index=True)
    df_full["Helix Number"] = range(1, len(df_full)+1)
    print(df_full)
    return df_full
df_full = Mirror_dataframe_to_get_full()
#%% Parameter_set_2_calculation()
def Parameter_set_2_calculation():
    global Max_length_subunit, Height_of_subunit, Total_staple_per_subunit, Scaffold_length_per_subunit, Scaffold_not_used
    Max_length_subunit = max(df_full.iloc[:, 1])
    Height_of_subunit = ((Max_length_subunit-1)*Rise_nm_bp)+(Interhelical_distance/tan_angle)
    Total_staple_per_subunit = sum(df_full.iloc[:, 2])
    Scaffold_length_per_subunit = round(7249 / Num_subunits)
    Scaffold_not_used = Scaffold_length_per_subunit-Total_staple_per_subunit
    return
Parameter_set_2_calculation()
#%% Function to calculate scaffold coordinates
def calculate_scaffold_coordinates(i, j, Interhelical_distance = Interhelical_distance, Num_helix = Num_helix, Height_of_subunit = Height_of_subunit):
    """
    Calculate the coordinates for the scaffold nucleotide in the jth helix.
    Parameters:
    i (int): Index of the scaffold nucleotide.
    j (int): Index of the helix.
    Returns:
    tuple: (x, y, z) coordinates of the scaffold nucleotide.
    """
    x = -Height_of_subunit + 0.34 * (i - 1)
    if j % 2 == 1:  # j is odd
        y = (Interhelical_distance - 2) / 2 + (Num_helix / 2 - j) * Interhelical_distance + (1 - math.cos(math.radians(360 * (i - 1) / 10.44 + 15)))
        z = -math.sin(math.radians(360 * (i - 1) / 10.44 + 15))
    else:  # j is even
        y = -(Interhelical_distance - 2) / 2 + (Num_helix / 2 - j + 1) * Interhelical_distance - (1 - math.cos(math.radians(360 * (i - 1) / 10.44 - 15)))
        z = math.sin(math.radians(360 * (i - 1) / 10.44 - 15))
    
    return (round(x, 2), round(y, 2), round(z, 2))

#%% Function to calculate staple coordinates
def calculate_staple_coordinates(i, j, Interhelical_distance = Interhelical_distance, Num_helix = Num_helix, Height_of_subunit = Height_of_subunit):
    """
    Calculate the coordinates for the staple nucleotide in the jth helix.
    Parameters:
    i (int): Index of the staple nucleotide.
    j (int): Index of the helix.
    Returns:
    tuple: (x, y, z) coordinates of the staple nucleotide.
    """
    x = -Height_of_subunit + 0.34 * (i - 1)
    if j % 2 == 1:  # j is odd
        y = (Interhelical_distance - 2) / 2 + (Num_helix / 2 - j) * Interhelical_distance + (1 - math.cos(math.radians(360 * (i - 1) / 10.44 + 150 + 15)))
        z = -math.sin(math.radians(360 * (i - 1) / 10.44 + 150 + 15))
    else:  # j is even
        y = -(Interhelical_distance - 2) / 2 + (Num_helix / 2 - j + 1) * Interhelical_distance - (1 - math.cos(math.radians(360 * (i - 1) / 10.44 - 150 - 15)))
        z = math.sin(math.radians(360 * (i - 1) / 10.44 - 150 - 15))
    
    return (round(x, 2), round(y, 2), round(z, 2))
#%% list calculation
def Coordination_list_calculation():
    indexes = []
    scaffold_nucleotide_coords_list = []
    staple_nucleotide_coords_list = []

    for j, i in df_full.iloc[:, [0, 2]].itertuples(index=False, name=None):  
        scaffold_nucleotide_coords = calculate_scaffold_coordinates(i, j)  # i = staple nucleotide length, j = helix bundle number  
        staple_nucleotide_coords = calculate_staple_coordinates(i, j)

        indexes.append(f"i = {i}, j = {j}")
        scaffold_nucleotide_coords_list.append(scaffold_nucleotide_coords)
        staple_nucleotide_coords_list.append(staple_nucleotide_coords)
    df = pd.DataFrame({
        'Index': indexes,
        'Scaffold Nucleotide Coordinates': scaffold_nucleotide_coords_list,
        'Staple Nucleotide Coordinates': staple_nucleotide_coords_list
    })
    print(df)
    return df
Coordination_list = Coordination_list_calculation()
#%%
def Coordination_list_transforamtion():
    # Define the rotation angle and convert to radians
    theta_deg = 360/Num_subunits * (Num_subunits - 1) # here because pentagon is formed by 5 triangle, the acute angle thus will be 360/5 = 72 degrees. Rotate that triangle by 4×72 will align the triangle to its right.
    theta_rad = math.radians(theta_deg)

    # Define the rotation matrix
    rotation_matrix = [
        [math.cos(theta_rad), -math.sin(theta_rad)],
        [math.sin(theta_rad), math.cos(theta_rad)]
    ]

    def rotate_coordinates(x, y, rotation_matrix):
        """
        Rotate the (x, y) coordinates using the provided rotation matrix.
        """
        x_rot = rotation_matrix[0][0] * x + rotation_matrix[0][1] * y
        y_rot = rotation_matrix[1][0] * x + rotation_matrix[1][1] * y
        return (round(x_rot, 2), round(y_rot, 2))

    # Initialize lists for transformed coordinates
    scaffold_nucleotide_transformed = []
    staple_nucleotide_transformed = []

    # Apply transformation to each row in the DataFrame
    for _, row in Coordination_list.iterrows():
        x_nuc, y_nuc, z_nuc = row['Scaffold Nucleotide Coordinates']
        x_staple, y_staple, z_staple = row['Staple Nucleotide Coordinates']
        
        nuc_transformed = rotate_coordinates(x_nuc, y_nuc, rotation_matrix)
        staple_transformed = rotate_coordinates(x_staple, y_staple, rotation_matrix)
        
        scaffold_nucleotide_transformed.append((nuc_transformed[0], nuc_transformed[1], z_nuc))
        staple_nucleotide_transformed.append((staple_transformed[0], staple_transformed[1], z_staple))

    # Create df_transform DataFrame
    df_transform = pd.DataFrame({
        'Index': Coordination_list['Index'],
        'Transformed Scaffold Nucleotide Coordinates': scaffold_nucleotide_transformed,
        'Transformed Staple Nucleotide Coordinates': staple_nucleotide_transformed
    })

    df_transform_reindexed = df_transform.iloc[::-1].reset_index(drop=True)
    print(df_transform_reindexed)
    return df_transform_reindexed
Coordination_list_transformed = Coordination_list_transforamtion()
#%%
def calculate_linker_length(x1, y1, z1, x2, y2, z2):
    """
    Calculate the linker length between two points in 3D space. (euclidean distance)
    
    Parameters:
    x1, y1, z1 (float): Coordinates of the first point.
    x2, y2, z2 (float): Coordinates of the second point.
    
    Returns:
    float: The calculated linker length.
    """
    distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    linker_length = (distance / 0.4) - 1
    return round(linker_length, 2)

#%%
# Function to calculate scaffold loop lengths between consecutive rows

"""Have not tested with an odd number of helix in the subunit, dataframe dropping maynot work"""
df_scaffold_loop_linker_calculation = pd.DataFrame()
def calculate_scaffold_loop_lengths(coords_list):
    lengths = []
    for idx in range(len(coords_list) - 1):
        x1, y1, z1 = coords_list[idx]
        x2, y2, z2 = coords_list[idx + 1]
        length = calculate_linker_length(x1, y1, z1, x2, y2, z2)
        lengths.append(round(length))
    lengths.append(None)  # Append None for the last row which has no next row
    return lengths

# Calculate scaffold loop lengths
df_scaffold_loop_linker_calculation['Scaffold Loop Length'] = calculate_scaffold_loop_lengths(Coordination_list.iloc[:, 1])
# Drop even-numbered indices (including 0)
df_scaffold_loop_linker_calculation = df_scaffold_loop_linker_calculation.iloc[1::2].reset_index(drop=True)
# Create "Helix" column
df_scaffold_loop_linker_calculation.insert(
    0, 
    "Helix", 
    [f"{i} to {i+1}" for i in range(1, len(df_scaffold_loop_linker_calculation) * 2, 2)]
)
df_scaffold_loop_linker_calculation.insert(0, "Helix", df_scaffold_loop_linker_calculation.pop("Helix"))
df_scaffold_loop_linker_calculation = df_scaffold_loop_linker_calculation.dropna()
Total_scaffold_loop_length_subunit = sum(df_scaffold_loop_linker_calculation.iloc[:, 1])
print(df_scaffold_loop_linker_calculation)

#%%

"""Have not tested with an odd number of helix in the subunit, dataframe dropping maynot work"""
def calculate_staple_bridge_length(df, df_transform):
    # Ensure DataFrames are of the same size
    assert df.shape[0] == df_transform.shape[0], "DataFrames must have the same number of rows"

    lengths = []

    for idx in range(df.shape[0]):
        # Extract coordinates from the original DataFrame
        x1, y1, z1 = df.iloc[idx, 2]
        # Extract coordinates from the transformed DataFrame
        x2, y2, z2 = df_transform.iloc[idx, 2]
        
        # Calculate linker length
        length = calculate_linker_length(x1, y1, z1, x2, y2, z2)
        lengths.append(round(length))
    
    # Create DataFrame for linker lengths
    df_staple_bridge_calculation = pd.DataFrame({
        'Staple Bridge Length': lengths
    })
    
    return df_staple_bridge_calculation

# Calculate linker lengths
df_staple_bridge_calculation = calculate_staple_bridge_length(df = Coordination_list, df_transform = Coordination_list_transformed)
# Calculate the max index and half of it
max_index = df_staple_bridge_calculation.index.max()
half_max_index = max_index // 2
# Shrink DataFrame by keeping rows with indices <= half of the max index
df_staple_bridge_calculation = df_staple_bridge_calculation[df_staple_bridge_calculation.index <= half_max_index].reset_index(drop=True)
# Drop odd-numbered indices
#df_staple_bridge_calculation = df_staple_bridge_calculation.iloc[::2].reset_index(drop=True)
# Create "Helix" column
df_staple_bridge_calculation["Helix"] = [f"{i + 1} to {Num_helix - i}" for i in range(len(df_staple_bridge_calculation))]
df_staple_bridge_calculation.insert(0, "Helix", df_staple_bridge_calculation.pop("Helix"))
df_staple_bridge_calculation = df_staple_bridge_calculation.drop(index=0).reset_index(drop=True)

print(df_staple_bridge_calculation)
#%% Summary
def create_summary():
    # Collect results into a list of tuples
    result = [
        ("01. Number of Subunits", Num_subunits),
        ("02. Rise per Base Pair (nm)", Rise_nm_bp),
        ("03. Interhelical Distance (nm)", Interhelical_distance),
        ("04. Minimum Helix Length (nm)", Min_length_subunit),
        ("05. Number of Helices", Num_helix),
        ("06. Subunit Angle (degrees)", Subunit_angles),
        ("07. Tangent angle ratio", tan_angle),
        ("08. Additional base pairs (bp)", Additional_bases),
        ("09. Length of the shortest helix (nm)", length_shortest_helix),
        ("10. Maximum Helix Length (nm)", Max_length_subunit),
        ("11. Height of the subunit (nm)", Height_of_subunit),
        ("12. Total number of staples per subunit (nt)", Total_staple_per_subunit),
        ("13. Scaffold length per subunit (nt)", Scaffold_length_per_subunit),
        ("14. Excess scaffold length", Scaffold_not_used),
        ("15. Scaffold leftover (after calculation of scaffold loop)", Scaffold_not_used-Total_scaffold_loop_length_subunit)
        ]
    global ALL_CALCULATION_SUMMARY
    ALL_CALCULATION_SUMMARY = dict(result)
    print(ALL_CALCULATION_SUMMARY)
    return
create_summary()
#%%
ALL_FINAL_OUTCOME_FOR_CADNANO = {
    "01. Staples length in bp": df_full,
    "02. Scaffold loop length in nt": df_scaffold_loop_linker_calculation,
    "03. Bridge staples length in nt": df_staple_bridge_calculation
}
