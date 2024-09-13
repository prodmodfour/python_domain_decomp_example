import numpy as np
import os


def read_number_atoms(filename):
    with open(filename) as file:
        lines = file.readlines()
        number_atoms = int(lines[0])
        return number_atoms
    
def initialise_arrays(filename, positions, velocities, forces, number_atoms, number_impact_atoms):

    with open(filename) as file:
        # We skip the first two lines as the first atom starts on line 3
        lines = file.readlines()[2:]

        for line in range(number_atoms):
            converted_line = np.zeros(4, dtype=float)

            split_line = lines[line].split()
            for index in range(4):
                if index == 0:
                    continue
                converted_element = float(split_line[index])
                converted_line[index] = converted_element
            
            positions[line, 0] = line
            velocities[line, 0] = line
            forces[line, 0] = line
            for index in range(1, 4):
                positions[line, index] = converted_line[index]

        for i in range(number_impact_atoms):
            positions[number_atoms + i, 0] = number_atoms + i
            velocities[number_atoms + i, 0] = number_atoms + i
            forces[number_atoms + i, 0] = number_atoms + i

def reset_origin(positions_np_array):
    lowest_values = positions_np_array.min(axis=0)

    if (lowest_values[1] < 0):
        positions_np_array = subtract_from_column(positions_np_array, 1, lowest_values[1])
    if (lowest_values[2] < 0):
        positions_np_array = subtract_from_column(positions_np_array, 2, lowest_values[2])
    if (lowest_values[3] < 0):
        positions_np_array = subtract_from_column(positions_np_array, 3, lowest_values[3])
    
    return positions_np_array

def subtract_from_column(array, column_index, value):
    new_array = array.copy()
    new_array[:,column_index] -= value
    return new_array

def initialise_trajectory_file( positions_np_array, number_atoms, trajectory_filename):
        number_of_atoms = str(number_atoms)
        os.makedirs(os.path.dirname('output/'), exist_ok=True)

        with open(trajectory_filename, "w") as file:
            file.write(number_of_atoms)

        with open(trajectory_filename, 'a') as file:
            file.write('\ntimestep 0') # Title
            write_to_trajectory_file(file, positions_np_array, number_atoms)

def determine_output_filename(filename):
    available_filename_found = False
    filename_index = 1
    while available_filename_found is False:

        if filename_index == 1:
            new_filename = f'output/{filename}.xyz'
        else:
            new_filename = f'output/{filename}_{filename_index}.xyz'
        
        if os.path.exists(new_filename) is True:
            filename_index +=1
        elif os.path.exists(new_filename) is False:
            available_filename_found = True


    return new_filename

def write_to_trajectory_file(file, positions_np_array, number_atoms):
    for index in range(number_atoms):
        atom_key = 'Cu'

        position = positions_np_array[index]
        i, x, y, z, cell_index = position

        # Change sig figs
        x = round(x, 4)
        y = round(y, 4)
        z = round(z, 4)
        line = f"\n{atom_key} {x} {y} {z}"
        file.write(line)

def write_history(positions_np_array, number_atoms,  timestep, filename):
    number_of_atoms = str(number_atoms)


    with open(filename, 'a') as file:
        file.write(f"\n{number_of_atoms}")
        file.write(f'\ntimestep {timestep}')

        write_to_trajectory_file(file, positions_np_array, number_atoms)

def initialise_energies_file(filename, number_timesteps):
        number_timesteps = str(number_timesteps)
        os.makedirs(os.path.dirname('output/'), exist_ok=True)

        with open(filename, "w") as file:
            file.write("Timestep TPE TKE Total_Energy")

def write_energies_to_file(filename, timestep, tpe, tke):
    with open(filename, 'a') as file:
        # Change sig figs
        # tpe = round(tpe, 4)
        # tke = round(tke, 4)

        file.write(f"\n{timestep} {tpe} {tke} {tpe + tke}")

def write_pair_checks_to_file(checked_pairs, calculated_pairs):
    filename = "checked_pairs"
    filename = determine_output_filename(filename)

    with open(filename, "w") as file:
        file.write("Timestep Checked Calculated Efficiency")
        file.write(f"\nSum {sum(checked_pairs)} {sum(calculated_pairs)} {(sum(checked_pairs) / sum(calculated_pairs)) * 100}%")
        for i in range(len(checked_pairs)):
            file.write(f"\n{i} {checked_pairs[i]} {calculated_pairs[i]} ")