import math

def determine_domain_parameters(domain_parameters, positions, r_cutoff):
    lowest_values = positions.min(axis=0)
    highest_values = positions.max(axis=0)
    
    domain_parameters['origin'] = lowest_values
    domain_parameters['x_length'] = float(highest_values[1] - lowest_values[1])
    domain_parameters['y_length'] = float(highest_values[2] - lowest_values[2])
    domain_parameters['z_length'] = float(highest_values[3] - lowest_values[3])


def determine_cell_parameters(domain_parameters, cell_parameters, r_cutoff):
    cell_parameters['cells_in_x'] = math.ceil(domain_parameters['x_length'] / r_cutoff)
    cell_parameters['cells_in_y'] = math.ceil(domain_parameters['y_length'] / r_cutoff)    
    cell_parameters['cells_in_z'] = math.ceil(domain_parameters['z_length'] / r_cutoff)

    cell_parameters['number_cells'] = cell_parameters['cells_in_x'] * cell_parameters['cells_in_y'] * cell_parameters['cells_in_z']

def assign_atoms_to_cells(cell_parameters, domain_parameters,  positions):
    number_atoms, row_length = positions.shape
    width = cell_parameters['cells_in_x']
    depth = cell_parameters['cells_in_y']
    origin = domain_parameters['origin']

    max_x_index = cell_parameters['cells_in_x'] - 1
    max_y_index = cell_parameters['cells_in_y'] - 1
    max_z_index = cell_parameters['cells_in_z'] - 1

    for i in range(number_atoms):
        relative_x = positions[i, 1] - origin[1]
        x_index = int((relative_x / domain_parameters['x_length']) * max_x_index)
        
        relative_y = positions[i, 2] - origin[2]
        y_index = int((relative_y / domain_parameters['y_length']) * max_y_index)
        



        relative_z = positions[i, 3] - origin[3]
        z_index = int((relative_z / domain_parameters['z_length']) * max_z_index) 
        


        assigned_cell = (z_index * depth * width) + (y_index * width) + x_index
        positions[i, 4] = assigned_cell


def limit_to_range(value, floor, ceiling):
    if value < floor:
        return floor
    elif value > ceiling:
        return ceiling
    else:
        return value



def check_if_cells_already_paired(cell_i, cell_j, i_iterator, j_iterators):
    pair = {cell_i, cell_j}

    for i in range(len(i_iterator)):
        j_iterator = j_iterators[i]
        if j_iterator is None:
            continue
        cell_i2 = i_iterator[i]

        
        for j in range(len(j_iterator)):
            cell_j2 = j_iterator[j]
            pair2 = {cell_i2, cell_j2}
            if pair == pair2:
                return True
    
    return False


def find_indexes(cell_index, cell_parameters):
    depth = cell_parameters['cells_in_y']
    width = cell_parameters['cells_in_x']

    for x_index in range(cell_parameters['cells_in_x']):
        for y_index in range(cell_parameters['cells_in_y']):
            for z_index in range(cell_parameters['cells_in_z']):
                cell = (z_index * depth * width) + (y_index * width) + x_index

                if cell_index == cell:
                    return x_index, y_index, z_index

def create_cell_position_dictionary_and_list(cell_parameters):
    cell_dictionary = dict()
    cell_positions = list()

    cell = 0
    for z in range(cell_parameters['cells_in_z']):
        for y in range(cell_parameters['cells_in_y']):
            for x in range(cell_parameters['cells_in_x']):
                cell_position = f"{x} {y} {z}"
                cell_dictionary[cell_position] = cell
                cell_positions.append([x, y, z]) 
                cell += 1

    return cell_dictionary, cell_positions

def create_cell_pairs_list(cell_parameters, positions):
    cell_dictionary, cell_positions = create_cell_position_dictionary_and_list(cell_parameters)
    valid_cells = set(positions[:,4])

    cell_pairs = list()
    for i in range(cell_parameters['number_cells']):
        if i not in valid_cells:
            continue
        cell_i_indexes = cell_positions[i]

        for z_offset in [0, 1, -1]:
            for y_offset in [0, 1, -1]:
                for x_offset in [0, 1, -1]:
                    cell_j_indexes = [cell_i_indexes[0] + x_offset, cell_i_indexes[1] + y_offset, cell_i_indexes[2] + z_offset]
                    cell_j_position = f"{cell_j_indexes[0]} {cell_j_indexes[1]} {cell_j_indexes[2]}"
                    
                    try: 
                        j = cell_dictionary[cell_j_position]
                    except KeyError:
                        continue

                    if j not in valid_cells:
                        continue

                    cell_pair = {i, j}
                    if cell_pair in cell_pairs:
                        continue
                    else:
                        cell_pairs.append(cell_pair)

    return cell_pairs