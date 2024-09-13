import numpy as np
import math

EV_TO_J_PER_MOLE =  96400.0
J_PER_MOLE_TO_EV =  1.037e-5
cu_mass = 63.546 # Atomic weight in g/mol

def zero_forces(forces, number_atoms):
    forces[:,1] = np.zeros(number_atoms, dtype= float)
    forces[:,2] = np.zeros(number_atoms, dtype= float)
    forces[:,3] = np.zeros(number_atoms, dtype= float)

def evaluate_forces(epsilon, sigma, number_atoms, positions, forces, r_cutoff):
    epsilon4 = 4 * epsilon
    epsilon24 = 24 * epsilon # energy factor
    sigma_squared = sigma**2
    total_potential = 0
    r_cutoff_squared = r_cutoff * r_cutoff

    checked_pairs = 0
    calculated_pairs = 0

    # Loop over all unique pairs of atom i - j
    for i in range(number_atoms):
        xi = positions[i, 1]
        yi = positions[i, 2]
        zi = positions[i, 3]

        for j in range(number_atoms):
            if i >= j:
                continue
            
            xj = positions[j, 1]
            yj = positions[j, 2]
            zj = positions[j, 3]

            checked_pairs += 1
            dx = xi - xj
            if dx > r_cutoff:
                continue
            dy = yi - yj
            if dy > r_cutoff:
                continue
            dz = zi - zj
            if dz > r_cutoff:
                continue
            
            squared_distance = dx*dx + dy*dy + dz*dz
            if squared_distance > r_cutoff_squared:
                continue
            distance = math.sqrt(squared_distance)
            calculated_pairs += 1

            # Calculate force and potential between i and j
            sr2 = sigma_squared / squared_distance
            sr6 = sr2 * sr2 * sr2
            sr12 = sr6 * sr6

            reciprocal_distance = 1 / distance
            unit_vector_x = dx * reciprocal_distance
            unit_vector_y = dy * reciprocal_distance
            unit_vector_z = dz * reciprocal_distance

            potential_ij = sr12 - sr6
            total_potential += potential_ij

            fij = ((2 * sr12) - sr6) * reciprocal_distance
            fxij = fij * unit_vector_x
            fyij = fij * unit_vector_y
            fzij = fij * unit_vector_z

            # Add the force on i due to j
            forces[i, 1] += fxij
            forces[i, 2] += fyij
            forces[i, 3] += fzij

            # Add the force on j to due i in ij direction
            forces[j, 1] -= fxij
            forces[j, 2] -= fyij
            forces[j, 3] -= fzij


    # Multiply all forces by energy factor 24epsilon

    forces[:,1] = forces[:,1] * epsilon24
    forces[:,2] = forces[:,2] * epsilon24
    forces[:,3] = forces[:,3] * epsilon24

    # Multiply by energy factor 4epsilon
    total_potential *= epsilon4
    total_potential *= J_PER_MOLE_TO_EV
    return total_potential, checked_pairs, calculated_pairs


def calculate_kinetic_energy(sum_v_squared):
    kinetic_energy = 0.5 * cu_mass * sum_v_squared
    kinetic_energy *= J_PER_MOLE_TO_EV
    return kinetic_energy

def calculate_pair_interactions_between_cells(i, j, positions, forces, epsilon, sigma, r_cutoff):
    epsilon4 = 4 * epsilon
    epsilon24 = 24 * epsilon # energy factor
    sigma_squared = sigma**2
    r_cutoff_squared = r_cutoff * r_cutoff
    cell_i = positions[np.in1d(positions[:,4], i)]
    cell_j = positions[np.in1d(positions[:,4], j)]

    number_atoms_i, columns = cell_i.shape
    number_atoms_j, columns = cell_j.shape

    potential_energy = 0
    checked_pairs = 0
    calculated_pairs = 0

    for n in range(number_atoms_i):
        a = int(cell_i[n, 0])
        xi = cell_i[n, 1]
        yi = cell_i[n, 2]
        zi = cell_i[n, 3]

        for m in range(number_atoms_j):
            if i == j:
                if n >= m:
                    continue
            
            b = int(cell_j[m, 0])
            xj = cell_j[m, 1]
            yj = cell_j[m, 2]
            zj = cell_j[m, 3]

            checked_pairs += 1
            dx = xi - xj
            if dx > r_cutoff:
                continue
            dy = yi - yj
            if dy > r_cutoff:
                continue
            dz = zi - zj
            if dz > r_cutoff:
                continue
            
            squared_distance = dx*dx + dy*dy + dz*dz
            if squared_distance > r_cutoff_squared:
                continue
            distance = math.sqrt(squared_distance)
            calculated_pairs += 1

            # Calculate force and potential between i and j
            sr2 = sigma_squared / squared_distance
            sr6 = sr2 * sr2 * sr2
            sr12 = sr6 * sr6

            reciprocal_distance = 1 / distance
            unit_vector_x = dx * reciprocal_distance
            unit_vector_y = dy * reciprocal_distance
            unit_vector_z = dz * reciprocal_distance

            potential_ij = sr12 - sr6
            potential_energy += potential_ij

            fij = ((2 * sr12) - sr6) * reciprocal_distance * epsilon24
            fxij = fij * unit_vector_x
            fyij = fij * unit_vector_y
            fzij = fij * unit_vector_z

            # Add the force on i due to j
            forces[a, 1] += fxij
            forces[a, 2] += fyij
            forces[a, 3] += fzij

            # Add the force on j to due i in ij direction
            forces[b, 1] -= fxij
            forces[b, 2] -= fyij
            forces[b, 3] -= fzij

            
        # Multiply by energy factor 4epsilon
        potential_energy *= epsilon4
        potential_energy *= J_PER_MOLE_TO_EV

        return potential_energy, checked_pairs , calculated_pairs


def evaluate_cell_pair(pair, positions, forces, epsilon, sigma, r_cutoff):
    epsilon4 = 4 * epsilon
    epsilon24 = 24 * epsilon # energy factor
    sigma_squared = sigma**2
    r_cutoff_squared = r_cutoff * r_cutoff

    if len(pair) == 1:
        i = list(pair)[0]
        j = list(pair)[0]
    else:
        i, j = pair

    cell_i = positions[np.in1d(positions[:,4], i)]
    cell_j = positions[np.in1d(positions[:,4], j)]

    number_atoms_i, columns = cell_i.shape
    number_atoms_j, columns = cell_j.shape
    
    potential_energy = 0
    checked_pairs = 0
    calculated_pairs = 0
    

    for n in range(number_atoms_i):
        global_index_i = int(cell_i[n, 0])
        xi = cell_i[n, 1]
        yi = cell_i[n, 2]
        zi = cell_i[n, 3]

        for m in range(number_atoms_j):
            global_index_j = int(cell_j[m, 0])
            if i == j:
                if global_index_i >= global_index_j:
                    continue
        
            xj = cell_j[m, 1]
            yj = cell_j[m, 2]
            zj = cell_j[m, 3]

            checked_pairs += 1
            dx = xi - xj
            if dx > r_cutoff:
                continue
            dy = yi - yj
            if dy > r_cutoff:
                continue
            dz = zi - zj
            if dz > r_cutoff:
                continue
            
            squared_distance = dx*dx + dy*dy + dz*dz
            if squared_distance > r_cutoff_squared:
                continue
            distance = math.sqrt(squared_distance)
            calculated_pairs += 1

            # Calculate force and potential between i and j
            sr2 = sigma_squared / squared_distance
            sr6 = sr2 * sr2 * sr2
            sr12 = sr6 * sr6

            reciprocal_distance = 1 / distance
            unit_vector_x = dx * reciprocal_distance
            unit_vector_y = dy * reciprocal_distance
            unit_vector_z = dz * reciprocal_distance

            potential_ij = sr12 - sr6
            potential_energy += potential_ij

            fij = ((2 * sr12) - sr6) * reciprocal_distance
            fxij = fij * unit_vector_x
            fyij = fij * unit_vector_y
            fzij = fij * unit_vector_z

            # Add the force on i due to j
            forces[global_index_i, 1] += fxij
            forces[global_index_i, 2] += fyij
            forces[global_index_i, 3] += fzij

            # Add the force on j to due i in ij direction
            forces[global_index_j, 1] -= fxij
            forces[global_index_j, 2] -= fyij
            forces[global_index_j, 3] -= fzij

 
    # Multiply by energy factor 4epsilon
    potential_energy *= epsilon4
    potential_energy *= J_PER_MOLE_TO_EV

   
    return potential_energy, checked_pairs , calculated_pairs
