import file_functions#
import numpy as np
import math
import time
import physics
import domain_decomposition

EV_TO_J_PER_MOLE =  96400.0
J_PER_MOLE_TO_EV =  1.037e-5
cu_mass = 63.546 # Atomic weight in g/mol

def dd_method(input_filename):
    start = time.time()
    scaling_factor = 0.01

    timestep_size = 0.001 # Picoseconds
    epsilon = 0.4802 * EV_TO_J_PER_MOLE # j/mol
    sigma = 2.285 # Angstrom
    r_cutoff = 2.5 * sigma # Angstrom

    # We read the atom positions from a file and then add in impact atom(s) 
    number_atoms = file_functions.read_number_atoms(input_filename)

    number_impact_atoms = 1
    # [Index, x, y, z, cell_index]
    positions = np.zeros((number_atoms + number_impact_atoms,  5), dtype=float)
    velocities = np.zeros((number_atoms + number_impact_atoms, 5), dtype=float)
    forces = np.zeros((number_atoms + number_impact_atoms,  5), dtype=float)
    file_functions.initialise_arrays(input_filename, positions, velocities, forces, number_atoms, number_impact_atoms)
    number_atoms += number_impact_atoms
    

    # Find centre of top surface
    lowest_values = positions.min(axis=0)
    highest_values = positions.max(axis=0)
    cx = (highest_values[1] + lowest_values[1]) / 2
    cy = (highest_values[2] + lowest_values[2]) / 2
    cz = highest_values[3]


    # Add in impact atom
    applied_energy = 1000 #eV
    applied_energy *= EV_TO_J_PER_MOLE
    impact_atom_speed = -math.sqrt((2 * applied_energy) / cu_mass)
    impact_atom = [number_atoms - 1, cx, cy, cz + 1, 0]

    positions[number_atoms - 1] = impact_atom
    velocities[number_atoms - 1, 0] = number_atoms - 1
    velocities[number_atoms - 1, 3] = impact_atom_speed
    forces[number_atoms - 1, 0] = number_atoms - 1

    number_timesteps = 1000
    trajectory_filename = file_functions.determine_output_filename("trajectories")
    file_functions.initialise_trajectory_file(positions, number_atoms, trajectory_filename)
    energies_filename = file_functions.determine_output_filename("energies")
    file_functions.initialise_energies_file(energies_filename, number_timesteps)

    history_interval = 20
    velocity_scale = (scaling_factor * timestep_size) / cu_mass

    domain_parameters = {
        'origin': [0, 0, 0],
        'x_length': 0,
        'y_length': 0,
        'z_length': 0,
    }

    cell_parameters = {
        'number_cells': 0,
        'cells_in_x': 0,
        'cells_in_y': 0,
        'cells_in_z': 0,
    }

    all_checked_pairs = list()
    all_calculated_pairs = list()
    

    for timestep in range(number_timesteps):
        physics.zero_forces(forces, number_atoms)

        # DOMAIN DECOMPOSITION
        domain_decomposition.determine_domain_parameters(domain_parameters, positions, r_cutoff)
        domain_decomposition.determine_cell_parameters(domain_parameters, cell_parameters, r_cutoff)
        domain_decomposition.assign_atoms_to_cells(cell_parameters, domain_parameters, positions)

        cell_pairs = domain_decomposition.create_cell_pairs_list(cell_parameters, positions)

        total_potential_energy = 0
        for pair in cell_pairs:
            potential_energy, checked_pairs, calculated_pairs = physics.evaluate_cell_pair(pair, positions, forces, epsilon, sigma, r_cutoff)
            total_potential_energy += potential_energy
            all_checked_pairs.append(checked_pairs)
            all_calculated_pairs.append(calculated_pairs)

        epsilon24 = 24 * epsilon # energy factor
        # Multiply all forces by energy factor 24epsilon
        forces[:,1] = forces[:,1] * epsilon24
        forces[:,2] = forces[:,2] * epsilon24
        forces[:,3] = forces[:,3] * epsilon24

        # Update velocities and positions
        sum_v_squared = 0
        for i in range(number_atoms):
            # Calculate velocity V(t + 0.5dt)
            vxi = velocities[i, 1]
            vyi = velocities[i, 2]
            vzi = velocities[i, 3]

            fxi = forces[i, 1]
            fyi = forces[i, 2]
            fzi = forces[i, 3]

            delta_vxi = fxi * velocity_scale
            delta_vyi = fyi * velocity_scale
            delta_vzi = fzi * velocity_scale

            vxi2 = vxi + delta_vxi
            vyi2 = vyi + delta_vyi
            vzi2 = vzi + delta_vzi

            # Update positions
            positions[i, 1] += vxi2 * timestep_size
            positions[i, 2] += vyi2 * timestep_size
            positions[i, 3] += vzi2 * timestep_size

            # Calculate actual velocity at time t
            # For kinetic energy calculations only
            vxi3 = (vxi + vxi2) / 2
            vyi3 = (vyi + vyi2) / 2
            vzi3 = (vzi + vzi2) / 2
            sum_v_squared += vxi3**2 + vyi3**2 + vzi3**2

            # Update velocities
            velocities[i, 1] = vxi2
            velocities[i, 2] = vyi2
            velocities[i, 3] = vzi2

        if timestep % history_interval == 0 and timestep != 0:
            file_functions.write_history(positions, number_atoms, timestep, trajectory_filename) 
        
        total_kinetic_energy = physics.calculate_kinetic_energy(sum_v_squared)
        print(f"Timestep: {timestep} Potential Energy: {total_potential_energy} Kinetic Energy: {total_kinetic_energy} Total Energy: {total_potential_energy + total_kinetic_energy}")
        file_functions.write_energies_to_file(energies_filename, timestep, total_potential_energy, total_kinetic_energy)

    
    end = time.time()
    print(f"Run Time: {end - start}")
    file_functions.write_pair_checks_to_file(all_checked_pairs, all_calculated_pairs)


dd_method("case2.xyz")

