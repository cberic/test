# Julia Script: Data Extraction and High-Pressure Equilibrium Geometry Calculations
# Usage:

# run the script with julia, passing as commaned line arguments the two file names 
# (gasphase freq calculation and high pressure force calculation) and the 
# indexes of totally_symmetric_modes of the gas phase freq calculation (separated
# by comma and no space allowed)

# julia xppcm-Qeq-v1.jl freq.log force.log 1,5,17


using Printf

# ---------------------------
# Data Extraction Functions
# ---------------------------
#= example frequency and normal mode section of Gaussian output
"""
 Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering
 activities (A**4/AMU), depolarization ratios for plane and unpolarized
 incident light, reduced masses (AMU), force constants (mDyne/A),
 and normal coordinates:
                      1                      2                      3
                     A2U                    A1U                    EU
 Frequencies --  -1816.4345                38.4199               322.6374
 Red. masses --      1.0962                 1.0078                 1.0116
 Frc consts  --      2.1309                 0.0009                 0.0620
 IR Inten    --    154.6514                 0.0000                 1.6377
  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z
     1   1    -0.00  -0.00   0.97    -0.00  -0.00  -0.00    -0.00   0.44   0.00
     2   6     0.00  -0.00  -0.06    -0.00  -0.00  -0.00    -0.00   0.01  -0.00
     3   6    -0.00  -0.00  -0.06    -0.00   0.00   0.00    -0.00   0.01  -0.00
     4   1     0.00  -0.03   0.09     0.41   0.00  -0.00     0.00  -0.12   0.49
     5   1     0.00  -0.03   0.09     0.41  -0.00   0.00     0.00  -0.12   0.49
     6   1     0.03   0.01   0.09    -0.20   0.35   0.00     0.00  -0.13  -0.24
     7   1    -0.03   0.01   0.09    -0.20  -0.35  -0.00    -0.00  -0.13  -0.24
     8   1    -0.03   0.01   0.09    -0.20  -0.35  -0.00    -0.00  -0.13  -0.24
     9   1     0.03   0.01   0.09    -0.20   0.35   0.00     0.00  -0.13  -0.24
                      4                      5                      6
                     EU                     A1G                    EG
 Frequencies --    322.6379               525.8519               701.1675
 Red. masses --      1.0116                 3.3246                 1.1595
 Frc consts  --      0.0620                 0.5416                 0.3359
 IR Inten    --      1.6375                 0.0000                 0.0000
  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z
     1   1     0.44   0.00   0.00     0.00   0.00   0.00     0.00  -0.00   0.00
     2   6     0.01   0.00  -0.00    -0.00   0.00   0.32    -0.00   0.08  -0.00
     3   6     0.01   0.00  -0.00     0.00  -0.00  -0.32     0.00  -0.08   0.00
     4   1    -0.13  -0.00   0.00     0.00  -0.01   0.36     0.00  -0.07   0.56
     5   1    -0.13  -0.00   0.00     0.00   0.01  -0.36    -0.00   0.07  -0.56
     6   1    -0.12   0.00  -0.42     0.01   0.00   0.36     0.02  -0.10  -0.28
     7   1    -0.12  -0.00   0.42     0.01  -0.00  -0.36     0.02   0.10   0.28
     8   1    -0.12  -0.00   0.42    -0.01   0.00   0.36    -0.02  -0.10  -0.28
     9   1    -0.12   0.00  -0.42    -0.01  -0.00  -0.36    -0.02   0.10   0.28
"""
=#
function extract_normal_modes(filename)
    # Initialize output variables
    displacements = Float64[]
    symmetry_labels = String[]
    frequencies = Float64[]
    force_constants = Float64[]

    open(filename, "r") do io
        frequency_section_flag = false
        for line in eachline(io)
            if frequency_section_flag == false # i.e., not in the right section
                # Look for the start of the harmonic frequencies section and reset the line counter to 1
                if occursin("Harmonic frequencies (cm**-1)", line)
                    frequency_section_flag = true
                end
                continue  # skip the rest and go to the next line
            end
            # Extract Cartesian displacements using regex matching, 2 integers + 3 float numbers
            if occursin(r"^\s*\d\s+\d\s+[-+]?\d\.\d\d(?:\s+[-+]?\d\.\d\d){8}$", line)
                displacement_strings = split(line)[3:end]
                displacement_values = parse.(Float64, displacement_strings)
                append!(displacements, displacement_values)
            # Extract the symmetry labels using regex matching, 3 symmetry labels with 2 or 3 (upper letters or number)
            elseif occursin(r"^\s*[A-Z][A-Z0-9]{0,2}(?:\s+[A-Z][A-Z0-9]{0,2}){2}$", line)
                symmetry_values = split(line)
                append!(symmetry_labels, symmetry_values)
            # Extract frequencies
            elseif occursin("Frequencies --", line)
                # Skip the "Frequencies" and "--" tokens
                freq_strings = split(line)[3:end]
                freq_values = parse.(Float64, freq_strings)
                append!(frequencies, freq_values)
            # Extract force constant
            elseif occursin("Frc consts  --", line)
                # Skip the "Frc", "consts" and "--" tokens
                force_constants_strings = split(line)[4:end]
                force_constants_values = parse.(Float64, force_constants_strings)
                append!(force_constants, force_constants_values)
            # Exit the `for line in eachline(io)` loop; i.e., stop reading file
            elseif occursin("- Thermochemistry -", line)
                break
            end
        end
    end
    return displacements, frequencies, symmetry_labels, force_constants
end

function transform_normal_modes(displacements_vector, num_atoms, num_modes)
    formatted_matrix = zeros(3, num_atoms, num_modes)
    N = reshape(displacements_vector, 9, :)
    for i in 1:num_modes
        formatted_matrix[:,:,i] = N[((i-1)%3*3+1):((i-1)%3*3+3), ((ceil(Int, i/3)-1)*num_atoms+1):ceil(Int, i/3)*num_atoms]
    end
    return formatted_matrix
end

# Function to extract the initial geometry as a num_atoms × 3 matrix
function extract_initial_geometry(filename)
    atom_numbers = Int64[]
    geometry_data = Float64[]

    open(filename, "r") do io
        lines = readlines(io)
        # Find all indices where "Standard orientation" or "Input orientation" occurs
        orientation_indices = []
        orientation_type = ""

        for (i, line) in enumerate(lines)
            if occursin("Standard orientation", line)
                push!(orientation_indices, i)
                orientation_type = "Standard orientation"
            elseif occursin("Input orientation", line)
                if isempty(orientation_indices)
                    push!(orientation_indices, i)
                    orientation_type = "Input orientation"
                end
            end
        end

        if isempty(orientation_indices)
            error("No 'Standard orientation' or 'Input orientation' section found in the file.")
        end

        # Start reading from the last orientation section found
        start_index = orientation_indices[end] + 5  # Skip header lines
        i = start_index
        while i <= length(lines)
            line = strip(lines[i])
            if isempty(line) || occursin("----", line)
                # End of geometry section
                break
            else
                data = split(replace(line, r"\s+" => " "))
                if length(data) >= 6
                    atom_num = parse(Int, data[2])
                    x = parse(Float64, data[4])
                    y = parse(Float64, data[5])
                    z = parse(Float64, data[6])
                    push!(atom_numbers, atom_num)
                    # Collect the coordinates in a flat array
                    push!(geometry_data, x)
                    push!(geometry_data, y)
                    push!(geometry_data, z)
                end
                i += 1
            end
        end
    end

    num_atoms = length(atom_numbers)
    # Convert the flat array into a 2D matrix of size num_atoms × 3
    geometry_matrix = reshape(geometry_data, (3, num_atoms))'
    # The transpose ensures the shape is num_atoms × 3

    return atom_numbers, geometry_matrix, num_atoms
end

# Function to extract pressures and volume gradients from a different output file
function extract_pressures_and_volume_gradients(filename, num_atoms)
    pressure_values = Float64[]
    DRV_list = []

    open(filename, "r") do io
        lines = readlines(io)
        i = 1
        while i <= length(lines)
            line = lines[i]
            if occursin("p(au)/p(GPa)=", line)
                # Extract pressure value
                pressure_line = line
                pressure_value = nothing
                # Extract the value after '=' and before '/'
                pressure_parts = split(pressure_line, "=")
                if length(pressure_parts) >= 2
                    pressure_str = strip(split(pressure_parts[2], "/")[1])
                    # Replace D with E in case of scientific notation
                    pressure_str = replace(pressure_str, "D" => "E")
                    pressure_value = parse(Float64, pressure_str)
                    push!(pressure_values, pressure_value)
                else
                    error("Cannot parse pressure value in line: $line")
                end

                # Now read Vx, Vy, Vz values for each atom
                volume_gradients = Float64[]
                read_lines = 0
                while read_lines < num_atoms * 3 && i + 1 <= length(lines)
                    i += 1
                    vg_line = strip(lines[i])
                    if isempty(vg_line)
                        continue
                    end
                    # Each line should be like 'Vx    value' or 'Vy    value' or 'Vz    value'
                    # Split the line
                    data = split(vg_line)
                    if length(data) >= 2
                        value_str = data[2]
                        # Replace D with E if necessary
                        value_str = replace(value_str, "D" => "E")
                        value = parse(Float64, value_str)
                        push!(volume_gradients, value)
                        read_lines += 1
                    else
                        error("Cannot parse volume gradient in line: $vg_line")
                    end
                end
                # Append volume_gradients to DRV_list
                push!(DRV_list, volume_gradients)
            else
                i += 1
            end
        end
    end

    # Now, convert DRV_list to a 2D array
    num_pressures = length(pressure_values)
    num_coords = num_atoms * 3
    DRV = zeros(num_coords, num_pressures)
    for j in 1:num_pressures
        DRV[:, j] = DRV_list[j]
    end

    return pressure_values, DRV
end
#=
# Function to extract the point group from Gaussian output
function extract_point_group(filename)
    point_group = ""
    open(filename, "r") do io
        for line in eachline(io)
            if occursin("Full point group", line)
                parts = split(line)
                # The point group is the 4th element (assuming line format: ' Full point group           D2h     NOp   8')
                if length(parts) >= 4
                    point_group = parts[4]
                    break
                end
            end
        end
    end
    if point_group == ""
        error("Could not find point group in Gaussian output.")
    end
    return point_group
end

# Mapping of point groups to their totally symmetric representations
const totally_symmetric_repr = Dict(
    "C1" => "A",
    "Cs" => "A'",
    "Ci" => "AG",
    "C2" => "A",
    "C2v" => "A1",
    "C3v" => "A1",
    "D2h" => "A1G",
    # Add more point groups as needed
)
=#
# Function to get atom labels from atomic numbers
function get_atom_labels(atom_numbers)
    atomic_symbols = Dict(
        1 => "H",    2 => "He",   3 => "Li",   4 => "Be",   5 => "B",    6 => "C",    7 => "N",
        8 => "O",    9 => "F",    10 => "Ne",  11 => "Na",  12 => "Mg",  13 => "Al",  14 => "Si",
        15 => "P",   16 => "S",   17 => "Cl",  18 => "Ar",  19 => "K",   20 => "Ca",  21 => "Sc",
        22 => "Ti",  23 => "V",   24 => "Cr",  25 => "Mn",  26 => "Fe",  27 => "Co",  28 => "Ni",
        29 => "Cu",  30 => "Zn",  31 => "Ga",  32 => "Ge",  33 => "As",  34 => "Se",  35 => "Br",
        36 => "Kr",  37 => "Rb",  38 => "Sr",  39 => "Y",   40 => "Zr",  41 => "Nb",  42 => "Mo",
        43 => "Tc",  44 => "Ru",  45 => "Rh",  46 => "Pd",  47 => "Ag",  48 => "Cd",  49 => "In",
        50 => "Sn",  51 => "Sb",  52 => "Te",  53 => "I",   54 => "Xe",  55 => "Cs",  56 => "Ba",
        57 => "La",  58 => "Ce",  59 => "Pr",  60 => "Nd",  61 => "Pm",  62 => "Sm",  63 => "Eu",
        64 => "Gd",  65 => "Tb",  66 => "Dy",  67 => "Ho",  68 => "Er",  69 => "Tm",  70 => "Yb",
        71 => "Lu",  72 => "Hf",  73 => "Ta",  74 => "W",   75 => "Re",  76 => "Os",  77 => "Ir",
        78 => "Pt",  79 => "Au",  80 => "Hg",  81 => "Tl",  82 => "Pb",  83 => "Bi",  84 => "Po",
        85 => "At",  86 => "Rn",  87 => "Fr",  88 => "Ra",  89 => "Ac",  90 => "Th",  91 => "Pa",
        92 => "U",   93 => "Np",  94 => "Pu",  95 => "Am",  96 => "Cm",  97 => "Bk",  98 => "Cf",
        99 => "Es", 100 => "Fm", 101 => "Md", 102 => "No", 103 => "Lr", 104 => "Rf", 105 => "Db",
       106 => "Sg", 107 => "Bh", 108 => "Hs", 109 => "Mt", 110 => "Ds", 111 => "Rg", 112 => "Cn",
       113 => "Nh", 114 => "Fl", 115 => "Mc", 116 => "Lv", 117 => "Ts", 118 => "Og"
    )
    atom_labels = [atomic_symbols[i] for i in atom_numbers]
    return atom_labels
end

# ---------------------------
# Main Computational Code
# ---------------------------

function main()
    # ---------------------------
    # Input
    # ---------------------------

    # Replace 'your_gaussian_output.log' with the path to your Gaussian output file
    freq_calculation_output_file = ARGS[1]

    # Replace 'your_pressure_output.log' with the path to the pressure output file
    force_calculation_output_file = ARGS[2]

    # List the indexes of the totally symmetric normal modes in the gas phase calculation
    totally_symmetric_mode_indexes = eval(Meta.parse(ARGS[3]))

    # ---------------------------
    # Data Extraction
    # ---------------------------

    # Extract initial geometry as a matrix
    atom_numbers, initial_geometry, num_atoms = extract_initial_geometry(freq_calculation_output_file)

    # Number of Cartesian coordinates
    num_coords = num_atoms * 3

    # Number of vibrational normal modes
    num_modes = num_coords - 6  # For nonlinear molecules

    # Extract freq and normal modes
    displacements, frequencies, symmetry_labels, force_constants_mdyne_per_angstrom = extract_normal_modes(freq_calculation_output_file)

    normal_modes = transform_normal_modes(displacements, num_atoms, num_modes)

    # Extract the point group
    #point_group = extract_point_group(freq_calculation_output_file)

    # Get the totally symmetric representation
    #if haskey(totally_symmetric_repr, point_group)
    #    ts_repr = totally_symmetric_repr[point_group]
    #else
    #    error("Point group $point_group not recognized in totally symmetric representation mapping.")
    #end

    # Identify the totally symmetric modes
    #is_totally_symmetric = [sym_label == ts_repr for sym_label in symmetry_labels]

    # Extract pressures and volume gradients from pressure output file
    # DRV is a 2D matrix with dimensions num_coords * num_pressures
    pressure_values, DRV = extract_pressures_and_volume_gradients(force_calculation_output_file, num_atoms)
    num_pressures = length(pressure_values)

    # ---------------------------
    # Data Reporting
    # ---------------------------

    # Number of atoms
    println("Number of atoms: $num_atoms\n")

    # Initial geometry
    atom_labels = get_atom_labels(atom_numbers)
    println("Initial geometry (Å):")
    println("Atom       x (Å)        y (Å)        z (Å)")
    println("---------------------------------------------")
    for i in 1:num_atoms
        @printf("%-6s %12.6f %12.6f %12.6f\n", atom_labels[i], initial_geometry[i, 1], initial_geometry[i, 2], initial_geometry[i, 3])
    end
    println("\n")

    # Force constants
    println("Force constants (mDyne/Å):")
    for (i, k) in enumerate(force_constants_mdyne_per_angstrom)
        @printf("Mode %3d: %10.6f mDyne/Å\n", i, k)
    end
    println("\n")

    # Normal modes
    println("Normal modes:")
    for i in axes(normal_modes, 3)
        println("Mode $i:")
        println("Frequency: $(frequencies[i]) cm⁻¹")
        println("Symmetry label: $(symmetry_labels[i])")
        println("Displacement matrix:")
        mode_matrix = normal_modes[:,:,i]'
        println("Atom     Δx         Δy         Δz")
        println("------------------------------------")
        for j in 1:num_atoms
            @printf("%-6s %9.6f %9.6f %9.6f\n", atom_labels[j], mode_matrix[j, 1], mode_matrix[j, 2], mode_matrix[j, 3])
        end
        println("\n")
    end

    # Volume gradients
    println("Volume gradients (a₀²):")
    for J in 1:num_pressures
        println("Pressure p[$J]: $(pressure_values[J]) Eh/a₀³")
        DRV_J = DRV[:, J]
        DRV_matrix = reshape(DRV_J, 3, num_atoms)'
        println("Atom     Vx          Vy          Vz")
        println("----------------------------------------")
        for i in 1:num_atoms
            @printf("%-6s %11.6f %11.6f %11.6f\n", atom_labels[i], DRV_matrix[i, 1], DRV_matrix[i, 2], DRV_matrix[i, 3])
        end
        println("\n")
    end

    # ---------------------------
    # Computational Steps
    # ---------------------------

    # Check that we have read the correct number of normal modes
    if size(normal_modes)[3] != num_modes
        error("Number of normal modes read does not match expected number.")
    end

    # Construct the transformation matrix U (num_coords x num_modes)
    U = zeros(num_coords, num_modes)

    for i in 1:num_modes
        # Each normal mode is a (num_atoms x 3) matrix
        mode_matrix = normal_modes[:,:,i]
        # Check dimensions
        if size(mode_matrix) != (3, num_atoms)
            error("Normal mode $i has incorrect dimensions.")
        end
        # Flatten the mode_matrix into a column vector of length num_coords
        # The ordering is [x1, y1, z1, x2, y2, z2, ..., xN, yN, zN]
        mode_vector = vec(mode_matrix)
        # Assign to the transformation matrix U
        U[:, i] = mode_vector
    end

    # ---------------------------
    # Force Constants
    # ---------------------------

    # Check length
    if length(force_constants_mdyne_per_angstrom) != num_modes
        error("Force constants data has incorrect length.")
    end

    # Check zero value
    if 0.0 in force_constants_mdyne_per_angstrom
        error("Force constants contain zero value(s).")
    end

    # Convert force constants to atomic units (Eh/a₀²)
    # Transformation factor: 1 mDyne/Å = 0.064231 Eh/a₀²
    force_constants_au = force_constants_mdyne_per_angstrom * 0.064231

    # ---------------------------
    # Compute DQV
    # ---------------------------

    # Compute the normal modes gradient vector DQV of the compression cavity volume
    # DQV has dimensions (num_modes x num_pressures)
    DQV = transpose(U) * DRV  # U^T * DRV

    # ---------------------------
    # Adjust for Totally Symmetric Modes
    # ---------------------------

    # Set the volume gradient to zero for modes that are not totally symmetric
    for j in 1:num_pressures
        for i in 1:num_modes
            if !(i in totally_symmetric_mode_indexes)
            #if !is_totally_symmetric[i]
                DQV[i, j] = 0.0
            end
        end
    end

    # ---------------------------
    # Compute Compression Propensity Vector
    # ---------------------------

    # Compute the normal modes compression propensity vector DQDP (Eh·a₀⁴)
    #compression_propensity_vector = zeros(num_modes, num_pressures)

    compression_propensity_vector = -DQV ./ force_constants_au

    # ---------------------------
    # Compute Normal Mode Displacements
    # ---------------------------

    # Compute the compression along the normal modes DQ at various pressures (a₀)
    normal_mode_displacements = zeros(num_modes, num_pressures)
    for J in 1:num_pressures
        normal_mode_displacements[:, J] = compression_propensity_vector[:, J] * pressure_values[J]
    end

    # ---------------------------
    # Compute Atomic Displacements
    # ---------------------------

    # Compute the compression along the Cartesian coordinates of atoms at the pressures (Å)
    atomic_displacements_angstrom = zeros(num_coords, num_pressures)
    for J in 1:num_pressures
        atomic_displacements_angstrom[:, J] = U * normal_mode_displacements[:, J] * 0.529177  # Convert a₀ to Å
    end

    # ---------------------------
    # Output Results
    # ---------------------------
    println("Compresseion matrices and new equiliurm geometries:")
    for J in 1:num_pressures
        #println("\n-----------------------------------------")
        @printf("Pressure p[%d] (Eh/a₀³): %.6e\n", J, pressure_values[J])
        println("Compression matrix at p[$J] (Å):")
        DR_ang_matrix = reshape(atomic_displacements_angstrom[:, J], 3, num_atoms)'
        # Print header
        println("Atom     Δx (Å)      Δy (Å)      Δz (Å)")
        println("-----------------------------------------")
        for i in 1:num_atoms
            @printf("%-6s %10.6f %10.6f %10.6f\n", atom_labels[i], DR_ang_matrix[i, 1], DR_ang_matrix[i, 2], DR_ang_matrix[i, 3])
        end
        println("\nEquilibrium geometry at p[$J] (Å):")
        new_geometry = initial_geometry + DR_ang_matrix
        # Print header
        println("Atom       x (Å)        y (Å)        z (Å)")
        println("-----------------------------------------")
        for i in 1:num_atoms
            @printf("%-6s %12.6f %12.6f %12.6f\n", atom_labels[i], new_geometry[i, 1], new_geometry[i, 2], new_geometry[i, 3])
        end
        println();println()
    end
    #println("\n-----------------------------------------\n")
end

# Call the main function
main()