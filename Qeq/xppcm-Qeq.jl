# Julia Script: Data Extraction and High-Pressure Equilibrium Geometry Calculations
# Usage:
# julia xppcm-Qeq.jl 

using Printf

include("input-Qeq.jl")
#filename_without_extension = "CH4"
#include(ARGS[1])
#filename_without_extension = replace(ARGS[1], ".jl" => "")  # remove the ".jl" extension

#------------------------------------------------------------------------------
# Solvents.jl
#------------------------------------------------------------------------------
struct Solvent <: Real #FieldVector{5, Real} #
	𝜀::Float64  # dielectric
	𝜌::Float64  # solvent density
	𝑀::Float64  # molar mass
	𝑛::Int64    # number of valence electrons
	𝑟::Float64  # molecular radius
end

function Solvent(s::String)
	if s == "cyclohexane"
		Solvent(2.0165, 0.7781, 84.1595, 36, 2.815)
	elseif s == "benzene"
		Solvent(2.2706, 0.8756, 78.1118, 30, 2.63 )
	elseif s == "argon"
		Solvent(1.43  , 1.3954, 39.948 ,  8, 1.705)
	else
		error("solvent not implemented. Try cyclohexane, benzene, or argon.")
	end
end
	#implemented = ("cyclohexane", "benzene", "argon")
	#cyclohexane = Solvent(2.0165, 0.7781, 84.1595, 36, 2.815)
	#benzene     = Solvent(2.2706, 0.8756, 78.1118, 30, 2.63 )
	#argon       = Solvent(1.43  , 1.3954, 39.948 ,  8, 1.705)

function get_sol_params(s::String = solvent)
	sol = Solvent(s)
	if @isdefined(dielectric)
		Solvent(dielectric, sol.𝜌, sol.𝑀, sol.𝑛, sol.𝑟)
	else
		sol
	end
end

#------------------------------------------------------------------------------
# AtomicRadii.jl
#------------------------------------------------------------------------------
bondi = Base.ImmutableDict(
    #dummy atom
    "X" => 0.0,
    #1s
    "H" => 1.20,    "1" => 1.20,
    "He"=> 1.40,    "2" => 1.40,
    #2s
    "Li"=> 1.82,    "3" => 1.82,
    "Be"=> 1.53,    "4" => 1.53,
    #2p
    "B" => 1.92,    "5" => 1.92,
    "C" => 1.70,    "6" => 1.70,
    "N" => 1.55,    "7" => 1.55,
    "O" => 1.52,    "8" => 1.52,
    "F" => 1.47,    "9" => 1.47,
    "Ne"=> 1.54,    "10"=> 1.54,
    #3s
    "Na"=> 2.27,    "11"=> 2.27,
    "Mg"=> 1.73,    "12"=> 1.73,
    #3p
    "Al"=> 1.84,    "13"=> 1.84,
    "Si"=> 2.10,    "14"=> 2.10,
    "P" => 1.80,    "15"=> 1.80,
    "S" => 1.80,    "16"=> 1.80,
    "Cl"=> 1.75,    "17"=> 1.75,
    "Ar"=> 1.88,    "18"=> 1.88,
    "K" => 2.75,    "19"=> 2.75,
    "Ca"=> 2.31,    "20"=> 2.31,
    #4p
    "Ga"=> 1.87,    "31"=> 1.87,
    "Ge"=> 2.11,    "32"=> 2.11,
    "As"=> 1.85,    "33"=> 1.85,
    "Se"=> 1.90,    "34"=> 1.90,
    "Br"=> 1.85,    "35"=> 1.85,
    "Kr"=> 2.02,    "36"=> 2.02,
    #5p
    "In"=> 1.93,    "49"=> 1.93,
    "Sn"=> 2.17,    "50"=> 2.17,
    "Sb"=> 2.06,    "51"=> 2.06,
    "Te"=> 2.06,    "52"=> 2.06,
    "I" => 1.98,    "53"=> 1.98,
    "Xe"=> 2.16,    "54"=> 2.16,
    #6p
    "Tl"=> 1.96,    "81"=> 1.96,
    "Pb"=> 2.02,    "82"=> 2.02,
    "Bi"=> 2.07,    "83"=> 2.07,
    "Po"=> 1.97,    "84"=> 1.97,
    "At"=> 2.02,    "85"=> 2.02,
    "Rn"=> 2.20,    "86"=> 2.20 )
# add more if needed from https://en.wikipedia.org/wiki/Van_der_Waals_radius

rahm = Base.ImmutableDict(
    #dummy atom
    "X" => 0.0,
    #1s
    "H" => 1.54,    "1" => 1.54,
    "He"=> 1.34,    "2" => 1.34,
    #2s
    "Li"=> 2.20,    "3" => 2.20,
    "Be"=> 2.19,    "4" => 2.19,
    #2p
    "B" => 2.05,    "5" => 2.05,
    "C" => 1.90,    "6" => 1.90,
    "N" => 1.79,    "7" => 1.79,
    "O" => 1.71,    "8" => 1.71,
    "F" => 1.63,    "9" => 1.63,
    "Ne"=> 1.56,    "10"=> 1.56,
    #3s
    "Na"=> 2.25,    "11"=> 2.25,
    "Mg"=> 2.40,    "12"=> 2.40,
    #3p
    "Al"=> 2.39,    "13"=> 2.39,
    "Si"=> 2.32,    "14"=> 2.32,
    "P" => 2.23,    "15"=> 2.23,
    "S" => 2.14,    "16"=> 2.14,
    "Cl"=> 2.06,    "17"=> 2.06,
    "Ar"=> 1.97,    "18"=> 1.97,
    #4s
    "K" => 2.34,    "19"=> 2.34,
    "Ca"=> 2.70,    "20"=> 2.70,
    #3d
    "Sc"=> 2.63,    "21"=> 2.63,
    "Ti"=> 2.57,    "22"=> 2.57,
    "V" => 2.51,    "23"=> 2.51,
    "Cr"=> 2.32,    "24"=> 2.32,
    "Mn"=> 2.42,    "25"=> 2.42,
    "Fe"=> 2.37,    "26"=> 2.37,
    "Co"=> 2.33,    "27"=> 2.33,
    "Ni"=> 2.19,    "28"=> 2.19,
    "Cu"=> 2.16,    "29"=> 2.16,
    "Zn"=> 2.22,    "30"=> 2.22,
    #4p
    "Ga"=> 2.33,    "31"=> 2.33,
    "Ge"=> 2.34,    "32"=> 2.34,
    "As"=> 2.29,    "33"=> 2.29,
    "Se"=> 2.25,    "34"=> 2.25,
    "Br"=> 2.19,    "35"=> 2.19,
    "Kr"=> 2.12,    "36"=> 2.12,
    #5s
    "Rb"=> 2.40,    "37"=> 2.40,
    "Sr"=> 2.79,    "38"=> 2.79,
    #4d
    "Y" => 2.74,    "39"=> 2.74,
    "Zr"=> 2.69,    "40"=> 2.69,
    "Nb"=> 2.51,    "41"=> 2.51,
    "Mo"=> 2.44,    "42"=> 2.44,
    "Tc"=> 2.52,    "43"=> 2.52,
    "Ru"=> 2.37,    "44"=> 2.37,
    "Rh"=> 2.33,    "45"=> 2.33,
    "Pd"=> 2.15,    "46"=> 2.15,
    "Ag"=> 2.25,    "47"=> 2.25,
    "Cd"=> 2.38,    "48"=> 2.38,
    #5p
    "In"=> 2.46,    "49"=> 2.46,
    "Sn"=> 2.48,    "50"=> 2.48,
    "Sb"=> 2.46,    "51"=> 2.46,
    "Te"=> 2.42,    "52"=> 2.42,
    "I" => 2.38,    "53"=> 2.38,
    "Xe"=> 2.32,    "54"=> 2.32,
    #6p
    "Tl"=> 2.42,    "81"=> 2.42,
    "Pb"=> 2.49,    "82"=> 2.49,
    "Bi"=> 2.50,    "83"=> 2.50,
    "Po"=> 2.50,    "84"=> 2.50,
    "At"=> 2.47,    "85"=> 2.47,
    "Rn"=> 2.43,    "86"=> 2.43 )
# add more if needed from https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/chem.201700610

function get_atom_radius(atom::AbstractString, 
                         type::String = @isdefined(radiustype) ? radiustype : "bondi",
                         r::Base.ImmutableDict{String,Float64} = eval(Symbol(type)))
    r[atom]
end

#------------------------------------------------------------------------------
# Geometry.jl
#------------------------------------------------------------------------------
# Get number of atomss
function calc_num_atoms_cartesian(s::String = cartesian)
    #count(r"\n", s) + 1
    count(r"\n", strip(s)) + 1
end

function calc_num_atoms_zmatrix(s::String = atomlist)
    length(split(s))
end

function calc_num_atoms()
    if @isdefined(cartesian)
        calc_num_atoms_cartesian()
    elseif @isdefined(zmatrix)
        calc_num_atoms_zmatrix()
    end
end

# Get the atomic label of the `i`th atom
function get_atomlabel_cartesian(i_atom::Int64, s::String = cartesian)
    index = (i_atom-1)*4 + 1
    split(s)[index]
end

function get_atomlabel_zmatrix(i_atom::Int64, s::String = atomlist)
    split(s)[i_atom]
end

#------------------------------------------------------------------------------
# Gaussian.jl
#------------------------------------------------------------------------------

# Gaussian input sections are explained here: https://gaussian.com/input/?tabid=0
# Link0, Route, Title, Molecule Specification, etc. sections
function print_link0(io::IO, jobtype::String, i_𝑓::Int64, np::Int64 = nproc, mem::String = mem)
    if jobtype == "Vc"
        println(io, "%kjob l301")
        println(io, "%nproc=1")
        println(io, "%mem=1gb")
    elseif jobtype == "Ger"
        println(io, "%subst l301 $exedir")
        println(io, "%subst l502 $exedir")
        println(io, "%subst l701 $exedir")
        println(io, i_𝑓 == 1 ? "" : "%kjob l502\n", "%chk=Ger.chk")
        println(io, "%nproc=",np)
        println(io, "%mem=",mem)
    elseif jobtype == "force"
        println(io, "%subst l301 $exedir")
        println(io, "%subst l502 $exedir")
        println(io, "%subst l701 $exedir")
        println(io, "%chk=force.chk")
        println(io, "%nproc=",np)
        println(io, "%mem=",mem)
    end
end

function print_route(io::IO, jobtype::String, i_𝑓::Int64, kws::String = keywords, sol::String = solvent)
    println(io, "#p ", kws, ((jobtype == "Ger" || jobtype == "force") && i_𝑓 > 1) ? " guess=read" : "")
    if jobtype == "force"
        #if @isdefined(cartesian)
            println(io, "#p force")
        #elseif @isdefined(zmatrix)
        #    println(io, "#p force=(z-matrix)")
        #end
    end
    println(io, "#p scrf=(iefpcm,solvent=",sol,",read) nosymm 6d 10f")
end

function print_title(io::IO, jobtype::String, i_𝑓::Int64, 𝑓list = scalingfactors)
    println(io, jobtype, " calculation with scalingfactor = ", 𝑓list[i_𝑓])
end

function print_mol_spec(io::IO, chrg::Int64 = charge, mulplct::Int64 = multiplicity)
    println(io, chrg, " ", mulplct)
    if @isdefined(cartesian)
        s = split(cartesian)
        for i_atom in 1:calc_num_atoms_cartesian()
            println(io, s[4*(i_atom-1)+1], "    ", s[4*(i_atom-1)+2], "  ", s[4*(i_atom-1)+3], "  ", s[4*(i_atom-1)+4])
        end
    elseif @isdefined(zmatrix)
        println(io, strip(zmatrix))
    end
end

function print_pcm_spec(io::IO, jobtype::String, i_𝑓::Int64, cav::String = cavity, sp::Solvent = get_sol_params())
    # determine the number of spheres
    if cav == "custom"
        nsfe = size(spherespec, 1)
    else
        nsfe = calc_num_atoms()
    end
    # determine whether to use surface charge smoothing
    if @isdefined(tesserae) # no smoothing
        smoothing = "g03defaults tsare=$tesserae"
    elseif @isdefined(pdens) # York-Karplus smoothing
        smoothing = "pdens=$pdens"
    end
    # determine if addsph is needed and what to write on the nsfe line
    if cav in ("vdw", "custom") 
        nsfeline = "nsfe=$nsfe noaddsph"
    elseif cav == "ses"
        nsfeline = "nsfe=$nsfe addsph rsolv=$(sp.𝑟)"
    end
    # print for different jobtypes
    if jobtype == "Vc"
        println(io, "qrep pcmdoc geomview nodis nocav ", smoothing)
    elseif jobtype in ("Ger", "force")
        println(io, "qrep pcmdoc geomview nodis nocav ", smoothing)
        println(io, "nvesolv=", sp.𝑛, " solvmw=", sp.𝑀)
        #println(io, "eps=", 𝜀[i_𝑓], " rhos=", 𝜌[i_𝑓])  # 𝜀 and 𝜌 are global variables of 1D array of length nosf
        println(io, "eps=2.0165 rhos=0.7781")  # 𝜀 and 𝜌 at f=1.2
        println(io, "sten=", float(𝜂))
        if jobtype == "Ger"
            println(io, "cmf=0")
        else # i.e. jobtype == "force"
            println(io, "cmf=100")
            #println(io, "dsten=", 𝑝[i_𝑓])
            println(io, "dsten=0.00003398927294") # 1GPa in a.u. as a random pressure
            println(io, "tce=1.0") # tce=1.0 using Cavity step function theory for volume gradients
        end
    end
    println(io, nsfeline)
end

function print_sphere_spec(io::IO, i_𝑓::Int64, 𝑓list = scalingfactors)
    if cavity == "custom"
        for i_sph in axes(spherespec, 1) # axes() gives 1:number_of_spheres
            println(io, Int(spherespec[i_sph,1]), "    ", spherespec[i_sph,2], "    ", 𝑓list[i_𝑓])
        end
    else # i.e. cavity in ("vdw", "ses")
        if @isdefined(cartesian)
            for i_atom in 1:calc_num_atoms_cartesian()
                radius = get_atom_radius(get_atomlabel_cartesian(i_atom))
                println(io, i_atom, "    ", radius, "    ", 𝑓list[i_𝑓])
            end
        elseif @isdefined(zmatrix)
            for i_atom in 1:calc_num_atoms_zmatrix()
                radius = get_atom_radius(get_atomlabel_zmatrix(i_atom))
                println(io, i_atom, "    ", radius, "    ", 𝑓list[i_𝑓])
            end
        end
    end
end

# combining the above pieces
function print_content(io::IO, jobtype::String, i_𝑓::Int64, nosf::Int64 = length(scalingfactors))
    print_link0(io, jobtype, i_𝑓)
    print_route(io, jobtype, i_𝑓)
    println(io)
    print_title(io, jobtype, i_𝑓)
    println(io)
    print_mol_spec(io)
    println(io)
    print_pcm_spec(io, jobtype, i_𝑓)
    println(io)
    print_sphere_spec(io, i_𝑓)
    println(io)
    if i_𝑓 != nosf; println(io, "--link1--"); end
end

function write_gjf(jobtype::String, nosf::Int64 = length(scalingfactors))
    open("$jobtype.gjf", "w") do file
        for i_𝑓 in 1:nosf
            print_content(file, jobtype, i_𝑓)
        end
    end
end

# check if g16 or g09 is installed and loaded
function get_gau_ver()
    if occursin("atlas", gethostname())
        return "g16"
    elseif typeof(Sys.which("g16")) === String
        return "g16"
    elseif typeof(Sys.which("g09")) === String
        return "g09"
    else
        error("Command `g16` or `g09` not found.")
    end
end

function run_gaussian(jobtype::String)
    gau = get_gau_ver()
    run(`$gau $jobtype.gjf`)
end

# extract data from Gaussian .log files
function get_data(jobtype::String, searchstring::String, fieldnum::Int64, nosf::Int64 = length(scalingfactors))
    data = Vector{Float64}(undef, nosf)    # 1D array
    i_𝑓 = 1    # i_𝑓 ranges from 1:nosf
    open("$jobtype.log", "r") do file
        for line in eachline(file)
            if occursin(searchstring, line)
                data[i_𝑓] = parse(Float64, split(line)[fieldnum])
                i_𝑓 += 1
            end
        end
    end
    data
end

# Input error checking
function chk_input_error()
    if @isdefined(cartesian) && @isdefined(zmatrix)
        error("Both `cartesian` and `zmatrix` are found; use one only.")
    end
    if @isdefined(zmatrix) && !@isdefined(atomlist)
        error("No atom list found; provide `atomlist` for the zmatrix.")
    end
    if cavity == "custom" && !@isdefined(spherespec)
        error("No custom sphere specification found; provide `spherespec` for the custom cavity.")
    end
    if @isdefined(tesserae) && @isdefined(pdens)
        error("Both `tesserae` and `pdens` are found; use one only.")
    end
    if !@isdefined(tesserae) && !@isdefined(pdens)
        error("No `tesserae` or `pdens` found; define one.")
    end
end

#--------------------
# xppcm-Qeq-v1.jl
#--------------------
function extract_normal_modes(filename, num_atoms, num_modes)
    # Initialize output variables
    displacements = Float64[]
    symmetry_labels = String[]
    frequencies = Float64[]
    force_constants = Float64[]
    formatted_displacements_matrix = zeros(3, num_atoms, num_modes)

    open(filename, "r") do io
        frequency_section_flag = false
        hpmodes_flag = false
        mode_num_counter = 0
        for line in eachline(io)
            if hpmodes_flag == false # i.e., normal freq output format
                if occursin(r"hpmodes"i, line)
                    hpmodes_flag = true
                    continue  # skip the rest and go to the next line
                end
            end
            if frequency_section_flag == false # i.e., not in the right section
                # Look for the start of the harmonic frequencies section and reset the line counter to 1
                if occursin("Harmonic frequencies (cm**-1)", line)
                    frequency_section_flag = true
                end
                continue  # skip the rest and go to the next line
            end
            if hpmodes_flag == true
                # Extract Cartesian displacements using regex matching, 2 integers + 3 float numbers
                if occursin(r"^\s*[1-3]\s+\d+\s+\d+\s+[-+]?\d\.\d\d\d\d\d(?:\s+[-+]?\d\.\d\d\d\d\d){0,4}$", line)
                    coord_component_string = split(line)[1]
                    atom_num_string = split(line)[2]
                    displacement_strings = split(line)[4:end]
                    coord_component_value = parse(Int64, coord_component_string)
                    atom_num_value = parse(Int64, atom_num_string)
                    displacement_values = parse.(Float64, displacement_strings)
                    for i in eachindex(displacement_values)
                        formatted_displacements_matrix[coord_component_value, atom_num_value, mode_num_counter+i] = displacement_values[i]
                    end
                    if atom_num_value == num_atoms && coord_component_value == 3 # the last line of the displacement for this mode block
                        mode_num_counter += 5
                    end
                # Extract the symmetry labels using regex matching, 3 symmetry labels with 2 or 3 (upper letters or number)
                elseif occursin(r"^\s*[A-Z][A-Z0-9]{0,2}(?:\s+[A-Z][A-Z0-9]{0,2}){0,4}$", line)
                    symmetry_values = split(line)
                    append!(symmetry_labels, symmetry_values)
                # Extract frequencies
                elseif occursin("Frequencies ---", line)
                    # Skip the "Frequencies" and "---" tokens
                    freq_strings = split(line)[3:end]
                    freq_values = parse.(Float64, freq_strings)
                    append!(frequencies, freq_values)
                # Extract force constant
                elseif occursin("Force constants ---", line)
                    # Skip the "Force", "constants" and "---" tokens
                    force_constants_strings = split(line)[4:end]
                    force_constants_values = parse.(Float64, force_constants_strings)
                    append!(force_constants, force_constants_values)
                # Exit the `for line in eachline(io)` loop; i.e., stop reading file
                elseif occursin("Harmonic frequencies (cm**-1)", line) # stop when reading the string 2nd time.
                    break
                end
            else # hpmodes is off
                # Extract Cartesian displacements using regex matching, 2 integers + 3 float numbers
                if occursin(r"^\s*\d+\s+\d+\s+[-+]?\d\.\d\d(?:\s+[-+]?\d\.\d\d){8}$", line)
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
        if hpmodes_flag == false
            N = reshape(displacements, 9, :)
            for i in 1:num_modes
                formatted_displacements_matrix[:,:,i] = N[((i-1)%3*3+1):((i-1)%3*3+3), ((ceil(Int, i/3)-1)*num_atoms+1):ceil(Int, i/3)*num_atoms]
            end
        end
    end
    # change force constants to negative for imaginary frequencies
    imaginary_freq_indices = findall(x -> x < 0, frequencies)
    for i in imaginary_freq_indices
        force_constants[i] *= -1
    end
    return formatted_displacements_matrix, frequencies, symmetry_labels, force_constants
end

#=
function transform_normal_modes(displacements_vector, num_atoms, num_modes)
    formatted_matrix = zeros(3, num_atoms, num_modes)
    N = reshape(displacements_vector, 9, :)
    for i in 1:num_modes
        formatted_matrix[:,:,i] = N[((i-1)%3*3+1):((i-1)%3*3+3), ((ceil(Int, i/3)-1)*num_atoms+1):ceil(Int, i/3)*num_atoms]
    end
    return formatted_matrix
end
=#
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

#= 
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

# Always use the volume gradients at the first pressure (or cavity) for all pressures (i.e., cavities)
function extract_pressures_and_volume_gradients_v2(filename, num_atoms)
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
        # Use the volume gradients at the first pressure (or cavity) for all pressures (i.e., cavities)
        DRV[:, j] = DRV_list[1]
    end

    return pressure_values, DRV
end
 =#

# Function to extract the pressures from force.log output
function extract_pressures(filename)
    pressure_values = Float64[]
    for line in eachline(filename) 
        if occursin("p(au)/p(GPa)=", line)
            # Extract pressure value
            pressure_line = line
            pressure_value = nothing
            # Extract the value after '=' and before '/'
            pressure_parts = split(pressure_line, "=")
            if length(pressure_parts) >= 2
                pressure_str = strip(split(pressure_parts[2], "/")[1])
                pressure_value = parse(Float64, pressure_str)
                push!(pressure_values, pressure_value)
            else
                error("Cannot parse pressure value in line: $line")
            end
        end
    end
    return pressure_values
end

# Function to extract pressures and volume gradients from a different output file
function extract_volume_gradients(filename, num_atoms)
    DRV_list = Float64[]

    open(filename, "r") do io
        lines = readlines(io)
        i = 1
        while i <= length(lines)
            line = lines[i]
            if occursin("Cavity step function theory", line)
                # Now read Vx, Vy, Vz values for each atom
                #volume_gradients = Float64[]
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
                        push!(DRV_list, value)
                        read_lines += 1
                    else
                        error("Cannot parse volume gradient in line: $vg_line")
                    end
                end
                # Append volume_gradients to DRV_list
                #push!(DRV_list, volume_gradients)
            else
                i += 1
            end
        end
    end

    # Now, convert DRV_list to a 2D array
    #num_coords = num_atoms * 3
    #DRV = zeros(num_coords, num_pressures)
    #for j in 1:num_pressures
        # Use the volume gradients at the first pressure (or cavity) for all pressures (i.e., cavities)
        #DRV[:, j] = DRV_list[1]
    #end

    return DRV_list
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
#--------------------
# end of xppcm-Qeq-v1.jl
#--------------------





#------------------------------------------------------------------------------
# main.jl
#------------------------------------------------------------------------------
function main()
    chk_input_error()

    # Step 1: cavity volume 𝑉𝑐(𝑓) and solvent property calculations
    #write_gjf("Vc")
    #run_gaussian("Vc")
    #𝑉𝑐 = get_data("Vc", "Cavity volume", 5)

    #𝑠 = @. ∛(𝑉𝑐/𝑉𝑐[1])  # linear scaling factor 𝑠, as the cubic root of the volume ratio

    #𝜀₀ = get_sol_params().𝜀
    #global 𝜀 = @. 1 + (𝜀₀ - 1) / 𝑠^3  # dielectric permitivity 𝜀 = 1 + (𝜀₀-1)/𝑠̄³

    #𝜌₀ = get_sol_params().𝜌
    #global 𝜌 = @. 𝜌₀ / 𝑠^(3+𝜂)  # solvent density 𝜌 = 𝜌₀/𝑠̄⁽³⁺𝜂⁾

    # Step 2: electronic structure Gaussian jobs and pressure calculations
    #write_gjf("Ger")
    #run_gaussian("Ger")
    #𝐺𝑒𝑟 = get_data("Ger", "SCF Done", 5)
    #global 𝑝 = get_data("Ger", "-dG/dV", 3)  # in a.u.; 1 Ha/bohr³ = 29421.0471 GPa

    # Step 3: volume gradients, force calculation
    write_gjf("force")
    run_gaussian("force")





    # Step 4: equilibrium geometry calculations
    #--------------------
    # xppcm-Qeq-v1.jl
    #--------------------
    # ---------------------------
    # Input
    # ---------------------------

    # Replace 'your_gaussian_output.log' with the path to your Gaussian output file
    #req_calculation_output_file = ARGS[1]

    # Replace 'your_pressure_output.log' with the path to the pressure output file
    force_calculation_output_file = "force.log"

    # List the indices of the totally symmetric normal modes in the gas phase calculation
    #totally_symmetric_mode_indices = eval(Meta.parse(ARGS[3]))

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
    normal_modes, frequencies, symmetry_labels, force_constants_mdyne_per_angstrom = extract_normal_modes(freq_calculation_output_file, num_atoms, num_modes)

    #normal_modes = transform_normal_modes(displacements, num_atoms, num_modes)

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
    #pressure_values, DRV = extract_pressures_and_volume_gradients(force_calculation_output_file, num_atoms)
    #pressure_values, DRV = extract_pressures_and_volume_gradients_v2(force_calculation_output_file, num_atoms)
    #pressure_values = extract_pressures(force_calculation_output_file)
    pressure_values = pressure_values_GPa/29421.0471
    num_pressures = length(pressure_values)
    DRV = extract_volume_gradients(force_calculation_output_file, num_atoms)

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
            @printf("%-6s %9.5f %9.5f %9.5f\n", atom_labels[j], mode_matrix[j, 1], mode_matrix[j, 2], mode_matrix[j, 3])
        end
        println("\n")
    end

    # pressures
    println("Pressures:")
    for J in 1:num_pressures
        @printf("p[%s]:  %10.3f GPa    %.8e Eh/a₀³\n", J, pressure_values[J]*29421.0471, pressure_values[J])
    end
    println("\n")

    # Volume gradients
    println("Volume gradients (a₀²):")
    #println(DRV)
    #for J in 1:num_pressures
        #println("Pressure p[$J]: $(pressure_values[J]) Eh/a₀³")
        #DRV_J = DRV[:, J]
        #DRV_matrix = reshape(DRV_J, 3, num_atoms)'
        DRV_matrix = reshape(DRV, 3, num_atoms)'
        println("Atom     Vx          Vy          Vz")
        println("----------------------------------------")
        for i in 1:num_atoms
            @printf("%-6s %11.6f %11.6f %11.6f\n", atom_labels[i], DRV_matrix[i, 1], DRV_matrix[i, 2], DRV_matrix[i, 3])
        end
        println("\n")
    #end

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
    # DQV is 1D vector of length num_modes
    DQV = transpose(U) * DRV  # U^T * DRV
    compression_propensity_vector = -DQV ./ force_constants_au


    # Print compression_propensity_vector before adjustment for totally symmetric modes
    println("-DQV/k (sorted and in a₀⁴/Eh):")
    println("* indicates the modes used in Qeq calculation")
    if length(compression_propensity_vector) == length(symmetry_labels)
        sorted_indices = sortperm(compression_propensity_vector, by=abs, rev=true)  # Sort indices by absolute value in descending order
        @printf("              -DQV/k       DQV       k      symm_label\n")
        for i in sorted_indices
            @printf("Mode %3d: %10.3f %10.3f %10.3f      %s%s\n", i, compression_propensity_vector[i], DQV[i], force_constants_au[i], symmetry_labels[i], i in totally_symmetric_mode_indices ? "   *" : "")
        end
        #for (i, k) in enumerate(compression_propensity_vector)
        #    @printf("Mode %3d: %10.3f   %s\n", i, k, symmetry_labels[i])
        #end
        println("Modes used in Qeq calculation: $totally_symmetric_mode_indices")
        println("\n")
    end

    # ---------------------------
    # Adjust for Totally Symmetric Modes
    # ---------------------------

    # Set the volume gradient to zero for modes that are not totally symmetric
    #for j in 1:num_pressures
        for i in 1:num_modes
            if !(i in totally_symmetric_mode_indices)
            #if !is_totally_symmetric[i]
                #DQV[i, j] = 0.0
                DQV[i] = 0.0
            end
        end
    #end

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
        normal_mode_displacements[:, J] = compression_propensity_vector * pressure_values[J]
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
        println("Pressure p[$J]: $(pressure_values_GPa[J]) GPa")
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
        # Print equilibrium geometries in xyz format
        open("p$J.xyz", "w") do file
            println(file, num_atoms)
            @printf(file, "Pressure:  %.3f GPa    %.8e Eh/a₀³\n", pressure_values_GPa[J], pressure_values[J])
            for i in 1:num_atoms
                @printf(file, "%-6s %12.6f %12.6f %12.6f\n", atom_labels[i], new_geometry[i, 1], new_geometry[i, 2], new_geometry[i, 3])
            end
        end
    end
    #--------------------
    # end of xppcm-Qeq-v1.jl
    #--------------------

end  # function main

main()
