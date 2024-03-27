#using Statistics
#using Printf
#using LsqFit

#include("CH4.jl")
#filename_without_extension = "CH4"
include(ARGS[1])
filename_without_extension = replace(ARGS[1], ".jl" => "")  # remove the ".jl" extension

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
	@isdefined(dielectric)
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
function calc_num_atoms_cartesian(s::String = cartesian)
    #count(r"\n", s) + 1
    count(r"\n", strip(s)) + 1
end

function calc_num_atoms_zmatrix(s::String = atomlist)
    length(split(atomlist))
end

# Get the atomic label of the `i`th atom
function get_atomlabel_cartesian(i_atom::Int64, s::String = cartesian)
    index = (i_atom-1)*4 + 1
    split(s)[index]
end

function get_atomlabel_zmatrix(i_atom::Int64, s::String = atomlist)
    split(s)[i_atom]
end

#= # Get the xyz coordinates of the `i`th atom
function get_atomcoor(i_atom::Int64, s::String = cartesian)
    index = (i_atom-1)*4 + 1
    tidygeom = strip(s)
    tidygeom[index+1], tidygeom[index+2], tidygeom[index+3]
end

function print_line(io::IO, type::String, i_atom::Int64, 𝑓::Float64)
    if @isdefined(cartesian)
        atomlable = get_atomlabel(i_atom)
        atomcoor = get_atomcoor(i_atom)
        radius = get_atom_radius(atomlable)
        if type == "structure"
            println(io, atomlable, "    ", atomcoor[1], " ", atomcoor[2], " ", atomcoor[3])
        elseif type in ("Vc", "Ger")
            println(io, atomcoor[1], " ", atomcoor[2], " ", atomcoor[3], "    ", radius, "    ", 𝑓)
        elseif type == "Gcav"
            if sphere == "hard"
                𝑓 = scalingfactors[1]
            end
            println(io, n, "    ", radius * 𝑓, "    1.0")
        elseif type == "opt"
            println(io, i_atom, "    ", radius, "    ", 𝑓)
        end
    end
    if @isdefined(zmatrix)
        println(io, i_atom, "    ", get_atom_radius(split(atomlist)[i_atom]), "    ", 𝑓)
    end
end
# @time print_line(stdout, "structure", 1, 1, 1.2)

function print_structure(io::IO, type::String, 𝑓::Float64, noa::Int64 = numatoms)
    for i_atom in 1:noa
        print_line(io, type, i_atom, 𝑓)
    end
end =#

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
        println(io, i_𝑓 == 1 ? "" : "%kjob l502\n", "%chk=$filename_without_extension-Ger.chk")
        println(io, "%nproc=",np)
        println(io, "%mem=",mem)
    elseif jobtype == "opt"
        println(io, "%subst l301 $exedir")
        println(io, "%subst l502 $exedir")
        println(io, "%subst l701 $exedir")
        println(io, "%chk=$filename_without_extension-opt.chk")
        println(io, "%nproc=",np)
        println(io, "%mem=",mem)
    end
end

function print_route(io::IO, jobtype::String, i_𝑓::Int64, kws::String = keywords, sol::String = solvent)
    println(io, "#p ", kws, ((jobtype == "Ger" || jobtype == "opt") && i_𝑓 > 1) ? " guess=read" : "")
    if jobtype == "opt"
        if @isdefined(cartesian)
            println(io, "#p opt")
        elseif @isdefined(zmatrix)
            println(io, "#p popt=(z-matrix)")
        end
    end
    println(io, "#p scrf=(iefpcm,solvent=",sol,",read) nosymm 6d 10f")
end

function print_title(io::IO, jobtype::String, i_𝑓::Int64, 𝑓list = scalingfactors)
    println(io, jobtype, " calculation with scalingfactor = ", 𝑓list[i_𝑓])
end

function print_mol_spec(io::IO, chrg::Int64 = charge, mulplct::Int64 = multiplicity, noa::Int64 = numatoms)
    println(io, chrg, " ", mulplct)
    if @isdefined(cartesian)
        s = split(cartesian)
        for i_atom in 1:noa
            println(io, s[4*(i_atom-1)+1], "    ", s[4*(i_atom-1)+2], "  ", s[4*(i_atom-1)+3], "  ", s[4*(i_atom-1)+4])
        end
    elseif @isdefined(zmatrix)
        println(io, strip(zmatrix))
    end
end

function print_pcm_spec(io::IO, jobtype::String, i_𝑓::Int64, noa::Int64 = numatoms, tsare::Float64 = tesserae, cav::String = cavity, sp::Solvent = get_sol_params())
    # determine the number of spheres
    if cavity == "custom"
        nsfe = size(spherespec, 1)
    else
        nsfe = noa
    end
    # determine if addsph is needed and what to write on the nsfe line
    if cav in ("vdw", "custom") 
        nsfeline = "nsfe=$nsfe noaddsph"
    elseif cav == "ses"
        nsfeline = "nsfe=$nsfe addsph rsolv=$(sp.𝑟)"
    end
    # print for differnt jobtype
    if jobtype == "Vc"
        println(io, "pcmdoc geomview g03defaults tsare=", tsare)
    elseif jobtype in ("Ger", "opt")
        println(io, "qrep pcmdoc geomview nodis nocav g03defaults tsare=", tsare)
        println(io, "nvesolv=", sp.𝑛, " solvmw=", sp.𝑀)
        println(io, "eps=", 𝜀[i_𝑓], " rhos=", 𝜌[i_𝑓])  # 𝜀 and 𝜌 are global variables of 1D array of length nosf
        println(io, "sten=", float(𝜂))
        if jobtype == "Ger"
            println(io, "cmf=0")
        else # i.e. jobtype == "opt"
            println(io, "cmf=100")
            println(io, "dsten=", 𝑝[i_𝑓])
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
        for i_atom in 1:numatoms
            radius = get_atom_radius(atomlist_array[i_atom])
            println(io, i_atom, "    ", radius, "    ", 𝑓list[i_𝑓])
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
    open("$filename_without_extension-$jobtype.gjf", "w") do file
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
    run(`$gau $filename_without_extension-$jobtype.gjf`)
end

# extract data from Gaussian .log files
function get_data(jobtype::String, searchstring::String, fieldnum::Int64, nosf::Int64 = length(scalingfactors))
    data = Vector{Float64}(undef, nosf)    # 1D array
    i_𝑓 = 1    # i_𝑓 ranges from 1:nosf
    open("$filename_without_extension-$jobtype.log", "r") do file
        for line in eachline(file)
            if occursin(searchstring, line)
                data[i_𝑓] = parse(Float64, split(line)[fieldnum])
                i_𝑓 += 1
            end
        end
    end
    data
end

#------------------------------------------------------------------------------
# main.jl
#------------------------------------------------------------------------------
#function main()
# Input error checking
if @isdefined(cartesian) && @isdefined(zmatrix)
    error("Both `cartesian` and `zmatrix` are found; use one only.")
end
if @isdefined(zmatrix) && !@isdefined(atomlist)
    error("No atom list found; provide `atomlist` for the zmatrix.")
end
if cavity == "custom" && !@isdefined(spherespec)
    error("No custom sphere specification found; provide `spherespec` for the custom cavity.")
end

if @isdefined(cartesian)
    numatoms = calc_num_atoms_cartesian(cartesian)
    atomlist_array = [get_atomlabel_cartesian(i_atom) for i_atom in 1:numatoms]
elseif @isdefined(zmatrix)
    numatoms = calc_num_atoms_zmatrix(atomlist)
    atomlist_array = [get_atomlabel_zmatrix(i_atom) for i_atom in 1:numatoms] 
end

# Step 1: cavity volume 𝑉𝑐(𝑓) and solvent property calculations
write_gjf("Vc")
run_gaussian("Vc")
𝑉𝑐 = get_data("Vc", "Cavity volume", 5)

𝑠 = @. ∛(𝑉𝑐/𝑉𝑐[1])  # linear scaling factor 𝑠, as the cubic root of the volume ratio

𝜀₀ = get_sol_params().𝜀
𝜀 = @. 1 + (𝜀₀ - 1) / 𝑠^3  # dielectric permitivity 𝜀 = 1 + (𝜀₀-1)/𝑠̄³

𝜌₀ = get_sol_params().𝜌
𝜌 = @. 𝜌₀ / 𝑠^(3+𝜂)  # solvent density 𝜌 = 𝜌₀/𝑠̄⁽³⁺𝜂⁾

# Step 2: electronic structure Gaussian jobs and pressure calculations
write_gjf("Ger")
run_gaussian("Ger")
#𝐺𝑒𝑟 = get_data("Ger", "SCF Done", 5)
𝑝 = get_data("Ger", "-dG/dV", 3)  # in a.u.; 1 Ha/bohr³ = 29421.0471 GPa

# Step 3: geometry opt (Ger + pVc) at constant pressure
write_gjf("opt")
#run_gaussian("opt")

#end  # function main

#main()
