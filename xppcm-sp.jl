# =============================================================================
#                                  xppcm.jl
# -----------------------------------------------------------------------------
#                                Roberto Cammi
#                   Department of Chemical Science (SCVSA)
#                           University of Parma, ITALY

#                                   Bo Chen
#                  Donostia International Physics Center, SPAIN
# -----------------------------------------------------------------------------
#
# xppcm.jl is a Julia script to conduct XP-PCM calculations for
# studying chemical reactions under pressure [1,2].
# References:
#         [1] R. Cammi, J. Comp. Chem., 36, 2246-2259 (2015)
#         [2] B. Chen. R. Hoffmann and R. Cammi,
#            Angew. Chem. Int. Ed. 56, 11126-11142 (2017) and related
#                       Supporting Information (SI)
#         [3] R. Cammi, B. Chen, M. Rahm, J. Comp. Chem. 39, 2243 (2018)
#         [4] R. Cammi, J. Chem. Phys.,150,164122 (2019)
#         [5] M. Rahm, R. Cammi, N.W. Ashcroft. R. Hoffman, JACS,   (2020)
#
# Note: The presents script is under development and its use is confidential.
# =============================================================================

using Statistics
using Printf
using LsqFit
include("input.jl")

#------------------------------------------------------------------------------
# Solvents.jl
#------------------------------------------------------------------------------
struct Solvent <: Real #FieldVector{5, Real} #
	𝜀::Float64  # dielectric
	𝜌::Float64  # valence electron density
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
#@time Solvent("benzene")

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
# Geometries.jl
#------------------------------------------------------------------------------
function calc_num_structs(s::String = geometries)
    count(r"\n+\s*\n+", strip(s)) + 1
end
#@time calc_num_structs()

function split_geoms(s::String = geometries)
    split(strip(s), r"\n+\s*\n+", keepempty=false)
end

function calc_num_atoms(s::AbstractString)
    count(r"\n", s) + 1
end

geomdata = split(geometries)
numatoms_vector = calc_num_atoms.(split_geoms()) # 1D array collecting numatom of each structure

# Get the atomic label of the `iₐₜₒₘ`th atom in the `iₛₜᵣᵤ`th structure from `geomdata`
function get_atomlabel(iₐₜₒₘ::Int64, iₛₜᵣᵤ::Int64, geom::Vector{SubString{String}} = geomdata, noa::Vector{Int64} = numatoms_vector)
    #nos = length(noa)
    index = sum(noa[i]*4 for i in 1:iₛₜᵣᵤ) - noa[iₛₜᵣᵤ]*4 + (iₐₜₒₘ-1)*4 + 1
    geom[index]
end
# @time get_atomlabel(1,1)

# Get the xyz coordinates of the `iₐₜₒₘ`th atom in the `iₛₜᵣᵤ`th structure from `geomdata`
function get_atomcoor(iₐₜₒₘ::Int64, iₛₜᵣᵤ::Int64, geom::Vector{SubString{String}} = geomdata, noa::Vector{Int64} = numatoms_vector)
    #nos = length(noa)
    index = sum(noa[i]*4 for i in 1:iₛₜᵣᵤ) - noa[iₛₜᵣᵤ]*4 + (iₐₜₒₘ-1)*4 + 1  # the index of the atomlabel 
    geom[index+1], geom[index+2], geom[index+3]  # the following three elements in the array are the xyz coordinates
end
# @time get_atomcoor.(1:18,103)

function print_line(io::IO, type::String, iₐₜₒₘ::Int64, iₛₜᵣᵤ::Int64, 𝑓::Float64, Gcav_sph::String = Gcav_sphere)
    atomlable = get_atomlabel(iₐₜₒₘ, iₛₜᵣᵤ)
    atomcoor = get_atomcoor(iₐₜₒₘ, iₛₜᵣᵤ)
    atomradius = get_atom_radius(atomlable)
    if type == "structure"
        println(io, atomlable, "    ", atomcoor[1], "  ", atomcoor[2], "  ", atomcoor[3])
    elseif type in ("Vc", "Ger")
        println(io, iₐₜₒₘ, "    ", atomradius, "    ", 𝑓)
        #println(io, atomcoor[1], "  ", atomcoor[2], "  ", atomcoor[3], "    ", atomradius, "    ", 𝑓)
    elseif type == "Gcav"
        if Gcav_sph == "hard" # alwasys use the first scaling factor, otherwise (i.e. sphere=="soft"), use appropriate scaling factors
            𝑓 = scalingfactors[1]
        end
        println(io, iₐₜₒₘ, "    ", atomradius * 𝑓, "    1.0")
        # In Gcav calculation, the following line separating the radius and scaling factor does not work.
        # The scaling factor seems to be ignored.
        #println(io, iₐₜₒₘ, "    ", atomradius, "    ", 𝑓)
    end
end
# @time print_line(stdout, "structure", 1, 1, 1.2)

function print_structure(io::IO, type::String, iₛₜᵣᵤ::Int64, 𝑓::Float64, noa::Vector{Int64} = numatoms_vector)
    for iₐₜₒₘ in 1:noa[iₛₜᵣᵤ]
        print_line(io, type, iₐₜₒₘ, iₛₜᵣᵤ, 𝑓)
    end
end

#------------------------------------------------------------------------------
# Gaussian.jl
#------------------------------------------------------------------------------
function calc_num_scalingfactors(𝑓list = scalingfactors)
    length(𝑓list)
end

# Gaussian input sections are explained here: https://gaussian.com/input/?tabid=0
# Link0, Route, Title, Molecule Specification, etc. sections
function print_link0(io::IO, jobtype::String, iₛₜᵣᵤ::Int64, i_𝑓::Int64, np::Int64 = nproc, mem::String = mem)
    if jobtype !== "Ger"  # jobtype == "Vc" or "Gcav"
        println(io, "%kjob l301")
        println(io, "%nproc=1")
        println(io, "%mem=1gb")
    else  # jobtype == "Ger"
        println(io, i_𝑓 == 1 ? "" : "%kjob l502\n", "%chk=structure-", iₛₜᵣᵤ, "-Ger.chk")
        println(io, "%nproc=", np)
        println(io, "%mem=", mem)
    end
end

function print_route(io::IO, jobtype::String, i_𝑓::Int64, kws::String = keywords, sol::String = solvent)
    println(io, "#p ", kws, (jobtype == "Ger" && i_𝑓 > 1) ? " guess=read" : "")
    println(io, "#p scrf=(iefpcm,solvent=", sol, ",read) nosymm 6d 10f")
end

function print_title(io::IO, jobtype::String, i_𝑓::Int64, 𝑓list = scalingfactors)
    println(io, jobtype, " calculation with scalingfactor = ", 𝑓list[i_𝑓])
end

function print_mol_spec(io::IO, iₛₜᵣᵤ::Int64, i_𝑓::Int64, chrg::Int64 = charge, mulplct::Int64 = multiplicity, 𝑓list = scalingfactors)
    println(io, chrg," ",mulplct)
    print_structure(io, "structure", iₛₜᵣᵤ, 𝑓list[i_𝑓])
end

function print_pcm_spec(io::IO, jobtype::String, iₛₜᵣᵤ::Int64, i_𝑓::Int64, tsare::Float64 = tesserae, noa::Vector{Int64} = numatoms_vector, cav::String = cavity)
    sp = get_sol_params()
    if jobtype == "Vc"
        println(io, "pcmdoc geomview g03defaults tsare=",tsare)
        println(io, "nsfe=",noa[iₛₜᵣᵤ], cav == "vdw" ? " noaddsph" : " rsolv=$(sp.𝑟)")
    elseif jobtype == "Ger"
        #𝜀 = calc_𝜀()    # data of 𝜀 and 𝜌 needed for the Ger gjf files
        #𝜌 = calc_𝜌()
        println(io, "qrep pcmdoc geomview nodis nocav g03defaults tsare=",tsare)
        println(io, "nsfe=",noa[iₛₜᵣᵤ], cav == "vdw" ? " noaddsph" : " rsolv=$(sp.𝑟)")
        println(io, "nvesolv=",sp.𝑛," solvmw=",sp.𝑀)
        println(io, "eps=",𝜀[i_𝑓]," rhos=",𝜌[i_𝑓])
    elseif jobtype == "Gcav"
        #𝑉ₘ = calc_𝑉ₘ()    # molar volume 𝑉ₘ of the solvent
        println(io, "norep nodis cav g03defaults tsare=",tsare)
        println(io, "nsfe=",noa[iₛₜᵣᵤ], cav == "vdw" ? " noaddsph" : " rsolv=$(sp.𝑟)")
        println(io, "Vmol=",𝑉ₘ[i_𝑓])
    end
end

function print_sphere_spec(io::IO, jobtype::String, iₛₜᵣᵤ::Int64, i_𝑓::Int64, 𝑓list = scalingfactors)
    print_structure(io, jobtype, iₛₜᵣᵤ, 𝑓list[i_𝑓])
end

# combining the above pieces
function print_content(io::IO, jobtype::String, iₛₜᵣᵤ::Int64, i_𝑓::Int64)
    nosf = calc_num_scalingfactors()
    print_link0(io, jobtype, iₛₜᵣᵤ, i_𝑓)
    print_route(io, jobtype, i_𝑓)
    println(io)
    print_title(io, jobtype, i_𝑓)
    println(io)
    print_mol_spec(io, iₛₜᵣᵤ, i_𝑓)
    println(io)
    print_pcm_spec(io, jobtype, iₛₜᵣᵤ, i_𝑓)
    println(io)
    print_sphere_spec(io, jobtype, iₛₜᵣᵤ, i_𝑓)
    println(io)
    if i_𝑓 != nosf println(io, "--link1--") end
end

function write_gjf(jobtype::String)
    nos = calc_num_structs()
    nosf = calc_num_scalingfactors()
    Threads.@threads for iₛₜᵣᵤ in 1:nos  # use multithreading
        open("tmp/structure-$iₛₜᵣᵤ-$jobtype.gjf", "w") do file
            for i_𝑓 in 1:nosf
                print_content(file, jobtype, iₛₜᵣᵤ, i_𝑓)
            end
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

function submit_job(gau::String, i_job::Int64, jobtype::String)
    run(`$gau structure-$i_job-$jobtype.gjf`)
end

function submit_jobs(gau::String, jobnums::Union{Vector{Int64},UnitRange{Int64}}, jobtype::String)
    Threads.@threads for i_job in jobnums
        run(`$gau structure-$i_job-$jobtype.gjf`)
    end
end

function run_gaussian(jobtype::String, jobnums::Union{Vector{Int64},UnitRange{Int64}} = 1:calc_num_structs(), ismulthrd::Bool = ismultithreading)
    gau = get_gau_ver()
    cd("tmp")
    # Only when ismultithreading = false and only for Ger jobs, do not use multithreading 
    if !ismulthrd #&& jobtype == "Ger" 
        submit_job.(gau, jobnums, jobtype)
    else  # is multithreading
        submit_jobs(gau, jobnums, jobtype)
    end
    cd("..")
end

function find_unfinished_jobs(jobtype::String)
    nos = calc_num_structs()
    nosf = calc_num_scalingfactors()
    searchstringdict = Dict("Vc" => "Cavity volume", "Ger" => "After PCM corrections", "Gcav" => "PCM non-electrostatic energy")
    unfinished = Int64[]   # empty array to collect the unfinished job numbers
    Threads.@threads for i_job in 1:nos
        try # try open file
            file = read("tmp/structure-$i_job-$jobtype.log", String)
            n = count(searchstringdict[jobtype], file)
            if n != nosf
                push!(unfinished,i_job)  # collect i_job if the number of energies found is != nosf
            end
        catch  # file not found; collect the file number
            push!(unfinished,i_job)
        end
    end
    unfinished
end

function restart_jobs(jobtype::String)
    unfinished = find_unfinished_jobs(jobtype)
    if !isempty(unfinished)
        run_gaussian(jobtype, unfinished)
    end
end

# extract data from Gaussian .log files
function get_data(jobtype::String, searchstring::String, fieldnum::Int64)
    nos = calc_num_structs()
    nosf = calc_num_scalingfactors()
    data = Matrix{Float64}(undef, nos,nosf)    # nos * nosf 2D array
    Threads.@threads for iₛₜᵣᵤ in 1:nos
        i_𝑓 = 1    # i_𝑓 ranges from 1:nosf
        open("tmp/structure-$iₛₜᵣᵤ-$jobtype.log", "r") do file
            for line in eachline(file)
                if occursin(searchstring, line)
                    data[iₛₜᵣᵤ,i_𝑓] = parse(Float64, split(line)[fieldnum])
                    i_𝑓 += 1
                end 
            end
        end
    end
    data
end

#𝑉𝑐 = get_data("Vc", "Cavity volume", 5)
#𝐺𝑒𝑟 = get_data("Ger", "After PCM corrections", 7)
#𝑉𝑐𝑎𝑣 = get_data("Gcav", "Cavity volume", 5) # 𝑉𝑐𝑎𝑣 could be different from 𝑉𝑐 in Gcav calculation using hard sphere
#𝐺𝑛𝑜𝑛_𝑝𝑉 = get_data("Gcav", "PCM non-electrostatic energy", 5)


#------------------------------------------------------------------------------
# write_properties.jl
#------------------------------------------------------------------------------
function write_properties(𝑓list, 
                        𝑉𝑐::Matrix{Float64}, 
                        𝑠::Matrix{Float64}, 
                        𝑠̄::Matrix{Float64}, 
                        𝜀::Matrix{Float64}, 
                        𝜌::Matrix{Float64}, 
                        𝑉ₘ::Matrix{Float64}, 
                        𝐺𝑒𝑟::Matrix{Float64}, 
                        𝑝::Matrix{Float64}, 
                        𝑝̄::Matrix{Float64}, 
                        𝑉𝑐𝑎𝑣::Matrix{Float64}, 
                        𝐺𝑛𝑜𝑛_𝑝𝑉::Matrix{Float64}, 
                        𝐺𝑐𝑎𝑣::Matrix{Float64}, 
                        𝐺𝑡𝑜𝑡::Matrix{Float64}, 
                        Δ𝐺𝑡𝑜𝑡::Matrix{Float64})
    nos = calc_num_structs()
    nosf = calc_num_scalingfactors()
    #𝑠 = calc_𝑠()
    #𝑠̄ = mean(𝑠, dims=2)
    #𝜀 = calc_𝜀()
    #𝜌 = calc_𝜌()
    #𝑉ₘ = calc_𝑉ₘ()
    #𝑝 = calc_𝑝()
    #𝑝̄ = mean(𝑝, dims=2)
    #𝐺𝑐𝑎𝑣 = calc_𝐺𝑐𝑎𝑣()
    #𝐺𝑡𝑜𝑡 = calc_𝐺𝑡𝑜𝑡()
    #Δ𝐺𝑡𝑜𝑡 = calc_Δ𝐺𝑡𝑜𝑡()
    𝐸𝑔𝑎𝑠 = get_data("Ger", "<psi(f)|   H    |psi(f)>", 6) # gas phase energy with the solvation field wave function
    𝐸ₑₗₑₛₜₐₜ = get_data("Ger", "(Polarized solute)-Solvent", 5)
    𝐸ₚₐᵤₗᵢ = get_data("Ger", "Quantum repulsion energy", 6)
    open("properties.dat", "w") do file
        for i in 1:nos
            println(file, "structure $i")
            println(file, "#      𝑓       𝑉𝑐       𝑠       𝑠̄       𝜀     𝜌ₛₒₗ       𝑉ₘ            𝐸𝑔𝑎𝑠  𝐸ₑₗₑₛₜₐₜ    𝐸ₚₐᵤₗᵢ             𝐺𝑒𝑟        𝑝        𝑝̄     𝑉𝑐𝑎𝑣        𝑝̄𝑉𝑐𝑎𝑣     𝐺𝑛𝑜𝑛_𝑝𝑉         𝐺𝑐𝑎𝑣            𝐺𝑡𝑜𝑡     Δ𝐺𝑡𝑜𝑡")
            println(file, "               Å³                                                        Eₕ  kcal/mol  kcal/mol              Eₕ      GPa      GPa       Å³           Eₕ          Eₕ           Eₕ              Eₕ  kcal/mol")
            for j in 1:nosf
                pv = 𝑝̄[j] * 𝑉𝑐𝑎𝑣[i,j] * 2.293712569e-4
                @printf(file, "%-2d %5.3f  %7.3f  %6.4f  %6.4f  %6.4f  %7.4f  %7.3f  %14.8f    %6.2f    %6.2f  %14.8f  %7.3f  %7.3f  %7.3f  %11.8f  %10.8f  %11.8f  %14.8f   %7.2f\n", 
                        j, 𝑓list[j], 𝑉𝑐[i,j], 𝑠[i,j], 𝑠̄[j], 𝜀[j], 𝜌[j], 𝑉ₘ[j], 𝐸𝑔𝑎𝑠[i,j], 𝐸ₑₗₑₛₜₐₜ[i,j], 𝐸ₚₐᵤₗᵢ[i,j], 𝐺𝑒𝑟[i,j], 𝑝[i,j], 𝑝̄[j], 𝑉𝑐𝑎𝑣[i,j], pv, 𝐺𝑛𝑜𝑛_𝑝𝑉[i,j], 𝐺𝑐𝑎𝑣[i,j], 𝐺𝑡𝑜𝑡[i,j], Δ𝐺𝑡𝑜𝑡[i,j])
            end
            println(file)
        end
    end
end

#------------------------------------------------------------------------------
# Algebra.jl
#------------------------------------------------------------------------------
# # linear scaling factor 𝑠, as the cubic root of the volume scaling
# 𝑠 = @. ∛(𝑉𝑐/𝑉𝑐[:,1])  # nos * nosf 2D array

# # average of 𝑠 over all structures at the same scalingfactor 𝑓
# 𝑠̄ = mean(𝑠, dims=1)  # 1D array of length nosf

# # dielectric permitivity 𝜀 = 1 + (𝜀₀-1)/𝑠̄³
# 𝜀₀ = get_sol_params().𝜀
# 𝜀 = @. 1 + (𝜀₀ - 1) / 𝑠̄^3  # 1D array of length nosf

# # solvent density 𝜌 = 𝜌₀/𝑠̄⁽³⁺𝜂⁾
# 𝜌₀ = get_sol_params().𝜌
# 𝜌 = @. 𝜌₀ / 𝑠̄^(3+𝜂)  # 1D array of length nosf

# # molar volume of solvent 𝑉ₘ = (𝑀/𝜌₀) * 𝑠̄³
# 𝑀 = get_sol_params().𝑀
# 𝑉ₘ = @. (𝑀/𝜌₀) * 𝑠̄^3  # 1D array of length nosf

# Murnaghan equation of state fitting for pressure 𝑝 calculation
# using LsqFit
function eos_fitting(𝑉𝑐, 𝐺𝑒𝑟)
    nos = calc_num_structs()
    abc_parameters = Array{Float64}(undef, nos,3)
    Threads.@threads for iₛₜᵣᵤ in 1:nos
    # python: y = (a/b)*(1/x)**b+(a-c)*x; y is Ger, x is Vc
    # mathematica: a*x ((1/b)*(t[[1, 1]]/x)^(b + 1) + 1) - c*x
    # LsqFit: a=p[1], b=p[2], c=p[3], x is Vc
        @. model(x, p) = (p[1]/p[2])*x^(-p[2]) + (p[1]-p[3])*x
        xdata = 𝑉𝑐[iₛₜᵣᵤ,:] ./ 𝑉𝑐[iₛₜᵣᵤ,1]
        ydata = 𝐺𝑒𝑟[iₛₜᵣᵤ,:] .- 𝐺𝑒𝑟[iₛₜᵣᵤ,1]
        p0 = [0.0, 5.0, 0.0]
        fit = curve_fit(model, xdata, ydata, p0)
        abc_parameters[iₛₜᵣᵤ,1] = fit.param[1]/𝑉𝑐[iₛₜᵣᵤ,1]
        abc_parameters[iₛₜᵣᵤ,2] = fit.param[2]
        abc_parameters[iₛₜᵣᵤ,3] = fit.param[3]/𝑉𝑐[iₛₜᵣᵤ,1]
    end
    return abc_parameters    # nos * 3 2D array
end

function calc_𝑝(𝑉𝑐, 𝐺𝑒𝑟)
    abc = eos_fitting(𝑉𝑐, 𝐺𝑒𝑟)    # nos * 3 2D array
    𝑎 = abc[:,1]   # 1D array of length nosf
    𝑏 = abc[:,2]
    𝑐 = abc[:,3]
    # nos * nosf 2D array; 1 hartree/Å³ = 4359.74417 GPa
    return @. (𝑎 * ( (𝑉𝑐[:,1]/𝑉𝑐)^(𝑏+1) - 1 ) + 𝑐) * 4359.74417
end

function calc_Δ𝐺𝑡𝑜𝑡(𝐺𝑡𝑜𝑡, mol = molecularity)
    Δ𝐺𝑡𝑜𝑡 = Array{Float64}(undef, size(𝐺𝑡𝑜𝑡))  # nos * nosf 2D array
    if mol == "uni"
        for i in 1:length(𝐺𝑡𝑜𝑡[1,:])
            @. Δ𝐺𝑡𝑜𝑡[:,i] = (𝐺𝑡𝑜𝑡[:,i] - 𝐺𝑡𝑜𝑡[1,i]) * 627.509  # 1 hartree = 627.509 kcal/mol
        end
    elseif mol == "bi"
        Δ𝐺𝑡𝑜𝑡[1,:] .= Δ𝐺𝑡𝑜𝑡[2,:] .= 0.0
        for i in 1:length(𝐺𝑡𝑜𝑡[1,:])
            @. Δ𝐺𝑡𝑜𝑡[3:end,i] = (𝐺𝑡𝑜𝑡[3:end,i] - 𝐺𝑡𝑜𝑡[1,i] - 𝐺𝑡𝑜𝑡[2,i]) * 627.509
        end
    end
    return Δ𝐺𝑡𝑜𝑡
end

#= need to figure out how to calculate barrier; use the TS structure or the maximum?
function calculateΔ𝑉activation()
    𝑝̄ = average𝑝()    # 1 * nosf 2D array
    Δ𝐺𝑡𝑜𝑡 = calculateΔ𝐺𝑡𝑜𝑡()   # nos * nosf 2D array
    Δ𝐺𝑡𝑜𝑡activation = Δ𝐺𝑡𝑜𝑡[50,:]    # 1D array of nosf
    slope, intercept = [ones(length(𝑝̄)) vec(𝑝̄)] \ Δ𝐺𝑡𝑜𝑡activation
    return slope * 4.184    # 1 kcal mol⁻¹ / GPa = 4.184 cm³/mol; 4.184 * 10^3 / 10^9 * 10^6
end
=#

#------------------------------------------------------------------------------
# main.jl
#------------------------------------------------------------------------------
#function main()
    # Step 1: cavity volume 𝑉𝑐(𝑓) and solvent property calculations
    if restart  # restart "Vc" jobs
        restart_jobs("Vc")
    else  #new jobs
        mkpath("tmp")    # creat a tmp folder in current directory
        write_gjf("Vc")
        run_gaussian("Vc")
    end

    𝑉𝑐 = get_data("Vc", "Cavity volume", 5)
    # linear scaling factor 𝑠, as the cubic root of the volume scaling
    𝑠 = @. ∛(𝑉𝑐/𝑉𝑐[:,1])  # nos * nosf 2D array
    # average of 𝑠 over all structures at the same scalingfactor 𝑓
    𝑠̄ = mean(𝑠, dims=1)  # 1D array of length nosf
    # dielectric permitivity 𝜀 = 1 + (𝜀₀-1)/𝑠̄³
    𝜀₀ = get_sol_params().𝜀
    global 𝜀 = @. 1 + (𝜀₀ - 1) / 𝑠̄^3  # 1D array of length nosf
    # solvent density 𝜌 = 𝜌₀/𝑠̄⁽³⁺𝜂⁾
    𝜌₀ = get_sol_params().𝜌
    global 𝜌 = @. 𝜌₀ / 𝑠̄^(3+𝜂)  # 1D array of length nosf
    # molar volume of solvent 𝑉ₘ = (𝑀/𝜌₀) * 𝑠̄³
    𝑀 = get_sol_params().𝑀
    global 𝑉ₘ = @. (𝑀/𝜌₀) * 𝑠̄^3  # 1D array of length nosf

    # Step 2: electronic structure Gaussian jobs and pressure calculations
    if restart  # restart Ger jobs
        restart_jobs("Ger")
    else
        write_gjf("Ger")
        run_gaussian("Ger")
    end
    𝐺𝑒𝑟 = get_data("Ger", "After PCM corrections", 7)
    𝑝 = calc_𝑝(𝑉𝑐, 𝐺𝑒𝑟)
    # average of 𝑝 over all structures at the same scalingfactor 𝑓
    𝑝̄ = mean(𝑝, dims=1)   # 1 * nosf 2D array

    # Step 3: cavitation energy Gaussian jobs
    if restart  # restart Gcav jobs
        restart_jobs("Gcav")
    else
        write_gjf("Gcav")
        run_gaussian("Gcav")
    end
    𝑉𝑐𝑎𝑣 = get_data("Gcav", "Cavity volume", 5) # 𝑉𝑐𝑎𝑣 could be different from 𝑉𝑐 in Gcav calculation using hard sphere
    𝐺𝑛𝑜𝑛_𝑝𝑉 = get_data("Gcav", "PCM non-electrostatic energy", 5)
    #𝐺𝑛𝑜𝑛_𝑝𝑉 and 𝑉𝑐𝑎𝑣 are nos * nosf 2D arrays; 1 GPa*Å³ = 2.293712569e-4 Hartree
    𝐺𝑐𝑎𝑣 = @. 𝐺𝑛𝑜𝑛_𝑝𝑉 + 𝑝̄ * 𝑉𝑐𝑎𝑣 * 2.293712569e-4    # nos * nosf 2D array
    𝐺𝑡𝑜𝑡 = 𝐺𝑒𝑟 .+ 𝐺𝑐𝑎𝑣    # nos * nosf 2D array
    Δ𝐺𝑡𝑜𝑡 = calc_Δ𝐺𝑡𝑜𝑡(𝐺𝑡𝑜𝑡)

    # print results to properties.dat file
    #𝑓 = scalingfactors
    write_properties(scalingfactors, 𝑉𝑐, 𝑠, 𝑠̄, 𝜀, 𝜌, 𝑉ₘ, 𝐺𝑒𝑟, 𝑝, 𝑝̄, 𝑉𝑐𝑎𝑣, 𝐺𝑛𝑜𝑛_𝑝𝑉, 𝐺𝑐𝑎𝑣, 𝐺𝑡𝑜𝑡, Δ𝐺𝑡𝑜𝑡)
#end

#main()
