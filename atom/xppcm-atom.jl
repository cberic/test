#using Statistics
using Printf
using LsqFit
#include("Li.jl")
#filename = "Li"
include(ARGS[1])
filename = replace(ARGS[1], ".jl" => "")  # remove the ".jl" extension

#------------------------------------------------------------------------------
# Solvents.jl
#------------------------------------------------------------------------------
struct Solvent <: Real #FieldVector{5, Real} #
	𝜀::Float64  # dielectric constant
	𝜌::Float64  # effective solvent density
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


#------------------------------------------------------------------------------
# Gaussian.jl
#------------------------------------------------------------------------------
# Gaussian input sections are explained here: https://gaussian.com/input/?tabid=0
# Link0, Route, Title, Molecule Specification, etc. sections
function print_content(io::IO, i_𝑓::Int64, np::Int64 = nproc, mem::String = mem, kws::String = keywords, sol::String = solvent, 𝑓list = scalingfactors, chrg::Int64 = charge, mulplct::Int64 = multiplicity, geom::String = geometry, tsare::Float64 = tesserae, sp = get_sol_params(), 𝜀 = 𝜀, 𝜌 = 𝜌, 𝑟₀ = 𝑟₀)
    # link0
    println(io, "%chk=$filename-Ger.chk")
    println(io, "%nproc=", np)
    println(io, "%mem=", mem)
    # route
    println(io, "#p ", kws, i_𝑓 > 1 ? " guess=read" : "")
    println(io, "#p scrf=(iefpcm,solvent=", sol, ",read) nosymm 6d 10f")
    println(io)
    #title
    println(io, "Ger calculation with scalingfactor = ", 𝑓list[i_𝑓])
    println(io)
    # molecule specification
    println(io, chrg," ",mulplct)
    println(io, strip(geom))
    println(io)
    # PCM keywords specification
    println(io, "qrep pcmdoc geomview nodis nocav g03defaults tsare=",tsare)
    println(io, "nvesolv=",sp.𝑛," solvmw=",sp.𝑀)
    println(io, "eps=",𝜀[i_𝑓]," rhos=",𝜌[i_𝑓])
    println(io, "nsfe=1 noaddsph")
    println(io)
    # PCM sphere specification
    println(io, "1    ", 𝑟₀, "    ", 𝑓list[i_𝑓])
    println(io)
end

function write_gjf(𝑓list = scalingfactors)
    open("$filename-Ger.gjf", "w") do file
        for i_𝑓 in eachindex(𝑓list)
            print_content(file, i_𝑓)
            if i_𝑓 != lastindex(𝑓list) println(file, "--link1--") end
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

function run_gaussian()
    gau = get_gau_ver()
    run(`$gau $filename-Ger.gjf`)
end

# extract data from Gaussian .log files
function get_data(searchstring::String, fieldnum::Int64, 𝑓list = scalingfactors)
    data = Vector{Float64}(undef, lastindex(𝑓list))  # 1D array
    i_𝑓 = 1    # i_𝑓 ranges from 1 to lastindex(𝑓list)
    open("$filename-Ger.log", "r") do file
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
# write_properties.jl
#------------------------------------------------------------------------------
function write_properties(𝑓list, 
                        𝑉𝑐::Vector{Float64}, 
                        𝑠::Vector{Float64}, 
                        𝜀::Vector{Float64}, 
                        𝜌::Vector{Float64}, 
                        𝐺𝑒𝑟::Vector{Float64}, 
                        𝑝::Vector{Float64})
    𝐸𝑔𝑎𝑠 = get_data("<psi(f)|   H    |psi(f)>", 6) # gas phase energy with the solvation field wave function
    𝐸ₑₗₑₛₜₐₜ = get_data("(Polarized solute)-Solvent", 5)
    𝐸ₚₐᵤₗᵢ = get_data("Quantum repulsion energy", 6)
    open("$filename-properties.dat", "w") do file
        println(file, "#      𝑓       𝑉𝑐       𝑠       𝜀     𝜌ₛₒₗ           𝐸𝑔𝑎𝑠  𝐸ₑₗₑₛₜₐₜ    𝐸ₚₐᵤₗᵢ             𝐺𝑒𝑟        𝑝")
        println(file, "               Å³                                      Eₕ  kcal/mol  kcal/mol              Eₕ      GPa")
        for j in eachindex(𝑓list)
            @printf(file, "%-2d %5.3f  %7.3f  %6.4f  %6.4f  %7.4f %14.8f    %6.2f    %6.2f  %14.8f  %7.3f\n", j, 𝑓list[j], 𝑉𝑐[j], 𝑠[j], 𝜀[j], 𝜌[j], 𝐸𝑔𝑎𝑠[j], 𝐸ₑₗₑₛₜₐₜ[j], 𝐸ₚₐᵤₗᵢ[j], 𝐺𝑒𝑟[j], 𝑝[j])
        end
        #println(file)
    end
end

#------------------------------------------------------------------------------
# Murnaghan equation of state fitting for pressure 𝑝 calculation
#------------------------------------------------------------------------------
# using LsqFit
function eos_fitting(𝑉𝑐, 𝐺𝑒𝑟)
    abc_parameters = Vector{Float64}(undef, 3)
    # python: y = (a/b)*(1/x)**b+(a-c)*x; y is Ger-Ger(Vc_0), x is Vc/Vc_0
    # mathematica: a*x ((1/b)*(t[[1, 1]]/x)^(b + 1) + 1) - c*x
    # LsqFit: a=p[1], b=p[2], c=p[3], x is Vc/Vc_0
    @. model(x, p) = (p[1]/p[2])*x^(-p[2]) + (p[1]-p[3])*x - p[1] - p[1]/p[2] + p[3]
    xdata = 𝑉𝑐 / 𝑉𝑐[1]
    ydata = 𝐺𝑒𝑟 .- 𝐺𝑒𝑟[1]
    p0 = [0.1, 5.0, 0.1]
    fit = curve_fit(model, xdata, ydata, p0)
    abc_parameters[1] = fit.param[1]/𝑉𝑐[1]
    abc_parameters[2] = fit.param[2]
    abc_parameters[3] = fit.param[3]/𝑉𝑐[1]
    return abc_parameters    # 1D array
end

function calc_𝑝(𝑉𝑐, 𝐺𝑒𝑟)
    abc = eos_fitting(𝑉𝑐, 𝐺𝑒𝑟)  # 1D array
    𝑎 = abc[1]
    𝑏 = abc[2]
    𝑐 = abc[3]
    # 1D array; 1 hartree/Å³ = 4359.74417 GPa
    return @. (𝑎 * ( (𝑉𝑐[1]/𝑉𝑐)^(𝑏+1) - 1 ) + 𝑐) * 4359.74417
end

#------------------------------------------------------------------------------
# main
#------------------------------------------------------------------------------
#function main()
    # Step 1: cavity volume 𝑉𝑐(𝑓) and solvent property calculations
    # Because we are dealing with single atoms in this project, the 
    # cavity volume can be computed directly without asking Gaussian
    # to do it.
    𝑟₀ = get_atom_radius(split(geometry)[1])
    𝑓 = scalingfactors  # 1D array
    𝑉𝑐 = 4/3 * pi * (𝑓 * 𝑟₀).^3  # 1D array
    # linear scaling factor 𝑠, as the cubic root of the volume scaling
    # In the case of a single atom, 𝑠 has a simple relation to 𝑓
    𝑠 = 𝑓/𝑓[1]  # 1D array
    # dielectric permitivity 𝜀 = 1 + (𝜀₀-1)/𝑠³
    𝜀₀ = get_sol_params().𝜀
    𝜀 = @. 1 + (𝜀₀ - 1) / 𝑠^3  # 1D array
    # effective solvent density 𝜌 = 𝜌₀/𝑠⁽³⁺𝜂⁾
    𝜌₀ = get_sol_params().𝜌
    𝜌 = @. 𝜌₀ / 𝑠^(3+𝜂)  # 1D array
    # molar volume of solvent 𝑉ₘ = (𝑀/𝜌₀) * 𝑠³
    #𝑀 = get_sol_params().𝑀
    #global 𝑉ₘ = @. (𝑀/𝜌₀) * 𝑠^3  # 1D array

    # Step 2: electronic structure Gaussian jobs and pressure calculations
    write_gjf()
    run_gaussian()
    𝐺𝑒𝑟 = get_data("After PCM corrections", 7)
    𝑝 = calc_𝑝(𝑉𝑐, 𝐺𝑒𝑟)  # 1D array

    # print results to properties.dat file
    write_properties(𝑓, 𝑉𝑐, 𝑠, 𝜀, 𝜌, 𝐺𝑒𝑟, 𝑝)
#end

#main()
