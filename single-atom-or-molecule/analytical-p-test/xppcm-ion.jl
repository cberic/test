using Printf
using LsqFit
using DelimitedFiles

#include("Cl.jl")
#filename_without_extension = "Cl"
include(ARGS[1])
filename_without_extension = replace(ARGS[1], ".jl" => "")  # remove the ".jl" extension

#------------------------------------------------------------------------------
# solvent
#------------------------------------------------------------------------------
# returns a 5-element tuple of (dielectric constant 𝜀, solvent density 𝜌,
# molar mass 𝑀, number of valence electrons, molecular radius (Ang.) )
function solventparameters(s = solvent)
    if s == "cyclohexane"
        # if dielectric is defined in input.jl, use it, otherwise use default
        𝜀 = @isdefined(dielectric) ? dielectric : 2.0165
        return (𝜀, 0.7781, 84.1595, 36, 2.815)
    elseif s == "benzene"
        𝜀 = @isdefined(dielectric) ? dielectric : 2.2706
        return (𝜀, 0.8756, 78.1118, 30, 2.63)
    elseif s == "argon"
        𝜀 = @isdefined(dielectric) ? dielectric : 1.43
        return (𝜀, 1.3954, 39.948, 8, 1.705)
    elseif s == "xenon" # CsCl compression chamber
        𝜀 = @isdefined(dielectric) ? dielectric : 1.6403^2
        return (𝜀, 3.99, 168.36, 16, 2.49)
    elseif s == "krypton" # NaCl compression chamber
        𝜀 = @isdefined(dielectric) ? dielectric : 1.544^2
        return (𝜀, 2.16, 58.44, 16, 2.49)
    else 
        error("solvent not implemented. Try cyclohexane, benzene, or argon")
    end
end

#------------------------------------------------------------------------------
# atomicradii
#------------------------------------------------------------------------------
function atomicradii(type = radiustype)
    if type == "bondi"
        return Dict(
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
            "Rn"=> 2.20,    "86"=> 2.20)
        # add more if needed from https://en.wikipedia.org/wiki/Van_der_Waals_radius
    elseif type == "rahm"
        return Dict(
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
            "Rn"=> 2.43,    "86"=> 2.43)
        # add more if needed from https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/chem.201700610
    elseif type == "ionic"
        return Dict(
		    #1s

            #2s

            #2p
			"F" => 1.92,    "9" => 1.92, # From RA16
            #3s
            #"Na"=> 0.95,    "11"=> 0.95, # From collin2005, hard ionic radii from Pauling
			#"Na"=> 1.20,    "11"=> 1.20, # From Alvaro
			"Na"=> 1.33,    "11"=> 1.33, # From RA16
            #3p
			"Cl"=> 2.29,    "17"=> 2.29, # From RA16
			#"Cl"=> 2.001,    "17"=> 2.001, # From Alvaro
            #4s
            #"K" => 2.34,    "19"=> 2.34,
            #"Ca"=> 2.70,    "20"=> 2.70,
            #3d

            #4p
			"Br"=> 2.41,    "35"=> 2.41, #RA16
            #5s

			#6s

            #4d

            #5p
			"I" => 2.59,    "53"=> 2.59, #From RA16
            #6p
            #"Tl"=> 2.42,    "81"=> 2.42, Keeping these as examples for now
            #"Pb"=> 2.49,    "82"=> 2.49,
            #"Bi"=> 2.50,    "83"=> 2.50,
            #"Po"=> 2.50,    "84"=> 2.50,
            #"At"=> 2.47,    "85"=> 2.47,
            #"Rn"=> 2.43,    "86"=> 2.43
            )
        # add more if needed from https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/chem.201700610
    else
        error("radiustype not supported. Try bondi or rahm")
    end
end

#------------------------------------------------------------------------------
# geometries
#------------------------------------------------------------------------------
function tidygeometries(geom = geometries)
    # remove leading and trailing spaces/blank lines 
    # and add a space to the beginning
    # "*" is the string concatenation operator
    a = " " * strip(geom)
    
    # remove leading and trailing spaces of each coordinate line 
    # and add a space to the begininig of coordinate line
    return replace(a, r"\s*\n\s*" => "\n ")
end


function structurelines(geom = geometries)
    g = tidygeometries(geom)
    return split(g, "\n", keepempty=false)
end


function numberofatoms(geom = geometries)
    return length(structurelines(geom))
end


function atomlist(geom = geometries)
    lines = structurelines(geom)
    noa = numberofatoms(geom)
    atoms = Array{String}(undef, noa)
    for i in 1:noa
        atoms[i] = split(lines[i], limit=2)[1]
    end
    return atoms
end


function coordinatelines(geom = geometries)
    lines = structurelines(geom)
    noa = numberofatoms(geom)
    coordlines = Array{String}(undef, noa)
    for i in 1:noa
        # get_ the 2nd-last fields of each coordinate line & add a leading space
        coordlines[i] = " " * split(lines[i], limit=2)[2]
    end
    return coordlines
end


#------------------------------------------------------------------------------
# charge sphere
#------------------------------------------------------------------------------

function charge_sphere(t = t_new, r = r_new)
    #Defining neccessary arrays that will be used in the coming calculations
    phi = zeros(0)
    theta = zeros(0)
    x = zeros(0)
    y = zeros(0)
    z = zeros(0)

    #Specified values
    num_pts = 1000
    indices = 0.5:999.5 #Temporary solution until I understand exactly how to handle these arrays in Julia
    if model =="basic"; q = 0; end
    if model =="pointcharges"; q= (-1)*charge; end
    point_charge = q/num_pts
    point_charges = fill(point_charge, num_pts)
    point_charges = string.(point_charges)

    #Determine the angles used to create the sphere based on the number of points and indices
    for k in 1:length(indices)
        phi_new = acos(1.0 - 2.0*indices[k]/num_pts)
        theta_new = pi * (1.0 + 5.0^(0.5)) * indices[k]
        append!(phi, phi_new)
        append!(theta, theta_new)
    end

    #calculate x, y, z positions for the charge sphere
    for k in 1:length(indices)
        x_new = t * r * cos(theta[k]) * sin(phi[k])
        y_new = t * r * sin(theta[k]) * sin(phi[k])
        z_new = t * r * cos(phi[k])
        append!(x, x_new)
        append!(y, y_new)
        append!(z, z_new)
    end

    #This is an inelegant solution but what is done in the final parts of the function is that we turn all values
    #into strings and then concate each one of them in the proper form that Gaussian wants
    x_string = string.(x)
    y_string = string.(y)
    z_string = string.(z)

    sphere = ""
    old_string = ""
    
    for m in 1:length(x_string)
        sphere = old_string * " " * x_string[m] * " " * y_string[m] * " " * z_string[m] * " " * point_charges[m] * "\n"
        old_string = sphere
    end
    return sphere
end


#------------------------------------------------------------------------------
# get electric field modulus and Potential energy of the charges-electron interaction
#------------------------------------------------------------------------------

function get_EFM(𝑓 = scalingfactors)
    #some initial parameters that is needed in the script
    if model =="basic"; totalcharge = 0; end
    if model =="pointcharges"; totalcharge= (-1)*charge; end
    N = 1000
    pq = totalcharge/N
    L = length(𝑓)
    EFM = Array{Float64}(undef, L, N)
    EPot = Array{Float64}(undef, L, N)
    sum_EFM = Array{Float64}(undef, L)
    sum_EPot = Array{Float64}(undef, L)
    j = 1 #j is used for iterative processes
    fprimc = readlines("Ger.log") #read Gaussian output file
    t = length(fprimc)
    p = 1

    for i in 1:t
        # first, identify where the electrostatic terms from prop keyword are in the output
        if occursin("Electrostatic Properties (Atomic Units)", fprimc[i])
            ion = i+6 # line which is the atom on, then iterate over number of charge points after this line
            Ec = Array{Float64}(undef, N)
            Potc = Array{Float64}(undef, N)
            for m in ion+1:ion+N
                a = parse(Float64, split(fprimc[m])[3])
                b = parse(Float64, split(fprimc[m])[4])
                c = parse(Float64, split(fprimc[m])[5])
                d = (a*a + b*b + c*c)^(0.5) #modulus of the electric field
                e = parse(Float64, split(fprimc[m])[2]) #potential
                Ec[j] = d
                Potc[j] = e
                j += 1
            end
            for q in 1:N
                EFM[p,q] = (Ec[q])*pq
                EPot[p,q] = (Potc[q])*pq
            end
            # resetting j for iterative purposes
            j = 1
            p += 1
        end
    end

    # calculating the wanted terms of electric field modulus and potential
    for i in 1:size(EPot,1)
        sum_EFM[i] = sum(EFM[i,:])
        sum_EPot[i] = sum(EPot[i,:])
    end
    return sum_EFM, sum_EPot
end


#------------------------------------------------------------------------------
# io
#------------------------------------------------------------------------------
function writegjf(jobtype)
    if jobtype == "Vc"        # Cavitation volume jobs
        gjfvc()
        gjfvc2()
    elseif jobtype == "Ger"   # Electronic energy SCRF jobs
        gjfgeranalytical()
    end
end


function custombasis(file = "gen")
    if occursin(" gen", lowercase(keywords)) || occursin("/gen", lowercase(keywords))
        a = read(file, String)
        return replace(a, r"\s*\n+\s*\n+\s*" => "")
    end
end

# write a single gjf file for Vc jobs
function gjfvc(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    g = tidygeometries(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    coordlines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    gen = custombasis(gen_filename)
    
    for i in 1:a
        # writing mode for the first and appending mode for other 𝑓
        open("Vc.gjf", "$(i == 1 ? "w" : "a")") do file
            write(file, """
                %kjob l301
                %nproc=$nproc
                %mem=$mem
                #p $keywords
                # scrf=(iefpcm,solvent=$solvent,read) nosym guess=only pop=none
                
                scaling factor = $(𝑓[i])
                
                $charge $multiplicity
                $g
                $(gen===nothing ? "" : "\n$gen\n")
                qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae
                nsfe=$noa
                nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])
                eps=$(sp[1]) rhos=$(sp[2])
                
                """)
                
            for j in 1:noa
                write(file, " $(coordlines[j])    $(𝑟ₐ[atoms[j]])    $(𝑓[i])\n")
            end
            
            write(file, "\n")
            
            if i != a  # do not write --link1-- for the last scaling factor
                write(file, "--link1--\n")
            end
        end
    end
end


# separate gjf files for Vc jobs for each scaling factors
function gjfvc2(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    g = tidygeometries(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    coordlines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    gen = custombasis(gen_filename)
    
    for i in 1:a
        open("Vc-$(𝑓[i]).gjf", "w") do file
            write(file, """
                %kjob l301
                %nproc=$nproc
                %mem=$mem
                #p $keywords
                # scrf=(iefpcm,solvent=$solvent,read) nosym guess=only pop=none
                
                scaling factor = $(𝑓[i])
                
                $charge $multiplicity
                $g
                $(gen===nothing ? "" : "\n$gen\n")
                qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae
                nsfe=$noa
                nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])
                eps=$(sp[1]) rhos=$(sp[2])
                
                """)
                
            for j in 1:noa
                write(file, " $(coordlines[j])    $(𝑟ₐ[atoms[j]])    $(𝑓[i])\n")
            end
            write(file, "\n")
        end
    end
end


# write gjf files for Ger calculations for analytical p
function gjfgeranalytical(𝜌 = calc_𝜌(), geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    g = tidygeometries(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    coordlines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    𝜀 = calc_𝜀()    # data of 𝜀 and 𝜌 needed for the gjf files
    #𝜌 = calc_𝜌()
    𝑛𝑡𝑠 = get_numberoftesserae()
    gen = custombasis(gen_filename)
    r_new = 𝑟ₐ[atoms[1]] #Jonatan: Note that this is a temporary solution which will only work if you only use single ions.  

        # writing mode for the first and appending mode for other 𝑓
    open("Ger.gjf", "w") do file
        for j in 1:a
            if model == "pointcharges" || model == "basic"
                t_new = 𝑓[j]
                charge_sp = charge_sphere(t_new, r_new)
            end
            write(file, """
                $(j == 1 ? "" : "%oldchk=$(𝑓[j-1]).chk\n")%chk=$(𝑓[j]).chk
                %nproc=$nproc
                %mem=$mem
                #p $keywords $(j == 1 ? "" : "guess=read") $(model == "pointcharges" ? "charge" : "")
                # scrf=(iefpcm,solvent=$solvent,read) iop(5/33=1) iop(6/17=2) prop(efg,read,grid) nosym 6d 10f

                scaling factor = $(𝑓[j])

                $charge $multiplicity
                $g
                $(model == "pointcharges" ? "\n$charge_sp" : "")$(gen===nothing ? "" : "\n$gen\n")
                qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae
                nsfe=$noa
                nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])
                eps=$(𝜀[j]) rhos=$(𝜌[j])

                """)

            for k in 1:noa
                write(file, " $(coordlines[k])    $(𝑟ₐ[atoms[k]])    $(𝑓[j])\n")
            end
            write(file, "\n")
            write(file, "$charge_sp\n")
            write(file, "$(𝑛𝑡𝑠[j]), 1, $(round(Int, 𝑓[j]*10000)), $(round(Int, 𝑓[j]*10000+1))\n")
            write(file, "\n")
            # do not write "--link1--" for the last scaling factor
            if j != a
                write(file, "--link1--\n")
            end
        end
    end
end


#input for the Ger job of the first scaling factore
function gjfger_1st_scalingfactor(𝜌, geom = geometries, 𝑓 = scalingfactors)
    g = tidygeometries(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    coordlines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    𝜀 = calc_𝜀()    # data of 𝜀 and 𝜌 needed for the gjf files
    #𝜌 = calc_𝜌()
    𝑛𝑡𝑠 = get_numberoftesserae()
    gen = custombasis(gen_filename)
    r_new = 𝑟ₐ[atoms[1]] #Jonatan: Note that this is a temporary solution which will only work if you only use single ions.  

    open("Ger.gjf", "w") do file
        if model == "pointcharges"|| model == "basic"
            t_new = 𝑓[1]
            charge_sp = charge_sphere(t_new, r_new)
        end
        write(file, """
            %chk=$(𝑓[1]).chk
            %nproc=$nproc
            %mem=$mem
            #p $keywords $(model == "pointcharges" ? "charge" : "")
            # scrf=(iefpcm,solvent=$solvent,read) iop(5/33=1) iop(6/17=2) prop(efg,read,grid) nosym 6d 10f

            scaling factor = $(𝑓[1])

            $charge $multiplicity
            $g
            $(model == "pointcharges" ? "\n$charge_sp" : "")$(gen===nothing ? "" : "\n$gen\n")
            qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae
            nsfe=$noa
            nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])
            eps=$(𝜀[1]) rhos=$𝜌

            """)

        for k in 1:noa
            write(file, " $(coordlines[k])    $(𝑟ₐ[atoms[k]])    $(𝑓[1])\n")
        end
        write(file, "\n")
        write(file, "$charge_sp\n")
        write(file, "$(𝑛𝑡𝑠[1]), 1, $(round(Int, 𝑓[1]*10000)), $(round(Int, 𝑓[1]*10000+1))\n")
        write(file, "\n")
    end
end

# extract data from Gaussian .log write_properties files
function get_data(filename::String, searchstring::String, fieldnum::Int64, 𝑓 = scalingfactors)
    a = length(𝑓)
    #data = Array{Float64}(undef, a)
    data = zeros(a)
    j = 1    # j ranges from 1:length(𝑓)
    open(filename, "r") do file
        for line in eachline(file)
            if occursin(searchstring, line)
                data[j] = parse(Float64, split(line)[fieldnum])
                j += 1
            end
        end
    end
    data
end

# extract the total number of tesserae from tesserae-𝑓.off
function get_numberoftesserae(𝑓 = scalingfactors)
    a = length(𝑓)
    𝑛𝑡𝑠 = Array{Int64}(undef, a)
    for j in 1:a
        open("tesserae-$(𝑓[j]).off", "r") do file
            readline(file)
            secondline=readline(file)
            # the 2nd number on the 2nd line is the number of tesserae
            𝑛𝑡𝑠[j] = parse(Int64, split(secondline)[2])
        end
    end
    return 𝑛𝑡𝑠
end


# extract and calculate the xyz coordinates of the centers of the tesserae
# from tesserae-𝑓.off 
function writetesseragrid(geom = geometries, 𝑓 = scalingfactors)
    𝑛𝑡𝑠 = get_numberoftesserae()
    atoms = atomlist(geom)
    a = length(𝑓)
    𝑟ₐ = atomicradii()
    for j in 1:a
        tesseraecoordinates = Array{Float64}(undef, 𝑛𝑡𝑠[j],3)  # n*3 2D array
        x1 = x2 = x3 = y1 = y2 = y3 = z1 = z2 = z3 = 0.0
        linecount = 1
        i = 1   # i ranges from 1:n
        open("tesserae-$(𝑓[j]).off", "r") do file
            for line in eachline(file)
                if linecount >= 3 && linecount <= 3 * 𝑛𝑡𝑠[j] + 2
                    if linecount % 3 == 0
                        # xyz coordinates of the 1st vertex of a tessera
                        (x1, y1, z1) = [parse(Float64, s) for s in split(line)[1:3]]
                    elseif linecount % 3 == 1
                        # xyz coordinates of the 2nd vertex of a tessera
                        (x2, y2, z2) = [parse(Float64, s) for s in split(line)[1:3]]
                    elseif linecount % 3 == 2
                        # xyz coordinates of the 3rd vertex of a tessera
                        (x3, y3, z3) = [parse(Float64, s) for s in split(line)[1:3]]
                        # sum of three vectors, and normalized
                        (x, y, z) = (x1 + x2 + x3, y1 + y2 + y3, z1 + z2 + z3)
                        tesseraecoordinates[i,1] = x / sqrt(x^2 + y^2 + z^2)
                        tesseraecoordinates[i,2] = y / sqrt(x^2 + y^2 + z^2)
                        tesseraecoordinates[i,3] = z / sqrt(x^2 + y^2 + z^2)
                        i += 1  # i ranges from 1:n
                    end
                end
                linecount += 1 
            end
        end
        #writedlm("$(𝑓[j]).tsrcoord", tesseraecoordinates * 𝑟ₐ[atoms[1]] * 𝑓[j])
        writedlm("fort.$(round(Int, 𝑓[j]*10000))", tesseraecoordinates * 𝑟ₐ[atoms[1]] * 𝑓[j])
    end
end


# extract Electric Field Gradie𝑛𝑡𝑠 from fort.1201 etc files
function get_𝑒𝑓𝑔╱𝑛𝑡𝑠(𝑓 = scalingfactors)
    𝑛𝑡𝑠 = get_numberoftesserae()
    a = length(𝑓)
    𝑒𝑓𝑔sum = zeros(a)
    for j in 1:a
        linecount = 1
        open("fort.$(round(Int, 𝑓[j]*10000+1))", "r") do file
            for line in eachline(file)
                if linecount % 4 == 2
                    # zz component of the 𝑒𝑓𝑔
                    𝑒𝑓𝑔sum[j] -= parse(Float64, split(line)[4])
                elseif linecount % 4 == 3
                    # yy and xx compone𝑛𝑡𝑠 of the 𝑒𝑓𝑔
                    𝑒𝑓𝑔sum[j] -= sum([parse(Float64, s) for s in split(line)[1:2]])
                end
                linecount += 1
            end
        end
    end
    return @. 𝑒𝑓𝑔sum / 4 / pi / 𝑛𝑡𝑠
end

# extract electronic energy 𝐺ₑᵣ data from gaussian output files
function get_𝐸ₚₐᵤₗᵢ(𝑓 = scalingfactors)
    a = length(𝑓)
    𝐸ₚₐᵤₗᵢ = Array{Float64}(undef, a)
    j = 1    # j ranges from 1:length(𝑓)
    open("Ger.log", "r") do file
        for line in eachline(file)
            if occursin("QRepSI", line)
                𝐸ₚₐᵤₗᵢ[j] = parse(Float64, split(line)[6]) / 627.503 # 1 hartree = 627.503 kcal/mol
            elseif occursin("SCF Done", line)
                j += 1    # j ranges from 1:length(𝑓)
            end
        end
    end
    return 𝐸ₚₐᵤₗᵢ
end


function get_𝑊ₚₒₗ′(𝑓 = scalingfactors)
    a = length(𝑓)
    𝑊ₚₒₗ′ = Array{Float64}(undef, a)
    j = 1    # j ranges from 1:length(𝑓)
    open("Ger.log", "r") do file
        for line in eachline(file)
            if occursin("Polarized solute", line)
                array = split(line)
                if length(array) == 5
                    𝑊ₚₒₗ′[j] = parse(Float64, array[5]) / 627.503 # 1 hartree = 627.503 kcal/mol
                elseif length(array) == 4  # remove the equal sign in front of the number, e.g. array[4] == "=-1798336.01"
                    𝑊ₚₒₗ′[j] = parse(Float64, chop(array[4], head = 1)) / 627.503
                end
                j += 1    # j ranges from 1:length(𝑓)
            end
        end
    end
    return 𝑊ₚₒₗ′
end


function get_orbitalenergy(𝑓 = scalingfactors)
    a = length(𝑓)
    orbital = Array{String}(undef, a)
    j = 0
    open("Ger.log", "r") do file
        for line in eachline(file)
            if occursin("Population analysis", line)
                j += 1    # j ranges from 1:length(𝑓)
            end
            if occursin("eigenvalues", line)
                isassigned(orbital, j) ? orbital[j] *= (line * "\n") : orbital[j] = line * "\n"
            end
        end
    end
    return orbital
end


# use the Printf package to write the properties.dat file
function writeproperties(𝑉𝑐 = 𝑉𝑐, 𝐺ₑᵣ = 𝐺ₑᵣ, 𝑓 = scalingfactors)
    a = length(𝑓)
    𝑠 = calc_𝑠()
    𝜀 = calc_𝜀()
    𝜌 = calc_𝜌()
    𝑝n = calc_numerical𝑝()
    𝑝a = calc_analytical𝑝()
    Eorbital = get_orbitalenergy()
    open("$filename_without_extension-properties.dat", "w") do file
        write(file, "#     𝑓         𝑉𝑐(𝑓) Å³   𝑠(𝑓)        𝜀(𝑠)        𝜌(𝑠)        𝐺ₑᵣ(𝑓) a.u.    𝑝(𝑓)-numeric. -analyt.GPa\n")
        for j in 1:a
            @printf(file, "%-2d    %.3f     %7.3f    %.6f    %.6f    %9.6f    %.8f    %6.3f    %6.3f\n", 
                            j,   𝑓[j],   𝑉𝑐[j],   𝑠[j],    𝜀[j],   𝜌[j],   𝐺ₑᵣ[j],  𝑝n[j],  𝑝a[j])
        end
        write(file, "\n")
        for j in 1:a
            @printf(file, "𝑓 = %.3f    𝑝 = %6.3f GPa ----orbital energies in a.u.----\n", 𝑓[j], 𝑝a[j])
            write(file, Eorbital[j])
            write(file, "\n")
        end
    end
end


function writeproperties2(𝑉𝑐 = 𝑉𝑐, 𝑓 = scalingfactors)
    a = length(𝑓)
    #𝑠 = calc_𝑠()
    𝜀 = calc_𝜀()
    #𝜌 = calc_𝜌()
    #𝐺ₑᵣ = get_data("Ger.log", "SCF Done", 5)
    Eorbital = get_orbitalenergy()
    open("$filename_without_extension-properties.dat", "w") do file
        write(file, "#     𝑓         𝑉𝑐      𝑠       𝜀     𝜌ₛₒₗ        𝒵         𝑊ₑ       𝑊ₚₒₗ    𝐸ₚₐᵤₗᵢ           𝐺ₑᵣ     𝑊ₗ=𝑊ₑ+𝐺ₑᵣ        𝑝  𝑉_cell\n")
        write(file, "#               Å³                    g/ml   mol/ml         Eₕ         Eₕ        Eₕ            Eₕ            Eₕ      GPa      Å³\n")
        for j in 1:a
            @printf(file, "%-2d  %5.3f  %7.3f  %5.3f  %6.4f  %7.4f  %7.4f  %9.6f  %9.6f  %8.6f  %12.6f  %12.6f  %7.3f  %6.3f\n", 
                            j, 𝑓[j], 𝑉𝑐[j], 𝑠[j], 𝜀[j], 𝜌[j], 𝒵[j], 𝑊ₑ[j], 𝑊ₚₒₗ[j], 𝐸ₚₐᵤₗᵢ[j], 𝐺ₑᵣ[j], 𝑊ₗ[j], 𝑝[j], 𝑉_cell[j])
        end
        write(file, "\n")
        for j in 1:a
            @printf(file, "𝑓 = %.3f    𝑝 = %6.3f GPa ----orbital energies in a.u.----\n", 𝑓[j], 𝑝[j])
            write(file, Eorbital[j])
            write(file, "\n")
        end
    end
end


function writeproperties3(𝑉𝑐 = 𝑉𝑐, 𝑓 = scalingfactors)
    a = length(𝑓)
    #𝑠 = calc_𝑠()
    #𝜀 = calc_𝜀()
    #𝜌 = calc_𝜌()
    #𝐺ₑᵣ = get_data("Ger.log", "SCF Done", 5)
    Eorbital = get_orbitalenergy()
    open("$filename_without_extension-properties.dat", "w") do file
        write(file, "# atom=$(atomlist()[1]), charge=$charge, multiplicity=$multiplicity, radius=$(atomicradii()[atomlist()[1]]) Å\n")
        write(file, "# solvent=$solvent, dielectric=$dielectric, 𝜂=$𝜂, tesserae=$tesserae, xppcm-model=$model\n")
        write(file, "# $keywords; basis-set=$gen_filename\n\n")
        write(file, "#      𝑓         𝑉𝑐      𝑠       𝜀     𝜌ₛₒₗ        𝒵   𝐸(nu-ch)   𝑑𝐸(nu-ch)╱𝑑𝑠   𝐸(ch-ch)   𝑑𝐸(ch-ch)╱𝑑𝑠   𝐸(el-ch)   𝑑𝐸(el-ch)╱𝑑𝑠     𝐸ₚₒₗₐᵣ   𝑑𝐸ₚₒₗₐᵣ╱𝑑𝑠     𝐸ₚₐᵤₗᵢ   𝑑𝐸ₚₐᵤₗᵢ╱𝑑𝑠           𝐺ₑᵣ       𝑝ₐ       𝑝ₙ      𝑝ₙ′     𝑉/𝑉₀\n")
        write(file, "#                Å³                    g/ml                  Eₕ             Eₕ         Eₕ             Eₕ         Eₕ             Eₕ         Eₕ           Eₕ         Eₕ           Eₕ            Eₕ      GPa      GPa      GPa         \n")
        for j in 1:a
            @printf(file, "%-3d  %5.3f  %7.3f  %5.3f  %6.4f  %7.4f  %7.4f  %9.6f  %13.6f  %9.6f  %13.6f  %9.6f  %13.6f  %9.6f  %11.6f  %9.6f  %11.6f  %12.6f  %7.3f  %7.3f  %7.3f  %7.3f\n", 
                            j, 𝑓[j], 𝑉𝑐[j], 𝑠[j], 𝜀[j], 𝜌[j], 𝒵[j], 𝐸_nuclei_charges[j], 𝑑𝐸_nuclei_charges╱𝑑𝑠[j], 𝐸_charges_charges[j], 𝑑𝐸_charges_charges╱𝑑𝑠[j], 𝐸_electrons_charges[j], 𝑑𝐸_electrons_charges╱𝑑𝑠[j], 𝑊ₚₒₗ′[j], 𝑑𝑊ₚₒₗ′╱𝑑𝑠[j], 𝐸ₚₐᵤₗᵢ[j], 𝑑𝐸ᵣ╱𝑑𝑠[j], 𝐺ₑᵣ[j], 𝑝ₐ[j], 𝑝ₙ[j], 𝑝ₙ′[j], 𝑉𝑐[j]/𝑉𝑐[1])
        end
        write(file, "\n")
        write(file, "#      𝑓        𝑠       𝐸(nu-ch)      --- 𝑑𝐸(nu-ch)╱𝑑𝑠 ---       𝐸(ch-ch)      --- 𝑑𝐸(ch-ch)╱𝑑𝑠 ---       𝐸(el-ch)      --- 𝑑𝐸(el-ch)╱𝑑𝑠 ---      ------ 𝐸ₚₒₗₐᵣ ------      ---- 𝑑𝐸ₚₒₗₐᵣ╱𝑑𝑠 ----        𝐸ₚₐᵤₗᵢ      ---- 𝑑𝐸ₚₐᵤₗᵢ╱𝑑𝑠 ----               𝐺ₑᵣ       ----------- 𝑑𝐺ₑᵣ╱𝑑𝑠 ----------\n")
        write(file, "#                                         numer     analyt                         numer     analyt                         numer     analyt       Gaussian       Born          numer     analyt                        numer     analyt                            numer                      \n")
        write(file, "#                              A              B          C              D              E          F              G              H          I              J          K              L          M             N              O          P                 Q              R      CFILO      CFILP\n")
        for j in 1:a
            @printf(file, "%-3d  %5.3f  %5.3f      %9.6f      %9.6f  %9.6f      %9.6f      %9.6f  %9.6f      %9.6f      %9.6f  %9.6f      %9.6f  %9.6f      %9.6f  %9.6f     %9.6f      %9.6f  %9.6f      %12.6f      %9.6f  %9.6f  %9.6f\n", 
                            j, 𝑓[j], 𝑠[j], 𝐸_nuclei_charges[j], 𝑑𝐸_nuclei_charges╱𝑑𝑠ₙ[j], 𝑑𝐸_nuclei_charges╱𝑑𝑠[j], 𝐸_charges_charges[j], 𝑑𝐸_charges_charges╱𝑑𝑠ₙ[j], 𝑑𝐸_charges_charges╱𝑑𝑠[j], 𝐸_electrons_charges[j], 𝑑𝐸_electrons_charges╱𝑑𝑠ₙ[j], 𝑑𝐸_electrons_charges╱𝑑𝑠[j], 𝑊ₚₒₗ′[j], 𝑊ₚₒₗ[j], 𝑑𝑊ₚₒₗ′╱𝑑𝑠[j], 𝑑𝑊ₚₒₗ╱𝑑𝑠[j], 𝐸ₚₐᵤₗᵢ[j], 𝑑𝐸ₚₐᵤₗᵢ╱𝑑𝑠ₙ[j], 𝑑𝐸ᵣ╱𝑑𝑠[j], 𝐺ₑᵣ[j], 𝑑𝐺ₑᵣ╱𝑑𝑠ₙ[j], 
                            𝑑𝐸_nuclei_charges╱𝑑𝑠[j]+𝑑𝐸_charges_charges╱𝑑𝑠[j]+𝑑𝐸_electrons_charges╱𝑑𝑠[j]+𝑑𝑊ₚₒₗ′╱𝑑𝑠[j]+𝑑𝐸ₚₐᵤₗᵢ╱𝑑𝑠ₙ[j],
                            𝑑𝐸_nuclei_charges╱𝑑𝑠[j]+𝑑𝐸_charges_charges╱𝑑𝑠[j]+𝑑𝐸_electrons_charges╱𝑑𝑠[j]+𝑑𝑊ₚₒₗ′╱𝑑𝑠[j]+𝑑𝐸ᵣ╱𝑑𝑠[j])
        end
        write(file, "\n")
        for j in 1:a
            @printf(file, "𝑓 = %.3f    𝑝 = %6.3f GPa ----orbital energies in a.u.----\n", 𝑓[j], 𝑝ₐ[j])
            write(file, Eorbital[j])
            write(file, "\n")
        end
    end
end


function debug(𝑉𝑐 = 𝑉𝑐, 𝐺ₑᵣ = 𝐺ₑᵣ, 𝑓 = scalingfactors)
    a = length(𝑓)
    𝑠 = calc_𝑠()
    𝜀 = calc_𝜀()
    𝜌 = calc_𝜌()
    𝐸ₚₐᵤₗᵢ = get_𝐸ₚₐᵤₗᵢ()
    𝑒𝑓𝑔╱𝑛𝑡𝑠 = get_𝑒𝑓𝑔╱𝑛𝑡𝑠()
    𝒵 = calc_𝒵()
    #𝑝n = calc_numerical𝑝()
    𝑝a = calc_analytical𝑝()
    #Eorbital = get_orbitalenergy()
    open("$filename_without_extension-debug.dat", "w") do file
        write(file, "#    𝑓       𝑉𝑐(𝑓) Å³   𝑠(𝑓)         𝜀(𝑠)        𝜌ₛₒₗ(𝑠)     𝐺ₑᵣ(𝑓) a.u.     𝑝a(𝑓) GPa      PauliE(𝑓)     𝑒𝑓𝑔/𝑛𝑡𝑠     𝒵(𝑓)\n")
        for j in 1:a
            @printf(file, "%-2d    %.3f     %7.3f    %.6f    %.6f    %9.6f    %.8f    %6.3f    %9.6f    %9.6f    %9.6f\n", 
                            j,   𝑓[j],   𝑉𝑐[j],   𝑠[j],    𝜀[j],   𝜌[j],   𝐺ₑᵣ[j],  𝑝a[j],  𝐸ₚₐᵤₗᵢ[j], 𝑒𝑓𝑔╱𝑛𝑡𝑠[j], 𝒵[j])
        end
    end
end


function debug2(𝑉𝑐 = 𝑉𝑐, 𝑓 = scalingfactors)
    a = length(𝑓)
    #𝑠 = calc_𝑠()
    𝜀 = calc_𝜀()
    #𝜌 = calc_𝜌()
    open("$filename_without_extension-debug.dat", "w") do file
        write(file, "#     𝑓         𝑉𝑐(𝑓) Å³   𝑠(𝑓)        𝜀(𝑠)        𝜌ₛₒₗ(𝑠)    𝒵(𝑠)      𝑒𝑓𝑔╱𝑛𝑡𝑠(𝑠)   𝑊ₑ(𝑠)      𝑊ₚₒₗ(𝑠)      𝑊ₚₒₗ′(𝑠)      𝐸ₚₐᵤₗᵢ(𝑠)    𝑊ₗ(𝑠)        𝑝(𝑠) GPa\n")
        for j in 1:a
            @printf(file, "%-2d    %.3f     %7.3f    %.6f    %.6f    %.4f    %.4f    %.6f    %.6f    %.6f    %.6f    %.6f    %.6f    %6.3f\n", 
                            j,   𝑓[j],   𝑉𝑐[j],   𝑠[j],    𝜀[j],   𝜌[j],   𝒵[j], 𝑒𝑓𝑔╱𝑛𝑡𝑠[j],𝑊ₑ[j], 𝑊ₚₒₗ[j],𝑊ₚₒₗ′[j], 𝐸ₚₐᵤₗᵢ[j], 𝑊ₗ[j], 𝑝[j])
        end
    end
end

#------------------------------------------------------------------------------
# rungaussian
#------------------------------------------------------------------------------
# check if g16 or g09 is installed and loaded
function gaussianversion()
    if gethostname() == "atlas-fdr-login-01" || gethostname() == "atlas-fdr-login-02"
        return "g16"
    elseif typeof(Sys.which("g16")) === String
        return "g16"
    elseif typeof(Sys.which("g09")) === String
        return "g09"
    else
        error("g16 and g09 not found. Make sure g16 or g09 correctly calls Gaussian16 or Gaussian09")
    end
end


function rungaussian(jobtype)
    gau = gaussianversion()
    run(`$gau $jobtype.gjf`)
end


#------------------------------------------------------------------------------
# algebra
#------------------------------------------------------------------------------
# linear scaling factor 𝑠, as the cubic root of the volume scaling
function calc_𝑠(𝑉𝑐 = 𝑉𝑐)  # 𝑉𝑐 is a global variable
    # @. is a macro for vectoried operation/broadcasting
    return @. cbrt(𝑉𝑐/𝑉𝑐[1])
end


# dielectric permitivity 𝜀 = 1 + (𝜀₀-1)/𝑠̄³
function calc_𝜀()
    𝜀₀ = solventparameters()[1]
    𝑠 = calc_𝑠()
    return @. 1 + (𝜀₀ - 1) / 𝑠^3    # 1D array of length a
end


# Pauli repulson barrier 𝜌 = 𝜌₀/𝑠̄⁽³⁺𝜂⁾, where 𝜌₀ = 𝜌₀
function calc_𝜌(𝜂=𝜂)
    𝜌₀ = solventparameters()[2]
    𝑠 = calc_𝑠()
    return @. 𝜌₀ / 𝑠^(3+𝜂)    # 1D array of length a
end


# Murnaghan equation of state fitting for pressure 𝑝 calculation using LsqFit
function eosfitting(𝑉, 𝐸)
    abc_parameters = Array{Float64}(undef, 3)
    # y = (a/b)*(1/x)**b+(a-c)*x; y is Ger, x is Vc
    # a=p[1], b=p[2], c=p[3], x is Vc
    @. model(x, p) = (p[1]/p[2])*x^(-p[2]) + (p[1]-p[3])*x
    xdata = 𝑉 ./ 𝑉[1]
    ydata = 𝐸 .- 𝐸[1]
    p0 = [0.0, 5.0, 0.0]
    fit = curve_fit(model, xdata, ydata, p0)
    abc_parameters[1] = fit.param[1]/𝑉[1]
    abc_parameters[2] = fit.param[2]
    abc_parameters[3] = fit.param[3]/𝑉[1]
    return abc_parameters
end


function calc_numerical𝑝(𝑉 = 𝑉𝑐, 𝐸 = 𝐺ₑᵣ)
    𝑎, 𝑏, 𝑐 = eosfitting(𝑉, 𝐸)
    # 1 hartree/Å³ = 4359.74417 GPa
    return @. (𝑎 * ( (𝑉𝑐[1]/𝑉𝑐)^(𝑏+1) - 1 ) + 𝑐) * 4359.74417 
end


# eq (9) in DOI:10.1002/jcc.25544
# In this eq, 𝒵 is called Alpha in Gaussian
function calc_𝒵()
    sp = solventparameters()
    # 𝒵 = fA*RhoS*NVES/SolvMW  ; fA=0.063d0
    # 𝒵0 = 0.063 * sp[2] * sp[4] / sp[3]
    𝜌 = calc_𝜌()
    return @. 0.063 * 𝜌 * sp[4] / sp[3]
end

function calc_𝒵_new(𝒵, 𝑅𝑟𝑒𝑓, 𝑓=[scalingfactors[1]])
    #sp = solventparameters()
    #𝜌 = calc_𝜌(𝜂)
    𝐸ₚₐᵤₗᵢ = get_𝐸ₚₐᵤₗᵢ(𝑓)[1]
    𝑛𝑡𝑠 = get_numberoftesserae(𝑓)[1]
    𝑒𝑓𝑔╱𝑛𝑡𝑠 = get_𝑒𝑓𝑔╱𝑛𝑡𝑠(𝑓)[1]
    𝑒𝑓𝑔 = 𝑛𝑡𝑠 * 𝑒𝑓𝑔╱𝑛𝑡𝑠
    𝐼₁ = 𝐸ₚₐᵤₗᵢ / 𝒵
    𝐼₂ = 4π * 𝑅𝑟𝑒𝑓^3 * 𝑒𝑓𝑔╱𝑛𝑡𝑠 #-𝑅𝑟𝑒𝑓 * (4π * 𝑅𝑟𝑒𝑓^2 / 𝑛𝑡𝑠) * 𝑒𝑓𝑔
    denominator = (3 + 𝜂) * 𝐼₁ + 𝐼₂
    numerator =  𝛼ᵣ * abs(charge)^2 / 𝑟₀ + 0.5(1 - 1/dielectric) * abs(charge)^2 / 𝑅𝑟𝑒𝑓 * (1 + 3/dielectric)
    𝒵_new =  numerator / denominator

    open("iterativeZ.dat", "a") do file
        println(file, #"𝜌_sol ", 𝜌, 
            " 𝒵 ", 𝒵, 
            " 𝐸ₚₐᵤₗᵢ ", 𝐸ₚₐᵤₗᵢ, 
            " 𝑛𝑡𝑠 ", 𝑛𝑡𝑠, 
            " 𝑒𝑓𝑔/𝑛𝑡𝑠 ", 𝑒𝑓𝑔╱𝑛𝑡𝑠,
            " 𝑒𝑓𝑔 ", 𝑒𝑓𝑔, 
            " 𝐼₁ ", 𝐼₁, 
            " 𝐼₂ ", 𝐼₂, 
            " numerator ", numerator,
            " denominator ", denominator,
            " 𝒵_new ", 𝒵_new)
    end

    return 𝒵_new
end

function calc_𝒵_new_pointcharges(𝒵, 𝑅𝑟𝑒𝑓, 𝑓=[scalingfactors[1]])
    #sp = solventparameters()
    #𝜌 = calc_𝜌(𝜂)
    𝐸_nuclei_charges = get_data("Ger.log", "Nuclei-charges interaction", 4, 𝑓)[1]
    𝐸_charges_charges = get_data("Ger.log", "Self energy", 7, 𝑓)[1]
    𝐸ₚₐᵤₗᵢ = get_𝐸ₚₐᵤₗᵢ(𝑓)[1]
    𝑛𝑡𝑠 = get_numberoftesserae(𝑓)[1]
    𝑒𝑓𝑔╱𝑛𝑡𝑠 = get_𝑒𝑓𝑔╱𝑛𝑡𝑠(𝑓)[1]
    𝑒𝑓𝑔 = 𝑛𝑡𝑠 * 𝑒𝑓𝑔╱𝑛𝑡𝑠
    𝐼₁ = 𝐸ₚₐᵤₗᵢ / 𝒵
    𝐼₂ = 4π * 𝑅𝑟𝑒𝑓^3 * 𝑒𝑓𝑔╱𝑛𝑡𝑠 #-𝑅𝑟𝑒𝑓 * (4π * 𝑅𝑟𝑒𝑓^2 / 𝑛𝑡𝑠) * 𝑒𝑓𝑔
    denominator = (3 + 𝜂) * 𝐼₁ + 𝐼₂
    numerator =  0.5(1 - 1/dielectric) * abs(charge)^2 / 𝑅𝑟𝑒𝑓 * (1 + 3/dielectric) -𝐸_nuclei_charges - 𝐸_charges_charges + 𝑅𝑟𝑒𝑓 * get_EFM()[1][1]
    𝒵_new =  numerator / denominator

    open("iterativeZ.dat", "a") do file
        println(file, 
            "𝜌_sol ", 𝒵 * sp[3] / sp[4] / 0.063, 
            " 𝒵 ", 𝒵, 
            " 𝐸_nuclei_charges ", 𝐸_nuclei_charges, 
            " 𝐸_charges_charges ", 𝐸_charges_charges, 
            " 𝐸_electrons_charges ", get_EFM()[2][1],
            " 𝐸ₚₐᵤₗᵢ ", 𝐸ₚₐᵤₗᵢ, 
            " 𝑛𝑡𝑠 ", 𝑛𝑡𝑠, 
            " 𝑒𝑓𝑔/𝑛𝑡𝑠 ", 𝑒𝑓𝑔╱𝑛𝑡𝑠,
            " 𝑒𝑓𝑔 ", 𝑒𝑓𝑔, 
            " 𝐼₁ ", 𝐼₁, 
            " 𝐼₂ ", 𝐼₂, 
            " numerator ", numerator,
            " denominator ", denominator,
            " 𝒵_new ", 𝒵_new)
    end

    return 𝒵_new
end


# eq (24) in DOI:10.1002/jcc.25544
# for atoms only; lattice energy (Born model) or ion-charges energy (point charges model) not considered
function calc_analytical𝑝(𝜂 = 𝜂, 𝑉𝑐 = 𝑉𝑐)
    𝐸ₚₐᵤₗᵢ = get_𝐸ₚₐᵤₗᵢ()
    𝒵 = calc_𝒵()
    𝑒𝑓𝑔╱𝑛𝑡𝑠 = get_𝑒𝑓𝑔╱𝑛𝑡𝑠()
    # 1 angstrom = 1.88973 bohr; 1 hartree/Å³ = 4359.74417 GPa
    return @. ((3 + 𝜂)/(3 * 𝑉𝑐 * 1.88973^3) * 𝐸ₚₐᵤₗᵢ + 𝒵 * 𝑒𝑓𝑔╱𝑛𝑡𝑠) * 1.88973^3 * 4359.74417
end

#------------------------------------------------------------------------------
# main
#------------------------------------------------------------------------------
#function main(𝑓 = scalingfactors)
    # Step 1: cavity volume 𝑉𝑐(𝑓) Gaussian jobs and solvent property calculations
    writegjf("Vc")
    𝑓 = scalingfactors
    for j in 𝑓
        rungaussian("Vc-$j")
        # combine the output files of Vc jobs at different f into one file
        open("Vc.log", "$(j == first(𝑓) ? "w" : "a")") do file
            write(file, read("Vc-$j.log", String))
        end
        run(`cp tesserae.off tesserae-$j.off`)
    end

    global 𝑉𝑐 = get_data("Vc.log", "GePol: Cavity volume", 5)

    # Step 2: electronic structure Gaussian jobs
    writetesseragrid()

    # if radiustype !== "ionic"
    #     writegjf("Ger")
    #     rungaussian("Ger")
    #     global 𝐺ₑᵣ = get_data("Ger.log", "SCF Done", 5)
    #     writeproperties()
    #     #debug()
    # end

    # if radiustype == "ionic"
    #     if model == "Born"
    #         if impose_equilibrium == true # at the 1st scalingfactor so that p(f0)=0
    #         end
    #         if impose_equilibrium == false
    #         end
    #         # self-consistent calculation of 𝒵
    #         sp = solventparameters()
    #         ##RC-101221!
    #         𝜌_guess = 5.0
    #         𝒵_guess = 0.063 * 𝜌_guess * sp[4] / sp[3]
    #         gjfger_1st_scalingfactor(𝜌_guess)
    #         rungaussian("Ger")
    #         open("iterativeZ.dat", "w") do file end
    #         𝑟ₐ = atomicradii()
    #         atoms = atomlist()
    #         𝑅𝑟𝑒𝑓 = 𝑓[1] * 𝑟ₐ[atoms[1]] * 1.88973 # reference radius of Cl- in bohr
    #         𝒵_new = calc_𝒵_new(𝒵_guess, 𝑅𝑟𝑒𝑓)
    #         while !(0.999 < 𝒵_new/𝒵_guess < 1.001)
    #             global 𝜌_guess = 𝜌_guess * 𝒵_new / 𝒵_guess
    #             gjfger_1st_scalingfactor(𝜌_guess)
    #             rungaussian("Ger")
    #             global 𝒵_guess = 𝒵_new
    #             global 𝒵_new = calc_𝒵_new(𝒵_guess, 𝑅𝑟𝑒𝑓)
    #         ##RC-101221!
    #         end
    #         # lattice Coulomb energy
    #         𝑠 = calc_𝑠()
    #         𝑊ₑ = -𝛼ᵣ * abs(charge)^2 / 𝑟₀ ./ 𝑠
    #         𝑑𝑊ₑ╱𝑑𝑠 = -𝑊ₑ ./ 𝑠

    #         # lattice polarization energy
    #         𝜀 = calc_𝜀()
    #         𝛼ₚₒₗ = 0.5(1 .- 1 ./ 𝜀)
    #         𝑊ₚₒₗ = @. -𝛼ₚₒₗ * abs(charge)^2 / 𝑠 / 𝑅𝑟𝑒𝑓
    #         𝑑𝑊ₚₒₗ╱𝑑𝑠 = @. -𝑊ₚₒₗ / 𝑠 * (1 + 3/𝜀)

    #         # xp-pcm energy, 𝐺ₑᵣ with polarization contribution and 𝐸ᵣ without
    #         ##RC-101221!
    #         𝒵 = @. 𝒵_new / 𝑠^(3 + 𝜂)
    #         𝜌 = 𝒵 * sp[3] / sp[4] / 0.063
    #         gjfgeranalytical(𝜌)
    #         𝜌 = @. sp[2] / 𝑠^(3 + 𝜂)
    #         𝒵 = 𝜌*0.063* sp[4] / sp[3]
    #         writegjf("Ger")
    #         ##RC-101221!
    #         rungaussian("Ger")
    #         𝐺ₑᵣ = get_data("Ger.log", "SCF Done", 5)
    #         𝑊ₚₒₗ′ = get_𝑊ₚₒₗ′()
    #         𝐸ₚₐᵤₗᵢ = get_𝐸ₚₐᵤₗᵢ()
    #         𝑒𝑓𝑔╱𝑛𝑡𝑠 = get_𝑒𝑓𝑔╱𝑛𝑡𝑠()
    #         𝑅 = 𝑅𝑟𝑒𝑓 * 𝑠
    #         𝐼₂ = @. -4π * 𝑅𝑟𝑒𝑓 * 𝑅^2 * 𝑒𝑓𝑔╱𝑛𝑡𝑠
    #         𝑑𝐸ᵣ╱𝑑𝑠 = @. -(3 + 𝜂) * 𝐸ₚₐᵤₗᵢ / 𝑠 + 𝒵 * 𝐼₂
    #         # total lattice energy
    #         𝑊ₗ = 𝑊ₑ + 𝐺ₑᵣ
    #         𝑑𝑊ₗ╱𝑑𝑠 = 𝑑𝑊ₑ╱𝑑𝑠 + 𝑑𝑊ₚₒₗ╱𝑑𝑠 + 𝑑𝐸ᵣ╱𝑑𝑠

    #         # unit cell volume per formula unit
    #         if lattice == "NaCl"
    #             𝑎_cell = 2𝑟₀
    #             𝑉_cell = @. (𝑎_cell * 𝑠 * 0.529177)^3 / 4 # 1 bohr = 0.529177 Å; devided by 4 because there are 4 formula units of NaCl
    #             𝑑𝑉_cell╱𝑑𝑠 = @. 3𝑉_cell / 𝑠
    #         end
    #         if lattice == "CsCl"
    #             𝑎_cell = 2𝑟₀ / √3
    #             𝑉_cell = @. (𝑎_cell * 𝑠 * 0.529177)^3 # 1 bohr = 0.529177 Å
    #             𝑑𝑉_cell╱𝑑𝑠 = @. 3𝑉_cell / 𝑠
    #         end
    #         ##!!RC121121
    #         if lattice == "noLattice"
    #             𝑉_cell =  𝑉𝑐
    #             𝑑𝑉_cell╱𝑑𝑠 = @. 3𝑉_cell / 𝑠
    #         end
    #         ##!!RC121121
    #         # analytical pressure
    #         𝑝 = @. -𝑑𝑊ₗ╱𝑑𝑠 / 𝑑𝑉_cell╱𝑑𝑠 * 4359.7 # 1 hartree/bohr = 4359.7 GPa
    #         #𝑉_cell╱𝑉₀ = 𝑉_cell / 𝑉_cell[1]
    #     end

    #     if model == "pointcharges"

            𝑠 = calc_𝑠()
            sp = solventparameters()
            𝑟ₐ = atomicradii()
            atoms = atomlist()
            𝑅𝑟𝑒𝑓 = 𝑓[1] * 𝑟ₐ[atoms[1]] * 1.88973 # reference radius of Cl- in bohr

            if impose_equilibrium # at 1st scalingfactor so that p(f0)=0
                open("iterativeZ.dat", "w") do file end # erase the file content if exists
                𝜌_old = sp[2]  # initial solvent density
                𝒵_old = 0.063 * 𝜌_old * sp[4] / sp[3]
                𝒵_new = 𝒵_old * 1.1  # a guess of 𝒵_new
                while !(0.999 < 𝒵_new/𝒵_old < 1.001) # self-consistent calculation of 𝒵
                    global 𝜌_old = 𝜌_old * 𝒵_new / 𝒵_old
                    gjfger_1st_scalingfactor(𝜌_old)
                    rungaussian("Ger")
                    global 𝒵_old = 𝒵_new
                    global 𝒵_new = calc_𝒵_new_pointcharges(𝒵_old, 𝑅𝑟𝑒𝑓)
                end
                𝒵 = @. 𝒵_new / 𝑠^(3 + 𝜂)
                𝜌 = 𝒵 * sp[3] / sp[4] / 0.063
                gjfgeranalytical(𝜌)  # using self-consistently determined 𝜌 and 𝒵 for Ger calculation
            else  # not impose_equilibrium
                𝜌 = @. sp[2] / 𝑠^(3 + 𝜂)
                𝒵 = 𝜌 * 0.063 * sp[4] / sp[3]
                writegjf("Ger")  # using initial 𝜌 and 𝒵 for Ger calculation
            end

            rungaussian("Ger")

            function finitedifference(energy::Vector{Float64}, 𝑠=𝑠)
                derivative = zeros(length(energy))
                derivative[1] = (energy[2] - energy[1]) / (𝑠[2] - 𝑠[1])
                for j in 2:length(𝑓)-1
                    derivative[j] = (energy[j+1] - energy[j-1]) / (𝑠[j+1] - 𝑠[j-1])
                end
                derivative[end] = (energy[end] - energy[end-1]) / (𝑠[end] - 𝑠[end-1])
                derivative
            end

            # ion-medium polarization energy
            𝑊ₚₒₗ′ = get_𝑊ₚₒₗ′()
            𝑑𝑊ₚₒₗ′╱𝑑𝑠 = finitedifference(𝑊ₚₒₗ′)
            𝜀 = calc_𝜀()
            𝛼ₚₒₗ = 0.5(1 .- 1 ./ 𝜀)
            𝑊ₚₒₗ = @. -𝛼ₚₒₗ * abs(charge)^2 / 𝑠 / 𝑅𝑟𝑒𝑓  # 𝑊ₚₒₗ is an approximate to 𝑊ₚₒₗ′
            𝑑𝑊ₚₒₗ╱𝑑𝑠 = @. -𝑊ₚₒₗ′ / 𝑠 * (1 + 3/𝜀)

            # Pauli repulsion energy
            𝐸ₚₐᵤₗᵢ = get_𝐸ₚₐᵤₗᵢ()
            𝑑𝐸ₚₐᵤₗᵢ╱𝑑𝑠ₙ = finitedifference(𝐸ₚₐᵤₗᵢ)
            𝑒𝑓𝑔╱𝑛𝑡𝑠 = get_𝑒𝑓𝑔╱𝑛𝑡𝑠()
            𝑅 = 𝑅𝑟𝑒𝑓 * 𝑠
            𝐼₂ = @. -4π * 𝑅𝑟𝑒𝑓 * 𝑅^2 * 𝑒𝑓𝑔╱𝑛𝑡𝑠
            𝑑𝐸ᵣ╱𝑑𝑠 = @. -(3 + 𝜂) * 𝐸ₚₐᵤₗᵢ / 𝑠 + 𝒵 * 𝐼₂

            # nuclei-charges and charges-charges Coulomb energies
            𝐸_nuclei_charges = get_data("Ger.log", "Nuclei-charges interaction", 4)
            𝐸_charges_charges = get_data("Ger.log", "Self energy", 7)
            𝑑𝐸_nuclei_charges╱𝑑𝑠 = -𝐸_nuclei_charges ./ 𝑠
            𝑑𝐸_charges_charges╱𝑑𝑠 = -𝐸_charges_charges ./ 𝑠
            𝑑𝐸_nuclei_charges╱𝑑𝑠ₙ = finitedifference(𝐸_nuclei_charges)
            𝑑𝐸_charges_charges╱𝑑𝑠ₙ = finitedifference(𝐸_charges_charges)

            # electrons-charges Coulomb energy
            𝐸_electrons_charges = get_EFM()[2]
            𝑑𝐸_electrons_charges╱𝑑𝑠 = 𝑅𝑟𝑒𝑓 * get_EFM()[1]
            𝑑𝐸_electrons_charges╱𝑑𝑠ₙ = finitedifference(𝐸_electrons_charges)

            # total energy
            𝐺ₑᵣ = get_data("Ger.log", "SCF Done", 5)
            𝑑𝐺ₑᵣ╱𝑑𝑠ₙ = finitedifference(𝐺ₑᵣ)
            𝐸ₜₒₜ = 𝐺ₑᵣ
            𝑑𝐸ₜₒₜ╱𝑑𝑠 = 𝑑𝑊ₚₒₗ╱𝑑𝑠 + 𝑑𝐸ᵣ╱𝑑𝑠 + 𝑑𝐸_nuclei_charges╱𝑑𝑠 + 𝑑𝐸_charges_charges╱𝑑𝑠 + 𝑑𝐸_electrons_charges╱𝑑𝑠

            # analytical pressure
            𝑑𝑉𝑐╱𝑑𝑠 = @. 3𝑉𝑐 / 𝑠
            𝑝ₐ = @. -𝑑𝐸ₜₒₜ╱𝑑𝑠 / 𝑑𝑉𝑐╱𝑑𝑠 * 4359.7 # 1 hartree/bohr = 4359.7 GPa

            # numerical pressure
            𝑝ₙ = calc_numerical𝑝(𝑉𝑐, 𝐸ₜₒₜ)  # by Murnaghan ESO fitting
            𝑝ₙ′ = -finitedifference(𝐸ₜₒₜ) ./ 𝑑𝑉𝑐╱𝑑𝑠 * 4359.7 # by finite difference

        #end
        # print output
        writeproperties3()
        #debug2()
    #end
    write("1.sh", "rm -rf fort.* *.off Vc-*.gjf Vc-*.log")
    run(`bash 1.sh`)
    run(`rm -rf 1.sh`)
#end

#main()
