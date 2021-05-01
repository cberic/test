# 2 checkpoint files, old.chk
# custum basis set keyword gen 

using Printf
using LsqFit
using DelimitedFiles

include("input.jl")

# Check whether pressure is to be calculated analytically or numerically.
# isnumerical = true if pressurecalc is defined and its value is "numerical"
isnumerical = @isdefined(pressurecalc) && pressurecalc == "numerical"
# When isnumerical = true, load the LsqFit package
#isnumerical && using LsqFit
# if false, load DelimitedFiles
#!isnumerical && using DelimitedFiles

#------------------------------------------------------------------------------
# solvent
#------------------------------------------------------------------------------
# returns a 5-element tuple of (dielectric constant 𝜀, valence electron density 𝜌,
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
            "Rn"=> 2.20,    "86"=> 2.20
        # add more if needed from https://en.wikipedia.org/wiki/Van_der_Waals_radius
        )
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
            "K" => 2.34,    "19"=> 2.34,
            "Ca"=> 2.70,    "20"=> 2.70,
            #4p
            "Ga"=> 2.33,    "31"=> 2.33,
            "Ge"=> 2.34,    "32"=> 2.34,
            "As"=> 2.31,    "33"=> 2.31,
            "Se"=> 2.24,    "34"=> 2.24,
            "Br"=> 2.19,    "35"=> 2.19,
            "Kr"=> 2.12,    "36"=> 2.12,
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
            "Rn"=> 2.43,    "86"=> 2.43
        # add more if needed from https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/chem.201700610
        )
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
    # "*" is string concatenation operator
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
        # get the 2nd-last fields of each coordinate line & add a leading space
        coordlines[i] = " " * split(lines[i], limit=2)[2]
    end
    return coordlines
end


#------------------------------------------------------------------------------
# io
#------------------------------------------------------------------------------
function writegjf(jobtype)
    if jobtype == "Vc"        # Cavitation volume jobs
        isnumerical ? gjfvc() : gjfvc2()
    elseif jobtype == "Ger"   # Electronic energy SCRF jobs
        isnumerical ? gjfgernumerical() : gjfgeranalytical()
    end
end


# write gjf files for Vc jobs
function gjfvc(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    g = tidygeometries(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    coordlines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    
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


# write gjf files for Vc jobs in analytical pressure calculation
function gjfvc2(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    g = tidygeometries(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    coordlines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    
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


# write gjf files for Ger calculations for numerical p
function gjfgernumerical(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    g = tidygeometries(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    coordlines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    𝜀 = calculate𝜀()    # data of 𝜀 and 𝑍 needed for the gjf files
    𝑍 = calculate𝑍()

    for j in 1:a
        # writing mode for the first and appending mode for other 𝑓
        open("Ger.gjf", "$(j == 1 ? "w" : "a")") do file
            write(file, """
                %chk=Ger.chk
                %nproc=$nproc
                %mem=$mem
                #p $keywords $(j == 1 ? "" : "guess=read")
                # scrf=(iefpcm,solvent=$solvent,read) nosym 6d 10f
                
                scaling factor = $(𝑓[j])
                
                $charge $multiplicity
                $g
                
                qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae
                nsfe=$noa
                nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])
                eps=$(𝜀[j]) rhos=$(𝑍[j])
                
                """)
            
            for k in 1:noa
                write(file, " $(coordlines[k])    $(𝑟ₐ[atoms[k]])    $(𝑓[j])\n")
            end
            
            write(file, "\n")
            # do not write "--link1--" for the last scaling factor
            if j != a
                write(file, "--link1--\n")
            end
        end
    end
end


# write gjf files for Ger calculations for analytical p
function gjfgeranalytical(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    g = tidygeometries(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    coordlines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    𝜀 = calculate𝜀()    # data of 𝜀 and 𝑍 needed for the gjf files
    𝑍 = calculate𝑍()
    n = getnumberoftesserae()
    
        # writing mode for the first and appending mode for other 𝑓
    open("Ger.gjf", "w") do file
        for j in 1:a
            write(file, """
                %chk=$(𝑓[j]).chk
                %nproc=$nproc
                %mem=$mem
                #p $keywords
                # scrf=(iefpcm,solvent=$solvent,read) iop(5/33=1) prop(efg,grid) nosym 6d 10f
                
                scaling factor = $(𝑓[j])
                
                $charge $multiplicity
                $g
                
                qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae
                nsfe=$noa
                nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])
                eps=$(𝜀[j]) rhos=$(𝑍[j])
                
                """)
            
            for k in 1:noa
                write(file, " $(coordlines[k])    $(𝑟ₐ[atoms[k]])    $(𝑓[j])\n")
            end
            write(file, "\n")
            write(file, "$(n[j]), 1, $(round(Int, 𝑓[j]*1000)), $(round(Int, 𝑓[j]*1000+1))\n")
            write(file, "\n")
            # do not write "--link1--" for the last scaling factor
            if j != a
                write(file, "--link1--\n")
            end
        end
    end
end


# extract volume 𝑉𝑐 data from gaussian output files
function get𝑉𝑐(𝑓 = scalingfactors)
    a = length(𝑓)
    𝑉𝑐 = Array{Float64}(undef, a)
    j = 1    # j ranges from 1:length(𝑓)
    open("Vc.log", "r") do file
        for line in eachline(file)
            if occursin("GePol: Cavity volume", line)
                𝑉𝑐[j] = parse(Float64, split(line)[5])
                j += 1    # j ranges from 1:length(𝑓)
            end 
        end
    end
    return 𝑉𝑐
end


# extract the total number of tesserae from tesserae-𝑓.off
function getnumberoftesserae(𝑓 = scalingfactors)
    a = length(𝑓)
    n = Array{Int64}(undef, a)
    for j in 1:a
        open("tesserae-$(𝑓[j]).off", "r") do file
            readline(file)
            secondline=readline(file)
            # the 2nd number on the 2nd line is the number of tesserae
            n[j] = parse(Int64, split(secondline)[2])
        end
    end
    return n
end


# extract and calculate the xyz coordinates of the centers of the tesserae
# from tesserae-𝑓.off 
function writetesseragrid(geom = geometries, 𝑓 = scalingfactors)
    n = getnumberoftesserae()
    atoms = atomlist(geom)
    a = length(𝑓)
    𝑟ₐ = atomicradii()
    for j in 1:a
        tesseraecoordinates = Array{Float64}(undef, n[j],3)  # n*3 2D array
        x1 = x2 = x3 = y1 = y2 = y3 = z1 = z2 = z3 = 0.0
        linecount = 1
        i = 1   # i ranges from 1:n
        open("tesserae-$(𝑓[j]).off", "r") do file
            for line in eachline(file)
                if linecount >= 3 && linecount <= 3 * n[j] + 2
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
        writedlm("fort.$(round(Int, 𝑓[j]*1000))", tesseraecoordinates * 𝑟ₐ[atoms[1]] * 𝑓[j])
    end
end


# extract Electric Field Gradients from fort.1201 etc files
function getEFG(𝑓 = scalingfactors)
    n = getnumberoftesserae()
    a = length(𝑓)
    efgsum = zeros(a)
    for j in 1:a
        linecount = 1
        open("fort.$(round(Int, 𝑓[j]*1000+1))", "r") do file
            for line in eachline(file)
                if linecount % 4 == 2
                    # zz component of the EFG
                    efgsum[j] += parse(Float64, split(line)[4])
                elseif linecount % 4 == 3
                    # yy and xx components of the EFG
                    efgsum[j] += sum([parse(Float64, s) for s in split(line)[1:2]])
                end
                linecount += 1
            end
        end
    end
    return @. efgsum / 4 / pi / n
end

# extract electronic energy 𝐺𝑒𝑟 data from gaussian output files
function getPauli𝐸(𝑓 = scalingfactors)
    a = length(𝑓)
    Pauli𝐸 = Array{Float64}(undef, a)
    j = 1    # j ranges from 1:length(𝑓)
    open("Ger.log", "r") do file
        for line in eachline(file)
            if occursin("QRepSI", line)
                Pauli𝐸[j] = parse(Float64, split(line)[6]) / 627.503 # 1 hartree = 627.503 kcal/mol
            elseif occursin("SCF Done", line)
                j += 1    # j ranges from 1:length(𝑓)
            end
        end
    end
    return Pauli𝐸
end


# extract electronic energy 𝐺𝑒𝑟 data from gaussian output files
function get𝐺𝑒𝑟(𝑓 = scalingfactors)
    a = length(𝑓)
    𝐺𝑒𝑟 = Array{Float64}(undef, a)
    j = 1    # j ranges from 1:length(𝑓)
    open("Ger.log", "r") do file
        for line in eachline(file)
            if occursin("SCF Done", line)
                𝐺𝑒𝑟[j] = parse(Float64, split(line)[5])
                j += 1    # j ranges from 1:length(𝑓)
            end
        end
    end
    return 𝐺𝑒𝑟
end


function getorbitalenergy(𝑓 = scalingfactors)
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


# extract edensity from .cube file
function getedensity(𝑓 = scalingfactors)
    a = length(𝑓)
    edensity = zeros(a)
    for j in 1:a
        open("$(𝑓[j]).cube", "r") do file
            for line in eachline(file)
                # negative sign for a negative number of edensity
                edensity[j] -= parse(Float64, split(line)[4])
            end
        end
    end
    n = getnumberoftesserae()
    return @. edensity / n
end


# use the Printf package to write the properties.dat file
function writeproperties(𝑉𝑐 = 𝑉𝑐, 𝐺𝑒𝑟 = 𝐺𝑒𝑟, 𝑓 = scalingfactors)
    a = length(𝑓)
    𝑠 = calculate𝑠()
    𝜀 = calculate𝜀()
    𝑍 = calculate𝑍()
    # for numerical p, call calculatenumerical𝑝(); for analytical p, call calculateanalytical𝑝()
    #𝑝 = isnumerical ? calculatenumerical𝑝() : calculateanalytical𝑝() 
    𝑝n = calculatenumerical𝑝()
    𝑝a = calculateanalytical𝑝()
    Eorbital = getorbitalenergy()
    open("properties.dat", "w") do file
        write(file, "#    𝑓       𝑉𝑐(𝑓) Å³   𝑠(𝑓)         𝜀(𝑠)        𝑍(𝑠)        𝐺𝑒𝑟(𝑓) a.u.     𝑝(𝑓)-numeric. -analyt.GPa\n")
        for j in 1:a
            @printf(file, "%d    %.3f     %7.3f    %.6f    %.6f    %9.6f    %.8f    %6.3f    %6.3f\n", 
                            j,   𝑓[j],   𝑉𝑐[j],   𝑠[j],    𝜀[j],   𝑍[j],   𝐺𝑒𝑟[j],  𝑝n[j],  𝑝a[j])
        end
        write(file, "\n")
        for j in 1:a
            @printf(file, "𝑓 = %.3f    𝑝 = %6.3f GPa ----orbital energies in a.u.----\n", 𝑓[j], 𝑝a[j])
            write(file, Eorbital[j])
            write(file, "\n")
        end
    end
end


function debug(𝑉𝑐 = 𝑉𝑐, 𝐺𝑒𝑟 = 𝐺𝑒𝑟, 𝑓 = scalingfactors)
    a = length(𝑓)
    𝑠 = calculate𝑠()
    𝜀 = calculate𝜀()
    𝑍 = calculate𝑍()
    PauliE = getPauli𝐸()
    edensity1 = getedensity()
    edensity2 = getEFG()
    alpha = calculateAlpha()
    𝑝n = calculatenumerical𝑝()
    𝑝a = calculateanalytical𝑝()
    Eorbital = getorbitalenergy()
    open("debug.dat", "w") do file
        write(file, "#    𝑓       𝑉𝑐(𝑓) Å³   𝑠(𝑓)         𝜀(𝑠)        𝑍(𝑠)        𝐺𝑒𝑟(𝑓) a.u.     𝑝a(𝑓) GPa      PauliE(𝑓)     edensity     efg/nts     Alpha(𝑓)\n")
        for j in 1:a
            @printf(file, "%d    %.3f     %7.3f    %.6f    %.6f    %9.6f    %.8f    %6.3f    %9.6f    %9.6f    %9.6f    %9.6f\n", 
                            j,   𝑓[j],   𝑉𝑐[j],   𝑠[j],    𝜀[j],   𝑍[j],   𝐺𝑒𝑟[j],  𝑝a[j],  PauliE[j], edensity1[j], edensity2[j], alpha[j])
        end
    end
end


#------------------------------------------------------------------------------
# rungaussian
#------------------------------------------------------------------------------
# check if g16 or g09 is installed and loaded
function gaussianversion()
    if typeof(Sys.which("g16")) === String
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
function calculate𝑠(𝑉𝑐 = 𝑉𝑐)  # 𝑉𝑐 is a global variable
    # @. is a macro for vectoried operation/broadcasting
    return @. cbrt(𝑉𝑐/𝑉𝑐[1])
end


# dielectric permitivity 𝜀 = 1 + (𝜀₀-1)/𝑠̄³
function calculate𝜀()
    𝜀₀ = solventparameters()[1]
    𝑠 = calculate𝑠()
    return @. 1 + (𝜀₀ - 1) / 𝑠^3    # 1D array of length a
end


# Pauli repulson barrier 𝑍 = 𝑍₀/𝑠̄⁽³⁺𝜂⁾, where 𝑍₀ = 𝜌₀
function calculate𝑍(𝜂=𝜂)
    𝜌₀ = solventparameters()[2]
    𝑠 = calculate𝑠()
    return @. 𝜌₀ / 𝑠^(3+𝜂)    # 1D array of length a
end


# Murnaghan equation of state fitting for pressure 𝑝 calculation
# using LsqFit
function eosfitting(𝑉𝑐 = 𝑉𝑐, 𝐺𝑒𝑟 = 𝐺𝑒𝑟)
    abc_parameters = Array{Float64}(undef, 3)
    # y = (a/b)*(1/x)**b+(a-c)*x; y is Ger, x is Vc
    # a=p[1], b=p[2], c=p[3], x is Vc
    @. model(x, p) = (p[1]/p[2])*x^(-p[2]) + (p[1]-p[3])*x
    xdata = 𝑉𝑐 ./ 𝑉𝑐[1]
    ydata = 𝐺𝑒𝑟 .- 𝐺𝑒𝑟[1]
    p0 = [0.0, 5.0, 0.0]
    fit = curve_fit(model, xdata, ydata, p0)
    abc_parameters[1] = fit.param[1]/𝑉𝑐[1]
    abc_parameters[2] = fit.param[2]
    abc_parameters[3] = fit.param[3]/𝑉𝑐[1]
    return abc_parameters
end


function calculatenumerical𝑝(𝑉𝑐 = 𝑉𝑐)
    𝑎, 𝑏, 𝑐 = eosfitting()
    # 1 hartree/Å³ = 4359.74417 GPa
    return @. (𝑎 * ( (𝑉𝑐[1]/𝑉𝑐)^(𝑏+1) - 1 ) + 𝑐) * 4359.74417 
end


# eq (9) in DOI:10.1002/jcc.25544
# In this eq, 𝒵 is called Alpha in Gaussian
function calculateAlpha()
    sp = solventparameters()
    # Alpha = fA*RhoS*NVES/SolvMW  ; fA=0.063d0
    # Alpha0 = 0.063 * sp[2] * sp[4] / sp[3]
    𝜌 = calculate𝑍()
    # $alpha[$n_fact]=0.063*$rhos[$n_fact]*$nves/$mw;
    return @. 0.063 * 𝜌 * sp[4] / sp[3]
end


# eq (24) in DOI:10.1002/jcc.25544
function calculateanalytical𝑝(𝜂 = 𝜂, 𝑉𝑐 = 𝑉𝑐)
    Pauli𝐸 = getPauli𝐸()
    # $p0[$n_fact]=(3.0+$eta)/($three*($volume[$n_fact]*1.88973**3))*$qrep[$n_fact]; 
    #@. firstterm = (3 + 𝜂) / (3 * 𝑉𝑐 * 1.88973^3) * Pauli𝐸  # 1 angstrom = 1.88973 bohr; all in atomic units
    # $alpha[$n_fact]*$edensity[$n_fact]/$nts;  
    alpha = calculateAlpha()
    # if efg == true, use getEFG() to calculate electron density, otherwise use getedensity()
    @isdefined(efg) && efg == true ? edensity = getEFG() : edensity = getedensity()
    #@. secondterm = alpha * edensity * 1.88973^3
    return @. ((3 + 𝜂)/(3 * 𝑉𝑐 * 1.88973^3) * Pauli𝐸 - alpha * edensity) * 1.88973^3 * 4359.74417 # 1 hartree/Å³ = 4359.74417 GPa, 
end

#------------------------------------------------------------------------------
# main
#------------------------------------------------------------------------------
function main(𝑓 = scalingfactors)    
    # Step 1: cavity volume 𝑉𝑐(𝑓) Gaussian jobs and solvent property calculations
    writegjf("Vc")
#    if isnumerical
#        rungaussian("Vc") 
#    else 
        for j in 𝑓
            rungaussian("Vc-$j")
            open("Vc.log", "$(j == first(𝑓) ? "w" : "a")") do file
                write(file, read(`cat Vc-$j.log`, String))
            end
            run(`cp tesserae.off tesserae-$j.off`)
        end
#    end
    global 𝑉𝑐 = get𝑉𝑐()

    # Step 2: electronic structure Gaussian jobs
    writetesseragrid()
    writegjf("Ger")
    rungaussian("Ger")
    global 𝐺𝑒𝑟 = get𝐺𝑒𝑟()

    # Step 3: analytical pressure calculation
    # formchk K_xp-060.chk 42.fchk 
    # cubegen 0 density=scf 42.fchk 42.cube -5 < fort.42
#    if !isnumerical
        for j in 𝑓
            run(`formchk $j.chk $j.fchk`)
            #write("1.sh", "cubegen 0 density=scf $j.fchk $j.cube -5 < $j.tsrcoord")
            write("1.sh", "cubegen 0 density=scf $j.fchk $j.cube -5 < fort.$(round(Int, j*1000))")            
            run(`bash 1.sh`)
        end
        run(`rm -rf 1.sh`)
        debug()
#    end
    # print results to properties.dat file
    writeproperties()
end

main()
