using Statistics
using Printf
using DelimitedFiles

include("input.jl")


# function definitions below
#------------------------------------------------------------------------------
# solvent.jl
#------------------------------------------------------------------------------
# returns a 5-element tuple of solvent properties: 
# (dielectric constant 𝜀, density 𝜌, molar mass 𝑀, number of valence electrons, molecular radius (Ang.) )
function solventparameters(s = solvent)
    if s == "cyclohexane"
        return (2.0165, 0.7781, 84.1595, 36, 2.815)
    elseif s == "benzene"
        return (2.2706, 0.8756, 78.1118, 30, 2.63)
    elseif s == "argon"
        return (1.43, 1.3954, 39.948, 8, 1.705)
    else 
        println("Solvent not implemented. Try cyclohexane, benzene, or argon")
    end
end

#------------------------------------------------------------------------------
# atomicradii.jl
#------------------------------------------------------------------------------
function atomicradii()
    return Dict(
    "H" => 1.20,    "1" => 1.20,
    "B" => 1.92,    "5" => 1.92,  
    "C" => 1.70,    "6" => 1.70,
    "N" => 1.55,    "7" => 1.55,
    "O" => 1.52,    "8" => 1.52,
    "F" => 1.47,    "9" => 1.47)
end

#------------------------------------------------------------------------------
# geometries.jl
#------------------------------------------------------------------------------
function tidygeometries(geom = geometries)
    # remove leading and trailing spaces/blank lines 
    # and add a space to the begininig
    # "*" is string concatenation operator
    a = " " * strip(geom)

    # change block separators to "NEXT"
    # r"" represents a regular expression
    b = replace(a, r"\s*\n+\s*\n+\s*" => "NEXT")
    
    # remove leading and trailing spaces of each coordinate line 
    # and add a space to the begininig of coordinate line
    c = replace(b, r"\s*\n\s*" => "\n ")
    
    # change the block separators back to "\n\n" and add a space
    return replace(c, "NEXT" => "\n\n ")
end


#= the function below generates an array of blocks; one block is shown
 C    1.10712900   -1.48465400    0.00000000
 C   -0.00001100   -0.72877000    0.00000000
 C    0.00000000    0.72875800    0.00000000
 C   -1.10712300    1.48466300    0.00000000
 H    2.09998900   -1.03986800    0.00000000
 H    1.05995700   -2.56940600    0.00000000
 H   -0.97811200   -1.21104800    0.00000000
 H    0.97811000    1.21101800    0.00000000
 H   -2.09999700    1.03991100    0.00000000
 H   -1.05991600    2.56941400    2.00000000
=#
function structures(geom = geometries)
    return split(tidygeometries(geom), "\n\n", keepempty=false)
end


function numberofstructures(geom = geometries)
    return length(structures(geom))
end


#= the function below generates an array of arrays of lines
 C    1.10712900   -1.48465400    0.00000000
=#
function structurelines(geom = geometries)
    blocks = structures(geom)
    nos = numberofstructures(geom)
    # lines is an array of arrays of length nos
    # lines[i] is a 1D array of length numberofatoms
    lines = Array{Array{String}}(undef, nos)
    for i in 1:nos
        lines[i] = split(blocks[i], "\n", keepempty=false)
    end
    return lines
end


function numberofatoms(geom = geometries)
    return length(structurelines(geom)[1])
end

#= the function below generates an array of atoms (no coordinates)
 C
=#
function atomlist(geom = geometries)
    nos = numberofstructures(geom)
    lines = structurelines(geom)
    noa = numberofatoms(geom)
    atoms = Array{String}(undef, nos,noa)  # 2D array of length nos * noa
    for i in 1:nos
        for j in 1:noa
            atoms[i,j] = split(lines[i][j], keepempty=false, limit=2)[1]
        end
    end
    return atoms
end

# a slower function with much more memory allocations
#function atomlist2(geometries)
#    noa = numberofatoms(geometries)
#    atoms = Array{String}(undef, noa) 
#    for i in 1:noa
#        atoms[i] = split(geometries)[4*i-3]
#    end
#    return atoms
#end


#= the function below generates an array of arrays of coordinate lines (no atom label)
 1.10712900   -1.48465400    0.00000000
=#
function coordinatelines(geom = geometries)
    nos = numberofstructures(geom)
    lines = structurelines(geom)
    noa = numberofatoms(geom)
    for i in 1:nos
        for j in 1:noa
            lines[i][j] = " " * split(lines[i][j], keepempty=false, limit=2)[2]
        end
    end
    return lines    # array (of length nos) of arrays of length noa
end

# a varient of the above function using a 2D array, instead of array of arrays
#function coordinatelines2(geom = geometries)
    #nos = numberofstructures(geom)
    #lines = structurelines(geom)
    #noa = numberofatoms(geom)
    #array = Array{String}(undef, nos,noa)
    #for i in 1:nos
     #   for j in 1:noa
     #       array[i,j] = " " * split(lines[i][j], keepempty=false, limit=2)[2]
     #   end
    #end
    #return array
#end


function atomcoordinates(geom = geometries)
    # split using the default space dilimiter, generating a very long 1D array
    array = split(geom)
    nos = numberofstructures(geom)
    noa = numberofatoms(geom)

    # 3D array storing each coordinate
    # the last dimension contains 3 elements corresponding to x, y, and z coordinates
    coordinates = Array{String}(undef, nos,noa,3)
    for i in 1:nos
        for j in 1:noa
            for k in 1:3
                coordinates[i,j,k] = array[(i-1)*40 + (j-1)*4 + k+1]
            end
        end
    end
    return coordinates
end


#------------------------------------------------------------------------------
# io.jl
#------------------------------------------------------------------------------
function writegjf(jobtype)
    if jobtype == "Vc"
        gjfvc()
    elseif jobtype == "Ger"
        gjfger()
    elseif jobtype == "Gcav"
        gjfgcav()
    end
end


# write gjf files for Vc calculations
function gjfvc(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    st = structures(geom)
    nos = numberofstructures(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    lines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    
    run(`mkdir -p tmp`)  # make a tmp dir in the working folder

    # loop through structures; use nproc=1 and multithreading
    Threads.@threads for i in 1:nos
        for j in 1:a   # loop through scaling factors
            # use writing mode for the first scalingfactor 
            # and appending mode for the rest
            open("tmp/structure-$i-Vc.gjf", "$(j == 1 ? "w" : "a")") do file
                # job killed after L301 (%kjob l301);
                write(file, """
                %kjob l301
                %nproc=1
                %mem=1000mb
                #p $keywords
                # scrf=(iefpcm,solvent=$solvent,read) nosym guess=only pop=none

                title

                $charge $multiplicity
                $(st[i])

                qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae
                nsfe=$noa
                nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])
                eps=$(sp[1]) rhos=$(sp[2])                  

                """)

                for k in 1:noa
                    write(file, " $(lines[i][k])  $(𝑟ₐ[atoms[i,k]])  $(𝑓[j])\n")
                end

                write(file, "\n")

                if j != a  # do not write --link1-- for the last scaling factor
                    write(file, "--link1--\n")
                end
            end
        end
    end
end


# write gjf files for Ger calculations
function gjfger(geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    st = structures(geom)
    nos = numberofstructures(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    lines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    𝜀 = calculate𝜀()    # needs data of 𝜀 and 𝜌 at different scalingfactors
    𝑍 = calculate𝑍()
    
    Threads.@threads for i in 1:nos  # loop through structures and multithreading
        for j in 1:a   # loop through scaling factors
            # use writing mode for the first scalingfactor and appending mode for the rest
            open("tmp/structure-$i-Ger.gjf", "$(j == 1 ? "w" : "a")") do file
                # job killed after L502
                write(file, """
                %kjob l502
                %chk=structure-$i-Ger.chk
                %nproc=$nproc
                %mem=$mem
                #p $keywords
                # scrf=(iefpcm,solvent=$solvent,read) nosym 6d 10f $(j != 1 ? "guess=read" : "")

                title

                $charge $multiplicity
                $(st[i])

                qrep pcmdoc geomview nodis nocav g03defaults tsare=$tesserae
                nsfe=$noa
                nvesolv=$(sp[4]) solvmw=$(sp[3]) rsolv=$(sp[5])
                eps=$(𝜀[j]) rhos=$(𝑍[j])                

                """)

                for k in 1:noa
                    write(file, " $(lines[i][k])  $(𝑟ₐ[atoms[i,k]])  $(𝑓[j])\n")
                end

                write(file, "\n")
                # do not write "--link1--" for the last scaling factor
                if j != a
                    write(file, "--link1--\n")
                end
            end
        end
    end
end


# write gjf files for Gcav calculations
function gjfgcav(cav = cavity, geom = geometries, 𝑓 = scalingfactors)
    a = length(𝑓)
    st = structures(geom)
    nos = numberofstructures(geom)
    noa = numberofatoms(geom)
    atoms = atomlist(geom)
    #lines = coordinatelines(geom)
    sp = solventparameters()
    𝑟ₐ = atomicradii()
    𝑉ₘ = calculate𝑉ₘ()    # needs data of molar volume of solvent 𝑉ₘ
    
    Threads.@threads for i in 1:nos  # loop through structures and multithreading
        for j in 1:a   # loop through scaling factors
            # writing mode for the first scalingfactor and appending mode for the rest
            open("tmp/structure-$i-Gcav.gjf", "$(j == 1 ? "w" : "a")") do file
                # job killed after L301
                # also testing for vdw cavity to add the noaddsph keyword
                write(file, """
                %kjob l301
                %nproc=1
                %mem=1000mb
                #p $keywords
                # scrf=(iefpcm,solvent=$solvent,read) nosym guess=only pop=none

                title

                $charge $multiplicity
                $(st[i])

                norep pcmdoc geomview nodis cav g03defaults $(cav == "vdw" ? "noaddsph" : "" ) tsare=$tesserae
                nsfe=$noa
                Vmol=$(𝑉ₘ[j]) rsolv=$(sp[5])

                """)

                for k in 1:noa
                    write(file, " $k  $(𝑟ₐ[atoms[i,k]] * 𝑓[1])  1.0\n")
                end

                write(file, "\n")
                # do not write "--link1--" after the last structure
                if j != a
                    write(file, "--link1--\n")
                end
            end
        end
    end
end


# extract volume 𝑉𝑐 from gaussian output files 
# 𝑉𝑐 as a function of scaling factor 𝑓
function get𝑉𝑐(geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝑉𝑐 = Array{Float64}(undef, nos,a)    # 2D array with dimensions nos * a
    Threads.@threads for i in 1:nos
        j = 1    # j indexes the length(𝑓)
        open("tmp/structure-$i-Vc.log") do file
            for line in eachline(file)
                if occursin("GePol: Cavity volume", line)
                    𝑉𝑐[i,j] = parse(Float64, split(line)[5])
                    j += 1    # j indexes the length(𝑓)
                end 
            end
        end
    end
    return 𝑉𝑐
end


# extract electronic energy 𝐺𝑒𝑟 from gaussian output files
function get𝐺𝑒𝑟(geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝐺𝑒𝑟 = Array{Float64}(undef, nos,a)    # 2D array with dimensions nos * a
    Threads.@threads for i in 1:nos
        j = 1    # j should index the length(𝑓)
        open("tmp/structure-$i-Ger.log") do file
            for line in eachline(file)
                if occursin("SCF Done", line)
                    𝐺𝑒𝑟[i,j] = parse(Float64, split(line)[5])
                    j += 1    # j should index the length(𝑓)
                end
            end
        end
    end
    return 𝐺𝑒𝑟
end


# extract non-electrostatic cavitation energy 𝐸𝑐𝑎𝑣 from gaussian output files
function get𝐸𝑐𝑎𝑣(geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝐸𝑐𝑎𝑣 = Array{Float64}(undef, nos,a)    # 2D array with dimensions nos * a
    Threads.@threads for i in 1:nos
        j = 1    # j should index the length(𝑓)
        open("tmp/structure-$i-Gcav.log") do file
            for line in eachline(file)
                if occursin("PCM non-electrostatic energy", line)
                    𝐸𝑐𝑎𝑣[i,j] = parse(Float64, split(line)[5])
                    j += 1    # j should index the length(𝑓)
                end
            end
        end
    end
    return 𝐸𝑐𝑎𝑣
end


# write the properties.dat file; the Printf package is used
function writeproperties(𝑉𝑐 = 𝑉𝑐, 𝐺𝑒𝑟 = 𝐺𝑒𝑟, 𝐸𝑐𝑎𝑣 = 𝐸𝑐𝑎𝑣, geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝑠 = calculate𝑠()
    𝑠̄ = average𝑠()
    𝜀 = calculate𝜀()
    𝑍 = calculate𝑍()
    𝑉ₘ = calculate𝑉ₘ()
    𝑝 = calculate𝑝()
    𝑝̄ = average𝑝()
    𝐺𝑐𝑎𝑣 = calculate𝐺𝑐𝑎𝑣()
    𝐺𝑡𝑜𝑡 = calculate𝐺𝑡𝑜𝑡()
    open("properties.dat", "w") do file
        for i in 1:nos
            write(file, "structure $i\n")
            write(file, "#    𝑓       𝑉𝑐(𝑓)      𝑠(𝑓)         𝑠̄(𝑓,𝑛ₛ)      𝜀(𝑠̄)        𝑍(𝑠̄)        𝑉ₘ(𝑠̄)      𝐺𝑒𝑟(𝑓)            𝑝(𝑓)      𝑝̄(𝑠̄)      𝐸𝑐𝑎𝑣(𝑓)        𝐺𝑐𝑎𝑣(𝑓)        𝐺𝑡𝑜𝑡(𝑓)\n")
            for j in 1:a
                @printf(file, "%d    %.2f    %7.3f    %.6f    %.6f    %.6f    %.6f    %7.3f    %.8f    %7.3f    %6.3f    %.8f      %.8f      %.8f\n", 
                               j,    𝑓[j],   𝑉𝑐[i,j], 𝑠[i,j],  𝑠̄[j],   𝜀[j],   𝑍[j],  𝑉ₘ[j], 𝐺𝑒𝑟[i,j], 𝑝[i,j], 𝑝̄[j], 𝐸𝑐𝑎𝑣[i,j], 𝐺𝑐𝑎𝑣[i,j], 𝐺𝑡𝑜𝑡[i,j])
            end
            write(file, "\n")
        end
    end
end


#------------------------------------------------------------------------------
# rungaussian.jl
#------------------------------------------------------------------------------
# multithreading for short Vc and Cav calculations.
# add return value for output file "return "$jobtype finished at time()""
function rungaussian(jobtype; geom=geometries, multi = multithreading)
    nos = numberofstructures(geom)
    cd("tmp")
    if jobtype == "Vc" || jobtype == "Gcav"
        Threads.@threads for i in 1:nos
            run(`g16 structure-$i-$jobtype.gjf`)
        end
    elseif jobtype == "Ger" && multi == "on"
        Threads.@threads for i in 1:nos
            run(`g16 structure-$i-$jobtype.gjf`)
        end
    elseif jobtype == "Ger" && multi == "off"
        for i in 1:nos
            run(`g16 structure-$i-$jobtype.gjf`)
        end
    end
    cd("..")
end


#------------------------------------------------------------------------------
# arithmetic.jl
#------------------------------------------------------------------------------
# calculate linear scaling factor 𝑠, as the cubic root of volume scaling 𝑉𝑐/𝑉𝑐𝑓=1.2
function calculate𝑠(𝑉𝑐 = 𝑉𝑐, geom = geometries, 𝑓 = scalingfactors)
    #𝑉𝑐 = get𝑉𝑐()
    # @. is a maroc for vectoried operation
    return @. cbrt(𝑉𝑐/𝑉𝑐[:,1]) # 2D array of dimensions nos * a
end


# 𝑠̄ is the average of 𝑠 over all structures at the same scalingfactor 𝑓
function average𝑠()
    𝑠 = calculate𝑠()
    # use the mean function in the Staticstics package
    return mean(𝑠, dims=1)   # 1D array of length a
end


# calculate dielectric permitivity 𝜀 as a function of 𝑠̄
function calculate𝜀()
    𝜀₀ = solventparameters()[1]
    𝑠̄ = average𝑠()
    # 𝜀 = 1 + (𝜀₀-1)/𝑠̄³
    return @. 1 + (𝜀₀ - 1) / 𝑠̄^3    # using vectorized "dot" operator for 1D array
end


# calculate Pauli repulson barrier 𝑍 as a function of 𝑠̄ and 𝜂
function calculate𝑍(𝜂=𝜂)
    𝜌₀ = solventparameters()[2]
    𝑠̄ = average𝑠()
    # 𝑍 = 𝑍₀/𝑠̄⁽³⁺𝜂⁾   𝑍₀ = 𝜌₀
    return @. 𝜌₀ / 𝑠̄^(3+𝜂)    # using vectorized "dot" operator for 1D array
end


# calculate the molar volume of solvent 𝑉ₘ as a function of 𝑠̄
function calculate𝑉ₘ()
    (𝜌₀, 𝑀) = solventparameters()[2:3]
    𝑠̄ = average𝑠()
    # 𝑉ₘ = (𝑀/𝜌₀) * 𝑠̄³
    return @. (𝑀/𝜌₀) * 𝑠̄^3    # using vectorized "dot" operator for 1D array
end


# Murnaghan equation of state fitting to calculate pressure 𝑝;
# three methods - julia (default), python, or mathematica
function murnaghan(ft = fitting)
    if ft == "julia"
        return juliafitting()
    elseif ft =="python"
        return pythonfitting()
    elseif ft == "mathematica"
        return mathematicafitting()
    end
end


# load the LsqFit package; 
# needs to be installed first by "Pkg.add("LsqFit")"
#using LsqFit
function juliafitting()
    # to be implemented
end


# need to install the lmfit package `pip install lmfit`
# make sure to use the correct python verion
# note the python call with full path inside the function
function pythonfitting(𝑉𝑐 = 𝑉𝑐, 𝐺𝑒𝑟 = 𝐺𝑒𝑟, geom = geometries)
    nos = numberofstructures(geom)
    #𝑉𝑐 = get𝑉𝑐()
    #𝐺𝑒𝑟 = get𝐺𝑒𝑟()
    array = Array{Float64}(undef, nos,3)
    Threads.@threads for i in 1:nos
        script = """#!/usr/bin/env python
            #<examples/doc_model1.py>
            from numpy import sqrt, pi, exp, linspace, loadtxt, savetxt
            from lmfit import Model
            
            ### import matplotlib.pyplot as plt
            
            #data = loadtxt('Vc-Ger.dat')
            #x = data[:, 0]/data[0,0]
            #y = data[:, 1]-data[0,1]
            x = $(𝑉𝑐[i,:] ./ 𝑉𝑐[i,1])
            y = $(𝐺𝑒𝑟[i,:] .- 𝐺𝑒𝑟[i,1])
            
            def murnaghan(x, a, b, c):
                "1-d gaussian: murnaghan(x, a,b, c)"
                return (a/b)*(1/x)**b+(a-c)*x
                  
            gmod = Model(murnaghan)
            result = gmod.fit(y, x=x, a=0, b=5, c=0)
            
            print(result.fit_report())
            
            #<end examples/doc_model1.py>"""
        
        # read in fitting results and write to structure-$i-fitting.out files
        results = read(pipeline(`echo $script`, `/scratch/bochen/Python-2.7.18/bin/python`), String)
        #results = read(pipeline(`echo $script`, `python3`), String)
        open("tmp/structure-$i-fitting.out", "w") do file
            write(file, results)
        end

        # use an awk script to extract the a, b, c parameters of the equation of state.
        awkscript = raw"/a:/ {a=$2}; /b:/ {b=$2}; /c:/ {c=$2}; END {print a,b,c}"
        abc = read(`awk $awkscript "tmp/structure-$i-fitting.out"`, String)
        array[i,1] = parse(Float64, split(abc)[1])/𝑉𝑐[i,1]
        array[i,2] = parse(Float64, split(abc)[2])
        array[i,3] = parse(Float64, split(abc)[3])/𝑉𝑐[i,1]
    end
    return array   # 2D array of dimension nos * 3
end



function mathematicafitting(𝑉𝑐 = 𝑉𝑐, 𝐺𝑒𝑟 = 𝐺𝑒𝑟, geom = geometries)
    nos = numberofstructures(geom)
    #𝑉𝑐 = get𝑉𝑐()
    #𝐺𝑒𝑟 = get𝐺𝑒𝑟()
    array = Array{Float64}(undef, nos,3)
    Threads.@threads for i in 1:nos
        # write data to Vc-Ger.dat file which is then read by Mathematica
        # the DelimitedFiles package is used
        open("tmp/Vc-Ger.dat", "w") do file
            writedlm(file, [𝑉𝑐[i,:] 𝐺𝑒𝑟[i,:]])
        end
        script = """###!/usr/local/bin/WolframScript -script

                    (* generate high-precision samples of a mixed distribution *)

                    t = Import["tmp/Vc-Ger.dat", "Table"]

                    Print[murnaghan-eos]

                    << NonlinearRegression`

                    Print[NonlinearRegress[t, 
                    t[[1, 2]] + a*x ((1/b)*(t[[1, 1]]/x)^(b + 1) + 1) - c*x, {a, b, 
                    c}, x, RegressionReport -> {BestFitParameters, ParameterCITable, 
                    ANOVATable, BestFit, PredictedResponse}]]"""
        results = read(`math -script $script`, String)
        #results = read(pipeline(`echo $script`, `math`), String)
        open("tmp/structure-$i-fitting.out", "w") do file
            write(file, results)
        end
        # extraction of abc parameters to be implemented
    end
    rm("tmp/Vc-Ger.dat")
    return array    # 2D array of dimension nos * 3
end


function calculate𝑝(𝑉𝑐 = 𝑉𝑐)
    #𝑉𝑐 = get𝑉𝑐()    # 2D array of dimension nos * length(𝑓)
    abc = murnaghan("python")    # 2D array of dimension nos * 3 
    𝑎 = abc[:,1]   # 1D array of length nos
    𝑏 = abc[:,2]
    𝑐 = abc[:,3]
    #for i in 1:nos
     #   𝑎[i] = abc[i,1]
      #  𝑏[i] = abc[i,2]
       # 𝑐[i] = abc[i,3]
        #    for j in 1:a
            # $p[$n_fact]=($a * (($volume[1] / $volume[$n_fact]) ** ($b + 1) -1 ) + $c)*$Eh_o_Ang3_to_GPa;
            # 4359.74417 is the conversion factor from hartree/Å³ to GPa
            #𝑝[i,j] = (𝑎[i] * ((𝑉𝑐[i,1] / 𝑉𝑐[i,j]) ^ (𝑏[i] + 1) - 1) + 𝑐[i]) * 4359.74417
            #𝑝[i,j] = (𝑎[i] * (1 - (𝑉𝑐[i,1] / 𝑉𝑐[i,j]) ^ 𝑏[i] ) + 𝑐[i]) * 4359.74417
        #end
    #end
    # 4359.74417 is the conversion factor from hartree/Å³ to GPa
    return @. (𝑎 * ( (𝑉𝑐[:,1]/𝑉𝑐)^(𝑏+1) - 1 ) + 𝑐) * 4359.74417    # 2D array of dimension nos * length(𝑓)
    #return @. (𝑎 * (1 - (𝑉𝑐[:,1]/𝑉𝑐)^𝑏) + 𝑐) * 4359.74417    # 2D array of dimension nos * length(𝑓)
end


function average𝑝()
    𝑝 = calculate𝑝()         # 2D array of dimension nos*length(𝑓)
    return mean(𝑝, dims=1)   # 1D array of length(𝑓)
    # use the mean function in the Staticstics package
end


function calculate𝐺𝑐𝑎𝑣(𝐸𝑐𝑎𝑣 = 𝐸𝑐𝑎𝑣, 𝑉𝑐 = 𝑉𝑐)
    #𝐸𝑐𝑎𝑣 = get𝐸𝑐𝑎𝑣()    # 2D array of dimension nos*length(𝑓)
    𝑝̄ = average𝑝()      # 1D array of length(𝑓)
    #𝑉𝑐 = get𝑉𝑐()        # 2D array of dimension nos*length(𝑓)
    return @. 𝐸𝑐𝑎𝑣 + 𝑝̄ * 𝑉𝑐 * 2.293712569e-4    # 2D array of dimension nos*length(𝑓)
    # 2.293712569e-4 is the conversion factor from GPa*Å³ to Hartree
end

function calculate𝐺𝑡𝑜𝑡(𝐺𝑒𝑟 = 𝐺𝑒𝑟)
    return 𝐺𝑒𝑟 .+ calculate𝐺𝑐𝑎𝑣()    # 2D array of dimension nos*length(𝑓)
end


# assuming job stopped during electronic energy calculation jobs, i.e., Ger jobs
# The idea is to first generate a set of input filenames, and remove those finished
# from the set. Then restart Ger jobs for the rest of filenames in the set.
function restartger(geom=geometries, 𝑓 = scalingfactors, multi = multithreading)
    nos = numberofstructures(geom)
    a = length(𝑓)
    all = [1:nos;]    # Int64 array containing all job numbers
    
    cd("tmp")

    script = raw"for file in *Ger.log; do i=${file#*structure-}; i=${i%-Ger.log}; grep 'SCF Done' $file | wc -l; echo $i; done | paste - - | awk '/" * "$a" * raw"/ {print $2}'"
    open("restart.sh", "w") do file
        write(file, "$script")
    end

    string = read(`bash restart.sh`, String)
    finished = parse.(Int64, split(string))    # parse the string into an Int64 array
    unfinished = setdiff(all, finished)    # remove the finished job numbers from all
    
    if multi == "off"
        for i in unfinished
            run(`g16 structure-$i-Ger.gjf`)
        end
    elseif multi == "on"
        Threads.@threads for i in unfinished
            run(`g16 structure-$i-Ger.gjf`)
        end
    end

    cd("..")
end

#-------------------------------------
# main program
#-------------------------------------
if restart == "no"
    writegjf("Vc")
    rungaussian("Vc")
    const 𝑉𝑐 = get𝑉𝑐()
    writegjf("Ger")          # write .gjf files for cavity volume "Ger" calculation 
    rungaussian("Ger")
elseif restart == "yes"
    const 𝑉𝑐 = get𝑉𝑐()
    restartger()
else
    println("restart only accepts \"yes\" or \"no\"")
end

#------------------------------------------------------------------------------
# Step 1: cavity volume 𝑉𝑐(𝑓) Gaussian jobs and solvent property calculations
#------------------------------------------------------------------------------

#writegjf("Vc")          # write .gjf files for cavity volume "Vc" calculation 

#rungaussian("Vc")       # run Gaussian jobs

#const 𝑉𝑐 = get𝑉𝑐()       # extract cavity volume 𝑉𝑐 from Gaussian output

#calculate𝑠()            # calculate linear scaling 𝑠 from 𝑉𝑐 date

#average𝑠()              # calculated the average of 𝑠 over all structures

#calculate𝜀()            # calculate dielectric permitivity 𝜀

#calculate𝑍()            # calculate Pauli repulsion barrier 𝑍

#calculate𝑉ₘ()           # calculate the molar volume of solvent 𝑉ₘ

#------------------------------------------------------------------------------
# Step 2: electronic structure Gaussian jobs and pressure calculations 
#------------------------------------------------------------------------------

#    writegjf("Ger")      # write .gjf files for cavity volume "Ger" calculation

#    rungaussian("Ger")   # run Gaussian jobs

const 𝐺𝑒𝑟 = get𝐺𝑒𝑟()      # extract 𝐺𝑒𝑟 from Gaussian output

#calculate𝑝()             # calculate pressure 𝑝

#------------------------------------------------------------------------------
# Step 3: cavitation energy Gaussian jobs
#------------------------------------------------------------------------------

writegjf("Gcav")         # write .gjf files for cavitation energy "Gcav" calculation 

rungaussian("Gcav")      # run Gaussian jobs

const 𝐸𝑐𝑎𝑣 = get𝐸𝑐𝑎𝑣()    # extract 𝐺𝑐𝑎𝑣 from Gaussian output

#calculate𝐺𝑐𝑎𝑣()           # calculate cavitation energy 𝐺𝑐𝑎𝑣

#calculate𝐺𝑡𝑜𝑡()            # calculate total energy 𝐺𝑡𝑜𝑡

#------------------------------------------------------------------------------
# print results
#------------------------------------------------------------------------------

writeproperties()        # write properties.dat file
