# calculate linear scaling factor 𝑠, as the cubic root of volume scaling 𝑉𝑐/𝑉𝑐𝑓=1.2
function calculate𝑠(geom = geometries, 𝑓 = scalingfactors)
    nos = numberofstructures(geom)
    a = length(𝑓)
    𝑉𝑐 = get𝑉𝑐()
    𝑠 = Array{Float64}(undef, nos,a)
    for i in 1:nos
        for j in 1:a
            𝑠[i,j] = cbrt(𝑉𝑐[i,j]/𝑉𝑐[i,1])
        end
    end
    return 𝑠
end


# 𝑠̄ is the average of 𝑠 over all structures at the same scalingfactor 𝑓
using Statistics
function average𝑠()
    𝑠 = calculate𝑠()
    return mean(𝑠, dims=1)   # use the mean function in the Staticstics package
end


# calculate dielectric permitivity 𝜀 as a function of 𝑠̄
function calculate𝜀()
    𝜀₀ = solventparameters()[1]
    𝑠̄ = average𝑠()
    # 𝜀 = 1 + (𝜀₀-1)/𝑠̄³
    return 1 .+ (𝜀₀ .- 1) ./ 𝑠̄.^3    # using vectorized "dot" operator for 1D array
end


# calculate Pauli repulson barrier 𝑍 as a function of 𝑠̄ and 𝜂
function calculate𝑍(𝜂=𝜂)
    𝜌₀ = solventparameters()[2]
    𝑠̄ = average𝑠()
    # 𝑍 = 𝑍₀/𝑠̄⁽³⁺𝜂⁾   𝑍₀ = 𝜌₀
    return 𝜌₀ ./ 𝑠̄.^(3+𝜂)    # using vectorized "dot" operator for 1D array
end


# calculate the molar volume of solvent 𝑉ₘ as a function of 𝑠̄
function calculate𝑉ₘ()
    (𝜌₀, 𝑀) = solventparameters()[2:3]
    𝑠̄ = average𝑠()
    # 𝑉ₘ = (𝑀/𝜌₀) * 𝑠̄³
    return (𝑀/𝜌₀) .* 𝑠̄.^3    # using vectorized "dot" operator for 1D array
end


# Murnaghan equation of state fitting to calculate pressure 𝑝
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
end

# need to install the lmfit package `pip install lmfit`
# make sure to use the correct python verion
# note the python3 call inside the function
function pythonfitting(geom = geometries)
    nos = numberofstructures(geom)
    𝑉𝑐 = get𝑉𝑐()
    𝐺𝑒𝑟 = get𝐺𝑒𝑟()
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
        
        results = read(pipeline(`echo $script`, `python3`), String)
        open("tmp/structure-$i-fitting.out", "w") do file
            write(file, results)
        end

        awkscript = raw"/a:/ {a=$2}; /b:/ {b=$2}; /c:/ {c=$2}; END {print a,b,c}"
        abc = read(`awk $awkscript "tmp/structure-$i-fitting.out"`, String)
        for j = 1:3
            array[i,j] = parse(Float64, split(abc)[j])
        end
    end
    return array
end

# write data to Vc-Ger.dat file which is then read by Mathematica
using DelimitedFiles
function mathematicafitting(geom = geometries)
    nos = numberofstructures(geom)
    𝑉𝑐 = get𝑉𝑐()
    𝐺𝑒𝑟 = get𝐺𝑒𝑟()
    array = Array{Float64}(undef, nos,3)
    Threads.@threads for i in 1:nos
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
        open("tmp/structure-$i-fitting.out", "w") do file
            write(file, results)
        end
    end
    rm("tmp/Vc-Ger.dat")
    return array
    # add feature to extract abc data
    #results = read(pipeline(`echo $script`, `math`), String)
end


function calculate𝑝()
    𝑉𝑐 = get𝑉𝑐()    # 2D array of dimension nos*length(𝑓)
    abc = murnaghan("python")    # 2D array of dimension nos*3 
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
    return @. (𝑎 * (1 - (𝑉𝑐[:,1]/𝑉𝑐)^𝑏) + 𝑐) * 4359.74417    # 2D array of dimension nos*length(𝑓)
end


function average𝑝()
    𝑝 = calculate𝑝()         # 2D array of dimension nos*length(𝑓)
    return mean(𝑝, dims=1)   # 1D array of length(𝑓)
    # use the mean function in the Staticstics package
end


function calculate𝐺𝑐𝑎𝑣()
    𝐸𝑐𝑎𝑣 = get𝐸𝑐𝑎𝑣()    # 2D array of dimension nos*length(𝑓)
    𝑝̄ = average𝑝()      # 1D array of length(𝑓)
    𝑉𝑐 = get𝑉𝑐()        # 2D array of dimension nos*length(𝑓)
    return @. 𝐸𝑐𝑎𝑣 + 𝑝̄ * 𝑉𝑐 * 2.293712569e-4    # 2D array of dimension nos*length(𝑓)
    # 2.293712569e-4 is the conversion factor from GPa*Å³ to Hartree
end

function calculate𝐺𝑡𝑜𝑡()
    return @. get𝐺𝑒𝑟() + calculate𝐺𝑐𝑎𝑣()    # 2D array of dimension nos*length(𝑓)
end