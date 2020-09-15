# multithreading for short Vc and Cav calculations.
# add return value for output file "return "$jobtype finished at time()""
function rungaussian(jobtype; geom=geometries)
    nos = numberofstructures(geom)
    run(`cd tmp`)
    if jobtype == "Vc" || jobtype == "Cav"
        Threads.@threads for i in 1:nos
            run(`g16 structure-$i-$jobtype.gjf`)
        end
    elseif jobtype == "Ger"
        for i in 1:nos
            run(`g16 structure-$i-$jobtype.gjf`)
        end
    end
    run(`cd ..`)
end