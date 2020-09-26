# assuming job stopped during electronic energy calculation jobs, i.e., Ger jobs
# The idea is to first generate a set of input filenames, and remove those finished
# from the set. Then restart Ger jobs for the rest of filenames in the set.
function restart(geom=geometries, 𝑓 = scalingfactors)
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
    
    for i in unfinished
        run(`g16 structure-$i-Ger.gjf`)
    end

    cd("..")
end