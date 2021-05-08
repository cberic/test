# Downloading and installing Julia

https://julialang.org/downloads/

Use the latest version; choose the corresponding version for Linux/Mac/Windows. For Linux, run the following commands to download and extract the julia package:

``` bash
wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.1-linux-x86_64.tar.gz

tar -xzvf julia-1.6.1-linux-x86_64.tar.gz
```

Add julia to PATH (replace "pathtojulia" in the following command with the actual path on your machine):

``` bash
export PATH=/pathtojulia/julia-1.6.1/bin:$PATH
```

# Installing the LsqFit package

This package will be used for the least-square fitting of the Ger-V curve to compute the numberial pressure.
In your terminal, type `julia`, which will open the julia REPL, and in the PEPL, type

``` julia
using Pkg
Pkg.add("LsqFit")
```

The installation of the LsqFit package may take a while. After installing the LsqFit package, run the following command to test if it is working:

``` julia
using LsqFit  # this may take >1 minute to load

xdata = [ 1.0
        0.9410359787896052
        0.892747272010167
        0.8440028046803102
        0.7944607563872212
        0.768315877119944
        0.7444147421008809]

ydata = [ 0.0
        0.002892932000008841
        0.007156206999979986
        0.012909014999991086
        0.020529311000018424
        0.025222962999976062
        0.03096462100000963]

p0 = [0.0, 5.0, 0.0]

@. model(x, p) = (p[1]/p[2])*x^(-p[2]) + (p[1]-p[3])*x

fit = curve_fit(model, xdata, ydata, p0)

print(fit.param)
```

The last command should print an array of 3 floating numbers:

``` julia
[0.04396688697055344, 4.864867834353412, 0.0533354593802239]
```

You may exit the julia REPL using

``` julia
exit()
```

# Running XP-PCM calculations

Copy `xppcm-atom-analyticalp.jl` (the name may be different depending the specific version of the script) and `input.jl` to a working folder. The input file `input.jl` contains Julia codes that specify the calculation parameters. The script `xppcm-atom-analyticalp.jl` will read those parameters before running the XP-PCM calculation. Modify `input.jl` according to the instructions in the file. 

Make sure Gaussian09 or Gaussian16 is properly installed and the exact command `g09` or `g16` will actually call the program. The script will prefer `g16` over `g09` if both are installed. 

Then use the following command (4 cpu cores/threads are assigned) to run an XP-PCM calculation:

``` bash
julia --threads 4 xppcm-atom-analyticalp.jl
```

Below is an example PBS script for running XP-PCM calculations on a cluster. Modify the PBS script according your cluster specifications. Important thing is to let the computing node know the paths to `g09`/`g16` and `julia`.

``` pbs
#!/bin/bash
#PBS -q parallel
#PBS -l nodes=1:ppn=4
#PBS -l mem=8gb
#PBS -l cput=24:00:00 
#PBS -N xppcm

cd $PBS_O_WORKDIR

module load Gaussian/16

/scratch/user/julia-1.5.1/bin/julia --threads 4 xppcm-atom-analyticalp.jl
```

The output will be printed to the `properties.dat` file when the calculation is finished.

# The input file

Here we use the Hydrogen atom as an example. First, give the xyz coordinate in `input.jl`, as one would do for a normal Gaussian calculation, but here in the geometry block with triple quotation marks.

``` julia
# Keep the coordinates (in Angstrom) within the triple """ block.
# Do not include comments or other text in the """ block.
# Blank lines and spaces are ok.
geometries = """

H 0.0 0.0 0.0

"""
```

Then, modify the Gaussian parameters as needed, such as the functional/basis set, charge, and multiplicity. Custom basis set can be used by specifiying the "gen" keyword. The custom basis set should be provided by a separate file named `gen`.

``` julia
# Gaussian 09/16 parameters
# If "gen" is used, provide the custom basis set in a separate file named "gen"
nproc = 4
mem = "4gb"
keywords = "pbe1pbe scf=(Damp,fermi,conver=6) int=finegrid gen"
charge = 0
multiplicity = 2
```

In the `gen` file:
```
H 0
6-31g*
****
```

Finally, modify the XP-PCM parameters as necessary. If not sure, use the default.

``` julia
# Only "cyclohexane", "benzene", and "argon" are currently implemented
solvent = "cyclohexane"

# Change the dielectric permittivity of the solvent to be close to 1 for
# calculations on charged systems. Leave the below commented to use the default value.
#dielectric = 1.0025

# Choose either the Bondi or Rahm vdW radii and the corresponding set of 
# scaling factors
#radiustype = "bondi"
#scalingfactors = [1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]
radiustype = "rahm"
scalingfactors = [1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95, 0.925, 0.90]

# The mean area in Ang^2 of the tesserae by which the surfaces of the cavity
# is partitioned. Suggested value = 0.075.
tesserae = 0.075

# Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
# Larger 𝜂 leads to higher calculated pressures
𝜂 = 6
```

# Run the calculations for multiple atoms

Right now the input filename `input.jl` is hard-coded. If the filename is changed, the script will not work. 

So if we want to calculate multiple atoms, we may (1) use separate folders for each atom; copy the script and input file to each folder, and run the calculation separately, or (2) calculate the first atom, rename its output file, and then modify the input and calculate the second atom in the same folder.

The second approach is illustrated below.

``` julia
# first calculation for the H atom
julia xppcm-atom-analyticalp.jl
mv properties.dat H.dat
# modify the charge, multiplicity, and custom basis set etc. as necessary in input.jl for the He atom
julia xppcm-atom-analyticalp.jl
mv properties.dat He.dat
```

# Understand the output

The script and the XP-PCM calculation will generate several intermediate files, which will be deleted after the job finishes. The remainig files are Gaussian input files `Vc.gjf` and `Ger.gjf`, corresponding output files `Vc.log` and `Ger.log`, and checkpoint files of the Ger calculations for each scaling factor. The file `properties.dat` contains  xppcm results in a tabular form and the orbital energies extracted from `Ger.log` without additional munipulation.

- 𝑓 is the scaling factor of the vdW radii
- 𝑉𝑐 is the cavity volume
- 𝑠 is the cubic root of the volume ratio
- 𝜀(𝑠) is the dielectric permitivity of the solvent
- 𝑍(𝑠) is the density of the solvent
- 𝐺𝑒𝑟 is the SCRF energy in a.u. (hartree)
- 𝑝 is the pressure, numerical and analytical, in GPa
```
#     𝑓         𝑉𝑐(𝑓) Å³   𝑠(𝑓)        𝜀(𝑠)        𝑍(𝑠)        𝐺𝑒𝑟(𝑓) a.u.    𝑝(𝑓)-numeric. -analyt.GPa
 1    1.300      33.611    1.000000    2.016500     0.778100    -0.49886581     0.424     0.183
 2    1.250      29.880    0.961538    2.143426     1.107485    -0.49863442     0.666     0.386
 3    1.200      26.436    0.923077    2.292388     1.599162    -0.49818307     1.130     0.817
 4    1.150      23.267    0.884612    2.468414     2.345625    -0.49730137     2.042     1.739
 5    1.100      20.362    0.846148    2.677909     3.499600    -0.49557670     3.889     3.720
 6    1.050      17.710    0.807692    2.929169     5.318927    -0.49220341     7.752     7.980
 7    1.000      15.299    0.769238    3.233191     8.250696    -0.48563338    16.123    17.064
 8    0.975      14.180    0.750006    3.409420    10.362177    -0.48032872    23.624    24.828
 9    0.950      13.117    0.730776    3.604680    13.091101    -0.47302176    34.987    35.870
10    0.925      12.108    0.711536    3.821736    16.644183    -0.46305764    52.403    51.245
11    0.900      11.153    0.692314    4.063353    21.296322    -0.44969526    79.364    71.943

𝑓 = 1.300    𝑝 =  0.183 GPa ----orbital energies in a.u.----
 Alpha  occ. eigenvalues --   -0.32929
 Alpha virt. eigenvalues --    0.72788
  Beta virt. eigenvalues --    0.00949   0.99577

𝑓 = 1.250    𝑝 =  0.386 GPa ----orbital energies in a.u.----
 Alpha  occ. eigenvalues --   -0.32889
 Alpha virt. eigenvalues --    0.72873
  Beta virt. eigenvalues --    0.01009   0.99665

  ...

```
