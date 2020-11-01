# Downloading and installing Julia

https://julialang.org/downloads/

Use the latest version. Run the following commands to download and extract the julia package on Linux:

``` bash
wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.2-linux-x86_64.tar.gz

tar -xzvf julia-1.5.2-linux-x86_64.tar.gz
```

Add julia to PATH (replace "pathtojulia" in the following command with the actual path on your machine):

``` bash
export PATH=/pathtojulia/julia-1.5.2/bin:$PATH
```

# Installing the LsqFit package

In your terminal, type `julia`, which will open the julia REPL, and and in the PEPL, type

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

Copy `xppcm-molecule.jl` and `input.jl` to a working folder. The input file `input.jl`, in which the calculation parameters are specified, is in the format of the Julia programming language. The script `xppcm-molecule.jl` will read those parameters before running the XP-PCM calculation. Modify `input.jl` according to the instructions in the file. 

Make sure Gaussian09 or Gaussian16 is properly installed, and the command `g09` or `g16` will actually call the program. The script will use `g16` over `g09` if both are installed. 

Then use the following command (4 cpu cores/threads are assigned) to run an XP-PCM calculation:

``` bash
julia --threads 4 xppcm-molecule.jl
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

/scratch/user/julia-1.5.1/bin/julia --threads 4 xppcm-molecule.jl
```

The output will be printed to the `properties.dat` file when the calculation is finished.

# Calculating IP and EA at high pressure

## Prepare the input file

Let's use the Helium atom as an example. First, give the xyz coordinate in `input.jl`, as one would do for a normal Gaussian calculation, but here in the geometry block with triple quotation marks.

``` julia
# Keep the coordinates (in Angstrom) within the triple """ block.
# Do not include comments or other text in the """ block.
# Blank lines and spaces are ok.
geometries = """

He 0.0 0.0 0.0

"""
```

Then, modify the Gaussian parameters as needed, such as the functional/basis set, charge, and multiplicity

``` julia
# Gaussian 09/16 parameters
nproc = 1
mem = "1gb"
keywords = "b3lyp 6-31g* int=finegrid"    # Gaussian keywords; add more if needed
charge = 0
multiplicity = 1
```

Finally, modify the XP-PCM parameters as necessary. If not sure, use the default.

``` julia
# Only "cyclohexane", "benzene", and "argon" are currently implemented
solvent = "cyclohexane"

# The default dielectric permittivity of the solvent may be set close to 1
# for calculations on charged systems. Comment it out to use the default value.
dielectric = 1.002

# Choose either the Bondi or Rahm vdW radii and the corresponding set of 
# scaling factors
#radiustype = "bondi"
#scalingfactors = [1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]
radiustype = "rahm"
scalingfactors = [1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 1.0, 0.975, 0.95]

# The mean area in Ang^2 of the tesserae by which the surfaces of the cavity
# is partitioned. Suggested value = 0.075.
tesserae = 0.075

# Empirical Pauli repulsion parameter; recommended values: 3, 6 or 9
# Larger 𝜂 leads to higher calculated pressures
𝜂 = 6
```

## Run the calculations for the neutral and cation

Right now the input filename `input.jl` is hard-coded. If the filename is changed, the script will not work. 

So if we want to calculate both the neutral and the cation, we may 1) use separate folders for each species and copy the script and input file to each folder, or 2) calculate the neutral first, rename its output file, and then calculate the cation.

I will use the second approach here.

``` julia
# first calculation for the neutral
julia xppcm-molecule.jl
mv properties.dat He-neutral.dat
# modify the change and multiplicity in input.jl
julia xppcm-molecule.jl
mv properties.dat He-cation.dat
```

## Understanding the output

The XP-PCM calculation will generate two Gaussian input files (`Vc.gjf` and `Ger.gjf`) and correspondingly two Gaussian output files (`Vc.log` and `Ger.log`). A third output file is the `properties.dat` file that contains the xppcm results in a tabular form, and the orbital energies extracted from `Ger.log` without additional munipulation.

- 𝑓 is the scaling factor for the vdW radii
- 𝐺𝑒𝑟 is the energy in a.u. (or hartree)
- 𝑝 is the pressure in GPa

```
#    𝑓       𝑉𝑐(𝑓) Å³   𝑠(𝑓)         𝜀(𝑠)        𝑍(𝑠)        𝐺𝑒𝑟(𝑓) a.u.     𝑝(𝑓) GPa
1    1.30     22.143    1.000000    1.002000     0.778100    -2.90699560     0.143
2    1.25     19.685    0.961538    1.002250     1.107487    -2.90692950     0.260
3    1.20     17.416    0.923076    1.002543     1.599188    -2.90678300     0.515
4    1.15     15.328    0.884606    1.002889     2.345786    -2.90645949     1.091
5    1.10     13.415    0.846158    1.003301     3.499232    -2.90574716     2.442
6    1.05     11.667    0.807683    1.003796     5.319456    -2.90417935     5.739
7    1.00     10.079    0.769237    1.004394     8.250732    -2.90073320    14.119
8    0.97      9.342    0.750011    1.004741    10.361561    -2.89768272    22.548
9    0.95      8.641    0.730762    1.005125    13.093410    -2.89316146    36.486

𝑓 = 1.30    𝑝 =  0.143 GPa ----orbital energies in a.u.----
 Alpha  occ. eigenvalues --   -0.64975
 Alpha virt. eigenvalues --    1.11994

𝑓 = 1.25    𝑝 =  0.260 GPa ----orbital energies in a.u.----
 Alpha  occ. eigenvalues --   -0.64969
 Alpha virt. eigenvalues --    1.12015

𝑓 = 1.20    𝑝 =  0.515 GPa ----orbital energies in a.u.----
 Alpha  occ. eigenvalues --   -0.64954
 Alpha virt. eigenvalues --    1.12060

𝑓 = 1.15    𝑝 =  1.091 GPa ----orbital energies in a.u.----
 Alpha  occ. eigenvalues --   -0.64921
 Alpha virt. eigenvalues --    1.12160

𝑓 = 1.10    𝑝 =  2.442 GPa ----orbital energies in a.u.----
 Alpha  occ. eigenvalues --   -0.64848
 Alpha virt. eigenvalues --    1.12379

𝑓 = 1.05    𝑝 =  5.739 GPa ----orbital energies in a.u.----
 Alpha  occ. eigenvalues --   -0.64689
 Alpha virt. eigenvalues --    1.12863

𝑓 = 1.00    𝑝 = 14.119 GPa ----orbital energies in a.u.----
 Alpha  occ. eigenvalues --   -0.64338
 Alpha virt. eigenvalues --    1.13926

𝑓 = 0.97    𝑝 = 22.548 GPa ----orbital energies in a.u.----
 Alpha  occ. eigenvalues --   -0.64028
 Alpha virt. eigenvalues --    1.14870

𝑓 = 0.95    𝑝 = 36.486 GPa ----orbital energies in a.u.----
 Alpha  occ. eigenvalues --   -0.63567
 Alpha virt. eigenvalues --    1.16275
```

Let's compare the XP-PCM results of the neutral and the cation. The scaling factors 𝑓 and 𝑠, volume 𝑉𝑐, dielectric constant 𝜀, Pauli repulsion barrier 𝑍 are the same for both species, because these parameters are based on the same atomic radius, scaling factors and solvent properties.

The energy 𝐺𝑒𝑟 and pressure 𝑝 are different for the two species at the same scaling factor. The neutral helium atom has one more electron than the cation, and its electron density is more diffuse than the cation, which leads to a stronger repulsion exerted by the solvent continuum and thus a higher pressure for the neutral.

```
# neutral Helium
#    𝑓       𝑉𝑐(𝑓) Å³   𝑠(𝑓)         𝜀(𝑠)        𝑍(𝑠)        𝐺𝑒𝑟(𝑓) a.u.     𝑝(𝑓) GPa
1    1.30     22.143    1.000000    1.002000     0.778100    -2.90699560     0.143
2    1.25     19.685    0.961538    1.002250     1.107487    -2.90692950     0.260
3    1.20     17.416    0.923076    1.002543     1.599188    -2.90678300     0.515
4    1.15     15.328    0.884606    1.002889     2.345786    -2.90645949     1.091
5    1.10     13.415    0.846158    1.003301     3.499232    -2.90574716     2.442
6    1.05     11.667    0.807683    1.003796     5.319456    -2.90417935     5.739
7    1.00     10.079    0.769237    1.004394     8.250732    -2.90073320    14.119
8    0.97      9.342    0.750011    1.004741    10.361561    -2.89768272    22.548
9    0.95      8.641    0.730762    1.005125    13.093410    -2.89316146    36.486

# cation Helium
#    𝑓       𝑉𝑐(𝑓) Å³   𝑠(𝑓)         𝜀(𝑠)        𝑍(𝑠)        𝐺𝑒𝑟(𝑓) a.u.     𝑝(𝑓) GPa
1    1.30     22.143    1.000000    1.002000     0.778100    -1.99345028     0.020
2    1.25     19.685    0.961538    1.002250     1.107487    -1.99348860     0.026
3    1.20     17.416    0.923076    1.002543     1.599188    -1.99352264     0.041
4    1.15     15.328    0.884606    1.002889     2.345786    -1.99353653     0.086
5    1.10     13.415    0.846158    1.003301     3.499232    -1.99349310     0.221
6    1.05     11.667    0.807683    1.003796     5.319456    -1.99330558     0.652
7    1.00     10.079    0.769237    1.004394     8.250732    -1.99278113     2.103
8    0.97      9.342    0.750011    1.004741    10.361561    -1.99227228     3.883
9    0.95      8.641    0.730762    1.005125    13.093410    -1.99148485     7.301
```

The vertical IP can be calculated by taking the energy difference between the neutral and the cation Ger(cation)-Ger(neutral) at the same scaling factor.

```
# 𝑝(𝑓) GPa    IP in a.u.    IP in eV
 0.143       0.91354532     24.858
 0.260       0.91344089     24.856
 0.515       0.91326036     24.851
 1.091       0.91292296     24.841
 2.442       0.91225406     24.823
 5.739       0.91087377     24.786
14.119       0.90795206     24.706
22.548       0.90541044     24.637
36.486       0.90167661     24.535
```

IP decreases as the pressure increases. The reason is that at high pressure, the valence electrons in the neutral He atom experiences a stronger repulsion than at low pressure; consequently it is easier to remove one valence electron from He at high pressure than at low pressure.

Similarly, vertical EAs can be calculated by taking the energy differences between the neutral and the anion Ger(neutral)-Ger(anion) at the same scaling factors. Negative EA indicates destabilization when one electron is added to the neutral. Due to the lack of diffuse functions in the basis set, the calculated energies of the anions are not reliable, so are the EA values.


```
# 𝑝(𝑓) GPa   IP in eV     EA in eV
 0.143       24.858       -11.124
 0.260       24.856       -11.131
 0.515       24.851       -11.146
 1.091       24.841       -11.180
 2.442       24.823       -11.255
 5.739       24.786       -11.419
14.119       24.706       -11.781
22.548       24.637       -12.102
36.486       24.535       -12.579
```

# Testing calculations for H to Ar

Changes of Mulliken electronegativities of atoms at 50 GPa, calculated at the pbe1pbe/aug-cc-pvtz level. The energy units are in eV.

![](~/Dropbox/work/DIPC/Conceptual-DFT/electronegativity.png)