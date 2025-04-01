# Download and install Julia

https://julialang.org/downloads/

Use the latest version. For example, run the following commands to download and extract the julia package on Linux:

``` bash
wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.5-linux-x86_64.tar.gz

tar -xzvf julia-1.8.5-linux-x86_64.tar.gz
```

Add julia to PATH (replace "pathtojulia" in the following command with the actual path on your machine):

``` bash
export PATH=/pathtojulia/julia-1.8.5/bin:$PATH
```

# Install the LsqFit.jl package

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

# Run XP-PCM calculations

Copy `xppcm-sp.jl` and `input.jl` to a working folder. Modify `input.jl` according to the instructions in the file. 

Make sure Gaussian09 or Gaussian16 is properly installed, and `g09` or `g16` will actually call the program. The script will use `g16` over `g09` if both are installed. 

Then run the following command (4 cpu cores/threads are assigned) to start the XP-PCM calculation:

``` bash
julia --threads 4 xppcm-sp.jl
```

Below is an example PBS script for running XP-PCM calculations on a cluster. Modify the PBS script according your cluster specifications. Important thing is to let the computing node know the paths to `g09`/`g16` and `julia`.

``` pbs
#!/bin/bash
#PBS -q parallel
#PBS -l nodes=1:ppn=24
#PBS -l mem=48gb
#PBS -l cput=24:00:00 
#PBS -N xppcm

cd $PBS_O_WORKDIR

module load Gaussian/16

/scratch/user/julia-1.8.5/bin/julia --threads 24 xppcm-sp.jl
```

# Understand the output

When the calculation starts, a `tmp` folder will be created in the working directory and all Gaussian jobs will be run in that folder. Gaussian input and output files for different structures will be numbered sequentially according to the order of their coordinates given in the `input.jl` file. For each structure, three Gaussian calcuations will be performed sequentially:
- cavity volume (Vc)
- XP-PCM electronic energy (Ger)
- cavitation energy (Gcav)

When the calculation is done, the XP-PCM data will be printed to the `properties.dat` file in the working directory. Below shows the structure of a working folder

```
├── properties.dat             # XP-PCM output
├── xppcm-sp.jl                # XP-PCM script
├── input.jl                   # XP-PCM input
└── tmp                        # A fold created by the script for Gaussian jobs
    ├── tesserae.off
    ├── structure-3-Vc.log
    ├── structure-3-Vc.gjf
    ├── structure-3-Ger.log
    ├── structure-3-Ger.gjf
    ├── structure-3-Ger.chk
    ├── structure-3-Gcav.log
    ├── structure-3-Gcav.gjf
    ├── structure-2-Vc.log
    ├── structure-2-Vc.gjf
    ├── structure-2-Ger.log
    ├── structure-2-Ger.gjf
    ├── structure-2-Ger.chk
    ├── structure-2-Gcav.log
    ├── structure-2-Gcav.gjf
    ├── structure-1-Vc.log
    ├── structure-1-Vc.gjf
    ├── structure-1-Ger.log
    ├── structure-1-Ger.gjf
    ├── structure-1-Ger.chk
    ├── structure-1-Gcav.log
    ├── structure-1-Gcav.gjf
    └── charge.off
```

Below shows the `properties.dat` file of an XP-PCM calculations on three structures: H2 with the H---H distance of 0.74, 0.73 and 0.72 Å. If some of the symbols are not displayed correctly, consider using the JuliaMono font (https://juliamono.netlify.app) for your text editor.
- 𝑓 is the scaling factor for the vdW radii
- 𝑉𝑐 is the volume of the molecule
- 𝑠 is the linear scaling factor
- 𝑠̄  is the average of 𝑠 over all structures
- 𝜀, 𝜌ₛₒₗ and 𝑉ₘ are solvent properties
- 𝐸𝑔𝑎𝑠 is gasphase energy
- 𝐸ₑₗₑₛₜₐₜ is the electrostatic interaction energy between the molecule and the solvent of the PCM model
- 𝐸ₚₐᵤₗᵢ is the Pauli repulsion energy between the molecule and the solvent
- 𝐺𝑒𝑟 is the total XP-PCM electronic energy, the sum of previous three terms
- 𝑝 is the pressure
- 𝑝̄  is the average pressure over all structures
- 𝑉𝑐𝑎𝑣 is the cavity volume used in cavitation energy calculation
- 𝑝̄𝑉𝑐𝑎𝑣 is perssure times cavity volume in cavitation energy calculation
- 𝐸𝑐𝑎𝑣 is the non-pV term of the cavitation energy
- 𝐺𝑐𝑎𝑣 is the total cavitation energy, the sum of previous two terms
- 𝐺𝑡𝑜𝑡 is the total XP-PCM energy, the sum of 𝐺𝑒𝑟 and 𝐺𝑐𝑎𝑣
- Δ𝐺𝑡𝑜𝑡 is the relative energy of 𝐺𝑡𝑜𝑡. This is usually the quantity of interest for a high-pressure reaction.

```
structure 1
#      𝑓       𝑉𝑐       𝑠       𝑠̄       𝜀     𝜌ₛₒₗ       𝑉ₘ            𝐸𝑔𝑎𝑠  𝐸ₑₗₑₛₜₐₜ    𝐸ₚₐᵤₗᵢ             𝐺𝑒𝑟        𝑝        𝑝̄     𝑉𝑐𝑎𝑣        𝑝̄𝑉𝑐𝑎𝑣     𝐺𝑛𝑜𝑛_𝑝𝑉         𝐺𝑐𝑎𝑣            𝐺𝑡𝑜𝑡     Δ𝐺𝑡𝑜𝑡
               Å³                                                        Eₕ  kcal/mol  kcal/mol              Eₕ      GPa      GPa       Å³           Eₕ          Eₕ           Eₕ              Eₕ  kcal/mol
1  1.200   17.254  1.0000  1.0000  2.0165   0.7781  108.160     -1.17152200     -0.06      0.87     -1.17022012    1.142    1.149   17.254   0.00454688  0.00549992   0.01004680     -1.16017332      0.00
2  1.150   15.365  0.9621  0.9621  2.1416   0.9814   96.310     -1.17150600     -0.07      1.31     -1.16953841    2.097    2.068   17.254   0.00818511  0.00767044   0.01585555     -1.15368286      0.00
3  1.100   13.616  0.9241  0.9240  2.2883   1.2499   85.338     -1.17147300     -0.09      2.00     -1.16843428    3.652    3.586   17.254   0.01419048  0.01184523   0.02603572     -1.14239856      0.00
4  1.050   12.001  0.8860  0.8859  2.4620   1.6096   75.201     -1.17140200     -0.12      3.14     -1.16659183    6.240    6.150   17.254   0.02434089  0.02217362   0.04651451     -1.12007732      0.00
5  1.000   10.520  0.8480  0.8478  2.6679   2.0948   65.919     -1.17125300     -0.16      4.79     -1.16387697   10.628   10.562   17.254   0.04180031  0.06591840   0.10771871     -1.05615826      0.00
6  0.975    9.825  0.8289  0.8287  2.7860   2.4021   61.559     -1.17112300     -0.18      5.98     -1.16188538   13.922   13.914   17.254   0.05506760  0.19072389   0.24579149     -0.91609389      0.00

structure 2
#      𝑓       𝑉𝑐       𝑠       𝑠̄       𝜀     𝜌ₛₒₗ       𝑉ₘ            𝐸𝑔𝑎𝑠  𝐸ₑₗₑₛₜₐₜ    𝐸ₚₐᵤₗᵢ             𝐺𝑒𝑟        𝑝        𝑝̄     𝑉𝑐𝑎𝑣        𝑝̄𝑉𝑐𝑎𝑣     𝐺𝑛𝑜𝑛_𝑝𝑉         𝐺𝑐𝑎𝑣            𝐺𝑡𝑜𝑡     Δ𝐺𝑡𝑜𝑡
               Å³                                                        Eₕ  kcal/mol  kcal/mol              Eₕ      GPa      GPa       Å³           Eₕ          Eₕ           Eₕ              Eₕ  kcal/mol
1  1.200   17.193  1.0000  1.0000  2.0165   0.7781  108.160     -1.17119700     -0.06      0.86     -1.16992012    1.139    1.149   17.193   0.00453081  0.00548472   0.01001553     -1.15990458      0.17
2  1.150   15.309  0.9621  0.9621  2.1416   0.9814   96.310     -1.17118100     -0.07      1.29     -1.16924580    2.072    2.068   17.193   0.00815617  0.00764925   0.01580542     -1.15344037      0.15
3  1.100   13.565  0.9240  0.9240  2.2883   1.2499   85.338     -1.17114800     -0.09      1.97     -1.16815773    3.600    3.586   17.193   0.01414032  0.01181251   0.02595283     -1.14220490      0.12
4  1.050   11.955  0.8859  0.8859  2.4620   1.6096   75.201     -1.17107900     -0.12      3.08     -1.16635239    6.157    6.150   17.193   0.02425483  0.02211237   0.04636720     -1.11998518      0.06
5  1.000   10.478  0.8478  0.8478  2.6679   2.0948   65.919     -1.17093100     -0.16      4.71     -1.16366879   10.523   10.562   17.193   0.04165253  0.06573631   0.10738884     -1.05627995     -0.08
6  0.975    9.785  0.8287  0.8287  2.7860   2.4021   61.559     -1.17080400     -0.18      5.89     -1.16170557   13.816   13.914   17.193   0.05487292  0.19019703   0.24506994     -0.91663563     -0.34

structure 3
#      𝑓       𝑉𝑐       𝑠       𝑠̄       𝜀     𝜌ₛₒₗ       𝑉ₘ            𝐸𝑔𝑎𝑠  𝐸ₑₗₑₛₜₐₜ    𝐸ₚₐᵤₗᵢ             𝐺𝑒𝑟        𝑝        𝑝̄     𝑉𝑐𝑎𝑣        𝑝̄𝑉𝑐𝑎𝑣     𝐺𝑛𝑜𝑛_𝑝𝑉         𝐺𝑐𝑎𝑣            𝐺𝑡𝑜𝑡     Δ𝐺𝑡𝑜𝑡
               Å³                                                        Eₕ  kcal/mol  kcal/mol              Eₕ      GPa      GPa       Å³           Eₕ          Eₕ           Eₕ              Eₕ  kcal/mol
1  1.200   17.130  1.0000  1.0000  2.0165   0.7781  108.160     -1.17071400     -0.05      0.84     -1.16946703    1.166    1.149   17.130   0.00451421  0.00546953   0.00998374     -1.15948329      0.43
2  1.150   15.252  0.9620  0.9621  2.1416   0.9814   96.310     -1.17069900     -0.07      1.26     -1.16879569    2.035    2.068   17.130   0.00812629  0.00762806   0.01575435     -1.15304134      0.40
3  1.100   13.513  0.9240  0.9240  2.2883   1.2499   85.338     -1.17066600     -0.09      1.94     -1.16772370    3.505    3.586   17.130   0.01408850  0.01177979   0.02586829     -1.14185541      0.34
4  1.050   11.904  0.8858  0.8859  2.4620   1.6096   75.201     -1.17059800     -0.12      3.00     -1.16600145    6.055    6.150   17.130   0.02416596  0.02205112   0.04621707     -1.11978437      0.18
5  1.000   10.436  0.8477  0.8478  2.6679   2.0948   65.919     -1.17045200     -0.15      4.64     -1.16330136   10.534   10.562   17.130   0.04149990  0.06555421   0.10705412     -1.05624725     -0.06
6  0.975    9.745  0.8286  0.8287  2.7860   2.4021   61.559     -1.17032600     -0.17      5.80     -1.16136624   14.006   13.914   17.130   0.05467185  0.18967016   0.24434201     -0.91702423     -0.58
```
![alt text](image.png)