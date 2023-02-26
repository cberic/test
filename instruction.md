# Download and install Julia

https://julialang.org/downloads/

Use the latest version. Run the following commands to download and extract the julia package on Linux:

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

Copy `xppcm.jl` and `input.jl` to a working folder. Modify `input.jl` according to the instructions in the file. 

Make sure Gaussian09 or Gaussian16 is properly installed, and `g09` or `g16` will actually call the program. The script will use `g16` over `g09` if both are installed. 

Then run the following command (4 cpu cores/threads are assigned) to start the XP-PCM calculation:

``` bash
julia --threads 4 xppcm.jl
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

/scratch/user/julia-1.8.5/bin/julia --threads 24 xppcm-test.jl
```

# Understand the output

When the calculation starts, a `tmp` folder will be created in the working directory and all Gaussian jobs will be run in that folder. Gaussian input and output files for different structures will be numbered sequentially according to the order of their coordinates given in the `input.jl` file. For each structure, three Gaussian calcuations will be performed sequentially:
- cavity volume (Vc)
- XP-PCM electronic energy (Ger)
- cavitation energy (Gcav)

When the calculation is done, the XP-PCM data will be printed to the `properties.dat` file in the working directory. Below shows the structure of a working folder

```
├── properties.dat             # XP-PCM output
├── xppcm.jl                   # XP-PCM script
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
#      𝑓       𝑉𝑐       𝑠       𝑠̄       𝜀     𝜌ₛₒₗ       𝑉ₘ            𝐸𝑔𝑎𝑠  𝐸ₑₗₑₛₜₐₜ    𝐸ₚₐᵤₗᵢ             𝐺𝑒𝑟        𝑝        𝑝̄     𝑉𝑐𝑎𝑣        𝑝̄𝑉𝑐𝑎𝑣       𝐸𝑐𝑎𝑣          𝐺𝑐𝑎𝑣            𝐺𝑡𝑜𝑡     Δ𝐺𝑡𝑜𝑡
               Å³                                                        Eₕ  kcal/mol  kcal/mol              Eₕ      GPa      GPa       Å³           Eₕ          Eₕ           Eₕ              Eₕ  kcal/mol
1  1.200   17.224  1.0000  1.0000  2.0165   0.7781  108.160     -1.16000500     -0.06      0.90     -1.15866013    1.504    1.468   17.224   0.00579789  0.00549992   0.01129781     -1.14736232      0.00
2  1.150   15.320  0.9617  0.9617  2.1427   0.9833   96.214     -1.15999000     -0.07      1.37     -1.15792596    2.293    2.246   17.224   0.00887175  0.00769465   0.01656640     -1.14135956      0.00
3  1.100   13.569  0.9236  0.9235  2.2905   1.2541   85.196     -1.15995600     -0.09      2.06     -1.15681496    3.682    3.626   17.224   0.01432338  0.01192676   0.02625014     -1.13056482      0.00
4  1.050   11.963  0.8856  0.8855  2.4639   1.6137   75.107     -1.15988600     -0.12      3.23     -1.15492812    6.183    6.116   17.224   0.02416326  0.02234271   0.04650598     -1.10842215      0.00
5  1.000   10.472  0.8472  0.8470  2.6727   2.1069   65.729     -1.15973400     -0.16      4.92     -1.15214245   10.887   10.831   17.224   0.04278991  0.06819110   0.11098101     -1.04116144      0.00
6  0.975    9.787  0.8283  0.8280  2.7905   2.4142   61.404     -1.15960500     -0.18      6.14     -1.15011080   14.582   14.566   17.224   0.05754589  0.20139153   0.25893742     -0.89117338      0.00
7  0.950    9.107  0.8086  0.8086  2.9228   2.7841   57.180     -1.15941200     -0.21      7.66     -1.14753072   19.948   19.909   17.224   0.07865309  5.59331140   5.67196449      4.52443377      0.00

structure 2
#      𝑓       𝑉𝑐       𝑠       𝑠̄       𝜀     𝜌ₛₒₗ       𝑉ₘ            𝐸𝑔𝑎𝑠  𝐸ₑₗₑₛₜₐₜ    𝐸ₚₐᵤₗᵢ             𝐺𝑒𝑟        𝑝        𝑝̄     𝑉𝑐𝑎𝑣        𝑝̄𝑉𝑐𝑎𝑣       𝐸𝑐𝑎𝑣          𝐺𝑐𝑎𝑣            𝐺𝑡𝑜𝑡     Δ𝐺𝑡𝑜𝑡
               Å³                                                        Eₕ  kcal/mol  kcal/mol              Eₕ      GPa      GPa       Å³           Eₕ          Eₕ           Eₕ              Eₕ  kcal/mol
1  1.200   17.163  1.0000  1.0000  2.0165   0.7781  108.160     -1.15954400     -0.06      0.88     -1.15822223    1.471    1.468   17.163   0.00577736  0.00548472   0.01126208     -1.14696015      0.25
2  1.150   15.268  0.9618  0.9617  2.1427   0.9833   96.214     -1.15952800     -0.07      1.34     -1.15750042    2.249    2.246   17.163   0.00884033  0.00767339   0.01651372     -1.14098670      0.23
3  1.100   13.519  0.9235  0.9235  2.2905   1.2541   85.196     -1.15949500     -0.09      2.03     -1.15639905    3.629    3.626   17.163   0.01427265  0.01189381   0.02616647     -1.13023258      0.21
4  1.050   11.918  0.8855  0.8855  2.4639   1.6137   75.107     -1.15942600     -0.12      3.18     -1.15454794    6.118    6.116   17.163   0.02407769  0.02228099   0.04635868     -1.10818926      0.15
5  1.000   10.429  0.8470  0.8470  2.6727   2.1069   65.729     -1.15927600     -0.15      4.85     -1.15178992   10.831   10.831   17.163   0.04263837  0.06800272   0.11064109     -1.04114882      0.01
6  0.975    9.743  0.8280  0.8280  2.7905   2.4142   61.404     -1.15914800     -0.18      6.05     -1.14978486   14.560   14.566   17.163   0.05734209  0.20083520   0.25817728     -0.89160757     -0.27
7  0.950    9.074  0.8086  0.8086  2.9228   2.7841   57.180     -1.15895800     -0.20      7.57     -1.14721997   19.881   19.909   17.163   0.07837453  5.57786027   5.65623480      4.50901483     -9.68

structure 3
#      𝑓       𝑉𝑐       𝑠       𝑠̄       𝜀     𝜌ₛₒₗ       𝑉ₘ            𝐸𝑔𝑎𝑠  𝐸ₑₗₑₛₜₐₜ    𝐸ₚₐᵤₗᵢ             𝐺𝑒𝑟        𝑝        𝑝̄     𝑉𝑐𝑎𝑣        𝑝̄𝑉𝑐𝑎𝑣       𝐸𝑐𝑎𝑣          𝐺𝑐𝑎𝑣            𝐺𝑡𝑜𝑡     Δ𝐺𝑡𝑜𝑡
               Å³                                                        Eₕ  kcal/mol  kcal/mol              Eₕ      GPa      GPa       Å³           Eₕ          Eₕ           Eₕ              Eₕ  kcal/mol
1  1.200   17.102  1.0000  1.0000  2.0165   0.7781  108.160     -1.15892400     -0.05      0.87     -1.15763068    1.428    1.468   17.102   0.00575683  0.00546953   0.01122636     -1.14640432      0.60
2  1.150   15.214  0.9618  0.9617  2.1427   0.9833   96.214     -1.15890900     -0.07      1.33     -1.15690310    2.195    2.246   17.102   0.00880891  0.00765214   0.01646105     -1.14044205      0.58
3  1.100   13.469  0.9235  0.9235  2.2905   1.2541   85.196     -1.15887600     -0.09      2.00     -1.15582517    3.565    3.626   17.102   0.01422193  0.01186087   0.02608279     -1.12974238      0.52
4  1.050   11.873  0.8855  0.8855  2.4639   1.6137   75.107     -1.15880700     -0.11      3.10     -1.15405558    6.048    6.116   17.102   0.02399211  0.02221927   0.04621138     -1.10784419      0.36
5  1.000   10.389  0.8469  0.8470  2.6727   2.1069   65.729     -1.15866000     -0.15      4.78     -1.15127780   10.774   10.831   17.102   0.04248683  0.06781435   0.11030118     -1.04097662      0.12
6  0.975    9.701  0.8278  0.8280  2.7905   2.4142   61.404     -1.15853300     -0.17      5.97     -1.14929910   14.557   14.566   17.102   0.05713828  0.20027887   0.25741715     -0.89188195     -0.44
7  0.950    9.039  0.8085  0.8086  2.9228   2.7841   57.180     -1.15834500     -0.20      7.47     -1.14675159   19.897   19.909   17.102   0.07809598  5.56240913   5.64050511      4.49375352    -19.25
```
