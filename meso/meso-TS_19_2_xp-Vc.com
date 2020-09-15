%chk=meso-TS_19_2_xp-12.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120849  -0.437571   0.247983 
  6  -0.595896  -0.213392  -1.011977 
  1  -4.084599  -0.548890   0.731234 
  1   0.369499  -0.121508  -1.489233 
  6   1.095787   1.409118   0.700116 
  6   0.595896   0.213392   1.011977 
  1   0.513469   2.310071   0.845290 
  1  -0.369499   0.121508   1.489233 
  6   1.478063  -1.000287   1.004227 
  1   1.940410  -1.028087   1.998902 
  1   0.895946  -1.915783   0.915280 
  6   2.423651   1.516981   0.110100 
  1   2.832064   2.505394  -0.064213 
  6   3.120849   0.437571  -0.247983 
  1   4.084599   0.548890  -0.731234 
  6   2.590846  -0.955955  -0.045910 
  1   2.243218  -1.349244  -1.009243 
  1   3.400776  -1.620708   0.261854 
  6  -1.095787  -1.409118  -0.700116 
  1  -0.513469  -2.310071  -0.845290 
  6  -2.423651  -1.516981  -0.110100 
  1  -2.832064  -2.505394   0.064213 
  6  -1.478063   1.000287  -1.004227 
  1  -0.895946   1.915783  -0.915280 
  1  -1.940410   1.028087  -1.998902 
  6  -2.590846   0.955955   0.045910 
  1  -3.400776   1.620708  -0.261855 
  1  -2.243218   1.349244   1.009243 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120849    -0.437571     0.247983     1.700000     1.200000
   -0.595896    -0.213392    -1.011977     1.700000     1.200000
   -4.084599    -0.548890     0.731234     1.200000     1.200000
    0.369499    -0.121508    -1.489233     1.200000     1.200000
    1.095787     1.409118     0.700116     1.700000     1.200000
    0.595896     0.213392     1.011977     1.700000     1.200000
    0.513469     2.310071     0.845290     1.200000     1.200000
   -0.369499     0.121508     1.489233     1.200000     1.200000
    1.478063    -1.000287     1.004227     1.700000     1.200000
    1.940410    -1.028087     1.998902     1.200000     1.200000
    0.895946    -1.915783     0.915280     1.200000     1.200000
    2.423651     1.516981     0.110100     1.700000     1.200000
    2.832064     2.505394    -0.064213     1.200000     1.200000
    3.120849     0.437571    -0.247983     1.700000     1.200000
    4.084599     0.548890    -0.731234     1.200000     1.200000
    2.590846    -0.955955    -0.045910     1.700000     1.200000
    2.243218    -1.349244    -1.009243     1.200000     1.200000
    3.400776    -1.620708     0.261854     1.200000     1.200000
   -1.095787    -1.409118    -0.700116     1.700000     1.200000
   -0.513469    -2.310071    -0.845290     1.200000     1.200000
   -2.423651    -1.516981    -0.110100     1.700000     1.200000
   -2.832064    -2.505394     0.064213     1.200000     1.200000
   -1.478063     1.000287    -1.004227     1.700000     1.200000
   -0.895946     1.915783    -0.915280     1.200000     1.200000
   -1.940410     1.028087    -1.998902     1.200000     1.200000
   -2.590846     0.955955     0.045910     1.700000     1.200000
   -3.400776     1.620708    -0.261855     1.200000     1.200000
   -2.243218     1.349244     1.009243     1.200000     1.200000

--link1--
%chk=meso-TS_19_2_xp-115.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120849  -0.437571   0.247983 
  6  -0.595896  -0.213392  -1.011977 
  1  -4.084599  -0.548890   0.731234 
  1   0.369499  -0.121508  -1.489233 
  6   1.095787   1.409118   0.700116 
  6   0.595896   0.213392   1.011977 
  1   0.513469   2.310071   0.845290 
  1  -0.369499   0.121508   1.489233 
  6   1.478063  -1.000287   1.004227 
  1   1.940410  -1.028087   1.998902 
  1   0.895946  -1.915783   0.915280 
  6   2.423651   1.516981   0.110100 
  1   2.832064   2.505394  -0.064213 
  6   3.120849   0.437571  -0.247983 
  1   4.084599   0.548890  -0.731234 
  6   2.590846  -0.955955  -0.045910 
  1   2.243218  -1.349244  -1.009243 
  1   3.400776  -1.620708   0.261854 
  6  -1.095787  -1.409118  -0.700116 
  1  -0.513469  -2.310071  -0.845290 
  6  -2.423651  -1.516981  -0.110100 
  1  -2.832064  -2.505394   0.064213 
  6  -1.478063   1.000287  -1.004227 
  1  -0.895946   1.915783  -0.915280 
  1  -1.940410   1.028087  -1.998902 
  6  -2.590846   0.955955   0.045910 
  1  -3.400776   1.620708  -0.261855 
  1  -2.243218   1.349244   1.009243 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120849    -0.437571     0.247983     1.700000     1.150000
   -0.595896    -0.213392    -1.011977     1.700000     1.150000
   -4.084599    -0.548890     0.731234     1.200000     1.150000
    0.369499    -0.121508    -1.489233     1.200000     1.150000
    1.095787     1.409118     0.700116     1.700000     1.150000
    0.595896     0.213392     1.011977     1.700000     1.150000
    0.513469     2.310071     0.845290     1.200000     1.150000
   -0.369499     0.121508     1.489233     1.200000     1.150000
    1.478063    -1.000287     1.004227     1.700000     1.150000
    1.940410    -1.028087     1.998902     1.200000     1.150000
    0.895946    -1.915783     0.915280     1.200000     1.150000
    2.423651     1.516981     0.110100     1.700000     1.150000
    2.832064     2.505394    -0.064213     1.200000     1.150000
    3.120849     0.437571    -0.247983     1.700000     1.150000
    4.084599     0.548890    -0.731234     1.200000     1.150000
    2.590846    -0.955955    -0.045910     1.700000     1.150000
    2.243218    -1.349244    -1.009243     1.200000     1.150000
    3.400776    -1.620708     0.261854     1.200000     1.150000
   -1.095787    -1.409118    -0.700116     1.700000     1.150000
   -0.513469    -2.310071    -0.845290     1.200000     1.150000
   -2.423651    -1.516981    -0.110100     1.700000     1.150000
   -2.832064    -2.505394     0.064213     1.200000     1.150000
   -1.478063     1.000287    -1.004227     1.700000     1.150000
   -0.895946     1.915783    -0.915280     1.200000     1.150000
   -1.940410     1.028087    -1.998902     1.200000     1.150000
   -2.590846     0.955955     0.045910     1.700000     1.150000
   -3.400776     1.620708    -0.261855     1.200000     1.150000
   -2.243218     1.349244     1.009243     1.200000     1.150000

--link1--
%chk=meso-TS_19_2_xp-110.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120849  -0.437571   0.247983 
  6  -0.595896  -0.213392  -1.011977 
  1  -4.084599  -0.548890   0.731234 
  1   0.369499  -0.121508  -1.489233 
  6   1.095787   1.409118   0.700116 
  6   0.595896   0.213392   1.011977 
  1   0.513469   2.310071   0.845290 
  1  -0.369499   0.121508   1.489233 
  6   1.478063  -1.000287   1.004227 
  1   1.940410  -1.028087   1.998902 
  1   0.895946  -1.915783   0.915280 
  6   2.423651   1.516981   0.110100 
  1   2.832064   2.505394  -0.064213 
  6   3.120849   0.437571  -0.247983 
  1   4.084599   0.548890  -0.731234 
  6   2.590846  -0.955955  -0.045910 
  1   2.243218  -1.349244  -1.009243 
  1   3.400776  -1.620708   0.261854 
  6  -1.095787  -1.409118  -0.700116 
  1  -0.513469  -2.310071  -0.845290 
  6  -2.423651  -1.516981  -0.110100 
  1  -2.832064  -2.505394   0.064213 
  6  -1.478063   1.000287  -1.004227 
  1  -0.895946   1.915783  -0.915280 
  1  -1.940410   1.028087  -1.998902 
  6  -2.590846   0.955955   0.045910 
  1  -3.400776   1.620708  -0.261855 
  1  -2.243218   1.349244   1.009243 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120849    -0.437571     0.247983     1.700000     1.100000
   -0.595896    -0.213392    -1.011977     1.700000     1.100000
   -4.084599    -0.548890     0.731234     1.200000     1.100000
    0.369499    -0.121508    -1.489233     1.200000     1.100000
    1.095787     1.409118     0.700116     1.700000     1.100000
    0.595896     0.213392     1.011977     1.700000     1.100000
    0.513469     2.310071     0.845290     1.200000     1.100000
   -0.369499     0.121508     1.489233     1.200000     1.100000
    1.478063    -1.000287     1.004227     1.700000     1.100000
    1.940410    -1.028087     1.998902     1.200000     1.100000
    0.895946    -1.915783     0.915280     1.200000     1.100000
    2.423651     1.516981     0.110100     1.700000     1.100000
    2.832064     2.505394    -0.064213     1.200000     1.100000
    3.120849     0.437571    -0.247983     1.700000     1.100000
    4.084599     0.548890    -0.731234     1.200000     1.100000
    2.590846    -0.955955    -0.045910     1.700000     1.100000
    2.243218    -1.349244    -1.009243     1.200000     1.100000
    3.400776    -1.620708     0.261854     1.200000     1.100000
   -1.095787    -1.409118    -0.700116     1.700000     1.100000
   -0.513469    -2.310071    -0.845290     1.200000     1.100000
   -2.423651    -1.516981    -0.110100     1.700000     1.100000
   -2.832064    -2.505394     0.064213     1.200000     1.100000
   -1.478063     1.000287    -1.004227     1.700000     1.100000
   -0.895946     1.915783    -0.915280     1.200000     1.100000
   -1.940410     1.028087    -1.998902     1.200000     1.100000
   -2.590846     0.955955     0.045910     1.700000     1.100000
   -3.400776     1.620708    -0.261855     1.200000     1.100000
   -2.243218     1.349244     1.009243     1.200000     1.100000

--link1--
%chk=meso-TS_19_2_xp-105.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120849  -0.437571   0.247983 
  6  -0.595896  -0.213392  -1.011977 
  1  -4.084599  -0.548890   0.731234 
  1   0.369499  -0.121508  -1.489233 
  6   1.095787   1.409118   0.700116 
  6   0.595896   0.213392   1.011977 
  1   0.513469   2.310071   0.845290 
  1  -0.369499   0.121508   1.489233 
  6   1.478063  -1.000287   1.004227 
  1   1.940410  -1.028087   1.998902 
  1   0.895946  -1.915783   0.915280 
  6   2.423651   1.516981   0.110100 
  1   2.832064   2.505394  -0.064213 
  6   3.120849   0.437571  -0.247983 
  1   4.084599   0.548890  -0.731234 
  6   2.590846  -0.955955  -0.045910 
  1   2.243218  -1.349244  -1.009243 
  1   3.400776  -1.620708   0.261854 
  6  -1.095787  -1.409118  -0.700116 
  1  -0.513469  -2.310071  -0.845290 
  6  -2.423651  -1.516981  -0.110100 
  1  -2.832064  -2.505394   0.064213 
  6  -1.478063   1.000287  -1.004227 
  1  -0.895946   1.915783  -0.915280 
  1  -1.940410   1.028087  -1.998902 
  6  -2.590846   0.955955   0.045910 
  1  -3.400776   1.620708  -0.261855 
  1  -2.243218   1.349244   1.009243 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120849    -0.437571     0.247983     1.700000     1.050000
   -0.595896    -0.213392    -1.011977     1.700000     1.050000
   -4.084599    -0.548890     0.731234     1.200000     1.050000
    0.369499    -0.121508    -1.489233     1.200000     1.050000
    1.095787     1.409118     0.700116     1.700000     1.050000
    0.595896     0.213392     1.011977     1.700000     1.050000
    0.513469     2.310071     0.845290     1.200000     1.050000
   -0.369499     0.121508     1.489233     1.200000     1.050000
    1.478063    -1.000287     1.004227     1.700000     1.050000
    1.940410    -1.028087     1.998902     1.200000     1.050000
    0.895946    -1.915783     0.915280     1.200000     1.050000
    2.423651     1.516981     0.110100     1.700000     1.050000
    2.832064     2.505394    -0.064213     1.200000     1.050000
    3.120849     0.437571    -0.247983     1.700000     1.050000
    4.084599     0.548890    -0.731234     1.200000     1.050000
    2.590846    -0.955955    -0.045910     1.700000     1.050000
    2.243218    -1.349244    -1.009243     1.200000     1.050000
    3.400776    -1.620708     0.261854     1.200000     1.050000
   -1.095787    -1.409118    -0.700116     1.700000     1.050000
   -0.513469    -2.310071    -0.845290     1.200000     1.050000
   -2.423651    -1.516981    -0.110100     1.700000     1.050000
   -2.832064    -2.505394     0.064213     1.200000     1.050000
   -1.478063     1.000287    -1.004227     1.700000     1.050000
   -0.895946     1.915783    -0.915280     1.200000     1.050000
   -1.940410     1.028087    -1.998902     1.200000     1.050000
   -2.590846     0.955955     0.045910     1.700000     1.050000
   -3.400776     1.620708    -0.261855     1.200000     1.050000
   -2.243218     1.349244     1.009243     1.200000     1.050000

--link1--
%chk=meso-TS_19_2_xp-100.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120849  -0.437571   0.247983 
  6  -0.595896  -0.213392  -1.011977 
  1  -4.084599  -0.548890   0.731234 
  1   0.369499  -0.121508  -1.489233 
  6   1.095787   1.409118   0.700116 
  6   0.595896   0.213392   1.011977 
  1   0.513469   2.310071   0.845290 
  1  -0.369499   0.121508   1.489233 
  6   1.478063  -1.000287   1.004227 
  1   1.940410  -1.028087   1.998902 
  1   0.895946  -1.915783   0.915280 
  6   2.423651   1.516981   0.110100 
  1   2.832064   2.505394  -0.064213 
  6   3.120849   0.437571  -0.247983 
  1   4.084599   0.548890  -0.731234 
  6   2.590846  -0.955955  -0.045910 
  1   2.243218  -1.349244  -1.009243 
  1   3.400776  -1.620708   0.261854 
  6  -1.095787  -1.409118  -0.700116 
  1  -0.513469  -2.310071  -0.845290 
  6  -2.423651  -1.516981  -0.110100 
  1  -2.832064  -2.505394   0.064213 
  6  -1.478063   1.000287  -1.004227 
  1  -0.895946   1.915783  -0.915280 
  1  -1.940410   1.028087  -1.998902 
  6  -2.590846   0.955955   0.045910 
  1  -3.400776   1.620708  -0.261855 
  1  -2.243218   1.349244   1.009243 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120849    -0.437571     0.247983     1.700000     1.000000
   -0.595896    -0.213392    -1.011977     1.700000     1.000000
   -4.084599    -0.548890     0.731234     1.200000     1.000000
    0.369499    -0.121508    -1.489233     1.200000     1.000000
    1.095787     1.409118     0.700116     1.700000     1.000000
    0.595896     0.213392     1.011977     1.700000     1.000000
    0.513469     2.310071     0.845290     1.200000     1.000000
   -0.369499     0.121508     1.489233     1.200000     1.000000
    1.478063    -1.000287     1.004227     1.700000     1.000000
    1.940410    -1.028087     1.998902     1.200000     1.000000
    0.895946    -1.915783     0.915280     1.200000     1.000000
    2.423651     1.516981     0.110100     1.700000     1.000000
    2.832064     2.505394    -0.064213     1.200000     1.000000
    3.120849     0.437571    -0.247983     1.700000     1.000000
    4.084599     0.548890    -0.731234     1.200000     1.000000
    2.590846    -0.955955    -0.045910     1.700000     1.000000
    2.243218    -1.349244    -1.009243     1.200000     1.000000
    3.400776    -1.620708     0.261854     1.200000     1.000000
   -1.095787    -1.409118    -0.700116     1.700000     1.000000
   -0.513469    -2.310071    -0.845290     1.200000     1.000000
   -2.423651    -1.516981    -0.110100     1.700000     1.000000
   -2.832064    -2.505394     0.064213     1.200000     1.000000
   -1.478063     1.000287    -1.004227     1.700000     1.000000
   -0.895946     1.915783    -0.915280     1.200000     1.000000
   -1.940410     1.028087    -1.998902     1.200000     1.000000
   -2.590846     0.955955     0.045910     1.700000     1.000000
   -3.400776     1.620708    -0.261855     1.200000     1.000000
   -2.243218     1.349244     1.009243     1.200000     1.000000

--link1--
%chk=meso-TS_19_2_xp-0975.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120849  -0.437571   0.247983 
  6  -0.595896  -0.213392  -1.011977 
  1  -4.084599  -0.548890   0.731234 
  1   0.369499  -0.121508  -1.489233 
  6   1.095787   1.409118   0.700116 
  6   0.595896   0.213392   1.011977 
  1   0.513469   2.310071   0.845290 
  1  -0.369499   0.121508   1.489233 
  6   1.478063  -1.000287   1.004227 
  1   1.940410  -1.028087   1.998902 
  1   0.895946  -1.915783   0.915280 
  6   2.423651   1.516981   0.110100 
  1   2.832064   2.505394  -0.064213 
  6   3.120849   0.437571  -0.247983 
  1   4.084599   0.548890  -0.731234 
  6   2.590846  -0.955955  -0.045910 
  1   2.243218  -1.349244  -1.009243 
  1   3.400776  -1.620708   0.261854 
  6  -1.095787  -1.409118  -0.700116 
  1  -0.513469  -2.310071  -0.845290 
  6  -2.423651  -1.516981  -0.110100 
  1  -2.832064  -2.505394   0.064213 
  6  -1.478063   1.000287  -1.004227 
  1  -0.895946   1.915783  -0.915280 
  1  -1.940410   1.028087  -1.998902 
  6  -2.590846   0.955955   0.045910 
  1  -3.400776   1.620708  -0.261855 
  1  -2.243218   1.349244   1.009243 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120849    -0.437571     0.247983     1.700000     0.975000
   -0.595896    -0.213392    -1.011977     1.700000     0.975000
   -4.084599    -0.548890     0.731234     1.200000     0.975000
    0.369499    -0.121508    -1.489233     1.200000     0.975000
    1.095787     1.409118     0.700116     1.700000     0.975000
    0.595896     0.213392     1.011977     1.700000     0.975000
    0.513469     2.310071     0.845290     1.200000     0.975000
   -0.369499     0.121508     1.489233     1.200000     0.975000
    1.478063    -1.000287     1.004227     1.700000     0.975000
    1.940410    -1.028087     1.998902     1.200000     0.975000
    0.895946    -1.915783     0.915280     1.200000     0.975000
    2.423651     1.516981     0.110100     1.700000     0.975000
    2.832064     2.505394    -0.064213     1.200000     0.975000
    3.120849     0.437571    -0.247983     1.700000     0.975000
    4.084599     0.548890    -0.731234     1.200000     0.975000
    2.590846    -0.955955    -0.045910     1.700000     0.975000
    2.243218    -1.349244    -1.009243     1.200000     0.975000
    3.400776    -1.620708     0.261854     1.200000     0.975000
   -1.095787    -1.409118    -0.700116     1.700000     0.975000
   -0.513469    -2.310071    -0.845290     1.200000     0.975000
   -2.423651    -1.516981    -0.110100     1.700000     0.975000
   -2.832064    -2.505394     0.064213     1.200000     0.975000
   -1.478063     1.000287    -1.004227     1.700000     0.975000
   -0.895946     1.915783    -0.915280     1.200000     0.975000
   -1.940410     1.028087    -1.998902     1.200000     0.975000
   -2.590846     0.955955     0.045910     1.700000     0.975000
   -3.400776     1.620708    -0.261855     1.200000     0.975000
   -2.243218     1.349244     1.009243     1.200000     0.975000

--link1--
%chk=meso-TS_19_2_xp-095.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120849  -0.437571   0.247983 
  6  -0.595896  -0.213392  -1.011977 
  1  -4.084599  -0.548890   0.731234 
  1   0.369499  -0.121508  -1.489233 
  6   1.095787   1.409118   0.700116 
  6   0.595896   0.213392   1.011977 
  1   0.513469   2.310071   0.845290 
  1  -0.369499   0.121508   1.489233 
  6   1.478063  -1.000287   1.004227 
  1   1.940410  -1.028087   1.998902 
  1   0.895946  -1.915783   0.915280 
  6   2.423651   1.516981   0.110100 
  1   2.832064   2.505394  -0.064213 
  6   3.120849   0.437571  -0.247983 
  1   4.084599   0.548890  -0.731234 
  6   2.590846  -0.955955  -0.045910 
  1   2.243218  -1.349244  -1.009243 
  1   3.400776  -1.620708   0.261854 
  6  -1.095787  -1.409118  -0.700116 
  1  -0.513469  -2.310071  -0.845290 
  6  -2.423651  -1.516981  -0.110100 
  1  -2.832064  -2.505394   0.064213 
  6  -1.478063   1.000287  -1.004227 
  1  -0.895946   1.915783  -0.915280 
  1  -1.940410   1.028087  -1.998902 
  6  -2.590846   0.955955   0.045910 
  1  -3.400776   1.620708  -0.261855 
  1  -2.243218   1.349244   1.009243 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120849    -0.437571     0.247983     1.700000     0.950000
   -0.595896    -0.213392    -1.011977     1.700000     0.950000
   -4.084599    -0.548890     0.731234     1.200000     0.950000
    0.369499    -0.121508    -1.489233     1.200000     0.950000
    1.095787     1.409118     0.700116     1.700000     0.950000
    0.595896     0.213392     1.011977     1.700000     0.950000
    0.513469     2.310071     0.845290     1.200000     0.950000
   -0.369499     0.121508     1.489233     1.200000     0.950000
    1.478063    -1.000287     1.004227     1.700000     0.950000
    1.940410    -1.028087     1.998902     1.200000     0.950000
    0.895946    -1.915783     0.915280     1.200000     0.950000
    2.423651     1.516981     0.110100     1.700000     0.950000
    2.832064     2.505394    -0.064213     1.200000     0.950000
    3.120849     0.437571    -0.247983     1.700000     0.950000
    4.084599     0.548890    -0.731234     1.200000     0.950000
    2.590846    -0.955955    -0.045910     1.700000     0.950000
    2.243218    -1.349244    -1.009243     1.200000     0.950000
    3.400776    -1.620708     0.261854     1.200000     0.950000
   -1.095787    -1.409118    -0.700116     1.700000     0.950000
   -0.513469    -2.310071    -0.845290     1.200000     0.950000
   -2.423651    -1.516981    -0.110100     1.700000     0.950000
   -2.832064    -2.505394     0.064213     1.200000     0.950000
   -1.478063     1.000287    -1.004227     1.700000     0.950000
   -0.895946     1.915783    -0.915280     1.200000     0.950000
   -1.940410     1.028087    -1.998902     1.200000     0.950000
   -2.590846     0.955955     0.045910     1.700000     0.950000
   -3.400776     1.620708    -0.261855     1.200000     0.950000
   -2.243218     1.349244     1.009243     1.200000     0.950000

