%chk=meso-TS_16_2_xp-12.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120053  -0.437296   0.247601 
  6  -0.582765  -0.208122  -0.986031 
  1  -4.084035  -0.549421   0.730161 
  1   0.376360  -0.119236  -1.476138 
  6   1.092928   1.407235   0.695695 
  6   0.582765   0.208122   0.986031 
  1   0.514464   2.309116   0.849982 
  1  -0.376360   0.119237   1.476138 
  6   1.473441  -0.999882   0.998114 
  1   1.930236  -1.007810   1.995578 
  1   0.900886  -1.922279   0.917856 
  6   2.421602   1.516304   0.111065 
  1   2.830156   2.505119  -0.060741 
  6   3.120053   0.437296  -0.247601 
  1   4.084035   0.549421  -0.730161 
  6   2.591259  -0.956820  -0.045733 
  1   2.248403  -1.352514  -1.009798 
  1   3.400569  -1.620098   0.266906 
  6  -1.092928  -1.407235  -0.695695 
  1  -0.514464  -2.309116  -0.849982 
  6  -2.421602  -1.516304  -0.111065 
  1  -2.830156  -2.505120   0.060741 
  6  -1.473441   0.999882  -0.998114 
  1  -0.900886   1.922279  -0.917856 
  1  -1.930236   1.007810  -1.995578 
  6  -2.591259   0.956820   0.045733 
  1  -3.400569   1.620098  -0.266906 
  1  -2.248403   1.352515   1.009798 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120053    -0.437296     0.247601     1.700000     1.200000
   -0.582765    -0.208122    -0.986031     1.700000     1.200000
   -4.084035    -0.549421     0.730161     1.200000     1.200000
    0.376360    -0.119236    -1.476138     1.200000     1.200000
    1.092928     1.407235     0.695695     1.700000     1.200000
    0.582765     0.208122     0.986031     1.700000     1.200000
    0.514464     2.309116     0.849982     1.200000     1.200000
   -0.376360     0.119237     1.476138     1.200000     1.200000
    1.473441    -0.999882     0.998114     1.700000     1.200000
    1.930236    -1.007810     1.995578     1.200000     1.200000
    0.900886    -1.922279     0.917856     1.200000     1.200000
    2.421602     1.516304     0.111065     1.700000     1.200000
    2.830156     2.505119    -0.060741     1.200000     1.200000
    3.120053     0.437296    -0.247601     1.700000     1.200000
    4.084035     0.549421    -0.730161     1.200000     1.200000
    2.591259    -0.956820    -0.045733     1.700000     1.200000
    2.248403    -1.352514    -1.009798     1.200000     1.200000
    3.400569    -1.620098     0.266906     1.200000     1.200000
   -1.092928    -1.407235    -0.695695     1.700000     1.200000
   -0.514464    -2.309116    -0.849982     1.200000     1.200000
   -2.421602    -1.516304    -0.111065     1.700000     1.200000
   -2.830156    -2.505120     0.060741     1.200000     1.200000
   -1.473441     0.999882    -0.998114     1.700000     1.200000
   -0.900886     1.922279    -0.917856     1.200000     1.200000
   -1.930236     1.007810    -1.995578     1.200000     1.200000
   -2.591259     0.956820     0.045733     1.700000     1.200000
   -3.400569     1.620098    -0.266906     1.200000     1.200000
   -2.248403     1.352515     1.009798     1.200000     1.200000

--link1--
%chk=meso-TS_16_2_xp-115.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120053  -0.437296   0.247601 
  6  -0.582765  -0.208122  -0.986031 
  1  -4.084035  -0.549421   0.730161 
  1   0.376360  -0.119236  -1.476138 
  6   1.092928   1.407235   0.695695 
  6   0.582765   0.208122   0.986031 
  1   0.514464   2.309116   0.849982 
  1  -0.376360   0.119237   1.476138 
  6   1.473441  -0.999882   0.998114 
  1   1.930236  -1.007810   1.995578 
  1   0.900886  -1.922279   0.917856 
  6   2.421602   1.516304   0.111065 
  1   2.830156   2.505119  -0.060741 
  6   3.120053   0.437296  -0.247601 
  1   4.084035   0.549421  -0.730161 
  6   2.591259  -0.956820  -0.045733 
  1   2.248403  -1.352514  -1.009798 
  1   3.400569  -1.620098   0.266906 
  6  -1.092928  -1.407235  -0.695695 
  1  -0.514464  -2.309116  -0.849982 
  6  -2.421602  -1.516304  -0.111065 
  1  -2.830156  -2.505120   0.060741 
  6  -1.473441   0.999882  -0.998114 
  1  -0.900886   1.922279  -0.917856 
  1  -1.930236   1.007810  -1.995578 
  6  -2.591259   0.956820   0.045733 
  1  -3.400569   1.620098  -0.266906 
  1  -2.248403   1.352515   1.009798 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120053    -0.437296     0.247601     1.700000     1.150000
   -0.582765    -0.208122    -0.986031     1.700000     1.150000
   -4.084035    -0.549421     0.730161     1.200000     1.150000
    0.376360    -0.119236    -1.476138     1.200000     1.150000
    1.092928     1.407235     0.695695     1.700000     1.150000
    0.582765     0.208122     0.986031     1.700000     1.150000
    0.514464     2.309116     0.849982     1.200000     1.150000
   -0.376360     0.119237     1.476138     1.200000     1.150000
    1.473441    -0.999882     0.998114     1.700000     1.150000
    1.930236    -1.007810     1.995578     1.200000     1.150000
    0.900886    -1.922279     0.917856     1.200000     1.150000
    2.421602     1.516304     0.111065     1.700000     1.150000
    2.830156     2.505119    -0.060741     1.200000     1.150000
    3.120053     0.437296    -0.247601     1.700000     1.150000
    4.084035     0.549421    -0.730161     1.200000     1.150000
    2.591259    -0.956820    -0.045733     1.700000     1.150000
    2.248403    -1.352514    -1.009798     1.200000     1.150000
    3.400569    -1.620098     0.266906     1.200000     1.150000
   -1.092928    -1.407235    -0.695695     1.700000     1.150000
   -0.514464    -2.309116    -0.849982     1.200000     1.150000
   -2.421602    -1.516304    -0.111065     1.700000     1.150000
   -2.830156    -2.505120     0.060741     1.200000     1.150000
   -1.473441     0.999882    -0.998114     1.700000     1.150000
   -0.900886     1.922279    -0.917856     1.200000     1.150000
   -1.930236     1.007810    -1.995578     1.200000     1.150000
   -2.591259     0.956820     0.045733     1.700000     1.150000
   -3.400569     1.620098    -0.266906     1.200000     1.150000
   -2.248403     1.352515     1.009798     1.200000     1.150000

--link1--
%chk=meso-TS_16_2_xp-110.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120053  -0.437296   0.247601 
  6  -0.582765  -0.208122  -0.986031 
  1  -4.084035  -0.549421   0.730161 
  1   0.376360  -0.119236  -1.476138 
  6   1.092928   1.407235   0.695695 
  6   0.582765   0.208122   0.986031 
  1   0.514464   2.309116   0.849982 
  1  -0.376360   0.119237   1.476138 
  6   1.473441  -0.999882   0.998114 
  1   1.930236  -1.007810   1.995578 
  1   0.900886  -1.922279   0.917856 
  6   2.421602   1.516304   0.111065 
  1   2.830156   2.505119  -0.060741 
  6   3.120053   0.437296  -0.247601 
  1   4.084035   0.549421  -0.730161 
  6   2.591259  -0.956820  -0.045733 
  1   2.248403  -1.352514  -1.009798 
  1   3.400569  -1.620098   0.266906 
  6  -1.092928  -1.407235  -0.695695 
  1  -0.514464  -2.309116  -0.849982 
  6  -2.421602  -1.516304  -0.111065 
  1  -2.830156  -2.505120   0.060741 
  6  -1.473441   0.999882  -0.998114 
  1  -0.900886   1.922279  -0.917856 
  1  -1.930236   1.007810  -1.995578 
  6  -2.591259   0.956820   0.045733 
  1  -3.400569   1.620098  -0.266906 
  1  -2.248403   1.352515   1.009798 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120053    -0.437296     0.247601     1.700000     1.100000
   -0.582765    -0.208122    -0.986031     1.700000     1.100000
   -4.084035    -0.549421     0.730161     1.200000     1.100000
    0.376360    -0.119236    -1.476138     1.200000     1.100000
    1.092928     1.407235     0.695695     1.700000     1.100000
    0.582765     0.208122     0.986031     1.700000     1.100000
    0.514464     2.309116     0.849982     1.200000     1.100000
   -0.376360     0.119237     1.476138     1.200000     1.100000
    1.473441    -0.999882     0.998114     1.700000     1.100000
    1.930236    -1.007810     1.995578     1.200000     1.100000
    0.900886    -1.922279     0.917856     1.200000     1.100000
    2.421602     1.516304     0.111065     1.700000     1.100000
    2.830156     2.505119    -0.060741     1.200000     1.100000
    3.120053     0.437296    -0.247601     1.700000     1.100000
    4.084035     0.549421    -0.730161     1.200000     1.100000
    2.591259    -0.956820    -0.045733     1.700000     1.100000
    2.248403    -1.352514    -1.009798     1.200000     1.100000
    3.400569    -1.620098     0.266906     1.200000     1.100000
   -1.092928    -1.407235    -0.695695     1.700000     1.100000
   -0.514464    -2.309116    -0.849982     1.200000     1.100000
   -2.421602    -1.516304    -0.111065     1.700000     1.100000
   -2.830156    -2.505120     0.060741     1.200000     1.100000
   -1.473441     0.999882    -0.998114     1.700000     1.100000
   -0.900886     1.922279    -0.917856     1.200000     1.100000
   -1.930236     1.007810    -1.995578     1.200000     1.100000
   -2.591259     0.956820     0.045733     1.700000     1.100000
   -3.400569     1.620098    -0.266906     1.200000     1.100000
   -2.248403     1.352515     1.009798     1.200000     1.100000

--link1--
%chk=meso-TS_16_2_xp-105.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120053  -0.437296   0.247601 
  6  -0.582765  -0.208122  -0.986031 
  1  -4.084035  -0.549421   0.730161 
  1   0.376360  -0.119236  -1.476138 
  6   1.092928   1.407235   0.695695 
  6   0.582765   0.208122   0.986031 
  1   0.514464   2.309116   0.849982 
  1  -0.376360   0.119237   1.476138 
  6   1.473441  -0.999882   0.998114 
  1   1.930236  -1.007810   1.995578 
  1   0.900886  -1.922279   0.917856 
  6   2.421602   1.516304   0.111065 
  1   2.830156   2.505119  -0.060741 
  6   3.120053   0.437296  -0.247601 
  1   4.084035   0.549421  -0.730161 
  6   2.591259  -0.956820  -0.045733 
  1   2.248403  -1.352514  -1.009798 
  1   3.400569  -1.620098   0.266906 
  6  -1.092928  -1.407235  -0.695695 
  1  -0.514464  -2.309116  -0.849982 
  6  -2.421602  -1.516304  -0.111065 
  1  -2.830156  -2.505120   0.060741 
  6  -1.473441   0.999882  -0.998114 
  1  -0.900886   1.922279  -0.917856 
  1  -1.930236   1.007810  -1.995578 
  6  -2.591259   0.956820   0.045733 
  1  -3.400569   1.620098  -0.266906 
  1  -2.248403   1.352515   1.009798 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120053    -0.437296     0.247601     1.700000     1.050000
   -0.582765    -0.208122    -0.986031     1.700000     1.050000
   -4.084035    -0.549421     0.730161     1.200000     1.050000
    0.376360    -0.119236    -1.476138     1.200000     1.050000
    1.092928     1.407235     0.695695     1.700000     1.050000
    0.582765     0.208122     0.986031     1.700000     1.050000
    0.514464     2.309116     0.849982     1.200000     1.050000
   -0.376360     0.119237     1.476138     1.200000     1.050000
    1.473441    -0.999882     0.998114     1.700000     1.050000
    1.930236    -1.007810     1.995578     1.200000     1.050000
    0.900886    -1.922279     0.917856     1.200000     1.050000
    2.421602     1.516304     0.111065     1.700000     1.050000
    2.830156     2.505119    -0.060741     1.200000     1.050000
    3.120053     0.437296    -0.247601     1.700000     1.050000
    4.084035     0.549421    -0.730161     1.200000     1.050000
    2.591259    -0.956820    -0.045733     1.700000     1.050000
    2.248403    -1.352514    -1.009798     1.200000     1.050000
    3.400569    -1.620098     0.266906     1.200000     1.050000
   -1.092928    -1.407235    -0.695695     1.700000     1.050000
   -0.514464    -2.309116    -0.849982     1.200000     1.050000
   -2.421602    -1.516304    -0.111065     1.700000     1.050000
   -2.830156    -2.505120     0.060741     1.200000     1.050000
   -1.473441     0.999882    -0.998114     1.700000     1.050000
   -0.900886     1.922279    -0.917856     1.200000     1.050000
   -1.930236     1.007810    -1.995578     1.200000     1.050000
   -2.591259     0.956820     0.045733     1.700000     1.050000
   -3.400569     1.620098    -0.266906     1.200000     1.050000
   -2.248403     1.352515     1.009798     1.200000     1.050000

--link1--
%chk=meso-TS_16_2_xp-100.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120053  -0.437296   0.247601 
  6  -0.582765  -0.208122  -0.986031 
  1  -4.084035  -0.549421   0.730161 
  1   0.376360  -0.119236  -1.476138 
  6   1.092928   1.407235   0.695695 
  6   0.582765   0.208122   0.986031 
  1   0.514464   2.309116   0.849982 
  1  -0.376360   0.119237   1.476138 
  6   1.473441  -0.999882   0.998114 
  1   1.930236  -1.007810   1.995578 
  1   0.900886  -1.922279   0.917856 
  6   2.421602   1.516304   0.111065 
  1   2.830156   2.505119  -0.060741 
  6   3.120053   0.437296  -0.247601 
  1   4.084035   0.549421  -0.730161 
  6   2.591259  -0.956820  -0.045733 
  1   2.248403  -1.352514  -1.009798 
  1   3.400569  -1.620098   0.266906 
  6  -1.092928  -1.407235  -0.695695 
  1  -0.514464  -2.309116  -0.849982 
  6  -2.421602  -1.516304  -0.111065 
  1  -2.830156  -2.505120   0.060741 
  6  -1.473441   0.999882  -0.998114 
  1  -0.900886   1.922279  -0.917856 
  1  -1.930236   1.007810  -1.995578 
  6  -2.591259   0.956820   0.045733 
  1  -3.400569   1.620098  -0.266906 
  1  -2.248403   1.352515   1.009798 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120053    -0.437296     0.247601     1.700000     1.000000
   -0.582765    -0.208122    -0.986031     1.700000     1.000000
   -4.084035    -0.549421     0.730161     1.200000     1.000000
    0.376360    -0.119236    -1.476138     1.200000     1.000000
    1.092928     1.407235     0.695695     1.700000     1.000000
    0.582765     0.208122     0.986031     1.700000     1.000000
    0.514464     2.309116     0.849982     1.200000     1.000000
   -0.376360     0.119237     1.476138     1.200000     1.000000
    1.473441    -0.999882     0.998114     1.700000     1.000000
    1.930236    -1.007810     1.995578     1.200000     1.000000
    0.900886    -1.922279     0.917856     1.200000     1.000000
    2.421602     1.516304     0.111065     1.700000     1.000000
    2.830156     2.505119    -0.060741     1.200000     1.000000
    3.120053     0.437296    -0.247601     1.700000     1.000000
    4.084035     0.549421    -0.730161     1.200000     1.000000
    2.591259    -0.956820    -0.045733     1.700000     1.000000
    2.248403    -1.352514    -1.009798     1.200000     1.000000
    3.400569    -1.620098     0.266906     1.200000     1.000000
   -1.092928    -1.407235    -0.695695     1.700000     1.000000
   -0.514464    -2.309116    -0.849982     1.200000     1.000000
   -2.421602    -1.516304    -0.111065     1.700000     1.000000
   -2.830156    -2.505120     0.060741     1.200000     1.000000
   -1.473441     0.999882    -0.998114     1.700000     1.000000
   -0.900886     1.922279    -0.917856     1.200000     1.000000
   -1.930236     1.007810    -1.995578     1.200000     1.000000
   -2.591259     0.956820     0.045733     1.700000     1.000000
   -3.400569     1.620098    -0.266906     1.200000     1.000000
   -2.248403     1.352515     1.009798     1.200000     1.000000

--link1--
%chk=meso-TS_16_2_xp-0975.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120053  -0.437296   0.247601 
  6  -0.582765  -0.208122  -0.986031 
  1  -4.084035  -0.549421   0.730161 
  1   0.376360  -0.119236  -1.476138 
  6   1.092928   1.407235   0.695695 
  6   0.582765   0.208122   0.986031 
  1   0.514464   2.309116   0.849982 
  1  -0.376360   0.119237   1.476138 
  6   1.473441  -0.999882   0.998114 
  1   1.930236  -1.007810   1.995578 
  1   0.900886  -1.922279   0.917856 
  6   2.421602   1.516304   0.111065 
  1   2.830156   2.505119  -0.060741 
  6   3.120053   0.437296  -0.247601 
  1   4.084035   0.549421  -0.730161 
  6   2.591259  -0.956820  -0.045733 
  1   2.248403  -1.352514  -1.009798 
  1   3.400569  -1.620098   0.266906 
  6  -1.092928  -1.407235  -0.695695 
  1  -0.514464  -2.309116  -0.849982 
  6  -2.421602  -1.516304  -0.111065 
  1  -2.830156  -2.505120   0.060741 
  6  -1.473441   0.999882  -0.998114 
  1  -0.900886   1.922279  -0.917856 
  1  -1.930236   1.007810  -1.995578 
  6  -2.591259   0.956820   0.045733 
  1  -3.400569   1.620098  -0.266906 
  1  -2.248403   1.352515   1.009798 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120053    -0.437296     0.247601     1.700000     0.975000
   -0.582765    -0.208122    -0.986031     1.700000     0.975000
   -4.084035    -0.549421     0.730161     1.200000     0.975000
    0.376360    -0.119236    -1.476138     1.200000     0.975000
    1.092928     1.407235     0.695695     1.700000     0.975000
    0.582765     0.208122     0.986031     1.700000     0.975000
    0.514464     2.309116     0.849982     1.200000     0.975000
   -0.376360     0.119237     1.476138     1.200000     0.975000
    1.473441    -0.999882     0.998114     1.700000     0.975000
    1.930236    -1.007810     1.995578     1.200000     0.975000
    0.900886    -1.922279     0.917856     1.200000     0.975000
    2.421602     1.516304     0.111065     1.700000     0.975000
    2.830156     2.505119    -0.060741     1.200000     0.975000
    3.120053     0.437296    -0.247601     1.700000     0.975000
    4.084035     0.549421    -0.730161     1.200000     0.975000
    2.591259    -0.956820    -0.045733     1.700000     0.975000
    2.248403    -1.352514    -1.009798     1.200000     0.975000
    3.400569    -1.620098     0.266906     1.200000     0.975000
   -1.092928    -1.407235    -0.695695     1.700000     0.975000
   -0.514464    -2.309116    -0.849982     1.200000     0.975000
   -2.421602    -1.516304    -0.111065     1.700000     0.975000
   -2.830156    -2.505120     0.060741     1.200000     0.975000
   -1.473441     0.999882    -0.998114     1.700000     0.975000
   -0.900886     1.922279    -0.917856     1.200000     0.975000
   -1.930236     1.007810    -1.995578     1.200000     0.975000
   -2.591259     0.956820     0.045733     1.700000     0.975000
   -3.400569     1.620098    -0.266906     1.200000     0.975000
   -2.248403     1.352515     1.009798     1.200000     0.975000

--link1--
%chk=meso-TS_16_2_xp-095.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.120053  -0.437296   0.247601 
  6  -0.582765  -0.208122  -0.986031 
  1  -4.084035  -0.549421   0.730161 
  1   0.376360  -0.119236  -1.476138 
  6   1.092928   1.407235   0.695695 
  6   0.582765   0.208122   0.986031 
  1   0.514464   2.309116   0.849982 
  1  -0.376360   0.119237   1.476138 
  6   1.473441  -0.999882   0.998114 
  1   1.930236  -1.007810   1.995578 
  1   0.900886  -1.922279   0.917856 
  6   2.421602   1.516304   0.111065 
  1   2.830156   2.505119  -0.060741 
  6   3.120053   0.437296  -0.247601 
  1   4.084035   0.549421  -0.730161 
  6   2.591259  -0.956820  -0.045733 
  1   2.248403  -1.352514  -1.009798 
  1   3.400569  -1.620098   0.266906 
  6  -1.092928  -1.407235  -0.695695 
  1  -0.514464  -2.309116  -0.849982 
  6  -2.421602  -1.516304  -0.111065 
  1  -2.830156  -2.505120   0.060741 
  6  -1.473441   0.999882  -0.998114 
  1  -0.900886   1.922279  -0.917856 
  1  -1.930236   1.007810  -1.995578 
  6  -2.591259   0.956820   0.045733 
  1  -3.400569   1.620098  -0.266906 
  1  -2.248403   1.352515   1.009798 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.120053    -0.437296     0.247601     1.700000     0.950000
   -0.582765    -0.208122    -0.986031     1.700000     0.950000
   -4.084035    -0.549421     0.730161     1.200000     0.950000
    0.376360    -0.119236    -1.476138     1.200000     0.950000
    1.092928     1.407235     0.695695     1.700000     0.950000
    0.582765     0.208122     0.986031     1.700000     0.950000
    0.514464     2.309116     0.849982     1.200000     0.950000
   -0.376360     0.119237     1.476138     1.200000     0.950000
    1.473441    -0.999882     0.998114     1.700000     0.950000
    1.930236    -1.007810     1.995578     1.200000     0.950000
    0.900886    -1.922279     0.917856     1.200000     0.950000
    2.421602     1.516304     0.111065     1.700000     0.950000
    2.830156     2.505119    -0.060741     1.200000     0.950000
    3.120053     0.437296    -0.247601     1.700000     0.950000
    4.084035     0.549421    -0.730161     1.200000     0.950000
    2.591259    -0.956820    -0.045733     1.700000     0.950000
    2.248403    -1.352514    -1.009798     1.200000     0.950000
    3.400569    -1.620098     0.266906     1.200000     0.950000
   -1.092928    -1.407235    -0.695695     1.700000     0.950000
   -0.514464    -2.309116    -0.849982     1.200000     0.950000
   -2.421602    -1.516304    -0.111065     1.700000     0.950000
   -2.830156    -2.505120     0.060741     1.200000     0.950000
   -1.473441     0.999882    -0.998114     1.700000     0.950000
   -0.900886     1.922279    -0.917856     1.200000     0.950000
   -1.930236     1.007810    -1.995578     1.200000     0.950000
   -2.591259     0.956820     0.045733     1.700000     0.950000
   -3.400569     1.620098    -0.266906     1.200000     0.950000
   -2.248403     1.352515     1.009798     1.200000     0.950000

