%chk=meso-TS_1_1_xp-12.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.118978  -0.431679   0.247428 
  6  -0.498946  -0.160979  -0.842137 
  1  -4.084002  -0.548806   0.725908 
  1   0.353974  -0.115713  -1.511513 
  6   1.104740   1.419096   0.680914 
  6   0.498946   0.160979   0.842137 
  1   0.544380   2.322265   0.886395 
  1  -0.353974   0.115713   1.511513 
  6   1.462618  -0.999815   0.987483 
  1   1.898715  -0.937376   1.990350 
  1   0.929624  -1.949747   0.936184 
  6   2.406372   1.522548   0.118104 
  1   2.824451   2.509227  -0.045547 
  6   3.118978   0.431679  -0.247428 
  1   4.084002   0.548806  -0.725908 
  6   2.592276  -0.958579  -0.043389 
  1   2.248579  -1.357021  -1.006512 
  1   3.401237  -1.622533   0.271000 
  6  -1.104740  -1.419096  -0.680914 
  1  -0.544380  -2.322265  -0.886395 
  6  -2.406372  -1.522548  -0.118104 
  1  -2.824451  -2.509227   0.045547 
  6  -1.462618   0.999815  -0.987483 
  1  -0.929624   1.949747  -0.936184 
  1  -1.898715   0.937376  -1.990350 
  6  -2.592276   0.958579   0.043389 
  1  -3.401236   1.622533  -0.271000 
  1  -2.248578   1.357021   1.006512 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.118978    -0.431679     0.247428     1.700000     1.200000
   -0.498946    -0.160979    -0.842137     1.700000     1.200000
   -4.084002    -0.548806     0.725908     1.200000     1.200000
    0.353974    -0.115713    -1.511513     1.200000     1.200000
    1.104740     1.419096     0.680914     1.700000     1.200000
    0.498946     0.160979     0.842137     1.700000     1.200000
    0.544380     2.322265     0.886395     1.200000     1.200000
   -0.353974     0.115713     1.511513     1.200000     1.200000
    1.462618    -0.999815     0.987483     1.700000     1.200000
    1.898715    -0.937376     1.990350     1.200000     1.200000
    0.929624    -1.949747     0.936184     1.200000     1.200000
    2.406372     1.522548     0.118104     1.700000     1.200000
    2.824451     2.509227    -0.045547     1.200000     1.200000
    3.118978     0.431679    -0.247428     1.700000     1.200000
    4.084002     0.548806    -0.725908     1.200000     1.200000
    2.592276    -0.958579    -0.043389     1.700000     1.200000
    2.248579    -1.357021    -1.006512     1.200000     1.200000
    3.401237    -1.622533     0.271000     1.200000     1.200000
   -1.104740    -1.419096    -0.680914     1.700000     1.200000
   -0.544380    -2.322265    -0.886395     1.200000     1.200000
   -2.406372    -1.522548    -0.118104     1.700000     1.200000
   -2.824451    -2.509227     0.045547     1.200000     1.200000
   -1.462618     0.999815    -0.987483     1.700000     1.200000
   -0.929624     1.949747    -0.936184     1.200000     1.200000
   -1.898715     0.937376    -1.990350     1.200000     1.200000
   -2.592276     0.958579     0.043389     1.700000     1.200000
   -3.401236     1.622533    -0.271000     1.200000     1.200000
   -2.248578     1.357021     1.006512     1.200000     1.200000

--link1--
%chk=meso-TS_1_1_xp-115.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.118978  -0.431679   0.247428 
  6  -0.498946  -0.160979  -0.842137 
  1  -4.084002  -0.548806   0.725908 
  1   0.353974  -0.115713  -1.511513 
  6   1.104740   1.419096   0.680914 
  6   0.498946   0.160979   0.842137 
  1   0.544380   2.322265   0.886395 
  1  -0.353974   0.115713   1.511513 
  6   1.462618  -0.999815   0.987483 
  1   1.898715  -0.937376   1.990350 
  1   0.929624  -1.949747   0.936184 
  6   2.406372   1.522548   0.118104 
  1   2.824451   2.509227  -0.045547 
  6   3.118978   0.431679  -0.247428 
  1   4.084002   0.548806  -0.725908 
  6   2.592276  -0.958579  -0.043389 
  1   2.248579  -1.357021  -1.006512 
  1   3.401237  -1.622533   0.271000 
  6  -1.104740  -1.419096  -0.680914 
  1  -0.544380  -2.322265  -0.886395 
  6  -2.406372  -1.522548  -0.118104 
  1  -2.824451  -2.509227   0.045547 
  6  -1.462618   0.999815  -0.987483 
  1  -0.929624   1.949747  -0.936184 
  1  -1.898715   0.937376  -1.990350 
  6  -2.592276   0.958579   0.043389 
  1  -3.401236   1.622533  -0.271000 
  1  -2.248578   1.357021   1.006512 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.118978    -0.431679     0.247428     1.700000     1.150000
   -0.498946    -0.160979    -0.842137     1.700000     1.150000
   -4.084002    -0.548806     0.725908     1.200000     1.150000
    0.353974    -0.115713    -1.511513     1.200000     1.150000
    1.104740     1.419096     0.680914     1.700000     1.150000
    0.498946     0.160979     0.842137     1.700000     1.150000
    0.544380     2.322265     0.886395     1.200000     1.150000
   -0.353974     0.115713     1.511513     1.200000     1.150000
    1.462618    -0.999815     0.987483     1.700000     1.150000
    1.898715    -0.937376     1.990350     1.200000     1.150000
    0.929624    -1.949747     0.936184     1.200000     1.150000
    2.406372     1.522548     0.118104     1.700000     1.150000
    2.824451     2.509227    -0.045547     1.200000     1.150000
    3.118978     0.431679    -0.247428     1.700000     1.150000
    4.084002     0.548806    -0.725908     1.200000     1.150000
    2.592276    -0.958579    -0.043389     1.700000     1.150000
    2.248579    -1.357021    -1.006512     1.200000     1.150000
    3.401237    -1.622533     0.271000     1.200000     1.150000
   -1.104740    -1.419096    -0.680914     1.700000     1.150000
   -0.544380    -2.322265    -0.886395     1.200000     1.150000
   -2.406372    -1.522548    -0.118104     1.700000     1.150000
   -2.824451    -2.509227     0.045547     1.200000     1.150000
   -1.462618     0.999815    -0.987483     1.700000     1.150000
   -0.929624     1.949747    -0.936184     1.200000     1.150000
   -1.898715     0.937376    -1.990350     1.200000     1.150000
   -2.592276     0.958579     0.043389     1.700000     1.150000
   -3.401236     1.622533    -0.271000     1.200000     1.150000
   -2.248578     1.357021     1.006512     1.200000     1.150000

--link1--
%chk=meso-TS_1_1_xp-110.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.118978  -0.431679   0.247428 
  6  -0.498946  -0.160979  -0.842137 
  1  -4.084002  -0.548806   0.725908 
  1   0.353974  -0.115713  -1.511513 
  6   1.104740   1.419096   0.680914 
  6   0.498946   0.160979   0.842137 
  1   0.544380   2.322265   0.886395 
  1  -0.353974   0.115713   1.511513 
  6   1.462618  -0.999815   0.987483 
  1   1.898715  -0.937376   1.990350 
  1   0.929624  -1.949747   0.936184 
  6   2.406372   1.522548   0.118104 
  1   2.824451   2.509227  -0.045547 
  6   3.118978   0.431679  -0.247428 
  1   4.084002   0.548806  -0.725908 
  6   2.592276  -0.958579  -0.043389 
  1   2.248579  -1.357021  -1.006512 
  1   3.401237  -1.622533   0.271000 
  6  -1.104740  -1.419096  -0.680914 
  1  -0.544380  -2.322265  -0.886395 
  6  -2.406372  -1.522548  -0.118104 
  1  -2.824451  -2.509227   0.045547 
  6  -1.462618   0.999815  -0.987483 
  1  -0.929624   1.949747  -0.936184 
  1  -1.898715   0.937376  -1.990350 
  6  -2.592276   0.958579   0.043389 
  1  -3.401236   1.622533  -0.271000 
  1  -2.248578   1.357021   1.006512 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.118978    -0.431679     0.247428     1.700000     1.100000
   -0.498946    -0.160979    -0.842137     1.700000     1.100000
   -4.084002    -0.548806     0.725908     1.200000     1.100000
    0.353974    -0.115713    -1.511513     1.200000     1.100000
    1.104740     1.419096     0.680914     1.700000     1.100000
    0.498946     0.160979     0.842137     1.700000     1.100000
    0.544380     2.322265     0.886395     1.200000     1.100000
   -0.353974     0.115713     1.511513     1.200000     1.100000
    1.462618    -0.999815     0.987483     1.700000     1.100000
    1.898715    -0.937376     1.990350     1.200000     1.100000
    0.929624    -1.949747     0.936184     1.200000     1.100000
    2.406372     1.522548     0.118104     1.700000     1.100000
    2.824451     2.509227    -0.045547     1.200000     1.100000
    3.118978     0.431679    -0.247428     1.700000     1.100000
    4.084002     0.548806    -0.725908     1.200000     1.100000
    2.592276    -0.958579    -0.043389     1.700000     1.100000
    2.248579    -1.357021    -1.006512     1.200000     1.100000
    3.401237    -1.622533     0.271000     1.200000     1.100000
   -1.104740    -1.419096    -0.680914     1.700000     1.100000
   -0.544380    -2.322265    -0.886395     1.200000     1.100000
   -2.406372    -1.522548    -0.118104     1.700000     1.100000
   -2.824451    -2.509227     0.045547     1.200000     1.100000
   -1.462618     0.999815    -0.987483     1.700000     1.100000
   -0.929624     1.949747    -0.936184     1.200000     1.100000
   -1.898715     0.937376    -1.990350     1.200000     1.100000
   -2.592276     0.958579     0.043389     1.700000     1.100000
   -3.401236     1.622533    -0.271000     1.200000     1.100000
   -2.248578     1.357021     1.006512     1.200000     1.100000

--link1--
%chk=meso-TS_1_1_xp-105.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.118978  -0.431679   0.247428 
  6  -0.498946  -0.160979  -0.842137 
  1  -4.084002  -0.548806   0.725908 
  1   0.353974  -0.115713  -1.511513 
  6   1.104740   1.419096   0.680914 
  6   0.498946   0.160979   0.842137 
  1   0.544380   2.322265   0.886395 
  1  -0.353974   0.115713   1.511513 
  6   1.462618  -0.999815   0.987483 
  1   1.898715  -0.937376   1.990350 
  1   0.929624  -1.949747   0.936184 
  6   2.406372   1.522548   0.118104 
  1   2.824451   2.509227  -0.045547 
  6   3.118978   0.431679  -0.247428 
  1   4.084002   0.548806  -0.725908 
  6   2.592276  -0.958579  -0.043389 
  1   2.248579  -1.357021  -1.006512 
  1   3.401237  -1.622533   0.271000 
  6  -1.104740  -1.419096  -0.680914 
  1  -0.544380  -2.322265  -0.886395 
  6  -2.406372  -1.522548  -0.118104 
  1  -2.824451  -2.509227   0.045547 
  6  -1.462618   0.999815  -0.987483 
  1  -0.929624   1.949747  -0.936184 
  1  -1.898715   0.937376  -1.990350 
  6  -2.592276   0.958579   0.043389 
  1  -3.401236   1.622533  -0.271000 
  1  -2.248578   1.357021   1.006512 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.118978    -0.431679     0.247428     1.700000     1.050000
   -0.498946    -0.160979    -0.842137     1.700000     1.050000
   -4.084002    -0.548806     0.725908     1.200000     1.050000
    0.353974    -0.115713    -1.511513     1.200000     1.050000
    1.104740     1.419096     0.680914     1.700000     1.050000
    0.498946     0.160979     0.842137     1.700000     1.050000
    0.544380     2.322265     0.886395     1.200000     1.050000
   -0.353974     0.115713     1.511513     1.200000     1.050000
    1.462618    -0.999815     0.987483     1.700000     1.050000
    1.898715    -0.937376     1.990350     1.200000     1.050000
    0.929624    -1.949747     0.936184     1.200000     1.050000
    2.406372     1.522548     0.118104     1.700000     1.050000
    2.824451     2.509227    -0.045547     1.200000     1.050000
    3.118978     0.431679    -0.247428     1.700000     1.050000
    4.084002     0.548806    -0.725908     1.200000     1.050000
    2.592276    -0.958579    -0.043389     1.700000     1.050000
    2.248579    -1.357021    -1.006512     1.200000     1.050000
    3.401237    -1.622533     0.271000     1.200000     1.050000
   -1.104740    -1.419096    -0.680914     1.700000     1.050000
   -0.544380    -2.322265    -0.886395     1.200000     1.050000
   -2.406372    -1.522548    -0.118104     1.700000     1.050000
   -2.824451    -2.509227     0.045547     1.200000     1.050000
   -1.462618     0.999815    -0.987483     1.700000     1.050000
   -0.929624     1.949747    -0.936184     1.200000     1.050000
   -1.898715     0.937376    -1.990350     1.200000     1.050000
   -2.592276     0.958579     0.043389     1.700000     1.050000
   -3.401236     1.622533    -0.271000     1.200000     1.050000
   -2.248578     1.357021     1.006512     1.200000     1.050000

--link1--
%chk=meso-TS_1_1_xp-100.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.118978  -0.431679   0.247428 
  6  -0.498946  -0.160979  -0.842137 
  1  -4.084002  -0.548806   0.725908 
  1   0.353974  -0.115713  -1.511513 
  6   1.104740   1.419096   0.680914 
  6   0.498946   0.160979   0.842137 
  1   0.544380   2.322265   0.886395 
  1  -0.353974   0.115713   1.511513 
  6   1.462618  -0.999815   0.987483 
  1   1.898715  -0.937376   1.990350 
  1   0.929624  -1.949747   0.936184 
  6   2.406372   1.522548   0.118104 
  1   2.824451   2.509227  -0.045547 
  6   3.118978   0.431679  -0.247428 
  1   4.084002   0.548806  -0.725908 
  6   2.592276  -0.958579  -0.043389 
  1   2.248579  -1.357021  -1.006512 
  1   3.401237  -1.622533   0.271000 
  6  -1.104740  -1.419096  -0.680914 
  1  -0.544380  -2.322265  -0.886395 
  6  -2.406372  -1.522548  -0.118104 
  1  -2.824451  -2.509227   0.045547 
  6  -1.462618   0.999815  -0.987483 
  1  -0.929624   1.949747  -0.936184 
  1  -1.898715   0.937376  -1.990350 
  6  -2.592276   0.958579   0.043389 
  1  -3.401236   1.622533  -0.271000 
  1  -2.248578   1.357021   1.006512 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.118978    -0.431679     0.247428     1.700000     1.000000
   -0.498946    -0.160979    -0.842137     1.700000     1.000000
   -4.084002    -0.548806     0.725908     1.200000     1.000000
    0.353974    -0.115713    -1.511513     1.200000     1.000000
    1.104740     1.419096     0.680914     1.700000     1.000000
    0.498946     0.160979     0.842137     1.700000     1.000000
    0.544380     2.322265     0.886395     1.200000     1.000000
   -0.353974     0.115713     1.511513     1.200000     1.000000
    1.462618    -0.999815     0.987483     1.700000     1.000000
    1.898715    -0.937376     1.990350     1.200000     1.000000
    0.929624    -1.949747     0.936184     1.200000     1.000000
    2.406372     1.522548     0.118104     1.700000     1.000000
    2.824451     2.509227    -0.045547     1.200000     1.000000
    3.118978     0.431679    -0.247428     1.700000     1.000000
    4.084002     0.548806    -0.725908     1.200000     1.000000
    2.592276    -0.958579    -0.043389     1.700000     1.000000
    2.248579    -1.357021    -1.006512     1.200000     1.000000
    3.401237    -1.622533     0.271000     1.200000     1.000000
   -1.104740    -1.419096    -0.680914     1.700000     1.000000
   -0.544380    -2.322265    -0.886395     1.200000     1.000000
   -2.406372    -1.522548    -0.118104     1.700000     1.000000
   -2.824451    -2.509227     0.045547     1.200000     1.000000
   -1.462618     0.999815    -0.987483     1.700000     1.000000
   -0.929624     1.949747    -0.936184     1.200000     1.000000
   -1.898715     0.937376    -1.990350     1.200000     1.000000
   -2.592276     0.958579     0.043389     1.700000     1.000000
   -3.401236     1.622533    -0.271000     1.200000     1.000000
   -2.248578     1.357021     1.006512     1.200000     1.000000

--link1--
%chk=meso-TS_1_1_xp-0975.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.118978  -0.431679   0.247428 
  6  -0.498946  -0.160979  -0.842137 
  1  -4.084002  -0.548806   0.725908 
  1   0.353974  -0.115713  -1.511513 
  6   1.104740   1.419096   0.680914 
  6   0.498946   0.160979   0.842137 
  1   0.544380   2.322265   0.886395 
  1  -0.353974   0.115713   1.511513 
  6   1.462618  -0.999815   0.987483 
  1   1.898715  -0.937376   1.990350 
  1   0.929624  -1.949747   0.936184 
  6   2.406372   1.522548   0.118104 
  1   2.824451   2.509227  -0.045547 
  6   3.118978   0.431679  -0.247428 
  1   4.084002   0.548806  -0.725908 
  6   2.592276  -0.958579  -0.043389 
  1   2.248579  -1.357021  -1.006512 
  1   3.401237  -1.622533   0.271000 
  6  -1.104740  -1.419096  -0.680914 
  1  -0.544380  -2.322265  -0.886395 
  6  -2.406372  -1.522548  -0.118104 
  1  -2.824451  -2.509227   0.045547 
  6  -1.462618   0.999815  -0.987483 
  1  -0.929624   1.949747  -0.936184 
  1  -1.898715   0.937376  -1.990350 
  6  -2.592276   0.958579   0.043389 
  1  -3.401236   1.622533  -0.271000 
  1  -2.248578   1.357021   1.006512 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.118978    -0.431679     0.247428     1.700000     0.975000
   -0.498946    -0.160979    -0.842137     1.700000     0.975000
   -4.084002    -0.548806     0.725908     1.200000     0.975000
    0.353974    -0.115713    -1.511513     1.200000     0.975000
    1.104740     1.419096     0.680914     1.700000     0.975000
    0.498946     0.160979     0.842137     1.700000     0.975000
    0.544380     2.322265     0.886395     1.200000     0.975000
   -0.353974     0.115713     1.511513     1.200000     0.975000
    1.462618    -0.999815     0.987483     1.700000     0.975000
    1.898715    -0.937376     1.990350     1.200000     0.975000
    0.929624    -1.949747     0.936184     1.200000     0.975000
    2.406372     1.522548     0.118104     1.700000     0.975000
    2.824451     2.509227    -0.045547     1.200000     0.975000
    3.118978     0.431679    -0.247428     1.700000     0.975000
    4.084002     0.548806    -0.725908     1.200000     0.975000
    2.592276    -0.958579    -0.043389     1.700000     0.975000
    2.248579    -1.357021    -1.006512     1.200000     0.975000
    3.401237    -1.622533     0.271000     1.200000     0.975000
   -1.104740    -1.419096    -0.680914     1.700000     0.975000
   -0.544380    -2.322265    -0.886395     1.200000     0.975000
   -2.406372    -1.522548    -0.118104     1.700000     0.975000
   -2.824451    -2.509227     0.045547     1.200000     0.975000
   -1.462618     0.999815    -0.987483     1.700000     0.975000
   -0.929624     1.949747    -0.936184     1.200000     0.975000
   -1.898715     0.937376    -1.990350     1.200000     0.975000
   -2.592276     0.958579     0.043389     1.700000     0.975000
   -3.401236     1.622533    -0.271000     1.200000     0.975000
   -2.248578     1.357021     1.006512     1.200000     0.975000

--link1--
%chk=meso-TS_1_1_xp-095.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.118978  -0.431679   0.247428 
  6  -0.498946  -0.160979  -0.842137 
  1  -4.084002  -0.548806   0.725908 
  1   0.353974  -0.115713  -1.511513 
  6   1.104740   1.419096   0.680914 
  6   0.498946   0.160979   0.842137 
  1   0.544380   2.322265   0.886395 
  1  -0.353974   0.115713   1.511513 
  6   1.462618  -0.999815   0.987483 
  1   1.898715  -0.937376   1.990350 
  1   0.929624  -1.949747   0.936184 
  6   2.406372   1.522548   0.118104 
  1   2.824451   2.509227  -0.045547 
  6   3.118978   0.431679  -0.247428 
  1   4.084002   0.548806  -0.725908 
  6   2.592276  -0.958579  -0.043389 
  1   2.248579  -1.357021  -1.006512 
  1   3.401237  -1.622533   0.271000 
  6  -1.104740  -1.419096  -0.680914 
  1  -0.544380  -2.322265  -0.886395 
  6  -2.406372  -1.522548  -0.118104 
  1  -2.824451  -2.509227   0.045547 
  6  -1.462618   0.999815  -0.987483 
  1  -0.929624   1.949747  -0.936184 
  1  -1.898715   0.937376  -1.990350 
  6  -2.592276   0.958579   0.043389 
  1  -3.401236   1.622533  -0.271000 
  1  -2.248578   1.357021   1.006512 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.118978    -0.431679     0.247428     1.700000     0.950000
   -0.498946    -0.160979    -0.842137     1.700000     0.950000
   -4.084002    -0.548806     0.725908     1.200000     0.950000
    0.353974    -0.115713    -1.511513     1.200000     0.950000
    1.104740     1.419096     0.680914     1.700000     0.950000
    0.498946     0.160979     0.842137     1.700000     0.950000
    0.544380     2.322265     0.886395     1.200000     0.950000
   -0.353974     0.115713     1.511513     1.200000     0.950000
    1.462618    -0.999815     0.987483     1.700000     0.950000
    1.898715    -0.937376     1.990350     1.200000     0.950000
    0.929624    -1.949747     0.936184     1.200000     0.950000
    2.406372     1.522548     0.118104     1.700000     0.950000
    2.824451     2.509227    -0.045547     1.200000     0.950000
    3.118978     0.431679    -0.247428     1.700000     0.950000
    4.084002     0.548806    -0.725908     1.200000     0.950000
    2.592276    -0.958579    -0.043389     1.700000     0.950000
    2.248579    -1.357021    -1.006512     1.200000     0.950000
    3.401237    -1.622533     0.271000     1.200000     0.950000
   -1.104740    -1.419096    -0.680914     1.700000     0.950000
   -0.544380    -2.322265    -0.886395     1.200000     0.950000
   -2.406372    -1.522548    -0.118104     1.700000     0.950000
   -2.824451    -2.509227     0.045547     1.200000     0.950000
   -1.462618     0.999815    -0.987483     1.700000     0.950000
   -0.929624     1.949747    -0.936184     1.200000     0.950000
   -1.898715     0.937376    -1.990350     1.200000     0.950000
   -2.592276     0.958579     0.043389     1.700000     0.950000
   -3.401236     1.622533    -0.271000     1.200000     0.950000
   -2.248578     1.357021     1.006512     1.200000     0.950000

