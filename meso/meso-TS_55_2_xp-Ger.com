%chk=meso-TS_55_2_xp-12.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.127160  -0.441387   0.265808 
  6  -0.742781  -0.253682  -1.253790 
  1  -4.066590  -0.531634   0.798854 
  1   0.204971  -0.172502  -1.771631 
  6   1.175378   1.427467   0.803019 
  6   0.742781   0.253682   1.253790 
  1   0.582914   2.324068   0.940997 
  1  -0.204971   0.172502   1.771631 
  6   1.574588  -0.989808   1.127258 
  1   2.116105  -1.122154   2.072691 
  1   0.927535  -1.861790   1.025285 
  6   2.460878   1.530950   0.111768 
  1   2.860246   2.516008  -0.098202 
  6   3.127160   0.441387  -0.265808 
  1   4.066590   0.531634  -0.798854 
  6   2.572730  -0.937949  -0.034385 
  1   2.091291  -1.273998  -0.960989 
  1   3.384763  -1.644565   0.147485 
  6  -1.175378  -1.427467  -0.803019 
  1  -0.582914  -2.324068  -0.940997 
  6  -2.460878  -1.530950  -0.111768 
  1  -2.860246  -2.516008   0.098202 
  6  -1.574588   0.989808  -1.127258 
  1  -0.927535   1.861790  -1.025285 
  1  -2.116105   1.122154  -2.072691 
  6  -2.572730   0.937949   0.034385 
  1  -3.384763   1.644565  -0.147485 
  1  -2.091291   1.273998   0.960989 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.778099999999998

   -3.127160    -0.441387     0.265808     1.700000     1.200000
   -0.742781    -0.253682    -1.253790     1.700000     1.200000
   -4.066590    -0.531634     0.798854     1.200000     1.200000
    0.204971    -0.172502    -1.771631     1.200000     1.200000
    1.175378     1.427467     0.803019     1.700000     1.200000
    0.742781     0.253682     1.253790     1.700000     1.200000
    0.582914     2.324068     0.940997     1.200000     1.200000
   -0.204971     0.172502     1.771631     1.200000     1.200000
    1.574588    -0.989808     1.127258     1.700000     1.200000
    2.116105    -1.122154     2.072691     1.200000     1.200000
    0.927535    -1.861790     1.025285     1.200000     1.200000
    2.460878     1.530950     0.111768     1.700000     1.200000
    2.860246     2.516008    -0.098202     1.200000     1.200000
    3.127160     0.441387    -0.265808     1.700000     1.200000
    4.066590     0.531634    -0.798854     1.200000     1.200000
    2.572730    -0.937949    -0.034385     1.700000     1.200000
    2.091291    -1.273998    -0.960989     1.200000     1.200000
    3.384763    -1.644565     0.147485     1.200000     1.200000
   -1.175378    -1.427467    -0.803019     1.700000     1.200000
   -0.582914    -2.324068    -0.940997     1.200000     1.200000
   -2.460878    -1.530950    -0.111768     1.700000     1.200000
   -2.860246    -2.516008     0.098202     1.200000     1.200000
   -1.574588     0.989808    -1.127258     1.700000     1.200000
   -0.927535     1.861790    -1.025285     1.200000     1.200000
   -2.116105     1.122154    -2.072691     1.200000     1.200000
   -2.572730     0.937949     0.034385     1.700000     1.200000
   -3.384763     1.644565    -0.147485     1.200000     1.200000
   -2.091291     1.273998     0.960989     1.200000     1.200000

--link1--
%chk=meso-TS_55_2_xp-115.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.127160  -0.441387   0.265808 
  6  -0.742781  -0.253682  -1.253790 
  1  -4.066590  -0.531634   0.798854 
  1   0.204971  -0.172502  -1.771631 
  6   1.175378   1.427467   0.803019 
  6   0.742781   0.253682   1.253790 
  1   0.582914   2.324068   0.940997 
  1  -0.204971   0.172502   1.771631 
  6   1.574588  -0.989808   1.127258 
  1   2.116105  -1.122154   2.072691 
  1   0.927535  -1.861790   1.025285 
  6   2.460878   1.530950   0.111768 
  1   2.860246   2.516008  -0.098202 
  6   3.127160   0.441387  -0.265808 
  1   4.066590   0.531634  -0.798854 
  6   2.572730  -0.937949  -0.034385 
  1   2.091291  -1.273998  -0.960989 
  1   3.384763  -1.644565   0.147485 
  6  -1.175378  -1.427467  -0.803019 
  1  -0.582914  -2.324068  -0.940997 
  6  -2.460878  -1.530950  -0.111768 
  1  -2.860246  -2.516008   0.098202 
  6  -1.574588   0.989808  -1.127258 
  1  -0.927535   1.861790  -1.025285 
  1  -2.116105   1.122154  -2.072691 
  6  -2.572730   0.937949   0.034385 
  1  -3.384763   1.644565  -0.147485 
  1  -2.091291   1.273998   0.960989 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.08074995319988 RhoS=0.87957140140931

   -3.127160    -0.441387     0.265808     1.700000     1.150000
   -0.742781    -0.253682    -1.253790     1.700000     1.150000
   -4.066590    -0.531634     0.798854     1.200000     1.150000
    0.204971    -0.172502    -1.771631     1.200000     1.150000
    1.175378     1.427467     0.803019     1.700000     1.150000
    0.742781     0.253682     1.253790     1.700000     1.150000
    0.582914     2.324068     0.940997     1.200000     1.150000
   -0.204971     0.172502     1.771631     1.200000     1.150000
    1.574588    -0.989808     1.127258     1.700000     1.150000
    2.116105    -1.122154     2.072691     1.200000     1.150000
    0.927535    -1.861790     1.025285     1.200000     1.150000
    2.460878     1.530950     0.111768     1.700000     1.150000
    2.860246     2.516008    -0.098202     1.200000     1.150000
    3.127160     0.441387    -0.265808     1.700000     1.150000
    4.066590     0.531634    -0.798854     1.200000     1.150000
    2.572730    -0.937949    -0.034385     1.700000     1.150000
    2.091291    -1.273998    -0.960989     1.200000     1.150000
    3.384763    -1.644565     0.147485     1.200000     1.150000
   -1.175378    -1.427467    -0.803019     1.700000     1.150000
   -0.582914    -2.324068    -0.940997     1.200000     1.150000
   -2.460878    -1.530950    -0.111768     1.700000     1.150000
   -2.860246    -2.516008     0.098202     1.200000     1.150000
   -1.574588     0.989808    -1.127258     1.700000     1.150000
   -0.927535     1.861790    -1.025285     1.200000     1.150000
   -2.116105     1.122154    -2.072691     1.200000     1.150000
   -2.572730     0.937949     0.034385     1.700000     1.150000
   -3.384763     1.644565    -0.147485     1.200000     1.150000
   -2.091291     1.273998     0.960989     1.200000     1.150000

--link1--
%chk=meso-TS_55_2_xp-110.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.127160  -0.441387   0.265808 
  6  -0.742781  -0.253682  -1.253790 
  1  -4.066590  -0.531634   0.798854 
  1   0.204971  -0.172502  -1.771631 
  6   1.175378   1.427467   0.803019 
  6   0.742781   0.253682   1.253790 
  1   0.582914   2.324068   0.940997 
  1  -0.204971   0.172502   1.771631 
  6   1.574588  -0.989808   1.127258 
  1   2.116105  -1.122154   2.072691 
  1   0.927535  -1.861790   1.025285 
  6   2.460878   1.530950   0.111768 
  1   2.860246   2.516008  -0.098202 
  6   3.127160   0.441387  -0.265808 
  1   4.066590   0.531634  -0.798854 
  6   2.572730  -0.937949  -0.034385 
  1   2.091291  -1.273998  -0.960989 
  1   3.384763  -1.644565   0.147485 
  6  -1.175378  -1.427467  -0.803019 
  1  -0.582914  -2.324068  -0.940997 
  6  -2.460878  -1.530950  -0.111768 
  1  -2.860246  -2.516008   0.098202 
  6  -1.574588   0.989808  -1.127258 
  1  -0.927535   1.861790  -1.025285 
  1  -2.116105   1.122154  -2.072691 
  6  -2.572730   0.937949   0.034385 
  1  -3.384763   1.644565  -0.147485 
  1  -2.091291   1.273998   0.960989 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.14797859684703 RhoS=0.992403372535132

   -3.127160    -0.441387     0.265808     1.700000     1.100000
   -0.742781    -0.253682    -1.253790     1.700000     1.100000
   -4.066590    -0.531634     0.798854     1.200000     1.100000
    0.204971    -0.172502    -1.771631     1.200000     1.100000
    1.175378     1.427467     0.803019     1.700000     1.100000
    0.742781     0.253682     1.253790     1.700000     1.100000
    0.582914     2.324068     0.940997     1.200000     1.100000
   -0.204971     0.172502     1.771631     1.200000     1.100000
    1.574588    -0.989808     1.127258     1.700000     1.100000
    2.116105    -1.122154     2.072691     1.200000     1.100000
    0.927535    -1.861790     1.025285     1.200000     1.100000
    2.460878     1.530950     0.111768     1.700000     1.100000
    2.860246     2.516008    -0.098202     1.200000     1.100000
    3.127160     0.441387    -0.265808     1.700000     1.100000
    4.066590     0.531634    -0.798854     1.200000     1.100000
    2.572730    -0.937949    -0.034385     1.700000     1.100000
    2.091291    -1.273998    -0.960989     1.200000     1.100000
    3.384763    -1.644565     0.147485     1.200000     1.100000
   -1.175378    -1.427467    -0.803019     1.700000     1.100000
   -0.582914    -2.324068    -0.940997     1.200000     1.100000
   -2.460878    -1.530950    -0.111768     1.700000     1.100000
   -2.860246    -2.516008     0.098202     1.200000     1.100000
   -1.574588     0.989808    -1.127258     1.700000     1.100000
   -0.927535     1.861790    -1.025285     1.200000     1.100000
   -2.116105     1.122154    -2.072691     1.200000     1.100000
   -2.572730     0.937949     0.034385     1.700000     1.100000
   -3.384763     1.644565    -0.147485     1.200000     1.100000
   -2.091291     1.273998     0.960989     1.200000     1.100000

--link1--
%chk=meso-TS_55_2_xp-105.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.127160  -0.441387   0.265808 
  6  -0.742781  -0.253682  -1.253790 
  1  -4.066590  -0.531634   0.798854 
  1   0.204971  -0.172502  -1.771631 
  6   1.175378   1.427467   0.803019 
  6   0.742781   0.253682   1.253790 
  1   0.582914   2.324068   0.940997 
  1  -0.204971   0.172502   1.771631 
  6   1.574588  -0.989808   1.127258 
  1   2.116105  -1.122154   2.072691 
  1   0.927535  -1.861790   1.025285 
  6   2.460878   1.530950   0.111768 
  1   2.860246   2.516008  -0.098202 
  6   3.127160   0.441387  -0.265808 
  1   4.066590   0.531634  -0.798854 
  6   2.572730  -0.937949  -0.034385 
  1   2.091291  -1.273998  -0.960989 
  1   3.384763  -1.644565   0.147485 
  6  -1.175378  -1.427467  -0.803019 
  1  -0.582914  -2.324068  -0.940997 
  6  -2.460878  -1.530950  -0.111768 
  1  -2.860246  -2.516008   0.098202 
  6  -1.574588   0.989808  -1.127258 
  1  -0.927535   1.861790  -1.025285 
  1  -2.116105   1.122154  -2.072691 
  6  -2.572730   0.937949   0.034385 
  1  -3.384763   1.644565  -0.147485 
  1  -2.091291   1.273998   0.960989 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.22376880117094 RhoS=1.12776705486046

   -3.127160    -0.441387     0.265808     1.700000     1.050000
   -0.742781    -0.253682    -1.253790     1.700000     1.050000
   -4.066590    -0.531634     0.798854     1.200000     1.050000
    0.204971    -0.172502    -1.771631     1.200000     1.050000
    1.175378     1.427467     0.803019     1.700000     1.050000
    0.742781     0.253682     1.253790     1.700000     1.050000
    0.582914     2.324068     0.940997     1.200000     1.050000
   -0.204971     0.172502     1.771631     1.200000     1.050000
    1.574588    -0.989808     1.127258     1.700000     1.050000
    2.116105    -1.122154     2.072691     1.200000     1.050000
    0.927535    -1.861790     1.025285     1.200000     1.050000
    2.460878     1.530950     0.111768     1.700000     1.050000
    2.860246     2.516008    -0.098202     1.200000     1.050000
    3.127160     0.441387    -0.265808     1.700000     1.050000
    4.066590     0.531634    -0.798854     1.200000     1.050000
    2.572730    -0.937949    -0.034385     1.700000     1.050000
    2.091291    -1.273998    -0.960989     1.200000     1.050000
    3.384763    -1.644565     0.147485     1.200000     1.050000
   -1.175378    -1.427467    -0.803019     1.700000     1.050000
   -0.582914    -2.324068    -0.940997     1.200000     1.050000
   -2.460878    -1.530950    -0.111768     1.700000     1.050000
   -2.860246    -2.516008     0.098202     1.200000     1.050000
   -1.574588     0.989808    -1.127258     1.700000     1.050000
   -0.927535     1.861790    -1.025285     1.200000     1.050000
   -2.116105     1.122154    -2.072691     1.200000     1.050000
   -2.572730     0.937949     0.034385     1.700000     1.050000
   -3.384763     1.644565    -0.147485     1.200000     1.050000
   -2.091291     1.273998     0.960989     1.200000     1.050000

--link1--
%chk=meso-TS_55_2_xp-100.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.127160  -0.441387   0.265808 
  6  -0.742781  -0.253682  -1.253790 
  1  -4.066590  -0.531634   0.798854 
  1   0.204971  -0.172502  -1.771631 
  6   1.175378   1.427467   0.803019 
  6   0.742781   0.253682   1.253790 
  1   0.582914   2.324068   0.940997 
  1  -0.204971   0.172502   1.771631 
  6   1.574588  -0.989808   1.127258 
  1   2.116105  -1.122154   2.072691 
  1   0.927535  -1.861790   1.025285 
  6   2.460878   1.530950   0.111768 
  1   2.860246   2.516008  -0.098202 
  6   3.127160   0.441387  -0.265808 
  1   4.066590   0.531634  -0.798854 
  6   2.572730  -0.937949  -0.034385 
  1   2.091291  -1.273998  -0.960989 
  1   3.384763  -1.644565   0.147485 
  6  -1.175378  -1.427467  -0.803019 
  1  -0.582914  -2.324068  -0.940997 
  6  -2.460878  -1.530950  -0.111768 
  1  -2.860246  -2.516008   0.098202 
  6  -1.574588   0.989808  -1.127258 
  1  -0.927535   1.861790  -1.025285 
  1  -2.116105   1.122154  -2.072691 
  6  -2.572730   0.937949   0.034385 
  1  -3.384763   1.644565  -0.147485 
  1  -2.091291   1.273998   0.960989 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.30480546666792 RhoS=1.28207131826017

   -3.127160    -0.441387     0.265808     1.700000     1.000000
   -0.742781    -0.253682    -1.253790     1.700000     1.000000
   -4.066590    -0.531634     0.798854     1.200000     1.000000
    0.204971    -0.172502    -1.771631     1.200000     1.000000
    1.175378     1.427467     0.803019     1.700000     1.000000
    0.742781     0.253682     1.253790     1.700000     1.000000
    0.582914     2.324068     0.940997     1.200000     1.000000
   -0.204971     0.172502     1.771631     1.200000     1.000000
    1.574588    -0.989808     1.127258     1.700000     1.000000
    2.116105    -1.122154     2.072691     1.200000     1.000000
    0.927535    -1.861790     1.025285     1.200000     1.000000
    2.460878     1.530950     0.111768     1.700000     1.000000
    2.860246     2.516008    -0.098202     1.200000     1.000000
    3.127160     0.441387    -0.265808     1.700000     1.000000
    4.066590     0.531634    -0.798854     1.200000     1.000000
    2.572730    -0.937949    -0.034385     1.700000     1.000000
    2.091291    -1.273998    -0.960989     1.200000     1.000000
    3.384763    -1.644565     0.147485     1.200000     1.000000
   -1.175378    -1.427467    -0.803019     1.700000     1.000000
   -0.582914    -2.324068    -0.940997     1.200000     1.000000
   -2.460878    -1.530950    -0.111768     1.700000     1.000000
   -2.860246    -2.516008     0.098202     1.200000     1.000000
   -1.574588     0.989808    -1.127258     1.700000     1.000000
   -0.927535     1.861790    -1.025285     1.200000     1.000000
   -2.116105     1.122154    -2.072691     1.200000     1.000000
   -2.572730     0.937949     0.034385     1.700000     1.000000
   -3.384763     1.644565    -0.147485     1.200000     1.000000
   -2.091291     1.273998     0.960989     1.200000     1.000000

--link1--
%chk=meso-TS_55_2_xp-0975.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.127160  -0.441387   0.265808 
  6  -0.742781  -0.253682  -1.253790 
  1  -4.066590  -0.531634   0.798854 
  1   0.204971  -0.172502  -1.771631 
  6   1.175378   1.427467   0.803019 
  6   0.742781   0.253682   1.253790 
  1   0.582914   2.324068   0.940997 
  1  -0.204971   0.172502   1.771631 
  6   1.574588  -0.989808   1.127258 
  1   2.116105  -1.122154   2.072691 
  1   0.927535  -1.861790   1.025285 
  6   2.460878   1.530950   0.111768 
  1   2.860246   2.516008  -0.098202 
  6   3.127160   0.441387  -0.265808 
  1   4.066590   0.531634  -0.798854 
  6   2.572730  -0.937949  -0.034385 
  1   2.091291  -1.273998  -0.960989 
  1   3.384763  -1.644565   0.147485 
  6  -1.175378  -1.427467  -0.803019 
  1  -0.582914  -2.324068  -0.940997 
  6  -2.460878  -1.530950  -0.111768 
  1  -2.860246  -2.516008   0.098202 
  6  -1.574588   0.989808  -1.127258 
  1  -0.927535   1.861790  -1.025285 
  1  -2.116105   1.122154  -2.072691 
  6  -2.572730   0.937949   0.034385 
  1  -3.384763   1.644565  -0.147485 
  1  -2.091291   1.273998   0.960989 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.34745273026226 RhoS=1.36724935387803

   -3.127160    -0.441387     0.265808     1.700000     0.975000
   -0.742781    -0.253682    -1.253790     1.700000     0.975000
   -4.066590    -0.531634     0.798854     1.200000     0.975000
    0.204971    -0.172502    -1.771631     1.200000     0.975000
    1.175378     1.427467     0.803019     1.700000     0.975000
    0.742781     0.253682     1.253790     1.700000     0.975000
    0.582914     2.324068     0.940997     1.200000     0.975000
   -0.204971     0.172502     1.771631     1.200000     0.975000
    1.574588    -0.989808     1.127258     1.700000     0.975000
    2.116105    -1.122154     2.072691     1.200000     0.975000
    0.927535    -1.861790     1.025285     1.200000     0.975000
    2.460878     1.530950     0.111768     1.700000     0.975000
    2.860246     2.516008    -0.098202     1.200000     0.975000
    3.127160     0.441387    -0.265808     1.700000     0.975000
    4.066590     0.531634    -0.798854     1.200000     0.975000
    2.572730    -0.937949    -0.034385     1.700000     0.975000
    2.091291    -1.273998    -0.960989     1.200000     0.975000
    3.384763    -1.644565     0.147485     1.200000     0.975000
   -1.175378    -1.427467    -0.803019     1.700000     0.975000
   -0.582914    -2.324068    -0.940997     1.200000     0.975000
   -2.460878    -1.530950    -0.111768     1.700000     0.975000
   -2.860246    -2.516008     0.098202     1.200000     0.975000
   -1.574588     0.989808    -1.127258     1.700000     0.975000
   -0.927535     1.861790    -1.025285     1.200000     0.975000
   -2.116105     1.122154    -2.072691     1.200000     0.975000
   -2.572730     0.937949     0.034385     1.700000     0.975000
   -3.384763     1.644565    -0.147485     1.200000     0.975000
   -2.091291     1.273998     0.960989     1.200000     0.975000

--link1--
%chk=meso-TS_55_2_xp-095.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym 6d 10f 

endo-2
xp-pcm for eta=3

0 1
  6  -3.127160  -0.441387   0.265808 
  6  -0.742781  -0.253682  -1.253790 
  1  -4.066590  -0.531634   0.798854 
  1   0.204971  -0.172502  -1.771631 
  6   1.175378   1.427467   0.803019 
  6   0.742781   0.253682   1.253790 
  1   0.582914   2.324068   0.940997 
  1  -0.204971   0.172502   1.771631 
  6   1.574588  -0.989808   1.127258 
  1   2.116105  -1.122154   2.072691 
  1   0.927535  -1.861790   1.025285 
  6   2.460878   1.530950   0.111768 
  1   2.860246   2.516008  -0.098202 
  6   3.127160   0.441387  -0.265808 
  1   4.066590   0.531634  -0.798854 
  6   2.572730  -0.937949  -0.034385 
  1   2.091291  -1.273998  -0.960989 
  1   3.384763  -1.644565   0.147485 
  6  -1.175378  -1.427467  -0.803019 
  1  -0.582914  -2.324068  -0.940997 
  6  -2.460878  -1.530950  -0.111768 
  1  -2.860246  -2.516008   0.098202 
  6  -1.574588   0.989808  -1.127258 
  1  -0.927535   1.861790  -1.025285 
  1  -2.116105   1.122154  -2.072691 
  6  -2.572730   0.937949   0.034385 
  1  -3.384763   1.644565  -0.147485 
  1  -2.091291   1.273998   0.960989 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.39339884323873 RhoS=1.4620813611373

   -3.127160    -0.441387     0.265808     1.700000     0.950000
   -0.742781    -0.253682    -1.253790     1.700000     0.950000
   -4.066590    -0.531634     0.798854     1.200000     0.950000
    0.204971    -0.172502    -1.771631     1.200000     0.950000
    1.175378     1.427467     0.803019     1.700000     0.950000
    0.742781     0.253682     1.253790     1.700000     0.950000
    0.582914     2.324068     0.940997     1.200000     0.950000
   -0.204971     0.172502     1.771631     1.200000     0.950000
    1.574588    -0.989808     1.127258     1.700000     0.950000
    2.116105    -1.122154     2.072691     1.200000     0.950000
    0.927535    -1.861790     1.025285     1.200000     0.950000
    2.460878     1.530950     0.111768     1.700000     0.950000
    2.860246     2.516008    -0.098202     1.200000     0.950000
    3.127160     0.441387    -0.265808     1.700000     0.950000
    4.066590     0.531634    -0.798854     1.200000     0.950000
    2.572730    -0.937949    -0.034385     1.700000     0.950000
    2.091291    -1.273998    -0.960989     1.200000     0.950000
    3.384763    -1.644565     0.147485     1.200000     0.950000
   -1.175378    -1.427467    -0.803019     1.700000     0.950000
   -0.582914    -2.324068    -0.940997     1.200000     0.950000
   -2.460878    -1.530950    -0.111768     1.700000     0.950000
   -2.860246    -2.516008     0.098202     1.200000     0.950000
   -1.574588     0.989808    -1.127258     1.700000     0.950000
   -0.927535     1.861790    -1.025285     1.200000     0.950000
   -2.116105     1.122154    -2.072691     1.200000     0.950000
   -2.572730     0.937949     0.034385     1.700000     0.950000
   -3.384763     1.644565    -0.147485     1.200000     0.950000
   -2.091291     1.273998     0.960989     1.200000     0.950000

