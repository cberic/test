%chk=meso-TS_19_1_xp-12.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.119127  -0.424214   0.248750 
  6  -0.411716  -0.117621  -0.690842 
  1  -4.089584  -0.546101   0.713998 
  1   0.331104  -0.106248  -1.492216 
  6   1.113220   1.433745   0.656146 
  6   0.411716   0.117621   0.690842 
  1   0.575541   2.332289   0.931979 
  1  -0.331104   0.106248   1.492216 
  6   1.443013  -0.993893   0.965601 
  1   1.846597  -0.847218   1.971019 
  1   0.957769  -1.971415   0.956603 
  6   2.393256   1.531202   0.129580 
  1   2.834891   2.513727   0.002918 
  6   3.119127   0.424214  -0.248750 
  1   4.089584   0.546101  -0.713998 
  6   2.592435  -0.960873  -0.043075 
  1   2.254014  -1.365839  -1.006400 
  1   3.399024  -1.624515   0.280295 
  6  -1.113220  -1.433745  -0.656146 
  1  -0.575541  -2.332289  -0.931979 
  6  -2.393256  -1.531202  -0.129580 
  1  -2.834891  -2.513727  -0.002918 
  6  -1.443013   0.993893  -0.965601 
  1  -0.957769   1.971415  -0.956603 
  1  -1.846597   0.847218  -1.971019 
  6  -2.592435   0.960873   0.043075 
  1  -3.399023   1.624515  -0.280295 
  1  -2.254014   1.365839   1.006400 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.119127    -0.424214     0.248750     1.700000     1.200000
   -0.411716    -0.117621    -0.690842     1.700000     1.200000
   -4.089584    -0.546101     0.713998     1.200000     1.200000
    0.331104    -0.106248    -1.492216     1.200000     1.200000
    1.113220     1.433745     0.656146     1.700000     1.200000
    0.411716     0.117621     0.690842     1.700000     1.200000
    0.575541     2.332289     0.931979     1.200000     1.200000
   -0.331104     0.106248     1.492216     1.200000     1.200000
    1.443013    -0.993893     0.965601     1.700000     1.200000
    1.846597    -0.847218     1.971019     1.200000     1.200000
    0.957769    -1.971415     0.956603     1.200000     1.200000
    2.393256     1.531202     0.129580     1.700000     1.200000
    2.834891     2.513727     0.002918     1.200000     1.200000
    3.119127     0.424214    -0.248750     1.700000     1.200000
    4.089584     0.546101    -0.713998     1.200000     1.200000
    2.592435    -0.960873    -0.043075     1.700000     1.200000
    2.254014    -1.365839    -1.006400     1.200000     1.200000
    3.399024    -1.624515     0.280295     1.200000     1.200000
   -1.113220    -1.433745    -0.656146     1.700000     1.200000
   -0.575541    -2.332289    -0.931979     1.200000     1.200000
   -2.393256    -1.531202    -0.129580     1.700000     1.200000
   -2.834891    -2.513727    -0.002918     1.200000     1.200000
   -1.443013     0.993893    -0.965601     1.700000     1.200000
   -0.957769     1.971415    -0.956603     1.200000     1.200000
   -1.846597     0.847218    -1.971019     1.200000     1.200000
   -2.592435     0.960873     0.043075     1.700000     1.200000
   -3.399023     1.624515    -0.280295     1.200000     1.200000
   -2.254014     1.365839     1.006400     1.200000     1.200000

--link1--
%chk=meso-TS_19_1_xp-115.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.119127  -0.424214   0.248750 
  6  -0.411716  -0.117621  -0.690842 
  1  -4.089584  -0.546101   0.713998 
  1   0.331104  -0.106248  -1.492216 
  6   1.113220   1.433745   0.656146 
  6   0.411716   0.117621   0.690842 
  1   0.575541   2.332289   0.931979 
  1  -0.331104   0.106248   1.492216 
  6   1.443013  -0.993893   0.965601 
  1   1.846597  -0.847218   1.971019 
  1   0.957769  -1.971415   0.956603 
  6   2.393256   1.531202   0.129580 
  1   2.834891   2.513727   0.002918 
  6   3.119127   0.424214  -0.248750 
  1   4.089584   0.546101  -0.713998 
  6   2.592435  -0.960873  -0.043075 
  1   2.254014  -1.365839  -1.006400 
  1   3.399024  -1.624515   0.280295 
  6  -1.113220  -1.433745  -0.656146 
  1  -0.575541  -2.332289  -0.931979 
  6  -2.393256  -1.531202  -0.129580 
  1  -2.834891  -2.513727  -0.002918 
  6  -1.443013   0.993893  -0.965601 
  1  -0.957769   1.971415  -0.956603 
  1  -1.846597   0.847218  -1.971019 
  6  -2.592435   0.960873   0.043075 
  1  -3.399023   1.624515  -0.280295 
  1  -2.254014   1.365839   1.006400 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.119127    -0.424214     0.248750     1.700000     1.150000
   -0.411716    -0.117621    -0.690842     1.700000     1.150000
   -4.089584    -0.546101     0.713998     1.200000     1.150000
    0.331104    -0.106248    -1.492216     1.200000     1.150000
    1.113220     1.433745     0.656146     1.700000     1.150000
    0.411716     0.117621     0.690842     1.700000     1.150000
    0.575541     2.332289     0.931979     1.200000     1.150000
   -0.331104     0.106248     1.492216     1.200000     1.150000
    1.443013    -0.993893     0.965601     1.700000     1.150000
    1.846597    -0.847218     1.971019     1.200000     1.150000
    0.957769    -1.971415     0.956603     1.200000     1.150000
    2.393256     1.531202     0.129580     1.700000     1.150000
    2.834891     2.513727     0.002918     1.200000     1.150000
    3.119127     0.424214    -0.248750     1.700000     1.150000
    4.089584     0.546101    -0.713998     1.200000     1.150000
    2.592435    -0.960873    -0.043075     1.700000     1.150000
    2.254014    -1.365839    -1.006400     1.200000     1.150000
    3.399024    -1.624515     0.280295     1.200000     1.150000
   -1.113220    -1.433745    -0.656146     1.700000     1.150000
   -0.575541    -2.332289    -0.931979     1.200000     1.150000
   -2.393256    -1.531202    -0.129580     1.700000     1.150000
   -2.834891    -2.513727    -0.002918     1.200000     1.150000
   -1.443013     0.993893    -0.965601     1.700000     1.150000
   -0.957769     1.971415    -0.956603     1.200000     1.150000
   -1.846597     0.847218    -1.971019     1.200000     1.150000
   -2.592435     0.960873     0.043075     1.700000     1.150000
   -3.399023     1.624515    -0.280295     1.200000     1.150000
   -2.254014     1.365839     1.006400     1.200000     1.150000

--link1--
%chk=meso-TS_19_1_xp-110.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.119127  -0.424214   0.248750 
  6  -0.411716  -0.117621  -0.690842 
  1  -4.089584  -0.546101   0.713998 
  1   0.331104  -0.106248  -1.492216 
  6   1.113220   1.433745   0.656146 
  6   0.411716   0.117621   0.690842 
  1   0.575541   2.332289   0.931979 
  1  -0.331104   0.106248   1.492216 
  6   1.443013  -0.993893   0.965601 
  1   1.846597  -0.847218   1.971019 
  1   0.957769  -1.971415   0.956603 
  6   2.393256   1.531202   0.129580 
  1   2.834891   2.513727   0.002918 
  6   3.119127   0.424214  -0.248750 
  1   4.089584   0.546101  -0.713998 
  6   2.592435  -0.960873  -0.043075 
  1   2.254014  -1.365839  -1.006400 
  1   3.399024  -1.624515   0.280295 
  6  -1.113220  -1.433745  -0.656146 
  1  -0.575541  -2.332289  -0.931979 
  6  -2.393256  -1.531202  -0.129580 
  1  -2.834891  -2.513727  -0.002918 
  6  -1.443013   0.993893  -0.965601 
  1  -0.957769   1.971415  -0.956603 
  1  -1.846597   0.847218  -1.971019 
  6  -2.592435   0.960873   0.043075 
  1  -3.399023   1.624515  -0.280295 
  1  -2.254014   1.365839   1.006400 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.119127    -0.424214     0.248750     1.700000     1.100000
   -0.411716    -0.117621    -0.690842     1.700000     1.100000
   -4.089584    -0.546101     0.713998     1.200000     1.100000
    0.331104    -0.106248    -1.492216     1.200000     1.100000
    1.113220     1.433745     0.656146     1.700000     1.100000
    0.411716     0.117621     0.690842     1.700000     1.100000
    0.575541     2.332289     0.931979     1.200000     1.100000
   -0.331104     0.106248     1.492216     1.200000     1.100000
    1.443013    -0.993893     0.965601     1.700000     1.100000
    1.846597    -0.847218     1.971019     1.200000     1.100000
    0.957769    -1.971415     0.956603     1.200000     1.100000
    2.393256     1.531202     0.129580     1.700000     1.100000
    2.834891     2.513727     0.002918     1.200000     1.100000
    3.119127     0.424214    -0.248750     1.700000     1.100000
    4.089584     0.546101    -0.713998     1.200000     1.100000
    2.592435    -0.960873    -0.043075     1.700000     1.100000
    2.254014    -1.365839    -1.006400     1.200000     1.100000
    3.399024    -1.624515     0.280295     1.200000     1.100000
   -1.113220    -1.433745    -0.656146     1.700000     1.100000
   -0.575541    -2.332289    -0.931979     1.200000     1.100000
   -2.393256    -1.531202    -0.129580     1.700000     1.100000
   -2.834891    -2.513727    -0.002918     1.200000     1.100000
   -1.443013     0.993893    -0.965601     1.700000     1.100000
   -0.957769     1.971415    -0.956603     1.200000     1.100000
   -1.846597     0.847218    -1.971019     1.200000     1.100000
   -2.592435     0.960873     0.043075     1.700000     1.100000
   -3.399023     1.624515    -0.280295     1.200000     1.100000
   -2.254014     1.365839     1.006400     1.200000     1.100000

--link1--
%chk=meso-TS_19_1_xp-105.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.119127  -0.424214   0.248750 
  6  -0.411716  -0.117621  -0.690842 
  1  -4.089584  -0.546101   0.713998 
  1   0.331104  -0.106248  -1.492216 
  6   1.113220   1.433745   0.656146 
  6   0.411716   0.117621   0.690842 
  1   0.575541   2.332289   0.931979 
  1  -0.331104   0.106248   1.492216 
  6   1.443013  -0.993893   0.965601 
  1   1.846597  -0.847218   1.971019 
  1   0.957769  -1.971415   0.956603 
  6   2.393256   1.531202   0.129580 
  1   2.834891   2.513727   0.002918 
  6   3.119127   0.424214  -0.248750 
  1   4.089584   0.546101  -0.713998 
  6   2.592435  -0.960873  -0.043075 
  1   2.254014  -1.365839  -1.006400 
  1   3.399024  -1.624515   0.280295 
  6  -1.113220  -1.433745  -0.656146 
  1  -0.575541  -2.332289  -0.931979 
  6  -2.393256  -1.531202  -0.129580 
  1  -2.834891  -2.513727  -0.002918 
  6  -1.443013   0.993893  -0.965601 
  1  -0.957769   1.971415  -0.956603 
  1  -1.846597   0.847218  -1.971019 
  6  -2.592435   0.960873   0.043075 
  1  -3.399023   1.624515  -0.280295 
  1  -2.254014   1.365839   1.006400 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.119127    -0.424214     0.248750     1.700000     1.050000
   -0.411716    -0.117621    -0.690842     1.700000     1.050000
   -4.089584    -0.546101     0.713998     1.200000     1.050000
    0.331104    -0.106248    -1.492216     1.200000     1.050000
    1.113220     1.433745     0.656146     1.700000     1.050000
    0.411716     0.117621     0.690842     1.700000     1.050000
    0.575541     2.332289     0.931979     1.200000     1.050000
   -0.331104     0.106248     1.492216     1.200000     1.050000
    1.443013    -0.993893     0.965601     1.700000     1.050000
    1.846597    -0.847218     1.971019     1.200000     1.050000
    0.957769    -1.971415     0.956603     1.200000     1.050000
    2.393256     1.531202     0.129580     1.700000     1.050000
    2.834891     2.513727     0.002918     1.200000     1.050000
    3.119127     0.424214    -0.248750     1.700000     1.050000
    4.089584     0.546101    -0.713998     1.200000     1.050000
    2.592435    -0.960873    -0.043075     1.700000     1.050000
    2.254014    -1.365839    -1.006400     1.200000     1.050000
    3.399024    -1.624515     0.280295     1.200000     1.050000
   -1.113220    -1.433745    -0.656146     1.700000     1.050000
   -0.575541    -2.332289    -0.931979     1.200000     1.050000
   -2.393256    -1.531202    -0.129580     1.700000     1.050000
   -2.834891    -2.513727    -0.002918     1.200000     1.050000
   -1.443013     0.993893    -0.965601     1.700000     1.050000
   -0.957769     1.971415    -0.956603     1.200000     1.050000
   -1.846597     0.847218    -1.971019     1.200000     1.050000
   -2.592435     0.960873     0.043075     1.700000     1.050000
   -3.399023     1.624515    -0.280295     1.200000     1.050000
   -2.254014     1.365839     1.006400     1.200000     1.050000

--link1--
%chk=meso-TS_19_1_xp-100.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.119127  -0.424214   0.248750 
  6  -0.411716  -0.117621  -0.690842 
  1  -4.089584  -0.546101   0.713998 
  1   0.331104  -0.106248  -1.492216 
  6   1.113220   1.433745   0.656146 
  6   0.411716   0.117621   0.690842 
  1   0.575541   2.332289   0.931979 
  1  -0.331104   0.106248   1.492216 
  6   1.443013  -0.993893   0.965601 
  1   1.846597  -0.847218   1.971019 
  1   0.957769  -1.971415   0.956603 
  6   2.393256   1.531202   0.129580 
  1   2.834891   2.513727   0.002918 
  6   3.119127   0.424214  -0.248750 
  1   4.089584   0.546101  -0.713998 
  6   2.592435  -0.960873  -0.043075 
  1   2.254014  -1.365839  -1.006400 
  1   3.399024  -1.624515   0.280295 
  6  -1.113220  -1.433745  -0.656146 
  1  -0.575541  -2.332289  -0.931979 
  6  -2.393256  -1.531202  -0.129580 
  1  -2.834891  -2.513727  -0.002918 
  6  -1.443013   0.993893  -0.965601 
  1  -0.957769   1.971415  -0.956603 
  1  -1.846597   0.847218  -1.971019 
  6  -2.592435   0.960873   0.043075 
  1  -3.399023   1.624515  -0.280295 
  1  -2.254014   1.365839   1.006400 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.119127    -0.424214     0.248750     1.700000     1.000000
   -0.411716    -0.117621    -0.690842     1.700000     1.000000
   -4.089584    -0.546101     0.713998     1.200000     1.000000
    0.331104    -0.106248    -1.492216     1.200000     1.000000
    1.113220     1.433745     0.656146     1.700000     1.000000
    0.411716     0.117621     0.690842     1.700000     1.000000
    0.575541     2.332289     0.931979     1.200000     1.000000
   -0.331104     0.106248     1.492216     1.200000     1.000000
    1.443013    -0.993893     0.965601     1.700000     1.000000
    1.846597    -0.847218     1.971019     1.200000     1.000000
    0.957769    -1.971415     0.956603     1.200000     1.000000
    2.393256     1.531202     0.129580     1.700000     1.000000
    2.834891     2.513727     0.002918     1.200000     1.000000
    3.119127     0.424214    -0.248750     1.700000     1.000000
    4.089584     0.546101    -0.713998     1.200000     1.000000
    2.592435    -0.960873    -0.043075     1.700000     1.000000
    2.254014    -1.365839    -1.006400     1.200000     1.000000
    3.399024    -1.624515     0.280295     1.200000     1.000000
   -1.113220    -1.433745    -0.656146     1.700000     1.000000
   -0.575541    -2.332289    -0.931979     1.200000     1.000000
   -2.393256    -1.531202    -0.129580     1.700000     1.000000
   -2.834891    -2.513727    -0.002918     1.200000     1.000000
   -1.443013     0.993893    -0.965601     1.700000     1.000000
   -0.957769     1.971415    -0.956603     1.200000     1.000000
   -1.846597     0.847218    -1.971019     1.200000     1.000000
   -2.592435     0.960873     0.043075     1.700000     1.000000
   -3.399023     1.624515    -0.280295     1.200000     1.000000
   -2.254014     1.365839     1.006400     1.200000     1.000000

--link1--
%chk=meso-TS_19_1_xp-0975.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.119127  -0.424214   0.248750 
  6  -0.411716  -0.117621  -0.690842 
  1  -4.089584  -0.546101   0.713998 
  1   0.331104  -0.106248  -1.492216 
  6   1.113220   1.433745   0.656146 
  6   0.411716   0.117621   0.690842 
  1   0.575541   2.332289   0.931979 
  1  -0.331104   0.106248   1.492216 
  6   1.443013  -0.993893   0.965601 
  1   1.846597  -0.847218   1.971019 
  1   0.957769  -1.971415   0.956603 
  6   2.393256   1.531202   0.129580 
  1   2.834891   2.513727   0.002918 
  6   3.119127   0.424214  -0.248750 
  1   4.089584   0.546101  -0.713998 
  6   2.592435  -0.960873  -0.043075 
  1   2.254014  -1.365839  -1.006400 
  1   3.399024  -1.624515   0.280295 
  6  -1.113220  -1.433745  -0.656146 
  1  -0.575541  -2.332289  -0.931979 
  6  -2.393256  -1.531202  -0.129580 
  1  -2.834891  -2.513727  -0.002918 
  6  -1.443013   0.993893  -0.965601 
  1  -0.957769   1.971415  -0.956603 
  1  -1.846597   0.847218  -1.971019 
  6  -2.592435   0.960873   0.043075 
  1  -3.399023   1.624515  -0.280295 
  1  -2.254014   1.365839   1.006400 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.119127    -0.424214     0.248750     1.700000     0.975000
   -0.411716    -0.117621    -0.690842     1.700000     0.975000
   -4.089584    -0.546101     0.713998     1.200000     0.975000
    0.331104    -0.106248    -1.492216     1.200000     0.975000
    1.113220     1.433745     0.656146     1.700000     0.975000
    0.411716     0.117621     0.690842     1.700000     0.975000
    0.575541     2.332289     0.931979     1.200000     0.975000
   -0.331104     0.106248     1.492216     1.200000     0.975000
    1.443013    -0.993893     0.965601     1.700000     0.975000
    1.846597    -0.847218     1.971019     1.200000     0.975000
    0.957769    -1.971415     0.956603     1.200000     0.975000
    2.393256     1.531202     0.129580     1.700000     0.975000
    2.834891     2.513727     0.002918     1.200000     0.975000
    3.119127     0.424214    -0.248750     1.700000     0.975000
    4.089584     0.546101    -0.713998     1.200000     0.975000
    2.592435    -0.960873    -0.043075     1.700000     0.975000
    2.254014    -1.365839    -1.006400     1.200000     0.975000
    3.399024    -1.624515     0.280295     1.200000     0.975000
   -1.113220    -1.433745    -0.656146     1.700000     0.975000
   -0.575541    -2.332289    -0.931979     1.200000     0.975000
   -2.393256    -1.531202    -0.129580     1.700000     0.975000
   -2.834891    -2.513727    -0.002918     1.200000     0.975000
   -1.443013     0.993893    -0.965601     1.700000     0.975000
   -0.957769     1.971415    -0.956603     1.200000     0.975000
   -1.846597     0.847218    -1.971019     1.200000     0.975000
   -2.592435     0.960873     0.043075     1.700000     0.975000
   -3.399023     1.624515    -0.280295     1.200000     0.975000
   -2.254014     1.365839     1.006400     1.200000     0.975000

--link1--
%chk=meso-TS_19_1_xp-095.chk
%nproc=24
%mem=48GB
#p uwb97xd guess=mix def2tzvp 
# scrf=(iefpcm,solvent=Cyclohexane,read) nosym guess=only pop=none 

endo-2
xp-pcm for eta=3

0 1
  6  -3.119127  -0.424214   0.248750 
  6  -0.411716  -0.117621  -0.690842 
  1  -4.089584  -0.546101   0.713998 
  1   0.331104  -0.106248  -1.492216 
  6   1.113220   1.433745   0.656146 
  6   0.411716   0.117621   0.690842 
  1   0.575541   2.332289   0.931979 
  1  -0.331104   0.106248   1.492216 
  6   1.443013  -0.993893   0.965601 
  1   1.846597  -0.847218   1.971019 
  1   0.957769  -1.971415   0.956603 
  6   2.393256   1.531202   0.129580 
  1   2.834891   2.513727   0.002918 
  6   3.119127   0.424214  -0.248750 
  1   4.089584   0.546101  -0.713998 
  6   2.592435  -0.960873  -0.043075 
  1   2.254014  -1.365839  -1.006400 
  1   3.399024  -1.624515   0.280295 
  6  -1.113220  -1.433745  -0.656146 
  1  -0.575541  -2.332289  -0.931979 
  6  -2.393256  -1.531202  -0.129580 
  1  -2.834891  -2.513727  -0.002918 
  6  -1.443013   0.993893  -0.965601 
  1  -0.957769   1.971415  -0.956603 
  1  -1.846597   0.847218  -1.971019 
  6  -2.592435   0.960873   0.043075 
  1  -3.399023   1.624515  -0.280295 
  1  -2.254014   1.365839   1.006400 

qrep pcmdoc geomview  nodis nocav g03defaults tsare=0.075 
nsfe=28
NVESolv=36  SolvMW=84.1595   Rsolv=2.815 
eps=2.0165 RhoS=0.7781

   -3.119127    -0.424214     0.248750     1.700000     0.950000
   -0.411716    -0.117621    -0.690842     1.700000     0.950000
   -4.089584    -0.546101     0.713998     1.200000     0.950000
    0.331104    -0.106248    -1.492216     1.200000     0.950000
    1.113220     1.433745     0.656146     1.700000     0.950000
    0.411716     0.117621     0.690842     1.700000     0.950000
    0.575541     2.332289     0.931979     1.200000     0.950000
   -0.331104     0.106248     1.492216     1.200000     0.950000
    1.443013    -0.993893     0.965601     1.700000     0.950000
    1.846597    -0.847218     1.971019     1.200000     0.950000
    0.957769    -1.971415     0.956603     1.200000     0.950000
    2.393256     1.531202     0.129580     1.700000     0.950000
    2.834891     2.513727     0.002918     1.200000     0.950000
    3.119127     0.424214    -0.248750     1.700000     0.950000
    4.089584     0.546101    -0.713998     1.200000     0.950000
    2.592435    -0.960873    -0.043075     1.700000     0.950000
    2.254014    -1.365839    -1.006400     1.200000     0.950000
    3.399024    -1.624515     0.280295     1.200000     0.950000
   -1.113220    -1.433745    -0.656146     1.700000     0.950000
   -0.575541    -2.332289    -0.931979     1.200000     0.950000
   -2.393256    -1.531202    -0.129580     1.700000     0.950000
   -2.834891    -2.513727    -0.002918     1.200000     0.950000
   -1.443013     0.993893    -0.965601     1.700000     0.950000
   -0.957769     1.971415    -0.956603     1.200000     0.950000
   -1.846597     0.847218    -1.971019     1.200000     0.950000
   -2.592435     0.960873     0.043075     1.700000     0.950000
   -3.399023     1.624515    -0.280295     1.200000     0.950000
   -2.254014     1.365839     1.006400     1.200000     0.950000

