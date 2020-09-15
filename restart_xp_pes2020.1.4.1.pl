#------------------------------------------------------------------------------
#!/usr/bin/perl
#
# ======================================================================
#                     >>>>>>>>>>>>>><<<<<<<<<<<<<<
#                     >>>>>>>>>>        <<<<<<<<<<
#                     >>>>>>>>> XP-PCM17 <<<<<<<<<
#                     >>>>>>>>>>        <<<<<<<<<<
#                     >>>>>>>>>>>>>><<<<<<<<<<<<<<
#                                 
#                            Roberto Cammi, 
#                 Department of Chemical Science (SCVSA)
#                         University of Parma
#======================================================================
#
# XP-PCM17: a Perl script to conduct XP-PCM calculations for
# studying chemical reactions under pressure [1,2]. 
# XP-PCM17 requires Gaussian-16 and Mathematica/Python [4]/[5] as 
# external packages.
# XP-PCM7 allows both mono- or bimolecular  reactive
# processes in a specified external non-polar medium (the solvent).
# Note: The presents script is still under development and its use is
#       confidential. 
#
#     
# References:
#         [1] R. Cammi, J. Comp. Chem., 36, 2246-2259 (2015)
#         [2] B. Chen. R. Hoffmann and R. Cammi, 
#            Angew. Chem. Int. Ed. 56, 11126-11142 (2017) and related
#                       Supporting Information (SI)
#         [3] R. Cammi, B. Chen, M. Rahm, J. Comp. Chem. 39, 2243 (2018)
#         [4] R. Cammi, J. Chem. Phys.,150,164122 (2019)
#         [5] M. Rahm, R. Cammi, N.W. Ashcroft. R. Hoffman, JACS,   (2020)
#   
# Note: The presents script is still under development and its use is for 
# now  confidential.    
# 
#------------------------------------------------------------------------------                     
# Usage: 
#       restart_xp_pes2020.1.4.1.pl  ns fn_1.com fn_2.com ... fn_ns.com  \
#                   f0 f1 ... fn \
#                   solvent  eta cavity fcav tsare molec mrng rest   
#
# The classes of input arguments:
# 1)  ns: 
#     ns .... the number of the  selected geometries of the reactive system 
# 2)  {fn_1. com,... fn_ns.com}: 
#     A set of ns Gaussian-16 input files. 
#     The input files must be as for a ground state (HF/DFT) calculation in the 
#     gas phase. Furthermore, they must have: 
#     (a) The route section must start with the string '#p', to require 
#     additional printing in the output files of Gaussian 16. 
#     (b) The geometry must be specified in Cartesian  coordinates, with at 
#     least one blank space before the atomic number/symbol of the nuclei.
# 
# 3)  {f_1, f_2, ..., f_nf}, with f_{n+1}<f_{n} and nf=> 6:
#     A set of decreasing cavity scaling factors for the shrinking of the 
#     molecular cavities (of Solvende Excluded Surface type) used in the calculation 
#     of the electronic energy, G_ef({i};f_n).
#     The cavity scaling factors modulate the van der Waals atomic sphere's radii 
#     centered on the atomic nuclei of the reactive system: R_i=f*R_i^vdw. The reference 
#     van der Waals atomic radii {R_i^vdw} are taken from Bondi (R_h=1.2 Ang.;
#     R_C=1.7 Ang; R_O=1.52 Ang.).  
#     The first element of the set (f_1) must have a fixed value of f_1=1.2, and the 
#     last value is bounded to f_{nf}=>0.90.  The XP-PCM calculations can span a pressure
#     range of around 1-15/20 GPa.
#
# 4) 'solvent': 
#    'solvent' is a string of characters that specifies the external medium.
#     Values for 'solvent' are:
#               C6H12, for cyclohexane
#               C6H6, for benzene
#               Ar, for argon
# 5) eta:  
#    eta specifies the value for the 'semi-empirical' Pauli repulsion parameter. 
#    The reference value is eta=6.
#                             
# 6) cavity: 
#    cavity is a string of character defining the typology of the cavity used 
#    for the calculation of the cavitation energy, G_cav. 
#    Values for cavity are:
#               vdw, for a van der Waals type cavity
#               ses, for a surface excluded type solvent cavity.
#               Reference value: vdw
# 
# 7) fcav: 
#    fcav ... the fixed scaling factor for the calculation of G_cav. 
#    Reference value fcav=1.2
#
# 8) tsare: 
#    tsare specifies the mean area (Ang^2) of the of tesserae in which the 
#    surfaces of all the molecular cavities are partitioned.
#    Reference value= 0.075
#
# 9) molec:
#    molec = 0 ... mono-molecular process
#          = 1 ... bi-molecular process (A+B-> AB).
#                  Where the separated reactant A and B are the first two 
#                  structures (fn_1.com,fn_2.com) in set of Gaussian-16 input files
# 
# 10) mrgn:
#     mrgn = 'math' Murnaghan fitting with Mathematica
#           = 'py'   Murnaghan fitting with Python
#  
# 11) rest:
#     rest ... specfies the point of restart of a XP-PCM job that was suddenly stopped. 
#              Assuming the job was interrupted during the step of calculation of the 
#              electronic energy of the selected molecular structures, the point of restart is 
#              the serial number  of the current molecular structure at the time of stop.
#              'rest' can be determined from the Vp.dat file of the interrupted job.  This 
#              file presents the serial numbers of all the molecular structures whose electronic 
#              energy calculations have been completed, and the value of 'rest' must 
#              correspond to the serial number after that of the last completed structure.
#              (For example, if Vp.dat looks like with three completed structures calculations we have to set rest = 4)     
#              For the restart to be successful, all output files *Ger.log of the electronic energy calculations of 
#              the completed structures are required. 
#                   
#
# The main output of XP-PES script is the the potential energy profile
# of the reaction path along the selected geometries. The output file
# has the fixed name xp_pes.log
#
# Output files: 
# {fn*_xp-Ger.com/log}: gaussian 16 calculations of the electronic energy
#
# {fn*_xp-Gcav.comlog}: gaussian 16 calculations of the cavitation energy
#
# xp_pes.log: The main xp-pcm outcomes collected in two tables: 
#
# Table I:  The electronic (G_er;i,p), cavitation (G_cav;i,p) and total 
#           [(G_tot;i,p)= (G_er;i,p)+(G_cav;i,p)] energies of the reactive 
#           system as a function of the pressure.
# Table II: The potential energy profile (DG_tot;p) as a function of the 
#           pressure. The first molecular geometry is used as reference state.    
#-----------------------------------------------------------------------------
#
#$n_par=6;
#$n_par=8;
#RC120820
$n_par=9;
print '$#ARGV=',$#ARGV,"\n".' $n_par=',$n_par,"\n";
$n_struct= $ARGV[0]; #print '$n_struct=',$n_struct;
#
# read cavity scaling factors.
#
#RC250217--
$nfm1=$#ARGV-($n_par+$n_struct); print '$nfm1=', $nfm1, "\n";
#for($i=0; $i<=$#ARGV-(6+$n_struct);$i++){$factor[$i]=$ARGV[$i+(1+$n_struct)]};
for($i=0; $i<=$nfm1;$i++){$factor[$i]=$ARGV[$i+(1+$n_struct)]};
##RC250217--
$factor0=$factor[0];
#
# Reading Solvent,Eta,Cavity type, Radii factor for Gcav,Tsare
#
$npointer=$n_struct+$nfm1+2;
#$solv=$ARGV[$#ARGV-4];print $solv,"\n";
#$eta=$ARGV[$#ARGV-3]; print $eta,"\n";
#$cavity=$ARGV[$#ARGV-2]; print $cavity,"\n";
#$frad=$ARGV[$#ARGV-1]; print $frad,"\n";
$solv=$ARGV[$npointer];print $solv,"\n";$npointer++;
$eta=$ARGV[$npointer]; print $eta,"\n";$npointer++;
$cavity=$ARGV[$npointer]; print $cavity,"\n";$npointer++;
$frad=$ARGV[$npointer]; print $frad,"\n";$npointer++;
if( $cavity =~ /vdw/){$noaddsph='noaddsph '};
if( $cavity =~ /ses/){$noaddsph=' '};
print $noaddsph," \n";
$tsare=$ARGV[$npointer]; print $tsare,"\n";$npointer++;
$nmolec=$ARGV[$npointer];$npointer++;
if($nmolec == 0){print 'mono-molecular process: A_1->A_2->....',"\n"};
if($nmolec == 1){print 'bi-molecular process: A+B-> AB_1->....',"\n"};
$mrgn=$ARGV[$npointer]; $npointer++;
if($mrgn =~ /math/){print 'Mathematica',"\n"};
if($mrgn =~ /py/){print 'Python',"\n"};
$rest = $ARGV[$npointer];
#
# Physical constants and conversion factors
#
$avogadro = 0.6022140857;
$cal2J = 4.184;
$Eh2J = 4.35974417;
$Eh_o_Ang3_to_GPa =  $Eh2J*10**3;
$Eh2kcal_mol=$Eh2J/$cal2J*$avogadro*10**3 ;
print $avogadro, " \n";
print $Eh_o_Ang3_to_GPa, " \n";
print $Eh2kcal_mol, " \n";
#
# Physical parameters of the solvent
#
# Cyclohexane 
#
if($solv =~ /C6H12/){
$solvent='Cyclohexane';
$eps_0 = 2.0165; 	# dielectric constant
$rho_0 = 0.7781; 	# density (gr/cm^3)
$mw = 84.1595; 	 	# molar mass (gr/mol)
$nves = 36;  	 	# number of valence electrons 
$rsolv = 2.815; 	# molecular radius (Ang.)
};
#
# Benzene
# 
if($solv =~ /C6H6/){
$solvent='Benzene';
$eps_0 = 2.2706;		# dielectric constant
$rho_0 = 0.8756;	# density (gr/cm^3)
$mw = 78.1118;		# molar mass (gr/mol)
$nves = 30;  		# number of valence electrons 
$rsolv = 2.63; 	# molecular radius (Ang.)
};
#
# Argon
# 
if($solv =~ /Ar/){
$solvent='Argon';
$eps_0 = 1.43; 	# dielectric constant
$rho_0 = 1.3954;# density (gr/cm^3) at b.p. https://en.wikipedia.org/wiki/Argon
$mw = 39.948; 	# molar mass (gr/mol)
$nves = 8;  	# number of valence electrons 
#$rsolv = 1.875;	# molecular radius (Ang.)
$rsolv = 1.705;	# molecular radius (Ang.)
};
print $solvent, " \n";
Print_1();
#$OUT= 'xp_pes.log';     
#
#close(OUT,$OUT);
#
#----------------------------------------------------
#
# Defining intermediate files
#
#----------------------------------------------------
if($rest == 0) {open(FNI3, ">properties.dat");
close(FNI3, "properties.dat");
open(FNI3, ">Vp.dat");
close(FNI3, "Vp.dat")}
#------------------------------------------------------------------------------
#
# Step 1: Computinthe cavity (SES) volume V_{c}(f) as a function of 
#         the cavity scaling factor f: V_{c}(f0)...V_c(fn) 
#
#-----------------------------------------------------------------------------
Eps_Rhos();
#-----------------------------------------------------------------------------
#
# Step 3: Computing the electronic energy G_er(f) as a function of the 
#         cavity scaling factor f. 
#
#-----------------------------------------------------------------------------
if ($rest == 0){$kstart = 1}
 else {$kstart = $rest;
my $cmd = 'rm Vp.dat ';
system($cmd)};
ElectronicEnergy();
#------------------------------------------------------------------------------
#
#  Step 5: Computing the mean value of the pressure 
#
#------------------------------------------------------------------------------
MeanPressure();
#-----------------------------------------------------------------------------
#
# Step 6 and 7: Cavitation, G_cav (p), and total G_tot(p) energies.
#
#-----------------------------------------------------------------------------
CavitationEnergy();
#------------------------------------------------------------------------------
#
# Step 8: Computing the reaction energy profile, DG_tot(ns;p), as a function
#         of the pressure p:
#
#         DG_tot(ns;p)= G_(tot)(ns;p)-[G_tot(ref;p)+G_tot(ref;p)]
# 
# 
#-----------------------------------------------------------------------------
PotentialEnergyProfile();
#
#
#
print "\n rm *.chk";
my $cmd = 'rm *.chk ';
system($cmd);
my $cmd = 'rm *scrfile* ';
system("rm *scrfile*");
# 
#
sub PotentialEnergyProfile{
 $i_fact = $nfm1+1;
if($nmolec == 0){
  $n_ref =  1;
#
for($i=1; $i<=$i_fact;$i++){ 
  for($j=1; $j<=$n_struct;$j++){
    $dgtot[$i][$j] = ($gtot_mat[$j][$i]-$gtot_mat[$n_ref][$i])*$Eh2kcal_mol;
    $dgcav[$i][$j] = $cav_mat[$j][$i]-$cav_mat[$n_ref][$i];
    $dger[$i][$j] = ($ger_mat[$j][$i]-$ger_mat[$n_ref][$i])*$Eh2kcal_mol;
                               }
                            }
#
# Reaction volume profile
#
$dV[1]=0.0;
#
for($j=2; $j<=$n_struct;$j++)
   { open(OUT3, ">pGtot.dat");
     for($i=1; $i<=$i_fact;$i++){print OUT3 $p[$i]," ",$cal2J*$dgtot[$i][$j],"\n"};
     close(OUT3);
#     system("math -script LinReg.m > lr.out");
if($mrgn =~ /math/){system("math -script LinReg.m > lr.out")};
if($mrgn =~ /py/){system("python LinReg.py > lr.out")};
     open(IN,'lr.out'); 
#     print 'structure: ',$j,"\n";
while($line = <IN>) {
        if($mrgn =~ /math/){
          @array1= split(/\s+/,$line);
          print 'math=',$array1[1];
          @array2= split(/\*/,$array1[2]);
#          $string = join( " ", $array1[1], $array2[0]);
         $string = $array1[1].$array2[0];
#          $dV[$j] = $array2[0]};
          $dV[$j] = $string};  
        if($mrgn =~ /py/){ if($line=~ /Slope/){
          @array1= split(/\,/,$line);
          $dV[$j] = $array1[1]}              
                         }
                     }
    }
#
# Reference pressure, p0,  for the analysis of the effect of the pressure DDG_tot(x;p)= DG_(tot)(x;p)-DG_(tot)(x;p0)
# 1 -> p at f=1.2
#
  $n_p0 = 1;
#
for($i=1; $i<=$i_fact;$i++){ 
  for($j=1; $j<=$n_struct;$j++){
    $ddgtot[$i][$j] = $dgtot[$i][$j]-$dgtot[$n_p0][$j];
    $ddgcav[$i][$j] = $dgcav[$i][$j]-$dgcav[$n_p0][$j];
    $ddger[$i][$j] = $dger[$i][$j]-$dger[$n_p0][$j];
                               }
                            }
#
print OUT "\n\n";
print OUT "----------------------------------------------------------------\n\n";
print OUT "Table II:
Potential energy profile DG_tot, (kcal/mol) as a function of the pressure, p (GPa). 
DG_er and DG_cav are the  electronic and cavitation contributions, respectively.\n\n";
for($i=1; $i<=$i_fact;$i++){ print OUT "p[GPa]= "; printf  OUT "%8.2f \n", $p[$i];
 printf OUT " %8s %18s %18s %18s\r\n", 'Structure','DG_tot[kcal/mol] ','DG_er[kcal/mol] ','DG_cav[kcal/mol] ';
  for($j=1; $j<=$n_struct;$j++){
  @array1=split(/\./,$ARGV[$j]);
#  printf OUT "%5d %18.1f %18.1f %18.1f %18.1f %18.1f %18.1f\r\n", $j,$dgtot[$i][$j],$dger[$i][$j],$dgcav[$i][$j],$ddgtot[$i][$j],$ddger[$i][$j],$ddgcav[$i][$j];
  printf OUT "%8s %18.1f %18.1f %18.1f \r\n", $array1[0],$dgtot[$i][$j],$dger[$i][$j],$dgcav[$i][$j];
                              }
                              };
print OUT "\n\n";
print OUT "----------------------------------------------------------------\n\n";
print OUT "Table III:
Analysis of DG_tot as a function of the pressure, p:
DDG_tot({i},p)=DG_tot({i},p)-DG_tot({i},p=";
printf  OUT "%4.1f", $p[1];
print OUT "GPa)
DDG_er and DDG_cav are the electronic and cavitation contributions, respectively.\n\n";
for($i=1; $i<=$i_fact;$i++){ print OUT "p[GPa]= "; printf  OUT "%6.1f \n", $p[$i];
 printf OUT " %8s %18s %18s %18s \r\n", 'Structure','DDG_tot[kcal/mol]','DDG_er[kcal/mol]','DDG_cav[kcal/mol] ';
  for($j=1; $j<=$n_struct;$j++){
  @array1=split(/\./,$ARGV[$j]);
#  printf OUT "%5d %18.1f %18.1f %18.1f %18.1f %18.1f %18.1f\r\n", $j,$dgtot[$i][$j],$dger[$i][$j],$dgcav[$i][$j],$ddgtot[$i][$j],$ddger[$i][$j],$ddgcav[$i][$j];
  printf OUT "%8s %18.1f %18.1f %18.1f \r\n", $array1[0],$ddgtot[$i][$j],$ddger[$i][$j],$ddgcav[$i][$j];
                              }
                              };
print OUT "\n\n";
print OUT "----------------------------------------------------------------\n\n";
print OUT "Table IV:
Reaction volume profile DV({i}) (cm^3/mol)\n";
 printf OUT " %16s %18s \r\n", 'Structure','DV_tot[cm^3/mol] ';
@array1=split(/\./,$ARGV[1]);
printf OUT "%16s %18.1f  \r\n", $array1[0],$dV[1];
  for($j=2; $j<=$n_struct;$j++){
  @array1=split(/\./,$ARGV[$j]);
  printf OUT "%16s %18.1f \r\n", $array1[0],$dV[$j];
                              }
 };
#----------------------------------------------------------------------------------
if($nmolec == 1){
  $n_ref =  1;
#
  for($i=1; $i<=$i_fact;$i++){ 
    $dgtot[$i][2] = (($gtot_mat[$n_ref][$i]+$gtot_mat[$n_ref+1][$i])-($gtot_mat[$n_ref][$i]+$gtot_mat[$n_ref+1][$i]))*$Eh2kcal_mol;
    $dgcav[$i][2] = ($cav_mat[$n_ref][$i]+$cav_mat[$n_ref+1][$i])-($cav_mat[$n_ref][$i]+$cav_mat[$n_ref+1][$i]);
    $dger[$i][2] = (($ger_mat[$n_ref][$i]+$ger_mat[$n_ref+1][$i])-($ger_mat[$n_ref][$i]+$ger_mat[$n_ref+1][$i]))*$Eh2kcal_mol;
  for($j=3; $j<=$n_struct;$j++){
    $dgtot[$i][$j] = ($gtot_mat[$j][$i]-($gtot_mat[$n_ref][$i]+$gtot_mat[$n_ref+1][$i]))*$Eh2kcal_mol;
    $dgcav[$i][$j] = $cav_mat[$j][$i]-($cav_mat[$n_ref][$i]+$cav_mat[$n_ref+1][$i]);
    $dger[$i][$j] = ($ger_mat[$j][$i]-($ger_mat[$n_ref][$i]+$ger_mat[$n_ref+1][$i]))*$Eh2kcal_mol;
                               }
                            }
#
#
# Reaction volume profile
#
$dV[2]=0.0;
#
for($j=3; $j<=$n_struct;$j++)
   { open(OUT3, ">pGtot.dat");
     for($i=1; $i<=$i_fact;$i++){print OUT3 $p[$i]," ",$cal2J*$dgtot[$i][$j],"\n"};
     close(OUT3);
#     system("math -script LinReg.m > lr.out");
if($mrgn =~ /math/){system("math -script LinReg.m > lr.out")};
if($mrgn =~ /py/){system("python LinReg.py > lr.out")};
     open(IN,'lr.out'); 
#     print 'structure: ',$j,"\n";
while($line = <IN>) {
        if($mrgn =~ /math/){
          @array1= split(/\s+/,$line);
          @array2= split(/\*/,$array1[2]);
          $string = $array1[1].$array2[0];
          $dV[$j] = $string}; 
        if($mrgn =~ /py/){ if($line=~ /Slope/){
          @array1= split(/\,/,$line);
          $dV[$j] = $array1[1]}              
                         }
                     }
    }
#
# Reference pressure, p0,  for the analysis of the effect of the pressure DDG_tot(x;p)= DG_(tot)(x;p)-DG_(tot)(x;p0)
# 1 -> p at f=1.2
#
  $n_p0 = 1;
#
for($i=1; $i<=$i_fact;$i++){ 
#  for($j=1; $j<=$n_struct;$j++){
  for($j=2; $j<=$n_struct;$j++){
    $ddgtot[$i][$j] = $dgtot[$i][$j]-$dgtot[$n_p0][$j];
    $ddgcav[$i][$j] = $dgcav[$i][$j]-$dgcav[$n_p0][$j];
    $ddger[$i][$j] = $dger[$i][$j]-$dger[$n_p0][$j];
                               }
                            }
#
print OUT "\n\n";
print OUT "----------------------------------------------------------------\n\n";
print OUT "Table II:
Potential energy profile DG_tot({i},p), (kcal/mol) as a function of the pressure, p (GPa).
DG_er and DG_cav are the  electronic and cavitation contributions, respectively.\n\n";
for($i=1; $i<=$i_fact;$i++){ print OUT "p[GPa]= "; printf  OUT "%6.1f \n", $p[$i];
 printf OUT " %16s %18s %18s %18s\r\n", 'Structure','DG_tot[kcal/mol] ','DG_er[kcal/mol] ','DG_cav[kcal/mol] ';
 $j=2;
@array1=split(/\./,$ARGV[1]);
@array2=split(/\./,$ARGV[2]);
$string = join( "+", $array1[0], $array2[0]);
 printf OUT "%16s %18.1f %18.1f %18.1f \r\n", $string,$dgtot[$i][$j],$dger[$i][$j],$dgcav[$i][$j];
#  for($j=1; $j<=$n_struct;$j++){
  for($j=3; $j<=$n_struct;$j++){
  @array1=split(/\./,$ARGV[$j]);
  printf OUT "%16s %18.1f %18.1f %18.1f \r\n", $array1[0],$dgtot[$i][$j],$dger[$i][$j],$dgcav[$i][$j];
                              }
                              };
print OUT "\n\n";
print OUT "----------------------------------------------------------------\n\n";
print OUT "Table III:
Analysis of DG_tot as a function of the pressure, p:
DDG_tot({i},p)=DG_tot({i},p)-DG_tot({i},p=";
printf  OUT "%4.1f", $p[1];
print OUT "GPa)
DDG_er and DDG_cav are the electronic and cavitation contributions, respectively.\n\n";
for($i=2; $i<=$i_fact;$i++){ print OUT "p[GPa]= "; printf  OUT "%6.1f \n", $p[$i];
 printf OUT " %16s %18s %18s %18s\r\n", 'Structure','DDG_tot[kcal/mol] ','DDG_er[kcal/mol] ','DDG_cav[kcal/mol] ';
 $j=2;
@array1=split(/\./,$ARGV[1]);
@array2=split(/\./,$ARGV[2]);
$string = join( "+", $array1[0], $array2[0]);
 printf OUT "%16s %18.1f %18.1f %18.1f \r\n", $string,$ddgtot[$i][$j],$ddger[$i][$j],$ddgcav[$i][$j];
#  for($j=1; $j<=$n_struct;$j++){
  for($j=3; $j<=$n_struct;$j++){
  @array1=split(/\./,$ARGV[$j]);
  printf OUT "%16s %18.1f %18.1f %18.1f \r\n", $array1[0],$ddgtot[$i][$j],$ddger[$i][$j],$ddgcav[$i][$j];
                              }
                              };
print OUT "\n\n";
print OUT "----------------------------------------------------------------\n\n";
print OUT "Table IV:
Reaction volume profile DV({i}) (cm^3/mol)\n";
 printf OUT " %16s %18s \r\n", 'Structure','DV_tot[cm^3/mol] ';
@array1=split(/\./,$ARGV[1]);
@array2=split(/\./,$ARGV[2]);
$string = join( "+", $array1[0], $array2[0]);
printf OUT "%16s %18.1f  \r\n", $string,$dV[2];
  for($j=3; $j<=$n_struct;$j++){
  @array1=split(/\./,$ARGV[$j]);
  printf OUT "%16s %18.1f \r\n", $array1[0],$dV[$j];
                              }
};
#----------------------------------------------------------------------------------
#
print OUT "\n----------------------------------------------------------------\n";
$datestring = localtime();
print OUT "$datestring\n\n";
}
#
#
#
sub Print_1 {
$OUT= 'xp_pes.log';     
open(OUT,">$OUT"); 
print OUT "\n
================================================================\n
                  >>>>>>>>>>>>>><<<<<<<<<<<<<<
                  >>>>>>>>>>        <<<<<<<<<<
                  >>>>>>>>> XP-PCM17 <<<<<<<<<
                  >>>>>>>>>>        <<<<<<<<<<
                  >>>>>>>>>>>>>><<<<<<<<<<<<<<

                          Roberto Cammi
                    University of Parma,ITALY

The XP-PCM method for studying reactions under pressure [1,2].
XP-PCM17 uses Gaussian 16 [3] and Mathematica(or Python)[4,5] 
as external packages.

References:
         [1] R. Cammi, J. Comp. Chem., 36, 2246-2259 (2015)
         [2] B. Chen. R. Hoffmann and R. Cammi, 
            Angew. Chem. Int. Ed. 56, 11126-11142 (2017);
         [3] Gaussian 16, M. J. Frisch et al, Gaussian, 
             Inc., Wallingford CT, 2016.
         [4] Wolfram Research, Inc., Mathematica, Version 12.0, 
             Champaign, IL (2020).
         [5] http://www.python.org\n
================================================================\n\n";                               
#
$datestring = localtime();
print OUT "$datestring\n\n";
print OUT "----------------------------------------------------------------\n";
if($nmolec == 0){print OUT "XP-PCM energy profile for monomolecular reactions\n\n"};
if($nmolec == 1){print OUT "XP-PCM energy profile for bimolecular reactions\n\n"};
print OUT "Number of molecular structures: ",$n_struct,"\n\n";
print OUT "Structures: ";
for($k=1; $k<=$n_struct;$k++){
#
$oldfile = $ARGV[$k];
@array1=split(/\./,$oldfile);
print OUT $array1[0],' '}; 
#
print OUT "

XP-PCM parameters:\n",
"eta=",$eta,"\n",
"tsare=",$tsare,"\n",
"Cavity(Ger)=(ses; f=",$factor[0],"-",$factor[$nfm1],")\n",
"Cavity(Gcav)=(",$cavity,"; f=",$frad,")\n",
"Solvent=",$solvent,"\n",
"         dielectric constant=",$eps_0,"\n",
"         density (gr/cm^3)=",$rho_0,"\n",
"         molar mass (gr/mol)=",$mw,"\n",
"         number of valence electrons=",$nves,"\n",
"         molecular radius(Ang.)=", $rsolv,"\n";
print OUT "----------------------------------------------------------------\n";
close(OUT,$OUT);   
}
#
sub Eps_Rhos {
#----------------------------------------------------
#
# Template input files for xp-pcm 
#
#----------------------------------------------------
for($k=1; $k<=$n_struct;$k++){
#
$oldfile = $ARGV[$k];
@array1=split(/\./,$oldfile);
$scrfile = $array1[0].'_scrfile'; 
$n_blanks = 0;
$n_atoms = 0;
open(OF, $oldfile);
open(SCR, ">$scrfile");
#
print SCR "%chk=\n";
# Reading a gas phase input file
$n_blanks = 0;
$n_atoms = 0;
$n_mult = 0;
#RC250217
$one = 1.0;
#RC250217
while ($line = <OF>) {
  if($line=~ /\S/) 
    {if(!($line=~/chk/)){print SCR $line};
      if($line =~ /#/) {print SCR "# scrf=(iefpcm,solvent=",$solvent,",read) nosym guess=only pop=none \n"};
      if($n_blanks == 1) {print SCR "xp-pcm for eta=",$eta,"\n"};
      if($n_blanks == 2)
        {if($n_mult == 0){$n_mult++}
         else 
          {chop($line);
           $n_atoms++; 
           @entries = split(/\s+/,$line);
           $atoms[$n_atoms] = $entries[1];
           $x[$n_atoms] = $entries[2];
           $y[$n_atoms] = $entries[3];
           $z[$n_atoms] = $entries[4];
           if ( $atoms[$n_atoms] == 1) {$radius[$n_atoms] = 1.2};
           if ( $atoms[$n_atoms] =~  /H/) {$radius[$n_atoms] = 1.2};
           if ( $atoms[$n_atoms] == 6) {$radius[$n_atoms] = 1.7};
           if ( $atoms[$n_atoms] =~  /C/) {$radius[$n_atoms] = 1.7};
           if ( $atoms[$n_atoms] == 7) {$radius[$n_atoms] = 1.55};
           if ( $atoms[$n_atoms] =~  /N/) {$radius[$n_atoms] = 1.55};
           if ( $atoms[$n_atoms] == 8) {$radius[$n_atoms] = 1.52};
           if ( $atoms[$n_atoms] =~  /O/) {$radius[$n_atoms] = 1.52};
           if ( $atoms[$n_atoms] == 9) {$radius[$n_atoms] = 1.47};
           if ( $atoms[$n_atoms] =~  /F/) {$radius[$n_atoms] = 1.47};
#           if ( $atoms[$n_atoms] == 17) {$radius[$n_atoms] = 1.75};
#           if ( $atoms[$n_atoms] =~  /Cl/) {$radius[$n_atoms] = 1.75};
          } 
        }
     }
   else
       {print SCR $line;$n_blanks++;
        if($n_blanks == 3){$n_sfe=$n_atoms;
                           print SCR "qrep pcmdoc geomview  nodis nocav g03defaults tsare=",$tsare," \n";
                           print SCR "nsfe=",$n_sfe,"\n";
                           print SCR "NVESolv=",$nves," "," SolvMW=",$mw," ","  Rsolv=",$rsolv," \n";
                           print SCR "eps=$eps_0 RhoS=$rho_0\n";
                           print SCR "\n";
                           for($i=1; $i<=$n_sfe;$i++)
                            {printf SCR "%5d %12.6f %12.6f\r\n",
                             $i,$radius[$i]*$frad,$one
                            };
                           print SCR "\n"}}
  }
#}
#
#sub Eps_Rhos{
@array1=split(/\./,$oldfile);
$newfile = $array1[0].'_xp-Vc'.'.com'; 
open(NF, ">$newfile");
#
#for($j=0; $j<=$#ARGV-(7+$n_struct);$j++){
for($j=0; $j<=($nfm1-1);$j++){
open(SCR, $scrfile);
@array2=split(/\./,$factor[$j]);
$chk= $array1[0].'_xp-'.$array2[0].$array2[1].'.chk';
$n_blanks=0;
while ($line = <SCR>) {
  if($line=~ /\S/) { 
      if(!($n_blanks >= 4)) {if($line=~/chk/){chop($line); print NF $line,$chk,"\n"}; if(!($line=~/chk/)){print NF $line}};
      if($n_blanks == 4){ for($i=1; $i<=$n_sfe;$i++)
                            {printf NF "%12.6f %12.6f %12.6f %12.6f %12.6f\r\n",
                            $x[$i], $y[$i], $z[$i],$radius[$i],$factor[$j]
                            };print NF "\n--link1--";$n_blanks++}};
       if(!($line=~ /\S/)) {print NF $line;$n_blanks++}
}
#
close(SCR, $scrfile)};
#
open(SCR, $scrfile);
@array2=split(/\./,$factor[$j]);
$chk= $array1[0].'_xp-'.$array2[0].$array2[1].'.chk';
$n_blanks=0;
while ($line = <SCR>) {
  if($line=~ /\S/) 
    { if(!($n_blanks >= 4)) {if($line=~/chk/){chop($line); print NF $line,$chk,"\n"}; if(!($line=~/chk/)) {print NF $line}};
      if($n_blanks == 4){ for($i=1; $i<=$n_sfe;$i++)
                            {printf NF "%12.6f %12.6f %12.6f %12.6f %12.6f\r\n",
                             $x[$i], $y[$i], $z[$i],$radius[$i],$factor[$j]
                            };$n_blanks++}};
       if(!($line=~ /\S/)) {print NF $line;$n_blanks++}
};
close(SCR, $scrfile);
#
# Running Gaussian 09 
#
my $cmd = 'g16 '.$newfile;
  system($cmd);
print "\nend g16";
#
# Computing parameter s(f) of Eq. 7 in SI in Ref. (2) for the structure $k
#
$flog = $array1[0].'_xp-Vc'.'.log';
open(FL, $flog);
open(FNI3, ">>properties.dat");
print FNI3 "geom ",$k,"\n";
#
$n_fact = 0;
    while($line = <FL> ) {
#          if($line =~ /ISph/){
#             $line = <FL>; 
#             if(!($line =~ /--/)){$n_fact++;chop($line);
#              @simpson = split(/\s+/,$line);$f[$n_fact] = $simpson[3]}};
          if($line =~ / GePol: Cavity volume/){
              $n_fact++;$f[$n_fact]=$factor[$n_fact-1];
              chop($line);
              @simpson = split(/\s+/,$line);
              $volume[$n_fact] = $simpson[5];print '$n_fact= ',$n_fact,'$volume[$n_fact]=',$volume[$n_fact];
              $s[$n_fact] = ($volume[$n_fact]/$volume[1])**(1/3); 
              print FNI3 $n_fact," ",$f[$n_fact]," ",$volume[$n_fact]," ",$s[$n_fact],"\n"};
                         }
# end loop over structures  
};
close(FNI3, "properties.dat");
#------------------------------------------------------------------------------
#
#  Step 2: Computing the permittivity \eps(f) and density RhoS(f) of the 
#          external medium (see Eq.s (9) and (10) of SI in Ref. (2))
#
#------------------------------------------------------------------------------
open(FNI3, "properties.dat");
open(FNI4, ">eps_rhos.dat");
#
$k = 0;
    while($line = <FNI3> ) {
           if($line =~ /geom/){$k++;$n_fact = 0}; 

          if(!($line =~ /geom/)){$n_fact++;if($k == 1){$s[$n_fact] = 0};
              chop($line);
              @simpson = split(/\s+/,$line);
              $s[$n_fact] =  $s[$n_fact]+$simpson[3]/$n_struct; 
              if ($k == $n_struct){$eps[$n_fact]= 1+($eps_0-1)/$s[$n_fact]**3;
              $rhos[$n_fact]= $rho_0/$s[$n_fact]**(3+$eta);
              $vmol[$n_fact]= ($mw/$rho_0)*$s[$n_fact]**3;
                                   print FNI4 $n_fact," ",$f[$n_fact]," ",$s[$n_fact]," ",$eps[$n_fact]," ",$rhos[$n_fact]," ", $vmol[$n_fact],"\n"}
                               };
	                      };
}
#
# Subroutine ElectronicEnergy
#
sub ElectronicEnergy {
#Start Ger
for($k=1; $k<=$n_struct;$k++){
#
# Running Gaussian   for the electronic energy G_{er}(f). 
#
$oldfile = $ARGV[$k];
@array1=split(/\./,$oldfile);
if ($k >= $kstart ){
$ger_fn= $array1[0].'_xp-Ger'.'.com';
open(FNI1, ">$ger_fn"); 
$newfile = $array1[0].'_xp-Vc'.'.com'; 
open(FI,$newfile);
$i_fact = 0;
    while($line = <FI> ) {
          if(!($line =~ /# scrf=/)){
            if(!($line =~ /eps=/)){print FNI1 $line};
                                   }
          if($line =~ /# scrf=/) {print FNI1 "# scrf=(iefpcm,solvent=",$solvent,",read) nosym 6d 10f \n"};
          if($line =~ /eps=/){$i_fact++;
                               print FNI1 "eps=",$eps[$i_fact]," RhoS=",$rhos[$i_fact],"\n";
                              }
                          }
close(FI);
my $cmd = 'g16 '.$ger_fn;
  system($cmd);
print "\nend g16";
}
#
# Saving G_{er}(f) vs V_c(f) . 
#
$flog = $array1[0].'_xp-Ger'.'.log'; 
$fout = $array1[0].'_VcGer'.'.dat';
open(FL, $flog);
open(OUT, ">$fout");
open(OUT2, ">VGer.dat");
print OUT "#Vc_Ger(",$solvent," eta=",$eta,")\n";
printf OUT "%1s %11s %14s\r\n", '#','Vc','G_{er}';
printf OUT "%1s %11s %14s\r\n", '#','Ang^3','a.u.';
$n_fact=0;
    while($line = <FL> ) {
              if($line =~ / GePol: Cavity volume/){$n_fact++;
                chop($line);
                @simpson = split(/\s+/,$line);
                $volume[$n_fact] = $simpson[5]};
              if($line =~ /SCF Done/){
                chop($line);
                @simpson = split(/\s+/,$line);$e[$n_fact] = $simpson[5];               
                $ger_mat[$k][$n_fact] = $e[$n_fact]; 
                 printf OUT "%12.6f %14.8f\r\n", $volume[$n_fact], $e[$n_fact];
                 print OUT2  $volume[$n_fact]," ", $e[$n_fact],"\n"};
	                     };
close(OUT, $fout);
#------------------------------------------------------------------------------
#
# Step 4: Computing the pressure p=-dG_er/dVc
#
#------------------------------------------------------------------------------
 Pressure()}
}
#End Ger----------------------------------------------------------------------------------------
#
# Subroutine Pressure
#
sub Pressure {
#
#  ========== Mathematica ==========
#
if($mrgn =~ /math/){
system("math -script murnaghan.m > t.out");
open(IN,'t.out'); 
open(OUT,'>murnaghan.dat');
print OUT "Murnaghan fitting (",$solvent," eta=",$eta,")\n";
while($line = <IN>) {
      if($line =~ /BestFitParameters/){
        @array1= split(/\s+/,$line);
        @array2= split(/\,/,$array1[4]);
        @array3= split(/\,/,$array1[7]);
        @array4= split(/\}/,$array1[10]);
        $a = $array2[0]; print OUT 'a=',$a,"\n";        
        $b = $array3[0]; print OUT 'b=',$b,"\n"; 
        $c = $array4[0]; print OUT 'c=',$c,"\n";            
                                       }
                     }
}
#-------------------------------------------
#  ========== Phyton ==========
#
if($mrgn =~ /py/){
system("python murnaghan.py > t.out");
open(IN,'t.out'); 
open(IN2,'VGer.dat');
$line = <IN2>; @array0= split(/\s+/,$line);
$V0=$array0[0];
open(OUT,'>murnaghan.dat');
print OUT "Murnaghan fitting (",$solvent," eta=",$eta,")\n";
while($line = <IN>) {
      if($line =~ /a:/){
        @array1= split(/\s+/,$line);
        $a = $array1[2]/$V0; print OUT 'a=',$a,"\n"; 
        $line = <IN>;  @array1= split(/\s+/,$line);
        $b = $array1[2]; print OUT 'b=',$b,"\n"; 
        $line = <IN>;  @array1= split(/\s+/,$line);       
        $c = $array1[2]/$V0; print OUT 'c=',$c,"\n";              
                                       }
                     }
}
##
@array1=split(/\./,$oldfile);
$out2 = $array1[0].'_Vp'.'.dat';
open(FL, "VGer.dat");
open(OUT, ">>Vp.dat");
#open(OUT, ">Vp.dat");
print OUT "Structure ",$k,"\n";
open(OUT2, ">$out2");
print OUT2 "p-Vc (",$solvent," eta=",$eta,")\n";
printf OUT2 "%8.4s %12.6s %10.3s\r\n", 'f','Vc','p';
printf OUT2 "%8.4s %12.6s %10.3s\r\n", '','Ang^3','GPa';
$n_fact = 0;
#
while($line = <FL> ) {$n_fact++;
                chop($line);
                @simpson = split(/\s+/,$line);
                $volume[$n_fact] = $simpson[0];
                $p[$n_fact]=($a*(($volume[1]/$volume[$n_fact])**($b+1)-1)+$c)*$Eh_o_Ang3_to_GPa;
                printf OUT2 "%8.4f %12.6f %10.3f\r\n", $f[$n_fact],$volume[$n_fact],$p[$n_fact];
                print OUT  $f[$n_fact]," ",$volume[$n_fact]," ",$p[$n_fact],"\n"}
}
#End Pressure-----------------------------------------------------------------------------------
#
# Subroutine MeanPressure
#
sub MeanPressure {
close(OUT, "Vp.dat");
open(OUT, "Vp.dat");
open(FNI3, "Vp.dat");
open(FNI4, ">PES_Vp.dat");
print FNI4 "p values averaged over ",$n_struct," structures (",$solvent," eta=",$eta,")\n";
printf FNI4 "%8s %10.3s\r\n", 'f','p'; 
printf FNI4 "%8.4s %10.3s\r\n", '','GPa';
#
    while($line = <FNI3> ) {
          if($line =~ /Stru/){$n_fact = 0;chop($line);@simpson = split(/\s+/,$line);
                $k = $simpson[1]}; 
          if(!($line =~ /Stru/)){$n_fact++;if($k == 1){$p[$n_fact] = 0};
              chop($line);
              @simpson = split(/\s+/,$line);
              $f[$n_fact] = $simpson[0];
              $p[$n_fact] =  $p[$n_fact]+$simpson[2]/$n_struct};
               if (($k == $n_struct)&&($n_fact > 0)){printf FNI4 "%8.4f %10.3f\r\n", $f[$n_fact],$p[$n_fact]};
}
}
#
# Subroutine CavitationEnergy
#
sub CavitationEnergy{
$OUT= 'xp_pes.log';   
$OUT2= 'xp_Gcav.log';   
open(OUT,">>$OUT"); 
print OUT "\nTable I: 
Cavity scaling factor,f, pressure, p (GPa), Electronic (G_er), cavitation (G_cav), and total (G_tot) 
energies";
open(OUT2,">$OUT2");print OUT2 "XP_Gcav results: for solvent=",$solvent," eta=",$eta," cavity=",$cavity," frad=",$frad," tsare=",$tsare,"\n";
close(OUT,$OUT);
close(OUT2,$OUT2);
for($k=1; $k<=$n_struct;$k++){  #Start loop for structures
#
# Creating the input files 
# 
$oldfile = $ARGV[$k];
@array1=split(/\./,$oldfile);
$scrfile = $array1[0].'_scrfile';
#
@array1=split(/\./,$oldfile);
$newfile = $array1[0].'_xp-Gcav'.'.com'; 
open(NF, ">$newfile");
#
#for($j=0; $j<=$#ARGV-(7+$n_struct);$j++){ #Start loop for factors f
for($j=0; $j<=($nfm1-1);$j++){ #Start loop for factors f
open(SCR, $scrfile);
$chk= 'scr.chk';
$n_blanks=0;
while ($line = <SCR>) {if(!($line=~ /\S/)){$n_blanks++;};
                       if(!($line=~ /chk/)&&!($line=~ /xp-pcm/)&&!($line=~ /qrep /)&&!($line=~ /NVE/)&&!($line=~ /eps/)) {print NF $line};
                       if($line=~ /chk/) {print NF '%chk=',$chk,"\n"};
                       if($line=~ /xp-pcm/) {print NF "xp-pcm with Vmol=",$vmol[$j+1]," cm^3/mol ","\n"};
                       if($line=~ /qrep /) {print NF "norep  pcmdoc geomview  nodis cav g03defaults ", $noaddsph," tsare=",$tsare," \n"};
                       if($line=~ /NVE/) {print NF "Vmol=",$vmol[$j+1], " Rsolv=", $rsolv," \n"};

                       if($n_blanks == 5){print NF "--link1--\n"};
                       } #end while
                                        };#End loop for factors f
# Last  f
open(SCR, $scrfile);
$n_blanks=0;
while ($line = <SCR>) {if(!($line=~ /\S/)){$n_blanks++};
                       if(!($line=~ /chk/)&&!($line=~ /xp-pcm/)&&!($line=~ /qrep /)&&!($line=~ /NVE/)&&!($line=~ /eps/)) {print NF $line};
                       if($line=~ /chk/) {print NF '%chk=',$chk,"\n"};
                       if($line=~ /xp-pcm/) {print NF "xp-pcm with Vmol=",$vmol[$j+1]," cm^3/mol ","\n"};
                       if($line=~ /qrep /) {print NF "norep  pcmdoc geomview  nodis cav g03defaults ", $noaddsph," tsare=",$tsare," \n"};
                       if($line=~ /NVE/) {print NF "Vmol=",$vmol[$j+1], " Rsolv=", $rsolv," \n"};
                       }; #end while
#
# Running Guassian 16 for the cavitation energy
#
my $cmd = 'g16 '.$newfile;
  system($cmd);
print "\nend g16";
#
$flog= $array1[0].'_xp-Gcav'.'.log';
open(FL, $flog); 
open(OUT,">>$OUT");
open(OUT2,">>$OUT2");
#print OUT "\n\n", "Structure ", $k,"\n";
print OUT "\n\n", "Structure ", $array1[0],"\n";
#printf OUT " %5s %8s %12s %12s %12s %12s %12s %12s\r\n", 'f','p[GPa]','G_er[au] ','G_cav[au]','G_tot[au]','G_er[kcal/mol] ','G_cav[kcal/mol]','G_tot[kcal/mol]';
#printf OUT " %5s %6s %10s %10s %10s %10s %10s %10s\r\n", 'f','p[GPa]','G_er[au] ','G_cav[au]','G_tot[au]','G_er[kcal/mol] ','G_cav[kcal/mol]','G_tot[kcal/mol]';
printf OUT " %5s %6s %10s %8s %10s %10s %10s %10s\r\n", 'f  ',' p   ','G_er   ','G_cav  ','G_tot  ','G_er   ','G_cav  ','G_tot   ';
printf OUT " %5s %5s %10s %8s %10s %10s %10s %10s\r\n", ' ',' [GPa]', '[Eh]   ','[Eh]   ','[Eh]   ','[kcal/mol] ','[kcal/mol] ','[kcal/mol]     ';
#
#print OUT2 "\n\n", "Structure ", $k,"\n";
print OUT2 "\n\n", "Structure ", $array1[0],"\n";
printf OUT2 " %5s %8s %10s %12s %12s %12s\r\n", 'f','p[GPa]','V_cav[Ang^3]','G_cav-pV[kcal/mol]','pV[kcal/mol]','G_cav[kcal/mol]';
$n_fact = 0;  
 while($line = <FL> ) {
              if($line =~ / GePol: Cavity volume/){$n_fact++;
                chop($line);
                @simpson = split(/\s+/,$line);
                $volume[$n_fact] = $simpson[5]};
              if($line =~ /PCM non-electrostatic energy/){
                chop($line);
                @simpson = split(/\s+/,$line);$ecav[$n_fact] = $simpson[5];
                $pV[$n_fact]=$p[$n_fact]*$volume[$n_fact]*(1.0/$cal2J)*$avogadro;
                $gcav[$n_fact]=$ecav[$n_fact]*$Eh2kcal_mol+$pV[$n_fact];
                $cav_mat[$k][$n_fact]= $gcav[$n_fact];
                $gtot[$n_fact]=$e[$n_fact]+$gcav[$n_fact]/$Eh2kcal_mol;
                $gtot_mat[$k][$n_fact] =$ger_mat[$k][$n_fact]+ $cav_mat[$k][$n_fact]/$Eh2kcal_mol;
                printf OUT2 "%5.4f %8.3f %10.3f %12.2f %12.2f %12.2f \r\n", $f[$n_fact],$p[$n_fact],$volume[$n_fact],$ecav[$n_fact]*$Eh2kcal_mol, $pV[$n_fact],  $cav_mat[$k][$n_fact];
#                printf OUT "%5.4f %8.3f %12.6f %12.6f %12.6f %12.2f %12.2f %12.2f\r\n", $f[$n_fact],$p[$n_fact],$ger_mat[$k][$n_fact],$cav_mat[$k][$n_fact]/$Eh2kcal_mol,$gtot_mat[$k][$n_fact], $ger_mat[$k][$n_fact]*$Eh2kcal_mol,$cav_mat[$k][$n_fact], $gtot_mat[$k][$n_fact]*$Eh2kcal_mol      
#                printf OUT "%5.3f %6.1f %10.4f %10.4f %10.4f %11.1f %11.1f %11.1f\r\n", $f[$n_fact],$p[$n_fact],$ger_mat[$k][$n_fact],$cav_mat[$k][$n_fact]/$Eh2kcal_mol,$gtot_mat[$k][$n_fact], $ger_mat[$k][$n_fact]*$Eh2kcal_mol,$cav_mat[$k][$n_fact], $gtot_mat[$k][$n_fact]*$Eh2kcal_mol  
printf OUT "%5.3f %5.1f %10.4f %8.4f %10.4f %10.1f %10.1f %10.1f\r\n", $f[$n_fact],$p[$n_fact],$ger_mat[$k][$n_fact],$cav_mat[$k][$n_fact]/$Eh2kcal_mol,$gtot_mat[$k][$n_fact], $ger_mat[$k][$n_fact]*$Eh2kcal_mol,$cav_mat[$k][$n_fact], $gtot_mat[$k][$n_fact]*$Eh2kcal_mol                                
                          };
	                  } # end while
                                  } #End loop for structures
}

