#!/usr/bin/perl -w
# use Parallel::ForkManager;
# my $forks = 3;
# # shift or die "Usage: $2 N\n";
# my $pm = Parallel::ForkManager->new($forks);



$GMXPATH="/usr/local/Cellar/gromacs/2021.2/bin";
$CONTFILE="/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/157_1/con.txt";
$CUTOFF=1.2;
$SKIPFRAMES=1;
open(CONaa,"$CONTFILE")  or die "no CA contacts file";

$CONNUMaa=0;
while(<CONaa>){
	$ConI=$_;
	chomp($ConI);
	$CONNUMaa ++;

	@TMP=split(" ",$ConI);
	$Iaa[$CONNUMaa]=$TMP[0];
	$Jaa[$CONNUMaa]=$TMP[1];
	# this gets the native distances from the 6-12 parameters and converts to Ang.
	$Raa[$CONNUMaa]=(6/5*$TMP[3]/$TMP[2])**(1/2.0)*10.0;
	
}

close(CONaa);


sub Q{
$TPR=$_[0];
chomp($TPR);
# go through all of the xtc files
$XTCfile=$_[1];
chomp($XTCfile);
print "\n\n\n\n\n";
print "This script will calculate the number of native contacts for frames n*$SKIPFRAMES (n starts at 0) for trajectory file $XTCfile.\n";
print "IMPORTANT: A contact is defined as any native pair listed in $CONTFILE that is withing $CUTOFF times the distance in $CONTFILE\n";
print "\n\n\n\n\n";
# you can change the flags sent to gromacs if you want to reduce the Q sampling
# i.e. every other frame is analyzed.  The trajectory file is converted to a pdb
# and then this script analyzes the pdb.
`echo 0 | $GMXPATH/gmx trjconv  -skip $SKIPFRAMES -s $TPR -o tmp.pdb -f  $XTCfile`;
open(CAQ,">$XTCfile.CA.Q") or die "couldn't open Q file";
open(CAQI,">$XTCfile.CA.Qi") or die "couldn't open Qi file";
open(PDB,"tmp.pdb") or die "pdb file missing somehow";
#GRAB THE ATOM residue numbers
while(<PDB>){
     	$LINE=$_;
       	chomp($LINE);
        if(substr($LINE,0,5) eq "MODEL"){
        	# start the residue storing
               	$ATOMNUM=0;
        }
        if(substr($LINE,0,4) eq "ATOM"){
           	$ATOMNUM=$ATOMNUM+1;
        }
        if(substr($LINE,0,3) eq "TER"){
		last;
	}
}
close PDB;
open(PDB,"tmp.pdb") or die "pdb file missing somehow";
###$SAMPLES=0; #I think this is obsolete
	#GRAB THE ATOM COORDS
while(<PDB>){
	$LINE=$_;
	chomp($LINE);
	if(substr($LINE,0,5) eq "MODEL"){
		$ATOMNUM=0;
	}
        if(substr($LINE,0,4) eq "ATOM"){
	        # store positions, index and residue number
		$ATOMNUM=$ATOMNUM+1;
	        $X[$ATOMNUM]=substr($LINE,30,8);
                $Y[$ATOMNUM]=substr($LINE,38,8);
               	$Z[$ATOMNUM]=substr($LINE,46,8);

        }
        if(substr($LINE,0,3) eq "TER"){			
	# do contact calculations
		# determine the atom-atom contacts
			# for($k = 0;$k<=6930; $k +=110){ ###-------> uncomment this block for crowded cases and comment the following block 
               	# $Qaa=0;
         		# for($i = 1;$i <= $CONNUMaa; $i +=1){
			    # 	$R=sqrt( ($X[$Iaa[$i]+$k]-$X[$Jaa[$i]+$k])**2.0 + 
			    #        	 ($Y[$Iaa[$i]+$k]-$Y[$Jaa[$i]+$k])**2.0 + 
  			 	#         ($Z[$Iaa[$i]+$k]-$Z[$Jaa[$i]+$k])**2.0);
 			    # 	if($R < $CUTOFF*$Raa[$i]){
				#     	$Qaa ++;
				#     	print CAQI "$i\n";
			    # 	}  
		    	# }
	            # print CAQ "$Qaa\n";
            # }	###--------> uncomment this for crowded cases
			 	$Qaa=0;

		for($i = 1;$i <= $CONNUMaa; $i +=1){
			$R=sqrt( ($X[$Iaa[$i]]-$X[$Jaa[$i]])**2.0 + 
			       	 ($Y[$Iaa[$i]]-$Y[$Jaa[$i]])**2.0 + 
  			 	 ($Z[$Iaa[$i]]-$Z[$Jaa[$i]])**2.0);
			 print "$R\n";		
 			if($R < $CUTOFF*$Raa[$i]){
				$Qaa ++;
				 print CAQI "$i\n";
			}
		}
		# sleep(60);
	         print CAQ "$Qaa\n";
         } 
}
close(CAQ);
close(CAQI);

`rm tmp.pdb`;
}



# @a = ('130_1', '140_1', '150_1', '155_1', '158_1', '159_1', '160_1', '162_1', '163_1', '164_1', '165_1', '170_1', '180_1', '190_1', '200_1', '210_1', '220_1'); 
# @a = ( '150_1', '155_1', '156_1',  '158_1', '159_1', '160_1', '170_1'); 
#  @a = ( '120_3', '140_3', '150_3', '155_3', '156_3', '157_3', '158_3', '159_3', '160_3', '180_3' ); 
@a = ('170_3');
for $i (@a){
    # my $pid = $pm->start and next;

    $mc_c="/Users/atrayee/Desktop/SOD_project/Data/WT_Charged/phic_0/$i/";
     $mc_u="/Users/atrayee/Desktop/SOD_project/Data/WT_Uncharged/phic_0/$i/";
     $mc_m="/Users/atrayee/Desktop/SOD_project/Data/G41D_Mutant/phic_0/$i/";
    chdir($mc_c) or die "$!";
    Q('topol.tpr','traj_comp_pbcmolcenter_WT_c.xtc');
    chdir($mc_u) or die "$!";
    Q('topol.tpr','traj_comp_pbcmolcenter_WT_u.xtc');
    chdir($mc_m) or die "$!";
    Q('topol.tpr','traj_comp_pbcmolcenter_WT_m.xtc');

    # $pm->finish;

}
# $pm->wait_all_children;