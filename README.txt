*This is README file for Residue Contact Frequency described in the following paper

Binding interface prediction by combining protein–protein docking results
H Hwang, T Vreven, Z Weng - Proteins: Structure, Function, and Bioinformatics, 2014

RCF.pl calculates Residue Contact Frequency (RCF) of each residue for a pair of protein structures using ZDOCK 3.0 output file.
It requires A ZDOCK 3.0 output file and ZDOCK input protein structures that were used for ZDOCK run.
The ZDOCK input protein structure files and ZDOCK output file must be at the current directory with RCF.pl.

*Construct executables

The source codes are written in c++ and perl.
The c++ source codes are under the source directory and executables can be created by following command.

>make

*To run RCF.pl

Put the ZDOCK 3.0 output file and its input protein structures in the current directory and run RCF.pl as the following 
command. (The ZDOCK 3.0 CG output file <1A2K.zd3.0.cg.out> and its two ZDOCK input protein structures <1A2K_r_u.pdb and 
1A2K_l_u.pdb> are included as an example)

>perl RCF.pl

Then it will show how to use the program as the followings.

>perl RCF.pl zdock_out_file rotational_sampling_mode

For rotational_sampling_mode type FG for 54K fine sampling mode
Or type CG for 3600 coarse sampling mode

For an example,

>perl RCF.pl 1A2K.zd3.0.cg.out CG

It will create two files with .RCF extension which contain RCF for each residue as the followings.

1A2K_r_u.pdb.CG.RCF
1A2K_l_u.pdb.CG.RCF