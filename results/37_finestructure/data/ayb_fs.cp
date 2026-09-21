section:fsproject
# This file contains the information about the project. It has the following sections:
############## stage0 (data processing to create the required input files):
############## stage1 (chromopainter parameter estimation via EM)
############## stage2 (chromopainter full painting)
############## stage3 (finestructure mcmc)
############## stage4 (finestructure tree)
#
# You may modify this file by hand, but it is possible to break it if you are careless.  It is advisable instead to set parameters via the "fs" interface; either from the command line or by placing the set of parameters to set in a settings file and using -import.
# Each section is separated by "section:<sectionname>".  Within each section the parameters are defined in the format:
# paramname:value
# A '#' symbol indicates the end of the readable line, everything afterwards is treated as a comment.
# Everything from the colon to the end of the line (or first #) is treated as content. Parameters that are empty are assigned their default values as described in "fs -h".

############## Command history:
section:history
CMDLINE: -idfile ayb.ids -phasefiles Ss01.phase,Ss02.phase,Ss03.phase,Ss04.phase,Ss05.phase,Ss06.phase,Ss07.phase,Ss08.phase,Ss09.phase,Ss10.phase,Ss11.phase -recombfiles Ss01.recombfile,Ss02.recombfile,Ss03.recombfile,Ss04.recombfile,Ss05.recombfile,Ss06.recombfile,Ss07.recombfile,Ss08.recombfile,Ss09.recombfile,Ss10.recombfile,Ss11.recombfile -hpc 0 -ploidy 2 -go
CMD:-idfile ayb.ids
CMD:-phasefiles Ss01.phase Ss02.phase Ss03.phase Ss04.phase Ss05.phase Ss06.phase Ss07.phase Ss08.phase Ss09.phase Ss10.phase Ss11.phase
CMD:-recombfiles Ss01.recombfile Ss02.recombfile Ss03.recombfile Ss04.recombfile Ss05.recombfile Ss06.recombfile Ss07.recombfile Ss08.recombfile Ss09.recombfile Ss10.recombfile Ss11.recombfile
CMD:-hpc 0
CMD:-ploidy 2
CMD:-countdata
CMD:-makes1
CMD:-dos1
CMD:-combines1
CMD:-makes2
CMD:-dos2
CMD:-combines2
CMD:-makes3
CMD:-dos3
CMD:-makes3
CMD:-dos3
CMD:-combines3
CMD:-combines3
CMD:-makes4
CMD:-dos4
CMD:-combines4
CMD:-go
CMDLINE: -go
CMD:-go
CMDLINE: -go
CMD:-go
CMDLINE: -go
CMD:-go
CMDLINE: -go
CMD:-go
CMDLINE: -go
CMD:-go
CMDLINE: -go
CMD:-go
CMDLINE: -go
CMD:-go
CMDLINE: -go
CMD:-go
CMDLINE: -go
CMD:-go
CMDLINE: -go
CMD:-go
CMDLINE: -go
CMD:-go
CMDLINE: -go 
CMD:-go

section:parameters

###################
stage:0
## Data preparation and conversion. Not yet implemented

validatedoutput:1,1,1,1,1  # Derived. Whether we have validated output from each stage of the analysis (0-4)

###################
stage:1
## Universal Stage1-4 properties of the analysis

exec:fs  # Finestructure command line. Set this to be able to use a specific version of this software. (default: fs)
hpc:0  # THIS IS IMPORTANT FOR BIG DATASETS! Set hpc mode. 0: Commands are run 'inline' (see 'numthreads' to control how many CPU's to use). 1: Stop computation for an external batch process, creating a file containing commands to generate the results of each stage. 2: Run commands inline, but create the commands for reference. (default: 0.)
numthreads:0  # Maximum parallel threads in 'hpc=0' mode. Default: 0, meaning all available CPUs.
ploidy:2  # Haplotypes per individual. =1 if haploid, 2 if diploid. (default: 2)
linkagemode:linked  # unlinked/linked. Whether we use the linked model. default: unlinked / linked if recombination files provided.
indsperproc:0  # Desired number of individuals per process (default: 0, meaning autocalculate: use 1 In HPC mode, ceiling(N/numthreads) otherwise. Try to choose it such that you get a sensible number of commands compared to the number of cores you have available.
outputlogfiles:1  # 1=Commands are written to file with redirection to log files. 0: no redirection. (default:1)
allowdep:1  # Whether dependency resolution is allowed. 0=no, 1=yes. Main use is for pipelining. (default:1).
## ChromoPainter Stage1-2 generic properties

s12inputtype:phase  # What type of data input (currently only "phase" supported)
idfile:ayb.ids  # IDfile location, containing the labels of each individual. REQUIRED, no default (unless -createids is used).
s12args:  # arguments to be passed to Chromopainter (default: empty)
## Quantities observed from data

ninds:95  # Derived. number of individuals observed in the idfile
nindsUsed:95  # Derived. number of individuals retained for processing from the idfile
nsnps:1454  # Derived. number of SNPs in total, over all files
## ChromoPainter Stage1 (EM) properties

s1args:-in -iM --emfilesonly  # Arguments passed to stage1 (default:-in -iM --emfilesonly)
s1emits:10  # Number of EM iterations (chromopainter -i <n>, default: 10)
s1minsnps:10000  # Minimum number of SNPs for EM estimation (for chromopainter -e, default: 10000)
s1snpfrac:0.1  # fraction of genome to use for EM estimation. (default: 0.1)
s1indfrac:1  # fraction of individuals to use for EM estimation. (default: 1.0)
s1outputroot:ayb_fs_tmp_EM_linked  # output file for stage 1 (default is autoconstructed from filename)
## ChromoPainter Stage2 properties inferred from Stage1

Neinf:2039.17  # Derived. Inferred `Effective population size Ne' (chromopainter -n).
muinf:0.00231611  # Derived. Inferred Mutation rate mu (chromopainter -M)

###################
stage:2
## ChromoPainter Stage2 (main run) properties

s2chunksperregion:-1  # number of chunks in a "region" (-ve: use default of 100 for linked, nsnps/100 for unlinked)
s2samples:0  # number of samples of the painting to obtain per recipient haplotype, for examining the details of the painting. (Populates <root>.samples.out; default 0. Warning: these file can get large)
s2args:  # Additional arguments for stage 2 (default: none, "")
s2outputroot:ayb_fs_tmp_mainrun.linked  # Output file name for stage 2 (default: autoconstructed).
s2combineargs:  # Additional arguments for stage 2 combine (fs combine; default: none, "")
## FineSTRUCTURE Stage3 properties inferred from Stage2

cval:0.0275982  # Derived. 'c' as inferred using chromopainter. This is only used for sanity checking. See s34 args for setting it manually.
cproot:ayb_fs_linked  # The name of the final chromopainter output. (Default: <filename>, the project file name)
cpchunkcounts:ayb_fs_linked.chunkcounts.out  # the finestructure input file, derived name of the chunkcounts file from cproot.

###################
stage:3
## FineSTRUCTURE Stage3-4 generic properties

fsroot:ayb_fs_linked  # The name of the finestructure output (Default: <filename>, the project file name).
s34args:  # Additional arguments to both finestructure mcmc and tree steps. Add "-c <val>" to manually override 'c'.
## FineSTRUCTURE Stage3 MCMC inference

s3iters:-1  # Number of TOTAL iterations to use for MCMC. By default we assign half to burnin and half to sampling. (default: 100000)
s3iterssample:100000  # Number of iterations to use for MCMC (default: -ve, meaning derive from s3iters)
s3itersburnin:0  # Number of iterations to use for MCMC burnin (default: -ve, meaning derive from s3iters)
numskip:200  # Number of mcmc iterations per retained sample; (default: -ve, meaning derive from maxretained)
maxretained:500  # Maximum number of samples to retain when numskip -ve. (default: 500)
nummcmcruns:2  # Number of *independent* mcmc runs. (default: 2)
fsmcmcoutput:ayb_fs/stage3/ayb_fs_linked_mcmc  # Filename to use for mcmc output (default: autogenerated)
mcmcGR:1.00125,1.03938,0.999371,1.01237,0.999847  # Derived. Gelman-Rubin diagnostics obtained from combining MCMC runs, for log-posterior, K,log-beta,delta,f respectively
threshGR:1.3  # Threshold for the Gelman-Rubin statistic to allow moving on to the tree building stage. We always move on if thresGR<0. (Default: 1.3)

###################
stage:4
## FineSTRUCTURE Stage4 tree inference

s4args:  # Extra arguments to the tree building step. (default: none, "")
s4iters:100000  # Number of maximization steps when finding the best state from which the tree is built. (default: 100000)
fstreeoutput:ayb_fs/stage4/ayb_fs_linked_tree  # Filename to use for finestructure tree output. (default: autogenerated)

###################
stage:1
## Vector quantities placed at the end for readability

phasefiles:Ss01.phase,Ss02.phase,Ss03.phase,Ss04.phase,Ss05.phase,Ss06.phase,Ss07.phase,Ss08.phase,Ss09.phase,Ss10.phase,Ss11.phase  # Comma or space separated list of all 'phase' files containing the (phased) SNP details for each haplotype. Required. Must be sorted alphanumerically to ensure chromosomes are correctly ordered. So don't use *.phase, use file{1..22}.phase. Override this with upper case -PHASEFILES.
recombfiles:Ss01.recombfile,Ss02.recombfile,Ss03.recombfile,Ss04.recombfile,Ss05.recombfile,Ss06.recombfile,Ss07.recombfile,Ss08.recombfile,Ss09.recombfile,Ss10.recombfile,Ss11.recombfile  # Comma or space separated list of all recombination map files containing the recombination distance between SNPs. If provided, a linked analysis is performed. Otherwise an 'unlinked' analysis is performed. Note that linkage is very important for dense markers!
nsnpsvec:126,119,222,208,106,111,77,157,120,134,74  # Derived. Comma separated list of the number of SNPs in each phase file.
s1outputrootvec:ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind1-8,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind9-16,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind17-24,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind25-32,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind33-40,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind41-48,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind49-56,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind57-64,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind65-72,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind73-80,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind81-88,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file1_ind89-95,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind1-8,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind9-16,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind17-24,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind25-32,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind33-40,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind41-48,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind49-56,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind57-64,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind65-72,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind73-80,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind81-88,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file2_ind89-95,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind1-8,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind9-16,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind17-24,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind25-32,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind33-40,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind41-48,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind49-56,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind57-64,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind65-72,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind73-80,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind81-88,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file3_ind89-95,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind1-8,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind9-16,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind17-24,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind25-32,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind33-40,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind41-48,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind49-56,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind57-64,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind65-72,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind73-80,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind81-88,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file4_ind89-95,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind1-8,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind9-16,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind17-24,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind25-32,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind33-40,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind41-48,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind49-56,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind57-64,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind65-72,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind73-80,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind81-88,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file5_ind89-95,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind1-8,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind9-16,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind17-24,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind25-32,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind33-40,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind41-48,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind49-56,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind57-64,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind65-72,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind73-80,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind81-88,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file6_ind89-95,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind1-8,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind9-16,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind17-24,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind25-32,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind33-40,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind41-48,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind49-56,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind57-64,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind65-72,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind73-80,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind81-88,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file7_ind89-95,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind1-8,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind9-16,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind17-24,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind25-32,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind33-40,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind41-48,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind49-56,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind57-64,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind65-72,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind73-80,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind81-88,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file8_ind89-95,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind1-8,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind9-16,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind17-24,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind25-32,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind33-40,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind41-48,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind49-56,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind57-64,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind65-72,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind73-80,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind81-88,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file9_ind89-95,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind1-8,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind9-16,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind17-24,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind25-32,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind33-40,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind41-48,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind49-56,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind57-64,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind65-72,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind73-80,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind81-88,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file10_ind89-95,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind1-8,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind9-16,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind17-24,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind25-32,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind33-40,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind41-48,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind49-56,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind57-64,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind65-72,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind73-80,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind81-88,ayb_fs/stage1/ayb_fs_tmp_EM_linked_file11_ind89-95  # Derived. Comma separated list of the stage 1 output files names.

###################
stage:2
## 

s2outputrootvec:ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind1-8,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind9-16,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind17-24,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind25-32,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind33-40,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind41-48,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind49-56,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind57-64,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind65-72,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind73-80,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind81-88,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file1_ind89-95,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind1-8,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind9-16,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind17-24,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind25-32,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind33-40,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind41-48,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind49-56,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind57-64,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind65-72,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind73-80,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind81-88,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file2_ind89-95,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind1-8,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind9-16,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind17-24,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind25-32,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind33-40,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind41-48,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind49-56,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind57-64,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind65-72,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind73-80,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind81-88,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file3_ind89-95,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind1-8,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind9-16,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind17-24,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind25-32,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind33-40,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind41-48,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind49-56,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind57-64,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind65-72,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind73-80,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind81-88,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file4_ind89-95,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind1-8,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind9-16,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind17-24,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind25-32,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind33-40,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind41-48,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind49-56,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind57-64,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind65-72,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind73-80,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind81-88,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file5_ind89-95,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind1-8,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind9-16,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind17-24,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind25-32,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind33-40,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind41-48,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind49-56,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind57-64,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind65-72,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind73-80,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind81-88,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file6_ind89-95,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind1-8,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind9-16,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind17-24,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind25-32,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind33-40,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind41-48,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind49-56,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind57-64,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind65-72,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind73-80,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind81-88,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file7_ind89-95,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind1-8,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind9-16,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind17-24,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind25-32,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind33-40,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind41-48,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind49-56,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind57-64,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind65-72,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind73-80,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind81-88,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file8_ind89-95,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind1-8,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind9-16,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind17-24,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind25-32,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind33-40,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind41-48,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind49-56,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind57-64,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind65-72,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind73-80,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind81-88,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file9_ind89-95,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind1-8,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind9-16,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind17-24,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind25-32,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind33-40,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind41-48,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind49-56,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind57-64,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind65-72,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind73-80,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind81-88,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file10_ind89-95,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind1-8,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind9-16,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind17-24,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind25-32,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind33-40,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind41-48,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind49-56,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind57-64,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind65-72,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind73-80,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind81-88,ayb_fs/stage2/ayb_fs_tmp_mainrun.linked_file11_ind89-95  # Derived. Comma separated list of the stage 2 output files names.

###################
stage:3
## 

fsmcmcoutputvec:ayb_fs_linked_mcmc.xml,ayb_fs/stage3/ayb_fs_linked_mcmc_run1.xml  # Derived. Comma separated list of the stage 3 output files names.
old_fsmcmcoutputvec:ayb_fs/stage3/ayb_fs_linked_mcmc_x50000_y50000_z100.xml,ayb_fs/stage3/ayb_fs_linked_mcmc_run1_x50000_y50000_z100.xml  # Derived. Comma separated list of the stage 3 output files names, if we need to continue a too-short MCMC run.

###################
stage:4
## 

fstreeoutputvec:ayb_fs_linked_tree.xml,ayb_fs/stage4/ayb_fs_linked_tree_run1.xml  # Derived. Comma separated list of the stage 4 output files names.

section:fsprojectend
stage:5  # Derived. Don't mess with this! The internal measure of which stage of processing we've reached. Change it via -reset or -duplicate.
