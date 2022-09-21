# DNAfold
This simulation code was developed to predicting 3D structure and stability for double/single-stranded DNAs in ion solutions based on the coarse-grained (CG) model developed by Ya-Zhou Shi, Zi-Chun Mu, etc. This CG model was extended from the RNA CG model introduced by Ya-Zhou Shi, Lei Jin, Chen-Jie Feng, Xunxun Wang, Ya-Lan Tan and Zhi-Jie Tan at Tan's group, to predict 3D structure, stability and salt effect for RNAs from their sequences.
# Usage under Linux
1. Keep the following files (DNA.c, initial_cof.dat, and P.dat (option)) in the same directory;
3. Compile: gcc -Wall DNA.c -o DNA -lm
4. Run: ./DNA, and typing from the keyboard according to the hit;
     or ./DNA <P.dat, Ensuring that the P.dat is under the directory and modifying the file based on its description (see below). 
# Descriptions of input and output files
a. Input files
1. initial_conf.dat: initial conformation as input;
2. P.dat: parameters for the simulation
     (1). Folding Anealing or at constant T? 1. Folding (annealing) 2. Folding (constant temperature)
        --> 1 (using the SA algorithm for structure)   ***highly Recommend***
        --> 2 (using constant temperature for structure)
     (2). The number of chain:
        --> 2 (for dsDNA)
     (3). Input the 1-th chain length (nt):
        --> *** (Input the first chain length, e.g.11)
     (4). Input the 1-th chain sequence:
	      --> *********** (Input the first chain sequence, e.g.CGGACAAGAAG, capital letters of ATCG)
     (5). Input the 2-th chain length (nt):
   	    --> *** (Input the second chain length, e.g.11)
     (6). Input the 2-th chain sequence:
	      --> *********** (Input the first chain sequence, e.g.CTTCTTGTCCG, capital letters)
     (7). Including salt? yes or no 0. no 1. yes
	      --> 0 or 1
	   (8). if salt=0, input t0 (no salt)
  	 	  --> ** (Initial temperature (℃)), 
          if salt=1, input CNa CMg t0 (with salt)
  	   	--> ** ** ** ([Na+] (mM), [Mg2+] (mM), Initial temperature)
     (9). Input the chain concentration (mM)
	      --> *** (chain concentration)
(b) Output files
1. conf.dat: predicted conformations at different MC steps;
2. BP.dat: predicted base pairs at different MC steps;
3. U.dat: energy of predicted conformations at different MC steps;
4. jg.dat: radius of gyration, end-to-end distance, and persistent length of predicted conformations at different MC steps;
5. para.dat: some information for this running, including start time, output frequency, etc.
# Examples
Two examples for single- and double-stranded DNAs, respectively.
# References
1.	Shi, Y.Z., Wang, F.H., Wu, Y.Y. and Tan, Z.J. (2014) A coarse-grained model with implicit salt for RNAs: predicting 3D structure, stability and salt effect. J. Chem. Phys., 141, 105102.
2.	Shi, Y.Z., Jin, L., Wang, F.H., Zhu, X.L. and Tan, Z.J. (2015) Predicting 3D structure, flexibility, and stability of RNA hairpins in monovalent and divalent ion solutions. Biophys. J., 109, 2654-2665.
3.	Jin, L., Tan, Y.L., Wu, Y., Wang, X., Shi, Y.Z. and Tan, Z.J. (2019) Structure folding of RNA kissing complexes in salt solutions: predicting 3D structure, stability, and folding pathway. RNA, 25, 1532-1548.
4.	Shi, Y.Z., Jin, L., Feng, C.J., Tan, Y.L. and Tan, Z.J. (2018) Predicting 3D structure and stability of RNA pseudoknots in monovalent and divalent ion solutions. PLoS Comput. Biol., 14, e1006222.
5.	Jin, L., Shi, Y.Z., Feng, C.J., Tan, Y.L. and Tan, Z.J. (2018) Modeling structure, stability, and flexibility of double-stranded RNAs in salt solutions. Biophys. J., 115, 1403-1416.
