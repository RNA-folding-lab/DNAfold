# DNAfold
This simulation code was developed to predict 3D structure and stability for double/single-stranded DNAs in ion solutions based on the coarse-grained (CG) model developed by Ya-Zhou Shi, Zi-Chun Mu, etc. This CG model was extended from the RNA CG model introduced by Ya-Zhou Shi, Lei Jin, Chen-Jie Feng, Xunxun Wang, Ya-Lan Tan and Zhi-Jie Tan at Tan's group, to predict 3D structure, stability and salt effect for RNAs from their sequences.
# Usage under Linux
1. Keep the following files (DNA.c, initial_cof.dat, and P.dat (option)) in the same directory;
2. Create a new folder named as "results" in the directory;
3. Compile: gcc -Wall DNA.c -o DNA -lm
4. Run: ./DNA, and typing from the keyboard according to the hit;
##or ./DNA <P.dat, Ensuring that the P.dat is under the directory and modifying the file based on its description (see below). 
# Descriptions of input and output files
(a). Input files
1. initial_conf.dat: initial conformation as input;
2. P.dat: parameters for the simulation

   I. Folding Anealing or at constant T? 1. Folding (annealing) 2. Folding (constant temperature)
   --> 1 (using the SA algorithm for structure)   ***highly Recommend***
   --> 2 (using constant temperature for structure.
   II. The number of chain:
   --> 2 (for dsDNA)
   III. Input the 1-th chain length (nt):
   --> *** (Input the first chain length, e.g.11)
   IV. Input the 1-th chain sequence:
   --> *********** (Input the first chain sequence, e.g.CGGACAAGAAG, capital letters of ATCG)
   V. Input the 2-th chain length (nt):
   --> *** (Input the second chain length, e.g.11)
   VI. Input the 2-th chain sequence:
   --> *********** (Input the first chain sequence, e.g.CTTCTTGTCCG, capital letters)
   VII. Including salt? yes or no 0. no 1. yes
   --> 0 or 1
   VIII. if salt=0, input t0 (no salt)
   --> ** (Initial temperature (℃)), 
   if salt=1, input CNa CMg t0 (with salt)
   --> ** ** ** ([Na+] (mM), [Mg2+] (mM), Initial temperature)
   IX. Input the chain concentration (mM) 
   --> *** (chain concentration)

(b) Output files in results/
1. conf.dat: predicted conformations at different MC steps;
2. BP.dat: predicted base pairs at different MC steps;
3. U.dat: energy of predicted conformations at different MC steps;
4. jg.dat: radius of gyration, end-to-end distance, and persistent length of predicted conformations at different MC steps;
5. Secondary.dat: secondary structure of predicted conformations at different MC steps;
6. tt.dat: average of number of base pairs as function of temperature;
7. para.dat: some information for this running, including start time, output frequency, etc.
# Examples
Two examples for single- and double-stranded DNAs, respectively.
1. ssDNA: 1JVE-27nt hairpin (5'-CCTAATTATAACGAAGTTATAATTAGG-3'); 
2. dsDNA: 1QSK-duplex with bulge loop (strand1: 17nt, 5'-GCATCGAAAAAGCTACG-3'; strand2: 12nt, 5'-CGTAGCCGATGC-3')
3. annealing from 100℃ to 20℃, slat=0
4. run: ./DNA<P.dat
# References
1.	Shi, Y.Z., Wang, F.H., Wu, Y.Y. and Tan, Z.J. (2014) A coarse-grained model with implicit salt for RNAs: predicting 3D structure, stability and salt effect. J. Chem. Phys., 141, 105102.
2.	Shi, Y.Z., Jin, L., Wang, F.H., Zhu, X.L. and Tan, Z.J. (2015) Predicting 3D structure, flexibility, and stability of RNA hairpins in monovalent and divalent ion solutions. Biophys. J., 109, 2654-2665.
3.	Jin, L., Tan, Y.L., Wu, Y., Wang, X., Shi, Y.Z. and Tan, Z.J. (2019) Structure folding of RNA kissing complexes in salt solutions: predicting 3D structure, stability, and folding pathway. RNA, 25, 1532-1548.
4.	Shi, Y.Z., Jin, L., Feng, C.J., Tan, Y.L. and Tan, Z.J. (2018) Predicting 3D structure and stability of RNA pseudoknots in monovalent and divalent ion solutions. PLoS Comput. Biol., 14, e1006222.
5.	Jin, L., Shi, Y.Z., Feng, C.J., Tan, Y.L. and Tan, Z.J. (2018) Modeling structure, stability, and flexibility of double-stranded RNAs in salt solutions. Biophys. J., 115, 1403-1416.
