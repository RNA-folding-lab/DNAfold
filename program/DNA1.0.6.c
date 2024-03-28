/*This program developed by Ya-Zhou Shi & Zi-Chun Mu  (2022.2.15)
  to predict 3D structure and stability for double-stranded DNAs in ion solutions, 
  based on the program for RNAs present by Ya-Zhou Shi, Lei Jin & Zhi-Jie Tan;
  named as DNA.c  (Simulated Annealing)
  Here is the new version DNA1.0.1 (2023.7.18 for code optimization)
$ DNA1.0.3: ss-stacking only for unpairing bases
$ DNA1.0.5: Modify the distance between stacking N-N & entire program;
$ DNA1.0.6: modify the ss_stacking potential
$ gcc -Wall -o3 DNA1.0.6.c -o DNA -lm
$ ./DNA <P.dat or ./DNA
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#define pi 3.1415
#define step_pivot .5      //Translation step in Pivot move in Folding
#define step_3nt .8        //Translation step in 3nt fragment move in Folding
#define step_angle 0.5     //Angle step in Pivot_angle move in Folding
#define total1 10000000    //Max steps in Folding annealing (Structure predicton)
#define total2 200000000   //Max steps in Folding at constant T (melting)         
#define Alpha 0.98         //Annealing rate
#define Anneal 1           //Annealing shecdle
#define jt 10              //Cooling with constant T (e.g., jt=10C) 
#define tran 10            //rotational frequency
#define gamma 0.10         //If one base-pairing energy Ubp0<gamma*A, then base pair is formed
#define tconf  10000       //The frequency of conformational output
#define tenergy 1000       //The frequency of energy calculation
#define tprint 10000       //The frequency of screen output
#define t_si 1000          //The output frequency of structure information
#define p_move_3nt 0.0     //The relative probality of move_3nt
#define p_move_pivot_angle 0.5
#define max_size 500
#define chain_maxsize 3
#define Nstem_maxsize 40
#define Lstem_maxsize 30
float T,Kd,I,Ek,q4,CNa,CMg,fNa,D,Lbox,lt;
int   N,N0,total,Nbp,t,l,salt,nm,Folding,Energy,flag;
int   s[chain_maxsize][chain_maxsize][max_size][max_size],ss[chain_maxsize][chain_maxsize][max_size][max_size],a[chain_maxsize][max_size],c[chain_maxsize][max_size],bp,BP;
double  rand01;
float phi,theta,rm,step;
float rgp,rep,rlp1,rlp2,rlp3,PBP,pbp,G,GG,t0,M,m;
float x[chain_maxsize][max_size],y[chain_maxsize][max_size],z[chain_maxsize][max_size],xx[chain_maxsize][max_size],yy[chain_maxsize][max_size],zz[chain_maxsize][max_size],xxx[max_size],yyy[max_size],zzz[max_size],Q[chain_maxsize][max_size],R[chain_maxsize][max_size],f[chain_maxsize][max_size];
float U,Umin,Up,Up0,du,uu,u,ulj,uulj,uc,uuc,ub,uub,ue,uue,ud,uud,uN,uuN,us,uus,uco,uuco,uN1,uuN1,uss,uuss;
char  type[chain_maxsize][max_size],sequence[chain_maxsize][max_size];
int   Nchain,Lchain[chain_maxsize],L0_chain[chain_maxsize],chain_i0;
int   l_stem[Nstem_maxsize],stem_i[Nstem_maxsize][Lstem_maxsize],stem_j[Nstem_maxsize][Lstem_maxsize],N_stem,sstem_i[Nstem_maxsize][Lstem_maxsize],sstem_j[Nstem_maxsize][Lstem_maxsize],ll_stem[Nstem_maxsize],NN_stem,stem_i_c[Nstem_maxsize],stem_j_c[Nstem_maxsize],sstem_i_c[Nstem_maxsize],sstem_j_c[Nstem_maxsize];
/********The parameters (including bonded & nonbonded) of the CG Model************/
float e=0.26,A=-4.5,B0=-13.7,Beta=0.67,Beta1=0.0,SB0=-1.0,kss=11.0;
float kpc=98.2,lpc=3.95,kcp=70.5,lcp=3.95,kcN=45.8,lcN=3.55;                                      //Bond
float kpcp=9.8,Apcp=2.1,kcpc=8.6,Acpc=1.8,kpcN=3.2,ApcN=1.69,kNcp=7.1,ANcp=1.7;            //Angle
float kpcpc=0.28*1.3,dpcpc=2.5,kcpcp=0.28*4.0,dcpcp=-2.95,kcpcN=0.2*1.6,dcpcN=-1.3,kNcpc=0.2*0.8,dNcpc=0.9;  //Dihedral
float dNN1=8.6,dNN2=9.4,kNN=2.66,rNN=8.9,kCN=1.37,rCN=12.1,kPN=0.464,rPN=14.05;                //base-pairing
float sigma=4.53,dco=4.9,dkco=1.0,kdco=2.5;                                                     //base/coaxial-stacking
float Kter0=-0.3,Kter1=0.15,dter=8.8,dHmis=-3.435,dSmis=-8.24,Kmis=-1.0,ATend=0.25;  // Triplex & Mismatch & End penalty
float bl=5.45,IMg=3.0;
float ksNN=0.8,rsNN=4.8,ks=1.0;//ksCN=1.0,rsCN1=5.9,rsCN2=7.1,ksPC=1.0,rsPC=10.0;
/***********************************************************/
FILE *fp,*fp1,*fp2,*fp3,*fp5,*fp50,*fp6,*fp7;
int main()
{
 int   i,j,duo1,duo2,jj,Temp_len;
 float x00,y00,z00,R0,Q0,f0,x0[chain_maxsize][max_size],y0[chain_maxsize][max_size],z0[chain_maxsize][max_size];
 char  type_0;
 void  Set_Files(),MC_Annealing(),OutputPara(),Set_Box();

 Set_Files();         //defined the input parameters and out file names;
 OutputPara();        //the output of the parametes;

 i=1;j=1; jj=1; Temp_len=0;
 while(!feof(fp)) {   //Read the inital conformation (ch.dat) & distinguish each chain
      fscanf(fp,"%d %d %s %f %f %f %f %f %f\n",&duo1,&duo2,&type_0,&x00,&y00,&z00,&R0,&Q0,&f0);
      if((i-2)/3==strlen(sequence[j])) {
         Temp_len+=3*strlen(sequence[j])+1;
         if(L0_chain[j]!=i-1) printf("Error: \n"); 
         j++; i=1;}
      if(j>Nchain) break;
      if(jj>3*strlen(sequence[1])&&(jj==Temp_len+1||jj==Temp_len+2))  { }
      else {
         if(fmod(i,3)==0) {type_0=sequence[j][i/3-1];}
         if(fmod(i-1,3)==0) {Q0=-1.0;}
         x0[j][i]=x00; y0[j][i]=y00; z0[j][i]=z00; 
         type[j][i]=type_0; R[j][i]=R0; Q[j][i]=Q0;f[j][i]=f0;
//printf("%d %d %c %f %f %f %f %f %f\n",t,i,type[j][i],x0[j][i],y0[j][i],z0[j][i],R[j][i],Q[j][i],f[j][i]);
         i++;}
       jj++;
 }
 fclose(fp);        //input of the initial conformation 
        /*for (j=1;j<=Nchain;j++)  {
             for (i=1;i<=L0_chain[j];i++) {
                  printf("%d %d %c %f %f %f %f %f %f\n",t,i,type[j][i],x0[j][i],y0[j][i],z0[j][i],R[j][i],Q[j][i],f[j][i]);}
        }*/ 
 N0=3*N+j-1;          //N0:Total number of CG beads; N: Num. of nt
/****Add the simulation box, let the centre P of 1-th chain move to the centre of box*****/
 int N10,N1,outbox=0;
 float mv_x,mv_y,mv_z;
 N10=floor(L0_chain[1]/2)+1; 
 if (fmod(N10,3)==0) {N1=N10-2;}
 else if (fmod(N10+1,3)==0) {N1=N10-1;}
 else {N1=N10;}                         //N1:Fixed the P in centre of chain
 printf("N %d N0 %d N1 %d\n",N,N0,N1);  //The number of nucleotides & atoms & the actionless atom
 mv_x=Lbox/2.-x0[1][N1]; mv_y=Lbox/2.-y0[1][N1]; mv_z=Lbox/2.-z0[1][N1];
 for(j=1;j<=Nchain;j++) {
     for(i=1;i<=L0_chain[j];i++) { 
         x0[j][i]+=mv_x; y0[j][i]+=mv_y; z0[j][i]+=mv_z;
         if ((x0[j][i]<0.||x0[j][i]>Lbox)||(y0[j][i]<0.||y0[j][i]>Lbox)||(z0[j][i]<0.||z0[j][i]>Lbox)) {
            outbox=i;}
     }
     if (outbox>0) {printf("Warning: chain %d (%d) is out of the box\n",j,outbox);}
 }
/************** MC simulation  **************/
 MC_Annealing(x0,y0,z0);   
/********************************************/
 fclose(fp1);fclose(fp2);fclose(fp3); fclose(fp5);fclose(fp50); fclose(fp6); 
 return 0;
}  
/* &%$#@!~&%$#@!~&%$#@!~   Some functions or modules for move and calculation    &%$#@!~&%$#@!~&%$#@!~ */
void Set_Files (void) 
{
 printf("This program is designed for ssDNA/dsDNA. Please note the input\n");
 printf("Folding Anealing or at constant T?\n1. Folding (annealing)\n2. Folding (constant temperature)\n"); 
 scanf("%d",&Folding);            //Define that the Program is used for Folding or Optimizing?
/******Input the number of chain & length for each chain*******/
 printf("The number of chain:\n"); 
 scanf("%d",&Nchain);

 for (int i=1;i<=Nchain;i++) {
      printf("Input the %d-th chain length (nt):\n",i);
      scanf("%d",&Lchain[i]);
      N+=Lchain[i]; L0_chain[i]=3*Lchain[i]+1;
      printf("Input the %d-th chain sequence:\n",i);
      scanf("%s",sequence[i]); printf("%d %s\n",i,sequence[i]);
      if(Lchain[i]!=strlen(sequence[i])) {printf("Warning: the input length is not equal to the sequence length\n"); exit(0);}
 }
 printf("Including salt? yes or no\n0. no\n1. yes\n"); 
 scanf("%d",&salt); 
 if (salt==1) {
     printf("Input CNa CMg t0\n"); 
     scanf("%f %f %f",&CNa,&CMg,&t0);}
 else   {printf("Input t0\n"); scanf("%f",&t0);}     //Input the initial parametes & salt conditions.
 if (Nchain==1) {Nbp=0.5*Lchain[1];}
 else {Nbp=Lchain[1]; A=-2.6; B0=-9.6;}
 //N:number of nucleotides; Nbp:number of basepairs in Native structure; CNa/CMg:Na/Mg concentration; t0:initial Temperature; 
 float C_chain,Vbox;
 if (Nchain==1) {Lbox=500.0; C_chain=1.0;}
 else {
    printf("Input the chain concentration (mM)\n"); //Input the concentration of chain
    scanf("%f",&C_chain);
    Vbox=Nchain/(0.000602*C_chain); Lbox=floor(pow(Vbox,1.0/3.0)*10.); // The size of simulation box;
    printf("Vbox: %f; Lbox: %f (%f)\n",Vbox,Lbox,floor(pow(Nchain/(0.602*C_chain),1./3.)*100));
 }
 if (Folding==1)      {total=total1; lt=20.0;}      //Steps in Folding annealing
 else                 {total=total2; lt=t0;}        //Steps in Folding at constant temperature

 const char folder[20]="results/";  //file folder for simulated results
 if (access(folder, 0)!=0) {
    mkdir(folder,0755);  //if the results file is not exit, create it
    printf("File of results/ is created\n");}
 fp =fopen("initial_conf.dat","r+");               //initial conformation
 fp1=fopen("results/conf.dat","w+");             //output of conformations
 fp2=fopen("results/jg.dat","w+");               //Structure Information: Rg, Ree, Lp etc.
 fp3=fopen("results/BP.dat","w+");               //Num. of base-pairs in structures within simulations
 fp5=fopen("results/Secondary.dat","w+");         //Details of base-pairs 
 fp50=fopen("results/tt.dat","w+");              //The last mean information of the system
 fp6=fopen("results/U.dat","w+");                //The Energy of the conformations 
 fp7=fopen("results/para.dat","w+");             //Just like the log file for outputing the important paramaters. 
}
//******************************************//
void Parameters_T(float t01)
{
   float tt,qq4,qqq=0.0;
   tt=0.0+t01; T=273.15+tt*1.0; D=T*2.0*pow(10,-3);          //D: KT
   if (salt==1) {
      fNa=0.001*CNa/(0.001*CNa+(8.1-32.4/(N*0.5))*(5.2-log(0.001*CNa))*0.001*CMg); //The percentage of Na+ in Mixture in TBI_Helix
      Ek=87.740-0.4008*tt+9.398*1e-4*tt*tt-1.41*1e-6*tt*tt*tt;  //Permittivity
      I=CNa+IMg*CMg; Kd=sqrt((0.396*Ek*T)/I);                   //Ionic strength & Debye length
      qq4=5.998*1e-6*bl*Ek*T*0.5*(fNa+1);
      if (N*5.5<Kd) {
         qqq=log(Kd/bl)/log(N); if (qqq>1.) q4=qq4*qqq; else q4=qq4;} 
      else {q4=qq4;}           // q4=b/lB
   }
   printf("Temp: %f; I: %f; fNa: %f\n",tt,I,fNa); 
   printf("A: %f; B0: %f;\n",A,B0); 
   fprintf(fp7,"Temp %f T %f Kd %f I %f Ek %f q4 %f\n",tt,T,Kd,I,Ek,q4); 
   fflush(fp7);  //The output of important para.
}
/***********Initialize the parameters & calculate the information for initial conformation(ch.dat)****************/
/****************I don't know why this Ucs module cannot be moved to the end and just declared here***************/
float UCoS(int jc,int jc1,int jc2,int jc3,int i1,int j1,int k1,int k2,float x[chain_maxsize][max_size],float y[chain_maxsize][max_size],float z[chain_maxsize][max_size],float dco_1,float dco_2,float kdco_1,float kdco_2)
{
 float dco1,Ucos,kco,dco2,Ucos1=0.0,Ucos2=0.0;
 float BSt();
 kco=BSt(jc,jc1,i1,j1,k1,k2);
 dco1=sqrt((x[jc][i1]-x[jc2][k1])*(x[jc][i1]-x[jc2][k1])+(y[jc][i1]-y[jc2][k1])*(y[jc][i1]-y[jc2][k1])+(z[jc][i1]-z[jc2][k1])*(z[jc][i1]-z[jc2][k1])); 
 dco2=sqrt((x[jc3][k2]-x[jc1][j1])*(x[jc3][k2]-x[jc1][j1])+(y[jc3][k2]-y[jc1][j1])*(y[jc3][k2]-y[jc1][j1])+(z[jc3][k2]-z[jc1][j1])*(z[jc3][k2]-z[jc1][j1])); 
 Ucos1=-0.5*(kco-dkco)*((1-exp(-(dco1-dco_1)/kdco_1))*(1-exp(-(dco1-dco_1)/kdco_1))-1);
 Ucos2=-0.5*(kco-dkco)*((1-exp(-(dco2-dco_2)/kdco_2))*(1-exp(-(dco2-dco_2)/kdco_2))-1);
 Ucos=Ucos1+Ucos2; 
 return Ucos;
}
void Initialize_Para_T(float x0[chain_maxsize][max_size],float y0[chain_maxsize][max_size],float z0[chain_maxsize][max_size])
{ 
   rep=0.0;rgp=0.0; rlp1=0.0; rlp2=0.0; rlp3=0.0; //para for average of Rg & Ree 
   G=0.0;m=0;pbp=0.0;M=0;PBP=0.0;GG=0.0; //para for average of base-pairs
   Umin=1000.;Up0=0.0;Up=0.0;            //para for Energy
/****** Initialize the conformation of each temperature *******/
   int   i,j,jc,jc1,k1=0,Lij1,Lij2,Lij3,Lij4;
   float HB(),St();
   for (j=1;j<=Nchain;j++) {    //The initconf.=last one at previous T 
       for (i=1;i<=L0_chain[j];i++) {x[j][i]=x0[j][i];y[j][i]=y0[j][i];z[j][i]=z0[j][i];c[j][i]=0;}
   } 
/***********The energy of uN/us/uco for initial conformation************/
   bp=0; uN=0.0; us=0.0;  uco=0.0;
   for (jc=1;jc<=Nchain;jc++)  {
        for(i=1;i<=L0_chain[jc];i++) {
            if (fmod(i,3)==0) {
                for(jc1=jc;jc1<=Nchain;jc1++) {
                    for(j=1;j<=L0_chain[jc1];j++)  {
                        s[jc][jc1][i][j]=0;ss[jc][jc1][i][j]=0; 
                        if(((jc==jc1&&j>=i+12)||(jc<jc1))&&fmod(j,3)==0) {
                           if ((type[jc][i]=='G'&&type[jc1][j]=='C')||(type[jc][i]=='C'&&type[jc1][j]=='G')
                             ||(type[jc][i]=='A'&&type[jc1][j]=='T')||(type[jc][i]=='T'&&type[jc1][j]=='A')) {
                             uN1=0.0;
                             if (c[jc][i]==0&&c[jc1][j]==0)  {
                                 uN1=HB(jc,jc1,i,j,x,y,z);
                                 if (uN1!=0)  {c[jc][i]=1;c[jc1][j]=1; s[jc][jc1][i][j]=1; bp++; uN+=uN1;}
                             }
                          }                           
                        }
                    }
                }
             }
        }
   }
   fprintf(fp7,"The information for the initial conformation at Temperature (%f ℃ - %f K)\n",t0,T); fflush(fp7);  
   for (i=0;i<Nstem_maxsize;i++)     {
       l_stem[i]=1; stem_i_c[i]=0; stem_j_c[i]=0;                      //l_stem[i]:length of i-th stem
       for (j=0;j<Lstem_maxsize;j++) {stem_i[i][j]=0; stem_j[i][j]=0; } //The j-th base-pair (stem_i-stem_j) in stem i;
   } 
   for(jc=1;jc<=Nchain;jc++)  { 
      for(i=1;i<=L0_chain[jc];i++) { 
          if(fmod(i,3)==0) { 
             for(jc1=jc;jc1<=Nchain;jc1++) { 
                 for(j=1;j<=L0_chain[jc1];j++)  { 
                     if(((jc==jc1&&j>=i+12)||(jc<jc1))&&fmod(j,3)==0&&s[jc][jc1][i][j]==1) {  
                          if (s[jc][jc1][i-3][j+3]==1)  {l_stem[k1]++;}
                          else  {k1++;}
                          stem_i[k1][l_stem[k1]]=i;  stem_j[k1][l_stem[k1]]=j;  //nucleotides in stem[k1], stem_i paired with stem_j 
                          stem_i_c[k1]=jc; stem_j_c[k1]=jc1;   //chian's No. of the stem[k1]
                          fprintf(fp7,"stem_%d %d %d %d %d %c %d %c\n",k1,l_stem[k1],stem_i_c[k1],stem_j_c[k1],stem_i[k1][l_stem[k1]]/3,type[jc][i],stem_j[k1][l_stem[k1]]/3,type[jc1][j]); fflush(fp7); 
                            break;

                     }
                 }
              }
          }
      }
   }
   N_stem=k1;
   for(i=1;i<=N_stem;i++) {
       if (l_stem[i]>=2) {    //base-stacking
           for (j=1;j<l_stem[i];j++) {us+=St(stem_i_c[i],stem_j_c[i],stem_i[i][j],stem_j[i][j],x,y,z);}
           if (type[stem_i_c[i]][stem_i[i][1]]=='T'||type[stem_j_c[i]][stem_j[i][1]]=='T') {us+=ATend;}
           if (type[stem_i_c[i]][stem_i[i][l_stem[i]]]=='T'||type[stem_j_c[i]][stem_j[i][l_stem[i]]]=='T') {us+=ATend;}
       }
       if (l_stem[i]>=3) {   //coaxial-stacking
           for (j=i+1;j<=N_stem;j++)  {
                if (l_stem[j]>=3)     {
                    Lij1=(stem_i[j][1]-stem_i[i][l_stem[i]])/3; Lij2=(stem_j[j][1]-stem_j[i][l_stem[i]])/3;
                    Lij3=(stem_i[j][1]-stem_j[i][1])/3;   Lij4=(stem_j[i][l_stem[i]]-stem_i[j][l_stem[j]])/3;
                    if (stem_i_c[i]==stem_i_c[j]&&stem_j_c[i]==stem_j_c[j]) {
                        if (Lij1<5||Lij2<5)  {
                            if (Lij1==1||Lij2==1) {  //Bulge
                                uco+=UCoS(stem_i_c[i],stem_j_c[i],stem_i_c[j],stem_j_c[j],stem_i[i][l_stem[i]],stem_j[i][l_stem[i]],stem_i[j][1],stem_j[j][1],x,y,z,dco,dco,kdco,kdco);
                            }
                            else if (Lij2>=5) {  //Internal or Juntion
                                uco+=UCoS(stem_i_c[i],stem_j_c[i],stem_i_c[j],stem_j_c[j],stem_i[i][l_stem[i]],stem_j[i][l_stem[i]],stem_i[j][1],stem_j[j][1],x,y,z,Lij1*dco,Lij1*dco,kdco,kdco); 
                            }      
                            else if (Lij1>=5) {  //Internal or Juntion
                                uco+=UCoS(stem_i_c[i],stem_j_c[i],stem_i_c[j],stem_j_c[j],stem_i[i][l_stem[i]],stem_j[i][l_stem[i]],stem_i[j][1],stem_j[j][1],x,y,z,Lij2*dco,Lij2*dco,kdco,kdco); 
                            }  
                            else  {
                                uco+=UCoS(stem_i_c[i],stem_j_c[i],stem_i_c[j],stem_j_c[j],stem_i[i][l_stem[i]],stem_j[i][l_stem[i]],stem_i[j][1],stem_j[j][1],x,y,z,Lij1*dco,Lij2*dco,kdco,kdco); 
                            }               
                        }
                        else if (Lij3<5) {  //Juntion
                                uco+=UCoS(stem_i_c[i],stem_j_c[i],stem_i_c[j],stem_j_c[j],stem_i[i][1],stem_j[i][1],stem_j[j][1],stem_i[j][1],x,y,z,Lij3*dco,Lij3*dco,kdco,kdco);                            
                        }
                        else if (Lij4<3) {  //Pseudoknot: should be further improved!!!
                            uco+=UCoS(stem_i_c[i],stem_j_c[i],stem_i_c[j],stem_j_c[j],stem_i[i][l_stem[i]],stem_j[i][l_stem[i]],stem_j[j][l_stem[j]],stem_i[j][l_stem[j]],x,y,z,Lij4*dco,Lij4*dco,kdco,kdco);
                        }  
                        else {uco=uco;}
                    }
                    else if (stem_i_c[i]==stem_i_c[j]&&Lij1<5) {
                        uco+=UCoS(stem_i_c[i],stem_j_c[i],stem_i_c[j],stem_j_c[j],stem_i[i][l_stem[i]],stem_j[i][l_stem[i]],stem_i[j][1],stem_j[j][1],x,y,z,Lij1*dco,Lij1*dco,kdco,kdco);                     
                    }
                    else if (stem_j_c[i]==stem_j_c[j]) {
                        uco+=UCoS(stem_i_c[i],stem_j_c[i],stem_i_c[j],stem_j_c[j],stem_i[i][l_stem[i]],stem_j[i][l_stem[i]],stem_i[j][1],stem_j[j][1],x,y,z,Lij2*dco,Lij2*dco,kdco,kdco);                     
                    }
                    else {uco=uco;}
                }
           }
       }
  }
  fprintf(fp7,"Num of bp: %d; Basepairing: %f; BaseStacking: %f; Coaxial-Stacking: %f\n",bp,uN,us,uco); fflush(fp7); 
}
/****************Monte Carlo simulated annealing*********************/
void MC_Annealing(float MC_x0[chain_maxsize][max_size],float MC_y0[chain_maxsize][max_size],float MC_z0[chain_maxsize][max_size])
{
 int k,i,j;
 void MC_T();
 l=0;                   //MC steps independent of T
 for (k=0;k<=100;k++)   //Temperature cycle mechanism
{
  Parameters_T(t0);     //parameters at any t0: steps, T, Debye length, ionic strength, ion fraction etc.
  Initialize_Para_T(MC_x0,MC_y0,MC_z0);  //Parameters initialized at T, e.g., mean value or minimal value;   
/*****Do simulations at current temperature*****/
  MC_T();
  fprintf(fp50,"%f %f %f %f %f %f\n",t0,pbp,G,PBP,GG,Up); fflush(fp50);
/** Cool down the system and take the last conformation as the input for the next temperature **/
  for (j=1;j<=Nchain;j++) {
       for (i=1;i<=L0_chain[j];i++) {
            MC_x0[j][i]=x[j][i];MC_y0[j][i]=y[j][i];MC_z0[j][i]=z[j][i];}
  }    //Take the last conf. as the initial conf. for next T
  if (Anneal==0) {t0=273.15*(Alpha-1)+t0*Alpha;}  //Cool down the system as exp. speed
  if (Anneal==1) {t0=t0-jt;}                      //Cool down the system as constant T=jt;
  if (t0<lt) break;
 }     
}
/*****************Monte Carlo simulation at given Temperature******************************/
void MC_T()
{
 int i,ii,j;
 void MC_Each_Step(),StructureInformation(),ENERGY();
 float accept_ratio,p_accept_ratio;
 //srand((unsigned)time(NULL));           //Random number changes over time
 accept_ratio=0.0; p_accept_ratio=0.0;
 for(t=1;t<=total;t++)   {    
    if (fmod((t-1),100)==0) {nm=0;}     //nm: the Num. of acceptance each 100 steps;
    for(ii=1;ii<=N0;ii++)   {MC_Each_Step();}    //Do simulation at each step                                 
    if (fmod(t,tconf)==0) {    //Output of conformations in Folding
        for (j=1;j<=Nchain;j++)  {
             for (i=1;i<=L0_chain[j];i++) {
                  fprintf(fp1,"%d %d %c %f %f %f %f %f %f\n",t,i,type[j][i],x[j][i],y[j][i],z[j][i],R[j][i],Q[j][i],f[j][i]);}
        }
        fflush(fp1); 
    }
/***********Calculate energy & structure information (Ree, Rg and lp) of one conformation*******************/
    if (fmod(t,tenergy)==0) {ENERGY();}                //Calculate energy of one conformation        
    if (fmod(t,t_si)==0)    {StructureInformation();}  //Calculate some structural information (Ree, Rg and lp)
    if (fmod(t,100)==0)     {accept_ratio+=nm*1.0/(100.0*N0); p_accept_ratio=accept_ratio/(t/(100.0));}
    if (fmod(t,tprint)==0||t==1||t==10000)  {printf("%d %f %f %f %f %f %f (%f)\n",t,rgp,rep,rlp1,rlp2,rlp3,Up,p_accept_ratio);}	                      
   }     //End of the step t cycle; 
}
/*****************Monte Carlo for each step******************************/
void MC_Each_Step()
{
   int  atom_i0;
   void FOLD(),Metropolis(); 
   u=0.0;uu=0.0;                                    //Initialize the energies;
   ulj=0.0;ub=0.0;ue=0.0;ud=0.0;uc=0.0;uss=0.0; /*uN=0.0;us=0.0;uco=0.0;*/
   uulj=0.0;uub=0.0;uue=0.0;uud=0.0;uuN=0.0;uus=0.0;uuc=0.0;uuco=0.0;uuss=0.0;
   Energy=0;                      //if Energy=1, fuction ENERGY() is running 
   re:;
   rand01=rand()/(RAND_MAX+1.); chain_i0=floor(rand01*Nchain)+1;            //Random choose one chain;
   rand01=rand()/(RAND_MAX+1.); atom_i0=floor(rand01*L0_chain[chain_i0])+1;  //Random choose one atom at chain_i0;
   if(chain_i0>Nchain||atom_i0>L0_chain[chain_i0]) {printf("warning: %d %d %d %lf %lf\n",t,chain_i0,atom_i0,rand01,rand01*Nchain); goto re;} 
/***********   Folding   ***********/
   FOLD(atom_i0);                 //Folding Process
   Metropolis();
 }
/*********************************************/
void Metropolis(void)
{
   int i,j,jc,jc1;
   float p=0.0;
   //u=ulj+uN+us+uc+(ub+ue+ud)*0.5963+uco+kss*uss;           //Energy of conf. before changed
   //uu=uulj+uuN+uus+uuc+(uub+uue+uud)*0.5963+uuco+kss*uuss;  //Energy of conf. after changed
   u=ulj+uc+(ub+ue+ud)*0.5963+kss*uss;           //Energy of conf. before changed
   uu=uulj+uuc+(uub+uue+uud)*0.5963+kss*uuss;  //Energy of conf. after changed
   du=uu-u;  //du: Energy changes before & after moves;
//if(t==1) printf("~~ %d %f %f %f %f %f\n",t,ub+ue+ud,uub+uue+uud,uss,uuss,du);
   if(du<=0.0)   {     
      for (jc=1;jc<=Nchain;jc++) {
          for (i=1;i<=L0_chain[jc];i++) {
               x[jc][i]=xx[jc][i];y[jc][i]=yy[jc][i];z[jc][i]=zz[jc][i];
               c[jc][i]=a[jc][i];  
               for (jc1=1;jc1<=Nchain;jc1++) {
                    for(j=1;j<=L0_chain[jc1];j++) {s[jc][jc1][i][j]=ss[jc][jc1][i][j];} }
          }
      }  //To update the conf.               
      nm=nm+1; 
      uN=uuN; us=uus; uco=uuco; N_stem=NN_stem;   //Update the next initial energy for base-pairing/stacking 
      for (i=1;i<=N_stem;i++)  {                  //Update the next initial 2D structure information 
           l_stem[i]=ll_stem[i]; stem_i_c[i]=sstem_i_c[i]; stem_j_c[i]=sstem_j_c[i];
           for (j=1;j<=l_stem[i];j++) {stem_i[i][j]=sstem_i[i][j]; stem_j[i][j]=sstem_j[i][j];} }         
   }
   else { 
        rand01=rand()/(RAND_MAX+1.);  p=exp((-1)*du/D);
        if(rand01<=p)   {
           for (jc=1;jc<=Nchain;jc++) {
                for (i=1;i<=L0_chain[jc];i++) {
                     x[jc][i]=xx[jc][i];y[jc][i]=yy[jc][i];z[jc][i]=zz[jc][i];
                     c[jc][i]=a[jc][i];  
                     for (jc1=1;jc1<=Nchain;jc1++) {
                          for(j=1;j<=L0_chain[jc1];j++) {s[jc][jc1][i][j]=ss[jc][jc1][i][j];} }
                }
           }  //To update the conf.               
           nm=nm+1; 
           uN=uuN; us=uus; uco=uuco; N_stem=NN_stem;   //Update the next initial energy for base-pairing/stacking 
           for (i=1;i<=N_stem;i++)  {                  //Update the next initial 2D structure information 
                l_stem[i]=ll_stem[i]; stem_i_c[i]=sstem_i_c[i]; stem_j_c[i]=sstem_j_c[i];
                for (j=1;j<=l_stem[i];j++) {stem_i[i][j]=sstem_i[i][j]; stem_j[i][j]=sstem_j[i][j];} }    
        }
        else     {
           for (jc=1;jc<=Nchain;jc++) {
                for (i=1;i<=L0_chain[jc];i++) {
                     x[jc][i]=x[jc][i];y[jc][i]=y[jc][i];z[jc][i]=z[jc][i];
                     c[jc][i]=c[jc][i];  
                     for (jc1=1;jc1<=Nchain;jc1++) {
                          for(j=1;j<=L0_chain[jc1];j++) {s[jc][jc1][i][j]=s[jc][jc1][i][j];} }
                }
           }  //To update the conf.               
           uN=uN; us=us; uco=uco; N_stem=N_stem;   //Update the next initial energy for base-pairing/stacking 
           for (i=1;i<=N_stem;i++)  {                  //Update the next initial 2D structure information 
                l_stem[i]=l_stem[i]; stem_i_c[i]=stem_i_c[i]; stem_j_c[i]=stem_j_c[i];
                for (j=1;j<=l_stem[i];j++) {stem_i[i][j]=stem_i[i][j]; stem_j[i][j]=stem_j[i][j];} }    
        }
   }
 }
//********************  Move for translation & Pivot  ***************/
void Rand01()       //Generate rodom Euler angle;
{
 rand01=rand()/(RAND_MAX+1.);phi=rand01*pi;  
 rand01=rand()/(RAND_MAX+1.);theta=rand01*2.*pi; 
 rand01=rand()/(RAND_MAX+1.);rm=rand01*2.*pi;
 rand01=rand()/(RAND_MAX+1.);
}
/**********************Translation***********************************/
void Translate(int jc,int i1) //Translation of one atom;
{
 xxx[i1]=x[jc][i1]+step*rand01*sin(phi)*cos(theta);
 yyy[i1]=y[jc][i1]+step*rand01*sin(phi)*sin(theta);
 zzz[i1]=z[jc][i1]+step*rand01*cos(phi);
}
/**********************Pivot Move************************************/
void Pivot(int jc,int i1,int i10)  //Pivot moves for one segment;
{   //Euler rotation Z-X-Z;
 xx[jc][i1]=(xxx[i1]-xxx[i10])*(cos(theta)*cos(rm)-cos(phi)*sin(theta)*sin(rm))+(yyy[i1]-yyy[i10])*(sin(theta)*cos(rm)+cos(phi)*cos(theta)*sin(rm))+(zzz[i1]-zzz[i10])*sin(phi)*sin(rm)+xxx[i10];
 yy[jc][i1]=-(xxx[i1]-xxx[i10])*(cos(theta)*sin(rm)+cos(phi)*sin(theta)*cos(rm))+(yyy[i1]-yyy[i10])*(cos(phi)*cos(theta)*cos(rm)-sin(theta)*sin(rm))+(zzz[i1]-zzz[i10])*sin(phi)*cos(rm)+yyy[i10];
 zz[jc][i1]=(xxx[i1]-xxx[i10])*sin(phi)*sin(theta)-sin(phi)*cos(theta)*(yyy[i1]-yyy[i10])+(zzz[i1]-zzz[i10])*cos(phi)+zzz[i10];
}
/*********************Side-Chian move=translate&rotate*****/
void MoveN(int jc,int i1) //Movement of each base;
{
  Rand01();  Translate(jc,i1);
  Rand01();  Pivot(jc,i1,i1-1);
}
/***************GaussianNum()-2020.6.1******************/
float GaussianNum()
{
 float v1,v2,v3,v;
 do {
    v1=2.*rand()/(RAND_MAX+1.)-1.; v2=2.*rand()/(RAND_MAX+1.)-1.;
    v3=v1*v1+v2*v2; 
 }  while(v3>=1||v3==0);
 v=v1*sqrt(-2.0*log(v3)/v3);       
 return v;
}
/***************Rotate by delt_angle-2020.6.1******************/
void Pivot_angle(int jc,int i1,int i2,int i3,float delt_angle)
{   //vector i1-i3 rotated by vector i1-i2 with delt_angle
 float norm,nx,ny,nz,tan_A,f_A,Ru,Rv,Rw,R00,R01,R02,R10,R11,R12,R20,R21,R22,d_A=delt_angle;
 float r_x,r_y,r_z;
 nx=x[jc][i2]-x[jc][i1];  ny=y[jc][i2]-y[jc][i1];  nz=z[jc][i2]-z[jc][i1];
 norm=1./sqrt(nx*nx+ny*ny+nz*nz);
 nx*=norm;  ny*=norm;  nz*=norm;  
 tan_A=tan(d_A*0.5);
 Ru=tan_A*nx; Rv=tan_A*ny; Rw=tan_A*nz; f_A=2./(1+tan_A*tan_A);
 R00=0.5*f_A*(1+Ru*Ru-Rv*Rv-Rw*Rw); R01=f_A*(Ru*Rv-Rw); R02=f_A*(Ru*Rw+Rv);
 R10=f_A*(Ru*Rv+Rw);  R11=0.5*f_A*(1-Ru*Ru+Rv*Rv-Rw*Rw); R12= f_A*(Rv*Rw-Ru);
 R20=f_A*(Ru*Rw-Rv);  R21=f_A*(Rv*Rw+Ru);  R22=0.5*f_A*(1-Ru*Ru-Rv*Rv+Rw*Rw); 
 r_x=x[jc][i3]-x[jc][i1];  
 r_y=y[jc][i3]-y[jc][i1];  
 r_z=z[jc][i3]-z[jc][i1];
 xx[jc][i3]=R00*r_x+R01*r_y+R02*r_z+x[jc][i1]; 
 yy[jc][i3]=R10*r_x+R11*r_y+R12*r_z+y[jc][i1]; 
 zz[jc][i3]=R20*r_x+R21*r_y+R22*r_z+z[jc][i1];
}
/****************************Functions********************************/
/* ~~~~~~~~~~ Details of bonded potential calculation ~~~~~~~~~~~~~~ */
float PC(int i1,int jc,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float d,ul;
 d=sqrt((x1[jc][i1-2]-x1[jc][i1-1])*(x1[jc][i1-2]-x1[jc][i1-1])+(y1[jc][i1-2]-y1[jc][i1-1])*(y1[jc][i1-2]-y1[jc][i1-1])+(z1[jc][i1-2]-z1[jc][i1-1])*(z1[jc][i1-2]-z1[jc][i1-1]));
 ul=kpc*(d-lpc)*(d-lpc);
 return ul;
}
float CP(int i1,int jc,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float d,ul;
 d=sqrt((x1[jc][i1+1]-x1[jc][i1-1])*(x1[jc][i1+1]-x1[jc][i1-1])+(y1[jc][i1+1]-y1[jc][i1-1])*(y1[jc][i1+1]-y1[jc][i1-1])+(z1[jc][i1+1]-z1[jc][i1-1])*(z1[jc][i1+1]-z1[jc][i1-1]));
 ul=kcp*(d-lcp)*(d-lcp);
 return ul;
}
float CN(int i1,int jc,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float d,ul;
 d=sqrt((x1[jc][i1]-x1[jc][i1-1])*(x1[jc][i1]-x1[jc][i1-1])+(y1[jc][i1]-y1[jc][i1-1])*(y1[jc][i1]-y1[jc][i1-1])+(z1[jc][i1]-z1[jc][i1-1])*(z1[jc][i1]-z1[jc][i1-1]));
 ul=kcN*(d-lcN)*(d-lcN);
 return ul;
}
float PCP(int i1,int jc,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float d1,d2,d3,w,a1,ue0;
 d1=sqrt((x1[jc][i1-2]-x1[jc][i1-1])*(x1[jc][i1-2]-x1[jc][i1-1])+(y1[jc][i1-2]-y1[jc][i1-1])*(y1[jc][i1-2]-y1[jc][i1-1])+(z1[jc][i1-2]-z1[jc][i1-1])*(z1[jc][i1-2]-z1[jc][i1-1]));
 d2=sqrt((x1[jc][i1-1]-x1[jc][i1+1])*(x1[jc][i1-1]-x1[jc][i1+1])+(y1[jc][i1-1]-y1[jc][i1+1])*(y1[jc][i1-1]-y1[jc][i1+1])+(z1[jc][i1-1]-z1[jc][i1+1])*(z1[jc][i1-1]-z1[jc][i1+1]));
 d3=sqrt((x1[jc][i1+1]-x1[jc][i1-2])*(x1[jc][i1+1]-x1[jc][i1-2])+(y1[jc][i1+1]-y1[jc][i1-2])*(y1[jc][i1+1]-y1[jc][i1-2])+(z1[jc][i1+1]-z1[jc][i1-2])*(z1[jc][i1+1]-z1[jc][i1-2]));
 w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
 if (w<=-1.0) {a1=3.14;}
 else if (w>=1.0) {a1=0.;}
 else  {a1=acos(w);}
 ue0=kpcp*(a1-Apcp)*(a1-Apcp);  
 return ue0;
}
float CPC(int i1,int jc,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float d1,d2,d3,w,a1,ue0;
 if(i1+2>N0)  ue0=0.0;
 else {
 d1=sqrt((x1[jc][i1-1]-x1[jc][i1+1])*(x1[jc][i1-1]-x1[jc][i1+1])+(y1[jc][i1-1]-y1[jc][i1+1])*(y1[jc][i1-1]-y1[jc][i1+1])+(z1[jc][i1-1]-z1[jc][i1+1])*(z1[jc][i1-1]-z1[jc][i1+1]));
 d2=sqrt((x1[jc][i1+1]-x1[jc][i1+2])*(x1[jc][i1+1]-x1[jc][i1+2])+(y1[jc][i1+1]-y1[jc][i1+2])*(y1[jc][i1+1]-y1[jc][i1+2])+(z1[jc][i1+1]-z1[jc][i1+2])*(z1[jc][i1+1]-z1[jc][i1+2]));
 d3=sqrt((x1[jc][i1+2]-x1[jc][i1-1])*(x1[jc][i1+2]-x1[jc][i1-1])+(y1[jc][i1+2]-y1[jc][i1-1])*(y1[jc][i1+2]-y1[jc][i1-1])+(z1[jc][i1+2]-z1[jc][i1-1])*(z1[jc][i1+2]-z1[jc][i1-1]));
 w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
 if (w<=-1.0) {a1=3.14;}
 else if (w>=1.0) {a1=0.;}
 else  {a1=acos(w);}
 ue0=kcpc*(a1-Acpc)*(a1-Acpc);  
 } 
 return ue0;
}
float PCN(int i1,int jc,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float d1,d2,d3,w,a1,ue0;
 d1=sqrt((x1[jc][i1-1]-x1[jc][i1-2])*(x1[jc][i1-1]-x1[jc][i1-2])+(y1[jc][i1-1]-y1[jc][i1-2])*(y1[jc][i1-1]-y1[jc][i1-2])+(z1[jc][i1-1]-z1[jc][i1-2])*(z1[jc][i1-1]-z1[jc][i1-2]));
 d2=sqrt((x1[jc][i1-1]-x1[jc][i1])*(x1[jc][i1-1]-x1[jc][i1])+(y1[jc][i1-1]-y1[jc][i1])*(y1[jc][i1-1]-y1[jc][i1])+(z1[jc][i1-1]-z1[jc][i1])*(z1[jc][i1-1]-z1[jc][i1]));
 d3=sqrt((x1[jc][i1-2]-x1[jc][i1])*(x1[jc][i1-2]-x1[jc][i1])+(y1[jc][i1-2]-y1[jc][i1])*(y1[jc][i1-2]-y1[jc][i1])+(z1[jc][i1-2]-z1[jc][i1])*(z1[jc][i1-2]-z1[jc][i1]));
 w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
 if (w<=-1.0) {a1=3.14;}
 else if (w>=1.0) {a1=0.;}
 else  {a1=acos(w);}
 ue0=kpcN*(a1-ApcN)*(a1-ApcN);
 return ue0;
}
float NCP(int i1,int jc,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float d1,d2,d3,w,a1,ue0;
 d1=sqrt((x1[jc][i1-1]-x1[jc][i1])*(x1[jc][i1-1]-x1[jc][i1])+(y1[jc][i1-1]-y1[jc][i1])*(y1[jc][i1-1]-y1[jc][i1])+(z1[jc][i1-1]-z1[jc][i1])*(z1[jc][i1-1]-z1[jc][i1]));
 d2=sqrt((x1[jc][i1-1]-x1[jc][i1+1])*(x1[jc][i1-1]-x1[jc][i1+1])+(y1[jc][i1-1]-y1[jc][i1+1])*(y1[jc][i1-1]-y1[jc][i1+1])+(z1[jc][i1-1]-z1[jc][i1+1])*(z1[jc][i1-1]-z1[jc][i1+1]));
 d3=sqrt((x1[jc][i1+1]-x1[jc][i1])*(x1[jc][i1+1]-x1[jc][i1])+(y1[jc][i1+1]-y1[jc][i1])*(y1[jc][i1+1]-y1[jc][i1])+(z1[jc][i1+1]-z1[jc][i1])*(z1[jc][i1+1]-z1[jc][i1]));
 w=(d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
 if (w<=-1.0) {a1=3.14;}
 else if (w>=1.0) {a1=0.;}
 else  {a1=acos(w);}
 ue0=kNcp*(a1-ANcp)*(a1-ANcp);
 return ue0;
}
float PCPC(int i1,int jc,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
 if(i1+2>L0_chain[jc]) ud0=0.0;
 else  {
 c1=((y1[jc][i1-2]-y1[jc][i1-1])*(z1[jc][i1-1]-z1[jc][i1+1])-(z1[jc][i1-2]-z1[jc][i1-1])*(y1[jc][i1-1]-y1[jc][i1+1]));
 c2=((z1[jc][i1-2]-z1[jc][i1-1])*(x1[jc][i1-1]-x1[jc][i1+1])-(x1[jc][i1-2]-x1[jc][i1-1])*(z1[jc][i1-1]-z1[jc][i1+1]));
 c3=((x1[jc][i1-2]-x1[jc][i1-1])*(y1[jc][i1-1]-y1[jc][i1+1])-(y1[jc][i1-2]-y1[jc][i1-1])*(x1[jc][i1-1]-x1[jc][i1+1]));
 p1=((y1[jc][i1-1]-y1[jc][i1+1])*(z1[jc][i1+1]-z1[jc][i1+2])-(z1[jc][i1-1]-z1[jc][i1+1])*(y1[jc][i1+1]-y1[jc][i1+2]));
 p2=((z1[jc][i1-1]-z1[jc][i1+1])*(x1[jc][i1+1]-x1[jc][i1+2])-(x1[jc][i1-1]-x1[jc][i1+1])*(z1[jc][i1+1]-z1[jc][i1+2]));
 p3=((x1[jc][i1-1]-x1[jc][i1+1])*(y1[jc][i1+1]-y1[jc][i1+2])-(y1[jc][i1-1]-y1[jc][i1+1])*(x1[jc][i1+1]-x1[jc][i1+2]));
 e1=sqrt(c1*c1+c2*c2+c3*c3); f1=sqrt(p1*p1+p2*p2+p3*p3);
 pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
 g1=(x1[jc][i1-2]-x1[jc][i1+2]); g2=(y1[jc][i1-2]-y1[jc][i1+2]); g3=(z1[jc][i1-2]-z1[jc][i1+2]);
 gg1=sqrt(g1*g1+g2*g2+g3*g3); hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
 if (pp1<=-1.0) {di=-3.14;}
 else if (pp1>=1.0) {di=0.;}
 else if (hh1>=0.) {di=acos(pp1);}
 else {di=-acos(pp1);}
 ud0=kpcpc*((1-cos(di-dpcpc))+0.5*(1-cos(3.*(di-dpcpc)))); }
 return ud0;
}
float CPCP(int i1,int jc,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
 if(i1+2>L0_chain[jc]) ud0=0.0;
 else  {
 c1=((y1[jc][i1-1]-y1[jc][i1+1])*(z1[jc][i1+1]-z1[jc][i1+2])-(z1[jc][i1-1]-z1[jc][i1+1])*(y1[jc][i1+1]-y1[jc][i1+2]));
 c2=((z1[jc][i1-1]-z1[jc][i1+1])*(x1[jc][i1+1]-x1[jc][i1+2])-(x1[jc][i1-1]-x1[jc][i1+1])*(z1[jc][i1+1]-z1[jc][i1+2]));
 c3=((x1[jc][i1-1]-x1[jc][i1+1])*(y1[jc][i1+1]-y1[jc][i1+2])-(y1[jc][i1-1]-y1[jc][i1+1])*(x1[jc][i1+1]-x1[jc][i1+2]));
 p1=((y1[jc][i1+1]-y1[jc][i1+2])*(z1[jc][i1+2]-z1[jc][i1+4])-(z1[jc][i1+1]-z1[jc][i1+2])*(y1[jc][i1+2]-y1[jc][i1+4]));
 p2=((z1[jc][i1+1]-z1[jc][i1+2])*(x1[jc][i1+2]-x1[jc][i1+4])-(x1[jc][i1+1]-x1[jc][i1+2])*(z1[jc][i1+2]-z1[jc][i1+4]));
 p3=((x1[jc][i1+1]-x1[jc][i1+2])*(y1[jc][i1+2]-y1[jc][i1+4])-(y1[jc][i1+1]-y1[jc][i1+2])*(x1[jc][i1+2]-x1[jc][i1+4]));
 e1=sqrt(c1*c1+c2*c2+c3*c3); f1=sqrt(p1*p1+p2*p2+p3*p3);
 pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
 g1=(x1[jc][i1-1]-x1[jc][i1+4]); g2=(y1[jc][i1-1]-y1[jc][i1+4]); g3=(z1[jc][i1-1]-z1[jc][i1+4]);
 gg1=sqrt(g1*g1+g2*g2+g3*g3); hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
 if (pp1<=-1.0) {di=-3.14;}
 else if (pp1>=1.0) {di=0.;}
 else if (hh1>=0.) {di=acos(pp1);}
 else {di=-acos(pp1);}
 ud0=kcpcp*((1-cos(di-dcpcp))+0.5*(1-cos(3.*(di-dcpcp))));  }
 return ud0;
}
float CPCN(int i1,int jc,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
 if(i1-4<0) ud0=0.0;
 else  {
 c1=((y1[jc][i1-4]-y1[jc][i1-2])*(z1[jc][i1-2]-z1[jc][i1-1])-(z1[jc][i1-4]-z1[jc][i1-2])*(y1[jc][i1-2]-y1[jc][i1-1]));
 c2=((z1[jc][i1-4]-z1[jc][i1-2])*(x1[jc][i1-2]-x1[jc][i1-1])-(x1[jc][i1-4]-x1[jc][i1-2])*(z1[jc][i1-2]-z1[jc][i1-1]));
 c3=((x1[jc][i1-4]-x1[jc][i1-2])*(y1[jc][i1-2]-y1[jc][i1-1])-(y1[jc][i1-4]-y1[jc][i1-2])*(x1[jc][i1-2]-x1[jc][i1-1]));
 p1=((y1[jc][i1-2]-y1[jc][i1-1])*(z1[jc][i1-1]-z1[jc][i1])-(z1[jc][i1-2]-z1[jc][i1-1])*(y1[jc][i1-1]-y1[jc][i1]));
 p2=((z1[jc][i1-2]-z1[jc][i1-1])*(x1[jc][i1-1]-x1[jc][i1])-(x1[jc][i1-2]-x1[jc][i1-1])*(z1[jc][i1-1]-z1[jc][i1]));
 p3=((x1[jc][i1-2]-x1[jc][i1-1])*(y1[jc][i1-1]-y1[jc][i1])-(y1[jc][i1-2]-y1[jc][i1-1])*(x1[jc][i1-1]-x1[jc][i1]));
 e1=sqrt(c1*c1+c2*c2+c3*c3); f1=sqrt(p1*p1+p2*p2+p3*p3);
 pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
 g1=(x1[jc][i1-4]-x1[jc][i1]); g2=(y1[jc][i1-4]-y1[jc][i1]); g3=(z1[jc][i1-4]-z1[jc][i1]);
 gg1=sqrt(g1*g1+g2*g2+g3*g3); hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
 if (pp1<=-1.0) {di=-3.14;}
 else if (pp1>=1.0) {di=0.;}
 else if (hh1>=0.) {di=acos(pp1);}
 else {di=-acos(pp1);}
 ud0=kcpcN*((1-cos(di-dcpcN))+0.5*(1-cos(3.*(di-dcpcN)))); }
 return ud0;
}
float NCPC(int i1,int jc,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float c1,c2,c3,p1,p2,p3,e1,f1,pp1,g1,g2,g3,gg1,hh1,di,ud0=0.0;
 if(i1+2>L0_chain[jc]) ud0=0.0;
 else  {
 c1=((y1[jc][i1]-y1[jc][i1-1])*(z1[jc][i1-1]-z1[jc][i1+1])-(z1[jc][i1]-z1[jc][i1-1])*(y1[jc][i1-1]-y1[jc][i1+1]));
 c2=((z1[jc][i1]-z1[jc][i1-1])*(x1[jc][i1-1]-x1[jc][i1+1])-(x1[jc][i1]-x1[jc][i1-1])*(z1[jc][i1-1]-z1[jc][i1+1]));
 c3=((x1[jc][i1]-x1[jc][i1-1])*(y1[jc][i1-1]-y1[jc][i1+1])-(y1[jc][i1]-y1[jc][i1-1])*(x1[jc][i1-1]-x1[jc][i1+1]));
 p1=((y1[jc][i1-1]-y1[jc][i1+1])*(z1[jc][i1+1]-z1[jc][i1+2])-(z1[jc][i1-1]-z1[jc][i1+1])*(y1[jc][i1+1]-y1[jc][i1+2]));
 p2=((z1[jc][i1-1]-z1[jc][i1+1])*(x1[jc][i1+1]-x1[jc][i1+2])-(x1[jc][i1-1]-x1[jc][i1+1])*(z1[jc][i1+1]-z1[jc][i1+2]));
 p3=((x1[jc][i1-1]-x1[jc][i1+1])*(y1[jc][i1+1]-y1[jc][i1+2])-(y1[jc][i1-1]-y1[jc][i1+1])*(x1[jc][i1+1]-x1[jc][i1+2]));
 e1=sqrt(c1*c1+c2*c2+c3*c3);
 f1=sqrt(p1*p1+p2*p2+p3*p3);
 pp1=(c1*p1+c2*p2+c3*p3)/(e1*f1);
 g1=(x1[jc][i1]-x1[jc][i1+2]); g2=(y1[jc][i1]-y1[jc][i1+2]); g3=(z1[jc][i1]-z1[jc][i1+2]);
 gg1=sqrt(g1*g1+g2*g2+g3*g3); hh1=(p1*g1+p2*g2+p3*g3)/(f1*gg1);
 if (pp1<=-1.0) {di=-3.14;}
 else if (pp1>=1.0) {di=0.;}
 else if (hh1>=0.) {di=acos(pp1);}
 else {di=-acos(pp1);}
 ud0=kNcpc*((1-cos(di-dNcpc))+0.5*(1-cos(3.*(di-dNcpc)))); }
 return ud0;
}
/*~~~~~~~~~~~~~~base pairing between two complementary bases (AU,GC,and GU)~~~~~~~~~~~~~~~*/
float HB(int jc1,int jc2,int i1,int j1,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
 float d,hb,UHB,d0,d01,d1,d11;
 d=sqrt((x1[jc1][i1]-x1[jc2][j1])*(x1[jc1][i1]-x1[jc2][j1])+(y1[jc1][i1]-y1[jc2][j1])*(y1[jc1][i1]-y1[jc2][j1])+(z1[jc1][i1]-z1[jc2][j1])*(z1[jc1][i1]-z1[jc2][j1])); //NN
 if (d>=dNN1&&d<=dNN2)  {  //The pairing formation condition
     d0=sqrt((x1[jc1][i1]-x1[jc2][j1-1])*(x1[jc1][i1]-x1[jc2][j1-1])+(y1[jc1][i1]-y1[jc2][j1-1])*(y1[jc1][i1]-y1[jc2][j1-1])+(z1[jc1][i1]-z1[jc2][j1-1])*(z1[jc1][i1]-z1[jc2][j1-1]));  //NiCj
     d01=sqrt((x1[jc2][j1]-x1[jc1][i1-1])*(x1[jc2][j1]-x1[jc1][i1-1])+(y1[jc2][j1]-y1[jc1][i1-1])*(y1[jc2][j1]-y1[jc1][i1-1])+(z1[jc2][j1]-z1[jc1][i1-1])*(z1[jc2][j1]-z1[jc1][i1-1]));  //CiNj
     d1=sqrt((x1[jc1][i1-2]-x1[jc2][j1])*(x1[jc1][i1-2]-x1[jc2][j1])+(y1[jc1][i1-2]-y1[jc2][j1])*(y1[jc1][i1-2]-y1[jc2][j1])+(z1[jc1][i1-2]-z1[jc2][j1])*(z1[jc1][i1-2]-z1[jc2][j1]));    //PiNj
     d11=sqrt((x1[jc2][j1-2]-x1[jc1][i1])*(x1[jc2][j1-2]-x1[jc1][i1])+(y1[jc2][j1-2]-y1[jc1][i1])*(y1[jc2][j1-2]-y1[jc1][i1])+(z1[jc2][j1-2]-z1[jc1][i1])*(z1[jc2][j1-2]-z1[jc1][i1])); //NiPj
     if ((type[jc1][i1]=='G'&&type[jc2][j1]=='C')||(type[jc1][i1]=='C'&&type[jc2][j1]=='G'))  {hb=A;}
else if ((type[jc1][i1]=='A'&&type[jc2][j1]=='T')||(type[jc1][i1]=='T'&&type[jc2][j1]=='A'))  {hb=Beta*A;} 
else {hb=Beta1*A;}   //mismatches,Beta1=0--without
     UHB=hb/(1+kNN*(d-rNN)*(d-rNN)+kCN*((d0-rCN)*(d0-rCN)+(d01-rCN)*(d01-rCN))+kPN*((d1-rPN)*(d1-rPN)+(d11-rPN)*(d11-rPN))); 
 }
 else UHB=0.0;
 return UHB;
}
 void BasePairing()
{
    int i,j,jc,jc1;
    /*for (jc=1;jc<=Nchain;jc++)  {
        for(i=1;i<=L0_chain[jc];i++) {a[jc][i]=0;
            for(jc1=1;jc1<=Nchain;jc1++) {
                for(j=1;j<=L0_chain[jc1];j++)  {ss[jc][jc1][i][j]=0;}}}
    }*/
    memset(a,0,sizeof(a)); memset(ss,0,sizeof(ss)); 
    for (jc=1;jc<=Nchain;jc++)  {
        for(i=1;i<=L0_chain[jc];i++) {
            if (fmod(i,3)==0) {
                for(jc1=1;jc1<=Nchain;jc1++) {
                    for(j=1;j<=L0_chain[jc1];j++)  {
                        if(((jc==jc1&&j>=i+12)||(jc<jc1))&&fmod(j,3)==0) {
                           if ((type[jc][i]=='G'&&type[jc1][j]=='C')||(type[jc][i]=='C'&&type[jc1][j]=='G')
                            ||(type[jc][i]=='A'&&type[jc1][j]=='T')||(type[jc][i]=='T'&&type[jc1][j]=='A')) {
                              uuN1=0.0;   
                              if (a[jc][i]==0&&a[jc1][j]==0) {
                                 uuN1=HB(jc,jc1,i,j,xx,yy,zz); 
                                 if (uuN1!=0) {a[jc][i]=1;a[jc1][j]=1; ss[jc][jc1][i][j]=1; }
                                 uuN+=uuN1; 
                             }
                          } 
                          else   {uuN=uuN;}                           
                        }
                   }
                }
             }
        }
    }
 }

// ~~~~~~~~~~~~~~~~Calculation of base stacking~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/
/*********BSt:experimental thermo-parameters for stacking*****************/
 float BSt(int jc,int jc1,int i1,int j1,int k1,int k2)  
{
   float B=0.0;   // 5'-i1-k1- 3'    SantaLucia and Hicks, Annu. Rev. Biophys. Struct. 2004, 33:415-440;
                  // 3'-j1/k2- 5'    SantaLucia et al. Biochemistry, 1996, 35:3555-3562.
     if (((type[jc][i1]=='A'&&type[jc][k1]=='A')&&(type[jc1][j1]=='T'&&type[jc1][k2]=='T'))
     ||((type[jc][i1]=='T'&&type[jc][k1]=='T')&&(type[jc1][j1]=='A'&&type[jc1][k2]=='A'))) {B=-7.6-T*0.001*(-21.3-B0);}   // AA/TT  (TT/AA)
else if (((type[jc][i1]=='A'&&type[jc][k1]=='C')&&(type[jc1][j1]=='T'&&type[jc1][k2]=='G')) 
     ||((type[jc][i1]=='G'&&type[jc][k1]=='T')&&(type[jc1][j1]=='C'&&type[jc1][k2]=='A'))) {B=-8.4-T*0.001*(-22.4-B0);}   // AC/TG  (GT/CA)
else if (((type[jc][i1]=='A'&&type[jc][k1]=='G')&&(type[jc1][j1]=='T'&&type[jc1][k2]=='C'))
     ||((type[jc][i1]=='C'&&type[jc][k1]=='T')&&(type[jc1][j1]=='G'&&type[jc1][k2]=='A'))) {B=-7.8-T*0.001*(-21.0-B0);}   // AG/TC  (CT/GA)
else if ((type[jc][i1]=='A'&&type[jc][k1]=='T')&&(type[jc1][j1]=='T'&&type[jc1][k2]=='A')) {B=-7.2-T*0.001*(-20.4-B0);}   // AT/TA
else if (((type[jc][i1]=='T'&&type[jc][k1]=='G')&&(type[jc1][j1]=='A'&&type[jc1][k2]=='C')) 
     ||((type[jc][i1]=='C'&&type[jc][k1]=='A')&&(type[jc1][j1]=='G'&&type[jc1][k2]=='T'))) {B=-8.5-T*0.001*(-22.7-B0);}   // TG/AC  (CA/GT)
else if (((type[jc][i1]=='G'&&type[jc][k1]=='G')&&(type[jc1][j1]=='C'&&type[jc1][k2]=='C'))
     ||((type[jc][i1]=='C'&&type[jc][k1]=='C')&&(type[jc1][j1]=='G'&&type[jc1][k2]=='G'))) {B=-8.0-T*0.001*(-19.9-B0);}   // GG/CC
else if ((type[jc][i1]=='C'&&type[jc][k1]=='G')&&(type[jc1][j1]=='G'&&type[jc1][k2]=='C')) {B=-10.6-T*0.001*(-27.2-B0);}  // CG/GC
else if (((type[jc][i1]=='T'&&type[jc][k1]=='C')&&(type[jc1][j1]=='A'&&type[jc1][k2]=='G'))
     ||((type[jc][i1]=='G'&&type[jc][k1]=='A')&&(type[jc1][j1]=='C'&&type[jc1][k2]=='T'))) {B=-8.2-T*0.001*(-22.2-B0);}   // TC/AG  (GA/CT)
else if ((type[jc][i1]=='G'&&type[jc][k1]=='C')&&(type[jc1][j1]=='C'&&type[jc1][k2]=='G')) {B=-9.8-T*0.001*(-24.4-B0);}   // GC/CG
else if ((type[jc][i1]=='T'&&type[jc][k1]=='A')&&(type[jc1][j1]=='A'&&type[jc1][k2]=='T')) {B=-7.2-T*0.001*(-21.3-B0);}   // TA/AT
else    {B=0.0;}   //Mismatched stacking
   return B;
}
/************************Potential for base-pair stacking*********************************/
 float St(int jc,int jc1,int i1,int j1,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
  float d,d0,kst,ULJ,d1,d5,d6,d10,d12,d01,d05,d06,d010,d012;
  d=sqrt((x1[jc][i1]-x1[jc][i1+3])*(x1[jc][i1]-x1[jc][i1+3])+(y1[jc][i1]-y1[jc][i1+3])*(y1[jc][i1]-y1[jc][i1+3])+(z1[jc][i1]-z1[jc][i1+3])*(z1[jc][i1]-z1[jc][i1+3])); 
  d1=sigma/d; d5=d1*d1*d1*d1*d1;d6=d1*d5;d10=d5*d5;d12=d6*d6; 
  d0=sqrt((x1[jc1][j1]-x1[jc1][j1-3])*(x1[jc1][j1]-x1[jc1][j1-3])+(y1[jc1][j1]-y1[jc1][j1-3])*(y1[jc1][j1]-y1[jc1][j1-3])+(z1[jc1][j1]-z1[jc1][j1-3])*(z1[jc1][j1]-z1[jc1][j1-3]));   
  d01=sigma/d0;d05=d01*d01*d01*d01*d01;d06=d01*d05;d010=d05*d05;d012=d06*d06;
  kst=BSt(jc,jc1,i1,j1,i1+3,j1-3);
  if (kst>=0) {ULJ=0.0;}
  else {ULJ=-0.5*kst*((5*d12-6*d10)+(5*d012-6*d010));}
  return ULJ;
} 
 void BaseStacking()
{
  int i,j;
  for (i=1;i<=NN_stem;i++) {
       if (ll_stem[i]>=2) {  
           for (j=1;j<ll_stem[i];j++) {uus+=St(sstem_i_c[i],sstem_j_c[i],sstem_i[i][j],sstem_j[i][j],xx,yy,zz);}
           if (type[stem_i_c[i]][sstem_i[i][1]]=='T'||type[stem_j_c[i]][sstem_j[i][1]]=='T') {uus+=ATend;}
           if (type[stem_i_c[i]][sstem_i[i][ll_stem[i]]]=='T'||type[stem_j_c[i]][sstem_j[i][ll_stem[i]]]=='T') {uus+=ATend;}
       }
  }
 }
/*********SBSt:experimental thermo-parameters for Single-base stacking*****************/
float ssBSt(int jc,int i1,int k1)
{
     float SB=0.0;  // 3'-i1-k1-5'     Thirumalai et al. J. Chem. Theory Comput., 2018, 14:3763-3779.
     if (type[jc][i1]=='A'&&type[jc][k1]=='A')     {SB=-3.53-T*0.001*(-10.96-SB0);}  // A-A
else if ((type[jc][i1]=='A'&&type[jc][k1]=='T')
      ||(type[jc][i1]=='T'&&type[jc][k1]=='A'))    {SB=-3.06-T*0.001*(-10.43-SB0);}  // A-T  (T-A)
else if ((type[jc][i1]=='A'&&type[jc][k1]=='G')
      ||(type[jc][i1]=='G'&&type[jc][k1]=='A'))    {SB=-3.76-T*0.001*(-11.26-SB0);}  // A-G  (G-A)
else if ((type[jc][i1]=='A'&&type[jc][k1]=='C')
      ||(type[jc][i1]=='C'&&type[jc][k1]=='A'))    {SB=-3.06-T*0.001*(-10.43-SB0);}  // A-C  (C-A)
else if (type[jc][i1]=='G'&&type[jc][k1]=='G')     {SB=-3.39-T*0.001*(-9.56-SB0);}   // G-G
else if ((type[jc][i1]=='G'&&type[jc][k1]=='C')
      ||(type[jc][i1]=='C'&&type[jc][k1]=='G'))    {SB=-4.28-T*0.001*(-12.9-SB0);}   // G-C  (C-G)
else if ((type[jc][i1]=='G'&&type[jc][k1]=='T')
      ||(type[jc][i1]=='T'&&type[jc][k1]=='G'))    {SB=-4.03-T*0.001*(-12.13-SB0);}  // G-T  (T-G)
else if (type[jc][i1]=='C'&&type[jc][k1]=='C')     {SB=-2.98-T*0.001*(-10.33-SB0);}  // C-C
else if ((type[jc][i1]=='C'&&type[jc][k1]=='T')
      ||(type[jc][i1]=='T'&&type[jc][k1]=='C'))    {SB=-2.98-T*0.001*(-10.33-SB0);}  // C-T  (T-C)
else if (type[jc][i1]=='T'&&type[jc][k1]=='T')     {SB=-2.98-T*0.001*(-10.33-SB0);}  // T-T
else    {SB=0.0;}   //Unknow Single-base stacking
     if (SB>0) {SB=0.0;}
     return SB;
}
float ssST(int jc,int i1,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size],int mark[chain_maxsize][max_size])
{
 float d,sst0,sst;//d0,d1,d2;
 //if ((mark[jc][i1]==1&&mark[jc][i1-3]==1)||i1<6) {sst=0.0;}
 //else {
    sst0=ssBSt(jc,i1,i1-3);
    if (sst0!=0) {
 d=sqrt((x1[jc][i1]-x1[jc][i1-3])*(x1[jc][i1]-x1[jc][i1-3])+(y1[jc][i1]-y1[jc][i1-3])*(y1[jc][i1]-y1[jc][i1-3])+(z1[jc][i1]-z1[jc][i1-3])*(z1[jc][i1]-z1[jc][i1-3])); //NN
/* d0=sqrt((x1[jc][i1]-x1[jc][j1-1])*(x1[jc][i1]-x1[jc][j1-1])+(y1[jc][i1]-y1[jc][j1-1])*(y1[jc][i1]-y1[jc][j1-1])+(z1[jc][i1]-z1[jc][j1-1])*(z1[jc][i1]-z1[jc][j1-1]));  //NiCj 
 d1=sqrt((x1[jc][j1]-x1[jc][i1-1])*(x1[jc][j1]-x1[jc][i1-1])+(y1[jc][j1]-y1[jc][i1-1])*(y1[jc][j1]-y1[jc][i1-1])+(z1[jc][j1]-z1[jc][i1-1])*(z1[jc][j1]-z1[jc][i1-1]));  //CiNj
 d2=sqrt((x1[jc][i1-2]-x1[jc][j1-1])*(x1[jc][i1-2]-x1[jc][j1-1])+(y1[jc][i1-2]-y1[jc][j1-1])*(y1[jc][i1-2]-y1[jc][j1-1])+(z1[jc][i1-2]-z1[jc][j1-1])*(z1[jc][i1-2]-z1[jc][j1-1]));  //NiCj 
       sst=sst0/(1+ksNN*(d-rsNN)*(d-rsNN)+ksCN*((d0-rsCN1)*(d0-rsCN1)+(d1-rsCN2)*(d1-rsCN2))+ksPC*(d2-rsPC)*(d2-rsPC));   */
      // sst=sst0/(1+ksNN*(d-rsNN)*(d-rsNN)+ks*(PCPC(i1-3,jc,x1,y1,z1)+NCPC(i1-3,jc,x1,y1,z1)+CPCP(i1,jc,x1,y1,z1)+CPCN(i1,jc,x1,y1,z1)));
       if (d<=4.0) {sst=0.0;}
       else {
       sst=sst0/(1+ksNN*(d-rsNN)*(d-rsNN)+ks*(5*PCPC(i1-3,jc,x1,y1,z1)+NCPC(i1-3,jc,x1,y1,z1)+5*CPCP(i1,jc,x1,y1,z1)+CPCN(i1,jc,x1,y1,z1)));}
      //if(t==1) printf("$$$ %d %d %f %f %f %f %f %f\n",t,i1,sst0,d,PCPC(i1-3,jc,x1,y1,z1),NCPC(i1-3,jc,x1,y1,z1),CPCP(i1,jc,x1,y1,z1),CPCN(i1,jc,x1,y1,z1));
    }
    else {sst=0.0;}
 //}
 return sst;
}
void SStacking(int i0,int jc,int mv)
{
 if (mv==1) {
    uss=ssST(jc,i0+2,x,y,z,c)+ssST(jc,i0+5,x,y,z,c)+ssST(jc,i0+8,x,y,z,c)+ssST(jc,i0+11,x,y,z,c);
    uuss=ssST(jc,i0+2,xx,yy,zz,a)+ssST(jc,i0+5,xx,yy,zz,a)+ssST(jc,i0+8,xx,yy,zz,a)+ssST(jc,i0+11,xx,yy,zz,a);
 }
 else if (mv==2) {
    if (fmod(i0+1,3)==0) {
       uss=ssST(jc,i0+1,x,y,z,c)+ssST(jc,i0+4,x,y,z,c);
       uuss=ssST(jc,i0+1,xx,yy,zz,a)+ssST(jc,i0+4,xx,yy,zz,a);}
    else if (fmod(i0+2,3)==0) {
       uss=ssST(jc,i0+2,x,y,z,c)+ssST(jc,i0-1,x,y,z,c);  //i0-1=3,how?
       uuss=ssST(jc,i0+2,xx,yy,zz,a)+ssST(jc,i0-1,xx,yy,zz,a);}
    else {
       if (i0==3) {uss=ssST(jc,i0+3,x,y,z,c);uuss=ssST(jc,i0+3,xx,yy,zz,a);}
       else if (i0==L0_chain[jc]-1) {uss=ssST(jc,i0,x,y,z,c);uuss=ssST(jc,i0,xx,yy,zz,a);}
       else {
          uss=ssST(jc,i0,x,y,z,c)+ssST(jc,i0+3,x,y,z,c);
          uuss=ssST(jc,i0,xx,yy,zz,a)+ssST(jc,i0+3,xx,yy,zz,a);}
    }
  }
  else {uss=0.0;uuss=0.0;}
}
/*********Potential of Exculded Volume between any two beads*******/
float LJ0(int jc,int jc1,int i1,int j1,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
  float d,d1,d6,d12,r1,r2,ULJ;
  r1=0.0;r2=0.0;
  d=sqrt((x1[jc][i1]-x1[jc1][j1])*(x1[jc][i1]-x1[jc1][j1])+(y1[jc][i1]-y1[jc1][j1])*(y1[jc][i1]-y1[jc1][j1])+(z1[jc][i1]-z1[jc1][j1])*(z1[jc][i1]-z1[jc1][j1])); 
  if (fmod(i1,3)==0&&fmod(j1,3)==0) {r1=2.15;r2=2.15;}  //if N-N, radii are changed to 2.15;
  else {r1=R[jc][i1];r2=R[jc1][j1];}  //radii: P 1.9; C 1.7; N 1.6;
  if (d<=(r1+r2))  {
      d1=(r1+r2)/(1.09*d); d6=d1*d1*d1*d1*d1*d1; d12=d6*d6;
      ULJ=(4.0*e*(d12-d6)+e);   //LJ(12-6) potential;
  }
  else ULJ=0.0;
  return ULJ;
} 
/*************Energy calculation of Exculded Volume for changed atoms*************/
 void ExcludedVolume(int i0,int i01,int i01_end)
{
  int i,j,jc;    //i01-i01_end：atoms of move;
  for(jc=1;jc<=Nchain;jc++) {      
      if(jc==chain_i0) {
         if(flag!=1) {
            for (i=i01;i<=i01_end;i++)   {
                 for (j=1;j<=L0_chain[jc];j++)  {        
                     if(j<i01||j>i01_end) { 
                          ulj+=LJ0(jc,jc,i,j,x,y,z);                             
                          uulj+=LJ0(jc,jc,i,j,xx,yy,zz);
                     } 
                 }   
            }
         }
         else {ulj=ulj; uulj=uulj;}        
       }
       else {
          for (i=i01;i<=i01_end;i++)   {
              for (j=1;j<=L0_chain[jc];j++)  {        
                   ulj+=LJ0(chain_i0,jc,i,j,x,y,z);                                
                   uulj+=LJ0(chain_i0,jc,i,j,xx,yy,zz);    
              } 
            }         
       } 
  }
}
// ~~~~~~~Calculation of elsectrostatic between two bead: Debye combining with CC theroy~~~~~~~~~~~~
float QQ(int jc,int jc1,int i1,int j1,float x1[chain_maxsize][max_size],float y1[chain_maxsize][max_size],float z1[chain_maxsize][max_size])
{
  float d,UQQ=0.0;
  d=sqrt((x1[jc][i1]-x1[jc1][j1])*(x1[jc][i1]-x1[jc1][j1])+(y1[jc][i1]-y1[jc1][j1])*(y1[jc][i1]-y1[jc1][j1])+(z1[jc][i1]-z1[jc1][j1])*(z1[jc][i1]-z1[jc1][j1])); 
  UQQ=q4*q4*(330.9*Q[jc][i1]*Q[jc1][j1]/(Ek*d))*exp(-d/Kd);
  return UQQ;
}
 void Electrostatic(int i0,int i01,int i01_end)
{
  int i,j,i02,i02_end,jc;
  if (salt==0||fmod(i0,3)==0) {uc=0.0;uuc=0.0;}
  else  {
      i02=4; i02_end=L0_chain[chain_i0]-3;
      for(jc=1;jc<=Nchain;jc++) {      
          if(jc==chain_i0) {
             if(flag!=1) {
                for(i=i01;i<=i01_end;i++) {
                    if(fmod(i+2,3)==0)   { 
                       for (j=i02;j<=i02_end;j++)  {        
                             if((j<i01||j>i01_end)&&fmod(j+2,3)==0) {  
                                uc+=QQ(chain_i0,jc,i,j,x,y,z);    uuc+=QQ(chain_i0,jc,i,j,xx,yy,zz); 
                             } 
                       } 
                    }    
                }
              }
           }
           else {
              for(i=i01;i<=i01_end;i++) {
                  if(fmod(i+2,3)==0)   { 
                     for (j=4;j<=L0_chain[jc];j++)  {        
                          if(fmod(j+2,3)==0) {  
                             uc+=QQ(chain_i0,jc,i,j,x,y,z);  uuc+=QQ(chain_i0,jc,i,j,xx,yy,zz); 
                          } 
                     } 
                  }    
              }             
           }
      }
  }
 }
// Bonded potential before and after conformational change: bond, angle and dihedral
 void Bonded_Move_3nt(int i0,int jc)
{
     if (fmod(i0+2,3)==0)  {
        ub= PC(i0+2,jc,x,y,z)+CP(i0+8,jc,x,y,z); uub=PC(i0+2,jc,xx,yy,zz)+CP(i0+8,jc,xx,yy,zz);
        ue= CPC(i0-1,jc,x,y,z)+PCP(i0+2,jc,x,y,z)+PCN(i0+2,jc,x,y,z)+PCP(i0+8,jc,x,y,z)+CPC(i0+8,jc,x,y,z)+NCP(i0+8,jc,x,y,z);
        uue=CPC(i0-1,jc,xx,yy,zz)+PCP(i0+2,jc,xx,yy,zz)+PCN(i0+2,jc,xx,yy,zz)+PCP(i0+8,jc,xx,yy,zz)+CPC(i0+8,jc,xx,yy,zz)+NCP(i0+8,jc,xx,yy,zz);
        ud= PCPC(i0-1,jc,x,y,z)+CPCP(i0-1,jc,x,y,z)+PCPC(i0+2,jc,x,y,z)+NCPC(i0-1,jc,x,y,z)+CPCN(i0+2,jc,x,y,z)+CPCP(i0+5,jc,x,y,z)+PCPC(i0+8,jc,x,y,z)+CPCP(i0+8,jc,x,y,z)+NCPC(i0+8,jc,x,y,z)+CPCN(i0+11,jc,x,y,z);
        uud=PCPC(i0-1,jc,xx,yy,zz)+CPCP(i0-1,jc,xx,yy,zz)+PCPC(i0+2,jc,xx,yy,zz)+NCPC(i0-1,jc,xx,yy,zz)+CPCN(i0+2,jc,xx,yy,zz)+CPCP(i0+5,jc,xx,yy,zz)+PCPC(i0+8,jc,xx,yy,zz)+CPCP(i0+8,jc,xx,yy,zz)+NCPC(i0+8,jc,xx,yy,zz)+CPCN(i0+11,jc,xx,yy,zz); 
     }
     else if (fmod(i0+1,3)==0)  {
        ub= CP(i0+1,jc,x,y,z)+CN(i0+1,jc,x,y,z)+PC(i0+10,jc,x,y,z); uub=CP(i0+1,jc,xx,yy,zz)+CN(i0+1,jc,xx,yy,zz)+PC(i0+10,jc,xx,yy,zz);
        ue= PCP(i0+1,jc,x,y,z)+CPC(i0+1,jc,x,y,z)+PCN(i0+1,jc,x,y,z)+NCP(i0+1,jc,x,y,z)+CPC(i0+7,jc,x,y,z)+PCP(i0+10,jc,x,y,z)+PCN(i0+10,jc,x,y,z);
        uue=PCP(i0+1,jc,xx,yy,zz)+CPC(i0+1,jc,xx,yy,zz)+PCN(i0+1,jc,xx,yy,zz)+NCP(i0+1,jc,xx,yy,zz)+CPC(i0+7,jc,xx,yy,zz)+PCP(i0+10,jc,xx,yy,zz)+PCN(i0+10,jc,xx,yy,zz);
        ud= CPCP(i0-2,jc,x,y,z)+PCPC(i0+1,jc,x,y,z)+CPCP(i0+1,jc,x,y,z)+NCPC(i0+1,jc,x,y,z)+CPCN(i0+1,jc,x,y,z)+CPCN(i0+4,jc,x,y,z)+CPCP(i0+7,jc,x,y,z)+PCPC(i0+7,jc,x,y,z)+PCPC(i0+10,jc,x,y,z)+NCPC(i0+7,jc,x,y,z)+CPCN(i0+10,jc,x,y,z);
        uud=CPCP(i0-2,jc,xx,yy,zz)+PCPC(i0+1,jc,xx,yy,zz)+CPCP(i0+1,jc,xx,yy,zz)+NCPC(i0+1,jc,xx,yy,zz)+CPCN(i0+1,jc,xx,yy,zz)+CPCN(i0+4,jc,xx,yy,zz)+CPCP(i0+7,jc,xx,yy,zz)+PCPC(i0+7,jc,xx,yy,zz)+PCPC(i0+10,jc,xx,yy,zz)+NCPC(i0+7,jc,xx,yy,zz)+CPCN(i0+10,jc,xx,yy,zz); 
     }
     else {ub=0.0; ue=0.0; ud=0.0; uub=0.0; uue=0.0; uud=0.0;}
}
 void Bonded_Pivot_angle(int i0,int jc)
{
     if (fmod(i0+2,3)==0)  {
         ud= CPCP(i0-1,jc,x,y,z)+PCPC(i0+2,jc,x,y,z)+CPCN(i0+2,jc,x,y,z);
         uud=CPCP(i0-1,jc,xx,yy,zz)+PCPC(i0+2,jc,xx,yy,zz)+CPCN(i0+2,jc,xx,yy,zz); 
     }
     else if (fmod(i0+1,3)==0)  {
         ue= PCN(i0+1,jc,x,y,z); uue=PCN(i0+1,jc,xx,yy,zz);
         ud= PCPC(i0+1,jc,x,y,z)+CPCP(i0+1,jc,x,y,z)+NCPC(i0+1,jc,x,y,z)+CPCN(i0+1,jc,x,y,z)+CPCN(i0+4,jc,x,y,z);
         uud=PCPC(i0+1,jc,xx,yy,zz)+CPCP(i0+1,jc,xx,yy,zz)+NCPC(i0+1,jc,xx,yy,zz)+CPCN(i0+1,jc,xx,yy,zz)+CPCN(i0+4,jc,xx,yy,zz); 
     }
     else {ue=0.0; uue=0.0; ud=0.0; uud=0.0;}          
}
 void Bonded_Rot_3nt(int i0,int jc)
{
     ue=CPC(i0-1,jc,x,y,z)+CPC(i0+8,jc,x,y,z); uue=CPC(i0-1,jc,xx,yy,zz)+CPC(i0+8,jc,xx,yy,zz);
     ud=PCPC(i0-1,jc,x,y,z)+CPCP(i0-1,jc,x,y,z)+NCPC(i0-1,jc,x,y,z)+CPCN(i0+2,jc,x,y,z)+PCPC(i0+8,jc,x,y,z)+CPCP(i0+8,jc,x,y,z)+NCPC(i0+8,jc,x,y,z)+CPCN(i0+11,jc,x,y,z);
     uud=PCPC(i0-1,jc,xx,yy,zz)+CPCP(i0-1,jc,xx,yy,zz)+NCPC(i0-1,jc,xx,yy,zz)+CPCN(i0+2,jc,xx,yy,zz)+PCPC(i0+8,jc,xx,yy,zz)+CPCP(i0+8,jc,xx,yy,zz)+NCPC(i0+8,jc,xx,yy,zz)+CPCN(i0+11,jc,xx,yy,zz);
}
 void Bonded_Pivot(int i0,int jc)
{
     if (fmod((i0+1),3)==0)  {       // Choose C atoms  
         ub=PC(i0+1,jc,x,y,z);           uub=PC(i0+1,jc,xx,yy,zz);
         if (i0==2)  {
             ue=PCN(i0+1,jc,x,y,z)+PCP(i0+1,jc,x,y,z);  
             uue=PCN(i0+1,jc,xx,yy,zz)+PCP(i0+1,jc,xx,yy,zz);             
             ud=PCPC(i0+1,jc,x,y,z)+CPCP(i0+1,jc,x,y,z)+NCPC(i0+1,jc,x,y,z)+CPCN(i0+4,jc,x,y,z);
             uud=PCPC(i0+1,jc,xx,yy,zz)+CPCP(i0+1,jc,xx,yy,zz)+NCPC(i0+1,jc,xx,yy,zz)+CPCN(i0+4,jc,xx,yy,zz); } 
         else if (i0==(L0_chain[chain_i0]-2)) {
             ue=PCN(i0+1,jc,x,y,z)+PCP(i0+1,jc,x,y,z)+CPC(i0-2,jc,x,y,z);                               
             uue=PCN(i0+1,jc,xx,yy,zz)+PCP(i0+1,jc,xx,yy,zz)+CPC(i0-2,jc,xx,yy,zz);  
             ud=PCPC(i0-2,jc,x,y,z)+CPCP(i0-2,jc,x,y,z)+NCPC(i0-2,jc,x,y,z)+CPCN(i0+1,jc,x,y,z);
             uud=PCPC(i0-2,jc,xx,yy,zz)+CPCP(i0-2,jc,xx,yy,zz)+NCPC(i0-2,jc,xx,yy,zz)+CPCN(i0+1,jc,xx,yy,zz); }   
         else  {
             ue=PCN(i0+1,jc,x,y,z)+PCP(i0+1,jc,x,y,z)+CPC(i0-2,jc,x,y,z);                               
             uue=PCN(i0+1,jc,xx,yy,zz)+PCP(i0+1,jc,xx,yy,zz)+CPC(i0-2,jc,xx,yy,zz);  
             ud=PCPC(i0-2,jc,x,y,z)+PCPC(i0+1,jc,x,y,z)+CPCP(i0-2,jc,x,y,z)+NCPC(i0-2,jc,x,y,z)+CPCN(i0+1,jc,x,y,z);
             uud=PCPC(i0-2,jc,xx,yy,zz)+PCPC(i0+1,jc,xx,yy,zz)+CPCP(i0-2,jc,xx,yy,zz)+NCPC(i0-2,jc,xx,yy,zz)+CPCN(i0+1,jc,xx,yy,zz); }   
      }
      else if (fmod((i0+2),3)==0)  {      // Choose P atoms 
         if (i0!=1) {ub=CP(i0-1,jc,x,y,z);uub=CP(i0-1,jc,xx,yy,zz);}
         if (i0==1) {
             ue=PCN(i0+2,jc,x,y,z)+PCP(i0+2,jc,x,y,z);    uue=PCN(i0+2,jc,xx,yy,zz)+PCP(i0+2,jc,xx,yy,zz); 
             ud=PCPC(i0+2,jc,x,y,z);                   uud=PCPC(i0+2,jc,xx,yy,zz); }
         else if (i0==L0_chain[chain_i0])  {
             ue=NCP(i0-1,jc,x,y,z)+PCP(i0-1,jc,x,y,z); uue=NCP(i0-1,jc,xx,yy,zz)+PCP(i0-1,jc,xx,yy,zz);
             ud=CPCP(i0-4,jc,x,y,z);                uud=CPCP(i0-4,jc,xx,yy,zz); } 
         else  {
             ue=NCP(i0-1,jc,x,y,z)+PCP(i0-1,jc,x,y,z)+CPC(i0-1,jc,x,y,z);              
             uue=NCP(i0-1,jc,xx,yy,zz)+PCP(i0-1,jc,xx,yy,zz)+CPC(i0-1,jc,xx,yy,zz);
             ud=CPCP(i0-4,jc,x,y,z)+PCPC(i0-1,jc,x,y,z)+CPCP(i0-1,jc,x,y,z)+NCPC(i0-1,jc,x,y,z)+CPCN(i0+2,jc,x,y,z);
             uud=CPCP(i0-4,jc,xx,yy,zz)+PCPC(i0-1,jc,xx,yy,zz)+CPCP(i0-1,jc,xx,yy,zz)+NCPC(i0-1,jc,xx,yy,zz)+CPCN(i0+2,jc,xx,yy,zz); }   
      }
      else     {           // Choose N atoms
         ub=CN(i0,jc,x,y,z);                 uub=CN(i0,jc,xx,yy,zz);
         ue=NCP(i0,jc,x,y,z)+PCN(i0,jc,x,y,z);  uue=NCP(i0,jc,xx,yy,zz)+PCN(i0,jc,xx,yy,zz); 
         if (i0==3) {ud=NCPC(i0,jc,x,y,z); uud=NCPC(i0,jc,xx,yy,zz);}
         else if (i0==L0_chain[chain_i0]-1) {ud=CPCN(i0,jc,x,y,z); uud=CPCN(i0,jc,xx,yy,zz);}
         else       {ud=NCPC(i0,jc,x,y,z)+CPCN(i0,jc,x,y,z); uud=NCPC(i0,jc,xx,yy,zz)+CPCN(i0,jc,xx,yy,zz);}
      }
}

 void CoaxialStacking()
{
   int i,j,LLij1,LLij2,LLij3,LLij4;
   for(i=1;i<=NN_stem;i++) {
       if (ll_stem[i]>=3) {   //coaxial-stacking
           for (j=i+1;j<=NN_stem;j++)  {
                if (ll_stem[j]>=3)     {
                    LLij1=(sstem_i[j][1]-sstem_i[i][ll_stem[i]])/3; LLij2=(sstem_j[j][1]-sstem_j[i][ll_stem[i]])/3;
                    LLij3=(sstem_i[j][1]-sstem_j[i][1])/3;   LLij4=(sstem_j[i][ll_stem[i]]-sstem_i[j][ll_stem[j]])/3;
                    if (sstem_i_c[i]==sstem_i_c[j]&&sstem_j_c[i]==sstem_j_c[j]) {
                        if (LLij1<5||LLij2<5)  {
                            if (LLij1==1||LLij2==1) {  //Bulge
                                uuco+=UCoS(sstem_i_c[i],sstem_j_c[i],sstem_i_c[j],sstem_j_c[j],sstem_i[i][ll_stem[i]],sstem_j[i][ll_stem[i]],sstem_i[j][1],sstem_j[j][1],xx,yy,zz,dco,dco,kdco,kdco);
                            }
                            else if (LLij2>=5) {  //Internal or Juntion
                                uuco+=UCoS(sstem_i_c[i],sstem_j_c[i],sstem_i_c[j],sstem_j_c[j],sstem_i[i][ll_stem[i]],sstem_j[i][ll_stem[i]],sstem_i[j][1],sstem_j[j][1],xx,yy,zz,LLij1*dco,LLij1*dco,kdco,kdco); 
                            }
                            else if (LLij1>=5) {  //Internal or Juntion
                                uuco+=UCoS(sstem_i_c[i],sstem_j_c[i],sstem_i_c[j],sstem_j_c[j],sstem_i[i][ll_stem[i]],sstem_j[i][ll_stem[i]],sstem_i[j][1],sstem_j[j][1],xx,yy,zz,LLij2*dco,LLij2*dco,kdco,kdco); 
                            }  
                            else  {
                                uuco+=UCoS(sstem_i_c[i],sstem_j_c[i],sstem_i_c[j],sstem_j_c[j],sstem_i[i][ll_stem[i]],sstem_j[i][ll_stem[i]],sstem_i[j][1],sstem_j[j][1],xx,yy,zz,LLij1*dco,LLij2*dco,kdco,kdco); 
                            }               
                        }
                        else if (LLij3<5) {  //Juntion
                                uuco+=UCoS(sstem_i_c[i],sstem_j_c[i],sstem_i_c[j],sstem_j_c[j],sstem_i[i][1],sstem_j[i][1],sstem_j[j][1],sstem_i[j][1],xx,yy,zz,LLij3*dco,LLij3*dco,kdco,kdco);                            
                        }
                        else if (LLij4<3) {  //Pseudoknot: should be further improved!!!
                            uuco+=UCoS(sstem_i_c[i],sstem_j_c[i],sstem_i_c[j],sstem_j_c[j],sstem_i[i][ll_stem[i]],sstem_j[i][ll_stem[i]],sstem_j[j][ll_stem[j]],sstem_i[j][ll_stem[j]],xx,yy,zz,LLij4*dco,LLij4*dco,kdco,kdco);
                        }  
                        else {uuco=uuco;}
                    }
                    else if (sstem_i_c[i]==sstem_i_c[j]&&LLij1<5) {
                        uuco+=UCoS(sstem_i_c[i],sstem_j_c[i],sstem_i_c[j],sstem_j_c[j],sstem_i[i][ll_stem[i]],sstem_j[i][ll_stem[i]],sstem_i[j][1],sstem_j[j][1],xx,yy,zz,LLij1*dco,LLij1*dco,kdco,kdco);                     
                    }
                    else if (sstem_j_c[i]==sstem_j_c[j]) {
                        uuco+=UCoS(sstem_i_c[i],sstem_j_c[i],sstem_i_c[j],sstem_j_c[j],sstem_i[i][ll_stem[i]],sstem_j[i][ll_stem[i]],sstem_i[j][1],sstem_j[j][1],xx,yy,zz,LLij2*dco,LLij2*dco,kdco,kdco);                     
                    }
                    else {uuco=uuco;}
                }
           }
       }
  }
}
void SecondaryStruct()
{
  int i,j,jc,jc1,k1=0;
   for (i=0;i<Nstem_maxsize;i++)     {
       ll_stem[i]=1; sstem_i_c[i]=0; sstem_j_c[i]=0;                      //l_stem[i]:length of i-th stem
       for (j=0;j<Lstem_maxsize;j++) {sstem_i[i][j]=0; sstem_j[i][j]=0; } //The j-th base-pair (stem_i-stem_j) in stem i;
   } 
   for(jc=1;jc<=Nchain;jc++)  {
      for(i=1;i<=L0_chain[jc];i++) {
          if(fmod(i,3)==0) {
             for(jc1=1;jc1<=Nchain;jc1++) {
                 for(j=1;j<=L0_chain[jc1];j++)  {
                     if(((jc==jc1&&j>=i+12)||(jc<jc1))&&fmod(j,3)==0&&ss[jc][jc1][i][j]==1) { 
                          if (ss[jc][jc1][i-3][j+3]==1)  {ll_stem[k1]++;}
                          else  {k1++;}
                          sstem_i[k1][ll_stem[k1]]=i;  sstem_j[k1][ll_stem[k1]]=j;  //nucleotides in stem[k1], stem_i paired with stem_j 
                          sstem_i_c[k1]=jc; sstem_j_c[k1]=jc1;   //chian's No. of the stem[k1]                           
                          break;
                     }
                 }
              }
          }
      }
   }
   NN_stem=k1;
}
// Energy of one confromation
 void ENERGY()
 {
     int i,j,jc,jc1;
     float u0=0.0,uN_Temp=0.0;
     Energy=1;
     ulj=0.0;ub=0.0;uc=0.0; U=0.0; bp=0; BP=0; uss=0.0;
/***********************************/
     for(jc=1;jc<=Nchain;jc++)  {
         for(i=1;i<=L0_chain[jc];i++)  {
             for(j=i+1;j<=L0_chain[jc];j++) {
                 ulj+=LJ0(jc,jc,i,j,x,y,z);  
                 if(salt!=0&&fmod(i+2,3)==0&&fmod(j+2,3)==0)  {uc+=QQ(jc,jc,i,j,x,y,z);}  
                 if(j>=i+12&&fmod(i,3)==0&&fmod(j,3)==0&&s[jc][jc][i][j]==1)  {
                    bp++; uN1=HB(jc,jc,i,j,x,y,z); uN_Temp+=uN1;
                    if ((type[jc][i]=='C'||type[jc][j]=='C')&&uN1<=gamma*A)  {BP++;}
                    if ((type[jc][i]=='T'||type[jc][j]=='T')&&uN1<=gamma*Beta*A)  {BP++;}
                 }                         
             } 
             for(jc1=jc+1;jc1<=Nchain;jc1++)   {
                 for(j=1;j<=L0_chain[jc1];j++) {
                     ulj+=LJ0(jc,jc1,i,j,x,y,z);  
                     if(salt!=0&&fmod(i+2,3)==0&&fmod(j+2,3)==0)  {uc+=QQ(jc,jc1,i,j,x,y,z);}  
                     if(fmod(i,3)==0&&fmod(j,3)==0&&s[jc][jc1][i][j]==1)  {
                        bp++; uN1=HB(jc,jc1,i,j,x,y,z); uN_Temp+=uN1;
                        if ((type[jc][i]=='C'||type[jc1][j]=='C')&&uN1<=gamma*A)  {BP++;}
                        if ((type[jc][i]=='T'||type[jc1][j]=='T')&&uN1<=gamma*Beta*A)  {BP++;}
                     }                      
                }

             }
             u0=0.0;
             if(fmod(i,3)==0)  {
                if(i==3)   {u0=PC(i,jc,x,y,z)+CP(i,jc,x,y,z)+CN(i,jc,x,y,z)+PCP(i,jc,x,y,z)+CPC(i,jc,x,y,z)+PCN(i,jc,x,y,z)+NCP(i,jc,x,y,z)+PCPC(i,jc,x,y,z)+CPCP(i,jc,x,y,z)+NCPC(i,jc,x,y,z);}
                else if(i==(L0_chain[jc]-1))  {u0=PC(i,jc,x,y,z)+CP(i,jc,x,y,z)+CN(i,jc,x,y,z)+PCP(i,jc,x,y,z)+PCN(i,jc,x,y,z)+NCP(i,jc,x,y,z)+CPCN(i,jc,x,y,z);} 
                else {u0=PC(i,jc,x,y,z)+CP(i,jc,x,y,z)+CN(i,jc,x,y,z)+PCP(i,jc,x,y,z)+CPC(i,jc,x,y,z)+PCN(i,jc,x,y,z)+NCP(i,jc,x,y,z)+PCPC(i,jc,x,y,z)+CPCP(i,jc,x,y,z)+CPCN(i,jc,x,y,z)+NCPC(i,jc,x,y,z);}
                ub+=u0*0.5963; 
                if(i>3)  {uss+=ssST(jc,i,x,y,z,c);}   
             }     
         }                           
    }
    if (uN_Temp!=uN) {printf("Error: Energy of bp is wrong!!\n");}  
    if (salt==0) {uc=0.0;uuc=0.0;}

    U=ulj+ub+uN+us+uc+uco+kss*uss; l++;  
    if (U<=Umin) {Umin=U;}
    else {Umin=Umin;}
    Up0+=U; Up=Up0/(t/tenergy);
    fprintf(fp6,"%d %d %f %f %f %f %f %f %f %f %f %f\n",t,l,U,Up,Umin,ulj,ub,uN,us,uc,uco,kss*uss);  fflush(fp6);
   //Statistics of the number of base pairs     
    m+=bp; pbp=1.0*m/(t/tenergy+1); G=1-pbp/Nbp; 
    M+=BP; PBP=1.0*M/(t/tenergy+1); GG=1-PBP/Nbp;
    fprintf(fp3,"%d %d %d %f %f %d %f %f\n",t,l,bp,pbp,G,BP,PBP,GG);  fflush(fp3);
    for (i=1;i<=N_stem;i++) {
         for (j=1;j<=l_stem[i];j++) {
              fprintf(fp5,"%d stem_%d %d %d %d %d %c %d %c\n",t,i,j,stem_i_c[i],stem_j_c[i],stem_i[i][j]/3,type[stem_i_c[i]][stem_i[i][j]],stem_j[i][j]/3,type[stem_j_c[i]][stem_j[i][j]]); 
              fflush(fp5);
         } 
    }
         fprintf(fp5,"\n");   fflush(fp5);
 }
// Structure Information, such as end-to-end distance, raduis of gyration and persistence length
 void StructureInformation()
{
    int   i,tk,j=0,jc;
    float dis=0.0,dis0=0.0,Xc=0.0,Yc=0.0,Zc=0.0,Rg,end,lp1,b_proj=0.0,lp2,lp3;
       for(jc=1;jc<=Nchain;jc++)  {
         for(i=1;i<=L0_chain[jc];i++)  {             
             if (fmod(i+1,3)==0) {
                Xc+=x[jc][i]; Yc+=y[jc][i]; Zc+=z[jc][i];               
                if (i>4) {
                   dis+=sqrt((x[jc][i]-x[jc][i-3])*(x[jc][i]-x[jc][i-3])+(y[jc][i]-y[jc][i-3])*(y[jc][i]-y[jc][i-3])+(z[jc][i]-z[jc][i-3])*(z[jc][i]-z[jc][i-3]));
                   b_proj+=(x[jc][5]-x[jc][2])*(x[jc][i]-x[jc][i-3])+(y[jc][5]-y[jc][2])*(y[jc][i]-y[jc][i-3])+(z[jc][5]-z[jc][2])*(z[jc][i]-z[jc][i-3]);
                }
                j++;
             }
         }
         lp3+=(x[jc][5]-x[jc][2])*(x[jc][L0_chain[jc]-1]-x[jc][2])+(y[jc][5]-y[jc][2])*(y[jc][L0_chain[jc]-1]-y[jc][2])+(z[jc][5]-z[jc][2])*(z[jc][L0_chain[jc]-1]-z[jc][2]);
       }
       Xc=Xc/(j*1.); Yc=Yc/(j*1.); Zc=Zc/(j*1.);  //Center of structure 
       for(jc=1;jc<=Nchain;jc++)  {
         for(i=1;i<=L0_chain[jc];i++)  {
             if (fmod(i+1,3)==0) {
                dis0+=(x[jc][i]-Xc)*(x[jc][i]-Xc)+(y[jc][i]-Yc)*(y[jc][i]-Yc)+(z[jc][i]-Zc)*(z[jc][i]-Zc);}
         }
       }
       Rg=sqrt(dis0/(j*1.));
       end=sqrt((x[1][1]-x[Nchain][L0_chain[Nchain]])*(x[1][1]-x[Nchain][L0_chain[Nchain]])+(y[1][1]-y[Nchain][L0_chain[Nchain]])*(y[1][1]-y[Nchain][L0_chain[Nchain]])+(z[1][1]-z[Nchain][L0_chain[Nchain]])*(z[1][1]-z[Nchain][L0_chain[Nchain]]));
       lp1=end*end/(2*dis);
       float rb=dis/((j-1)*1.0); //mean of bond Ci-Ci
       lp2=b_proj/rb;
       lp3=lp3/rb;
       if(t==1) {rgp=Rg;rep=end;rlp1=lp1;rlp2=lp2;}      
       tk=t/t_si;
       rgp=(rgp*(tk-1)+Rg)/(tk);   
       rep=(rep*(tk-1)+end)/(tk); 
       rlp1=(rlp1*(tk-1)+lp1)/(tk); 
       rlp2=(rlp2*(tk-1)+lp2)/(tk);
       rlp3=(rlp3*(tk-1)+lp3)/(tk);  
       fprintf(fp2,"%d %d %f %f %f %f %f %f %f %f %f %f %d\n",t,tk,Rg,end,lp1,lp2,lp3,rgp,rep,rlp1,rlp2,rlp3,nm); fflush(fp2);
 }

//Out put some structural parameters and input value
void OutputPara()
{
 time_t timep; time(&timep);fprintf(fp7,"This program starts at :%s\n",ctime(&timep)); //Output time;
 fprintf(fp7,"*************************\nSequence information: \nChain number: %d\n",Nchain);
 for (int i=1;i<=Nchain;i++) {fprintf(fp7,"The %d-th chain: Length = %d & Sequence = %s\n",i,Lchain[i],sequence[i]);}
 fprintf(fp7,"*************************\nCondition:\nInitial Temperature = %f  &  Steps = %d  & Annealing way = %d\n",t0,total,Anneal);
 fprintf(fp7,"Folding = %d & salt= %d\n",Folding,salt);
 fprintf(fp7,"Output parameter:  tconf = %d  &  tenergy = %d  &  t_si = %d  &  tprintf = %d\n",tconf,tenergy,t_si,tprint);
 fprintf(fp7,"Ion condition:     CNa = %f  &  CMg = %f\n***************\n",CNa,CMg);
 fprintf(fp7,"*************************\nMove sets: p_Move_3nt = %f (step = %f) & p_move_pivot_angle = %f (step = %f)  & step_pivot = %f\n",p_move_3nt,step_3nt,p_move_pivot_angle,step_angle,step_pivot);
 fprintf(fp7,"Energy parameters: base-pairing = %f  &  base-stacking = %f  & A(G)U(T)/GC = %f/\n",A,B0,Beta);
 fflush(fp7); 
}
void FOLD(int i0)
{
 int   i,j,i0_1,i0_end,move_way=0; 
 float mv_x,mv_y,mv_z;
 float GaussianNum();
 for(j=1;j<=Nchain;j++)  {
    if (j!=chain_i0) { //keep the other chain fixed,only the chain_i0-th chain move;
        for (i=1;i<=L0_chain[j];i++) {xx[j][i]=x[j][i];yy[j][i]=y[j][i];zz[j][i]=z[j][i];} 
        move_way=0;
    }
    else { 
        if (rand()/(RAND_MAX+1.)<p_move_3nt&&i0>=4&&i0<L0_chain[chain_i0]-12&&fmod(i0+2,3)==0) {    //3nt fragment translation;
            move_way=1;
            step=step_3nt; 
            i0_1=i0+1; i0_end=i0+8;      
            rand01=rand()/(RAND_MAX+1.);
            for(i=1;i<=L0_chain[j];i++)  {      //i0+1-->i0+8:translated
                if (i>i0&&i-i0<=8) {Translate(j,i);xx[j][i]=xxx[i];yy[j][i]=yyy[i];zz[j][i]=zzz[i];}
                else  {xx[j][i]=x[j][i];yy[j][i]=y[j][i];zz[j][i]=z[j][i];}  }  
            Bonded_Move_3nt(i0,j);              
        }      
        else  {      //Pivot_moves  
            move_way=2;
            if (i0!=1) {step=step_pivot;} else {step=20*step_pivot;} //else {step=-19*step_pivot*t/(total-1)+20*step_pivot+19*step_pivot/(total-1);}
            if (fmod(i0,3)!=0)    {           //Backbone                
                i0_1=i0; i0_end=L0_chain[j]; 
                Rand01();
                for (i=1;i<=L0_chain[j];i++)   {
                     if (i<i0) {xx[j][i]=x[j][i];yy[j][i]=y[j][i];zz[j][i]=z[j][i];}
                     else      {
                         Translate(j,i);                      
                         if (i==i0) {xx[j][i]=xxx[i];yy[j][i]=yyy[i];zz[j][i]=zzz[i];} }
                }  
                if (i0!=L0_chain[j])   {
                   for (i=i0+1;i<=L0_chain[j];i++)  Pivot(j,i,i0);}   
            }
            else  {                     //N(side-chain) is chosen                
                i0_1=i0; i0_end=i0;
                for(i=1;i<=L0_chain[j];i++)  {
                    if (i!=i0) {xx[j][i]=x[j][i];yy[j][i]=y[j][i];zz[j][i]=z[j][i];}
                    else  {MoveN(j,i0);}
                }            
            }
            Bonded_Pivot(i0,j); 
        }  
        flag=0;
        if(fmod(i0,3)!=0) {
           for(i=i0_1;i<=i0_end;i++) {
               if((xx[j][i]<0||xx[j][i]>Lbox)||(yy[j][i]<0||yy[j][i]>Lbox)||(zz[j][i]<0||zz[j][i]>Lbox)) {
                  flag=1;             
                  mv_x=Lbox/2.-x[j][i0]; mv_y=Lbox/2.-y[j][i0]; mv_z=Lbox/2.-z[j][i0];
                  int k;
                  for(k=1;k<=L0_chain[j];k++) {
                      xx[j][k]=x[j][k]+mv_x; yy[j][k]=y[j][k]+mv_y; zz[j][k]=z[j][k]+mv_z;
                  }
               break;
               }
           }
        } 
        if (flag==1) {i0_1=1;i0_end=L0_chain[chain_i0];}           
    }
  }  
/*******Energy of conformations before and after moves;********/
        ExcludedVolume(i0,i0_1,i0_end);
        //Electrostatic(i0,i0_1,i0_end);      
        //BasePairing(); 
        //SecondaryStruct();
        //BaseStacking();  
        //CoaxialStacking();  
        SStacking(i0,chain_i0,move_way);

}
