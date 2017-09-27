#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <cblas.h>
#include <unistd.h>
#include <float.h>


/*_________________________________________________________*/
/**
 * @author Nicolas Alcala;  
 */

void shuffle(int *array, int n){
  if (n > 1){
    int i;
    for (i = 0; i < n-1; i++){
      int j = i + rand() / (RAND_MAX / (n - i) + 1);
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
}

//call: compPePc(Pe, Pc, D, piall, etmp, pCtmp, M, n, nstates, j);

void compPePc(double * Pee, double * Pcc, double * D,  int ** piall, double e, double * pC, double ** M, int n, long nstates, int j){
  //printf("(e=%lf, n=%d, nstates=%d, j=%d\n)",e,n,nstates,j);
  //compute probability transition from state piold to all states in piall
  int i = 0;
  double etmp = e;
  if(etmp>1) etmp = 1;
  int * pitmp = piall[j];
  for(i=0;i<nstates;i++){// for all possible states at time point k
    //printf("(comp state %d)\t",i);
    int * piold = piall[i];

    double pijc = 0; //1;
    int s1 = 0; //1 -> 0
    int s2 = 0; //1 -> 1
    int k;
    
    int ispos = 1;
    for(k=0; k<n; k++){
      if(pitmp[k]*(1-piold[k])>0){
	ispos = 0;
	break;//0 -> 1 during extinction phase: impossible!
      }
      s1 += (1-pitmp[k])*piold[k];
      s2 += pitmp[k]*piold[k];
      //pijc *= pitmp[k]*piold[k] + (1-pitmp[k])*( (1-piold[k])*(1-pC[k]) + piold[k]*pC[k] ); //note: piold is considered the new state here
      if(pitmp[k] ==0 ) pijc += log((1-piold[k])*(1-pC[k]) + piold[k]*pC[k] ); //note: piold is considered the new state here
    }
    if(ispos==1){
      //printf("Pee: %d; Pcc: %d\t",i*nstates + j,j*nstates + i);
      Pcc[ j*nstates + i]  = exp(pijc);
      //printf("Pcc done\n");
      Pee[ i*nstates + j] = pow(etmp,s1)*pow(1-etmp,s2);
      //printf("Pee done\n");
    }else{
      //printf("Pee: %d; Pcc: %d\t",i*nstates + j,j*nstates + i);
      Pcc[ j*nstates + i]  = 0;
      //printf("Pcc done\n");
      Pee[ i*nstates + j] = 0;
      //printf("Pee done\n");
    }//ispos
    
  }//i
}

double ylik(int ** Y, int ** J, int j,int n,double p, int * pi){
  double lres = 0;
  int i;
  if(p==0){
    for(i=0;i<n;i++){
      if(Y[j][i]>0) return 0;
    }
  }else{
    if(p==1){
      for(i=0;i<n;i++){
	if(Y[j][i]<J[j][i]) return 0;
	if( (Y[j][i]>0)&&(pi[i]==0)) return 0;
      }
    }else{
      for(i=0;i<n;i++){
    //printf("",pi);
    //if(Y[j][i]>= pi[i] ){
	if( Y[j][i]>0 ){
	  if(pi[i]==0){  return 0;} //printf("p=0,Y>0\n");
	  else lres += log(p)*Y[j][i] + log(1-p)*(J[j][i]-Y[j][i]);
	}else{
	  if(pi[i]!=0) lres += log(1-p)*J[j][i];
	}
	/*}else{
      return 0;
      }*/
      }
    }
  }
  //for(i=0;i<n;i++) printf("%d, %d/%d\t",pi[i],Y[j][i],J[j][i]);
  //printf("\t lres=%lf\n",lres);
  //printf("lres=%lf, e^lres=%lf\n",lres,exp(lres) );//"log(%lf)=%lf, log(1-%lf)=%lf",p,log(p),p,log(1-p));
  return exp(lres);
}

int pi_to_id(int* pi,int n){
  //converts state pi to identifier (in binary)
  int ip;
  int res=0;
  for(ip=0;ip<n;ip++) res += pi[ip]*pow(2,n-ip-1);
  return res;
}


int main(int argc, char ** argv)//takes the path to an input file as argument
{

  printf("sizeof size_t = %zx, SIZE_MAX = %zx\n", sizeof(size_t), SIZE_MAX);
  //initialize random number generator
  //main loop
  //records the time of execution
  time_t start,end;         //beginning and ending times
  double dif;               //difference between beginning and ending

  time (&start);//records the beginning time
  srand (time(NULL));
  
  printf("------ MIDASPOM, beta version ------\n-> N. Alcala, E. M. Cole, and N. A. Rosenberg <-\n");
    
  /*-------------- model and simulation Parameters--------------------*/
  float prioroc = 0.5; // prior occupancy for missing data the first year
  char* Mfname = NULL;
  char* fname = "input.txt";
  char* fout  = "posterior.txt";
  char* Jname = "J.txt";
  double d    = 100;
  double win  = 0.01;
  int nstepe = 101;
  int nstepc = 101;
  double a    = 1.0/400; //dispersal distance
  double emin = 0.0;
  double emax = 1.0;
  double cmin = 0.0;
  double cmax = 1.0;
  double p    = 1; //proba of observation given presence
  double C    = 100000;//scaling constant
  long nstates = 0; //pow(2,nvar); //number of possible states
  long maxnstates = 0; //pow(2,nvar); //number of possible states
  double thres = 0;
  
  int c;
  
  opterr = 0;
  while ((c = getopt (argc, argv, "m:M:p:d:i:J:o:s:l:h:u:v:f:t:n:C:")) != -1)
    switch (c)
      {
      case 'm':
        a = 1.0/atof(optarg);
        break;
      case 'M':
        Mfname = optarg;
        break;
      case 'p':
        prioroc = atof(optarg);
        break;
      case 'd':
        d = atof(optarg);
        break;
      case 'i':
        fname = optarg;
        break;
      case 'J':
        Jname = optarg;
        break;
      case 'o':
        fout = optarg;
        break;
      case 's':
        win = atof(optarg);
	break;
      case 'l':
        emin = atof(optarg);
        break;
      case 'h':
        cmin = atof(optarg);
        break;
      case 'u':
        emax = atof(optarg);
        break;
      case 'v':
        cmax = atof(optarg);
        break;
      case 't':
        thres = atof(optarg);
        break;
      case 'n':
        maxnstates = atoi(optarg);
        break;
      case 'C':
        C = atof(optarg);
        break;
      case 'f':
        p = atof(optarg);
        break;
      case '?':
        if (optopt == 'c')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }
  
  nstepe = (int) ((emax-emin)/win + 1);
  nstepc = (int) ((cmax-cmin)/win + 1);
  
  printf("Parameters for numerical approximation of the posterior density:\n\tWindow size=%lf, number of steps=%dx%d\n",win,nstepe,nstepc );
   
  //model parameters
  int i,j,k,l,m;

  //read observed data
  FILE * fo;
  printf("Reading observations from file %s... ",fname);
  int tmax=0;
  fo = fopen(fname,"rb");   //
  int n=1;
  while ( (c=fgetc(fo)) != EOF ) {
    if ( c == '\n' )
      tmax++;
    if(tmax==0){
      if( (c==' ')|(c=='\t') ) n++;
    }
  }
  fclose(fo);
  
  //int** piobs = malloc(tmax*sizeof(int*)); //piobs[tmax][n]; 
  //for(i=0; i<tmax; i++) piobs[i] = malloc(n*sizeof(int));

  int** Y     = malloc(tmax*sizeof(int*)); //piobs[tmax][n]; 
  for(i=0; i<tmax; i++) Y[i] = malloc(n*sizeof(int));
  
  fo = fopen(fname,"rb");   //
  int out;
  //int var[n];
  //for(j=0; j<n; j++) var[j]=0;
  for(i=0; i<tmax; i++){
    for(j=0; j<n; j++){
      out = fscanf(fo,"%d",&(Y[i][j]));
      /*if(Y[i][j]!=0){
	var[j]=1;
	//piobs[i][j] = 1;
      }else{
	var[j]=1;
	//piobs[i][j] = 0;
	}*/
    }
  }
  fclose(fo);
  printf("done\n");

  printf("Number of habitat patches: %d\nNumber of sampled years: %d\n",n,tmax);

  FILE * fJ;
  printf("Reading survey number from file %s... ",Jname);
  fJ = fopen(Jname,"rb");   //

  int** J = malloc(tmax*sizeof(int*)); //piobs[tmax][n]; 
  for(i=0; i<tmax; i++) J[i] = malloc(n*sizeof(int));
  
  for(i=0; i<tmax; i++){
    for(j=0; j<n; j++){
      out = fscanf(fJ,"%d",&(J[i][j]));
      if(J[i][j]==0){
	//piobs[i][j] = -1;
	//var[j] = 1;
      }
    }
  }
  fclose(fJ);
  printf("done\n");
  
  int nvar=n; //number of variable segments
  //for(j=0; j<n; j++) if(var[j]==1) nvar++;

  long nstatestot = (long) pow(2,n);
  if(maxnstates==0) maxnstates = nstatestot;
  printf("nstatestot=%ld\n",maxnstates); //ok
  
  double **M = malloc(n*sizeof(double*));
  for(i = 0; i < n; i++) M[i] = malloc(n*sizeof(double));

  if(Mfname==NULL){
    for(i=0; i<n; i++){
      for(j=i; j<n; j++){
	if(i==j) M[i][j]=0.0;
	else{
	  M[i][j]=exp(-a*(j-i)*d);
	  M[j][i]=exp(-a*(j-i)*d);
	}
      }
    }
  }else{
    printf("Reading migration rates from file %s... ",Jname);
    FILE * fM;
    fM = fopen(Mfname,"rb");   //
    for(i=0; i<n; i++){
      for(j=0; j<n; j++){
	out = fscanf(fM,"%lf",&(M[i][j]));
      }
    }
    fclose(fM);
    printf("done\n");
  }
  printf("Dispersal matrix:\n");
  for(i=0; i < n; i++){
    for(j=0; j<n; j++){
      printf("%.3f ",M[i][j]);
    }
    printf("\n");
  }
  
  //all possible states
  int maxnstatestmp = 0;
  int** piall  = malloc(10*tmax*maxnstates*sizeof(int*));
  int idpiall[10*tmax*maxnstates];
  for(k=0; k<10*tmax*maxnstates;k++){
    piall[k] = malloc(n*sizeof(int));
    idpiall[k] = 0;
  }
  double** prmin = malloc(tmax*sizeof(double*)); //list of maxnstates states with top proba /t value
  for(i=0;i<tmax;i++){
    prmin[i] = malloc(maxnstates*sizeof(double));
    for(j=0;j<maxnstates;j++) prmin[i][j] = -1;
  }
  double prtmp; //list of temporary proba /t value
  int sprtmp[tmax];   //list of current nb of val in prmin /t value
  for(i=0;i<tmax;i++) sprtmp[i] = 0;
  int pialltmp[n];
  nstates=0;
  for(i=0; i<nstatestot; i++){
    //printf("i=%d\t",i);
    int jtmp = 0;
    for(j=0; j<n;j++){
      int a = pow(2,n-jtmp-1);
      pialltmp[j] = i/a%2;
      jtmp++;
    }
    //sprtmp = 0;//ylik(Y,J,0,n,p,piall[i]);
    //for(k=0;k<n;k++) printf("%d ", pialltmp[k] );
    for(j=0; j<tmax; j++){
      //boolean added = false;
      prtmp = ylik(Y,J,j,n,p,pialltmp)*C; //tmp proba  int ** Y, int ** J, int j,int n,double p, int * pi
      //printf(", prtmp[%d]=%.20lf\n",j,prtmp);
      k=0;
      while(prtmp>=prmin[j][k]){
	if( (k>0)&&(sprtmp[j]==maxnstates) ){//if full 
	  prmin[j][k-1] = prmin[j][k]; //move all to the left to make room for new val
	  for(l=0; l<n;l++) piall[k-1 + j*maxnstates ][l] = piall[k + j*maxnstates ][l];
	  idpiall[k-1 + j*maxnstates ] = idpiall[k + j*maxnstates ];
	}
	if( ((sprtmp[j]<maxnstates)&&(k==sprtmp[j]))||((sprtmp[j]==maxnstates)&&(k==(maxnstates-1))) ){//if not full and reached val or full and reached end
	  //for(l=0; l<n;l++) piall[k + j*maxnstates ][l] = pialltmp[l];
	  break;
	}
	k++;
      }//while
      //printf("k=%d\n",k);
      if(sprtmp[j]==maxnstates){//full
	if( (k>0)||(prtmp>=prmin[j][0]) ){
	  if(prtmp<prmin[j][k]) k--;
	  prmin[j][k] = prtmp;
	  for(l=0; l<n;l++) piall[k + j*maxnstates ][l] = pialltmp[l];
	  idpiall[k + j*maxnstates ] = i;
	  //maxnstatestmp++;
	  //break;
	}
      }else{//not full: take everyone
	for(l=sprtmp[j];l>k;l--){//move everything to the right
	  prmin[j][l] = prmin[j][l-1];
	  for(m=0; m<n;m++) piall[l + j*maxnstates ][m] = piall[l-1 + j*maxnstates ][m];
	  idpiall[l + j*maxnstates ] = idpiall[l-1 + j*maxnstates ];
	}
	prmin[j][k] = prtmp;
	for(l=0; l<n;l++) piall[k + j*maxnstates ][l] = pialltmp[l];
	idpiall[k + j*maxnstates ] = i;
	//maxnstatestmp++;
	//break;
      }
      if(sprtmp[j]<maxnstates) sprtmp[j]++;
    }//j   
  }

  for(j=0; j<tmax; j++){
    for(l=0;l<maxnstates;l++){
      printf("pr[%d][%d]=%lf\t",j,l,log(prmin[j][l]) );
      for(k=0;k<n;k++) printf("%d ", piall[l+j*maxnstates][k] );
    }
    printf("\n");
  }
  printf("done\n");

  //find duplicate states
 
  /*int** pialle  = malloc(9*tmax*maxnstates*sizeof(int*)); //10 extinction/possible state
  for(k=0; k<maxnstatestmp; k++){
    pialle[k] = malloc(n*sizeof(int));
    }*/


  /*for(k=0; k<2*tmax*maxnstates;k++){
    printf("k=%d/%d\t",k,2*tmax*maxnstates);
    piall2[k] = malloc(n*sizeof(int));
    }*/

  printf("start\n");
  double pp;
  for(k=1;k<10;k++){
    printf("k=%d\t",k);
    for(i=0;i<tmax;i++){
      for(j=0;j<maxnstates;j++){
	idpiall[j + i*maxnstates + k*tmax*maxnstates] = 0;
	for(l=0;l<n;l++){
	  if(piall[j + i*maxnstates][l]==1){
	    pp = (double)rand()/(double)(RAND_MAX);
	    if(pp>(0.1*k)){
	      piall[j + i*maxnstates + k*tmax*maxnstates][l] = 1;
	      idpiall[j + i*maxnstates + k*tmax*maxnstates] += pow(2,n-l-1);
	    }else{
	      piall[j + i*maxnstates + k*tmax*maxnstates][l] = 0;
	    }
	  }else{
	    piall[j + i*maxnstates + k*tmax*maxnstates][l] = 0;
	  }
	}
      }//j
    }//i
  }//k
  
  int nstatestmp = 10*maxnstates*tmax;
  int idkeep[nstatestmp];
  for(j=0;j<nstatestmp;j++) idkeep[j] = 1;
  for(j=0;j<10*maxnstates*tmax-1;j++){
    for(k=j+1;k<10*maxnstates*tmax;k++){
      if(idkeep[j]){
	if(idpiall[j]==idpiall[k]){
	  idkeep[k] = 0;
	  nstatestmp--;
	}
      }
    }
  }

  
  /*for(k=0;k<10*maxnstates*tmax;k++){
    printf("k=%d, idkeep=%d, idpiall=%d\t",k,idkeep[k],idpiall[k]);
    for(l=0;l<n;l++) printf("%d\t",piall[k][l]);
    printf("\n");
    }*/

  int** piall2  = malloc(nstatestmp*sizeof(int*));
  for(k=0; k<nstatestmp;k++){
    piall2[k] = malloc(n*sizeof(int));
  }
  
  l=0;
  for(j=0;j<10*maxnstates*tmax;j++){
    if(idkeep[j]){
      piall2[l] = piall[j];
      l++;
    }
  }
  
  //for(j=0;j<tmax*maxnstates;j++) free(piall[j]);
  piall = piall2;
  
  //memcpy(&Pold[0], &Ptmp[0], sizeof(double)*nstates);
  
  for(j=0;j<nstatestmp;j++){
    for(k=0;k<n;k++) printf("%d ", piall[j][k] );
    printf("\n");
  }
  
   
  nstates=nstatestmp;
  printf("Number of states to compute: %d; total matrix size: %d\n",nstates,nstates*nstates);
  
  //number of possible states at each time step
  //int npstates[tmax]; 
  //int nnpstates = 0;
  //int** pstates = malloc(tmax*sizeof(int*)); //list of state id in all states list
  //int** simppstates = malloc(tmax*sizeof(int*)); //list of state id in states short list
  float* priorst;//prior probability of each initial state
  //int s1;
  //int st1;
  //int nextid = 0;
  //for(i=0; i<tmax; i++){
  //s1 = n;//0;
    /*for(j=0; j<n; j++){
      if(piobs[i][j]==-1){
	s1++; //number of segments with missing data
      }
      }*/
    //npstates[i] = pow(2,s1); //number of possible states
    //nnpstates += npstates[i]; //total number of possible states
    //if(i==0){
  priorst = malloc( nstates*sizeof(float));
  for(k=0; k<nstates; k++){
    priorst[k]=1;
    for(i=0;i<n;i++) priorst[k]*= piall[k][i]*prioroc + (1-piall[k][i])*(1-prioroc);
  }
     //print input
  printf("Input occupancy data:\n");
  for(i=0; i<tmax; i++){
    printf("Year %d: ",i);
    for(j=0; j<n; j++) printf("%d/%d\t",Y[i][j],J[i][j]); //,piobs[i][j]);
    printf("\n");
  }
  
  printf("Number of states to compute: %d\n",nstates);

  
  
  //create input tables
  double ellik[nstepe];
  double cllik[nstepc];
  for(i=0;i<nstepe-1;i++){
    ellik[i] = ((double)i)*win + emin;
  }
  for(i=0;i<nstepc-1;i++){
    cllik[i] = ((double)i)*win + cmin; 
  }
  ellik[nstepe-1] = emax;
  cllik[nstepc-1] = cmax;
  //create output tables
  double Lik[nstepe][nstepc];
  for(i=0;i<nstepe;i++){
    for(j=0;j<nstepc;j++){
      Lik[i][j] = 0;
    }
  }
  setbuf(stdout, NULL);
  
  int ie,ic;
  double etmp,ctmp;//Ltmp
  
  // start loop
  int stloop,endloop;
  stloop  = 0; 
  endloop = nstepe;
  double *P  = malloc(nstates*nstates*sizeof(double));
  double *D  = malloc(nstates*tmax*sizeof(double)); 
  double *pCtmp = malloc(n*sizeof(double)); //proba of colonization

  for(ie=stloop; ie<endloop; ie++){
    etmp = ellik[ie];
    for(ic=0; ic<nstepc; ic++){
      printf("%d %d\n",ie,ic);
      ctmp = cllik[ic];

      long long n2 = nstates*nstates;
      double * Pe = calloc(n2,sizeof(double));
      double * Pc = calloc(n2,sizeof(double));
      printf("Allocated %zu bytes from %p to %p\n", n2, &Pe, &Pe + n2);
      /*for(j=0;j<n2;j++){
	printf("%d:\t",j);
	Pe[j]=0;
	Pc[j]=0;
	printf("Pe[%d]=%lf, Pc[%d]=%lf\n",j,Pe[j],j,Pc[j]);
	}*/
      
      int snp;
      printf("Comp PePC states ");
      for(j=0;j<nstates;j++){
	//printf("%d\t",j);
	for(k=0; k<n; k++){//for all patches
	  double s1 = 0;
	  for(l=0; l<n; l++){
	    if(l!=k) s1 += M[l][k]*piall[j][l];
	  }//l
	  pCtmp[k] = ctmp*s1;
	  if(pCtmp[k]>1) pCtmp[k] = 1;
	}//k
	//printf("pCtmp done\t");
	//printf("Run compPePc, e=%lf, c=%lf, n=%d, nstates=%d, j=%d\n",etmp,ctmp,n,nstates,j);
	compPePc(Pe, Pc, D, piall, etmp, pCtmp, M, n, nstates, j);
	//printf("PePc done\n");
      }//j
      /*double SPe=0,SPc=0;
      for(j=0;j<nstates;j++){
	for(k=0;k<nstates;k++){
	  SPe += Pe[ j*nstates + k];
	  SPc += Pc[ k*nstates + j];
	}
      }
      printf("ie=%d,ic=%d: SPe=%lf, SPc=%lf\t",ie,ic,SPe,SPc);//OK
      */
      //printf("PePc totally done\n");
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nstates, nstates, nstates, 1.0, Pe, nstates, Pc, nstates, 0, P, nstates); //multiply Pe and Pc using cblas routine
      double SPec=0;
      //printf("\n\tP:");
      for(j=0;j<nstates;j++){
	for(k=0;k<nstates;k++){
	  SPec += P[ j*nstates + k];
	  //printf("%e\t",P[j*nstates + k]);//SPe OK, SPc weird
	  //if(P[j*nstates + k]==0) printf("error!");
	}
	//printf("\n");
      }
      printf("SPec=%lf\t",SPec);
      
      //printf("P done\n");
      free(Pe);
      free(Pc);

      double *Pold = malloc(nstates*sizeof(double));
      for(k=0;k<nstates;k++) Pold[k] = C*C*C*ylik(Y,J,0,n,p,piall[k])*priorst[k]; //Pold = Phi_0*D(p_1)
      double SP=0;
      for(l=0;l<nstates;l++){
	SP += Pold[ l];
      }
      printf("SPold0=%lf\t",log(SP));
      
      snp = 0;
      for(j=1; j<tmax; j++){//for each year
	//multiply by diagonal matrix indicating proba of observations
	double *Pcur = calloc(nstates*nstates,sizeof(double));   
	/*for(k=0;k<nstates;k++){//for each possible state at time j-1
	  for(l=0;l<nstates;l++){//for each possible state at time j
	    Pcur[k*nstates + l] = P[ simppstates[j-1][k]*nstates + simppstates[j][l] ]; //copy from precomputed list of transitions
	    //printf("%d: Pcur[%d]=%e\t",j, k*nstates + l,Pcur[k*nstates + l] ); //not always 1
	  }//l
	  }//k*/

	for(k=0; k<nstates; k++) Pcur[k*nstates + k] = C*ylik(Y,J,j,n,p,piall[k]);// D(p_t)
	
	SP=0;
	for(l=0;l<nstates;l++){
	  for(k=0;k<nstates;k++){
	    SP += Pcur[ l*nstates + k];
	  }
	}
	printf("SPcur=%lf\t",log(SP));

	double *Ptmp2 = malloc(nstates*nstates*sizeof(double));
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nstates, nstates, nstates, C, P, nstates, Pcur, nstates, 0, Ptmp2, nstates);//Ptmp2 = Phi_t*D(p_t)
	SP=0;
	for(l=0;l<nstates;l++){
	  for(k=0;k<nstates;k++){
	    SP += Ptmp2[ l*nstates + k];
	  }
	}
	printf("SPtmp2=%lf\t",log(SP));

	double *Ptmp = malloc(nstates*sizeof(double));
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 1, nstates, nstates, C, Pold, nstates, Ptmp2, nstates, 0, Ptmp, nstates);//Ptmp = Pold*Ptmp2
	
	SP=0;
	for(l=0;l<nstates;l++){
	  SP += Pold[ l];
	}
	printf("SPold=%lf\t",log(SP));
	
	free(Ptmp2);
	//printf("Pold*Pcur done\n");
	free(Pcur);
	memcpy(&Pold[0], &Ptmp[0], sizeof(double)*nstates);
	free(Ptmp);
	snp +=nstates;

	double SPtmp=0;
	for(l=0;l<nstates;l++){
	  //for(k=0;k<nstates;k++){
	  SPtmp += Pold[l]; //+ l*nstates];
	    //if(Pold[ l*nstates + k]==0) nSPcur ++;
	    //}
	}
	printf("SPtmp=%lf\n",log(SPtmp));
	if(SPtmp==0)  break;
	//printf("P[%d]=%lf\t", 0,Pold[0] ); //always 1?
	//printf("bef-endloop-Pold[%d]=%lf\t", 0,log(Pold[0]) ); //always 1?
      }//j

      //printf("af-endloop-Pold[%d]=%lf\t", 0,log(Pold[0]) ); //always 1?
      for(k=0;k<nstates;k++){
	//for(l=0;l<nstates;l++){
	    //if(!isinf(log(Pold[k*nstates + l]))) printf("logLik-Pold[%d]=%.10lf\t", k*nstates+l,log(Pold[k*nstates + l]) ); //always 1?
	    //if(!isinf(log(Pold[k*nstates + l]))) printf("logylik[%d]=%.10lf\t", k*nstates+l,log(ylik(Y,J,0,n,p,piall[k]) )); //always 1?
	    /*if(ylik(Y,J,0,n,p,piall[k])!=0){
	      int ii;
	      for(ii=0;ii<n;ii++) printf("%d/%d,%d\t", Y[0][ii],J[0][ii],  piall[k][ii]); //always 1?
	      printf("\n");
	      }*/
	Lik[ie][ic] += exp( log(Pold[k]) + log(C));
	    //if(Pold[k*nstates + l]>0) printf("Lik[%d][%d]temp=%.30e\n",ie,ic,Pold[k*nstates + l]); //=> always -Inf?
	    //}//l
      }//k
      printf("L done\n");
      free(Pold);
      printf("LLik[%d][%d]=%lf\t",ie,ic,log(Lik[ie][ic]));
      Lik[ie][ic] = log(Lik[ie][ic]); //compute log likelihood to simplify display of output data
      printf("Lik[%d][%d]=%.30e\n",ie,ic,Lik[ie][ic]); 
    }//ic
    printf("%.2f%% done\n",((float)ie+1-stloop)*100.0/(endloop-stloop));
  }//ie
  
  //free memory
  free(priorst);
  for(j=0;j<10*tmax*maxnstates;j++) free(piall[j]);
  free(piall);
  //free(piall2);
  for(i = 0; i < n; i++) free(M[i]);
  free(M);
  free(pCtmp);
  //free(all2short);
  //for(j=0;j<tmax;j++) free(simppstates[j]);
  //free(simppstates);
  free(P);
  for(i=0; i<tmax; i++){
    free(J[i]);
    free(Y[i]);
    free(prmin[i]);
  }
  for(k=0; k<nstatestmp;k++) free(piall2[k]);
  free(piall2);
  free(J);
  free(Y);
  free(D);
  free(prmin); 
  
  //start sending data
  //total Likelihood
  double Ltot=0;
  double coef;
  for(k=0; k<nstepe; k++){
    for(l=0; l<nstepc; l++){
      coef = 1;
      if( (k==0)|(k==(nstepe-1))) coef *= 0.5;
      if( (l==0)|(l==(nstepc-1))) coef *= 0.5;
      Ltot += exp(Lik[k][l])*coef;
    }
  }
  Ltot = 2*log(win) + log(Ltot);
  printf("Total log-likelihood=%.5lf\n",Ltot);
  //print results
  FILE * fe;
  printf("Writing output in file %s... ",fout);
  fe=fopen(fout,"wb");  
  for(i=0; i<nstepe; i++){
    for(j=0; j<nstepc; j++){
      if(Ltot==0) fprintf(fe,"%.20lf\t", Lik[i][j] );
      else fprintf(fe,"%.20lf\t", Lik[i][j] ); //-Ltot) ); //do not normalize if different values of p are to be compared
    }
    fprintf(fe,"\n");
  }
  fclose(fe);
  time (&end);//records ending time
  dif = difftime (end,start);
  printf ("done\n Total running time: %.2lf min\n", dif/60.0 );
  
  return 0;
}
