#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <cblas.h>
#include <mpi.h>
#include <unistd.h>
#include <float.h>


/*_________________________________________________________*/
/**
 * @author Nicolas Alcala;  
 */

#define send_data_tag 2001
#define return_data_tag 2002

void compPePc(double * Pee, double * Pcc, double * D,  unsigned int ** piall, double e, double * pC, double ** M, unsigned int n, unsigned int nstates, unsigned int j){
  //compute probability transition from state piold to all states in piall
  unsigned int i = 0;
  double etmp = e;
  if(etmp>1) etmp = 1;
  unsigned int * pitmp = piall[j];
  for(i=0;i<nstates;i++){// for all possible states at time point k
    unsigned int * piold = piall[i];

    double pijc = 0; //1;
    unsigned int s1 = 0; //1 -> 0
    unsigned int s2 = 0; //1 -> 1
    unsigned int k;
    
    unsigned int ispos = 1;
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
      Pee[ i*nstates + j] = pow(etmp,s1)*pow(1-etmp,s2);
      Pcc[ j*nstates + i]  = exp(pijc);
    }else{
      Pee[ i*nstates + j] = 0;
      Pcc[ j*nstates + i]  = 0;
    }//ispos
  }//i
}

double ylik(int ** Y, int ** J, int j,unsigned int n,double p, unsigned int * pi){
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

unsigned int pi_to_id(int* pi,unsigned int n){
  //converts state pi to identifier (in binary)
  unsigned int ip;
  unsigned int res=0;
  for(ip=0;ip<n;ip++) res += pi[ip]*pow(2,n-ip-1);
  return res;
}


int main(int argc, char ** argv)//takes the path to an input file as argument
{
  /*--------------MPI Parameters-------------------*/
  MPI_Status status;
  int my_id, root_process, ierr, num_procs, an_id, avg_rows_per_process;
  /* Creation of parallel processes */
  ierr = MPI_Init(&argc, &argv);
  root_process = 0;

  /* find out process ID and total number of processes */
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  if(my_id == root_process)
    printf("------ MIDASPOM, beta MPI version ------\n-> N. Alcala, E. M. Cole, and N. A. Rosenberg <-\n");
    
  /*-------------- model and simulation Parameters--------------------*/
  float prioroc = 0.5; // prior occupancy for missing data the first year
  char* fname = "input.txt";
  char* fout  = "posterior.txt";
  char* Jname  = "J.txt";
  double d = 100;
  unsigned int nstep = 101;
  double a = 1.0/400; //dispersal distance
  double ecmin = 0.0;
  double ecmax = 1.0;
  double p     = 1; //proba of observation given presence
  double C     = 100000;//scaling constant
  
  int c;
  
  opterr = 0;
  while ((c = getopt (argc, argv, "m:p:d:i:J:o:s:l:u:f:")) != -1)
    switch (c)
      {
      case 'm':
        a = 1.0/atof(optarg);
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
        nstep = atoi(optarg);
	break;
      case 'l':
        ecmin = atof(optarg);
        break;
      case 'u':
        ecmax = atof(optarg);
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

  double win = (ecmax-ecmin)/(nstep-1);

  if(my_id == root_process) printf("Parameters for numerical approximation of the posterior density:\n\tWindow size=%lf, number of steps=%d\n",win,nstep );
 
  //records the time of execution
  time_t start,end;         //beginning and ending times
  double dif;               //difference between beginning and ending
  
  //model parameters
  unsigned int i,j,k,l;

  //read observed data
  FILE * fo;
  if(my_id == root_process) printf("Reading observations from file %s... ",fname);
  unsigned int tmax=0;
  fo = fopen(fname,"rb");   //
  unsigned int n=1;
  while ( (c=fgetc(fo)) != EOF ) {
    if ( c == '\n' )
      tmax++;
    if(tmax==0){
      if( (c==' ')|(c=='\t') ) n++;
    }
  }
  fclose(fo);
  
  int** piobs = malloc(tmax*sizeof(int*)); //piobs[tmax][n]; 
  for(i=0; i<tmax; i++) piobs[i] = malloc(n*sizeof(int));

  int** Y     = malloc(tmax*sizeof(int*)); //piobs[tmax][n]; 
  for(i=0; i<tmax; i++) Y[i] = malloc(n*sizeof(int));
  
  fo = fopen(fname,"rb");   //
  int out;
  unsigned int var[n];
  for(j=0; j<n; j++) var[j]=0;
  for(i=0; i<tmax; i++){
    for(j=0; j<n; j++){
      out = fscanf(fo,"%d",&(Y[i][j]));
      if(Y[i][j]!=0){
	var[j]=1;
	piobs[i][j] = 1;
      }else{
	var[j]=1;
	piobs[i][j] = 0;
      }
    }
  }
  fclose(fo);
  if(my_id == root_process) printf("done\n");

  if(my_id == root_process) printf("Number of habitat patches: %d\nNumber of sampled years: %d\n",n,tmax);

  FILE * fJ;
  if(my_id == root_process) printf("Reading survey number from file %s... ",Jname);
  fJ = fopen(Jname,"rb");   //

  int** J = malloc(tmax*sizeof(int*)); //piobs[tmax][n]; 
  for(i=0; i<tmax; i++) J[i] = malloc(n*sizeof(int));
  
  for(i=0; i<tmax; i++){
    for(j=0; j<n; j++){
      out = fscanf(fJ,"%d",&(J[i][j]));
      if(J[i][j]==0){
	piobs[i][j] = -1;
	var[j] = 1;
      }
    }
  }
  fclose(fJ);
  if(my_id == root_process) printf("done\n");
  
  unsigned int nvar=0; //number of variable segments
  for(j=0; j<n; j++) if(var[j]==1) nvar++;

  unsigned int nstates = pow(2,nvar); //number of possible states
  
  avg_rows_per_process = nstep/num_procs;
  double **M = malloc(n*sizeof(double*));
  for(i = 0; i < n; i++) M[i] = malloc(n*sizeof(double));
  
  for(i=0; i<n; i++){
    for(j=i; j<n; j++){
      if(i==j) M[i][j]=0.0;
      else{
	M[i][j]=exp(-a*(j-i)*d);
	M[j][i]=exp(-a*(j-i)*d);
      }
    }
  }
  if(my_id == root_process){
    printf("Dispersal matrix:\n");
    for(i=0; i < n; i++){
      for(j=0; j<n; j++){
	printf("%.3f ",M[i][j]);
      }
      printf("\n");
    }
  }

  //all possible states
  unsigned int** piall = malloc(nstates*sizeof(unsigned int*));
  for(i=0; i<nstates; i++) piall[i] = malloc(n*sizeof(unsigned int));//piall[nstates][n]; 
  for(i=0; i<nstates; i++){
    unsigned int jtmp = 0;
    for(j=0; j<n;j++){
      if(var[j]==0){
	piall[i][j] = 0;
      }else{ 
	unsigned int a = pow(2,nvar-jtmp-1);
	piall[i][j] = i/a%2;
	jtmp++;
      }
    }
  }
       
  //number of possible states at each time step
  unsigned int npstates[tmax]; 
  unsigned int nnpstates = 0;
  unsigned int** pstates = malloc(tmax*sizeof(unsigned int*)); //list of state id in all states list
  unsigned int** simppstates = malloc(tmax*sizeof(unsigned int*)); //list of state id in states short list
  float* priorst;//prior probability of each initial state
  unsigned int s1;
  unsigned int st1;
  unsigned int nextid = 0;
  for(i=0; i<tmax; i++){
    s1 = n;//0;
    /*for(j=0; j<n; j++){
      if(piobs[i][j]==-1){
	s1++; //number of segments with missing data
      }
      }*/
    npstates[i] = pow(2,s1); //number of possible states
    nnpstates += npstates[i]; //total number of possible states
    if(i==0){
      priorst = malloc( npstates[0]*sizeof(float));
      for(k=0; k<npstates[0]; k++) priorst[k]=1;
    }
    pstates[i]     = malloc(npstates[i]*sizeof(unsigned int));
    simppstates[i] = malloc(npstates[i]*sizeof(unsigned int));
    for(k=0; k<npstates[i]; k++) pstates[i][k]=0;
    s1 = 0;
    unsigned int jtmp = 0;
    for(j=0; j<n; j++){
      if(piobs[i][j]==-1) s1 ++;
      for(k=0; k<npstates[i]; k++){
	if(piobs[i][j]>-1){
	  pstates[i][k] += (unsigned int) (piobs[i][j]*pow(2,nvar-jtmp-1)); // to change if some sites always at 0
	  //printf("j=%d: pstates[%d][%d]=%d\t",j,i,k,pstates[i][k] ); 
	}else{
	  st1 = npstates[i]/pow(2,s1);
	  pstates[i][k] += (unsigned int) ((k/st1%2)*pow(2,nvar-jtmp-1));
	  if(i==0){
	    priorst[k] *= (k/st1%2)*prioroc + (1-k/st1%2)*(1-prioroc);
	  }
	}  
      }//for state k
      if(var[j]>0) jtmp++;
    }//for patch j
    
    if(i==0){
      for(k=0; k<npstates[0]; k++){
	simppstates[0][k] = nextid;
	nextid++;
      }
    }else{
      for(k=0; k<npstates[i]; k++){
	//printf("pstates[%d][%d]=%d\t",i,k,pstates[i][k]); // OK
	int ispresent = 0;
	for(j=0; j<i; j++){
	  for(l=0; l<npstates[j]; l++){
	    if(pstates[i][k]==pstates[j][l]){
	      simppstates[i][k] = simppstates[j][l];
	      ispresent = 1;
	    }
	  }//l
	}//j
	if(ispresent == 0){
	  simppstates[i][k] = nextid;
	  nextid++;
	}
	
      }//k
    }
  }//for time point i

  //printf("nextid=%d\n",nextid);
  
  unsigned int* all2short = malloc(nextid*sizeof(unsigned int)); //correspondance between table ID piall and short table
  for(i=0; i<tmax; i++){
    for(k=0; k<npstates[i]; k++){
      all2short[simppstates[i][k]] = pstates[i][k];
    }
  }

  for(i=0; i<tmax; i++) free(pstates[i]);
  free(pstates);

  //print input
  if(my_id == root_process){
    printf("Input occupancy data:\n");
    for(i=0; i<tmax; i++){
      printf("Year %d: ",i);
      for(j=0; j<n; j++) printf("%d/%d=%d\t",Y[i][j],J[i][j],piobs[i][j]);
      printf("\n");
    }

    printf("Number of possible states per year:\n");
    for(i=0; i<tmax; i++){
      printf("Year %d: %d\n",i,npstates[i]);
    }
  }
  if(my_id == root_process) printf("Number of states to compute: %d\n",nstates);
  
  unsigned int st = pi_to_id(piobs[0],n);

  for(i=0; i<tmax; i++) free(piobs[i]);
  free(piobs);
  
     
  //create input tables
  double cllik[nstep];
  double ellik[nstep];
  for(i=0;i<nstep-1;i++){
    ellik[i] = ((double)i)*win + ecmin; 
    cllik[i] = ((double)i)*win + ecmin; 
  }
  ellik[nstep-1] = ecmax;
  cllik[nstep-1] = ecmax;
  //create output tables
  double Lik[nstep][nstep];
  for(i=0;i<nstep;i++){
    for(j=0;j<nstep;j++){
      Lik[i][j] = 0;
    }
  }
  setbuf(stdout, NULL);

  //main loop
  time (&start);//records the beginning time
  srand (time(NULL));
  
  printf("Starting parallel likelihood computation process %i/%i\n", my_id+1, num_procs);
  unsigned int ie,ic;
  double etmp,ctmp;//Ltmp
  
  // start loop
  unsigned int stloop,endloop;
  if(my_id==root_process){
    stloop  = 0; 
    endloop = avg_rows_per_process + nstep%num_procs;
  }else{
    stloop  = my_id*avg_rows_per_process + nstep%num_procs;
    endloop = (my_id+1)*avg_rows_per_process + nstep%num_procs;
  }
  double *P  = malloc(nstates*nstates*sizeof(double));
  double *D  = malloc(nstates*tmax*sizeof(double)); 
  double *pCtmp = malloc(n*sizeof(double)); //proba of colonization

  for(ie=stloop; ie<endloop; ie++){
    etmp = ellik[ie];
    for(ic=0; ic<nstep; ic++){
      //printf("%d %d\n",ie,ic);
      ctmp = cllik[ic];

      double *Pe = calloc(nstates*nstates,sizeof(double));
      double *Pc = calloc(nstates*nstates,sizeof(double));
      
      unsigned int snp;
      for(j=0;j<nstates;j++){
	for(k=0; k<n; k++){//for all patches
	  double s1 = 0;
	  for(l=0; l<n; l++){
	    if(l!=k) s1 += M[l][k]*piall[j][l];
	  }//l
	  pCtmp[k] = ctmp*s1;
	  if(pCtmp[k]>1) pCtmp[k] = 1;
	}//k
	
	compPePc(Pe, Pc, D, piall, etmp, pCtmp, M, n, nstates, j);
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
      //printf("PePc done\n");
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nstates, nstates, nstates, 1.0, Pe, nstates, Pc, nstates, 0, P, nstates); //multiply Pe and Pc using cblas routine
      /*double SPec=0;
      //printf("\n\tP:");
      for(j=0;j<nstates;j++){
	for(k=0;k<nstates;k++){
	  SPec += P[ j*nstates + k];
	  //printf("%e\t",P[j*nstates + k]);//SPe OK, SPc weird
	  //if(P[j*nstates + k]==0) printf("error!");
	}
	//printf("\n");
      }
      //printf("SPec=%lf\t",SPec);//SPe OK, SPc weird
      */
      //printf("P done\n");
      free(Pe);
      free(Pc);

      double *Pold = malloc(nstates*nstates*sizeof(double));
      for(k=0;k<nstates;k++) Pold[k] = C*ylik(Y,J,0,n,p,piall[k])*priorst[k];
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

	for(k=0; k<nstates; k++) Pcur[k*nstates + k] = C*ylik(Y,J,j,n,p,piall[k]);
	
	/*double SP=0;
	for(l=0;l<nstates;l++){
	  for(k=0;k<nstates;k++){
	    SP += Pcur[ l*nstates + k];
	  }
	}
	printf("SPcur=%lf\t",log(SP));*/

	double *Ptmp2 = malloc(nstates*nstates*sizeof(double));
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nstates, nstates, nstates, C, P, nstates, Pcur, nstates, 0, Ptmp2, nstates);
	/*SP=0;
	for(l=0;l<nstates;l++){
	  for(k=0;k<nstates;k++){
	    SP += Ptmp2[ l*nstates + k];
	  }
	}
	printf("SPtmp2=%lf\t",log(SP));*/

	double *Ptmp = malloc(nstates*nstates*sizeof(double));
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nstates, nstates, nstates, C, Pold, nstates, Ptmp2, nstates, 0, Ptmp, nstates);
	/*SP=0;
	for(l=0;l<nstates;l++){
	  for(k=0;k<nstates;k++){
	    SP += Ptmp[ l*nstates + k];
	  }
	}
	printf("SPtmp=%lf\t",log(SP));*/
	
	free(Ptmp2);
	//printf("Pold*Pcur done\n");
	free(Pcur);
	memcpy(&Pold[0], &Ptmp[0], sizeof(double)*nstates*nstates);
	free(Ptmp);
	snp +=nstates;

	/*double SPcur=0;
	int nSPcur=0;
	for(l=0;l<nstates;l++){
	  for(k=0;k<nstates;k++){
	    SPcur += Pold[ l*nstates + k];
	    if(Pold[ l*nstates + k]==0) nSPcur ++;
	  }
	  }*/
	//printf("SPcur=%lf, nSPcur=%d\t",log(SPcur),nSPcur);
	//printf("P[%d]=%lf\t", 0,Pold[0] ); //always 1?
	//printf("bef-endloop-Pold[%d]=%lf\t", 0,log(Pold[0]) ); //always 1?
      }//j

      //printf("af-endloop-Pold[%d]=%lf\t", 0,log(Pold[0]) ); //always 1?
      for(k=0;k<nstates;k++){
	  for(l=0;l<nstates;l++){
	    //if(!isinf(log(Pold[k*nstates + l]))) printf("logLik-Pold[%d]=%.10lf\t", k*nstates+l,log(Pold[k*nstates + l]) ); //always 1?
	    //if(!isinf(log(Pold[k*nstates + l]))) printf("logylik[%d]=%.10lf\t", k*nstates+l,log(ylik(Y,J,0,n,p,piall[k]) )); //always 1?
	    /*if(ylik(Y,J,0,n,p,piall[k])!=0){
	      int ii;
	      for(ii=0;ii<n;ii++) printf("%d/%d,%d\t", Y[0][ii],J[0][ii],  piall[k][ii]); //always 1?
	      printf("\n");
	      }*/
	    Lik[ie][ic] += exp( log(Pold[k*nstates + l]) + log(C));
	    //if(Pold[k*nstates + l]>0) printf("Lik[%d][%d]temp=%.30e\n",ie,ic,Pold[k*nstates + l]); //=> always -Inf?
	  }//l
      }//k
      //printf("L done\n");
      free(Pold);
      printf("LLik[%d][%d]=%lf\t",ie,ic,log(Lik[ie][ic])); //=> always 0?
      Lik[ie][ic] = log(Lik[ie][ic]); //compute log likelihood to simplify display of output data
      //printf("Lik[%d][%d]=%.30e\n",ie,ic,Lik[ie][ic]); //=> always -Inf?
    }//ic
    if(my_id==root_process) printf("%.2f%% done\n",((float)ie+1-stloop)*100.0/(endloop-stloop));
  }//ie
  
  //free memory
  free(priorst);
  for(j=0;j<nstates;j++) free(piall[j]);
  free(piall);
  for(i = 0; i < n; i++) free(M[i]);
  free(M);
  free(pCtmp);
  free(all2short);
  for(j=0;j<tmax;j++) free(simppstates[j]);
  free(simppstates);
  free(P);
  for(i=0; i<tmax; i++){
    free(J[i]);
    free(Y[i]);
  }
  free(J);
  free(Y);
  free(D);

  printf("end likelihood computation process %i/%i\n", my_id+1, num_procs);
  //end root loop

  //start sending data
  if(my_id!=root_process){
    printf("Sending data (proc %d)... ",my_id);
    double MSGsend[2+avg_rows_per_process*nstep];
    MSGsend[0] = my_id*avg_rows_per_process + nstep%num_procs ;
    MSGsend[1] = (my_id+1)*avg_rows_per_process + nstep%num_procs ;
    for(k=MSGsend[0],j=0; k<MSGsend[1]; k++,j++){
      for(l=0; l<nstep; l++){
	MSGsend[2 + j*nstep + l] = Lik[k][l];
      }
    }
    ierr = MPI_Send(MSGsend, 2+avg_rows_per_process*nstep, MPI_DOUBLE, root_process, return_data_tag, MPI_COMM_WORLD); 
    printf("done\n");
  }else{
    printf("Gathering data from %d proc... ",num_procs-1);
    for(an_id = 1; an_id < num_procs; an_id++) {//for each proc
      double MSGrec[2+avg_rows_per_process*nstep];
      ierr = MPI_Recv(MSGrec, 2+avg_rows_per_process*nstep, MPI_DOUBLE, MPI_ANY_SOURCE, return_data_tag, MPI_COMM_WORLD, &status); //receive data
      for(k=MSGrec[0],j=0; k<MSGrec[1]; k++,j++){
	for(l=0; l<nstep; l++){
	  Lik[k][l] = MSGrec[2 + j*nstep + l]; //store data
	}
      }
    }
    printf("done\n");

    //total Likelihood
    double Ltot=0;
    double coef;
    for(k=0; k<nstep; k++){
      for(l=0; l<nstep; l++){
	coef = 1;
	if( (k==0)|(k==(nstep-1))) coef *= 0.5;
	if( (l==0)|(l==(nstep-1))) coef *= 0.5;
	Ltot += exp(Lik[k][l])*coef;
      }
    }
    Ltot = 2*log(win) + log(Ltot);
    printf("Total log-likelihood=%.5lf\n",Ltot);
    //print results
    FILE * fe;
    printf("Writing output in file %s... ",fout);
    fe=fopen(fout,"wb");  
    for(i=0; i<nstep; i++){
      for(j=0; j<nstep; j++){
	if(Ltot==0) fprintf(fe,"%.20lf\t", Lik[i][j] );
	else fprintf(fe,"%.20lf\t", Lik[i][j] ); //-Ltot) ); //do not normalize if different values of p are to be compared
      }
      fprintf(fe,"\n");
    }
    fclose(fe);
    time (&end);//records ending time
    dif = difftime (end,start);
    printf ("done\n Total running time: %.2lf min\n", dif/60.0 );
  }
  //end MPI
  ierr = MPI_Finalize();
  
  return 0;
}
