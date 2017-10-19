#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <cblas.h>
#include <mpi.h>
#include <unistd.h>

/*_________________________________________________________*/
/**
 * @author Nicolas Alcala;  
 */

#define send_data_tag 2001
#define return_data_tag 2002

void matpow(double *x, int n, int k, double *z){
  //matrix power
  if (k == 0) { /* return identity matrix */
    int i, j;
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
	z[i * n + j] = (i == j) ? 1.0 : 0.0;
    return;
  }else if (k < 0) {
    error("power must be a positive integer; use solve() directly for negative powers");
  }else { /* k >= 1 */
    int nSqr = n * n;
    double *tmp = malloc(nSqr*sizeof(double)); /* temporary matrix */
    /* Take powers in multiples of 2 until there is only one product left to make. That is, if k = 5, compute (x * x), then ((x * x) * (x * x)) and finally ((x * x) * (x * x)) * x. */
    memcpy(&z[0], &x[0], sizeof(double)*n*n);
    k--;
    while (k > 0) {
      if (k & 1) { /* z := z * x */
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n, n, n, 1.0, x, n, z, n, 0, tmp, n);
	memcpy(&z[0], &tmp[0], sizeof(double)*n*n);
      }//if k&1
      if(k == 1){
	break;
      }
      k >>= 1; /* efficient division by 2; now have k >= 1 */
      /* x := x * x */
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n, n, n, 1.0, x, n, x, n, 0.0, tmp, n);
      memcpy(&x[0], &tmp[0], sizeof(double)*n*n);
    }
    free(tmp);
    return;
  }
}

double pije(int* piold, int* pitmp, double e, double K, int n){
  //extinction phase
  int s1 = 0; //1 -> 0
  int s2 = 0; //1 -> 1
  int k;
  double E = e/K;
  if(E>1) E=1;
  for(k=0; k<n; k++){
    if(pitmp[k]*(1-piold[k])>0) return 0;//0 -> 1 impossible!
    s1 += (1-pitmp[k])*piold[k];
    s2 += pitmp[k]*piold[k];
  }
  return pow(E,s1)*pow(1-E,s2);
}

double pijc(int* pitmp, int* pinew, double c, double K, double** M, int n){
  //colonization phase
  double pCi[n]; 
  double res=1;
  int k;
  for(k=0; k<n; k++){
    if( (pitmp[k]*(1-pinew[k]))>0 ) return 0;
    double s1 = 0;
    int l;
    for(l=0; l<n; l++){
      if(l!=k) s1 += M[l][k]*pitmp[l];
    }
    pCi[k] = c*s1*K;
    if(pCi[k]>1) pCi[k] = 1;
    res *= pitmp[k] + (1-pitmp[k])*(1-pinew[k])*(1-pCi[k]) + (1-pitmp[k])*(pinew[k])*pCi[k];
  }
  return res;
}

unsigned int simE(unsigned int* piold, unsigned int* pitmp, double e, double K, unsigned int n){
  //extinction phase
  unsigned int k;
  double E = e/K;
  double pp;
  if(E>1) E=1;
  for(k=0; k<n; k++){
    if(piold[k]==1){
      pp = (double)rand()/(double)(RAND_MAX);
      if(pp>E) pitmp[k] = 1;
      else pitmp[k] = 0;
    }else{
      pitmp[k] = 0;
    }
    //printf("%d ",pitmp[k]);
  }
  return 1;
}

unsigned int simC(unsigned int* pitmp, unsigned int* pinew, double c, double K, double** M, unsigned int n){
  unsigned int k;
  double pp;
  double pCi; 
  unsigned int res=0;
  for(k=0; k<n; k++){
    if(pitmp[k]==0){
      double s1 = 0;
      unsigned int l;
      for(l=0; l<n; l++){
	if(l!=k) s1 += pitmp[l]*K*M[l][k];//*;
      }
      pCi = c*s1;
      if(pCi>1) pCi = 1;
      pp = (double)rand()/(double)(RAND_MAX);
      if(pp<pCi) pinew[k] = 1;
      else pinew[k] = 0;
    }else{
      pinew[k] = 1;
    }
  }
  return res;
}


int pi_to_id(int* pi,int n){
  int ip;
  int res=0;
  for(ip=0;ip<n;ip++) res += pi[ip]*pow(2,n-ip-1);
  return res;
}

int main(int argc, char ** argv)//takes the path to an input file as argument
{
  /*--------------MPI Parameters-------------------*/
  MPI_Status status;
  int my_id, root_process, ierr, num_procs, an_id, avg_rows_per_process, sender;
  /* Creation of parallel processes */
  ierr = MPI_Init(&argc, &argv);
  root_process = 0;

  /* find out process ID and total number of processes */
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  if(my_id == root_process)
    printf("------ MIDASPOM, in situ die-off hypothesis, beta MPI version -------\n-> N. Alcala, E. M. Cole, N. A. Rosenberg <-\n");
    
  /*----------------Parameters--------------------*/
  //records the time of execution
  time_t start,end;         //beginning and ending times
  double dif;               //difference between beginning and ending
  

  //model parameters
  //int n; //size of the Creek
  int ts=20; //number of time steps before the event increasing die-off
  int tdis; //number of time steps after the event

  unsigned int nstep=151;
  double eB; //extinction rate
  double cB; //colonization rate
  double Kmin = 0.1;
  double Kmax = 100.0;
  double a = 1.0/400.0; //dispersal distance
  float prioroc = 0.5; // prior occupancy for missing data the first year
  char* fname = "input.txt";
  char* fout  = "lh_dieoff.txt";
  char* fstname = "piall.txt";
  double d = 200;
  int readpi = 0;
  //float win = 0.01;
      
  int i,j,k,l;
  unsigned int index;
  int c;
  
  opterr = 0;
  while ((c = getopt (argc, argv, "b:a:e:c:m:p:d:i:o:s:l:r:u:w:")) != -1)
    switch (c)
      {
      case 'b':
        ts = atoi(optarg);
        break;
      case 'a':
        tdis = atoi(optarg);
        break;
      case 'e':
        eB = atof(optarg);
        break;
      case 'c':
        cB = atof(optarg);
        break;
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
      case 'o':
        fout = optarg;
        break;
      case 's':
        nstep = atoi(optarg);
	break;
      case 'l':
        Kmin = atof(optarg);
        break;
      case 'u':
        Kmax = atof(optarg);
        break;
      case 'r':
	readpi = 1;
	fstname = optarg;
	break;
      case 'w':
	fstname = optarg;
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
 
  //nstep = (unsigned int)(round(ceil(1/win) +1));

  avg_rows_per_process = nstep/num_procs;

  if(my_id == root_process) printf("%d years before increased die-off, %d years after increased die-off\n",ts,tdis);

  //read observed data
  FILE * fo;
  if(my_id == root_process) printf("Reading observations from file %s... ",fname);
  fo = fopen(fname,"rb");   //
  unsigned int n=1;
  while ( (c=fgetc(fo)) != EOF ) {
    if ( c == '\n' )
      break;
    if( (c==' ')|(c=='\t') ) n++;
  }
  fclose(fo);

  int nstates = pow(2,n); //number of possible states
  if(my_id == root_process) printf("%d patches\n",n);
  
  double* pendp   = malloc(n*sizeof(double));
  int nstart = n;
  fo = fopen(fname,"rb");   //
  int out;
  for(j=0; j<n; j++){//scan only initial state
    out = fscanf(fo,"%lf",&(pendp[j]));
    if(my_id == root_process) printf("%lf\t",pendp[j]); //for time point i
  }
  fclose(fo);
  if(my_id == root_process) printf("\n");


  double Kllik[nstep];
  for(i=0;i<nstep;i++){
    Kllik[i] = pow(10.0,((double)i)/(nstep-1)*(log10(Kmax)-log10(Kmin))+log10(Kmin)); //((double)i)*Kmax/(nstep);
    //if(my_id == root_process) printf("K[%d]=%lf\t",i,Kllik[i]);
  }
  
  
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
    printf("Migration matrix:\n");
    for(i=0; i < n; i++){
      for(j=0; j<n; j++){
	printf("%.3f ",M[i][j]);
      }
      printf("\n");
    }
  }
  
  
  //all possible states
  int** piall;
  if(!readpi){
    int thin = 10;
    int tmax = tdis+ts;
    nstates = 2*(nstep/thin+1)*tmax+1;
    printf("nstates=%d\n",nstates);
    piall = malloc((nstates+1)*sizeof(int*));
    for(i=0; i<nstates; i++){
      //printf("i=%d\t",i);
      piall[i] = malloc(n*sizeof(int));//piall[nstates][n];
    }
    int idpiall[nstates];
    /*for(i=0; i<nstates; i++){
      for(j=0; j<n;j++){
      int a = pow(2,n-1)/pow(2,j);
      piall[i][j] = i/a%2;
      }
      }*/
    
    //SIMULATE STATES
    k=0;
    double pp;
    int kt;
    for(kt=0;kt<nstep;kt++){// for all values, simulate 2X maxnstates/time point
      if(kt%thin==0){
	printf("kt=%d, k=%d\t",kt,k);
	//initial value
	for(l=0; l<n; l++){
	  pp = (double)rand()/(double)(RAND_MAX);
	  idpiall[k*tmax*2] = 0;
	  if(pp>0.5){
	    piall[k*tmax*2][l] = 1;
	    idpiall[k*tmax*2] += pow(2,n-l-1);
	  }else{
	    piall[k*tmax*2][l] = 0;
	  }
	}
	for(j=0; j<ts; j++){
	  idpiall[2*k*tmax+2*j+1] = 0;
	  idpiall[2*k*tmax+2*j+2] = 0;
	  simE(piall[2*k*tmax+2*j],   piall[2*k*tmax+2*j+1], eB, Kllik[kt],n);
	  simC(piall[2*k*tmax+2*j+1], piall[2*k*tmax+2*j+2], cB, Kllik[kt],M,n);
	  for(l=0; l<n; l++){
	    if(piall[2*k*tmax+2*j+1][l]==1) idpiall[2*k*tmax+2*j+1] += pow(2,n-l-1);
	    if(piall[2*k*tmax+2*j+2][l]==1) idpiall[2*k*tmax+2*j+2] += pow(2,n-l-1);
	  }
	  //printf("2*k*tmax+2*j+2=%d\t",2*k*tmax+2*j+2);
	}
	for(j=ts; j<tmax; j++){
	  idpiall[2*k*tmax+2*j+1] = 0;
	  idpiall[2*k*tmax+2*j+2] = 0;
	  simE(piall[2*k*tmax+2*j],   piall[2*k*tmax+2*j+1], eB, 1,n);
	  simC(piall[2*k*tmax+2*j+1], piall[2*k*tmax+2*j+2], cB, 1,M,n);
	  //printf("2*k*tmax+2*j+2=%d\t",2*k*tmax+2*j+2);
	}
	
	k++;
      }
      // add state corresponding to observed data
      idpiall[nstates-1] = 0;
      for(l=0; l<n; l++){
	pp = (double)rand()/(double)(RAND_MAX);
	if( pp < pendp[l] ){
	  piall[nstates-1][l] = 1;
	  idpiall[nstates-1] += pow(2,n-l-1);
	}else{
	  piall[nstates-1][l] = 0;
	}
      }
    }

    //find duplicate states
    printf("start\n");    
    int nstatestmp = nstates;
    int idkeep[nstatestmp];
    for(j=0;j<nstatestmp;j++) idkeep[j] = 1;
    for(j=0;j<nstates-1;j++){
      for(k=j+1;k<nstates;k++){
	if(idkeep[j]){
	  if(idpiall[j]==idpiall[k]){
	    idkeep[k] = 0;
	    nstatestmp--;
	  }
	}
      }
    }
    
    l=0;
    for(j=0;j<nstates;j++){
      if(idkeep[j]){
	if(l!=j) piall[l] = piall[j];
	l++;
      }
    }
     
    for(j=0;j<nstatestmp;j++){
      for(k=0;k<n;k++) printf("%d ", piall[j][k] );
      printf("\n");
    }
    
    nstates=nstatestmp;
    printf("Number of states to compute: %d; total matrix size: %d\n",nstates,nstates*nstates);

    
    //
    FILE * fpi;
    printf("Writing list of occupancy states in file %s\n",fstname);
    fpi=fopen(fstname,"wb");  
    for(i=0; i<nstates; i++){
      for(j=0; j<n; j++){
	fprintf(fpi,"%d\t", piall[i][j] ); //-Ltot) ); //do not normalize if different values of p are to be compared
      }
      fprintf(fpi,"\n");
    }
    fclose(fpi);

  }else{
    nstates = 0;
    FILE * fst;
    printf("Reading observations from file %s... ",fstname);
    fst = fopen(fstname,"rb");   //
    while ( (c=fgetc(fst)) != EOF ) {
      if ( c == '\n' )
	nstates++;
    }
    fclose(fst);
    printf("nstates=%d\n",nstates);
    piall  = malloc(nstates*sizeof(int*));
    for(i=0; i<nstates; i++) piall[i] = malloc(n*sizeof(int));
    
    fst = fopen(fstname,"rb");  
    for(i=0; i<nstates; i++){
      for(j=0; j<n; j++){
	out = fscanf(fst,"%d",&(piall[i][j]));
      }
    }
    fclose(fst);
    
    
    for(i=0; i<nstates; i++){
      for(j=0; j<n; j++){
	printf("%d\t",piall[i][j]);
      }
      printf("\n");
    }
  }
  
  
  //prior probability of each initial state
  printf("Start computing proba of observation in first year\n");
  float* priorst;
  int s1=0;
  int st1;
  for(j=0; j<n; j++){
    if(pendp[j]==-1){
      s1++; //number of occupied segments
    }
  }
  priorst = malloc( nstates*sizeof(float));
  for(k=0; k<nstates; k++){
    priorst[k]=1;
  }
  s1 = 0;
  for(j=0; j<n; j++){
    for(k=0; k<nstates; k++){
      //st1 = nstates/pow(2,s1);
      priorst[k] *= piall[k][j]*pendp[j] + (1-piall[k][j])*(1-pendp[j]);//(k/st1%2)*prioroc + (1-k/st1%2)*(1-prioroc);
      //printf("priorst[k]=%lf\n",priorst[k]);
    }//for state k
  }//for patch j

  free(pendp);  
  
  //create input tables
  //if(my_id == root_process) printf("Creating %dx%d input tables\n",nstates,nstates);
  double *P = malloc(nstates*nstates*sizeof(double));
  double *Pe = malloc(nstates*nstates*sizeof(double));
  double *Pc = malloc(nstates*nstates*sizeof(double));
  double *PK = malloc(nstates*nstates*sizeof(double));
  double *PKe = malloc(nstates*nstates*sizeof(double));
  double *PKc = malloc(nstates*nstates*sizeof(double));
  //if(my_id == root_process) printf("\n");
  
  //create output tables
  double Lik[nstep];
  for(i=0;i<nstep;i++){
    Lik[i] = 0;
  }
    
  //main loop
  time (&start);//records the beginning time
  srand (time(NULL));
  
  printf("Starting parallel likelihood computation process %i/%i\n", my_id+1, num_procs);
  double LLtmp,Ltmp2;
  // start root loop
  int stloop,endloop,t;
  if(my_id==root_process){
    stloop  = 0; 
    endloop = avg_rows_per_process + nstep%num_procs;
  }else{
    stloop  = my_id*avg_rows_per_process + nstep%num_procs;
    endloop = (my_id+1)*avg_rows_per_process + nstep%num_procs;
  }
  double *Ppow    = (double *) malloc(nstates*nstates*sizeof(double));
  double *PKpow   = (double *) malloc(nstates*nstates*sizeof(double));
  double *Ptotpow = (double *) malloc(nstates*nstates*sizeof(double)); //for product Ppow*PKpow
  
  // transition matrix after disease
  for(i=0;i<nstates;i++){
    for(j=0;j<nstates;j++){
      Pe[i*nstates + j] = pije(piall[i],piall[j],eB,1.0,n);
      Pc[i*nstates + j] = pijc(piall[i],piall[j],cB,1.0,M,n);
    }
  }
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nstates, nstates, nstates, 1.0, Pe, nstates, Pc, nstates, 0, P, nstates);
  free(Pe);
  free(Pc);
  matpow(P,nstates,tdis,Ppow);
  free(P);

  // transition before disease
  int iK;
  double Ktmp;
  for(iK=stloop; iK<endloop; iK++){//nstep; ie++){
    Ktmp = Kllik[iK];
    //printf("iK=%d\n",iK);
    LLtmp = 0;
    for(i=0;i<nstates;i++){
      for(j=0;j<nstates;j++){
	PKe[i*nstates + j] = pije(piall[i],piall[j],eB,Ktmp,n);
	PKc[i*nstates + j] = pijc(piall[i],piall[j],cB,Ktmp,M,n);
      }
    }
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nstates, nstates, nstates, 1.0, PKe, nstates, PKc, nstates, 0, PK, nstates);
    matpow(PK,nstates,ts,PKpow);
    double sumPKe=0,sumPKc=0,sumPK=0,sumPKpow=0;
    for(i=0;i<nstates;i++){
      for(j=0;j<nstates;j++){
	sumPKe += PKe[i*nstates + j];
	sumPKc += PKc[i*nstates + j];
	sumPK += PK[i*nstates + j];
	sumPKpow += PKpow[i*nstates + j];
      }
    }
        
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nstates, nstates, nstates, 1.0, PKpow, nstates, Ppow, nstates, 0, Ptotpow, nstates);
    
    for(i=0;i<nstates;i++){
      for(j=0; j<nstates; j++){
	LLtmp += Ptotpow[i*nstates + j]*priorst[j];
      }
    }
    Lik[iK] = LLtmp;
  }

  //free memory
  free(PK);
  free(PKe);
  free(PKc);
  free(Ppow);
  free(PKpow);
  free(Ptotpow);
  free(priorst);
  for(i=0; i<nstates; i++) free(piall[i]);
  free(piall);
  for(i = 0; i < n; i++) free(M[i]);
  free(M);
  
  printf("end likelihood computation process %i/%i\n", my_id+1, num_procs);
  

  if(my_id!=root_process){
    printf("Sending data (proc %d)... ",my_id);
    double MSGsend[2+avg_rows_per_process];
    MSGsend[0] = my_id*avg_rows_per_process + nstep%num_procs ;
    MSGsend[1] = (my_id+1)*avg_rows_per_process + nstep%num_procs ;
    for(k=MSGsend[0],j=0; k<MSGsend[1]; k++,j++){
      MSGsend[2 + j] = Lik[k];
    }
    ierr = MPI_Send(MSGsend, 2+avg_rows_per_process, MPI_DOUBLE, root_process, return_data_tag, MPI_COMM_WORLD); 
    printf("done\n");
  }else{
    printf("Gathering data from %d proc... ",num_procs-1);
    for(an_id = 1; an_id < num_procs; an_id++) {
      double MSGrec[2+avg_rows_per_process];
      ierr = MPI_Recv(MSGrec, 2+avg_rows_per_process, MPI_DOUBLE, MPI_ANY_SOURCE, return_data_tag, MPI_COMM_WORLD, &status);
      for(k=MSGrec[0],j=0; k<MSGrec[1]; k++,j++){
	Lik[k] = MSGrec[2 + j];
      }
      printf("proc %d done\n",an_id);
    }
    //print results
    printf("Writing on file %s... ",fout);
    FILE * fe;
    fe=fopen(fout,"wb");  
    for(i=0; i<nstep; i++){
      fprintf(fe,"%.20lf\t",log(Lik[i]) ); 
    }
    fclose(fe);
    printf("done\n");
    time (&end);//records ending time
    dif = difftime (end,start);
    printf ("Finished. It took  %.2lf min\n", dif/60.0 );
  }
  //end MPI
  ierr = MPI_Finalize();
  
  return 0;
}
