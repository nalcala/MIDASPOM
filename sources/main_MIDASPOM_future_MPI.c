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
    double *tmp = (double *) malloc(nSqr*sizeof(double)); /* temporary matrix */
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

void shuffle(int *array, int n){
  if(n > 1){
    int i;
    for(i=0; i<n-1; i++){
      int j = i + rand()/(RAND_MAX/(n-i)+1 );
      int tmp = array[j];
      array[j] = array[i];
      array[i] = tmp;
    }
  }
}

 
int simpij(int* piold, int* pinew, double e, double c, double K, double Ksource, double** M, int n){
  //extinction phase
  int k;
  double E = e/K;
  double Es= e/Ksource;
  double pp;
  int pitmp[n];
  /*printf("Old:\n");
  for(k=0;k<n;k++) printf("%d ",piold[k]);
  printf("\nExt (E=%.2lf):\n",E);*/
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
  //printf("\nColo (C=%.2lf):\n",c*K);
  //colonization phase
  double pCi[n]; 
  int res=0;
  for(k=0; k<n; k++){
    double s1 = 0;
    int l;
    for(l=0; l<n; l++){
      if(l!=k) s1 += M[l][k]*pitmp[l]*K;
    }
    s1 += M[n][k]*Ksource;
    pCi[k] = c*s1;
    if(pCi[k]>1) pCi[k] = 1;
    if(pitmp[k]==0){
      pp = (double)rand()/(double)(RAND_MAX);
      if(pp<pCi[k]) pinew[k] = 1;
      else pinew[k] = 0;
    }else{
      pinew[k] = 1;
    }
    res += pinew[k];
    //printf("%d ",pinew[k]);
  }
  //printf("\n");
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
    printf("------MIDASPOM, beta MPI version -------\n-> N. Alcala, E. M. Cole, N. A. Rosenberg  <-\n");
    
  /*----------------Parameters--------------------*/
  //records the time of execution
  time_t start,end;         //beginning and ending times
  double dif;               //difference between beginning and ending
  
  //model parameters
  int tfut  = 50; //number of time steps without disease 
  double KS = 0;
  double dS = 200;
  double KD = 1;
  int nsimul  = 10000;
  double a = 1.0/400; //dispersal distance
  double d = 200;
  float prioroc = 0.5; //prior occupancy for missing data the first year
  char* finame = "posterior.txt"; //input file with joint posterior distribution of e and c
  char* fname = "input.txt";
  char* fout  = "pext_future.txt";

  int i,j,k,l,tmp;
  unsigned int index;
  int c;
  
  opterr = 0;
  while ((c = getopt (argc, argv, "n:a:m:p:q:d:i:o:S:s:D:")) != -1)
    switch (c)
      {
      case 'n':
        nsimul = atoi(optarg);
        break;
      case 'a':
        tfut = atoi(optarg);
        break;
      case 'm':
        a = 1.0/atof(optarg);
        break;
      case 'p':
        prioroc = atof(optarg);
        break;
      case 'q':
        finame = optarg;
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
      case 'S':
        KS = atof(optarg);
        break;
      case 's':
        dS = atof(optarg);
        break;
      case 'D':
        KD = atof(optarg);
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
 
  avg_rows_per_process = nsimul/num_procs;

   if(my_id == root_process) printf("%d years in the future\n",tfut);
  
  //read observed data
  FILE * fo;

  int tmax=0;
  if(my_id == root_process) printf("Reading observations from file %s... ",fname);
  fo = fopen(fname,"rb");   //
  int n=1;
  while ( (c=fgetc(fo)) != EOF ) {
    if ( c == '\n' )
      tmax++;
    if(tmax==0){
      //printf("%c",c);
      if( (c==' ')|(c=='\t') ) n++;
    }
  }
  fclose(fo);
  int nstates = pow(2,n); //number of possible states
  
  int* pend = (int *) malloc(n*sizeof(int));
    
  fo = fopen(fname,"rb");   //
  int out;
  int var[n];
  for(j=0; j<n; j++) var[j]=0;
  for(i=0; i<tmax; i++){
    for(j=0; j<n; j++){
      out = fscanf(fo,"%d",&(pend[j]));
      if(pend[j]!=0){
	var[j]=1;
      }
    }
  }
  fclose(fo);


  if(my_id == root_process){
    printf("Last occupancy survey:\n");
    for(j=0; j<n; j++){
      printf("%d ",pend[j]);
    }
    printf("\n");
  }
  if(my_id == root_process) printf("\n done\n");

  if(my_id == root_process){
    printf("Number of habitat patches: %d\nNumber of sampled years: %d\n",n,tmax);
  }
  
  //read posterior distrib of e and c
  FILE * fi2;
  if(my_id == root_process) printf("Reading posterior distribution from file %s... ",finame);
  fi2 = fopen(finame,"rb");   //
  int necstep = 0;
  while ( (c=fgetc(fi2)) != EOF ) {
    if ( c == '\n' )
      break;
    if( (c==' ')|(c=='\t') ) necstep++;
  }
  fclose(fi2);
  if(my_id == root_process) printf("%dX%d posterior distribution\n",necstep,necstep);
  
  //int necstep = 101;
  double post[necstep][necstep];
  fi2 = fopen(finame,"rb");   //
  double totpost = 0;
  for(i=0; i<necstep; i++){
    for(j=0; j<necstep; j++){
      out = fscanf(fi2,"%lf",&(post[i][j]));
      //printf("%lf\t",post[i][j]);
      totpost += post[i][j];
    }
    //printf("\n");
  }
  fclose(fi2);

 
  // build migration matrix
  double **M = (double **) malloc((n+1)*sizeof(double*));
  for(i = 0; i < n+1; i++) M[i] = (double *) malloc(n*sizeof(double));
  for(i=0; i<n; i++){
    for(j=i; j<n; j++){
      if(i==j) M[i][j]=0.0;
      else{
	M[i][j]=exp(-a*(j-i)*d);
	M[j][i]=exp(-a*(j-i)*d);
      }
    }
  }
  for(j=0; j<n; j++) M[n][j]=exp( -a*(j+1)*dS ); // distance to source
  if(my_id == root_process){
    printf("Migration matrix:\n");
    for(i=0; i < n+1; i++){
      for(j=0; j<n; j++){
	printf("%.3f ",M[i][j]);
      }
      printf("\n");
    }
  }
  
  int npstates; 
  int s1=0; //number of segments with missing data
  int st1;
  int pmin = 0; //min number of occupied segments
  
  
  for(j=0; j<n; j++){
    if(pend[j]==1) pmin++;
    if(pend[j]==-1){
      s1++; 
    }
  }
  npstates = pow(2,s1); //number of possible states
  printf("npstates = %d\n",npstates);
  int** pstates = (int**) malloc(npstates*sizeof(int*)); 
  for(i=0; i<npstates; i++) pstates[i] = (int*) malloc(n*sizeof(int));
  
  double* priorst = (double*) malloc( npstates*sizeof(double)); //prior probability of each initial state
  for(k=0; k<npstates; k++){
    priorst[k]=1;
    for(l=0; l<n; l++){
      pstates[k][l]=0;
    }
  }

  s1 = 0;
  for(j=0; j<n; j++){
    if(pend[j]==-1) s1 ++;
    for(k=0; k<npstates; k++){
      if(pend[j]>-1){
	pstates[k][j] = pend[j]; // to change if some sites always at 0
      }else{
	st1 = npstates/pow(2,s1);
	pstates[k][j] = (k/st1%2);
	priorst[k] *= (k/st1%2)*prioroc + (1-k/st1%2)*(1-prioroc);
      }  
    }//for state k
  }//for patch j
 
  free(pend);

  if(my_id == root_process){
    printf("Last occupancy survey:\n");
    for(k=0; k<npstates; k++){
      printf("\t");
      for(j=0; j<n;j++) printf("%d ",pstates[k][j]);
      printf("; pr=%lf\n",priorst[k]);
    }//for time point i
  }

  free(priorst);

  //create output tables
  int* Lik = (int *) malloc( tfut*sizeof(int));
  
  for(tmp=0;tmp<tfut;tmp++){
    Lik[tmp] = 0;
  }
  //printf("end Lik\n");
  
  
  //main loop
  time (&start);//records the beginning time
  srand (time(NULL));
  
  printf("Starting parallel likelihood computation process %i/%i\n", my_id+1, num_procs);
  double LLtmp,Ltmp2;
  int p0[n];
  // start root loop
  int stloop,endloop,t;
  if(my_id==root_process){
    stloop  = 0; 
    endloop = avg_rows_per_process + nsimul%num_procs;
  }else{
    stloop  = my_id*avg_rows_per_process + nsimul%num_procs;
    endloop = (my_id+1)*avg_rows_per_process + nsimul%num_procs;
  }
  
   // simul for each set of parameter values
  int jtmp,ie,ic;
  double pec,etmp=0,ctmp=0;
  int* pcur = (int *) malloc(sizeof(int)*n); 
  int* ptmp = (int *) malloc(sizeof(int)*n); 
  
  for(i=stloop; i<endloop; i++){//nstep; ie++){
    // sample e and c from prior distrib
    pec = (necstep-1)*(necstep-1)*(double)rand()/(double)(RAND_MAX);//rand()%(necstep*necstep);
    double pcum = 0;
    for(ie=0;ie<necstep;ie++ ){
      for(ic=0;ic<necstep;ic++ ){
	double w = 1.0;
	if( (ie==0)||(ie==necstep-1) ) w*= 0.5;
	if( (ic==0)||(ic==necstep-1) ) w*= 0.5;
	pcum += w*post[ie][ic];
	if(pec<pcum){
	  etmp = ie*0.01;
	  ctmp = ic*0.01;
	  goto end_nested_loop;
	}
      }
    }
  end_nested_loop:
    1;
    int init = rand()%npstates; //randomly choose an initial state
    memcpy(&pcur[0], &pstates[init][0], sizeof(int)*n);
    //printf("end copy\n");
    for(tmp=0;tmp<tfut;tmp++){
      jtmp = simpij(pcur,ptmp, etmp, ctmp, KD,  KS, M, n); 
      if(jtmp==0) Lik[tmp] += 1; // extinct 
      memcpy(&pcur[0], &ptmp[0], sizeof(int)*n);
    }
  }
  
  //free memory
  for(i = 0; i < n+1; i++) free(M[i]);
  free(M);
  free(pcur);
  free(ptmp);
  for(j=0;j<npstates;j++) free(pstates[j]);
  free(pstates);
  printf("end likelihood computation process %i/%i\n", my_id+1, num_procs);
  
  if(my_id!=root_process){
    printf("Sending data (proc %d)... ",my_id);
    double MSGsend[tfut];
    for(tmp=0; tmp<tfut; tmp++){
      MSGsend[tmp] = Lik[tmp];
    }
    ierr = MPI_Send(MSGsend, tfut, MPI_DOUBLE, root_process, return_data_tag, MPI_COMM_WORLD); 
    printf("done\n");
  }else{
    printf("Gathering data from %d proc... ",num_procs-1);
    for(an_id = 1; an_id < num_procs; an_id++) {
      double MSGrec[tfut];
      ierr = MPI_Recv(MSGrec, tfut, MPI_DOUBLE, MPI_ANY_SOURCE, return_data_tag, MPI_COMM_WORLD, &status);
      for(tmp=0; tmp<tfut; tmp++){
	Lik[tmp] += MSGrec[tmp];
      }
      printf("done\n");
    }
    //print results
    printf("Writing on file %s... ",fout);
    FILE * fe;
    fe=fopen(fout,"wb");  
    for(tmp=0; tmp<tfut; tmp++){
      fprintf(fe,"%d\t",Lik[tmp]);
    }
    fclose(fe);
    printf("done\n");
    time (&end);//records ending time
    dif = difftime (end,start);
    printf ("Finished. It took  %.2lf min\n", dif/60.0 );
  }
  
  //end MPI
  ierr = MPI_Finalize();
  
  //free memory
  free(Lik);
  
  return 0;
}
