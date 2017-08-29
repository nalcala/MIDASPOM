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

unsigned int simpij(unsigned int* piold, unsigned int* pinew, double e, double c, double K, double Ksource, double** M, unsigned int n){
  //extinction phase
  unsigned int k;
  double E = e/K;
  //double Es= e/Ksource;
  double pp;
  unsigned int pitmp[n];
  //printf("Old:\n");
  //for(k=0;k<n;k++) printf("%d ",piold[k]);
  //printf("\nExt (E=%.2lf):\n",E);
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
  double pCi; 
  unsigned int res=0;
  for(k=0; k<n; k++){
    //printf("segment %d\t",k);
    if(pitmp[k]==0){
      double s1 = 0;
      unsigned int l;
      for(l=0; l<n; l++){
	//printf("segment %d\t",l);
	if(l!=k) s1 += pitmp[l]*K*M[l][k];//*;
      }
      //printf("s1=%lf\t",s1);
      //if(Ksource>0) s1 += M[n][k]*Ksource;
      pCi = c*s1;
      if(pCi>1) pCi = 1;
      //printf("C=%lf",pCi);
      pp = (double)rand()/(double)(RAND_MAX);
      if(pp<pCi) pinew[k] = 1;
      else pinew[k] = 0;
      //printf(" end C");
    }else{
      //printf("no C %d\t",pinew[k]);
      pinew[k] = 1;
    }
    //res += pinew[k];
    //printf("-> %d \n",pinew[k]);
  }
  //printf("\n");
  return res;
}

//compPePc(Pe, Pc, pcur,ptmp, etmp, pCtmp, M, n, nsim, j);
void compPePc(double * Pee, double * Pcc, unsigned int ** p0,  unsigned int ** p1, double e, double * pC, double ** M, unsigned int n, unsigned int nsim){
  //compute probability transition from state piold to all states in piall
  unsigned int i,j;
  double etmp = e;
  if(etmp>1) etmp = 1;
  for(j=0;j<nsim;j++){// for all possible states at time point k
  unsigned int * pitmp = p1[j];
  for(i=0;i<nsim;i++){// for all possible states at time point k
    unsigned int * piold = p0[i];

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
      Pee[ i*nsim + j] = pow(etmp,s1)*pow(1-etmp,s2);
      Pcc[ j*nsim + i] = exp(pijc);
    }else{
      Pee[ i*nsim + j] = 0;
      Pcc[ j*nsim + i] = 0;
    }//ispos
  }//i
  }//j
}

double ylik(int ** Y, int ** J, unsigned int j,unsigned int n,double p, unsigned int * pi){
  double lres = 0;
  unsigned int i;
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
  char* Jname = "J.txt";
  double d    = 100;
  double win  = 0.01;
  unsigned int nstepe = 101;
  unsigned int nstepc = 101;
  double a    = 1.0/400; //dispersal distance
  double emin = 0.0;
  double emax = 1.0;
  double cmin = 0.0;
  double cmax = 1.0;
  double p    = 1; //proba of observation given presence
  double C    = 100000;//scaling constant
  unsigned int nsim        = 10000;
  
  int c;
  
  opterr = 0;
  while ((c = getopt (argc, argv, "m:p:d:i:J:o:s:l:h:u:v:f:n:")) != -1)
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
      case 'f':
        p = atof(optarg);
        break;
      case 'n':
        nsim = atoi(optarg);
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

  nstepe = (unsigned int) ((emax-emin)/win + 1);
  nstepc = (unsigned int) ((cmax-cmin)/win + 1);

  if(my_id == root_process) printf("Parameters for numerical approximation of the posterior density:\n\tWindow size=%lf, number of steps=%dx%d\n",win,nstepe,nstepc );
 
  //records the time of execution
  time_t start,end;         //beginning and ending times
  double dif;               //difference between beginning and ending
  
  //model parameters
  unsigned int i,j,k,l,tmp;

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
  
  unsigned int** piobs = malloc(tmax*sizeof(unsigned int*)); //piobs[tmax][n]; 
  for(i=0; i<tmax; i++) piobs[i] = malloc(n*sizeof(unsigned int));

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

  double* pendp = (double *) malloc(n*sizeof(double));
  
  unsigned int nvar=0; //number of variable segments
  for(j=0; j<n; j++){
    if(var[j]==1) nvar++;
    if(Y[0][j]>0){
      pendp[j] = 1;
    }else{
      pendp[j] = prioroc*pow(p,Y[0][j])*pow(1-p,J[0][j]-Y[0][j])/( prioroc*pow(p,Y[0][j])*pow(1-p,J[0][j]-Y[0][j]) + 1-prioroc );
    }
  }

  unsigned int nstates = pow(2,nvar); //number of possible states
  
  avg_rows_per_process = nstepe/num_procs;
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
  /*unsigned int** piall = malloc(nstates*sizeof(unsigned int*));
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
    }*/
       
  //number of possible states at each time step
  //unsigned int npstates[tmax]; 
  //unsigned int nnpstates = 0;
  //unsigned int** pstates = malloc(tmax*sizeof(unsigned int*)); //list of state id in all states list
  //unsigned int** simppstates = malloc(tmax*sizeof(unsigned int*)); //list of state id in states short list
  //double priorst = 1; //(float*) malloc( n*sizeof(float));;//prior probability of initial patch state
  //unsigned int s1;
  //unsigned int st1;
  //unsigned int nextid = 0;
 
  //printf("nextid=%d\n",nextid);
  
  /*unsigned int* all2short = malloc(nextid*sizeof(unsigned int)); //correspondance between table ID piall and short table
  for(i=0; i<tmax; i++){
    for(k=0; k<npstates[i]; k++){
      all2short[simppstates[i][k]] = pstates[i][k];
    }
    }*/

  /*for(i=0; i<tmax; i++) free(pstates[i]);
    free(pstates);*/
  
  //print input
  if(my_id == root_process){
    printf("Input occupancy data:\n");
    for(i=0; i<tmax; i++){
      printf("Year %d: ",i);
      for(j=0; j<n; j++) printf("%d/%d=%d\t",Y[i][j],J[i][j],piobs[i][j]);
      printf("\n");
    }

  }
  if(my_id == root_process) printf("Number of states to compute: %d\n",nstates);
  
  //unsigned int st = pi_to_id(piobs[0],n);

  for(i=0; i<tmax; i++) free(piobs[i]);
  free(piobs);
  
     
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
      Lik[i][j] = 0; //1
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
    endloop = avg_rows_per_process + nstepe%num_procs;
  }else{
    stloop  = my_id*avg_rows_per_process + nstepe%num_procs;
    endloop = (my_id+1)*avg_rows_per_process + nstepe%num_procs;
  }
  double *P  = malloc(nsim*nsim*sizeof(double));
  //double *D  = malloc(nstates*tmax*sizeof(double)); 
  double *pCtmp = malloc(n*sizeof(double)); //proba of colonization
  
  unsigned int* pcur = (unsigned int *) malloc(sizeof(unsigned int)*n); 
  unsigned int* ptmp = (unsigned int *) malloc(sizeof(unsigned int)*n);
  /*for(j=0; j<nsim;j++){
    pcur[j] = (unsigned int *) malloc(sizeof(unsigned int)*n); 
    ptmp[j] = (unsigned int *) malloc(sizeof(unsigned int)*n);
    }*/
  printf("Create occupancy tables DONE\n");
  double Liktmp;
  int nlik;
  
  for(ie=stloop; ie<endloop; ie++){
    etmp = ellik[ie];
    for(ic=0; ic<nstepc; ic++){
      //printf("%d %d\n",ie,ic);
      ctmp = cllik[ic];

      // t=0
      //double *Pold = malloc(nsim*sizeof(double));
      //for(k=0;k<nsim;k++) Pold[k] = C*ylik(Y,J,0,n,p,pcur[k]);
      //printf("Computed Pold:\n");
      /*for(k=0;k<nsim;k++){
	printf("Pr(");
	for(l=0;l<n;l++) printf("%d\t",pcur[k][l]);
	printf(")=%lf\n",ylik(Y,J,0,n,p,pcur[k]));
      }
      printf("Done\n");*/
      //double *Pe = calloc(nsim*nsim,sizeof(double));
      //double *Pc = calloc(nsim*nsim,sizeof(double));
      nlik = 0;
      for(k=0;k<nsim;k++){
	Liktmp = C*ylik(Y,J,0,n,p,pcur);
	
	for(j=0;j<n;j++){
	  double pr = (double)rand()/(double)(RAND_MAX);
	  if( pr<=pendp[j]){
	    pcur[j] = 1;
	  }else{
	    pcur[j] = 0;
	  }
	}
      
	for(tmp=0;tmp<tmax-1;tmp++){
	//printf("t=%d\n",tmp);
	//for(k=0;k<nsim;k++) jtmp = simpij(pcur[k],ptmp[k], etmp, ctmp, 1,  0, M, n);//compute ptmp values
	//printf("Sim done\n");
	//unsigned int snp;
	/*for(j=0;j<nsim;j++){
	  for(k=0; k<n; k++){//for all patches
	    double s1 = 0;
	    for(l=0; l<n; l++){
	      if(l!=k) s1 += M[l][k]*pcur[j][l];
	    }//l
	    pCtmp[k] = ctmp*s1;
	    if(pCtmp[k]>1) pCtmp[k] = 1;
	  }//k
	  //compPePc(Pe, Pc, pcur,ptmp, etmp, pCtmp, M, n, nsim, j);

	  }//j*/
	//cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nsim, nsim, nsim, 1.0, Pe, nsim, Pc, nsim, 0, P, nsim);
	//free(Pe);
	//free(Pc);

	//printf("end copy\n");
	  //printf("sim %d\n",k);
	  unsigned int jtmp = simpij(pcur,ptmp, etmp, ctmp, 1,  0, M, n);
	  //printf("occ done\n");
	  Liktmp *= C*ylik(Y,J,tmp+1,n,p,ptmp);// extinct
	  //printf("Lik done\n");
	  memcpy(&pcur[0], &ptmp[0], sizeof(unsigned int)*n);
	}//tmp
	if(Liktmp>0) nlik ++;
	
	Lik[ie][ic] += Liktmp;//*priorst; //exp( log(Pold[k*nstates + l]) + log(C));
	//}
      
	//for(j=1; j<tmax; j++){//for each year
	//multiply by diagonal matrix indicating proba of observations
	/*double *Pcur = calloc(nsim*nsim,sizeof(double));   
	for(k=0; k<nsim; k++) Pcur[k*nsim + k] = C*ylik(Y,J,j,n,p,ptmp[k]);
	double *Ptmp2 = malloc(nsim*nsim*sizeof(double));
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nsim, nsim, nsim, C, P, nsim, Pcur, nsim, 0, Ptmp2, nsim);
	
	double *Ptmp = malloc(nsim*nsim*sizeof(double));
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, nsim, nsim, nsim, C, Pold, nsim, Ptmp2, nsim, 0, Ptmp, nsim);
	free(Ptmp2);
	//printf("Pold*Pcur done\n");
	free(Pcur);
	memcpy(&Pold[0], &Ptmp[0], sizeof(double)*nsim*nsim);
	free(Ptmp);
	snp +=nsim;*/
      }//k*/

      /*for(k=0;k<nsim;k++){
	for(l=0;l<nsim;l++){
	  Lik[ie][ic] += ;//exp( log(Pold[k*nsim + l]) + log(C));
	}//l
	}//k*/
      
      //free(Pold);
      //printf("LLik[%d][%d]=%lf (%d)\t",ie,ic,log(Lik[ie][ic]),nlik); //=> always 0?
      Lik[ie][ic] = log(Lik[ie][ic]); //compute log likelihood to simplify display of output data
    }//ic
    if(my_id==root_process) printf("%.2f%% done\n",((float)ie+1-stloop)*100.0/(endloop-stloop));
  }//ie
  
  //free memory
  //free(priorst);
  /*for(j=0;j<nstates;j++) free(piall[j]);
    free(piall);*/
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
  }
  free(J);
  free(Y);
  //free(D);

  printf("end likelihood computation process %i/%i\n", my_id+1, num_procs);
  //end root loop

  //start sending data
  if(my_id!=root_process){
    printf("Sending data (proc %d)... ",my_id);
    double MSGsend[2+avg_rows_per_process*nstepc];
    MSGsend[0] = my_id*avg_rows_per_process + nstepe%num_procs ;
    MSGsend[1] = (my_id+1)*avg_rows_per_process + nstepe%num_procs ;
    for(k=MSGsend[0],j=0; k<MSGsend[1]; k++,j++){
      for(l=0; l<nstepc; l++){
	MSGsend[2 + j*nstepc + l] = Lik[k][l];
      }
    }
    ierr = MPI_Send(MSGsend, 2+avg_rows_per_process*nstepc, MPI_DOUBLE, root_process, return_data_tag, MPI_COMM_WORLD); 
    printf("done\n");
  }else{
    printf("Gathering data from %d proc... ",num_procs-1);
    for(an_id = 1; an_id < num_procs; an_id++) {//for each proc
      double MSGrec[2+avg_rows_per_process*nstepc];
      ierr = MPI_Recv(MSGrec, 2+avg_rows_per_process*nstepc, MPI_DOUBLE, MPI_ANY_SOURCE, return_data_tag, MPI_COMM_WORLD, &status); //receive data
      for(k=MSGrec[0],j=0; k<MSGrec[1]; k++,j++){
	for(l=0; l<nstepc; l++){
	  Lik[k][l] = MSGrec[2 + j*nstepc + l]; //store data
	}
      }
    }
    printf("done\n");

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
  }
  //end MPI
  ierr = MPI_Finalize();
  
  return 0;
}
