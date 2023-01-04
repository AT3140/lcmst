//commit 1
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#define MAX 50
#define FILENAME "./inst/lcmste50.14"
#define TARGET_FILE "501.txt"

const int count_in=(int)((float)MAX*(0.1));  //count of internal nodes
void initgraph(float *g, float totalCost[MAX], FILE *fp1, float x[MAX], float y[MAX]);
float computeTreeCost(int *T, float *g);
static int compare (const void * x, const void * y);
void cpu_time(double *t); 
void totxt(int *T,float*g);

int minKey(float key[], int mstSet[], int nodes[])
{
    // Initialize min value
    float min = INT_MAX;
    int min_index;
    for (int vi=0; vi<count_in;vi++){
      int v=nodes[vi];;
      if (mstSet[v] == 0 && key[vi] <= min){          
        min=key[vi]; min_index=vi;   
      }
    }

    //return nodes[min_index];
    return min_index;
}

float prims_mst(int nodes[], const int size, float g[], int* T){
  float key[size];
  int pn[size];
  for(int i=0;i<size;i++){
    key[i]=INT_MAX;
    pn[i]=INT_MAX;
  }
  key[0]=0; 
  pn[0]=0;
  int mstset[MAX];
  for(int i=0;i<MAX;i++){mstset[i]=0;}

  for(int count=0;count<size-1;count++){
    //int u=minKey(key,mstset,nodes);
    int ui=minKey(key,mstset,nodes);
    int u=nodes[ui];
    mstset[u]=1; 
    for(int vi=0;vi<size;vi++){
      int v=nodes[vi];
      if(g[MAX*u+v]<key[vi] && mstset[v]==0){
        key[vi]=g[MAX*u+v];
        pn[vi]=ui;
      }
    }
  } 
  
  float wt=0;
  for(int i=0;i<size;i++){
    wt+=key[i];
    *(T+MAX*nodes[i]+nodes[pn[i]])=*(T+MAX*nodes[pn[i]]+nodes[i])=1;
  }
  return wt;
}

float tc(int ni,int count_nb,int *taken_in,float* g){
  float tc=0;
  int taken[MAX]={0};
  taken[ni]=1;

  if(taken_in){
    for(int i=0;i<MAX;i++){
      taken[i]=taken_in[i];
    }
    taken_in[ni]=1;
  }

  for(int i=0;i<count_nb;i++){
    int min=INT_MAX;
    int mi;
    for(int j=0;j<MAX;j++){
      if(taken[j]==0 && *(g+MAX*ni+j)<min){
        min=*(g+MAX*ni+j);
        mi=j;
      }
    }
    tc+=*(g+MAX*ni+mi);
    taken[mi]=1; //node mi is a neighbour of ni
    if(taken_in)
      taken_in[mi]=1;
  }
  return tc;
}

void select_in(int inds[],float* g){
  float costs[MAX]; //costs of node i to count_nb nearest nodes
  int count_nb=(MAX-count_in)/count_in;
  for(int i=0;i<MAX;i++){
    costs[i]=tc(i,count_nb,NULL,g);
  }
  int taken_in[MAX]={0};
  //fill inds using costs
  for(int j=0;j<count_in;j++){
    int min=INT_MAX;
    int mon;
    for(int i=0;i<MAX;i++){
     if(taken_in[i]==0 && costs[i]<min){
      min=costs[i];
      mon=i;
     }
    }
    costs[mon]=tc(mon,count_nb,taken_in,g);
    inds[j]=mon;
  }
}

//finds euclidean distance between points (x1,x2) and (y1,y2)
float distance(float x1,float y1, float x2, float y2){
  return sqrt( ((x1-x2)*(x1-x2)) + ((y1-y2)*(y1-y2)) );
}

float algo(int* T, float* g, float* totalCost, float* x, float* y){
  int root[MAX]; 
  int inds[count_in]; //stands for internal nodes
  select_in(inds,g);
  //implement Prims on inds and obtain an mst
  float wt=prims_mst(inds,count_in,g,T);
  
  for(int lf=0;lf<MAX;lf++){//for all leaves
    float min=INT_MAX;
    int nin;//nearest internal node
    for(int i=0; i<count_in; i++){
      int in=inds[i];
      if(*(g+MAX*lf+in)<min){
        nin=in;
        min=*(g+MAX*lf+in);
      }
    }
    wt+=*(g+MAX*lf+nin);
    *(T+MAX*lf+nin)=*(T+MAX*nin+lf)=1;
  }
  //printf("Original wt: %f\n",wt);
  return wt;
}

int main(){
  float x[MAX],y[MAX], totalCost[MAX];
  float g[MAX*MAX];
  int T[MAX*MAX]; 
  FILE *fp1;  
  double *t=(double*)malloc(sizeof(double));
  cpu_time(NULL);
  initgraph(g,totalCost,fp1,x,y);

  printf("n=%d and l=%d\n",MAX,MAX-count_in);

  algo(T,g,totalCost,x,y);

  printf("%f",computeTreeCost(T,g)); 
  cpu_time(t);

  printf("\nTime Taken (ms): %lf",*t);

  totxt(T,g);

  return 0;
}

void totxt(int* T, float* g){
  FILE *fp;
  fp=fopen(TARGET_FILE,"w");
  for(int i=0;i<MAX;i++){
    for(int j=0;j<MAX;j++){
      if(*(T+MAX*i+j)!=0.)
        fprintf(fp,"%f ",*(g+MAX*i+j));
      else fprintf(fp,"%f ",0.);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}

void initgraph(float *g, float totalCost[MAX], FILE *fp1, float x[MAX], float y[MAX]) {
  float val, nodes;
  int i,j,k;

//    FILE *fp1;	
  if ((fp1 = fopen(FILENAME,"r")) == NULL) {
    printf("\n%s not found\n",FILENAME);
    exit(0);
  }
  for (i=0; i<MAX; totalCost[i++]=0.0);

// Euclidean Instances
///*
  //fscanf(fp1,"%f",&nodes); //uncomment this for lcmstr files
  for (i=0;i<MAX;i++)
    fscanf(fp1,"%f %f",&x[i], &y[i]);
  fclose(fp1);

  for (i=0;i<MAX;i++)
    for (j=i;j<MAX;j++) {
      if (i==j)
        *(g+MAX*i+j) = *(g+MAX*j+i)=0.0;
      else
        *(g+MAX*i+j) = *(g+MAX*j+i) = sqrt( ((x[j]-x[i])*(x[j]-x[i])) + ((y[j]-y[i])*(y[j]-y[i])) );
    } // Each g(i,j) is assigned the euclidean distance between nodes i and j
//*/

// Random Instances
/*
     fscanf(fp1,"%d",&nodes);

     for (i=0;i<MAX;i++) {
         for (j=0;j<MAX;j++) {

             fscanf(fp1,"%f",&val);
             *(g+MAX*i+j) = val;
         }
     }
*/

  for (i=0;i<MAX;i++)
    for (j=0;j<MAX;j++) {
      totalCost[i] += *(g+MAX*i+j);
    } // total cost of joining all other nodes to each node is stored here

    /* 
	for (k=0;k<MAX;k++) printf("totalCost[%d]=%f\n",k,totalCost[k]);
	*/
}

float computeTreeCost(int *T, float *g) {
  int i,j;
  float cost=0.0;

  for (i=0; i<MAX; i++)			// if there is an edge between two nodes
    for (j=i; j<MAX; j++)  		// then add the cost of that edge
      if (*(T+MAX*i+j) == 1)
        cost = cost + *(g+MAX*i+j);

  return cost;
}

static int compare (const void * x, const void * y) {

  if (*(float *)(x)<*(float *)(y)) return -1;
  if (*(float *)(x)==*(float *)(y)) return 0;
  if (*(float *)(x)>*(float *)(y)) return 1;

}

//cpu_time(NULL) to start time
//cpu_time(&double_var) to stop and record
void cpu_time(double *t) //milliseconds
{
  // static double start,stop;
  static clock_t start,stop;
  if (t == NULL) {
    //start=time(NULL);
    start=clock();
  } else {
    //stop=time(NULL);
    stop=clock();
    //printf("\n%f\n%f\n",(double)stop,(double)start);
    *t = ((double) (stop - start) / (CLOCKS_PER_SEC/1000));
  }
}
// void cpu_time(double *t) //microseconds
// {
//   static struct timeval start,stop;
//   if (t == NULL) {
//     gettimeofday(&start,NULL);
//   } else {
//     gettimeofday(&stop,NULL);
//     *t = (double) ((stop.tv_sec - start.tv_sec)*(int)1e6 + (stop.tv_usec - start.tv_usec));
//   }
// }

// void cpu_time(double *t) //nanoseconds
// {
//   static struct timespec start,stop;
//   if (t == NULL) {
//     clock_gettime(CLOCK_MONOTONIC, &start);
//     ios_base::sync_with_stdio(false); 
//   } else {
//     clock_gettime(CLOCK_MONOTONIC, &stop);
//     *t = (double) ((stop.tv_sec - start.tv_sec)*(int)1e9 + (stop.tv_nsec - start.tv_nsec)*1e-9);
//   }
// }