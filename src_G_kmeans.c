// cd "d:\lcmst\" ; if ($?) { gcc src_G_kmeans.c -o src_G_kmeans } ; if ($?) { .\src_G_kmeans }
// python -u "d:\lcmst\kmeans_G.py"
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <conio.h>

#define rep(i,a,b) for(int i=a;i<b;i++)

#define MAX 250
#define FILENAME "./inst/lcmste250.15"
#define TARGET_FILE "2501.txt"

#define KMeansClusterCenters "km_centres.txt"
#define KMeansLabels "k_labels.txt"
#define CURR_GRAPH "curr_graph.txt"
#define CURR_INDS "curr_inds.txt"
#define KM_INDS "kmeans_inds.txt"
#define SGMT './segment.txt'
#define CNTRDS "./centroids.txt"
#define segment_no 0

const int count_in=(int)((float)MAX*(0.1));  //count of internal nodes
void initgraph(float *g, float totalCost[MAX], FILE *fp1, float x[MAX], float y[MAX]);
float computeTreeCost(int *T, float *g);
static int compare (const void * x, const void * y);
void cpu_time(double *t); 
void totxt(int *T,float*g);
void xytotxt(const char*,float x[],float y[],int n);
void nodestotxt(const char *s,int nodes[], float x[],float y[],int n);

typedef struct segment{
  int n; //no of elements in the array pointed to by ptr
  float ox,oy,ux,uy;
} segment;

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

//finds euclidean distance between points (x1,x2) and (y1,y2)
float distance(float x1,float y1, float x2, float y2){
  return sqrt( ((x1-x2)*(x1-x2)) + ((y1-y2)*(y1-y2)) );
}

//select internal nodes
void select_in(int inds[], segment s[], int root[], float x[], float y[]){
  //for(int i=0; i<count_in;i++){printf("%d ",s[i].n);}
  //calc avg of all nodes of a segment
  float ax[count_in]; //avg over x coord
  float ay[count_in]; //avg over y coord
  for(int i=0;i<count_in;i++){
    ax[i]=0; ay[i]=0;
  }
  for(int i=0;i<MAX;i++){
    int si=root[i];
    ax[si]+=x[i];
    ay[si]+=y[i];
  }
  
  for(int i=0;i<count_in;i++){
    ax[i]/=s[i].n; 
    ay[i]/=s[i].n;
  }
  //printing axy in txt
  printf("\nSegment Centroids to txt");
  xytotxt("centroids.txt",ax,ay,count_in);

  //printing segment in txt
  //int segment_no=0;
  float sx[s[segment_no].n],sy[s[segment_no].n];
  int p=0;
  rep(i,0,MAX){
    if(root[i]==segment_no){
      sx[p]=x[i];
      sy[p]=y[i]; p++;
    }
  }
  printf("\nNodes from segment %d to txt",segment_no);
  xytotxt("segment.txt",sx,sy,s[segment_no].n); 

  //find nearest nodes to each segment's central point (as determined via averaging)
  float da[count_in];
  for(int i=0; i<count_in; i++)
    da[i]=INT_MAX; //stores 
  for(int i=0;i<MAX;i++){
    float dtc=distance(ax[root[i]],ay[root[i]],x[i],y[i]); //dtc: distance to centermost point of the same segment
    if(dtc<da[root[i]]){
      da[root[i]]=dtc;
      inds[root[i]]=i;
    }
  }
}

//is node i eligible to be part of the segment s
int is_eligible(segment s, int i,float x[],float y[]){
  //is node i eligible to be part of the segment i.e. are the coordinates of node i within the boundaries of segment s
  
  // if(!(x[i]>=s.ox && x[i]<=s.ox+s.ux)){return 0;}
  // else if(!(y[i]>=s.oy && y[i]<=s.oy+s.uy)){return 0;}
  // else return 1;
  
  int c1=0,c2=0;
  if(x[i]>=s.ox && x[i]<=(s.ox+s.ux)){
    c1=1;
  }
  if(y[i]>=s.oy && y[i]<=(s.oy+s.uy)){
    c2=1;
  }
  if(c1==1 && c2==1){
    return 1;
  }
  else return 0;
}

//bubbling a segment into two
void bubble(segment s[],int ri, int ci, int root[], float x[], float y[]){ //ci: child index ri: richest seg index
  //specify new boundaries
  int hor=0; //(hor==1) ? bubble hor : b vert
  if(s[ri].ux<s[ri].uy){hor=1;}
  if(hor==0){
    s[ci].ox=s[ri].ox+s[ri].ux/2;
    s[ci].oy=s[ri].oy;
    s[ci].ux=s[ri].ux/2;
    s[ci].uy=s[ri].uy;
    s[ri].ux/=2;
  }
  else{
    s[ci].oy=s[ri].oy+s[ri].uy/2;
    s[ci].ox=s[ri].ox;
    s[ci].uy=s[ri].uy/2;
    s[ci].ux=s[ri].ux;
    s[ri].uy/=2;
  }
  //modify root and segment array (n) accordingly)
  for(int i=0; i<MAX; i++){
    if(is_eligible(s[ci],i,x,y)==1){
      //printf("d");//test
      s[root[i]].n--;
      s[ci].n++;
      root[i]=ci;
    }
  }
}

//partition the unit square into segments equal to number of internal nodes 
void partition(segment s[], int root[],float x[], float y[]){  
  for(int i=0 ;i<MAX; i++){root[i]=0;}
  s[0].n=MAX;
  s[0].ox=0;s[0].oy=0;s[0].ux=1;s[0].uy=1;
  for(int i=1;i<count_in;i++){
    s[i].n=0;
    s[i].ox=0;s[i].oy=0;s[i].ux=0;s[i].uy=0;
  }
  for(int ci=1; ci<count_in; ci++){
    //determine rich segment
    int ri; //richest segment index
    int max=INT_MIN;
    for(int i=0; i<count_in; i++){
      if(s[i].n>max){
        ri=i;
        max=s[i].n;
      }
    }
    //bubble from rich segment to current segment
    bubble(s,ri,ci,root,x,y);
  }
}

void inds_from_kmeans(int inds[], int root[], float* x, float* y){
  float ax[count_in], ay[count_in]; 
  //Store kmeans-cluster-centers in ax and ay
  FILE *f=fopen(KMeansClusterCenters,"r");
  rep(i,0,count_in){
    fscanf(f,"%f",ax+i);
    fscanf(f,"%f",ay+i);
  }
  fclose(f);
  //kmeans labels in root[]
  f=fopen(KMeansLabels,"r");
  rep(i,0,MAX){
    fscanf(f,"%d",root+i);
  }
  fclose(f);

  //find nearest nodes to each kmeans-cluster-centers
  //internal nodes in inds[]
  float da[count_in]; //smallest distances obtained from cluster centers
  for(int i=0; i<count_in; i++)
    da[i]=INT_MAX; //stores 
  for(int i=0;i<MAX;i++){
    float dtc=distance(ax[root[i]],ay[root[i]],x[i],y[i]); //dtc: distance to centermost point of the same segment
    if(dtc<da[root[i]]){
      da[root[i]]=dtc;
      inds[root[i]]=i;
    }
  }
}

float algo(int* T, float* g, float* totalCost, float* x, float* y){
  segment s[count_in];
  FILE *f;
  int root[MAX]; 
  int inds[count_in]; //stands for internal nodes
  partition(s,root,x,y); //nodes distributed among these segments
  printf("Partitioned");
  select_in(inds,s,root,x,y); //indexes of internal nodes stored in inds array
  
  //input coord of inds in txt file
  printf("\ninds to txt");
  nodestotxt(CURR_INDS,inds,x,y,count_in);
  //input coord of all nodes in txt file
  printf("\ngraph to txt");
  xytotxt(CURR_GRAPH,x,y,MAX);

  printf("\nPress any key to continue..");getch();
  
  inds_from_kmeans(inds,root,x,y);
  printf("\ninds updated\ninds to txt");

  nodestotxt(KM_INDS,inds,x,y,count_in);

  //implement Prims on inds and obtain an mst
  printf("\nApplying Prims");
  float wt=prims_mst(inds,count_in,g,T);
  
  printf("\nJoining nodes");
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
  printf("\nMST obtained successfully");
  printf("\nTree wt: %f",wt);
  printf("\nTree Cost: %f",computeTreeCost(T,g));

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

  printf("\n%f",computeTreeCost(T,g)); 
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

void xytotxt(const char *s,float x[],float y[],int n){
  FILE *f;
  f=fopen(s,"w");
  for(int i=0;i<n;i++){
    fprintf(f,"%f %f\n",x[i],y[i]);
  }
  fclose(f);
}

void nodestotxt(const char *s,int nodes[], float x[],float y[],int n){
  FILE *f;
  f=fopen(s,"w");
  for(int i=0;i<n;i++)
    fprintf(f,"%f %f\n",*(x+nodes[i]),*(y+nodes[i]));
  fclose(f);
}