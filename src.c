#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#define MAX 250
#define FILENAME "./inst/lcmste250.10"
#define TARGET_FILE "2501.txt"

const int count_in=(int)((float)MAX*(0.1));  //count of internal nodes
void initgraph(float *g, float totalCost[MAX], FILE *fp1, float x[MAX], float y[MAX]);
float computeTreeCost(int *T, float *g);
static int compare (const void * x, const void * y);
void cpu_time(double *t); 
void totxt(int *T,float*g);

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
  //printf("y:%f uy:%f\n",y[i],s.oy);
  if(!(x[i]>=s.ox && x[i]<=s.ox+s.ux)){return 0;}
  else if(!(y[i]>=s.oy && y[i]<=s.oy+s.uy)){return 0;}
  else return 1;
}

//bubbling a segment into two
void bubble(segment s[],int ri, int ci, int root[], float x[], float y[]){ //ci: child index ri: richest seg index
  //specify new boundaries
  int hor=0; //(hor==1) ? bubble vertically : b horizontally
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

void rep_graph(int*T,float*g,int inds[]){
  //reset T
  for(int i=0;i<MAX*MAX; i++){
    *(T+i)=0;
  }
  //reconnect internal nodes
  prims_mst(inds,count_in,g,T);

  //join leaves
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
    *(T+MAX*lf+nin)=*(T+MAX*nin+lf)=1;
  }  
}

float hill_climb(int inds[],float wt,int*T,float* g,int root[]){ 
  float bw=computeTreeCost(T,g);

  const int count_iter=MAX;
  int lasti=0;

  //k th iteration for randomized hill climb
  for(int k=0;k<count_iter;k++){
    int r_guide[count_in];
    int taken[count_in];
    for(int i=0;i<count_in;i++) taken[i]=0;

    //prepare randomization guide array r_guide[]
    for(int i=0;i<count_in;i++){
      int nin; //inds index
      while(taken[nin=(int)rand()%count_in]!=0);
      taken[nin]=1;
      r_guide[i]=nin;
    }

    float ciw=bw;  //current iteration initial w
    //run hillclimb on r_guide
    for(int i=0;i<count_in;i++){
      int index=r_guide[i];
      int bn=inds[index];
      int cn=inds[index];

      for(int j=0;j<MAX;j++){
        if(root[j]==root[cn] && cn!=j){
          inds[index]=j;
          rep_graph(T,g,inds);
          float cw=computeTreeCost(T,g);
          if(cw>bw){ //changes make no improvement => restore
            inds[index]=bn; 
            rep_graph(T,g,inds);
          }
          else{ //changes accepted
            bn=j;
            bw=computeTreeCost(T,g);
          }
        }
      }

    }

    if(bw<ciw)
      lasti=k+1;

  }

  printf("Last improvement at iteration number @ %d\n",lasti);

  return bw;
}

float algo(int* T, float* g, float* totalCost, float* x, float* y){
  segment s[count_in];
  int root[MAX]; 
  int inds[count_in]; //stands for internal nodes
  partition(s,root,x,y); //nodes distributed among these segments
  select_in(inds,s,root,x,y); //indexes of internal nodes stored in inds array
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
  wt=hill_climb(inds,wt,T,g,root);
  printf("After Hill_Climb: %f\n",wt);//test
  return wt;
}

int main(){
  srand(time(NULL));
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