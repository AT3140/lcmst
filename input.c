
void initgraph(float *g, float totalCost[MAX], FILE *fp1, float x[MAX], float y[MAX]) {

      float val, nodes;

      int i,j,k;

//    FILE *fp1;
/*	
      if ((fp1 = fopen(FILENAME,"r")) == NULL) {
         printf("\n%s not found\n",FILENAME);
         exit(0);
      }
*/

      for (i=0; i<MAX; totalCost[i++]=0.0);

// Euclidean Instances
///*
      fscanf(fp1,"%f",&nodes);

      for (i=0;i<MAX;i++)
          fscanf(fp1,"%f %f",&x[i], &y[i]);

      //fclose(fp1);

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


 // Sorting using the inbuilt qsort routine uses compare ()
 // Syntax: qsort (&floatarrayname, MAX, sizeof(float), compare);
 static int compare (const void * x, const void * y) {

    if (*(float *)(x)<*(float *)(y)) return -1;
    if (*(float *)(x)==*(float *)(y)) return 0;
    if (*(float *)(x)>*(float *)(y)) return 1;

 }


/* ======================================================================
                                   cpu_time
   ====================================================================== */

/* The following routine is used for measuring the time used
 *   cpu_time(NULL)    - start of time measuring
 *   cpu_time(&t)      - returns the elapsed cpu-time in variable t (double).
 */

#include<sys/times.h>
void cpu_time(double *t)
{
  static struct tms start, stop;

  if (t == NULL) {
    times(&start);
  } else {
    times(&stop);
    *t = (double) (stop.tms_utime - start.tms_utime) / sysconf(_SC_CLK_TCK);
  }
}
