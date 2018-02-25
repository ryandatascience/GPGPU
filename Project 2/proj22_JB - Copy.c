/* ==================================================================
	Programmer: Yicheng Tu (ytu@cse.usf.edu)
	The basic SDH algorithm implementation for 3D data
	To compile: nvcc SDH.c -o SDH in the rc machines
   ==================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>


#define BOX_SIZE	23000 /* size of the data box on one dimension            */

/* descriptors for single atom in the tree */
typedef struct atomdesc {
	double x_pos;
	double y_pos;
	double z_pos;
} atom;

typedef struct hist_entry{
	//float min;
	//float max;
	unsigned long long int d_cnt;   /* need a long long type as the count might be huge */
} bucket;


bucket * histogram, * d_histogram,*fin_histogram, *dif_histogram;		/* list of all buckets in the histogram   */
long long	PDH_acnt;	/* total number of data points            */
int num_buckets;		/* total number of buckets in the histogram */
double   PDH_res;		/* value of w                             */
atom * atom_list, *d_atom_list;		/* list of all data points                */

/* These are for an old way of tracking time */
struct timezone Idunno;	
struct timeval startTime, endTime;


/* 
	distance of two points in the atom_list 
*/
  double p2p_distance(int ind1, int ind2) {
	
	double x1 = atom_list[ind1].x_pos;
	double x2 = atom_list[ind2].x_pos;
	double y1 = atom_list[ind1].y_pos;
	double y2 = atom_list[ind2].y_pos;
	double z1 = atom_list[ind1].z_pos;
	double z2 = atom_list[ind2].z_pos;
		
	return sqrt((x1 - x2)*(x1-x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}

__device__ double distance(int ind1, int ind2, atom *atom_list) {

        double x1 = atom_list[ind1].x_pos;
        double x2 = atom_list[ind2].x_pos;
        double y1 = atom_list[ind1].y_pos;
        double y2 = atom_list[ind2].y_pos;
        double z1 = atom_list[ind1].z_pos;
        double z2 = atom_list[ind2].z_pos;

        return sqrt((x1 - x2)*(x1-x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}



/* 
	brute-force SDH solution in a single CPU thread 
*/
  int PDH_baseline() {
	int i, j, h_pos;
	double dist;
	
	for(i = 0; i < PDH_acnt; i++) {
		for(j = i+1; j < PDH_acnt; j++) {
			dist = p2p_distance(i,j);
			h_pos = (int) (dist / PDH_res);
			histogram[h_pos].d_cnt++;
		} 
	}
	return 0;
}

__global__ void compute(long long n, atom *x, bucket *y, double res)
{
	int i=blockIdx.x*blockDim.x +threadIdx.x;
	int j;
	if(i<n/2){
	for(j=i+1; j<(n+1)/2+i+1;j++)
	{
	double dist=distance(i,j,x);
	int h_pos=(int)(dist/res);
	atomicAdd(&y[h_pos].d_cnt,1);
	//y[h_pos].d_cnt++;
	}
	}

	if(i>n/2-1&&i<n-1)
	{
		for(j=i+1; j<n;j++)
        	{
        	double dist=distance(i,j,x);
        	int h_pos=(int)(dist/res);
        	atomicAdd(&y[h_pos].d_cnt,1);
        	//y[h_pos].d_cnt++;
        	}
		
		if(n%2==0)
		{  int k=n-2-i;
		for(j=n/2+1+k; j<n;j++)
                   {
                double dist=distance(k,j,x);
                int h_pos=(int)(dist/res);
                atomicAdd(&y[h_pos].d_cnt,1);
                //y[h_pos].d_cnt++;
                   }
		}
		else{
		   int k=n-2-i;
                   for(j=n/2+2+k; j<n;j++)
                   {
                   double dist=distance(k,j,x);
                   int h_pos=(int)(dist/res);
                   atomicAdd(&y[h_pos].d_cnt,1);
                   //y[h_pos].d_cnt++;
                   }
		}

	}
}


/* 
	set a checkpoint and show the (natural) running time in seconds 
*/
double report_running_time() {
	long sec_diff, usec_diff;
	gettimeofday(&endTime, &Idunno);
	sec_diff = endTime.tv_sec - startTime.tv_sec;
	usec_diff= endTime.tv_usec-startTime.tv_usec;
	if(usec_diff < 0) {
		sec_diff --;
		usec_diff += 1000000;
	}
	printf("Running time for cuda version: %ld.%06ld\n", sec_diff, usec_diff);
	return (double)(sec_diff*1.0 + usec_diff/1000000.0);
}


/* 
	print the counts in all buckets of the histogram 
*/
void output_histogram(bucket *histogram){
	int i; 
	long long total_cnt = 0;
	for(i=0; i< num_buckets; i++) {
		if(i%5 == 0) /* we print 5 buckets in a row */
			printf("\n%02d: ", i);
		printf("%15lld ", histogram[i].d_cnt);
		total_cnt += histogram[i].d_cnt;
	  	/* we also want to make sure the total distance count is correct */
		if(i == num_buckets - 1)	
			printf("\n T:%lld \n", total_cnt);
		else printf("| ");
	}
}


int main(int argc, char **argv)
{
	int i;
	int thr_num;
	PDH_acnt = atoi(argv[1]);
	PDH_res	 = atof(argv[2]);
	thr_num=atoi(argv[3]);
//printf("args are %d and %f\n", PDH_acnt, PDH_res);

	num_buckets = (int)(BOX_SIZE * 1.732 / PDH_res) + 1;
	histogram = (bucket *)malloc(sizeof(bucket)*num_buckets);
	dif_histogram=(bucket *)malloc(sizeof(bucket)*num_buckets);
	fin_histogram=(bucket*)malloc(sizeof(bucket)*num_buckets);
	atom_list = (atom *)malloc(sizeof(atom)*PDH_acnt);
	memset(fin_histogram,0,sizeof(bucket)*num_buckets);
		
	srand(1);
	/* generate data following a uniform distribution */
	for(i = 0;  i < PDH_acnt; i++) {
		atom_list[i].x_pos = ((double)(rand()) / RAND_MAX) * BOX_SIZE;
		atom_list[i].y_pos = ((double)(rand()) / RAND_MAX) * BOX_SIZE;
		atom_list[i].z_pos = ((double)(rand()) / RAND_MAX) * BOX_SIZE;
	}
	cudaMalloc(&d_atom_list,PDH_acnt*sizeof(atom) );
	cudaMalloc(&d_histogram,num_buckets*sizeof(bucket));
	cudaMemset(d_histogram,0,sizeof(bucket)*num_buckets);
	cudaMemcpy(d_atom_list,atom_list,PDH_acnt*sizeof(atom),cudaMemcpyHostToDevice);
	
	// measure time
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);
	
	compute<<<(PDH_acnt+thr_num-1)/thr_num,thr_num>>>(PDH_acnt,d_atom_list,d_histogram,PDH_res);
	
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime,start,stop);
	printf("Time to generate: %0.5f ms\n",elapsedTime);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	cudaMemcpy(fin_histogram,d_histogram,num_buckets*sizeof(bucket),cudaMemcpyDeviceToHost);
	
	/* start counting time */
	//gettimeofday(&startTime, &Idunno);
	
	/* call CPU single thread version to compute the histogram */
	//PDH_baseline();
	/*for(int k=0; k<num_buckets; k++)
	{
		dif_histogram[k].d_cnt=histogram[k].d_cnt-fin_histogram[k].d_cnt;
	}*/
	/* check the total running time */ 
	//report_running_time();
	
	/* print out the histogram */
	//printf("\n CPU: \n");
	//output_histogram(histogram);
	printf("\n GPU: \n");
	output_histogram(fin_histogram);
	//printf("\n difference:\n");
	//output_histogram(dif_histogram);
	cudaFree(d_atom_list);
	cudaFree(d_histogram);
	return 0;
}


