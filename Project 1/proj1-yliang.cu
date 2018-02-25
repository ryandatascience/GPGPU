/* ==================================================================
	Programmer: Yu Liang (yliang@usf.edu)
	The basic SDH algorithm implementation for 3D data
	To compile: nvcc SDH.c -o SDH in the rc machines
	Based on Project 1 code from Dr. Yicheng Tu from USF
            6/6/2017
   ==================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <cuda.h>


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
	int d_cnt;   /* need a long long type as the count might be huge */
} bucket;


bucket * histogram;/* list of all buckets in the histogram   */
bucket * d_histogram; /* list of all buckets in the histogram in device  */
bucket * cuda_histogram; /* list of all buckets in the histogram from GPU  */
bucket * diff_histogram; /*difference between two histogram*/

long long	PDH_acnt;	/* total number of data points            */
int num_buckets;		/* total number of buckets in the histogram */
double   PDH_res;		/* value of w                             */
atom * atom_list;		/* list of all data points                */
atom * d_atom_list; /* list of all data points in device            */

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
	printf("Running time for CPU version: %ld.%06ld\n", sec_diff, usec_diff);
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

/*kernel function*/
__global__
void SDHKernel( atom *x, bucket *y,long long n, double PDH_res)
    {
	int i = threadIdx.x + blockDim.x * blockIdx.x;
	int j, h_pos;
	double dist;

	for(j = i + 1; j < n; j++)
        {
        double x1 = x[i].x_pos;        double x2 = x[j].x_pos;
        double y1 = x[i].y_pos;        double y2 = x[j].y_pos;
        double z1 = x[i].z_pos;        double z2 = x[j].z_pos;

	dist = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));

	h_pos = (int) (dist / PDH_res);

	atomicAdd(&y[h_pos].d_cnt, 1);
	/*histogram[h_pos].d_cnt++;*/
        }
    }

int main(int argc, char **argv)
{
	int i;
    int j;
	PDH_acnt = atoi(argv[1]);
	PDH_res	 = atof(argv[2]);
//printf("args are %d and %f\n", PDH_acnt, PDH_res);

	num_buckets = (int)(BOX_SIZE * 1.732 / PDH_res) + 1;
	histogram = (bucket *)malloc(sizeof(bucket)*num_buckets);
	cuda_histogram=(bucket*)malloc(sizeof(bucket)*num_buckets);
	diff_histogram=(bucket*)malloc(sizeof(bucket)*num_buckets);
	atom_list = (atom *)malloc(sizeof(atom)*PDH_acnt);


	srand(1);
	/* generate data following a uniform distribution */
	for(i = 0;  i < PDH_acnt; i++) {
		atom_list[i].x_pos = ((double)(rand()) / RAND_MAX) * BOX_SIZE;
		atom_list[i].y_pos = ((double)(rand()) / RAND_MAX) * BOX_SIZE;
		atom_list[i].z_pos = ((double)(rand()) / RAND_MAX) * BOX_SIZE;
	}

	/*GPU calculation*/
	cudaMalloc(&d_atom_list, PDH_acnt*sizeof(atom) );
	cudaMalloc(&d_histogram, num_buckets*sizeof(bucket));


	cudaMemset(d_histogram, 0, num_buckets*sizeof(bucket)); /*memory reset*/

	cudaMemcpy(d_atom_list,atom_list,PDH_acnt*sizeof(atom),cudaMemcpyHostToDevice); /*copy data from host to device*/

	/*run kernel*/
	SDHKernel<<<ceil(PDH_acnt/256.0), 256>>>(d_atom_list,d_histogram,PDH_acnt,PDH_res);

	cudaMemcpy(cuda_histogram,d_histogram,num_buckets*sizeof(bucket),cudaMemcpyDeviceToHost);/*copy result from device to host*/

	/* start counting time */
	gettimeofday(&startTime, &Idunno);

	/* call CPU single thread version to compute the histogram */
	PDH_baseline();

	/* check the total running time */
	report_running_time();

    	for(j = 0; j < num_buckets; j++) {
		diff_histogram[j].d_cnt =histogram[j].d_cnt-cuda_histogram[j].d_cnt;
		}

	/* print out the histogram */
	output_histogram(histogram);
	output_histogram(cuda_histogram);
	output_histogram(diff_histogram);

	/*free cuda memory*/
	cudaFree(d_atom_list);
	cudaFree(d_histogram);

	return 0;
}

