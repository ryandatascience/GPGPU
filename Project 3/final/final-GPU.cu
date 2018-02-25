#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include <semaphore.h>
#include <math.h>
#include <sys/time.h>

#define Max_VERTEX_NUM 1160000
int counta[10];
int countb[10];
int nodenuma[10];

int *atom,*d_atom;
//typedef enum{unvisited, visited}VisitIf;

typedef struct EBox
{
        
        int vex; //two vertexs connect to this
        struct EBox *next;
        
}EBox;

typedef struct VexBox
{

        EBox *firstedge;

}VexBox;

typedef struct vector
{
	int *nod;
	int *nod1;
}vector; 


vector *node, *d_node;
VexBox adjmulist[Max_VERTEX_NUM];
VexBox adjmulist1[Max_VERTEX_NUM];
int vernum[Max_VERTEX_NUM];
int vernum1[Max_VERTEX_NUM];
//void counttriangle(int[],int[][Max_VERTEX_NUM],int);
int max1(int a,int b)
{
	if(a>b) return a;
	else return b;
}

int min1(int a,int b)
{
	if(a<b) return a;
	else return b;
}

__device__ int max2(int a,int b)
{
        if(a>b) return a;
        else return b;
}

__device__ int min2(int a,int b)
{
        if(a<b) return a;
        else return b;
}



void createlist()
{	
	
	char const *fil;
	
        fil="com-youtube.txt";
        
	
	//printf("%d\n",thread_num);
	int i;
	unsigned a,b;
	//VexBox adjmulist[Max_VERTEX_NUM];
	FILE *fl=fopen(fil,"r");
	//fscanf(fl,"%d %d",&a,&b);
	//for(i=0;i<1000;i++)
	for( i=0;!feof(fl);i++)
	{	
		fscanf(fl,"%d %d",&a,&b);
		int times=0;
		EBox *search=adjmulist[min1(a,b)].firstedge;
		while(search!=NULL)
		{
		   if(search->vex==max1(a,b))
			{times=1; break;}
			search=search->next;
		}
	    if(times==0)
		{
		EBox *vertex=(EBox*)malloc(sizeof(EBox));
		vertex->vex=max1(a,b);
		vertex->next=adjmulist[min1(a,b)].firstedge;
		adjmulist[min1(a,b)].firstedge=vertex;
		vernum[min1(a,b)]++;
		
		EBox *vertex2=(EBox*)malloc(sizeof(EBox));
                vertex2->vex=max1(a,b);
                vertex2->next=adjmulist1[min1(a,b)].firstedge;
                adjmulist1[min1(a,b)].firstedge=vertex2;
                vernum1[min1(a,b)]++;

                EBox *vertex1=(EBox*)malloc(sizeof(EBox));
                vertex1->vex=min1(a,b);
                vertex1->next=adjmulist1[max1(a,b)].firstedge;
                adjmulist1[max1(a,b)].firstedge=vertex1;
                vernum1[max1(a,b)]++;
		
		}
		//fscanf(fl,"%d %d",&a,&b);	
	}
	
	
	
}

__global__ void compute(int *atom, int *num, int *closed, int *index, int total)
{
	int i=blockIdx.x*blockDim.x +threadIdx.x;
	//int j;
	if(i<total)
	{
	if(num[i]!=0)
	{   int j=index[i];
	    for(int k=0;k<num[i];k++)
	    {	int a=atom[j+k];
		for(int m=k+1;m<num[i];m++)
		{
		   int b=atom[j+m];
		   int q=min2(a,b);
		   int p=max2(a,b);
		   if(num[q]!=0){
		   for(int x=0;x<num[q];x++)
		   {
		      if(p==atom[index[q]+x])
			{
			 closed[i]++;    
			}
		   }}
		}
	    }
	}
	}
}



int main()
{	struct timeval start;
	struct timeval end;
	gettimeofday(&start,NULL);
	int i;
	long long count=0,count1=0,nodenum=0;
	int *d_vernum, *d_closed, *d_numindex;
	int closed[Max_VERTEX_NUM], numindex[Max_VERTEX_NUM];
	int number=Max_VERTEX_NUM;
	for( i=0;i<Max_VERTEX_NUM;i++)
        {
	  closed[i]=0;
          vernum[i]=0;
          adjmulist[i].firstedge=NULL;
	  vernum1[i]=0;
          adjmulist1[i].firstedge=NULL;
        }
	
	
	createlist();
	// cuda
	/*
	node=(vector *)malloc(sizeof(vector)*Max_VERTEX_NUM);
	for(i=0;i<Max_VERTEX_NUM;i++)
	{
	     if(vernum[i]!=0)
	    {
	
		int k=0;
	     node[i].nod=(int *)malloc(sizeof(int)*vernum[i]);
	     EBox *search=adjmulist[i].firstedge;
	     while(search!=NULL)
		{
		  node[i].nod[k]=search->vex;
		  search=search->next;
		  k++;
		}
	    } 
	}
	*/
	int c=0,d=0;
	for(i=0;i<Max_VERTEX_NUM;i++)
	{
	  //c=max1(vernum[i],c);
	  //d=max1(vernum1[i],d);
	  if(vernum1[i]!=0)
	  {nodenum++;}
	  c=vernum[i]+c;
	  d=vernum1[i]+d;  
	}
	//printf("vernum: %d  vernum1: %d\n",c,d);
	
	atom=(int *)malloc(sizeof(int)*c);
	memset(atom,0,sizeof(int)*c);
	
	//printf("vernum: %d  vernum1: %d\n",c,d);
	int k=0;
	for(i=0;i<Max_VERTEX_NUM;i++)
        {	
	    
             if(vernum[i]!=0)
            {

                //int k=i*c;
             numindex[i]=k;
             EBox *search=adjmulist[i].firstedge;
             while(search!=NULL)
                {
                  atom[k]=search->vex;
                  search=search->next;
                  k++;
                }
            }
	    	
        }
	//printf("vernum: %d  vernum[4]: %d\n",k,vernum[4]);

	cudaMalloc(&d_vernum,sizeof(int)*Max_VERTEX_NUM);
	cudaMemset(d_vernum,0,sizeof(int)*Max_VERTEX_NUM);
	
	cudaMalloc(&d_numindex,sizeof(int)*Max_VERTEX_NUM);
        cudaMemset(d_numindex,0,sizeof(int)*Max_VERTEX_NUM);

	cudaMalloc(&d_closed,sizeof(int)*Max_VERTEX_NUM);
        cudaMemset(d_closed,0,sizeof(int)*Max_VERTEX_NUM);

	cudaMalloc(&d_atom,sizeof(int)*c);
	cudaMemset(d_atom,0,sizeof(int)*c);

	cudaMemcpy(d_numindex,numindex,sizeof(int)*Max_VERTEX_NUM,cudaMemcpyHostToDevice);
	cudaMemcpy(d_atom,atom,sizeof(int)*c,cudaMemcpyHostToDevice);
	cudaMemcpy(d_vernum,vernum,sizeof(int)*Max_VERTEX_NUM,cudaMemcpyHostToDevice);

	//compute<<<(Max_VERTEX_NUM+255)/256,256>>>(d_atom, d_vernum,d_closed,d_numindex,number);
	compute<<<1,2>>>(d_atom, d_vernum,d_closed,d_numindex,number);
	cudaMemcpy(closed,d_closed,Max_VERTEX_NUM*sizeof(int),cudaMemcpyDeviceToHost);

	//printf("vernum: %d  vernum1: %d\n",c,d);
	
 	for(i=0;i<Max_VERTEX_NUM;i++)
	{
	   count=counta[i]+count;
	   //nodenum=nodenuma[i]+nodenum;
	   count1=closed[i]+count1;
	}
	gettimeofday(&end,NULL);
	long total=end.tv_sec-start.tv_sec;
	printf("the number of node: %lld\n",nodenum);
	printf("the number of closed triangle: %lld\n",count1);
        //printf("the number of unclosed triangles: %lld\n",count-3*count1);
	printf("the time: %lds\n",total);

	
	return 0;  	
}


