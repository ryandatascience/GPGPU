#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include <semaphore.h>
#include <math.h>
#include <sys/time.h>
#define NUM_THREAD 1024
#define Max_VERTEX_NUM 1160000

int counta[NUM_THREAD];
long long countb[NUM_THREAD];
int nodenum[NUM_THREAD];

int *atom;  // store the node


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



VexBox adjmulist[Max_VERTEX_NUM];
VexBox adjmulist1[Max_VERTEX_NUM];
int vernum[Max_VERTEX_NUM];
int vernum1[Max_VERTEX_NUM];
int numindex[Max_VERTEX_NUM];

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
/*
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
*/


void createlist()
{	
	
	char const *fil;
	
	//return;
        fil="com-youtube.txt"; // data 
        
	//printf("%d\n",thread_num);
	int i;
	unsigned a,b;
	//VexBox adjmulist[Max_VERTEX_NUM];
	FILE *fl=fopen(fil,"r");  
	
	// read data 
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
			
	}
	
	//return 0;
	
}

 void *compute(void *pthread_num)
{
	int thread_num=*(int*)pthread_num;
	int i,k,m,x;
	//int i=blockIdx.x*blockDim.x +threadIdx.x;
	int numthread=(Max_VERTEX_NUM+NUM_THREAD-1)/NUM_THREAD;
	int first=min1(thread_num*numthread,Max_VERTEX_NUM);
	int last=min1((thread_num+1)*numthread,Max_VERTEX_NUM);
	if(first<Max_VERTEX_NUM)
	{
	for( i=first;i<last;i++)
	{
	// compute the number of closed triangles
	if(vernum[i]!=0)
	{   int j=numindex[i];
	    for( k=0;k<vernum[i];k++)
	    {	int a=atom[j+k];
		for( m=k+1;m<vernum[i];m++)
		{
		   int b=atom[j+m];
		   int q=min1(a,b);
		   int p=max1(a,b);
		   if(vernum[q]!=0){
		   for( x=0;x<vernum[q];x++)
		   {
		      if(p==atom[numindex[q]+x])
			{
			 counta[thread_num]++;
			 break;    
			}
		   }}
		}
	    }
	}
	// compute the number of unclosed triangles
	if(vernum1[i]!=0)
	{
		nodenum[thread_num]++;
		countb[thread_num]=vernum1[i]*(vernum1[i]-1)/2+countb[thread_num];
	}

	}}
	//printf("thread num:%d the number of node: %d\n",thread_num,nodenum[thread_num]);
	//printf("thread num:%d the number of closed triangle: %d\n",thread_num,counta[thread_num]);
	//printf("thread num:%d the number of unclosed triangle: %lld\n",thread_num,countb[thread_num]);
}



int main()
{	// declare the time  function
	struct timeval start;
	struct timeval end;
	gettimeofday(&start,NULL);
	int i,j;
	int nodenum1=0;
	long long count=0;
	//int *d_vernum, *d_closed, *d_numindex;
	
	int number=Max_VERTEX_NUM;
	for( i=0;i<Max_VERTEX_NUM;i++)
        {
	  
          vernum[i]=0; // directed number 
          adjmulist[i].firstedge=NULL;
	  vernum1[i]=0;// undirected number
          adjmulist1[i].firstedge=NULL;
        }
	
	
	
	
	createlist();// read the data and store in linkedlist
	

	int c=0,d=0;
	for(i=0;i<Max_VERTEX_NUM;i++)
	{
	  //c=max1(vernum[i],c);
	  //d=max1(vernum1[i],d);
	  if(vernum1[i]!=0)
	  {nodenum1++;}
	  c=vernum[i]+c;  // compute all significance digit
	  d=vernum1[i]+d;  
	}
	
	
	atom=(int *)malloc(sizeof(int)*c);
	memset(atom,0,sizeof(int)*c);
	
	//printf("vernum: %d  vernum1: %d\n",c,d);
	int k=0; long long counting=0;
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
	   /* if(vernum1[i]!=0)
        {
                
                counting=vernum1[i]*(vernum1[i]-1)/2+counting;
        }*/
	    	
        }
	
	// create multithread
	int *thread2_id=(int *)malloc(NUM_THREAD*sizeof(int));
        pthread_t tt[NUM_THREAD];
        for(i=0;i<NUM_THREAD;i++)
        {
           counta[i]=0;
           countb[i]=0;
           nodenum[i]=0;
           thread2_id[i]=i;
        }
         for(i=0;i<NUM_THREAD;i++)
         {int ret2;
        
         ret2=pthread_create(&tt[i],NULL,compute,&thread2_id[i]);
         if(ret2!=0)
         {printf("Create 2set pthread error!\n");}
         }
        
         for(j=0;j<NUM_THREAD;j++)
         {
          pthread_join(tt[j],NULL);
         }
	

	// compute the total number from threads
	int count1=0, countnode=0;
 	for(i=0;i<NUM_THREAD;i++)
	{
	   count=countb[i]+count;
	   countnode=nodenum[i]+countnode;
	   count1=counta[i]+count1;
	}
	gettimeofday(&end,NULL);
	long total=end.tv_sec-start.tv_sec;
	int s=NUM_THREAD;
	long long cot=count-3*count1;
	printf("the number of thread: %d\n",s);
	printf("the number of node: %lld\n",nodenum1);
	//printf("the number of countnode: %d\n",countnode);
	printf("the number of closed triangle: %d\n",count1);
        printf("the number of unclosed triangles: %lld\n",cot);
	//printf("the number of unclosed counting triangles: %lld\n",counting);
	printf("the time: %lds\n",total);

	
	return 0;  	
}


