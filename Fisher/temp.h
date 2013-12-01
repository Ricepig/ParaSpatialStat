#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"
 
#define intsize sizeof(int)
#define floatsize sizeof(float)
#define A(x,y) A[x*N+y]
#define Q(x,y) Q[x*N+y]
#define a(x,y) a[x*N+y]
#define f(x) f[x]
float *A,*Q;
float *a,*f;
int N,v,m;
int p;
int myid;
MPI_Status status;
FILE *dataFile;
double starttime,endtime,time1;
 
void readData()
{
       int i,j;
       starttime = MPI_Wtime();
       dataFile = fopen("dataIn.txt","r");
       fscanf(dataFile,"%d",&N);
       A = (float *)malloc(floatsize*N*N);
       for(i = 0; i < N; i++)
       {
              for(j = 0; j < N; j++)
              {
                     fscanf(dataFile,"%f",&A(i,j));
              }
       }
       fclose(dataFile);
       printf("Input of file \"dataIn.txt\"\n");
       printf("%d\n",N);
       for(i = 0; i < N; i++)
       {
              for(j = 0; j < N; j++)
              {
                     printf("%f\t",A(i,j));
              }
              printf("\n");
       }
       Q = (float *)malloc(floatsize*N*N);
}
 
void printResult()
{
       int i,j;
       printf("\nOutput of Matrix Q\n");
       for(i = 0; i < N; i++)
       {
              for(j = 0; j < N; j++)
              {
                     printf("%f\t",Q(i,j));
              }
              printf("\n");
       }
       endtime = MPI_Wtime();
       printf("\n");
       printf("Whole running time    = %f seconds\n",endtime-starttime);
       printf("Distribute data time  = %f seconds\n",time1-starttime);
       printf("Parallel compute time = %f seconds\n",endtime-time1);
 
       dataFile=fopen("dataOut.txt", "w");
       fprintf(dataFile, "Input of file \"dataIn.txt\"\n");
       fprintf(dataFile, "%d\n",N);
       for(i=0;i<N;i++)
       {
              for(j=0; j<N; j++)
              {
                     fprintf(dataFile, "%f\t", A(i,j));
              }
              fprintf(dataFile, "\n");
        }
       fprintf(dataFile, "\nOutput of Matrix A's inversion\n");
       for(i=0; i<N; i++)
       {
              for(j=0; j<N; j++)
              {
                     fprintf(dataFile, "%f\t", Q(i,j));
              }
              fprintf(dataFile, "\n");
       }
        /* 0号进程将时间统计写入目标文件 */
       fprintf(dataFile, "\n");
       fprintf(dataFile, "Whole running time    = %f seconds\n",endtime-starttime);
       fprintf(dataFile, "Distribute data time  = %f seconds\n",time1-starttime);
       fprintf(dataFile, "Parallel compute time = %f seconds\n",endtime-time1);
       fprintf(dataFile, "Parallel Process number = %d\n", p);
       fclose(dataFile);
 
}
 
void broadcast(int i,int j,int v)
{
       int k;
       if(myid == j)
       {
              a(i,v) = (float)(1 / a(i,v));                               
              for(k = 0; k < N; k++)
              {
                     if(k != v)
                     {
                            a(i,k) = a(i,k) * a(i,v);
                     }
              }
 
              for(k = 0; k < N; k++)
              {
                     f(k) = a(i,k);
              }
 
              MPI_Bcast(&f(0),N,MPI_FLOAT,myid,MPI_COMM_WORLD);
       }
 
       else
       {
              MPI_Bcast(f,N,MPI_FLOAT,j,MPI_COMM_WORLD);
       }
}
 
void transform(int i,int j,int v)
{
    int k,w;
    if(myid != j)
       {
              for(k = 0; k < m; k++)
              {
                     for(w = 0; w < N; w++)
                     {
                            if(w != v)
                            {
                                   a(k,w) = a(k,w) - f(w) * a(k,v);
                            }
                     }
              }
 
              for(k = 0; k < m; k++)
              {
                     a(k,v) = -f(v) * a(k,v);
              }
       }
 
       if(myid == j)
       {
              for(k = 0; k < m; k++)
              {
                     if(k != i)
                     {
                            for(w = 0; w < N; w++)
                            {
                                   if(w != v)
                                   {
                                          a(k,w) = a(k,w) - f(w) * a(k,v);
                                   }
                            }
                     }
              }
 
              for(k = 0; k < m; k++)
              {
                     if(k != i)
                     {
                            a(k,v) = -f(v) * a(k,v);
                     }
              }
       }
}
 
int main2(int argc,char **argv)
{
       int i,j,k,w,group_size;
       MPI_Init(&argc,&argv);
       MPI_Comm_size(MPI_COMM_WORLD,&group_size);
       MPI_Comm_rank(MPI_COMM_WORLD,&myid);
       p = group_size;
 
       if(myid == 0)
       {
              readData();
       } 
       m = N/p; 
       if(myid == 0)
       {
              for(i = 1; i < p; i++)
              {
                     MPI_Send(&N,1,MPI_INT,i,i,MPI_COMM_WORLD);
                     MPI_Send(&m,1,MPI_INT,i,i,MPI_COMM_WORLD);
              }
       }
       if(myid != 0)
       {
              MPI_Recv(&N,1,MPI_INT,0,myid,MPI_COMM_WORLD,&status);
              MPI_Recv(&m,1,MPI_INT,0,myid,MPI_COMM_WORLD,&status);
       }
 
       a = (float *)malloc(floatsize*m*N);
       f = (float *)malloc(floatsize*N);
       if(a == NULL || f == NULL)
       {
              printf("Allocate space for a or f fail!");
       }
 
       if(myid == 0)
       {
              for(i = 0; i < m; i++)
              {
                     for(j = 0; j < N; j++)
                     {
                            a(i,j) = A(p*i,j);
                     }
              }
 
              for(k = 1; k < p; k++)
              {
                     for(i = 0; i < m; i++)
                     {
                         MPI_Send(&A(p*i+k,0),N,MPI_FLOAT,k,k,MPI_COMM_WORLD);
                     }
              }
       }
 
       if((myid > 0)&&(myid < p))
       {
              MPI_Recv(&a(0,0),m*N,MPI_FLOAT,0,myid,MPI_COMM_WORLD,&status);
       }
 
       if(myid < p)
       {
              for(i = 0; i < m; i++)
              {
                     for(j = 0; j < p; j++)
                     {
                            v = i * p + j;
                            broadcast(i,j,v);
                            transform(i,j,v);
                     }
              }
       }
      
       if((myid > 0)&&(myid < p))
       {
              MPI_Send(a,m*N,MPI_FLOAT,0,myid,MPI_COMM_WORLD);
       }
 
       if(myid == 0)
       {
              for(i = 0; i < m; i++)
              {
                     for(j = 0; j < N; j++)
                     {
                            Q(p*i,j)=a(i,j);
                     }
              }
 
              for(k = 1; k < p; k++)
              {
                     for(i = 0; i < m; i++)
                     {
                           MPI_Recv(&Q(p*i+k,0),N,MPI_FLOAT,k,k,MPI_COMM_WORLD,&status);
                     }
              }
              printResult();
              free(Q);
       }
 
       MPI_Finalize();
       if(myid < p)
       {
              free(a);
              free(f);
       }
       return(0);
}