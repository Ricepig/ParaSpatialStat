#include <cstdlib>
#include "stdio.h"
#include "mpi.h"
#include "paramatrix.h"

void singlematrix(double *A, double *B, double *C, int rA, int cA, int cB)
{
    for(int row=0;row<rA;row++)
    {
        for(int col=0;col<cB;col++)
        {
            double sum = 0;
            for(int i = 0;i<cA;i++)
            {
                sum+= A[row*cA + i] * B[i*cB + col];
            }
            C[row*cB + col] = sum;
        }
    }
}

void paramatrix(double *A, double *B, double* C, int rA, int cA, int cB, int numprocs, int myid)
{
    int rB = cA;
    MPI_Status status;

    int line = rA/numprocs;//将数据分为(进程数)个块,主进程也要处理数据
    
    //缓存大小大于等于要处理的数据大小，大于时只需关注实际数据那部分
    double * buffer = (double*)malloc(sizeof(double)*cA*line);//数据分组大小
    double *ans = (double*)malloc(sizeof(double)*cB*line);//保存数据块结算的结果
    int * sizes = (int*)malloc(sizeof(int)*3);
    //主进程对矩阵赋初值，并将矩阵N广播到各进程,将矩阵M分组广播到各进程
    if (myid==0)
    {
        sizes[0] = rA;
        sizes[1] = cA;
        sizes[2] = cB;
        
        MPI_Bcast(sizes, 4, MPI_INT, myid, MPI_COMM_WORLD);
        
        //将矩阵N发送给其他从进程
        MPI_Bcast(B,cB*rB,MPI_FLOAT, myid,MPI_COMM_WORLD);
        //依次将M的各行发送给各从进程
        for (int m=0;m<numprocs-1;m++)
        {
            MPI_Send(A+m*cA*line,cA*line,MPI_FLOAT,m+1,1,MPI_COMM_WORLD);
        } 
        
        //计算M剩下的数据
        for (int i=(numprocs-1)*line;i<rA;i++)
        {
            for (int j=0;j<cB;j++)
            {
                double temp=0.0;
                for (int k=0;k<cA;k++)
                    temp += A[i*cA+k] * B[k*cB+j];
                C[i*cB+j]=temp;
            }
        }
        
        int row;
        //接收从进程计算的结果
        for (int k=0;k<numprocs-1;k++)
        {
            MPI_Recv(ans,line*cB,MPI_FLOAT,k+1,3,MPI_COMM_WORLD,&status);
            
            //将结果传递给数组P
            for (int i=0;i<line;i++)
            {
                row = k+i;
                for (int j=0;j<cB;j++)
                {
                    C[row*cB+ j] = ans[i*cB + j];
                }
            }
        }
    }
    //其他进程接收数据，计算结果后，发送给主进程
    else
    {
        MPI_Bcast(sizes, 4, MPI_INT, 0, MPI_COMM_WORLD);
        rA = sizes[0];
        cA = rB = sizes[1];
        cB = sizes[2];
        B = (double*)malloc(sizeof(double)*cB*rB);
        //接收广播的数据(矩阵N)
        MPI_Bcast(B,cB*rB,MPI_FLOAT,0,MPI_COMM_WORLD);
        //MPI_Recv(N,Width*Width,MPI_FLOAT,0,0,MPI_COMM_WORLD,&status);
        
        MPI_Recv(buffer,cA*line,MPI_FLOAT,0,1,MPI_COMM_WORLD,&status);
        //计算乘积结果，并将结果发送给主进程
        for (int i=0;i<line;i++)
        {
            for (int j=0;j<cB;j++)
            {
                double temp=0.0;
                for(int k=0;k<cA;k++)
                    temp += buffer[i*cA+k]*B[k*cB+j];
                ans[i*cB+j]=temp;
            }
        }
        //将计算结果传送给主进程
        MPI_Send(ans,line*cB,MPI_FLOAT,0,3,MPI_COMM_WORLD);
        free(B);
    }
    free(buffer);
    free(ans);
    free(sizes);
}

