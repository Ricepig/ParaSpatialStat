#include <cstdlib>
#include "stdio.h"
#include "mpi.h"
#include "paramatrix.h"
#include "parainvert.h"
#include "gdal_priv.h"
#include "cpl_conv.h"

void fatal(const char *message) {
    printf("%s\n",message);
    exit(1); 
}

void testinvert(int argc, char **argv){
    int my_rank,group_size; 
    
    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&group_size); 
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank); 
    
    double *A,*B,*C;
    int cA = 4, rA = 2, cB = 3;
    A = (double*)malloc(sizeof(double)*cA*rA);
    B = (double*)malloc(sizeof(double)*cA*cB);
    C = (double*)malloc(sizeof(double)*rA*cB);
    A[0] = 2;
    A[1] = 4;
    A[2] = 1;
    A[3] = 0;
    A[4] = 3;
    A[5] = -1;
    A[6] = 0;
    A[7] = 3;
    
    B[0] = 1;
    B[1] = 3;
    B[2] = 2;
    B[3] = -2;
    B[4] = 0;
    B[5] = 1;
    B[6] = 3;
    B[7] = -1;
    B[8] = 5;
    B[9] = 0;
    B[10] = 2;
    B[11] = 4;
    paramatrix(A,B,C, rA,cA,cB, group_size, my_rank); 
    MPI_Finalize(); 
    free(A);
    free(B);
    free(C);
}

void testmultiply(int argc, char **argv)
{
    int my_rank,group_size; 
    
    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&group_size); 
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank); 
    double * matrix;
    int order = 3;
    
    if (my_rank==0) {
        matrix = (double*)malloc(sizeof(double)*order*order);
        matrix[0] = 1;
        matrix[1] = 2;
        matrix[2] = 3;
        matrix[3] = 2;
        matrix[4] = 1;
        matrix[5] = 2;
        matrix[6] = 1;
        matrix[7] = 3;
        matrix[8] = 3;
    }
    invert(matrix, order, group_size, my_rank);
    
    MPI_Finalize(); 
    
    if(my_rank==0)
    {
        printMatrix(matrix, order,order);
        free(matrix);
    }    
}

void print_error(const char * message)
{
    
}

int fisher_training(const char * classFile, const char * dataFile, int* bandIndices, int bandCount)
{
    int myRank,rankSize; 
    
    MPI_Comm_size(MPI_COMM_WORLD,&rankSize); 
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank); 
    
    GDALAllRegister();
    GDALDataset * pCateDataset;
    GDALDataset * poDataset;

    pCateDataset = (GDALDataset *) GDALOpen(classFile, GA_ReadOnly);
    poDataset = (GDALDataset *) GDALOpen(dataFile, GA_ReadOnly );    
    if( poDataset == NULL || pCateDataset == NULL )
    {
        print_error("Can't open the data file or the classification file.");        
        exit(1);
    }
    
    GDALRasterBand *catBand = pCateDataset->GetRasterBand(1);
    GDALRasterBand **bands = (GDALRasterBand**)malloc(sizeof(GDALRasterBand*)*bandCount);
    
    for(int i=0;i<bandCount;i++)
    {
        bands[i] = poDataset->GetRasterBand(bandIndices[i]); 
    }
    
    int xSize = catBand->GetXSize();
    int ySize = catBand->GetYSize();    
    int blocksize = (ySize+rankSize-1) / rankSize;
    int localsize = blocksize*(myRank+1)>ySize?ySize-blocksize*myRank:blocksize;
    
    float *scanline;    
    scanline = (float *) CPLMalloc(sizeof(float)*xSize);
    float * cat = (float *)CPLMalloc(sizeof(float)*xSize);
    int * counts = (int*)malloc(sizeof(int)*4);
    
    memset(counts, 0 , sizeof(int)*4);
    for(int i=myRank*blocksize;i<myRank*blocksize+localsize;i++)
    {
        catBand->RasterIO(GF_Read,0, i, xSize, 1, cat, xSize, 1, GDT_Float32, 0, 0 );
        for(int k=0;k<xSize;k++)
        {
            int sig = cat[k]>0?0:1;            
            counts[sig]++;
        }
    }
    MPI_Allreduce(counts, counts+2, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);    
    
    
    double * totals = (double*)malloc(sizeof(double)*bandCount*4);
    double * means = totals+bandCount*2; 
    memset(totals, 0 , sizeof(double)*bandCount*4);
    // Calc Mi
    for(int i = myRank*blocksize;i<myRank*blocksize+localsize;i++)
    {
        catBand->RasterIO(GF_Read,0, i, xSize, 1, cat, xSize, 1, GDT_Float32, 0, 0 );
        for(int j=0;j<bandCount;j++)
        {
            bands[j]->RasterIO(GF_Read, 0, i, xSize, 1, scanline, xSize, 1, GDT_Float32, 0, 0);
            for(int k=0;k<xSize;k++)
            {
                int sig = cat[k]>0?0:1;
                if(scanline[k]>0)
                        totals[bandCount*sig+j] += scanline[k];                
            }
        }
    }
    
    MPI_Allreduce(totals, means, bandCount*2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for(int i=0;i<bandCount*2;i++)
        means[i] /= counts[2+i/bandCount];
    
    // Calc Sw
    int scount = bandCount*bandCount;
    
    double *localS = (double*)malloc(sizeof(double)*scount*2);
    double *S = localS+scount;
    memset(localS, 0, sizeof(double)*scount*2);
    
    int sscount=scount*xSize;
    double *ss = (double*)malloc(sizeof(double)*sscount);
    
    for(int i = myRank*blocksize;i<myRank*blocksize+localsize;i++)
    {
        catBand->RasterIO(GF_Read,0, i, xSize, 1, cat, xSize, 1, GDT_Float32, 0, 0 );
        memset(ss, 0, sizeof(double)*sscount);
        for(int j=0;j<bandCount;j++)
        {
            bands[j]->RasterIO(GF_Read, 0, i, xSize, 1, scanline, xSize, 1, GDT_Float32, 0, 0);
            for(int k=0;k<xSize;k++)
            {
                int sig = cat[k]>0?0:1;
                double avg = means[sig*bandCount+j];
                double* sss = ss+k*scount;
                double value = (scanline[k]<0?0:scanline[k])-avg;
                for(int m=0;m<bandCount;m++)
                {
                    if(j>m)
                    {
                        sss[j*bandCount+m] = sss[m*bandCount+j] = sss[j*bandCount+m] * value;
                    }
                    else if(j==m)
                    {  
                        sss[j*bandCount+m] = value*value;
                    }
                    else
                    {
                        sss[j*bandCount+m] = sss[m*bandCount+j] = value;
                    }
                }
            }
        }
        
        for(int k=0;k<xSize;k++)
            for(int i=0;i<bandCount;i++)
                for(int j=0;j<bandCount;j++)
                    localS[i*bandCount+j] += ss[k*scount+i*bandCount+j];
    }
    
    free(ss);
    
    MPI_Reduce(localS, S, scount, MPI_DOUBLE,  MPI_SUM, 0, MPI_COMM_WORLD);
    
    if(bandCount<rankSize*10)
        invert(S, bandCount, 1, 0);
    else
        invert(S, bandCount, rankSize, myRank);
    
    /*if(myRank==0)
        printMatrix(S, bandCount, bandCount);*/
    
    double* m12 = NULL;
    double* sum12 = NULL;
    double* w = NULL;
    if(myRank ==0)
    {
        m12 = (double*)malloc(sizeof(double)*bandCount);
        sum12 = (double*)malloc(sizeof(double)*bandCount);
        for(int i=0;i<bandCount;i++)
        {
            m12[i] =  means[i]-means[bandCount+i];
            sum12[i] = means[i]+means[bandCount+i];
        }
        w = (double*)malloc(sizeof(double)*bandCount);
    }
    
    if(bandCount<rankSize*10){
        if(myRank==0)
        singlematrix(S, m12, w, bandCount, bandCount, 1);
    }else{
        paramatrix(S, m12, w, bandCount, bandCount, 1, rankSize, myRank);
    }
    
    if(myRank==0)
    {
        double sum = 0;
        for(int i=0;i<bandCount;i++)
        {
            sum+=w[i]*means[i]*counts[2] + w[i]*means[bandCount+i]*counts[3]; 
        }
        sum/=(counts[2]+counts[3]);
        
        printf("Final Output:\nPhi=%f\nDiscriminant function:\n", sum);
        for(int i=0;i<bandCount;i++) 
        {
            if(i>0)
                printf(" + ");
            printf("%f*X%d",w[i],i+1);
        }
        
        free(w);
        free(m12);
    }
    
    free(localS);
    free(totals);
    free(counts);
    free(cat);
    free(scanline);
    
    GDALClose(poDataset);
    free(bands);
    
    return 1;
}

void printUsage()
{
    printf("Usage:\nfisher [classification file] [data file] [index of band1] [index of band2] ...\rn");
}
                
int main(int argc, char **argv) {
    testmultiply(argc, argv);
    return 1;
    GDALAllRegister();
    MPI_Init(&argc,&argv);
    if(strcmp("test", argv[0])){
        int indices[7];
        indices[0] = 1;
        indices[1] = 2;
        indices[2] = 3;
        indices[3] = 4;
        indices[4] = 5;
        indices[5] = 6;
        indices[6] = 7;
        fisher_training("/home/ricepig/class1.tif", "/home/ricepig/data1.tif", indices, 3);
    }else{
        if(argc<4){
            printUsage();
            
        }else{
            int start = 3;
            int * ind = (int*)malloc(sizeof(int)*(argc-start));
            for(int i = start;i<argc;i++)
                ind[i-start] = atoi(argv[i]);
            
            fisher_training(argv[1], argv[2],ind, argc-start);
        }        
    }
    MPI_Finalize();
    return 1;
}

    