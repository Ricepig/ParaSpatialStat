#include <cstdlib>
#include "mpi.h"
#include "paramatrix.h"
#include "parainvert.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "common.h"

using namespace std;
/**
 * print usage of this program
 */
void print_usage()
{
	printf("Fisher Discriminance Program (Training) Usage:\n");
	printf("fisher [classification raster] [raster 1] [raster 2] ... [raster n]\n");
}

/**
 * check arguments
 * @param argc
 * @param argv
 * @return true if arguments are valid
 */
bool check_args(int argc, char **argv)
{
    for(int i=1;i<argc;i++){
        if(strcmp(argv[i], "-h") || strcmp(argv[i], "--help")){
            return false;
        }
    }
    if(argc<4){
        return false;
    }
    return true;
}

void fatal(const char *message) {
    printf("%s\n",message);
    exit(1); 
}

int fisher_training(const char * classFile, char ** dataFiles, int bandCount)
{
    int myRank,rankSize; 
    
    MPI_Comm_size(MPI_COMM_WORLD,&rankSize); 
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank); 
    
    GDALAllRegister();
    double t1 = MPI_Wtime();
    double t2,t3,iotime = 0;
    // Open rasters
    GDALDataset * pCateDataset;
    pCateDataset = (GDALDataset *) GDALOpen(classFile, GA_ReadOnly);
    if(myRank==0)
        cout<<"[DEBUG] [OPTIONS] class raster:"<<classFile<<endl;
    
    if(pCateDataset == NULL){
        cout<<"[ERROR] Can't open the classification file: "<< classFile << endl;
        exit(1);
    }
    
    GDALDataset ** datasets = new GDALDataset*[bandCount];
    for(int i=0;i<bandCount;i++){
        datasets[i] = (GDALDataset *) GDALOpen(dataFiles[i], GA_ReadOnly );
        if(myRank==0)
            cout<<"[DEBUG] [OPTIONS] data raster:"<<dataFiles[i]<<endl;
        if(datasets[i] == NULL){
            cout<<"[ERROR] Can't open raster file: "<< dataFiles[i] << endl;
            exit(1);
        }
    }
    
    GDALRasterBand *catBand = pCateDataset->GetRasterBand(1);
    GDALRasterBand **bands = new GDALRasterBand*[bandCount];
    for(int i=0;i<bandCount;i++){
        bands[i] = datasets[i]->GetRasterBand(1); 
    }
    
    // Stage 1: Calc item count for each categories
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
        t2 = MPI_Wtime();
        catBand->RasterIO(GF_Read,0, i, xSize, 1, cat, xSize, 1, GDT_Float32, 0, 0 );
        t3 = MPI_Wtime();
        iotime += (t3-t2);
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
    
    // Stage 2: Calc Mi
    for(int i = myRank*blocksize;i<myRank*blocksize+localsize;i++)
    {
        t2 = MPI_Wtime();
        catBand->RasterIO(GF_Read,0, i, xSize, 1, cat, xSize, 1, GDT_Float32, 0, 0 );
        t3 = MPI_Wtime();
        iotime += (t3-t2);
        for(int j=0;j<bandCount;j++)
        {
            t2 = MPI_Wtime();
            bands[j]->RasterIO(GF_Read, 0, i, xSize, 1, scanline, xSize, 1, GDT_Float32, 0, 0);
            t3 = MPI_Wtime();
            iotime += (t3-t2);
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
    
    // Stage 3: Calc Sw
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
    
    // Stage 4: Calc reverse matrix of Sw and final result
    if(bandCount<rankSize*10)
        invert(S, bandCount, 1, 0);
    else
        invert(S, bandCount, rankSize, myRank);
   
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
        
        double t4 = MPI_Wtime();
        cout << "[OUTPUT] Phi: " << sum << endl;
        cout << "[OUTPUT] Discriminant function: ";
        for(int i=0;i<bandCount;i++) 
        {
            if(i>0)
                cout<<" + ";
            cout<<w[i]<<"x"<<i+1;
        }
        cout<<endl;
        
        cout << "[DEBUG] [TIMESPAN] [IO] " << iotime << endl;
        cout << "[DEBUG] [TIMESPAN] [COMPUTING] " << t4-t1-iotime << endl;
        cout << "[DEBUG] [TIMESPAN] [TOTAL] " << t4-t1 << endl;
        free(w);
        free(m12);
    }
    
    free(localS);
    free(totals);
    free(counts);
    free(cat);
    free(scanline);
    
    delete[] bands;
    GDALClose(pCateDataset);
    for(int i=0;i<bandCount;i++){
        GDALClose(datasets[i]);
    }
    delete[] datasets;
    
    return 0;
}
        
int main(int argc, char **argv) {
    if(check_args(argc, argv)==false){
        print_usage();
        return 0;
    }
    
    GDALAllRegister();
    MPI_Init(&argc,&argv);
    fisher_training(argv[1], argv+2, argc-2);
    MPI_Finalize();
    return 0;
}


/*
void testmultiply(int argc, char **argv){
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

void testinvert(int argc, char **argv)
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
*/
    
