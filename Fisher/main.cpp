#include <cstdlib>
#include "mpi.h"
#include "paramatrix.h"
#include "parainvert.h"
#include "string.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include "common.h"

using namespace std;

struct Timespans {
	double io;
	double computing;
};

/**
 * print usage of this program
 */
void print_usage()
{
	printf("Fisher Discriminance Program (Training) Usage:\n");
	printf("fisher [classification raster] [raster 1] [raster 2] ... [raster n]\n");
}

char** split_tokens(char * arg, int& count){
	const char signal = ',';
	
	count = 0;
	
	for(int i=0;i<strlen(arg);i++){
		if(arg[i] == signal)
			count++;
	}
	
	count++;
	
	char ** rt = (char**)malloc(sizeof(char*)*count);
	char * token = strtok(arg, ",");
	rt[0] = token;
	int index = 1;
	while((token = strtok(NULL, ","))){
		if(index >= count)
			return NULL;
		rt[index] = token;
		index++;
	}
	if(index != count)
		return NULL;
	return rt;
}

/**
 * parse arguments
 * @param argc
 * @param argv
 * @return true if arguments are valid
 */
bool parse_args(int argc, char **argv, char *** pLearnFiles, char *** pInferFiles, int * pBandCount, char ** pOutputFile, char** pClassFile)
{
	bool hasClass = false;
	bool hasOutput = false;
	bool hasLearn = false;
	bool hasInfer = false;
	int count1, count2;
	
    for(int i=1;i<argc;i++){
		
		if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0){
            return false;
        }
		
		if(strcmp(argv[i], "-pc") == 0 && i+1<argc){
			*pClassFile = argv[i+1];
			hasClass = true;
		} else if(strcmp(argv[i], "-o") == 0 && i+1<argc){
			*pOutputFile = argv[i+1];
			hasOutput = true;
		} else if(strcmp(argv[i], "-pi") == 0 && i+1<argc){
			*pInferFiles = split_tokens(argv[i+1], count1);
			if(pInferFiles != NULL)
				hasInfer = true;
		} else if(strcmp(argv[i], "-pl") == 0 && i+1<argc){
			*pLearnFiles = split_tokens(argv[i+1], count2);
			if(pLearnFiles != NULL)
				hasLearn = true;
		}	
    }
	
	if(!hasOutput || !hasInfer || !hasLearn || !hasClass)
		return false;
	
	if(count1!=count2)
		return false;
	
	*pBandCount = count1;
	return true;
}

void fatal(const char *message) {
    printf("%s\n",message);
    exit(1); 
}

int create_raster(const char* filename, double left, double top, int nXSize, int nYSize, double pixelSize, const char* spatialRefWkt )
{
	const char *pszFormat = "GTiff";
	
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if( poDriver == NULL ){
		printf("[ERROR] Can't find the driver for writing GeoTiff.\n");
		return 1;
	}
    
	GDALDataset *poDstDS;       
    char **papszOptions = NULL;
    poDstDS = poDriver->Create( filename, nXSize, nYSize, 1, GDT_Byte, 
                                papszOptions );  
	
	if(poDstDS == NULL){
		printf("[ERROR] Can't create the raster file as output.\n");
		return 1;
	}
	
	double adfGeoTransform[6] = {left, pixelSize, 0, top, 0, -pixelSize};
	poDstDS->SetGeoTransform(adfGeoTransform);
	poDstDS->SetProjection(spatialRefWkt);
	
	GDALClose((GDALDatasetH)poDstDS);
	return 0;
}

int open_raster(const char* filename, GDALAccess eAccess, GDALDataset ** pDS, GDALRasterBand** pBand)
{
	*pDS = (GDALDataset*)GDALOpen(filename, eAccess);
	if(*pDS == NULL){
		printf("[ERROR] Can't open the output file.\n");
		return 1;
	}
	*pBand = (*pDS)->GetRasterBand(1);
	return 0;
}

int close_raster(GDALDataset *pDS)
{
	GDALClose((GDALDatasetH)pDS);
	return 0;
}

int fisher_discriminating(char ** dataFiles, int bandCount, const char * outputFile, double * coefficients, double pivot, int myRank, int rankSize, struct Timespans &benchmark)
{
	double t1,t2,t3,t4;
	
	t1 = MPI_Wtime();
	GDALDataset ** datasets = new GDALDataset*[bandCount];
	GDALRasterBand **bands = new GDALRasterBand*[bandCount];
	
    for(int i=0;i<bandCount;i++){
		if(open_raster(dataFiles[i], GA_ReadOnly, datasets+i, bands+i)!=0)
			return 1;
    }
    if(myRank == 0){
		double adfGeoTransform[6];
		datasets[0]->GetGeoTransform( adfGeoTransform );
		if(create_raster(outputFile, adfGeoTransform[0], adfGeoTransform[3],
			datasets[0]->GetRasterXSize(), datasets[0]->GetRasterYSize(), adfGeoTransform[1], datasets[0]->GetProjectionRef()) != 0)
			return 1;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	double calc_duration = 0;
	
	GDALDataset * destset;
	GDALRasterBand * destband;
	if(open_raster(outputFile, GA_Update, &destset, &destband) != 0)
		return 1;
	
	int width = destband->GetXSize();
	int height = destband->GetYSize();
	float* buffer = (float*)malloc(sizeof(float)*width); 
	double* results = (double*)malloc(sizeof(double)*width);
	char* buffer2 = (char*)malloc(sizeof(char)*width);
	int count = 0;
	for(int i=0;i<height+rankSize-1;i+=rankSize){
		int row = i+myRank;
		
		if(row < height){
			count++;
			for(int j=0;j<bandCount;j++){
				bands[j]->RasterIO(GF_Read, 0, row, width, 1, buffer, width, 1, GDT_Float32, 0, 0);
				t2 = MPI_Wtime();
				for(int k=0;k<width;k++){
					results[k] += (buffer[k] * coefficients[j]);
				}
				t3 = MPI_Wtime();
				calc_duration += t3-t2;
			}
			t2 = MPI_Wtime();
			for(int k=0;k<width;k++){
				buffer2[k] = results[k]>pivot?0:1;
			}
			t3 = MPI_Wtime();
			calc_duration += t3-t2;
			
			destband->RasterIO(GF_Write, 0, row, width, 1, buffer2, width, 1, GDT_Byte, 0, 0);
			memset((void*)results, 0, sizeof(double)*width);
		}
	}
	
	printf("thread %d, comp time: %f, rows: %d\n", myRank, calc_duration,count );
	
	double * temp = new double[2];
	*temp = calc_duration;
	MPI_Reduce(temp, temp+1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	free(buffer2);
	free(results);
	free(buffer);
	
    close_raster(destset);
	for(int j=0;j<bandCount;j++)
		close_raster(datasets[j]);
	delete[] bands;
	delete[] datasets;
	t4 = MPI_Wtime();
	if(myRank == 0){
		benchmark.computing = temp[1];
		benchmark.io = t4-t1-temp[1];
	}
	delete[] temp;
	return 0;
}

int fisher_training(const char * classFile, char ** dataFiles, int bandCount, int myRank, int rankSize, double ** coefficients, double &pivot, struct Timespans &benchmark)
{
    double t1 = MPI_Wtime();
    double t2,t3,iotime = 0;
    // Open rasters
    GDALDataset * pCateDataset;
    pCateDataset = (GDALDataset *) GDALOpen(classFile, GA_ReadOnly);
	if(myRank==0)
		cout<< "class:" << classFile<<endl;
    if(pCateDataset == NULL){
        cout<<"[ERROR] Can't open the classification file: "<< classFile << endl;
        exit(1);
    }
    
    GDALDataset ** datasets = new GDALDataset*[bandCount];
    for(int i=0;i<bandCount;i++){
        datasets[i] = (GDALDataset *) GDALOpen(dataFiles[i], GA_ReadOnly );
		if(myRank==0)
			cout<<"band:"<<dataFiles[i]<<endl;
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
    for(int i=myRank*blocksize;i<myRank*blocksize+localsize;i++){
        t2 = MPI_Wtime();
        catBand->RasterIO(GF_Read,0, i, xSize, 1, cat, xSize, 1, GDT_Float32, 0, 0 );
        t3 = MPI_Wtime();
        iotime += (t3-t2);
        for(int k=0;k<xSize;k++){
            int sig = cat[k]>0?0:1;            
            counts[sig]++;
        }
    }
    MPI_Allreduce(counts, counts+2, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);    
    
    double * totals = (double*)malloc(sizeof(double)*bandCount*4);
    double * means = totals+bandCount*2; 
    memset(totals, 0 , sizeof(double)*bandCount*4);
    
    // Stage 2: Calc Mi
    for(int i = myRank*blocksize;i<myRank*blocksize+localsize;i++){
        t2 = MPI_Wtime();
        catBand->RasterIO(GF_Read,0, i, xSize, 1, cat, xSize, 1, GDT_Float32, 0, 0 );
        t3 = MPI_Wtime();
        iotime += (t3-t2);
        for(int j=0;j<bandCount;j++){
            t2 = MPI_Wtime();
            bands[j]->RasterIO(GF_Read, 0, i, xSize, 1, scanline, xSize, 1, GDT_Float32, 0, 0);
            t3 = MPI_Wtime();
            iotime += (t3-t2);
            for(int k=0;k<xSize;k++){
                int sig = cat[k]>0?0:1;
                if(scanline[k]>-9998)
                        totals[bandCount*sig+j] += scanline[k];                
            }
        }
    }
    
    MPI_Allreduce(totals, means, bandCount*2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for(int i=0;i<bandCount*2;i++){
        means[i] /= counts[2+i/bandCount];
	}
	
    // Stage 3: Calc Sw
    int scount = bandCount*bandCount;
    
    double *localS = (double*)malloc(sizeof(double)*scount*2);
    double *S = localS+scount;
    memset(localS, 0, sizeof(double)*scount*2);
    
    int sscount=scount*xSize;
    double *ss = (double*)malloc(sizeof(double)*sscount);
    
    for(int i = myRank*blocksize;i<myRank*blocksize+localsize;i++){
        catBand->RasterIO(GF_Read,0, i, xSize, 1, cat, xSize, 1, GDT_Float32, 0, 0 );
        memset(ss, 0, sizeof(double)*sscount);
        for(int j=0;j<bandCount;j++){
            bands[j]->RasterIO(GF_Read, 0, i, xSize, 1, scanline, xSize, 1, GDT_Float32, 0, 0);
            for(int k=0;k<xSize;k++){
                int sig = cat[k]>0?0:1;
                double avg = means[sig*bandCount+j];
                double* sss = ss+k*scount;
                double value = (scanline[k]<-9998?0:scanline[k])-avg;
                for(int m=0;m<bandCount;m++){
                    if(j>m){
                        sss[j*bandCount+m] = sss[m*bandCount+j] = sss[j*bandCount+m] * value;
                    } else if(j==m) {  
                        sss[j*bandCount+m] = value*value;
                    } else {
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
    
    invert(S, bandCount, 1, 0);    
    double* m12 = NULL;
    double* sum12 = NULL;
    double* w = NULL;
    if(myRank ==0){
        m12 = (double*)malloc(sizeof(double)*bandCount);
        sum12 = (double*)malloc(sizeof(double)*bandCount);
        for(int i=0;i<bandCount;i++){
            m12[i] =  means[i]-means[bandCount+i];
            sum12[i] = means[i]+means[bandCount+i];
        }
        w = (double*)malloc(sizeof(double)*bandCount);
    }
    
    if(myRank==0)
        singlematrix(S, m12, w, bandCount, bandCount, 1);
    
    if(myRank==0){
        double sum = 0;
        for(int i=0;i<bandCount;i++){
            sum+=w[i]*means[i]*counts[2] + w[i]*means[bandCount+i]*counts[3]; 
        }
        sum/=(counts[2]+counts[3]);
        
        double t4 = MPI_Wtime();
		pivot = sum;
        
        for(int i=0;i<bandCount;i++) 
			(*coefficients)[i] = w[i];
        cout<<endl;
		benchmark.io = iotime;
		benchmark.computing = t4-t1-iotime;
       
        free(w);
		free(sum12);
        free(m12);
    }
    
    free(localS);
    free(totals);
    free(counts);
    free(cat);
    free(scanline);
    
    GDALClose(pCateDataset);
    delete[] bands;
	for(int i=0;i<bandCount;i++){
        GDALClose(datasets[i]);
    }
    delete[] datasets;
    
    return 0;
}
        
int main(int argc, char **argv) {
	char ** pLearnFiles;
	char * pClassFile;
	char ** pInferFiles;
	char * pOutputFile;
	int bandCount;
	
	if(!parse_args(argc, argv, &pLearnFiles, &pInferFiles, &bandCount, &pOutputFile, &pClassFile)){
        print_usage();
        return 0;
    }
    
    GDALAllRegister();
    MPI_Init(&argc,&argv);
	int myRank,rankSize; 
    MPI_Comm_size(MPI_COMM_WORLD,&rankSize); 
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
	struct Timespans benchmark1, benchmark2;
	
	double pivot = 0;
	double* coefficients = new double[bandCount];
	
    fisher_training(pClassFile, pLearnFiles, bandCount, myRank, rankSize, &coefficients, pivot, benchmark1);
	fisher_discriminating(pInferFiles, bandCount, pOutputFile, coefficients, pivot, myRank, rankSize, benchmark2);
	if(myRank==0){
		cout << "[OUTPUT] Phi: " << pivot << endl;
		cout << "[OUTPUT] Discriminant function: ";
		for(int i=0;i<bandCount;i++){
			if(i>0)
				cout<< " + ";
			cout<<coefficients[i]<<"x"<<"X"<<(i+1);
		}
		cout << endl;
		cout << "[DEBUG] [TIMESPAN] [IO] " << benchmark1.io + benchmark2.io << endl;
		cout << "[DEBUG] [TIMESPAN] [COMPUTING] " << benchmark1.computing + benchmark2.computing << endl;
		cout << "[DEBUG] [TIMESPAN] [TOTAL] " << benchmark1.io + benchmark2.io + benchmark1.computing + benchmark2.computing << endl;
    }
	delete[] coefficients;
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
    
