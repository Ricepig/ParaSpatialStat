/* 
 * File:   main.cpp
 * Author: ricepig
 *
 */

#include <cstdlib>
#include "gdal_priv.h"
#include "cpl_conv.h"

using namespace std;

/**
 * print usage of this program
 */
void print_usage()
{
	printf("Fisher Discriminance Program (Training) Usage:\n");
	printf("fisher [classification raster] [raster 1] [raster 2] ... [raster n]\n");
}

bool create_raster(const char* filename, double left, double top, int nXSize, int nYSize, double pixelSize, char* spatialRefWkt )
{
	const char *pszFormat = "GTiff";
	
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if( poDriver == NULL )
	{
		printf("[ERROR] Can't find the driver for writing GeoTiff.\n");
		return false;
	}
    
	GDALDataset *poDstDS;       
    char **papszOptions = NULL;
    poDstDS = poDriver->Create( filename, nXSize, nYSize, 1, GDT_Float32, 
                                papszOptions );  
	
	if(poDstDS == NULL)
	{
		printf("[ERROR] Can't create the raster file as output.\n");
		return false;
	}
	
	double adfGeoTransform[6] = {left, pixelSize, 0, top, 0, -pixelSize};
	poDstDS->SetGeoTransform(adfGeoTransform);
	
	poDstDS->SetProjection(spatialRefWkt);
	CPLFree(spatialRefWkt);
	GDALClose((GDALDatasetH)poDstDS);
}

/*
 * 
 */
int main(int argc, char** argv) {

    return 0;
    
    GDALAllRegister();
    
    GDALDataset  *poDataset;

    

    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
    if( poDataset == NULL )
    {
        ...;
    }
}

