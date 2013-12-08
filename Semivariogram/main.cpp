/* 
 * File:   main.cpp
 * Author: ricepig
 *
 * Created on 2013年11月30日, 下午11:13
 */

#include <cstdlib>
#include <string>
#include "ogrsf_frmts.h"
#include "ElementContainer.h"
#include "LagContainer.h"
#include "math.h"
#include "mpi.h"

#define EPSILON 0.0000001

using namespace std;

void fatal(const char* msg)
{
    printf("Fatal error: %s\n", msg);
    exit(0);
}

struct ElementContainer loadData(const char * shapefile, int fieldIndex, int rankSize, int myRank, int *pcount)
{
    OGRDataSource *poDS = OGRSFDriverRegistrar::Open(shapefile, FALSE);
    if(poDS == NULL)
        fatal("Can't open the data source.");
    string path(shapefile);
    int pos = path.find_last_of("\\");
    OGRLayer  *poLayer = poDS->GetLayerByName( path.substr(pos+1, path.length()-pos-5) );
    if(poLayer == NULL)
        fatal("Can't open the data source.");
    
    struct ElementContainer *ec = (struct ElementContainer*)malloc(sizeof(struct ElementContainer));
    ECInit(ec);
    
    OGRFeature *poFeature;
    int count = 0;
    poLayer->ResetReading();
    while( (poFeature = poLayer->GetNextFeature()) != NULL )
    {
        if(count % rankSize == myRank)
        {
            OGRGeometry *poGeometry = poFeature->GetGeometryRef();
            if(poGeometry != NULL){
                
                if(wkbFlatten(poGeometry->getGeometryType()) =! wkbPoint)
                    fatal("The layer type is restricted to point.");
            
                OGRPoint *poPoint = (OGRPoint*)poGeometry;
                
                ECAdd(ec, poPoint->getX(), poPoint->getY(), (float)poFeature->GetFieldAsDouble(fieldIndex));
            }
            
        }
        OGRFeature::DestroyFeature( poFeature );
        count++;
    }

    OGRDataSource::DestroyDataSource( poDS );
    *pcount = count;
    return ec;
}

void classify(struct ElementContainer* ec, int lag, int lagCount)
{
    
}


inline void inner_classify(struct Element* poE1, struct Element* poE2, struct LagContainer *lc, double sqLag, double range)
{
    double dist = (poE1->X - poE2->X) * (poE1->X - poE2->X) + (poE1->Y - poE2->Y) * (poE1->Y - poE2->Y);
    if(dist <= range)
    {
        int index = (dist + sqLag - EPSILON)/sqLag;
        (lc->counts[index])++;
        (lc->sums[index])+= abs(poE1->Value - poE2->Value);
    }
}

void classify(struct ElementContainer* ec1, struct ElementContainer* ec2, struct LagContainer *lc, double sqLag, int lagCount)
{
    double range = sqLag * lagCount;
    struct Element* poE1 = ec1->Head;
    
    for(size_t i=0;i<ec1->Length;i++, poE1++){
        struct Element* poE2 = ec2->Head;
        for(size_t j=0;j<ec1->Length;j++, poE2++)
            inner_classify(poE1, poE2, lc, sqLag, range);
    }
        
}

void classify_self(struct ElementContainer* ec, struct LagContainer *lc, double sqLag, int lagCount)
{
    double range = sqLag * lagCount;
    struct Element* poE1 = ec->Head;
    for(size_t i=0;i<ec->Length;i++, poE1++){
        struct Element* poE2 = ec->Head+i+1;
        for(size_t j=i+1;j<ec->Length;j++, poE2++)
            inner_classify(poE1, poE2, lc, sqLag, range);
    }
}

int variogram(const char* shapefile, int fieldIndex, double lag, int lagCount, int rankSize, int myRank)
{
    int count;
    struct ElementContainer * ec = loadData(shapefile, fieldIndex, rankSize, myRank, &count);
    
    struct LagContainer * lc = (struct LagContainer)malloc(sizeof(struct Element));
    LCInit(lc, lagCount);
    
    double sqLag = lag*lag;
    
    int blockSize = (count + rankSize - 1)/rankSize;
    struct ElementContainer * ec2 = (struct ElementContainer *)malloc(sizeof(struct ElementContainer));
    ECInitWithSize(ec2, blockSize);
    
    for(int i=0;i<rankSize;i++)
    {
        int currentSize = blockSize*(rankSize-1) + i <= count?blockSize:blockSize-1;
        MPI_Bcast(ec2->Head, currentSize*sizeof(struct Element), MPI_BYTE, i, MPI_COMM_WORLD);
        ec2->Length = currentSize;
        if(i==myRank){
            classify_self(ec, lc, sqLag, lagCount);
        }else{
            if(i<myRank)
                classify(ec, ec2, lc, sqLag, lagCount);
        }
    }
    
    struct LagContainer * lc2 = (struct LagContainer)malloc(sizeof(struct Element));
    LCInit(lc2, lagCount);
    MPI_Allreduce(lc->counts, lc2->counts, lc->size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(lc->sums, lc2->sums, lc->size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    LCCalcAverage(lc);
    
    
    
}


int variogram2(const char * shapefile, int fieldIndex, int lag, int lagCount, int rankSize, int myRank)
{
    struct ElementContainer * ec = loadData(shapefile, fieldIndex, rankSize, myRank);
    int dest, src;
    for(size_t stride = 1;stride<rankSize;stride++)
    {
        // phase 1: send data from the first half stride to the second half
        if(myRank / stride % 2 == 0){
            //sender
            dest = myRank + stride;
            if(dest<rankSize){
                //send msg to dest;
            }
        }else{
            //receiver
            src = myRank - stride;
        }
        
        // phase 2: send data from the second half stride to the first half of the next stride
        if(myRank / stride %2 == 0 ){
            //receiver
            src = myRank - stride;
            if(src >= 0){
                //receive
            }
        } else {
            dest = myRank + stride;
            if(dest<rankSize){
            
            }
        }
            
    }
    
    
    
}

int main(int argc, char** argv) {
    
    OGRRegisterAll();
    OGRDataSource *poDS;

    poDS = OGRSFDriverRegistrar::Open("point.shp", FALSE);
    if( poDS == NULL )
    {
        printf( "Open failed.\n" );
        exit( 1 );
    }

    OGRLayer  *poLayer;
    poLayer = poDS->GetLayerByName( "point" );

    OGRFeature *poFeature;

    poLayer->ResetReading();
    while( (poFeature = poLayer->GetNextFeature()) != NULL )
    {
        OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
        int iField;

        for( iField = 0; iField < poFDefn->GetFieldCount(); iField++ )
        {
            OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn( iField );

            if( poFieldDefn->GetType() == OFTInteger )
                printf( "%d,", poFeature->GetFieldAsInteger( iField ) );
            else if( poFieldDefn->GetType() == OFTReal )
                printf( "%.3f,", poFeature->GetFieldAsDouble(iField) );
            else if( poFieldDefn->GetType() == OFTString )
                printf( "%s,", poFeature->GetFieldAsString(iField) );
            else
                printf( "%s,", poFeature->GetFieldAsString(iField) );
        }

        OGRGeometry *poGeometry;

        poGeometry = poFeature->GetGeometryRef();
        if( poGeometry != NULL 
            && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint )
        {
            OGRPoint *poPoint = (OGRPoint *) poGeometry;

            printf( "%.3f,%3.f\n", poPoint->getX(), poPoint->getY() );
        }
        else
        {
            printf( "no point geometry\n" );
        }       
        OGRFeature::DestroyFeature( poFeature );
    }

    OGRDataSource::DestroyDataSource( poDS );
    
    
    return 0;
}



