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

using namespace std;

void fatal(const char* msg)
{
    printf("Fatal error: %s\n", msg);
    exit(0);
}

struct ElementContainer loadData(const char * shapefile, int fieldIndex, int rankSize, int myRank)
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
    return ec;
}

int variogram(const char * shapefile, int fieldIndex, int lag, int lagCount, int rankSize, int myRank)
{
    struct ElementContainer * ec = loadData(shapefile, fieldIndex, rankSize, myRank);
    
    
    
    
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



