/* 
 * File:   pararevert.h
 * Author: Ricepig
 *
 * Created on 2013年11月12日, 上午3:56
 */

#ifndef PARAREVERT_H
#define	PARAREVERT_H

#ifdef	__cplusplus
extern "C" {
#endif

void invert(double* matrix, int order, int group_size, int my_rank);

void printMatrix(double* matrix, int sizex, int sizey);

#ifdef	__cplusplus
}
#endif

#endif	/* PARAREVERT_H */

