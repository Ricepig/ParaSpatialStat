/* 
 * File:   paramatrix.h
 * Author: Ricepig
 *
 * Created on 2013年11月12日, 上午3:48
 */

#ifndef PARAMATRIX_H
#define	PARAMATRIX_H

#ifdef	__cplusplus
extern "C" {
#endif

void paramatrix(double *A, double *B, double* C, int rA, int cA, int cB, int numprocs, int myid);
void singlematrix(double *A, double *B, double *C, int rA, int cA, int cB);

#ifdef	__cplusplus
}
#endif

#endif	/* PARAMATRIX_H */

