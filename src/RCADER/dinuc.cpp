// Copyright 2014 Hamed S. Najafabadi

/********************************************************************

This file is part of RCOpt.

RCOpt is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RCOpt is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RCOpt.  If not, see <http://www.gnu.org/licenses/>.

********************************************************************/


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>

#include "declarations.h"

//////////////////////////////////////////////////////////

const char *__nuc1 = "AC";
const char *__nuc2 = "TG";

const char *__dinuc1[] = {
					"AA", "AC", "AG", "AT",
					"CA", "CC", "CG",
					"GA", "GC",
					"TA" };
const char *__dinuc2[] = {
					"TT", "GT", "CT", "AT",
					"TG", "GG", "CG",
					"TC", "GC",
					"TA" };

///////////////////////////////////////////////////////////////////////////////////////////
void _calc_nuc(	s_seq *seq )
{
	int i, j;
	for( j = 0; j < 2; j ++ )
	{
		seq ->nuc[ j ] = 0;
		
		for( i = 0; i < seq ->seq_length; i ++ )
			if( seq ->seq[ i ] == __nuc1[ j ] ||
				seq ->seq[ i ] == __nuc2[ j ] )
					seq ->nuc[ j ] ++;
					
		seq ->nuc[ j ] /= seq ->seq_length;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
void _calc_dinuc( s_seq *seq )
{
	int i, j;
	for( j = 0; j < 10; j ++ )
	{
		seq ->dinuc[ j ] = 0;
		
		for( i = 0; i < seq ->seq_length - 1; i ++ )
			if( ( seq ->seq[ i ] == __dinuc1[ j ][ 0 ] && seq ->seq[ i+1 ] == __dinuc1[ j ][ 1 ] ) ||
				( seq ->seq[ i ] == __dinuc2[ j ][ 0 ] && seq ->seq[ i+1 ] == __dinuc2[ j ][ 1 ] ) )
					seq ->dinuc[ j ] ++;
					
		seq ->dinuc[ j ] /= (seq ->seq_length-1);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
void calc_nuc_dinuc(
	s_seq *seqs[],
	int num_seqs )
{
	int i;
	for( i = 0; i < num_seqs; i ++ )
	{
		_calc_nuc( seqs[ i ] );
		_calc_dinuc( seqs[ i ] );
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
void correct_for_dinuc(
	s_seq *seqs[],
	int num_seqs )
{
	// initialize the objects for GSL regression
	int max_num_cov = 2+10 + 1; // nuc + dinuc + intercept
	gsl_matrix *X = gsl_matrix_calloc( num_seqs, max_num_cov );
	gsl_vector *Y = gsl_vector_alloc( num_seqs );
	gsl_vector *beta = gsl_vector_alloc( max_num_cov );
	gsl_matrix *cov = gsl_matrix_alloc( max_num_cov, max_num_cov );
	gsl_multifit_linear_workspace *wspc = gsl_multifit_linear_alloc( num_seqs, max_num_cov );
	
	int i;
	for( i = 0; i < num_seqs; i ++ )
	{
		gsl_vector_set( Y, i, seqs[ i ] ->residual );
		
		// set the elements of the covariates
		int j;
		for( j = 0; j < 2; j ++ )
			gsl_matrix_set( X, i, j, seqs[ i ] ->nuc[ j ] + drand48()/1000000 );
		for( j = 0; j < 10; j ++ )
			gsl_matrix_set( X, i, 2+j, seqs[ i ] ->dinuc[ j ] + drand48()/1000000 );
		gsl_matrix_set( X, i, 12, 1 );
	}
	
	// solve the least-square regression
	double chisq;
	gsl_multifit_linear( X, Y, beta, cov, &chisq, wspc );

	// get the residuals, and also calculate the correlation between the fit and the observed values
	double EX = 0, EX2 = 0, EY = 0, EY2 = 0, EXY = 0;
	for( i = 0; i < num_seqs; i ++ )
	{
		int j;
		double x = 0;
		for( j = 0; j < 13; j ++ )
			x += gsl_vector_get( beta, j ) * gsl_matrix_get( X, i, j );
			
		double y = seqs[ i ] ->residual;
		
		EX += x;
		EX2 += (x*x);
		EY += y;
		EY2 += (y*y);
		EXY += (x*y);
		
		// calculate the residual
		seqs[ i ] ->target -= x;
		seqs[ i ] ->residual -= x;
	}
	
  	EX /= num_seqs;
	EY /= num_seqs;
	EX2 /= num_seqs;
	EY2 /= num_seqs;
	EXY /= num_seqs;
	
	double sX = sqrt( EX2 - EX * EX );
	double sY = sqrt( EY2 - EY * EY );
	
	double rXY = ( EXY - EX * EY ) / ( sX * sY );
	
	cout << "PCC of nucleotide/dinucleotide composition with target scores: " << rXY << endl << endl;
	
	// release the memory allocated to 
	gsl_matrix_free (X);
	gsl_vector_free (Y);
	gsl_vector_free (beta);
	gsl_matrix_free (cov);
  
}
