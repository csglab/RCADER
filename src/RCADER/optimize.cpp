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

#define MAX_ITER			1000
#define MAX_UNCHANGED		20
#define CONVERGENCE_RATE	0.1

////////////////////////////////////////////////////////////////////////////////
int _compare_motif_scores( const void *a, const void *b )
{
	s_motif *motif1 = *( (s_motif **) a );
	s_motif *motif2 = *( (s_motif **) b );

	if( motif1 ->score > motif2 ->score )
		return -1;
	else if( motif1 ->score < motif2 ->score )
		return 1;
	
	return 0;
}

//////////////////////////////////////////////////////////
int _n_index( char nucleotide );

//////////////////////////////////////////////////////////
void _update_PWM(
	s_seq *seqs[],
	int num_seqs,
	double *PWM[],
	int PWM_width,
	int max_PWM_width )
{
	// initialize the objects for GSL regression
	int max_num_cov = 4*max_PWM_width + 1;
	static gsl_matrix *X = gsl_matrix_calloc( num_seqs, max_num_cov );
	static gsl_vector *Y = gsl_vector_alloc( num_seqs );
	static gsl_vector *beta = gsl_vector_alloc( max_num_cov );
	static gsl_matrix *cov = gsl_matrix_alloc( max_num_cov, max_num_cov );
	static gsl_multifit_linear_workspace *wspc = gsl_multifit_linear_alloc( num_seqs, max_num_cov );

	// keep track of the number of instances each covariate is used (i.e. each base is represented in the regression)
	static int *sums = new int[ max_num_cov ];
	memset( sums, 0, sizeof(int)*max_num_cov );

	// set the values for the regression covariates

	int i, n, x;
	for( i = 0; i < num_seqs; i ++ )
	{
		gsl_vector_set( Y, i, seqs[ i ] ->residual );
		
		// preset all the elements of the covariate matrix to (almost) zero, except for the last one, which is set to 1.0
		int j;
		for( j = 0; j < max_num_cov-1; j ++ )
			gsl_matrix_set( X, i, j, drand48()/1000000 );
		gsl_matrix_set( X, i, j, 1 );

		// set the relevant covariates based on the sequence at the previous motif matches
		
		for( x = 0; x < PWM_width; x ++ )
		{
			n = _n_index( seqs[ i ] ->seq[ seqs[ i ] ->max_pos + x ] );
			
			if( n > 0 ) // this base is known
			{
				if( seqs[ i ] ->max_dir == FORWARD_ONLY )
				{
					gsl_matrix_set( X, i,
						x * NUM_N_LETTERS + n,
						1 );
					sums[ x * NUM_N_LETTERS + n ] ++;
				}
				else
				{
					gsl_matrix_set( X, i,
						(PWM_width-x-1) * NUM_N_LETTERS + (NUM_N_LETTERS-n-1),
						1 );
					sums[ (PWM_width-x-1) * NUM_N_LETTERS + (NUM_N_LETTERS-n-1) ] ++;
				}
			}
		}
	}
	
	//for( x = 0; x < max_num_cov; x ++ )
	//	cout << char(9) << sums[x];
	//cout << endl;
	
	// solve the least-square regression
	double chisq;
	gsl_multifit_linear( X, Y, beta, cov, &chisq, wspc );
		
	// update the PWM based on the regression covariates
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		for( x = 0; x < PWM_width; x ++ )
		{
			if( sums[ x * NUM_N_LETTERS + n ] > 4 ) // at least five instances of this nucleotide at this position must have been seen by the regressor
				PWM[ n ][ x ] +=
					gsl_vector_get( beta, x * NUM_N_LETTERS + n ) * CONVERGENCE_RATE;
				
			// cout << gsl_vector_get( beta, x * NUM_N_LETTERS + n ) << endl;
		}
		
	//for( x = 0; x < max_num_cov; x ++ )
	//	cout << char(9) << gsl_vector_get( beta, x );
	//cout << endl;


}

//////////////////////////////////////////////////////////
void _calc_cor(
	s_seq *seqs[],
	int num_seqs,
	double *cor,
	double *p,
	double *PWM[],
	int PWM_width )
{
	// find the regression statistics between motif PWM scores and target values
	// update the residual values
	// scale the PWM according to regression slope
	
	double EX = 0, EX2 = 0, EY = 0, EY2 = 0, EXY = 0;
	int i;
	
	// find the regression statistics and Pearson correlation
	for( i = 0; i < num_seqs; i ++ )
	{
		double x = seqs[ i ] ->pwm_score;
		double y = seqs[ i ] ->target;
		
		EX += x;
		EX2 += (x*x);
		EY += y;
		EY2 += (y*y);
		EXY += (x*y);
	}
	
	EX /= num_seqs;
	EY /= num_seqs;
	EX2 /= num_seqs;
	EY2 /= num_seqs;
	EXY /= num_seqs;
	
	double sX = sqrt( EX2 - EX * EX );
	double sY = sqrt( EY2 - EY * EY );
	
	double rXY = ( EXY - EX * EY ) / ( sX * sY );
	
	double b = rXY * sY / sX; // the slope of regression
	double a = EY - b * EX; // the intercept of regression
	
	// calculate regression residuals
	for( i = 0; i < num_seqs; i ++ )
	{
		double x = seqs[ i ] ->pwm_score;
		double y = seqs[ i ] ->target;

		seqs[ i ] ->residual = y - ( a + b*x );
	}
	
	// scale the PWM
	int n, x;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		for( x = 0; x < PWM_width; x ++ )
			PWM[ n ][ x ] *= b;
			
	// store the statistics required for selecting the best motif
	*cor = rXY;
	*p = normal_cum_p( 0.5 * log( (1+rXY)/(1-rXY) ) * sqrt( num_seqs-3 ) );
}

//////////////////////////////////////////////////////////
int _optimize_PWM(
	s_seq *seqs[],
	int num_seqs,
	double *ini_PWM[],
	double *opt_PWM[],
	int PWM_width,
	int max_PWM_width,
	double *ini_cor,
	double *ini_p,
	double *opt_cor,
	double *opt_p )
// returns 1 if optimized, 0 otherwise
{
	*ini_cor = -1;
	*ini_p = 1;
	*opt_cor = -1;
	*opt_p = 1;

	// initialize internal variables and buffers
	double best_cor = -1;
	double best_cor_p = -1;
	double best_score = -1;
	static double *best_PWM[ NUM_N_LETTERS ];
	static double *PWM[ NUM_N_LETTERS ];
	int n, x;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
	{
		if( !best_PWM[ n ] )
			best_PWM[ n ] = new double[ max_PWM_width ];
		memset( best_PWM[ n ], 0, sizeof(double) * max_PWM_width  );
		if( !PWM[ n ] )
			PWM[ n ] = new double[ max_PWM_width ];
		memset( PWM[ n ], 0, sizeof(double) * max_PWM_width  );
	}
	
	// initialize the PWM as the initial PWM of this motif
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		memcpy( PWM[ n ] , ini_PWM[ n ], sizeof(double)*PWM_width );
		
	// initialize the residual values to be the same as target values
	int i;
	for( i = 0; i < num_seqs; i ++ )
		seqs[ i ] ->residual = seqs[ i ] ->target;
		
	cout << "Cor(optimized_PWM,target)" << char(9) << "P(optimized_PWM,target)" << char(9)
		<< "Cor(optimized_PWM,initial_PWM)" << char(9) << "P(optimized_PWM,initial_PWM)" << endl;

	int counter = 0;
	for( i = 0; i < MAX_ITER; i ++ )
	{
		// scan and the sequences
		scan_sequences( seqs, num_seqs, PWM, PWM_width, BOTH_STRANDS );
		
		// calculate PCC and related statistics (including residuals of regression)
		// also scale the PWM according to the slope of regression
		double cor = -1;
		double p_cor = 1;
		_calc_cor( seqs, num_seqs, &cor, &p_cor, PWM, PWM_width );

		// calculate the correlation of this motif with the initial motif
		double mean_correl = 0, stdev_correl = 1;
		double correl = get_correlation( ini_PWM, PWM, PWM_width, &mean_correl, &stdev_correl );
		double z_correl = ( correl - mean_correl ) / stdev_correl;
		double p_correl = normal_cum_p( z_correl );
		
		cout << cor << char(9) << p_cor << char(9) << correl << char(9) << p_correl << endl;
		
		double score = cor;

		// store the initial PCC and its p-value
		if( i == 0 )
		{
			*ini_cor = cor;
			*ini_p = p_cor;
		}
		
		double prev_best_score = best_score;
		if( i == 0 || best_score < score ) // this is the best PWM so far
		{
			best_score = score;
			best_cor = cor;
			best_cor_p = p_cor;
			for( n = 0; n < NUM_N_LETTERS; n ++ )
				memcpy( best_PWM[ n ], PWM[ n ], sizeof(double)*PWM_width );
		}
			
		if( p_cor > 0.001 )
		{
			cout << "Pearson correlation is insignificant. Optimization halted." << endl;
			return 0;
		}

		if( correl < 0.3 )
		{
			cout << "Motif diverged from the seed. Optimization halted." << endl;
			return 0;
		}
		
		// check whether any improvement is achieved
		if( score < prev_best_score + 0.001 ) // there must be at least a minuscule improvement in score
		{
			counter ++;
			
			if( counter >= MAX_UNCHANGED ) // for this many times in a row, the score did not improve significantly
				break;
		}
		else
			counter = 0;
		
		if( i < MAX_ITER - 1 ) // there will be at least one more iteration
		// update the PWM based on residual scores
			_update_PWM( seqs, num_seqs, PWM, PWM_width, max_PWM_width );
	}
	
	// store the best PFM
	*opt_cor = best_cor;
	*opt_p = best_cor_p;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		memcpy( opt_PWM[ n ], best_PWM[ n ], sizeof(double)*PWM_width );

	return 1;
}


//////////////////////////////////////////////////////////
void _optimize_PWMs(
	s_seq *seqs[],
	int num_seqs,
	s_motif *motifs[],
	int num_motifs )
{
	// find the maximum number of zinc fingers per motif
	int max_PWM_width = -1;
	int i;
	for( i = 0; i < num_motifs; i ++ )
		if( max_PWM_width < motifs[ i ] ->PWM_width )
			max_PWM_width = motifs[ i ] ->PWM_width;
			
	// initialize the PWM
	double *PWM[ NUM_N_LETTERS ];
	int n;
	for( n = 0; n < NUM_N_LETTERS; n ++ )
		PWM[ n ] = new double[ max_PWM_width ];
		
	// optimize each PWM
	for( i = 0; i < num_motifs; i ++ )
	{
		cout << "Optimizing motif " << motifs[ i ] ->name << " ..." << endl;
		motifs[ i ] ->score = -2;

		// optimize this PWM
		motifs[ i ] ->PWM_optimized =
			_optimize_PWM( seqs, num_seqs,
				motifs[ i ] ->PWM, PWM, motifs[ i ] ->PWM_width,
				max_PWM_width,
				&motifs[ i ] ->ini_cor, &motifs[ i ] ->ini_p,
				&motifs[ i ] ->opt_cor, &motifs[ i ] ->opt_p );
		
		// store the optimized PWM, and calculate correlation with the original PWM
		if( motifs[ i ] ->PWM_optimized )
		{
			for( n = 0; n < NUM_N_LETTERS; n ++ )
			{
				motifs[ i ] ->opt_PWM[ n ] = new double[ motifs[ i ] ->PWM_width ];
				memcpy( motifs[ i ] ->opt_PWM[ n ] , PWM[ n ], sizeof(double)*motifs[ i ] ->PWM_width );
			}
			
			double mean_correl = 0;
			double stdev_correl = 1;
			motifs[ i ] ->correl = get_correlation( motifs[ i ] ->PWM, PWM,
				motifs[ i ] ->PWM_width, &mean_correl, &stdev_correl );
			double z_correl = ( motifs[ i ] ->correl - mean_correl ) / stdev_correl;
			motifs[ i ] ->p_correl = normal_cum_p( z_correl );
			
			motifs[ i ] ->score = motifs[ i ] ->opt_cor;

		}
		
		cout << endl << endl;
	}
	
	// sort the motifs by descending order of their scores
	qsort( motifs, num_motifs, sizeof(s_motif*), _compare_motif_scores );
	

	///// prepare the target values for the next round of optimization
	
	if( motifs[ 0 ] ->PWM_optimized ) // there is at least one optimized motif
	{
		// scan all the sequences with the best motif
		scan_sequences( seqs, num_seqs, motifs[ 0 ] ->opt_PWM, motifs[ 0 ] ->PWM_width, BOTH_STRANDS );
		
		// calculate the residuals
		double cor = -1;
		double p_cor = 1;
		_calc_cor( seqs, num_seqs, &cor, &p_cor, motifs[ 0 ] ->opt_PWM, motifs[ 0 ] ->PWM_width );
		
		cout << endl << ".........." << endl;
		cout << "The best optimized motif is " << motifs[ 0 ] ->name
			<< ", PCC=" << cor << ", P=" << p_cor << endl;
		cout << endl << ".........." << endl << endl;
		
		// replace the target values with residuals, for the next round of optimization
		for( i = 0; i < num_seqs; i ++ )
			seqs[ i ] ->target = seqs[ i ] ->residual;
	}
}

//////////////////////////////////////////////////////////
void optimize_PWMs(
	s_seq *seqs[],
	int num_seqs,
	s_motif *motifs[],
	int num_motifs )
{
	int i;
	for( i = 0; i < num_motifs; i ++ )
	{
		cout << endl << endl << "..... Optimization round " << i+1 << " ....." << endl;
		_optimize_PWMs( seqs, num_seqs, motifs + i, num_motifs - i );

		if( !motifs[ i ] ->PWM_optimized ) // no motif was optimized
			break;
	}
}
