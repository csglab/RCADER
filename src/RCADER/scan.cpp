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

#include "declarations.h"

//////////////////////////////////////////////////////////
int _n_index( char nucleotide )
// returns the index of the specified nucleotide
{
	int i;
	for( i = 0; i < NUM_N_LETTERS; i ++ )
		if( nucleotide == __n_letters[ i ] )
			return i;
			
	return -1;
}

//////////////////////////////////////////////////////////
double _get_binding_score(
	char *seq,
	double *PWM[],
	int PWM_width,
	int direction,
	int *max_dir )
// This function returns the score of a particular sequence stretch relative to the given PWM
{
	double score_f = 0;
	double score_r = 0;
	
	int x;
	for( x = 0; x < PWM_width; x ++ )
	{
		int pos_f = x;
		int pos_r = PWM_width - x - 1;
		
		double this_score_f = 0;
		double this_score_r = 0;
		
		int n = _n_index( seq[ x ] );

		if( n < 0 ) // the nucleotide is unknown; thus, the scores will be averaged over this position of the PWM
			for( n = 0; n < NUM_N_LETTERS; n ++ )
			{
				this_score_f += PWM[ n ][ pos_f ] / double(NUM_N_LETTERS);
				this_score_r += PWM[ n ][ pos_r ] / double(NUM_N_LETTERS);
			}
		else // the nucleotide is known
		{
			this_score_f = PWM[ n ][ pos_f ];
			this_score_r = PWM[ NUM_N_LETTERS - n - 1 ][ pos_r ];
		}
		
		// since we're dealing with PWMs, the likelihoods are added
		score_f += this_score_f;
		score_r += this_score_r;
	}

	if( direction == FORWARD_ONLY || // only the forward direction is requested,
		( direction == BOTH_STRANDS && score_f > score_r ) )  // or both strands are requested and the forward score is larger
	{
		*max_dir = FORWARD_ONLY;
		return score_f;
	}
	// else, only the reverse direction is requested, or both directions are requested and the reverse score is larger
	*max_dir = REVERSE_ONLY;
	return score_r;
}

///////////////////////////////////////////////////////////////////////////////////////////
void _scan_sequence(
	s_seq *seq,
	double *PWM[],
	int PWM_width,
	int direction )
// This function scans the given sequence for the instances of the motif
{
	// calculate the max of the scores for this motif over the specified sequence
	// also, store the position of the window that has the maximum score
	double max_score = 0;
	int max_pos = -1;
	int max_dir = -1;
	int i;
	for( i = 0; i <= seq ->seq_length - PWM_width; i ++ )
	{
		// get the score of this position
		int dir = -1;
		double score = _get_binding_score( seq ->seq + i, PWM, PWM_width, direction, &dir );
		
		// update the max score
		if( i == 0 || // either this is the first score
			max_score < score ) // or this is the best score so far
		{
			max_score = score;
			max_pos = i;
			max_dir = dir;
		}
	}
	
	seq ->pwm_score = max_score;
	seq ->max_pos = max_pos;
	seq ->max_dir = max_dir;
}

///////////////////////////////////////////////////////////////////////////////////////////
void scan_sequences(
	s_seq *seqs[],
	int num_seqs,
	double *PWM[],
	int PWM_width,
	int direction )
{
	// scan the sequences
	int i;
	for( i = 0; i < num_seqs; i ++ )
		_scan_sequence( seqs[ i ], PWM, PWM_width, direction );
}
