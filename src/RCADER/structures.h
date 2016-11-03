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

#ifndef _H_STRUCTURES // this is to make sure that this file will not be included twice
//*******************************************************************
#define _H_STRUCTURES // this macro indicates that this file is now included

#ifndef _NAMESPACE
using namespace std;
#define _NAMSPACE
#endif


#include <string.h>
#include <stdlib.h>

// ************************************************ macro definitions

// constants
#define MAX_MOTIFS				100000
#define MAX_LINE_LENGTH			100000
#define MAX_STRING_LENGTH		1000
#define MAX_ZF_PER_GENE			100
#define MAX_SEQS				500000
#define MAX_SEQ_LENGTH			10000

#define FORWARD_ONLY			0
#define REVERSE_ONLY			1
#define BOTH_STRANDS			2

#define _NA						(-10000)

// macros
#define RELEASE(x)		do{ if(x) delete (x); }while(false) // this do-while loop is called only once
#define _CALL(x)		do{ if(!(x)) return false; }while(false)
#define _COPY_STR(x,y)	do{ x=new char[strlen(y)+1]; strcpy((x),(y)); }while(false)
#define MAX(x,y)		(((x)>(y))?(x):(y))
#define MIN(x,y)		(((x)<(y))?(x):(y))

// ********************************************** typedef definitions

typedef unsigned char BYTE;

// ********************************************** local variables
static char __aa_letters[] = "ACDEFGHIKLMNPQRSTVWY-";
#define NUM_AA_LETTERS	21

static char __n_letters[] = "ACGT";
#define NUM_N_LETTERS	4


// ******************************************** structure definitions

// this structure will hold the information for the c2h2 zinc fingers
struct s_c2h2
{
	s_c2h2()
	{
		int i;
		for( i = 0; i < NUM_N_LETTERS; i ++ )
			memset( PWM, 0, sizeof(double) * 4 );
			
		weight = 1;
	}

	double PWM[ NUM_N_LETTERS ][ 4 ];
	
	double weight; // the weight that this zf is given when computing the PWM of the whole protein
};

// this structure will hold the information for the input motifs
struct s_motif
{
	s_motif()
	{
		// initialize all the variables, setting them to zero
		
		name = NULL; // no memory allocated yet
		
		protein_name = NULL;
		first_zf_in_protein = -1;
		last_zf_in_protein = -1;

		memset( zfs, 0, sizeof(s_c2h2*) * MAX_ZF_PER_GENE ); // no memory allocated yet
		num_zfs = 0;
		
		memset( PWM, 0, sizeof(double*) * NUM_N_LETTERS );

		PWM_optimized = 0;
		memset( opt_PWM, 0, sizeof(double*) * NUM_N_LETTERS );
		ini_cor = -1;
		ini_p = 1;
		opt_cor = -1;
		opt_p = 1;
		
		correl = 0;
		p_correl = 1;
	}
	~s_motif()
	{
		// release all the allocated memory
	
		RELEASE( name );
		RELEASE( protein_name );
		
		int n;
		for( n = 0; n < NUM_N_LETTERS; n ++ )
		{
			RELEASE( PWM[ n ] );
			RELEASE( opt_PWM[ n ] );
		}
	}

	char *name; // the name of this motif
	
	char *protein_name; // the name of the protein that this motif belongs to
	int first_zf_in_protein; // the index of the first ZF in the original protein
	int last_zf_in_protein; // the index of the last ZF in the original protein

	s_c2h2 *zfs[ MAX_ZF_PER_GENE ]; // the zinc fingers included in this motif
	int num_zfs; // the number of zinc fingers
	
	double *PWM[ NUM_N_LETTERS ];
	int PWM_width;
	
	int PWM_optimized; // whether the PWM was optimized for this motif
	double *opt_PWM[ NUM_N_LETTERS ]; // the PWM after optimization
	double ini_cor; // the Pearson correlation with scores before optimization
	double ini_p; // the p-value of the above correlation
	double opt_cor; // the Pearson correlation with scores after optimization
	double opt_p; // the p-value of the above correlation
	
	double correl; // the correlation between the initial and optimized PWMs
	double p_correl; // the p-value association with the above correlation
	
	double score; // the final score for this motif (same as opt_cor)
	
};

// this structure will hold the information for the input sequences
struct s_seq
{
	s_seq()
	{
		// initialize all the variables, setting them to zero
		
		name = NULL; // no memory allocated yet
		seq = NULL; // no memory allocated yet
		seq_length = 0; // the initial length of the sequence is zero
		
		target = _NA;
		residual = _NA;
		
		pwm_score = 0;
		max_pos = -1;
		max_dir = -1;
	}
	~s_seq()
	{
		// release all the allocated memory
	
		RELEASE( name );
		RELEASE( seq );
	}

	char *name; // the name of this sequence, as read from the input FASTA file

	char *seq; // the sequence, as read from the input FASTA file
	int seq_length; // the length of the nucleotide sequence
	
	double nuc[ 2 ]; // the nucleotide frequencies of this sequence
	double dinuc[ 10 ]; // the dinucleotide frequencies of this sequence
	
	double target; // the target value that needs to be predicted by the motifs
	double residual; // the residual value after regression of motif scores vs. target value

	double target_orig; // a copy of the original target score, as read from the file
	
	double pwm_score; // the score for scanning this sequence with the active PFM
	int max_pos; // the position for the best hit
	int max_dir; // the direction for the best hit
};


//*******************************************************************
#endif // this is to make sure that this file will not be included twice
