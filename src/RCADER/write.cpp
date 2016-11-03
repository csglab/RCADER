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
#include "graphics.h"

#define _PTRCHECK(x)	if( !x ){ cout << "ERROR: Unexpected end of string." << endl; return false; }
static char _orsign[] = { '|', 0 };

void _scan_sequence( s_seq *seq, double *PWM[], int PWM_width, int direction );

extern const char *__nuc1;
extern const char *__nuc2;
extern const char *__dinuc1[];
extern const char *__dinuc2[];


///////////////////////////////////////////////////////////////////////////////////////////
bool write_scores(
	ofstream &ofs,
	s_seq *seqs[],
	int num_seqs,
	s_motif *motifs[],
	int num_motifs )
{
	// write the header
	ofs << "SequenceID" << char(9) << "TargetScore";
	int i;
	for( i = 0; i < 2; i ++ )
		ofs << char(9) << __nuc1[ i ] << "/" << __nuc2[ i ];
	for( i = 0; i < 10; i ++ )
		ofs << char(9) << __dinuc1[ i ] << "/" << __dinuc2[ i ];
	
	for( i = 0; i < num_motifs; i ++ )
		if( motifs[ i ] ->PWM_optimized ) // this motif has an optimized version
			ofs << char(9) << motifs[ i ] ->name << "|opt|pwm_score"
				<< char(9) << motifs[ i ] ->name << "|opt|pwm_position"
				<< char(9) << motifs[ i ] ->name << "|opt|pwm_direction";
	ofs << endl;
	
	for( i = 0; i < num_seqs; i ++ )
	{
		ofs << seqs[ i ] ->name << char(9) << seqs[ i ] ->target_orig;

		int j;
		for( j = 0; j < 2; j ++ )
			ofs << char(9) << seqs[ i ] ->nuc[ j ];
		for( j = 0; j < 10; j ++ )
			ofs << char(9) << seqs[ i ] ->dinuc[ j ];

		for( j = 0; j < num_motifs; j ++ )
			if( motifs[ j ] ->PWM_optimized ) // this motif has an optimized version
			{
				_scan_sequence( seqs[ i ], motifs[ j ] ->opt_PWM, motifs[ j ] ->PWM_width, BOTH_STRANDS );
				ofs << char(9) << seqs[ i ] ->pwm_score
					<< char(9) << seqs[ i ] ->max_pos
					<< char(9) << ( ( seqs[ i ] ->max_dir == FORWARD_ONLY )? "forward":"reverse" );
			}
		ofs << endl;
	}
	
	return true;
}


///////////////////////////////////////////////////////////////////////////////////////////
bool write_PFMs(
	ofstream &ofs,
	s_motif *motifs[],
	int num_motifs,
	const char *experiment_name,
	bool opt_only )
{
	// write the PFM for each motif
	int i;
	for( i = 0; i < num_motifs; i ++ )
		if( motifs[ i ] ->PWM_optimized ) // this motif has an optimized version
		{
			////// write the initial PFM
			if( !opt_only )
			{
				ofs << "TF" << char(9) << "Unknown" << endl
					<< "TF Name" << char(9) << "Unknown" << endl
					<< "Gene" << char(9) << motifs[ i ] ->protein_name << endl
					<< "Motif" << char(9) << motifs[ i ] ->name << endl
					<< "Family" << char(9) << "C2H2 ZF" << endl
					<< "Species" << char(9) << "Unknown" << endl;
	
				ofs << "Pos";
				int n;
				for( n = 0; n < NUM_N_LETTERS; n ++ )
					ofs << char(9) << __n_letters[ n ];
				ofs << endl;

				int j;
				for( j = 0; j < motifs[ i ] ->PWM_width; j ++ )
				{
					ofs << j+1;
					for( n = 0; n < NUM_N_LETTERS; n ++ )
						ofs << char(9) << pow( 10, motifs[ i ] ->PWM[ n ][ j ] ) / NUM_N_LETTERS;
					ofs << endl;
				}

				ofs << endl << endl;
			}

			////// write the optimized PFM
	
			ofs << "TF" << char(9) << "Unknown" << endl
				<< "TF Name" << char(9) << "Unknown" << endl
				<< "Gene" << char(9) << motifs[ i ] ->protein_name << endl
				<< "Motif" << char(9) << motifs[ i ] ->name << "|opt" << experiment_name << endl
				<< "Family" << char(9) << "C2H2 ZF" << endl
				<< "Species" << char(9) << "Unknown" << endl;
		
			ofs << "Pos";
			int n;
			for( n = 0; n < NUM_N_LETTERS; n ++ )
				ofs << char(9) << __n_letters[ n ];
			ofs << endl;
	
			int j;
			for( j = 0; j < motifs[ i ] ->PWM_width; j ++ )
			{
				ofs << j+1;
				for( n = 0; n < NUM_N_LETTERS; n ++ )
					ofs << char(9) << pow( 10, motifs[ i ] ->opt_PWM[ n ][ j ] ) / NUM_N_LETTERS;
				ofs << endl;
			}
	
			ofs << endl << endl;
		}
	
	return true;
}


///////////////////////////////////////////////////////////////////////////////////////////
void write_opt_MEME(
	ofstream &ofs,
	s_motif *motif,
	const char *experiment_name )
{
	if( !motif ->PWM_optimized )
		return; // if this motif is not optimized, do not write it
	
	// write the PWM for this motif in a format that is recognized by the MEME suite
	
	ofs << "MEME version 4" << endl
		<< endl
		<< "ALPHABET= ACGT" << endl
		<< endl
		<< "strands: + -" << endl
		<< endl;
	
	ofs << "MOTIF " << motif ->name << "|opt" << experiment_name << endl
		<< "letter-probability matrix: alength= " << NUM_N_LETTERS
			<< " w= " << motif ->PWM_width
			<< " nsites= 10000" << endl;

	int i;
	for( i = 0; i < motif ->PWM_width; i ++ )
	{
		int n;
		for( n = 0; n < NUM_N_LETTERS; n ++ )
		{
			if( n > 0 )
				ofs << char(9);
			ofs << pow( 10, motif ->opt_PWM[ n ][ i ] ) / NUM_N_LETTERS;
		}
		ofs << endl;
	}
			
	ofs << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
void write_report(
	ofstream &ofs,
	s_motif *motifs[],
	int num_motifs,
	const char *experiment_name )
{
	
	// produce a report
	ofs << "EXPERIMENT" << char(9) << "MOTIF" << char(9) << "OPTIMIZED" << char(9)
		<< "INITIAL_AUC" << char(9) << "INITIAL_P" << char(9)
		<< "OPTIMIZED_AUC" << char(9) << "OPTIMIZED_P" << char(9)
		<< "CORRELATION" << char(9) << "CORRELATION_P" << endl;
		
	int i;
	for( i = 0; i < num_motifs; i ++ )
		ofs << experiment_name << char(9)
			<< motifs[ i ] ->name << char(9)
			<< motifs[ i ] ->PWM_optimized << char(9)
			<< motifs[ i ] ->ini_cor << char(9)
			<< motifs[ i ] ->ini_p << char(9)
			<< motifs[ i ] ->opt_cor << char(9)
			<< motifs[ i ] ->opt_p << char(9)
			<< motifs[ i ] ->correl << char(9)
			<< motifs[ i ] ->p_correl << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
void write_graphics(
	ofstream &ofs,
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

	// draw the logos
	char string[ MAX_STRING_LENGTH + 1 ];
	for( i = 0; i < num_motifs * 2; i ++ )
	{
		// copy the information for the relevant PWM
		int motif_index = i / 2;
		if( !motifs[ motif_index ] ->PWM_optimized )
			break;
		int PWM_width = motifs[ motif_index ] ->PWM_width;
		int n, x;
		for( n = 0; n < NUM_N_LETTERS; n ++ )
			for( x = 0; x < PWM_width; x ++ )
				if( i%2 == 0 )
					PWM[ n ][ x ] = motifs[ motif_index ] ->PWM[ n ][ x ];
				else
					PWM[ n ][ x ] = motifs[ motif_index ] ->opt_PWM[ n ][ x ];

		// if this is the beginning of the page, add margin					
		if( i%10 == 0 )
			//ps_translate( ofs, LEFT_MARGIN, PAGE_HEIGHT - TOP_MARGIN - LOGO_HEIGHT );
			ps_translate( ofs, PAGE_WIDTH / 5, PAGE_HEIGHT - TOP_MARGIN - LOGO_HEIGHT );

		// write the motif name
		ps_translate( ofs, 0, LOGO_HEIGHT / 2 + 10 );
		ps_font( ofs, "Helvetica-Bold", BASE_FONT_SIZE * 0.75 );
		if( i%2 == 0 )
			ps_write_string( ofs, motifs[ motif_index ] ->name, 0, 0, 0 );
		ps_translate( ofs, 0, -LOGO_HEIGHT / 2 - 10 );

		// write the motif
		if( i%2 == 0 )
			ps_translate( ofs, 0, -LOGO_HEIGHT * 0.6 );
		write_motif( ofs, PWM, PWM_width, BASE_FONT_SIZE / 2, 0 );

		// draw a box, with the color depending on whether this is the seed or optimized motif
		if( i%2 == 0 )
			ps_box_path( ofs,
				-COLUMN_WIDTH / 2, -LOGO_HEIGHT * 0.05,
				PWM_width * COLUMN_WIDTH + COLUMN_WIDTH,
				LOGO_HEIGHT * 1.2,
				0.7, 0.7, 0.7 );
		else
			ps_box_path( ofs,
				-COLUMN_WIDTH / 2, -LOGO_HEIGHT * 0.05,
				PWM_width * COLUMN_WIDTH + COLUMN_WIDTH,
				LOGO_HEIGHT * 1.2,
				0, 0, 0 );
				
		// draw a line connecting the logo pairs, and write the correlation statistics
		if( i%2 == 0 )
		{
			ps_line( ofs,
				PWM_width * COLUMN_WIDTH + 5*COLUMN_WIDTH, LOGO_HEIGHT /2,
				PWM_width * COLUMN_WIDTH + 5*COLUMN_WIDTH, LOGO_HEIGHT /2 - BASE_FONT_SIZE * FONT_UNIT_HEIGHT * 0.867 - LOGO_HEIGHT - 5,
				0.5, 0, 0, 0 );
			ps_line( ofs,
				PWM_width * COLUMN_WIDTH + 5*COLUMN_WIDTH, LOGO_HEIGHT /2,
				PWM_width * COLUMN_WIDTH + 4*COLUMN_WIDTH, LOGO_HEIGHT /2,
				0.5, 0, 0, 0 );
			ps_line( ofs,
				PWM_width * COLUMN_WIDTH + 5*COLUMN_WIDTH, LOGO_HEIGHT /2 - BASE_FONT_SIZE * FONT_UNIT_HEIGHT * 0.867 - LOGO_HEIGHT - 5,
				PWM_width * COLUMN_WIDTH + 4*COLUMN_WIDTH, LOGO_HEIGHT /2 - BASE_FONT_SIZE * FONT_UNIT_HEIGHT * 0.867 - LOGO_HEIGHT - 5,
				0.5, 0, 0, 0 );

			// Pearson correlation
			ps_translate( ofs, PWM_width * COLUMN_WIDTH + 6*COLUMN_WIDTH,
				LOGO_HEIGHT /2 + (- BASE_FONT_SIZE * FONT_UNIT_HEIGHT * 0.867 - LOGO_HEIGHT - 5 ) /2 + LOGO_HEIGHT/6 );
			ps_font( ofs, "Helvetica", BASE_FONT_SIZE * 0.75 );
			sprintf( string, "Pearson correlation: %0.2f", motifs[ motif_index ] ->correl );
			ps_write_string( ofs, string, 0, 0, 0 );

			ps_translate( ofs, 0, -LOGO_HEIGHT/3 );
			ps_font( ofs, "Helvetica", BASE_FONT_SIZE * 0.75 );
			sprintf( string, "p-value: %0.0e", motifs[ motif_index ] ->p_correl );
			ps_write_string( ofs, string, 0, 0, 0 );
			ps_translate( ofs, 0, LOGO_HEIGHT/3 );
			
			ps_translate( ofs, - (PWM_width * COLUMN_WIDTH + 6*COLUMN_WIDTH ),
				- ( LOGO_HEIGHT /2 + (- BASE_FONT_SIZE * FONT_UNIT_HEIGHT * 0.867 - LOGO_HEIGHT - 5 ) /2 + LOGO_HEIGHT/6 ) );
		}
				
		////////// write the motif info and statistics
		
		// whether it's a seed or optimized motif
		ps_translate( ofs, -PAGE_WIDTH/5 + LEFT_MARGIN, LOGO_HEIGHT*0.8 );
		ps_font( ofs, "Helvetica-Bold", BASE_FONT_SIZE * 0.75 );
		if( i%2 == 0 )
			ps_write_string( ofs, "Seed", 0, 0, 0 );
		else
			ps_write_string( ofs, "Optimized", 0, 0, 0 );

		// Motif AUC
		ps_translate( ofs, 0, -LOGO_HEIGHT/3 );
		ps_font( ofs, "Helvetica", BASE_FONT_SIZE * 0.75 );
		if( i%2 == 0 )
			sprintf( string, "PCC: %0.2f", motifs[ motif_index ] ->ini_cor );
		else
			sprintf( string, "PCC: %0.2f", motifs[ motif_index ] ->opt_cor );

		ps_write_string( ofs, string, 0, 0, 0 );

		// Motif p-value
		ps_translate( ofs, 0, -LOGO_HEIGHT/3 );
		ps_font( ofs, "Helvetica", BASE_FONT_SIZE * 0.75 );
		if( i%2 == 0 )
			sprintf( string, "p-value: %0.0e", motifs[ motif_index ] ->ini_p );
		else
			sprintf( string, "p-value: %0.0e", motifs[ motif_index ] ->opt_p );

		ps_write_string( ofs, string, 0, 0, 0 );

		ps_translate( ofs, +PAGE_WIDTH/5 - LEFT_MARGIN, -LOGO_HEIGHT*0.8 + 2*LOGO_HEIGHT/3 );
		

		// make room for the next logo
		ps_translate( ofs, 0, -BASE_FONT_SIZE * FONT_UNIT_HEIGHT * 0.867 - LOGO_HEIGHT - 5 );

		// if this is the optimized PWM, make room for the title of the next page
		if( i%2 == 1 )
			ps_translate( ofs, 0, -15 );
	
		if( i%10 == 9 )
			ofs << "showpage" << endl;
	}

	if( i%10 != 0 )
		ofs << "showpage" << endl;
}
