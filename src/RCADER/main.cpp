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


// main.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "declarations.h"

extern const char *__rf_file;
extern const char *__filter;
extern const char *__score_file;
extern const char *__fasta_file;
extern const double *__enrichment_fold;
extern const int *__correct_dinuc;
extern const char *__experiment;
extern const char *__output_file;


void welcome_message()
{
	cout << endl
		<< "****************************** RCADER version 1.0 *****************************" << endl
		<< "*.............................................................................*" << endl
		<< "* Copyright 2016 Hamed S. Najafabadi .........................................*" << endl
		<< "*.............................................................................*" << endl
		<< "*.............................................................................*" << endl
		<< endl;
}


int main( int argc, char* argv[] )
{	
	welcome_message();
	
	srand ( time(NULL) );

	if( argc <= 1 )
	// no argument is profided
	{
		print_commandline_format();
		return 0;
	}
	
	if( !read_arguments( argc, argv ) )
		return 1;
		
	print_arguments();

	//******************* open output

	ofstream ofs_PWM;
	if( !open_output( ofs_PWM, __output_file, ".PWM.txt" ) )
		return 1;

	//ofstream ofs_opt_PFM;
	//if( !open_output( ofs_opt_PFM, __output_file, ".opt.PFM.txt" ) )
	//	return 1;

	//ofstream ofs_opt_PFM_MEME;
	//if( !open_output( ofs_opt_PFM_MEME, __output_file, ".opt.PFM.meme.txt" ) )
	//	return 1;

	ofstream ofs_ps;
	if( !open_output( ofs_ps, __output_file, ".ps" ) )
		return 1;

	//ofstream ofs_opt_ps;
	//if( !open_output( ofs_opt_ps, __output_file, ".opt.ps" ) )
	//	return 1;

	ofstream ofs_report;
	if( !open_output( ofs_report, __output_file, ".report.txt" ) )
		return 1;

	ofstream ofs_scores;
	if( !open_output( ofs_scores, __output_file, ".scores.txt" ) )
		return 1;

	//******************* open the random forest file

	s_motif *motifs[ MAX_MOTIFS ];
	int num_motifs = 0;

	cout << "Opening the randomForest output file..." << endl;
	if( !open_rndforest( __rf_file, motifs, &num_motifs, __filter ) )
		return 1;
		
	//******************* open the FASTA file

	s_seq *seqs[ MAX_SEQS ];
	int num_seqs = 0;
	cout << "Opening the FASTA file..." << endl;
	if( !open_FASTA( __fasta_file, seqs, &num_seqs ) )
		return 1;
		
	//******************* open the score file

	cout << "Opening the score file..." << endl;
	if( !open_scores( __score_file, seqs, &num_seqs ) )
		return 1;
	
	//******************* calculations
	cout << "Generating initial PWMs..." << endl;
	if( !generate_PWMs( motifs, num_motifs ) )
		return 1;
	
	cout << "Calculating sequence nucleotide compositions..." << endl << endl;
	calc_nuc_dinuc( seqs, num_seqs );

	if( *__correct_dinuc )
	{
		cout << "Correcting the effect of nucleotide compositions..." << endl;
		correct_for_dinuc( seqs, num_seqs );
	}

	cout << "Optimizing PWMs..." << endl;
	optimize_PWMs( seqs, num_seqs, motifs, num_motifs );
	
	if( *__enrichment_fold > 0 )
	{
		cout << "Normalizing PWMs..." << endl;
		convert_PWMs_to_likelihood( motifs, num_motifs, *__enrichment_fold );
	}			
	
	//******************* write the output files
	cout << "Writing to output..." << endl;
	
	// write all optimized PFMs
	if( !write_scores( ofs_scores, seqs, num_seqs, motifs, num_motifs ) )
		return 1;

	// write all optimized PFMs
	if( !write_PWMs( ofs_PWM, motifs, num_motifs, __experiment, true ) )
		return 1;

	// write the optimized PFM for the best-scoring motif
	//if( !write_PFMs( ofs_opt_PFM, motifs, 1, __experiment, true ) )
	//	return 1;

	// write the optimized PFM for the best-scoring motif in MEME format
	//if( num_motifs )
	//	write_opt_MEME( ofs_opt_PFM_MEME, motifs[ 0 ], __experiment );

	// write the postscript
	write_graphics( ofs_ps, motifs, num_motifs );
	//write_graphics( ofs_opt_ps, motifs, 1 );

	// write report
	write_report( ofs_report, motifs, num_motifs, __experiment );

	cout << endl << "Job finished successfully." << endl;

	return 0;
}
