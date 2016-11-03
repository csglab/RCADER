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

using namespace std;

#include "structures.h"
#include "declarations.h"

#define _EOFCHECK(x)	if( !x || x.eof() ){ cout << "ERROR: Unexpected end of file." << endl << endl; return false; }
#define _PTRCHECK(x)	if( !x ){ cout << "ERROR: Unexpected end of line." << endl; return false; }

static int _line_number; // this local variable will hold the number of lines that are read so far

static char _delimiters[] = { char(9), 0 };

extern const int *__mode;


///////////////////////////////////////////////////////////////////////////////////////////
bool read_scores(
	ifstream &ifs,
	s_seq *seqs[],
	int *num_seqs,
	char this_EOL )
// returns true if successful, false otherwise
{
	// A line cannot be longer than MAX_LINE_LENGTH characters
	char string[ MAX_LINE_LENGTH + 1 ];
	
	// read the file line by line
	_line_number = 0; // reset the line number
	while( true )
	{
		// read the line
		ifs.getline( string, MAX_LINE_LENGTH, this_EOL );
		if( !ifs || ifs.eof() ) // the end of the file has reached
			break;

		_line_number ++; // update the line number that was just read
			
		if( string[ 0 ] == '#' || string[ 0 ] == 0 )
		// this is either a comment line or an empty line
			continue; // ignore the line
		
		// get the sequence name	
		char seq_name[ MAX_STRING_LENGTH ];
		_PTRCHECK( extract_phrase( string, 0, seq_name, _delimiters ) );
		int seq_index = find_seq( seq_name, seqs, *num_seqs );
		if( seq_index >= *num_seqs )
		{
			cout << "WARNING at line " << _line_number << ": Sequence " << seq_name
				<< " was not found in the sequence list." << endl;
			continue;
		}
		
		// get the score
		char score[ MAX_STRING_LENGTH ];
		_PTRCHECK( extract_phrase( string, 1, score, _delimiters ) );
		
		if( *__mode == 1 )
			seqs[ seq_index ] ->target_orig =
			seqs[ seq_index ] ->target =
			seqs[ seq_index ] ->residual = log10( atof( score ) );
		else if( *__mode == 2 )
			seqs[ seq_index ] ->target_orig =
			seqs[ seq_index ] ->target =
			seqs[ seq_index ] ->residual = atof( score );
		else
		{
			cout << "ERROR: Unknown score transformation mode." << endl;
			return false;
		}
	}
	
	bool warning = false;
	// check if all sequences have target values
	int i;
	for( i = 0; i < *num_seqs; i ++ )
		if( seqs[ i ] ->target == _NA )
		{
			if( !warning )
			{
				cout << "WARNING: One or more sequences do not have associated scores and will be removed." << endl;
				warning = true;
			}			
			
			// delete seqs[ i ], shift all other seq objects in the list, and retry the new seqs[ i ]
			
			RELEASE( seqs[ i ] );
			(*num_seqs) --;

			int j;
			for( j = i; j < *num_seqs; j ++ )
				seqs[ j ] = seqs[ j+1 ];
				
			i --;
		}
		
	return true;
}
