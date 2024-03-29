/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
	This file is part of foam-extend.

	foam-extend is free software: you can redistribute it and/or modify it
	under the terms of the GNU General Public License as published by the
	Free Software Foundation, either version 3 of the License, or (at your
	option) any later version.

	foam-extend is distributed in the hope that it will be useful, but
	WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
	crAddressing

Description
	Compressed row matrix addressing

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
	crAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef crAddressing_H
#define crAddressing_H

#include "labelList.H"
#include "refCount.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class crAddressing;
Ostream& operator<<(Ostream&, const crAddressing&);



class crAddressing
:
	public refCount
{
	// Private data

		//- Number of rows
		label nRows_;

		//- Number of columns
		label nCols_;

		//- Row array.
		//  Provides start index for each row, dimensioned to nRows + 1
		labelList rowStart_;

		//- Column index array
		labelList column_;


	// Private Member Functions

		//- Set row sizes
		void setRowSizes(const labelList& count);


public:

	// Constructors

		//- Construct null
		crAddressing();

		//- Construct given size.  Column and coefficients set later
		crAddressing
		(
			const label nRows,
			const label nCols,
			const labelList& nEntries
		);

		//- Construct from components
		crAddressing
		(
			const label nRows,
			const label nCols,
			const labelList& rowStart,
			const labelList& column
		);

		//- Construct as copy
		crAddressing(const crAddressing&);

		//- Construct from Istream
		crAddressing(Istream&);


	// Destructor - default


	// Member Functions

		// Access

			//- Return number of rows
			label nRows() const
			{
				return nRows_;
			}

			//- Return number of columns
			label nCols() const
			{
				return nCols_;
			}

			//- Return number of coefficients
			label nEntries() const
			{
				return column_.size();
			}

			//- Return row array
			const labelList& rowStart() const
			{
				return rowStart_;
			}

			//- Return column array
			const labelList& column() const
			{
				return column_;
			}


		// Edit

			//- Return column array
			labelList& column()
			{
				return column_;
			}


		// Operations

			//- Return transpose addressing
			tmp<crAddressing> T() const;


	// Member Operators

		void operator=(const crAddressing&);


	// IOstream Operators

		friend Ostream& operator<<(Ostream&, const crAddressing&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
