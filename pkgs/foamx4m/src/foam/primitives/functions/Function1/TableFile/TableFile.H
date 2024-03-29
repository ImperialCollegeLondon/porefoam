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
	Foam::Function1Types::TableFile

Description
	Templated table container data entry where data is read from file.

	\verbatim
		<entryName> tableFile;
		<entryName>Coeffs
		{
			fileName            dataFile;    // name of data file
			outOfBounds         clamp;       // optional out-of-bounds handling
			interpolationScheme linear;      // optional interpolation method
		}
	\endverbatim

	Items are stored in a list of Tuple2's. First column is always stored as
	scalar entries.  Data is read in the form, e.g. for an entry \<entryName\>
	that is (scalar, vector):
	\verbatim
		(
			(0.0 (1 2 3))
			(1.0 (4 5 6))
		);
	\endverbatim


SourceFiles
	TableFile.C

\*---------------------------------------------------------------------------*/

#ifndef TableFile_H
#define TableFile_H

#include "Function1.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{


template<class Type>
class TableFile
:
	public TableBase<Type>
{
	// Private data

		//- File name for csv table (optional)
		fileName fName_;


	// Private Member Functions

		//- Disallow default bitwise assignment
		void operator=(const TableFile<Type>&);


public:

	//- Runtime type information
	TypeName("tableFile");


	// Constructors

		//- Construct from entry name and Istream
		TableFile(const word& entryName, const dictionary& dict);

		//- Copy constructor
		TableFile(const TableFile<Type>& tbl);

		//- Construct and return a clone
		virtual tmp<Function1<Type> > clone() const
		{
			return tmp<Function1<Type> >(new TableFile<Type>(*this));
		}


	//- Destructor
	virtual ~TableFile();


	// I/O

		//- Write in dictionary format
		virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
	#include "TableFile.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
