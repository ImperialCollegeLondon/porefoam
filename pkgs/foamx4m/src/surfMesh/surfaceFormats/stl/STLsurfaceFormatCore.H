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
	Foam::fileFormats::STLsurfaceFormatCore

Description
	Internal class used by the STLsurfaceFormat

SourceFiles
	STLsurfaceFormatCore.C
	STLsurfaceFormatASCII.L

\*---------------------------------------------------------------------------*/

#ifndef STLsurfaceFormatCore_H
#define STLsurfaceFormatCore_H

#include "STLtriangle.H"
#include "triFace.H"
#include "IFstream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{


class STLsurfaceFormatCore
{
	// Private Data

		bool sorted_;

		//- The points supporting the facets
		pointField points_;

		//- The zones associated with the faces
		labelList zoneIds_;

		//- The solid names, in the order of their first appearance
		wordList names_;

		//- The solid count, in the order of their first appearance
		labelList sizes_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		STLsurfaceFormatCore(const STLsurfaceFormatCore&);

		//- Disallow default bitwise assignment
		void operator=(const STLsurfaceFormatCore&);

		//- Determine the file type
		static int detectBINARY(const fileName&);

		//- Read ASCII
		bool readASCII(istream&, const off_t)  {  FatalErrorIn("STLsurfaceFormatCore::readSTLASCII deactivated") << exit(FatalError); return false;  };

		//- Read BINARY
		bool readBINARY(istream&, const off_t);


public:

	// Static Data

		//- The number of bytes in the STL binary header
		static const unsigned int headerSize = 80;


	// Static Member Functions

		//- Write "STL binary file" and number of triangles to stream
		static void writeHeaderBINARY(ostream&, unsigned int);


	// Constructors

		//- Read from file, filling in the information
		STLsurfaceFormatCore(const fileName&);


	// Destructor

		~STLsurfaceFormatCore();


	// Member Functions

		//- File read was already sorted
		bool sorted() const
		{
			return sorted_;
		}

		//- Flush all values
		void clear()
		{
			sorted_ = true;
			points_.clear();
			zoneIds_.clear();
			names_.clear();
			sizes_.clear();
		}

		//- Return full access to the points
		pointField& points()
		{
			return points_;
		}

		//- Return full access to the zoneIds
		labelList& zoneIds()
		{
			return zoneIds_;
		}

		//- The list of solid names in the order of their first appearance
		wordList& names()
		{
			return names_;
		}

		//- The list of solid sizes in the order of their first appearance
		labelList& sizes()
		{
			return sizes_;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
