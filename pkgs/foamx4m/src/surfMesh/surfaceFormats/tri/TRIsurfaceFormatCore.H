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
	Foam::fileFormats::TRIsurfaceFormatCore

Description
	Internal class used by the TRIsurfaceFormat

SourceFiles
	TRIsurfaceFormatCore.C

\*---------------------------------------------------------------------------*/

#ifndef TRIsurfaceFormatCore_H
#define TRIsurfaceFormatCore_H

#include "surfaceFormatsCore.H"
#include "triFace.H"

#include "IFstream.H"
#include "Ostream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{


class TRIsurfaceFormatCore
:
	public surfaceFormatsCore
{
	// Private Data

		bool sorted_;

		//- The points supporting the facets
		pointField points_;

		//- The zones associated with the faces
		labelList zoneIds_;

		//- The solid count, in the order of their first appearance
		labelList sizes_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		TRIsurfaceFormatCore(const TRIsurfaceFormatCore&);

		//- Disallow default bitwise assignment
		void operator=(const TRIsurfaceFormatCore&);

		bool read(const fileName&);


public:

	// Constructors

		//- Read from file, filling in the information
		TRIsurfaceFormatCore(const fileName&);


	// Destructor

		~TRIsurfaceFormatCore();


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
			sizes_.clear();
		}

		//- Return full access to the points
		pointField& points()
		{
			return points_;
		}

		//- Return full access to the zones
		labelList& zoneIds()
		{
			return zoneIds_;
		}

		//- The list of zone sizes in the order of their first appearance
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
