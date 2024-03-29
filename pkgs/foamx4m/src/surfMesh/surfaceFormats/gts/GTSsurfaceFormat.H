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
	Foam::fileFormats::GTSsurfaceFormat

Description
	Provide a means of reading/writing GTS format.
	The output is never sorted by zone and is only written if it consists
	entirely of triangles.

SourceFiles
	GTSsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef GTSsurfaceFormat_H
#define GTSsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{


template<class Face>
class GTSsurfaceFormat
:
	public UnsortedMeshedSurface<Face>
{
	// Private Member Functions

		//- Disallow default bitwise copy construct
		GTSsurfaceFormat(const GTSsurfaceFormat<Face>&);

		//- Disallow default bitwise assignment
		void operator=(const GTSsurfaceFormat<Face>&);


public:

	// Constructors

		//- Construct from file name
		GTSsurfaceFormat(const fileName&);


	// Selectors

		//- Read file and return surface
		static autoPtr<UnsortedMeshedSurface<Face> > New(const fileName& name)
		{
			return autoPtr<UnsortedMeshedSurface<Face> >
			(
				new GTSsurfaceFormat<Face>(name)
			);
		}


	//- Destructor
	virtual ~GTSsurfaceFormat()
	{}


	// Member Functions

		//- Write MeshedSurface
		static void write(const fileName&, const MeshedSurface<Face>&);

		//- Write UnsortedMeshedSurface, the output remains unsorted
		static void write(const fileName&, const UnsortedMeshedSurface<Face>&);

		//- Read from file
		virtual bool read(const fileName&);

		//- Write object
		virtual void write(const fileName& name) const
		{
			write(name, *this);
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "GTSsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
