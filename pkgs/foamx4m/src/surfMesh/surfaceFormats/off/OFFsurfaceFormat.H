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
	Foam::fileFormats::OFFsurfaceFormat

Description
	Provide a means of reading/writing Geomview OFF polyList format.


See Also
	The <a href="http://www.geoview.org">Geoview</a>
	file format information:
	http://www.geomview.org/docs/html/OFF.html#OFF

Note
	When reading, the optional @a colorspec is ignored.
	When writing, it is set to the zone number (integer).

SourceFiles
	OFFsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef OFFsurfaceFormat_H
#define OFFsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{


template<class Face>
class OFFsurfaceFormat
:
	public MeshedSurface<Face>
{
	// Private Member Functions

		//- Disallow default bitwise copy construct
		OFFsurfaceFormat(const OFFsurfaceFormat&);

		//- Disallow default bitwise assignment
		void operator=(const OFFsurfaceFormat&);


public:

	// Constructors

		//- Construct from file name
		OFFsurfaceFormat(const fileName&);


	// Selectors

		//- Read file and return surface
		static autoPtr<MeshedSurface<Face> > New(const fileName& name)
		{
			return autoPtr<MeshedSurface<Face> >
			(
				new OFFsurfaceFormat(name)
			);
		}


	//- Destructor
	virtual ~OFFsurfaceFormat()
	{}


	// Member Functions

		//- Write surface mesh components by proxy
		static void write(const fileName&, const MeshedSurfaceProxy<Face>&);

		//- Read from file
		virtual bool read(const fileName&);

		//- Write object
		virtual void write(const fileName& name) const
		{
			write(name, MeshedSurfaceProxy<Face>(*this));
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#	include "OFFsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
