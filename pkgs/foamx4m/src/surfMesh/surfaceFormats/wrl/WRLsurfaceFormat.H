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
	Foam::fileFormats::WRLsurfaceFormat

Description
	Provide a means of writing VRML97 (wrl) format.

SourceFiles
	WRLsurfaceFormat.C

\*---------------------------------------------------------------------------*/

#ifndef WRLsurfaceFormat_H
#define WRLsurfaceFormat_H

#include "MeshedSurface.H"
#include "MeshedSurfaceProxy.H"
#include "UnsortedMeshedSurface.H"
#include "WRLsurfaceFormatCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{


template<class Face>
class WRLsurfaceFormat
:
	public MeshedSurface<Face>,
	public WRLsurfaceFormatCore

{
	// Private Member Functions

		//- Disallow default bitwise copy construct
		WRLsurfaceFormat(const WRLsurfaceFormat<Face>&);

		//- Disallow default bitwise assignment
		void operator=(const WRLsurfaceFormat<Face>&);

public:

	// Constructors

		//- Construct null
		WRLsurfaceFormat();


	// Destructor

		virtual ~WRLsurfaceFormat()
		{}


	// Member Functions

		//- Write surface mesh components by proxy
		static void write(const fileName&, const MeshedSurfaceProxy<Face>&);

		//- Write object file
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
#	include "WRLsurfaceFormat.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
