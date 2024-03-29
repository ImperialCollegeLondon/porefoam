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
	Foam::MeshedSurfaceIOAllocator

Description
	A helper class for storing points, faces and zones with IO capabilities.

SourceFiles
	MeshedSurfaceIOAllocator.C

\*---------------------------------------------------------------------------*/

#ifndef MeshedSurfaceIOAllocator_H
#define MeshedSurfaceIOAllocator_H

#include "pointIOField.H"
#include "faceIOList.H"
#include "surfZoneIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


//- A helper class for storing points, faces and zones
class MeshedSurfaceIOAllocator
{
	// Private data

		//- Points
		pointIOField   points_;

		//- Faces
		faceIOList	 faces_;

		//- Surface zones
		surfZoneIOList zones_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		MeshedSurfaceIOAllocator(const MeshedSurfaceIOAllocator&);

		//- Disallow default bitwise assignment
		void operator=(const MeshedSurfaceIOAllocator&);


public:

	// Constructors

		//- Read construct from IOobjects
		MeshedSurfaceIOAllocator
		(
			const IOobject& ioPoints,
			const IOobject& ioFaces,
			const IOobject& ioZones
		);

		//- Construct from IOobjects, copying components
		MeshedSurfaceIOAllocator
		(
			const IOobject& ioPoints,
			const pointField& points,
			const IOobject& ioFaces,
			const faceList& faces,
			const IOobject& ioZones,
			const surfZoneList& zones
		);

		//- Construct from IOobjects, possibly transferring components
		MeshedSurfaceIOAllocator
		(
			const IOobject& ioPoints,
			const Xfer< pointField >& points,
			const IOobject& ioFaces,
			const Xfer< faceList >& faces,
			const IOobject& ioZones,
			const Xfer< surfZoneList >& zones
		);

	// Destructor


	// Member Functions

		// Access

			//- Non-const access to the points
			pointIOField& storedIOPoints()
			{
				return points_;
			}

			//- Non-const access to the faces
			faceIOList& storedIOFaces()
			{
				return faces_;
			}

			//- Non-const access to the zones
			surfZoneIOList& storedIOZones()
			{
				return zones_;
			}

			//- Const access to the points
			const pointIOField& storedIOPoints() const
			{
				return points_;
			}

			//- Const access to the faces
			const faceIOList& storedIOFaces() const
			{
				return faces_;
			}

			//- Const access to the zones
			const surfZoneIOList& storedIOZones() const
			{
				return zones_;
			}


		//  Storage management

			//- Clear primitive data (points, faces and zones)
			void clear();

			//- Reset primitive data (points, faces and zones)
			//  Note, optimized to avoid overwriting data (with Xfer::null)
			void resetFaces
			(
				const Xfer< faceList >& faces,
				const Xfer< surfZoneList >& zones
			);

			//- Reset primitive data (points, faces and zones)
			//  Note, optimized to avoid overwriting data (with Xfer::null)
			void reset
			(
				const Xfer< pointField >& points,
				const Xfer< faceList >& faces,
				const Xfer< surfZoneList >& zones
			);

			//- Reset primitive data (points, faces and zones)
			//  Note, optimized to avoid overwriting data (with Xfer::null)
			void reset
			(
				const Xfer< List<point> >& points,
				const Xfer< faceList >& faces,
				const Xfer< surfZoneList >& zones
			);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
