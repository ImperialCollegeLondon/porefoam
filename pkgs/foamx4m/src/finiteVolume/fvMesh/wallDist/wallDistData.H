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
	Foam::wallDistData

Description
	Wall distance calculation. Like wallDist but also transports extra
	data (template argument).

	Used for e.g reflection vector calculation or vanDriest damping.

	Templated on two parameters:
		- TransferType: type of overall data transported
		  (e.g. wallPointData\<vector\>)

SourceFiles
	wallDistData.C

\*---------------------------------------------------------------------------*/

#ifndef wallDistData_H
#define wallDistData_H

#include "cellDistFuncs.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class TransferType>
class wallDistData
:
	public volScalarField,
	public cellDistFuncs
{

private:

	typedef typename TransferType::dataType Type;


	// Private Member Data

		//- Reference to field whose data to use (on walls) and update
		//  (every cell and non-wall face)
		GeometricField<Type, fvPatchField, volMesh>& field_;

		//- Do accurate distance calculation for near-wall cells.
		bool correctWalls_;

		//- Number of unset cells and faces.
		label nUnset_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		wallDistData(const wallDistData&);

		//- Disallow default bitwise assignment
		void operator=(const wallDistData&);


public:

	// Constructors

		//- Construct from mesh and flag whether or not to correct wall.
		//  Calculate for all cells. correctWalls : correct wall (face&point)
		//  cells for correct distance, searching neighbours.
		wallDistData
		(
			const fvMesh& mesh,
			GeometricField<Type, fvPatchField, volMesh>&,
			bool correctWalls = true
		);


	// Destructor

		virtual ~wallDistData();


	// Member Functions

		const volScalarField& y() const
		{
			return *this;
		}

		label nUnset() const
		{
			return nUnset_;
		}

		//- Access field
		const GeometricField<Type, fvPatchField, volMesh>& data() const
		{
			return field_;
		}

		//- Correct for mesh geom/topo changes
		virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "wallDistData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
