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
	Foam::MeshWave

Description
	FaceCellWave plus data

SourceFiles
	MeshWave.C

\*---------------------------------------------------------------------------*/

#ifndef MeshWave_H
#define MeshWave_H

#include "FaceCellWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


TemplateName(MeshWave);



template <class Type>
class MeshWave
:
	public MeshWaveName
{
	// Private data

		//- Wall information for all faces
		List<Type> allFaceInfo_;

		//- Wall information for all cells
		List<Type> allCellInfo_;

		//- Wave calculation engine.
		FaceCellWave<Type> calc_;

	// Private Member Functions

		//- Disallow default bitwise copy construct
		MeshWave(const MeshWave&);

		//- Disallow default bitwise assignment
		void operator=(const MeshWave&);

public:

	// Constructors

		//- Construct from mesh and list of changed faces with the Type
		//  for these faces. Iterates until nothing changes or maxIter reached.
		//  (maxIter can be 0)
		MeshWave
		(
			const polyMesh& mesh,
			const labelList& initialChangedFaces,
			const List<Type>& changedFacesInfo,
			const label maxIter
		);

		//- Construct from mesh, list of changed faces with the Type
		//  for these faces and initial field.
		//  Iterates until nothing changes or maxIter reached.
		//  (maxIter can be 0)
		MeshWave
		(
			const polyMesh& mesh,
			const labelList& initialChangedFaces,
			const List<Type>& changedFacesInfo,
			const List<Type>& allCellInfo,
			const label maxIter
		);


	// Member Functions

		//- Get allFaceInfo
		const List<Type>& allFaceInfo() const
		{
			return allFaceInfo_;
		}

		//- Get allCellInfo
		const List<Type>& allCellInfo() const
		{
			return allCellInfo_;
		}


		//- Iterate until no changes or maxIter reached. Returns number of
		//  unset cells (see getUnsetCells)
		label iterate(const label maxIter)
		{
			return calc_.iterate(maxIter);
		}

		//- Get number of unvisited cells, i.e. cells that were not (yet)
		//  reached from walking across mesh. This can happen from
		//  - not enough iterations done
		//  - a disconnected mesh
		//  - a mesh without walls in it
		label getUnsetCells() const
		{
			return calc_.getUnsetCells();
		}

		//- Get number of unvisited faces
		label getUnsetFaces() const
		{
			return calc_.getUnsetFaces();
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "MeshWave.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
