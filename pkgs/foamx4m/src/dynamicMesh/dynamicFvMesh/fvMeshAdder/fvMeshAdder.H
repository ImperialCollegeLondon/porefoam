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
	Foam::fvMeshAdder

Description
	Adds two fvMeshes without using any polyMesh morphing.
	Uses fvMeshAdder.

SourceFiles
	fvMeshAdder.C
	fvMeshAdderTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshAdder_H
#define fvMeshAdder_H

#include "polyMeshAdder.H"
#include "fvPatchFieldsFwd.H"
#include "fvsPatchFieldsFwd.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class IOobject;
class faceCoupleInfo;
class IOobjectList;
class fvMesh;
class volMesh;
class surfaceMesh;
class mapAddedPolyMesh;


class fvMeshAdder
:
	public polyMeshAdder
{
private:

	// Private class

	class directFvPatchFieldMapper
	:
		public fvPatchFieldMapper
	{
		// Private data

			//- Size before mapping
			label sizeBeforeMapping_;

			//- Addressing
			const labelList& directAddressing_;

	public:

		// Constructors

			//- Construct from components
			directFvPatchFieldMapper
			(
				const label sizeBeforeMapping,
				const labelList& directAddressing
			)
			:
				fvPatchFieldMapper(),
				sizeBeforeMapping_(sizeBeforeMapping),
				directAddressing_(directAddressing)
			{}


		// Destructor

			virtual ~directFvPatchFieldMapper()
			{}


		// Member Functions

			label size() const
			{
				return directAddressing_.size();
			}

			label sizeBeforeMapping() const
			{
				return sizeBeforeMapping_;
			}

			bool direct() const
			{
				return true;
			}

			const unallocLabelList& directAddressing() const
			{
				return directAddressing_;
			}
	};


	// Private Member Functions

		//- Calculate map from new patch faces to old patch faces. -1 where
		//  could not map.
		static labelList calcPatchMap
		(
			const label oldStart,
			const label oldSize,
			const labelList& oldToNew,
			const polyPatch& newPatch,
			const label unmappedIndex
		);

		//- Map from old to new according to map. Handles map = -1.
		template<class Type>
		static void map
		(
			const Field<Type>&,
			const labelList& map,
			Field<Type>&
		);

		//- Update single volField.
		template<class Type>
		static void MapVolField
		(
			const mapAddedPolyMesh& meshMap,

			GeometricField<Type, fvPatchField, volMesh>& fld,
			const GeometricField<Type, fvPatchField, volMesh>& fldToAdd
		);

		//- Update single surfaceField.
		template<class Type>
		static void MapSurfaceField
		(
			const mapAddedPolyMesh& meshMap,

			GeometricField<Type, fvsPatchField, surfaceMesh>& fld,
			const GeometricField<Type, fvsPatchField, surfaceMesh>& fldToAdd
		);

public:

	// Member Functions

		//- Inplace add mesh to fvMesh. Maps all stored fields. Returns map.
		static autoPtr<mapAddedPolyMesh> add
		(
			fvMesh& mesh0,
			const fvMesh& mesh1,
			const faceCoupleInfo& coupleInfo,
			const bool validBoundary = true
		);

		//- Map all volFields of Type
		template<class Type>
		static void MapVolFields
		(
			const mapAddedPolyMesh&,
			const fvMesh& mesh,
			const fvMesh& meshToAdd
		);

		//- Map all surfaceFields of Type
		template<class Type>
		static void MapSurfaceFields
		(
			const mapAddedPolyMesh&,
			const fvMesh& mesh,
			const fvMesh& meshToAdd
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "fvMeshAdderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
