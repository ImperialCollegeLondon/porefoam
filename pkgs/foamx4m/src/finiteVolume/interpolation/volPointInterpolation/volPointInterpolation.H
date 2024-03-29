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
	Foam::volPointInterpolation

Description
	Foam::volPointInterpolation

SourceFiles
	volPointInterpolation.C
	volPointInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef volPointInterpolation_H
#define volPointInterpolation_H

#include "MeshObject.H"
#include "pointPatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class pointMesh;


class volPointInterpolation
:
	public MeshObject<fvMesh, volPointInterpolation>
{
	// Private data

		//- Boundary interpolation engine
		pointPatchInterpolation boundaryInterpolator_;

		//- Interpolation scheme weighting factor array using lazy evaluation
		//  Updated for MeshObject.  HJ, 30/Aug/2010
		mutable scalarListList* pointWeightsPtr_;


	// Private member functions

		//- Construct point weighting factors
		void makeWeights() const;

		//- Disallow default bitwise copy construct
		volPointInterpolation(const volPointInterpolation&);

		//- Disallow default bitwise assignment
		void operator=(const volPointInterpolation&);


public:

	// Declare name of the class and its debug switch
	TypeName("volPointInterpolation");


	// Constructors

		//- Constructor given fvMesh.  pointMesh will be created or
		//  looked up from objectRegistry
		explicit volPointInterpolation(const fvMesh&);


	// Destructor

		virtual ~volPointInterpolation();


	// Member functions

		// Access

			//- Const reference to fvMesh
			const fvMesh& mesh() const
			{
				return boundaryInterpolator_.mesh();
			}

			//- Return interpolation weights (lazy evaluation)
			const scalarListList& pointWeights() const;


		// Edit

			//- Correct weighting factors for moving mesh.
			//  Updated for MeshObject.  HJ, 30/Aug/2010
			virtual bool movePoints() const;

			//- Update mesh topology using the morph engine
			//  Updated for MeshObject.  HJ, 30/Aug/2010
			virtual bool updateMesh(const mapPolyMesh&) const;


	// Interpolation functions

		//- Interpolate internal field from volField to pointField
		//  using inverse distance weighting
		template<class Type>
		void interpolateInternalField
		(
			const GeometricField<Type, fvPatchField, volMesh>&,
			GeometricField<Type, pointPatchField, pointMesh>&
		) const;

		//- Interpolate from volField to pointField
		//  using inverse distance weighting
		template<class Type>
		void interpolate
		(
			const GeometricField<Type, fvPatchField, volMesh>&,
			GeometricField<Type, pointPatchField, pointMesh>&
		) const;

		//- Interpolate volField using inverse distance weighting
		//  returning pointField with the same patchField types
		template<class Type>
		tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
		(
			const GeometricField<Type, fvPatchField, volMesh>&,
			const wordList& patchFieldTypes
		) const;

		//- Interpolate tmp<volField> using inverse distance weighting
		//  returning pointField with the same patchField types
		template<class Type>
		tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
		(
			const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
			const wordList& patchFieldTypes
		) const;

		//- Interpolate volField using inverse distance weighting
		//  returning pointField
		template<class Type>
		tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) const;

		//- Interpolate tmp<volField> using inverse distance weighting
		//  returning pointField
		template<class Type>
		tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
		(
			const tmp<GeometricField<Type, fvPatchField, volMesh> >&
		) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "volPointInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
