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
	Foam::pointPatchInterpolation

Description
	Foam::pointPatchInterpolation

SourceFiles
	pointPatchInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef pointPatchInterpolation_H
#define pointPatchInterpolation_H

#include "primitivePatchInterpolation.H"
#include "PtrList.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "scalarList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class pointMesh;


class pointPatchInterpolation
{
	// Private data

		//- Reference to volume mesh
		const fvMesh& fvMesh_;

		//- Primitive patch interpolators
		mutable PtrList<primitivePatchInterpolation> patchInterpolators_;

		//- List of patch-patch edge points that require special treatement
		mutable labelList* patchPatchPointsPtr_;

		//- List of constraint points
		mutable labelList* patchPatchPointConstraintPointsPtr_;

		//- List of constraint tensors for constrained points
		mutable tensorField* patchPatchPointConstraintTensorsPtr_;

		//- Weights for patch-patch boundary points
		mutable scalarListList* patchPatchPointWeightsPtr_;


	// Private member functions

		//- Construct addressing for patch-patch boundary points
		void makePatchPatchAddressing() const;

		//- Construct weights for patch-patch boundary points
		void makePatchPatchWeights() const;

		//- Clear all demand driven data
		void clearOut() const;

		//- Disallow default bitwise copy construct
		pointPatchInterpolation(const pointPatchInterpolation&);

		//- Disallow default bitwise assignment
		void operator=(const pointPatchInterpolation&);


public:

	// Declare name of the class and its debug switch
	ClassName("pointPatchInterpolation");


	// Constructors

		//- Constructor given fvMesh and pointMesh.
		pointPatchInterpolation(const fvMesh&);


	// Destructor

		~pointPatchInterpolation();


	// Member functions

		// Access

			//- Return const reference to fvMesh
			const fvMesh& mesh() const
			{
				return fvMesh_;
			}

			//- Return primitive patch interpolators
			const PtrList<primitivePatchInterpolation>& patchInterpolators
			() const;

			//- Return patch-patch edge points that require special treatment
			const labelList& patchPatchPoints() const;

			//- Return weights for patch-patch boundary points
			const scalarListList& patchPatchPointWeights() const;

			//- Return constraint points
			const labelList& patchPatchPointConstraintPoints() const;

			//- Return constraint tensors
			const tensorField& patchPatchPointConstraintTensors() const;


		// Edit

			//- Update mesh topology using the morph engine
			void updateMesh();

			//- Correct weighting factors for moving mesh.
			bool movePoints();


	// Interpolation functions

		template<class Type>
		void interpolate
		(
			const GeometricField<Type, fvPatchField, volMesh>&,
			GeometricField<Type, pointPatchField, pointMesh>&,
			bool overrideFixedValue
		) const;

		template<class Type>
		void applyCornerConstraints
		(
			GeometricField<Type, pointPatchField, pointMesh>& pf
		) const;
};


template<>
void pointPatchInterpolation::applyCornerConstraints<scalar>
(
	GeometricField<scalar, pointPatchField, pointMesh>& pf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "pointPatchInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
