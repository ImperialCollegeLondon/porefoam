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
	Foam::localMin

Description
	LocalMin-mean differencing scheme class.

	This scheme interpolates 1/field using a scheme specified at run-time
	and return the reciprocal of the interpolate.

SourceFiles
	localMin.C

\*---------------------------------------------------------------------------*/

#ifndef localMin_H
#define localMin_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class localMin
:
	public surfaceInterpolationScheme<Type>
{
	// Private Member Functions

		//- Disallow default bitwise assignment
		void operator=(const localMin&);


public:

	//- Runtime type information
	TypeName("localMin");


	// Constructors

		//- Construct from mesh
		localMin(const fvMesh& mesh)
		:
			surfaceInterpolationScheme<Type>(mesh)
		{}

		//- Construct from Istream.
		//  The name of the flux field is read from the Istream and looked-up
		//  from the mesh objectRegistry
		localMin
		(
			const fvMesh& mesh,
			Istream& is
		)
		:
			surfaceInterpolationScheme<Type>(mesh)
		{}

		//- Construct from faceFlux and Istream
		localMin
		(
			const fvMesh& mesh,
			const surfaceScalarField& faceFlux,
			Istream& is
		)
		:
			surfaceInterpolationScheme<Type>(mesh)
		{}


	// Member Functions

		//- Return the interpolation weighting factors
		virtual tmp<surfaceScalarField> weights
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) const
		{
			notImplemented
			(
				"localMin::weights"
				"(const GeometricField<Type, fvPatchField, volMesh>&)"
			);

			return tmp<surfaceScalarField>(nullptr);
		}

		//- Return the face-interpolate of the given cell field
		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		interpolate
		(
			const GeometricField<Type, fvPatchField, volMesh>& vf
		) const
		{
			const fvMesh& mesh = vf.mesh();

			tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tvff
			(
				new GeometricField<Type, fvsPatchField, surfaceMesh>
				(
					IOobject
					(
					    vf.name(),
					    mesh.time().timeName(),
					    mesh
					),
					mesh,
					vf.dimensions()
				)
			);
			GeometricField<Type, fvsPatchField, surfaceMesh>& vff = tvff();

			forAll(vff.boundaryField(), patchi)
			{
				vff.boundaryField()[patchi] = vf.boundaryField()[patchi];
			}

			const unallocLabelList& own = mesh.owner();
			const unallocLabelList& nei = mesh.neighbour();

			forAll(vff, facei)
			{
				vff[facei] = minMod(vf[own[facei]], vf[nei[facei]]);
			}

			return tvff;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
