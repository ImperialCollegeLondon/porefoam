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
	Foam::localBlended

Description
	Two-scheme localBlended differencing scheme.

SourceFiles
	localBlended.C

\*---------------------------------------------------------------------------*/

#ifndef localBlended_H
#define localBlended_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class localBlended
:
	public surfaceInterpolationScheme<Type>
{
	// Private Member Functions

		//- Scheme 1
		tmp<surfaceInterpolationScheme<Type> > tScheme1_;

		//- Scheme 2
		tmp<surfaceInterpolationScheme<Type> > tScheme2_;


		//- Disallow default bitwise copy construct
		localBlended(const localBlended&);

		//- Disallow default bitwise assignment
		void operator=(const localBlended&);


public:

	//- Runtime type information
	TypeName("localBlended");


	// Constructors

		//- Construct from mesh and Istream
		//  The name of the flux field is read from the Istream and looked-up
		//  from the mesh objectRegistry
		localBlended
		(
			const fvMesh& mesh,
			Istream& is
		)
		:
			surfaceInterpolationScheme<Type>(mesh),
			tScheme1_
			(
				surfaceInterpolationScheme<Type>::New(mesh, is)
			),
			tScheme2_
			(
				surfaceInterpolationScheme<Type>::New(mesh, is)
			)
		{}

		//- Construct from mesh, faceFlux and Istream
		localBlended
		(
			const fvMesh& mesh,
			const surfaceScalarField& faceFlux,
			Istream& is
		)
		:
			surfaceInterpolationScheme<Type>(mesh),
			tScheme1_
			(
				surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
			),
			tScheme2_
			(
				surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
			)
		{}


	// Member Functions

		//- Return the interpolation weighting factors
		tmp<surfaceScalarField> weights
		(
			const GeometricField<Type, fvPatchField, volMesh>& vf
		) const
		{
			const surfaceScalarField& blendingFactor =
				this->mesh().objectRegistry::template
				lookupObject<const surfaceScalarField>
				(
					word(vf.name() + "BlendingFactor")
				);

			return
				blendingFactor*tScheme1_().weights(vf)
			  + (scalar(1) - blendingFactor)*tScheme2_().weights(vf);
		}

		//- Return the face-interpolate of the given cell field
		//  with explicit correction
		tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		interpolate(const GeometricField<Type, fvPatchField, volMesh>& vf) const
		{
			const surfaceScalarField& blendingFactor =
			(
				this->mesh().objectRegistry::template
				lookupObject<const surfaceScalarField>
				(
					word(vf.name() + "BlendingFactor")
				)
			);

			return
				blendingFactor*tScheme1_().interpolate(vf)
			  + (scalar(1) - blendingFactor)*tScheme2_().interpolate(vf);
		}


		//- Return true if this scheme uses an explicit correction
		virtual bool corrected() const
		{
			return tScheme1_().corrected() || tScheme2_().corrected();
		}


		//- Return the explicit correction to the face-interpolate
		//  for the given field
		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		correction
		(
			const GeometricField<Type, fvPatchField, volMesh>& vf
		) const
		{
			const surfaceScalarField& blendingFactor =
				this->mesh().objectRegistry::template
				lookupObject<const surfaceScalarField>
				(
					word(vf.name() + "BlendingFactor")
				);

			if (tScheme1_().corrected())
			{
				if (tScheme2_().corrected())
				{
					return
					(
					    blendingFactor
					  * tScheme1_().correction(vf)
					  + (scalar(1.0) - blendingFactor)
					  * tScheme2_().correction(vf)
					);
				}
				else
				{
					return
					(
					    blendingFactor
					  * tScheme1_().correction(vf)
					);
				}
			}
			else if (tScheme2_().corrected())
			{
				return
				(
					(scalar(1.0) - blendingFactor)
				  * tScheme2_().correction(vf)
				);
			}
			else
			{
				return tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
				(
					nullptr
				);
			}
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
