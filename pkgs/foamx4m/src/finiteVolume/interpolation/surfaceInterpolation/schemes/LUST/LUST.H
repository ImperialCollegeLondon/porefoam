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
	Foam::LUST

Description
	LUST: Linear-upwind stabilised transport.

	Interpolation scheme class derived from linearUpwind which returns blended
	linear/linear-upwind weighting factors and also applies a explicit
	gradient-based correction obtained from the linearUpwind scheme.  The
	blending-factor is set to 0.75 linear which optimises the balance between
	accuracy and stability on a range of LES cases with a range of mesh quality.

SourceFiles
	LUST.C

\*---------------------------------------------------------------------------*/

#ifndef LUST_H
#define LUST_H

#include "linearUpwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class LUST
:
	public linearUpwind<Type>
{
	// Private Member Functions

		//- Disallow default bitwise copy construct
		LUST(const LUST&);

		//- Disallow default bitwise assignment
		void operator=(const LUST&);


public:

	//- Runtime type information
	TypeName("LUST");


	// Constructors

		//- Construct from mesh and Istream
		LUST
		(
			const fvMesh& mesh,
			Istream& schemeData
		)
		:
			linearUpwind<Type>(mesh, schemeData)
		{}

		//- Construct from mesh, faceFlux and Istream
		LUST
		(
			const fvMesh& mesh,
			const surfaceScalarField& faceFlux,
			Istream& schemeData
		)
		:
			linearUpwind<Type>(mesh, faceFlux, schemeData)
		{}


	// Member Functions

		//- Return the interpolation weighting factors
		virtual tmp<surfaceScalarField> weights
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) const
		{
			return
				0.75*this->mesh().surfaceInterpolation::weights()
			  + 0.25*linearUpwind<Type>::weights();
		}

		//- Return true if this scheme uses an explicit correction
		virtual bool corrected() const
		{
			return true;
		}

		//- Return the explicit correction to the face-interpolate
		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
		correction
		(
			const GeometricField<Type, fvPatchField, volMesh>& vf
		) const
		{
			return 0.25*linearUpwind<Type>::correction(vf);
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
