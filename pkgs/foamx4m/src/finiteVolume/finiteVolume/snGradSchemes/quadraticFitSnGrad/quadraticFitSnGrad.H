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
	quadraticFitSnGrad

Description
	Simple central-difference snGrad scheme with quadratic fit correction from
	a larger stencil.

SourceFiles
	quadraticFitSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticFitSnGrad_H
#define quadraticFitSnGrad_H

#include "snGradScheme.H"
#include "quadraticFitSnGradData.H"
#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


template<class Type>
class quadraticFitSnGrad
:
	public snGradScheme<Type>
{
	// Private Data
		//- weights for central stencil
		const scalar centralWeight_;

	// Private Member Functions

		//- Disallow default bitwise assignment
		void operator=(const quadraticFitSnGrad&);


public:

	//- Runtime type information
	TypeName("quadraticFit");


	// Constructors

		//- Construct from mesh and scheme data
		quadraticFitSnGrad
		(
			const fvMesh& mesh,
			const scalar centralWeight
		)
		:
			snGradScheme<Type>(mesh),
			centralWeight_(centralWeight)
		{}


		//- Construct from mesh and data stream
		quadraticFitSnGrad(const fvMesh& mesh, Istream& is)
		:
			snGradScheme<Type>(mesh),
			centralWeight_(readScalar(is))
		{}


	// Destructor

		virtual ~quadraticFitSnGrad() {}


	// Member Functions

		//- Return the interpolation weighting factors for the given field
		virtual tmp<surfaceScalarField> deltaCoeffs
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) const
		{
			return this->mesh().deltaCoeffs();
		}

		//- Return true if this scheme uses an explicit correction
		virtual bool corrected() const
		{
			return true;
		}

		//- Return the explicit correction to the quadraticFitSnGrad
		//  for the given field
		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		correction(const GeometricField<Type, fvPatchField, volMesh>& vf) const
		{
			const fvMesh& mesh = this->mesh();

			const quadraticFitSnGradData& qfd = quadraticFitSnGradData::New
			(
				mesh,
				centralWeight_
			);

			const extendedCellToFaceStencil& stencil = qfd.stencil();
			const List<scalarList>& f = qfd.fit();

			tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > sft
				= stencil.weightedSum(vf, f);

			sft().dimensions() /= dimLength;

			return sft;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
