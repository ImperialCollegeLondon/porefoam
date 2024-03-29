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
	Foam::fv::gaussGrad

Description
	Basic second-order gradient scheme using face-interpolation
	and Gauss' theorem.

SourceFiles
	gaussGrad.C

\*---------------------------------------------------------------------------*/

#ifndef gaussGrad_H
#define gaussGrad_H

#include "gradScheme.H"
#include "surfaceInterpolationScheme.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


template<class Type>
class gaussGrad
:
	public fv::gradScheme<Type>
{
	// Private data

		tmp<surfaceInterpolationScheme<Type> > tinterpScheme_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		gaussGrad(const gaussGrad&);

		//- Disallow default bitwise assignment
		void operator=(const gaussGrad&);


public:

	//- Runtime type information
	TypeName("Gauss");


	// Constructors

		//- Construct from mesh
		gaussGrad(const fvMesh& mesh)
		:
			gradScheme<Type>(mesh),
			tinterpScheme_(new linear<Type>(mesh))
		{}

		//- Construct from mesh and Istream
		gaussGrad(const fvMesh& mesh, Istream& is)
		:
			gradScheme<Type>(mesh),
			tinterpScheme_(nullptr)
		{
			if (is.eof())
			{
				tinterpScheme_ =
					tmp<surfaceInterpolationScheme<Type> >
					(
					    new linear<Type>(mesh)
					);
			}
			else
			{
				tinterpScheme_ =
					tmp<surfaceInterpolationScheme<Type> >
					(
					    surfaceInterpolationScheme<Type>::New(mesh, is)
					);
			}
		}


	// Member Functions

		//- Return the gradient of the given field
		//  calculated using Gauss' theorem on the given surface field
		static
		tmp
		<
			GeometricField
			<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
		> gradf
		(
			const GeometricField<Type, fvsPatchField, surfaceMesh>&,
			const word& name
		);

		//- Return the gradient of the given field calculated
		//  using Gauss' theorem on the interpolated field
		virtual tmp
		<
			GeometricField
			<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
		> calcGrad
		(
			const GeometricField<Type, fvPatchField, volMesh>& vsf,
			const word& name
		) const;

		//- Return the BlockLduSystem corresponding to the implicit grad
		// discretization. For block coupled systems.
		tmp
		<
			BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
		> fvmGrad
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) const;

		// Correct boundary conditions moved to base class
		// HJ, 14/Jun/2013
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalarGaussGrad.H"

#ifdef NoRepository
#	include "gaussGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
