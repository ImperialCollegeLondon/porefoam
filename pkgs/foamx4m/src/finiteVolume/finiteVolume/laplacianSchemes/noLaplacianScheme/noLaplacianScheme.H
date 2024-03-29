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
	Foam::fv::noLaplacianScheme

Description
	No laplacian scheme.  Used to switch off diffusion without recompilation.

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
	noLaplacianScheme.C

\*---------------------------------------------------------------------------*/

#ifndef noLaplacianScheme_H
#define noLaplacianScheme_H

#include "laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


template<class Type, class GType>
class noLaplacianScheme
:
	public fv::laplacianScheme<Type, GType>
{
	// Private Member Functions

		//- Disallow default bitwise copy construct
		noLaplacianScheme(const noLaplacianScheme&);

		//- Disallow default bitwise assignment
		void operator=(const noLaplacianScheme&);


public:

	//- Runtime type information
	TypeName("off");


	// Constructors

		//- Construct null
		noLaplacianScheme(const fvMesh& mesh)
		:
			laplacianScheme<Type, GType>(mesh)
		{}

		//- Construct from Istream
		noLaplacianScheme(const fvMesh& mesh, Istream& is)
		:
			laplacianScheme<Type, GType>(mesh, is)
		{}

		//- Construct from mesh, interpolation and snGradScheme schemes
		noLaplacianScheme
		(
			const fvMesh& mesh,
			const tmp<surfaceInterpolationScheme<GType> >& igs,
			const tmp<snGradScheme<Type> >& sngs
		)
		:
			laplacianScheme<Type, GType>(mesh, igs, sngs)
		{}


	// Destructor

		virtual ~noLaplacianScheme()
		{}


	// Member Functions

		tmp<GeometricField<Type, fvPatchField, volMesh> > fvcLaplacian
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		);

		tmp<fvMatrix<Type> > fvmLaplacian
		(
			const GeometricField<GType, fvsPatchField, surfaceMesh>&,
			const GeometricField<Type, fvPatchField, volMesh>&
		);

		tmp<GeometricField<Type, fvPatchField, volMesh> > fvcLaplacian
		(
			const GeometricField<GType, fvsPatchField, surfaceMesh>&,
			const GeometricField<Type, fvPatchField, volMesh>&
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "noLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
