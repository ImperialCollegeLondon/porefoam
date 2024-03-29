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
	Foam::fv::fourthGrad

Description
	Second-order gradient scheme using least-squares.

SourceFiles
	fourthGrad.C

\*---------------------------------------------------------------------------*/

#ifndef fourthGrad_H
#define fourthGrad_H

#include "gradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


template<class Type>
class fourthGrad
:
	public fv::gradScheme<Type>
{
	// Private Member Functions

		//- Disallow default bitwise copy construct
		fourthGrad(const fourthGrad&);

		//- Disallow default bitwise assignment
		void operator=(const fourthGrad&);


public:

	//- Runtime type information
	TypeName("fourth");


	// Constructors

		//- Construct from Istream
		fourthGrad(const fvMesh& mesh, Istream&)
		:
			gradScheme<Type>(mesh)
		{}


	// Member Functions

		tmp
		<
			GeometricField
			<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
		> calcGrad
		(
			const GeometricField<Type, fvPatchField, volMesh>& vsf,
			const word& name
		) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "fourthGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
