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

Namespace
	Foam::SRF

Description
	Namespace for single rotating frame (SRF) models

Class
	Foam::SRF::SRFModel

Description
	Top level model for single rotating frame
	- Steady state only - no time derivatives included

SourceFiles
	SRFModel.C

\*---------------------------------------------------------------------------*/

#ifndef SRFModel_H
#define SRFModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "volFields.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace SRF
{


class SRFModel
:
	public IOdictionary
{
protected:

	// Protected data

		//- Reference to the relative velocity field
		const volVectorField& Urel_;

		//- Reference to the mesh
		const fvMesh& mesh_;

		//- Axis of rotation
		vector axis_;

		//- SRF model coeficients dictionary
		dictionary SRFModelCoeffs_;

		//- Angular velocity of the frame (rad/s)
		dimensionedVector omega_;


private:

	// Private Member Functions

		//- Disallow default bitwise copy construct
		SRFModel(const SRFModel&);

		//- Disallow default bitwise assignment
		void operator=(const SRFModel&);


public:

	//- Runtime type information
	TypeName("SRFModel");


	// Declare runtime constructor selection table

		 declareRunTimeSelectionTable
		 (
			 autoPtr,
			 SRFModel,
			 dictionary,
			 (
				 const volVectorField& Urel
			 ),
			 (Urel)
		 );


	// Constructors

		//- Construct from components
		SRFModel
		(
			const word& type,
			const volVectorField& Urel
		);


	// Selectors

		 //- Return a reference to the selected SRF model
		 static autoPtr<SRFModel> New
		 (
			 const volVectorField& Urel
		 );


	// Destructor

		virtual ~SRFModel();


	// Member Functions

		// Edit

			//- Read radiationProperties dictionary
			virtual bool read();


		// Access

			//- Return the axis of rotation
			const vector& axis() const;

			//- Return the angular velocity field [rad/s]
			const dimensionedVector& omega() const;

			//- Return the coriolis force
			tmp<DimensionedField<vector, volMesh> > Fcoriolis() const;

			//- Return the centrifugal force
			tmp<DimensionedField<vector, volMesh> > Fcentrifugal() const;

			//- Source term component for momentum equation
			tmp<DimensionedField<vector, volMesh> > Su() const;

			//- Return velocity vector from positions
			vectorField velocity(const vectorField& positions) const;

			//- Return rotational velocity of SRF for complete mesh
			tmp<volVectorField> U() const;

			//- Return rotational face velocity of SRF for complete mesh
			tmp<surfaceVectorField> faceU() const;

			//- Return absolute velocity for complete mesh
			tmp<volVectorField> Uabs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SRF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

