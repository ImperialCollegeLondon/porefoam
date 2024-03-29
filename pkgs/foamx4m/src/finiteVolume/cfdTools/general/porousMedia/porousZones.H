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
	Foam::porousZones

Description
	A centralized porousZone collection.

	Container class for a set of porousZones with the porousZone member
	functions implemented to loop over the functions for each porousZone.

	The input file @c constant/porousZone is implemented as
	IOPtrList\<porousZone\> and contains the following type of data:

	@verbatim
	1
	(
	cat1
	{
		coordinateSystem    system_10;
		porosity    0.781;
		Darcy
		{
			d   d [0 -2 0 0 0]  (-1000 -1000 0.50753e+08);
			f   f [0 -1 0 0 0]  (-1000 -1000 12.83);
		}
		heatTransfer
		{
			Maux = Mass flow of aux fluid                         // [kg/s]
			Qaux = Heat to be eliminated                          // [W]
			Caux = Specific heat capacity of aux fluid            // [J / kg K]
			Taux = initial guess for aux fluid inlet temperature  // [K]
		}
	}
	)
	@endverbatim

SourceFiles
	porousZones.C

\*---------------------------------------------------------------------------*/

#ifndef porousZones_H
#define porousZones_H

#include "porousZone.H"
#include "IOPtrList.H"

#include "volFieldsFwd.H"
#include "fvMatrices.H"
#include "geometricOneField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class porousZones
:
	public IOPtrList<porousZone>
{
	// Private data

		//- Reference to the finite volume mesh this zone is part of
		const fvMesh& mesh_;

	// Private Member Functions

		//- Disallow default bitwise copy construct
		porousZones(const porousZones&);

		//- Disallow default bitwise assignment
		void operator=(const porousZones&);


		//- Modify time derivative elements
		template<class Type>
		void modifyDdt(fvMatrix<Type>&) const;

public:

	// Constructors

		//- Construct from fvMesh
		//  with automatically constructed coordinate systems list
		porousZones(const fvMesh&);


	// Member Functions

		//- mirror fvm::ddt with porosity
		template<class Type>
		tmp<fvMatrix<Type> > ddt
		(
			GeometricField<Type, fvPatchField, volMesh>&
		);

		//- mirror fvm::ddt with porosity
		template<class Type>
		tmp<fvMatrix<Type> > ddt
		(
			const geometricOneField&,
			GeometricField<Type, fvPatchField, volMesh>&
		);

		//- mirror fvm::ddt with porosity
		template<class Type>
		tmp<fvMatrix<Type> > ddt
		(
			const dimensionedScalar&,
			GeometricField<Type, fvPatchField, volMesh>&
		);

		//- mirror fvm::ddt with porosity
		template<class Type>
		tmp<fvMatrix<Type> > ddt
		(
			const volScalarField&,
			GeometricField<Type, fvPatchField, volMesh>&
		);

		//- Add the viscous and inertial resistance force contribution
		//  to the momentum equation
		void addResistance(fvVectorMatrix& UEqn) const;

		//- Add the viscous and inertial resistance force contribution
		//  to the tensorial diagonal
		void addResistance
		(
			const fvVectorMatrix& UEqn,
			volTensorField& AU
		) const;

		//- Add the heat source contribution to enalphy equation (energy
		// equation)
		void addHeatResistance
		(
			fvScalarMatrix& hTEqn,
			const volScalarField& T,
			volScalarField& Taux,
			volScalarField& Qaux,
			const volVectorField& U,
			const volScalarField& Macro,
			const volScalarField& posFlux
		) const;

		//- Order cells for the Dual Stream model
		void macroCellOrder
		(
			volScalarField& Taux,
			volScalarField& Macro,
			volScalarField& posFlux,
			const surfaceScalarField& phi
		) const;

		//- Read modified data
		virtual bool readData(Istream&);

		//- Write data
		bool writeData(Ostream&, bool subDict = true) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "porousZonesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
