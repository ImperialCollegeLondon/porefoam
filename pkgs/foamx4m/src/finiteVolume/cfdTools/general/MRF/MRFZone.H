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
	Foam::MRFZone

Description
	MRF zone definition based on cell zone and parameters
	obtained from a control dictionary constructed from the given stream.

	The rotation of the MRF region is defined by an origin and axis of
	rotation and an angular speed in rad/s.

SourceFiles
	MRFZone.C

\*---------------------------------------------------------------------------*/

#ifndef MRFZone_H
#define MRFZone_H

#include "dictionary.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;


class MRFZone
{
	// Private data

		//- Name of the MRF zone
		const word name_;

		//- Mesh reference
		const fvMesh& mesh_;

		//- Reference to dictionary
		const dictionary dict_;

		//- MRF cell zone ID: cells in MRF region
		label cellZoneID_;

		//- List of patch names touching the MRF zone which do not move
		const wordList excludedPatchNames_;

		//- List of patch labels touching the MRF zone which do not move
		labelList excludedPatchLabels_;

		//- Internal faces that are part of MRF
		labelList internalFaces_;

		//- Outside faces (per patch) that move with the MRF
		labelListList includedFaces_;

		//- Excluded faces (per patch) that do not move with the MRF
		labelListList excludedFaces_;

		//- Origin of rotation
		const dimensionedVector origin_;

		//- Axis of rotation
		dimensionedVector axis_;

		//- Rotational velocity in rad/s
		dimensionedScalar omega_;

		//- Ramping time scale
		scalar rampTime_;

		//- Mesh velocity
		//  Note: for precise calculation of the mesh velocity on bad meshes
		//  the velocity is calculated from the swept volume
		mutable surfaceScalarField* meshVelocityPtr_;

		//- Time for which mesh velocity is calculated
		mutable scalar meshVelTime_;


	// Private Static Data

		//- Rotational angle in degrees used to determine the time-step size
		static const debug::tolerancesSwitch rotAngle_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		MRFZone(const MRFZone&);

		//- Disallow default bitwise assignment
		void operator=(const MRFZone&);


		//- Divide faces in frame according to patch
		void setMRFFaces();

		//- Calculate mesh velocity
		void calcMeshVelocity() const;


		//- Make the given absolute mass/vol flux relative within MRF region
		template<class RhoFieldType>
		void relativeRhoFlux
		(
			const RhoFieldType& rho,
			surfaceScalarField& phi
		) const;

		//- Make the given relative mass/vol flux absolute within MRF region
		template<class RhoFieldType>
		void absoluteRhoFlux
		(
			const RhoFieldType& rho,
			surfaceScalarField& phi
		) const;


public:

	// Declare name of the class and its debug switch
	ClassName("MRFZone");


	// Constructors

		//- Construct from fvMesh and Istream
		MRFZone(const fvMesh& mesh, Istream& is);

		//- Return clone
		autoPtr<MRFZone> clone() const
		{
			notImplemented("autoPtr<MRFZone> clone() const");
			return autoPtr<MRFZone>(nullptr);
		}

		//- Return a pointer to a new MRFZone created on freestore
		//  from Istream
		class iNew
		{
			const fvMesh& mesh_;

		public:

			iNew(const fvMesh& mesh)
			:
				mesh_(mesh)
			{}

			autoPtr<MRFZone> operator()(Istream& is) const
			{
				return autoPtr<MRFZone>(new MRFZone(mesh_, is));
			}
		};

	//- Destructor
	~MRFZone();


	// Member Functions

		// Access Functions

			//- Return the MRFZone name
			const word& name() const
			{
				return name_;
			}

			//- Return cellZone number
			label zoneId() const
			{
				return cellZoneID_;
			}

			//- Return cellZone
			const cellZone& zone() const
			{
				return mesh_.cellZones()[cellZoneID_];
			}

			//- Return the MRFZone dictionary
			const dictionary& dict() const
			{
				return dict_;
			}

			//- Return the MRFZone axis of rotation
			const dimensionedVector& axis() const
			{
				return axis_;
			}

			//- Return rotational vector
			vector Omega() const;

			//- Return the MRFZone omega
			const dimensionedScalar& omega() const
			{
				return omega_;
			}

			//- Return mesh velocity for all faces
			const surfaceScalarField& meshVelocity() const;


	// Member functions used in a top level solver

		//- Update the mesh corresponding to given map
		void updateMesh(const mapPolyMesh& mpm)
		{
			// Clear mesh velocity
			deleteDemandDrivenData(meshVelocityPtr_);

			// Only updates face addressing
			setMRFFaces();
		}

		//- Add the Coriolis force contribution to the momentum equation
		void addCoriolis(fvVectorMatrix& UEqn) const;

		//- Add the Coriolis force contribution to the momentum equation
		void addCoriolis
		(
			const volScalarField& rho,
			fvVectorMatrix& UEqn
		) const;

		//- Add omega within the MRF region
		void addOmega(volVectorField& omg) const;

		//- Make the given absolute velocity relative within the MRF region
		void relativeVelocity(volVectorField& U) const;

		//- Make the given relative velocity absolute within the MRF region
		void absoluteVelocity(volVectorField& U) const;

		//- Make the given absolute flux relative within the MRF region
		void relativeFlux(surfaceScalarField& phi) const;

		//- Make the given absolute mass-flux relative within the MRF region
		void relativeFlux
		(
			const surfaceScalarField& rho,
			surfaceScalarField& phi
		) const;

		//- Make the given relative flux absolute within the MRF region
		void absoluteFlux(surfaceScalarField& phi) const;

		//- Make the given relative mass-flux absolute within the MRF region
		void absoluteFlux
		(
			const surfaceScalarField& rho,
			surfaceScalarField& phi
		) const;

		//- Compute pseudo mesh phi for the MRF region
		void meshPhi(surfaceScalarField& phi) const;

		//- Correct the boundary velocity for the roation of the MRF region
		void correctBoundaryVelocity(volVectorField& U) const;

		//- Compute source term for volScalarFields in each cellZone
		void Su
		(
			const volScalarField& phi,
			const volVectorField& gradPhi,
			volScalarField& source
		) const;

		//- Compute source term for volVectorFields in each cellZone
		void Su
		(
			const volVectorField& phi,
			const volTensorField& gradPhi,
			volVectorField& source
		) const;

	// Ostream Operator

		friend Ostream& operator<<(Ostream& os, const MRFZone&)
		{
			notImplemented("Ostream& operator<<(Ostream& os, const MRFZone&)");
			return os;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "MRFZoneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
