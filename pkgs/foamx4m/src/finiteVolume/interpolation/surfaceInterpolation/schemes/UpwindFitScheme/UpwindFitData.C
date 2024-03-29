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

\*---------------------------------------------------------------------------*/

#include "UpwindFitData.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "SVD.H"
#include "syncTools.H"
#include "extendedUpwindCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Polynomial>
Foam::UpwindFitData<Polynomial>::UpwindFitData
(
	const fvMesh& mesh,
	const extendedUpwindCellToFaceStencil& stencil,
	const bool linearCorrection,
	const scalar linearLimitFactor,
	const scalar centralWeight
)
:
	FitData
	<
		UpwindFitData<Polynomial>,
		extendedUpwindCellToFaceStencil,
		Polynomial
	>
	(
		mesh, stencil, linearCorrection, linearLimitFactor, centralWeight
	),
	owncoeffs_(mesh.nFaces()),
	neicoeffs_(mesh.nFaces())
{
	if (debug)
	{
		Info<< "Contructing UpwindFitData<Polynomial>" << endl;
	}

	calcFit();

	if (debug)
	{
		Info<< "UpwindFitData<Polynomial>::UpwindFitData() :"
			<< "Finished constructing polynomialFit data"
			<< endl;
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Polynomial>
void Foam::UpwindFitData<Polynomial>::calcFit()
{
	const fvMesh& mesh = this->mesh();

	const surfaceScalarField& w = mesh.surfaceInterpolation::weights();
	const surfaceScalarField::Boundary& bw = w.boundaryField();

	// Owner stencil weights
	// ~~~~~~~~~~~~~~~~~~~~~

	// Get the cell/face centres in stencil order.
	List<List<point> > stencilPoints(mesh.nFaces());
	this->stencil().collectData
	(
		this->stencil().ownMap(),
		this->stencil().ownStencil(),
		mesh.C(),
		stencilPoints
	);

	// find the fit coefficients for every owner

	//Pout<< "-- Owner --" << endl;
	for(label facei = 0; facei < mesh.nInternalFaces(); facei++)
	{
		FitData
		<
			UpwindFitData<Polynomial>,
			extendedUpwindCellToFaceStencil,
			Polynomial
		>::calcFit(owncoeffs_[facei], stencilPoints[facei], w[facei], facei);

		//Pout<< "    facei:" << facei
		//    << " at:" << mesh.faceCentres()[facei] << endl;
		//forAll(owncoeffs_[facei], i)
		//{
		//    Pout<< "    point:" << stencilPoints[facei][i]
		//        << "\tweight:" << owncoeffs_[facei][i]
		//        << endl;
		//}
	}

	forAll(bw, patchi)
	{
		const fvsPatchScalarField& pw = bw[patchi];

		if (pw.coupled())
		{
			label facei = pw.patch().patch().start();

			forAll(pw, i)
			{
				FitData
				<
					UpwindFitData<Polynomial>,
					extendedUpwindCellToFaceStencil,
					Polynomial
				>::calcFit
				(
					owncoeffs_[facei], stencilPoints[facei], pw[i], facei
				);
				facei++;
			}
		}
	}


	// Neighbour stencil weights
	// ~~~~~~~~~~~~~~~~~~~~~~~~~

	// Note:reuse stencilPoints since is major storage
	this->stencil().collectData
	(
		this->stencil().neiMap(),
		this->stencil().neiStencil(),
		mesh.C(),
		stencilPoints
	);

	// find the fit coefficients for every neighbour

	//Pout<< "-- Neighbour --" << endl;
	for(label facei = 0; facei < mesh.nInternalFaces(); facei++)
	{
		FitData
		<
			UpwindFitData<Polynomial>,
			extendedUpwindCellToFaceStencil,
			Polynomial
		>::calcFit(neicoeffs_[facei], stencilPoints[facei], w[facei], facei);

		//Pout<< "    facei:" << facei
		//    << " at:" << mesh.faceCentres()[facei] << endl;
		//forAll(neicoeffs_[facei], i)
		//{
		//    Pout<< "    point:" << stencilPoints[facei][i]
		//        << "\tweight:" << neicoeffs_[facei][i]
		//        << endl;
		//}
	}

	forAll(bw, patchi)
	{
		const fvsPatchScalarField& pw = bw[patchi];

		if (pw.coupled())
		{
			label facei = pw.patch().patch().start();

			forAll(pw, i)
			{
				FitData
				<
					UpwindFitData<Polynomial>,
					extendedUpwindCellToFaceStencil,
					Polynomial
				>::calcFit
				(
					neicoeffs_[facei], stencilPoints[facei], pw[i], facei
				);
				facei++;
			}
		}
	}
}


// ************************************************************************* //
