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

#include "MULES.H"
#include "upwind.H"
#include "uncorrectedSnGrad.H"
#include "gaussConvectionScheme.H"
#include "gaussLaplacianScheme.H"
#include "uncorrectedSnGrad.H"
#include "surfaceInterpolate.H"
#include "fvcSurfaceIntegrate.H"
#include "slicedSurfaceFields.H"
#include "syncTools.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RhoType, class SpType, class SuType>
void Foam::MULES::explicitSolve
(
	const RhoType& rho,
	volScalarField& psi,
	const surfaceScalarField& phi,
	surfaceScalarField& phiPsi,
	const SpType& Sp,
	const SuType& Su,
	const scalar psiMax,
	const scalar psiMin
)
{
	Info<< "MULES: Solving for " << psi.name() << endl;

	const fvMesh& mesh = psi.mesh();
	psi.correctBoundaryConditions();

	surfaceScalarField phiBD = upwind<scalar>(psi.mesh(), phi).flux(psi);

	surfaceScalarField& phiCorr = phiPsi;
	phiCorr -= phiBD;

	scalarField allLambda(mesh.nFaces(), 1.0);

	slicedSurfaceScalarField lambda
	(
		IOobject
		(
			"lambda",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE,
			false
		),
		mesh,
		dimless,
		allLambda,
		false   // Use slices for the couples
	);

	limiter
	(
		allLambda,
		rho,
		psi,
		phiBD,
		phiCorr,
		Sp.field(),
		Su.field(),
		psiMax,
		psiMin,
		3
	);

	phiPsi = phiBD + lambda*phiCorr;

	scalarField& psiIf = psi;
	const scalarField& psi0 = psi.oldTime();
	const scalar deltaT = mesh.time().deltaT().value();

	psiIf = 0.0;
	fvc::surfaceIntegrate(psiIf, phiPsi);

	if (mesh.moving())
	{
		psiIf =
		(
			mesh.Vsc0()*rho.oldTime()*psi0/(deltaT*mesh.Vsc())
		  + Su.field()
		  - psiIf
		)/(rho/deltaT - Sp.field());
	}
	else
	{
		psiIf =
		(
			rho.oldTime()*psi0/deltaT
		  + Su.field()
		  - psiIf
		)/(rho/deltaT - Sp.field());
	}

	psi.correctBoundaryConditions();
}


template<class RhoType, class SpType, class SuType>
void Foam::MULES::implicitSolve
(
	const RhoType& rho,
	volScalarField& psi,
	const surfaceScalarField& phi,
	surfaceScalarField& phiPsi,
	const SpType& Sp,
	const SuType& Su,
	const scalar psiMax,
	const scalar psiMin
)
{
	const fvMesh& mesh = psi.mesh();

	const dictionary& MULEScontrols =
		mesh.solutionDict().solverDict(psi.name());

	label maxIter
	(
		readLabel(MULEScontrols.lookup("maxIter"))
	);

	label nLimiterIter
	(
		readLabel(MULEScontrols.lookup("nLimiterIter"))
	);

	scalar maxUnboundedness
	(
		readScalar(MULEScontrols.lookup("maxUnboundedness"))
	);

	scalar CoCoeff
	(
		readScalar(MULEScontrols.lookup("CoCoeff"))
	);

	scalarField allCoLambda(mesh.nFaces());

	{
		surfaceScalarField Cof =
			mesh.time().deltaT()*mesh.surfaceInterpolation::deltaCoeffs()
		   *mag(phi)/mesh.magSf();

		slicedSurfaceScalarField CoLambda
		(
			IOobject
			(
				"CoLambda",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE,
				false
			),
			mesh,
			dimless,
			allCoLambda,
			false   // Use slices for the couples
		);

		CoLambda == 1.0/max(CoCoeff*Cof, scalar(1));
	}

	scalarField allLambda(allCoLambda);
	//scalarField allLambda(mesh.nFaces(), 1.0);

	slicedSurfaceScalarField lambda
	(
		IOobject
		(
			"lambda",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE,
			false
		),
		mesh,
		dimless,
		allLambda,
		false   // Use slices for the couples
	);

	linear<scalar> CDs(mesh);
	upwind<scalar> UDs(mesh, phi);
	//fv::uncorrectedSnGrad<scalar> snGrads(mesh);

	fvScalarMatrix psiConvectionDiffusion
	(
		fvm::ddt(rho, psi)
	  + fv::gaussConvectionScheme<scalar>(mesh, phi, UDs).fvmDiv(phi, psi)
		//- fv::gaussLaplacianScheme<scalar, scalar>(mesh, CDs, snGrads)
		//.fvmLaplacian(Dpsif, psi)
	  - fvm::Sp(Sp, psi)
	  - Su
	);

	surfaceScalarField phiBD = psiConvectionDiffusion.flux();

	surfaceScalarField& phiCorr = phiPsi;
	phiCorr -= phiBD;

	for (label i=0; i<maxIter; i++)
	{
		if (i != 0 && i < 4)
		{
			allLambda = allCoLambda;
		}

		limiter
		(
			allLambda,
			rho,
			psi,
			phiBD,
			phiCorr,
			Sp.field(),
			Su.field(),
			psiMax,
			psiMin,
			nLimiterIter
		);

		solve
		(
			psiConvectionDiffusion + fvc::div(lambda*phiCorr),
			MULEScontrols
		);

		scalar maxPsiM1 = gMax(psi.internalField()) - 1.0;
		scalar minPsi = gMin(psi.internalField());

		scalar unboundedness = max(max(maxPsiM1, 0.0), -min(minPsi, 0.0));

		if (unboundedness < maxUnboundedness)
		{
			break;
		}
		else
		{
			Info<< "MULES: max(" << psi.name() << " - 1) = " << maxPsiM1
				<< " min(" << psi.name() << ") = " << minPsi << endl;

			phiBD = psiConvectionDiffusion.flux();

			/*
			word gammaScheme("div(phi,gamma)");
			word gammarScheme("div(phirb,gamma)");

			const surfaceScalarField& phir =
				mesh.lookupObject<surfaceScalarField>("phir");

			phiCorr =
				fvc::flux
				(
					phi,
					psi,
					gammaScheme
				)
			  + fvc::flux
				(
					-fvc::flux(-phir, scalar(1) - psi, gammarScheme),
					psi,
					gammarScheme
				)
				- phiBD;
			*/
		}
	}

	phiPsi = psiConvectionDiffusion.flux() + lambda*phiCorr;
}


template<class RhoType, class SpType, class SuType>
void Foam::MULES::limiter
(
	scalarField& allLambda,
	const RhoType& rho,
	const volScalarField& psi,
	const surfaceScalarField& phiBD,
	const surfaceScalarField& phiCorr,
	const SpType& Sp,
	const SuType& Su,
	const scalar psiMax,
	const scalar psiMin,
	const label nLimiterIter
)
{
	const scalarField& psiIf = psi.internalField();
	const volScalarField::Boundary& psiBf = psi.boundaryField();

	const scalarField& psi0 = psi.oldTime();

	const fvMesh& mesh = psi.mesh();

	const unallocLabelList& owner = mesh.owner();
	const unallocLabelList& neighb = mesh.neighbour();
	tmp<volScalarField::Internal> tVsc = mesh.Vsc();
	const scalarField& V = tVsc();
	const scalar deltaT = mesh.time().deltaT().value();

	const scalarField& phiBDIf = phiBD;
	const surfaceScalarField::Boundary& phiBDBf =
		phiBD.boundaryField();

	const scalarField& phiCorrIf = phiCorr;
	const surfaceScalarField::Boundary& phiCorrBf =
		phiCorr.boundaryField();

	slicedSurfaceScalarField lambda
	(
		IOobject
		(
			"lambda",
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE,
			false
		),
		mesh,
		dimless,
		allLambda,
		false   // Use slices for the couples
	);

	scalarField& lambdaIf = lambda;
	surfaceScalarField::Boundary& lambdaBf =
		lambda.boundaryField();

	scalarField psiMaxn(psiIf.size(), psiMin);
	scalarField psiMinn(psiIf.size(), psiMax);

	scalarField sumPhiBD(psiIf.size(), 0.0);

	scalarField sumPhip(psiIf.size(), VSMALL);
	scalarField mSumPhim(psiIf.size(), VSMALL);

	forAll (phiCorrIf, facei)
	{
		label own = owner[facei];
		label nei = neighb[facei];

		psiMaxn[own] = max(psiMaxn[own], psiIf[nei]);
		psiMinn[own] = min(psiMinn[own], psiIf[nei]);

		psiMaxn[nei] = max(psiMaxn[nei], psiIf[own]);
		psiMinn[nei] = min(psiMinn[nei], psiIf[own]);

		sumPhiBD[own] += phiBDIf[facei];
		sumPhiBD[nei] -= phiBDIf[facei];

		scalar phiCorrf = phiCorrIf[facei];

		if (phiCorrf > 0.0)
		{
			sumPhip[own] += phiCorrf;
			mSumPhim[nei] += phiCorrf;
		}
		else
		{
			mSumPhim[own] -= phiCorrf;
			sumPhip[nei] -= phiCorrf;
		}
	}

	forAll (phiCorrBf, patchi)
	{
		const fvPatchScalarField& psiPf = psiBf[patchi];
		const scalarField& phiBDPf = phiBDBf[patchi];
		const scalarField& phiCorrPf = phiCorrBf[patchi];

		const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

		if (psiPf.coupled())
		{
			scalarField psiPNf = psiPf.patchNeighbourField();

			forAll (phiCorrPf, pFacei)
			{
				label pfCelli = pFaceCells[pFacei];

				psiMaxn[pfCelli] = max(psiMaxn[pfCelli], psiPNf[pFacei]);
				psiMinn[pfCelli] = min(psiMinn[pfCelli], psiPNf[pFacei]);
			}
		}
		else
		{
			forAll (phiCorrPf, pFacei)
			{
				label pfCelli = pFaceCells[pFacei];

				psiMaxn[pfCelli] = max(psiMaxn[pfCelli], psiPf[pFacei]);
				psiMinn[pfCelli] = min(psiMinn[pfCelli], psiPf[pFacei]);
			}
		}

		forAll (phiCorrPf, pFacei)
		{
			label pfCelli = pFaceCells[pFacei];

			sumPhiBD[pfCelli] += phiBDPf[pFacei];

			scalar phiCorrf = phiCorrPf[pFacei];

			if (phiCorrf > 0.0)
			{
				sumPhip[pfCelli] += phiCorrf;
			}
			else
			{
				mSumPhim[pfCelli] -= phiCorrf;
			}
		}
	}

	psiMaxn = min(psiMaxn, psiMax);
	psiMinn = max(psiMinn, psiMin);

	//scalar smooth = 0.5;
	//psiMaxn = min((1.0 - smooth)*psiIf + smooth*psiMaxn, psiMax);
	//psiMinn = max((1.0 - smooth)*psiIf + smooth*psiMinn, psiMin);

	if (mesh.moving())
	{
		tmp<volScalarField::Internal> V0 = mesh.Vsc0();

		psiMaxn =
			V*((rho/deltaT - Sp)*psiMaxn - Su)
		  - (V0()/deltaT)*rho.oldTime()*psi0
		  + sumPhiBD;

		psiMinn =
			V*(Su - (rho/deltaT - Sp)*psiMinn)
		  + (V0/deltaT)*rho.oldTime()*psi0
		  - sumPhiBD;
	}
	else
	{
		psiMaxn =
			V*((rho/deltaT - Sp)*psiMaxn - (rho.oldTime()/deltaT)*psi0 - Su)
		  + sumPhiBD;

		psiMinn =
			V*((rho/deltaT)*psi0 - (rho.oldTime()/deltaT - Sp)*psiMinn + Su)
		  - sumPhiBD;
	}

	scalarField sumlPhip(psiIf.size());
	scalarField mSumlPhim(psiIf.size());

	for(int j=0; j<nLimiterIter; j++)
	{
		sumlPhip = 0.0;
		mSumlPhim = 0.0;

		forAll (lambdaIf, facei)
		{
			label own = owner[facei];
			label nei = neighb[facei];

			scalar lambdaPhiCorrf = lambdaIf[facei]*phiCorrIf[facei];

			if (lambdaPhiCorrf > 0.0)
			{
				sumlPhip[own] += lambdaPhiCorrf;
				mSumlPhim[nei] += lambdaPhiCorrf;
			}
			else
			{
				mSumlPhim[own] -= lambdaPhiCorrf;
				sumlPhip[nei] -= lambdaPhiCorrf;
			}
		}

		forAll (lambdaBf, patchi)
		{
			scalarField& lambdaPf = lambdaBf[patchi];
			const scalarField& phiCorrfPf = phiCorrBf[patchi];

			const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

			forAll (lambdaPf, pFacei)
			{
				label pfCelli = pFaceCells[pFacei];

				scalar lambdaPhiCorrf = lambdaPf[pFacei]*phiCorrfPf[pFacei];

				if (lambdaPhiCorrf > 0.0)
				{
					sumlPhip[pfCelli] += lambdaPhiCorrf;
				}
				else
				{
					mSumlPhim[pfCelli] -= lambdaPhiCorrf;
				}
			}
		}

		forAll (sumlPhip, celli)
		{
			sumlPhip[celli] =
				max(min
				(
					(sumlPhip[celli] + psiMaxn[celli])/mSumPhim[celli],
					1.0), 0.0
				);

			mSumlPhim[celli] =
				max(min
				(
					(mSumlPhim[celli] + psiMinn[celli])/sumPhip[celli],
					1.0), 0.0
				);
		}

		const scalarField& lambdam = sumlPhip;
		const scalarField& lambdap = mSumlPhim;

		forAll (lambdaIf, facei)
		{
			if (phiCorrIf[facei] > 0.0)
			{
				lambdaIf[facei] = min
				(
					lambdaIf[facei],
					min(lambdap[owner[facei]], lambdam[neighb[facei]])
				);
			}
			else
			{
				lambdaIf[facei] = min
				(
					lambdaIf[facei],
					min(lambdam[owner[facei]], lambdap[neighb[facei]])
				);
			}
		}


		forAll (lambdaBf, patchi)
		{
			fvsPatchScalarField& lambdaPf = lambdaBf[patchi];
			const scalarField& phiCorrfPf = phiCorrBf[patchi];

			const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

			forAll (lambdaPf, pFacei)
			{
				label pfCelli = pFaceCells[pFacei];

				if (phiCorrfPf[pFacei] > 0.0)
				{
					lambdaPf[pFacei] = min(lambdaPf[pFacei], lambdap[pfCelli]);
				}
				else
				{
					lambdaPf[pFacei] = min(lambdaPf[pFacei], lambdam[pfCelli]);
				}
			}
		}

		syncTools::syncFaceList(mesh, allLambda, minEqOp<scalar>(), false);
	}
}


// ************************************************************************* //
