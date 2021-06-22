/*-------------------------------------------------------------------------*\

 Correct interface at contact lines 

// Copyright (C) 2014-2020  Mosayeb Shams
// Copyright (C) 2017-2020  Ali Qaseminejad Raeini 

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
\*-------------------------------------------------------------------------*/


namespace Foam
{
	void correctForContactAngle
	(
		const volScalarField& alpha1,
		const volVectorField U,
		const surfaceVectorField nw_,
		pointVectorField& gradAlphaP,
		const labelList& interfPoints
	)
	{

		const polyBoundaryMesh& boundaryMesh = alpha1.mesh().boundaryMesh();
		const fvBoundaryMesh& patches = alpha1.mesh().boundary();
		forAll(boundaryMesh, bi)
		{
			const labelList& patchPoints = boundaryMesh[bi].meshPoints();
			if (patches[bi].coupled()) {}
			else	if (isA<alphaContactAngleFvPatchScalarField>(alpha1.boundaryField()[bi]))
			{
				primitivePatchInterpolation pinterpolator(boundaryMesh[bi]);
				const Foam::scalar convertToRad =  Foam::mathematicalConstant::pi/180.0;

				pointField nWp= pinterpolator.faceToPointInterpolate(nw_.boundaryField()[bi]);
				const  alphaContactAngleFvPatchScalarField& acap =  refCast<const alphaContactAngleFvPatchScalarField> ( alpha1.boundaryField()[bi] );
				vectorField Uprelative= pinterpolator.faceToPointInterpolate( U.internalField()[bi]-U.boundaryField()[bi] );
				vectorField nHatp=0.0*nWp;
				forAll(patchPoints, ppI) nHatp[ppI] = gradAlphaP[patchPoints[ppI]];
				nHatp /= mag(nHatp) + 1e-37;
				scalarField theta = convertToRad*acap.theta(Uprelative, nHatp, nWp);

				vectorField ss=nHatp-(nHatp & nWp)*nWp;
				ss/=mag(ss)+1e-37;

				forAll(patchPoints, ppI)
				{
					label pI=patchPoints[ppI];
					if (interfPoints[pI])
					{
						gradAlphaP[pI] = mag(gradAlphaP[pI]) * ( sin(theta[0])*ss[ppI]+cos(theta[0])*nWp[ppI] );
					}
				}
			}
		}
	}
}
