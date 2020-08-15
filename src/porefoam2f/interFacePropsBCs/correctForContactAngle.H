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
		forAll(boundaryMesh, patchi)
		{
			const labelList& patchPoints = boundaryMesh[patchi].meshPoints();
			if (patches[patchi].coupled()) {}
			else	if (isA<alphaContactAngleFvPatchScalarField>(alpha1.boundaryField()[patchi]))
			{
				primitivePatchInterpolation pinterpolator(boundaryMesh[patchi]);
				const Foam::scalar convertToRad =  Foam::mathematicalConstant::pi/180.0;

				pointField nWp= pinterpolator.faceToPointInterpolate(nw_.boundaryField()[patchi]);
				const  alphaContactAngleFvPatchScalarField& acap =  refCast<const alphaContactAngleFvPatchScalarField> ( alpha1.boundaryField()[patchi] );
				vectorField Uprelative= pinterpolator.faceToPointInterpolate( U.internalField()[patchi]-U.boundaryField()[patchi] );
				vectorField nHatp=0.0*nWp;
				forAll(patchPoints, ppI) nHatp[ppI] = gradAlphaP[patchPoints[ppI]];
				nHatp /= mag(nHatp) + 1.0e-37;
				scalarField theta = convertToRad*acap.theta(Uprelative, nHatp, nWp);

				vectorField ss=nHatp-(nHatp & nWp)*nWp;
				ss/=mag(ss)+1.0e-37;

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
