/*-------------------------------------------------------------------------*\

 Compute interface face areas

 Copyright (C) 2014-2020  Mosayeb Shams
 Copyright (C) 2017-2020  Ali Qaseminejad Raeini 

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
\*-------------------------------------------------------------------------*/

//! This file is part of calcCurvatureForce.C
 



surfaceVectorField interfaceSf
	(	IOobject (   "interfaceSf", timeName(), msh ),
		msh, dimensionedVector("interfaceSf", dimArea, vector::zero)
	);



// calculate interface normal vectors, interfaceSf
forAll(magDelS,fI)
{
	if (magDelS[fI] > CONTRAST_FACTOR)
	{
		// defining the right hand rule coeficient, such that interfaceSf_ always points from 0(Fluid1) to 1(Fluid2)
		scalar rightHandRuleCoef =deltaS_[fI]/magDelS[fI];

		// return the points of fI
		const face& facePoints = faces[fI];

		pointField 	xInterface(facePoints.size());

		forAll(facePoints,pointI)
		{
			 label ip1(facePoints[pointI]);
			 xInterface[pointI] = points[ip1] + distPointInterface_[ip1];
		}

		label nPoints = xInterface.size();

		if (nPoints == 3)
		{
			interfaceSf[fI] = rightHandRuleCoef*triPointRef
				(
					xInterface[0],
					xInterface[1],
					xInterface[2]
				).normal();
		}
		else
		{
			point centrePoint = point::zero;

			for(Foam::label pointI=0;pointI<nPoints;++pointI)
			{
				  centrePoint += xInterface[pointI];
			}

			centrePoint /= nPoints;

			for (Foam::label pointI=0; pointI<nPoints; ++pointI)
			{
				const point& currentPoint = xInterface[pointI];
				const point& nextPoint = xInterface[(pointI + 1) % nPoints];

					interfaceSf[fI] += rightHandRuleCoef*triPointRef
					  (
							currentPoint,
							nextPoint,
							centrePoint
					  ).normal();
			}
		}
	}
}

	forAll(patches, bI)
	{
		const scalarField & pMagDelS = magDelS.boundaryField()[bI];
		const scalarField & pSnGradAlphaSharp = deltaS_.boundaryField()[bI];
		vectorField & pInterfaceSf = interfaceSf.boundaryField()[bI];
		if (patches[bI].coupled())
		{
			forAll(pMagDelS, pfI)
			{
			  if (pMagDelS[pfI] > CONTRAST_FACTOR)
			  {
					// defining the right hand rule coeficient, such that interfaceSf_ always points from 0(Fluid1) to 1(Fluid2)
					scalar rightHandRuleCoef = pSnGradAlphaSharp[pfI]/pMagDelS[pfI];

					// return the points of face pfI
					const face& facePoints = faces[pfI+patches[bI].patch().start()];

					pointField 	xInterface(facePoints.size());

					forAll(facePoints,pointI)
					{
						 label ip1(facePoints[pointI]);
						 xInterface[pointI] = points[ip1] + distPointInterface_[ip1];
					}

					const label nPoints = xInterface.size();

					if (nPoints == 3)
					{
						 pInterfaceSf[pfI] = rightHandRuleCoef*triPointRef
						 (
							  xInterface[0],
							  xInterface[1],
							  xInterface[2]
						 ).normal();
					}
					else
					{
						 point centrePoint = point::zero;

						 for(Foam::label pointI=0;pointI<nPoints;++pointI)
						 {
							  centrePoint += xInterface[pointI];
						 }

						 centrePoint /= nPoints;

						 for (Foam::label pointI=0; pointI<nPoints; ++pointI)
						 {
							  const point& currentPoint = xInterface[pointI];
							  const point& nextPoint = xInterface[(pointI + 1) % nPoints];

							  pInterfaceSf[pfI] += rightHandRuleCoef*triPointRef
							  (
									currentPoint,
									nextPoint,
									centrePoint
							  ).normal();
						 }
					}
				}
			}

		}
		else		 {			 }
	}

	//interfaceSf=(interfaceSf&curvatureNSf)*curvatureNSf;
	//interfaceSf=msh.magSf()/(curvatureNSf&msh.Sf()+1.0e-36)*msh.Sf();


