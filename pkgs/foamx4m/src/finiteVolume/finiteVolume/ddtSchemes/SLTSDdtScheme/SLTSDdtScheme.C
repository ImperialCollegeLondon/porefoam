
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

#include "SLTSDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void SLTSDdtScheme<Type>::relaxedDiag
(
	scalarField& rD,
	const surfaceScalarField& phi
) const
{
	const unallocLabelList& owner = mesh().owner();
	const unallocLabelList& neighbour = mesh().neighbour();
	scalarField diag(rD.size(), 0.0);

	forAll(owner, faceI)
	{
		if (phi[faceI] > 0.0)
		{
			diag[owner[faceI]] += phi[faceI];
			rD[neighbour[faceI]] += phi[faceI];
		}
		else
		{
			diag[neighbour[faceI]] -= phi[faceI];
			rD[owner[faceI]] -= phi[faceI];
		}
	}

	forAll(phi.boundaryField(), patchi)
	{
		const fvsPatchScalarField& pphi = phi.boundaryField()[patchi];
		const unallocLabelList& faceCells = pphi.patch().patch().faceCells();

		forAll(pphi, patchFacei)
		{
			if (pphi[patchFacei] > 0.0)
			{
				diag[faceCells[patchFacei]] += pphi[patchFacei];
			}
			else
			{
				rD[faceCells[patchFacei]] -= pphi[patchFacei];
			}
		}
	}

	rD += (1.0/alpha_ - 2.0)*diag;
}


template<class Type>
tmp<volScalarField> SLTSDdtScheme<Type>::SLrDeltaT() const
{
	const surfaceScalarField& phi =
		mesh().objectRegistry::template
		lookupObject<surfaceScalarField>(phiName_);

	const dimensionedScalar& deltaT = mesh().time().deltaT();

	tmp<volScalarField> trDeltaT
	(
		new volScalarField
		(
			IOobject
			(
				"rDeltaT",
				phi.instance(),
				mesh()
			),
			mesh(),
			dimensionedScalar("rDeltaT", dimless/dimTime, 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);

	volScalarField& rDeltaT = trDeltaT();

	relaxedDiag(rDeltaT, phi);

	if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
	{
		rDeltaT.internalField() = max
		(
			rDeltaT.internalField()/mesh().V(),
			scalar(1)/deltaT.value()
		);
	}
	else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
	{
		const volScalarField& rho =
			mesh().objectRegistry::template
			lookupObject<volScalarField>(rhoName_)
		   .oldTime();

		rDeltaT.internalField() = max
		(
			rDeltaT.internalField()/(rho.internalField()*mesh().V()),
			scalar(1)/deltaT.value()
		);
	}
	else
	{
		FatalErrorIn("SLTSDdtScheme<Type>::CofrDeltaT() const")
			<< "Incorrect dimensions of phi: " << phi.dimensions()
			<< abort(FatalError);
	}

	rDeltaT.correctBoundaryConditions();

	return trDeltaT;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
SLTSDdtScheme<Type>::fvcDdt
(
	const dimensioned<Type>& dt
)
{
	volScalarField rDeltaT = SLrDeltaT();

	IOobject ddtIOobject
	(
		"ddt("+dt.name()+')',
		mesh().time().timeName(),
		mesh()
	);

	if (mesh().moving())
	{
		tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt
		(
			new GeometricField<Type, fvPatchField, volMesh>
			(
				ddtIOobject,
				mesh(),
				dimensioned<Type>
				(
					"0",
					dt.dimensions()/dimTime,
					pTraits<Type>::zero
				)
			)
		);

		tdtdt().internalField() =
			rDeltaT.internalField()*dt.value()*(1.0 - mesh().V0()/mesh().V());

		return tdtdt;
	}
	else
	{
		return tmp<GeometricField<Type, fvPatchField, volMesh> >
		(
			new GeometricField<Type, fvPatchField, volMesh>
			(
				ddtIOobject,
				mesh(),
				dimensioned<Type>
				(
					"0",
					dt.dimensions()/dimTime,
					pTraits<Type>::zero
				),
				calculatedFvPatchField<Type>::typeName
			)
		);
	}
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
SLTSDdtScheme<Type>::fvcDdt
(
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	volScalarField rDeltaT = SLrDeltaT();

	IOobject ddtIOobject
	(
		"ddt("+vf.name()+')',
		mesh().time().timeName(),
		mesh()
	);

	if (mesh().moving())
	{
		return tmp<GeometricField<Type, fvPatchField, volMesh> >
		(
			new GeometricField<Type, fvPatchField, volMesh>
			(
				ddtIOobject,
				mesh(),
				rDeltaT.dimensions()*vf.dimensions(),
				rDeltaT.internalField()*
				(
					vf.internalField()
				  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
				),
				rDeltaT.boundaryField()*
				(
					vf.boundaryField() - vf.oldTime().boundaryField()
				)
			)
		);
	}
	else
	{
		return tmp<GeometricField<Type, fvPatchField, volMesh> >
		(
			new GeometricField<Type, fvPatchField, volMesh>
			(
				ddtIOobject,
				rDeltaT*(vf - vf.oldTime())
			)
		);
	}
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
SLTSDdtScheme<Type>::fvcDdt
(
	const dimensionedScalar& rho,
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	volScalarField rDeltaT = SLrDeltaT();

	IOobject ddtIOobject
	(
		"ddt("+rho.name()+','+vf.name()+')',
		mesh().time().timeName(),
		mesh()
	);

	if (mesh().moving())
	{
		return tmp<GeometricField<Type, fvPatchField, volMesh> >
		(
			new GeometricField<Type, fvPatchField, volMesh>
			(
				ddtIOobject,
				mesh(),
				rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
				rDeltaT.internalField()*rho.value()*
				(
					vf.internalField()
				  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
				),
				rDeltaT.boundaryField()*rho.value()*
				(
					vf.boundaryField() - vf.oldTime().boundaryField()
				)
			)
		);
	}
	else
	{
		return tmp<GeometricField<Type, fvPatchField, volMesh> >
		(
			new GeometricField<Type, fvPatchField, volMesh>
			(
				ddtIOobject,
				rDeltaT*rho*(vf - vf.oldTime())
			)
		);
	}
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
SLTSDdtScheme<Type>::fvcDdt
(
	const volScalarField& rho,
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	volScalarField rDeltaT = SLrDeltaT();

	IOobject ddtIOobject
	(
		"ddt("+rho.name()+','+vf.name()+')',
		mesh().time().timeName(),
		mesh()
	);

	if (mesh().moving())
	{
		return tmp<GeometricField<Type, fvPatchField, volMesh> >
		(
			new GeometricField<Type, fvPatchField, volMesh>
			(
				ddtIOobject,
				mesh(),
				rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
				rDeltaT.internalField()*
				(
					rho.internalField()*vf.internalField()
				  - rho.oldTime().internalField()
				   *vf.oldTime().internalField()*mesh().V0()/mesh().V()
				),
				rDeltaT.boundaryField()*
				(
					rho.boundaryField()*vf.boundaryField()
				  - rho.oldTime().boundaryField()
				   *vf.oldTime().boundaryField()
				)
			)
		);
	}
	else
	{
		return tmp<GeometricField<Type, fvPatchField, volMesh> >
		(
			new GeometricField<Type, fvPatchField, volMesh>
			(
				ddtIOobject,
				rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
			)
		);
	}
}


template<class Type>
tmp<fvMatrix<Type> >
SLTSDdtScheme<Type>::fvmDdt
(
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	tmp<fvMatrix<Type> > tfvm
	(
		new fvMatrix<Type>
		(
			vf,
			vf.dimensions()*dimVol/dimTime
		)
	);

	fvMatrix<Type>& fvm = tfvm();

	scalarField rDeltaT = SLrDeltaT()().internalField();

	Info<< "max/min rDeltaT " << gMax(rDeltaT) << " " << gMin(rDeltaT) << endl;

	fvm.diag() = rDeltaT*mesh().V();

	if (mesh().moving())
	{
		fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V0();
	}
	else
	{
		fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V();
	}

	return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
SLTSDdtScheme<Type>::fvmDdt
(
	const dimensionedScalar& rho,
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	tmp<fvMatrix<Type> > tfvm
	(
		new fvMatrix<Type>
		(
			vf,
			rho.dimensions()*vf.dimensions()*dimVol/dimTime
		)
	);
	fvMatrix<Type>& fvm = tfvm();

	scalarField rDeltaT = SLrDeltaT()().internalField();

	fvm.diag() = rDeltaT*rho.value()*mesh().V();

	if (mesh().moving())
	{
		fvm.source() = rDeltaT
			*rho.value()*vf.oldTime().internalField()*mesh().V0();
	}
	else
	{
		fvm.source() = rDeltaT
			*rho.value()*vf.oldTime().internalField()*mesh().V();
	}

	return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
SLTSDdtScheme<Type>::fvmDdt
(
	const volScalarField& rho,
	const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
	tmp<fvMatrix<Type> > tfvm
	(
		new fvMatrix<Type>
		(
			vf,
			rho.dimensions()*vf.dimensions()*dimVol/dimTime
		)
	);
	fvMatrix<Type>& fvm = tfvm();

	scalarField rDeltaT = SLrDeltaT()().internalField();

	fvm.diag() = rDeltaT*rho.internalField()*mesh().V();

	if (mesh().moving())
	{
		fvm.source() = rDeltaT
			*rho.oldTime().internalField()
			*vf.oldTime().internalField()*mesh().V0();
	}
	else
	{
		fvm.source() = rDeltaT
			*rho.oldTime().internalField()
			*vf.oldTime().internalField()*mesh().V();
	}

	return tfvm;
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtPhiCorr
(
	const volScalarField& rA,
	const GeometricField<Type, fvPatchField, volMesh>& U,
	const fluxFieldType& phi
)
{
	IOobject ddtIOobject
	(
		"ddtPhiCorr(" + rA.name() + ',' + U.name() + ',' + phi.name() + ')',
		mesh().time().timeName(),
		mesh()
	);

	if (mesh().moving())
	{
		return tmp<fluxFieldType>
		(
			new fluxFieldType
			(
				ddtIOobject,
				mesh(),
				dimensioned<typename flux<Type>::type>
				(
					"0",
					rA.dimensions()*phi.dimensions()/dimTime,
					pTraits<typename flux<Type>::type>::zero
				)
			)
		);
	}
	else
	{
		volScalarField rDeltaT = SLrDeltaT();

		return tmp<fluxFieldType>
		(
			new fluxFieldType
			(
				ddtIOobject,
				this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())*
				(
					fvc::interpolate(rDeltaT*rA)*phi.oldTime()
				  - (fvc::interpolate(rDeltaT*rA*U.oldTime()) & mesh().Sf())
				)
			)
		);
	}
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtPhiCorr
(
	const volScalarField& rA,
	const volScalarField& rho,
	const GeometricField<Type, fvPatchField, volMesh>& U,
	const fluxFieldType& phi
)
{
	IOobject ddtIOobject
	(
		"ddtPhiCorr("
	  + rA.name() + ',' + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
		mesh().time().timeName(),
		mesh()
	);

	if (mesh().moving())
	{
		return tmp<fluxFieldType>
		(
			new fluxFieldType
			(
				ddtIOobject,
				mesh(),
				dimensioned<typename flux<Type>::type>
				(
					"0",
					rA.dimensions()*phi.dimensions()/dimTime,
					pTraits<typename flux<Type>::type>::zero
				)
			)
		);
	}
	else
	{
		volScalarField rDeltaT = SLrDeltaT();

		if
		(
			U.dimensions() == dimVelocity
		 && phi.dimensions() == dimVelocity*dimArea
		)
		{
			return tmp<fluxFieldType>
			(
				new fluxFieldType
				(
					ddtIOobject,
					this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
				   *(
					    fvc::interpolate(rDeltaT*rA*rho.oldTime())*phi.oldTime()
					  - (fvc::interpolate(rDeltaT*rA*rho.oldTime()*U.oldTime())
					  & mesh().Sf())
					)
				)
			);
		}
		else if
		(
			U.dimensions() == dimVelocity
		 && phi.dimensions() == dimDensity*dimVelocity*dimArea
		)
		{
			return tmp<fluxFieldType>
			(
				new fluxFieldType
				(
					ddtIOobject,
					this->fvcDdtPhiCoeff
					(
					    U.oldTime(),
					    phi.oldTime()/fvc::interpolate(rho.oldTime())
					)
				   *(
					    fvc::interpolate(rDeltaT*rA*rho.oldTime())
					   *phi.oldTime()/fvc::interpolate(rho.oldTime())
					  - (
					        fvc::interpolate
					        (
					            rDeltaT*rA*rho.oldTime()*U.oldTime()
					        ) & mesh().Sf()
					    )
					)
				)
			);
		}
		else if
		(
			U.dimensions() == dimDensity*dimVelocity
		 && phi.dimensions() == dimDensity*dimVelocity*dimArea
		)
		{
			return tmp<fluxFieldType>
			(
				new fluxFieldType
				(
					ddtIOobject,
					this->fvcDdtPhiCoeff(rho.oldTime(), U.oldTime(), phi.oldTime())
				   *(
					    fvc::interpolate(rDeltaT*rA)*phi.oldTime()
					  - (
					        fvc::interpolate(rDeltaT*rA*U.oldTime())&mesh().Sf()
					    )
					)
				)
			);
		}
		else
		{
			FatalErrorIn
			(
				"SLTSDdtScheme<Type>::fvcDdtPhiCorr"
			)   << "dimensions of phi are not correct"
				<< abort(FatalError);

			return fluxFieldType::null();
		}
	}
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtConsistentPhiCorr
(
	const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
	const GeometricField<Type, fvPatchField, volMesh>& U,
	const surfaceScalarField& rAUf
)
{
	tmp<fluxFieldType> toldTimeFlux =
		(mesh().Sf() & faceU.oldTime())*rAUf*fvc::interpolate(SLrDeltaT());

	if (mesh().moving())
	{
		// Mesh is moving, need to take into account the ratio between old and
		// current cell volumes
		volScalarField V0ByV
		(
			IOobject
			(
				"V0ByV",
				mesh().time().timeName(),
				mesh(),
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh(),
			dimensionedScalar("one", dimless, 1.0),
			zeroGradientFvPatchScalarField::typeName
		);
		V0ByV.internalField() = mesh().V0()/mesh().V();
		V0ByV.correctBoundaryConditions();

		// Correct the flux with interpolated volume ratio
		toldTimeFlux() *= fvc::interpolate(V0ByV);
	}

	return toldTimeFlux;
}


template<class Type>
tmp<surfaceScalarField> SLTSDdtScheme<Type>::meshPhi
(
	const GeometricField<Type, fvPatchField, volMesh>&
)
{
	return tmp<surfaceScalarField>
	(
		new surfaceScalarField
		(
			IOobject
			(
				"meshPhi",
				mesh().time().timeName(),
				mesh()
			),
			mesh(),
			dimensionedScalar("0", dimVolume/dimTime, 0.0)
		)
	);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
