

#include "fvcDDt.H"
#include "fvcDiv.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
DDt
(
	const surfaceScalarField& phi,
	const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
	tmp<GeometricField<Type, fvPatchField, volMesh> > ddtDivPhiPsi
		= fvc::ddt(psi) + fvc::div(phi, psi);

	if (phi.mesh().moving())
	{
		return ddtDivPhiPsi - fvc::div(phi + phi.mesh().phi())*psi;
	}
	else
	{
		return ddtDivPhiPsi - fvc::div(phi)*psi;
	}
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
DDt
(
	const tmp<surfaceScalarField>& tphi,
	const GeometricField<Type, fvPatchField, volMesh>& psi
)
{
	tmp<GeometricField<Type, fvPatchField, volMesh> > DDtPsi
	(
		fvc::DDt(tphi(), psi)
	);
	tphi.clear();
	return DDtPsi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
