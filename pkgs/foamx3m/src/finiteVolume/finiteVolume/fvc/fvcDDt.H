


#ifndef fvcDDt_H
#define fvcDDt_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


namespace fvc
{
	template<class Type>
	tmp<GeometricField<Type, fvPatchField, volMesh> > DDt
	(
		const surfaceScalarField&,
		const GeometricField<Type, fvPatchField, volMesh>&
	);

	template<class Type>
	tmp<GeometricField<Type, fvPatchField, volMesh> > DDt
	(
		const tmp<surfaceScalarField>&,
		const GeometricField<Type, fvPatchField, volMesh>&
	);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "fvcDDt.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
