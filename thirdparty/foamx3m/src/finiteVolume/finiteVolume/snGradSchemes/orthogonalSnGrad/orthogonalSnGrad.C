
#include "orthogonalSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
orthogonalSnGrad<Type>::~orthogonalSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
orthogonalSnGrad<Type>::correction
(
	const GeometricField<Type, fvPatchField, volMesh>&
) const
{
	notImplemented
	(
		"orthogonalSnGrad<Type>::correction"
		"(const GeometricField<Type, fvPatchField, volMesh>&)"
	);
	return tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >(NULL);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
