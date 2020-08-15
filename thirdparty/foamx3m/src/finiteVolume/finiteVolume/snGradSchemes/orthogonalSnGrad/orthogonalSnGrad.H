
#ifndef orthogonalSnGrad_H
#define orthogonalSnGrad_H

#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


template<class Type>
class orthogonalSnGrad
:
	public snGradScheme<Type>
{
	// Private Member Functions

		//- Disallow default bitwise assignment
		void operator=(const orthogonalSnGrad&);


public:

	//- Runtime type information
	TypeName("orthogonal");


	// Constructors

		//- Construct from mesh
		orthogonalSnGrad(const fvMesh& mesh)
		:
			snGradScheme<Type>(mesh)
		{}


		//- Construct from mesh and data stream
		orthogonalSnGrad(const fvMesh& mesh, Istream&)
		:
			snGradScheme<Type>(mesh)
		{}


	//- Destructor
	virtual ~orthogonalSnGrad();


	// Member Functions

		//- Return the interpolation weighting factors for the given field
		virtual tmp<surfaceScalarField> deltaCoeffs
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) const
		{
			return this->mesh().deltaCoeffs();
		}

		//- Return true if this scheme uses an explicit correction
		virtual bool corrected() const
		{
			return false;
		}

		//- Return the explicit correction to the orthogonalSnGrad
		//  for the given field
		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "orthogonalSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
