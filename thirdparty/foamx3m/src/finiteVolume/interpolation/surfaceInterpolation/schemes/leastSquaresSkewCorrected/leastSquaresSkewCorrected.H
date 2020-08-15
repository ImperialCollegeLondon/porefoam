
#ifndef leastSquaresSkewCorrected_H
#define leastSquaresSkewCorrected_H

#include "surfaceInterpolationScheme.H"
#include "skewCorrectionVectors.H"
#include "linear.H"
#include "gaussGrad.H"
#include "leastSquaresGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class leastSquaresSkewCorrected
:
	public surfaceInterpolationScheme<Type>
{
	// Private member data

		tmp<surfaceInterpolationScheme<Type> > tScheme_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		leastSquaresSkewCorrected(const leastSquaresSkewCorrected&);

		//- Disallow default bitwise assignment
		void operator=(const leastSquaresSkewCorrected&);


public:

	//- Runtime type information
	TypeName("leastSquaresSkewCorrected");


	// Constructors

		//- Construct from mesh and Istream
		leastSquaresSkewCorrected
		(
			const fvMesh& mesh,
			Istream& is
		)
		:
			surfaceInterpolationScheme<Type>(mesh),
			tScheme_
			(
				surfaceInterpolationScheme<Type>::New(mesh, is)
			)
		{}


		//- Construct from mesh, faceFlux and Istream
		leastSquaresSkewCorrected
		(
			const fvMesh& mesh,
			const surfaceScalarField& faceFlux,
			Istream& is
		)
		:
			surfaceInterpolationScheme<Type>(mesh),
			tScheme_
			(
				surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
			)
		{}


	// Member Functions

		//- Return the interpolation weighting factors
		tmp<surfaceScalarField> weights
		(
			const GeometricField<Type, fvPatchField, volMesh>& vf
		) const
		{
			return tScheme_().weights(vf);
		}

		//- Return true if this scheme uses an explicit correction
		virtual bool corrected() const
		{
			return
				tScheme_().corrected()
			 || skewCorrectionVectors::New(this->mesh()).skew();
		}

		tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		skewCorrection
		(
			const GeometricField<Type, fvPatchField, volMesh>& vf
		) const
		{
			const fvMesh& mesh = this->mesh();

			const skewCorrectionVectors& scv =
				skewCorrectionVectors::New(mesh);

			tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr
			(
				new GeometricField<Type, fvsPatchField, surfaceMesh>
				(
					IOobject
					(
						vf.name(),
						mesh.time().timeName(),
						mesh
					),
					mesh,
					dimensioned<Type>
					(
						vf.name(),
						vf.dimensions(),
						pTraits<Type>::zero
					)
				)
			);

			for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
			{
				tsfCorr().replace
				(
					cmpt,
					scv() & linear
					<
						typename outerProduct
						<
							vector,
							typename pTraits<Type>::cmptType
						>::type
					> (mesh).interpolate
					(
						fv::leastSquaresGrad<typename pTraits<Type>::cmptType>
						(mesh).grad(vf.component(cmpt))
					)
				);
			}

			return tsfCorr;
		}


		//- Return the explicit correction to the face-interpolate
		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		correction
		(
			const GeometricField<Type, fvPatchField, volMesh>& vf
		) const
		{
			if
			(
				tScheme_().corrected()
			 && skewCorrectionVectors::New(this->mesh()).skew()
			)
			{
				return tScheme_().correction(vf) + skewCorrection(vf);
			}
			else if (tScheme_().corrected())
			{
				return tScheme_().correction(vf);
			}
			else if (skewCorrectionVectors::New(this->mesh()).skew())
			{
				return skewCorrection(vf);
			}
			else
			{
				return
					tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >(NULL);
			}
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
