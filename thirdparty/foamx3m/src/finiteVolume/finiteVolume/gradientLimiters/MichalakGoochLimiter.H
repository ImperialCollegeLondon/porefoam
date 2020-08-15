
#ifndef MichalakGoochLimiter_H
#define MichalakGoochLimiter_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class MichalakGoochLimiter
{
	// Private data

		//- Coefficient
		const scalar k_;


public:

	// Constructors

		//- Construct null
		MichalakGoochLimiter()
		:
			k_(1.5)
		{}


	// Destructor - default


	// Member functions

		//- Set scalar limiter value
		inline void limiter
		(
			scalar& lim,
			const scalar& cellVolume,
			const scalar& deltaOneMax,
			const scalar& deltaOneMin,
			const scalar& extrapolate
		)
		{
			// If there is no span or extrapolation, do not limit.  Note that
			// deltaOneMax and deltaOneMin are bound by zero and min is
			// negative
			if
			(
				deltaOneMax - deltaOneMin < SMALL
			 || mag(extrapolate) < SMALL
			)
			{
				return;
			}

			scalar y;

			if (extrapolate > 0)
			{
				y = deltaOneMax/extrapolate;
			}
			else
			{
				y = deltaOneMin/extrapolate;
			}

			if (y < k_)
			{
				const scalar C2 = (3 - 2*k_)/sqr(k_);
				const scalar C3 = -1/(3*sqr(k_)) - 2/(3*k_)*C2;

				lim = min(lim, y + C2*sqr(y) + C3*pow3(y));
			}
		}

		//- Set Type limiter
		template<class Type>
		inline void limiter
		(
			Type& lim,
			const scalar& cellVolume,
			const Type& deltaOneMax,
			const Type& deltaOneMin,
			const Type& extrapolate
		)
		{
			for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
			{
				limiter
				(
					lim.component(cmpt),
					cellVolume,
					deltaOneMax.component(cmpt),
					deltaOneMin.component(cmpt),
					extrapolate.component(cmpt)
				);
			}
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
