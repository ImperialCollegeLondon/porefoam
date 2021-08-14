
#ifndef WangLimiter_H
#define WangLimiter_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class WangLimiter
{
	// Private data

		//- Limiter value
		scalar epsilonPrime_;


public:

	// Constructor

		//- Construct null
		WangLimiter()
		:
			epsilonPrime_(0.05)
		{}

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
			scalar epsilonSquare =
				sqr(epsilonPrime_*(deltaOneMax - deltaOneMin));

			if (extrapolate > VSMALL)
			{
				lim = max
				(
					0,
					min
					(
						(
							(sqr(deltaOneMax) + epsilonSquare)*extrapolate
						  + 2*sqr(extrapolate)*deltaOneMax
						)/
						stabilise
						(
							extrapolate*
							(
								sqr(deltaOneMax)
							  + 2*sqr(extrapolate)
							  + deltaOneMax*extrapolate
							  + epsilonSquare
							),
							SMALL
						),
						lim
					)
				);
			}
			else if (extrapolate < VSMALL)
			{
				lim = max
				(
					0,
					min
					(
						(
							(sqr(deltaOneMin) + epsilonSquare)*extrapolate
						  + 2*sqr(extrapolate)*deltaOneMin
						)/
						stabilise
						(
							extrapolate*
							(
								sqr(deltaOneMin)
							  + 2*sqr(extrapolate)
							  + deltaOneMin*extrapolate
							  + epsilonSquare
							),
							SMALL
						),
						lim
					)
				);
			}
			else
			{
				// Limiter remains unchanged
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
