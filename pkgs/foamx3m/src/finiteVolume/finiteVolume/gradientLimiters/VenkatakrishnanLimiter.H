
#ifndef VenkatakrishnanLimiter_H
#define VenkatakrishnanLimiter_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class VenkatakrishnanLimiter
{
	// Private data

		//- Limiter value
		scalar k_;


public:

	// Constructor

		//- Construct null
		VenkatakrishnanLimiter()
		:
			k_(5)
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
			scalar epsilonSquare = pow3(k_)*cellVolume;

			if (mag(extrapolate) < SMALL)
			{
				// Limiter remains unchanged
				return;
			}
			else if (extrapolate > VSMALL)
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
