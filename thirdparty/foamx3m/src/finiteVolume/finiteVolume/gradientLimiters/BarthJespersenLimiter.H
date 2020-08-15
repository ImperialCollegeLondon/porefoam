
#ifndef BarthJespersenLimiter_H
#define BarthJespersenLimiter_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class BarthJespersenLimiter
{
public:

	// Constructors

		//- Construct null
		BarthJespersenLimiter()
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
			if (mag(extrapolate) < SMALL)
			{
				return;
			}

			if (extrapolate - deltaOneMax > SMALL)
			{
				lim = min(lim, deltaOneMax/extrapolate);
			}
			else if (extrapolate - deltaOneMin < -SMALL)
			{
				lim = min(lim, deltaOneMin/extrapolate);
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
