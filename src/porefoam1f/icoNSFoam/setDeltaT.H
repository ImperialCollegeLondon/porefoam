
//Description
//    Reset the timestep to maintain a constant maximum courant Number.
//    Reduction of time-step is immediate, but increase is damped to avoid
//    unstable oscillations.


if (adjustTimeStep)
{
	if(!(CoNum<10000.)){Info<<"Error CoNum="<<CoNum<<endl;  exit(-1);}

    deltaTFact =  min(maxCo/(CoNum + SMALL),1.2);



    if(( (deltaTFact-1.)>0.05 )||( (1.-deltaTFact)>0.02 ))
    {
        runTime.setDeltaT
        (
            min
            (
                deltaTFact*runTime.deltaT().value(),
                maxDeltaT
            )
        );
    }
    
    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}

// ************************************************************************* //
