

if (adjustTimeStep)
{


    scalar maxDeltaTFactRe = maxCo/(CoNum + SMALL);
    
    scalar maxDeltaTFactCa=maxAlphaCo/(alphaCoNum + SMALL);

    scalar maxDeltaTFact = min(maxDeltaTFactCa, maxDeltaTFactRe);
        
    
    
    
    
    if( maxDeltaTFact>1.2 && (CoNum < CoNumOld*deltaTFact || maxDeltaTFact > 1.5) )
    {
        deltaTFact = min(min(maxDeltaTFact, 1. + 0.1*(maxDeltaTFact-1.)), 1.3);
        runTime.setDeltaT
        (
            min
            (
                deltaTFact*runTime.deltaT().value(),
                maxDeltaT
            )
        );
    }
    else if ( maxDeltaTFact< 0.95 )
    {
        deltaTFact = min(maxDeltaTFact, 0.81);

        runTime.setDeltaT
        (
            min
            (
                deltaTFact*runTime.deltaT().value(),
                maxDeltaT
            )
        );
    }
    else
    {
        deltaTFact = 1.;
    }
        
    CoNumOld=CoNum;   


    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}

// ************************************************************************* //
