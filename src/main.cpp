#include "CDomain.h" //ds grid
#include "Timer.h"   //ds time measurement
#include <stdlib.h>  //ds atoi
#include <math.h>    //ds sqrt, etc.
#include <iostream>  //ds cout

int main( int argc, char** argv )
{
    //ds check simple input arguments - CAUTION: the implementation expects real numbers, the simulation will be corrupted if invalid values are entered
    if( 4 != argc )
    {
        //ds inform
        std::cout << "usage: diffusion_pse [Number of grid points 1D] [Number of time steps] [Performance mode: 1(yes)/0(no)]" << std::endl;
        return 0;
    }

    //ds start timing
    Timer tmTimer; tmTimer.start( );

    //ds get parameters
    const double dDiffusionCoefficient( 1.0 );
    const std::pair< double, double > prBoundaries( 0.0, 1.0 );
    const double dBoundarySize = fabs( prBoundaries.second - prBoundaries.first );
    const unsigned int uNumberOfGridPoints1D( atoi( argv[1] ) );
    const unsigned int uNumberOfParticles( uNumberOfGridPoints1D*uNumberOfGridPoints1D );
    const double dGridPointSpacing( dBoundarySize/( uNumberOfGridPoints1D - 1 ) );
    const unsigned int uNumberOfTimeSteps( atoi( argv[2] ) );
    const double dTimeStepSize( 0.5*dGridPointSpacing*dGridPointSpacing/dDiffusionCoefficient );

    //ds user information
    std::cout << "\n---------------------------- DIFFUSION PSE SETUP ----------------------------" << std::endl;
    std::cout << "   Diffusion Coefficient: "  << dDiffusionCoefficient << std::endl;
    std::cout << "           Boundary (2D): [" << prBoundaries.first << ", " << prBoundaries.second << "]" << std::endl;
    std::cout << "Number of Grid Points 1D: "  << uNumberOfGridPoints1D << std::endl;
    std::cout << "     Number of Particles: "  << uNumberOfParticles << std::endl;
    std::cout << "      Grid Point Spacing: "  << dGridPointSpacing << std::endl;
    std::cout << "    Number of Time Steps: "  << uNumberOfTimeSteps << std::endl;
    std::cout << "          Time Step Size: "  << dTimeStepSize << std::endl;
    std::cout << "------------------------------------------------------------------------------" << std::endl;

    //ds allocate domain (automatically creates initial density distribution)
    Diffusion::CDomain cDomain( dDiffusionCoefficient, prBoundaries, dBoundarySize, uNumberOfGridPoints1D, uNumberOfParticles, dGridPointSpacing, dTimeStepSize );

    //ds information
    std::cout << "               Status:  0% done - current time: 0";

    //ds get the mode mode
    const unsigned int uMode( atoi( argv[3] ) );

    //ds check for performance run (no streaming)
    if( 1 == uMode )
    {
        //ds start simulation
        for( unsigned int uCurrentTimeStep = 1; uCurrentTimeStep < uNumberOfTimeSteps+1; ++uCurrentTimeStep )
        {
            //ds calculate percentage done
            const double dPercentageDone( 100.0*uCurrentTimeStep/uNumberOfTimeSteps );

            //ds and time
            const double dCurrentTime( uCurrentTimeStep*dTimeStepSize );

            //ds get a formatted string -> 100% -> 3 digits
            char chBuffer[4];

            //ds fill the buffer
            std::snprintf( chBuffer, 4, "%3.0f", dPercentageDone );

            //ds print info
            std::cout << '\xd';
            std::cout << "               Status: " << chBuffer << "% done - current time: " << dCurrentTime;

            //ds update domain
            cDomain.updateHeatDistributionNumerical( );
        }
    }
    else
    {
        //ds start simulation
        for( unsigned int uCurrentTimeStep = 1; uCurrentTimeStep < uNumberOfTimeSteps+1; ++uCurrentTimeStep )
        {
            //ds calculate percentage done
            const double dPercentageDone( 100.0*uCurrentTimeStep/uNumberOfTimeSteps );

            //ds and time
            const double dCurrentTime( uCurrentTimeStep*dTimeStepSize );

            //ds get a formatted string -> 100% -> 3 digits
            char chBuffer[4];

            //ds fill the buffer
            std::snprintf( chBuffer, 4, "%3.0f", dPercentageDone );

            //ds print info
            std::cout << '\xd';
            std::cout << "               Status: " << chBuffer << "% done - current time: " << dCurrentTime;

            //ds update domain
            cDomain.updateHeatDistributionNumerical( );

            //ds streaming
            cDomain.saveHeatGridToStream( );
            cDomain.saveNormsToStream( dCurrentTime );
        }

        //ds save the streams to a file
        cDomain.writeHeatGridToFile( "bin/simulation.txt", uNumberOfTimeSteps );
        cDomain.writeNormsToFile( "bin/norms.txt", uNumberOfTimeSteps, dTimeStepSize );
    }

    //ds stop timing
    const double dDurationSeconds( tmTimer.stop( ) );

    //ds cause an output ostream
    std::cout << std::endl;
    std::cout << "     Computation time: " << dDurationSeconds << std::endl;
    std::cout << "------------------------------------------------------------------------------\n" << std::endl;

    return 0;
}
