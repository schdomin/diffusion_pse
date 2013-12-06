#include "CDomain.h"
#include <math.h>    //ds fabs, etc.
#include <iostream>



namespace Diffusion
{

//ds ctor/dtor
CDomain::CDomain( const double& p_dDiffusionCoefficient,
             const std::pair< double, double >& p_prBoundaries,
             const double& p_dBoundarySize,
             const unsigned int& p_uNumberOfGridPoints1D,
             const double& p_dGridPointSpacing,
             const double& p_dTimeStepSize ) : m_dDiffusionCoefficient( p_dDiffusionCoefficient ),
                                               m_prBoundaries( p_prBoundaries ),
                                               m_dBoundarySize( p_dBoundarySize ),
                                               m_uNumberOfGridPoints1D( p_uNumberOfGridPoints1D ),
                                               m_dGridPointSpacing( p_dGridPointSpacing ),
                                               m_dTimeStepSize( p_dTimeStepSize ),
                                               m_dDiffusionFactor( p_dTimeStepSize*p_dDiffusionCoefficient/( 2*( p_dGridPointSpacing*p_dGridPointSpacing ) ) ),
                                               m_dEpsilon( 2*p_dGridPointSpacing ),
                                               m_dVolP( p_dGridPointSpacing*p_dGridPointSpacing ),
                                               m_strLogHeatDistribution( "" ),
                                               m_strLogNorms( "" )

{
    //ds allocate memory for the data structure
    m_gridHeat = new double*[m_uNumberOfGridPoints1D];

    //ds for each element
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        m_gridHeat[u] = new double[m_uNumberOfGridPoints1D];
    }

    //ds initialize grid
    setInitialHeatDistribution( );
};

CDomain::~CDomain( )
{
    //ds deallocate heat structure
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        delete[] m_gridHeat[u];
    }

    delete[] m_gridHeat;
};

//ds accessors
void CDomain::updateHeatDistributionNumerical( )
{
    //ds for all grid points
    for( unsigned int uk = 0; uk < m_uNumberOfGridPoints1D; ++uk )
    {
        for( unsigned int vk = 0; vk < m_uNumberOfGridPoints1D; ++vk )
        {
            //ds get 2d vector of current coordinates
            const double dXk[2] = { uk*m_dGridPointSpacing, vk*m_dGridPointSpacing };

            //ds inner sum of formula
            double dInnerSum( 0.0 );

            //ds loop 20x20
            for( unsigned int u = 0; u < 20; ++u )
            {
                for( unsigned int v = 0; v < 20; ++v )
                {

                }
            }

            //ds add final part of formula
            m_gridHeat[uk][vk] = m_dTimeStepSize*m_dDiffusionCoefficient/( m_dEpsilon*m_dEpsilon )*m_dVolP*dInnerSum + m_gridHeat[uk][vk];
        }
    }
}

void CDomain::updateHeatDistributionAnalytical( const double& p_dCurrentTime )
{
    //ds for all grid points
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            //ds get the exact heat at this point
            m_gridHeat[u][v] = getHeatAnalytical( u*m_dGridPointSpacing, v*m_dGridPointSpacing, p_dCurrentTime );
        }
    }
}

void CDomain::saveHeatGridToStream( )
{
    //ds buffer for snprintf
    char chBuffer[16];

    //ds add each element
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            //ds get the integrals stream
            std::snprintf( chBuffer, 16, "%f", m_gridHeat[u][v] );

            //ds add buffer and space to string
            m_strLogHeatDistribution += chBuffer;
            m_strLogHeatDistribution += " ";
        }

        //ds add new line for new row
        m_strLogHeatDistribution += '\n';
    }
}

void CDomain::saveNormsToStream( const double& p_dCurrentTime )
{
    //ds get total heat
    double dTotalHeat( 0.0 );

    //ds norms
    double dLInfinity( 0.0 );
    double dL1( 0.0 );
    double dL2( 0.0 );

    //ds for all grid points
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            //ds add the current heat value
            dTotalHeat += m_gridHeat[u][v];

            //ds calculate the error (numerical - analytical) - the absolute value is already here taken since it does not change L
            const double dErrorAbsolute( fabs( m_gridHeat[u][v] - getHeatAnalytical( u*m_dGridPointSpacing, v*m_dGridPointSpacing, p_dCurrentTime ) ) );

            //ds check for LInf
            if( dLInfinity < dErrorAbsolute )
            {
                dLInfinity = dErrorAbsolute;
            }

            //ds L1, L2
            dL1 += dErrorAbsolute;
            dL2 += dErrorAbsolute*dErrorAbsolute;
        }
    }

    //ds scale L1, L2
    dL1 /= m_uNumberOfGridPoints1D;
    dL2 /= m_uNumberOfGridPoints1D;
    dL2 = sqrt( dL2 );

    //ds buffer for snprintf
    char chBuffer[64];

    //ds get the norms stream: E Linf L1 L2
    std::snprintf( chBuffer, 64, "%f %f %f %f", dTotalHeat, dLInfinity, dL1, dL2 );

    //ds append the buffer to our string
    m_strLogNorms += chBuffer;
    m_strLogNorms += "\n";
}

void CDomain::writeHeatGridToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps ) const
{
    //ds ofstream object
    std::ofstream ofsFile;

    //ds open the file for writing
    ofsFile.open( p_strFilename.c_str( ), std::ofstream::out );

    //ds if it worked
    if( ofsFile.is_open( ) )
    {
        //ds first dump setup information number of points and time steps
        ofsFile << m_uNumberOfGridPoints1D << " " << p_uNumberOfTimeSteps << "\n" << m_strLogHeatDistribution;
    }

    //ds close the file
    ofsFile.close( );
}

void CDomain::writeNormsToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps, const double& p_dTimeStepSize ) const
{
    //ds ofstream object
    std::ofstream ofsFile;

    //ds open the file for writing
    ofsFile.open( p_strFilename.c_str( ), std::ofstream::out );

    //ds if it worked
    if( ofsFile.is_open( ) )
    {
        //ds dump information to file
        ofsFile << p_uNumberOfTimeSteps << " " << p_dTimeStepSize << "\n" << m_strLogNorms;
    }

    //ds close the file
    ofsFile.close( );
}

//ds helpers
void CDomain::setInitialHeatDistribution( )
{
    //ds loop over all indexi
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        for( unsigned int v = 0; v < m_uNumberOfGridPoints1D; ++v )
        {
            //ds set initial heat value
            m_gridHeat[u][v] = sin( u*m_dGridPointSpacing*2*M_PI )*sin( v*m_dGridPointSpacing*2*M_PI );
        }
    }

    //ds set boundary values
    for( unsigned int u = 0; u < m_uNumberOfGridPoints1D; ++u )
    {
        m_gridHeat[u][0]                         = 0.0;
        m_gridHeat[u][m_uNumberOfGridPoints1D-1] = 0.0;
        m_gridHeat[0][u]                         = 0.0;
        m_gridHeat[m_uNumberOfGridPoints1D-1][u] = 0.0;
    }
}

double CDomain::getHeatAnalytical( const double& p_dX, const double& p_dY, const double& p_dT ) const
{
    //ds formula
    return sin( p_dX*2*M_PI )*sin( p_dY*2*M_PI )*exp( -8*m_dDiffusionCoefficient*M_PI*M_PI*p_dT );
}

double CDomain::getKernel( const double p_dXp[2], const double p_dXk[2] ) const
{
    //ds compute distance (using pow for readability)
    const double dDistance = sqrt( ( p_dXp[1] - p_dXk[1] )*( p_dXp[1] - p_dXk[1] ) + ( p_dXp[2] - p_dXk[2] )*( p_dXp[2] - p_dXk[2] ) );

    //ds if we are out of the cutoff
    if( 5*m_dEpsilon < dDistance )
    {
        return 0.0;
    }
    else
    {
        //return kernel function
        return 16.0/( M_PI*M_PI*( pow( dDistance, 8 ) + 1.0 ) );
    }
}

} //namespace Diffusion
