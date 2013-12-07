#ifndef CDOMAIN_H_
#define CDOMAIN_H_

#include <fstream> //ds streaming to file for visualization



namespace Diffusion
{

class CDomain
{

//ds ctor/dtor
public:

    CDomain( const double& p_dDiffusionCoefficient,
             const std::pair< double, double >& p_prBoundaries,
             const double& p_dBoundarySize,
             const unsigned int& p_uNumberOfGridPoints1D,
             const unsigned int& p_uNumberOfParticles,
             const double& p_dGridPointSpacing,
             const double& p_dTimeStepSize );

    ~CDomain( );

//ds attributes
private:

    //ds heat structure - we use a standard structure since the constant size is known at allocation and we benefit from access speed
    double** m_gridHeat;

    //domain properties
    const double m_dDiffusionCoefficient;
    const std::pair< double, double > m_prBoundaries;
    const double m_dBoundarySize;
    const double m_dGridPointSpacing;
    const unsigned int m_uNumberOfGridPoints1D;
    const unsigned int m_uNumberOfParticles;
    const double m_dTimeStepSize;

    //ds 2d-kernel
    const double m_dEpsilon;
    const double m_dVolP;
    const double m_PSEFactor;

    //ds stream for offline data - needed for the movie and graphs
    std::string m_strLogHeatDistribution;
    std::string m_strLogNorms;

//ds accessors
public:

    void updateHeatDistributionNumerical( );
    void updateHeatDistributionAnalytical( const double& p_dCurrentTime );
    void saveHeatGridToStream( );
    void saveNormsToStream( const double& p_dCurrentTime );
    void writeHeatGridToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps ) const;
    void writeNormsToFile( const std::string& p_strFilename, const unsigned int& p_uNumberOfTimeSteps, const double& p_dTimeStepSize ) const;

//ds helpers
private:

    void setInitialHeatDistribution( );
    double getHeatAnalytical( const double& p_dX, const double& p_dY, const double& p_dT ) const;

    //kernel
    double getKernel( const double p_dXp[2], const double p_dXk[2] ) const;

}; //class CDomain

} //namespace Diffusion



#endif //CDOMAIN_H_
