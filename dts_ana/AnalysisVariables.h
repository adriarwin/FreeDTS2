#ifndef ANALYSISVARIABLES_H
#define ANALYSISVARIABLES_H

#include <vector>
#include <string>
#include "SimDef.h"
#include "Nfunction.h"

class AnalysisVariables {
public:

    //Folder variables
    inline static std::string GetInputFolderName() {return "FolderName";}
    void SetFolderName(std::string input) {m_pFolderName=input;}
    bool OpenFolder() {return  Nfunction::OpenFolder(m_pFolderName);}
    std::string GetFolderName() {return m_pFolderName;}

    inline static std::string GetInputTopologyName() {return "Topology";}
    void SetTopology(std::string topology){m_pTopology=topology;}
    std::string GetTopology(){return m_pTopology;}


    inline static std::string GetFluctuationSpectrumName() {return "FluctuationSpectrum";}
    void SetFluctationSpectrumActive() {m_pFluctuationSpectrumActive=true;}
    bool GetFluctuationSpectrumActive() {return m_pFluctuationSpectrumActive;}
    std::string GetNameQVectorFile()    {return m_NameQVectorFile;}
    std::string GetNameHQVectorFile()    {return m_NameHQVectorFile;}
    std::string GetNameHPQVectorFile()    {return m_NameHPQVectorFile;}
    std::string GetNamePQVectorFile()    {return m_NamePQVectorFile;}


    inline static std::string GetAreaName() {return "Area";}
    void SetAreaCalculationActive() {m_pAreaCalculationActive=true;}
    bool GetAreaCalculationActive() {return m_pAreaCalculationActive;}

    inline static std::string GetVolumeName() {return "Volume";}
    void SetVolumeCalculationActive() {m_pVolumeCalculationActive=true;}
    bool GetVolumeCalculationActive() {return m_pVolumeCalculationActive;}

    inline static std::string GetEnergyName() {return "Energy";}
    void SetEnergyCalculationActive() {m_pEnergyCalculationActive=true;}
    bool GetEnergyCalculationActive() {return m_pEnergyCalculationActive;}


    inline static std::string GetProjectedAreaName() {return "ProjectedArea";}
    void SetProjectedAreaCalculationActive() {m_pProjectedAreaCalculationActive=true;}
    bool GetProjectedAreaCalculationActive() {return m_pProjectedAreaCalculationActive;}

    inline static std::string GetMeanCurvatureName() {return "MeanCurvature";}
    void SetMeanCurvatureCalculationActive() {m_pMeanCurvatureCalculationActive=true;}
    bool GetMeanCurvatureCalculationActive() {return m_pMeanCurvatureCalculationActive;}

    inline static std::string GetGaussianCurvatureName() {return "GaussianCurvature";}
    void SetGaussianCurvatureCalculationActive() {m_pGaussianCurvatureCalculationActive=true;}
    bool GetGaussianCurvatureCalculationActive() {return m_pGaussianCurvatureCalculationActive;}

    inline static std::string GetThicknessName() {return "Thickness";}
    void SetThicknessCalculationActive() {m_pThicknessCalculationActive=true;}
    bool GetThicknessCalculationActive() {return m_pThicknessCalculationActive;}

    inline static std::string GetInclusionName() {return "Inclusion";}
    void SetInclusionCalculationActive() {m_pInclusionCalculationActive=true;}
    void SetInclusionCalculationOff() {m_pInclusionCalculationActive=false;}
    bool GetInclusionCalculationActive() {return m_pInclusionCalculationActive;}

    inline static std::string GetInclusionClusterName() {return "InclusionCluster";}
    void SetInclusionClusterCalculationActive() {m_pInclusionClusterCalculationActive=true;}
    bool GetInclusionClusterCalculationActive() {return m_pInclusionClusterCalculationActive;}
    std::string GetNameInclusionClusterFile() {return m_NameInclusionClusterFile;}


    inline static std::string GetVisualizationName() {return "Visualization";}
    void SetVisualizationActive() {m_pVisualizationActive=true;}
    bool GetVisualizationActive() {return m_pVisualizationActive;}

    std::string GetNameGeneralAnalysisFile() {return m_NameGeneralAnalysisFile;}
    

private:

    bool m_pFluctuationSpectrumActive=false;
    bool m_pAreaCalculationActive=false;
    bool m_pVolumeCalculationActive=false;
    bool m_pEnergyCalculationActive=false;
    bool m_pMeanCurvatureCalculationActive=false;
    bool m_pProjectedAreaCalculationActive=false;
    bool m_pGaussianCurvatureCalculationActive=false;
    bool m_pThicknessCalculationActive=false;
    bool m_pVisualizationActive=false;
    bool m_pInclusionCalculationActive=false;
    bool m_pInclusionClusterCalculationActive=false;

    std::string m_pTopology;

    std::string m_pFolderName;
    std::string m_NameQVectorFile="qvector.ana";
    std::string m_NameHQVectorFile="hqvector.ana";
    std::string m_NameHPQVectorFile="hpqvector.ana";
    std::string m_NamePQVectorFile="pqvector.ana";
    std::string m_NameGeneralAnalysisFile="dts.ana";
    std::string m_NameInclusionClusterFile="inclusioncluster.ana";
};

#endif // READTRAJTSIFILES_H