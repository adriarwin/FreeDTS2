#include <stdio.h>
#include "AnnealedImportanceSamplingMove.h"
#include "Nfunction.h"
#include "vertex.h"
#include "State.h"
#include "Voxelization.h"
#include <sstream>
#include <iostream>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <unordered_set>
/*
===============================================================================================================
 Last update Aug 2023 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class changes the box in x and y direction to minimize the energy. It is called in input file by Frame_Tension  = on 0 2, where the first flag should be on to start this command and the second argument is frame tension and the third argument is the priodic of the operating.

The way it works is based on the changing the box in x and y direction by a small size of dr by calling "ChangeBoxSize" function.
=================================================================================================================
*/
AnnealedImportanceSamplingMove::AnnealedImportanceSamplingMove(State *pState, std::string period_file, std::string temperature_file) :
        m_pState(pState),
        m_pPeriodFile(period_file),
        m_pTemperatureFile(temperature_file)
{
}
//-- convert the type string into the direction vector

AnnealedImportanceSamplingMove::~AnnealedImportanceSamplingMove() {

}

bool AnnealedImportanceSamplingMove::AnnealedImportanceSamplingMoveOn(){
    return true;
}

void AnnealedImportanceSamplingMove::Initialize() {


    m_pCounter = 0;


    std::cout<<"---> the algorithm for Annealed Importance Sampling involves applying this: "<< GetBaseDefaultReadName()<<" \n";
    std::vector<double> localBetaVec;
    localBetaVec = ReadTemperatures();
    m_pBetaVector=localBetaVec;
    

    std::vector<int> localPeriodVec;
    localPeriodVec = ReadPeriods();
    m_pPeriodVector=localPeriodVec;

    

    

    //Checking that both vectors have the same size
    if (m_pBetaVector.size() != m_pPeriodVector.size()) {
        std::ostringstream errorMsg;
        errorMsg <<  "Vector sizes do not match. "<<
                  "m_pBetaVector size: " << m_pBetaVector.size() << ", "
                 << "m_pPeriodVector size: " << m_pPeriodVector.size();
        std::cerr << errorMsg.str() << std::endl;
    }





    


    int ini,fi;
    double beta;


    int totalSteps = 0;
    for (const auto& period : m_pPeriodVector) {
        totalSteps += period;
    }

    // Set initial and final steps
    if (m_pRestart==false){
        ini = 0;
        fi = totalSteps-1;
        m_pPeriod=m_pPeriodVector[0];
        beta=m_pBetaVector[0];
    }
    else{
        ini = m_pState->GetSimulation()->GetInitialStep();
        fi = totalSteps-1;
        m_pPeriod = m_pPeriodVector[0];
        m_pCounter = 0;
        int i = 0;
        while (ini > m_pPeriod) {
            i++;
            m_pPeriod += m_pPeriodVector[i];
        }
        
        beta=m_pBetaVector[i];

        m_pCounter=i;
    }




    

    m_pState->GetSimulation()->SetBeta(beta, 0);
    m_pState->GetSimulation()->UpdateInitialStep(ini);
    m_pState->GetSimulation()->UpdateFinalStep(fi);



    if (m_pRestart==false){
        m_TimeSeriesFile.open(GetOutputFileName(), std::ios_base::out);
        m_TimeSeriesFile << "Weights: "<<std::endl;
            // Call the new function to write m_RankAtTempID to m_TimeSeriesFile
    }
    else{
        m_TimeSeriesFile.open(GetOutputFileName(), std::ios_base::app);
            // Call the new function to write m_RankAtTempID to m_TimeSeriesFile

    }


}
void AnnealedImportanceSamplingMove::WriteWeight(double weight) {

    if (!m_TimeSeriesFile.is_open()) {
        std::cerr << "Error: Time series file is not open!" << std::endl;
        return;
    }

    m_TimeSeriesFile << weight << std::endl;

}


bool AnnealedImportanceSamplingMove::EvolveOneStep(int step) {
    /**
     * @brief a call function to exchange temperature between replicas.
     *
     * This function attempts a temperature swap between different replicas.
     * The algorithm is a decentralized one, in the sense that each couple of
     * neighbouring replicas that attempt the exchange do not need of the other ones
     * to continue with the simulation.
     * 
     * @param step The current step in the simulation.
     * @return true if the box size was changed, false otherwise.
     */
    //---> if does not match the period, return false
    if (step <m_pPeriod)
        return false;


    m_pState->GetRestart()->UpdateRestartStateNoCondition(step, m_pState->GetVertexPositionUpdate()->GetDR(), m_pState->GetDynamicBox()->GetDR());
    m_pCounter++;
    m_pPeriod+=m_pPeriodVector[m_pCounter];

    double local_energy=m_pState->GetEnergyCalculator()->GetEnergy();

    double weight=(-1)*(m_pBetaVector[m_pCounter]-m_pBetaVector[m_pCounter-1])*(local_energy);

    
    WriteWeight(weight);
    

    //1. Broadcast new ranks to each rank (each rank will know what its next configuration will be)
    //2. Use restart?
    m_pState->GetSimulation()->SetBeta(m_pBetaVector[m_pCounter], 0);

    m_pState->GetEnergyCalculator()->UpdateTotalEnergy(m_pState->GetEnergyCalculator()->CalculateAllLocalEnergy());
        
   return true;

    

    
}

void AnnealedImportanceSamplingMove::SetRestart() {
    m_pRestart=true;
}




std::vector<int> AnnealedImportanceSamplingMove::ReadPeriods() {
    std::vector<int> periods;
    std::ifstream file(m_pPeriodFile);
    
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open period file: " + m_pPeriodFile);
    }

    int period;
    while (file >> period) {
        periods.push_back(period);
    }

    if (file.fail() && !file.eof()) {
        throw std::runtime_error("Non-integer value found in period file");
    }

    if (periods.empty()) {
        throw std::runtime_error("No valid periods found in file");
    }

    return periods;
}

std::vector<double> AnnealedImportanceSamplingMove::ReadTemperatures() {

    Nfunction f;
    std::string filename = m_pTemperatureFile;

    if (f.FileExist(filename) == false) {
        std::cout << "-----> Error: the temperature file name with the name " << filename << " does not exist" << std::endl;
        exit(0);
    }

    std::ifstream TemperatureFile(filename);

    // Check if the file was opened successfully
    if (!TemperatureFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(0);
    }

    std::vector<double> betavec;

    // File processing code goes here

    std::string line;
    while (std::getline(TemperatureFile, line)) {
        // Skip empty or whitespace-only lines
        if (is_line_empty(line)) {
            continue;
        }

        std::stringstream ss(line);  // Use a stringstream to parse the line
        double value;
        // Read the double value from the line
        if (ss >> value) {
            // Ensure there is no extra data in the line
            if (ss >> std::ws && !ss.eof()) {
                std::cerr << "Error: More than one value or extra data in line: \"" << line << "\"" << std::endl;
                exit(0);  // Exit with error code
            }

            betavec.push_back(value);  // Add each double to the data vector
        } else {
            std::cerr << "Error: Could not read a double value from line: \"" << line << "\"" << std::endl;
            exit(0);  // Exit with error code
        }
    }

    // Close the file
    TemperatureFile.close();

    for (size_t i = 1; i < betavec.size(); ++i) {
        if (betavec[i] > betavec[i - 1]) {
            std::cerr << "Error: Vector is not sorted in descending order." << std::endl;
            exit(0);  // Exit with error code
        }
    }

    for (size_t i = 0; i < betavec.size(); ++i) {
        betavec[i] = 1 / betavec[i];
    }

    return betavec;
}

bool AnnealedImportanceSamplingMove::is_line_empty(const std::string& line) {
    // Check if the line is empty or contains only whitespace
    return line.empty() || std::all_of(line.begin(), line.end(), [](unsigned char c) { return std::isspace(c); });
}


// Function to check if a string is a valid integer
bool AnnealedImportanceSamplingMove::isInteger(const std::string& s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

// Function to split a string by spaces and convert to integers
std::vector<int> AnnealedImportanceSamplingMove::splitStringToIntVector(const std::string& line) {
    std::vector<int> result;
    std::stringstream ss(line);
    std::string temp;

    while (ss >> temp) {
        if (isInteger(temp)) {
            result.push_back(std::stoi(temp));
        } else {
            // If any part is not an integer, return an empty vector
            return {};
        }
    }

    return result;
}



std::string AnnealedImportanceSamplingMove::CurrentState() {

    std::string state = GetBaseDefaultReadName() + " = " + this->GetDerivedDefaultReadName();
    //What is this used for?
    //state = state +" "+ Nfunction::D2S(m_Period)+" "+Nfunction::D2S(m_SigmaP)+" "+  m_Type;
    return state;
}

