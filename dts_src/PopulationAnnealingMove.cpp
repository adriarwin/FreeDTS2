#include <stdio.h>
#include "PopulationAnnealingMove.h"
#include "Nfunction.h"
#include "vertex.h"
#include "State.h"
#include "Voxelization.h"
#ifdef MPI_DETECTED
#include <mpi.h>
#include <sstream>
#include <iostream>
#include <filesystem>

/*
===============================================================================================================
 Last update Aug 2023 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class changes the box in x and y direction to minimize the energy. It is called in input file by Frame_Tension  = on 0 2, where the first flag should be on to start this command and the second argument is frame tension and the third argument is the priodic of the operating.

The way it works is based on the changing the box in x and y direction by a small size of dr by calling "ChangeBoxSize" function.
=================================================================================================================
*/
PopulationAnnealingMove::PopulationAnnealingMove(State *pState, std::string period_file, std::string temperature_file, std::string topology_files, int n_processors) :
        m_pState(pState),
        m_pPeriodFile(period_file),
        m_pTemperatureFile(temperature_file),
        m_pTopologyFile(topology_files),
        m_pInputSize(n_processors)
{
}
//-- convert the type string into the direction vector

PopulationAnnealingMove::~PopulationAnnealingMove() {

}

bool PopulationAnnealingMove::PopulationAnnealingMoveOn(){
    return true;
}

void PopulationAnnealingMove::Initialize() {

#ifdef MPI_DETECTED
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //std::cout<<"Rank: "<<rank<<std::endl;
    m_pRank = rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //std::cout<<"Size: "<<size<<std::endl;
    m_pSize = size;
    m_pCounter = 0;

    if(m_pRank==0){
    std::cout<<"---> the algorithm for Populated Annealing involves applying this: "<< GetBaseDefaultReadName()<<" \n";}
    std::vector<double> localBetaVec;
    if (rank == 0) {
        localBetaVec = ReadTemperatures();
    }

    // Broadcast the size of the vector to all ranks
    int vecSize = (rank == 0) ? localBetaVec.size() : 0;
    MPI_Bcast(&vecSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resize the local vector based on the size received
    std::vector<double> globalBetaVec(vecSize);
    if (rank == 0) {
        // Ensure localBetaVec.data() is the correct pointer for broadcasting
        MPI_Bcast(localBetaVec.data(), vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Bcast(globalBetaVec.data(), vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Now each rank has a copy of the vector
    m_pBetaVector = globalBetaVec;
    if (rank==0){m_pBetaVector=localBetaVec;}

    std::vector<int> localPeriodVec;
    if (rank == 0) {
        localPeriodVec = ReadPeriods();
    }
    
    // Broadcast the size of the vector to all ranks
    int vecSizep = (rank == 0) ? localPeriodVec.size() : 0;
    MPI_Bcast(&vecSizep, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resize the local vector based on the size received
    std::vector<int> globalPeriodVec(vecSizep);
    if (rank == 0) {
        // Ensure localBetaVec.data() is the correct pointer for broadcasting
        MPI_Bcast(localPeriodVec.data(), vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Bcast(globalPeriodVec.data(), vecSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Now each rank has a copy of the vector
    m_pPeriodVector = globalPeriodVec;
    if (rank==0){m_pPeriodVector=localPeriodVec;}

    //Checking that both vectors have the same size
    if (m_pBetaVector.size() != m_pPeriodVector.size()) {
        std::ostringstream errorMsg;
        errorMsg << "Error on rank " << m_pRank << ": Vector sizes do not match. "
                 << "m_pBetaVector size: " << m_pBetaVector.size() << ", "
                 << "m_pPeriodVector size: " << m_pPeriodVector.size();
        std::cerr << errorMsg.str() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);  // Abort all MPI processes
    }

    // Read input files
    if (rank == 0) {
        std::ifstream file(m_pTopologyFile);
        std::string line;
        while (std::getline(file, line)) {
            m_pTopologyFilesVector.push_back(line);
        }

        // Check if number of files matches number of processors
        if (m_pTopologyFilesVector.size() != m_pSize) {
            std::cerr << "Error: Number of input files (" << m_pTopologyFilesVector.size() 
                      << ") does not match number of processors (" << size << ")" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Check if all files exist
        for (const auto& file : m_pTopologyFilesVector) {
            if (!std::filesystem::exists(file)) {
                std::cerr << "Error: File " << file << " does not exist" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    // Broadcast the number of files
    int numFiles = m_pTopologyFilesVector.size();
    MPI_Bcast(&numFiles, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate space for each rank's file
    std::string myFile;

    if (rank == 0) {
        // Send each filename to the corresponding rank
        for (int i = 1; i < size; ++i) {
            MPI_Send(m_pTopologyFilesVector[i].c_str(), m_pTopologyFilesVector[i].size() + 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }
        myFile = m_pTopologyFilesVector[0];
    } else {
        // Receive the filename
        MPI_Status status;
        MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
        int count;
        MPI_Get_count(&status, MPI_CHAR, &count);
        char* buffer = new char[count];
        MPI_Recv(buffer, count, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        myFile = std::string(buffer);
        delete[] buffer;
    }

    std::cout << "Rank " << rank << " will use file: " << myFile << std::endl;

    // Store the filename for later use
    m_pTopologyFileString = myFile;

    std::string newFolderName = "_" + Nfunction::Int_to_String(m_pRank);

    m_pState->GetNonbinaryTrajectory()->SetFolderName(m_pState->GetNonbinaryTrajectory()->GetOriginalFolderName() +"_" + Nfunction::Int_to_String(m_pRank));
    m_pState->GetTimeSeriesDataOutput()->SetCustomFileName(m_pState ->GetRunTag() + "_" +Nfunction::Int_to_String(m_pRank)+TimeSeriDataExt);
    m_pState->GetRestart()->SetUniqueRestartFileName("_" + Nfunction::Int_to_String(m_pRank));
    m_pState->GetVisualization()->ChangeFolderName(newFolderName);


    int totalSteps = 0;
    for (const auto& period : m_pPeriodVector) {
        totalSteps += period;
    }

    // Set initial and final steps
    int ini = 0;
    int fi = totalSteps-1;
    m_pPeriod=m_pPeriodVector[0];


    m_pState->GetSimulation()->SetBeta(m_pBetaVector[m_pCounter], 0);
    m_pState->GetSimulation()->UpdateInitialStep(ini);
    m_pState->GetSimulation()->UpdateFinalStep(fi);

    if (m_pRank==0){
        m_pEnergyVector.resize(m_pSize);
    }

#endif
}

bool PopulationAnnealingMove::EvolveOneStep(int step) {
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

    MPI_Barrier(MPI_COMM_WORLD);

    double local_energy=m_pState->GetEnergyCalculator()->GetEnergy();

    MPI_Gather(&local_energy, 1, MPI_DOUBLE, 
               m_pEnergyVector.data(), 1, MPI_DOUBLE, 
               0, MPI_COMM_WORLD);

    
    if (m_pRank == 0) {
        // Now energy_vector contains energies from all ranks
        // You can process or print the energies here
        std::cout << "Energies from all ranks:" << std::endl;
        for (int i = 0; i < m_pSize; ++i) {
            std::cout << "Rank " << i << ": " << m_pEnergyVector[i] << std::endl;
        }
    }

    
    //Update the counter
    m_pCounter++;
    m_pPeriod+=m_pPeriodVector[m_pCounter];

    //Now it comes the exchange part.


    m_pState->GetSimulation()->SetBeta(m_pBetaVector[m_pCounter], 0);

    //std::cout << "Rank " << m_pRank << ": Counter updated to: " << m_pCounter << std::endl;
    //std::cout << "Rank " << m_pRank << ": New period set to: " << m_pPeriod << std::endl;

    

    return true;
}


std::vector<int> PopulationAnnealingMove::ReadPeriods() {
    std::vector<int> periods;
    std::ifstream file(m_pPeriodFile);
    int period;
    while (file >> period) {
        periods.push_back(period);
    }
    return periods;
}

std::vector<double> PopulationAnnealingMove::ReadTemperatures() {

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

bool PopulationAnnealingMove::is_line_empty(const std::string& line) {
    // Check if the line is empty or contains only whitespace
    return line.empty() || std::all_of(line.begin(), line.end(), [](unsigned char c) { return std::isspace(c); });
}


// Function to check if a string is a valid integer
bool PopulationAnnealingMove::isInteger(const std::string& s) {
    return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

// Function to split a string by spaces and convert to integers
std::vector<int> PopulationAnnealingMove::splitStringToIntVector(const std::string& line) {
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

// Function to read the last line of a file and parse it into a vector
std::vector<int> PopulationAnnealingMove::ReadLastLineOutput(const std::string& output, int N) {
    std::ifstream file(output);
    std::string lastLine;
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << output << std::endl;
        // Return a default vector if the file cannot be opened
        std::vector<int> defaultVec;
        for (int i = 0; i < N; ++i) {
            defaultVec.push_back(i);
        }
        return defaultVec;
    }
    std::cout << "Here" << std::endl;
    // Read the file to find the last non-empty line
    while (std::getline(file, line)) {
        if (!line.empty()) {
            lastLine = line;
        }
    }
    std::cout << "Here 2" << std::endl;

    file.close();

    // Try to convert the last line to a vector of integers
    std::vector<int> result = splitStringToIntVector(lastLine);
    if (!result.empty() && result.size() == N) {
        return result;
    }

    // If the last line is not a valid integer vector, return a default vector
    std::vector<int> defaultVec;
    for (int i = 0; i < N; ++i) {
        defaultVec.push_back(i);
    }
    return defaultVec;
}

std::string PopulationAnnealingMove::CurrentState() {

    std::string state = GetBaseDefaultReadName() + " = " + this->GetDerivedDefaultReadName();
    //What is this used for?
    //state = state +" "+ Nfunction::D2S(m_Period)+" "+Nfunction::D2S(m_SigmaP)+" "+  m_Type;
    return state;
}

#endif