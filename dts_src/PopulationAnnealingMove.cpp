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

/*
===============================================================================================================
 Last update Aug 2023 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class changes the box in x and y direction to minimize the energy. It is called in input file by Frame_Tension  = on 0 2, where the first flag should be on to start this command and the second argument is frame tension and the third argument is the priodic of the operating.

The way it works is based on the changing the box in x and y direction by a small size of dr by calling "ChangeBoxSize" function.
=================================================================================================================
*/
PopulationAnnealingMove::PopulationAnnealingMove(State *pState, std::string period_file, std::string temperature_file, std::string input_files, int n_processors) :
        m_pState(pState),
        m_pPeriodFile(period_file),
        m_pTemperatureFile(temperature_file),
        m_pInputFiles(input_files),
        m_pInputSize(n_processors)
{
}
//-- convert the type string into the direction vector

PopulationAnnealingMove::~PopulationAnnealingMove() {

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
    if (step % m_pPeriod != 0)
        return false;

    MPI_Barrier(MPI_COMM_WORLD);

    return true;
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
        if (betavec[i] < betavec[i - 1]) {
            std::cerr << "Error: Vector is not sorted in ascending order." << std::endl;
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