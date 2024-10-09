
#include <complex>
#include <cmath>
#include <time.h>
#include "FluctuationSpectrum.h"
#include "State.h"
#include <algorithm> // For std::sort

using Complex = std::complex<double>;
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Energy of a single vertex
 Energy of a link (when connected vertices has inclusions)
 Energy of the whole system
 */
FluctuationSpectrum::FluctuationSpectrum(State* pState,int nx,int ny){
    m_pState = pState;
    m_Nx = nx;
    m_Ny = ny;
    GenerateVectorOrder();
    m_SpectrumSize=m_VectorOrder.size();


    GenerateZeroAndNonZeroVectorOrder();

    


}
FluctuationSpectrum::~FluctuationSpectrum() {
    
}

void FluctuationSpectrum::OpenOutputStreams() {

    std::string hqfilename=m_pState->GetAnalysisVariables()->GetFolderName() + '/' + m_pState->GetAnalysisVariables()->GetNameHQVectorFile();
    std::string qfilename=m_pState->GetAnalysisVariables()->GetFolderName() + '/' + m_pState->GetAnalysisVariables()->GetNameQVectorFile();
    m_QVector.open(qfilename);
    m_HQVector.open(hqfilename);
}

void FluctuationSpectrum::CloseOutputStreams() {


    if (m_QVector.is_open()) {
        m_QVector.flush(); 
        m_QVector.close();
     }

    if (m_HQVector.is_open()) {
        m_HQVector.flush(); 
        m_HQVector.close();
     }

}

void FluctuationSpectrum::GenerateVectorOrder(){
    std::vector<int> input(2,0);
    for(int i = 0; i < m_Ny; i++){
        for(int j = 0; j <= i; j++){
            input[1]=j;
            input[0]=i;
            m_VectorOrder.push_back(input);
        }
    }

    auto distance = [](const std::vector<int>& v) {
        return std::sqrt(v[0] * v[0] + v[1] * v[1]);
    };

    // Sort m_VectorOrder based on the calculated distance
    std::sort(m_VectorOrder.begin(), m_VectorOrder.end(),
              [&distance](const std::vector<int>& a, const std::vector<int>& b) {
                  return distance(a) < distance(b);
              });
}

void FluctuationSpectrum::GenerateZeroAndNonZeroVectorOrder(){

    for (size_t i = 0; i < m_VectorOrder.size(); ++i) {
        std::vector<std::vector<int>> MatrixOrder(4, std::vector<int>(2, 0));


        if(m_VectorOrder[i][0]==m_VectorOrder[i][1]){
            MatrixOrder[0]=m_VectorOrder[i];
            MatrixOrder[1] = {m_VectorOrder[i][0], -m_VectorOrder[i][1]};
            m_MatrixOrder.push_back(MatrixOrder);}
        else if (m_VectorOrder[i][0]==0 || m_VectorOrder[i][1]==0){
            MatrixOrder[0]=m_VectorOrder[i];
            MatrixOrder[1] = {m_VectorOrder[i][1], m_VectorOrder[i][0]};
            m_MatrixOrder.push_back(MatrixOrder);
        }
        
        else{
            MatrixOrder[0]=m_VectorOrder[i];
            MatrixOrder[1] = {m_VectorOrder[i][1], m_VectorOrder[i][0]};
            MatrixOrder[2] = {m_VectorOrder[i][0], -m_VectorOrder[i][1]};
            MatrixOrder[3] = {m_VectorOrder[i][1], -m_VectorOrder[i][0]};
            m_MatrixOrder.push_back(MatrixOrder);
        }

    }
}

double FluctuationSpectrum::FourierTransform(std::vector<std::vector<double>> qvector){

    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    std::vector<Complex> sum_p(qvector.size(),Complex(0,0));

    double x;
    double y;
    double z;
    double result=0;

    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
            x=(*it)->GetVXPos();
            y=(*it)->GetVYPos();
            z=(*it)->GetVZPos();
            for (size_t i = 0; i < qvector.size(); ++i) {
            sum_p[i]=sum_p[i]+Complex(z-m_AverageHeight,0)*std::exp(Complex(0, -qvector[i][0]*x-qvector[i][1]*y));
            }
    }

    for (size_t i = 0; i < qvector.size(); ++i) {
        result+=(sum_p[i]*std::conj(sum_p[i])).real();
        //std::cout<<"Fourier result: "<<result<<std::endl;
    }

    result=result/static_cast<double>(qvector.size());

    return result;
}

void FluctuationSpectrum::AverageHeight(){
    
    double average_height=0;
    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetActiveV();
    double number_vertex=0;

    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin() ; it != pAllVertices.end(); ++it) {
            average_height+=(*it)->GetVZPos();
            number_vertex+=1;
    }

    average_height/=number_vertex;
    m_AverageHeight=average_height;
}

void FluctuationSpectrum::CalculateSpectrum(){
    //Generate two vectors, one that gives the modulus of q and another one that gives the value of hq
    //First, get the size of the box. Lx*Ly

    double Lx=(*(m_pState->GetMesh()->GetBox()))(0);
    double Ly=(*(m_pState->GetMesh()->GetBox()))(1);
    AverageHeight();

    std::vector<double> qvector(m_VectorOrder.size(),0);
    std::vector<double> hqvector(m_VectorOrder.size(),0);

    for (size_t i = 0; i < m_MatrixOrder.size(); ++i) {

        std::vector<double> q={m_VectorOrder[i][0]*2*PI/Lx,m_VectorOrder[i][1]*2*PI/Ly};
        qvector[i]=std::sqrt(std::pow(q[0], 2) + std::pow(q[1], 2));

        if((m_VectorOrder[i][0]==m_VectorOrder[i][1]) || (m_VectorOrder[i][0]==0 || m_VectorOrder[i][1]==0)){
            std::vector<std::vector<double>> qmatrix(2,std::vector<double>(2,0));
            for (size_t j = 0; j < 2; ++j) {
                qmatrix[j][0]=static_cast<double>(m_MatrixOrder[i][j][0])*2*PI/Lx;
                qmatrix[j][1]=static_cast<double>(m_MatrixOrder[i][j][1])*2*PI/Ly;
            }

            hqvector[i]=FourierTransform(qmatrix);
        }   
        else{
            std::vector<std::vector<double>> qmatrix(4,std::vector<double>(2,0));
            for (size_t j = 0; j < 4; ++j) {
                qmatrix[j][0]=static_cast<double>(m_MatrixOrder[i][j][0])*2*PI/Lx;
                qmatrix[j][1]=static_cast<double>(m_MatrixOrder[i][j][1])*2*PI/Ly;
            }

            hqvector[i]=FourierTransform(qmatrix);  
            }

    }

    
    // Write the qvector and hqvector to files
    

    if (m_QVector.is_open()) {
        for (size_t i = 0; i < qvector.size(); ++i) {
            m_QVector << qvector[i] << " ";
        }
        m_QVector << "\n";
    } else {
        std::cerr << "Unable to open qvector.txt for writing." << std::endl;
    }

    if (m_HQVector.is_open()) {
        for (size_t i = 0; i < hqvector.size(); ++i) {
            m_HQVector << hqvector[i] << " ";
        }
        m_HQVector << "\n";
    } else {
        std::cerr << "Unable to open hqvector.txt for writing." << std::endl;
    }

    
}


