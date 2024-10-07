
#include <complex>
#include <cmath>
#include <time.h>
#include "FluctuationSpectrum.h"
#include "State.h"

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
    std::cout<<"Spectrum size: "<<m_SpectrumSize<<std::endl;

    for (size_t i = 0; i < m_VectorOrder.size(); ++i) {
        for (size_t j = 0; j < m_VectorOrder[i].size(); ++j) {
            std::cout << m_VectorOrder[i][j] << " ";
        }
        std::cout << std::endl;  // Newline after each inner vector
    }

    GenerateZeroAndNonZeroVectorOrder();

    


}
FluctuationSpectrum::~FluctuationSpectrum() {
    
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
}

void FluctuationSpectrum::GenerateZeroAndNonZeroVectorOrder(){

    for (size_t i = 0; i < m_VectorOrder.size(); ++i) {
        std::vector<std::vector<int>> MatrixOrder(4, std::vector<int>(2, 0));

        std::cout << "m_VectorOrder[" << i << "]: ["
                  << m_VectorOrder[i][0] << ", " << m_VectorOrder[i][1] << "]" << std::endl;


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
            //Do something
        }

        std::cout << "m_MatrixOrder[" << i << "]: " << std::endl;
        for (size_t j = 0; j < m_MatrixOrder[i].size(); ++j) {
            std::cout << "[";
            for (size_t k = 0; k < m_MatrixOrder[i][j].size(); ++k) {
                std::cout << m_MatrixOrder[i][j][k];
                if (k < m_MatrixOrder[i][j].size() - 1) std::cout << ", ";
            }
            std::cout << "]" << std::endl;
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
                std::cout<<m_MatrixOrder[i][j][0]<<m_MatrixOrder[i][j][1]<<std::endl;
                qmatrix[j][0]=static_cast<double>(m_MatrixOrder[i][j][0])*2*PI/Lx;
                qmatrix[j][1]=static_cast<double>(m_MatrixOrder[i][j][1])*2*PI/Ly;
            }

            hqvector[i]=FourierTransform(qmatrix);
        }   
        else{
            std::vector<std::vector<double>> qmatrix(4,std::vector<double>(2,0));
            for (size_t j = 0; j < 4; ++j) {
                std::cout<<m_MatrixOrder[i][j][0]<<m_MatrixOrder[i][j][1]<<std::endl;
                qmatrix[j][0]=static_cast<double>(m_MatrixOrder[i][j][0])*2*PI/Lx;
                qmatrix[j][1]=static_cast<double>(m_MatrixOrder[i][j][1])*2*PI/Ly;
            }

            hqvector[i]=FourierTransform(qmatrix);  
            }

    }

    std::cout << "Matrix order: ";
    for (size_t i = 0; i < qvector.size(); ++i) {
        std::cout << m_MatrixOrder[i][0][0] <<m_MatrixOrder[i][0][1] << " ";
    }
    std::cout << std::endl; // New line for better readability
    std::cout << std::endl; // New line for better readability
    std::cout << "qvector: ";
    for (size_t i = 0; i < qvector.size(); ++i) {
        std::cout << qvector[i] << " ";
    }
    std::cout << std::endl; // New line for better readability

    std::cout << "hqvector: ";
    for (size_t i = 0; i < hqvector.size(); ++i) {
        std::cout << hqvector[i] << " ";
    }
    std::cout << std::endl; // New line for better readability


}


