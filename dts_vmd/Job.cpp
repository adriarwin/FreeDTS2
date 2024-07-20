

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com && weria.pezeshkian@nbi.ku.dk)
 Copyright (c) Weria Pezeshkian
 This class finds checks the name of the executable, Not a very important task for the current version
 */

#include "Job.h"
#include "State.h"
#include "Nfunction.h"
#include "Analyze.h"
#include "RNG.h"


Job::Job(std::vector <std::string> argument)
{
std::string Exe=argument.at(0);
    std::string ExeName = "ANA";
    Nfunction f;

    if (Exe.size()>3 && Exe.at(Exe.size()-4)!='/')
    {
        std::cout<<" Error:  (b) unrecognized exacutable name --->"<<Exe<<" :( "<<std::endl;
    }
    else
    {
        char L1 = Exe.at(Exe.size()-1);
        char L2 = Exe.at(Exe.size()-2);
        char L3 = Exe.at(Exe.size()-3);
        if(L1==ExeName.at(2) && L2==ExeName.at(1) && L3==ExeName.at(0))
        {
            State tem_S(argument);
            Analyze ANA(&tem_S);
        }
    }
}
Job::~Job()
{
    
}

