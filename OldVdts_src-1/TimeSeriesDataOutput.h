#if !defined(AFX_TimeSeriesDataOutput_H_INCLUDED_)
#define AFX_TimeSeriesDataOutput_H_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This for reading and writing the TimeSeriesDataOutput file (a binary file)
 */
#include <fstream>

class State;
class TimeSeriesDataOutput
{
public:
    TimeSeriesDataOutput();
	TimeSeriesDataOutput(State* pstate);
	 ~TimeSeriesDataOutput();
    
        inline int GetPeriodic()                 const      {return m_Periodic;}

public:
    
    bool WrireTimeSeriesDataOutput(int step);
    bool OpenTimeSeriesDataOutput();
    bool FlushTimeSeriesFile(); // the eneergy file should be flushed first
    void UpdatePeriod(int period);
    bool CheckTimeSeriesFile(int endstep); // this checks if the energy file matches the restart end step
private:
    State *m_pState;
    int m_Periodic;
    std::ofstream m_TimeSeriesFile;  // Member variable to hold the file stream



};

#endif
