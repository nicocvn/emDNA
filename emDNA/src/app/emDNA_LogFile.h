// emDNA_LogFile class
// Nicolas Clauvelin


// class for log file management


#ifndef emDNA_emDNA_LogFile_h
#define emDNA_emDNA_LogFile_h


#include <emDNA_Includes.h>
class emDNA_CommandLine;
class BpCollection;
class MinimResults;


class emDNA_LogFile {


public:

    // constructors
    emDNA_LogFile();
    emDNA_LogFile(const emDNA_LogFile& log_file);
    ~emDNA_LogFile();

    // copy operator
    emDNA_LogFile& operator=(const emDNA_LogFile& log_file);

    // log file creation method
    void create_log_file(const std::string& filename);

    // logging methods
    void log_emDNA_CommandLine(const emDNA_CommandLine& cmd_line);
    void log_BpCollection(const BpCollection& bp_collection);
    void log_MinimResults(const MinimizationResults& minim_res);
    void log_EnergyComposition(const MinimizationResults& minim_res);


private:

    // log file
    OutputFileHandler m_log_file;


};


#endif  // emDNA_emDNA_LogFile_h
