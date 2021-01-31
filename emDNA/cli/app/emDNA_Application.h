// emDNA_Application class
// Nicolas Clauvelin


// application class


#ifndef emDNA_emDNA_Application_h
#define emDNA_emDNA_Application_h


#include <emDNA_CommandLine.h>
#include <emDNA_LogFile.h>
#include <BpCollection.h>
#include <BpCollection_Interface.h>
struct AlglibMinSettings;
using AlglibMinSettings_Ptr = std::unique_ptr<const AlglibMinSettings>;


// emDNA_Application class
class emDNA_Application {


public:

    // constructors
    emDNA_Application();
    ~emDNA_Application();

    // application command line parsing
    void parse_command_line(int argc, char* argv[]);

    // application setup method
    void setup_application();

    // minimization method
    void perform_minimization();

    // results output method
    void output_optimized_bp_collection();


private:

    // setup methods
    void setup_bp_collection();
    void setup_log_file(const std::string& filename);

    // message methods
    void application_setup_message() const;
    void minimization_pre_flight_message(const AlglibMinSettings_Ptr&
                                         m_set) const;
    void minimization_post_flight_message(const MinimizationResults& res) const;
    void apllication_tear_down_message() const;

    // static method for minimization type description
    static std::string minimization_type_description(const CollectionType&
                                                     coll_type);

    // application data
    emDNA_CommandLine m_cmd_line;       // command line parser
    emDNA_LogFile m_log_file;           // log file
    BpCollection m_bp_collection;       // bp collection
    std::shared_ptr<BpCollection_Interface> m_bp_intf_ptr;


};


// entry point
int main(int argc, char* argv[]);


#endif  // emDNA_emDNA_Application_h
