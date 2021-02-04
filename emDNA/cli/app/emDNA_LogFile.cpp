// emDNA_LogFile class
// Nicolas Clauvelin


#include <emDNA_CommandLine.h>
//#include <BpCollection.h>
//#include <IndexManager.h>
//#include <MinimizerAgent.h>
//#include <BpCollection_x3DNA.h>
//#include <BpCollectionElasticEnergy.h>
#include <emDNA_LogFile.h>


// constructors
emDNA_LogFile::emDNA_LogFile() : m_log_file() {};


// class constructor by copy
emDNA_LogFile::emDNA_LogFile(const emDNA_LogFile& log_file) :
m_log_file() {
    m_log_file.set_file_name(log_file.m_log_file.file_name());
    m_log_file.set_file_path(log_file.m_log_file.file_path());
};


// class destructor
emDNA_LogFile::~emDNA_LogFile() {
    if (m_log_file.is_open())
        m_log_file.close();
};

// copy operator
emDNA_LogFile& emDNA_LogFile::operator=(const emDNA_LogFile& log_file) {
    m_log_file.set_file_name(log_file.m_log_file.file_name());
    m_log_file.set_file_path(log_file.m_log_file.file_path());
    return *this;
};


// log file creation method
void emDNA_LogFile::create_log_file(const std::string& filename) {

    // file name
    m_log_file.set_file_name(filename);
    m_log_file.set_file_path(".");

    // date
    std::string timestamp("log file timestamp: ");
    timestamp += DateParser::get_date();
    timestamp += "\n\n";

    // file creation
    m_log_file.open();
    m_log_file.write_line(timestamp);
    m_log_file.close();

};


// command line input logging method
void emDNA_LogFile::log_emDNA_CommandLine(const emDNA_CommandLine& cmd_line) {

    // stream
    std::stringstream s;
    s << ">> command line input\n";
    cmd_line.report_command_line_input_data(s);
    s << "\n";

    // logging
    m_log_file.open_in_append_mode();
    m_log_file.write_stream(s);
    m_log_file.close();

};


// bp collection logging method
void emDNA_LogFile::log_BpCollection(const BpCollection& bp_collection) {

    // stream
    std::stringstream s;

    // header
    s << ">> initial bp collection\n";

    // bp collection data
    s << "number of base pairs: " << bp_collection.n_of_base_pairs() << "\n";
    s << "number of bp steps: " << bp_collection.n_of_bp_steps() << "\n";
    s << "sequence: " << bp_collection.collection_sequence() << "\n";
    s << "step parameters:\n";
    BpCollection_x3DNA::format_as_step_parameters_file(s, bp_collection);
    s << "\n";

    // logging
    m_log_file.open_in_append_mode();
    m_log_file.write_stream(s);
    m_log_file.close();

};


// minimization results
void emDNA_LogFile::log_MinimResults(const MinimizationResults& minim_res) {

    // stream
    std::stringstream s;

    // header
    s << ">> minimization results\n";

    // results data
    s << std::fixed << std::setprecision(REAL_WIDTH);
    s << "minimization description: " << minim_res._minim_desc << "\n";
    s << "initial energy: " << minim_res._initial_E << "\n";
    s << "final energy: " << minim_res._final_E << "\n";
    s << "free dofs gradient norm: "
        << minim_res._free_dofs_gradient.norm() << "\n";
    s << "# iterations: " << minim_res._n_iterations << "\n";
    s << "return code: " << minim_res._return_code << "\n";
    s << "\n";

    // logging
    m_log_file.open_in_append_mode();
    m_log_file.write_stream(s);
    m_log_file.close();

};


// energy composition logging method
void emDNA_LogFile::log_EnergyComposition(const MinimizationResults& minim_res)
{

    // optimized bp collection
    const BpCollection& bp_collection_opt = minim_res._optimized_bp_collection;

    // index manager
    IndexManager idx_mgr;
    idx_mgr.set_n_of_bp_step(bp_collection_opt.n_of_bp_steps());
    idx_mgr.set_frozen_steps_indexes(bp_collection_opt.frozen_steps_domains());

    // container
    const MatrixN energy_split =
    BpCollectionElasticEnergy::elastic_energy_splits(bp_collection_opt,
                                                     idx_mgr);

    // stream
    std::stringstream s;

    // energy contribs header
    s << ">> energy contributions details\n";
    for (Size i=0; i<minim_res._energy_contribs.size(); ++i) {
        s << minim_res._energy_contribs[i].first;
        s << " -> ";
        s << minim_res._energy_contribs[i].second;
        s << "\n";
    };
    s << "\n";

    // elastic energy splits
    s << ">> elastic energy contributions per step parameter pairs\n";
    s << std::fixed << std::setprecision(REAL_WIDTH);
    for (Size i=0; i<emDNAConstants::StepParametersDim; ++i) {
        s << energy_split.row(i) << "\n";
    };
    s << "\n";

    // logging
    m_log_file.open_in_append_mode();
    m_log_file.write_stream(s);
    m_log_file.close();

};

