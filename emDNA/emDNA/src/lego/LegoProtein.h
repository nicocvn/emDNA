// LegoProtein class
// Nicolas Clauvelin


// base implementation of a LEGO protein

// each LEGO protein is implemented as a collection of step parameters
//
// a single vector of intrinsic step parameters is required and additional
// parameters depending on the specfic implementation


#ifndef emDNA_LegoProtein_h
#define emDNA_LegoProtein_h


#include <emDNA_Includes.h>
#include <BpCollection.h>


class LegoProtein {


public:

    // virtual destructor
    virtual ~LegoProtein() = default;

    // bp step parameters LEGO method
    const std::vector<BpStepParams>& bp_step_params() const {
        return m_lego_coll.bp_step_params();
    };

    // base pairs LEGO method
    std::vector<BasePair> base_pairs() const {
        return m_lego_coll.base_pairs();
    };


protected:

    // constructors
    LegoProtein() = default;
    LegoProtein(const LegoProtein& lego_protein) = default;
    LegoProtein(LegoProtein&& lego_protein) = default;

    // copy and move operators
    LegoProtein& operator=(const LegoProtein& lego_protein) = default;
    LegoProtein& operator=(LegoProtein&& lego_protein) = default;

    // step parameters computation method
    virtual void compute_bp_collection() =0;

    // LEGO step parameters data
    BpCollection m_lego_coll;


};


#endif  // emDNA_LegoProtein_h
