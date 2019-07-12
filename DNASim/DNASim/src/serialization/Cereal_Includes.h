// Cereal Includes header file
// Nicolas Clauvelin


// project-wide cereal include directives


#ifndef DNASim_Cereal_Includes_h
#define DNASim_Cereal_Includes_h


// suppress warnings
#pragma clang diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wsign-compare"


// archive headers
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>

// archive types
using CerealBinaryOutput = cereal::PortableBinaryOutputArchive;
using CerealBinaryInput = cereal::PortableBinaryInputArchive;
using CerealJSONOutput = cereal::JSONOutputArchive;
using CerealJSONInput = cereal::JSONInputArchive;
using CerealXMLOutput = cereal::XMLOutputArchive;
using CerealXMLInput = cereal::XMLInputArchive;

// stl types serialization
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>

// inheritance and polymorphism
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>


#endif  // DNASim_Cereal_Includes_h
