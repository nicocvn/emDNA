// DNASim library header.
// Nicolas Clauvelin (n.clauvelin@gmail.com)


#ifndef EMDNA_DNASIM_H
#define EMDNA_DNASIM_H


#include "DNASim_Includes.h"


// dna
#include "dna/BpGeometry.h"
#include "dna/ForceConstantsDB.h"
#include "dna/ForceConstants_AnisoDNA.h"
#include "dna/ForceConstants_IdealDNA.h"
#include "dna/ForceConstants_Olson1998.h"
#include "dna/ForceConstants_Cohen2017_tet.h" //Added by Zoe Wefers (McGill University, July 2021, DIMACS REU)
#include "dna/ForceConstants_Cohen2017_dim.h" //Added by Zoe Wefers (McGill University, July 2021, DIMACS REU)
#include "dna/Parser_x3DNA.h"
#include "dna/Sequence.h"
#include "dna/SequenceDepenceModels.h"
#include "dna/TetramerDepenceModels.h" //Added by Zoe Wefers (McGill University, July 2021, DIMACS REU)
#include "dna/StepArray.h"
#include "dna/StepParameters.h"
#include "dna/StepParametersDB.h"
#include "dna/StepParameters_AnisoDNA.h"
#include "dna/StepParameters_AnisoDNA_304.h"
#include "dna/StepParameters_IdealDNA.h"
#include "dna/StepParameters_IdealDNA_304.h"
#include "dna/StepParameters_Olson1998.h"
#include "dna/StepParameters_Cohen2017_tet.h" //Added by Zoe Wefers (McGill University, July 2021, DIMACS REU)
#include "dna/StepParameters_Cohen2017_dim.h" //Added by Zoe Wefers (McGill University, July 2021, DIMACS REU)

// file_io
#include "file_io/EnhancedString.h"
#include "file_io/FileException.h"
#include "file_io/FileHandler.h"
#include "file_io/InputFileHandler.h"
#include "file_io/OutputFileHandler.h"

// geometry
#include "geometry/BpStepTwistDensity.h"
#include "geometry/CurveTopology.h"
#include "geometry/CylinderShape.h"
#include "geometry/DiscreteRibbon.h"
#include "geometry/Segment.h"
#include "geometry/Shape.h"
#include "geometry/SphereShape.h"

// maths
#include "maths/AffineTransformation.h"
#include "maths/Eigen_Includes.h"
#include "maths/Matrix3.h"
#include "maths/Matrix4.h"
#include "maths/MatrixN.h"
#include "maths/SparseMatrixEntries.h"
#include "maths/SparseMatrixN.h"
#include "maths/Triad.h"
#include "maths/Vector3.h"
#include "maths/VectorN.h"

// pdb
#include "pdb/PDBAtom.h"
#include "pdb/PDBAtomRecord.h"
#include "pdb/PDBChain.h"
#include "pdb/PDBParser.h"
#include "pdb/PDBResidue.h"
#include "pdb/PDBStructure.h"

// prn-generators
#include "prn-generators/PRN_Integer.h"
#include "prn-generators/PRN_Real.h"
#include "prn-generators/Seed.h"
#include "prn-generators/dsUUID.h"

// serialization
#include "serialization/Cereal_Includes.h"
#include "serialization/DNASimArchive.h"
#include "serialization/DNASimSerialization.h"

// simulations
#include "simulations/DateParser.h"
#include "simulations/OptionsManager.h"
#include "simulations/PointCloud.h"
#include "simulations/PointCloudAdaptor.h"
#include "simulations/PointCloudKdTree.h"
#include "simulations/SparseIndex.h"
#include "simulations/Trigger.h"


// Register polymorphic types for Cereal serialization.
CEREAL_REGISTER_TYPE(DNASim::SphereShape);
CEREAL_REGISTER_TYPE(DNASim::CylinderShape);
CEREAL_REGISTER_POLYMORPHIC_RELATION(DNASim::Shape, DNASim::SphereShape)
CEREAL_REGISTER_POLYMORPHIC_RELATION(DNASim::Shape, DNASim::CylinderShape)


#endif    // EMDNA_DNASIM_H
