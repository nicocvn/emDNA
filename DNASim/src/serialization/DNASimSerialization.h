// DNASimSerialization class
// Nicolas Clauvelin


// cereal serialization load/save functions

// to use serialization:
// - include DNASimSerialization header to get load/save functions
// - include DNASimArchive header
// - declare polymorphic types if needed

// for polymorphic types registration is needed:
// CEREAL_REGISTER_TYPE(Type);
// required for: Shape, SphereShape, CylinderShape


#ifndef DNASim_DNASimSerialization_h
#define DNASim_DNASimSerialization_h


#include "serialization/Cereal_Includes.h"
#include "file_io/EnhancedString.h"
#include "maths/Vector3.h"
#include "maths/VectorN.h"
#include "maths/Matrix3.h"
#include "maths/Matrix4.h"
#include "maths/MatrixN.h"
#include "maths/Triad.h"
#include "maths/AffineTransformation.h"
#include "dna/StepParameters.h"
#include "dna/StepArray.h"
#include "geometry/SphereShape.h"
#include "geometry/CylinderShape.h"
#include "simulations/OptionsManager.h"
#include "dna/Sequence.h"
#include "dna/StepParametersDB.h"
#include "dna/ForceConstantsDB.h"


namespace DNASim {


    // EnhancedString
    template <class Archive>
    void save(Archive& archive, const EnhancedString& es) {
        archive(cereal::make_nvp("string", std::string(es)));
    };
    template <class Archive>
    void load(Archive& archive, EnhancedString& es) {
        std::string es_str;
        archive(es_str);
        es = EnhancedString(es_str);
    };


    // Vector3
    template <class Archive>
    void save(Archive& archive, const Vector3& v) {
        archive(cereal::make_nvp("X", v[X]),
                cereal::make_nvp("Y", v[Y]),
                cereal::make_nvp("Z", v[Z]));
    };
    template <class Archive>
    void load(Archive& archive, Vector3& v) {
        archive(v[X], v[Y], v[Z]);
    };


    // VectorN
    template <class Archive>
    void save(Archive& archive, const VectorN& v) {
        const std::vector<Real> std_v = v.std_vector();
        archive(cereal::make_nvp("vector", std_v));
    };
    template <class Archive>
    void load(Archive& archive, VectorN& v) {
        std::vector<Real> std_vec;
        archive(std_vec);
        v = VectorN(std_vec);
    };


    // Matrix3
    template <class Archive>
    void save(Archive& archive, const Matrix3& m) {
        archive(CEREAL_NVP(m(0,0)),
                CEREAL_NVP(m(0,1)),
                CEREAL_NVP(m(0,2)),
                CEREAL_NVP(m(1,0)),
                CEREAL_NVP(m(1,1)),
                CEREAL_NVP(m(1,2)),
                CEREAL_NVP(m(2,0)),
                CEREAL_NVP(m(2,1)),
                CEREAL_NVP(m(2,2)));
    };
    template <class Archive>
    void load(Archive& archive, Matrix3& m) {
        archive(m(0,0), m(0,1), m(0,2),
                m(1,0), m(1,1), m(1,2),
                m(2,0), m(2,1), m(2,2));
    };


    // Matrix4
    template <class Archive>
    void save(Archive& archive, const Matrix4& m) {
        archive(CEREAL_NVP(m(0,0)),
                CEREAL_NVP(m(0,1)),
                CEREAL_NVP(m(0,2)),
                CEREAL_NVP(m(0,3)),
                CEREAL_NVP(m(1,0)),
                CEREAL_NVP(m(1,1)),
                CEREAL_NVP(m(1,2)),
                CEREAL_NVP(m(1,3)),
                CEREAL_NVP(m(2,0)),
                CEREAL_NVP(m(2,1)),
                CEREAL_NVP(m(2,2)),
                CEREAL_NVP(m(2,3)),
                CEREAL_NVP(m(3,0)),
                CEREAL_NVP(m(3,1)),
                CEREAL_NVP(m(3,2)),
                CEREAL_NVP(m(3,3)));
    };
    template <class Archive>
    void load(Archive& archive, Matrix4& m) {
        archive(m(0,0), m(0,1), m(0,2), m(0,3),
                m(1,0), m(1,1), m(1,2), m(1,3),
                m(2,0), m(2,1), m(2,2), m(2,3),
                m(3,0), m(3,1), m(3,2), m(3,3));
    };


    // MatrixN
    template <class Archive>
    void save(Archive& archive, const MatrixN& m) {

        // sizing data
        const Size n_rows = m.n_rows();
        const Size n_cols = m.n_cols();

        // values
        std::vector<Real> mat_values;
        mat_values.reserve(n_rows*n_cols);
        for (Size i=0; i<n_rows; ++i)
            for (Size j=0; j<n_cols; ++j)
                mat_values.push_back(m(i,j));

        // archive
        archive(CEREAL_NVP(n_rows),
                CEREAL_NVP(n_cols),
                cereal::make_nvp("inline_values", mat_values));

    };
    template <class Archive>
    void load(Archive& archive, MatrixN& m) {

        // sizing data
        Size n_rows, n_cols;

        // values
        std::vector<Real> mat_values;

        // loading
        archive(n_rows, n_cols, mat_values);

        // filling
        m = MatrixN(n_rows, n_cols);
        for (Size i=0; i<n_rows; ++i)
            for (Size j=0; j<n_cols; ++j)
                m(i,j) = mat_values[i*n_rows+j];

    };


    // Triad
    template <class Archive>
    void save(Archive& archive, const Triad& frame) {
        archive(cereal::make_nvp("matrix", frame.matrix_representation()));
    };
    template <class Archive>
    void load(Archive& archive, Triad& frame) {
        Matrix4 frame_mat;
        archive(frame_mat);
        frame = Triad(frame_mat);
    };


    // AffineTransformation
    template <class Archive>
    void save(Archive& archive, const AffineTransformation& aff_t) {
        archive(cereal::make_nvp("transformation_matrix", aff_t.matrix()));
    };
    template <class Archive>
    void load(Archive& archive, AffineTransformation& aff_t) {
        Matrix4 aff_t_mat;
        archive(aff_t_mat);
        aff_t = AffineTransformation(aff_t_mat);
    };


    // StepParameters
    template <class Archive>
    void save(Archive& archive, const StepParameters& p) {
        archive(cereal::make_nvp("tilt", p[TILT]),
                cereal::make_nvp("roll", p[ROLL]),
                cereal::make_nvp("twist", p[TWIST]),
                cereal::make_nvp("shift", p[SHIFT]),
                cereal::make_nvp("slide", p[SLIDE]),
                cereal::make_nvp("rise", p[RISE]));
    };
    template <class Archive>
    void load(Archive& archive, StepParameters& p) {
        archive(p[TILT], p[ROLL], p[TWIST],
                p[SHIFT], p[SLIDE], p[RISE]);
    };


    // StepArray
    template <class Archive>
    void save(Archive& archive, const StepArray& step_array) {
        const Size n_steps = step_array.n_step();
        archive(CEREAL_NVP(n_steps));
        for (Size i=0; i<n_steps; ++i)
            archive(cereal::make_nvp("step_" +
                                     EnhancedString::convert_to_string(i),
                                     step_array.parameters(i)));
    };
    template <class Archive>
    void load(Archive& archive, StepArray& step_array) {
        step_array.clear();
        Size n_steps;
        archive(n_steps);
        for (Size i=0; i<n_steps; ++i) {
            StepParameters p;
            archive(p);
            step_array.append_step_parameters(p);
        };
    };


    // Shape
    template <class Archive>
    void save(Archive& archive, const DNASim::Shape& shape) {
        archive(cereal::make_nvp("shape_frame", shape.frame()));
    };
    template <class Archive>
    void load(Archive& archive, DNASim::Shape& shape) {
        Triad shape_frame;
        archive(shape_frame);
        shape.set_frame(shape_frame);
    };


    // SphereShape
    template <class Archive>
    void save(Archive& archive, const SphereShape& sphere) {
        archive(cereal::make_nvp("base_data",
                                 cereal::base_class<Shape>(&sphere)),
                cereal::make_nvp("radius", sphere.radius()));
    };
    template <class Archive>
    void load(Archive& archive, DNASim::SphereShape& sphere) {
        Real sphere_radius;
        archive(cereal::virtual_base_class<Shape>(&sphere), sphere_radius);
        sphere.set_radius(sphere_radius);
    };


    // CylinderShape
    template <class Archive>
    void save(Archive& archive, const CylinderShape& cylinder) {
        archive(cereal::make_nvp("base_data",
                                 cereal::base_class<Shape>(&cylinder)),
                cereal::make_nvp("radius", cylinder.radius()),
                cereal::make_nvp("height", cylinder.height()));
    };
    template <class Archive>
    void load(Archive& archive, CylinderShape& cylinder) {
        Real cylinder_radius;
        Real cylinder_height;
        archive(cereal::virtual_base_class<Shape>(&cylinder),
                cylinder_radius,
                cylinder_height);
        cylinder.set_radius(cylinder_radius);
        cylinder.set_height(cylinder_height);
    };


    // OptionsManager
    template <class Archive>
    void save(Archive& archive, const OptionsManager& opt_mgr) {
        std::map<std::string,std::string> opts(opt_mgr.begin(), opt_mgr.end());
        archive(cereal::make_nvp("options_map", opts));
    };
    template <class Archive>
    void load(Archive& archive, OptionsManager& opt_mgr) {
        std::map<std::string,std::string> opts;
        archive(opts);
        opt_mgr = OptionsManager(opts);
    };


    // StepSequence
    template <class Archive>
    void save(Archive& archive, const StepSequence& seq) {
        archive(cereal::make_nvp("first_base",
                                 static_cast<Size>(seq.first_base())));
        archive(cereal::make_nvp("last_base",
                                 static_cast<Size>(seq.last_base())));
    };
    template <class Archive>
    void load(Archive& archive, StepSequence& seq) {
        Size first, last;
        archive(first);
        archive(last);
        seq = StepSequence(static_cast<BaseSymbol>(first),
                           static_cast<BaseSymbol>(last));
    };


    // Sequence
    template <class Archive>
    void save(Archive& archive, const Sequence& sequence) {
        archive(cereal::make_nvp("linear_sequence",
                                 sequence.linear_sequence()));
    };
    template <class Archive>
    void load(Archive& archive, Sequence& sequence) {
        std::string linear_sequence;
        archive(linear_sequence);
        sequence = Sequence(linear_sequence);
    };


    // StepParametersDB
#define SEQ_DIM 4

    //Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    template <class ArchiveType>
    void save(ArchiveType& archive, const StepParametersDB& steps_db) {

        // model name
        archive(cereal::make_nvp("model_name", steps_db.model_name()));

        // intrinsic step parameters
        for (Size i=0; i<SEQ_DIM+1; ++i)
            for (Size j=0; j<SEQ_DIM; ++j)
                for (Size k=0; k<SEQ_DIM; ++k)
                    for (Size l=0; l<SEQ_DIM+1; ++l) {

                        // nucleotides
                        const BaseSymbol first = static_cast<BaseSymbol>(i);
                        const BaseSymbol second = static_cast<BaseSymbol>(j);
                        const BaseSymbol third = static_cast<BaseSymbol>(k);
                        const BaseSymbol fourth = static_cast<BaseSymbol>(l);
                        const TetramerSequence tetra_seq(first, second, third, fourth);

                        // save
                        archive(cereal::make_nvp("first_base",
                                                    Base::str(first)),
                                cereal::make_nvp("second_base",
                                                    Base::str(second)),
                                cereal::make_nvp("third_base",
                                                    Base::str(third)),
                                cereal::make_nvp("fourth_base",
                                                    Base::str(fourth)),
                                cereal::make_nvp("intrinsic_step",
                                                    steps_db.
                                                    intrinsic_bp_step_params(tetra_seq)));
                    }

    };

    //Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    template <class ArchiveType>
    void load(ArchiveType& archive, StepParametersDB& steps_db) {

        // model name
        std::string model_name;
        archive(model_name);

        // intrinsic step parameters
        StepParameters db_data[SEQ_DIM+1][SEQ_DIM][SEQ_DIM][SEQ_DIM+1];
        for (Size i=0; i<SEQ_DIM+1; ++i)
            for (Size j=0; j<SEQ_DIM; ++j)
                for (Size k=0; k<SEQ_DIM; ++k)
                    for (Size l=0; l<SEQ_DIM+1; ++l) {
                        
                        // load
                        std::string first_str, second_str, third_str, fourth_str;
                        StepParameters step_params;
                        archive(first_str, second_str, third_str, fourth_str, step_params);

                        // nucleotide types
                        const char s1 = first_str.c_str()[0];
                        const char s2 = second_str.c_str()[0];
                        const char s3 = third_str.c_str()[0];
                        const char s4 = fourth_str.c_str()[0];
                        BaseSymbol first = Base::base_symbol_from_char(s1);
                        BaseSymbol second = Base::base_symbol_from_char(s2);
                        BaseSymbol third = Base::base_symbol_from_char(s3);
                        BaseSymbol fourth = Base::base_symbol_from_char(s4);
                        
                        // array filling
                        db_data[(Size)first][(Size)second][(Size)third][(Size)fourth] = step_params;             
                    }
        // setup
        steps_db = StepParametersDB(model_name, db_data);
        
    };
#undef SEQ_DIM


    // ForceConstantsDB
#define SEQ_DIM 4

    //Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    template <class ArchiveType>
    void save(ArchiveType& archive, const ForceConstantsDB& fmat_db) {

        // model name
        archive(cereal::make_nvp("model_name", fmat_db.model_name()));

        // force constants
        for (Size i=0; i<SEQ_DIM+1; ++i)
            for (Size j=0; j<SEQ_DIM; ++j)
                for (Size k=0; k<SEQ_DIM; ++k)
                    for (Size l=0; l<SEQ_DIM+1; ++l) {

                        // nucleotides
                        const BaseSymbol first = static_cast<BaseSymbol>(i);
                        const BaseSymbol second = static_cast<BaseSymbol>(j);
                        const BaseSymbol third = static_cast<BaseSymbol>(k);
                        const BaseSymbol fourth = static_cast<BaseSymbol>(l);
                        const TetramerSequence tetra_seq(first, second, third, fourth);

                        // save
                        archive(cereal::make_nvp("first_base",
                                                    Base::str(first)),
                                cereal::make_nvp("second_base",
                                                    Base::str(second)),
                                cereal::make_nvp("third_base",
                                                    Base::str(third)),
                                cereal::make_nvp("fourth_base",
                                                    Base::str(fourth)),
                                cereal::make_nvp("force_constants",
                                         fmat_db.force_constants(tetra_seq)));
                    }

    };

    //Changed by Zoe Wefers (McGill University, June 2021, DIMACS REU)
    template <class ArchiveType>
    void load(ArchiveType& archive, ForceConstantsDB& fmat_db) {

        // model name
        std::string model_name;
        archive(model_name);

        // intrinsic step parameters
        MatrixN db_data[SEQ_DIM+1][SEQ_DIM][SEQ_DIM][SEQ_DIM+1];
        for (Size i=0; i<SEQ_DIM+1; ++i)
            for (Size j=0; j<SEQ_DIM; ++j)
                for (Size k=0; k<SEQ_DIM; ++k)
                    for (Size l=0; l<SEQ_DIM+1; ++l) {

                        // load
                        std::string first_str, second_str, third_str, fourth_str;
                        MatrixN fmat_vals;
                        archive(first_str, second_str, third_str, fourth_str, fmat_vals);

                        // nucleotide types
                        const char s1 = first_str.c_str()[0];
                        const char s2 = second_str.c_str()[0];
                        const char s3 = third_str.c_str()[0];
                        const char s4 = fourth_str.c_str()[0];
                        BaseSymbol first = Base::base_symbol_from_char(s1);
                        BaseSymbol second = Base::base_symbol_from_char(s2);
                        BaseSymbol third = Base::base_symbol_from_char(s3);
                        BaseSymbol fourth = Base::base_symbol_from_char(s4);
                        
                        // array filling
                        db_data[(Size)first][(Size)second][(Size)third][(Size)fourth] = fmat_vals;
                        
                    }
        
        // setup
        fmat_db = ForceConstantsDB(model_name, db_data);
        
    };
#undef SEQ_DIM


}


#endif  // DNASim_DNASimSerialization_h
