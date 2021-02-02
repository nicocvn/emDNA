// test_DNASimArchive class
// Nicolas Clauvelin


#include "test_DNASimArchive.h"


namespace {


    // Vector3SerializationTest
    TEST_F(DNASimArchiveTest, Vector3SerializationTest) {

        // output
        {
            OutputArchive<CerealBinaryOutput> bin_archive("test.cereal", false);
            bin_archive.save(v3);
        }

        // input
        Vector3 read_v;
        InputArchive<CerealBinaryInput> bin_archive("test.cereal");
        bin_archive.load(read_v);

        EXPECT_NEAR((read_v-v3).norm(), FLOAT_INIT, FLOAT_INIT);

    };


    // VectorNSerializationTest
    TEST_F(DNASimArchiveTest, VectorNSerializationTest) {

        // output
        {
            OutputArchive<CerealBinaryOutput> bin_archive("test.cereal", false);
            bin_archive.save(vN);
        }

        // input
        VectorN read_v;
        InputArchive<CerealBinaryInput> bin_archive("test.cereal");
        bin_archive.load(read_v);

        EXPECT_NEAR((read_v-vN).norm(), FLOAT_INIT, FLOAT_INIT);

    };


    // TriadSerializationTest
    TEST_F(DNASimArchiveTest, TriadSerializationTest) {

        // output
        {
            OutputArchive<CerealBinaryOutput> bin_archive("test.cereal", false);
            bin_archive.save(f);
        }

        // input
        Triad read_f;
        InputArchive<CerealBinaryInput> bin_archive("test.cereal");
        bin_archive.load(read_f);

        EXPECT_NEAR((read_f.origin()-f.origin()).norm(), FLOAT_INIT,
                    FLOAT_INIT);
        EXPECT_NEAR((read_f.axis(I)-f.axis(I)).norm(), FLOAT_INIT,
                    FLOAT_INIT);
        EXPECT_NEAR((read_f.axis(J)-f.axis(J)).norm(), FLOAT_INIT,
                    FLOAT_INIT);
        EXPECT_NEAR((read_f.axis(K)-f.axis(K)).norm(), FLOAT_INIT,
                    FLOAT_INIT);

    };


    // StepParametersSerializationTest
    TEST_F(DNASimArchiveTest, StepParametersSerializationTest) {

        // output
        {
            OutputArchive<CerealBinaryOutput> bin_archive("test.cereal", false);
            bin_archive.save(p);
        }

        // input
        StepParameters read_p;
        InputArchive<CerealBinaryInput> bin_archive("test.cereal");
        bin_archive.load(read_p);

        EXPECT_NEAR((read_p.inline_vector()-p.inline_vector()).norm(),
                    FLOAT_INIT,
                    FLOAT_INIT);
        
    };


    // ShapeSerializationTest
    TEST_F(DNASimArchiveTest, ShapeSerializationTest) {

        // shape frame
        Triad f_shape;
        f_shape.set_origin(Vector3(Real(-0.626757), Real(-0.883724),
                                   Real(3.0627)));
        f_shape.set_axis(I, Vector3(Real(0.843184), Real(0.537622),
                                    Real(-0.00158442)));
        f_shape.set_axis(J, Vector3(Real(-0.537064), Real(0.842436),
                                    Real(0.043173)));
        f_shape.set_axis(K, Vector3(Real(0.0245455), Real(-0.0355518),
                                    Real(0.999066)));
        f_shape.orthogonalize();

        // dimensions
        Real sphere_radius = Real(1.27);
        Real cylinder_radius = Real(3.45);
        Real cylinder_height = Real(4.1);

        // sphere and cylinder
        std::shared_ptr<SphereShape> sphere =
        std::make_shared<SphereShape>(f_shape.origin(), sphere_radius);
        std::shared_ptr<CylinderShape> cylinder =
        std::make_shared<CylinderShape>(f_shape,
                                        cylinder_radius, cylinder_height);

        // output
        {
            OutputArchive<CerealBinaryOutput> bin_archive("test.cereal", false);
            bin_archive.save(sphere);
            bin_archive.save(cylinder);
        }

        // input
        std::shared_ptr<SphereShape> sphere_ptr;
        std::shared_ptr<CylinderShape> cylinder_ptr;
        InputArchive<CerealBinaryInput> bin_archive("test.cereal");
        bin_archive.load(sphere_ptr);
        bin_archive.load(cylinder_ptr);

        EXPECT_NEAR((sphere_ptr->frame().origin()-f_shape.origin()).norm(),
                    FLOAT_INIT,
                    FLOAT_INIT);
        EXPECT_NEAR((cylinder_ptr->frame().axis(J)-f_shape.axis(J)).norm(),
                    FLOAT_INIT,
                    FLOAT_INIT);
        EXPECT_NEAR((cylinder_ptr->height()-cylinder_height),
                    FLOAT_INIT,
                    FLOAT_INIT);

    };


    // StepArraySerializationTest
    TEST_F(DNASimArchiveTest, StepArraySerializationTest) {

        StepParameters p;
		p[TILT]=Real(2.35);
		p[ROLL]=Real(0.78);
		p[TWIST]=Real(32.52);
		p[SHIFT]=Real(-0.87);
		p[SLIDE]=Real(-0.61);
		p[RISE]=Real(3.07);
		StepArray steps;
		steps.append_step_parameters(p);
		steps.append_step_parameters(p);
		steps.append_step_parameters(p);

        // output
        {
            OutputArchive<CerealBinaryOutput> bin_archive("test.cereal", false);
            bin_archive.save(steps);
        }

        // input
        InputArchive<CerealBinaryInput> bin_archive("test.cereal");
        StepArray read_steps;
        bin_archive.load(read_steps);

        EXPECT_EQ(steps.n_step(), read_steps.n_step());
        EXPECT_NEAR((steps.parameters(2).inline_vector()-
                     read_steps.parameters(2).inline_vector()).norm(),
                    FLOAT_INIT,
                    FLOAT_INIT);

    };


    // SequenceSerializationTest
    TEST_F(DNASimArchiveTest, SequenceSerializationTest) {

        StepSequence seq(BaseSymbol::A,BaseSymbol::T);
        Sequence full_seq("ACTGCTC");

        // output
        {
            OutputArchive<CerealBinaryOutput> bin_archive("test.cereal", false);
            bin_archive.save(seq);
            bin_archive.save(full_seq);
        }

        // input
        InputArchive<CerealBinaryInput> bin_archive("test.cereal");
        StepSequence read_seq;
        Sequence read_full_seq;
        bin_archive.load(read_seq);
        bin_archive.load(read_full_seq);

        EXPECT_EQ(seq.first_base(), read_seq.first_base());
        EXPECT_EQ(seq.last_base(), read_seq.last_base());
        EXPECT_EQ(full_seq.linear_sequence(),
                  read_full_seq.linear_sequence());

    };


}

