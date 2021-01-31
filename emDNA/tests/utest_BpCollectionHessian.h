// BpCollectionHessian unit tests
// Nicolas Clauvelin


#ifndef emDNA_utest_BpCollectionHessian_h
#define emDNA_utest_BpCollectionHessian_h


#include <gtest/gtest.h>
#include <BpStepParams.h>


namespace {


    // testing class for BpCollectionHessianTest
    class BpCollectionHessianTest : public ::testing::Test {

    protected:

        BpCollectionHessianTest() {};
        virtual ~BpCollectionHessianTest() {};

        // set up before test
        // called after the constructor
        virtual void SetUp() {

            // step parameters
            p1 = VectorN({
                4.645880454195819, 8.39030823575434, 35.615435525650994,
                -0.024477835167500128, -0.018101457881504412, 3.846066481936653
            });
            p2 = VectorN({
                -9.904259329614288, 8.347924523680547, 37.80333336368872,
                0.2705645954420093, 0.165139557061885, 3.623383368039957
            });
            p3 = VectorN({
                -0.2117433900087775, 8.53375613491388, 43.052642685792605,
                0.6647932714787816, -0.4549898068199303, 2.508112987277481
            });
            test_params = {
                BpStepParams(p1),
                BpStepParams(p2),
                BpStepParams(p3)
            };

            // single step Hessian from mathematica
            // for the first step (diagonal terms)
            Real mm_single_step_H[36];
            mm_single_step_H[0] = 612.8393602981668;
            mm_single_step_H[1] = -13.439898275174178;
            mm_single_step_H[2] = -61.854316201144556;
            mm_single_step_H[3] = -10.3481636741989;
            mm_single_step_H[4] = 32.32577797376071;
            mm_single_step_H[5] = 1.3768517390203745;
            mm_single_step_H[6] = -13.439898275174066;
            mm_single_step_H[7] = 600.0110129603108;
            mm_single_step_H[8] = -81.90419343521025;
            mm_single_step_H[9] = -32.26513145538375;
            mm_single_step_H[10] = -10.399306318639917;
            mm_single_step_H[11] = 2.4865492341461835;
            mm_single_step_H[12] = -61.854316201144556;
            mm_single_step_H[13] = -81.90419343521025;
            mm_single_step_H[14] = 226.33874665107288;
            mm_single_step_H[15] = -0.5504389783688967;
            mm_single_step_H[16] = -2.78848789421231;
            mm_single_step_H[17] = FLOAT_INIT;
            mm_single_step_H[18] = -10.348163674198899;
            mm_single_step_H[19] = -32.26513145538375;
            mm_single_step_H[20] = -0.5504389783688969;
            mm_single_step_H[21] = 19.999999999999996;
            mm_single_step_H[22] = FLOAT_INIT;
            mm_single_step_H[23] = FLOAT_INIT;
            mm_single_step_H[24] = 32.32577797376071;
            mm_single_step_H[25] = -10.399306318639916;
            mm_single_step_H[26] = -2.7884878942123104;
            mm_single_step_H[27] = FLOAT_INIT;
            mm_single_step_H[28] = 19.999999999999996;
            mm_single_step_H[29] = FLOAT_INIT;
            mm_single_step_H[30] = 1.3768517390203745;
            mm_single_step_H[31] = 2.4865492341461835;
            mm_single_step_H[32] = FLOAT_INIT;
            mm_single_step_H[33] = FLOAT_INIT;
            mm_single_step_H[34] = FLOAT_INIT;
            mm_single_step_H[35] = 19.999999999999996;
            mm_single_step_Hessian.resize(6, 6);
            mm_single_step_Hessian << Array<Real>(36, mm_single_step_H);

            // full Hessian results for complete test
            // including sequence dependence
            Real free_h_elems[324];
            free_h_elems[0] = 723.9950689946296;
            free_h_elems[1] = 52.572727177016205;
            free_h_elems[2] = -35.830041748951224;
            free_h_elems[3] = -12.553866751252919;
            free_h_elems[4] = 11.477548467314989;
            free_h_elems[5] = 4.583516179691038;
            free_h_elems[6] = 78.66505047775863;
            free_h_elems[7] = -54.44343871429379;
            free_h_elems[8] = -8.560586655606773;
            free_h_elems[9] = -20.589887335298243;
            free_h_elems[10] = 19.890325427237723;
            free_h_elems[11] = -2.5780032155333563;
            free_h_elems[12] = -14.68337244110298;
            free_h_elems[13] = -29.325916999927458;
            free_h_elems[14] = 11.872393114110757;
            free_h_elems[15] = -10.571764852572858;
            free_h_elems[16] = 15.134707977724915;
            free_h_elems[17] = -17.2545189108269;
            free_h_elems[18] = 52.572727177016255;
            free_h_elems[19] = 334.1317548110185;
            free_h_elems[20] = 42.30685007463355;
            free_h_elems[21] = 0.4950620719902217;
            free_h_elems[22] = -3.372656570960172;
            free_h_elems[23] = -1.7417379036769456;
            free_h_elems[24] = 114.39728575432017;
            free_h_elems[25] = 78.44824454241811;
            free_h_elems[26] = -49.18572062171813;
            free_h_elems[27] = -20.784258734721426;
            free_h_elems[28] = 7.05597930828675;
            free_h_elems[29] = 4.158446571129871;
            free_h_elems[30] = 51.28890665096421;
            free_h_elems[31] = 3.1175439826641265;
            free_h_elems[32] = -35.586008186738894;
            free_h_elems[33] = -33.68382468567331;
            free_h_elems[34] = -3.9054542521216766;
            free_h_elems[35] = 16.840693833415777;
            free_h_elems[36] = -35.830041748951224;
            free_h_elems[37] = 42.306850074633545;
            free_h_elems[38] = 254.34291733711208;
            free_h_elems[39] = -2.0758996193370693;
            free_h_elems[40] = -6.493447186870016;
            free_h_elems[41] = -23.688056681056114;
            free_h_elems[42] = -57.77652932427928;
            free_h_elems[43] = -24.18756174060944;
            free_h_elems[44] = 19.92375541782394;
            free_h_elems[45] = 2.9127053252055157;
            free_h_elems[46] = -1.9668926856142126;
            free_h_elems[47] = -0.20733738853759848;
            free_h_elems[48] = -21.515926334149405;
            free_h_elems[49] = 9.217750955834923;
            free_h_elems[50] = 14.304892304092526;
            free_h_elems[51] = 15.84662692868695;
            free_h_elems[52] = 2.208134129745611;
            free_h_elems[53] = -4.351390652199881;
            free_h_elems[54] = -12.553866751252922;
            free_h_elems[55] = 0.4950620719902217;
            free_h_elems[56] = -2.07589961933707;
            free_h_elems[57] = 3.994090484025335;
            free_h_elems[58] = -2.166819843152103;
            free_h_elems[59] = 1.590008038491086;
            free_h_elems[60] = 0.;
            free_h_elems[61] = 0.;
            free_h_elems[62] = 0.;
            free_h_elems[63] = 0.;
            free_h_elems[64] = 0.;
            free_h_elems[65] = 0.;
            free_h_elems[66] = 0.;
            free_h_elems[67] = 0.;
            free_h_elems[68] = 0.;
            free_h_elems[69] = 0.;
            free_h_elems[70] = 0.;
            free_h_elems[71] = 0.;
            free_h_elems[72] = 11.477548467314985;
            free_h_elems[73] = -3.372656570960172;
            free_h_elems[74] = -6.493447186870014;
            free_h_elems[75] = -2.166819843152103;
            free_h_elems[76] = 9.976636084005502;
            free_h_elems[77] = 0.3905097152221981;
            free_h_elems[78] = 0.;
            free_h_elems[79] = 0.;
            free_h_elems[80] = 0.;
            free_h_elems[81] = 0.;
            free_h_elems[82] = 0.;
            free_h_elems[83] = 0.;
            free_h_elems[84] = 0.;
            free_h_elems[85] = 0.;
            free_h_elems[86] = 0.;
            free_h_elems[87] = 0.;
            free_h_elems[88] = 0.;
            free_h_elems[89] = 0.;
            free_h_elems[90] = 4.583516179691038;
            free_h_elems[91] = -1.7417379036769456;
            free_h_elems[92] = -23.688056681056118;
            free_h_elems[93] = 1.5900080384910855;
            free_h_elems[94] = 0.39050971522219813;
            free_h_elems[95] = 25.44227343196916;
            free_h_elems[96] = 0.;
            free_h_elems[97] = 0.;
            free_h_elems[98] = 0.;
            free_h_elems[99] = 0.;
            free_h_elems[100] = 0.;
            free_h_elems[101] = 0.;
            free_h_elems[102] = 0.;
            free_h_elems[103] = 0.;
            free_h_elems[104] = 0.;
            free_h_elems[105] = 0.;
            free_h_elems[106] = 0.;
            free_h_elems[107] = 0.;
            free_h_elems[108] = 78.66505047775863;
            free_h_elems[109] = 114.39728575432017;
            free_h_elems[110] = -57.77652932427928;
            free_h_elems[111] = 0.;
            free_h_elems[112] = 0.;
            free_h_elems[113] = 0.;
            free_h_elems[114] = 318.6981958035228;
            free_h_elems[115] = 33.98034750873533;
            free_h_elems[116] = -17.573908408833052;
            free_h_elems[117] = -2.0686358002820717;
            free_h_elems[118] = -6.898335618833821;
            free_h_elems[119] = 27.806254433348126;
            free_h_elems[120] = -7.959021525947669;
            free_h_elems[121] = -25.82988688807746;
            free_h_elems[122] = -12.273673250849512;
            free_h_elems[123] = -13.891549650395504;
            free_h_elems[124] = 12.768428470617359;
            free_h_elems[125] = -7.09713346508701;
            free_h_elems[126] = -54.44343871429379;
            free_h_elems[127] = 78.44824454241811;
            free_h_elems[128] = -24.18756174060944;
            free_h_elems[129] = 0.;
            free_h_elems[130] = 0.;
            free_h_elems[131] = 0.;
            free_h_elems[132] = 33.98034750873532;
            free_h_elems[133] = 256.75156231830704;
            free_h_elems[134] = 28.51949309282005;
            free_h_elems[135] = -2.6691021064969966;
            free_h_elems[136] = -5.9666812820214155;
            free_h_elems[137] = 23.62912334463799;
            free_h_elems[138] = 53.00561038005279;
            free_h_elems[139] = 48.4956768881012;
            free_h_elems[140] = -18.030109411027013;
            free_h_elems[141] = -10.89117153830384;
            free_h_elems[142] = -21.215996085728897;
            free_h_elems[143] = 22.18534799733751;
            free_h_elems[144] = -8.560586655606773;
            free_h_elems[145] = -49.18572062171813;
            free_h_elems[146] = 19.92375541782394;
            free_h_elems[147] = 0.;
            free_h_elems[148] = 0.;
            free_h_elems[149] = 0.;
            free_h_elems[150] = -17.57390840883304;
            free_h_elems[151] = 28.51949309282006;
            free_h_elems[152] = 235.67895858219836;
            free_h_elems[153] = -2.2622510746620588;
            free_h_elems[154] = -6.511039186024259;
            free_h_elems[155] = -22.82316975180133;
            free_h_elems[156] = -5.603523744265056;
            free_h_elems[157] = 3.617934590127035;
            free_h_elems[158] = 8.969013141392116;
            free_h_elems[159] = 3.056255798863724;
            free_h_elems[160] = 3.132285828919436;
            free_h_elems[161] = -2.3438687453258162;
            free_h_elems[162] = -20.589887335298243;
            free_h_elems[163] = -20.784258734721426;
            free_h_elems[164] = 2.9127053252055157;
            free_h_elems[165] = 0.;
            free_h_elems[166] = 0.;
            free_h_elems[167] = 0.;
            free_h_elems[168] = -2.068635800282071;
            free_h_elems[169] = -2.6691021064969958;
            free_h_elems[170] = -2.2622510746620588;
            free_h_elems[171] = 5.986746838202515;
            free_h_elems[172] = 2.499955059906813;
            free_h_elems[173] = 2.244388168716383;
            free_h_elems[174] = 0.;
            free_h_elems[175] = 0.;
            free_h_elems[176] = 0.;
            free_h_elems[177] = 0.;
            free_h_elems[178] = 0.;
            free_h_elems[179] = 0.;
            free_h_elems[180] = 19.890325427237723;
            free_h_elems[181] = 7.05597930828675;
            free_h_elems[182] = -1.9668926856142126;
            free_h_elems[183] = 0.;
            free_h_elems[184] = 0.;
            free_h_elems[185] = 0.;
            free_h_elems[186] = -6.898335618833821;
            free_h_elems[187] = -5.966681282021414;
            free_h_elems[188] = -6.51103918602426;
            free_h_elems[189] = 2.4999550599068137;
            free_h_elems[190] = 3.7070439954547934;
            free_h_elems[191] = 0.8699589726841307;
            free_h_elems[192] = 0.;
            free_h_elems[193] = 0.;
            free_h_elems[194] = 0.;
            free_h_elems[195] = 0.;
            free_h_elems[196] = 0.;
            free_h_elems[197] = 0.;
            free_h_elems[198] = -2.5780032155333563;
            free_h_elems[199] = 4.158446571129871;
            free_h_elems[200] = -0.20733738853759848;
            free_h_elems[201] = 0.;
            free_h_elems[202] = 0.;
            free_h_elems[203] = 0.;
            free_h_elems[204] = 27.806254433348126;
            free_h_elems[205] = 23.62912334463799;
            free_h_elems[206] = -22.823169751801323;
            free_h_elems[207] = 2.244388168716383;
            free_h_elems[208] = 0.8699589726841308;
            free_h_elems[209] = 22.448209166342686;
            free_h_elems[210] = 0.;
            free_h_elems[211] = 0.;
            free_h_elems[212] = 0.;
            free_h_elems[213] = 0.;
            free_h_elems[214] = 0.;
            free_h_elems[215] = 0.;
            free_h_elems[216] = -14.68337244110298;
            free_h_elems[217] = 51.28890665096421;
            free_h_elems[218] = -21.515926334149405;
            free_h_elems[219] = 0.;
            free_h_elems[220] = 0.;
            free_h_elems[221] = 0.;
            free_h_elems[222] = -7.959021525947669;
            free_h_elems[223] = 53.00561038005279;
            free_h_elems[224] = -5.603523744265056;
            free_h_elems[225] = 0.;
            free_h_elems[226] = 0.;
            free_h_elems[227] = 0.;
            free_h_elems[228] = 236.55153587041846;
            free_h_elems[229] = 10.742590673498132;
            free_h_elems[230] = -0.8091337938933123;
            free_h_elems[231] = -9.153657446784322;
            free_h_elems[232] = -8.034149626195846;
            free_h_elems[233] = 4.7343410817730405;
            free_h_elems[234] = -29.325916999927458;
            free_h_elems[235] = 3.1175439826641265;
            free_h_elems[236] = 9.217750955834923;
            free_h_elems[237] = 0.;
            free_h_elems[238] = 0.;
            free_h_elems[239] = 0.;
            free_h_elems[240] = -25.82988688807746;
            free_h_elems[241] = 48.4956768881012;
            free_h_elems[242] = 3.617934590127035;
            free_h_elems[243] = 0.;
            free_h_elems[244] = 0.;
            free_h_elems[245] = 0.;
            free_h_elems[246] = 10.742590673498118;
            free_h_elems[247] = 175.0455503438168;
            free_h_elems[248] = 77.34036572182082;
            free_h_elems[249] = 0.9226971595391072;
            free_h_elems[250] = -7.13886398969299;
            free_h_elems[251] = 4.961969399707001;
            free_h_elems[252] = 11.872393114110757;
            free_h_elems[253] = -35.586008186738894;
            free_h_elems[254] = 14.304892304092526;
            free_h_elems[255] = 0.;
            free_h_elems[256] = 0.;
            free_h_elems[257] = 0.;
            free_h_elems[258] = -12.273673250849512;
            free_h_elems[259] = -18.030109411027013;
            free_h_elems[260] = 8.969013141392116;
            free_h_elems[261] = 0.;
            free_h_elems[262] = 0.;
            free_h_elems[263] = 0.;
            free_h_elems[264] = -0.8091337938933139;
            free_h_elems[265] = 77.34036572182082;
            free_h_elems[266] = 152.12456940303315;
            free_h_elems[267] = -2.932218867493803;
            free_h_elems[268] = -1.1810119900651301;
            free_h_elems[269] = -9.247338224367944;
            free_h_elems[270] = -10.571764852572858;
            free_h_elems[271] = -33.68382468567331;
            free_h_elems[272] = 15.84662692868695;
            free_h_elems[273] = 0.;
            free_h_elems[274] = 0.;
            free_h_elems[275] = 0.;
            free_h_elems[276] = -13.891549650395504;
            free_h_elems[277] = -10.89117153830384;
            free_h_elems[278] = 3.056255798863724;
            free_h_elems[279] = 0.;
            free_h_elems[280] = 0.;
            free_h_elems[281] = 0.;
            free_h_elems[282] = -9.153657446784322;
            free_h_elems[283] = 0.9226971595391072;
            free_h_elems[284] = -2.932218867493803;
            free_h_elems[285] = 3.3375717276371963;
            free_h_elems[286] = 0.4116403528666656;
            free_h_elems[287] = 0.7565867268251006;
            free_h_elems[288] = 15.134707977724915;
            free_h_elems[289] = -3.9054542521216766;
            free_h_elems[290] = 2.208134129745611;
            free_h_elems[291] = 0.;
            free_h_elems[292] = 0.;
            free_h_elems[293] = 0.;
            free_h_elems[294] = 12.768428470617359;
            free_h_elems[295] = -21.215996085728897;
            free_h_elems[296] = 3.132285828919436;
            free_h_elems[297] = 0.;
            free_h_elems[298] = 0.;
            free_h_elems[299] = 0.;
            free_h_elems[300] = -8.034149626195848;
            free_h_elems[301] = -7.13886398969299;
            free_h_elems[302] = -1.1810119900651301;
            free_h_elems[303] = 0.41164035286666567;
            free_h_elems[304] = 2.425631918162899;
            free_h_elems[305] = 3.08476136023078;
            free_h_elems[306] = -17.2545189108269;
            free_h_elems[307] = 16.840693833415777;
            free_h_elems[308] = -4.351390652199881;
            free_h_elems[309] = 0.;
            free_h_elems[310] = 0.;
            free_h_elems[311] = 0.;
            free_h_elems[312] = -7.09713346508701;
            free_h_elems[313] = 22.18534799733751;
            free_h_elems[314] = -2.3438687453258162;
            free_h_elems[315] = 0.;
            free_h_elems[316] = 0.;
            free_h_elems[317] = 0.;
            free_h_elems[318] = 4.7343410817730405;
            free_h_elems[319] = 4.961969399707002;
            free_h_elems[320] = -9.247338224367944;
            free_h_elems[321] = 0.7565867268251008;
            free_h_elems[322] = 3.08476136023078;
            free_h_elems[323] = 13.287796354199905;
            mm_free_Hessian.resize(3*6,3*6);
            mm_free_Hessian << Array<Real>(3*6*3*6, free_h_elems);

        };

        // clean up before test
        // called after the destructor
        virtual void TearDown() {};

        // declare objects that will be used by all tests
        VectorN p1, p2, p3;
        std::vector<BpStepParams> test_params;
        MatrixN mm_single_step_Hessian, mm_free_Hessian;
        
    };
    
    
}


#endif  // emDNA_utest_BpCollectionHessian_h

