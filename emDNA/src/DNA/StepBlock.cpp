// StepBlock structure
// Nicolas Clauvelin


#include "DNA/StepBlock.h"


// constructor with initialization from four 3x3 matrices
StepBlock::StepBlock(const Matrix3& rot_rot,
                     const Matrix3& rot_trans,
                     const Matrix3& trans_rot,
                     const Matrix3& trans_trans) :
_rot_rot(rot_rot),
_rot_trans(rot_trans),
_trans_rot(trans_rot),
_trans_trans(trans_trans) {};


// constructor with initialization from a matrix
StepBlock::StepBlock(const MatrixN& step_matrix) :
_rot_rot(), _rot_trans(), _trans_rot(), _trans_trans() {

    // size checking
    DS_ASSERT(step_matrix.n_cols()==emDNAConstants::StepParametersDim&&
              step_matrix.n_rows()==emDNAConstants::StepParametersDim,
              "wrong-sized matrix to build a step block");

};


// full matrix method
MatrixN StepBlock::full_matrix() const {

    MatrixN full_mat(emDNAConstants::StepParametersDim);
    full_mat.set_block(0, 0, _rot_rot.block_matrix());
    full_mat.set_block(0, 3, _rot_trans.block_matrix());
    full_mat.set_block(3, 0, _trans_rot.block_matrix());
    full_mat.set_block(3, 3, _trans_trans.block_matrix());

    return full_mat;

};


// transpose method
void StepBlock::transpose() {

    // diagonal terms
    _rot_rot.transpose();
    _trans_trans.transpose();

    // off diagonal terms
    Matrix3 rtT = _rot_trans;
    rtT.transpose();
    Matrix3 trT = _trans_rot;
    trT.transpose();
    _rot_trans = trT;
    _trans_rot = rtT;

};
StepBlock StepBlock::get_transpose() const {

    return
    StepBlock(_rot_rot.get_transpose(), _trans_rot.get_transpose(),
              _rot_trans.get_transpose(), _trans_trans.get_transpose());

};


// multiplication operator
StepBlock& StepBlock::operator*=(const StepBlock& h) {

    StepBlock cp_this(*this);

    _rot_rot =
    cp_this._rot_rot*h._rot_rot+cp_this._rot_trans*h._trans_rot;
    _rot_trans =
    cp_this._rot_rot*h._rot_trans+cp_this._rot_trans*h._trans_trans;
    _trans_rot =
    cp_this._trans_rot*h._rot_rot+cp_this._trans_trans*h._trans_rot;
    _trans_trans =
    cp_this._trans_rot*h._rot_trans+cp_this._trans_trans*h._trans_trans;

    return *this;

};
StepBlock operator*(const StepBlock& h1, StepBlock h2) {

    StepBlock mul_h = h1;
    mul_h *= h2;

    return mul_h;

};


// addition operator
StepBlock& StepBlock::operator+=(const StepBlock& h) {

    StepBlock cp_this(*this);

    _rot_rot = cp_this._rot_rot + h._rot_rot;
    _rot_trans = cp_this._rot_trans + h._rot_trans;
    _trans_rot = cp_this._trans_rot + h._trans_rot;
    _trans_trans = cp_this._trans_trans + h._trans_trans;

    return *this;

};
StepBlock operator+(const StepBlock& h1, StepBlock h2) {

    StepBlock add_h = h1;
    add_h += h2;

    return add_h;

};


// constructor with size initialization
StepBlockArray2D::StepBlockArray2D(const Size& n_steps) :
m_n_steps(n_steps),
m_block_matrix(n_steps*n_steps, StepBlock()) {};


// accessor/modifier operator
const StepBlock& StepBlockArray2D::operator()(const Size& i, const Size& j)
const {
    return m_block_matrix[i*m_n_steps+j];
};
StepBlock& StepBlockArray2D::operator()(const Size& i, const Size& j) {
    return m_block_matrix[i*m_n_steps+j];
};


// regular matrix accessor
MatrixN StepBlockArray2D::regular_matrix() const {

    // container
    MatrixN mat(m_n_steps*emDNAConstants::StepParametersDim);

    // filling
    for (Size i=0; i<m_n_steps; ++i) {
        for (Size j=0; j<m_n_steps; ++j) {
            mat.set_block(i*emDNAConstants::StepParametersDim,
                          j*emDNAConstants::StepParametersDim,
                          m_block_matrix[i*m_n_steps+j].
                          full_matrix().block_matrix());
        };
    };

    return mat;

};


