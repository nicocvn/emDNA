// StepBlock structure
// Nicolas Clauvelin


#ifndef emDNA_StepBlock_h
#define emDNA_StepBlock_h


#include <emDNA_Includes.h>


// step hessian block structure
// the hessian is stored as the four dofs quadrants
struct StepBlock {

    Matrix3 _rot_rot;
    Matrix3 _rot_trans;
    Matrix3 _trans_rot;
    Matrix3 _trans_trans;
    StepBlock() = default;
    StepBlock(const Matrix3& rot_rot,
              const Matrix3& rot_trans,
              const Matrix3& trans_rot,
              const Matrix3& trans_trans);
    StepBlock(const MatrixN& step_matrix);
    ~StepBlock() = default;
    StepBlock(const StepBlock& sh) = default;
    StepBlock(StepBlock&& sh) = default;
    StepBlock& operator=(const StepBlock& sh) = default;
    StepBlock& operator=(StepBlock&& sh) = default;
    MatrixN full_matrix() const;
    friend std::ostream& operator<<(std::ostream& os,
                                    const StepBlock& h) {
        os << h.full_matrix();
        return os;
    };

    // transpose method
    void transpose();
    StepBlock get_transpose() const;

    // multiplication operator
    StepBlock& operator*=(const StepBlock& h);
    friend StepBlock operator*(const StepBlock& h1, StepBlock h2);

    // addition operator
    StepBlock& operator+=(const StepBlock& h);
    friend StepBlock operator+(const StepBlock& h1, StepBlock h2);

};
using StepBlockVec = std::vector<StepBlock>;


// quick and dirty implementation of a 2d array of step blocks
class StepBlockArray2D {


public:

    // constructors
    StepBlockArray2D() = default;
    StepBlockArray2D(const Size& n_steps);
    StepBlockArray2D(const StepBlockArray2D& array_2d) = default;
    StepBlockArray2D(StepBlockArray2D&& array_2d) = default;
    ~StepBlockArray2D() = default;

    // copy and move operators
    StepBlockArray2D& operator=(const StepBlockArray2D& array_2d) = default;
    StepBlockArray2D& operator=(StepBlockArray2D&& array_2d) = default;

    // accessor/modifier operator
    const StepBlock& operator()(const Size& i, const Size& j) const;
    StepBlock& operator()(const Size& i, const Size& j);

    // regular matrix accessor
    MatrixN regular_matrix() const;


private:

    Size m_n_steps;
    StepBlockVec m_block_matrix;


};


#endif  // emDNA_StepBlock_h
