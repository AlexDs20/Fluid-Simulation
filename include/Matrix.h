#ifndef MATRIX_H
#define MATRIX_H
using namespace std;

#define Real float

class Matrix{
  private:
    int width, height;
    Real inVal;
    Real* val;

  public:
    Matrix(int width, int height);
    ~Matrix(void);

    int idx(int i, int j);
    void setVal(Real inVal);
    Real getVal(int i, int j);
};

#endif
