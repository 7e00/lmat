#ifndef NARUTOACM_MATRIX_H_
#define NARUTOACM_MATRIX_H_

#include <algorithm>
#include <assert.h>
#include "matbase.h"

#ifndef ROWORDER
#define ROWORDER 0
#endif
#ifndef COLORDER
#define COLORDER 1
#endif
#ifndef MATRIXORDER
#define MATRIXORDER COLORDER
#endif

namespace narutoacm
{

template <typename eT, int M, int N>
class Matrix;

template <typename eT, int M, int N>
class Matrix : public SimpleMatrixBase<Matrix<eT,M,N>>
{
public:
    Matrix() = default;
    Matrix(eT initval)
    {
        std::fill(data_, data_ + Elems(), initval);
    }
    Matrix(const Matrix &) = default;
    Matrix &operator=(const Matrix &) = default;

    template <typename ExprT>
    Matrix(const MatrixBase<ExprT> &matexpr);
    template <typename ExprT>
    Matrix &operator=(const MatrixBase<ExprT> &matexpr);

    int Elems() const
    {
        return M*N;
    }
    int Rows() const
    {
        return M;
    }
    int Cols() const
    {
        return N;
    }
    int Capacity() const
    {
        return Elems();
    }

    eT &operator()(int idx)
    {
        return data_[idx];
    }
    const eT &operator()(int idx) const
    {
        return data_[idx];
    }
    eT &operator()(int r, int c)
    {
#if MATRIXORDER == ROWORDER
        return data_[r*N+c];
#else
        return data_[r+c*M];
#endif
    }
    const eT &operator()(int r, int c) const
    {
#if MATRIXORDER == ROWORDER
        return data_[r*N+c];
#else
        return data_[r+c*M];
#endif
    }

    void SetShape(int r, int c)
    {
        assert(M == r && N == c);
    }

private:
    eT data_[M*N];
};

template <typename eT>
class Matrix<eT,0,0> : public SimpleMatrixBase<Matrix<eT,0,0>>
{
public:
    Matrix()
        : data_(nullptr), capacity_(0), row_(0), col_(0)
    {
    }
    Matrix(int row, int col)
        : data_(nullptr)
    {
        capacity_ = row * col;
        data_ = new eT[capacity_];
        row_ = row;
        col_ = col;
    }
    Matrix(int row, int col, eT initval)
        : data_(nullptr)
    {
        capacity_ = row * col;
        data_ = new eT[capacity_];
        row_ = row;
        col_ = col;
        std::fill(data_, data_ + capacity_, initval);
    }
    Matrix(const Matrix &mat)
        : data_(nullptr)
    {
        capacity_ = mat.row_ * mat.col_;
        data_ = new eT[capacity_];
        memcpy(data_, mat.data_, capacity_ * sizeof(eT));
        row_ = mat.row_;
        col_ = mat.col_;
    }
    Matrix(Matrix &&mat)
        : data_(mat.data_), capacity_(mat.capacity_), row_(mat.row_), col_(mat.col_)
    {
        mat.data_ = nullptr;
        mat.capacity_ = 0;
    }
    Matrix &operator=(const Matrix &mat)
    {
        if (this != &mat)
        {
            SetShape(mat.row_, mat.col_);
            memcpy(data_, mat.data_, this->Elems() * sizeof(eT));
        }
    }
    Matrix &operator=(Matrix &&mat)
    {
        dump();
        data_ = mat.data_;
        capacity_ = mat.capacity_;
        row_ = mat.row_;
        col_ = mat.col_;
        mat.data_ = nullptr;
        mat.capacity_ = 0;
    }

    template <typename ExprT>
    Matrix(const MatrixBase<ExprT> &matexpr);
    template <typename ExprT>
    Matrix &operator=(const MatrixBase<ExprT> &matexpr);

    int Elems() const
    {
        return row_ * col_;
    }
    int Rows() const
    {
        return row_;
    }
    int Cols() const
    {
        return col_;
    }
    int Capacity() const
    {
        return capacity_;
    }

    eT &operator()(int idx)
    {
        return data_[idx];
    }
    const eT &operator()(int idx) const
    {
        return data_[idx];
    }
    eT &operator()(int r, int c)
    {
#if MATRIXORDER == ROWORDER
        return data_[r*col_+c];
#else
        return data_[r+c*row_];
#endif
    }
    const eT &operator()(int r, int c) const
    {
#if MATRIXORDER == ROWORDER
        return data_[r*col_+c];
#else
        return data_[r+c*row_];
#endif
    }

    void SetShape(int row, int col)
    {
        auto cap = row * col;
        if (cap > capacity_)
        {
            dump();
            capacity_ = cap;
            data_ = new eT[capacity_];
        }
        row_ = row;
        col_ = col;
    }
    void SetRows(int row)
    {
        SetShape(row, col_);
    }
    void SetCols(int col)
    {
        SetShape(row_, col);
    }

    eT *Data()
    {
        return data_;
    }
    const eT *Data() const
    {
        return data_;
    }

protected:
    void dump()
    {
        delete [] data_;
    }

private:
    eT *data_;
    int capacity_;
    int row_;
    int col_;
};

} // end of namespace narutoacm

#endif

