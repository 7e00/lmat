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
class Matrix : public EntityMatrixBase<eT,Matrix<eT,M,N>>
{
public:
    Matrix() = default;
    Matrix(eT initval)
    {
        std::fill(data_, data_ + Elems(), initval);
    }
    Matrix(const Matrix &) = default;
    Matrix &operator=(const Matrix &) = default;

    template <typename expreT, typename ExprT>
    Matrix(const MatrixBase<expreT,ExprT> &matexpr);
    template <typename expreT, typename ExprT>
    Matrix &operator=(const MatrixBase<expreT, ExprT> &matexpr);

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
    eT *Data()
    {
        return data_;
    }
    const eT *Data() const
    {
        return data_;
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
        return data_[r*Cols()+c];
#else
        return data_[r+c*Rows()];
#endif
    }
    const eT &operator()(int r, int c) const
    {
#if MATRIXORDER == ROWORDER
        return data_[r*Cols()+c];
#else
        return data_[r+c*Rows()];
#endif
    }

    void SetShape(int row, int col)
    {
        assert(M == row && N == col);
    }
    void SetRows(int row)
    {
        assert(M == row);
    }
    void SetCols(int col)
    {
        assert(N == col);
    }

private:
    eT data_[M*N];
};

template <typename eT>
class Matrix<eT,0,0> : public EntityMatrixBase<eT,Matrix<eT,0,0>>
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
        mat.dump();
    }
    ~Matrix()
    {
        dump();
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
        mat.dump();
    }

    template <typename expreT, typename ExprT>
    Matrix(const MatrixBase<expreT, ExprT> &matexpr);
    template <typename expreT, typename ExprT>
    Matrix &operator=(const MatrixBase<expreT, ExprT> &matexpr);

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
    eT *Data()
    {
        return data_;
    }
    const eT *Data() const
    {
        return data_;
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
        return data_[r*Cols()+c];
#else
        return data_[r+c*Rows()];
#endif
    }
    const eT &operator()(int r, int c) const
    {
#if MATRIXORDER == ROWORDER
        return data_[r*Cols()+c];
#else
        return data_[r+c*Rows()];
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


protected:
    void dump()
    {
        delete [] data_;
        data_ = nullptr;
        capacity_ = 0;
    }

private:
    eT *data_;
    int capacity_;
    int row_;
    int col_;
};

template <typename eT, int M>
class Matrix<eT,M,0> : public EntityMatrixBase<eT,Matrix<eT,M,0>>
{
public:
    Matrix()
        : data_(nullptr), capacity_(0), col_(0)
    {
    }
    Matrix(int col)
        : data_(nullptr)
    {
        capacity_ = M * col;
        data_ = new eT[capacity_];
        col_ = col;
    }
    Matrix(int col, eT initval)
        : data_(nullptr)
    {
        capacity_ = M * col;
        data_ = new eT[capacity_];
        col_ = col;
        std::fill(data_, data_ + capacity_, initval);
    }
    Matrix(const Matrix &mat)
        : data_(nullptr)
    {
        capacity_ = M * mat.col_;
        data_ = new eT[capacity_];
        memcpy(data_, mat.data_, capacity_ * sizeof(eT));
        col_ = mat.col_;
    }
    Matrix(Matrix &&mat)
        : data_(mat.data_), capacity_(mat.capacity_), col_(mat.col_)
    {
        mat.dump();
    }
    ~Matrix()
    {
        dump();
    }
    Matrix &operator=(const Matrix &mat)
    {
        if (this != &mat)
        {
            SetCols(mat.col_);
            memcpy(data_, mat.data_, this->Elems() * sizeof(eT));
        }
    }
    Matrix &operator=(Matrix &&mat)
    {
        dump();
        data_ = mat.data_;
        capacity_ = mat.capacity_;
        col_ = mat.col_;
        mat.dump();
    }

    template <typename expreT, typename ExprT>
    Matrix(const MatrixBase<expreT,ExprT> &matexpr);
    template <typename expreT, typename ExprT>
    Matrix &operator=(const MatrixBase<expreT,ExprT> &matexpr);

    int Elems() const
    {
        return M * col_;
    }
    int Rows() const
    {
        return M;
    }
    int Cols() const
    {
        return col_;
    }
    int Capacity() const
    {
        return capacity_;
    }
    eT *Data()
    {
        return data_;
    }
    const eT *Data() const
    {
        return data_;
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
        return data_[r*Cols()+c];
#else
        return data_[r+c*Rows()];
#endif
    }
    const eT &operator()(int r, int c) const
    {
#if MATRIXORDER == ROWORDER
        return data_[r*Cols()+c];
#else
        return data_[r+c*Rows()];
#endif
    }

    void SetShape(int row, int col)
    {
        assert(M == row);
        SetCols(col);
    }
    void SetRows(int row)
    {
        assert(M == row);
    }
    void SetCols(int col)
    {
        auto cap = M * col;
        if (cap > capacity_)
        {
            dump();
            capacity_ = cap;
            data_ = new eT[capacity_];
        }
        col_ = col;
    }

protected:
    void dump()
    {
        delete [] data_;
        data_ = nullptr;
        capacity_ = 0;
    }

private:
    eT *data_;
    int capacity_;
    int col_;
};

template <typename eT, int N>
class Matrix<eT,0,N> : public EntityMatrixBase<eT,Matrix<eT,0,N>>
{
public:
    Matrix()
        : data_(nullptr), capacity_(0), row_(0)
    {
    }
    Matrix(int row)
        : data_(nullptr)
    {
        capacity_ = row * N;
        data_ = new eT[capacity_];
        row_ = row;
    }
    Matrix(int row, eT initval)
        : data_(nullptr)
    {
        capacity_ = row * N;
        data_ = new eT[capacity_];
        row_ = row;
        std::fill(data_, data_ + capacity_, initval);
    }
    Matrix(const Matrix &mat)
        : data_(nullptr)
    {
        capacity_ = mat.row_ * N;
        data_ = new eT[capacity_];
        memcpy(data_, mat.data_, capacity_ * sizeof(eT));
        row_ = mat.row_;
    }
    Matrix(Matrix &&mat)
        : data_(mat.data_), capacity_(mat.capacity_), row_(mat.row_)
    {
        mat.dump();
    }
    ~Matrix()
    {
        dump();
    }
    Matrix &operator=(const Matrix &mat)
    {
        if (this != &mat)
        {
            SetRows(mat.row_);
            memcpy(data_, mat.data_, this->Elems() * sizeof(eT));
        }
    }
    Matrix &operator=(Matrix &&mat)
    {
        dump();
        data_ = mat.data_;
        capacity_ = mat.capacity_;
        row_ = mat.row_;
        mat.dump();
    }

    template <typename expreT, typename ExprT>
    Matrix(const MatrixBase<expreT,ExprT> &matexpr);
    template <typename expreT, typename ExprT>
    Matrix &operator=(const MatrixBase<expreT,ExprT> &matexpr);

    int Elems() const
    {
        return row_ * N;
    }
    int Rows() const
    {
        return row_;
    }
    int Cols() const
    {
        return N;
    }
    int Capacity() const
    {
        return capacity_;
    }
    eT *Data()
    {
        return data_;
    }
    const eT *Data() const
    {
        return data_;
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
        return data_[r*Cols()+c];
#else
        return data_[r+c*Rows()];
#endif
    }
    const eT &operator()(int r, int c) const
    {
#if MATRIXORDER == ROWORDER
        return data_[r*Cols()+c];
#else
        return data_[r+c*Rows()];
#endif
    }

    void SetShape(int row, int col)
    {
        assert(N == col);
        SetRows(row);
    }
    void SetRows(int row)
    {
        auto cap = row * N;
        if (cap > capacity_)
        {
            dump();
            capacity_ = cap;
            data_ = new eT[capacity_];
        }
        row_ = row;
    }
    void SetCols(int col)
    {
        assert(N == col);
    }

protected:
    void dump()
    {
        delete [] data_;
        data_ = nullptr;
        capacity_ = 0;
    }

private:
    eT *data_;
    int capacity_;
    int row_;
};

} // end of namespace narutoacm

#endif

