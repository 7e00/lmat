#ifndef NARUTOACM_MATSHELL_H_
#define NARUTOACM_MATSHELL_H_

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
class MatrixShell;

template <typename eT, int M, int N>
class MatrixShell : public EntityMatrixBase<eT,MatrixShell<eT,M,N>>
{
public:
    MatrixShell()
        : data_(nullptr)
    {
    }
    MatrixShell(eT *data, int cap)
        : data_(nullptr), capacity_(0)
    {
        SetData(data, cap);
    }
    MatrixShell(const MatrixShell &) = default;
    MatrixShell &operator=(const MatrixShell &) = default;

    void SetData(eT *data, int cap)
    {
        assert(cap >= M*N);
        data_ = data;
        capacity_ = cap;
    }

    template <typename expreT, typename ExprT>
    MatrixShell(const MatrixBase<expreT,ExprT> &matexpr);
    template <typename expreT, typename ExprT>
    MatrixShell &operator=(const MatrixBase<expreT, ExprT> &matexpr);

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
    eT *data_;
    int capacity_;
};

template <typename eT>
class MatrixShell<eT,0,0> : public EntityMatrixBase<eT,MatrixShell<eT,0,0>>
{
public:
    MatrixShell()
        : data_(nullptr), capacity_(0), row_(0), col_(0)
    {
    }
    MatrixShell(eT *data, int cap, int row = -1, int col = -1)
        : data_(nullptr), capacity_(0)
    {
        SetData(data, cap, row, col);
    }
    MatrixShell(const MatrixShell &mat) = default;
    MatrixShell &operator=(const MatrixShell &mat) = default;

    template <typename expreT, typename ExprT>
    MatrixShell(const MatrixBase<expreT, ExprT> &matexpr);
    template <typename expreT, typename ExprT>
    MatrixShell &operator=(const MatrixBase<expreT, ExprT> &matexpr);

    void SetData(eT *data, int cap, int row = -1, int col = -1)
    {
        assert(row * col <= cap);
        data_ = data;
        capacity_ = cap;
        row_ = (row <= 0 ? cap : row);
        col_ = (col <= 0 ? cap / row_ : col);
    }

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
        assert(cap <= capacity_);
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

private:
    eT *data_;
    int capacity_;
    int row_;
    int col_;
};

template <typename eT, int M>
class MatrixShell<eT,M,0> : public EntityMatrixBase<eT,MatrixShell<eT,M,0>>
{
public:
    MatrixShell()
        : data_(nullptr), capacity_(0), col_(0)
    {
    }
    MatrixShell(eT *data, int cap, int col = -1)
        : data_(nullptr), capacity_(0)
    {
        SetData(data, cap, col);
    }
    MatrixShell(const MatrixShell &mat) = default;
    MatrixShell &operator=(const MatrixShell &mat) = default;

    void SetData(eT *data, int cap, int col = -1)
    {
        assert(M * col <= cap);
        data_ = data;
        capacity_ = cap;
        col_ = (col <= 0 ? cap / M : col);
    }

    template <typename expreT, typename ExprT>
    MatrixShell(const MatrixBase<expreT,ExprT> &matexpr);
    template <typename expreT, typename ExprT>
    MatrixShell &operator=(const MatrixBase<expreT,ExprT> &matexpr);

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
        assert(cap <= capacity_);
        col_ = col;
    }

private:
    eT *data_;
    int capacity_;
    int col_;
};

template <typename eT, int N>
class MatrixShell<eT,0,N> : public EntityMatrixBase<eT,MatrixShell<eT,0,N>>
{
public:
    MatrixShell()
        : data_(nullptr), capacity_(0), row_(0)
    {
    }
    MatrixShell(eT *data, int cap, int row = -1)
        : data_(nullptr), capacity_(0)
    {
        SetData(data, cap, row);
    }
    MatrixShell(const MatrixShell &mat) = default;
    MatrixShell &operator=(const MatrixShell &mat) = default;

    void SetData(eT *data, int cap, int row = -1)
    {
        assert(row * N <= cap);
        data_ = data;
        capacity_ = cap;
        row_ = (row <= 0 ? cap / N : row);
    }

    template <typename expreT, typename ExprT>
    MatrixShell(const MatrixBase<expreT,ExprT> &matexpr);
    template <typename expreT, typename ExprT>
    MatrixShell &operator=(const MatrixBase<expreT,ExprT> &matexpr);

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
        assert(cap <= capacity_);
        row_ = row;
    }
    void SetCols(int col)
    {
        assert(N == col);
    }

private:
    eT *data_;
    int capacity_;
    int row_;
};

} // end of namespace narutoacm

#endif

