#ifndef NARUTOACM_MAT_H_
#define NARUTOACM_MAT_H_

#include <assert.h>

namespace narutoacm
{
// every matrix like object should herited this template base class
template<typename eT, typename ExprT>
class MatrixBase
{
public:
    inline const ExprT& Derived() const
    {
        return static_cast<const ExprT&>(*this);
    }
    inline ExprT& Derived()
    {
        return static_cast<ExprT&>(*this);
    }

    inline int Elems() const
    {
        return this->Derived().Elems();
    }
    inline int Rows() const
    {
        return this->Derived().Rows();
    }
    inline int Cols() const
    {
        return this->Derived().Cols();
    }
    
    inline const eT operator()(int idx) const
    {
        return this->Derived()(idx);
    }
    inline const eT operator()(int r, int c) const
    {
        return this->Derived()(r, c);
    }
};

// matrix, or matrix expression which can be access every element simplely
// e.g. a+b.t() is "Simple Matrix"
template<typename eT, typename ExprT>
class SimpleMatrixBase : public MatrixBase<eT,ExprT>
{
};

// Matrix or MatrixShell
template <typename eT, typename ExprT>
class EntityMatrixBase : public SimpleMatrixBase<eT,ExprT>
{
public:
    inline int Capacity() const
    {
        return this->Derived().Capacity();
    }
    inline eT *Data()
    {
        return this->Derived().Data();
    }
    inline const eT *Data() const
    {
        return this->Derived().Data();
    }
    inline void SetShape(int row, int col)
    {
        this->Derived().SetShape(row, col);
    }
    inline void SetRows(int row)
    {
        this->Derived().SetRows(row);
    }
    inline void SetCols(int col)
    {
        this->Derived().SetCols(col);
    }
    inline eT &operator()(int idx)
    {
        return this->Derived()(idx);
    }
    inline const eT &operator()(int idx) const
    {
        return this->Derived()(idx);
    }
    inline eT &operator()(int r, int c)
    {
        return this->Derived()(r, c);
    }
    inline const eT &operator()(int r, int c) const
    {
        return this->Derived()(r, c);
    }
};


} // end of namespace narutoacm

#endif

