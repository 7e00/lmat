#ifndef NARUTOACM_MAT_H_
#define NARUTOACM_MAT_H_

namespace narutoacm
{
// every matrix like object should herited this template base class
template<typename derivedT>
class MatrixBase
{
public:
    inline const derivedT& Derived() const
    {
        return static_cast<const derivedT&>(*this);
    }

    /*
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
    */
};

// matrix, or matrix expression which can be access every element simplely
// e.g. a+b.t() is "Simple Matrix"
template<typename derivedT>
class SimpleMatrixBase : public MatrixBase<derivedT>
{
};

} // end of namespace narutoacm

#endif

