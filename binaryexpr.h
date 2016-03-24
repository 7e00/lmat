#ifndef NARUTOACM_BINARYEXPR_H_
#define NARUTOACM_BINARYEXPR_H_

#include <memory>
#include "matbase.h"
#include "matmeta.h"

namespace narutoacm
{

template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
class BinaryExpr;

template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2, int type>
class BinaryExprHelper;

template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
class BinaryExprHelper<eT,M,N,opT,ExprT1,ExprT2,0> : public MatrixBase<eT,BinaryExpr<eT,M,N,opT,ExprT1,ExprT2>>
{
protected:
    BinaryExprHelper() = default;

public:
    eT operator()(int idx) const
    {
        const_cast<BinaryExprHelper<eT,M,N,opT,ExprT1,ExprT2,0> *>(this)->calres();
        return (*res_)(idx);
    } 
    eT operator()(int r, int c) const
    {
        const_cast<BinaryExprHelper<eT,M,N,opT,ExprT1,ExprT2,0> *>(this)->calres();
        return (*res_)(r, c);
    }

private:
    void calres()
    {
        if (!res_)
        {
            res_ = std::unique_ptr<Matrix<eT,M,N>>(new Matrix<eT,M,N>());
            opT op(*res_, this->Derived());
        }
    }

protected:
    std::unique_ptr<Matrix<eT,M,N>> res_;
};

template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
class BinaryExprHelper<eT,M,N,opT,ExprT1,ExprT2,1> : public SimpleMatrixBase<eT,BinaryExpr<eT,M,N,opT,ExprT1,ExprT2>>
{
protected:
    BinaryExprHelper() = default;

public:
    eT operator()(int idx) const
    {
        return opT::eval(this->Derived(), idx);
    }
    eT operator()(int r, int c) const
    {
        return opT::eval(this->Derived(), r, c);
    }
};

using meta::helper::is_same_type;
using meta::helper::select_type;
using meta::element_type;
using meta::matrix_expr_type;

template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
class BinaryExpr : public BinaryExprHelper<eT,M,N,opT,ExprT1,ExprT2,matrix_expr_type<BinaryExpr<eT,M,N,opT,ExprT1,ExprT2>>::value>
{
    friend opT;
    friend class Op;
    template <typename feT, int fM, int fN, typename fopT, typename fexpreT1, typename fExprT1, typename fexpreT2, typename fExprT2>
    friend BinaryExpr<feT,fM,fN,fopT,fExprT1,fExprT2> GenBinaryExpr(const MatrixBase<fexpreT1,fExprT1> &expr1, const MatrixBase<fexpreT2,fExprT2> &expr2);

protected:
    BinaryExpr(const ExprT1 &expr1, const ExprT2 &expr2)
        : expr1_(expr1), expr2_(expr2)
    {
    }

public:
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
    

protected:
    typename select_type<is_same_type<ExprT1,Scalar<typename element_type<ExprT1>::type>>::result, const ExprT1, const ExprT1&>::result_type expr1_;
    typename select_type<is_same_type<ExprT2,Scalar<typename element_type<ExprT2>::type>>::result, const ExprT2, const ExprT2&>::result_type expr2_;
};

template <typename eT, int M, typename opT, typename ExprT1, typename ExprT2>
class BinaryExpr<eT,M,0,opT,ExprT1,ExprT2> : public BinaryExprHelper<eT,M,0,opT,ExprT1,ExprT2,matrix_expr_type<BinaryExpr<eT,M,0,opT,ExprT1,ExprT2>>::value>
{
    friend opT;
    friend class Op;
    template <typename feT, int fM, int fN, typename fopT, typename fexpreT1, typename fExprT1, typename fexpreT2, typename fExprT2>
    friend BinaryExpr<feT,fM,fN,fopT,fExprT1,fExprT2> GenBinaryExpr(const MatrixBase<fexpreT1,fExprT1> &expr1, const MatrixBase<fexpreT2,fExprT2> &expr2);

protected:
    BinaryExpr(const ExprT1 &expr1, const ExprT2 &expr2)
        : expr1_(expr1), expr2_(expr2)
    {
        col_ = opT::cols(*this);
    }

public:
    int Elems() const
    {
        return M*col_;
    }
    int Rows() const
    {
        return M;
    }
    int Cols() const
    {
        return col_;
    }

protected:
    typename select_type<is_same_type<ExprT1,Scalar<typename element_type<ExprT1>::type>>::result, const ExprT1, const ExprT1&>::result_type expr1_;
    typename select_type<is_same_type<ExprT2,Scalar<typename element_type<ExprT2>::type>>::result, const ExprT2, const ExprT2&>::result_type expr2_;
    int col_;
};

template <typename eT, int N, typename opT, typename ExprT1, typename ExprT2>
class BinaryExpr<eT,0,N,opT,ExprT1,ExprT2> : public BinaryExprHelper<eT,0,N,opT,ExprT1,ExprT2,matrix_expr_type<BinaryExpr<eT,0,N,opT,ExprT1,ExprT2>>::value>
{
    friend opT;
    friend class Op;
    template <typename feT, int fM, int fN, typename fopT, typename fexpreT1, typename fExprT1, typename fexpreT2, typename fExprT2>
    friend BinaryExpr<feT,fM,fN,fopT,fExprT1,fExprT2> GenBinaryExpr(const MatrixBase<fexpreT1,fExprT1> &expr1, const MatrixBase<fexpreT2,fExprT2> &expr2);

protected:
    BinaryExpr(const ExprT1 &expr1, const ExprT2 &expr2)
        : expr1_(expr1), expr2_(expr2)
    {
        row_ = opT::rows(*this);
    }

public:
    int Elems() const
    {
        return row_*N;
    }
    int Rows() const
    {
        return row_;
    }
    int Cols() const
    {
        return N;
    }

protected:
    typename select_type<is_same_type<ExprT1,Scalar<typename element_type<ExprT1>::type>>::result, const ExprT1, const ExprT1&>::result_type expr1_;
    typename select_type<is_same_type<ExprT2,Scalar<typename element_type<ExprT2>::type>>::result, const ExprT2, const ExprT2&>::result_type expr2_;
    int row_;
};

template <typename eT, typename opT, typename ExprT1, typename ExprT2>
class BinaryExpr<eT,0,0,opT,ExprT1,ExprT2> : public BinaryExprHelper<eT,0,0,opT,ExprT1,ExprT2,matrix_expr_type<BinaryExpr<eT,0,0,opT,ExprT1,ExprT2>>::value>
{
    friend opT;
    friend class Op;
    template <typename feT, int fM, int fN, typename fopT, typename fexpreT1, typename fExprT1, typename fexpreT2, typename fExprT2>
    friend BinaryExpr<feT,fM,fN,fopT,fExprT1,fExprT2> GenBinaryExpr(const MatrixBase<fexpreT1,fExprT1> &expr1, const MatrixBase<fexpreT2,fExprT2> &expr2);

protected:
    BinaryExpr(const ExprT1 &expr1, const ExprT2 &expr2)
        : expr1_(expr1), expr2_(expr2)
    {
        row_ = opT::rows(*this);
        col_ = opT::cols(*this);
    }

public:
    int Elems() const
    {
        return row_*col_;
    }
    int Rows() const
    {
        return row_;
    }
    int Cols() const
    {
        return col_;
    }

protected:
    typename select_type<is_same_type<ExprT1,Scalar<typename element_type<ExprT1>::type>>::result, const ExprT1, const ExprT1&>::result_type expr1_;
    typename select_type<is_same_type<ExprT2,Scalar<typename element_type<ExprT2>::type>>::result, const ExprT2, const ExprT2&>::result_type expr2_;
    int row_;
    int col_;
};

template <typename eT, int M, int N, typename opT, typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
BinaryExpr<eT,M,N,opT,ExprT1,ExprT2> GenBinaryExpr(const MatrixBase<expreT1,ExprT1> &expr1, const MatrixBase<expreT2,ExprT2> &expr2)
{
    return BinaryExpr<eT,M,N,opT,ExprT1,ExprT2>(expr1.Derived(), expr2.Derived());
}

} // end of namespace narutoacm
#endif

