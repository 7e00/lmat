#ifndef NARUTOACM_BINARYEXPR_H_
#define NARUTOACM_BINARYEXPR_H_

#include <memory>
#include "matbase.h"
#include "matmeta.h"

namespace narutoacm
{
namespace meta
{
template <typename ExprT>
struct expr_type
{
    typedef const ExprT& result_type;
};
template <typename eT>
struct expr_type<Scalar<eT>>
{
    typedef const Scalar<eT> result_type;
};
} // end of namespace narutoacm::meta

template <typename eT, typename opT, typename ExprT1, typename ExprT2>
class BinaryExpr;

template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
class BinaryExprBase
{
protected:
    BinaryExprBase(const ExprT1 &expr1, const ExprT2 &expr2)
        :expr1_(expr1), expr2_(expr2)
    {
    }

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
    typename meta::expr_type<ExprT1>::result_type expr1_;
    typename meta::expr_type<ExprT2>::result_type expr2_;
};

template <typename eT, typename opT, typename ExprT1, typename ExprT2>
class BinaryExprBase<eT, 0, 0, opT, ExprT1, ExprT2>
{
protected:
    BinaryExprBase(const ExprT1 &expr1, const ExprT2 &expr2)
        :expr1_(expr1), expr2_(expr2)
    {
        row_ = opT::rows(*reinterpret_cast<BinaryExpr<eT,opT,ExprT1,ExprT2> *>(this));
        col_ = opT::cols(*reinterpret_cast<BinaryExpr<eT,opT,ExprT1,ExprT2> *>(this));
    }

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
    typename meta::expr_type<ExprT1>::result_type expr1_;
    typename meta::expr_type<ExprT2>::result_type expr2_;
    int row_;
    int col_;
};
template <typename eT, int M, typename opT, typename ExprT1, typename ExprT2>
class BinaryExprBase<eT, M, 0, opT, ExprT1, ExprT2>
{
protected:
    BinaryExprBase(const ExprT1 &expr1, const ExprT2 &expr2)
        :expr1_(expr1), expr2_(expr2)
    {
        col_ = opT::cols(static_cast<BinaryExpr<eT,opT,ExprT1,ExprT2>>(*this));
    }

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
    typename meta::expr_type<ExprT1>::result_type expr1_;
    typename meta::expr_type<ExprT2>::result_type expr2_;
    int col_;
};
template <typename eT, int N, typename opT, typename ExprT1, typename ExprT2>
class BinaryExprBase<eT, 0, N, opT, ExprT1, ExprT2>
{
protected:
    BinaryExprBase(const ExprT1 &expr1, const ExprT2 &expr2)
        :expr1_(expr1), expr2_(expr2)
    {
        row_ = opT::rows(static_cast<BinaryExpr<eT,opT,ExprT1,ExprT2>>(*this));
    }

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
    typename meta::expr_type<ExprT1>::result_type expr1_;
    typename meta::expr_type<ExprT2>::result_type expr2_;
    int row_;
};


template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2, int type>
class BinaryExprHelper;

template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
class BinaryExprHelper<eT, M, N, opT, ExprT1, ExprT2, 0> : public BinaryExprBase<eT,M,N,opT,ExprT1,ExprT2>,
                                                     public MatrixBase<BinaryExpr<eT,opT,ExprT1,ExprT2>>
{
protected:
    BinaryExprHelper(const ExprT1 &expr1, const ExprT2 &expr2)
        : BinaryExprBase<eT,M,N,opT,ExprT1,ExprT2>(expr1, expr2)
    {
    }

public:
    int Elems() const
    {
        return BinaryExprBase<eT,M,N,opT,ExprT1,ExprT2>::Elems();
    }
    int Rows() const
    {
        return BinaryExprBase<eT,M,N,opT,ExprT1,ExprT2>::Rows();
    }
    int Cols() const
    {
        return BinaryExprBase<eT,M,N,opT,ExprT1,ExprT2>::Cols();
    }
    eT operator()(int idx) const
    {
        calres();
        return (*res_)(idx);
    } 
    eT operator()(int r, int c) const
    {
        calres();
        return (*res_)(r, c);
    }

private:
    void calres() const
    {
        if (!res_)
        {
            res_ = std::unique_ptr<Matrix<eT,0,0>>(new Matrix<eT,0,0>());
            opT op(*res_, *this);
        }
    }

protected:
    std::unique_ptr<Matrix<eT,0,0>> res_;
};

template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
class BinaryExprHelper<eT, M, N, opT, ExprT1, ExprT2, 1> : public BinaryExprBase<eT,M,N,opT,ExprT1,ExprT2>,
                                                     public SimpleMatrixBase<BinaryExpr<eT,opT,ExprT1,ExprT2>>
{
protected:
    BinaryExprHelper(const ExprT1 &expr1, const ExprT2 &expr2)
        : BinaryExprBase<eT,M,N,opT,ExprT1,ExprT2>(expr1, expr2)
    {
    }

public:
    int Elems() const
    {
        return BinaryExprBase<eT,M,N,opT,ExprT1,ExprT2>::Elems();
    }
    int Rows() const
    {
        return BinaryExprBase<eT,M,N,opT,ExprT1,ExprT2>::Rows();
    }
    int Cols() const
    {
        return BinaryExprBase<eT,M,N,opT,ExprT1,ExprT2>::Cols();
    }
    eT operator()(int idx) const
    {
        return opT::eval(*this, idx);
    } 
    eT operator()(int r, int c) const
    {
        return opT::eval(this, r, c);
    }

};

template <typename eT, typename opT, typename ExprT1, typename ExprT2>
class BinaryExpr : public BinaryExprHelper<eT, meta::rows<BinaryExpr<eT,opT,ExprT1,ExprT2>>::value, meta::cols<BinaryExpr<eT,opT,ExprT1,ExprT2>>::value, opT, ExprT1, ExprT2, meta::matrix_expr_type<BinaryExpr<eT,opT,ExprT1,ExprT2>>::value>
{
    friend opT;
    friend class Op;
protected:
    BinaryExpr(const ExprT1 &expr1, const ExprT2 &expr2)
        : BinaryExprHelper<eT,meta::rows<BinaryExpr<eT,opT,ExprT1,ExprT2>>::value,meta::cols<BinaryExpr<eT,opT,ExprT1,ExprT2>>::value,opT,ExprT1,ExprT2,meta::matrix_expr_type<BinaryExpr<eT,opT,ExprT1,ExprT2>>::value>(expr1, expr2)
    {
    }
};

} // end of namespace narutoacm
#endif

