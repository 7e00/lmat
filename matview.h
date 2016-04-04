#ifndef NARUTOACM_MATVIEW_H_
#define NARUTOACM_MATVIEW_H_

#include <assert.h>
#include "matbase.h"
#include "matmeta.h"

namespace narutoacm
{

template <typename eT, int M, int N, typename ExprT>
class MatrixView;

template <typename eT, int M, int N, typename ExprT, int type>
class MatrixViewHelper;

using meta::element_type;

template <typename eT, int M, int N, typename ExprT>
class MatrixViewHelper<eT,M,N,ExprT,0> : public MatrixBase<eT,MatrixView<eT,M,N,ExprT>>
{
protected:
    MatrixViewHelper(const ExprT &expr, int sr, int er, int rstep, int sc, int ec, int cstep)
        : expr_(expr), startrow_(sr), startcol_(sc), rowstep_(rstep), colstep_(cstep)
    {
        row_ = (er - sr + rstep - 1) / rstep;
        col_ = (ec - sc + cstep - 1) / cstep;
        assert((M == 0 || M == row_) && (N == 0 || N == col_));
    }

public:
    eT operator()(int idx) const
    {
#if MATRIXORDER == ROWORDER
        return this->operator()(idx/this->col_, idx%this->col_);
#else
        return this->operator()(idx%this->row_, idx/this->row_);
#endif
    }
    eT operator()(int r, int c) const
    {
        return expr_(startrow_+r*rowstep_, startcol_+c*colstep_);
    }
    
protected:
    const ExprT &expr_;
    int startrow_;
    int startcol_;
    int rowstep_;
    int colstep_;
    int row_;
    int col_;
};

template <typename eT, int M, int N, typename ExprT>
class MatrixViewHelper<eT,M,N,ExprT,1> : public SimpleMatrixBase<eT,MatrixView<eT,M,N,ExprT>>
{
protected:
    MatrixViewHelper(const ExprT &expr, int sr, int er, int rstep, int sc, int ec, int cstep)
        : expr_(expr), startrow_(sr), startcol_(sc), rowstep_(rstep), colstep_(cstep)
    {
        row_ = (er - sr + rstep - 1) / rstep;
        col_ = (ec - sc + cstep - 1) / cstep;
        assert((M == 0 || M == row_) && (N == 0 || N == col_));
    }
    
public:
    eT operator()(int idx) const
    {
#if MATRIXORDER == ROWORDER
        return this->operator()(idx/this->col_, idx%this->col_);
#else
        return this->operator()(idx%this->row_, idx/this->row_);
#endif
    }
    eT operator()(int r, int c) const
    {
        return expr_(startrow_+r*rowstep_, startcol_+c*colstep_);
    }
    
protected:
    const ExprT &expr_;
    int startrow_;
    int startcol_;
    int rowstep_;
    int colstep_;
    int row_;
    int col_;
};

template <typename eT, int M, int N, typename ExprT>
class MatrixViewHelper<eT,M,N,ExprT,2> : public MutableMatrixBase<eT,MatrixView<eT,M,N,ExprT>>
{
protected:
    MatrixViewHelper(ExprT &expr, int sr, int er, int rstep, int sc, int ec, int cstep)
        : expr_(expr), startrow_(sr), startcol_(sc), rowstep_(rstep), colstep_(cstep)
    {
        row_ = (er - sr + rstep - 1) / rstep;
        col_ = (ec - sc + cstep - 1) / cstep;
        assert((M == 0 || M == row_) && (N == 0 || N == col_));
    }
    
public:
    eT &operator()(int idx)
    {
#if MATRIXORDER == ROWORDER
        return this->operator()(idx/this->col_, idx%this->col_);
#else
        return this->operator()(idx%this->row_, idx/this->row_);
#endif
    }
    const eT &operator()(int idx) const
    {
#if MATRIXORDER == ROWORDER
        return this->operator()(idx/this->col_, idx%this->col_);
#else
        return this->operator()(idx%this->row_, idx/this->row_);
#endif
    }
    eT &operator()(int r, int c)
    {
        return expr_(startrow_+r*rowstep_, startcol_+c*colstep_);
    }
    const eT &operator()(int r, int c) const
    {
        return expr_(startrow_+r*rowstep_, startcol_+c*colstep_);
    }
    
protected:
    ExprT &expr_;
    int startrow_;
    int startcol_;
    int rowstep_;
    int colstep_;
    int row_;
    int col_;
};

using meta::matrix_expr_type;
using meta::helper::select_type;

template <typename eT, int M, int N, typename ExprT>
class MatrixView : public MatrixViewHelper<eT,M,N,ExprT,matrix_expr_type<MatrixView<eT,M,N,ExprT>>::value>
{
    friend class MatOp;

protected:
    MatrixView(typename select_type<matrix_expr_type<MatrixView<eT,M,N,ExprT>>::value == 2, ExprT, const ExprT>::result_type &expr, int sr, int er, int rstep, int sc, int ec, int cstep)
        : MatrixViewHelper<eT,M,N,ExprT,matrix_expr_type<MatrixView<eT,M,N,ExprT>>::value>(expr,sr,er,rstep,sc,ec,cstep)
    {
    }

public:
    MatrixView &operator=(const MatrixView &matexpr);
    template <typename expreT, typename ExprT1>
    MatrixView &operator=(const MatrixBase<expreT, ExprT1> &matexpr);
    
    int Elems() const
    {
        return this->row_*this->col_;
    }
    int Rows() const
    {
        return this->row_;
    }
    int Cols() const
    {
        return this->col_;
    }
};



} // end of namespace narutoacm

#endif
