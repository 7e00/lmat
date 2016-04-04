#ifndef NARUTOACM_MATOP_H_
#define NARUTOACM_MATOP_H_

#include <iostream>
#include <algorithm>
#include <vector>
#include <assert.h>
#include "scalar.h"
#include "matshell.h"
#include "matrix.h"
#include "binaryexpr.h"
#include "unaryexpr.h"
#include "matview.h"
#include "sharedpool.h"

namespace narutoacm
{
using meta::element_type;
using meta::helper::assert;
using meta::rows;
using meta::cols;
using util::SharedPool;
//using std::min;
//using std::max;

namespace
{
int exgcd(int a, int b, int &x, int &y)
{
    if(b == 0)
    {
        x = 1;
        y = 0;
        return a;
    }
    int d = exgcd(b, a%b, y, x);
    y = y - a / b * x;
    return d;
}
} // end of namespace narutoacm::

class MatOp
{
    template <typename eT, int M, int N>
    friend class Matrix;
    template <typename eT, int M, int N>
    friend class MatrixShell;
    template <typename eT, int M, int N, typename ExprT>
    friend class MatrixView;

protected:
    class Cache
    {
    public:
        struct Block
        {
            Block(unsigned long sz)
                : ptr(nullptr), size(sz)
            {
                ptr = new unsigned char[sz];
            }
            ~Block()
            {
                delete [] ptr;
            }
            unsigned char *ptr;
            unsigned long size;
        };
        Cache(unsigned long blocksize)
            : blocksize_(blocksize)
        {
        }
        
        using ptr_type = SharedPool<Block>::ptr_type;
        
        ptr_type Acquire()
        {
            if (pool_.Empty())
            {
                std::cout << "add a block with size " << blocksize_ << std::endl;
                pool_.Add(std::unique_ptr<Block>(new Block(blocksize_)));
            }
            return pool_.Acquire();
        }

        ptr_type GetResult()
        {
            return std::move(result_);
        }

        void StoreResult(ptr_type &a)
        {
            result_.swap(a);
        }
        
        size_t Size() const
        {
            return pool_.Size();
        }
        
    private:
        unsigned long blocksize_;
        SharedPool<Block> pool_;
        ptr_type result_;
    };

protected:
    template <typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static bool is_ref_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const EntityMatrixBase<eT2,ExprT2> &expr)
    {
        return (void *)(mat.Data()) == (void *)(expr.Data());
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename opT, typename ExprT2, typename ExprT3>
    static bool is_ref_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,opT,ExprT2,ExprT3> &expr)
    {
        return is_ref_to_mat(mat, expr.expr1_) || is_ref_to_mat(mat, expr.expr2_);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename SeT>
    static bool is_ref_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,MatOpSMul,ExprT2,Scalar<SeT>> &expr)
    {
        return is_ref_to_mat(mat, expr.expr1_);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename SeT>
    static bool is_ref_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,MatOpSMul,Scalar<SeT>,ExprT2> &expr)
    {
        return is_ref_to_mat(mat, expr.expr2_);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2>
    static bool is_ref_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const MatrixView<eT2,M,N,ExprT2> &expr)
    {
        return is_ref_to_mat(mat, expr.expr_);
    }

    // 0: no ref, 1: ref but ordered to mat, 2: ref and not ordered to mat
    template <typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const EntityMatrixBase<eT2,ExprT2> &expr)
    {
        return (is_ref_to_mat(mat, expr) ? (sizeof(eT1)<=sizeof(eT2)?1:2) : 0);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename opT, typename ExprT2, typename ExprT3>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,opT,ExprT2,ExprT3> &expr)
    {
        int ar = relation_to_mat(mat, expr.expr1_);
        return (ar == 2 ? ar : std::max(ar, relation_to_mat(mat, expr.expr2_)));
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename SeT>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,MatOpSMul,ExprT2,Scalar<SeT>> &expr)
    {
        return relation_to_mat(mat, expr.expr1_);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename SeT>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,MatOpSMul,Scalar<SeT>,ExprT2> &expr)
    {
        return relation_to_mat(mat, expr.expr2_);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename ExprT3>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,MatOpMul,ExprT2,ExprT3> &expr)
    {
        return (is_ref_to_mat(mat, expr) ? 2 : 0);
    }
    template <typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const MatrixBase<eT2,ExprT2> &expr,
                                int row, int col, int sr, int rstep, int sc, int cstep)
    {
        int r = relation_to_mat(mat, expr.Derived());
        if (r != 1)
            return r;
        if (rstep > 0 && cstep > 0)
            return 1;
        return 2;
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const MatrixView<eT2,M,N,ExprT2> &expr,
                                int row, int col, int sr, int rstep, int sc, int cstep)
    {
        return relation_to_mat(mat, expr.expr_, row, col, expr.startrow_ + sr * expr.rowstep_, rstep * expr.rowstep_,
                                expr.startcol_ + sc * expr.colstep_, cstep * expr.colstep_);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const MatrixView<eT2,M,N,ExprT2> &expr)
    {
        return relation_to_mat(mat, expr.expr_, expr.row_, expr.col_, expr.startrow_, expr.rowstep_, expr.startcol_, expr.colstep_);
    }

    template <typename eT1, typename eT2, typename ExprT1, typename ExprT2>
    static int relation_to_mat_part(const EntityMatrixBase<eT1,ExprT1> &mat, const MatrixBase<eT2,ExprT2> &expr, int row, int col,
                                    int sr1, int rstep1, int sc1, int cstep1, int sr2, int rstep2, int sc2, int cstep2)
    {
        int r = relation_to_mat(mat, expr.Derived());
        if (r != 1)
            return r;
        if (sizeof(eT1) < sizeof(eT2))
        {
            sr1 /= (sizeof(eT2)/sizeof(eT1));
            sc1 /= (sizeof(eT2)/sizeof(eT1));
            rstep1 /= (sizeof(eT2)/sizeof(eT1));
            cstep1 /= (sizeof(eT2)/sizeof(eT1));
        }
        else if (sizeof(eT2) < sizeof(eT1))
        {
            sr2 /= (sizeof(eT1)/sizeof(eT2));
            sc2 /= (sizeof(eT1)/sizeof(eT2));
            rstep2 /= (sizeof(eT1)/sizeof(eT2));
            cstep2 /= (sizeof(eT1)/sizeof(eT2));
        }

        // don't cross each other
        if ((sr1 - sr2 + (row-1)*rstep1) < 0 || (sr2 - sr1 + (row-1)*rstep2) < 0)
            return 1;
        if ((sc1 - sc2 + (col-1)*cstep1 < 0) || (sc2 - sc1 + (col-1)*cstep2) < 0)
            return 1;

        // whole overlap
        if (sr1 == sr2 && sc1 == sc2 && rstep1 == rstep2 && cstep1 == cstep2)
            return 1;

        // solve sr1+a*rstep1 = sr2+b*rstep2 && sc1+c*cstep1 = sc2+d*cstep2
        int r1, r2, c1, c2;
        int rd = exgcd(rstep1, rstep2, r1, r2);
        if ((sr2-sr1) % rd)
            return 1;
        int cd = exgcd(cstep1, cstep2, c1, c2);
        if ((sc2-sc1) % cd)
            return 1;

        sr1 /= rd;
        sr2 /= rd;
        rstep1 /= rd;
        rstep2 /= rd;
        r1 *= (sr2-sr1);
        r1 = (r1 % rstep2 + rstep2) % rstep2;
        if (r1 >= row)
            return 1;
        int flag = false;
        {
            int l = 0;
            int r = (row-1-r1)/rstep2;
            while (l <= r)
            {
                int mid = (l+r)>>1;
                int val = r1 + mid * rstep2;
                r2 = (val * rstep1 + sr1 - sr2) / rstep2;
                if (r2 >= 0 && r2 < row)
                {
                    flag = true;
                    break;
                }
                else if (r2 >= row)
                    r = mid - 1;
                else
                    l = mid + 1;
            }
        }
        if (!flag)
            return 1;

        sc1 /= cd;
        sc2 /= cd;
        cstep1 /= cd;
        cstep2 /= cd;
        c1 *= (sc2 - sc1);
        c1 = (c1 % cstep2 + cstep2) % cstep2;
        if (c1 >= col)
            return 1;
        flag = false;
        {
            int l = 0;
            int r = (col-1-c1)/cstep2;
            while (l <= r)
            {
                int mid = (l+r)>>1;
                int val = c1 + mid * cstep2;
                c2 = (val * cstep1 + sc1 - sc2) / cstep2;
                if (c2 >= 0 && c2 < col)
                {
                    flag = true;
                    break;
                }
                else if (c2 >= col)
                    r = mid - 1;
                else
                    l = mid + 1;
            }
        }
        return (flag ? 2 : 1);
    }
    template <typename eT1, typename eT2, typename ExprT1, typename ExprT2>
    static int relation_to_mat_part(const EntityMatrixBase<eT1,ExprT1> &mat, const MatrixBase<eT2,ExprT2> &expr, int row, int col,
                                    int sr1, int rstep1, int sc1, int cstep1)
    {
        return relation_to_mat_part(mat, expr.Derived(), row, col, sr1, rstep1, sc1, cstep1, 0, 1, 0, 1);
    }
    template <typename eT1, typename eT2, typename ExprT1, int M, int N, typename ExprT2>
    static int relation_to_mat_part(const EntityMatrixBase<eT1,ExprT1> &mat, const MatrixView<eT2,M,N,ExprT2> &expr, int row, int col,
                                    int sr1, int rstep1, int sc1, int cstep1, int sr2, int rstep2, int sc2, int cstep2)
    {
        return relation_to_mat_part(mat, expr.expr_, row, col, sr1, rstep1, sc1, cstep1,
                                    expr.startrow_ + sr2 * expr.rowstep_,
                                    expr.rowstep_ * rstep2,
                                    expr.startcol_ + sc2 * expr.colstep_,
                                    expr.colstep_ * cstep2);
    }
    template <typename eT1, int M, int N, typename ExprT1, typename eT2, typename ExprT2>
    static int relation_to_mat_part(const MutableMatrixBase<eT1,MatrixView<eT1,M,N,ExprT1>> &mat, const MatrixBase<eT2,ExprT2> &expr, int row, int col,
                                    int sr1, int rstep1, int sc1, int cstep1, int sr2, int rstep2, int sc2, int cstep2)
    {
        auto &m = mat.Derived();
        return relation_to_mat_part(m.expr_, expr.Derived(), row, col,
                                    m.startrow_ + sr1 * m.rowstep_,
                                    m.rowstep_ * rstep1,
                                    m.startcol_ + sc1 * m.colstep_,
                                    m.colstep_ * cstep1,
                                    sr2, rstep2, sc2, cstep2);
    }
    template <typename eT1, int M, int N, typename ExprT1, typename eT2, typename ExprT2>
    static int relation_to_mat_part(const MutableMatrixBase<eT1,MatrixView<eT1,M,N,ExprT1>> &mat, const MatrixBase<eT2,ExprT2> &expr)
    {
        auto &m = mat.Derived();
        return relation_to_mat_part(m.expr_, expr.Derived(), m.row_, m.col_, m.startrow_, m.rowstep_, m.startcol_, m.colstep_);
    }
    

    // relation < 0 means that we don't know the relation between res and expr
    template <typename eT, typename ExprT1, typename ExprT2>
    static void eval(EntityMatrixBase<eT,ExprT1> &res, const EntityMatrixBase<eT,ExprT2> &expr, Cache &cache, int relation = -1)
    {
        assert<(rows<ExprT1>::value == 0 || rows<ExprT2>::value == 0 || rows<ExprT1>::value == rows<ExprT2>::value) 
            && (cols<ExprT1>::value == 0 || cols<ExprT2>::value == 0 || cols<ExprT1>::value == cols<ExprT2>::value)>();
        
        relation = (relation >= 0 ? relation : relation_to_mat(res, expr));
        if (relation == 0)
        {
            res.SetShape(expr.Rows(), expr.Cols());
            memcpy(res.Data(), expr.Data(), expr.Elems() * sizeof(eT));
        }
    }
    template <typename eT1, typename eT2, typename ExprT1, typename ExprT2>
    static void eval(EntityMatrixBase<eT1,ExprT1> &res, const EntityMatrixBase<eT2,ExprT2> &expr, Cache &cache, int relation = -1)
    {
        assert<(rows<ExprT1>::value == 0 || rows<ExprT2>::value == 0 || rows<ExprT1>::value == rows<ExprT2>::value) 
            && (cols<ExprT1>::value == 0 || cols<ExprT2>::value == 0 || cols<ExprT1>::value == cols<ExprT2>::value)>();

        relation = (relation >= 0 ? relation : relation_to_mat(res, expr));

        if (relation == 2)
        {
            auto block = cache.Acquire();
            MatrixShell<eT2,rows<ExprT2>::value,cols<ExprT2>::value> ms((eT2 *)(block->ptr), block->size/sizeof(eT2));
            eval(ms, expr, cache, 0);
            eval(res, ms, cache, 0);
        }
        else
        {
            res.SetShape(expr.Rows(), expr.Cols());
            for (int i = 0; i < expr.Elems(); ++i)
                res(i) = expr(i);
        }
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename opT, typename ExprT2, typename ExprT3>
    static void eval(EntityMatrixBase<eT1,ExprT1> &res, const BinaryExpr<eT2,M,N,opT,ExprT2,ExprT3> &expr, Cache &cache, int relation = -1)
    {
        assert<(rows<ExprT1>::value == 0 || M == 0 || rows<ExprT1>::value == M) 
        && (cols<ExprT1>::value == 0 || N == 0 || cols<ExprT1>::value == N)>();
        
        opT::eval(res, expr, cache, relation);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2>
    static void eval(EntityMatrixBase<eT1,ExprT1> &res, const SimpleMatrixBase<eT2,MatrixView<eT2,M,N,ExprT2>> &expr, Cache &cache, int relation = -1)
    {
        assert<(rows<ExprT1>::value == 0 || M == 0 || rows<ExprT1>::value == M) 
        && (cols<ExprT1>::value == 0 || N == 0 || cols<ExprT1>::value == N)>();
        
        relation = relation = (relation >= 0 ? relation : relation_to_mat(res, expr.Derived()));

        auto &view = expr.Derived();
        
        if (relation == 2)
        {
            typedef typename element_type<ExprT2>::type T;
            eval<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>(view.expr_, cache);
            auto block = cache.GetResult();
            MatrixShell<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms((T *)(block->ptr), block->size/sizeof(T));
            ms.SetShape(view.expr_.Rows(), view.expr_.Cols());
            MatrixView<eT2,M,N,MatrixShell<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>> msview(ms,
                                                            view.startrow_, view.startrow_+view.row_*view.rowstep_, 
                                                            view.startcol_, view.startcol_+view.col_*view.colstep_,
                                                            view.rowstep_, view.colstep_);
            eval(res, msview, cache, 0);
        }
        else
        {
            res.SetShape(expr.Rows(), expr.Cols());
            for (int i = 0; i < expr.Elems(); ++i)
                res(i) = expr(i);
        }
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2>
    static void eval(EntityMatrixBase<eT1,ExprT1> &res, const MatrixBase<eT2,MatrixView<eT2,M,N,ExprT2>> &expr, Cache &cache, int relation = -1)
    {
        assert<(rows<ExprT1>::value == 0 || M == 0 || rows<ExprT1>::value == M) 
        && (cols<ExprT1>::value == 0 || N == 0 || cols<ExprT1>::value == N)>();
        
        auto &view = expr.Derived();

        typedef typename element_type<ExprT2>::type T;
        eval<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>(view.expr_, cache);
        auto block = cache.GetResult();
        MatrixShell<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms((T *)(block->ptr), block->size/sizeof(T));
        ms.SetShape(view.expr_.Rows(), view.expr_.Cols());
        MatrixView<eT2,M,N,MatrixShell<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>> msview(ms,
                                                        view.startrow_, view.startrow_+view.row_*view.rowstep_, 
                                                        view.startcol_, view.startcol_+view.col_*view.colstep_,
                                                        view.rowstep_, view.colstep_);
        eval(res, msview, cache, 0);
        
    }
    template <typename eT1, int M, int N, typename ExprT1, typename eT2, typename ExprT2>
    static void eval(MutableMatrixBase<eT1,MatrixView<eT1,M,N,ExprT1>> &res, const SimpleMatrixBase<eT2,ExprT2> &expr, Cache &cache, int relation = -1)
    {
        assert<(rows<ExprT2>::value == 0 || M == 0 || rows<ExprT2>::value == M) 
        && (cols<ExprT2>::value == 0 || N == 0 || cols<ExprT2>::value == N)>();

        auto &view = res.Derived();

        relation = (relation >= 0  ? relation : relation_to_mat_part(view, expr.Derived()));

        if (relation == 2)
        {
            eval<eT2,rows<ExprT2>::value,cols<ExprT2>::value>(expr.Derived(), cache);
            auto block = cache.GetResult();
            MatrixShell<eT2,rows<ExprT2>::value,cols<ExprT2>::value> ms((eT2 *)(block->ptr), block->size/sizeof(eT2));
            ms.SetShape(expr.Rows(), expr.Cols());
            eval(view, ms, cache, 0);
        }
        else
        {
            for (int i = 0; i < expr.Elems(); ++i)
                view(i) = expr(i);
        }
    }
    template <typename eT1, int M, int N, typename ExprT1, typename eT2, typename ExprT2>
    static void eval(MutableMatrixBase<eT1,MatrixView<eT1,M,N,ExprT1>> &res, const MatrixBase<eT2,ExprT2> &expr, Cache &cache, int relation = -1)
    {
        assert<(rows<ExprT2>::value == 0 || M == 0 || rows<ExprT2>::value == M) 
        && (cols<ExprT2>::value == 0 || N == 0 || cols<ExprT2>::value == N)>();

        auto &view = res.Derived();

        eval<eT2,rows<ExprT2>::value,cols<ExprT2>::value>(expr.Derived(), cache);
        auto block = cache.GetResult();
        MatrixShell<eT2,rows<ExprT2>::value,cols<ExprT2>::value> ms((eT2 *)(block->ptr), block->size/sizeof(eT2));
        ms.SetShape(expr.Rows(), expr.Cols());
        eval(view, ms, cache, 0);
    }
    template <typename eT1, int M1, int N1, typename ExprT1, typename eT2, int M2, int N2, typename ExprT2>
    static void eval(MutableMatrixBase<eT1,MatrixView<eT1,M1,N1,ExprT1>> &res, const SimpleMatrixBase<eT2,MatrixView<eT2,M2,N2,ExprT2>> &expr, Cache &cache, int relation = -1)
    {
        assert<(M1 == 0 || M2 == 0 || M1 == M2) 
        && (N1 == 0 || N2 == 0 || N1== N2)>();

        auto &view = res.Derived();
        auto &exprview = expr.Derived();

        relation = (relation >= 0  ? relation : relation_to_mat_part(view, exprview));

        if (relation == 2)
        {
            typedef typename element_type<ExprT2>::type T;
            eval<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>(exprview.expr_, cache);
            auto block = cache.GetResult();
            MatrixShell<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms((T *)(block->ptr), block->size/sizeof(T));
            ms.SetShape(exprview.expr_.Rows(), exprview.expr_.Cols());
            MatrixView<eT2,M2,N2,MatrixShell<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>> msview(ms,
                                                            exprview.startrow_, exprview.startrow_+exprview.row_*exprview.rowstep_, 
                                                            exprview.startcol_, exprview.startcol_+exprview.col_*exprview.colstep_,
                                                            exprview.rowstep_, exprview.colstep_);
            eval(view, msview, cache, 0);
        }
        else
        {
            for (int i = 0; i < expr.Elems(); ++i)
                view(i) = exprview(i);
        }
    }
    template <typename eT1, int M1, int N1, typename ExprT1, typename eT2, int M2, int N2, typename ExprT2>
    static void eval(MutableMatrixBase<eT1,MatrixView<eT1,M1,N1,ExprT1>> &res, const MatrixBase<eT2,MatrixView<eT2,M2,N2,ExprT2>> &expr, Cache &cache, int relation = -1)
    {
        assert<(M1 == 0 || M2 == 0 || M1 == M2) 
        && (N1 == 0 || N2 == 0 || N1== N2)>();

        auto &view = res.Derived();
        auto &exprview = expr.Derived();

        typedef typename element_type<ExprT2>::type T;
        eval<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>(exprview.expr_, cache);
        auto block = cache.GetResult();
        MatrixShell<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms((T *)(block->ptr), block->size/sizeof(T));
        ms.SetShape(exprview.expr_.Rows(), exprview.expr_.Cols());
        MatrixView<eT2,M2,N2,MatrixShell<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>> msview(ms,
                                                        exprview.startrow_, exprview.startrow_+exprview.row_*exprview.rowstep_, 
                                                        exprview.startcol_, exprview.startcol_+exprview.col_*exprview.colstep_,
                                                        exprview.rowstep_, exprview.colstep_);
        eval(view, msview, cache, 0);
    }
    template <typename eT, int M, int N, typename eT1, typename ExprT1>
    static void eval(const EntityMatrixBase<eT1,ExprT1> &expr, Cache &cache)
    {
        assert<(rows<ExprT1>::value == 0 || M == 0 || rows<ExprT1>::value == M) 
        && (cols<ExprT1>::value == 0 || N == 0 || cols<ExprT1>::value == N)>();

        auto block = cache.Acquire();
        MatrixShell<eT,M,N> ms((eT *)(block->ptr), block->size/sizeof(eT));
        eval(ms, expr, cache, 0);
        cache.StoreResult(block);
    }
    template <typename eT, int M, int N, typename eT1, int M1, int N1, typename opT, typename ExprT1, typename ExprT2>
    static void eval(const BinaryExpr<eT1,M1,N1,opT,ExprT1,ExprT2> &expr, Cache &cache)
    {
        assert<(M == 0 || M1 == 0 || M == M1) 
        && (N == 0 || N1 == 0 || N== N1)>();

        return opT::template eval<eT,M,N>(expr, cache);
    }
    template <typename eT, int M, int N, typename eT1, int M1, int N1, typename ExprT>
    static void eval(const SimpleMatrixBase<eT1,MatrixView<eT1,M1,N1,ExprT>> &expr, Cache &cache)
    {
        assert<(M == 0 || M1 == 0 || M == M1) 
        && (N == 0 || N1 == 0 || N== N1)>();

        auto block = cache.Acquire();
        MatrixShell<eT,M,N> ms((eT *)(block->ptr), block->size/sizeof(eT));
        eval(ms, expr.Derived(), cache, 0);
        cache.StoreResult(block);
    }
    template <typename eT, int M, int N, typename eT1, int M1, int N1, typename ExprT>
    static void eval(const MatrixBase<eT1,MatrixView<eT1,M1,N1,ExprT>> &expr, Cache &cache)
    {
        assert<(M == 0 || M1 == 0 || M == M1) 
        && (N == 0 || N1 == 0 || N== N1)>();

        auto &exprview = expr.Derived();

        typedef typename element_type<ExprT>::type T;
        eval<T,rows<ExprT>::value,cols<ExprT>::value>(exprview.expr_, cache);
        auto block = cache.GetResult();
        MatrixShell<T,rows<ExprT>::value,cols<ExprT>::value> ms1((T *)(block->ptr), block->size/sizeof(T));
        ms1.SetShape(exprview.expr_.Rows(), exprview.expr_.Cols());
        MatrixView<eT1,M1,N1,MatrixShell<T,meta::rows<ExprT>::value,meta::cols<ExprT>::value>> msview(ms1,
                                                        exprview.startrow_, exprview.startrow_+exprview.row_*exprview.rowstep_, 
                                                        exprview.startcol_, exprview.startcol_+exprview.col_*exprview.colstep_,
                                                        exprview.rowstep_, exprview.colstep_);
        MatrixShell<eT,M,N> ms((eT *)(block->ptr), block->size/sizeof(eT));
        eval(ms, msview, cache, -1);
        cache.StoreResult(block);
    }

    template <typename eT, typename ExprT>
    static size_t get_cache_block_size(const EntityMatrixBase<eT,ExprT> &expr)
    {
        size_t sz = sizeof(eT) * expr.Elems();
        sz = (sz + sizeof(long) - 1) / sizeof(long) * sizeof(long);
        return sz;
    }
    template <typename eT, int M, int N, typename ExprT1, typename SeT>
    static size_t get_cache_block_size(const BinaryExpr<eT,M,N,MatOpSMul,ExprT1,Scalar<SeT>> &expr)
    {
        return get_cache_block_size(expr.expr1_);
    }
    template <typename eT, int M, int N, typename ExprT1, typename SeT>
    static size_t get_cache_block_size(const BinaryExpr<eT,M,N,MatOpSMul,Scalar<SeT>,ExprT1> &expr)
    {
        return get_cache_block_size(expr.expr2_);
    }
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
    static size_t get_cache_block_size(const BinaryExpr<eT,M,N,opT,ExprT1,ExprT2> &expr)
    {
        size_t sz = std::max(sizeof(eT) * expr.Elems(), std::max(get_cache_block_size(expr.expr1_), get_cache_block_size(expr.expr2_)));
        sz = (sz + sizeof(long) - 1) / sizeof(long) * sizeof(long);
        return sz;
    }
    template <typename eT, int M, int N, typename ExprT>
    static size_t get_cache_block_size(const MatrixView<eT,M,N,ExprT> &expr)
    {
        return get_cache_block_size(expr.expr_);
    }

public:
    template <typename eT, int M, int N, typename expreT, typename ExprT>
    static const MatrixView<eT,M,N,ExprT> GenMatrixView(const MatrixBase<expreT,ExprT> &expr, int sr, int er, int rstep, int sc, int ec, int cstep)
    {
        return MatrixView<eT,M,N,ExprT>(expr.Derived(), sr, er, rstep, sc, ec, cstep);
    }
    template <typename eT, int M, int N, typename expreT, typename ExprT>
    static const MatrixView<eT,M,N,const ExprT> GenMatrixView(const EntityMatrixBase<expreT,ExprT> &expr, int sr, int er, int rstep, int sc, int ec, int cstep)
    {
        return MatrixView<eT,M,N,const ExprT>(expr.Derived(), sr, er, rstep, sc, ec, cstep);
    }
    template <typename eT, int M, int N, typename expreT, typename ExprT>
    static MatrixView<eT,M,N,ExprT> GenMatrixView(EntityMatrixBase<expreT,ExprT> &expr, int sr, int er, int rstep, int sc, int ec, int cstep)
    {
        return MatrixView<eT,M,N,ExprT>(expr.Derived(), sr, er, rstep, sc, ec, cstep);
    }
    template <typename eT, int M, int N, typename expreT, int M1, int N1, typename ExprT>
    static const MatrixView<eT,M,N,ExprT> GenMatrixView(const MatrixView<expreT,M1,N1,ExprT> &expr, int sr, int er, int rstep, int sc, int ec, int cstep)
    {
        return MatrixView<eT,0,0,ExprT>(expr.expr_, expr.startrow_ + sr * expr.rowstep_, expr.startrow_ + er * expr.rowstep_, rstep * expr.rowstep_,
                                        expr.startcol_ + sc * expr.colstep_, expr.startcol_ + ec * expr.colstep_, cstep * expr.colstep_);
    }
    template <typename eT, int M, int N, typename expreT, int M1, int N1, typename ExprT>
    static MatrixView<eT,M,N,ExprT> GenMatrixView(MatrixView<expreT,M1,N1,ExprT> &expr, int sr, int er, int rstep, int sc, int ec, int cstep)
    {
        return MatrixView<eT,0,0,ExprT>(expr.expr_, expr.startrow_ + sr * expr.rowstep_, expr.startrow_ + er * expr.rowstep_, rstep * expr.rowstep_,
                                        expr.startcol_ + sc * expr.colstep_, expr.startcol_ + ec * expr.colstep_, cstep * expr.colstep_);
    }

protected:
    template <typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    MatOp(MutableMatrixBase<eT1,ExprT1> &mat, const MatrixBase<eT2,ExprT2> &expr, int relation = -1)
        : cache(get_cache_block_size(expr.Derived()))
    {
        eval(mat.Derived(), expr.Derived(), cache, relation);
    }

private:
    Cache cache;
};

template <typename eT, int M, int N>
template <typename expreT, typename ExprT>
Matrix<eT,M,N>::Matrix(const MatrixBase<expreT,ExprT> &matexpr)
{
    MatOp op(*this, matexpr.Derived(), 0);
}
template <typename eT, int M, int N>
template <typename expreT, typename ExprT>
Matrix<eT,M,N> &Matrix<eT,M,N>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    MatOp op(*this, matexpr.Derived(), -1);
}

template <typename eT>
template <typename expreT, typename ExprT>
Matrix<eT,0,0>::Matrix(const MatrixBase<expreT,ExprT> &matexpr)
    : data_(nullptr), capacity_(0), row_(0), col_(0)
{
    MatOp op(*this, matexpr.Derived(), 0);
}

template <typename eT>
template <typename expreT, typename ExprT>
Matrix<eT,0,0> &Matrix<eT,0,0>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    MatOp op(*this, matexpr.Derived(), -1);
}

template <typename eT, int M>
template <typename expreT, typename ExprT>
Matrix<eT,M,0>::Matrix(const MatrixBase<expreT,ExprT> &matexpr)
    : data_(nullptr), capacity_(0), col_(0)
{
    MatOp op(*this, matexpr.Derived(), 0);
}

template <typename eT, int M>
template <typename expreT, typename ExprT>
Matrix<eT,M,0> &Matrix<eT,M,0>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    MatOp op(*this, matexpr.Derived(), -1);
}

template <typename eT, int N>
template <typename expreT, typename ExprT>
Matrix<eT,0,N>::Matrix(const MatrixBase<expreT,ExprT> &matexpr)
    : data_(nullptr), capacity_(0), row_(0)
{
    MatOp op(*this, matexpr.Derived(), 0);
}

template <typename eT, int N>
template <typename expreT, typename ExprT>
Matrix<eT,0,N> &Matrix<eT,0,N>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    MatOp op(*this, matexpr.Derived(), -1);
}

template <typename eT, int M, int N>
template <typename expreT, typename ExprT>
MatrixShell<eT,M,N> &MatrixShell<eT,M,N>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    MatOp op(*this, matexpr.Derived(), -1);
}

template <typename eT>
template <typename expreT, typename ExprT>
MatrixShell<eT,0,0> &MatrixShell<eT,0,0>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    MatOp op(*this, matexpr.Derived(), -1);
}

template <typename eT, int M>
template <typename expreT, typename ExprT>
MatrixShell<eT,M,0> &MatrixShell<eT,M,0>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    MatOp op(*this, matexpr.Derived(), -1);
}

template <typename eT, int N>
template <typename expreT, typename ExprT>
MatrixShell<eT,0,N> &MatrixShell<eT,0,N>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    MatOp op(*this, matexpr.Derived(), -1);
}

template <typename eT, int M, int N, typename ExprT>
MatrixView<eT,M,N,ExprT>& MatrixView<eT,M,N,ExprT>::operator=(const MatrixView &matexpr)
{
    MatOp op(*this, matexpr, -1);
}
template <typename eT, int M, int N, typename ExprT>
template <typename expreT, typename ExprT1>
MatrixView<eT,M,N,ExprT>& MatrixView<eT,M,N,ExprT>::operator=(const MatrixBase<expreT, ExprT1> &matexpr)
{
    MatOp op(*this, matexpr.Derived(), -1);
}

template <typename eT, typename ExprT>
auto SubMatrix(const MatrixBase<eT,ExprT> &expr, int sr, int er, int rstep, int sc, int ec, int cstep) -> decltype(MatOp::GenMatrixView<eT,0,0>(expr.Derived(),sr,er,rstep,sc,ec,cstep))
{
    return MatOp::GenMatrixView<eT,0,0>(expr.Derived(), sr, er, rstep, sc, ec, cstep);
}
template <typename eT, typename ExprT>
auto SubMatrix(MatrixBase<eT,ExprT> &expr, int sr, int er, int rstep, int sc, int ec, int cstep) -> decltype(MatOp::GenMatrixView<eT,0,0>(expr.Derived(),sr,er,rstep,sc,ec,cstep))
{
    return MatOp::GenMatrixView<eT,0,0>(expr.Derived(), sr, er, rstep, sc, ec, cstep);
}

} // end of namespace narutoacm



#endif

