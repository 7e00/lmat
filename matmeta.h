#ifndef NARUTOACM_MATMETA_H_
#define NARUTOACM_MATMETA_H_

namespace narutoacm
{
template <typename eT>
class Scalar;
template <typename eT, int M, int N>
class Matrix;
template <typename eT, typename opT, typename ExprT1, typename ExprT2>
class BinaryExpr;


class OpAdd;
class OpSub;
class OpEMul; // element-wise multiplication
class OpEDiv;
class OpSMul; // matrix * scalar
class OpMul; // matrix multiplication
class OpTran; // transpose


typedef double f64;
typedef float f32;
typedef long long s64;
typedef unsigned long long u64;
typedef int s32;
typedef unsigned int u32;
typedef short s16;
typedef unsigned short u16;
typedef signed char s8;
typedef unsigned char u8;


namespace meta
{
namespace helper
{
template <bool>
struct assert;
template<>
struct assert<true> {};

template <int a, int b>
struct max
{
    enum { result = (a > b ? a : b) };
};
template <int a, int b>
struct min
{
    enum { result = (a < b ? a : b) };
};

template <typename eT1,typename eT2>
struct is_same_type
{
    enum { result = false };
};
template <typename eT>
struct is_same_type<eT,eT>
{
    enum { result = true };
};

template <bool, typename eT1, typename eT2>
struct select_type;
template <typename eT1, typename eT2>
struct select_type<true,eT1,eT2>
{
    typedef eT1 result_type;
};
template <typename eT1, typename eT2>
struct select_type<false,eT1,eT2>
{
    typedef eT2 result_type;
};

} // end of namespace narutoacm::meta::helper

namespace types
{
template <typename eT1, typename eT2>
struct prio_type { typedef eT1 result_type; };
template<typename eT> struct prio_type<eT,f64> { typedef f64 result_type; };
template<> struct prio_type<s64,f32> { typedef f32 result_type; };
template<> struct prio_type<s32,f32> { typedef f32 result_type; };
template<> struct prio_type<s16,f32> { typedef f32 result_type; };
template<> struct prio_type<s8,f32> { typedef f32 result_type; };
template<> struct prio_type<u64,f32> { typedef f32 result_type; };
template<> struct prio_type<u32,f32> { typedef f32 result_type; };
template<> struct prio_type<u16,f32> { typedef f32 result_type; };
template<> struct prio_type<u8,f32> { typedef f32 result_type; };
template<> struct prio_type<u64,s64> { typedef s64 result_type; };
template<> struct prio_type<s32,s64> { typedef s64 result_type; };
template<> struct prio_type<s16,s64> { typedef s64 result_type; };
template<> struct prio_type<s8,s64> { typedef s64 result_type; };
template<> struct prio_type<u32,s64> { typedef s64 result_type; };
template<> struct prio_type<u16,s64> { typedef s64 result_type; };
template<> struct prio_type<u8,s64> { typedef s64 result_type; };
template<> struct prio_type<s32,u64> { typedef s64 result_type; };
template<> struct prio_type<s16,u64> { typedef s64 result_type; };
template<> struct prio_type<s8,u64> { typedef s64 result_type; };
template<> struct prio_type<u64,s32> { typedef s64 result_type; };
template<> struct prio_type<u64,s16> { typedef s64 result_type; };
template<> struct prio_type<u64,s8> { typedef s64 result_type; };
template<> struct prio_type<u32,u64> { typedef u64 result_type; };
template<> struct prio_type<u16,u64> { typedef u64 result_type; };
template<> struct prio_type<u8,u64> { typedef u64 result_type; };
template<> struct prio_type<u32,s32> { typedef s32 result_type; };
template<> struct prio_type<u16,s32> { typedef s32 result_type; };
template<> struct prio_type<u8,s32> { typedef s32 result_type; };
template<> struct prio_type<s16,s32> { typedef s32 result_type; };
template<> struct prio_type<s8,s32> { typedef s32 result_type; };
template<> struct prio_type<u32,s16> { typedef s32 result_type; };
template<> struct prio_type<u32,s8> { typedef s32 result_type; };
template<> struct prio_type<s16,u32> { typedef s32 result_type; };
template<> struct prio_type<s8,u32> { typedef s32 result_type; };
template<> struct prio_type<u16,u32> { typedef u32 result_type; };
template<> struct prio_type<u8,u32> { typedef u32 result_type; };
template<> struct prio_type<u16,s16> { typedef s16 result_type; };
template<> struct prio_type<u8,s16> { typedef s16 result_type; };
template<> struct prio_type<s8,s16> { typedef s16 result_type; };
template<> struct prio_type<u16,s8> { typedef s16 result_type; };
template<> struct prio_type<s8,u16> { typedef s16 result_type; };
template<> struct prio_type<u8,u16> { typedef u16 result_type; };
template<> struct prio_type<u8,s8> { typedef s8 result_type; };

using helper::is_same_type;
using helper::select_type;

template <typename eT1, typename eT2>
struct inferior_type
{
    typedef typename select_type<
        is_same_type<eT1,typename prio_type<eT1,eT2>::result_type>::result,
        eT2,eT1>::result_type result_type;
};
} // end of namespace narutoacm::meta::types

using helper::assert;
using helper::max;
using helper::min;

template <typename ExprT>
struct element_type;
template <typename eT, int M, int N>
struct element_type<Matrix<eT,M,N>>
{
    typedef eT type;
};
template <typename eT, typename opT, typename ExprT1, typename ExprT2>
struct element_type<BinaryExpr<eT,opT,ExprT1,ExprT2>>
{
    typedef eT type;
};

template <typename ExprT>
struct rows;
template <typename ExprT>
struct cols;
template <typename eT, int M, int N>
struct rows<Matrix<eT,M,N>>
{
    enum { value = M };
};
template <typename eT, int M, int N>
struct cols<Matrix<eT,M,N>>
{
    enum { value = N };
};
// mat+mat, mat-mat, mat .* mat, mat ./ mat
template <typename eT, typename opT, typename ExprT1, typename ExprT2>
struct rows<BinaryExpr<eT,opT,ExprT1,ExprT2>> : public assert<(rows<ExprT1>::value == 0 
                                                                || rows<ExprT2>::value == 0
                                                                || rows<ExprT1>::value == rows<ExprT2>::value)
                                                                && (cols<ExprT1>::value == 0
                                                                || cols<ExprT2>::value == 0
                                                                || cols<ExprT1>::value == cols<ExprT2>::value)>
{
    enum { value = max<rows<ExprT1>::value, rows<ExprT2>::value>::result };
};
template <typename eT, typename opT, typename ExprT1, typename ExprT2>
struct cols<BinaryExpr<eT,opT,ExprT1,ExprT2>> : public assert<(rows<ExprT1>::value == 0 
                                                                || rows<ExprT2>::value == 0
                                                                || rows<ExprT1>::value == rows<ExprT2>::value)
                                                                && (cols<ExprT1>::value == 0
                                                                || cols<ExprT2>::value == 0
                                                                || cols<ExprT1>::value == cols<ExprT2>::value)>
{
    enum { value = max<cols<ExprT1>::value, cols<ExprT2>::value>::result };
};
// mat * mat
template <typename eT, typename ExprT1, typename ExprT2>
struct rows<BinaryExpr<eT,OpMul,ExprT1,ExprT2>> : public assert<cols<ExprT1>::value == 0
                                                                || rows<ExprT2>::value == 0
                                                                || cols<ExprT1>::value == rows<ExprT2>::value>
{
    enum { value = rows<ExprT1>::value };
};
template <typename eT, typename ExprT1, typename ExprT2>
struct cols<BinaryExpr<eT,OpMul,ExprT1,ExprT2>> : public assert<cols<ExprT1>::value == 0
                                                                || rows<ExprT2>::value == 0
                                                                || cols<ExprT1>::value == rows<ExprT2>::value>
{
    enum { value = cols<ExprT2>::value };
};

template <typename ExprT>
struct matrix_expr_type;
template <typename eT, int M, int N>
struct matrix_expr_type<Matrix<eT,M,N>>
{
    enum { value = 1 };
};
template <typename eT, typename opT, typename ExprT1, typename ExprT2>
struct matrix_expr_type<BinaryExpr<eT,opT,ExprT1,ExprT2>>
{
    enum { value = min<matrix_expr_type<ExprT1>::value, matrix_expr_type<ExprT2>::value>::result };
};
template <typename eT, typename ExprT1, typename ExprT2>
struct matrix_expr_type<BinaryExpr<eT,OpMul,ExprT1,ExprT2>>
{
    enum { value = 0 };
};

} // end of namespace narutoacm::meta

} // end of namespace narutoacm

#endif

