#ifndef _XF_SOLVER_CHOLESKY_HPP_
#define _XF_SOLVER_CHOLESKY_HPP_

#ifndef FAST_SQRT_FIXED
#define FAST_SQRT_FIXED 0
#endif

// 添加兜底：若 tcl/Makefile 没传 -DSEL_ARCH，则默认用 ARCH0
#ifndef SEL_ARCH
#define SEL_ARCH 0
#endif

#include "ap_fixed.h"
#include "hls_x_complex.h"
#include <complex>
#include "utils/std_complex_utils.h"
#include "utils/x_matrix_utils.hpp"
#include "hls_stream.h"

namespace xf {
namespace solver {

// ===================================================================================================================
// Default traits
template <bool LowerTriangularL, int RowsColsA, typename InputType, typename OutputType>
struct choleskyTraits {
    typedef InputType PROD_T;
    typedef InputType ACCUM_T;
    typedef InputType ADD_T;
    typedef InputType DIAG_T;
    typedef InputType RECIP_DIAG_T;
    typedef InputType OFF_DIAG_T;
    typedef OutputType L_OUTPUT_T;
    static const int ARCH =SEL_ARCH;
    static const int INNER_II = 1;
    static const int UNROLL_FACTOR = 1;
    static const int UNROLL_DIM = (LowerTriangularL == true ? 1 : 2);
    static const int ARCH2_ZERO_LOOP = true;
};

template <bool LowerTriangularL, int RowsColsA, typename InputBaseType, typename OutputBaseType>
struct choleskyTraits<LowerTriangularL, RowsColsA, hls::x_complex<InputBaseType>, hls::x_complex<OutputBaseType> > {
    typedef hls::x_complex<InputBaseType> PROD_T;
    typedef hls::x_complex<InputBaseType> ACCUM_T;
    typedef hls::x_complex<InputBaseType> ADD_T;
    typedef hls::x_complex<InputBaseType> DIAG_T;
    typedef InputBaseType RECIP_DIAG_T;
    typedef hls::x_complex<InputBaseType> OFF_DIAG_T;
    typedef hls::x_complex<OutputBaseType> L_OUTPUT_T;
    static const int ARCH = SEL_ARCH;
    static const int INNER_II = 1;
    static const int UNROLL_FACTOR = 1;
    static const int UNROLL_DIM = (LowerTriangularL == true ? 1 : 2);
    static const int ARCH2_ZERO_LOOP = true;
};

template <bool LowerTriangularL, int RowsColsA, typename InputBaseType, typename OutputBaseType>
struct choleskyTraits<LowerTriangularL, RowsColsA, std::complex<InputBaseType>, std::complex<OutputBaseType> > {
    typedef std::complex<InputBaseType> PROD_T;
    typedef std::complex<InputBaseType> ACCUM_T;
    typedef std::complex<InputBaseType> ADD_T;
    typedef std::complex<InputBaseType> DIAG_T;
    typedef InputBaseType RECIP_DIAG_T;
    typedef std::complex<InputBaseType> OFF_DIAG_T;
    typedef std::complex<OutputBaseType> L_OUTPUT_T;
    static const int ARCH = SEL_ARCH;
    static const int INNER_II = 1;
    static const int UNROLL_FACTOR = 1;
    static const int UNROLL_DIM = (LowerTriangularL == true ? 1 : 2);
    static const int ARCH2_ZERO_LOOP = true;
};

template <bool LowerTriangularL,
          int RowsColsA,
          int W1, int I1, ap_q_mode Q1, ap_o_mode O1, int N1,
          int W2, int I2, ap_q_mode Q2, ap_o_mode O2, int N2>
struct choleskyTraits<LowerTriangularL, RowsColsA, ap_fixed<W1, I1, Q1, O1, N1>, ap_fixed<W2, I2, Q2, O2, N2> > {
    typedef ap_fixed<W1 + W1, I1 + I1, AP_RND_CONV, AP_SAT, 0> PROD_T;
    typedef ap_fixed<(W1 + W1) + BitWidth<RowsColsA>::Value,
                     (I1 + I1) + BitWidth<RowsColsA>::Value,
                     AP_RND_CONV, AP_SAT, 0> ACCUM_T;
    typedef ap_fixed<W1 + 1, I1 + 1, AP_RND_CONV, AP_SAT, 0> ADD_T;
    typedef ap_fixed<(W1 + 1) * 2, I1 + 1, AP_RND_CONV, AP_SAT, 0> DIAG_T;
    typedef ap_fixed<(W1 + 1) * 2, I1 + 1, AP_RND_CONV, AP_SAT, 0> OFF_DIAG_T;
    typedef ap_fixed<2 + (W2 - I2) + W2, 2 + (W2 - I2), AP_RND_CONV, AP_SAT, 0> RECIP_DIAG_T;
    typedef ap_fixed<W2, I2, AP_RND_CONV, AP_SAT, 0> L_OUTPUT_T;
    static const int ARCH = SEL_ARCH;
    static const int INNER_II = 1;
    static const int UNROLL_FACTOR = 1;
    static const int UNROLL_DIM = (LowerTriangularL == true ? 1 : 2);
    static const int ARCH2_ZERO_LOOP = true;
};

template <bool LowerTriangularL,
          int RowsColsA,
          int W1, int I1, ap_q_mode Q1, ap_o_mode O1, int N1,
          int W2, int I2, ap_q_mode Q2, ap_o_mode O2, int N2>
struct choleskyTraits<LowerTriangularL,
                      RowsColsA,
                      hls::x_complex<ap_fixed<W1, I1, Q1, O1, N1> >,
                      hls::x_complex<ap_fixed<W2, I2, Q2, O2, N2> > > {
    typedef hls::x_complex<ap_fixed<W1 + W1, I1 + I1, AP_RND_CONV, AP_SAT, 0> > PROD_T;
    typedef hls::x_complex<ap_fixed<(W1 + W1) + BitWidth<RowsColsA>::Value,
                                    (I1 + I1) + BitWidth<RowsColsA>::Value,
                                    AP_RND_CONV, AP_SAT, 0> > ACCUM_T;
    typedef hls::x_complex<ap_fixed<W1 + 1, I1 + 1, AP_RND_CONV, AP_SAT, 0> > ADD_T;
    typedef hls::x_complex<ap_fixed<(W1 + 1) * 2, I1 + 1, AP_RND_CONV, AP_SAT, 0> > DIAG_T;
    typedef hls::x_complex<ap_fixed<(W1 + 1) * 2, I1 + 1, AP_RND_CONV, AP_SAT, 0> > OFF_DIAG_T;
    typedef ap_fixed<2 + (W2 - I2) + W2, 2 + (W2 - I2), AP_RND_CONV, AP_SAT, 0> RECIP_DIAG_T;
    typedef hls::x_complex<ap_fixed<W2, I2, AP_RND_CONV, AP_SAT, 0> > L_OUTPUT_T;
    static const int ARCH = SEL_ARCH;
    static const int INNER_II = 1;
    static const int UNROLL_FACTOR = 1;
    static const int UNROLL_DIM = (LowerTriangularL == true ? 1 : 2);
    static const int ARCH2_ZERO_LOOP = true;
};

template <bool LowerTriangularL,
          int RowsColsA,
          int W1, int I1, ap_q_mode Q1, ap_o_mode O1, int N1,
          int W2, int I2, ap_q_mode Q2, ap_o_mode O2, int N2>
struct choleskyTraits<LowerTriangularL,
                      RowsColsA,
                      std::complex<ap_fixed<W1, I1, Q1, O1, N1> >,
                      std::complex<ap_fixed<W2, I2, Q2, O2, N2> > > {
    typedef std::complex<ap_fixed<W1 + W1, I1 + I1, AP_RND_CONV, AP_SAT, 0> > PROD_T;
    typedef std::complex<ap_fixed<(W1 + W1) + BitWidth<RowsColsA>::Value,
                                  (I1 + I1) + BitWidth<RowsColsA>::Value,
                                  AP_RND_CONV, AP_SAT, 0> > ACCUM_T;
    typedef std::complex<ap_fixed<W1 + 1, I1 + 1, AP_RND_CONV, AP_SAT, 0> > ADD_T;
    typedef std::complex<ap_fixed<(W1 + 1) * 2, I1 + 1, AP_RND_CONV, AP_SAT, 0> > DIAG_T;
    typedef std::complex<ap_fixed<(W1 + 1) * 2, I1 + 1, AP_RND_CONV, AP_SAT, 0> > OFF_DIAG_T;
    typedef ap_fixed<2 + (W2 - I2) + W2, 2 + (W2 - I2), AP_RND_CONV, AP_SAT, 0> RECIP_DIAG_T;
    typedef std::complex<ap_fixed<W2, I2, AP_RND_CONV, AP_SAT, 0> > L_OUTPUT_T;
    static const int ARCH = SEL_ARCH;
    static const int INNER_II = 1;
    static const int UNROLL_FACTOR = 1;
    static const int UNROLL_DIM = (LowerTriangularL == true ? 1 : 2);
    static const int ARCH2_ZERO_LOOP = true;
};

// ===================================================================================================================
// Helper functions
// fixed-point reciprocal NR (constant-time, narrow intermediates)
template <int W, int I, ap_q_mode Q, ap_o_mode O, int N>
inline void cholesky_reciprocal_nr(const ap_fixed<W,I,Q,O,N>& a,
                                   ap_fixed<W,I,Q,O,N>& y) {
#pragma HLS INLINE
    typedef ap_fixed<W,I,Q,O,N> FT;     // 与输出一致的位宽类型
    FT x = (FT)a;

    // 顺序归一化到 [0.5,1)
    int k = 0;
    if (x < (FT)0.0625) { x = (FT)(x * (FT)16.0); k += 4; }
    if (x < (FT)0.25)   { x = (FT)(x * (FT)4.0);  k += 2; }
    if (x < (FT)0.5)    { x = (FT)(x * (FT)2.0);  k += 1; }
    if (x > (FT)2.0)    { x = (FT)(x * (FT)0.25); k -= 2; }
    if (x > (FT)1.0)    { x = (FT)(x * (FT)0.5);  k -= 1; }

    // y0 = 2 - x，2 次 NR；每步乘法后显式截断到 FT
    const FT two = (FT)2.0;
    FT yn = (FT)(two - x);
    FT t  = (FT)(x * yn);
    yn    = (FT)(yn * (FT)(two - t));
    t     = (FT)(x * yn);
    yn    = (FT)(yn * (FT)(two - t));

    // 2^k 缩放用移位，避免乘法
    if (k > 0)      yn = (FT)(yn << k);
    else if (k < 0) yn = (FT)(yn >> (-k));

    y = yn;
}
// type-safe zero helpers to avoid int -> complex(ap_fixed) implicit conversion issues
template <typename T>
inline void set_zero(T& v) {
    v = 0;
}
template <typename B>
inline void set_zero(hls::x_complex<B>& v) {
    v.real(0);
    v.imag(0);
}
template <typename B>
inline void set_zero(std::complex<B>& v) {
    v.real(0);
    v.imag(0);
}
#if FAST_SQRT_FIXED
// 高精度定点 sqrt：把 a 归一化到 [0.5,1)，对 1/sqrt(x) 做 4 次 NR（更宽中间位宽），再回缩。
template<int WI, int II, ap_q_mode QI, ap_o_mode OI, int NI,
         int WO, int IO, ap_q_mode QO, ap_o_mode OO, int NO>
inline int cholesky_sqrt_op(ap_fixed<WI,II,QI,OI,NI> a,
                            ap_fixed<WO,IO,QO,OO,NO>& b) {
#pragma HLS INLINE
    // 负数/零保护
    const ap_fixed<WI,II,QI,OI,NI> ZERO = 0;
    if (a <= ZERO) { b = 0; return 0; }

    // 宽工作位宽（比输出宽一些，保证收敛精度）
    typedef ap_fixed<WO+8, IO+4, AP_RND_CONV, AP_SAT> WT;
    WT x = (WT)a;

    // 归一化到 [0.5, 1)
    int k = 0;
    if (x < (WT)0.0625) { x *= (WT)16.0; k += 4; }
    if (x < (WT)0.25  ) { x *= (WT)4.0;  k += 2; }
    if (x < (WT)0.5   ) { x *= (WT)2.0;  k += 1; }
    if (x > (WT)2.0   ) { x *= (WT)0.25; k -= 2; }
    if (x > (WT)1.0   ) { x *= (WT)0.5;  k -= 1; }

    // 初值（在 x≈1 处的一阶展开，硬件便宜，4 次 NR 足够 16bit 精度）
    WT y = (WT)1.0 - (WT)0.5 * (x - (WT)1.0);

    const WT half  = (WT)0.5;
    const WT onep5 = (WT)1.5;

    // 4 次 NR： y = y*(1.5 - 0.5*x*y*y)
    for (int it = 0; it < 4; ++it) {
#pragma HLS UNROLL
        WT t = x * y;
        y = y * (onep5 - half * t * y);
    }

    // sqrt(x) ≈ x * y
    WT s = x * y;

    // 回缩：sqrt(a) = s * 2^(k/2)
    static const WT RT2     = (WT)1.4142135623730950488; // sqrt(2)
    static const WT INV_RT2 = (WT)0.7071067811865475244; // 1/sqrt(2)
    if (k > 0) {
        int p = k >> 1;              // floor(k/2)
        if (p)   s = (WT)(s >> p);   // /2^p
        if (k&1) s = s * INV_RT2;    // /sqrt(2)
    } else if (k < 0) {
        int p = (-k) >> 1;
        if (p)     s = (WT)(s << p); // *2^p
        if ((-k)&1) s = s * RT2;     // *sqrt(2)
    }

    b = (ap_fixed<WO,IO,QO,OO,NO>)s;
    return 0;
}

// 复数：取实部开方，虚部清零（和库行为一致）
template<int WI, int II, ap_q_mode QI, ap_o_mode OI, int NI,
         int WO, int IO, ap_q_mode QO, ap_o_mode OO, int NO>
inline int cholesky_sqrt_op(hls::x_complex<ap_fixed<WI,II,QI,OI,NI> > din,
                            hls::x_complex<ap_fixed<WO,IO,QO,OO,NO> >& dout) {
#pragma HLS INLINE
    typedef ap_fixed<WI,II,QI,OI,NI> TI;
    typedef ap_fixed<WO,IO,QO,OO,NO> TO;
    const TI ZERO = 0;
    TI a = din.real();
    dout.imag((TO)0);
    if (a <= ZERO) { dout.real((TO)0); return 0; }
    TO s;
    cholesky_sqrt_op(a, s); // 调用上面的实数版本
    dout.real(s);
    return 0;
}
#endif // FAST_SQRT_FIXED
// Square root
template <typename T_IN, typename T_OUT>
int cholesky_sqrt_op(T_IN a, T_OUT& b) {
Function_cholesky_sqrt_op_real:;
    const T_IN ZERO = 0;
    if (a < ZERO) {
        b = ZERO;
        return (1);
    }
    b = x_sqrt(a);
    return (0);
}
template <typename T_IN, typename T_OUT>
int cholesky_sqrt_op(hls::x_complex<T_IN> din, hls::x_complex<T_OUT>& dout) {
Function_cholesky_sqrt_op_complex:;
    const T_IN ZERO = 0;
    T_IN a = din.real();
    dout.imag(ZERO);
    if (a < ZERO) {
        dout.real(ZERO);
        return (1);
    }
    dout.real(x_sqrt(a));
    return (0);
}
template <typename T_IN, typename T_OUT>
int cholesky_sqrt_op(std::complex<T_IN> din, std::complex<T_OUT>& dout) {
Function_cholesky_sqrt_op_complex:;
    const T_IN ZERO = 0;
    T_IN a = din.real();
    dout.imag(ZERO);
    if (a < ZERO) {
        dout.real(ZERO);
        return (1);
    }
    dout.real(x_sqrt(a));
    return (0);
}

// Reciprocal square root.
template <typename InputType, typename OutputType>
void cholesky_rsqrt(InputType x, OutputType& res) {
Function_cholesky_rsqrt_default:;
    res = x_rsqrt(x);
}
template <int W1, int I1, ap_q_mode Q1, ap_o_mode O1, int N1, int W2, int I2, ap_q_mode Q2, ap_o_mode O2, int N2>
void cholesky_rsqrt(ap_fixed<W1, I1, Q1, O1, N1> x, ap_fixed<W2, I2, Q2, O2, N2>& res) {
Function_cholesky_rsqrt_fixed:;
    ap_fixed<W2, I2, Q2, O2, N2> one = 1;
    ap_fixed<W1, I1, Q1, O1, N1> sqrt_res;
    ap_fixed<W2, I2, Q2, O2, N2> sqrt_res_cast;
    sqrt_res = x_sqrt(x);
    sqrt_res_cast = sqrt_res;
    res = one / sqrt_res_cast;
}

// Local multiplier
template <typename AType, typename BType, typename CType>
void cholesky_prod_sum_mult(AType A, BType B, CType& C) {
Function_cholesky_prod_sum_mult_real:;
    C = A * B;
}
template <typename AType, typename BType, typename CType>
void cholesky_prod_sum_mult(hls::x_complex<AType> A, BType B, hls::x_complex<CType>& C) {
Function_cholesky_prod_sum_mult_complex:;
    C.real(A.real() * B);
    C.imag(A.imag() * B);
}
template <typename AType, typename BType, typename CType>
void cholesky_prod_sum_mult(std::complex<AType> A, BType B, std::complex<CType>& C) {
Function_cholesky_prod_sum_mult_complex:;
    C.real(A.real() * B);
    C.imag(A.imag() * B);
}

// ===================================================================================================================
// ARCH0: choleskyBasic (安全版优化 + 类型安全)
template <bool LowerTriangularL, int RowsColsA, typename CholeskyTraits, class InputType, class OutputType>
int choleskyBasic(const InputType A[RowsColsA][RowsColsA], OutputType L[RowsColsA][RowsColsA]) {
#pragma HLS INLINE off
    int return_code = 0;

    typename CholeskyTraits::PROD_T prod;
    typename CholeskyTraits::ACCUM_T A_cast_to_sum;
    typename CholeskyTraits::ACCUM_T prod_cast_to_sum;

    typename CholeskyTraits::ADD_T A_minus_sum;
    typename CholeskyTraits::DIAG_T new_L_diag;
    typename CholeskyTraits::OFF_DIAG_T new_L_off_diag;

    typename CholeskyTraits::L_OUTPUT_T new_L;
    OutputType retrieved_L;

    OutputType L_internal[RowsColsA][RowsColsA];
#pragma HLS ARRAY_PARTITION variable=L_internal complete dim=0

col_loop:
    for (int j = 0; j < RowsColsA; j++) {
        typename CholeskyTraits::ACCUM_T sum_val; set_zero(sum_val);

    diag_loop:
        for (int k = 0; k < RowsColsA; k++) {
            if (k <= (j - 1)) {
                if (LowerTriangularL == true) {
                    retrieved_L = L_internal[j][k];
                } else {
                    retrieved_L = L_internal[k][j];
                }
                // 修正：
                sum_val += hls::x_conj(retrieved_L) * retrieved_L;
            }
        }

        A_cast_to_sum = A[j][j];
        A_minus_sum = A_cast_to_sum - sum_val;

        if (cholesky_sqrt_op(A_minus_sum, new_L_diag)) {
#ifndef __SYNTHESIS__
            printf("ERROR: Trying to find the square root of a negative number\n");
#endif
            return_code = 1;
        }

        new_L = new_L_diag;

        if (LowerTriangularL == true) {
            L_internal[j][j] = new_L;
            L[j][j] = new_L;
        } else {
            L_internal[j][j] = hls::x_conj(new_L);
            L[j][j] = hls::x_conj(new_L);
        }

        // 每列一次精确倒数
        typename CholeskyTraits::RECIP_DIAG_T inv_diag;
        {
            typename CholeskyTraits::RECIP_DIAG_T denom =
    (typename CholeskyTraits::RECIP_DIAG_T)hls::x_real(L_internal[j][j]);
cholesky_reciprocal_nr(denom, inv_diag);  // NR 版本的倒数
        }

    off_diag_loop:
        for (int i = j + 1; i < RowsColsA; i++) { // 仅计算需要的 i
            typename CholeskyTraits::ACCUM_T sum_prod;
            if (LowerTriangularL == true) {
                sum_prod = A[i][j];
            } else {
                sum_prod = hls::x_conj(A[j][i]);
            }

        sum_loop:
            for (int k = 0; k < RowsColsA; k++) {
#pragma HLS PIPELINE II = CholeskyTraits::INNER_II
                if (k <= (j - 1)) {
                    OutputType Ljk = (LowerTriangularL ? L_internal[j][k] : L_internal[k][j]);
                    auto Ljk_conj = hls::x_conj(Ljk);

                    if (LowerTriangularL == true) {
                        prod = -L_internal[i][k] * Ljk_conj;
                    } else {
                        prod = -hls::x_conj(L_internal[k][i]) * (Ljk);
                    }
                    prod_cast_to_sum = prod;
                    sum_prod += prod_cast_to_sum;
                }
            }

            cholesky_prod_sum_mult(sum_prod, inv_diag, new_L_off_diag);
            new_L = new_L_off_diag;

            if (LowerTriangularL == true) {
                L[i][j] = new_L;
                L_internal[i][j] = new_L;
            } else {
                L[j][i] = hls::x_conj(new_L);
                L_internal[j][i] = hls::x_conj(new_L);
            }
        }
    }

    // 统一零化另一半三角
    OutputType zero_c; set_zero(zero_c);
zero_rows_loop:
    for (int r = 0; r < RowsColsA - 1; r++) {
#pragma HLS PIPELINE
    zero_cols_loop:
        for (int c = r + 1; c < RowsColsA; c++) {
            if (LowerTriangularL == true) {
                L[r][c] = zero_c;
            } else {
                L[c][r] = zero_c;
            }
        }
    }

    return (return_code);
}

// ===================================================================================================================
// ARCH1 (原样)
template <bool LowerTriangularL, int RowsColsA, typename CholeskyTraits, class InputType, class OutputType>
int choleskyAlt(const InputType A[RowsColsA][RowsColsA], OutputType L[RowsColsA][RowsColsA]) {
    int return_code = 0;

    OutputType L_internal[(RowsColsA * RowsColsA - RowsColsA) / 2];
    typename CholeskyTraits::RECIP_DIAG_T diag_internal[RowsColsA];

    typename CholeskyTraits::ACCUM_T square_sum;
    typename CholeskyTraits::ACCUM_T A_cast_to_sum;
    typename CholeskyTraits::ADD_T A_minus_sum;
    typename CholeskyTraits::DIAG_T A_minus_sum_cast_diag;
    typename CholeskyTraits::DIAG_T new_L_diag;
    typename CholeskyTraits::RECIP_DIAG_T new_L_diag_recip;
    typename CholeskyTraits::PROD_T prod;
    typename CholeskyTraits::ACCUM_T prod_cast_to_sum;
    typename CholeskyTraits::ACCUM_T product_sum;
    typename CholeskyTraits::OFF_DIAG_T prod_cast_to_off_diag;
    typename CholeskyTraits::RECIP_DIAG_T L_diag_recip;
    typename CholeskyTraits::OFF_DIAG_T new_L_off_diag;
    typename CholeskyTraits::L_OUTPUT_T new_L;
    typename CholeskyTraits::L_OUTPUT_T new_L_recip;

row_loop:
    for (int i = 0; i < RowsColsA; i++) {
        int i_sub1 = i - 1;
        int i_off = ((i_sub1 * i_sub1 - i_sub1) / 2) + i_sub1;

        square_sum = 0;
    col_loop:
        for (int j = 0; j < i; j++) {
#pragma HLS loop_tripcount max = 1 + RowsColsA / 2
            int j_sub1 = j - 1;
            int j_off = ((j_sub1 * j_sub1 - j_sub1) / 2) + j_sub1;
            if (LowerTriangularL == true) {
                product_sum = A[i][j];
            } else {
                product_sum = hls::x_conj(A[j][i]);
            }
        sum_loop:
            for (int k = 0; k < j; k++) {
#pragma HLS loop_tripcount max = 1 + RowsColsA / 2
#pragma HLS PIPELINE II = CholeskyTraits::INNER_II
                prod = -L_internal[i_off + k] * hls::x_conj(L_internal[j_off + k]);
                prod_cast_to_sum = prod;
                product_sum += prod_cast_to_sum;
            }
            prod_cast_to_off_diag = product_sum;
            L_diag_recip = diag_internal[j];
            cholesky_prod_sum_mult(prod_cast_to_off_diag, L_diag_recip, new_L_off_diag);
            new_L = new_L_off_diag;
            square_sum += hls::x_conj(new_L) * new_L;
            L_internal[i_off + j] = new_L;
            if (LowerTriangularL == true) {
                L[i][j] = new_L;
                L[j][i] = 0;
            } else {
                L[j][i] = hls::x_conj(new_L);
                L[i][j] = 0;
            }
        }

        A_cast_to_sum = A[i][i];
        A_minus_sum = A_cast_to_sum - square_sum;
        if (cholesky_sqrt_op(A_minus_sum, new_L_diag)) {
#ifndef __SYNTHESIS__
            printf("ERROR: Trying to find the square root of a negative number\n");
#endif
            return_code = 1;
        }
        new_L = new_L_diag;
        A_minus_sum_cast_diag = A_minus_sum;
        cholesky_rsqrt(hls::x_real(A_minus_sum_cast_diag), new_L_diag_recip);
        diag_internal[i] = new_L_diag_recip;
        if (LowerTriangularL == true) {
            L[i][i] = new_L;
        } else {
            L[i][i] = hls::x_conj(new_L);
        }
    }
    return (return_code);
}

// ARCH2 (原样)
template <bool LowerTriangularL, int RowsColsA, typename CholeskyTraits, class InputType, class OutputType>
int choleskyAlt2(const InputType A[RowsColsA][RowsColsA], OutputType L[RowsColsA][RowsColsA]) {
    int return_code = 0;

    OutputType L_internal[RowsColsA][RowsColsA];
    OutputType prod_column_top;
    typename CholeskyTraits::ACCUM_T square_sum_array[RowsColsA];
    typename CholeskyTraits::ACCUM_T A_cast_to_sum;
    typename CholeskyTraits::ADD_T A_minus_sum;
    typename CholeskyTraits::DIAG_T A_minus_sum_cast_diag;
    typename CholeskyTraits::DIAG_T new_L_diag;
    typename CholeskyTraits::RECIP_DIAG_T new_L_diag_recip;
    typename CholeskyTraits::PROD_T prod;
    typename CholeskyTraits::ACCUM_T prod_cast_to_sum;
    typename CholeskyTraits::ACCUM_T product_sum;
    typename CholeskyTraits::ACCUM_T product_sum_array[RowsColsA];
    typename CholeskyTraits::OFF_DIAG_T prod_cast_to_off_diag;
    typename CholeskyTraits::OFF_DIAG_T new_L_off_diag;
    typename CholeskyTraits::L_OUTPUT_T new_L;

#pragma HLS ARRAY_PARTITION variable = A cyclic dim = CholeskyTraits::UNROLL_DIM factor = CholeskyTraits::UNROLL_FACTOR
#pragma HLS ARRAY_PARTITION variable = L cyclic dim = CholeskyTraits::UNROLL_DIM factor = CholeskyTraits::UNROLL_FACTOR
#pragma HLS ARRAY_PARTITION variable = L_internal cyclic dim = CholeskyTraits::UNROLL_DIM factor = \
    CholeskyTraits::UNROLL_FACTOR
#pragma HLS ARRAY_PARTITION variable = square_sum_array cyclic dim = 1 factor = CholeskyTraits::UNROLL_FACTOR
#pragma HLS ARRAY_PARTITION variable = product_sum_array cyclic dim = 1 factor = CholeskyTraits::UNROLL_FACTOR

col_loop:
    for (int j = 0; j < RowsColsA; j++) {
        A_cast_to_sum = A[j][j];
        if (j == 0) {
            A_minus_sum = A_cast_to_sum;
        } else {
            A_minus_sum = A_cast_to_sum - square_sum_array[j];
        }
        if (cholesky_sqrt_op(A_minus_sum, new_L_diag)) {
#ifndef __SYNTHESIS__
            printf("ERROR: Trying to find the square root of a negative number\n");
#endif
            return_code = 1;
        }
        new_L = new_L_diag;
        A_minus_sum_cast_diag = A_minus_sum;
        cholesky_rsqrt(hls::x_real(A_minus_sum_cast_diag), new_L_diag_recip);
        if (LowerTriangularL == true) {
            L[j][j] = new_L;
        } else {
            L[j][j] = hls::x_conj(new_L);
        }

    sum_loop:
        for (int k = 0; k <= j; k++) {
#pragma HLS loop_tripcount max = 1 + RowsColsA / 2
            prod_column_top = -hls::x_conj(L_internal[j][k]);

        row_loop:
            for (int i = 0; i < RowsColsA; i++) {
#pragma HLS LOOP_FLATTEN off
#pragma HLS PIPELINE II = CholeskyTraits::INNER_II
#pragma HLS UNROLL FACTOR = CholeskyTraits::UNROLL_FACTOR
                if (i > j) {
                    prod = L_internal[i][k] * prod_column_top;
                    prod_cast_to_sum = prod;

                    if (k == 0) {
                        if (LowerTriangularL == true) {
                            A_cast_to_sum = A[i][j];
                        } else {
                            A_cast_to_sum = hls::x_conj(A[j][i]);
                        }
                        product_sum = A_cast_to_sum;
                    } else {
                        product_sum = product_sum_array[i];
                    }

                    if (k < j) {
                        product_sum_array[i] = product_sum + prod_cast_to_sum;
                    } else {
                        prod_cast_to_off_diag = product_sum;
                        cholesky_prod_sum_mult(prod_cast_to_off_diag, new_L_diag_recip, new_L_off_diag);
                        new_L = new_L_off_diag;
                        square_sum_array[j] = hls::x_conj(new_L) * new_L;
                        L_internal[i][j] = new_L;
                        if (LowerTriangularL == true) {
                            L[i][j] = new_L;
                            if (!CholeskyTraits::ARCH2_ZERO_LOOP) L[j][i] = 0;
                        } else {
                            L[j][i] = hls::x_conj(new_L);
                            if (!CholeskyTraits::ARCH2_ZERO_LOOP) L[i][j] = 0;
                        }
                    }
                }
            }
        }
    }

    if (CholeskyTraits::ARCH2_ZERO_LOOP) {
    zero_rows_loop:
        for (int i = 0; i < RowsColsA - 1; i++) {
        zero_cols_loop:
            for (int j = i + 1; j < RowsColsA; j++) {
#pragma HLS loop_tripcount max = 1 + RowsColsA / 2
#pragma HLS PIPELINE
                if (LowerTriangularL == true) {
                    L[i][j] = 0;
                } else {
                    L[j][i] = 0;
                }
            }
        }
    }
    return (return_code);
}

// ===================================================================================================================
template <bool LowerTriangularL, int RowsColsA, typename CholeskyTraits, class InputType, class OutputType>
int choleskyTop(const InputType A[RowsColsA][RowsColsA], OutputType L[RowsColsA][RowsColsA]) {
    switch (CholeskyTraits::ARCH) {
        case 0:
            return choleskyBasic<LowerTriangularL, RowsColsA, CholeskyTraits, InputType, OutputType>(A, L);
        case 1:
            return choleskyAlt<LowerTriangularL, RowsColsA, CholeskyTraits, InputType, OutputType>(A, L);
        case 2:
            return choleskyAlt2<LowerTriangularL, RowsColsA, CholeskyTraits, InputType, OutputType>(A, L);
        default:
            return choleskyBasic<LowerTriangularL, RowsColsA, CholeskyTraits, InputType, OutputType>(A, L);
    }
}

/**
* @brief cholesky
*/
template <bool LowerTriangularL,
          int RowsColsA,
          class InputType,
          class OutputType,
          typename TRAITS = choleskyTraits<LowerTriangularL, RowsColsA, InputType, OutputType> >
int cholesky(hls::stream<InputType>& matrixAStrm, hls::stream<OutputType>& matrixLStrm) {
    InputType A[RowsColsA][RowsColsA];
    OutputType L[RowsColsA][RowsColsA];

#pragma HLS ARRAY_PARTITION variable=A complete dim=0
#pragma HLS ARRAY_PARTITION variable=L complete dim=0

    for (int r = 0; r < RowsColsA; r++) {
#pragma HLS PIPELINE
        for (int c = 0; c < RowsColsA; c++) {
            matrixAStrm.read(A[r][c]);
        }
    }

    int ret = 0;
    ret = choleskyTop<LowerTriangularL, RowsColsA, TRAITS, InputType, OutputType>(A, L);

    for (int r = 0; r < RowsColsA; r++) {
#pragma HLS PIPELINE
        for (int c = 0; c < RowsColsA; c++) {
            matrixLStrm.write(L[r][c]);
        }
    }
    return ret;
}

} // end namespace solver
} // end namespace xf
#endif