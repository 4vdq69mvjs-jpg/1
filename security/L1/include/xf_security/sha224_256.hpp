/* security/L1/include/xf_security/sha224_256.hpp */
#ifndef _XF_SECURITY_SHA224_256_HPP_
#define _XF_SECURITY_SHA224_256_HPP_

#include <ap_int.h>
#include <hls_stream.h>

#include "xf_security/types.hpp"
#include "xf_security/utils.hpp"

#ifndef __SYNTHESIS__
#include <cstdio>
#endif
#ifndef _DEBUG
#define _DEBUG (0)
#endif
#define _XF_SECURITY_VOID_CAST static_cast<void>
#define _XF_SECURITY_PRINT(msg...) \
    do {                           \
        if (_DEBUG) printf(msg);   \
    } while (0)

#define ROTR(n, x) ((x >> n) | (x << (32 - n)))
#define ROTL(n, x) ((x << n) | (x >> (32 - n)))
#define SHR(n, x) (x >> n)
#define CH(x, y, z) ((x & y) ^ ((~x) & z))
#define MAJ(x, y, z) ((x & y) ^ (x & z) ^ (y & z))
#define BSIG0(x) (ROTR(2, x) ^ ROTR(13, x) ^ ROTR(22, x))
#define BSIG1(x) (ROTR(6, x) ^ ROTR(11, x) ^ ROTR(25, x))
#define SSIG0(x) (ROTR(7, x) ^ ROTR(18, x) ^ SHR(3, x))
#define SSIG1(x) (ROTR(17, x) ^ ROTR(19, x) ^ SHR(10, x))

namespace xf {
namespace security {
namespace internal {

struct SHA256Block {
    uint32_t M[16];
};

template <bool do_sha224>
struct sha256_digest_config;
template <> struct sha256_digest_config<true> { static const short numH = 7; };
template <> struct sha256_digest_config<false> { static const short numH = 8; };

/********************
 * preProcessing(32b)
 ********************/
inline void preProcessing(hls::stream<ap_uint<32> >& msg_strm,
                          hls::stream<ap_uint<64> >& len_strm,
                          hls::stream<bool>& end_len_strm,
                          hls::stream<SHA256Block>& blk_strm,
                          hls::stream<uint64_t>& nblk_strm,
                          hls::stream<bool>& end_nblk_strm) {
LOOP_SHA256_GENENERATE_MAIN_32:
    for (bool end_flag = end_len_strm.read(); !end_flag; end_flag = end_len_strm.read()) {
        uint64_t len = len_strm.read();
        uint64_t L = 8 * len;
        uint64_t blk_num = (len >> 6) + 1 + ((len & 0x3f) > 55);
        nblk_strm.write(blk_num);
        end_nblk_strm.write(false);

    LOOP_SHA256_GEN_FULL_BLKS_32:
        for (uint64_t j = 0; j < uint64_t(len >> 6); ++j) {
#pragma HLS pipeline II = 16
#pragma HLS loop_tripcount min = 0 max = 1
            SHA256Block b0;
#pragma HLS array_partition variable = b0.M complete
        LOOP_SHA256_GEN_ONE_FULL_BLK_32:
            for (int i = 0; i < 16; ++i) {
#pragma HLS unroll
                uint32_t l = msg_strm.read();
                l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8) |
                    ((0xff000000UL & l) >> 24);
                b0.M[i] = l;
            }
            blk_strm.write(b0);
        }

        char left = (char)(len & 0x3fULL);
        if (left == 0) {
            SHA256Block b;
#pragma HLS array_partition variable = b.M complete
            b.M[0] = 0x80000000UL;
            for (int i = 1; i < 14; ++i) {
#pragma HLS unroll
                b.M[i] = 0;
            }
            b.M[14] = (uint32_t)(0xffffffffUL & (L >> 32));
            b.M[15] = (uint32_t)(0xffffffffUL & (L));
            blk_strm.write(b);
        } else if (left < 56) {
            SHA256Block b;
#pragma HLS array_partition variable = b.M complete
        LOOP_SHA256_GEN_COPY_TAIL_AND_ONE_32:
            for (int i = 0; i < 14; ++i) {
#pragma HLS pipeline
                if (i < (left >> 2)) {
                    uint32_t l = msg_strm.read();
                    l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8) |
                        ((0xff000000UL & l) >> 24);
                    b.M[i] = l;
                } else if (i > (left >> 2)) {
                    b.M[i] = 0UL;
                } else {
                    uint32_t e = left & 3L;
                    if (e == 0) {
                        b.M[i] = 0x80000000UL;
                    } else if (e == 1) {
                        uint32_t l = msg_strm.read();
                        l = ((0x000000ffUL & l) << 24);
                        b.M[i] = l | 0x00800000UL;
                    } else if (e == 2) {
                        uint32_t l = msg_strm.read();
                        l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8);
                        b.M[i] = l | 0x00008000UL;
                    } else {
                        uint32_t l = msg_strm.read();
                        l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8);
                        b.M[i] = l | 0x00000080UL;
                    }
                }
            }
            b.M[14] = (uint32_t)(0xffffffffUL & (L >> 32));
            b.M[15] = (uint32_t)(0xffffffffUL & (L));
            blk_strm.write(b);
        } else {
            SHA256Block b;
#pragma HLS array_partition variable = b.M complete
        LOOP_SHA256_GEN_COPY_TAIL_ONLY_32:
            for (int i = 0; i < 16; ++i) {
#pragma HLS unroll
                if (i < (left >> 2)) {
                    uint32_t l = msg_strm.read();
                    l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8) |
                        ((0xff000000UL & l) >> 24);
                    b.M[i] = l;
                } else if (i > (left >> 2)) {
                    b.M[i] = 0UL;
                } else {
                    uint32_t e = left & 3L;
                    if (e == 0) {
                        b.M[i] = 0x80000000UL;
                    } else if (e == 1) {
                        uint32_t l = msg_strm.read();
                        l = ((0x000000ffUL & l) << 24);
                        b.M[i] = l | 0x00800000UL;
                    } else if (e == 2) {
                        uint32_t l = msg_strm.read();
                        l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8);
                        b.M[i] = l | 0x00008000UL;
                    } else {
                        uint32_t l = msg_strm.read();
                        l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8);
                        b.M[i] = l | 0x00000080UL;
                    }
                }
            }
            blk_strm.write(b);

            SHA256Block b1;
#pragma HLS array_partition variable = b1.M complete
            for (int i = 0; i < 14; ++i) {
#pragma HLS unroll
                b1.M[i] = 0;
            }
            b1.M[14] = (uint32_t)(0xffffffffUL & (L >> 32));
            b1.M[15] = (uint32_t)(0xffffffffUL & (L));
            blk_strm.write(b1);
        }
    }
    end_nblk_strm.write(true);
}

/********************
 * preProcessing(64b)
 ********************/
inline void preProcessing(hls::stream<ap_uint<64> >& msg_strm,
                          hls::stream<ap_uint<64> >& len_strm,
                          hls::stream<bool>& end_len_strm,
                          hls::stream<SHA256Block>& blk_strm,
                          hls::stream<uint64_t>& nblk_strm,
                          hls::stream<bool>& end_nblk_strm) {
LOOP_SHA256_GENENERATE_MAIN_64:
    for (bool end_flag = end_len_strm.read(); !end_flag; end_flag = end_len_strm.read()) {
        uint64_t len = len_strm.read();
        uint64_t L = 8 * len;
        uint64_t blk_num = (len >> 6) + 1 + ((len & 0x3f) > 55);
        nblk_strm.write(blk_num);
        end_nblk_strm.write(false);

    LOOP_SHA256_GEN_FULL_BLKS_64:
        for (uint64_t j = 0; j < uint64_t(len >> 6); ++j) {
#pragma HLS pipeline II = 16
#pragma HLS loop_tripcount min = 0 max = 1
            SHA256Block b0;
#pragma HLS array_partition variable = b0.M complete
        LOOP_SHA256_GEN_ONE_FULL_BLK_64:
            for (int i = 0; i < 16; i += 2) {
#pragma HLS unroll
                uint64_t ll = msg_strm.read().to_uint64();
                uint32_t l = ll & 0xffffffffUL;
                l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8) |
                    ((0xff000000UL & l) >> 24);
                b0.M[i] = l;
                l = (ll >> 32) & 0xffffffffUL;
                l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8) |
                    ((0xff000000UL & l) >> 24);
                b0.M[i + 1] = l;
            }
            blk_strm.write(b0);
        }

        char left = (char)(len & 0x3fULL);
        if (left == 0) {
            SHA256Block b;
#pragma HLS array_partition variable = b.M complete
            b.M[0] = 0x80000000UL;
            for (int i = 1; i < 14; ++i) {
#pragma HLS unroll
                b.M[i] = 0;
            }
            b.M[14] = (uint32_t)(0xffffffffUL & (L >> 32));
            b.M[15] = (uint32_t)(0xffffffffUL & (L));
            blk_strm.write(b);
        } else {
            SHA256Block b;
#pragma HLS array_partition variable = b.M complete
        LOOP_SHA256_GEN_COPY_TAIL_PAD_ONE_64:
            for (int i = 0; i < ((left < 56) ? 7 : 8); ++i) {
#pragma HLS pipeline
                if (i < (left >> 3)) {
                    uint64_t ll = msg_strm.read().to_uint64();
                    uint32_t l = ll & 0xffffffffUL;
                    l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8) |
                        ((0xff000000UL & l) >> 24);
                    b.M[i * 2] = l;
                    l = (ll >> 32) & 0xffffffffUL;
                    l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8) |
                        ((0xff000000UL & l) >> 24);
                    b.M[i * 2 + 1] = l;
                } else if (i > (left >> 3)) {
                    b.M[i * 2] = 0UL;
                    b.M[i * 2 + 1] = 0UL;
                } else {
                    if ((left & 4) == 0) {
                        uint32_t e = left & 3L;
                        if (e == 0) {
                            b.M[i * 2] = 0x80000000UL;
                        } else if (e == 1) {
                            uint32_t l = msg_strm.read().to_uint64() & 0xffffffffUL;
                            l = ((0x000000ffUL & l) << 24);
                            b.M[i * 2] = l | 0x00800000UL;
                        } else if (e == 2) {
                            uint32_t l = msg_strm.read().to_uint64() & 0xffffffffUL;
                            l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8);
                            b.M[i * 2] = l | 0x00008000UL;
                        } else {
                            uint32_t l = msg_strm.read().to_uint64() & 0xffffffffUL;
                            l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8);
                            b.M[i * 2] = l | 0x00000080UL;
                        }
                        b.M[i * 2 + 1] = 0UL;
                    } else {
                        uint64_t ll = msg_strm.read().to_uint64();
                        uint32_t l = ll & 0xffffffffUL;
                        l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8) |
                            ((0xff000000UL & l) >> 24);
                        b.M[i * 2] = l;
                        l = (ll >> 32) & 0xffffffffUL;
                        uint32_t e = left & 3L;
                        if (e == 0) {
                            b.M[i * 2 + 1] = 0x80000000UL;
                        } else if (e == 1) {
                            l = ((0x000000ffUL & l) << 24);
                            b.M[i * 2 + 1] = l | 0x00800000UL;
                        } else if (e == 2) {
                            l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8);
                            b.M[i * 2 + 1] = l | 0x00008000UL;
                        } else {
                            l = ((0x000000ffUL & l) << 24) | ((0x0000ff00UL & l) << 8) | ((0x00ff0000UL & l) >> 8);
                            b.M[i * 2 + 1] = l | 0x00000080UL;
                        }
                    }
                }
            }

            if (left < 56) {
                b.M[14] = (uint32_t)(0xffffffffUL & (L >> 32));
                b.M[15] = (uint32_t)(0xffffffffUL & (L));
                blk_strm.write(b);
            } else {
                blk_strm.write(b);
                SHA256Block b1;
#pragma HLS array_partition variable = b1.M complete
            LOOP_SHA256_GEN_L_ONLY_BLK_64:
                for (int i = 0; i < 14; ++i) {
#pragma HLS unroll
                    b1.M[i] = 0;
                }
                b1.M[14] = (uint32_t)(0xffffffffUL & (L >> 32));
                b1.M[15] = (uint32_t)(0xffffffffUL & (L));
                blk_strm.write(b1);
            }
        }
    }
    end_nblk_strm.write(true);
}

/***************
 * 单轮迭代 + 就地生成 Wt
 ***************/
inline void sha256_iter1_onfly(uint32_t& a, uint32_t& b, uint32_t& c, uint32_t& d,
                               uint32_t& e, uint32_t& f, uint32_t& g, uint32_t& h,
                               uint32_t W[16], short t, const uint32_t K[]) {
#pragma HLS inline
    uint32_t Wt = (t < 16) ? W[t & 15]
                           : (SSIG1(W[(t - 2) & 15]) + W[(t - 7) & 15]
                              + SSIG0(W[(t - 15) & 15]) + W[(t - 16) & 15]);
    if (t >= 16) W[t & 15] = Wt;

    uint32_t s1  = BSIG1(e);
    uint32_t chv = CH(e, f, g);
    uint32_t s0  = BSIG0(a);
    uint32_t maj = MAJ(a, b, c);

    // 平衡加法树（不使用 bind_op，工具更稳）
    uint32_t sum_h_s1 = h + s1;
    uint32_t sum_ck   = chv + K[t & 63];
    uint32_t t1_pre   = sum_h_s1 + sum_ck;
    uint32_t T1       = t1_pre + Wt;
    uint32_t T2       = s0 + maj;

    h = g;
    g = f;
    f = e;
    e = d + T1;
    d = c;
    c = b;
    b = a;
    a = T1 + T2;
}

/****************
 * Digest on-the-fly
 ****************/
template <int h_width>
void sha256Digest_onfly(hls::stream<SHA256Block>& blk_strm,
                        hls::stream<uint64_t>& nblk_strm,
                        hls::stream<bool>& end_nblk_strm,
                        hls::stream<ap_uint<h_width> >& hash_strm,
                        hls::stream<bool>& end_hash_strm) {
    XF_SECURITY_STATIC_ASSERT((h_width == 256) || (h_width == 224),
                              "Unsupported hash stream width, must be 224 or 256");

    static const uint32_t K[64] = {
        0x428a2f98UL, 0x71374491UL, 0xb5c0fbcfUL, 0xe9b5dba5UL, 0x3956c25bUL, 0x59f111f1UL, 0x923f82a4UL, 0xab1c5ed5UL,
        0xd807aa98UL, 0x12835b01UL, 0x243185beUL, 0x550c7dc3UL, 0x72be5d74UL, 0x80deb1feUL, 0x9bdc06a7UL, 0xc19bf174UL,
        0xe49b69c1UL, 0xefbe4786UL, 0x0fc19dc6UL, 0x240ca1ccUL, 0x2de92c6fUL, 0x4a7484aaUL, 0x5cb0a9dcUL, 0x76f988daUL,
        0x983e5152UL, 0xa831c66dUL, 0xb00327c8UL, 0xbf597fc7UL, 0xc6e00bf3UL, 0xd5a79147UL, 0x06ca6351UL, 0x14292967UL,
        0x27b70a85UL, 0x2e1b2138UL, 0x4d2c6dfcUL, 0x53380d13UL, 0x650a7354UL, 0x766a0abbUL, 0x81c2c92eUL, 0x92722c85UL,
        0xa2bfe8a1UL, 0xa81a664bUL, 0xc24b8b70UL, 0xc76c51a3UL, 0xd192e819UL, 0xd6990624UL, 0xf40e3585UL, 0x106aa070UL,
        0x19a4c116UL, 0x1e376c08UL, 0x2748774cUL, 0x34b0bcb5UL, 0x391c0cb3UL, 0x4ed8aa4aUL, 0x5b9cca4fUL, 0x682e6ff3UL,
        0x748f82eeUL, 0x78a5636fUL, 0x84c87814UL, 0x8cc70208UL, 0x90befffaUL, 0xa4506cebUL, 0xbef9a3f7UL, 0xC67178F2UL};
#pragma HLS array_partition variable = K complete

LOOP_SHA256_DIGEST_MAIN:
    for (bool end_flag = end_nblk_strm.read(); !end_flag; end_flag = end_nblk_strm.read()) {
        uint64_t blk_num = nblk_strm.read();

        uint32_t H[8];
#pragma HLS array_partition variable = H complete
        if (h_width == 224) {
            H[0] = 0xc1059ed8UL; H[1] = 0x367cd507UL; H[2] = 0x3070dd17UL; H[3] = 0xf70e5939UL;
            H[4] = 0xffc00b31UL; H[5] = 0x68581511UL; H[6] = 0x64f98fa7UL; H[7] = 0xbefa4fa4UL;
        } else {
            H[0] = 0x6a09e667UL; H[1] = 0xbb67ae85UL; H[2] = 0x3c6ef372UL; H[3] = 0xa54ff53aUL;
            H[4] = 0x510e527fUL; H[5] = 0x9b05688cUL; H[6] = 0x1f83d9abUL; H[7] = 0x5be0cd19UL;
        }

    LOOP_SHA256_DIGEST_NBLK:
        for (uint64_t n = 0; n < blk_num; ++n) {
#pragma HLS loop_tripcount min = 1 max = 1

            SHA256Block blk = blk_strm.read();
#pragma HLS array_partition variable = blk.M complete

            uint32_t W[16];
#pragma HLS array_partition variable = W complete
        INIT_W_RING:
            for (int t = 0; t < 16; ++t) {
#pragma HLS unroll
                W[t] = blk.M[t];
            }

            uint32_t a = H[0], b = H[1], c = H[2], d = H[3];
            uint32_t e = H[4], f = H[5], g = H[6], h = H[7];

        LOOP_SHA256_UPDATE_64_ROUNDS_1PER:
            for (short t = 0; t < 64; ++t) {
#pragma HLS pipeline II = 1
#pragma HLS dependence variable=W inter false
                sha256_iter1_onfly(a, b, c, d, e, f, g, h, W, t, K);
            }

            H[0] = a + H[0];
            H[1] = b + H[1];
            H[2] = c + H[2];
            H[3] = d + H[3];
            H[4] = e + H[4];
            H[5] = f + H[5];
            H[6] = g + H[6];
            H[7] = h + H[7];
        }

        if (h_width == 224) {
            ap_uint<224> w224;
        LOOP_SHA256_EMIT_H224:
            for (short i = 0; i < sha256_digest_config<true>::numH; ++i) {
#pragma HLS unroll
                uint32_t l = H[i];
                uint8_t t0 = (((l) >> 24) & 0xff);
                uint8_t t1 = (((l) >> 16) & 0xff);
                uint8_t t2 = (((l) >> 8) & 0xff);
                uint8_t t3 = (((l)) & 0xff);
                uint32_t l_little =
                    ((uint32_t)t0) | (((uint32_t)t1) << 8) | (((uint32_t)t2) << 16) | (((uint32_t)t3) << 24);
                w224.range(32 * i + 31, 32 * i) = l_little;
            }
            hash_strm.write(w224);
        } else {
            ap_uint<256> w256;
        LOOP_SHA256_EMIT_H256:
            for (short i = 0; i < sha256_digest_config<false>::numH; ++i) {
#pragma HLS unroll
                uint32_t l = H[i];
                uint8_t t0 = (((l) >> 24) & 0xff);
                uint8_t t1 = (((l) >> 16) & 0xff);
                uint8_t t2 = (((l) >> 8) & 0xff);
                uint8_t t3 = (((l)) & 0xff);
                uint32_t l_little =
                    ((uint32_t)t0) | (((uint32_t)t1) << 8) | (((uint32_t)t2) << 16) | (((uint32_t)t3) << 24);
                w256.range(32 * i + 31, 32 * i) = l_little;
            }
            hash_strm.write(w256);
        }
        end_hash_strm.write(false);
    }
    end_hash_strm.write(true);
}

/***************
 * Top
 ***************/
template <int m_width, int h_width>
inline void sha256_top(hls::stream<ap_uint<m_width> >& msg_strm,
                       hls::stream<ap_uint<64> >& len_strm,
                       hls::stream<bool>& end_len_strm,
                       hls::stream<ap_uint<h_width> >& hash_strm,
                       hls::stream<bool>& end_hash_strm) {
#pragma HLS DATAFLOW
    hls::stream<SHA256Block> blk_strm("blk_strm");
#pragma HLS STREAM variable = blk_strm depth = 32
#pragma HLS RESOURCE variable = blk_strm core = FIFO_LUTRAM

    hls::stream<uint64_t> nblk_strm("nblk_strm");
#pragma HLS STREAM variable = nblk_strm depth = 32
#pragma HLS RESOURCE variable = nblk_strm core = FIFO_LUTRAM

    hls::stream<bool> end_nblk_strm("end_nblk_strm");
#pragma HLS STREAM variable = end_nblk_strm depth = 32
#pragma HLS RESOURCE variable = end_nblk_strm core = FIFO_LUTRAM

    preProcessing(msg_strm, len_strm, end_len_strm, blk_strm, nblk_strm, end_nblk_strm);
    sha256Digest_onfly<h_width>(blk_strm, nblk_strm, end_nblk_strm, hash_strm, end_hash_strm);
}

} // namespace internal

template <int m_width>
void sha224(hls::stream<ap_uint<m_width> >& msg_strm,
            hls::stream<ap_uint<64> >& len_strm,
            hls::stream<bool>& end_len_strm,
            hls::stream<ap_uint<224> >& hash_strm,
            hls::stream<bool>& end_hash_strm) {
    internal::sha256_top(msg_strm, len_strm, end_len_strm, hash_strm, end_hash_strm);
}

template <int m_width>
void sha256(hls::stream<ap_uint<m_width> >& msg_strm,
            hls::stream<ap_uint<64> >& len_strm,
            hls::stream<bool>& end_len_strm,
            hls::stream<ap_uint<256> >& hash_strm,
            hls::stream<bool>& end_hash_strm) {
    internal::sha256_top(msg_strm, len_strm, end_len_strm, hash_strm, end_hash_strm);
}

} // namespace security
} // namespace xf

#undef ROTR
#undef ROTL
#undef SHR
#undef CH
#undef MAJ
#undef BSIG0
#undef BSIG1
#undef SSIG0
#undef SSIG1

#undef _XF_SECURITY_PRINT
#undef _XF_SECURITY_VOID_CAST

#endif // _XF_SECURITY_SHA224_256_HPP_