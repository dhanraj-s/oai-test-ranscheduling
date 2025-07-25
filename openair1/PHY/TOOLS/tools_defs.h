/*
 * Licensed to the OpenAirInterface (OAI) Software Alliance under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The OpenAirInterface Software Alliance licenses this file to You under
 * the OAI Public License, Version 1.1  (the "License"); you may not use this file
 * except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.openairinterface.org/?page_id=698
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *-------------------------------------------------------------------------------
 * For more information about the OpenAirInterface (OAI) Software Alliance:
 *      contact@openairinterface.org
 */

#ifndef __PHY_TOOLS_DEFS__H__
#define __PHY_TOOLS_DEFS__H__

/** @addtogroup _PHY_DSP_TOOLS_


* @{

*/

#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include "PHY/sse_intrin.h"

#include "common/utils/assertions.h"
#include "common/utils/utils.h"
#include "common/utils/LOG/log.h"

#define simd_q15_t simde__m128i
#define shiftright_int16(a,shift) simde_mm_srai_epi16(a,shift)
#define set1_int16(a) simde_mm_set1_epi16(a)
#define mulhi_int16(a,b) simde_mm_mulhrs_epi16 (a,b)
#define mulhi_s1_int16(a,b) simde_mm_slli_epi16(simde_mm_mulhi_epi16(a,b),2)
#define adds_int16(a,b) simde_mm_adds_epi16(a,b)
#define mullo_int16(a,b) simde_mm_mullo_epi16(a,b)

#ifdef __cplusplus
extern "C" {
#endif

#define ALIGNARRAYSIZE(a, b) (((a + b - 1) / b) * b)
#define ALNARS_16_4(a) ALIGNARRAYSIZE(a, 4)

  typedef struct complexd {
    double r;
    double i;
  } cd_t;

  typedef struct complexf {
    float r;
    float i;
  } cf_t;

  typedef struct complex8 {
    int8_t r;
    int8_t i;
  } c8_t;

  typedef struct complex16 {
    int16_t r;
    int16_t i;
  } c16_t;

  typedef struct complex32 {
    int32_t r;
    int32_t i;
  } c32_t;

  typedef struct complex64 {
    int64_t r;
    int64_t i;
  } c64_t;

  typedef struct {
    int dim1;
    int dim2;
    int dim3;
    int dim4;
    uint8_t data[] __attribute__((aligned(32)));
  } fourDimArray_t;

  static inline fourDimArray_t *allocateFourDimArray(int elmtSz, int dim1, int dim2, int dim3, int dim4)
  {
    int sz = elmtSz;
    DevAssert(dim1 > 0);
    sz *= dim1;
    if (dim2) {
      sz *= dim2;
      if (dim3) {
        sz *= dim3;
        if (dim4)
          sz *= dim4;
      }
    }
    fourDimArray_t *tmp = (fourDimArray_t *)malloc16_clear(sizeof(*tmp) + sz);
    AssertFatal(tmp, "no more memory\n");
    *tmp = (fourDimArray_t){dim1, dim2, dim3, dim4};
    return tmp;
  }

#define CheckArrAllocated(workingVar, elementType, ArraY, diM1, diM2, diM3, diM4, resizeAllowed)                           \
  if (!(ArraY))                                                                                                            \
    ArraY = allocateFourDimArray(sizeof(elementType), diM1, diM2, diM3, diM4);                                             \
  else {                                                                                                                   \
    if ((resizeAllowed)                                                                                                    \
        && ((diM1) != (ArraY)->dim1 || (diM2) != (ArraY)->dim2 || (diM3) != (ArraY)->dim3 || (diM4) != (ArraY)->dim4)) {   \
      LOG_I(PHY, "resizing %s to %d/%d/%d/%d\n", #ArraY, (diM1), (diM2), (diM3), (diM4));                                  \
      free(ArraY);                                                                                                         \
      ArraY = allocateFourDimArray(sizeof(elementType), diM1, diM2, diM3, diM4);                                           \
    } else                                                                                                                 \
      DevAssert((diM1) == (ArraY)->dim1 && (diM2) == (ArraY)->dim2 && (diM3) == (ArraY)->dim3 && (diM4) == (ArraY)->dim4); \
  }

#define cast1Darray(workingVar, elementType, ArraY) elementType *workingVar = (elementType *)((ArraY)->data);

#define allocCast1D(workingVar, elementType, ArraY, dim1, resizeAllowed)           \
  CheckArrAllocated(workingVar, elementType, ArraY, dim1, 0, 0, 0, resizeAllowed); \
  cast1Darray(workingVar, elementType, ArraY);

#define cast2Darray(workingVar, elementType, ArraY) \
  elementType(*workingVar)[(ArraY)->dim2] = (elementType(*)[(ArraY)->dim2])((ArraY)->data);

#define allocCast2D(workingVar, elementType, ArraY, dim1, dim2, resizeAllowed)        \
  CheckArrAllocated(workingVar, elementType, ArraY, dim1, dim2, 0, 0, resizeAllowed); \
  cast2Darray(workingVar, elementType, ArraY);

#define cast3Darray(workingVar, elementType, ArraY) \
  elementType(*workingVar)[(ArraY)->dim2][(ArraY)->dim3] = (elementType(*)[(ArraY)->dim2][(ArraY)->dim3])((ArraY)->data);

#define allocCast3D(workingVar, elementType, ArraY, dim1, dim2, dim3, resizeAllowed)     \
  CheckArrAllocated(workingVar, elementType, ArraY, dim1, dim2, dim3, 0, resizeAllowed); \
  cast3Darray(workingVar, elementType, ArraY);

#define cast4Darray(workingVar, elementType, ArraY)                       \
  elementType(*workingVar)[(ArraY)->dim2][(ArraY)->dim3][(ArraY)->dim4] = \
      (elementType(*)[(ArraY)->dim2][(ArraY)->dim3][(ArraY)->dim4])((ArraY)->data);

#define allocCast4D(workingVar, elementType, ArraY, dim1, dim2, dim3, dim4, resizeAllowed)  \
  CheckArrAllocated(workingVar, elementType, ArraY, dim1, dim2, dim3, dim4, resizeAllowed); \
  cast4Darray(workingVar, elementType, ArraY);

#define clearArray(ArraY, elementType) \
  memset((ArraY)->data,                  \
         0,                            \
         sizeof(elementType) * (ArraY)->dim1 * max((ArraY)->dim2, 1) * max((ArraY)->dim3, 1) * max((ArraY)->dim4, 1))

#define squaredMod(a) ((a).r*(a).r + (a).i*(a).i)
#define csum(res, i1, i2) (res).r = (i1).r + (i2).r ; (res).i = (i1).i + (i2).i

  __attribute__((always_inline)) inline c16_t c16conj(const c16_t a) {
    return (c16_t) {
      .r =  a.r,
      .i = (int16_t)-a.i
    };
  }

  __attribute__((always_inline)) inline uint32_t c16toI32(const c16_t a) {
    return *((uint32_t*)&a);
  }

  __attribute__((always_inline)) inline c16_t c16swap(const c16_t a) {
    return (c16_t){
        .r = a.i,
        .i = a.r
    };
  }

  __attribute__((always_inline)) inline uint32_t c16amp2(const c16_t a) {
    return a.r * a.r + a.i * a.i;
  }

  __attribute__((always_inline)) inline c16_t c16add(const c16_t a, const c16_t b)
  {
    return (c16_t){.r = (int16_t)(a.r + b.r), .i = (int16_t)(a.i + b.i)};
  }

  __attribute__((always_inline)) inline c16_t c16sub(const c16_t a, const c16_t b) {
    return (c16_t) {
        .r = (int16_t) (a.r - b.r),
        .i = (int16_t) (a.i - b.i)
    };
  }

  __attribute__((always_inline)) inline c16_t c16Shift(const c16_t a, const int Shift) {
    return (c16_t) {
        .r = (int16_t)(a.r >> Shift),
        .i = (int16_t)(a.i >> Shift)
    };
  }

  __attribute__((always_inline)) inline c16_t c16addShift(const c16_t a, const c16_t b, const int Shift) {
    return (c16_t) {
        .r = (int16_t)((a.r + b.r) >> Shift),
        .i = (int16_t)((a.i + b.i) >> Shift)
    };
  }

  __attribute__((always_inline)) inline c16_t c16mulShift(const c16_t a, const c16_t b, const int Shift) {
    return (c16_t) {
      .r = (int16_t)((a.r * b.r - a.i * b.i) >> Shift),
      .i = (int16_t)((a.r * b.i + a.i * b.r) >> Shift)
    };
  }

  __attribute__((always_inline)) inline c16_t c16mulRealShift(const c16_t a, const int32_t b, const int Shift)
  {
    return (c16_t){.r = (int16_t)((a.r * b) >> Shift), .i = (int16_t)((a.i * b) >> Shift)};
  }
  __attribute__((always_inline)) inline c16_t c16MulConjShift(const c16_t a, const c16_t b, const int Shift)
  {
    return (c16_t) {
      .r = (int16_t)((a.r * b.r + a.i * b.i) >> Shift),
      .i = (int16_t)((a.r * b.i - a.i * b.r) >> Shift)
    };
  }

  __attribute__((always_inline)) inline c16_t c16maddShift(const c16_t a, const c16_t b, c16_t c, const int Shift) {
    return (c16_t) {
      .r = (int16_t)(((a.r * b.r - a.i * b.i ) >> Shift) + c.r),
      .i = (int16_t)(((a.r * b.i + a.i * b.r ) >> Shift) + c.i)
    };
  }

  __attribute__((always_inline)) inline c32_t c32x16mulShift(const c16_t a, const c16_t b, const int Shift) {
    return (c32_t) {
      .r = (a.r * b.r - a.i * b.i) >> Shift,
      .i = (a.r * b.i + a.i * b.r) >> Shift
    };
  }

  __attribute__((always_inline)) inline c16_t c16xmulConstShift(const c16_t a, const int b, const int Shift)
  {
    return (c16_t){.r = (int16_t)((a.r * b) >> Shift), .i = (int16_t)((a.i * b) >> Shift)};
  }

  __attribute__((always_inline)) inline c32_t c32x16maddShift(const c16_t a, const c16_t b, const c32_t c, const int Shift) {
    return (c32_t) {
      .r = ((a.r * b.r - a.i * b.i) >> Shift) + c.r,
      .i = ((a.r * b.i + a.i * b.r) >> Shift) + c.i
    };
  }

  __attribute__((always_inline)) inline c16_t c16x32div(const c32_t a, const int div) {
    return (c16_t) {
      .r = (int16_t)(a.r / div),
      .i = (int16_t)(a.i / div)
    };
  }

  __attribute__((always_inline)) inline cd_t cdMul(const cd_t a, const cd_t b)
  {
    return (cd_t) {
        .r = a.r * b.r - a.i * b.i,
        .i = a.r * b.i + a.i * b.r
    };
  }

  // On N complex numbers
  //   y.r += (x * alpha.r) >> 14
  //   y.i += (x * alpha.i) >> 14
  // See regular C implementation at the end
  static __attribute__((always_inline)) inline void c16multaddVectRealComplex(const int16_t *x,
                                                                              const c16_t *alpha,
                                                                              c16_t *y,
                                                                              const int N)
  {
    // Default implementation for x86
    const int8_t makePairs[32] __attribute__((aligned(32)))={
      0,1,0+16,1+16,
      2,3,2+16,3+16,
      4,5,4+16,5+16,
      6,7,6+16,7+16,
      8,9,8+16,9+16,
      10,11,10+16,11+16,
      12,13,12+16,13+16,
      14,15,14+16,15+16};
    
    simde__m256i alpha256= simde_mm256_set1_epi32(*(int32_t *)alpha);
    simde__m128i *x128=(simde__m128i *)x;
    simde__m128i *y128=(simde__m128i *)y;
    AssertFatal(N%8==0,"Not implemented\n");
    for (int i=0; i<N/8; i++) {
      const simde__m256i xduplicate=simde_mm256_broadcastsi128_si256(*x128);
      const simde__m256i x_duplicate_ordered=simde_mm256_shuffle_epi8(xduplicate,*(simde__m256i*)makePairs);
      const simde__m256i x_mul_alpha_shift15 =simde_mm256_mulhrs_epi16(alpha256, x_duplicate_ordered);
      // Existing multiplication normalization is weird, constant table in alpha need to be doubled
      const simde__m256i x_mul_alpha_x2= simde_mm256_adds_epi16(x_mul_alpha_shift15,x_mul_alpha_shift15);
      *y128 = simde_mm_adds_epi16(simde_mm256_extracti128_si256(x_mul_alpha_x2,0),*y128);
      y128++;
      *y128= simde_mm_adds_epi16(simde_mm256_extracti128_si256(x_mul_alpha_x2,1),*y128);
      y128++;
      x128++;
    }
  }

/*!\fn void multadd_real_vector_complex_scalar(int16_t *x,int16_t *alpha,int16_t *y,uint32_t N)
This function performs componentwise multiplication and accumulation of a complex scalar and a real vector.
@param x Vector input (Q1.15)
@param alpha Scalar input (Q1.15) in the format  |Re0 Im0|
@param y Output (Q1.15) in the format  |Re0  Im0 Re1 Im1|,......,|Re(N-1)  Im(N-1) Re(N-1) Im(N-1)|
@param N Length of x WARNING: N>=8

The function implemented is : \f$\mathbf{y} = y + \alpha\mathbf{x}\f$
*/

  static inline void multadd_real_vector_complex_scalar(const int16_t *x, const c16_t alpha, c16_t *y, const uint32_t N)
  {
    // do 8 multiplications at a time
    simd_q15_t *x_128 = (simd_q15_t *)x, *y_128 = (simd_q15_t *)y;

    //  printf("alpha = %d,%d\n",alpha[0],alpha[1]);
    const simd_q15_t alpha_r_128 = set1_int16(alpha.r);
    const simd_q15_t alpha_i_128 = set1_int16(alpha.i);
    for (uint32_t i = 0; i < N >> 3; i++) {
      const simd_q15_t yr = mulhi_s1_int16(alpha_r_128, x_128[i]);
      const simd_q15_t yi = mulhi_s1_int16(alpha_i_128, x_128[i]);
      const simd_q15_t tmp = simde_mm_loadu_si128(y_128);
      simde_mm_storeu_si128(y_128++, simde_mm_adds_epi16(tmp, simde_mm_unpacklo_epi16(yr, yi)));
      const simd_q15_t tmp2 = simde_mm_loadu_si128(y_128);
      simde_mm_storeu_si128(y_128++, simde_mm_adds_epi16(tmp2, simde_mm_unpackhi_epi16(yr, yi)));
    }
  }

static __attribute__((always_inline)) inline void multadd_real_four_symbols_vector_complex_scalar(const int16_t *x,
                                                                                           const c16_t alpha,
                                                                                           c16_t *y)
{
    // do 8 multiplications at a time
    const simd_q15_t alpha_r_128 = set1_int16(alpha.r);
    const simd_q15_t alpha_i_128 = set1_int16(alpha.i);

    const simd_q15_t *x_128 = (const simd_q15_t *)x;
    // fixme: this loads too much data after pointer x (8 * int16_t instead of 4 * int16_t)
    const simd_q15_t yr = mulhi_s1_int16(alpha_r_128, *x_128);
    const simd_q15_t yi = mulhi_s1_int16(alpha_i_128, *x_128);
    simd_q15_t y_128 = simde_mm_loadu_si128((simd_q15_t *)y);
    y_128 = simde_mm_adds_epi16(y_128, simde_mm_unpacklo_epi16(yr, yi));
    y_128 = simde_mm_adds_epi16(y_128, simde_mm_unpackhi_epi16(yr, yi));

    simde_mm_storeu_si128((simd_q15_t *)y, y_128);
}

// Multiply two vectors of complex int16 and take the most significant bits (shift by 15 in normal case)
// works only with little endian storage (for big endian, modify the srai/ssli at the end)
static __attribute__((always_inline)) inline void mult_complex_vectors(const c16_t *in1,
                                                                       const c16_t *in2,
                                                                       c16_t *out,
                                                                       const int size,
                                                                       const int shift)
{
  const simde__m256i complex_shuffle256 = simde_mm256_set_epi8(29,
                                                               28,
                                                               31,
                                                               30,
                                                               25,
                                                               24,
                                                               27,
                                                               26,
                                                               21,
                                                               20,
                                                               23,
                                                               22,
                                                               17,
                                                               16,
                                                               19,
                                                               18,
                                                               13,
                                                               12,
                                                               15,
                                                               14,
                                                               9,
                                                               8,
                                                               11,
                                                               10,
                                                               5,
                                                               4,
                                                               7,
                                                               6,
                                                               1,
                                                               0,
                                                               3,
                                                               2);
  const simde__m256i conj256 = simde_mm256_set_epi16(-1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1);
  int i;
  // do 8 multiplications at a time
  for (i = 0; i < size - 7; i += 8) {
    const simde__m256i i1 = simde_mm256_loadu_si256((simde__m256i *)(in1 + i));
    const simde__m256i i2 = simde_mm256_loadu_si256((simde__m256i *)(in2 + i));
    const simde__m256i i2swap = simde_mm256_shuffle_epi8(i2, complex_shuffle256);
    const simde__m256i i2conj = simde_mm256_sign_epi16(i2, conj256);
    const simde__m256i re = simde_mm256_madd_epi16(i1, i2conj);
    const simde__m256i im = simde_mm256_madd_epi16(i1, i2swap);
    simde_mm256_storeu_si256(
        (simde__m256i *)(out + i),
        simde_mm256_blend_epi16(simde_mm256_srai_epi32(re, shift), simde_mm256_slli_epi32(im, 16 - shift), 0xAA));
  }
  if (size - i > 4) {
    const simde__m128i i1 = simde_mm_loadu_si128((simde__m128i *)(in1 + i));
    const simde__m128i i2 = simde_mm_loadu_si128((simde__m128i *)(in2 + i));
    const simde__m128i i2swap = simde_mm_shuffle_epi8(i2, *(simde__m128i *)&complex_shuffle256);
    const simde__m128i i2conj = simde_mm_sign_epi16(i2, *(simde__m128i *)&conj256);
    const simde__m128i re = simde_mm_madd_epi16(i1, i2conj);
    const simde__m128i im = simde_mm_madd_epi16(i1, i2swap);
    simde_mm_storeu_si128((simde__m128i *)(out + i),
                          simde_mm_blend_epi16(simde_mm_srai_epi32(re, shift), simde_mm_slli_epi32(im, 16 - shift), 0xAA));
    i += 4;
  }
  for (; i < size; i++)
    out[i] = c16mulShift(in1[i], in2[i], shift);
}
/*!\fn void multadd_complex_vector_real_scalar(int16_t *x,int16_t alpha,int16_t *y,uint8_t zero_flag,uint32_t N)
This function performs componentwise multiplication and accumulation of a real scalar and a complex vector.
@param x Vector input (Q1.15) in the format |Re0 Im0|Re1 Im 1| ...
@param alpha Scalar input (Q1.15) in the format  |Re0|
@param y Output (Q1.15) in the format  |Re0  Im0 Re1 Im1|,......,|Re(N-1)  Im(N-1) Re(N-1) Im(N-1)|
@param zero_flag Set output (y) to zero prior to accumulation
@param N Length of x WARNING: N>=8

The function implemented is : \f$\mathbf{y} = y + \alpha\mathbf{x}\f$
*/
static inline void multadd_complex_vector_real_scalar(const c16_t *x, const int16_t alpha, c16_t *y, uint32_t N)
{
  simd_q15_t alpha_128, *x_128 = (simd_q15_t *)x, *y_128 = (simd_q15_t *)y;
  alpha_128 = set1_int16(alpha);
  const uint32_t num_simd_adds = N / 4;
  const uint32_t num_adds = N % 4;
  for (uint32_t n = 0; n < num_simd_adds; n++) {
    y_128[n] = adds_int16(y_128[n], mulhi_int16(x_128[n], alpha_128));
  }
  for (uint32_t n = 0; n < num_adds; n++) {
    const uint32_t offset = num_simd_adds * 4;
    y[offset + n].r += (x[offset + n].r * alpha) >> 16;
    y[offset + n].i += (x[offset + n].i * alpha) >> 16;
  }
}

static inline void mult_complex_vector_real_scalar(const c16_t *x, const int16_t alpha, c16_t *y, const uint32_t N)
{
  simd_q15_t alpha_128, *x_128 = (simd_q15_t *)x, *y_128 = (simd_q15_t *)y;

  alpha_128 = set1_int16(alpha);
  const uint32_t num_simd_adds = N / 4;
  const uint32_t num_adds = N % 4;
  for (uint32_t n = 0; n < num_simd_adds; n++) {
    y_128[n] = mulhi_int16(x_128[n], alpha_128);
  }
  for (uint32_t n = 0; n < num_adds; n++) {
    const uint32_t offset = num_simd_adds * 4;
    y[offset + n].r = (x[offset + n].r * alpha) >> 16;
    y[offset + n].i = (x[offset + n].i * alpha) >> 16;
  }
}

/*!\fn void init_fft(uint16_t size,uint8_t logsize,uint16_t *rev)
\brief Initialize the FFT engine for a given size
@param size Size of the FFT
@param logsize log2(size)
@param rev Pointer to bit-reversal permutation array
*/

/*!
  Multiply elementwise the complex conjugate of x1 with x2.
  @param x1       - input 1    in the format  |Re0 Im0 Re1 Im1|,......,|Re(N-2)  Im(N-2) Re(N-1) Im(N-1)|
              We assume x1 with a dinamic of 15 bit maximum
  @param x2       - input 2    in the format  |Re0 Im0 Re1 Im1|,......,|Re(N-2)  Im(N-2) Re(N-1) Im(N-1)|
              We assume x2 with a dinamic of 14 bit maximum
  @param y        - output     in the format  |Re0 Im0 Re1 Im1|,......,|Re(N-2)  Im(N-2) Re(N-1) Im(N-1)|
  @param N        - the size f the vectors (this function does N cpx mpy. WARNING: N%4==0;
  @param output_shift  - shift to be applied to generate output
*/
static inline void mult_cpx_conj_vector(const c16_t *x1, const c16_t *x2, c16_t *y, const uint32_t N, int const output_shift)
{
  const simde__m128i *x1_128 = (simde__m128i *)x1;
  const simde__m128i *x2_128 = (simde__m128i *)x2;
  simde__m128i *y_128 = (simde__m128i *)y;

  // SSE compute 4 cpx multiply for each loop
  for (uint32_t i = 0; i < (N >> 2); i++)
    y_128[i] = oai_mm_cpx_mult_conj(x1_128[i], x2_128[i], output_shift);
}

/*!
  Element-wise multiplication and accumulation of two complex vectors x1 and x2.
  @param x1       - input 1    in the format  |Re0 Im0 Re1 Im1|,......,|Re(N-2)  Im(N-2) Re(N-1) Im(N-1)|
              We assume x1 with a dinamic of 15 bit maximum
  @param x2       - input 2    in the format  |Re0 Im0 Re1 Im1|,......,|Re(N-2)  Im(N-2) Re(N-1) Im(N-1)|
              We assume x2 with a dinamic of 14 bit maximum
  @param y        - output     in the format  |Re0 Im0 Re1 Im1|,......,|Re(N-2)  Im(N-2) Re(N-1) Im(N-1)|
  @param zero_flag Set output (y) to zero prior to accumulation
  @param N        - the size f the vectors (this function does N cpx mpy. WARNING: N>=4;
  @param output_shift  - shift to be applied to generate output
*/
static inline void mult_cpx_vector(const c16_t *x1, // Q15
                                   const c16_t *x2, // Q13
                                   c16_t *y,
                                   const uint32_t N,
                                   const int output_shift)
{
  const simde__m128i *x1_128 = (simde__m128i *)x1;
  const simde__m128i *x2_128 = (simde__m128i *)x2;
  simde__m128i *y_128 = (simde__m128i *)y;

  // right shift by 13 while p_a * x0 and 15 while
  //  SSE compute 4 cpx multiply for each loop
  for (uint32_t i = 0; i < (N >> 2); i++) {
    y_128[i] = oai_mm_cpx_mult(x1_128[i], x2_128[i], output_shift);
  }
}

static inline void multadd_cpx_vector(const c16_t *x1, const c16_t *x2, c16_t *y, const uint32_t N, const int output_shift)
{
  const simde__m128i *x1_128 = (simde__m128i *)x1;
  const simde__m128i *x2_128 = (simde__m128i *)x2;
  simde__m128i *y_128 = (simde__m128i *)y;
  // SSE compute 4 cpx multiply for each loop
  for (uint32_t i = 0; i < (N >> 2); i++) {
    simde__m128i result = oai_mm_cpx_mult(x1_128[i], x2_128[i], output_shift);
    y_128[i] = simde_mm_adds_epi16(y_128[i], result);
  }
}

// lte_dfts.c
void init_fft(uint16_t size,
              uint8_t logsize,
              uint16_t *rev);

#define FOREACH_DFTSZ(SZ_DEF) \
  SZ_DEF(12)                  \
  SZ_DEF(24)                  \
  SZ_DEF(36)                  \
  SZ_DEF(48)                  \
  SZ_DEF(60)                  \
  SZ_DEF(64)                  \
  SZ_DEF(72)                  \
  SZ_DEF(96)                  \
  SZ_DEF(108)                 \
  SZ_DEF(120)                 \
  SZ_DEF(128)                 \
  SZ_DEF(144)                 \
  SZ_DEF(180)                 \
  SZ_DEF(192)                 \
  SZ_DEF(216)                 \
  SZ_DEF(240)                 \
  SZ_DEF(256)                 \
  SZ_DEF(288)                 \
  SZ_DEF(300)                 \
  SZ_DEF(324)                 \
  SZ_DEF(360)                 \
  SZ_DEF(384)                 \
  SZ_DEF(432)                 \
  SZ_DEF(480)                 \
  SZ_DEF(512)                 \
  SZ_DEF(540)                 \
  SZ_DEF(576)                 \
  SZ_DEF(600)                 \
  SZ_DEF(648)                 \
  SZ_DEF(720)                 \
  SZ_DEF(768)                 \
  SZ_DEF(864)                 \
  SZ_DEF(900)                 \
  SZ_DEF(960)                 \
  SZ_DEF(972)                 \
  SZ_DEF(1024)                \
  SZ_DEF(1080)                \
  SZ_DEF(1152)                \
  SZ_DEF(1200)                \
  SZ_DEF(1296)                \
  SZ_DEF(1440)                \
  SZ_DEF(1500)                \
  SZ_DEF(1536)                \
  SZ_DEF(1620)                \
  SZ_DEF(1728)                \
  SZ_DEF(1800)                \
  SZ_DEF(1920)                \
  SZ_DEF(1944)                \
  SZ_DEF(2048)                \
  SZ_DEF(2160)                \
  SZ_DEF(2304)                \
  SZ_DEF(2400)                \
  SZ_DEF(2592)                \
  SZ_DEF(2700)                \
  SZ_DEF(2880)                \
  SZ_DEF(2916)                \
  SZ_DEF(3000)                \
  SZ_DEF(3072)                \
  SZ_DEF(3240)                \
  SZ_DEF(4096)                \
  SZ_DEF(6144)                \
  SZ_DEF(8192)                \
  SZ_DEF(12288)               \
  SZ_DEF(18432)               \
  SZ_DEF(24576)               \
  SZ_DEF(36864)               \
  SZ_DEF(49152)               \
  SZ_DEF(98304)

#define FOREACH_IDFTSZ(SZ_DEF) \
  SZ_DEF(64)                   \
  SZ_DEF(128)                  \
  SZ_DEF(256)                  \
  SZ_DEF(512)                  \
  SZ_DEF(768)                  \
  SZ_DEF(1024)                 \
  SZ_DEF(1536)                 \
  SZ_DEF(2048)                 \
  SZ_DEF(3072)                 \
  SZ_DEF(4096)                 \
  SZ_DEF(6144)                 \
  SZ_DEF(8192)                 \
  SZ_DEF(12288)                \
  SZ_DEF(16384)                \
  SZ_DEF(18432)                \
  SZ_DEF(24576)                \
  SZ_DEF(32768)                \
  SZ_DEF(36864)                \
  SZ_DEF(49152)                \
  SZ_DEF(65536)                \
  SZ_DEF(98304)

typedef  void(*dftfunc_t)(uint8_t sizeidx,int16_t *sigF,int16_t *sig,unsigned char scale_flag);
typedef void (*idftfunc_t)(uint8_t sizeidx, int16_t *sigF, int16_t *sig, unsigned char scale_flag);
extern dftfunc_t dft;
extern idftfunc_t idft;
int load_dftslib(void);

#define SZ_ENUM(Sz) DFT_##Sz,
typedef enum dft_size_idx {
  FOREACH_DFTSZ(SZ_ENUM)
  DFT_SIZE_IDXTABLESIZE
}  dft_size_idx_t;

/*******************************************************************
*
* NAME :         get_dft
*
* PARAMETERS :   size of ofdm symbol
*
* RETURN :       function for discrete fourier transform
*
* DESCRIPTION :  get dft function depending of ofdm size
*
*********************************************************************/
#define FIND_ENUM(Sz) \
  case Sz:            \
    return DFT_##Sz;  \
    break;
static inline dft_size_idx_t get_dft(int size)
{
  switch (size) {
    FOREACH_DFTSZ(FIND_ENUM)
    default:
      LOG_E(UTIL, "function get_dft : unsupported DFT size %d\n", size);
      break;
  }
  return DFT_SIZE_IDXTABLESIZE;
}

#define SZ_iENUM(Sz) IDFT_##Sz,
typedef enum idft_size_idx {
  FOREACH_IDFTSZ(SZ_iENUM)
  IDFT_SIZE_IDXTABLESIZE
}  idft_size_idx_t;

#ifdef OAIDFTS_MAIN
typedef void (*adftfunc_t)(int16_t *sigF, int16_t *sig, unsigned char scale_flag);
typedef void (*aidftfunc_t)(int16_t *sigF, int16_t *sig, unsigned char scale_flag);

#define SZ_FUNC(Sz) void dft##Sz(int16_t *x, int16_t *y, uint8_t scale_flag);

FOREACH_DFTSZ(SZ_FUNC)

#define SZ_iFUNC(Sz) void idft##Sz(int16_t *x, int16_t *y, uint8_t scale_flag);

FOREACH_IDFTSZ(SZ_iFUNC)
#define SZ_PTR(Sz) {dft ## Sz,Sz},
struct {
  adftfunc_t func;
  int size;
} const dft_ftab[] = {FOREACH_DFTSZ(SZ_PTR)};

#define SZ_iPTR(Sz)  {idft ## Sz,Sz},
struct {
  adftfunc_t func;
  int size;
} const idft_ftab[] = {FOREACH_IDFTSZ(SZ_iPTR)};

#endif

/*******************************************************************
*
* NAME :         get_idft
*
* PARAMETERS :   size of ofdm symbol
*
* RETURN :       index pointing to the dft func in the dft library
*
* DESCRIPTION :  get idft function depending of ofdm size
*
*********************************************************************/
#define FIND_iENUM(iSz) \
  case iSz:             \
    return IDFT_##iSz;  \
    break;

static inline idft_size_idx_t get_idft(int size)
{
  switch (size) {
    FOREACH_IDFTSZ(FIND_iENUM)
    default:
      AssertFatal(false, "function get_idft : unsupported iDFT size %d\n", size);
      break;
  }
  return IDFT_SIZE_IDXTABLESIZE; // never reached and will trigger assertion in idft function
}


/*!\fn int32_t rotate_cpx_vector(c16_t *x,c16_t *alpha,c16_t *y,uint32_t N,uint16_t output_shift)
This function performs componentwise multiplication of a vector with a complex scalar.
@param x Vector input (Q1.15)  in the format  |Re0  Im0|,......,|Re(N-1) Im(N-1)|
@param alpha Scalar input (Q1.15) in the format  |Re0 Im0|
@param y Output (Q1.15) in the format  |Re0  Im0|,......,|Re(N-1) Im(N-1)|
@param N Length of x WARNING: N>=4
@param output_shift Number of bits to shift output down to Q1.15 (should be 15 for Q1.15 inputs) WARNING: log2_amp>0 can cause overflow!!

The function implemented is : \f$\mathbf{y} = \alpha\mathbf{x}\f$
*/
static inline void rotate_cpx_vector(const c16_t *const x, const c16_t *const alpha, c16_t *y, uint32_t N, uint16_t output_shift)
{
  // multiply a complex vector with a complex value (alpha)
  // stores result in y
  // N is the number of complex numbers
  // output_shift reduces the result of the multiplication by this number of bits
#if defined(__x86_64__) || defined(__i386__)
  if (__builtin_cpu_supports("avx2")) {
    // output is 32 bytes aligned, but not the input

    const c16_t for_re = {alpha->r, (int16_t)-alpha->i};
    const simde__m256i alpha_for_real = simde_mm256_set1_epi32(*(uint32_t *)&for_re);
    const c16_t for_im = {alpha->i, alpha->r};
    const simde__m256i alpha_for_im = simde_mm256_set1_epi32(*(uint32_t *)&for_im);
    const simde__m256i perm_mask = simde_mm256_set_epi8(31,
                                                        30,
                                                        23,
                                                        22,
                                                        29,
                                                        28,
                                                        21,
                                                        20,
                                                        27,
                                                        26,
                                                        19,
                                                        18,
                                                        25,
                                                        24,
                                                        17,
                                                        16,
                                                        15,
                                                        14,
                                                        7,
                                                        6,
                                                        13,
                                                        12,
                                                        5,
                                                        4,
                                                        11,
                                                        10,
                                                        3,
                                                        2,
                                                        9,
                                                        8,
                                                        1,
                                                        0);
    simde__m256i *xd = (simde__m256i *)x;
    const simde__m256i *end = xd + N / 8;
    for (simde__m256i *yd = (simde__m256i *)y; xd < end; yd++, xd++) {
      const simde__m256i y256 = simde_mm256_lddqu_si256(xd);
      const simde__m256i xre = simde_mm256_srai_epi32(simde_mm256_madd_epi16(y256, alpha_for_real), output_shift);
      const simde__m256i xim = simde_mm256_srai_epi32(simde_mm256_madd_epi16(y256, alpha_for_im), output_shift);
      // a bit faster than unpacklo+unpackhi+packs
      const simde__m256i tmp = simde_mm256_packs_epi32(xre, xim);
      simde_mm256_storeu_si256(yd, simde_mm256_shuffle_epi8(tmp, perm_mask));
    }
    c16_t *alpha16 = (c16_t *)alpha, *yLast;
    yLast = ((c16_t *)y) + (N / 8) * 8;
    for (c16_t *xTail = (c16_t *)end; xTail < ((c16_t *)x) + N; xTail++, yLast++) {
      *yLast = c16mulShift(*xTail, *alpha16, output_shift);
    }
  } else {
#endif
    // Multiply elementwise two complex vectors of N elements
    // x        - input 1    in the format  |Re0  Im0 |,......,|Re(N-1) Im(N-1)|
    //            We assume x1 with a dynamic of 15 bit maximum
    //
    // alpha      - input 2    in the format  |Re0 Im0|
    //            We assume x2 with a dynamic of 15 bit maximum
    //
    // y        - output     in the format  |Re0  Im0|,......,|Re(N-1) Im(N-1)|
    //
    // N        - the size f the vectors (this function does N cpx mpy. WARNING: N>=4;
    //
    // log2_amp - increase the output amplitude by a factor 2^log2_amp (default is 0)
    //            WARNING: log2_amp>0 can cause overflow!!

    uint32_t i; // loop counter

    simd_q15_t *y_128, alpha_128;
    int32_t *xd = (int32_t *)x;

    simde__m128i shift = simde_mm_cvtsi32_si128(output_shift);

    ((int16_t *)&alpha_128)[0] = alpha->r;
    ((int16_t *)&alpha_128)[1] = -alpha->i;
    ((int16_t *)&alpha_128)[2] = alpha->i;
    ((int16_t *)&alpha_128)[3] = alpha->r;
    ((int16_t *)&alpha_128)[4] = alpha->r;
    ((int16_t *)&alpha_128)[5] = -alpha->i;
    ((int16_t *)&alpha_128)[6] = alpha->i;
    ((int16_t *)&alpha_128)[7] = alpha->r;
    y_128 = (simd_q15_t *)y;

    for (i = 0; i < N >> 2; i++) {
      y_128[i] = simde_mm_packs_epi32( // pack in 16bit integers with saturation [re im re im re im re im]
          simde_mm_sra_epi32( // shift right by shift in order to  compensate for the input amplitude
              simde_mm_madd_epi16( // complex multiply. result is 32bit [Re Im Re Im]
                  simde_mm_setr_epi32(xd[0 + i * 4], xd[0 + i * 4], xd[1 + i * 4], xd[1 + i * 4]),
                  alpha_128),
              shift),
          simde_mm_sra_epi32( // shift right by shift in order to  compensate for the input amplitude
              simde_mm_madd_epi16( // complex multiply. result is 32bit [Re Im Re Im]
                  simde_mm_setr_epi32(xd[2 + i * 4], xd[2 + i * 4], xd[3 + i * 4], xd[3 + i * 4]),
                  alpha_128),
              shift));
      // print_ints("y_128[0]=", &y_128[0]);
    }
#if defined(__x86__) || defined(__x86_64__)
  }
#endif
}

/*!\fn int32_t sub_cpx_vector16(c16_t *x,c16_t *y, c16_t z, uint32_t N)
This function performs componentwise subsctraction  of complex vectors
@param x Vector input (Q1.15)  in the format  |Re0  Im0 Re1 Im1|,......,|Re(N-2)  Im(N-2) Re(N-1) Im(N-1)|
@param y Output (Q1.15) in the format  |Re0  Im0 Re1 Im1|,......,|Re(N-2)  Im(N-2) Re(N-1) Im(N-1)|
@param N Length of x WARNING: N>=4

The function implemented is : \f$\mathbf{y} = \alpha + \mathbf{x}\f$
*/

static inline int32_t sub_cpx_vector16(const c16_t *x, const c16_t *y, c16_t *z, uint32_t N)
{
  simde__m128i *x_128 = (simde__m128i *)x;
  simde__m128i *y_128 = (simde__m128i *)y;
  simde__m128i *z_128 = (simde__m128i *)z;

  for (uint32_t i = 0; i < (N >> 3); i++)
    z_128[i] = simde_mm_subs_epi16(x_128[i], y_128[i]);
  return (0);
}

/*!\fn int32_t signal_energy(int *,uint32_t);
\brief Computes the signal energy per subcarrier
*/
int32_t signal_energy(int32_t *,uint32_t);

/*!\fn uint32_t signal_energy_nodc(c16_t *,uint32_t);
\brief Computes the signal energy per subcarrier, without DC removal
*/
uint32_t signal_energy_nodc(const c16_t *input, uint32_t length);

int32_t signal_power(int32_t *,uint32_t);
int32_t interference_power(int32_t *,uint32_t);

/*!\fn double signal_energy_fp(double *s_re[2], double *s_im[2],uint32_t, uint32_t,uint32_t);
\brief Computes the signal energy per subcarrier
*/
double signal_energy_fp(double *s_re[2], double *s_im[2], uint32_t nb_antennas, uint32_t length,uint32_t offset);

/*!\fn double signal_energy_fp2(struct complex *, uint32_t);
\brief Computes the signal energy per subcarrier
*/
double signal_energy_fp2(struct complexd *s, uint32_t length);


int32_t iSqrt(int32_t value);
uint8_t log2_approx(uint32_t);
uint8_t log2_approx64(unsigned long long int x);
int16_t invSqrt(int16_t x);
uint32_t angle(struct complex16 perrror);

/// computes the number of factors 2 in x
unsigned char factor2(unsigned int x);

int8_t dB_fixed(uint32_t x);
uint8_t dB_fixed64(uint64_t x);
int8_t dB_fixed2(uint32_t x,uint32_t y);
int16_t dB_fixed_times10(uint32_t x);
int16_t dB_fixed_x10(uint32_t x);

c32_t dot_product(const c16_t *x,
                  const c16_t *y,
                  const uint32_t N, // must be a multiple of 8
                  const int output_shift);

/** @} */

void oai_mm_separate_real_imag_parts(simde__m128i *out_re, simde__m128i *out_im, simde__m128i in0, simde__m128i in1);
void oai_mm256_separate_real_imag_parts(simde__m256i *out_re, simde__m256i *out_im, simde__m256i in0, simde__m256i in1);

void InitSinLUT(void);

// fast calculation of e^{i*2*pi*phase}
// ret.r == cosinus << 14
// ret.i == sinus << 14
c16_t get_sin_cos(double phase);

#ifdef __cplusplus
}
#endif

#endif //__PHY_TOOLS_DEFS__H__
