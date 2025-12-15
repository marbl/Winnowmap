
#define _POSIX_C_SOURCE 200112L

#define HAVE_ALIGNED_ALLOC 1
#define HAVE_POSIX_MEMALIGN 1

#define VERSION "2.4.2"

#define SIZEOF_INT 4

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#define restrict __restrict

#define HAVE_ALTIVEC 0

#define HAVE_NEON 0

#define HAVE_XGETBV 0
#define HAVE_SSE2 0
#define HAVE_SSE41 0


#ifdef MACHINETYPE_pcc64el
#undef  HAVE_ALTIVEC
#define HAVE_ALTIVEC 1
#endif


#ifdef MACHINETYPE_arm64
#undef  HAVE_NEON
#define HAVE_NEON 1
#endif


#ifdef MACHINETYPE_amd64
#undef  HAVE_XGETBV
#define HAVE_XGETBV 1

#undef  HAVE_SSE2
#define HAVE_SSE2 1
#define HAVE_SSE2_MM_SET1_EPI64X 1
#define HAVE_SSE2_MM_SET_EPI64X 1

#undef  HAVE_SSE41
#define HAVE_SSE41 1
#define HAVE_SSE41_MM_EXTRACT_EPI64 1
#define HAVE_SSE41_MM_INSERT_EPI64 1

#undef  HAVE_AVX2
#define HAVE_AVX2 1
#define HAVE_AVX2_MM256_EXTRACT_EPI16 1
#define HAVE_AVX2_MM256_EXTRACT_EPI32 1
#define HAVE_AVX2_MM256_EXTRACT_EPI64 1
#define HAVE_AVX2_MM256_EXTRACT_EPI8 1
#define HAVE_AVX2_MM256_INSERT_EPI16 1
#define HAVE_AVX2_MM256_INSERT_EPI32 1
#define HAVE_AVX2_MM256_INSERT_EPI64 1
#define HAVE_AVX2_MM256_INSERT_EPI8 1
#define HAVE_AVX2_MM256_SET1_EPI64X 1
#define HAVE_AVX2_MM256_SET_EPI64X 1

#define HAVE_AVX512BW 1
#define HAVE_AVX512F 1
#define HAVE_AVX512VBMI 1
#endif

