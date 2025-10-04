/* Provide Declarations */
#include <stdarg.h>
#include <setjmp.h>
#include <limits.h>
#ifdef NEED_CBEAPINT
#include <autopilot_cbe.h>
#else
#define aesl_fopen fopen
#define aesl_freopen freopen
#define aesl_tmpfile tmpfile
#endif
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#ifdef __STRICT_ANSI__
#define inline __inline__
#define typeof __typeof__ 
#endif
#define __isoc99_fscanf fscanf
#define __isoc99_sscanf sscanf
#undef ferror
#undef feof
/* get a declaration for alloca */
#if defined(__CYGWIN__) || defined(__MINGW32__)
#define  alloca(x) __builtin_alloca((x))
#define _alloca(x) __builtin_alloca((x))
#elif defined(__APPLE__)
extern void *__builtin_alloca(unsigned long);
#define alloca(x) __builtin_alloca(x)
#define longjmp _longjmp
#define setjmp _setjmp
#elif defined(__sun__)
#if defined(__sparcv9)
extern void *__builtin_alloca(unsigned long);
#else
extern void *__builtin_alloca(unsigned int);
#endif
#define alloca(x) __builtin_alloca(x)
#elif defined(__FreeBSD__) || defined(__NetBSD__) || defined(__OpenBSD__) || defined(__DragonFly__) || defined(__arm__)
#define alloca(x) __builtin_alloca(x)
#elif defined(_MSC_VER)
#define inline _inline
#define alloca(x) _alloca(x)
#else
#include <alloca.h>
#endif

#ifndef __GNUC__  /* Can only support "linkonce" vars with GCC */
#define __attribute__(X)
#endif

#if defined(__GNUC__) && defined(__APPLE_CC__)
#define __EXTERNAL_WEAK__ __attribute__((weak_import))
#elif defined(__GNUC__)
#define __EXTERNAL_WEAK__ __attribute__((weak))
#else
#define __EXTERNAL_WEAK__
#endif

#if defined(__GNUC__) && (defined(__APPLE_CC__) || defined(__CYGWIN__) || defined(__MINGW32__))
#define __ATTRIBUTE_WEAK__
#elif defined(__GNUC__)
#define __ATTRIBUTE_WEAK__ __attribute__((weak))
#else
#define __ATTRIBUTE_WEAK__
#endif

#if defined(__GNUC__)
#define __HIDDEN__ __attribute__((visibility("hidden")))
#endif

#ifdef __GNUC__
#define LLVM_NAN(NanStr)   __builtin_nan(NanStr)   /* Double */
#define LLVM_NANF(NanStr)  __builtin_nanf(NanStr)  /* Float */
#define LLVM_NANS(NanStr)  __builtin_nans(NanStr)  /* Double */
#define LLVM_NANSF(NanStr) __builtin_nansf(NanStr) /* Float */
#define LLVM_INF           __builtin_inf()         /* Double */
#define LLVM_INFF          __builtin_inff()        /* Float */
#define LLVM_PREFETCH(addr,rw,locality) __builtin_prefetch(addr,rw,locality)
#define __ATTRIBUTE_CTOR__ __attribute__((constructor))
#define __ATTRIBUTE_DTOR__ __attribute__((destructor))
#define LLVM_ASM           __asm__
#else
#define LLVM_NAN(NanStr)   ((double)0.0)           /* Double */
#define LLVM_NANF(NanStr)  0.0F                    /* Float */
#define LLVM_NANS(NanStr)  ((double)0.0)           /* Double */
#define LLVM_NANSF(NanStr) 0.0F                    /* Float */
#define LLVM_INF           ((double)0.0)           /* Double */
#define LLVM_INFF          0.0F                    /* Float */
#define LLVM_PREFETCH(addr,rw,locality)            /* PREFETCH */
#define __ATTRIBUTE_CTOR__
#define __ATTRIBUTE_DTOR__
#define LLVM_ASM(X)
#endif

#if __GNUC__ < 4 /* Old GCC's, or compilers not GCC */ 
#define __builtin_stack_save() 0   /* not implemented */
#define __builtin_stack_restore(X) /* noop */
#endif

#if __GNUC__ && __LP64__ /* 128-bit integer types */
typedef int __attribute__((mode(TI))) llvmInt128;
typedef unsigned __attribute__((mode(TI))) llvmUInt128;
#endif

#define CODE_FOR_MAIN() /* Any target-specific code for main()*/

#ifndef __cplusplus
typedef unsigned char bool;
#endif


/* Support for floating point constants */
typedef unsigned long long ConstantDoubleTy;
typedef unsigned int        ConstantFloatTy;
typedef struct { unsigned long long f1; unsigned short f2; unsigned short pad[3]; } ConstantFP80Ty;
typedef struct { unsigned long long f1; unsigned long long f2; } ConstantFP128Ty;


/* Global Declarations */
/* Helper union for bitcasts */
typedef union {
  unsigned int Int32;
  unsigned long long Int64;
  float Float;
  double Double;
} llvmBitCastUnion;

/* Function Declarations */
double fmod(double, double);
float fmodf(float, float);
long double fmodl(long double, long double);
signed int fir(float *llvm_cbe_y_in, float llvm_cbe_mu, float *llvm_cbe_ref, signed int llvm_cbe_nbTrain, float *llvm_cbe_output, signed int llvm_cbe_totalNumber);
static void aesl_internal_shift_insertion(float *llvm_cbe_shift_reg, float llvm_cbe_new_value);
static float aesl_internal_dotProduct(float *llvm_cbe_row, float *llvm_cbe_column);
static void aesl_internal_weightUpdate(float *llvm_cbe_w, float llvm_cbe_mu, float llvm_cbe_error, float *llvm_cbe_shift_reg);


/* Function Bodies */
static inline int llvm_fcmp_ord(double X, double Y) { return X == X && Y == Y; }
static inline int llvm_fcmp_uno(double X, double Y) { return X != X || Y != Y; }
static inline int llvm_fcmp_ueq(double X, double Y) { return X == Y || llvm_fcmp_uno(X, Y); }
static inline int llvm_fcmp_une(double X, double Y) { return X != Y; }
static inline int llvm_fcmp_ult(double X, double Y) { return X <  Y || llvm_fcmp_uno(X, Y); }
static inline int llvm_fcmp_ugt(double X, double Y) { return X >  Y || llvm_fcmp_uno(X, Y); }
static inline int llvm_fcmp_ule(double X, double Y) { return X <= Y || llvm_fcmp_uno(X, Y); }
static inline int llvm_fcmp_uge(double X, double Y) { return X >= Y || llvm_fcmp_uno(X, Y); }
static inline int llvm_fcmp_oeq(double X, double Y) { return X == Y ; }
static inline int llvm_fcmp_one(double X, double Y) { return X != Y && llvm_fcmp_ord(X, Y); }
static inline int llvm_fcmp_olt(double X, double Y) { return X <  Y ; }
static inline int llvm_fcmp_ogt(double X, double Y) { return X >  Y ; }
static inline int llvm_fcmp_ole(double X, double Y) { return X <= Y ; }
static inline int llvm_fcmp_oge(double X, double Y) { return X >= Y ; }

signed int fir(float *llvm_cbe_y_in, float llvm_cbe_mu, float *llvm_cbe_ref, signed int llvm_cbe_nbTrain, float *llvm_cbe_output, signed int llvm_cbe_totalNumber) {
  static  unsigned long long aesl_llvm_cbe_w_count = 0;
  float llvm_cbe_w[5];    /* Address-exposed local */
  static  unsigned long long aesl_llvm_cbe_shift_reg_count = 0;
  float llvm_cbe_shift_reg[5];    /* Address-exposed local */
  static  unsigned long long aesl_llvm_cbe_1_count = 0;
  static  unsigned long long aesl_llvm_cbe_2_count = 0;
  static  unsigned long long aesl_llvm_cbe_3_count = 0;
  static  unsigned long long aesl_llvm_cbe_4_count = 0;
  static  unsigned long long aesl_llvm_cbe_5_count = 0;
  static  unsigned long long aesl_llvm_cbe_6_count = 0;
  static  unsigned long long aesl_llvm_cbe_7_count = 0;
  static  unsigned long long aesl_llvm_cbe_8_count = 0;
  static  unsigned long long aesl_llvm_cbe_9_count = 0;
  static  unsigned long long aesl_llvm_cbe_10_count = 0;
  static  unsigned long long aesl_llvm_cbe_11_count = 0;
  static  unsigned long long aesl_llvm_cbe_12_count = 0;
  static  unsigned long long aesl_llvm_cbe_13_count = 0;
  static  unsigned long long aesl_llvm_cbe_14_count = 0;
  static  unsigned long long aesl_llvm_cbe_15_count = 0;
  static  unsigned long long aesl_llvm_cbe_16_count = 0;
  static  unsigned long long aesl_llvm_cbe_17_count = 0;
  static  unsigned long long aesl_llvm_cbe_18_count = 0;
  static  unsigned long long aesl_llvm_cbe_19_count = 0;
  static  unsigned long long aesl_llvm_cbe_20_count = 0;
  static  unsigned long long aesl_llvm_cbe_21_count = 0;
  static  unsigned long long aesl_llvm_cbe_22_count = 0;
  static  unsigned long long aesl_llvm_cbe_23_count = 0;
  static  unsigned long long aesl_llvm_cbe_24_count = 0;
  static  unsigned long long aesl_llvm_cbe_25_count = 0;
  static  unsigned long long aesl_llvm_cbe_26_count = 0;
  static  unsigned long long aesl_llvm_cbe_27_count = 0;
  static  unsigned long long aesl_llvm_cbe_or_2e_cond_count = 0;
  bool llvm_cbe_or_2e_cond;
  static  unsigned long long aesl_llvm_cbe_28_count = 0;
  static  unsigned long long aesl_llvm_cbe_or_2e_cond2_count = 0;
  bool llvm_cbe_or_2e_cond2;
  static  unsigned long long aesl_llvm_cbe_29_count = 0;
  static  unsigned long long aesl_llvm_cbe_30_count = 0;
  static  unsigned long long aesl_llvm_cbe_31_count = 0;
   char *llvm_cbe_tmp__1;
  static  unsigned long long aesl_llvm_cbe_32_count = 0;
   char *llvm_cbe_tmp__2;
  static  unsigned long long aesl_llvm_cbe_33_count = 0;
  static  unsigned long long aesl_llvm_cbe_34_count = 0;
   char *llvm_cbe_tmp__3;
  static  unsigned long long aesl_llvm_cbe_35_count = 0;
   char *llvm_cbe_tmp__4;
  static  unsigned long long aesl_llvm_cbe_36_count = 0;
  static  unsigned long long aesl_llvm_cbe_37_count = 0;
  static  unsigned long long aesl_llvm_cbe_38_count = 0;
  static  unsigned long long aesl_llvm_cbe_39_count = 0;
  static  unsigned long long aesl_llvm_cbe_40_count = 0;
  static  unsigned long long aesl_llvm_cbe_41_count = 0;
  static  unsigned long long aesl_llvm_cbe_42_count = 0;
  static  unsigned long long aesl_llvm_cbe_43_count = 0;
  static  unsigned long long aesl_llvm_cbe_44_count = 0;
  static  unsigned long long aesl_llvm_cbe_45_count = 0;
  static  unsigned long long aesl_llvm_cbe_46_count = 0;
  static  unsigned long long aesl_llvm_cbe_47_count = 0;
  static  unsigned long long aesl_llvm_cbe_48_count = 0;
  float *llvm_cbe_tmp__5;
  static  unsigned long long aesl_llvm_cbe_49_count = 0;
  float *llvm_cbe_tmp__6;
  static  unsigned long long aesl_llvm_cbe_50_count = 0;
  static  unsigned long long aesl_llvm_cbe_51_count = 0;
  static  unsigned long long aesl_llvm_cbe_52_count = 0;
  static  unsigned long long aesl_llvm_cbe_53_count = 0;
  static  unsigned long long aesl_llvm_cbe_54_count = 0;
  static  unsigned long long aesl_llvm_cbe_55_count = 0;
  static  unsigned long long aesl_llvm_cbe_56_count = 0;
  static  unsigned long long aesl_llvm_cbe_57_count = 0;
  static  unsigned long long aesl_llvm_cbe_58_count = 0;
  static  unsigned long long aesl_llvm_cbe_59_count = 0;
  static  unsigned long long aesl_llvm_cbe_60_count = 0;
  float *llvm_cbe_tmp__7;
  static  unsigned long long aesl_llvm_cbe_61_count = 0;
  float *llvm_cbe_tmp__8;
  static  unsigned long long aesl_llvm_cbe_62_count = 0;
  static  unsigned long long aesl_llvm_cbe_storemerge4_count = 0;
  unsigned int llvm_cbe_storemerge4;
  unsigned int llvm_cbe_storemerge4__PHI_TEMPORARY;
  static  unsigned long long aesl_llvm_cbe_63_count = 0;
  unsigned long long llvm_cbe_tmp__9;
  static  unsigned long long aesl_llvm_cbe_64_count = 0;
  float *llvm_cbe_tmp__10;
  static  unsigned long long aesl_llvm_cbe_65_count = 0;
  float llvm_cbe_tmp__11;
  static  unsigned long long aesl_llvm_cbe_66_count = 0;
  static  unsigned long long aesl_llvm_cbe_67_count = 0;
  float llvm_cbe_tmp__12;
  static  unsigned long long aesl_llvm_cbe_68_count = 0;
  float *llvm_cbe_tmp__13;
  static  unsigned long long aesl_llvm_cbe_69_count = 0;
  static  unsigned long long aesl_llvm_cbe_70_count = 0;
  float *llvm_cbe_tmp__14;
  static  unsigned long long aesl_llvm_cbe_71_count = 0;
  float llvm_cbe_tmp__15;
  static  unsigned long long aesl_llvm_cbe_72_count = 0;
  float llvm_cbe_tmp__16;
  static  unsigned long long aesl_llvm_cbe_73_count = 0;
  static  unsigned long long aesl_llvm_cbe_74_count = 0;
  unsigned int llvm_cbe_tmp__17;
  static  unsigned long long aesl_llvm_cbe_75_count = 0;
  static  unsigned long long aesl_llvm_cbe_76_count = 0;
  static  unsigned long long aesl_llvm_cbe_77_count = 0;
  static  unsigned long long aesl_llvm_cbe_78_count = 0;
  static  unsigned long long aesl_llvm_cbe_79_count = 0;
  static  unsigned long long aesl_llvm_cbe_80_count = 0;
  static  unsigned long long aesl_llvm_cbe_81_count = 0;
  static  unsigned long long aesl_llvm_cbe_82_count = 0;
  static  unsigned long long aesl_llvm_cbe_83_count = 0;
  static  unsigned long long aesl_llvm_cbe_exitcond7_count = 0;
  static  unsigned long long aesl_llvm_cbe_84_count = 0;
  static  unsigned long long aesl_llvm_cbe_storemerge13_count = 0;
  unsigned int llvm_cbe_storemerge13;
  unsigned int llvm_cbe_storemerge13__PHI_TEMPORARY;
  static  unsigned long long aesl_llvm_cbe_85_count = 0;
  unsigned long long llvm_cbe_tmp__18;
  static  unsigned long long aesl_llvm_cbe_86_count = 0;
  float *llvm_cbe_tmp__19;
  static  unsigned long long aesl_llvm_cbe_87_count = 0;
  float llvm_cbe_tmp__20;
  static  unsigned long long aesl_llvm_cbe_88_count = 0;
  static  unsigned long long aesl_llvm_cbe_89_count = 0;
  float llvm_cbe_tmp__21;
  static  unsigned long long aesl_llvm_cbe_90_count = 0;
  float *llvm_cbe_tmp__22;
  static  unsigned long long aesl_llvm_cbe_91_count = 0;
  static  unsigned long long aesl_llvm_cbe_92_count = 0;
  unsigned int llvm_cbe_tmp__23;
  static  unsigned long long aesl_llvm_cbe_93_count = 0;
  static  unsigned long long aesl_llvm_cbe_94_count = 0;
  static  unsigned long long aesl_llvm_cbe_95_count = 0;
  static  unsigned long long aesl_llvm_cbe_96_count = 0;
  static  unsigned long long aesl_llvm_cbe_97_count = 0;
  static  unsigned long long aesl_llvm_cbe_98_count = 0;
  static  unsigned long long aesl_llvm_cbe_99_count = 0;
  static  unsigned long long aesl_llvm_cbe_100_count = 0;
  static  unsigned long long aesl_llvm_cbe_exitcond_count = 0;
  static  unsigned long long aesl_llvm_cbe_101_count = 0;
  static  unsigned long long aesl_llvm_cbe_102_count = 0;
  static  unsigned long long aesl_llvm_cbe_103_count = 0;
  static  unsigned long long aesl_llvm_cbe_104_count = 0;
  static  unsigned long long aesl_llvm_cbe_105_count = 0;
  unsigned int llvm_cbe_tmp__24;
  unsigned int llvm_cbe_tmp__24__PHI_TEMPORARY;
  static  unsigned long long aesl_llvm_cbe_106_count = 0;

const char* AESL_DEBUG_TRACE = getenv("DEBUG_TRACE");
if (AESL_DEBUG_TRACE)
printf("\n\{ BEGIN @fir\n");
if (AESL_DEBUG_TRACE)
printf("\n  %%or.cond = or i1 %%1, %%2, !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_or_2e_cond_count);
  llvm_cbe_or_2e_cond = (bool )((((llvm_cbe_y_in) == (((float *)/*NULL*/0))) | ((llvm_cbe_ref) == (((float *)/*NULL*/0))))&1);
if (AESL_DEBUG_TRACE)
printf("\nor.cond = 0x%X\n", llvm_cbe_or_2e_cond);
if (AESL_DEBUG_TRACE)
printf("\n  %%or.cond2 = or i1 %%or.cond, %%3, !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_or_2e_cond2_count);
  llvm_cbe_or_2e_cond2 = (bool )((llvm_cbe_or_2e_cond | ((llvm_cbe_output) == (((float *)/*NULL*/0))))&1);
if (AESL_DEBUG_TRACE)
printf("\nor.cond2 = 0x%X\n", llvm_cbe_or_2e_cond2);
  if (llvm_cbe_or_2e_cond2) {
    llvm_cbe_tmp__24__PHI_TEMPORARY = (unsigned int )4294967295u;   /* for PHI node */
    goto llvm_cbe_tmp__25;
  } else {
    goto llvm_cbe_tmp__26;
  }

llvm_cbe_tmp__26:
if (AESL_DEBUG_TRACE)
printf("\n  %%5 = bitcast [5 x float]* %%w to i8*, !dbg !9 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_31_count);
  llvm_cbe_tmp__1 = ( char *)(( char *)(&llvm_cbe_w));
if (AESL_DEBUG_TRACE)
printf("\n  %%6 = call i8* @memset(i8* %%5, i32 0, i64 20 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_32_count);
  ( char *)memset(( char *)llvm_cbe_tmp__1, 0u, 20ull);
if (AESL_DEBUG_TRACE) {
printf("\nArgument  = 0x%X",0u);
printf("\nArgument  = 0x%I64X",20ull);
printf("\nReturn  = 0x%X",llvm_cbe_tmp__2);
}
if (AESL_DEBUG_TRACE)
printf("\n  %%7 = bitcast [5 x float]* %%shift_reg to i8*, !dbg !9 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_34_count);
  llvm_cbe_tmp__3 = ( char *)(( char *)(&llvm_cbe_shift_reg));
if (AESL_DEBUG_TRACE)
printf("\n  %%8 = call i8* @memset(i8* %%7, i32 0, i64 20 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_35_count);
  ( char *)memset(( char *)llvm_cbe_tmp__3, 0u, 20ull);
if (AESL_DEBUG_TRACE) {
printf("\nArgument  = 0x%X",0u);
printf("\nArgument  = 0x%I64X",20ull);
printf("\nReturn  = 0x%X",llvm_cbe_tmp__4);
}
  if ((((signed int )llvm_cbe_nbTrain) > ((signed int )0u))) {
    goto llvm_cbe__2e_lr_2e_ph6;
  } else {
    goto llvm_cbe__2e_preheader;
  }

llvm_cbe__2e_lr_2e_ph6:
if (AESL_DEBUG_TRACE)
printf("\n  %%10 = getelementptr inbounds [5 x float]* %%shift_reg, i64 0, i64 0, !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_48_count);
  llvm_cbe_tmp__5 = (float *)(&llvm_cbe_shift_reg[(((signed long long )0ull))
#ifdef AESL_BC_SIM
 % 5
#endif
]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%11 = getelementptr inbounds [5 x float]* %%w, i64 0, i64 0, !dbg !7 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_49_count);
  llvm_cbe_tmp__6 = (float *)(&llvm_cbe_w[(((signed long long )0ull))
#ifdef AESL_BC_SIM
 % 5
#endif
]);
if (AESL_DEBUG_TRACE) {
}
  llvm_cbe_storemerge4__PHI_TEMPORARY = (unsigned int )0u;   /* for PHI node */
  goto llvm_cbe_tmp__27;

llvm_cbe__2e_preheader:
  if ((((signed int )llvm_cbe_nbTrain) < ((signed int )llvm_cbe_totalNumber))) {
    goto llvm_cbe__2e_lr_2e_ph;
  } else {
    goto llvm_cbe__2e__crit_edge;
  }

llvm_cbe__2e_lr_2e_ph:
if (AESL_DEBUG_TRACE)
printf("\n  %%13 = getelementptr inbounds [5 x float]* %%shift_reg, i64 0, i64 0, !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_60_count);
  llvm_cbe_tmp__7 = (float *)(&llvm_cbe_shift_reg[(((signed long long )0ull))
#ifdef AESL_BC_SIM
 % 5
#endif
]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%14 = getelementptr inbounds [5 x float]* %%w, i64 0, i64 0, !dbg !7 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_61_count);
  llvm_cbe_tmp__8 = (float *)(&llvm_cbe_w[(((signed long long )0ull))
#ifdef AESL_BC_SIM
 % 5
#endif
]);
if (AESL_DEBUG_TRACE) {
}
  llvm_cbe_storemerge13__PHI_TEMPORARY = (unsigned int )llvm_cbe_nbTrain;   /* for PHI node */
  goto llvm_cbe_tmp__28;

  do {     /* Syntactic loop '' to make GCC happy */
llvm_cbe_tmp__27:
if (AESL_DEBUG_TRACE)
printf("\n  %%storemerge4 = phi i32 [ 0, %%.lr.ph6 ], [ %%24, %%15  for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_storemerge4_count);
  llvm_cbe_storemerge4 = (unsigned int )llvm_cbe_storemerge4__PHI_TEMPORARY;
if (AESL_DEBUG_TRACE) {
printf("\nstoremerge4 = 0x%X",llvm_cbe_storemerge4);
printf("\n = 0x%X",0u);
printf("\n = 0x%X",llvm_cbe_tmp__17);
}
if (AESL_DEBUG_TRACE)
printf("\n  %%16 = sext i32 %%storemerge4 to i64, !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_63_count);
  llvm_cbe_tmp__9 = (unsigned long long )((signed long long )(signed int )llvm_cbe_storemerge4);
if (AESL_DEBUG_TRACE)
printf("\n = 0x%I64X\n", llvm_cbe_tmp__9);
if (AESL_DEBUG_TRACE)
printf("\n  %%17 = getelementptr inbounds float* %%y_in, i64 %%16, !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_64_count);
  llvm_cbe_tmp__10 = (float *)(&llvm_cbe_y_in[(((signed long long )llvm_cbe_tmp__9))]);
if (AESL_DEBUG_TRACE) {
printf("\n = 0x%I64X",((signed long long )llvm_cbe_tmp__9));
}
if (AESL_DEBUG_TRACE)
printf("\n  %%18 = load float* %%17, align 4, !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_65_count);
  llvm_cbe_tmp__11 = (float )*llvm_cbe_tmp__10;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__11, *(int*)(&llvm_cbe_tmp__11));
if (AESL_DEBUG_TRACE)
printf("\n  call fastcc void @aesl_internal_shift_insertion(float* %%10, float %%18), !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_66_count);
  aesl_internal_shift_insertion((float *)llvm_cbe_tmp__5, llvm_cbe_tmp__11);
if (AESL_DEBUG_TRACE) {
printf("\nArgument  = %f,  0x%x",llvm_cbe_tmp__11, *(int*)(&llvm_cbe_tmp__11));
}
if (AESL_DEBUG_TRACE)
printf("\n  %%19 = call fastcc float @aesl_internal_dotProduct(float* %%11, float* %%10), !dbg !7 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_67_count);
  llvm_cbe_tmp__12 = (float )aesl_internal_dotProduct((float *)llvm_cbe_tmp__6, (float *)llvm_cbe_tmp__5);
if (AESL_DEBUG_TRACE) {
printf("\nReturn  = %f",llvm_cbe_tmp__12);
}
if (AESL_DEBUG_TRACE)
printf("\n  %%20 = getelementptr inbounds float* %%output, i64 %%16, !dbg !7 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_68_count);
  llvm_cbe_tmp__13 = (float *)(&llvm_cbe_output[(((signed long long )llvm_cbe_tmp__9))]);
if (AESL_DEBUG_TRACE) {
printf("\n = 0x%I64X",((signed long long )llvm_cbe_tmp__9));
}
if (AESL_DEBUG_TRACE)
printf("\n  store float %%19, float* %%20, align 4, !dbg !7 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_69_count);
  *llvm_cbe_tmp__13 = llvm_cbe_tmp__12;
if (AESL_DEBUG_TRACE)
printf("\n = %f\n", llvm_cbe_tmp__12);
if (AESL_DEBUG_TRACE)
printf("\n  %%21 = getelementptr inbounds float* %%ref, i64 %%16, !dbg !6 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_70_count);
  llvm_cbe_tmp__14 = (float *)(&llvm_cbe_ref[(((signed long long )llvm_cbe_tmp__9))]);
if (AESL_DEBUG_TRACE) {
printf("\n = 0x%I64X",((signed long long )llvm_cbe_tmp__9));
}
if (AESL_DEBUG_TRACE)
printf("\n  %%22 = load float* %%21, align 4, !dbg !6 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_71_count);
  llvm_cbe_tmp__15 = (float )*llvm_cbe_tmp__14;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__15, *(int*)(&llvm_cbe_tmp__15));
if (AESL_DEBUG_TRACE)
printf("\n  %%23 = fsub float %%22, %%19, !dbg !6 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_72_count);
  llvm_cbe_tmp__16 = (float )((float )(llvm_cbe_tmp__15 - llvm_cbe_tmp__12));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__16, *(int*)(&llvm_cbe_tmp__16));
if (AESL_DEBUG_TRACE)
printf("\n  call fastcc void @aesl_internal_weightUpdate(float* %%11, float %%mu, float %%23, float* %%10), !dbg !6 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_73_count);
  aesl_internal_weightUpdate((float *)llvm_cbe_tmp__6, llvm_cbe_mu, llvm_cbe_tmp__16, (float *)llvm_cbe_tmp__5);
if (AESL_DEBUG_TRACE) {
printf("\nArgument mu = %f,  0x%x",llvm_cbe_mu, *(int*)(&llvm_cbe_mu));
printf("\nArgument  = %f,  0x%x",llvm_cbe_tmp__16, *(int*)(&llvm_cbe_tmp__16));
}
if (AESL_DEBUG_TRACE)
printf("\n  %%24 = add nsw i32 %%storemerge4, 1, !dbg !10 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_74_count);
  llvm_cbe_tmp__17 = (unsigned int )((unsigned int )(llvm_cbe_storemerge4&4294967295ull)) + ((unsigned int )(1u&4294967295ull));
if (AESL_DEBUG_TRACE)
printf("\n = 0x%X\n", ((unsigned int )(llvm_cbe_tmp__17&4294967295ull)));
  if (((llvm_cbe_tmp__17&4294967295U) == (llvm_cbe_nbTrain&4294967295U))) {
    goto llvm_cbe__2e_preheader;
  } else {
    llvm_cbe_storemerge4__PHI_TEMPORARY = (unsigned int )llvm_cbe_tmp__17;   /* for PHI node */
    goto llvm_cbe_tmp__27;
  }

  } while (1); /* end of syntactic loop '' */
  do {     /* Syntactic loop '' to make GCC happy */
llvm_cbe_tmp__28:
if (AESL_DEBUG_TRACE)
printf("\n  %%storemerge13 = phi i32 [ %%nbTrain, %%.lr.ph ], [ %%31, %%25  for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_storemerge13_count);
  llvm_cbe_storemerge13 = (unsigned int )llvm_cbe_storemerge13__PHI_TEMPORARY;
if (AESL_DEBUG_TRACE) {
printf("\nstoremerge13 = 0x%X",llvm_cbe_storemerge13);
printf("\nnbTrain = 0x%X",llvm_cbe_nbTrain);
printf("\n = 0x%X",llvm_cbe_tmp__23);
}
if (AESL_DEBUG_TRACE)
printf("\n  %%26 = sext i32 %%storemerge13 to i64, !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_85_count);
  llvm_cbe_tmp__18 = (unsigned long long )((signed long long )(signed int )llvm_cbe_storemerge13);
if (AESL_DEBUG_TRACE)
printf("\n = 0x%I64X\n", llvm_cbe_tmp__18);
if (AESL_DEBUG_TRACE)
printf("\n  %%27 = getelementptr inbounds float* %%y_in, i64 %%26, !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_86_count);
  llvm_cbe_tmp__19 = (float *)(&llvm_cbe_y_in[(((signed long long )llvm_cbe_tmp__18))]);
if (AESL_DEBUG_TRACE) {
printf("\n = 0x%I64X",((signed long long )llvm_cbe_tmp__18));
}
if (AESL_DEBUG_TRACE)
printf("\n  %%28 = load float* %%27, align 4, !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_87_count);
  llvm_cbe_tmp__20 = (float )*llvm_cbe_tmp__19;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__20, *(int*)(&llvm_cbe_tmp__20));
if (AESL_DEBUG_TRACE)
printf("\n  call fastcc void @aesl_internal_shift_insertion(float* %%13, float %%28), !dbg !5 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_88_count);
  aesl_internal_shift_insertion((float *)llvm_cbe_tmp__7, llvm_cbe_tmp__20);
if (AESL_DEBUG_TRACE) {
printf("\nArgument  = %f,  0x%x",llvm_cbe_tmp__20, *(int*)(&llvm_cbe_tmp__20));
}
if (AESL_DEBUG_TRACE)
printf("\n  %%29 = call fastcc float @aesl_internal_dotProduct(float* %%14, float* %%13), !dbg !7 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_89_count);
  llvm_cbe_tmp__21 = (float )aesl_internal_dotProduct((float *)llvm_cbe_tmp__8, (float *)llvm_cbe_tmp__7);
if (AESL_DEBUG_TRACE) {
printf("\nReturn  = %f",llvm_cbe_tmp__21);
}
if (AESL_DEBUG_TRACE)
printf("\n  %%30 = getelementptr inbounds float* %%output, i64 %%26, !dbg !7 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_90_count);
  llvm_cbe_tmp__22 = (float *)(&llvm_cbe_output[(((signed long long )llvm_cbe_tmp__18))]);
if (AESL_DEBUG_TRACE) {
printf("\n = 0x%I64X",((signed long long )llvm_cbe_tmp__18));
}
if (AESL_DEBUG_TRACE)
printf("\n  store float %%29, float* %%30, align 4, !dbg !7 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_91_count);
  *llvm_cbe_tmp__22 = llvm_cbe_tmp__21;
if (AESL_DEBUG_TRACE)
printf("\n = %f\n", llvm_cbe_tmp__21);
if (AESL_DEBUG_TRACE)
printf("\n  %%31 = add nsw i32 %%storemerge13, 1, !dbg !10 for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_92_count);
  llvm_cbe_tmp__23 = (unsigned int )((unsigned int )(llvm_cbe_storemerge13&4294967295ull)) + ((unsigned int )(1u&4294967295ull));
if (AESL_DEBUG_TRACE)
printf("\n = 0x%X\n", ((unsigned int )(llvm_cbe_tmp__23&4294967295ull)));
  if (((llvm_cbe_tmp__23&4294967295U) == (llvm_cbe_totalNumber&4294967295U))) {
    goto llvm_cbe__2e__crit_edge;
  } else {
    llvm_cbe_storemerge13__PHI_TEMPORARY = (unsigned int )llvm_cbe_tmp__23;   /* for PHI node */
    goto llvm_cbe_tmp__28;
  }

  } while (1); /* end of syntactic loop '' */
llvm_cbe__2e__crit_edge:
  llvm_cbe_tmp__24__PHI_TEMPORARY = (unsigned int )0u;   /* for PHI node */
  goto llvm_cbe_tmp__25;

llvm_cbe_tmp__25:
if (AESL_DEBUG_TRACE)
printf("\n  %%33 = phi i32 [ 0, %%._crit_edge ], [ -1, %%0  for 0x%I64xth hint within @fir  --> \n", ++aesl_llvm_cbe_105_count);
  llvm_cbe_tmp__24 = (unsigned int )llvm_cbe_tmp__24__PHI_TEMPORARY;
if (AESL_DEBUG_TRACE) {
printf("\n = 0x%X",llvm_cbe_tmp__24);
printf("\n = 0x%X",0u);
printf("\n = 0x%X",4294967295u);
}
  if (AESL_DEBUG_TRACE)
      printf("\nEND @fir}\n");
  return llvm_cbe_tmp__24;
}


static void aesl_internal_shift_insertion(float *llvm_cbe_shift_reg, float llvm_cbe_new_value) {
  static  unsigned long long aesl_llvm_cbe_107_count = 0;
  static  unsigned long long aesl_llvm_cbe_108_count = 0;
  static  unsigned long long aesl_llvm_cbe_109_count = 0;
  static  unsigned long long aesl_llvm_cbe_110_count = 0;
  static  unsigned long long aesl_llvm_cbe_111_count = 0;
  static  unsigned long long aesl_llvm_cbe_112_count = 0;
  static  unsigned long long aesl_llvm_cbe_113_count = 0;
  static  unsigned long long aesl_llvm_cbe_114_count = 0;
  static  unsigned long long aesl_llvm_cbe_115_count = 0;
  static  unsigned long long aesl_llvm_cbe_116_count = 0;
  static  unsigned long long aesl_llvm_cbe_117_count = 0;
  static  unsigned long long aesl_llvm_cbe_118_count = 0;
  static  unsigned long long aesl_llvm_cbe_119_count = 0;
  static  unsigned long long aesl_llvm_cbe_120_count = 0;
  float *llvm_cbe_tmp__29;
  static  unsigned long long aesl_llvm_cbe_121_count = 0;
  float llvm_cbe_tmp__30;
  static  unsigned long long aesl_llvm_cbe_122_count = 0;
  float *llvm_cbe_tmp__31;
  static  unsigned long long aesl_llvm_cbe_123_count = 0;
  static  unsigned long long aesl_llvm_cbe_124_count = 0;
  static  unsigned long long aesl_llvm_cbe_125_count = 0;
  static  unsigned long long aesl_llvm_cbe_126_count = 0;
  static  unsigned long long aesl_llvm_cbe_127_count = 0;
  static  unsigned long long aesl_llvm_cbe_128_count = 0;
  static  unsigned long long aesl_llvm_cbe_129_count = 0;
  float *llvm_cbe_tmp__32;
  static  unsigned long long aesl_llvm_cbe_130_count = 0;
  float llvm_cbe_tmp__33;
  static  unsigned long long aesl_llvm_cbe_131_count = 0;
  static  unsigned long long aesl_llvm_cbe_132_count = 0;
  static  unsigned long long aesl_llvm_cbe_133_count = 0;
  static  unsigned long long aesl_llvm_cbe_134_count = 0;
  static  unsigned long long aesl_llvm_cbe_135_count = 0;
  static  unsigned long long aesl_llvm_cbe_136_count = 0;
  static  unsigned long long aesl_llvm_cbe_137_count = 0;
  float *llvm_cbe_tmp__34;
  static  unsigned long long aesl_llvm_cbe_138_count = 0;
  float llvm_cbe_tmp__35;
  static  unsigned long long aesl_llvm_cbe_139_count = 0;
  static  unsigned long long aesl_llvm_cbe_140_count = 0;
  static  unsigned long long aesl_llvm_cbe_141_count = 0;
  static  unsigned long long aesl_llvm_cbe_142_count = 0;
  static  unsigned long long aesl_llvm_cbe_143_count = 0;
  static  unsigned long long aesl_llvm_cbe_144_count = 0;
  static  unsigned long long aesl_llvm_cbe_145_count = 0;
  float llvm_cbe_tmp__36;
  static  unsigned long long aesl_llvm_cbe_146_count = 0;
  static  unsigned long long aesl_llvm_cbe_147_count = 0;
  static  unsigned long long aesl_llvm_cbe_148_count = 0;
  static  unsigned long long aesl_llvm_cbe_149_count = 0;
  static  unsigned long long aesl_llvm_cbe_150_count = 0;
  static  unsigned long long aesl_llvm_cbe_151_count = 0;
  static  unsigned long long aesl_llvm_cbe_152_count = 0;
  static  unsigned long long aesl_llvm_cbe_153_count = 0;

const char* AESL_DEBUG_TRACE = getenv("DEBUG_TRACE");
if (AESL_DEBUG_TRACE)
printf("\n\{ BEGIN @aesl_internal_shift_insertion\n");
if (AESL_DEBUG_TRACE)
printf("\n  %%1 = getelementptr inbounds float* %%shift_reg, i64 3, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_120_count);
  llvm_cbe_tmp__29 = (float *)(&llvm_cbe_shift_reg[(((signed long long )3ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%2 = load float* %%1, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_121_count);
  llvm_cbe_tmp__30 = (float )*llvm_cbe_tmp__29;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__30, *(int*)(&llvm_cbe_tmp__30));
if (AESL_DEBUG_TRACE)
printf("\n  %%3 = getelementptr inbounds float* %%shift_reg, i64 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_122_count);
  llvm_cbe_tmp__31 = (float *)(&llvm_cbe_shift_reg[(((signed long long )4ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  store float %%2, float* %%3, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_123_count);
  *llvm_cbe_tmp__31 = llvm_cbe_tmp__30;
if (AESL_DEBUG_TRACE)
printf("\n = %f\n", llvm_cbe_tmp__30);
if (AESL_DEBUG_TRACE)
printf("\n  %%4 = getelementptr inbounds float* %%shift_reg, i64 2, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_129_count);
  llvm_cbe_tmp__32 = (float *)(&llvm_cbe_shift_reg[(((signed long long )2ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%5 = load float* %%4, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_130_count);
  llvm_cbe_tmp__33 = (float )*llvm_cbe_tmp__32;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__33, *(int*)(&llvm_cbe_tmp__33));
if (AESL_DEBUG_TRACE)
printf("\n  store float %%5, float* %%1, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_131_count);
  *llvm_cbe_tmp__29 = llvm_cbe_tmp__33;
if (AESL_DEBUG_TRACE)
printf("\n = %f\n", llvm_cbe_tmp__33);
if (AESL_DEBUG_TRACE)
printf("\n  %%6 = getelementptr inbounds float* %%shift_reg, i64 1, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_137_count);
  llvm_cbe_tmp__34 = (float *)(&llvm_cbe_shift_reg[(((signed long long )1ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%7 = load float* %%6, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_138_count);
  llvm_cbe_tmp__35 = (float )*llvm_cbe_tmp__34;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__35, *(int*)(&llvm_cbe_tmp__35));
if (AESL_DEBUG_TRACE)
printf("\n  store float %%7, float* %%4, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_139_count);
  *llvm_cbe_tmp__32 = llvm_cbe_tmp__35;
if (AESL_DEBUG_TRACE)
printf("\n = %f\n", llvm_cbe_tmp__35);
if (AESL_DEBUG_TRACE)
printf("\n  %%8 = load float* %%shift_reg, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_145_count);
  llvm_cbe_tmp__36 = (float )*llvm_cbe_shift_reg;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__36, *(int*)(&llvm_cbe_tmp__36));
if (AESL_DEBUG_TRACE)
printf("\n  store float %%8, float* %%6, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_146_count);
  *llvm_cbe_tmp__34 = llvm_cbe_tmp__36;
if (AESL_DEBUG_TRACE)
printf("\n = %f\n", llvm_cbe_tmp__36);
if (AESL_DEBUG_TRACE)
printf("\n  store float %%new_value, float* %%shift_reg, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_shift_insertion  --> \n", ++aesl_llvm_cbe_152_count);
  *llvm_cbe_shift_reg = llvm_cbe_new_value;
if (AESL_DEBUG_TRACE)
printf("\nnew_value = %f\n", llvm_cbe_new_value);
  if (AESL_DEBUG_TRACE)
      printf("\nEND @aesl_internal_shift_insertion}\n");
  return;
}


static float aesl_internal_dotProduct(float *llvm_cbe_row, float *llvm_cbe_column) {
  static  unsigned long long aesl_llvm_cbe_154_count = 0;
  static  unsigned long long aesl_llvm_cbe_155_count = 0;
  static  unsigned long long aesl_llvm_cbe_156_count = 0;
  static  unsigned long long aesl_llvm_cbe_157_count = 0;
  static  unsigned long long aesl_llvm_cbe_158_count = 0;
  static  unsigned long long aesl_llvm_cbe_159_count = 0;
  static  unsigned long long aesl_llvm_cbe_160_count = 0;
  static  unsigned long long aesl_llvm_cbe_161_count = 0;
  static  unsigned long long aesl_llvm_cbe_162_count = 0;
  static  unsigned long long aesl_llvm_cbe_163_count = 0;
  static  unsigned long long aesl_llvm_cbe_164_count = 0;
  static  unsigned long long aesl_llvm_cbe_165_count = 0;
  static  unsigned long long aesl_llvm_cbe_166_count = 0;
  static  unsigned long long aesl_llvm_cbe_167_count = 0;
  static  unsigned long long aesl_llvm_cbe_168_count = 0;
  float llvm_cbe_tmp__37;
  static  unsigned long long aesl_llvm_cbe_169_count = 0;
  float llvm_cbe_tmp__38;
  static  unsigned long long aesl_llvm_cbe_170_count = 0;
  float llvm_cbe_tmp__39;
  static  unsigned long long aesl_llvm_cbe_171_count = 0;
  float llvm_cbe_tmp__40;
  static  unsigned long long aesl_llvm_cbe_172_count = 0;
  static  unsigned long long aesl_llvm_cbe_173_count = 0;
  static  unsigned long long aesl_llvm_cbe_174_count = 0;
  static  unsigned long long aesl_llvm_cbe_175_count = 0;
  static  unsigned long long aesl_llvm_cbe_176_count = 0;
  static  unsigned long long aesl_llvm_cbe_177_count = 0;
  static  unsigned long long aesl_llvm_cbe_178_count = 0;
  static  unsigned long long aesl_llvm_cbe_179_count = 0;
  static  unsigned long long aesl_llvm_cbe_180_count = 0;
  float *llvm_cbe_tmp__41;
  static  unsigned long long aesl_llvm_cbe_181_count = 0;
  float llvm_cbe_tmp__42;
  static  unsigned long long aesl_llvm_cbe_182_count = 0;
  float *llvm_cbe_tmp__43;
  static  unsigned long long aesl_llvm_cbe_183_count = 0;
  float llvm_cbe_tmp__44;
  static  unsigned long long aesl_llvm_cbe_184_count = 0;
  float llvm_cbe_tmp__45;
  static  unsigned long long aesl_llvm_cbe_185_count = 0;
  float llvm_cbe_tmp__46;
  static  unsigned long long aesl_llvm_cbe_186_count = 0;
  static  unsigned long long aesl_llvm_cbe_187_count = 0;
  static  unsigned long long aesl_llvm_cbe_188_count = 0;
  static  unsigned long long aesl_llvm_cbe_189_count = 0;
  static  unsigned long long aesl_llvm_cbe_190_count = 0;
  static  unsigned long long aesl_llvm_cbe_191_count = 0;
  static  unsigned long long aesl_llvm_cbe_192_count = 0;
  static  unsigned long long aesl_llvm_cbe_193_count = 0;
  static  unsigned long long aesl_llvm_cbe_194_count = 0;
  float *llvm_cbe_tmp__47;
  static  unsigned long long aesl_llvm_cbe_195_count = 0;
  float llvm_cbe_tmp__48;
  static  unsigned long long aesl_llvm_cbe_196_count = 0;
  float *llvm_cbe_tmp__49;
  static  unsigned long long aesl_llvm_cbe_197_count = 0;
  float llvm_cbe_tmp__50;
  static  unsigned long long aesl_llvm_cbe_198_count = 0;
  float llvm_cbe_tmp__51;
  static  unsigned long long aesl_llvm_cbe_199_count = 0;
  float llvm_cbe_tmp__52;
  static  unsigned long long aesl_llvm_cbe_200_count = 0;
  static  unsigned long long aesl_llvm_cbe_201_count = 0;
  static  unsigned long long aesl_llvm_cbe_202_count = 0;
  static  unsigned long long aesl_llvm_cbe_203_count = 0;
  static  unsigned long long aesl_llvm_cbe_204_count = 0;
  static  unsigned long long aesl_llvm_cbe_205_count = 0;
  static  unsigned long long aesl_llvm_cbe_206_count = 0;
  static  unsigned long long aesl_llvm_cbe_207_count = 0;
  static  unsigned long long aesl_llvm_cbe_208_count = 0;
  float *llvm_cbe_tmp__53;
  static  unsigned long long aesl_llvm_cbe_209_count = 0;
  float llvm_cbe_tmp__54;
  static  unsigned long long aesl_llvm_cbe_210_count = 0;
  float *llvm_cbe_tmp__55;
  static  unsigned long long aesl_llvm_cbe_211_count = 0;
  float llvm_cbe_tmp__56;
  static  unsigned long long aesl_llvm_cbe_212_count = 0;
  float llvm_cbe_tmp__57;
  static  unsigned long long aesl_llvm_cbe_213_count = 0;
  float llvm_cbe_tmp__58;
  static  unsigned long long aesl_llvm_cbe_214_count = 0;
  static  unsigned long long aesl_llvm_cbe_215_count = 0;
  static  unsigned long long aesl_llvm_cbe_216_count = 0;
  static  unsigned long long aesl_llvm_cbe_217_count = 0;
  static  unsigned long long aesl_llvm_cbe_218_count = 0;
  static  unsigned long long aesl_llvm_cbe_219_count = 0;
  static  unsigned long long aesl_llvm_cbe_220_count = 0;
  static  unsigned long long aesl_llvm_cbe_221_count = 0;
  static  unsigned long long aesl_llvm_cbe_222_count = 0;
  float *llvm_cbe_tmp__59;
  static  unsigned long long aesl_llvm_cbe_223_count = 0;
  float llvm_cbe_tmp__60;
  static  unsigned long long aesl_llvm_cbe_224_count = 0;
  float *llvm_cbe_tmp__61;
  static  unsigned long long aesl_llvm_cbe_225_count = 0;
  float llvm_cbe_tmp__62;
  static  unsigned long long aesl_llvm_cbe_226_count = 0;
  float llvm_cbe_tmp__63;
  static  unsigned long long aesl_llvm_cbe_227_count = 0;
  float llvm_cbe_tmp__64;
  static  unsigned long long aesl_llvm_cbe_228_count = 0;
  static  unsigned long long aesl_llvm_cbe_229_count = 0;
  static  unsigned long long aesl_llvm_cbe_230_count = 0;
  static  unsigned long long aesl_llvm_cbe_231_count = 0;
  static  unsigned long long aesl_llvm_cbe_232_count = 0;
  static  unsigned long long aesl_llvm_cbe_233_count = 0;
  static  unsigned long long aesl_llvm_cbe_234_count = 0;
  static  unsigned long long aesl_llvm_cbe_235_count = 0;
  static  unsigned long long aesl_llvm_cbe_236_count = 0;

const char* AESL_DEBUG_TRACE = getenv("DEBUG_TRACE");
if (AESL_DEBUG_TRACE)
printf("\n\{ BEGIN @aesl_internal_dotProduct\n");
if (AESL_DEBUG_TRACE)
printf("\n  %%1 = load float* %%row, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_168_count);
  llvm_cbe_tmp__37 = (float )*llvm_cbe_row;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__37, *(int*)(&llvm_cbe_tmp__37));
if (AESL_DEBUG_TRACE)
printf("\n  %%2 = load float* %%column, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_169_count);
  llvm_cbe_tmp__38 = (float )*llvm_cbe_column;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__38, *(int*)(&llvm_cbe_tmp__38));
if (AESL_DEBUG_TRACE)
printf("\n  %%3 = fmul float %%1, %%2, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_170_count);
  llvm_cbe_tmp__39 = (float )((float )(llvm_cbe_tmp__37 * llvm_cbe_tmp__38));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__39, *(int*)(&llvm_cbe_tmp__39));
if (AESL_DEBUG_TRACE)
printf("\n  %%4 = fadd float %%3, 0.000000e+00, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_171_count);
  llvm_cbe_tmp__40 = (float )((float )(llvm_cbe_tmp__39 + 0x0p0));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__40, *(int*)(&llvm_cbe_tmp__40));
if (AESL_DEBUG_TRACE)
printf("\n  %%5 = getelementptr inbounds float* %%row, i64 1, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_180_count);
  llvm_cbe_tmp__41 = (float *)(&llvm_cbe_row[(((signed long long )1ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%6 = load float* %%5, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_181_count);
  llvm_cbe_tmp__42 = (float )*llvm_cbe_tmp__41;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__42, *(int*)(&llvm_cbe_tmp__42));
if (AESL_DEBUG_TRACE)
printf("\n  %%7 = getelementptr inbounds float* %%column, i64 1, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_182_count);
  llvm_cbe_tmp__43 = (float *)(&llvm_cbe_column[(((signed long long )1ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%8 = load float* %%7, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_183_count);
  llvm_cbe_tmp__44 = (float )*llvm_cbe_tmp__43;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__44, *(int*)(&llvm_cbe_tmp__44));
if (AESL_DEBUG_TRACE)
printf("\n  %%9 = fmul float %%6, %%8, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_184_count);
  llvm_cbe_tmp__45 = (float )((float )(llvm_cbe_tmp__42 * llvm_cbe_tmp__44));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__45, *(int*)(&llvm_cbe_tmp__45));
if (AESL_DEBUG_TRACE)
printf("\n  %%10 = fadd float %%4, %%9, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_185_count);
  llvm_cbe_tmp__46 = (float )((float )(llvm_cbe_tmp__40 + llvm_cbe_tmp__45));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__46, *(int*)(&llvm_cbe_tmp__46));
if (AESL_DEBUG_TRACE)
printf("\n  %%11 = getelementptr inbounds float* %%row, i64 2, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_194_count);
  llvm_cbe_tmp__47 = (float *)(&llvm_cbe_row[(((signed long long )2ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%12 = load float* %%11, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_195_count);
  llvm_cbe_tmp__48 = (float )*llvm_cbe_tmp__47;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__48, *(int*)(&llvm_cbe_tmp__48));
if (AESL_DEBUG_TRACE)
printf("\n  %%13 = getelementptr inbounds float* %%column, i64 2, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_196_count);
  llvm_cbe_tmp__49 = (float *)(&llvm_cbe_column[(((signed long long )2ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%14 = load float* %%13, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_197_count);
  llvm_cbe_tmp__50 = (float )*llvm_cbe_tmp__49;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__50, *(int*)(&llvm_cbe_tmp__50));
if (AESL_DEBUG_TRACE)
printf("\n  %%15 = fmul float %%12, %%14, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_198_count);
  llvm_cbe_tmp__51 = (float )((float )(llvm_cbe_tmp__48 * llvm_cbe_tmp__50));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__51, *(int*)(&llvm_cbe_tmp__51));
if (AESL_DEBUG_TRACE)
printf("\n  %%16 = fadd float %%10, %%15, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_199_count);
  llvm_cbe_tmp__52 = (float )((float )(llvm_cbe_tmp__46 + llvm_cbe_tmp__51));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__52, *(int*)(&llvm_cbe_tmp__52));
if (AESL_DEBUG_TRACE)
printf("\n  %%17 = getelementptr inbounds float* %%row, i64 3, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_208_count);
  llvm_cbe_tmp__53 = (float *)(&llvm_cbe_row[(((signed long long )3ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%18 = load float* %%17, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_209_count);
  llvm_cbe_tmp__54 = (float )*llvm_cbe_tmp__53;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__54, *(int*)(&llvm_cbe_tmp__54));
if (AESL_DEBUG_TRACE)
printf("\n  %%19 = getelementptr inbounds float* %%column, i64 3, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_210_count);
  llvm_cbe_tmp__55 = (float *)(&llvm_cbe_column[(((signed long long )3ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%20 = load float* %%19, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_211_count);
  llvm_cbe_tmp__56 = (float )*llvm_cbe_tmp__55;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__56, *(int*)(&llvm_cbe_tmp__56));
if (AESL_DEBUG_TRACE)
printf("\n  %%21 = fmul float %%18, %%20, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_212_count);
  llvm_cbe_tmp__57 = (float )((float )(llvm_cbe_tmp__54 * llvm_cbe_tmp__56));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__57, *(int*)(&llvm_cbe_tmp__57));
if (AESL_DEBUG_TRACE)
printf("\n  %%22 = fadd float %%16, %%21, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_213_count);
  llvm_cbe_tmp__58 = (float )((float )(llvm_cbe_tmp__52 + llvm_cbe_tmp__57));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__58, *(int*)(&llvm_cbe_tmp__58));
if (AESL_DEBUG_TRACE)
printf("\n  %%23 = getelementptr inbounds float* %%row, i64 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_222_count);
  llvm_cbe_tmp__59 = (float *)(&llvm_cbe_row[(((signed long long )4ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%24 = load float* %%23, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_223_count);
  llvm_cbe_tmp__60 = (float )*llvm_cbe_tmp__59;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__60, *(int*)(&llvm_cbe_tmp__60));
if (AESL_DEBUG_TRACE)
printf("\n  %%25 = getelementptr inbounds float* %%column, i64 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_224_count);
  llvm_cbe_tmp__61 = (float *)(&llvm_cbe_column[(((signed long long )4ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%26 = load float* %%25, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_225_count);
  llvm_cbe_tmp__62 = (float )*llvm_cbe_tmp__61;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__62, *(int*)(&llvm_cbe_tmp__62));
if (AESL_DEBUG_TRACE)
printf("\n  %%27 = fmul float %%24, %%26, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_226_count);
  llvm_cbe_tmp__63 = (float )((float )(llvm_cbe_tmp__60 * llvm_cbe_tmp__62));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__63, *(int*)(&llvm_cbe_tmp__63));
if (AESL_DEBUG_TRACE)
printf("\n  %%28 = fadd float %%22, %%27, !dbg !5 for 0x%I64xth hint within @aesl_internal_dotProduct  --> \n", ++aesl_llvm_cbe_227_count);
  llvm_cbe_tmp__64 = (float )((float )(llvm_cbe_tmp__58 + llvm_cbe_tmp__63));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__64, *(int*)(&llvm_cbe_tmp__64));
  if (AESL_DEBUG_TRACE)
      printf("\nEND @aesl_internal_dotProduct}\n");
  return llvm_cbe_tmp__64;
}


static void aesl_internal_weightUpdate(float *llvm_cbe_w, float llvm_cbe_mu, float llvm_cbe_error, float *llvm_cbe_shift_reg) {
  static  unsigned long long aesl_llvm_cbe_237_count = 0;
  static  unsigned long long aesl_llvm_cbe_238_count = 0;
  static  unsigned long long aesl_llvm_cbe_239_count = 0;
  static  unsigned long long aesl_llvm_cbe_240_count = 0;
  static  unsigned long long aesl_llvm_cbe_241_count = 0;
  static  unsigned long long aesl_llvm_cbe_242_count = 0;
  static  unsigned long long aesl_llvm_cbe_243_count = 0;
  static  unsigned long long aesl_llvm_cbe_244_count = 0;
  static  unsigned long long aesl_llvm_cbe_245_count = 0;
  static  unsigned long long aesl_llvm_cbe_246_count = 0;
  static  unsigned long long aesl_llvm_cbe_247_count = 0;
  static  unsigned long long aesl_llvm_cbe_248_count = 0;
  static  unsigned long long aesl_llvm_cbe_249_count = 0;
  static  unsigned long long aesl_llvm_cbe_250_count = 0;
  static  unsigned long long aesl_llvm_cbe_251_count = 0;
  static  unsigned long long aesl_llvm_cbe_252_count = 0;
  float llvm_cbe_tmp__65;
  static  unsigned long long aesl_llvm_cbe_253_count = 0;
  float llvm_cbe_tmp__66;
  static  unsigned long long aesl_llvm_cbe_254_count = 0;
  float llvm_cbe_tmp__67;
  static  unsigned long long aesl_llvm_cbe_255_count = 0;
  float llvm_cbe_tmp__68;
  static  unsigned long long aesl_llvm_cbe_256_count = 0;
  float llvm_cbe_tmp__69;
  static  unsigned long long aesl_llvm_cbe_257_count = 0;
  static  unsigned long long aesl_llvm_cbe_258_count = 0;
  static  unsigned long long aesl_llvm_cbe_259_count = 0;
  static  unsigned long long aesl_llvm_cbe_260_count = 0;
  static  unsigned long long aesl_llvm_cbe_261_count = 0;
  static  unsigned long long aesl_llvm_cbe_262_count = 0;
  static  unsigned long long aesl_llvm_cbe_263_count = 0;
  float *llvm_cbe_tmp__70;
  static  unsigned long long aesl_llvm_cbe_264_count = 0;
  float llvm_cbe_tmp__71;
  static  unsigned long long aesl_llvm_cbe_265_count = 0;
  float llvm_cbe_tmp__72;
  static  unsigned long long aesl_llvm_cbe_266_count = 0;
  float *llvm_cbe_tmp__73;
  static  unsigned long long aesl_llvm_cbe_267_count = 0;
  float llvm_cbe_tmp__74;
  static  unsigned long long aesl_llvm_cbe_268_count = 0;
  float llvm_cbe_tmp__75;
  static  unsigned long long aesl_llvm_cbe_269_count = 0;
  static  unsigned long long aesl_llvm_cbe_270_count = 0;
  static  unsigned long long aesl_llvm_cbe_271_count = 0;
  static  unsigned long long aesl_llvm_cbe_272_count = 0;
  static  unsigned long long aesl_llvm_cbe_273_count = 0;
  static  unsigned long long aesl_llvm_cbe_274_count = 0;
  static  unsigned long long aesl_llvm_cbe_275_count = 0;
  float *llvm_cbe_tmp__76;
  static  unsigned long long aesl_llvm_cbe_276_count = 0;
  float llvm_cbe_tmp__77;
  static  unsigned long long aesl_llvm_cbe_277_count = 0;
  float llvm_cbe_tmp__78;
  static  unsigned long long aesl_llvm_cbe_278_count = 0;
  float *llvm_cbe_tmp__79;
  static  unsigned long long aesl_llvm_cbe_279_count = 0;
  float llvm_cbe_tmp__80;
  static  unsigned long long aesl_llvm_cbe_280_count = 0;
  float llvm_cbe_tmp__81;
  static  unsigned long long aesl_llvm_cbe_281_count = 0;
  static  unsigned long long aesl_llvm_cbe_282_count = 0;
  static  unsigned long long aesl_llvm_cbe_283_count = 0;
  static  unsigned long long aesl_llvm_cbe_284_count = 0;
  static  unsigned long long aesl_llvm_cbe_285_count = 0;
  static  unsigned long long aesl_llvm_cbe_286_count = 0;
  static  unsigned long long aesl_llvm_cbe_287_count = 0;
  float *llvm_cbe_tmp__82;
  static  unsigned long long aesl_llvm_cbe_288_count = 0;
  float llvm_cbe_tmp__83;
  static  unsigned long long aesl_llvm_cbe_289_count = 0;
  float llvm_cbe_tmp__84;
  static  unsigned long long aesl_llvm_cbe_290_count = 0;
  float *llvm_cbe_tmp__85;
  static  unsigned long long aesl_llvm_cbe_291_count = 0;
  float llvm_cbe_tmp__86;
  static  unsigned long long aesl_llvm_cbe_292_count = 0;
  float llvm_cbe_tmp__87;
  static  unsigned long long aesl_llvm_cbe_293_count = 0;
  static  unsigned long long aesl_llvm_cbe_294_count = 0;
  static  unsigned long long aesl_llvm_cbe_295_count = 0;
  static  unsigned long long aesl_llvm_cbe_296_count = 0;
  static  unsigned long long aesl_llvm_cbe_297_count = 0;
  static  unsigned long long aesl_llvm_cbe_298_count = 0;
  static  unsigned long long aesl_llvm_cbe_299_count = 0;
  float *llvm_cbe_tmp__88;
  static  unsigned long long aesl_llvm_cbe_300_count = 0;
  float llvm_cbe_tmp__89;
  static  unsigned long long aesl_llvm_cbe_301_count = 0;
  float llvm_cbe_tmp__90;
  static  unsigned long long aesl_llvm_cbe_302_count = 0;
  float *llvm_cbe_tmp__91;
  static  unsigned long long aesl_llvm_cbe_303_count = 0;
  float llvm_cbe_tmp__92;
  static  unsigned long long aesl_llvm_cbe_304_count = 0;
  float llvm_cbe_tmp__93;
  static  unsigned long long aesl_llvm_cbe_305_count = 0;
  static  unsigned long long aesl_llvm_cbe_306_count = 0;
  static  unsigned long long aesl_llvm_cbe_307_count = 0;
  static  unsigned long long aesl_llvm_cbe_308_count = 0;
  static  unsigned long long aesl_llvm_cbe_309_count = 0;
  static  unsigned long long aesl_llvm_cbe_310_count = 0;
  static  unsigned long long aesl_llvm_cbe_311_count = 0;

const char* AESL_DEBUG_TRACE = getenv("DEBUG_TRACE");
if (AESL_DEBUG_TRACE)
printf("\n\{ BEGIN @aesl_internal_weightUpdate\n");
if (AESL_DEBUG_TRACE)
printf("\n  %%1 = fmul float %%mu, %%error, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_252_count);
  llvm_cbe_tmp__65 = (float )((float )(llvm_cbe_mu * llvm_cbe_error));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__65, *(int*)(&llvm_cbe_tmp__65));
if (AESL_DEBUG_TRACE)
printf("\n  %%2 = load float* %%shift_reg, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_253_count);
  llvm_cbe_tmp__66 = (float )*llvm_cbe_shift_reg;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__66, *(int*)(&llvm_cbe_tmp__66));
if (AESL_DEBUG_TRACE)
printf("\n  %%3 = fmul float %%1, %%2, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_254_count);
  llvm_cbe_tmp__67 = (float )((float )(llvm_cbe_tmp__65 * llvm_cbe_tmp__66));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__67, *(int*)(&llvm_cbe_tmp__67));
if (AESL_DEBUG_TRACE)
printf("\n  %%4 = load float* %%w, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_255_count);
  llvm_cbe_tmp__68 = (float )*llvm_cbe_w;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__68, *(int*)(&llvm_cbe_tmp__68));
if (AESL_DEBUG_TRACE)
printf("\n  %%5 = fadd float %%4, %%3, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_256_count);
  llvm_cbe_tmp__69 = (float )((float )(llvm_cbe_tmp__68 + llvm_cbe_tmp__67));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__69, *(int*)(&llvm_cbe_tmp__69));
if (AESL_DEBUG_TRACE)
printf("\n  store float %%5, float* %%w, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_257_count);
  *llvm_cbe_w = llvm_cbe_tmp__69;
if (AESL_DEBUG_TRACE)
printf("\n = %f\n", llvm_cbe_tmp__69);
if (AESL_DEBUG_TRACE)
printf("\n  %%6 = getelementptr inbounds float* %%shift_reg, i64 1, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_263_count);
  llvm_cbe_tmp__70 = (float *)(&llvm_cbe_shift_reg[(((signed long long )1ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%7 = load float* %%6, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_264_count);
  llvm_cbe_tmp__71 = (float )*llvm_cbe_tmp__70;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__71, *(int*)(&llvm_cbe_tmp__71));
if (AESL_DEBUG_TRACE)
printf("\n  %%8 = fmul float %%1, %%7, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_265_count);
  llvm_cbe_tmp__72 = (float )((float )(llvm_cbe_tmp__65 * llvm_cbe_tmp__71));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__72, *(int*)(&llvm_cbe_tmp__72));
if (AESL_DEBUG_TRACE)
printf("\n  %%9 = getelementptr inbounds float* %%w, i64 1, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_266_count);
  llvm_cbe_tmp__73 = (float *)(&llvm_cbe_w[(((signed long long )1ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%10 = load float* %%9, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_267_count);
  llvm_cbe_tmp__74 = (float )*llvm_cbe_tmp__73;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__74, *(int*)(&llvm_cbe_tmp__74));
if (AESL_DEBUG_TRACE)
printf("\n  %%11 = fadd float %%10, %%8, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_268_count);
  llvm_cbe_tmp__75 = (float )((float )(llvm_cbe_tmp__74 + llvm_cbe_tmp__72));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__75, *(int*)(&llvm_cbe_tmp__75));
if (AESL_DEBUG_TRACE)
printf("\n  store float %%11, float* %%9, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_269_count);
  *llvm_cbe_tmp__73 = llvm_cbe_tmp__75;
if (AESL_DEBUG_TRACE)
printf("\n = %f\n", llvm_cbe_tmp__75);
if (AESL_DEBUG_TRACE)
printf("\n  %%12 = getelementptr inbounds float* %%shift_reg, i64 2, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_275_count);
  llvm_cbe_tmp__76 = (float *)(&llvm_cbe_shift_reg[(((signed long long )2ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%13 = load float* %%12, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_276_count);
  llvm_cbe_tmp__77 = (float )*llvm_cbe_tmp__76;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__77, *(int*)(&llvm_cbe_tmp__77));
if (AESL_DEBUG_TRACE)
printf("\n  %%14 = fmul float %%1, %%13, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_277_count);
  llvm_cbe_tmp__78 = (float )((float )(llvm_cbe_tmp__65 * llvm_cbe_tmp__77));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__78, *(int*)(&llvm_cbe_tmp__78));
if (AESL_DEBUG_TRACE)
printf("\n  %%15 = getelementptr inbounds float* %%w, i64 2, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_278_count);
  llvm_cbe_tmp__79 = (float *)(&llvm_cbe_w[(((signed long long )2ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%16 = load float* %%15, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_279_count);
  llvm_cbe_tmp__80 = (float )*llvm_cbe_tmp__79;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__80, *(int*)(&llvm_cbe_tmp__80));
if (AESL_DEBUG_TRACE)
printf("\n  %%17 = fadd float %%16, %%14, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_280_count);
  llvm_cbe_tmp__81 = (float )((float )(llvm_cbe_tmp__80 + llvm_cbe_tmp__78));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__81, *(int*)(&llvm_cbe_tmp__81));
if (AESL_DEBUG_TRACE)
printf("\n  store float %%17, float* %%15, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_281_count);
  *llvm_cbe_tmp__79 = llvm_cbe_tmp__81;
if (AESL_DEBUG_TRACE)
printf("\n = %f\n", llvm_cbe_tmp__81);
if (AESL_DEBUG_TRACE)
printf("\n  %%18 = getelementptr inbounds float* %%shift_reg, i64 3, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_287_count);
  llvm_cbe_tmp__82 = (float *)(&llvm_cbe_shift_reg[(((signed long long )3ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%19 = load float* %%18, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_288_count);
  llvm_cbe_tmp__83 = (float )*llvm_cbe_tmp__82;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__83, *(int*)(&llvm_cbe_tmp__83));
if (AESL_DEBUG_TRACE)
printf("\n  %%20 = fmul float %%1, %%19, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_289_count);
  llvm_cbe_tmp__84 = (float )((float )(llvm_cbe_tmp__65 * llvm_cbe_tmp__83));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__84, *(int*)(&llvm_cbe_tmp__84));
if (AESL_DEBUG_TRACE)
printf("\n  %%21 = getelementptr inbounds float* %%w, i64 3, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_290_count);
  llvm_cbe_tmp__85 = (float *)(&llvm_cbe_w[(((signed long long )3ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%22 = load float* %%21, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_291_count);
  llvm_cbe_tmp__86 = (float )*llvm_cbe_tmp__85;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__86, *(int*)(&llvm_cbe_tmp__86));
if (AESL_DEBUG_TRACE)
printf("\n  %%23 = fadd float %%22, %%20, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_292_count);
  llvm_cbe_tmp__87 = (float )((float )(llvm_cbe_tmp__86 + llvm_cbe_tmp__84));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__87, *(int*)(&llvm_cbe_tmp__87));
if (AESL_DEBUG_TRACE)
printf("\n  store float %%23, float* %%21, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_293_count);
  *llvm_cbe_tmp__85 = llvm_cbe_tmp__87;
if (AESL_DEBUG_TRACE)
printf("\n = %f\n", llvm_cbe_tmp__87);
if (AESL_DEBUG_TRACE)
printf("\n  %%24 = getelementptr inbounds float* %%shift_reg, i64 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_299_count);
  llvm_cbe_tmp__88 = (float *)(&llvm_cbe_shift_reg[(((signed long long )4ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%25 = load float* %%24, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_300_count);
  llvm_cbe_tmp__89 = (float )*llvm_cbe_tmp__88;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__89, *(int*)(&llvm_cbe_tmp__89));
if (AESL_DEBUG_TRACE)
printf("\n  %%26 = fmul float %%1, %%25, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_301_count);
  llvm_cbe_tmp__90 = (float )((float )(llvm_cbe_tmp__65 * llvm_cbe_tmp__89));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__90, *(int*)(&llvm_cbe_tmp__90));
if (AESL_DEBUG_TRACE)
printf("\n  %%27 = getelementptr inbounds float* %%w, i64 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_302_count);
  llvm_cbe_tmp__91 = (float *)(&llvm_cbe_w[(((signed long long )4ull))]);
if (AESL_DEBUG_TRACE) {
}
if (AESL_DEBUG_TRACE)
printf("\n  %%28 = load float* %%27, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_303_count);
  llvm_cbe_tmp__92 = (float )*llvm_cbe_tmp__91;
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__92, *(int*)(&llvm_cbe_tmp__92));
if (AESL_DEBUG_TRACE)
printf("\n  %%29 = fadd float %%28, %%26, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_304_count);
  llvm_cbe_tmp__93 = (float )((float )(llvm_cbe_tmp__92 + llvm_cbe_tmp__90));
if (AESL_DEBUG_TRACE)
printf("\n = %f,  0x%x\n", llvm_cbe_tmp__93, *(int*)(&llvm_cbe_tmp__93));
if (AESL_DEBUG_TRACE)
printf("\n  store float %%29, float* %%27, align 4, !dbg !5 for 0x%I64xth hint within @aesl_internal_weightUpdate  --> \n", ++aesl_llvm_cbe_305_count);
  *llvm_cbe_tmp__91 = llvm_cbe_tmp__93;
if (AESL_DEBUG_TRACE)
printf("\n = %f\n", llvm_cbe_tmp__93);
  if (AESL_DEBUG_TRACE)
      printf("\nEND @aesl_internal_weightUpdate}\n");
  return;
}

