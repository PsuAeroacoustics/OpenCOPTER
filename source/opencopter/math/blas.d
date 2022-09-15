module opencopter.math.blas;


import core.stdc.config;
import core.stdc.stdarg: va_list;
static import core.simd;
static import std.conv;

struct Int128 { long lower; long upper; }
struct UInt128 { ulong lower; ulong upper; }

struct __locale_data { int dummy; }

alias _Bool = bool;
struct dpp {
    static struct Opaque(int N) {
        void[N] bytes;
    }

    static bool isEmpty(T)() {
        return T.tupleof.length == 0;
    }
    static struct Move(T) {
        T* ptr;
    }

    static auto move(T)(ref T value) {
        return Move!T(&value);
    }
    mixin template EnumD(string name, T, string prefix) if(is(T == enum)) {
        private static string _memberMixinStr(string member) {
            import std.conv: text;
            import std.array: replace;
            return text(` `, member.replace(prefix, ""), ` = `, T.stringof, `.`, member, `,`);
        }
        private static string _enumMixinStr() {
            import std.array: join;
            string[] ret;
            ret ~= "enum " ~ name ~ "{";
            static foreach(member; __traits(allMembers, T)) {
                ret ~= _memberMixinStr(member);
            }
            ret ~= "}";
            return ret.join("\n");
        }
        mixin(_enumMixinStr());
    }
}

extern(C)
{
    void cblas_zgeadd(const(CBLAS_ORDER), const(int), const(int), const(double)*, double*, const(int), const(double)*, double*, const(int)) @nogc nothrow;
    void cblas_cgeadd(const(CBLAS_ORDER), const(int), const(int), const(float)*, float*, const(int), const(float)*, float*, const(int)) @nogc nothrow;
    void cblas_dgeadd(const(CBLAS_ORDER), const(int), const(int), const(double), double*, const(int), const(double), double*, const(int)) @nogc nothrow;
    void cblas_sgeadd(const(CBLAS_ORDER), const(int), const(int), const(float), float*, const(int), const(float), float*, const(int)) @nogc nothrow;
    void cblas_zimatcopy(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(double)*, double*, const(int), const(int)) @nogc nothrow;
    void cblas_cimatcopy(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(float)*, float*, const(int), const(int)) @nogc nothrow;
    void cblas_dimatcopy(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(double), double*, const(int), const(int)) @nogc nothrow;
    void cblas_simatcopy(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(float), float*, const(int), const(int)) @nogc nothrow;
    void cblas_zomatcopy(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(double)*, const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_comatcopy(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(float)*, const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_domatcopy(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(double), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_somatcopy(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(float), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_zaxpby(const(int), const(void)*, const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_caxpby(const(int), const(void)*, const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_daxpby(const(int), const(double), const(double)*, const(int), const(double), double*, const(int)) @nogc nothrow;
    void cblas_saxpby(const(int), const(float), const(float)*, const(int), const(float), float*, const(int)) @nogc nothrow;
    void cblas_xerbla(int, char*, char*, ...) @nogc nothrow;
    void cblas_zher2k(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(double), void*, const(int)) @nogc nothrow;
    void cblas_cher2k(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(float), void*, const(int)) @nogc nothrow;
    void cblas_zherk(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(int), const(int), const(double), const(void)*, const(int), const(double), void*, const(int)) @nogc nothrow;
    void cblas_cherk(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(int), const(int), const(float), const(void)*, const(int), const(float), void*, const(int)) @nogc nothrow;
    void cblas_zhemm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_chemm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_ztrsm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(void)*, const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_ctrsm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(void)*, const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_dtrsm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(double), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_strsm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(float), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_ztrmm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(void)*, const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_ctrmm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(void)*, const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_dtrmm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(double), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_strmm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(float), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_zsyr2k(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_csyr2k(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_dsyr2k(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(int), const(int), const(double), const(double)*, const(int), const(double)*, const(int), const(double), double*, const(int)) @nogc nothrow;
    void cblas_ssyr2k(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(int), const(int), const(float), const(float)*, const(int), const(float)*, const(int), const(float), float*, const(int)) @nogc nothrow;
    void cblas_zsyrk(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_dsyrk(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(int), const(int), const(double), const(double)*, const(int), const(double), double*, const(int)) @nogc nothrow;
    void cblas_ssyrk(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(int), const(int), const(float), const(float)*, const(int), const(float), float*, const(int)) @nogc nothrow;
    void cblas_zsymm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_csymm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_dsymm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(int), const(int), const(double), const(double)*, const(int), const(double)*, const(int), const(double), double*, const(int)) @nogc nothrow;
    void openblas_set_num_threads(int) @nogc nothrow;
    //void goto_set_num_threads(int) @nogc nothrow;
    //int openblas_get_num_threads() @nogc nothrow;
    //int openblas_get_num_procs() @nogc nothrow;
    //char* openblas_get_config() @nogc nothrow;
    //char* openblas_get_corename() @nogc nothrow;
    //int openblas_get_parallel() @nogc nothrow;
    enum CBLAS_ORDER
    {
        CblasRowMajor = 101,
        CblasColMajor = 102,
    }
    enum CblasRowMajor = CBLAS_ORDER.CblasRowMajor;
    enum CblasColMajor = CBLAS_ORDER.CblasColMajor;
    enum CBLAS_TRANSPOSE
    {
        CblasNoTrans = 111,
        CblasTrans = 112,
        CblasConjTrans = 113,
        CblasConjNoTrans = 114,
    }
    enum CblasNoTrans = CBLAS_TRANSPOSE.CblasNoTrans;
    enum CblasTrans = CBLAS_TRANSPOSE.CblasTrans;
    enum CblasConjTrans = CBLAS_TRANSPOSE.CblasConjTrans;
    enum CblasConjNoTrans = CBLAS_TRANSPOSE.CblasConjNoTrans;
    enum CBLAS_UPLO
    {
        CblasUpper = 121,
        CblasLower = 122,
    }
    enum CblasUpper = CBLAS_UPLO.CblasUpper;
    enum CblasLower = CBLAS_UPLO.CblasLower;
    enum CBLAS_DIAG
    {
        CblasNonUnit = 131,
        CblasUnit = 132,
    }
    enum CblasNonUnit = CBLAS_DIAG.CblasNonUnit;
    enum CblasUnit = CBLAS_DIAG.CblasUnit;
    enum CBLAS_SIDE
    {
        CblasLeft = 141,
        CblasRight = 142,
    }
    enum CblasLeft = CBLAS_SIDE.CblasLeft;
    enum CblasRight = CBLAS_SIDE.CblasRight;
    alias CBLAS_LAYOUT = CBLAS_ORDER;
    float cblas_sdsdot(const(int), const(float), const(float)*, const(int), const(float)*, const(int)) @nogc nothrow;
    double cblas_dsdot(const(int), const(float)*, const(int), const(float)*, const(int)) @nogc nothrow;
    float cblas_sdot(const(int), const(float)*, const(int), const(float)*, const(int)) @nogc nothrow;
    double cblas_ddot(const(int), const(double)*, const(int), const(double)*, const(int)) @nogc nothrow;
    cfloat cblas_cdotu(const(int), const(void)*, const(int), const(void)*, const(int)) @nogc nothrow;
    cfloat cblas_cdotc(const(int), const(void)*, const(int), const(void)*, const(int)) @nogc nothrow;
    cdouble cblas_zdotu(const(int), const(void)*, const(int), const(void)*, const(int)) @nogc nothrow;
    cdouble cblas_zdotc(const(int), const(void)*, const(int), const(void)*, const(int)) @nogc nothrow;
    void cblas_cdotu_sub(const(int), const(void)*, const(int), const(void)*, const(int), void*) @nogc nothrow;
    void cblas_cdotc_sub(const(int), const(void)*, const(int), const(void)*, const(int), void*) @nogc nothrow;
    void cblas_zdotu_sub(const(int), const(void)*, const(int), const(void)*, const(int), void*) @nogc nothrow;
    void cblas_zdotc_sub(const(int), const(void)*, const(int), const(void)*, const(int), void*) @nogc nothrow;
    float cblas_sasum(const(int), const(float)*, const(int)) @nogc nothrow;
    double cblas_dasum(const(int), const(double)*, const(int)) @nogc nothrow;
    float cblas_scasum(const(int), const(void)*, const(int)) @nogc nothrow;
    double cblas_dzasum(const(int), const(void)*, const(int)) @nogc nothrow;
    float cblas_ssum(const(int), const(float)*, const(int)) @nogc nothrow;
    double cblas_dsum(const(int), const(double)*, const(int)) @nogc nothrow;
    float cblas_scsum(const(int), const(void)*, const(int)) @nogc nothrow;
    double cblas_dzsum(const(int), const(void)*, const(int)) @nogc nothrow;
    float cblas_snrm2(const(int), const(float)*, const(int)) @nogc nothrow;
    double cblas_dnrm2(const(int), const(double)*, const(int)) @nogc nothrow;
    float cblas_scnrm2(const(int), const(void)*, const(int)) @nogc nothrow;
    double cblas_dznrm2(const(int), const(void)*, const(int)) @nogc nothrow;
    c_ulong cblas_isamax(const(int), const(float)*, const(int)) @nogc nothrow;
    c_ulong cblas_idamax(const(int), const(double)*, const(int)) @nogc nothrow;
    c_ulong cblas_icamax(const(int), const(void)*, const(int)) @nogc nothrow;
    c_ulong cblas_izamax(const(int), const(void)*, const(int)) @nogc nothrow;
    c_ulong cblas_isamin(const(int), const(float)*, const(int)) @nogc nothrow;
    c_ulong cblas_idamin(const(int), const(double)*, const(int)) @nogc nothrow;
    c_ulong cblas_icamin(const(int), const(void)*, const(int)) @nogc nothrow;
    c_ulong cblas_izamin(const(int), const(void)*, const(int)) @nogc nothrow;
    c_ulong cblas_ismax(const(int), const(float)*, const(int)) @nogc nothrow;
    c_ulong cblas_idmax(const(int), const(double)*, const(int)) @nogc nothrow;
    c_ulong cblas_icmax(const(int), const(void)*, const(int)) @nogc nothrow;
    c_ulong cblas_izmax(const(int), const(void)*, const(int)) @nogc nothrow;
    c_ulong cblas_ismin(const(int), const(float)*, const(int)) @nogc nothrow;
    c_ulong cblas_idmin(const(int), const(double)*, const(int)) @nogc nothrow;
    c_ulong cblas_icmin(const(int), const(void)*, const(int)) @nogc nothrow;
    c_ulong cblas_izmin(const(int), const(void)*, const(int)) @nogc nothrow;
    void cblas_saxpy(const(int), const(float), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_daxpy(const(int), const(double), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_caxpy(const(int), const(void)*, const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_zaxpy(const(int), const(void)*, const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_scopy(const(int), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_dcopy(const(int), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_ccopy(const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_zcopy(const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_sswap(const(int), float*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_dswap(const(int), double*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_cswap(const(int), void*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_zswap(const(int), void*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_srot(const(int), float*, const(int), float*, const(int), const(float), const(float)) @nogc nothrow;
    void cblas_drot(const(int), double*, const(int), double*, const(int), const(double), const(double)) @nogc nothrow;
    void cblas_srotg(float*, float*, float*, float*) @nogc nothrow;
    void cblas_drotg(double*, double*, double*, double*) @nogc nothrow;
    void cblas_srotm(const(int), float*, const(int), float*, const(int), const(float)*) @nogc nothrow;
    void cblas_drotm(const(int), double*, const(int), double*, const(int), const(double)*) @nogc nothrow;
    void cblas_srotmg(float*, float*, float*, const(float), float*) @nogc nothrow;
    void cblas_drotmg(double*, double*, double*, const(double), double*) @nogc nothrow;
    void cblas_sscal(const(int), const(float), float*, const(int)) @nogc nothrow;
    void cblas_dscal(const(int), const(double), double*, const(int)) @nogc nothrow;
    void cblas_cscal(const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_zscal(const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_csscal(const(int), const(float), void*, const(int)) @nogc nothrow;
    void cblas_zdscal(const(int), const(double), void*, const(int)) @nogc nothrow;
    void cblas_sgemv(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(float), const(float)*, const(int), const(float)*, const(int), const(float), float*, const(int)) @nogc nothrow;
    void cblas_dgemv(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(double), const(double)*, const(int), const(double)*, const(int), const(double), double*, const(int)) @nogc nothrow;
    void cblas_cgemv(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_zgemv(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_sger(const(CBLAS_ORDER), const(int), const(int), const(float), const(float)*, const(int), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_dger(const(CBLAS_ORDER), const(int), const(int), const(double), const(double)*, const(int), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_cgeru(const(CBLAS_ORDER), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_cgerc(const(CBLAS_ORDER), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_zgeru(const(CBLAS_ORDER), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_zgerc(const(CBLAS_ORDER), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_strsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_dtrsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_ctrsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_ztrsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_strmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_dtrmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_ctrmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_ztrmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_ssyr(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(float), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_dsyr(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(double), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_cher(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(float), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_zher(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(double), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_ssyr2(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(float), const(float)*, const(int), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_dsyr2(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(double), const(double)*, const(int), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_cher2(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_zher2(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_sgbmv(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(int), const(int), const(float), const(float)*, const(int), const(float)*, const(int), const(float), float*, const(int)) @nogc nothrow;
    void cblas_dgbmv(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(int), const(int), const(double), const(double)*, const(int), const(double)*, const(int), const(double), double*, const(int)) @nogc nothrow;
    void cblas_cgbmv(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_zgbmv(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(int), const(int), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_ssbmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(int), const(float), const(float)*, const(int), const(float)*, const(int), const(float), float*, const(int)) @nogc nothrow;
    void cblas_dsbmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(int), const(double), const(double)*, const(int), const(double)*, const(int), const(double), double*, const(int)) @nogc nothrow;
    void cblas_stbmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_dtbmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_ctbmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_ztbmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_stbsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(float)*, const(int), float*, const(int)) @nogc nothrow;
    void cblas_dtbsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(double)*, const(int), double*, const(int)) @nogc nothrow;
    void cblas_ctbsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_ztbsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(int), const(void)*, const(int), void*, const(int)) @nogc nothrow;
    void cblas_stpmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(float)*, float*, const(int)) @nogc nothrow;
    void cblas_dtpmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(double)*, double*, const(int)) @nogc nothrow;
    void cblas_ctpmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_ztpmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_stpsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(float)*, float*, const(int)) @nogc nothrow;
    void cblas_dtpsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(double)*, double*, const(int)) @nogc nothrow;
    void cblas_ctpsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_ztpsv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(CBLAS_TRANSPOSE), const(CBLAS_DIAG), const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_ssymv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(float), const(float)*, const(int), const(float)*, const(int), const(float), float*, const(int)) @nogc nothrow;
    void cblas_dsymv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(double), const(double)*, const(int), const(double)*, const(int), const(double), double*, const(int)) @nogc nothrow;
    void cblas_chemv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_zhemv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_sspmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(float), const(float)*, const(float)*, const(int), const(float), float*, const(int)) @nogc nothrow;
    void cblas_dspmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(double), const(double)*, const(double)*, const(int), const(double), double*, const(int)) @nogc nothrow;
    void cblas_sspr(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(float), const(float)*, const(int), float*) @nogc nothrow;
    void cblas_dspr(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(double), const(double)*, const(int), double*) @nogc nothrow;
    void cblas_chpr(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(float), const(void)*, const(int), void*) @nogc nothrow;
    void cblas_zhpr(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(double), const(void)*, const(int), void*) @nogc nothrow;
    void cblas_sspr2(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(float), const(float)*, const(int), const(float)*, const(int), float*) @nogc nothrow;
    void cblas_dspr2(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(double), const(double)*, const(int), const(double)*, const(int), double*) @nogc nothrow;
    void cblas_chpr2(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), void*) @nogc nothrow;
    void cblas_zhpr2(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), void*) @nogc nothrow;
    void cblas_chbmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_zhbmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_chpmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(void)*, const(void)*, const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_zhpmv(const(CBLAS_ORDER), const(CBLAS_UPLO), const(int), const(void)*, const(void)*, const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_sgemm(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(CBLAS_TRANSPOSE), const(int), const(int), const(int), const(float), const(float)*, const(int), const(float)*, const(int), const(float), float*, const(int)) @nogc nothrow;
    void cblas_dgemm(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(CBLAS_TRANSPOSE), const(int), const(int), const(int), const(double), const(double)*, const(int), const(double)*, const(int), const(double), double*, const(int)) @nogc nothrow;
    void cblas_cgemm(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(CBLAS_TRANSPOSE), const(int), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_cgemm3m(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(CBLAS_TRANSPOSE), const(int), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_zgemm(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(CBLAS_TRANSPOSE), const(int), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_zgemm3m(const(CBLAS_ORDER), const(CBLAS_TRANSPOSE), const(CBLAS_TRANSPOSE), const(int), const(int), const(int), const(void)*, const(void)*, const(int), const(void)*, const(int), const(void)*, void*, const(int)) @nogc nothrow;
    void cblas_ssymm(const(CBLAS_ORDER), const(CBLAS_SIDE), const(CBLAS_UPLO), const(int), const(int), const(float), const(float)*, const(int), const(float)*, const(int), const(float), float*, const(int)) @nogc nothrow;
}


struct __va_list_tag;
