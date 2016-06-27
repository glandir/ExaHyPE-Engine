#!/bin/env python

import Backend



def generateImmediateDSCAL(i_scalar: float, i_inVectorName: str, i_outVectorName: str, i_vectorSize: int) -> str:
    l_code = '#pragma simd\n'                                \
             '  for(int i=0;i<'+str(i_vectorSize)+';i++)\n'  \
             '    '+i_outVectorName+'[i] = '+ str(i_scalar) +' * '+i_inVectorName+'[i];\n'
    return l_code


def generateDSCAL(i_scalarName: str, i_inVectorName: str, i_outVectorName:str, i_vectorSize: int) -> str:
    l_code = '#pragma simd\n'                                \
             '  for(int i=0;i<'+str(i_vectorSize)+';i++) \n' \
             '    '+i_outVectorName+'[i] = '+ i_scalarName+' * '+i_inVectorName+'[i];\n'
    return l_code


# --------------------------------------------------------------------------------------

def __scatter_avx(i_nVar: int, i_chunkSize: int, i_baseAddr_lhs: int, i_baseAddr_rhs: int, i_simdWidth: int) -> str:
    """
    Scatter for architectures with a SIMD width of 4 working
    on four input vectors

    Args:
        i_nVar:
            total number of variables
        i_chunkSize:
            padded number of dofs that separates two distinct variables
        i_baseAddr_lhs:
            start address of the output buffer
        i_baseAddr_rhs:
            start address of the input buffer
        i_simdWidth:
            the SIMD width used for computing offsets (may be a fraction of the true SIMD width)
    Returns:
        l_code:
            intrinsics implementation of the scatter operation
    """
    l_iters = int(i_nVar/i_simdWidth)
    l_remainder = i_nVar % i_simdWidth
    l_offset = Backend.getSizeWithPadding(i_nVar)
    # E.g. 9 Vars => 2*simdWidth + 1 => 2*packed + 1*remainder


    # fully packed
    l_startAddress = 0
    l_code = ''
    for it in range(0, l_iters):
        # aligned load
        l_code = l_code + 'v1 = _mm256_load_pd(&in_buf['+str(i_baseAddr_rhs+(0*l_offset+i_simdWidth*it))+']);\n'
        l_code = l_code + 'v2 = _mm256_load_pd(&in_buf['+str(i_baseAddr_rhs+(1*l_offset+i_simdWidth*it))+']);\n'
        l_code = l_code + 'v3 = _mm256_load_pd(&in_buf['+str(i_baseAddr_rhs+(2*l_offset+i_simdWidth*it))+']);\n'
        l_code = l_code + 'v4 = _mm256_load_pd(&in_buf['+str(i_baseAddr_rhs+(3*l_offset+i_simdWidth*it))+']);\n'
        # vunpckhpd, vunpcklpd
        l_code = l_code + 'perm1 = _mm256_unpacklo_pd(v1,v2);\n'
        l_code = l_code + 'perm2 = _mm256_unpackhi_pd(v1,v2);\n'
        l_code = l_code + 'perm3 = _mm256_unpacklo_pd(v3,v4);\n'
        l_code = l_code + 'perm4 = _mm256_unpackhi_pd(v3,v4);\n'
        # vperm2f128
        l_code = l_code + 'res1 = _mm256_permute2f128_pd(perm1,perm3,0b00100000);\n'
        l_code = l_code + 'res3 = _mm256_permute2f128_pd(perm1,perm3,0b00110001);\n'
        l_code = l_code + 'res2 = _mm256_permute2f128_pd(perm2,perm4,0b00100000);\n'
        l_code = l_code + 'res4 = _mm256_permute2f128_pd(perm2,perm4,0b00110001);\n'
        # unaligned store. Should later on become _mm256_store_pd()
        l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(i_baseAddr_lhs+l_startAddress)+'],res1);\n'
        l_startAddress = l_startAddress + i_chunkSize
        l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(i_baseAddr_lhs+l_startAddress)+'],res2);\n'
        l_startAddress = l_startAddress + i_chunkSize
        l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(i_baseAddr_lhs+l_startAddress)+'],res3);\n'
        l_startAddress = l_startAddress + i_chunkSize
        l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(i_baseAddr_lhs+l_startAddress)+'],res4);\n'
        l_startAddress = l_startAddress + i_chunkSize

    # we need some scalar instructions for the remaining variables
    if(l_remainder > 0):
        # aligned packed load - in_buf is padded
        l_startRemainder = i_baseAddr_rhs + (l_iters * i_simdWidth)
        l_code = l_code + 'v1 = _mm256_load_pd(&in_buf['+str(0*l_offset+l_startRemainder)+']);\n'
        l_code = l_code + 'v2 = _mm256_load_pd(&in_buf['+str(1*l_offset+l_startRemainder)+']);\n'
        l_code = l_code + 'v3 = _mm256_load_pd(&in_buf['+str(2*l_offset+l_startRemainder)+']);\n'
        l_code = l_code + 'v4 = _mm256_load_pd(&in_buf['+str(3*l_offset+l_startRemainder)+']);\n'
        # vunpckhpd, vunpcklpd
        l_code = l_code + 'perm1 = _mm256_unpacklo_pd(v1,v2);\n'
        l_code = l_code + 'perm2 = _mm256_unpackhi_pd(v1,v2);\n'
        l_code = l_code + 'perm3 = _mm256_unpacklo_pd(v3,v4);\n'
        l_code = l_code + 'perm4 = _mm256_unpackhi_pd(v3,v4);\n'
        # vperm2f128
        l_code = l_code + 'res1 = _mm256_permute2f128_pd(perm1,perm3,0b00100000);\n'
        l_code = l_code + 'res3 = _mm256_permute2f128_pd(perm1,perm3,0b00110001);\n'
        l_code = l_code + 'res2 = _mm256_permute2f128_pd(perm2,perm4,0b00100000);\n'
        l_code = l_code + 'res4 = _mm256_permute2f128_pd(perm2,perm4,0b00110001);\n'

        if(l_remainder > 0):
            l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(i_baseAddr_lhs+l_startAddress)+'],res1);\n'
            l_startAddress = l_startAddress + i_chunkSize
        if(l_remainder > 1):
            l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(i_baseAddr_lhs+l_startAddress)+'],res2);\n'
            l_startAddress = l_startAddress + i_chunkSize
        if(l_remainder > 2):
            l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(i_baseAddr_lhs+l_startAddress)+'],res3);\n'

    return l_code


def __scatter_sse2(i_nVar: int, i_chunkSize: int, i_baseAddr_lhs: int, i_baseAddr_rhs: int, i_simdWidth: int) -> str:
    """
    Scatter for architectures with a SIMD width of 2 working
    on two input vectors

    Args:
        i_nVar:
            total number of variables
        i_chunkSize:
            padded number of dofs that separates two distinct variables
        i_baseAddr_lhs:
            start address of the output buffer
        i_baseAddr_rhs:
            start address of the input buffer
        i_simdWidth:
            the SIMD width used for computing offsets (may be a fraction of the true SIMD width)
    Returns:
        l_code:
            intrinsics implementation of the scatter operation
    """
    l_iters = int(i_nVar/i_simdWidth)
    l_remainder = i_nVar % i_simdWidth
    l_offset = Backend.getSizeWithPadding(i_nVar)
    # E.g. 9 Vars => 4*simdWidth + 1 => 4*packed + 1*remainder

    # fully packed
    l_startAddress = 0
    l_code = ''
    for it in range(0, l_iters):
        # aligned load
        l_code = l_code + 'v1 = _mm_load_pd(&in_buf['+str(i_baseAddr_rhs+(0*l_offset+i_simdWidth*it))+']);\n'
        l_code = l_code + 'v2 = _mm_load_pd(&in_buf['+str(i_baseAddr_rhs+(1*l_offset+i_simdWidth*it))+']);\n'
        # unpcklpd
        l_code = l_code + 'res = _mm_unpacklo_pd(v1, v2);\n'
        # unaligned store. Should later on become _mm_store_pd()
        l_code = l_code + '_mm_storeu_pd(&out_buf['+str(i_baseAddr_lhs+l_startAddress)+'],res);\n'
        l_startAddress = l_startAddress + i_chunkSize
        # unpackhpd
        l_code = l_code + 'res = _mm_unpackhi_pd(v1,v2);\n'
        l_code = l_code + '_mm_storeu_pd(&out_buf['+str(i_baseAddr_lhs+l_startAddress)+'],res);\n'
        l_startAddress = l_startAddress + i_chunkSize

    # there is one variable remaining
    if(l_remainder > 0):
        l_startRemainder = i_baseAddr_rhs + (l_iters * i_simdWidth)
        l_code = l_code + 'out_buf['+str(i_baseAddr_lhs+l_startAddress)+'] = in_buf['+str(l_startRemainder)+'];\n'
        l_code = l_code + 'out_buf['+str(i_baseAddr_lhs+l_startAddress+1)+'] = in_buf['+str(l_startRemainder+l_offset)+'];\n'

    return l_code


def __scatter_scalar(i_nVar: int, i_chunkSize: int, i_startAddr_lhs: int, i_startAddr_rhs: int) -> str:
    """
    Scalar scatter operation working on one input vector

    Args:
        i_nVar:
            total number of variables
        i_chunkSize:
            padded number of dofs that separates two variables
        i_startAddr_lhs:
            start address of the output buffer
        i_startAddr_rhs:
            start address of the input buffer

    Returns:
        l_code:
            plain C implementation of the scatter operation
    """
    l_startAddress = 0
    l_code = ''
    for iVar in range(0, i_nVar):
        l_inAddr  = i_startAddr_rhs+iVar
        l_outAddr = i_startAddr_lhs+iVar*i_chunkSize
        l_code = l_code + 'out_buf['+str(l_outAddr)+'] = in_buf['+str(l_inAddr)+'];\n'

    return l_code


def generateScatter(i_nVar: int, i_nVectors: int, i_chunkSize: int):
    """
    Interface for the generation of a scatter routine. Internally
    decomposes the scatter operation into micro kernels for which
    the corresponding intrinsics are generated. Creates a file
    asm_scatter.c with an intrinsics implementation.

    Args:
        i_nVar:
            total number of variables (without padding)
        i_nVectors:
            total number of vectors (without padding)
        i_chunkSize:
            padded number of dofs that separates two sets of variables
    """

    # decompose the number of input vectors
    i_architecture = Backend.m_architecture
    i_simdWidth = Backend.m_simdWidth['DP'][i_architecture]
    l_nVectorGroup = int(i_nVectors/i_simdWidth)
    l_nVectorRest  = i_nVectors % i_simdWidth

    print("i_nVectors", i_nVectors)
    print("l_nVectorGroup", l_nVectorGroup)
    print("l_nVectorRest", l_nVectorRest)
    print("i_chunkSize", i_chunkSize)

    l_sourceFile = open("asm_scatter.c", 'a')
    l_sourceFile.write('#if defined( __SSE3__) || defined(__MIC__)\n'\
                       '#include <immintrin.h>\n'\
                       '#endif\n\n'\
                       'void scatter(const double* restrict const in_buf, double* out_buf) {\n'
                       )

    l_code = ''
    l_startAddr_rhs = 0
    l_startAddr_lhs = 0

    if(i_architecture == 'hsw' or i_architecture == 'snb'):
        # declaration of variables for macro kernel
        l_code = '__m256d v1, v2, v3, v4, perm1, perm2, perm3, perm4, res1, res2, res3, res4;\n'

        # full vector groups (micro kernel 1)
        for iVec in range(0, l_nVectorGroup):
            l_code = l_code + __scatter_avx(i_nVar, i_chunkSize, l_startAddr_lhs, l_startAddr_rhs, i_simdWidth)
            l_startAddr_rhs = l_startAddr_rhs + 4*Backend.getSizeWithPadding(i_nVar)
            l_startAddr_lhs = l_startAddr_lhs + 4

        # rest (micro kernel 2)
        if(l_nVectorRest == 1):
            l_code = l_code + __scatter_scalar(i_nVar, i_chunkSize, l_startAddr_lhs, l_startAddr_rhs)
        if(l_nVectorRest == 2):
            l_code = l_code + __scatter_sse2(i_nVar, i_chunkSize, l_startAddr_lhs, l_startAddr_rhs, 2)
            l_startAddr_rhs = l_startAddr_rhs + 2*Backend.getSizeWithPadding(i_nVar)
            l_startAddr_lhs = l_startAddr_lhs + 2
        if(l_nVectorRest == 3):
            l_code = l_code + __scatter_sse2(i_nVar, i_chunkSize, l_startAddr_lhs, l_startAddr_rhs, 2)
            l_startAddr_rhs = l_startAddr_rhs + 2*Backend.getSizeWithPadding(i_nVar)
            l_startAddr_lhs = l_startAddr_lhs + 2
            l_code = l_code + __scatter_scalar(i_nVar, i_chunkSize, l_startAddr_lhs, l_startAddr_rhs)

        l_sourceFile.write(l_code)

    if(i_architecture == 'wsm'):
        # declaration of variables for macro kernel
        l_code = '__m128d v1,v2,res;\n'

        # process full vector groups (micro kernel 1)
        for iVec in range(0, l_nVectorGroup):
            l_code = l_code + __scatter_sse2(i_nVar, i_chunkSize, l_startAddr_lhs, l_startAddr_rhs, i_simdWidth)
            l_startAddr_rhs = l_startAddr_rhs + 2*Backend.getSizeWithPadding(i_nVar)
            l_startAddr_lhs = l_startAddr_lhs + 2

        # rest (micro kernel 2)
        for iVec in range(0, l_nVectorRest):
            l_code = l_code + __scatter_scalar(i_nVar, i_chunkSize, l_startAddr_lhs, l_startAddr_rhs)


        l_sourceFile.write(l_code)


    if(i_architecture == 'noarch'):
        for iVec in range(0, i_nVectors):
            l_code = l_code + __scatter_scalar(i_nVar, i_chunkSize, l_startAddr_lhs, l_startAddr_rhs)
            l_startAddr_rhs = l_startAddr_rhs + Backend.getSizeWithPadding(i_nVar)
            l_startAddr_lhs = l_startAddr_lhs + 1

        l_sourceFile.write(l_code)


    l_sourceFile.write('}')
    l_sourceFile.close()















