/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

#define b_ 5
#define s_ 5
#define B_ 32
#define S_ 32

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */

#define blk32 8
#define blk64 4
#define blk61 16
#define blk67 16
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    if(M == 32 && N == 32){
        int i, j, m, n;
        for(m = 0; m < M; m += blk32){
            for(n = 0; n < N; n += blk32){
                for(i = 0; i < blk32; ++i){
                    for(j = 0; j <= i; ++j){
                        B[j+m][i+n] = A[i+n][j+m];
                    }
                }
                for(j = 0; j < blk32; ++j){
                    for(i = 0; i < j; ++i){
                        B[j+m][i+n] = A[i+n][j+m];
                    }
                }
            }
        }        
    } else if (M == 61 && N == 67) {
        int i, j, m, n;
        for(m = 0; m < M; m += blk61){
            for(n = 0; n < N; n += blk67){
                for(i = n; i < blk67+n && i < N; ++i){
                    for(j = m; j < blk61+m && j <= i && j < M; ++j){
                        B[j][i] = A[i][j];
                    }
                }
                for(j = m; j < blk61+m && j < M; ++j){
                    for(i = n; i < blk67+n && i < j && i < N; ++i){
                        B[j][i] = A[i][j];
                    }
                }
            }
        }
    } else if (M == 64 && N == 64) {
        int i, j, m, n;
        int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
        for(n = 0; n < N; n += blk64){//处理对角线处
            for(m = n / (2*blk64) * (2*blk64); m < n / (2*blk64) * (2*blk64) + 2*blk64; m += blk64){
                // for(i = 0; i < blk64; ++i){
                //     for(j = 0; j <= i; ++j){
                //         B[j+m][i+n] = A[i+n][j+m];
                //     }
                // }
                tmp1 = A[n][m+1];
                tmp2 = A[n][m+2];
                tmp3 = A[n][m+3];
                B[m][n] = A[n][m];

                tmp4 = A[n+1][m+2];
                tmp5 = A[n+1][m+3];
                B[m][n+1] = A[n+1][m];
                B[m+1][n+1] = A[n+1][m+1];

                tmp6 = A[n+2][m+3];
                B[m][n+2] = A[n+2][m];
                B[m+1][n+2] = A[n+2][m+1];
                B[m+2][n+2] = A[n+2][m+2];
                
                B[m][n+3] = A[n+3][m];
                B[m+1][n+3] = A[n+3][m+1];
                B[m+2][n+3] = A[n+3][m+2];
                B[m+3][n+3] = A[n+3][m+3];

                // for(j = 0; j < blk64; ++j){
                //     for(i = 0; i < j; ++i){
                //         B[j+m][i+n] = A[i+n][j+m];
                //     }
                // }
                B[m+1][n] = tmp1;
                B[m+2][n] = tmp2;
                B[m+3][n] = tmp3;
                B[m+2][n+1] = tmp4;
                B[m+3][n+1] = tmp5;
                B[m+3][n+2] = tmp6;
            }
        }
        for(n = 0; n < N; n += blk64 * 2){
            for(m = 0; m < M; m += blk64 * 2){
                if(m == n) continue;
                for(i = n; i < n + blk64; ++i){
                    for(j = m; j < m + blk64; ++j){
                        B[j][i] = A[i][j];
                    }
                }
                tmp1 = A[n][m+4];
                tmp2 = A[n][m+5];
                tmp3 = A[n][m+6];
                tmp4 = A[n][m+7];
                tmp5 = A[n+1][m+4];
                tmp6 = A[n+1][m+5];
                tmp7 = A[n+1][m+6];
                tmp8 = A[n+1][m+7];            
                for(i = n + blk64; i < n + blk64 * 2; ++i){
                    for(j = m; j < m + blk64; ++j){
                        B[j][i] = A[i][j];
                    }
                }
                for(i = n + blk64; i < n + blk64 * 2; ++i){
                    for(j = m + blk64; j < m + blk64 * 2; ++j){
                        B[j][i] = A[i][j];
                    }
                }
                
                B[m+4][n] = tmp1;
                B[m+5][n] = tmp2;
                B[m+6][n] = tmp3;
                B[m+7][n] = tmp4;
                B[m+4][n+1] = tmp5;
                B[m+5][n+1] = tmp6;
                B[m+6][n+1] = tmp7;
                B[m+7][n+1] = tmp8;
                for(i = n + blk64 / 2; i < n + blk64; ++i){
                    for(j = m + blk64; j < m + blk64*2; ++j){
                        B[j][i] = A[i][j];
                    }
                }
            }
        }        
    }
}

/* 
* You can define additional transpose functions below. We've defined
* a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }    

}

char trans_test64x64_desc[] = "attempt 64x64 matrix";
void trans_test64x64(int M, int N, int A[N][M], int B[M][N])
{
    if(M != 64 || N != 64) return;
    int i, j, m, n;
    for(n = 0; n < N; n+=blk64){
        for(m = 0; m < M; m+=blk64){
            for(i = 0; i < blk64; ++i){
                for(j = 0; j <= i; ++j){
                    B[j+m][i+n] = A[i+n][j+m];
                }
            }
            for(j = 0; j < blk64; ++j){
                for(i = 0; i < j; ++i){
                    B[j+m][i+n] = A[i+n][j+m];
                }
            }
        }
    }
}

char trans1_test64x64_desc[] = "attempt 64x64 matrix 1";
void trans1_test64x64(int M, int N, int A[N][M], int B[M][N])
{
    if(M != 64 || N != 64) return;
    int i, j, m, n;
    for(n = 0; n < N; n += blk64){
        for(m = 0; m < M; m += blk64){
            for(i = n; i < n + blk64; ++i){
                for(j = m; j < m + blk64 && j <= i; ++j){
                    B[j][i] = A[i][j];
                }
            }
            for(j = m; j < m + blk64; ++j){
                for(i = n; i < n + blk64 && i < j; ++i){
                    B[j][i] = A[i][j];
                }
            }
        }
    }
}

char trans2_test64x64_desc[] = "attempt 64x64 matrix 2";
void trans2_test64x64(int M, int N, int A[N][M], int B[M][N])
{
    if(M != 64 || N != 64) return;
    int i, j, m, n;
    int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;
    for(n = 0; n < N; n += blk64){//处理对角线处
        for(m = n / (2*blk64) * (2*blk64); m < n / (2*blk64) * (2*blk64) + 2*blk64; m += blk64){
            // for(i = 0; i < blk64; ++i){
            //     for(j = 0; j <= i; ++j){
            //         B[j+m][i+n] = A[i+n][j+m];
            //     }
            // }
            tmp1 = A[n][m+1];
            tmp2 = A[n][m+2];
            tmp3 = A[n][m+3];
            B[m][n] = A[n][m];

            tmp4 = A[n+1][m+2];
            tmp5 = A[n+1][m+3];
            B[m][n+1] = A[n+1][m];
            B[m+1][n+1] = A[n+1][m+1];

            tmp6 = A[n+2][m+3];
            B[m][n+2] = A[n+2][m];
            B[m+1][n+2] = A[n+2][m+1];
            B[m+2][n+2] = A[n+2][m+2];

            B[m][n+3] = A[n+3][m];
            B[m+1][n+3] = A[n+3][m+1];
            B[m+2][n+3] = A[n+3][m+2];
            B[m+3][n+3] = A[n+3][m+3];

            // for(j = 0; j < blk64; ++j){
            //     for(i = 0; i < j; ++i){
            //         B[j+m][i+n] = A[i+n][j+m];
            //     }
            // }
            B[m+1][n] = tmp1;
            B[m+2][n] = tmp2;
            B[m+3][n] = tmp3;
            B[m+2][n+1] = tmp4;
            B[m+3][n+1] = tmp5;
            B[m+3][n+2] = tmp6;
        }
    }
    for(n = 0; n < N; n += blk64 * 2){
        for(m = 0; m < M; m += blk64 * 2){
            if(m == n) continue;
            for(i = n; i < n + blk64; ++i){
                for(j = m; j < m + blk64; ++j){
                    B[j][i] = A[i][j];
                }
            }
            tmp1 = A[n][m+4];
            tmp2 = A[n][m+5];
            tmp3 = A[n][m+6];
            tmp4 = A[n][m+7];
            tmp5 = A[n+1][m+4];
            tmp6 = A[n+1][m+5];
            tmp7 = A[n+1][m+6];
            tmp8 = A[n+1][m+7];            
            for(i = n + blk64; i < n + blk64 * 2; ++i){
                for(j = m; j < m + blk64; ++j){
                    B[j][i] = A[i][j];
                }
            }
            for(i = n + blk64; i < n + blk64 * 2; ++i){
                for(j = m + blk64; j < m + blk64 * 2; ++j){
                    B[j][i] = A[i][j];
                }
            }
            
            B[m+4][n] = tmp1;
            B[m+5][n] = tmp2;
            B[m+6][n] = tmp3;
            B[m+7][n] = tmp4;
            B[m+4][n+1] = tmp5;
            B[m+5][n+1] = tmp6;
            B[m+6][n+1] = tmp7;
            B[m+7][n+1] = tmp8;
            for(i = n + blk64 / 2; i < n + blk64; ++i){
                for(j = m + blk64; j < m + blk64*2; ++j){
                    B[j][i] = A[i][j];
                }
            }
        }
    }
}


/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc); 
    registerTransFunction(trans_test64x64, trans_test64x64_desc);
    registerTransFunction(trans1_test64x64, trans1_test64x64_desc);
    registerTransFunction(trans2_test64x64, trans2_test64x64_desc);
}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

