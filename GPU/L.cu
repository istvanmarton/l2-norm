/*
compile:
  nvcc -gencode arch=compute_80,code=sm_80 -o L L.cu
use:
  ./L s7.mat
*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define luint unsigned long long int
#define elemt int
#define maxHeight 96
#define maxWidth 96
#define maxWidthB 100
#define offs 32
#define procNum 1536

__device__ elemt best[1];
__device__ luint counter[1];
__device__ __constant__ elemt mat[maxHeight * maxWidthB];

__global__ void func(elemt *matt, const int grX, const int endX, const luint blocks) {
  elemt v1, v2, v3, tmp, m;
  luint b = blockIdx.x;
  const int l = threadIdx.x;
  elemt* mat = matt + l;
  elemt* matB = mat + offs;
  elemt* matC = matB + offs;
  elemt* cvec = matt + maxWidth;
  int gr   = grX  * maxWidthB;
  int end  = endX * maxWidthB;
  int mtx;

beg:
  mtx = 0;
  m = best[0];
  if (mtx == end) goto end;

  v1 = mat[mtx];
  v2 = matB[mtx];
  v3 = matC[mtx];
  while (mtx < gr - maxWidthB) {
    mtx += maxWidthB;
    tmp = mat[mtx];
    if (b & 1) tmp = -tmp;
    v1 += tmp;
    tmp = matB[mtx];
    if (b & 1) tmp = -tmp;
    v2 += tmp;
    tmp = matC[mtx];
    if (b & 1) tmp = -tmp;
    v3 += tmp;
    b >>= 1;
  }

d:
  mtx += maxWidthB;
  tmp = __reduce_add_sync(0xffffffff, abs(v1) + abs(v2) + abs(v3));

  if (mtx == end) {
    m = max(m, tmp);
    goto u;
  }
  if (tmp + cvec[mtx] > m) {
    b <<= 1;
    v1 += mat[mtx];
    v2 += matB[mtx];
    v3 += matC[mtx];
    goto d;
  }

u:
  if (mtx == gr) goto end;
  mtx -= maxWidthB;
  if (b & 1) {
    b >>= 1;
    v1 += mat[mtx];
    v2 += matB[mtx];
    v3 += matC[mtx];
    goto u;
  }
  b++;
  v1 -= mat[mtx] << 1;
  v2 -= matB[mtx] << 1;
  v3 -= matC[mtx] << 1;
  goto d;

end:
  if (l == 0) {
    tmp = best[0];
    if (m > tmp) atomicMax(best, m);
    b = atomicAdd(counter, (luint)1);
  }
  b = __shfl_sync(0xffffffff, b, 0);
  if ((long long int)b < blocks) goto beg;
}

elemt comp(int ofs, int granM, int height, elemt guess) {
  elemt *mtx;
  cudaGetSymbolAddress((void**)&mtx, mat);
  int gran = min(max(11, granM) + 1, height) - 1;
  luint blocks = (luint)1 << gran;
//  luint blk = (luint)1 << min(gran, 11);
  luint blk = gran < 32 ? min((int)blocks, procNum) : procNum;
  cudaMemcpyToSymbol(best, &guess, sizeof(elemt));
  cudaMemcpyToSymbol(counter, &blk, sizeof(luint));
  func<<<blk, 32>>>(mtx + ofs * maxWidthB, gran + 1, height, blocks);
  cudaMemcpyFromSymbol(&guess, best, sizeof(elemt));
  cudaMemcpy(mtx + ofs * maxWidthB + maxWidth, &guess, sizeof(elemt), cudaMemcpyHostToDevice);
  return guess;
}

int main(int argc, char* argv[]){
  elemt mtx[maxHeight * maxWidthB];
  FILE* fp = fopen(argv[1], "r");
  char *line = NULL;
  size_t len = 0;
  int height = 0, width = 0;
  while (1) {
    ssize_t read = getline(&line, &len, fp);
    if (read <= 1) break;
    int offset = 0, inc, j = 0;
    while (1) {
      int res = sscanf(line + offset, "%d%n", &mtx[height * maxWidthB + j], &inc);
      if (res < 1) break;
      offset += inc;
      j++;
    }
    for (int k = j; k < maxWidth; k++) mtx[height * maxWidthB + k] = 0;
    mtx[height * maxWidthB + maxWidth] = 1 << 30;
    height++;
    width = max(width, j);
  }
  free(line);
  fclose(fp);
  printf("%d x %d matrix\n", height, width);

  cudaMemcpyToSymbol(mat, mtx, maxWidthB * height * sizeof(elemt));

  elemt l1;
  int lim = height / 2 + 2;
  for (int i = height; i >= 0; i--) {
    if (i > lim || i == 0) {
      l1 = comp(i, (height - i) / 2 + 2, height - i, l1);
//      printf("L1 norm at %d: %d\n", i, l1);
      }
  }
  printf("L norm: %d\n", l1);

  return 0;
}
