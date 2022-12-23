#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define length 4096
#define NUM_OF_THREADS 16384
#define NUM_OF_BITS 8 * sizeof(unsigned long long int)

__global__ void func(int* d_mtx_to_vec, unsigned long long int steps, unsigned long long int steps_remainder, int *d_L1_vector, int *d_L1_strategy, int iLonger, int iShorter){
	int i, l;

	int temp[length], vect[NUM_OF_BITS + 1], product = 0, L1 = 0;
	unsigned long long int number, index, iMax, iMin, logical;

	index = blockIdx.x * blockDim.x + threadIdx.x;
	iMax = (index + 1) *(steps) - 1;
	iMin = index * (steps);
	if(index < (steps_remainder) ) iMax += index + 1;
	else iMax += (steps_remainder);
	if(index <= (steps_remainder)) iMin += index;
	else iMin += (steps_remainder);
//	if(index == 1) printf("index: %llu, iMin: %llu, iMax: %llu\n", index, iMin, iMax);
		 
		 for(number = iMin; number <= iMax; number++){
		 product = 0;
		 for(l=0; l < iLonger; l++) {temp[l] = d_mtx_to_vec[(iShorter-1) * iLonger + l];}
			for(i = 0 ; (iShorter - 1) > i; i++){
				logical = number & ((unsigned long long int) 1 << i);
				logical = logical == 0 ? 0 : 1;
				vect[i] = (int) 2 * logical - 1;
				for(l=0; l < iLonger; l++){
					temp[l] += d_mtx_to_vec[i * iLonger + l] * vect[i];
				}
			}
			for(l= 0; l < iLonger; l++) {product += abs(temp[l]);}
			if(product > L1){
				L1 = product;
				for(l=0; l<(iShorter - 1); l++){d_L1_strategy[index * (iShorter - 1) + l] = vect[l];}
			}
		}

d_L1_vector[index] = L1;
}

int** mtx_read(int *iRows, int *iCols, char* fileName){
	int i = 0,j = 0, k = 0;
	int *row, **mtx, value;
	
	mtx = NULL;
	row = NULL;
	
	char g, cNum[256];
	
	FILE *fp;
	fp = fopen(fileName,"r");
	
	do{
		
		g = fgetc(fp);
		
//		
		if((((g - '0') < 10) && ((g - '0') >= 0)) || (g == 'e') || ( g == 'E') || (g == '.') || (g == '+') || (g == '-')) {cNum[i] = g; i++;}
		else {
			cNum[i] = '\0'; 
			if(cNum[0] != '\0') {sscanf(cNum, "%d", &value); j++; i = 0;  row = (int*) realloc(row, j * sizeof(int)); row[j-1] = value;}
			if( ((g == '\n') || (g == EOF)) && (j > 0)){ *iCols = j; j = 0; k++; mtx = (int**) realloc(mtx, k * sizeof(int*)); mtx[k-1] = row; row = NULL;}
		}

	}while(!feof(fp));
	*iRows = k;
printf("rows: %d, cols: %d\n",*iRows, *iCols); 
	fclose(fp);
return mtx;
}

void fileN(char *fileName, char** argv, int *argc){
	int r;
	FILE *fp;
	if((*argc) < 2) {
		do{
			printf("Please give me a filename: "); 
			r = scanf("%s",fileName);
		}while(r != 1);
	}
	else sprintf(fileName,"%s", argv[1]);

	fp = fopen(fileName, "r");
	if(fp == NULL) {
		do{
			printf("Please give me a filename that exist within this directory: ");
			r = scanf("%s",fileName);
			fp = fopen(fileName, "r");
		}while(fp == NULL);
	}
	fclose(fp);
}

int main(int argc, char *argv[]){
cudaDeviceProp devProp;
cudaGetDeviceProperties(&devProp, 0);
     char fileName[1024];
     fileN(fileName, argv, &argc);     
     int i, j, iMax, iRows, iCols, **mtx, *mtx_to_vec, *d_mtx_to_vec, iShorter, iLonger;

     mtx = mtx_read(&iRows, &iCols, fileName);
     mtx_to_vec = (int*)calloc(iRows * iCols, sizeof(int));
	if(iRows > iCols){
		for(j = 0; j < iCols; j++){
			for(i = 0; i < iRows; i++){
				mtx_to_vec[j * iRows + i] = mtx[i][j];
			}
		}
	}
	else{
		for(i = 0; i < iRows; i++){
			for(j = 0; j < iCols; j++){
				mtx_to_vec[i * iCols + j] = mtx[i][j];
			}
		}
	}

	if(iRows < iCols) {iShorter = iRows; iLonger = iCols;}
	else {iShorter = iCols; iLonger = iRows;}
	if(iShorter > (NUM_OF_BITS + 1)) {printf("Matrix is too big. The number of rows or columns can not be more than %lu.\n", NUM_OF_BITS + 1); return 0;}
	if(iLonger > length) {printf("Matrix is too big. The length variable %d should be bigger or equal than %d.\n", length, iLonger); return 0;}
	cudaMalloc((void**)&d_mtx_to_vec, iRows * iCols * sizeof(int));
     	cudaMalloc((void**)&d_mtx_to_vec, iRows * iCols * sizeof(int));
     	unsigned long long int steps, steps_remainder, Inner_num = (unsigned long long int) 1 << (iShorter - 1), copyNum;
	copyNum = NUM_OF_THREADS > Inner_num ? Inner_num : NUM_OF_THREADS;
	int *L1_vector, *d_L1_vector, L1_max = 0, *L1_strategy, *d_L1_strategy, num_ofBlock = (int) ceil((float)copyNum/(float)devProp.warpSize), num_ofThread = copyNum < devProp.warpSize ? copyNum : devProp.warpSize;

	steps=Inner_num/copyNum; steps_remainder = Inner_num % copyNum;

	L1_vector = (int*) malloc(copyNum * sizeof(int));
	L1_strategy = (int*) malloc(copyNum * (iShorter - 1) * sizeof(int));

	cudaMalloc((void**)&d_L1_vector, copyNum * sizeof(int));
	cudaMalloc((void**)&d_L1_strategy, copyNum * (iShorter - 1) * sizeof(int));
	
	cudaMemcpy(d_mtx_to_vec, mtx_to_vec, iRows * iCols * sizeof(int), cudaMemcpyHostToDevice);
printf("num_ofBlock: %d, num_ofThread: %d\n",num_ofBlock,num_ofThread);
	 func<<<num_ofBlock,num_ofThread>>>(d_mtx_to_vec, steps, steps_remainder, d_L1_vector, d_L1_strategy, iLonger, iShorter);
	cudaMemcpy(L1_vector, d_L1_vector, copyNum * sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(L1_strategy, d_L1_strategy, copyNum * (iShorter - 1) * sizeof(int), cudaMemcpyDeviceToHost);
	 
	for(i = 0; i < copyNum; i++){
		if(L1_max < L1_vector[i]) {L1_max = L1_vector[i]; iMax = i;}
	}
FILE *fp;
fp = fopen("strategy_L1.txt", "w");
	for(i=0; i<(iShorter - 1); i++) {fprintf(fp, "%d\n", L1_strategy[iMax * (iShorter - 1) + i]);}
	fprintf(fp,"1\n");
fclose(fp);

	printf("L1 is: %d\n",L1_max);

	free(L1_vector);
	free(L1_strategy);
	free(mtx_to_vec);

	cudaFree(d_L1_vector);
	cudaFree(d_L1_strategy);
	cudaFree(d_mtx_to_vec);

     return 0;     
}
