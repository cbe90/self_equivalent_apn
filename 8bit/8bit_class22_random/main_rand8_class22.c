/* Compile with: g++ -O2 main_rand8_class22.c -o rand8_class22
Authors: C. Beierle, G. Leander  -- Feb 2020

This program randomly searches for 8-bit APN permutations F s.t. AF = FB with A, B corresponding to class 22, i.e.,
A = B = Comp(X^4+X^3+X^2+X+1) \oplus Comp(X^4+X^3+X^2+X+1) 
up to affine equivalence. The program is optimized for that particular class.
*/

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include <random>

// dimension of the APN permutation to find
#define  N 8
#define RE_SHUFFLE 60	// defines the seconds after a fresh random state is generated

time_t start;
time_t runtime;

// the matrix A
int A[(1<<N)] = {0, 2, 4, 6, 8, 10, 12, 14, 15, 13, 11, 9, 7, 5, 3, 1, 32, 34, 36, 38, 40, 42, 44, 46, 47, 45, 43, 41, 39, 37, 35, 33, 64, 66, 68, 70, 72, 74, 76, 78, 79, 77, 75, 73, 71, 69, 67, 65, 96, 98, 100, 102, 104, 106, 108, 110, 111, 109, 107, 105, 103, 101, 99, 97, 128, 130, 132, 134, 136, 138, 140, 142, 143, 141, 139, 137, 135, 133, 131, 129, 160, 162, 164, 166, 168, 170, 172, 174, 175, 173, 171, 169, 167, 165, 163, 161, 192, 194, 196, 198, 200, 202, 204, 206, 207, 205, 203, 201, 199, 197, 195, 193, 224, 226, 228, 230, 232, 234, 236, 238, 239, 237, 235, 233, 231, 229, 227, 225, 240, 242, 244, 246, 248, 250, 252, 254, 255, 253, 251, 249, 247, 245, 243, 241, 208, 210, 212, 214, 216, 218, 220, 222, 223, 221, 219, 217, 215, 213, 211, 209, 176, 178, 180, 182, 184, 186, 188, 190, 191, 189, 187, 185, 183, 181, 179, 177, 144, 146, 148, 150, 152, 154, 156, 158, 159, 157, 155, 153, 151, 149, 147, 145, 112, 114, 116, 118, 120, 122, 124, 126, 127, 125, 123, 121, 119, 117, 115, 113, 80, 82, 84, 86, 88, 90, 92, 94, 95, 93, 91, 89, 87, 85, 83, 81, 48, 50, 52, 54, 56, 58, 60, 62, 63, 61, 59, 57, 55, 53, 51, 49, 16, 18, 20, 22, 24, 26, 28, 30, 31, 29, 27, 25, 23, 21, 19, 17};

// Inverse of A
int A_inv[(1<<N)] = {0, 15, 1, 14, 2, 13, 3, 12, 4, 11, 5, 10, 6, 9, 7, 8, 240, 255, 241, 254, 242, 253, 243, 252, 244, 251, 245, 250, 246, 249, 247, 248, 16, 31, 17, 30, 18, 29, 19, 28, 20, 27, 21, 26, 22, 25, 23, 24, 224, 239, 225, 238, 226, 237, 227, 236, 228, 235, 229, 234, 230, 233, 231, 232, 32, 47, 33, 46, 34, 45, 35, 44, 36, 43, 37, 42, 38, 41, 39, 40, 208, 223, 209, 222, 210, 221, 211, 220, 212, 219, 213, 218, 214, 217, 215, 216, 48, 63, 49, 62, 50, 61, 51, 60, 52, 59, 53, 58, 54, 57, 55, 56, 192, 207, 193, 206, 194, 205, 195, 204, 196, 203, 197, 202, 198, 201, 199, 200, 64, 79, 65, 78, 66, 77, 67, 76, 68, 75, 69, 74, 70, 73, 71, 72, 176, 191, 177, 190, 178, 189, 179, 188, 180, 187, 181, 186, 182, 185, 183, 184, 80, 95, 81, 94, 82, 93, 83, 92, 84, 91, 85, 90, 86, 89, 87, 88, 160, 175, 161, 174, 162, 173, 163, 172, 164, 171, 165, 170, 166, 169, 167, 168, 96, 111, 97, 110, 98, 109, 99, 108, 100, 107, 101, 106, 102, 105, 103, 104, 144, 159, 145, 158, 146, 157, 147, 156, 148, 155, 149, 154, 150, 153, 151, 152, 112, 127, 113, 126, 114, 125, 115, 124, 116, 123, 117, 122, 118, 121, 119, 120, 128, 143, 129, 142, 130, 141, 131, 140, 132, 139, 133, 138, 134, 137, 135, 136};

int pos[51] = {1,3,5,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95}; // one element of each cycle within the cycle structure of A (without 0). The next free position x will be selected from this list of positions

int repr[(1<<N)] = {0, 1, 1, 3, 1, 5, 3, 3, 1, 5, 5, 5, 3, 5, 3, 1, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 16, 31, 17, 30, 18, 29, 19, 28, 20, 27, 21, 26, 22, 25, 23, 24, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 16, 24, 31, 23, 17, 25, 30, 22, 18, 26, 29, 21, 19, 27, 28, 20, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 48, 63, 49, 62, 50, 61, 51, 60, 52, 59, 53, 58, 54, 57, 55, 56, 48, 52, 56, 60, 63, 59, 55, 51, 49, 53, 57, 61, 62, 58, 54, 50, 16, 20, 24, 28, 31, 27, 23, 19, 17, 21, 25, 29, 30, 26, 22, 18, 80, 84, 88, 92, 95, 91, 87, 83, 81, 85, 89, 93, 94, 90, 86, 82, 80, 95, 81, 94, 82, 93, 83, 92, 84, 91, 85, 90, 86, 89, 87, 88, 80, 88, 95, 87, 81, 89, 94, 86, 82, 90, 93, 85, 83, 91, 92, 84, 48, 56, 63, 55, 49, 57, 62, 54, 50, 58, 61, 53, 51, 59, 60, 52, 80, 82, 84, 86, 88, 90, 92, 94, 95, 93, 91, 89, 87, 85, 83, 81, 48, 50, 52, 54, 56, 58, 60, 62, 63, 61, 59, 57, 55, 53, 51, 49, 16, 18, 20, 22, 24, 26, 28, 30, 31, 29, 27, 25, 23, 21, 19, 17}; // the cycle representatives of each element from 0-255.

// numbers of even Hamming weight up to 255 (needed for more efficient APN check)
int evens[127] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126, 129, 130, 132, 135, 136, 139, 141, 142, 144, 147, 149, 150, 153, 154, 156, 159, 160, 163, 165, 166, 169, 170, 172, 175, 177, 178, 180, 183, 184, 187, 189, 190, 192, 195, 197, 198, 201, 202, 204, 207, 209, 210, 212, 215, 216, 219, 221, 222, 225, 226, 228, 231, 232, 235, 237, 238, 240, 243, 245, 246, 249, 250, 252, 255};

int P[(1<<N)];	// contains the random initial state

int solutions;
long long iterations;
int proceed;

int sbox[(1<<N)]; // it stores the (partial) S-box that we are going to construct (undefined values represented by -1)
int sbox_DDT[(1<<N)][(1<<N)]; // it stores the (partial) DDT of sbox

// random number generators to generate integers uniformly in [0,(1<<N))
std::random_device rd;
std::mt19937_64 gen(rd());
std::uniform_int_distribution<int> dis(0,(1<<N)-1);

// applys the Fisher-Yates shuffle to the array arr to get a uniform permutation
void shuffle_array(int* arr) {
    uint32_t j = 0;
    uint32_t tmp = 0;
    for (int i=(1<<N)-1; i!=0; i--) {
        decltype(dis.param()) nrange (0, i);
        dis.param(nrange);
        j = dis(gen);
        // swap
        tmp = arr[j];
        arr[j] = arr[i];
        arr[i] = tmp;
    }
}

void printArray(int A[1<<N]) {
	for (int i=0;i<(1<<N);i++)
		printf("%2x ",i);
	printf("\n");
	for (int i=0;i<(1<<N);i++)
		if (A[i]!=-1) printf("%2x ",A[i]);
			else printf("-- ");
	printf("\n\n");
}

void printSkippedArray(int A[128]) {
	for (int i=0;i<32;i++)
		printf("%d ",A[i]);
	printf("\n\n");
}

void printArray2(int A[1<<N]) {
	printf("{0x%02x",A[0]);
	for (int i=1;i<(1<<N);i++)
		printf(",0x%02x",A[i]);
	printf("};\n");
}

void fprintArray2(FILE *fp, int A[1<<N]) {
	fprintf(fp, "{0x%02x",A[0]);
	for (int i=1;i<(1<<N);i++)
		fprintf(fp, ",0x%02x",A[i]);
	fprintf(fp, "};\n");
}

void printDDT(int A[1<<N][1<<N]) {
    for (int i=0; i<(1<<N); i++) {
        for (int j=0; j<(1<<N); j++) {
            printf("%d ",A[i][j]);
        }
        printf("\n");
    }
}

// checks if A does not contain -1
int isComplete(int A[1<<N]) {
    int val = A[(1<<N)-1];
    int i;
    if (val==-1) return 0;
    A[(1<<N)-1] = -1;
    for (i=0; A[i] != -1; i++);
    A[(1<<N)-1] = val;
	return i==((1<<N)-1);

}

// return 1 if y is not yet in the image of sbox. returns 0 if it is already taken (i.e., check whether y is contained in sbox or not). It only checks whether the cycle representative is already contained. The function does not work for y=0, since it is not included in the pos array
int isNotTaken(int y) {
	for (int i=0; i<51; i++) {
		if (repr[sbox[pos[i]]] == repr[y])
			return 0;
	}
	return 1;
}

// continue building DDT for new set point sbox[c]. Returns 0 if the S-box cannot be APN anymore (because one entry is set to 4), returns 1 otherwise
int addDDTInformation(int c) {
    int a;
    // for all input differences
	for (int i=0;i<127;i++) {
		a = evens[i];	// we only need to check those input differences corresponding to vectors with even Hamming Weight
                // if sbox[c + a] is set, increase DDT entry by 2
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]+=2;
                    // if set above 2, return that it cannot be APN anymore
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]>2) return 0;
                }

	}
	return 1;
}

// reducing DDT entries for new removed point sbox[c]. It reduces all entries that were increased during the APN check of addDDTInformation(c).
void removeDDTInformation(int c) {
	int a;
	for (int i=0;i<127;i++) {
		a = evens[i];
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]-=2;
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]==2) break; // important because addDDTInformation(c) did not increment any values above that point.
                }

	}
}

int maxDepth=0;

// recursive construction of the sbox.
void nextValue(int depth, char* filename) {
	int xS,yS,j,y;
	iterations++;

	if (depth>maxDepth) maxDepth=depth; // to show maxDepth
	if ((time(NULL)-start)>=RE_SHUFFLE) { // display the status every RE_SHUFFLE seconds
		start=time(NULL);
		printf(" depth:%d\n",depth);
		printArray(sbox);
		printf("solutions so far:%d maxDepth:%d \n",solutions,maxDepth);
		printf("iterations:%lli \n",iterations);
		printf("running time: %li sec\n", time(NULL)-runtime);
		fflush(stdout);
		proceed = 0;
	}

	//complete
	if (isComplete(sbox)) {
			solutions++;
			printf("found a new apn permutation: #%d\n",solutions);
			printArray2(sbox);
			fflush(stdout);
			FILE *fp = fopen(filename, "a");
			if (fp == NULL)
            {
                printf("Error opening file!\n");
		fflush(stdout);
            }
            else
            {
                fprintArray2(fp, sbox);
            }
            fclose(fp);
		return;
	}

	//not complete

	int x=pos[depth];	// get next free position in the look-up table
	for (int z=0;z<(1<<N);z++) {		
		y = P[z];
		if (y!=0) { // do not include y=0 since this value is already set in the beginning. The function isNotTaken(y) would not work for 0.
		if (isNotTaken(y)) { // we know the order of each element is 5
			xS = x;
			yS = y;
			for (int i=0; i<4; i++) {	// we know that the order of x is 5
				// we set sbox[xS] = yS. Add this to the partial DDT information and check whether it can still be APN
				sbox[xS]=yS;
				if (!addDDTInformation(xS)) goto UNDO; // undo if it can't be APN anymore
				xS=A[xS];
                    		yS=A[yS];
			}
			sbox[xS]=yS;
			if (!addDDTInformation(xS)) goto UNDO; // undo if it can't be APN anymore
			
                    	if(proceed)	nextValue(depth+1, filename);
            
			//undo the changes (i.e., delete the set points of sbox from this step and undo the changes in the DDT)
			UNDO:
			// we remove the points in reverse order to undo the DDT changes correctly
			do { 
				removeDDTInformation(xS);
				sbox[xS] = -1;
				xS = A_inv[xS];
			} while (xS != A_inv[x]);
		}
		}
	}
}

void test_for_matrix(char* filename, int max_runtime) {
	if (N!=8) {
		printf("invalid N\n");
		exit(-1);
	}
    int xS, yS;
    solutions=0;
    iterations=0;
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
            printf("Error opening file!\n");
    }
    fclose(fp);
    printf("Search S-boxes invariant for matrix: \n");
    printArray2(A);
    printf("\n\n\n");
    fflush(stdout);
    for (int x=0;x<(1<<N);x++) {
	P[x] = x;
    }

    while(time(NULL)-runtime < max_runtime) {
	memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
	for (int x=0;x<(1<<N);x++) {
        	sbox[x]=-1;
    	}
	// w.l.o.g., we can set sbox[0]=0
    	sbox[0]=0;
    	if (!addDDTInformation(0)) exit(0);
    
    	// shuffle the state
    	shuffle_array(P);
    	printf("Re-shuffled intitial state ..\n");
    	fflush(stdout);
	
    	// start the recursive search
    	proceed = 1;
    	nextValue(0,filename);
    }
}


int main(int argc, char* argv[])
{
	if (argc != 3) {
        	printf("Usage: %s <outfile> <max_runtime> ...\n", argv[0]);
        	return -1;
    	}
	srand (time(NULL));
	start = time(NULL);
	runtime = time(NULL);
	
	test_for_matrix(argv[1],atoi(argv[2]));
	printf("\nFinished. Total running time: %li sec\n", time(NULL)-runtime);
	exit(0);
}
