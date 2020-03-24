/* Compile with: g++ -O2 main_rand7_class23.c -o rand7_class23
Authors: C. Beierle, G. Leander  -- Feb 2020

This program randomly searches for 7-bit APN permutations F s.t. AF = FB with A, B corresponding to class 23, i.e.,
A = B = Comp(X^2+X+1) \oplus Comp(X^2+X+1) \oplus Comp(X^3+1)
up to affine equivalence. It may return multiple representatives in one affine-equivalence class. 
*/


#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include <random>

// dimension of the APN permutation to find
#define  N 7
#define RE_SHUFFLE 60	// defines the seconds after a fresh random state is generated

time_t start;
time_t runtime;

// the matrix A
int A[(1<<N)] = {0, 2, 3, 1, 8, 10, 11, 9, 12, 14, 15, 13, 4, 6, 7, 5, 32, 34, 35, 33, 40, 42, 43, 41, 44, 46, 47, 45, 36, 38, 39, 37, 64, 66, 67, 65, 72, 74, 75, 73, 76, 78, 79, 77, 68, 70, 71, 69, 96, 98, 99, 97, 104, 106, 107, 105, 108, 110, 111, 109, 100, 102, 103, 101, 16, 18, 19, 17, 24, 26, 27, 25, 28, 30, 31, 29, 20, 22, 23, 21, 48, 50, 51, 49, 56, 58, 59, 57, 60, 62, 63, 61, 52, 54, 55, 53, 80, 82, 83, 81, 88, 90, 91, 89, 92, 94, 95, 93, 84, 86, 87, 85, 112, 114, 115, 113, 120, 122, 123, 121, 124, 126, 127, 125, 116, 118, 119, 117};

// Inverse of A
int A_inv[(1<<N)] = {0, 3, 1, 2, 12, 15, 13, 14, 4, 7, 5, 6, 8, 11, 9, 10, 64, 67, 65, 66, 76, 79, 77, 78, 68, 71, 69, 70, 72, 75, 73, 74, 16, 19, 17, 18, 28, 31, 29, 30, 20, 23, 21, 22, 24, 27, 25, 26, 80, 83, 81, 82, 92, 95, 93, 94, 84, 87, 85, 86, 88, 91, 89, 90, 32, 35, 33, 34, 44, 47, 45, 46, 36, 39, 37, 38, 40, 43, 41, 42, 96, 99, 97, 98, 108, 111, 109, 110, 100, 103, 101, 102, 104, 107, 105, 106, 48, 51, 49, 50, 60, 63, 61, 62, 52, 55, 53, 54, 56, 59, 57, 58, 112, 115, 113, 114, 124, 127, 125, 126, 116, 119, 117, 118, 120, 123, 121, 122};

int pos[42] = {1,4,5,6,7,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,113,116,117,118,119}; // one element of each cycle within the cycle structure of A (without 0 and 112). The next free position x will be selected from this list of positions (long cycles indexed first)

int repr[(1<<N)] = {0, 1, 1, 1, 4, 5, 6, 7, 4, 7, 5, 6, 4, 6, 7, 5, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 16, 19, 17, 18, 28, 31, 29, 30, 20, 23, 21, 22, 24, 27, 25, 26, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 16, 18, 19, 17, 24, 26, 27, 25, 28, 30, 31, 29, 20, 22, 23, 21, 48, 50, 51, 49, 56, 58, 59, 57, 60, 62, 63, 61, 52, 54, 55, 53, 48, 51, 49, 50, 60, 63, 61, 62, 52, 55, 53, 54, 56, 59, 57, 58, 112, 113, 113, 113, 116, 117, 118, 119, 116, 119, 117, 118, 116, 118, 119, 117}; // the cycle representatives of each element from 0-127.

// numbers of even Hamming weight up to 127 (needed for more efficient APN check)
int evens[63] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126};

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

// return 1 if y is not yet in the image of sbox. returns 0 if it is already taken (i.e., check whether y is contained in sbox or not). It only checks whether the cycle representative is already contained. The function does not work for the fixed points, since they are not included in the pos array
int isNotTaken(int y) {
	for (int i=0; i<42; i++) {
		if (repr[sbox[pos[i]]] == repr[y])
			return 0;
	}
	return 1;
}

// continue building DDT for new set point sbox[c]. Returns 0 if the S-box cannot be APN anymore (because one entry is set to 4), returns 1 otherwise
int addDDTInformation(int c) {
    int a;
    // for all input differences
	for (int i=0;i<63;i++) {
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
	for (int i=0;i<63;i++) {
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
		if ((y!=0) && (y!=112)) { // fixed points 0 and 112 are already set in the beginning
		if (isNotTaken(y)) { // we know the order of each element is 3
			xS = x;
			yS = y;
			for (int i=0; i<2; i++) {	// we know that the order of x is 3
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
	if (N!=7) {
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
	// w.l.o.g., we can set the sbox on the fixed points 0,112 to be the identity
    	sbox[0]=0;
    	if (!addDDTInformation(0)) exit(0);
    	sbox[112]=112;
    	if (!addDDTInformation(112)) exit(0);
    
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
