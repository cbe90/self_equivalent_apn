/* Compile with: g++ -O2 main_det6_class5.c -o det6_class5
Authors: C. Beierle, G. Leander  -- Feb 2020

This program constructs and returns all 6-bit APN permutations F s.t. AF = FB with A, B corresponding to class 5, i.e.,
A = B = 
0 0 0 0 0 1
1 0 0 0 0 1
0 1 0 0 0 1
0 0 1 0 0 1
0 0 0 1 0 1
0 0 0 0 1 0
up to affine equivalence. It may return multiple representatives in one affine-equivalence class. 
*/

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include "commuting_matrices_A.h"
#include "commuting_matrices_B.h"

// dimension of the APN permutation to find
#define  N 6

time_t start;
time_t runtime;

// the matrix A
int A[(1<<N)] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5, 3, 1};

// the matrix B
int B[(1<<N)] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 63, 61, 59, 57, 55, 53, 51, 49, 47, 45, 43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11, 9, 7, 5, 3, 1};

// inverse of B
int B_inv[(1<<N)] = {0, 63, 1, 62, 2, 61, 3, 60, 4, 59, 5, 58, 6, 57, 7, 56, 8, 55, 9, 54, 10, 53, 11, 52, 12, 51, 13, 50, 14, 49, 15, 48, 16, 47, 17, 46, 18, 45, 19, 44, 20, 43, 21, 42, 22, 41, 23, 40, 24, 39, 25, 38, 26, 37, 27, 36, 28, 35, 29, 34, 30, 33, 31, 32};

// an ordering of the indices mod 2^n. It orders cycle-wise with respect to the cycles of B
int lex[(1<<N)] = {0,1, 2, 4, 8, 16, 32, 63,3, 6, 12, 24, 48, 31, 62,5, 10, 20, 40, 47, 33, 61,7, 14, 28, 56, 15, 30, 60,9, 18, 36, 55, 17, 34, 59,11, 22, 44, 39, 49, 29, 58,13, 26, 52, 23, 46, 35, 57,19, 38, 51, 25, 50, 27, 54,21, 42, 43, 41, 45, 37, 53};

int pos[9] = {1,3,5,7,9,11,13,19,21}; // one element of each cycle within the cycle structure of B (without 0). The next free position x will be selected from this list of positions (long cycles indexed first)

int repr[(1<<N)] = {0, 1, 1, 3, 1, 5, 3, 7, 1, 9, 5, 11, 3, 13, 7, 7, 1, 9, 9, 19, 5, 21, 11, 13, 3, 19, 13, 19, 7, 11, 7, 3, 1, 5, 9, 13, 9, 21, 19, 11, 5, 21, 21, 21, 11, 21, 13, 5, 3, 11, 19, 19, 13, 21, 19, 9, 7, 13, 11, 9, 7, 5, 3, 1}; // the cycle representatives of each element of A from 0-64.

// numbers of even Hamming weight up to 63 (needed for more efficient APN check)
int evens[31] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63};

int solutions;
long long iterations;
int n_skipped;
int skipped[128];

int sbox[(1<<N)]; // it stores the (partial) S-box that we are going to construct (undefined values represented by -1)
int sbox_DDT[(1<<N)][(1<<N)]; // it stores the (partial) DDT of sbox

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


// returns 1 if A is the smallest permutation up to conjugation with matrices in commuting_matrices. We use the lexicographic order of the lookup tables.
int is_smallest_in_class(int A[1<<N]){
    int E;
    for (int i=0;i<N_COMMUTING_MATRICES_B;i++) {
        for (int j=0; j<N_COMMUTING_MATRICES_A;j++) {
            for (int x=0;x<(1<<N);x++){
                if (A[lex[x]]==-1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                if (A[commuting_matrices_B[i][lex[x]]] == -1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                else {
                    E = commuting_matrices_A[j][A[commuting_matrices_B[i][lex[x]]]];
                    if (E < A[lex[x]]) return 0; // E is smaller, so A cannot be the smallest in class
                    if (E > A[lex[x]]) goto NEXT_MAT; // A is smaller, so consider next matrix
                } // consider next position
            }
            NEXT_MAT:;
        }
    }
    // A is the smallest of all, so return 1
    return 1;
}

// return 1 if y is not yet in the image of sbox. returns 0 if it is already taken (i.e., check whether y is contained in sbox or not). It only checks whether the cycle representative is already contained. The function does not work for y=0, since those are not included in the pos array
int isNotTaken(int y) {
	for (int i=0; i<9; i++) {
		if (repr[sbox[pos[i]]] == repr[y])
			return 0;
	}
	return 1;
}

// continue building DDT for new set point sbox[c]. Returns 0 if the S-box cannot be APN anymore (because one entry is set to 4), returns 1 otherwise
int addDDTInformation(int c) {
    int a;
    // for all input differences
	for (int i=0;i<31;i++) {
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
	for (int i=0;i<31;i++) {
		a = evens[i];
                if (sbox[c^a]!= -1) {
                    sbox_DDT[a][sbox[c]^sbox[c^a]]-=2;
                    if (sbox_DDT[a][sbox[c]^sbox[c^a]]==2) break; // important because addDDTInformation(c) did not increment any values above that point.
                }

	}
}

int maxDepth=0;

// recursive construction of the sbox.
void nextValue_with_filter(int depth, char* filename) {
	int xS,yS,j;
	iterations++;

	if (depth>maxDepth) maxDepth=depth; // to show maxDepth
	if ((time(NULL)-start)>=300) { // display the status every 300 seconds
		start=time(NULL);
		printf(" depth:%d\n",depth);
		printArray(sbox);
		printf("solutions so far:%d maxDepth:%d \n",solutions,maxDepth);
		printf("iterations:%lli  skipped: %d \nDepth skipped: \n",iterations,n_skipped);
		printSkippedArray(skipped);
		printf("running time: %li sec\n", time(NULL)-runtime);
		fflush(stdout);
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
	for (int y=(1<<N)-1;y>=1;y--) {		// do not include y=0 since those values are already set in the beginning. The function isNotTaken(y) would not work for 0.
		if (isNotTaken(y)) { // we know the order of each element is 7
			xS = x;
			yS = y;
			for (int i=0; i<6; i++) { // we know the order of each element is 7
				// we set sbox[xS] = yS. Add this to the partial DDT information and check whether it can still be APN
				sbox[xS]=yS;
				if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
				xS=B[xS];
                    		yS=A[yS];
			}
			sbox[xS]=yS;
			if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
			
                	if (depth <= 9) { // if it is too deep, it is faster to omit the check for smallest representative.
                        	if (is_smallest_in_class(sbox)) nextValue_with_filter(depth+1, filename);
                        	else {n_skipped++; skipped[depth]++;}
                	}
                	else { // omit check for smallest representative
                    		nextValue_with_filter(depth+1, filename);
                	}
            
			//undo the changes (i.e., delete the set points of sbox from this step and undo the changes in the DDT)
			UNDO_WITH_FILTER:
			// we remove the points in reverse order to undo the DDT changes correctly
			do { 
				removeDDTInformation(xS);
				sbox[xS] = -1;
				xS = B_inv[xS];
			} while (xS != B_inv[x]);
		}
	}
}

void test_for_matrix(char* filename) {
	if (N!=6) {
		printf("invalid N\n");
		exit(-1);
	}
    int xS, yS;
    solutions=0;
    iterations=0;
    memset(sbox_DDT,0,(1<<N)*(1<<N)*sizeof(int));
    for (int i=0; i<128; i++) skipped[i]=0;
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
            printf("Error opening file!\n");
    }
    fclose(fp);
    printf("Search S-boxes invariant for matrices: \n");
    printArray2(B);
    printArray2(A);
    printf("\n\n\n");
    fflush(stdout);
    for (int x=0;x<(1<<N);x++) {
        sbox[x]=-1;
    }

    // w.l.o.g., we can set sbox[0]=0
    sbox[0]=0;
    if (!addDDTInformation(0)) exit(0);

    // start the recursive search
    nextValue_with_filter(0, filename);
}

int main(int argc, char* argv[])
{
	if (argc != 2) {
        	printf("Usage: %s <outfile> ...\n", argv[0]);
        	return -1;
    	}
	srand (time(NULL));
	start = time(NULL);
	runtime = time(NULL);
	
	test_for_matrix(argv[1]);
	printf("\nFinished. Total running time: %li sec\n", time(NULL)-runtime);
	exit(0);
}
