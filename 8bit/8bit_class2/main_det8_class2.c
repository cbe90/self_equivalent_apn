/* Compile with: g++ -O2 main_det8_class2.c -o det8_class2
Authors: C. Beierle, G. Leander  -- Feb 2020

This program constructs and returns all 8-bit APN permutations F s.t. AF = FA with A corresponding to class 2, i.e.,
A = B = Comp(X^8+X^7+X^6+X^4+X^2+X+1)
up to affine equivalence. It may return multiple representatives in one affine-equivalence class. The program is optimized for that particular class.
*/

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include "commuting_matrices_8_class2.h"

// dimension of the APN permutation to find
#define  N 8

time_t start;
time_t runtime;

// the matrix A
int A[(1<<N)] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142, 144, 146, 148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 168, 170, 172, 174, 176, 178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 206, 208, 210, 212, 214, 216, 218, 220, 222, 224, 226, 228, 230, 232, 234, 236, 238, 240, 242, 244, 246, 248, 250, 252, 254, 215, 213, 211, 209, 223, 221, 219, 217, 199, 197, 195, 193, 207, 205, 203, 201, 247, 245, 243, 241, 255, 253, 251, 249, 231, 229, 227, 225, 239, 237, 235, 233, 151, 149, 147, 145, 159, 157, 155, 153, 135, 133, 131, 129, 143, 141, 139, 137, 183, 181, 179, 177, 191, 189, 187, 185, 167, 165, 163, 161, 175, 173, 171, 169, 87, 85, 83, 81, 95, 93, 91, 89, 71, 69, 67, 65, 79, 77, 75, 73, 119, 117, 115, 113, 127, 125, 123, 121, 103, 101, 99, 97, 111, 109, 107, 105, 23, 21, 19, 17, 31, 29, 27, 25, 7, 5, 3, 1, 15, 13, 11, 9, 55, 53, 51, 49, 63, 61, 59, 57, 39, 37, 35, 33, 47, 45, 43, 41};

// Inverse of A
int A_inv[(1<<N)] = {0, 235, 1, 234, 2, 233, 3, 232, 4, 239, 5, 238, 6, 237, 7, 236, 8, 227, 9, 226, 10, 225, 11, 224, 12, 231, 13, 230, 14, 229, 15, 228, 16, 251, 17, 250, 18, 249, 19, 248, 20, 255, 21, 254, 22, 253, 23, 252, 24, 243, 25, 242, 26, 241, 27, 240, 28, 247, 29, 246, 30, 245, 31, 244, 32, 203, 33, 202, 34, 201, 35, 200, 36, 207, 37, 206, 38, 205, 39, 204, 40, 195, 41, 194, 42, 193, 43, 192, 44, 199, 45, 198, 46, 197, 47, 196, 48, 219, 49, 218, 50, 217, 51, 216, 52, 223, 53, 222, 54, 221, 55, 220, 56, 211, 57, 210, 58, 209, 59, 208, 60, 215, 61, 214, 62, 213, 63, 212, 64, 171, 65, 170, 66, 169, 67, 168, 68, 175, 69, 174, 70, 173, 71, 172, 72, 163, 73, 162, 74, 161, 75, 160, 76, 167, 77, 166, 78, 165, 79, 164, 80, 187, 81, 186, 82, 185, 83, 184, 84, 191, 85, 190, 86, 189, 87, 188, 88, 179, 89, 178, 90, 177, 91, 176, 92, 183, 93, 182, 94, 181, 95, 180, 96, 139, 97, 138, 98, 137, 99, 136, 100, 143, 101, 142, 102, 141, 103, 140, 104, 131, 105, 130, 106, 129, 107, 128, 108, 135, 109, 134, 110, 133, 111, 132, 112, 155, 113, 154, 114, 153, 115, 152, 116, 159, 117, 158, 118, 157, 119, 156, 120, 147, 121, 146, 122, 145, 123, 144, 124, 151, 125, 150, 126, 149, 127, 148};

// an ordering of the indices mod 2^n. It orders cycle-wise with respect to the cycles of A
int lex[(1<<N)] = {0,1, 2, 4, 8, 16, 32, 64, 128, 215, 121, 242, 51, 102, 204, 79, 158, 235,3, 6, 12, 24, 48, 96, 192, 87, 174, 139, 193, 85, 170, 131, 209, 117, 234,5, 10, 20, 40, 80, 160, 151, 249, 37, 74, 148, 255, 41, 82, 164, 159, 233,7, 14, 28, 56, 112, 224, 23, 46, 92, 184, 167, 153, 229, 29, 58, 116, 232,9, 18, 36, 72, 144, 247, 57, 114, 228, 31, 62, 124, 248, 39, 78, 156, 239,11, 22, 44, 88, 176, 183, 185, 165, 157, 237, 13, 26, 52, 104, 208, 119, 238,15, 30, 60, 120, 240, 55, 110, 220, 111, 222, 107, 214, 123, 246, 59, 118, 236,17, 34, 68, 136, 199, 89, 178, 179, 177, 181, 189, 173, 141, 205, 77, 154, 227,19, 38, 76, 152, 231, 25, 50, 100, 200, 71, 142, 203, 65, 130, 211, 113, 226,21, 42, 84, 168, 135, 217, 101, 202, 67, 134, 219, 97, 194, 83, 166, 155, 225,27, 54, 108, 216, 103, 206, 75, 150, 251, 33, 66, 132, 223, 105, 210, 115, 230,35, 70, 140, 207, 73, 146, 243, 49, 98, 196, 95, 190, 171, 129, 213, 125, 250,43, 86, 172, 143, 201, 69, 138, 195, 81, 162, 147, 241, 53, 106, 212, 127, 254,45, 90, 180, 191, 169, 133, 221, 109, 218, 99, 198, 91, 182, 187, 161, 149, 253,47, 94, 188, 175, 137, 197, 93, 186, 163, 145, 245, 61, 122, 244, 63, 126, 252};

int pos[15] = {1,3,5,7,9,11,15,17,19,21,27,35,43,45,47}; // one element of each cycle within the cycle structure of A (without 0). The next free position x will be selected from this list of positions

int repr[(1<<N)] = {0, 1, 1, 3, 1, 5, 3, 7, 1, 9, 5, 11, 3, 11, 7, 15, 1, 17, 9, 19, 5, 21, 11, 7, 3, 19, 11, 27, 7, 7, 15, 9, 1, 27, 17, 35, 9, 5, 19, 9, 5, 5, 21, 43, 11, 45, 7, 47, 3, 35, 19, 1, 11, 43, 27, 15, 7, 9, 7, 15, 15, 47, 9, 47, 1, 19, 27, 21, 17, 43, 35, 19, 9, 35, 5, 27, 19, 17, 9, 1, 5, 43, 5, 21, 21, 3, 43, 3, 11, 17, 45, 45, 7, 47, 47, 35, 3, 21, 35, 45, 19, 21, 1, 27, 11, 27, 43, 15, 27, 45, 15, 15, 7, 19, 9, 27, 7, 3, 15, 11, 15, 1, 47, 15, 9, 35, 47, 43, 1, 35, 19, 3, 27, 45, 21, 21, 17, 47, 43, 3, 35, 17, 19, 43, 9, 47, 35, 43, 5, 45, 27, 5, 19, 7, 17, 21, 9, 11, 1, 5, 5, 45, 43, 47, 5, 11, 21, 7, 21, 45, 3, 35, 43, 17, 3, 47, 11, 17, 17, 17, 45, 17, 45, 11, 7, 11, 47, 45, 47, 17, 35, 45, 3, 3, 21, 43, 35, 47, 45, 17, 19, 43, 21, 19, 1, 17, 27, 35, 11, 3, 27, 19, 43, 35, 15, 1, 27, 21, 45, 21, 15, 45, 15, 27, 7, 21, 19, 17, 9, 7, 27, 19, 7, 5, 3, 1, 15, 11, 11, 9, 15, 43, 1, 35, 47, 47, 15, 9, 9, 5, 35, 27, 47, 45, 43, 5}; // the cycle representatives of each element of A from 0-255.

// numbers of even Hamming weight up to 255 (needed for more efficient APN check)
int evens[127] = {3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24, 27, 29, 30, 33, 34, 36, 39, 40, 43, 45, 46, 48, 51, 53, 54, 57, 58, 60, 63, 65, 66, 68, 71, 72, 75, 77, 78, 80, 83, 85, 86, 89, 90, 92, 95, 96, 99, 101, 102, 105, 106, 108, 111, 113, 114, 116, 119, 120, 123, 125, 126, 129, 130, 132, 135, 136, 139, 141, 142, 144, 147, 149, 150, 153, 154, 156, 159, 160, 163, 165, 166, 169, 170, 172, 175, 177, 178, 180, 183, 184, 187, 189, 190, 192, 195, 197, 198, 201, 202, 204, 207, 209, 210, 212, 215, 216, 219, 221, 222, 225, 226, 228, 231, 232, 235, 237, 238, 240, 243, 245, 246, 249, 250, 252, 255};

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
    for (int i=0;i<N_COMMUTING_MATRICES;i++) {
        for (int j=0; j<N_COMMUTING_MATRICES;j++) {
            for (int x=0;x<(1<<N);x++){
                if (A[lex[x]]==-1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                if (A[commuting_matrices[i][lex[x]]] == -1) goto NEXT_MAT; // We cannot say whether A is smaller than current matrix because of the undefined value
                else {
                    E = commuting_matrices[j][A[commuting_matrices[i][lex[x]]]];
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

// return 1 if y is not yet in the image of sbox. returns 0 if it is already taken (i.e., check whether y is contained in sbox or not). It only checks whether the cycle representative is already contained. The function does not work for y=0, since it is not included in the pos array
int isNotTaken(int y) {
	for (int i=0; i<15; i++) {
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
	for (int y=(1<<N)-1;y>=1;y--) {		// do not include y=0 since this value is already set in the beginning. The function isNotTaken(y) would not work for 0.
		if (isNotTaken(y)) { // we know the order of each element is 17
			xS = x;
			yS = y;
			for (int i=0; i<16; i++) {	// we know that the order of x is 17
				// we set sbox[xS] = yS. Add this to the partial DDT information and check whether it can still be APN
				sbox[xS]=yS;
				if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
				xS=A[xS];
                    		yS=A[yS];
			}
			sbox[xS]=yS;
			if (!addDDTInformation(xS)) goto UNDO_WITH_FILTER; // undo if it can't be APN anymore
			
                	if (depth <= 2) { // if it is too deep, it is faster to omit the check for smallest representative.
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
				xS = A_inv[xS];
			} while (xS != A_inv[x]);
		}
	}
}

void test_for_matrix(char* filename, int first, int second) {
	if (N!=8) {
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
    printf("Search S-boxes invariant for matrix: \n");
    printArray2(A);
    printf("\n\n\n");
    fflush(stdout);
    for (int x=0;x<(1<<N);x++) {
        sbox[x]=-1;
    }

    // w.l.o.g., we can set sbox[0]=0
    sbox[0]=0;
    if (!addDDTInformation(0)) exit(0);

    sbox[1]=first;
    if (!addDDTInformation(1)) exit(0);
    xS=A[1];
    yS=A[first];
    // set the orbits accordingly
    for (int i=0;i<16;i++) {
                    sbox[xS]=yS;
                    if (!addDDTInformation(xS)) exit(0);
                    xS=A[xS];
                    yS=A[yS];
    }
    sbox[3]=second;
    if (!addDDTInformation(3)) exit(0);
    xS=A[3];
    yS=A[second];
    // set the orbits accordingly
    for (int i=0;i<16;i++) {
                    sbox[xS]=yS;
                    if (!addDDTInformation(xS)) exit(0);
                    xS=A[xS];
                    yS=A[yS];
    }
    if (!is_smallest_in_class(sbox)) {
	exit(0);
    }

    // start the recursive search
    nextValue_with_filter(2, filename);
}


int main(int argc, char* argv[])
{
	if (argc != 4) {
        	printf("Usage: %s <outfile> <SBox[1]> <SBox[3]> ...\n", argv[0]);
        	return -1;
    	}
	srand (time(NULL));
	start = time(NULL);
	runtime = time(NULL);
	
	
	if ((atoi(argv[2]) == 0) or (atoi(argv[3]) == 0)) {
		printf("No permutation ..\n");
		return -1;
	}
	if (repr[atoi(argv[2])] == repr[atoi(argv[3])]) {
		printf("No permutation ..\n");
		return -1;
	}
	test_for_matrix(argv[1],atoi(argv[2]),atoi(argv[3]));
	printf("\nFinished. Total running time: %li sec\n", time(NULL)-runtime);
	exit(0);
}
