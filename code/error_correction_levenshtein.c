#include <stdio.h>
#include <string.h>



#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))


// from https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance
int levenshtein(unsigned int s1len, unsigned int s2len, char *s1, char *s2) {
    unsigned int x, y, lastdiag, olddiag;
    unsigned int column[s1len + 1];
    for (y = 1; y <= s1len; y++)
        column[y] = y;
    for (x = 1; x <= s2len; x++) {
        column[0] = x;
        for (y = 1, lastdiag = x - 1; y <= s1len; y++) {
            olddiag = column[y];
            column[y] = MIN3(column[y] + 1, column[y - 1] + 1, lastdiag + (s1[y-1] == s2[x - 1] ? 0 : 1));
            lastdiag = olddiag;
        }
    }
    return column[s1len];
}


int hamming_distance(int l, char* a, char* b) {
    int d = 0; 
    int i = 0; 

    for(; i < l; i++)
        d += (__builtin_popcount(a[i] ^ b[i]) != 0) ? 1 : 0; 

    return d;
}

// arguments -> 
// 1: cc0 file
// 2: all file
// 3: out file
int main( int argc, char *argv[] ){
  unsigned int dist_hamm, dist;


  char * line_cc0 = NULL;
  size_t len_cc0 = 0;
  ssize_t read_cc0;
  

  FILE *file_cc0;
  file_cc0 = fopen(argv[1], "r");
  if (file_cc0 == NULL) {
    printf("Error opening file.\n");
  }
  
  char * line_all = NULL;
  size_t len_all = 0;
  ssize_t read_all;
  FILE *file_all;
  file_all = fopen(argv[2], "r");
  if (file_all == NULL) {
    printf("Error opening file.\n");
  }


  FILE *fp;
  fp = fopen(argv[3], "w");


  while ((read_cc0 = getline(&line_cc0, &len_cc0, file_cc0)) != -1) {

    while ((read_all = getline(&line_all, &len_all, file_all)) != -1) {
      line_all[strcspn(line_all, "\n")] = 0;
      line_cc0[strcspn(line_cc0, "\n")] = 0;
      dist_hamm = hamming_distance(20, line_all, line_cc0);
      if ((dist_hamm<10) && (dist_hamm!=0)) {
        dist = levenshtein(20, 20, line_all, line_cc0);
        if (dist<5) {
          fprintf(fp, "%s,%s,%u\n", line_all, line_cc0, dist);
        }
      }
    }
    rewind (file_all);
  }
  fclose(file_all);
  fclose(fp);
  fclose(file_cc0);
  
  return 1;
}  