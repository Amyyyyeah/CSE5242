/*
  CSE 5242 Project 2, Fall 2024

  See class project handout for more extensive documentation.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/time.h>
#include <asm/unistd.h>
#include <immintrin.h>

/* uncomment out the following DEBUG line for debug info, for experiment comment the DEBUG line  */
// #define DEBUG

/* This will define whether it will run the band_join or band_join_simd function. 
  If you want to run band_join_simd, comment out this line  */
#define BAND_JOIN_EXPERIMENT

/* compare two int64_t values - for use with qsort */
static int compare(const void *p1, const void *p2)
{
  int a,b;
  a = *(int64_t *)p1;
  b = *(int64_t *)p2;
  if (a<b) return -1;
  if (a==b) return 0;
  return 1;
}

/* initialize searches and data - data is sorted and searches is a random permutation of data */
int init(int64_t* data, int64_t* searches, int count)
{
  for(int64_t i=0; i<count; i++){
    searches[i] = random();
    data[i] = searches[i]+1;
  }
  qsort(data,count,sizeof(int64_t),compare);
}

/* initialize outer probes of band join */
int band_init(int64_t* outer, int64_t size)
{
  for(int64_t i=0; i<size; i++){
    outer[i] = random();
  }
}

inline int64_t simple_binary_search(int64_t* data, int64_t size, int64_t target)
{
  int64_t left=0;
  int64_t right=size;
  int64_t mid;

  while(left<=right) {
    mid = (left + right)/2;   /* ignore possibility of overflow of left+right */
    if (data[mid]==target) return mid;
    if (data[mid]<target) left=mid+1;
    else right = mid-1;
  }
  return -1; /* no match */
}

inline int64_t low_bin_search(int64_t* data, int64_t size, int64_t target)
{
  /* this binary search variant
     (a) does only one comparison in the inner loop
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
	 identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
  */
  int64_t left=0;
  int64_t right=size;
  int64_t mid;

  while(left<right) {
    mid = (left + right)/2;   /* ignore possibility of overflow of left+right */
    if (data[mid]>=target)
      right=mid;
    else
      left=mid+1;
  }
  return right;
}

/* (a) question*/
inline int64_t low_bin_nb_arithmetic(int64_t* data, int64_t size, int64_t target)
{
  /* this binary search variant
     (a) does no comparisons in the inner loop by using multiplication and addition to convert control dependencies
         to data dependencies
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
	 identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
  */
  int64_t left=0;
  int64_t right=size;
  int64_t mid;

  while(left<right) {

    mid = (left + right) / 2;
    int64_t is_greater_equal = (data[mid] >= target);
    right = mid * is_greater_equal + right * (1 - is_greater_equal); 
    left = left * is_greater_equal + (mid + 1) * (1 - is_greater_equal);
  }
  return right;
}

inline int64_t low_bin_nb_mask(int64_t* data, int64_t size, int64_t target)
{
  /* this binary search variant
     (a) does no comparisons in the inner loop by using bit masking operations to convert control dependencies
         to data dependencies
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
	 identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
  */
  int64_t left=0;
  int64_t right=size;
  int64_t mid;
  int64_t update_right;

  while(left<right) {
    mid = (left + right)/2;
    update_right = (data[mid] >= target);

    // If the data[mid] >= target the right value should be updated and left value should remain the same.
    // But if it is not, then the right value should remain the same and left value should be updated
    right = ((0 - update_right) & mid) | ((0 - (!update_right)) & right);
    left = ((0 - (!update_right)) & (mid + 1)) | ((0 - update_right) & left);
  }
  return right;
}

void low_bin_nb_4x(int64_t* data, int64_t size, int64_t* targets, int64_t* right)
{
  /* this binary search variant
     (a) does no comparisons in the inner loop by using bit masking operations instead
     (b) doesn't require an exact match; instead it returns the index of the first key >= the search key.
         That's good in a DB context where we might be doing a range search, and using binary search to
	 identify the first key in the range.
     (c) If the search key is bigger than all keys, it returns size.
     (d) does 4 searches at the same time in an interleaved fashion, so that an out-of-order processor
         can make progress even when some instructions are still waiting for their operands to be ready.

     Note that we're using the result array as the "right" elements in the search so no need for a return statement
  */

    /* YOUR CODE HERE */

  int64_t left[4] = {0, 0, 0, 0};
  int64_t update_right[4];
  int64_t mid[4];
  int64_t quit_iter[4] = {0, 0, 0, 0};
  int64_t i;

  right[0] = size;
  right[1] = size;
  right[2] = size;
  right[3] = size;

  while (quit_iter[0] + quit_iter[1] + quit_iter[2] + quit_iter[3] != 4) {
    for (i = 0; i < 4; i++) {
      mid[i] = (left[i] + right[i])/2;
      update_right[i] = (data[mid[i]] >= targets[i]);
      quit_iter[i] = (left[i] >= right[i]);

      right[i] = ((0 - (update_right[i] & !quit_iter[i])) & mid[i]) | ((0 - (!update_right[i] | quit_iter[i])) & right[i]);
      left[i] = ((0 - (!update_right[i] & !quit_iter[i])) & (mid[i] + 1)) | ((0 - (update_right[i] | quit_iter[i])) & left[i]);
    }
  }

}


/* The following union type is handy to output the contents of AVX512 data types */
union int8x4 {
  __m256i a;
  int64_t b[4];
};

void printavx(char* name, __m256i v) {
  union int8x4 n;

  n.a=v;
  printf("Value in %s is [%ld %ld %ld %ld ]\n",name,n.b[0],n.b[1],n.b[2],n.b[3]);
}

/*
 * Optinal for using AVX-512

  union int8x8 {
    __m512i a;
    int64_t b[8];
  };

  void printavx512(char* name, __m512i v) {
    union int8x4 n;

    n.a=v;
    printf("Value in %s is [%ld %ld %ld %ld %ld %ld %ld %ld ]\n",name,n.b[0],n.b[1],n.b[2],n.b[3]);
  }

 */

inline void low_bin_nb_simd(int64_t* data, int64_t size, __m256i target, __m256i* result) 
{
    __m256i left_pointer = _mm256_set1_epi64x(0);
    __m256i right_pointer = _mm256_set1_epi64x(size-1);
    __m256i one = _mm256_set1_epi64x(1);

    while (_mm256_movemask_epi8(_mm256_cmpgt_epi64(right_pointer, left_pointer))) {
        __m256i mid = _mm256_add_epi64(left_pointer, right_pointer);
        mid = _mm256_srli_epi64(mid, 1);  
        __m256i mid_values = _mm256_i64gather_epi64((const long long int *)data, mid, 8);
        __m256i is_greater = _mm256_cmpgt_epi64(target, mid_values);
      
        left_pointer = _mm256_blendv_epi8(left_pointer, _mm256_add_epi64(mid, one), is_greater);
        right_pointer = _mm256_blendv_epi8(right_pointer, mid, _mm256_xor_si256(is_greater, _mm256_set1_epi64x(-1)));
    }

    *result = left_pointer;
}


void bulk_bin_search(int64_t* data, int64_t size, int64_t* searchkeys, int64_t numsearches, int64_t* results, int repeats)
{
  for(int j=0; j<repeats; j++) {
    /* Function to test a large number of binary searches

       we might need repeats>1 to make sure the events we're measuring are not dominated by various
       overheads, particularly for small values of size and/or numsearches

       we assume that we want exactly "size" searches, where "size" is the length if the searchkeys array
     */
    for(int64_t i=0;i<numsearches; i++) {
#ifdef DEBUG
      printf("Searching for %ld...\n",searchkeys[i]);
#endif

      // Uncomment one of the following to measure it
      results[i] = low_bin_search(data,size,searchkeys[i]);
      //results[i] = low_bin_nb_arithmetic(data,size,searchkeys[i]);
      //results[i] = low_bin_nb_mask(data,size,searchkeys[i]);

#ifdef DEBUG
      printf("Result is %ld\n",results[i]);
#endif
    }
  }
}

void bulk_bin_search_4x(int64_t* data, int64_t size, int64_t* searchkeys, int64_t numsearches, int64_t* results, int repeats)
{
  register __m256i searchkey_4x;

  for(int j=0; j<repeats; j++) {
    /* Function to test a large number of binary searches using one of the 8x routines

       we might need repeats>1 to make sure the events we're measuring are not dominated by various
       overheads, particularly for small values of size and/or numsearches

       we assume that we want exactly "size" searches, where "size" is the length if the searchkeys array
     */
    int64_t extras = numsearches % 4;
    for(int64_t i=0;i<numsearches-extras; i+=4) {
#ifdef DEBUG
      printf("Searching for %ld %ld %ld %ld  ...\n",
	     searchkeys[i],searchkeys[i+1],searchkeys[i+2],searchkeys[i+3]);
#endif

      // Uncomment one of the following depending on which routine you want to profile

      // Algorithm A
      low_bin_nb_4x(data,size,&searchkeys[i],&results[i]);

      // Algorithm B
      // searchkey_4x = _mm256_loadu_si256((__m256i *)&searchkeys[i]);
      // low_bin_nb_simd(data,size,searchkey_4x,(__m256i *)&results[i]);

#ifdef DEBUG
      printf("Result is %ld %ld %ld %ld  ...\n",
	     results[i],results[i+1],results[i+2],results[i+3]);
#endif
    }
    /* a little bit more work if numsearches is not a multiple of 8 */
    for(int64_t i=numsearches-extras;i<numsearches; i++) {

      results[i] = low_bin_nb_mask(data,size,searchkeys[i]);

    }

  }
}


int64_t band_join(int64_t* inner, int64_t inner_size, int64_t* outer, int64_t outer_size, int64_t* inner_results, int64_t* outer_results, int64_t result_size, int64_t bound)
{
  /* In a band join we want matches within a range of values.  If p is the probe value from the outer table, then all
     records in the inner table with a key in the range [p-bound,p+bound] inclusive should be part of the result.

     Results are returned via two arrays. outer_results stores the index of the outer table row that matches, and
     inner_results stores the index of the inner table row that matches.  result_size tells you the size of the
     output array that has been allocated. You should make sure that you don't exceed this size.  If there are
     more results than can fit in the result arrays, then return early with just a prefix of the results in the result
     arrays. The return value of the function should be the number of output results.

  */

  /* YOUR CODE HERE */
  int64_t i, j, k;
  int64_t extras = result_size % 4;
  int ret;

  // This will contain the index for the lower and upper bound within the inner array
  // int64_t lower[4], upper[4], index_low[4], index_up[4]; 
  int64_t total_cnt = 0;
  int64_t start, end;

  i = 0;
  while ((total_cnt < result_size) && (i <= (outer_size - 4))) {
    int64_t lower[4], upper[4], index_low[4], index_up[4]; 

    lower[0] = outer[i] - bound;
    upper[0] = outer[i] + bound;

    lower[1] = outer[i + 1] - bound;
    upper[1] = outer[i + 1] + bound;

    lower[2] = outer[i + 2] - bound;
    upper[2] = outer[i + 2] + bound;

    lower[3] = outer[i + 3] - bound;
    upper[3] = outer[i + 3] + bound;

    low_bin_nb_4x(inner, inner_size, &lower[0], &index_low[0]);
    low_bin_nb_4x(inner, inner_size, &upper[0], &index_up[0]);

    for (j = 0; j < 4; j++) {
      for (k = index_low[j]; k <= index_up[j]; k++) {
        if (total_cnt >= result_size) return total_cnt;

        inner_results[total_cnt] = k;
        outer_results[total_cnt] = i + j;
        total_cnt++;
      }  
    }
    
    i += 4;
  }

  for(i = result_size - extras; (i < result_size) && (i < outer_size); i++) {
    start = low_bin_nb_mask(inner, inner_size, outer[i] - bound);
    end = low_bin_nb_mask(inner, inner_size, outer[i] + bound);

    for (j = start; j <= end; j++) {
        if (total_cnt >= result_size) return total_cnt;
        
        inner_results[total_cnt] = j;
        outer_results[total_cnt] = i;
        total_cnt++;
    }
  }

  return total_cnt;
}

int64_t band_join_simd(int64_t* inner, int64_t inner_size, int64_t* outer, int64_t outer_size, int64_t* inner_results, int64_t* outer_results, int64_t result_size, int64_t bound)
{
  int64_t match_count = 0;
  int64_t i, j;
  int64_t start, end, idx;
  int64_t lower, upper;

  __m256i outer_vals;
  __m256i lower_bound;
  __m256i upper_bound;
  __m256i lower_idx;
  __m256i upper_idx;
  __m256i bound_vec = _mm256_set1_epi64x(bound);

  for (i = 0; i <= outer_size - 4; i += 4) {
    outer_vals = _mm256_loadu_si256((__m256i*)&outer[i]); 
    lower_bound = _mm256_sub_epi64(outer_vals, bound_vec); 
    upper_bound = _mm256_add_epi64(outer_vals, bound_vec); 

    
    low_bin_nb_simd(inner, inner_size, lower_bound, &lower_idx);
    low_bin_nb_simd(inner, inner_size, upper_bound, &upper_idx);


    for (j = 0; j < 4; j++) {
      start = ((int64_t*)&lower_idx)[j];
      end = ((int64_t*)&upper_idx)[j];

      for (idx = start; idx <= end; idx++) {
          if (match_count >= result_size) return match_count;
          inner_results[match_count] = idx;
          outer_results[match_count] = i + j;
          match_count++;
      }
    }
  }

  for (i = outer_size - outer_size % 4; i < outer_size; i++) {
    lower = outer[i] - bound;
    upper = outer[i] + bound;
    start = low_bin_nb_mask(inner, inner_size, lower);
    end = low_bin_nb_mask(inner, inner_size, upper);

    for (idx = start; idx <= end; idx++) {
        if (match_count >= result_size) return match_count;
        inner_results[match_count] = idx;
        outer_results[match_count] = i;
        match_count++;
    }
  }

  return match_count;
}


int main(int argc, char *argv[])
{
	 long long counter;
	 int64_t arraysize, outer_size, result_size;
	 int64_t bound;
	 int64_t *data, *queries, *results;
	 int ret;
	 struct timeval before, after;
	 int repeats;
	 int64_t total_results;

	 // band-join arrays
	 int64_t *outer, *outer_results, *inner_results;


	 if (argc >= 5)
	   {
	     arraysize = atoi(argv[1]);
	     outer_size = atoi(argv[2]);
	     result_size = atoi(argv[3]);
	     bound = atoi(argv[4]);
	   }
	 else
	   {
	     fprintf(stderr, "Usage: db5242 inner_size outer_size result_size bound <repeats>\n");
	     exit(EXIT_FAILURE);
	   }

	 if (argc >= 6)
	   repeats = atoi(argv[5]);
	 else
	   {
	     repeats=1;
	   }

	 /* allocate the array and the queries for searching */
	 ret=posix_memalign((void**) &data,64,arraysize*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }
	 ret=posix_memalign((void**) &queries,64,arraysize*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }
	 ret=posix_memalign((void**) &results,64,arraysize*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }

	 /* allocate the outer array and output arrays for band-join */
	 ret=posix_memalign((void**) &outer,64,outer_size*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }
	 ret=posix_memalign((void**) &outer_results,64,result_size*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }
	 ret=posix_memalign((void**) &inner_results,64,result_size*sizeof(int64_t));
	 if (ret)
	 {
	   fprintf(stderr, "Memory allocation error.\n");
	   exit(EXIT_FAILURE);
	 }


	   /* code to initialize data structures goes here so that it is not included in the timing measurement */
	   init(data,queries,arraysize);
	   band_init(outer,outer_size);

#ifdef DEBUG
	   /* show the arrays */
	   printf("data: ");
	   for(int64_t i=0;i<arraysize;i++) printf("%ld ",data[i]);
	   printf("\n");
	   printf("queries: ");
	   for(int64_t i=0;i<arraysize;i++) printf("%ld ",queries[i]);
	   printf("\n");
	   printf("outer: ");
	   for(int64_t i=0;i<outer_size;i++) printf("%ld ",outer[i]);
	   printf("\n");
#endif


	   /* now measure... */

	   gettimeofday(&before,NULL);

	   /* the code that you want to measure goes here; make a function call */
	   bulk_bin_search(data,arraysize,queries,arraysize,results, repeats);

	   gettimeofday(&after,NULL);
	   printf("Time in bulk_bin_search loop is %ld microseconds or %f microseconds per search\n", (after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec), 1.0*((after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec))/arraysize/repeats);

	   gettimeofday(&before,NULL);

	   /* the code that you want to measure goes here; make a function call */
	   bulk_bin_search_4x(data,arraysize,queries,arraysize,results, repeats);

	   gettimeofday(&after,NULL);
	   printf("Time in bulk_bin_search_4x loop is %ld microseconds or %f microseconds per search\n", (after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec), 1.0*((after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec))/arraysize/repeats);

	   gettimeofday(&before,NULL);

	   /* the code that you want to measure goes here; make a function call */
#ifdef BAND_JOIN_EXPERIMENT
	   total_results=band_join(data, arraysize, outer, outer_size, inner_results, outer_results, result_size, bound);
#else
     total_results=band_join_simd(data, arraysize, outer, outer_size, inner_results, outer_results, result_size, bound);
#endif

	   gettimeofday(&after,NULL);
	   printf("Band join result size is %ld with an average of %f matches per output record\n",total_results, 1.0*total_results/(1.0+outer_results[total_results-1]));
	   printf("Time in band_join loop is %ld microseconds or %f microseconds per outer record\n", (after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec), 1.0*((after.tv_sec-before.tv_sec)*1000000+(after.tv_usec-before.tv_usec))/outer_size);

#ifdef DEBUG
	   /* show the band_join results */
	   printf("band_join results: ");
	   for(int64_t i=0;i<total_results;i++) printf("(%ld,%ld) ",outer_results[i],inner_results[i]);
	   printf("\n");
#endif

}