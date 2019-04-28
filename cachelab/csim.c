/* 
 * csim.c - A cache simulator that can replay traces from Valgrind
 *     and output statistics such as number of hits, misses, and
 *     evictions.  The replacement policy is LRU.
 *
 * Implementation and assumptions:
 *  1. Each load/store can cause at most one cache miss. (I examined the trace,
 *  the largest request I saw was for 8 bytes).
 *  2. Instruction loads (I) are ignored, since we are interested in evaluating
 *  trans.c in terms of its data cache performance.
 *  3. data modify (M) is treated as a load followed by a store to the same
 *  address. Hence, an M operation can result in two cache hits, or a miss and a
 *  hit plus an possible eviction.
 *
 * The function printSummary() is given to print output.
 * Please use this function to print the number of hits, misses and evictions.
 * This is crucial for the driver to evaluate your work. 
 */
#include <getopt.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include <string.h>
#include "cachelab.h"

//#define DEBUG_ON 
#define ADDRESS_LENGTH 64

/* Type: Memory address */
typedef unsigned long long int mem_addr_t;

/* Type: Cache line
   LRU is a counter used to implement LRU replacement policy  */
typedef struct cache_line {
    mem_addr_t tag;
    struct cache_line *next;
} cache_line_t;

typedef struct cache_set {
    cache_line_t * line;
    int size;
} cache_set_t;

/* Globals set by command line args */
int verbosity = 0; /* print trace if set */
int s = 0; /* set index bits */
int b = 0; /* block offset bits */
int E = 0; /* associativity */
char* trace_file = NULL;

/* Derived from command line args */
unsigned int S; /* number of sets */
unsigned int B; /* block size (bytes) */

/* Counters used to record cache statistics */
int miss_count = 0;
int hit_count = 0;
int eviction_count = 0;
// unsigned long long int lru_counter = 1;

/* The cache we are simulating */
cache_set_t *cache;  
mem_addr_t set_index_mask;

/* 
 * initCache - Allocate memory, write 0's for valid and tag and LRU
 * also computes the set_index_mask
 */
void initCache()
{
    int i;
    cache = (cache_set_t*) malloc(S * sizeof(cache_set_t));
    for (i=0; i<S; i++){
        cache[i].line = NULL;
        cache[i].size = 0;
    }

    /* Computes set index mask */
    set_index_mask = (mem_addr_t) (S - 1);
}


/* 
 * freeCache - free allocated memory
 */
void freeCache()
{
    int i;
    for(i = 0; i < S; ++i) {
        if(cache[i].line) free(cache[i].line);
    }
    free(cache);
}


/* 
 * accessData - Access data at memory address addr.
 *   If it is already in cache, increast hit_count
 *   If it is not in cache, bring it in cache, increase miss count.
 *   Also increase eviction_count if a line is evicted.
 */
void accessData(mem_addr_t addr)
{
    mem_addr_t block_addr = (addr >> (b+s));
    mem_addr_t set_index = (addr >> b) & (S-1);
    int flag = 0;
    for(cache_line_t **p = &cache[set_index].line; *p; p = &((*p)->next)){
        // printf("\n%pd\n", *p);
        if((*p)->tag == block_addr){
            flag = 1 | flag;
            if(verbosity) printf(" hit");
            hit_count++;
            if(*p == cache[set_index].line) break;//already the head
            cache_line_t *tmp = (*p)->next;
            (*p)->next = cache[set_index].line;
            cache[set_index].line = *p;
            *p = tmp;
            break;
        } else if((*p)->next == NULL && cache[set_index].size == E) {
            free(*p);
            *p = NULL;
            break;
        }
    }
    if(!flag){//not find
        if(verbosity) printf(" miss");
        miss_count++;
        cache_line_t *newline = (cache_line_t*)malloc(sizeof(cache_line_t) + B);
        newline->tag = block_addr;
        newline->next = cache[set_index].line;
        cache[set_index].line = newline;
        if(cache[set_index].size == E){//already overflow
            if(verbosity) printf(" eviction");
            eviction_count++;
        } else {
            cache[set_index].size++;
        }
    }
}


/*
 * replayTrace - replays the given trace file against the cache 
 */
void replayTrace(char* trace_fn)
{
    FILE* trace_fp = fopen(trace_fn, "r");
    char op;
    unsigned long long address;
    int sz;
    while(fscanf(trace_fp, " %c%llx,%d", &op, &address,&sz) == 3){
        if(op == 'I') continue;
        else{
            if(verbosity) 
                printf("%c %llx,%d", op, address, sz);
            if(op == 'L' || op == 'S'){
                accessData(address);
            }else if(op == 'M'){
                accessData(address);
                accessData(address);
            }
        }   
        if(verbosity) putchar('\n');
    }
    fclose(trace_fp);
}

/*
 * printUsage - Print usage info
 */
void printUsage(char* argv[])
{
    printf("Usage: %s [-hv] -s <num> -E <num> -b <num> -t <file>\n", argv[0]);
    printf("Options:\n");
    printf("  -h         Print this help message.\n");
    printf("  -v         Optional verbose flag.\n");
    printf("  -s <num>   Number of set index bits.\n");
    printf("  -E <num>   Number of lines per set.\n");
    printf("  -b <num>   Number of block offset bits.\n");
    printf("  -t <file>  Trace file.\n");
    printf("\nExamples:\n");
    printf("  linux>  %s -s 4 -E 1 -b 4 -t traces/yi.trace\n", argv[0]);
    printf("  linux>  %s -v -s 8 -E 2 -b 4 -t traces/yi.trace\n", argv[0]);
    exit(0);
}

/*
 * main - Main routine 
 */
int main(int argc, char* argv[])
{
    char c;

    while( (c=getopt(argc,argv,"s:E:b:t:vh")) != -1){
        switch(c){
        case 's':
            s = atoi(optarg);
            S = 2 << s;
            break;
        case 'E':
            E = atoi(optarg);
            break;
        case 'b':
            b = atoi(optarg);
            B = 2 << b;
            break;
        case 't':
            trace_file = optarg;
            break;
        case 'v':
            verbosity = 1;
            break;
        case 'h':
            printUsage(argv);
            exit(0);
        default:
            printUsage(argv);
            exit(1);
        }
    }

    /* Make sure that all required command line args were specified */
    if (s == 0 || E == 0 || b == 0 || trace_file == NULL) {
        printf("%s: Missing required command line argument\n", argv[0]);
        printUsage(argv);
        exit(1);
    }

    /* Compute S, E and B from command line args */
    S = (unsigned int) pow(2, s);
    B = (unsigned int) pow(2, b);
 
    /* Initialize cache */
    initCache();

#ifdef DEBUG_ON
    printf("DEBUG: S:%u E:%d B:%u trace:%s\n", S, E, B, trace_file);
    printf("DEBUG: set_index_mask: %llu\n", set_index_mask);
#endif
 
    replayTrace(trace_file);

    /* Free allocated memory */
    freeCache();

    /* Output the hit and miss statistics for the autograder */
    printSummary(hit_count, miss_count, eviction_count);
    return 0;
}
