//
// Created by BJGomes on 20/02/2021.
//

#ifndef SARS_COV2_MUTATION_MAIN_H
#define SARS_COV2_MUTATION_MAIN_H

#endif //SARS_COV2_MUTATION_MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <math.h>

#define MAX100    100
#define BUF_SIZE  4096
#define DICT_SIZE 11

typedef struct gene {
    char name[MAX100];
    char * seq;
}GENE;

typedef struct genome {
    char name[MAX100];
    GENE * genes;
    int n_genes;
    int size_genes;
    struct genome * pnext;
    struct genome * pprev;
}GENOME;

typedef struct genome_list {
    GENOME * phead;
    GENOME * ptail;
    long n_genomes;
}GENOME_LIST;

typedef struct gene_dict {
    char name[MAX100];
    char prot[MAX100];
    int rows;
    int cols;
    int padding;
}GENE_DICT;

typedef struct int_array {
    int * arr;
    int n;
    int size;
}INT_ARRAY;

typedef struct mutation {
    char genome_a[MAX100];
    char genome_b[MAX100];
    char gene[MAX100];
    INT_ARRAY seq_mutations;
}MUTATION;

typedef struct mutation_array {
    int n_mutations;
    int size_mutations;
    MUTATION * mutations;
}MUTATION_ARRAY;

GENE_DICT gd[] = {{"gene=ORF1ab", "ORF1ab", 115, 115, 7},
                  {"gene=S", "surface", 62, 62, 22},
                  {"gene=ORF3a", "ORF3a", 28, 30, 12},
                  {"gene=E", "envelope", 15, 16, 12},
                  {"gene=M", "membrane", 26, 27, 33},
                  {"gene=ORF6", "ORF6", 14, 14, 10},
                  {"gene=ORF7a", "ORF7a", 19, 20, 14},
                  {"gene=ORF7b", "ORF7b", 12, 12, 12},
                  {"gene=ORF8", "ORF8", 19, 20, 14},
                  {"gene=N", "nucleocapsid", 36, 36, 36},
                  {"gene=ORF10", "ORF10", 11, 11, 4}};

GENE * find_gene (GENOME * genome, char * gene_name);
void insert_int_array (INT_ARRAY * int_array, int element);
INT_ARRAY * gene_cmp (GENE g1, GENE g2);
void insert_mutation (MUTATION_ARRAY * mutation_array, char * genome_a, char * genome_b, char * gene, INT_ARRAY * gene_mut);
void genome_cmp (GENOME * genome, MUTATION_ARRAY * mutation_array);
void read_genomes (GENOME_LIST * gl, char * path);
void remove_white_spaces (char * str);
void insert_gene (GENOME * genome, GENE * gene);
char * find_protein_name (char * protein);
GENOME * find_genome (GENOME_LIST * gl, char * g_id);
void insert_genome (GENOME_LIST * gl, GENOME * g);
void print_genome (GENOME genome);
GENE_DICT * find_gene_dict (char * name);
int get_gene_padding (char * name);
GENE * create_gene (char * name, char * seq);
void parse_genome (GENOME_LIST * gl, char * cds);
long get_time_useconds(long * time_usec);
void save_mutation_array (MUTATION_ARRAY * mutation_array, char * path, int detail);
void free_mutations(MUTATION_ARRAY * mutation_array);
void handler();
ssize_t readn(int fd, void *ptr, size_t n);
ssize_t writen(int fd, const void *ptr, size_t n);
void save_mutation(int fd, char* genome_a, char * genome_b, char * gene, int n);
void save_to_file(MUTATION *aux, char * path);
void send_pipe_mutation_array (MUTATION_ARRAY * mutation_array, int *fds, int detail);
void send_socket_mutation_array (MUTATION_ARRAY * mutation_array, int * uds, int detail);