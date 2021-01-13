/*
 * Description:
 *     Solve the Max Subarray/Subsequence Problem using Parallel Algorithm
 *
 * Tested with:
 *     ldd (Debian GLIBC 2.28-10) 2.28
 *     gcc (Debian 8.3.0-6) 8.3.0
 *     OpenMP 4.5 (201511)
 *
 * Not copyrighted
 */

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdarg.h>
#include <limits.h>
#include <math.h>

// Parallel library
#include <omp.h>

#define LOG2(X) ((unsigned) (8 * sizeof (unsigned long long) - \
                                __builtin_clzll((X)) - 1))
#define POW2(X) (1 << (X))

// OpenMP function declaration
#ifdef _OPENMP
int omp_get_num_threads (void);
void omp_set_num_threads (int);
void omp_set_dynamic (int);
void omp_set_nested (int);
#endif

// enum _EMODE represents the reading direction of a struct tablo
typedef enum _EMODE
{
    PREFIX = 0,
    SUFFIX = 1
} EMODE;

// Array representation
struct tablo
{
    int    * tab;
    int    size;
};

// Program utils
#define DEBUG 0
void _debug (const char *, ...);

// Prefix/Suffix Sum of a struct tablo
void sum (EMODE, struct tablo *, struct tablo **);
void montee (EMODE, struct tablo *, struct tablo *);
void descente (struct tablo *,struct tablo *);
void final (struct tablo *, struct tablo *);

// Array utils
void generate_array (struct tablo **);
size_t new_size (int);
void print_array (struct tablo *);
struct tablo * init_tablo (size_t);
size_t file_size (FILE *);
void read_array (struct tablo **, char *);

// Max Prefix/Suffix of a struct tablo
void max_montee (EMODE, struct tablo *, struct tablo *);
void max_descente (struct tablo *, struct tablo *);
void max_final (struct tablo *,struct tablo *);
void max (EMODE, struct tablo *, struct tablo **);

// Max Subarray/Subsequence Problem
void find_max_subarray (struct tablo *, struct tablo *);

int
main (int  argc,
      char * argv [])
{
    int i;
    struct tablo * Q;
    struct tablo * PSUM;
    struct tablo * SSUM;
    struct tablo * SMAX;
    struct tablo * PMAX;
    struct tablo * Ms;
    struct tablo * Mp;
    struct tablo * M;

#ifdef _OPENMP
    omp_set_num_threads(100);
#endif

    if (argc == 1)
    {
        generate_array(&Q);
    } else if (argc == 2)
    {
        read_array(&Q, argv[1]);
    } else
    {
        fprintf(stderr, "[ERR] Too much arguments\n" \
                        "Usage: ./program [OPTIONAL]\n" \
                        "\n"
                        "[OPTIONAL]\n"
                            "\t<FILE>\tFilename containing an array\n");
        return EXIT_FAILURE;
    }

    Ms = init_tablo(Q->size);
    Mp = init_tablo(Q->size);
    M = init_tablo(Q->size);

    sum(PREFIX, Q, &PSUM);
    _debug("PSUM: "); print_array(PSUM);
    sum(SUFFIX, Q, &SSUM);

    _debug("SSUM: "); print_array(SSUM);

    max(SUFFIX, PSUM, &SMAX);
    max(PREFIX, SSUM, &PMAX);

    _debug("SMAX: "); print_array(SMAX);
    _debug("PMAX: "); print_array(PMAX);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 1; i < Q->size; i++)
    {
        Ms->tab[i] = PMAX->tab[i] - SSUM->tab[i] + Q->tab[i];
        _debug("\tMs->tab[%d] = %d - %d + %d\n",
            i,
            PMAX->tab[i],
            SSUM->tab[i],
            Q->tab[i]);

        Mp->tab[i] = SMAX->tab[i] - PSUM->tab[i] + Q->tab[i];
        _debug("\tMp->tab[%d] = %d - %d + %d\n",
            i,
            SMAX->tab[i],
            PSUM->tab[i],
            Q->tab[i]);

        _debug("M->tab[%d] = %d + %d - %d\n",
            i,
            Ms->tab[i],
            Mp->tab[i],
            Q->tab[i]);
        M->tab[i] = Ms->tab[i] + Mp->tab[i] - Q->tab[i];
    }

    print_array(M);

    find_max_subarray(Q, M);

    free(Q->tab);
    free(Q);
    free(PSUM->tab);
    free(PSUM);
    free(SSUM->tab);
    free(SSUM);
    free(PMAX->tab);
    free(PMAX);
    free(SMAX->tab);
    free(SMAX);
    free(Ms->tab);
    free(Ms);
    free(Mp->tab);
    free(Mp);
    free(M->tab);
    free(M);
    return EXIT_SUCCESS;
}

/*
 * Description:
 *     generate_array() will generate an array
 * Parameters
 *     arr     The "struct tablo *" destination
 */
void
generate_array (struct tablo ** arr)
{
    //size_t i;

    //construction d'un tableau pour tester
    *arr = calloc (sizeof(struct tablo), 1);
    // 524288 + 1 // 20
    (*arr)->size = (1 << 4) + 1;
    (*arr)->tab = calloc ((*arr)->size * sizeof(int), sizeof(int));
    (*arr)->tab[0] = 0; // useless
    /*
    for (i = 1; i < (*arr)->size; i++)
    {
        (*arr)->tab[i] = (unsigned) rand() % 50 - 25;
    }*/
    (*arr)->tab[0] = 0;
    (*arr)->tab[1] = 3;
    (*arr)->tab[2] = 2;
    (*arr)->tab[3] = -7;
    (*arr)->tab[4] = 11;
    (*arr)->tab[5] = 10;
    (*arr)->tab[6] = -6;
    (*arr)->tab[7] = 4;
    (*arr)->tab[8] = 9;
    (*arr)->tab[9] = -6;
    (*arr)->tab[10] = 1;
    (*arr)->tab[11] = -2;
    (*arr)->tab[12] = -3;
    (*arr)->tab[13] = 4;
    (*arr)->tab[14] = -3;
    (*arr)->tab[15] = 0;
    (*arr)->tab[16] = 2;
}

/*
 * Description:
 *     new_size() will generate the new size for the pseudo perfect binary-tree
 * Parameters
 *     size     The current size of the number list
 */
size_t
new_size (int size)
{
    unsigned int log2size;

    log2size = LOG2(size - 1);
    // - 1 car c'est la taille attendue
    // + 1 car on rajoute le 0 inutile
    return POW2(log2size + 1) - 1 + 1;
}

/*
 * Description:
 *     _debug() is a function which will print debug information on stdout if
 *     DEBUG macro is not 0
 * Parameters
 *     fmt     As printf
 *     ...     As printf
 */
void
_debug (const char * fmt,
        ...)
{
        if (DEBUG)
        {
                va_list arg;
                va_start(arg, fmt);
                vfprintf(stdout, fmt, arg);
                va_end(arg);
        }
}

/*
 * Description:
 *     print_array() is a debug function for printing a "struct tablo *" type
 * Parameters
 *     tmp     "struct tablo *" to display
 */
void
print_array (struct tablo * tmp)
{
    int size = tmp->size;
    int i;

    _debug("[%i", tmp->tab[1]);
    for (i = 2; i < size; i++)
    {
        _debug(", %i", tmp->tab[i]);
    }
    _debug("]\n");
}

/*
 * Description:
 *     init_tablo() will init a "struct tablo *" given an array size
 * Parameters
 *     size     The number of integers in the array
 */
struct tablo *
init_tablo (size_t size)
{
    struct tablo * tmp;

    tmp = (struct tablo *) calloc (sizeof(struct tablo), sizeof(int));
    tmp->size = size;
    tmp->tab = calloc (size * sizeof(int), sizeof(int));
    return tmp;
}

/*
 * Description:
 *     montee() is the first stage of the prefix/suffix sum, this function will
 *     compute the sum of 2 childs and redirect it to the parent
 * Parameters
 *     mode     PREFIX (0) or SUFFIX (1), this changes the reading direction
 *     source   "struct tablo *", where we need to do the sum
 *     dest     The binary tree, his size is source->size * 2
 */
void
montee (EMODE mode,
        struct tablo * source,
        struct tablo * destination)
{
    size_t l;
    size_t j;
    size_t m;
    size_t i;

    // remplissage du tableau destination de taille 2*m en
    // copiant les donnees du tableau source dans destination,
    // a la bonne position
    // on suppose que le malloc de destination a ete fait avant

    if (mode == PREFIX)
    {
        // From MID to RIGHT
        for (i = 0; i < source->size; i++)
        {
            destination->tab[(destination->size - source->size) + i] =
                source->tab[i];
        }
    } else if (mode == SUFFIX)
    {
        // from RIGHT to MID
        for (i = 0; i < source->size; i++)
        {
            destination->tab[(destination->size - source->size) + i] =
                source->tab[source->size - i];
        }
    }

    // Boucle de calcul pour la montee dans l'arbre/tableau
    m = LOG2(source->size);
    for (l = m; l > 0; l--)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (j = POW2(l - 1); j < POW2(l); j++)
        {
            /*_debug("destination->tab[%d] = %d = %d + %d\n",
                j,
                destination->tab[j],
                destination->tab[2 * j],
                destination->tab[(2 * j) + 1]);*/

            destination->tab[j] = destination->tab[2 * j]
                  + destination->tab[(2 * j) + 1];
        }
    }
}

/*
 * Description:
 *     descente() is the second stage of the prefix/suffix sum, this function
 *     will, if the index is even, take the parent's value, else it will
 *     sum the left node in montee() binary tree and the parent's value
 * Parameters
 *     a        montee() binary tree
 *     b        The destination binary tree
 */
void
descente (struct tablo * a,
          struct tablo * b)
{
    size_t l;
    size_t j;
    size_t m;

    b->tab[1] = 0;
    m = LOG2(b->size);
    for (l = 1; l < m + 1; l++)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (j = POW2(l); j < POW2(l + 1); j ++)
        {
            // if even(j)
            if (j % 2 == 0)
            {
                // Takes the parent's value
                b->tab[j] = b->tab[j / 2];
            }
            else
            {
                // sum the left node in montee() binary tree and the parent's
                // value
                b->tab[j] = b->tab[j / 2] + a->tab[j - 1];
            }
        }
    }
}

/*
 * Description:
 *     final() is the last stage of the prefix/suffix sum, this function
 *     will sum the montee() and descente() binary tree
 * Parameters
 *     a        montee() binary tree
 *     b        descente() binary tree
 */
void
final (struct tablo * a,
       struct tablo * b)
{
    size_t j;
    size_t m;

    m = LOG2(b->size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (j = POW2(m - 1); j < POW2(m); j ++)
    {
        b->tab[j] += a->tab[j];
    }
}

/*
 * Description:
 *     sum() will do the sum prefix/suffix of the array source and will output
 *     on the destination array
 * Parameters
 *     mode     PREFIX (0) or SUFFIX (1), this changes the reading direction
 *     source   "struct tablo *", where we need to do the sum
 *     dest     "struct tablo **", where we will redirect the output
 */
void
sum (EMODE mode,
     struct tablo * source,
     struct tablo ** dest)
{
    int i;
    struct tablo * b;
    struct tablo * a;
    struct tablo * res;

    a = init_tablo(new_size(source->size));
    montee(mode, source, a);
    //print_array(a);

    b = init_tablo(new_size(source->size));
    descente(a, b);
    //print_array(b);

    final(a, b);
    //print_array(b);
    res = init_tablo(source->size);

    if (mode == PREFIX)
    {
        for (i = b->size - res->size; i < b->size; i++)
        {
            res->tab[i - (b->size - res->size)] = b->tab[i];
        }
    } else if (mode == SUFFIX)
    {
        for (i = 0; i < res->size; i++)
        {
            res->tab[i] = b->tab[b->size - i];
        }
    }
    *dest = res;

    free(a->tab);
    free(a);
    free(b->tab);
    free(b);
}

/*
 * Description:
 *     max_montee() is the first stage of the max prefix/suffix, this function
 *     will compute the max of 2 childs and redirect it to the parent
 * Parameters
 *     source   "struct tablo *", where we need to do the max
 *     dest     The binary tree, his size is source->size * 2
 */
void
max_montee (EMODE mode,
            struct tablo * source,
            struct tablo * destination)
{
    size_t l;
    size_t j;
    size_t m;
    size_t i;

    // remplissage du tableau destination de taille 2*m en
    // copiant les donnees du tableau source dans destination,
    // a la bonne position
    // on suppose que le malloc de destination a ete fait avant

    if (mode == PREFIX)
    {
        // From MID to RIGHT
        for (i = 0; i < source->size; i++)
        {
            destination->tab[destination->size - source->size + i] =
                source->tab[i];
        }
    } else if (mode == SUFFIX)
    {
        // from RIGHT to MID
        for (i = 0; i < source->size; i++)
        {
            destination->tab[destination->size - source->size + i] =
                source->tab[source->size - i];
        }
    }

    // Boucle de calcul pour la montee dans l'arbre/tableau
    m = LOG2(source->size);
    for (l = m; l > 0; l--)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (j = POW2(l - 1); j < POW2(l); j++)
        {
            /*_debug("destination->tab[%d] = %d = fmax(%d, %d)\n",
                j,
                destination->tab[j],
                destination->tab[2 * j],
                destination->tab[(2 * j) + 1]);*/

            destination->tab[j] = fmax(destination->tab[2 * j],
                  destination->tab[(2 * j) + 1]);
        }
    }
}

/*
 * Description:
 *     max_descente() is the second stage of the max prefix/suffix, this
 *     function will, if the index is even, take the parent's value, else it
 *     will fmax the left node in max_montee() binary tree and the parent's
 *     value
 * Parameters
 *     a        max_montee() binary tree
 *     b        The destination binary tree
 */
void
max_descente (struct tablo * a,
              struct tablo * b)
{
    size_t l;
    size_t j;
    size_t m;

    b->tab[1] = INT_MIN;
    m = LOG2(a->size);
    for (l = 1; l < m + 1; l++)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (j = POW2(l); j < POW2(l + 1); j++)
        {
            // si pair(j)
            if (j % 2 == 0)
            {
                b->tab[j] = b->tab[j / 2];
            }
            else
            {
                b->tab[j] = fmax(b->tab[j / 2], a->tab[j - 1]);
            }
        }
    }
}

/*
 * Description:
 *     max_final() is the last stage of the max prefix/suffix, this function
 *     will fmax the max_montee() and max_descente() binary tree
 * Parameters
 *     a        max_montee() binary tree
 *     b        max_descente() binary tree
 */
void
max_final (struct tablo * a,
           struct tablo * b)
{
    size_t j;
    size_t m;

    m = LOG2(b->size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (j = POW2(m - 1); j < POW2(m); j ++)
    {
        b->tab[j] = fmax(b->tab[j], a->tab[j]);
    }
}

/*
 * Description:
 *     max() will do the max prefix/suffix of the array source and will output
 *     on the destination array
 * Parameters
 *     source   "struct tablo *", where we need to do the max
 *     dest     "struct tablo **", where we will redirect the output
 */
void
max (EMODE mode,
     struct tablo * source,
     struct tablo ** dest)
{
    int i;
    struct tablo * b;
    struct tablo * a;
    struct tablo * res;

    a = init_tablo(new_size(source->size));
    max_montee(mode, source, a);
    //print_array(a);

    b = init_tablo(new_size(source->size));
    max_descente(a, b);
    //print_array(b);

    max_final(a, b);
    //print_array(b);
    res = init_tablo(source->size);

    if (mode == PREFIX)
    {
        for (i = b->size - res->size; i < b->size; i++)
        {
            res->tab[i - (b->size - res->size)] = b->tab[i];
        }
    } else if (mode == SUFFIX)
    {
        for (i = 0; i < res->size; i++)
        {
            res->tab[i] = b->tab[b->size - i];
        }
    }
    *dest = res;

    free(a->tab);
    free(a);
    free(b->tab);
    free(b);
}

/*
 * Description:
 *     find_max_subarray() will resolve the max subarray problem, it will find
 *     the max in M and sum the indexes on the max subarray of M in Q
 * Parameters
 *     Q        The test/user array
 *     M        The array where the index of the max subarray is/are
 */
void
find_max_subarray (struct tablo * Q,
                   struct tablo * M)
{
    unsigned int imax_start, imax_length;
    int max_value;
    int MSQ;
    unsigned int i;

    max_value = M->tab[1];
    MSQ = Q->tab[1];
    imax_start = 1;
    imax_length = 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 2; i < M->size; i++)
    {
        _debug("{M[i]: %d, Q[i]: %d\n", M->tab[i], Q->tab[i]);
        if (M->tab[i] > max_value)
        {
            _debug("\t%d > %d\n", M->tab[i], max_value);
            imax_start = i;
            imax_length = 1;
            max_value = M->tab[i];
            MSQ = Q->tab[i];
        }
        else if (M->tab[i] == max_value)
        {
            imax_length++;
            MSQ += Q->tab[i];
            _debug("\t%d == %d, MSQ = %d\n", M->tab[i], max_value, MSQ);
        }
    }

    fprintf(stdout, "%d ", MSQ);
    for (i = imax_start; i < imax_start + imax_length - 1; i++)
    {
        fprintf(stdout, "%d ", Q->tab[i]);
    }
    fprintf(stdout, "%d\n", Q->tab[imax_start + imax_length - 1]);
}

/*
 * Description:
 *     file_size() will find the size of a file
 * Parameters
 *     f        The filestream
 */
size_t
file_size (FILE * f)
{
    size_t size;

    fseek(f, 0L, SEEK_END);
    size = ftell(f);
    fseek(f, 0L, SEEK_SET);
    return size;
}

/*
 * Description:
 *     read_array() will read an integer array in a file
 * Parameters
 *     Q        The destination "struct tablo *"
 *     filename The file name
 */
void
read_array (struct tablo ** Q,
            char * filename)
{
    FILE * file;
    size_t size;
    unsigned int i;
    int value_read;

    file = fopen(filename, "r");
    size = file_size(file);

    *Q = init_tablo(size / 2 + 1);
    fseek(file, 0L, SEEK_SET);

    (*Q)->tab[0] = 0;
    i = 1;
    while (fscanf(file, "%d", &value_read) == 1)
    {
        (*Q)->tab[i] = value_read;
        i++;
    }
    (*Q)->size = i;
    fclose(file);
}
