#include <sys/time.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "compute.h"
#include "math.h"

#include <stdio.h>

#define SQRT2   1.41421356237
#define WN      ( SQRT2 / ( SQRT2 + 1.0 ) )
#define WD      ( 1.0 / ( SQRT2 + 1.0 ) )

/**
 * @brief This function copies a two dimensional table to another.
 *
 * This function copies a contiguous subset of a table to a contiguous subset
 * of another table. Both tables need to point to allocated space, and be of
 * sufficient size for the operation to succeed. Both tables must point to
 * contiguous memory.
 *
 * @param in The table to copy from.
 * @param out The table to copy to.
 * @param N The number of rows to copy.
 * @param M The column size of both tables.
 * @param in_offset The row offset from where to copy.
 * @param out_offset The row offset where to write.
 */
void copy_table(const double *in, double *out, int N, int M, int in_offset, int out_offset) {
    in_offset = in_offset*M;
    out_offset = out_offset*M;

    for (int i = 0; i < N*M; i++)
        out[i+out_offset] = in[i+in_offset];
}

/**
 * @brief Copies a single row of data.
 *
 * This function copies M elements from in to out.
 *
 * @param in The memory address where to start copying.
 * @param out The memory address where to write.
 * @param M The number of elements to copy.
 */
void copy_line(const double *in, double *out, int M) {
    for (int i = 0; i < M; i++)
        out[i] = in[i];
}

/**
 * @brief This function updates the current results for the simulation.
 *
 * @param r The results struct where to write.
 * @param table_value A table element.
 * @param table_diff The difference between the table element and its previous iteration.
 */
inline void update_results(struct results *r, double table_value, double table_diff) {
    r->tmin =    ( r->tmin > table_value ) ? table_value : r->tmin;
    r->tmax =    ( r->tmax < table_value ) ? table_value : r->tmax;
    r->maxdiff = ( r->maxdiff < table_diff  ) ? table_diff : r->maxdiff;
    r->tavg += table_value;
}

double get_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL))
        return -1;
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

void do_compute(const struct parameters* p, struct results *r) {
    int M = p->M;
    int N = p->N;
    int iter = 0, max_iter = p->maxiter;
    double curr_diff;

    int k = p->period;

    // Create tables which include the unmutable edges
    double *table1 = malloc((N+2) * M * sizeof(double));
    double *table2 = malloc((N+2) * M * sizeof(double));
    double *prev_table, *curr_table, *temp;

    double begin, now;

    // Initialize main table
    copy_table(p->tinit,table1,N,M,0,1);
    // Initialize edges on both tables
    copy_line(p->tinit,table1,M);
    copy_line(p->tinit,table2,M);
    copy_line(&p->tinit[(N-1)*M],&table1[(N+1)*M],M);
    copy_line(&p->tinit[(N-1)*M],&table2[(N+1)*M],M);

    prev_table = table1;
    curr_table = table2;

    begin = get_time();

    int MN = M*N;
    int count = MN;
    int loopLimit = MN + M;
    
    while( count && iter < max_iter ){
        iter++;
        count = MN;
        int print_time = ! (iter % k);

        // This is the index of the first cell of the current line
        int curr_line = 0;
        // This is the next index, independent of line
        int next_index = 0;
        for(int i = M; i < loopLimit; i++){
            // This is the previous index, independent of line
            // Now it contains current index so we can use it for branches and avoid doing multiple modules
            int prev_index = i % M;
            if ( prev_index == (M-1) )
                next_index = -1;
            next_index++;

            if ( prev_index )
                prev_index--; 
            else {
                curr_line += M;     // If current index is 0, we are in a new line
                prev_index = M-1;   // And we set manually the previous index
            }
            double neighbors_conductivity = 1 - p->conductivity[i-M];
            double wn = neighbors_conductivity * WN / 4;
            double wd = neighbors_conductivity * WD / 4;
            curr_table[i] = prev_table[i] * p->conductivity[i-M] +
                wn * ( prev_table[i-M] + prev_table[i+M] +
                       prev_table[curr_line + prev_index] + prev_table[curr_line + next_index] ) +
                wd * ( prev_table[curr_line - M + prev_index] + prev_table[curr_line - M + next_index] +
                       prev_table[curr_line + M + prev_index] + prev_table[curr_line + M + next_index] );

            curr_diff = fabs(curr_table[i]-prev_table[i]);

            if ( curr_diff < p->threshold ) count--;
            // Only update if we need to
            if ( print_time ) {
                update_results(r, curr_table[i], curr_diff);
            }
        }
        // Is this our last iteration?
        int done = (! count) || (iter == max_iter);

        // This is called only once for our last iteration, and only if we didn't do it already.
        if ( done && !print_time ) {
            for(int i = M; i < loopLimit; i++) {
                curr_diff = fabs(curr_table[i]-prev_table[i]);
                update_results(r, curr_table[i], curr_diff);
            }
        }
        // Finalizes report struct, and if needed prints it
        if ( print_time || done ) {
            now = get_time();

            r->niter = iter;
            r->time = now - begin;

            r->tavg /= MN;

            if ( !done ) {
                report_results(p, r);
                // Reset report values for next report
                r->tmin += r->tmax;
                r->tmax = r->tmin - r->tmax;
                r->tmin -= r->tmax;
                r->maxdiff = 0.0;    

                r->tavg = 0;
            }
        }

        temp = prev_table;
        prev_table = curr_table;
        curr_table = temp;
        /*
        begin_picture(iter, M, N+2, -100, 100);
        for(int i = 0; i < (M*(N+2)); i++)
            draw_point( i % M, i / M, curr_table[i]);
        end_picture();
        */
    }

    free(table1);
    free(table2);
}
