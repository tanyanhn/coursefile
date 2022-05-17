void LAX_ddot_strided_blas(LA_Dvector_ddot_stride_binop *xy,
                           double *result)
{
    int    start, stride;
    double local_sum = 0.0;
    LA_Dvector *x, *y;

    /* Dereference the binary vector operands: */
    x = xy -> binop.x;  y = xy -> binop.y;

    /* Determine the stride based on type */
    start  = xy -> local_start; 
    stride = xy -> local_stride;

    /* Sum up my part */
    blas_ddot(&(x->v.m), &(x->v.data[start]), &stride,
	      &(y->v.data[start]), &stride);

    /* Get the sum of all parts */
    MPI_Allreduce(&local_sum, result, 1, MPI_DOUBLE, MPI_SUM,
		  xy->binop.comm);
    /* Return result */
    xy -> binop.error = 0;
}
