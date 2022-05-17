typedef struct
{
   LA_Grid_2d *grid; /* grid on which the distribution is based */
   LA_Mapping *row;  /* row mapping */
   LA_Mapping *col;  /* col mapping */
} LA_Distrib_2d;
