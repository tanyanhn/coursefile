#ifndef SEQDEF_H_DEFINED
#define SEQDEF_H_DEFINED 1
void seqBegin(MPI_Comm);
void seqEnd(MPI_Comm);
void seqChangeOrder(MPI_Comm, int, int);
#endif
