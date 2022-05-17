   program main
   use mpi_f08
   type(MPI_Comm)     :: comm
   type(MPI_Datatype) :: rtype
   type(MPI_Request)  :: req
   integer            :: myrank
   real, asynchronous :: a(100)   ! see text on asynchronous
!
   call mpi_init()
   call mpi_comm_dup(MPI_COMM_WORLD, comm)
   call mpi_type_contiguous(100, MPI_REAL, rtype)
   call mpi_type_commit(rtype)
   call mpi_comm_rank(comm, myrank)
   if (myrank .eq. 0) then
     call mpi_isend(a, 1, rtype, 1, 0, comm, req)
     call mpi_wait(req, MPI_STATUS_IGNORE)
   elseif (myrank .eq. 1) then
     call mpi_recv(a, 1, rtype, 0, 0, comm, MPI_STATUS_IGNORE)
   endif
   call mpi_comm_free(comm)
   call mpi_type_free(rtype)
   call mpi_finalize()
!
   end
