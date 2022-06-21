program hello
  use mpi
  integer rank, size, ierror, tag, strlen, status(MPI_STATUS_SIZE)
  character*(MPI_MAX_PROCESSOR_NAME) node
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
  call MPI_Get_processor_name(node, strlen, ierror)
  print*, node, rank, ': Hello world'
  call MPI_FINALIZE(ierror)
end program hello
