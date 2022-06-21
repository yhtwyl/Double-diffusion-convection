!
!  This file contains a routine for producing a decomposition of a 1-d array
!  when given a number of processors.  It may be used in "direct" product
!  decomposition.  The values returned assume a "global" domain in [1:n]
!
subroutine MPE_DECOMP1D( n, numprocs, myid__, s, e )
  integer n, numprocs, myid__, s, e
  integer nlocal
  integer deficit
  !
  nlocal  = n / numprocs
  s      = myid__ * nlocal + 1
  deficit = mod(n,numprocs)
  s      = s + min(myid__,deficit)
  if (myid__ .lt. deficit) then
     nlocal = nlocal + 1
  endif
  e = s + nlocal - 1
  if (e .gt. n .or. myid__ .eq. numprocs-1) e = n
  return
end subroutine MPE_DECOMP1D
      
