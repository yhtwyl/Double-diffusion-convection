subroutine tridag(a__,b__,c__,d__,i__,x__,mp,where)
  integer mp,i,i__,k
  real, dimension(mp) :: a__,b__,c__,d__,x__,al,bta
  character(len=*) :: where 
  al(1) = b__(1)
!  print *,'The guilty is = ',where
  do I=2,i__
!     print *, 'tridag',i, al(i-1)
     if(abs(al(I-1)) .le. 1.0E-30) write(*,*)'I1=',I,'al(I-1)=',al(i-1)
     al(i)=b__(i)-(a__(i)*c__(i-1)/al(i-1))
  enddo
  bta(1)=d__(1)/b__(1)
  do i=2,i__
     if( abs(al(i)).le. 1.0E-30 ) write(*,*)'I2=',i,'Al(I)=',al(i)
     bta(i)=(d__(i)- a__(i)*bta(i-1))/al(i)
  enddo
  x__(i__)=bta(i__)
  k=i__
!  print *,'value of i__ = ',i__
  do i=1,i__-1
     k = k-1
     If(Abs(al(k)).Le. 1.E-30) Write(*,*)'K3=',k,'Al(K)=',al(k)
     x__(k)=bta(k)-(c__(k)*x__(k+1)/al(k))
  enddo
  return
end subroutine tridag
