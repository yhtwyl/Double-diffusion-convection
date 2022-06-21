subroutine filename(Quantname,itau_, fname_)
  character(len=*), intent(in) :: Quantname
  character(len=20), intent(out) :: fname_
  integer itau_
  if (len(Quantname) .eq. 1) then
     write(fname_,20)Quantname,itau_
20   format(A1,I9.9,'.out')
  endif
  return
end subroutine filename
