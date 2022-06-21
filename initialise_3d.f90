subroutine Initialise(m,n,p,x,xdu,xdv,xdw,xdt,xveu,xvwu,xweu,xwwu,xdeu,xdwu,xdev,xdwv,xdew,xdww,xdet,xdwt, &
     xtc,y,ydu,ydv,ydw,ydt,yunv,yusv,ywnv,ywsv,ydnu,ydsu,ydnv,ydsv,ydnw,ydsw,ydnt,ydst, &
     ytc,z,zdu,zdv,zdw,zdt,zutw,zubw,zvtw,zvbw,zdtu,zdbu,zdtv,zdbv,zdtw,zdbw,zdtt,zdbt,ztc, &
     nx,ny,nz,sx,ex,sy,ey,sz,ez,arxz,aryz,thigh,tlow,tinbot,cinbot,t_,c_,u_,v_,w_)

  integer m,n,p
  real, dimension(m) :: x,xdu,xdv,xdw,xdt,xveu,xvwu,xweu,xwwu,xdeu,xdwu,xdev,xdwv,xdew,xdww
  real, dimension(m) :: xdet,xdwt,xtl,xtr,xtc
  real, dimension(n) :: y,ydu,ydv,ydw,ydt,yunv,yusv,ywnv,ywsv,ydnu,ydsu,ydnv,ydsv,ydnw,ydsw
  real, dimension(n) :: ydnt,ydst,ytn,yts,ytc
  real, dimension(p) :: z,zdu,zdv,zdw,zdt,zutw,zubw,zvtw,zvbw,zdtu,zdbu,zdtv,zdbv,zdtw,zdbw
  real, dimension(p) :: zdtt,zdbt,ztt,ztb,ztc
  integer i,j,k, nx,ny,nz,sx,ex,sy,ey,sz,ez
  real, dimension(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1) :: t_,c_,u_,v_,w_
  real thigh,tlow,tinbot,cinbot,dx,dy,dz,arxz,aryz
  dx = arxz/(nx-1);dy = aryz/(ny-1);dz = 1.0/(nz-1);
  do i  = 1,nx
     x(i) = dx*(i-1);
  enddo
  
  do j = 1,ny
     y(j) = dy*(j-1);
  enddo
  
  do k = 1,nz
     z(k) = dz*(k-1);
  enddo

  
  do i=2,nx
     xtc(i)=(x(i)+x(i-1))/2.0
     xdt(i)= x(i)-x(i-1)
  enddo
  xtc(1)=x(1)
  xtc(nx+1)=x(nx)
  
  do i=2,nx-1
     xtr(i)=(x(i)+x(i+1))/2.0
  enddo
  xtr(nx)=x(nx)
  xtl(2)=x(1)

  do i=3,nx
     xtl(i)=(x(i-1)+x(i-2))/2.0
  enddo

  ytc(1)=y(1)

  do j=2,ny
     ytc(j)=(y(j) + y(j-1))*0.5
     ydt(j)= y(j)-y(j-1)
  enddo
  ytc(ny+1)=y(ny)
  
  do j=2,ny-1
     ytn(j)=(y(j)+y(j+1))*0.5
  enddo
  
  ytn(ny) = y(ny)
  yts(2) = y(1)

  do j=3,ny
     yts(j)=(y(j-1)+y(j-2))/2.0
  enddo

  do i=2,nx
     xdet(i)=xtr(i)-xtc(i)
     xdwt(i)=xtc(i)-xtl(i)
  enddo

  do j=2,ny
     ydnt(j)=ytn(j)-ytc(j)
     ydst(j)=ytc(j)-yts(j)
  enddo

  do i=2,nx-1
     xdeu(i)=x(i+1)-x(i)
     xdwu(i)=x(i)-x(i-1)
     xdu(i)=xtc(i+1)-xtc(i)
  enddo

  do j=2,ny
     ydnu(j)=ytc(j+1)-ytc(j)
     ydsu(j)=ytc(j)-ytc(j-1)
     ydu(j)=y(j)-y(j-1)
  enddo

  do i=2,nx-1
     xveu(i)=xtc(i+1)-x(i)
     xvwu(i)=x(i)-xtc(i)
     xweu(i)=xtc(i+1)-x(i)
     xwwu(i)=x(i) - xtc(i)
  enddo

  do i=2,nx
     xdev(i)=xtc(i+1)-xtc(i)
     xdwv(i)=xtc(i)-xtc(i-1)
     xdew(i)=xtc(i+1)-xtc(i)
     xdww(i)=xtc(i)-xtc(i-1)
     xdv(i)=x(i)-x(i-1)
     xdw(i)=x(i)-x(i-1)
  enddo

  do j=2,ny-1
     ydnv(j)=y(j+1)-y(j)
     ydsv(j)=y(j)-y(j-1)
     ydnw(j)=ytc(j+1)-ytc(j)
     ydsw(j)=ytc(j) - ytc(j-1)
     ydv(j)=ytc(j+1)-ytc(j)
     ydw(j)=y(j) - y(j-1)
  enddo
  ydw(ny) = y(ny) - y(ny-1)
  ydnw(ny) = ytc(ny+1)-ytc(ny)
  ydsw(ny) = ytc(ny) - ytc(ny-1)
  
  do j=2,ny-1
     yunv(j)=ytc(j+1)-y(j)
     yusv(j)=y(j)-ytc(j)
     ywnv(j)=ytc(j+1) - y(j)
     ywsv(j)=y(j) - ytc(j)
  enddo
  
  
  ztc(1)=z(1)
  do k=2,nz
     ztc(k)=(z(k)+z(k-1))/2.
     zdt(k)= z(k)-z(k-1)
  enddo
  ztc(nz+1)=z(nz)
  do k=1,nz-1
     ztt(k)=(z(k)+z(k+1))/2.
  enddo
  ztt(nz)=z(nz)
  ztb(2)=z(1)
  do k=3,nz
     ztb(k)=(z(k-1)+z(k-2))/2.
  enddo

  do k=1,nz
     zdtt(k)=ztt(k)-ztc(k)
     zdbt(k)=ztc(k)-ztb(k)
  enddo
  
  do k=2,nz
     zdtu(k)=ztc(k+1)-ztc(k)
     zdbu(k)=ztc(k)-ztc(k-1)
     zdu(k)=z(k)-z(k-1)
  enddo
  
  do k=2,nz-1
     zdtv(k) = ztc(k+1) - ztc(k)
     zdbv(k) = ztc(k) - ztc(k-1)
     zdtw(k) = z(k+1) - z(k)
     zdbw(k) = z(k) - z(k-1)
     zubw(k) = z(k) - ztc(k)
     zutw(k) = ztc(k+1) - z(k)
     zvbw(k) = z(k) - ztc(k)
     zvtw(k) = ztc(k+1) - z(k)
     zdw(k) = ztc(k+1) - ztc(k)
     zdv(k) = z(k) - z(k-1)
  enddo
  zdv(nz) = z(nz) - z(nz-1)
  zdtv(nz) = ztc(nz+1) - ztc(nz)
  zdbv(nz) = ztc(nz) - ztc(nz-1)

  do k = sz,ez
     do j = sy-1,ey+1 
       do i = sx-1,ex+1
           if(ztc(k) .lt. 0.5) then
              t_(i,j,k) = amin1(( 1 -  ztc(k) )*tinbot,thigh)
              c_(i,j,k) = amin1(( 1 -  ztc(k) )*cinbot,thigh)
              !           u_(i,j) = 0.0
           else
              t_(i,j,k)= tlow
              c_(i,j,k)= tlow
              !           u_(i,j) = tlow
              !           if(xtc(i) .lt. 0.5*xtc(i)) u_(i,j)= tlow;
           endif
           !        t_(i,j) = j*1000 + i;
        enddo
     enddo
  enddo

  return
end subroutine Initialise
