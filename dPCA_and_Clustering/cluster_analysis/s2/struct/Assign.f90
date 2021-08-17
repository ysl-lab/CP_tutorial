program assign
  implicit none
  integer, parameter :: ng = 50*50*50
  integer, dimension (ng) :: asn
  double precision, dimension(ng):: gx, gy, gz 
  double precision:: x, y, z, xsize_half, ysize_half, zsize_half
  integer :: ng_kept,ig, ig_max, cl, ik
  integer :: ios
  character(len=1000):: input, output


  gx = 0.d0
  gy = 0.d0
  gz = 0.d0
  asn = 0.d0

  call get_command_argument(1,input)
  call get_command_argument(2,output)

  xsize_half = 0.5*(3.0+3.0)/50.0
  ysize_half = 0.5*(3.0+2.0)/50.0
  zsize_half = 0.5*(3.0+2.0)/50.0


! Read the assignments for the density grid
  open(12,file='../GRID_ASSIGNATION',action='read')
  ig = 1
  do 
    read(12,*, iostat=ios) gx(ig), gy(ig), gz(ig), asn(ig)
    if (ios /= 0) exit
    ig = ig + 1
  end do
  close(12)
  ig_max = ig

  write(6,*) 'Number of assigned grids: ', ig_max

  open(12,file=input, action='read')
  open(13,file=output, action='write')

  ik = 1
  do
    read(12,*,iostat=ios) x, y, z
    if (ios /= 0) exit
    cl = 0
    do ig = 1, ig_max
      if (abs(x-gx(ig)) < xsize_half .and. abs(y-gy(ig)) < ysize_half .and. abs(z-gz(ig)) < zsize_half) then
        cl = asn(ig) 
      end if
    end do
    write(13,'(i8,f10.5,f10.5,f10.5,i8)') ik, x, y, z, cl
    ik = ik+1
    if (mod(ik,1000) .eq. 0) then
      write(6,*) ik, " frames calculated ..."
    end if
  end do
  close(12)
  close(13)


end program assign
