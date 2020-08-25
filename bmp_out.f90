! got from the web: https://qiita.com/cure_honey/items/795bd0e048ffeadc3e63
! however, several errors to run in this machine. to_bmp was not allowed to have arguments,
! thus we have to change the all nx and ny values when changed. ( here, 1024)
    module m_bmp
      implicit none
      type :: t_bmp_file_header
        sequence                      ! 14bytes
        character(len = 2) :: bfType = 'BM' !!integer(2) :: bfType = transfer('BM', 0_2, 1) ! BitMap
        integer(4) :: bfSize          ! file size in bytes
        integer(2) :: bfReserved1 = 0 ! always 0
        integer(2) :: bfReserved2 = 0 ! always 0
        integer(4) :: bfOffBits
      end type t_bmp_file_header
      ! 
      type :: t_bmp_info_header 
        sequence
        integer(4) :: biSize     = Z'28' ! size of bmp_info_header ; 40bytes 
        integer(4) :: biWidth
        integer(4) :: biHeight
        integer(2) :: biPlanes   = 1 ! always 1
        integer(2) :: biBitCount
        integer(4) :: biCompression = 0 !0:nocompression,1:8bitRLE,2:4bitRLE,3:bitfield
        integer(4) :: biSizeImage
        integer(4) :: biXPelsPerMeter = 3780 ! 96dpi
        integer(4) :: biYPelsPerMeter = 3780 ! 96dpi 
        integer(4) :: biClrUsed      = 0
        integer(4) :: biClrImportant = 0 
      end type t_bmp_info_header  
      !
      type :: t_rgb
        sequence
        character :: b, g, r  ! order is b g r 
      end type t_rgb 
      !
      type :: t_bmp!(nx, ny)
        integer :: nx = 1024, ny = 1024
!!$        integer, len:: nx, ny  
        type(t_rgb) :: rgb(1024, 1024)
      contains 
        procedure :: wr => wr_bmp
        procedure :: pr_bmp
!!$        generic :: write(formatted) => pr_bmp
      end type
    contains   
      subroutine wr_bmp(bmp, fn)
        class(t_bmp), intent(in) :: bmp
        character(len = *), intent(in) :: fn
        type(t_bmp_file_header) :: bmp_file_header
        type(t_bmp_info_header) :: bmp_info_header
        bmp_file_header%bfSize      = 14 + bmp_info_header%biSize + 0 + bmp%nx * bmp%ny * 3
        bmp_file_header%bfOffBits   = 14 + bmp_info_header%biSize
        bmp_info_header%biWidth     = bmp%nx       ! nx shouold be a multiple of 4
        bmp_info_header%biHeight    = bmp%ny
        bmp_info_header%biBitCount  = 24           ! color depth 24bits
        bmp_info_header%biSizeImage = bmp%nx * bmp%ny * 3
        open(9, file = fn, form = 'binary', status = 'unknown')
        write(9) bmp_file_header
        write(9) bmp_info_header
        write(9) bmp%rgb
        close(9)
        return
      end subroutine wr_bmp
 ! convert to t_RGB    
      pure elemental type(t_rgb) function to_rgb(ir, ig, ib)
        integer, intent(in) :: ir, ig, ib
        to_rgb = t_rgb(achar(ib), achar(ig), achar(ir))
      end function to_rgb  

      subroutine pr_bmp(dtv, unit, iotype, vlist, io, iomsg)
        class(t_bmp), intent(in) :: dtv
        integer, intent(in) :: unit
        character(len = *), intent(in) :: iotype
        integer, intent(in) :: vlist(:)
        integer, intent(out) :: io
        character(len = *), intent(in out) :: iomsg
        character(len = 30) :: fmt
        if (iotype == 'LISTDIRECTED') then
          write(unit, *, iostat = io) 'nx =', dtv%nx, ', ny =', dtv%ny
        end if    
      end subroutine pr_bmp
    end module m_bmp

    module m_mandel
      implicit none
      integer, parameter :: maxiter = 254
    contains
      pure elemental integer function mandel(c)
        complex, intent(in) :: c
        complex :: z
        z = (0.0, 0.0)
        do mandel = 0, maxiter
          z = z * z + c
          if (abs(z) > 2.0) exit
        end do    
      end function mandel 
    end module m_mandel

    subroutine write_bmp(nz,mark,filehead)
      use m_bmp
      use m_mandel 
      implicit none
      integer, parameter :: nx = 1024, ny = 1024
      integer :: nz
      integer :: mark(nx,ny,nz), color(nx,ny), k, iy
      character(80) :: filehead
      character(100) :: filename
      character(3) str
      type(t_bmp) :: bmp 
      do k = 1, nz
         write(str,'(i3.3)')k
         filename = trim(filehead)//str//'.bmp'
         write(*,*) filename
!!$      integer :: ix, iy, iter, mandelbrot(nx, ny)
!!$      real    :: x(nx), y(ny)
!!$      complex :: c(nx, ny)

         forall (iy = 1:ny) color(1:nx,iy) = min(mod(mark(1:nx,ny-iy+1,k),255)+1,mark(1:nx,ny-iy+1,k))
!!$         color = transpose(color)
         bmp%rgb = to_rgb(color, color, color)
!!$         bmp%rgb = to_rgb(mark(:,:,k) , mark(:,:,k) ,mark(:,:,k))
         call bmp%wr(filename)
      end do
!!$      print *, 'BMP size: ', bmp
!!$      stop 
    end subroutine write_bmp
