! got from the web: https://qiita.com/cure_honey/items/795bd0e048ffeadc3e63
! however, several errors to run in this machine. to_bmp was not allowed to have arguments,
! thus we have to change the all nx and ny values when changed. ( here, 1024)
    module m_tif
      implicit none
      type :: t_tif_file_header
        sequence                      ! 8bytes
        integer(1) :: bfOrder0   = 73  ! Intel Type little endian
        integer(1) :: bfOrder1   = 73  ! Intel Type little endian
        integer(2) :: bfVersion  = 42  ! fix
        integer(4) :: bfPointer = 8    ! Pointer to IFD(info) data, unit is byte
      end type t_tif_file_header
      ! 
      type :: t_tif_info_header 
        sequence
        integer(2) :: biSize     = 1 ! IFD entry number
        integer(2) :: biCode     = 0
        integer(2) :: biDataType   = 1 ! 1:BYTE(unsigned int), 2]ASCII, 3:SHORT, 4:LONG< 5,RATIONAL, 6:SBYTE, 7:UNDEFINED, 8:SSHORT, 9:SLONG, 10:SRATIONAL, 11:FLOAT, 12:DOUBLE
        integer(4) :: biCountField = 1024*1024
        integer(4) :: biDataPointer= 26 ! 8 header + 18 info 
        integer(4) :: biIfdPointer = 0 ! no next IFD
      end type t_tif_info_header  
      !
      type :: t_rgb
        sequence
        character :: g  ! gray scale
      end type t_rgb
      !
      type :: t_tif!(nx, ny)
        integer :: nx = 1024, ny = 1024
!!$        integer, len:: nx, ny  
        type(t_rgb) :: rgb(1024, 1024)
      contains 
        procedure :: wr => wr_tif
        procedure :: pr_tif
!!$        generic :: write(formatted) => pr_tif
      end type
    contains   
      subroutine wr_tif(tif, fn)
        class(t_tif), intent(in) :: tif
        character(len = *), intent(in) :: fn
        type(t_tif_file_header) :: tif_file_header
        type(t_tif_info_header) :: tif_info_header
        open(9, file = fn//'.tif', form = 'binary', status = 'unknown')
        write(9) tif_file_header
        write(9) tif_info_header
        write(9) tif%rgb
        close(9)
        return
      end subroutine wr_tif
 ! convert to t_RGB    
      pure elemental type(t_rgb) function to_rgb(ig)
        integer, intent(in) :: ig
        to_rgb = t_rgb(achar(ig))
      end function to_rgb  

      subroutine pr_tif(dtv, unit, iotype, vlist, io, iomsg)
        class(t_tif), intent(in) :: dtv
        integer, intent(in) :: unit
        character(len = *), intent(in) :: iotype
        integer, intent(in) :: vlist(:)
        integer, intent(out) :: io
        character(len = *), intent(in out) :: iomsg
        character(len = 30) :: fmt
        if (iotype == 'LISTDIRECTED') then
          write(unit, *, iostat = io) 'nx =', dtv%nx, ', ny =', dtv%ny
        end if    
      end subroutine pr_tif
    end module m_tif

    subroutine write_tif(nz,mark,filehead)
      use m_tif
      use m_mandel 
      implicit none
      integer, parameter :: nx = 1024, ny = 1024
      integer :: nz
      integer :: mark(nx,ny,nz), color(nx,ny), k
      character(80) :: filehead
      character(100) :: filename
      character(3) str
      type(t_tif) :: tif 
      do k = 1, nz
         write(str,'(i3.3)')k
         filename = trim(filehead)//str
         write(*,*) filename
!!$      integer :: ix, iy, iter, mandelbrot(nx, ny)
!!$      real    :: x(nx), y(ny)
!!$      complex :: c(nx, ny)

         color(1:nx,1:ny) = mod(mark(1:nx,1:ny,k),255)
         color = transpose(color)
         tif%rgb = to_rgb(color)
!!$         tif%rgb = to_rgb(mark(:,:,k))
         call tif%wr(filename)
      end do
!!$      print *, 'TIF size: ', tif
!!$      stop 
    end subroutine write_tif
