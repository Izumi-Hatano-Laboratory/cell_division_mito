!* ---------------------------------------------------
!*       nucleus out program
!*       19.02.28
!*       Programmed by: Hiroaki Tanaka
!*       Modified by: A. H.
!*       Ver8: nuclear coloring is largely changed. Depth First Search
!*       Ver9: DFS is imposed to smarker and marker
!*       Ver10: Size check is implemented
!*       Ver11: Size check and thresholding, for nuclear. fileIO is simplified
!*       Ver12: discarded  
!*       Ver13: Cluster analysis, size and cluster file output
!*       Ver14: Larger mitochondria are sub-divided with new thresholds
!* ---------------------------------------------------
      implicit none

      integer, parameter :: width=1024, height=1024, depth=400
      integer :: color(width,height,depth), scolor(width,height,depth), color0(width,height,depth),scolor0(width,height,depth),ncolor(width,height,depth), nmark(width,height,depth),&
           smark(width,height,depth), mark(width,height,depth), mark0(width,height,depth), cluster(width,height,depth), mask(width,height,depth)
      integer :: nuc, depth2, size_thresh,nclu, nuc1, i, nsmall, x,y,z, j,k, count
      integer, allocatable :: size0(:), mito_location(:,:), size1(:), mito_location1(:,:),mit2clu(:), if_cut(:)
      integer :: cristae_thresh, scolor_thresh, color_thresh, nuc_thresh

      character(100) filehead, filename
!c#######################################################################

      filehead = './in_txt/' ! mitochondrial region learned
      call read_color(width, height, depth, color, filehead)

!      call write_bmp(depth, color, filehead)
!      stop

      filehead = './out_txt/'! mitochondrial contour learned
      call read_color(width, height, depth, scolor, filehead)

      color0  = color
      scolor0 =scolor
      cristae_thresh = 220
      scolor_thresh = 220
      color_thresh = 190
      nuc_thresh   = 240
      
      call make_scolor(width, height, depth, color, scolor, cristae_thresh, scolor_thresh, color_thresh)

! Subtraction, gausian filtering and extracting nuclear region (making ncolor, but here, 255 or 0)
      call nucleus_extraction(width, height, depth, nuc_thresh, scolor, ncolor, filename)
!c      filehead = './nucleus/'
!c      call write_color(width, height, depth, ncolor, filehead)
!c      call read_color(width, height, depth, ncolor, filehead)

! Make ncolor colorful (giving each mitocondrion a different color) 
      call nmarker(width, height, depth, ncolor, nmark, nuc, filename, 0)
      allocate(size0(nuc),mito_location(3,nuc))
      size_thresh = 100
      call size_check(nmark, width,height,depth, size_thresh, size0, 1, nuc, mito_location)

!      filehead = './nmark/'
!      call write_color(width, height, depth, nmark, filehead)
!c      call read_color(width, height, depth, nmark, filehead)

! Growing process to subtracted (mitochondria - edge)
      call marker(width, height, depth, scolor, nmark, smark)
!c      filehead = './smark/'
!c      call write_color(width, height, depth, smark, filehead)
!c      call read_color(width, height, depth, smark, filehead)

! Growing process 3D to mitochondrial region
      call marker(width, height, depth, color, smark, mark)
!c      filehead = './mark/'
!c      call write_color(width, height, depth, mark, filehead)
!c      call read_color(width, height, depth, mark, filehead)

      size_thresh = 2000
      call size_check(mark, width,height,depth, size_thresh, size0, 1, nuc, mito_location)
      call sort_by_size(mark, width,height,depth, size0, nuc, mito_location)

!       filehead = './mark/'
!       call write_color(width, height, depth, mark, filehead)
!c      call read_color(width, height, depth, mark, filehead)

!C# FOR Sub Division and Analysis
!      nuc = maxval(mark)
!      if(.not. allocated(size)) allocate(size0(nuc))
!      if(.not. allocated(size)) allocate(mito_location(3,nuc))
!      size_thresh = 500
!      call size_check(mark, width,height,depth, size_thresh, size, 1, nuc, mito_location)
!      call sort_by_size(mark, width,height,depth, size, nuc, mito_location)

      mark0  = mark ! save previous result
      color  = color0
      scolor = scolor0
      cristae_thresh = 200      ! lower, then more fragmented  # original: 220
      scolor_thresh = 220  ! higher,then more fragmented  # original: 220
!c      color_thresh = 190 ! Do not change, otherwize outer line will be changed.
      nuc_thresh = 247          ! 240 higher,then more fragmented  # original: 240
      write(*,*) 'thresholds cristae', cristae_thresh
      write(*,*) 'thresholds scolor ', scolor_thresh
      write(*,*) 'thresholds  color ',  color_thresh
      write(*,*) 'thresholds nuclear', nuc_thresh

      call make_scolor(width, height, depth, color, scolor, cristae_thresh, scolor_thresh, color_thresh)

! make mask 
! 1:large mito
! 0:small mito
      mask = 0
      nsmall = 0
      do i = 1, nuc
         if(size0(i) .le. 200000)then
            nsmall = nsmall + 1
         else
            x = mito_location(1,i)
            y = mito_location(2,i)
            z = mito_location(3,i)
            call DFSrenum(mark,mask,x,y,z,width,height,depth,i,1)
         end if
      end do
      write(*,*) 'over 200000 starts from',nsmall + 1
      write(*,*) 'Large mitochondria', nuc - nsmall


!color and scolor become only non-zero values when mito is large
      call mask_color(width, height, depth, color, scolor, mask)
      
! Subtraction, gausian filtering and extracting nuclear region (making ncolor, but here, 255 or 0)
      call nucleus_extraction(width, height, depth, nuc_thresh, scolor, ncolor, filename)
!      filehead = './nucleus/'
!      call write_color(width, height, depth, ncolor, filehead)
!      call write_bmp(depth, ncolor, filehead)
! c$$$c      call read_color(width, height, depth, ncolor, filehead)
! c$$$
! Make ncolor colorful (giving each mitocondrion a different color) 
      call nmarker(width, height, depth, ncolor, nmark, nuc1, filename, nsmall)
      write(*,*) 'Large mitochondria',nuc-nsmall,'was re devided into',nuc1-nsmall

      allocate(size1(nuc1),mito_location1(3,nuc1))
      size_thresh = 70
      size1(1:nsmall) = size0(1:nsmall)
      mito_location1(1:3,1:nsmall) = mito_location(1:3,1:nsmall) ! copy small sized mito
      call size_check(nmark, width,height,depth, size_thresh, size1, nsmall+1, nuc1, mito_location1)
! c$$$
!      filehead = './nmark/'
!      call write_bmp(depth, nmark, filehead)
! c$$$c      call read_color(width, height, depth, nmark, filehead)
! c$$$
! Growing process to subtracted (mitochondria - edge)
      call marker(width, height, depth, scolor, nmark, smark)
!c      filehead = './smark/'
!c      call write_color(width, height, depth, smark, filehead)
!c      call read_color(width, height, depth, smark, filehead)
      size_thresh = 1000
      call size_check(smark, width,height,depth, size_thresh, size1, nsmall+1, nuc1, mito_location1)

! Growing process 3D to mitochondrial region
      call marker(width, height, depth, color, smark, mark) 
      mask = 1 - mask ! invert the mask, 0:large, 1:small
      mark = mark + mask*mark0 ! superimpose large subdevided mito (mark) and already devided small mitochondria (mark0)
      call size_check(mark, width,height,depth, size_thresh, size1, 1, nuc1, mito_location1)

      filehead = './mark_subdivide/'
      call write_bmp(depth, mark, filehead)
      call write_color(width, height, depth, mark, filehead) ! 

!!$      nuc1 = 1106
!!$      allocate(size1(nuc1),mito_location1(3,nuc1))
!!$      call read_color(width, height, depth, mark, filehead)
!!$      call size_check(mark, width,height,depth, size_thresh, size1, 1, nuc1, mito_location1)
      write(*,*) 'max',maxval(mark),'nuc1',nuc1
      write(*,*) 'min',minval(mark) !for dbg
      count = 0
      ! clean pixcels
      do k = 1, depth
         do j = 1, height
            do i = 1, width
               if( mark(i, j, k) .gt. nuc1)then
                  write(*,*) 'error pixel', i,j,k,mark(i,j,k)
                  count = count + 1
                  mark(i,j,k) = 0
               end if
            end do
         end do
      end do
      write(*,*) 'error num', count

      ! check if the mito is cut by the volume edge, because cut mito should be excluded by the size analysis
      allocate(if_cut(0:nuc1))
      call check_if_cut(mark, width, height, depth, if_cut, nuc1)
      write(*,*) 'if_cut made'

!Cluster Analysis
      allocate(mit2clu(nuc1))
      call cluster_analysis(width,height,depth, cluster, nuc1,mark,mito_location1, size1, nclu,mit2clu, if_cut)

      filehead = './cluster/'
      call write_bmp(depth, cluster, filehead)
!      call write_color(width, height, depth, cluster, filehead)
!      call read_color(width, height, depth, mark, filehead)
      
      open(15,file='mito_size.out')
      write(15,'(10a10)') 'number','size','mit2clu','x','y','z'
      do i = 1, nuc1
        write(15,'(10i10)') i,size1(i),mit2clu(i),mito_location1(1:3,i),if_cut(i)
      end do
      close(15)

      stop
      end program

!################################################################################################################################
      subroutine check_if_cut(mark, width, height, depth, if_cut, nuc)
        implicit none
        integer :: width, height, depth, nuc, mark(width,height,depth), if_cut(0:nuc)
        integer :: i, j, k
        if_cut = 0
        do k = 1, depth
           do j = 1, height
              if_cut(mark(    1, j, k)) = 1
              if_cut(mark(width, j, k)) = 1
           end do
        end do
        do j = 1, height
           do i = 1, width
              if_cut(mark( i, j,    1)) = 1
              if_cut(mark( i, j,depth)) = 1
           end do
        end do
        do k = 1, depth
           do i = 1, width
              if_cut(mark( i,      1, k)) = 1
              if_cut(mark( i, height, k)) = 1
           end do
        end do
!!$      if_cut(mark(      1, 1:height, 1:depth)) = 1
!!$      if_cut(mark(  width, 1:height, 1:depth)) = 1
!!$      if_cut(mark(1:width,        1, 1:depth)) = 1
!!$      if_cut(mark(1:width,   height, 1:depth)) = 1
!!$      if_cut(mark(1:width, 1:height,       1)) = 1
!!$      if_cut(mark(1:width, 1:height,   depth)) = 1
        return
      end subroutine check_if_cut

!################################################################################################################################
      subroutine read_color(width, height, depth, color, filehead)

      implicit none
      integer :: width, height, depth, color(width,height,depth)
      integer :: i, j, k, num
      character(*) :: filehead
      character(100) :: filename
      character(3) str

      do k = 1, depth
         write(str,'(i3.3)')k
         write(*,*) str
         filename = trim(filehead)//str//'.txt'
         write(*,*)filename
         open(11,file=filename)
         do j = 1, height
            do i = 1, width
               read(11,*)num
               color(i,j,k) = num
            end do
         end do
         close(11)
      end do

      return
      end

!################################################################################################################################
      subroutine write_color(width, height, depth, color, filehead)

      implicit none

      integer :: width, height, depth, color(width,height,depth)
      integer :: i, j, k
      character(*) :: filehead
      character(100) :: filename
      character(3) str
      do k = 1, depth
         write(str,'(i3.3)')k
         filename = trim(filehead)//str//'.txt'
         write(*,*)filename
         open(13,file=filename)
         do j = 1, height
            do i = 1, width
               write(13,'(i8)')color(i,j,k)
            end do
         end do
         close(13)
      end do
      return
      end

!################################################################################################################################
      subroutine mask_color(width, height, depth, color,scolor,mask)

      implicit none
      integer :: width, height, depth, color(width,height,depth), scolor(width,height,depth), mask(width,height,depth)
      integer :: i, j, k, num
      character(100) :: filename
      character(3) str

      color = color*mask
      scolor = scolor*mask
!!$c$$$  do k = 1, depth
!!$c$$$         do j = 1, height
!!$c$$$            do i = 1, width
!!$c$$$               if(mask(i,j,k) .eq. 0)then
!!$c$$$                  color(i,j,k) = 0
!!$c$$$                  scolor(i,j,k) = 0
!!$c$$$               end if
!!$c$$$            end do
!!$c$$$         end do
!!$c$$$      end do
      return
      end
               
!################################################################################################################################
      subroutine make_scolor(width, height, depth, color,scolor, cristae_thresh, scolor_thresh, color_thresh)

      implicit none

      integer :: width, height, depth, color(width,height,depth), scolor(width,height,depth)
      integer :: i, j, k, num
      integer :: cristae_thresh, scolor_thresh, color_thresh
      character(100) :: filename1, filename2
      character(3) str

      do k = 1, depth
         do j = 1, height
            do i = 1, width
               ! contour space threshold: 220 (remove cristae membranes)
               if(scolor(i,j,k) .lt. cristae_thresh) scolor(i,j,k) = 0
               ! subtract
               scolor(i,j,k)= color(i,j,k) - scolor(i,j,k)
               ! scolor threshold: 220 ( minus numbers to be zero)
               if(scolor(i,j,k) .lt. scolor_thresh) scolor(i,j,k) = 0
               ! mitochondrial space threshold: 190 , and binarization
               if(color(i,j,k) .lt. color_thresh)then
                  color(i,j,k) = 0
               else
                  color(i,j,k) = 255
               end if
            end do
         end do
      end do

      write(*,*) 'subtraction done'
      return
      end


!################################################################################################################################
      subroutine nucleus_extraction(width, height, depth, nuc_thresh, scolor, ncolor, filename)

      implicit none

      integer :: width, height, depth, nuc_thresh, scolor(width,height,depth), ncolor(width,height,depth)
      integer,parameter :: xr=10, yr=10, zr=1
      real(8) :: gauss2(-xr:xr,-yr:yr),gauss1(-zr:zr), sum_color, sigma,sum_gauss

      integer, allocatable::sscolor(:,:,:)

      integer :: x, y, z, i, j, k, num, xi, yj, zk, n
      complex(16),parameter :: ci=cmplx(0d0,1d0)
      real(8),parameter :: pi=atan(1.d0)*4.d0
      character(*) :: filename
      character(3) str
      character(4) str2

      allocate(sscolor(-xr+1:width+xr, -yr+1:height+yr, -zr+1:depth+zr))

      ncolor = scolor
      ! enlarge image for filter
      sscolor(1:width,1:height,1:depth)=scolor(1:width,1:height,1:depth)
      do x = -xr+1, 0
         sscolor(x,1:height,1:depth) = scolor(1,1:height,1:depth)
      end do
      do x = width+1, width+xr
         sscolor(x,1:height,1:depth) = scolor(width,1:height,1:depth)
      end do
      do y = -yr+1,0
         sscolor(1-xr:width+xr,y,1:depth) = sscolor(1-xr:width+xr,1,1:depth)
      end do
      do y = height+1,height+yr
         sscolor(1-xr:width+xr,y,1:depth) = sscolor(1-xr:width+xr,height,1:depth)
      end do
      do z = -zr+1,0
         sscolor(1-xr:width+xr,1-yr:height+yr,z) = sscolor(1-xr:width+xr,1-yr:height+yr,1)
      end do
      do z = depth+1,depth+zr
         sscolor(1-xr:width+xr,1-yr:height+yr,z) = sscolor(1-xr:width+xr,1-yr:height+yr,depth)
      end do
         
      sigma = 1.9
! 2D Gauss function for xy direction
      sum_gauss = 0
      do y = -yr, yr
         do x = -xr, xr
            gauss2(x,y) = 1.D0/(2.D0*pi*sigma**2) * dexp(-(dble(x)**2+dble(y)**2)/(2.0*sigma**2))
!c     write(*,*)gauss(x,y,z)
            sum_gauss = sum_gauss + gauss2(x,y)
         end do
         write(*,*) 'gauss2',y,gauss2(-xr:xr,y)
      end do
      write(*,*)'sum_gauss2',sum_gauss

! 1D Gauss function for z direction
      sigma = 0.57
      sum_gauss = 0
      do z = -zr, zr
         gauss1(z) = 1/(sqrt(2*pi)*sigma) * dexp(-dble(z)**2/(2.0*sigma**2))
         sum_gauss = sum_gauss + gauss1(z)
      end do
      write(*,*)'sum_gauss1',sum_gauss
      write(*,*)'gauss1',gauss1(-zr:zr)

      do z = 1, depth
      do y = 1, height
      do x = 1, width
         sum_color = 0
!c            write(*,*)color(y,x,1:3)
         do k = -zr, zr
         do j = -yr, yr
         do i = -xr, xr
            sum_color = sum_color + sscolor(x+i,y+j,z+k) * gauss2(i,j) * gauss1(k)
         end do
         end do
         end do
! If you would like to adjust the threshold level (240), better to comment out  else sum=255 and visualize # previous 230
         if(sum_color .lt. nuc_thresh)then
            sum_color = 0
         else
            sum_color = 255
         end if
         ncolor(x,y,z) = sum_color
      end do
      end do
      end do

! 2D shlinc.  min pooling
!!$c$$$      xr = 3
!!$c$$$      yr = 3
!!$c$$$      sscolor = 255
!!$c$$$      sscolor(1:width,1:height,1:depth)=ncolor(1:width,1:height,1:depth)
!!$c$$$      do z = 1, depth
!!$c$$$      do y = 1, height
!!$c$$$      do x = 1, width
!!$c$$$         sum_color = sscolor(x,y,z)
!!$c$$$         do j = -yr, yr
!!$c$$$         do i = -xr, xr
!!$c$$$            sum_color = min(sum_color,sscolor(x+i,y+j,z+k))
!!$c$$$         end do
!!$c$$$         end do
!!$c$$$         end do
!!$c$$$         ncolor(x,y,z) = sum_color
!!$c$$$      end do
!!$c$$$      end do
!!$c$$$      end do
!!$c$$$      write(*,*)scolor(60,3,1)
!!$c$$$      write(*,*)ncolor(60,2,1)
      deallocate(sscolor)

      write(*,*) 'nucleus extraction done'
      return
      end


!################################################################################################################################
      subroutine nmarker(width, height, depth, ncolor, nmark, nuc,filename, nlast)

      implicit none

      integer :: width, height, depth, ncolor(width,height,depth), nmark(width,height,depth), nuc
      integer :: count, mx_nucleus, icheck(20), ilabel,isize, nlast
      integer :: x, y, z, i, j, k, n, m, im, jm, km, l,ll, i1, i0
      integer, allocatable ::renum(:), weight(:)
      double precision, allocatable :: zcenter(:)
      character(*) :: filename
      character(3) str

!Depth First Search
! Initially nmark contains parent pixel number
! 6 direction search
      nmark = -ncolor
      ilabel = nlast ! continues from the last count
      do z = 1, depth
      do y = 1, height
      do x = 1, width
         if(nmark(x,y,z) .lt. 0)then ! inner mitochondria and unlabelled is -255
            ilabel = ilabel+1
            isize = 0
            call DFSearch(nmark,x,y,z,width,height,depth,ilabel,isize)
            write(*,*) ilabel, isize
         end if
      end do
      end do
      end do
      nuc = ilabel
      
      return
      end

!################################################################################################################################
      RECURSIVE subroutine DFSearch(nmark,x,y,z,width,height,depth, ilabel,isize)
      implicit none
      integer :: width, height, depth, nmark(width,height,depth)
      integer :: ilabel,x,y,z, isize

      if(  x.le.0 .or. x.gt.width .or. &
           y.le.0 .or. y.gt.height .or. &
           z.le.0 .or. z.gt.depth)then ! out of boundary
         return
      else if(nmark(x,y,z) .ge. 0)then ! not a mitochondrial pixel, or already marked
         return
      else
         nmark(x,y,z) = ilabel
         isize = isize + 1
         call DFSearch(nmark,x+1,y,z,width, height, depth, ilabel,isize)
         call DFSearch(nmark,x-1,y,z,width, height, depth, ilabel,isize)
         call DFSearch(nmark,x,y+1,z,width, height, depth, ilabel,isize)
         call DFSearch(nmark,x,y-1,z,width, height, depth, ilabel,isize)
         call DFSearch(nmark,x,y,z+1,width, height, depth, ilabel,isize)
         call DFSearch(nmark,x,y,z-1,width, height, depth, ilabel,isize)
         return
      end if
      return
      end
      
!################################################################################################################################
      subroutine smarker(width, height, depth, scolor, nmark, smark)
      implicit none

      integer :: width, height, depth, scolor(width,height,depth), nmark(width,height,depth), smark(width,height,depth)

      integer :: x, y, z, i, j, k, n, m, ilabel, icount
      integer :: memory(1:width,1:height,1:depth)
      integer :: wsize = 1

      smark=nmark
      do n = 1, 50 ! maximum 50 pixels growing
         memory = 0
         icount = 0
         do z = 1, depth
         do y = 1, height
         do x = 1, width
            if(smark(x,y,z) .ne. 0 .and. memory(x,y,z).eq.0)then
               ilabel = smark(x,y,z)
               call DFSgrow2D(smark,scolor,memory,x,y,z, width,height,depth,ilabel,icount)
            end if
         end do
         end do
         end do
         smark = memory
         write(*,*)n,'growed pixel',icount
      end do

      return
      end

!################################################################################################################################
      subroutine marker(width, height, depth, color, smark, mark)
      implicit none
      integer :: width, height, depth, color(width,height,depth),smark(width,height,depth), mark(width,height,depth)
      integer :: x, y, z, i, j, k, n, m, ilabel, icount
      integer :: memory(1:width,1:height,1:depth)
      integer :: wsize = 1

      mark=smark
      do n = 1, 50 ! maximum 30 pixels growing
         memory = 0
         icount = 0
         if(mod(n,20) .eq. 1)then
            do z = 1, depth
            do y = 1, height
            do x = 1, width
               if(mark(x,y,z) .ne. 0 .and. memory(x,y,z).eq.0)then
                  ilabel = mark(x,y,z)
                  call DFSgrow3D(mark,color,memory,x,y,z, width,height,depth,ilabel,icount)
               end if
            end do
            end do
            end do
         else
            do z = 1, depth
            do y = 1, height
            do x = 1, width
               if(mark(x,y,z) .ne. 0 .and. memory(x,y,z).eq.0)then
                  ilabel = mark(x,y,z)
                  call DFSgrow2D(mark,color,memory,x,y,z,width,height,depth,ilabel,icount)
               end if
            end do
            end do
            end do
         end if
         mark = memory
         write(*,*)n,'growed pixel',icount
      end do
      return
      end

!################################################################################################################################
      RECURSIVE subroutine DFSgrow2D(mark,color,grown,x,y,z, w,h,d,ilabel,icount)
      implicit none
      integer :: w, h, d
      integer :: mark(w,h,d) ! input, colored
      integer :: color(w,h,d)! target, 0 or 255
      integer :: grown(w,h,d)! temporal, 1 pixel growed, and flag of visited
      integer :: ilabel,x,y,z,icount

      if(  x.le.0 .or. x.gt.w .or. &
           y.le.0 .or. y.gt.h .or. &
           z.le.0 .or. z.gt.d)then ! out of boundary
         return
      else if(grown(x,y,z) .ne. 0)then ! already visited
         return
      else if(mark(x,y,z) .eq. 0)then ! reached an unmarked pixel
         if(color(x,y,z) .ne. 0)then
            grown(x,y,z)=ilabel
            icount=icount+1
         end if
         return
      else if(mark(x,y,z) .eq. ilabel)then
        grown(x,y,z)=ilabel
        call DFSgrow2D(mark,color,grown,x+1,y,z,w,h,d,ilabel,icount)
        call DFSgrow2D(mark,color,grown,x-1,y,z,w,h,d,ilabel,icount)
        call DFSgrow2D(mark,color,grown,x,y+1,z,w,h,d,ilabel,icount)
        call DFSgrow2D(mark,color,grown,x,y-1,z,w,h,d,ilabel,icount)
        return
      else ! already marked
         return
      end if
      return
      end
      
!################################################################################################################################
      RECURSIVE subroutine DFSgrow3D(mark,color,grown,x,y,z, w,h,d,ilabel,icount)
      implicit none
      integer :: w, h, d
      integer :: mark(w,h,d) ! input, colored
      integer :: color(w,h,d)! target, 0 or 255
      integer :: grown(w,h,d)! temporal, 1 pixel growed, and flag of visited
      integer :: ilabel,x,y,z,icount

      if(  x.le.0 .or. x.gt.w .or. &
           y.le.0 .or. y.gt.h .or. &
           z.le.0 .or. z.gt.d)then ! out of boundary
         return
      else if(grown(x,y,z) .ne. 0)then ! already visited
         return
      else if(mark(x,y,z) .eq. 0)then ! reached an unmarked pixel
         if(color(x,y,z) .ne. 0)then
            grown(x,y,z)=ilabel
            icount=icount+1
         end if
         return
      else if(mark(x,y,z) .eq. ilabel)then
        grown(x,y,z)=ilabel
        call DFSgrow3D(mark,color,grown,x+1,y,z,w,h,d,ilabel,icount)
        call DFSgrow3D(mark,color,grown,x-1,y,z,w,h,d,ilabel,icount)
        call DFSgrow3D(mark,color,grown,x,y+1,z,w,h,d,ilabel,icount)
        call DFSgrow3D(mark,color,grown,x,y-1,z,w,h,d,ilabel,icount)
        call DFSgrow3D(mark,color,grown,x,y,z+1,w,h,d,ilabel,icount)
        call DFSgrow3D(mark,color,grown,x,y,z-1,w,h,d,ilabel,icount)
        return
      else ! already marked
         return
      end if
      return
      end
      
!################################################################################################################################
      subroutine size_check(mark, width,height,depth, size_thresh, size, istart, nuc, mito_location)
      implicit none

      ! istart: mito index is from istart to nuc, size check is done and sorted
      integer :: width, height, depth, mark(width,height,depth), istart, nuc, size(nuc), mito_location(3,nuc), size_thresh, temp(3,nuc)
      integer :: count, ilabel, ichk(width,height,depth),isize
      integer :: renum(nuc)
      integer :: x, y, z, i, j, k
      character(100) :: filename
      character(3) :: str

      ichk = 0
      do ilabel = istart, nuc
      do z = 1, depth
      do y = 1, height
      do x = 1, width
        if(mark(x,y,z) .eq. ilabel .and. ichk(x,y,z) .eq. 0)then ! find a pixel of the mitochondrion
           mito_location(1,ilabel)=x
           mito_location(2,ilabel)=y
           mito_location(3,ilabel)=z
           isize = 0
          call DFScount(mark,ichk,x,y,z,width,height,depth,ilabel,isize)
          write(*,*) ilabel, isize
          size(ilabel) = isize
        end if
      end do
      end do
      end do
      end do

      renum = 0
      count = istart-1
      do i = istart, nuc
         if(size(i) .ge. size_thresh)then
            count = count + 1
            renum(i) = count
         else
            isize=0
            x=mito_location(1,i)
            y=mito_location(2,i)
            z=mito_location(3,i)
            call DFSzero(mark,x,y,z,width,height,depth,i,isize)
            write(*,*) i,'nuc removed',isize,size(i)
         end if
      end do
      do z = 1, depth
      do y = 1, height
      do x = 1, width
         if(mark(x,y,z) .gt. 0)then
            mark(x,y,z) = renum(mark(x,y,z))
         end if
      end do
      end do
      end do
      do i = istart, nuc
         if(renum(i) .ne. 0)then
            mito_location(1:3,renum(i))=mito_location(1:3,i)
            size(renum(i))=size(i)
         end if
      end do
      nuc = count

      return
      end

      subroutine sort_by_size(mark, width,height,depth, size, nuc, mito_location)
      implicit none
      integer :: width, height, depth, mark(width,height,depth), nuc, size(nuc), mito_location(3,nuc), temp(3,nuc)
      integer :: renum(nuc),turn(nuc)
      integer :: x, y, z, i, j, k
!     SORT by SIZE
      renum = 0
      call heapsort2(nuc, size, turn)
      do i = 1, nuc
         renum(turn(i)) = i
      end do
      do z = 1, depth
      do y = 1, height
      do x = 1, width
         if(mark(x,y,z) .gt. 0)then
            mark(x,y,z) = renum(mark(x,y,z))
         end if
      end do
      end do
      end do
      temp = mito_location
      do i = 1, nuc
         mito_location(1:3,renum(i))=temp(1:3,i)
      end do

      x = mito_location(1,1)
      y = mito_location(2,1)
      z = mito_location(3,1)
      if(mark(x,y,z) .ne. 1)then
         write(*,*) 'error, heapsort failed 1'
      end if
      x = mito_location(1,100)
      y = mito_location(2,100)
      z = mito_location(3,100)
      if(mark(x,y,z) .ne. 100)then
         write(*,*) 'error, heapsort failed 100'
      end if

! CHECK ROUTINE
!!$c$$$      do i = 1, nuc
!!$c$$$         write(*,*) 'nuc',i,
!!$c$$$     $    mark(mito_location(1,i),mito_location(2,i),mito_location(3,i))
!!$c$$$     $        ,size(i)
!!$c$$$      end do
!!$c$$$      count = 0
!!$c$$$      do z = 1, depth
!!$c$$$      do y = 1, height
!!$c$$$      do x = 1, width
!!$c$$$        if(mark(x,y,z) .gt. 0 .and. ichk(x,y,z) .eq. 0)then ! find a error pixel
!!$c$$$           write(*,*) 'error pixel', x,y,z
!!$c$$$        end if
!!$c$$$        if(mark(x,y,z) .gt. count)then
!!$c$$$           count = mark(x,y,z)
!!$c$$$        end if
!!$c$$$      end do
!!$c$$$      end do
!!$c$$$      end do
!!$c$$$      write(*,*) 'count',count,'nuc',nuc, 'should be consistent'
!!$c$$$
!!$c$$$!  SIZE OUTPUT ~ 100000
!!$c$$$      ichk = 0
!!$c$$$      count = 0
!!$c$$$      do i = 1, nuc
!!$c$$$         if(size(i) .le. 100000)then
!!$c$$$            count = count + 1
!!$c$$$            x=mito_location(1,i)
!!$c$$$            y=mito_location(2,i)
!!$c$$$            z=mito_location(3,i)
!!$c$$$            call DFSrenum(mark,ichk,x,y,z,width,height,depth,i,count)
!!$c$$$            write(*,*) 'smaller than 100000', i, count, size(i)
!!$c$$$         end if
!!$c$$$      end do
!!$c$$$      do k = 1, depth
!!$c$$$         write(str,'(i3.3)')k
!!$c$$$         filename = 'size_0_100000/'//str//'.txt'
!!$c$$$         write(*,*)filename
!!$c$$$         open(15,file=filename)
!!$c$$$         do j = 1, height
!!$c$$$            do i = 1, width
!!$c$$$               write(15,'(i8)')ichk(i,j,k)
!!$c$$$            end do
!!$c$$$         end do
!!$c$$$         close(15)
!!$c$$$      end do
!!$c$$$!  SIZE OUTPUT ~ 200000
!!$c$$$      ichk = 0
!!$c$$$      count = 0
!!$c$$$      do i = 1, nuc
!!$c$$$         if(size(i) .gt. 100000 .and. size(i) .lt. 200000)then
!!$c$$$            count = count + 1
!!$c$$$            x=mito_location(1,i)
!!$c$$$            y=mito_location(2,i)
!!$c$$$            z=mito_location(3,i)
!!$c$$$            call DFSrenum(mark,ichk,x,y,z,width,height,depth,i,count)
!!$c$$$            write(*,*) '100000 < size <= 200000', i, count, size(i)
!!$c$$$         end if
!!$c$$$      end do
!!$c$$$      do k = 1, depth
!!$c$$$         write(str,'(i3.3)')k
!!$c$$$         filename = 'size_100000_200000/'//str//'.txt'
!!$c$$$         write(*,*)filename
!!$c$$$         open(15,file=filename)
!!$c$$$         do j = 1, height
!!$c$$$            do i = 1, width
!!$c$$$               write(15,'(i8)')ichk(i,j,k)
!!$c$$$            end do
!!$c$$$         end do
!!$c$$$         close(15)
!!$c$$$      end do
!!$c$$$!  SIZE OUTPUT 200000 ~ 
!!$c$$$      ichk = 0
!!$c$$$      count = 0
!!$c$$$      do i = 1, nuc
!!$c$$$         if(size(i) .gt. 200000)then
!!$c$$$            count = count + 1
!!$c$$$            x=mito_location(1,i)
!!$c$$$            y=mito_location(2,i)
!!$c$$$            z=mito_location(3,i)
!!$c$$$            call DFSrenum(mark,ichk,x,y,z,width,height,depth,i,count)
!!$c$$$            write(*,*) '100000 < size <= 200000', i, count, size(i)
!!$c$$$         end if
!!$c$$$      end do
!!$c$$$      do k = 1, depth
!!$c$$$         write(str,'(i3.3)')k
!!$c$$$         filename = 'size_200000_upper/'//str//'.txt'
!!$c$$$         write(*,*)filename
!!$c$$$         open(15,file=filename)
!!$c$$$         do j = 1, height
!!$c$$$            do i = 1, width
!!$c$$$               write(15,'(i8)')ichk(i,j,k)
!!$c$$$            end do
!!$c$$$         end do
!!$c$$$         close(15)
!!$c$$$      end do

      return
      end

!################################################################################################################################
      RECURSIVE subroutine DFScount(mark,ichk,x,y,z,w,h,d,ilabel,isize)
      implicit none
      integer :: w, h, d, mark(w,h,d)
      integer :: ilabel,x,y,z, isize,ichk(w,h,d)

      if(  x.le.0 .or. x.gt.w .or.&
           y.le.0 .or. y.gt.h .or.&
           z.le.0 .or. z.gt.d)then ! out of boundary
         return
      else if(ichk(x,y,z) .gt. 0)then ! already visited
         return
      else if(mark(x,y,z) .eq. ilabel)then ! not a mitochondrial pixel, or already marked
         ichk(x,y,z) = ilabel
         isize = isize + 1
         call DFScount(mark,ichk,x+1,y,z,w, h, d, ilabel,isize)
         call DFScount(mark,ichk,x-1,y,z,w, h, d, ilabel,isize)
         call DFScount(mark,ichk,x,y+1,z,w, h, d, ilabel,isize)
         call DFScount(mark,ichk,x,y-1,z,w, h, d, ilabel,isize)
         call DFScount(mark,ichk,x,y,z+1,w, h, d, ilabel,isize)
         call DFScount(mark,ichk,x,y,z-1,w, h, d, ilabel,isize)
         return
      end if
      return
      end
      
!################################################################################################################################
      RECURSIVE subroutine DFSzero(mark,x,y,z,w,h,d,ilabel,isize)
      implicit none
      integer :: w, h, d, mark(w,h,d)
      integer :: ilabel,x,y,z, isize

      if(  x.le.0 .or. x.gt.w .or. &
           y.le.0 .or. y.gt.h .or. &
           z.le.0 .or. z.gt.d)then ! out of boundary
         return
      else if(mark(x,y,z) .ne. ilabel)then ! already visited or out of mito boundary
         return
      else ! ilabel mitochondrial pixel
         mark(x,y,z) = 0
         isize = isize + 1
         call DFSzero(mark,x+1,y,z,w, h, d, ilabel,isize)
         call DFSzero(mark,x-1,y,z,w, h, d, ilabel,isize)
         call DFSzero(mark,x,y+1,z,w, h, d, ilabel,isize)
         call DFSzero(mark,x,y-1,z,w, h, d, ilabel,isize)
         call DFSzero(mark,x,y,z+1,w, h, d, ilabel,isize)
         call DFSzero(mark,x,y,z-1,w, h, d, ilabel,isize)
         return
      end if
      return
      end
      
!################################################################################################################################
      RECURSIVE subroutine DFSrenum(mark,ichk,x,y,z,w,h,d, ilabel,irenum)
      implicit none
      integer :: w, h, d, mark(w,h,d),ichk(w,h,d)
      integer :: ilabel,x,y,z, irenum

      if(  x.le.0 .or. x.gt.w .or. &
           y.le.0 .or. y.gt.h .or. &
           z.le.0 .or. z.gt.d)then ! out of boundary
         return
      else if(ichk(x,y,z) .gt. 0)then ! already visited
         return
      else if(mark(x,y,z) .eq. ilabel)then !if not visited the very mitochondrion
         ichk(x,y,z) = irenum
         mark(x,y,z) = 0
         call DFSrenum(mark,ichk,x+1,y,z,w, h, d, ilabel,irenum)
         call DFSrenum(mark,ichk,x-1,y,z,w, h, d, ilabel,irenum)
         call DFSrenum(mark,ichk,x,y+1,z,w, h, d, ilabel,irenum)
         call DFSrenum(mark,ichk,x,y-1,z,w, h, d, ilabel,irenum)
         call DFSrenum(mark,ichk,x,y,z+1,w, h, d, ilabel,irenum)
         call DFSrenum(mark,ichk,x,y,z-1,w, h, d, ilabel,irenum)
         return
      end if
      return
      end
      
!#######################################################################################################################################
      subroutine cluster_analysis(width,height,depth,cluster,nuc, mark,mit_loc, msize, nclu,mit2clu, if_cut)
      implicit none
      integer :: width, height, depth, nuc, mit_loc(3,nuc), mark(width,height,depth),msize(nuc)

      integer :: x,y,z, i,j,kmit, ilabel, nclu, icluster, mxnmit, if_cut(0:nuc)
      integer :: cluster(width,height,depth),isize, clu_size(nuc),mit2clu(nuc), clu_nmit(nuc)
      integer :: ave_all, ave_uncut, ave(nuc), sum_vol_all, sum_vol_uncut, tot_vol
      integer, allocatable :: clu2mit(:,:)

      cluster=-mark
      clu_size=0

      ilabel = 0
      do z = 1, depth
      do y = 1, height
      do x = 1, width
         if(cluster(x,y,z) .lt. 0)then ! inner mitochondria and unlabelled is minus value
            ilabel = ilabel+1
            isize = 0
            call DFSearch(cluster,x,y,z,width,height,depth,ilabel,isize)
            write(*,*) ilabel, isize
            clu_size(ilabel) = isize
         end if
      end do
      end do
      end do
      nclu = ilabel

      clu_nmit = 0
      do i = 1, nuc
         icluster = cluster(mit_loc(1,i),mit_loc(2,i),mit_loc(3,i))
         mit2clu(i) = icluster
         clu_nmit(icluster) = clu_nmit(icluster) + 1
      end do

      mxnmit = 0
      do i = 1, nclu
         mxnmit = max(mxnmit,clu_nmit(i))
      end do
      allocate(clu2mit(mxnmit,nclu))

      sum_vol_uncut = 0
      clu_nmit = 0
      do i = 1, nuc
         icluster = mit2clu(i)
         clu_nmit(icluster) = clu_nmit(icluster) + 1
         clu2mit(clu_nmit(icluster),icluster)=i
         if(if_cut(i) .eq. 0)then
            sum_vol_uncut = sum_vol_uncut + msize(i)
         end if
      end do

      sum_vol_all = 0
      ave = 0

      do i = 1, nclu
         sum_vol_all = sum_vol_all + clu_size(i)
         ave(i) = clu_size(i)/clu_nmit(i)
      end do
      ave_all = sum_vol_all/nuc
      ave_uncut = sum_vol_uncut / (nuc - sum(if_cut))
      tot_vol = width*height*depth

      write(*,*) 'mitochondria number', nuc
      write(*,*) 'mitochondria number uncut', nuc - sum(if_cut)
      write(*,*) 'cluster num', nclu
      write(*,*) 'averaged msize', ave_all
      write(*,*) 'averaged uncut msize', ave_uncut

      open(15,file='mito_size.out')
      write(15,'(10a10)') 'number','msize','mit2clu','x','y','z'
      do i = 1, nuc
        write(15,'(10i10)') i,msize(i),mit2clu(i),mit_loc(1:3,i)
      end do
      close(15)

      open(16,file='cluster.out')
      write(16,*) 'mitochondria number', nuc
      write(*,*) 'mitochondria number uncut', nuc - sum(if_cut)
      write(16,*) 'averaged msize', ave_all
      write(16,*) 'averaged uncut msize', ave_uncut
      write(16,*) 'total mit volume, fraction', sum_vol_all, sum_vol_all/tot_vol ! calculated by pixel number
      write(16,*) ''
      write(16,*) 'cluster analysis'
      write(16,*) 'cluster num', nclu
      write(16,*) 'mx mito per cluster', mxnmit
      write(16,'(10a10)') 'cluster,','size','nmito','ave_size'
      do i = 1, nclu
         write(16,'(10i10)') i,clu_size(i),clu_nmit(i),ave(i)
      end do
      close(16)
      
      return
      end

!################################################################################################################################
      subroutine marker_line(width, height, depth, mark, mark_line, filename)
      implicit none

      integer :: width, height, depth, nucleus, mark(width,height,depth), mark_line(width,height,depth)
      integer :: x, y, z, i, j, k, n, check, i_nuc
      character(*) :: filename
      character(3) str


      mark_line = mark
      check = 0
      do z = 1, depth
      do y = 2, height-1 !BORDERS will not be erased
      do x = 2, width-1
         i_nuc = mark(x,y,z) 
         if(i_nuc .ne. 0)then
            if( mark(x-1,y,z).eq.i_nuc .and. mark(x+1,y,z).eq.i_nuc.and.&
               mark(x,y-1,z).eq.i_nuc .and. mark(x,y+1,z).eq.i_nuc)then
               mark_line(x,y,z)=0
            end if
         end if
      end do
      end do
      end do

      do k = 1, depth
         write(str,'(i3.3)')k
         filename = 'mark_line'//str//'.txt'
         write(*,*)filename
         open(27,file=filename)
         do j = 1, height
            do i = 1, width
!!$c               if(mark_line(i,j,k) .eq. 0)then
!!$c                  write(27,'(i8)')0
!!$c               else
                  write(27,'(i8)') mark_line(i,j,k)
!!$c               end if
            end do
         end do
         close(27)
      end do
      write(*,*)'mark_line'

      return
      end



!################################################################################################################################
      subroutine out_particles(width, height, depth, mark, nucleus, filename)
      implicit none

      integer :: width, height, depth, mark(width,height,depth), nucleus
      integer :: x, y, z, i, j, k, n, check
      character(*) :: filename
      open(21,file=filename)
      n = 0
      check = 0
      do z = 1, depth
      do y = 1, height
      do x = 1, width
         if(mark(x,y,z) .ne. 0)then
            do i = 1, 3
               if(x+(i-2).ge.1 .and. x+(i-2).le.width)then
               if(mark(x+(i-2),y,z) .ne. mark(x,y,z))then
                  check = 1
               end if
               end if
            end do
            do j = 1, 3
               if(y+(j-2).ge.1 .and. y+(j-2).le.height)then
               if(mark(x,y+(j-2),z) .ne. mark(x,y,z))then
                  check = 1
               end if
               end if
            end do
            do k = 1, 3
               if(z+(k-2).ge.1 .and. z+(k-2).le.depth)then
               if(mark(x,y,z+(k-2)) .ne. mark(x,y,z))then
                  check = 1
               end if
               end if
            end do

            if(check .eq. 1)then
               write(21,'(i0,x,i0,x,i0,x,i0)')x, y, z, mark(x,y,z)
               n = n + 1
               check = 0
            end if
         end if
      end do
      end do
      end do
      close(21)
      write(*,*)n

      return
      end


!################################################################################################################################
      subroutine out_fld(width, height, depth, mark, nucleus, filename, filename2)
      implicit none

      integer :: width, height, depth, mark(width,height,depth), nucleus
      integer :: x, y, z, i, j, k, n, check
      character(*) :: filename, filename2

      open(19,file=filename2)
      n = 0
      check = 0
!c      do z = 1, depth
      do z = 1, 5
      do y = 1, height
      do x = 1, width
         if(mark(x,y,z) .ne. 0)then
            do i = 1, 3
               if(x+(i-2).ge.1 .and. x+(i-2).le.width)then
               if(mark(x+(i-2),y,z) .ne. mark(x,y,z))then
                  check = 1
               end if
               end if
            end do
            do j = 1, 3
               if(y+(j-2).ge.1 .and. y+(j-2).le.height)then
               if(mark(x,y+(j-2),z) .ne. mark(x,y,z))then
                  check = 1
               end if
               end if
            end do
            do k = 1, 3
               if(z+(k-2).ge.1 .and. z+(k-2).le.depth)then
               if(mark(x,y,z+(k-2)) .ne. mark(x,y,z))then
                  check = 1
               end if
               end if
            end do

            if(check .eq. 1)then
               write(19,'(i0,2x,i0,2x,i0,2x,i0)')mark(x,y,z), x, y, z
               n = n + 1
               check = 0
            end if
         end if
      end do
      end do
      end do
      close(19)
      write(*,*)n

      open(16,file=filename)
      write(16,'(a)')'# AVS field file '
      write(16,'(a)')'ndim = 1'
      write(16,'(a,i0)')'dim1 = ', n
      write(16,'(a)')'nspace = 3'
      write(16,'(a)')'veclen = 1'
      write(16,'(a)')'data = int'
      write(16,'(a)')'field = irregular'
      write(16,'(a)')'label = pressure'

      write(16,'(a,a)')'variable 1 file=./marker ', 'filetype=ascii offset=0 stride=4'
      write(16,'(a,a)')'coord 1 file=./marker ', 'filetype=ascii offset=1 stride=4'
      write(16,'(a,a)')'coord 2 file=./marker ', 'filetype=ascii offset=2 stride=4'
      write(16,'(a,a)')'coord 3 file=./marker ', 'filetype=ascii offset=3 stride=4'
      close(16)

      return
      end

!#########################################################
      subroutine heapsort2(n,array,turn)
        implicit none
        integer,intent(in)::n
        integer,intent(out)::turn(1:n) ! turn(1) = 5 means, the first member after sorted was 5th member in original array.
        integer,intent(inout)::array(1:n)
        
        integer::i,k,j,l,m
        integer::t
        
        if(n.le.0)then
           write(6,*)"Error, at heapsort"; stop
        endif
        if(n.eq.1)return

        do i=1,N
           turn(i)=i
        enddo

        l=n/2+1
        k=n
        do while(k.ne.1)
           if(l.gt.1)then
              l=l-1
              t=array(l)
              m=turn(l)
           else
              t=array(k)
              m=turn(k)
              array(k)=array(1)
              turn(k)=turn(1)
              k=k-1
              if(k.eq.1) then
                 array(1)=t
                 turn(1)=m
                 exit
              endif
           endif
           i=l
           j=l+l
           do while(j.le.k)
              if(j.lt.k)then
                 if(array(j).lt.array(j+1))j=j+1
              endif
              if (t.lt.array(j))then
                 array(i)=array(j)
                 turn(i)=turn(j)
                 i=j
                 j=j+j
              else
                 j=k+1
              endif
           enddo
           array(i)=t
           turn(i)=m
        enddo

        return
      end subroutine heapsort2

