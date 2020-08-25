module jacobi
  ! 参考：http://hooktail.org/computer/index.php?Jacobi%CB%A1
  implicit none
!
contains
!
  subroutine jacobi_calcEigen(eigenValues,eigenVecs,inputMat,N)
    real(8),parameter::EPS = 1.d-6
    ! input  !
    integer,intent(in)::N
    real(8),intent(in)::inputMat(N,N)
    ! output !
    real(8),intent(out)::eigenValues(N),eigenVecs(N,N)

    real(8)::aMat(N,N),a2Mat(N,N),app,apq,aqq,alpha,beta,gamma,s,c, &
         & pnMat(N,N),tempMat(N,N)
    real(8)::nowMax,beforeMax
    integer::i,maxRow,maxColumn,iTemp, turn(N)

    eigenVecs(:,:) = 0.d0
    do i = 1,N
       eigenVecs(i,i) = 1.d0
    end do
    aMat(:,:)      = inputMat(:,:)
    iTemp = 0
!    write(*,*)'in jacobi'
    do while(.true.)
       iTemp = iTemp + 1
       ! 非対角成分のうち、上側のもので絶対値最大のものを返す
!!$       write(*,'(2i8,3e14.6)') iTemp, 0, nowMax
       beforeMax = nowMax*10
       call setMaxRowAndMaxColumn(maxRow,maxColumn,N,aMat)
       nowMax = abs(aMat(maxRow,maxColumn))
       write(*,'(i8,3e14.6)') iTemp, nowMax, beforeMax*0.1D0, (beforeMax*0.1D0-nowMax)/nowMax
!!$       write(*,'(2i8,3e14.6)') iTemp, 1, nowMax
       if(nowMax < EPS)exit
       if(iTemp /= 1 .and. abs((beforeMax*0.1D0-nowMax)/nowMax) < EPS)exit
       if(iTemp > 4*N)exit

       ! cosθ,sinθの算出
       app = aMat(maxRow,maxRow)
       apq = aMat(maxRow,maxColumn)
       aqq = aMat(maxColumn,maxColumn)

       alpha = 0.5d0*(app-aqq)
       beta  = - apq
       gamma = abs(alpha)/sqrt(alpha**2.d0 + beta**2.d0)

       c = sqrt(0.5d0*(1.d0+gamma))
       s = sqrt(0.5d0*(1.d0-gamma))
       if( (alpha*beta) < 0.d0)s = - s

       ! A'の計算
       a2Mat(:,:) = aMat(:,:)
       do i = 1,N
          a2Mat(maxRow,i)    = c * aMat(maxRow,i) - s * aMat(maxColumn,i)
          a2Mat(maxColumn,i) = s * aMat(maxRow,i) + c * aMat(maxColumn,i)
          a2Mat(i,maxRow)    = c * aMat(i,maxRow) - s * aMat(i,maxColumn)
          a2Mat(i,maxColumn) = s * aMat(i,maxRow) + c * aMat(i,maxColumn)
       end do

       a2Mat(maxRow,maxRow)      = c*(c*app-s*apq) - s*(c*apq-s*aqq)
       a2Mat(maxRow,maxColumn)   = s*(c*app-s*apq) + c*(c*apq-s*aqq)
       a2Mat(maxColumn,maxRow)   = a2Mat(maxRow,maxColumn)
       a2Mat(maxColumn,maxColumn)= s*(s*app+c*apq) + c*(s*apq+c*aqq)
       aMat(:,:) = a2Mat(:,:)
       call checkSym(a2Mat,N)

       ! eigenVecs =P1*P2*....*Pn
       pnMat(:,:) = 0.d0
       do i = 1,N
          pnMat(i,i) = 1.d0
       end do
       pnMat(maxRow,maxRow)      = c
       pnMat(maxRow,maxColumn)   = s
       pnMat(maxColumn,maxRow)   = -s
       pnMat(maxColumn,maxColumn)= c
       eigenVecs = matmul(eigenVecs,pnMat)
!        write(*,*)'Pn',pnMat
!        write(*,'(A,3F8.3)')'eigen1',eigenVecs(1:3,1)
!        write(*,'(A,3F8.3)')'eigen2',eigenVecs(1:3,2)
!        write(*,'(A,3F8.3)')'eigen3',eigenVecs(1:3,3)
!!$       write(*,'(2i8,3e14.6)') iTemp, 2, nowMax
    end do

    do i = 1,N
       eigenValues(i) = aMat(i,i)
    end do
    call heapsort_real(N,eigenValues,turn)
    tempMat = eigenVecs
    do i = 1,N
       eigenVecs(1:N,i) = tempMat(1:N,turn(i))
    end do
!!$    write(*,*)
!!$    write(*,'(A,4F8.3)')'eigen1',eigenValues(1),eigenVecs(1:3,1)
!!$    write(*,'(A,4F8.3)')'eigen2',eigenValues(2),eigenVecs(1:3,2)
!!$    write(*,'(A,4F8.3)')'eigen3',eigenValues(3),eigenVecs(1:3,3)

    return
  end subroutine jacobi_calcEigen
!!!!!!!!!!
  subroutine setMaxRowAndMaxColumn(row,column,N,aMat)
    ! input  !
    integer,intent(in)::N
    real(8),intent(in)::aMat(N,N)
    ! output !
    integer,intent(out)::row,column
    
    integer::i,j
    real(8)::dtemp,max

    max = 0.d0
    do i = 2,N
    do j = 1,(i-1)
       dtemp = abs(aMat(j,i))
       if(dtemp > max)then
          max = dtemp
          column = i
          row    = j
       end if
    end do
    end do
    return
  end subroutine setMaxRowAndMaxColumn
!!!!!!!!!!!
  subroutine checkSym(aMat,N)
    ! input  !
    integer,intent(in)::N
    real(8),intent(in)::aMat(N,N)

    integer::i,j
    real(8),parameter::MIN = 1.d-8
    do j = 1,N
    do i = (j+1),N
       if(abs(aMat(i,j)-aMat(j,i))>MIN)then
          write(*,*)'NOT SYM!!! at',i,j
          write(*,*)aMat(i,j),aMat(j,i)
          stop
       end if
    end do
    end do
    return
  end subroutine checkSym
end module jacobi

!#########################################################
subroutine heapsort_real(n,array,turn)
        implicit none
        integer,intent(in)::n
        integer,intent(out)::turn(1:n) ! turn(1) = 5 means, the first member after sorted was 5th member in original array.
        real(8),intent(inout)::array(1:n)
        
        integer::i,k,j,l,m
        real(8)::t
        
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
                 if(array(j).gt.array(j+1))j=j+1
              endif
              if (t.gt.array(j))then
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
      end subroutine heapsort_real
