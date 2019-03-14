

      subroutine odr_wrap_linear( n,m,np,nq,x,y, beta, ifix, &
                 fnm,infofnm, ux,uy, stdev_beta, cov, &
                 ressum, info)
      use real_precision
      use odrpack95

      implicit none

      !input
      integer, intent(in) :: n,m, np,nq
      real(kind=r8), intent(in) :: y(n,m), x(n,nq)
      real(kind=r8), intent(in) :: beta(np)
      real(kind=r8), intent(in) :: ux, uy
      character(len=500), intent(in) ::  fnm
      integer :: infofnm
      integer, intent(in) :: ifix(np)



      ! parameter
      integer :: info, k, l
      integer:: ldwe = 1, ld2we =1, ldwd = 1, ld2wd = 1
      integer::  lwork, liwork, sdi, niteri, vcvi, rvari, idfi

      real(kind=r8), allocatable, dimension(:) :: lowerbound
      real(kind=r8), allocatable, dimension(:) :: upperbound
      real(kind=r8) :: zero

      ! variable 
      integer :: niter
      integer :: mu
      integer, pointer, dimension(:) :: iwork
      integer, allocatable, target :: tiwork(:)


      real(kind=r8), allocatable, dimension(:,:,:) :: we, wd
      real(kind=r8), pointer, dimension(:) :: work
      real(kind=r8), allocatable, target :: twork(:)
      real(kind=r8) :: res_var

      !output
      real(kind=r8), intent(out) :: stdev_beta(np)
      real(kind=r8), intent(out) :: cov(np)
      real(kind=r8) :: covariance(np,np)
      real(kind=r8), intent(out) :: ressum

      character*256 ::  fnmrpt, fnmerr
      integer :: iosrpt, ioserr
      integer :: fidrpt, fiderr


      lwork = 18 + 13*np + np**2 + m + m**2 + &
               4*n*nq + 6*n*m + 2*n*nq*np + &
               2*n*nq*m + nq**2 + 5*nq + nq*(np+m) + n*nq*nq

      liwork = 20 + 2* np + nq *(np+m)               
      

      allocate(lowerbound(np))
      allocate(upperbound(np))

      allocate(wd(ldwd, ld2wd, m))
      allocate(we(ldwe, ld2we, nq))
      allocate(twork(lwork))
      allocate(tiwork(liwork))

      work => twork
      iwork => tiwork

!      ux = 1
!      uy = 0.001
!      write(*,*) 'ux', ux, 'uy', uy


      if (abs(ux)< 1e-10) then
              wd(1,1,1) = 1
      else
              wd(1,1,1) = -1./ux**2
      end if
      if (abs(uy)< 1e-10) then
              we(1,1,1) = 1
      else
              we(1,1,1) = -1/uy**2
      end if
 
      !if no filename suppress output
      if (infofnm .eq. 0) then
        fidrpt = 0
        infofnm = infofnm +1
        fiderr = 0
        infofnm = infofnm +2
      else
        !create rerport and error filenames
        fnmrpt=fnm(:infofnm)//'.rpt'
        fnmerr=fnm(:infofnm)//'.err'

        infofnm = 0 ! save information about file handles

        ! try to open the files, if that's not possible, send output to
        ! stdout
        fidrpt = 99
        open(unit = fidrpt, file=fnmrpt, iostat=iosrpt)
        if ( iosrpt > 0)  then
           fidrpt = 6
           infofnm = infofnm +1
        end if 

        fiderr = 88
        open(unit = fiderr, file=fnmerr, iostat=ioserr)
        if ( ioserr > 0) then
          fiderr = 6
          infofnm = infofnm +2
        end if 
      end if  


!      write (*,*) 'ifix', ifix !ifix = 0 fixed, ifix=1 unfixed

!     write(*,*) 'ux', ux, 'uy', uy
!     write(*,*) 'wd', wd, 'we', we

!     write(*,*) 'beta', beta
      
!     JOB 000000 auto forward fin. diff., 00010 - auto central fin. diff.; 00020 - user provided + check
! long output, print evrything possible on screen iprint 6616, no output iprint 0,
! iprint -1 => output to files lunrpt lunerr only
      ! iprint 2212 print long outputs to lunrpt and lunerr, don't print
      ! anything to stdout        
!      call odr(fitfunc, n, m, np, nq, beta, y,x, we = we ,wd = wd,   job=00020,  &
!            maxit = 1000, iprint = 0, info=info, iwork = iwork, work = work) 
      call odr(fitfunc, &
               n, m, np, nq, &
               beta, &
               y,x, we = we , wd = wd, &
               ifixb = ifix,&
               job = 00020, &
               maxit = 1000, iprint = 2212, &
               lunerr = fiderr, lunrpt = fidrpt, &
               work = work, iwork = iwork,&
               info=info)!, & 
          !     lower = lowerbound), &
!               upper = upperbound )

      if (infofnm == 0 .or. infofnm == 2 )       close(fidrpt)

      if (infofnm == 0 .or. infofnm == 1 )      close(fiderr)
      ! number of iterations
      niteri = nq*np + nq*m + np + 15
      niter = iwork(niteri)

      ! stdev of fitted parameters
      sdi = 2*n*m+2*n*nq+1
      do k =1, np
          stdev_beta(k) = work(sdi-1 +k)
      enddo

      ! covariance matrix befor rescaling
      vcvi = 2*n*m + 2*n*nq + np + 1
      do k = 1,np
         do l = 1,np
             covariance(k,l) = work(vcvi -1 + k + (l-1)*np)
         enddo
      enddo

      !residual variance
      rvari = 2*n*m + 2*n*nq + np + np*np+ 1
      res_var = work(rvari)
!      write(*,*) 'rvari', rvari
!      write(*,*) 'rres_var', res_var

      !degrees of freedom
      idfi = nq*np + nq*m + np + 6
      mu = iwork(idfi)
!      write(*,*) 'idfi', idfi
!      write(*,*) 'mu', mu
!      write(*,*) 'n', n

      !sum of squared residuals
      ressum = res_var*mu
!      write(*,*) 'ressum', ressum
 
      !rescale covariance matrix
      covariance = covariance*res_var

      !off-diagonal terms of covariance matrix
      if (np == 3) then
        cov(1) = covariance(2,3)
        cov(2) = covariance(1,3)
        cov(3) = covariance(1,2)
      else 
          cov(1) = covariance(1,2)
      end if

	!write(*,*) 'covariance'
!	write(*,*) covariance(1,1), covariance(1,2), covariance(1,3)
!	write(*,*) covariance(2,1), covariance(2,2), covariance(2,3)
!	write(*,*) covariance(3,1), covariance(3,2), covariance(3,3)

      !write(*,*) stdev_beta(1), stdev_beta(2),stdev_beta(3)
!      write(*,*) covariance(1,2)/stdev_beta(1)/stdev_beta(2), &
!                    covariance(1,3)/stdev_beta(1)/stdev_beta(3), &
!                    covariance(2,3)/stdev_beta(2)/stdev_beta(3)

      info = info + infofnm*1d6

      deallocate(twork)
      deallocate(tiwork)
      deallocate(wd)
      deallocate(we)

      contains

      real(kind=r8) function fun(x,a,b)

      implicit none
      real(kind=r8), intent(in) :: x,a,b

      fun = a*x+b

      end function fun




      real(kind=r8) function dfun_da(x,a,b)

      implicit none
      real(kind=r8), intent(in) :: x,a,b

      dfun_da= x

      end function dfun_da




      real(kind=r8) function dfun_db(x,a,b)

      implicit none
      real(kind=r8), intent(in) :: x,a,b

      dfun_db= 1

      end function dfun_db


      real(kind=r8) function dfun_ddeltai(x,a,b)

      implicit none
      real(kind=r8), intent(in) :: x,a,b

      dfun_ddeltai= a

      end function dfun_ddeltai

     

      subroutine fitfunc(n,m,np,nq,ldn, ldm, ldnp, beta, xplusd, ifixb, ifixx, &
      &      ldifx, ideval, f, fjacb, fjacd, istop)

      implicit none
      
      integer,intent(in) :: n,m, np, nq
      integer, intent(in) :: ldifx, ideval
      integer, intent(in) ::ldn, ldm, ldnp
      integer,intent(in) :: ifixb(np), ifixx(ldifx,m)
      
      real(kind=r8),intent(in) :: beta(np)
      real(kind=r8),intent(in) ::  xplusd(ldn,m)


      integer ::istop
      real(kind=r8) :: f(ldn, nq), fjacb(ldn, ldnp, nq), fjacd(ldn, ldm, nq) 

      real(kind=r8) :: x,a,b
      integer :: l,i,j
        
      if (istop .ne. 0) then
              return
      end if         

      if (mod(ideval,10) .ge. 1) then
              do l=1,nq
                 do i=1,n
                     f(i,l) = fun(xplusd(i,l), beta(1),beta(2))
                 end do
              end do    
      end if         
      if (mod(ideval/10,10) .ge. 1) then
              do l=1,nq
                 do i=1,n
                     fjacb(i,1,l) = dfun_da(xplusd(i,l), beta(1),beta(2))
                     fjacb(i,2,l) = dfun_db(xplusd(i,l), beta(1),beta(2))
                 end do
              end do    
      end if         
      if (mod(ideval/100,10) .ge. 1) then
              do l=1,nq
                 do j=1,m
                     do i=1,n
                         fjacd(i,j,l) = dfun_ddeltai(xplusd(i,j), beta(1),beta(2))
                     end do
                 end do
              end do    
      end if         

      end subroutine fitfunc

      end subroutine odr_wrap_linear

