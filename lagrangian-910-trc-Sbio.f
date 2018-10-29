      program floats
      use mod_floats ! HYCOM post-processing synthetic float interface
      use mod_za     ! HYCOM array I/O interface
c
      implicit none
c
c --- post-processing lagrangian floats from a sequence of HYCOM 2.0
c --- archive files.
c
cmtb- modified from original to read from lagrangian.input file
      character*120    flnm1,flnm2
      logical          lexist,lfatal,mtb
      integer          iexpt,yrflag,kdmin
      integer          iinc,it,itr,jinc,jt,nntr,ntracr
      integer          iorig,jorig,idms,jdms,first,wbin
      integer          fbiotyp,addtra,tottr,fbiodims
      real             depthi(0:1,0:1,0:99),duk,dukm1,dvk,dvkm1
      real             fwbw,irad,idens,ikvis
      real             bioZmax,bioAr,bioIk,bioKn,groMax,Smort
      double precision time_ai,time_ao,time1,time2
      real, allocatable :: work1(:,:)

c      common/offln1/   idms,jdms,iorig,jorig
      common  /inerpars/  irad,idens,ikvis
      common  /ftrcpars/  fbiodims
      common  /sbiopars/  fbiotyp,ntracr,addtra,tottr,
     &                    bioZmax,bioAr,bioIk,bioKn,groMax,Smort
c
      call xcspmd
      call zaiost
c
      lp=6
      onecm=0.01*onem
      tencm=0.10*onem
c
      pi = 4.0*atan(1.0)
      radian=pi/180.0
      mtb=.FALSE.
c
c --- 'flnm_fii' = name of input floats input file
c --- 'flnm_fio' = name of output floats input file for next time segment 
c --- 'flnm_fo'  = name of floats output file containing path positions
c --- 'flnm_ai'  = name of first input archive file
c --- 'flnm_ao'  = name of final input archive file
c --- 'iexpt '   = experiment number x10  (000=from archive file)
c --- 'yrflag'   = days in year flag (0=360J16,1=366J16,2=366J01,3=actual)
c --- 'archfq'   = no. of hours between archive input
c --- 'archin'   = no. of hours between time-interpolated archive input
c --- 'idm   '   = longitudinal array size
c --- 'jdm   '   = latitudinal  array size
c --- 'kdm   '   = number of layers
c
      open(unit=98,file='lagrangian.input')
      read (98,'(a)') flnm_fii
      write (lp,'(2a)')
     & 'input floats input file:  ',flnm_fii(1:len_trim(flnm_fii))
      call flush(lp)
      read (98,'(a)') flnm_fio
      write (lp,'(2a)')
     & 'output floats input file: ',flnm_fio(1:len_trim(flnm_fio))
      call flush(lp)
      read (98,'(a)') flnm_fo
      write (lp,'(2a)')
     & 'floats output file:       ',flnm_fo(1:len_trim(flnm_fo))
      call flush(lp)
      read (98,'(a)') flnm_ai
      write (lp,'(2a)')
     & 'archive input file:       ',flnm_ai(1:len_trim(flnm_ai))
      call flush(lp)
      read (98,'(a)') flnm_ao
      write (lp,'(2a)')
     & 'archive output file:      ',flnm_ao(1:len_trim(flnm_ao))
      call flush(lp)

cmtb

c
cmtb      open(unit=98,file='lagrangian.input')
      call blkini(iexpt, 'iexpt ')
      call blkini(yrflag,'yrflag')
      call blkinr(archfq,
     &           'diagfq','("blkinr: ",a6," =",f11.4," days")')
      call blkinr(archin,
     &           'archin','("blkinr: ",a6," =",f11.4," hours")')

c      archfq=archfq/24.0
      archin=archin/24.0
      nintrp=0
      write (lp,*)
     & 'archfq',archfq,' archin ',archin
      call flush(lp)
      if (archfq.ne.archin) then
        nintrp=nint(archfq/archin)
        write (lp,*) 'nintrp',nintrp
        if (archfq/archin.ne.float(nintrp)) then
          stop 'archfq must be an integer multiple of archin'
        endif
      endif

      call blkini(ii,    'idm   ')
      call blkini(jj,    'jdm   ')
      call blkini(kk,    'kdm    ')
      if     (ii.ne.idm .or. jj.ne.jdm) then
        write(lp,*)
        write(lp,*) 'error - wrong idm or jdm (should be:',idm,jdm,')'
        write(lp,*)
        call flush(lp)
        stop 
      endif
      call blkinr(thbase,
     &           'thbase','("blkinr: ",a6," =",f11.4," days")')
c
c --- 'fltdbg' = FLOATS: float number to be debugged
c --- 'intpfl' = FLOATS: horiz. interp. (0=2nd order+bilinear; 1=bilinear)
c --- 'iturbv' = FLOATS: add horiz. turb. advection velocity (0=no; 1=yes)
c --- 'tbvar'  = FLOATS: horizontal turbulent velocity variance scale
c --- 'tdecri' = FLOATS: inverse decorrelation time scale
c 
      call blkini(nfl_debug,'fltdbg')
      call blkini(intpfl,'intpfl')
      call blkini(iturbv,'iturbv')
      call blkinr(tbvar ,'tbvar ','(a6," =",f10.8," m**2/s**2")')
      call blkinr(tdecri,'tdecri','(a6," =",f10.4," 1/day")')
      call blkini(iorig, 'iorig ')
      call blkini(jorig, 'jorig ')
      call blkini(idms,  'idms  ')
      call blkini(jdms,  'jdms  ')
      call blkinr(fwbw,  'fwbw  ','(a6," =",f10.4,"timestep")')
      call blkini(first,  'firs  ')
      call blkini(wbin,   'wbin  ')
      turbvel = iturbv.ge.1
      call blkinl(trcout, 'trcout')
c      trcout = .false.
      call blkini(ntracr,   'ntracr')
      call blkinr(irad,   'irad   ','(a6," =",f10.4,"m")')
      call blkinr(idens,  'idens  ','(a6," =",f10.4,"g/m**2")')
      call blkinr(ikvis,  'ikvis  ','(a6," =",f10.4,"g/m**2")')
      call blkini(fbiotyp, 'fbiotyp')
      call blkini(addtra,  'addtra ')


      allocate( work1(ii,jj))
c
c --- determine size of float array
c --- fbiotyp: selects which bio code to run
c ---          1= MBrooks Sargassum v. 2.0
c ---          2= MBrooks Sargassum v. 3.0
      if(ntracr.gt.0) then
        tottr=ntracr+addtra
      else
        tottr=0
      endif
      if(fbiotyp.eq.0) then
        fbiodim=10+tottr
      elseif(fbiotyp.eq.1) then
        fbiodim=10+tottr+4 !sets maximum float bio state variables
      elseif(fbiotyp.eq.2) then
        fbiodim=10+tottr+7
      else
        write(lp,*)
        write(lp,*) 'error - illegal biotyp'
        write(lp,*)
        call flush(lp)
        stop '(biotyp)'
      endif
      !=standard float vars (=10) + tottr + fbiovars

cmtb
c
c --- geopar array allocation
c
      call floats_galloc
c
c --- land masks.
c
      call geopar
c --- adjust geopar arrays to new subregion framework
c
      ii=idms
      jj=jdms
      lperiod=.false.
      if(ii.le.idm) then
         do j=1,jj+3
           do i=1,ii+3
             work1(i,j)=depths(iorig+i-1,jorig+j-1)
             depths(i,j)=work1(i,j)
             work1(i,j)=depthu(iorig+i-1,jorig+j-1)
             depthu(i,j)=work1(i,j)
             work1(i,j)=depthv(iorig+i-1,jorig+j-1)
             depthv(i,j)=work1(i,j)
             work1(i,j)=plon(iorig+i-1,jorig+j-1)
             plon(i,j)=work1(i,j)
             work1(i,j)=plat(iorig+i-1,jorig+j-1)
             plat(i,j)=work1(i,j)
             work1(i,j)=ulon(iorig+i-1,jorig+j-1)
             ulon(i,j)=work1(i,j)
             work1(i,j)=ulat(iorig+i-1,jorig+j-1)
             ulat(i,j)=work1(i,j)
             work1(i,j)=vlon(iorig+i-1,jorig+j-1)
             vlon(i,j)=work1(i,j)
             work1(i,j)=vlat(iorig+i-1,jorig+j-1)
             vlat(i,j)=work1(i,j)
             work1(i,j)=scpx(iorig+i-1,jorig+j-1)
             scpx(i,j)=work1(i,j)
             work1(i,j)=scpy(iorig+i-1,jorig+j-1)
             scpy(i,j)=work1(i,j)
             work1(i,j)=scux(iorig+i-1,jorig+j-1)
             scux(i,j)=work1(i,j)
             work1(i,j)=scuy(iorig+i-1,jorig+j-1)
             scuy(i,j)=work1(i,j)
             work1(i,j)=scvx(iorig+i-1,jorig+j-1)
             scvx(i,j)=work1(i,j)
             work1(i,j)=scvy(iorig+i-1,jorig+j-1)
             scvy(i,j)=work1(i,j)
             work1(i,j)=cori(iorig+i-1,jorig+j-1)
             cori(i,j)=work1(i,j)

             work1(i,j)=ip(iorig+i-1,jorig+j-1)
             ip(i,j)=work1(i,j)
             work1(i,j)=iu(iorig+i-1,jorig+j-1)
             iu(i,j)=work1(i,j)
             work1(i,j)=iv(iorig+i-1,jorig+j-1)
             iv(i,j)=work1(i,j)
             work1(i,j)=iq(iorig+i-1,jorig+j-1)
             iq(i,j)=work1(i,j)
             work1(i,j)=scu2(iorig+i-1,jorig+j-1)
             scu2(i,j)=work1(i,j)
             work1(i,j)=scv2(iorig+i-1,jorig+j-1)
             scv2(i,j)=work1(i,j)
           enddo
         enddo
      endif


c --- float array allocation
      call floats_alloc
      write(lp,*) 'floats_alloc'
      call flush(lp)
c
      open (unit=nt1,file=flnm_fii,form='formatted',
     &      status='old',action='read' )
      open (unit=nt3,file=flnm_fo ,form='formatted',
     &      status='new',action='write')
c
c --- first and last model days from the input archives.
c
      call getdat_time(flnm_ai,time_ai,iexpt)
      call getdat_time(flnm_ao,time_ao,iexpt)
c
      narch = nint(fwbw*(time_ao-time_ai)/archfq) + 1
c
      write(lp,900) flnm_ai,flnm_ao
  900 format('first and last archives files:'//a120/a120)
      write(lp,910) narch,time_ai,time_ao
  910 format(/'number of archives, start and end times:'/i6,
     &       1p,2e16.8)
c
c --- initialize the floats
c
      iar=1
      call getdat_name(flnm_ai,time_ai,yrflag)
      write(lp,'(/a,i4,f10.2,/a)') 'iar,time_ai,flnm_ai = ',
     &                        iar,time_ai,flnm_ai(1:len_trim(flnm_ai))
      call flush(lp)
c
      write(lp,*) flnm_ai,time_ai,iexpt,yrflag,trcout
	  call flush(lp)
c
c first =1 means read archive files otherwise read reduced binary files
      if (first .eq. 1) then
         call getdat_sub(flnm_ai,time_ai,iexpt,yrflag,trcout,
     &                 iorig,jorig,idms,jdms)
      else
         call getdat_bin(flnm_ai,time_ai,iexpt,yrflag,trcout,
     &                 iorig,jorig,idms,jdms)
      endif
      if(wbin .eq. 1) then
         call writebin(flnm_ai,time_ai,iexpt,yrflag,trcout,
     &                 iorig,jorig,idms,jdms)
      endif
c      write(lp,*)'222 plat,platf',plat(1,1),
c     &   plat(ii,1)
c      call flush(lp)

c
c      open(unit=nt3,file=flnm_fo,status='unknown',
c     &     form='formatted',position='append')
c      do j=1,jdms
c       do i=1,idms
c        write(nt3,991) i,j,saln(i,j,1)
c 991    format(i6,i6,f12.4)
c       enddo
c      enddo
c      close(nt3)
      do j= 1,jj
        do i= 1,ii
            ubaro(i,j)=fwbw*ubaro(i,j) !allows frontwards/backwards
            vbaro(i,j)=fwbw*vbaro(i,j) 
            do k= 1,kk
               u(i,j,k)=fwbw*u(i,j,k) 
               v(i,j,k)=fwbw*v(i,j,k) 
            enddo !k
        enddo !i
      enddo !j

      call floats_init
c ---     initial archive
c
c         calculate dpu and dpv at u, v points.
c

          do j= 1,jj
            do i= 1,ii
              duk           = 0.0
              dvk           = 0.0
              depthi(:,:,0) = 0.0
              do k= 1,kk
                depthi(1,1,k) =
     &          depthi(1,1,k-1) + dp(    i,     j,     k)
                depthi(0,1,k) =
     &          depthi(0,1,k-1) + dp(max(i-1,1),j,     k)
                if(lperiod.and.i.eq.1) then
                  depthi(0,1,k) =
     &            depthi(0,1,k-1) + dp(idm,     j,     k)
                endif
                depthi(1,0,k) =
     &          depthi(1,0,k-1) + dp(    i, max(j-1,1),k)
                dukm1 = duk
                duk   = min( depthu(i,j),
     &                       0.5*(depthi(0,1,k) + depthi(1,1,k)) )
                dvkm1 = dvk
                dvk   = min( depthv(i,j),
     &                       0.5*(depthi(1,0,k) + depthi(1,1,k)) )
                dpu(i,j,k) = max( 0.0, duk-dukm1 )
                dpv(i,j,k) = max( 0.0, dvk-dvkm1 )
              enddo !k
            enddo !i
          enddo !j
c
c --- initial call to floats_advect writes out float info. at initial time
c
      call floats_advect_back(fwbw)
c
      flnm1=flnm_ai
      flnm2=flnm_ai
c
      do iar=2,narch
c
c ---   loop through all input model archives
c
        if (nintrp.gt.1) then
          time1=time_ai+fwbw*(iar-2)*archfq
          call getdat_name(flnm1,time1,yrflag)
          write(lp,'(a,i4,f10.2,/a)')'iar,time1,flnm1 = ',
     &                                iar,time1,flnm1(1:len_trim(flnm1))
          call flush(lp)
        endif
        time2=time_ai+fwbw*(iar-1)*archfq
        call getdat_name(flnm2,time2,yrflag)
        write(lp,'(a,i4,f10.2,/a)') 'iar,time2,flnm2 = ',
     &                               iar,time2,flnm2(1:len_trim(flnm2))
        call flush(lp)
c
      if (first .eq. 1) then
         call getdat_sub(flnm2,time2,iexpt,yrflag,trcout,
     &                 iorig,jorig,idms,jdms)
      else
         call getdat_bin(flnm2,time2,iexpt,yrflag,trcout,
     &                 iorig,jorig,idms,jdms)
      endif
      if(wbin .eq. 1) then
         call writebin(flnm2,time2,iexpt,yrflag,trcout,
     &                 iorig,jorig,idms,jdms)
      endif
c
        if (nintrp.le.1) then
c
c ---     no temporal interpolation
c
c         calculate dpu and dpv at u, v points.
c
          do j= 1,jj
            do i= 1,ii
              duk           = 0.0
              dvk           = 0.0
              depthi(:,:,0) = 0.0
              do k= 1,kk
                depthi(1,1,k) =
     &          depthi(1,1,k-1) + dp(    i,     j,     k)
                depthi(0,1,k) =
     &          depthi(0,1,k-1) + dp(max(i-1,1),j,     k)
                if(lperiod.and.i.eq.1) then
                  depthi(0,1,k) =
     &            depthi(0,1,k-1) + dp(idm,     j,     k)
                endif
                depthi(1,0,k) =
     &          depthi(1,0,k-1) + dp(    i, max(j-1,1),k)
                dukm1 = duk
                duk   = min( depthu(i,j),
     &                       0.5*(depthi(0,1,k) + depthi(1,1,k)) )
                dvkm1 = dvk
                dvk   = min( depthv(i,j),
     &                       0.5*(depthi(1,0,k) + depthi(1,1,k)) )
                dpu(i,j,k) = max( 0.0, duk-dukm1 )
                dpv(i,j,k) = max( 0.0, dvk-dvkm1 )
              enddo !k
            enddo !i
          enddo !j

          
c
          call floats_advect_back(fwbw)
c
        else
c
c ---     interpolate archives in time
c
c         fill alternate arrays for second archive.
c
          do j= 1,jj
            do i= 1,ii
              do k=1,kk
                temp2(i,j,k)=temp(i,j,k)
                saln2(i,j,k)=saln(i,j,k)
                th3d2(i,j,k)=th3d(i,j,k)
                u2(i,j,k)=u(i,j,k)
                v2(i,j,k)=v(i,j,k)
                dp2(i,j,k)=dp(i,j,k)
                if(trcout) then
                 do ktr=1,ntracr
                  tracr2(i,j,k,ktr)=tracer(i,j,k,ktr)
                 enddo !ktr
                endif !trcout
              enddo !k
            enddo !i
          enddo !j
c
          do jar=1,nintrp
c
            q=float(jar)/float(nintrp)
            if (q.lt.1) then
c
              if (jar.eq.1) then
                if (first .eq. 1) then
                  call getdat_sub(flnm1,time1,iexpt,yrflag,trcout,
     &                 iorig,jorig,idms,jdms)
                else
                  call getdat_bin(flnm1,time1,iexpt,yrflag,trcout,
     &                 iorig,jorig,idms,jdms)
                endif
             if(wbin .eq. 1) then
                call writebin(flnm1,time1,iexpt,yrflag,trcout,
     &                 iorig,jorig,idms,jdms)
             endif

               do j= 1,jj
                 do i= 1,ii
                   do k=1,kk
                     temp1(i,j,k)=temp(i,j,k)
                     saln1(i,j,k)=saln(i,j,k)
                     th3d1(i,j,k)=th3d(i,j,k)
                     u1(i,j,k)=u(i,j,k)
                     v1(i,j,k)=v(i,j,k)
                     dp1(i,j,k)=dp(i,j,k)
                     if(trcout) then
                      do ktr=1,ntracr
                       tracr1(i,j,k,ktr)=tracer(i,j,k,ktr)
                      enddo !ktr
                     endif !trcout
                   enddo !k
                 enddo !i
               enddo !j
              endif
c
c ---         time interpolate
c
              do j=1,jj
                do i=1,ii
                  ubaro(i,j)=fwbw*ubaro(i,j)
                  vbaro(i,j)=fwbw*vbaro(i,j)
                  do k=1,kk
                    temp(i,j,k)=(1.0-q)*temp1(i,j,k)+q*temp2(i,j,k)
                    saln(i,j,k)=(1.0-q)*saln1(i,j,k)+q*saln2(i,j,k)
                    th3d(i,j,k)=(1.0-q)*th3d1(i,j,k)+q*temp2(i,j,k)
                    u(i,j,k)=fwbw*((1.0-q)*u1(i,j,k)+q*u2(i,j,k))
                    v(i,j,k)=fwbw*((1.0-q)*v1(i,j,k)+q*v2(i,j,k))
                    dp(i,j,k)=(1.0-q)*dp1(i,j,k)+q*dp2(i,j,k)
                    if(trcout) then
                     do ktr=1,ntracr
                     tracer(i,j,k,ktr)=(1.0-q)*tracr1(i,j,k,ktr)
     &                     +q*tracr2(i,j,k,ktr)
                     enddo !ktr
                    endif !trcout
                  enddo !k
                enddo !i
              enddo !j
            else
              do j=1,jj
                do i=1,ii
                  do k=1,kk
                    temp(i,j,k)=temp2(i,j,k)
                    saln(i,j,k)=saln2(i,j,k)
                    th3d(i,j,k)=th3d2(i,j,k)
                    u(i,j,k)=fwbw*u2(i,j,k)
                    v(i,j,k)=fwbw*v2(i,j,k)
                    dp(i,j,k)=dp2(i,j,k)
                    if(trcout) then
                     do ktr=1,ntracr
                     tracer(i,j,k,ktr)=tracr2(i,j,k,ktr)
                     enddo !ktr
                    endif !trcout
                  enddo !k
                enddo !i
              enddo !j
c
            endif !q
c
c           calculate dpu and dpv at u, v points.
c
            do j= 1,jj
              do i= 1,ii
                duk           = 0.0
                dvk           = 0.0
                depthi(:,:,0) = 0.0
                do k= 1,kk
                  depthi(1,1,k) = 
     &            depthi(1,1,k-1) + dp(    i,     j,     k)
                  depthi(0,1,k) = 
     &            depthi(0,1,k-1) + dp(max(i-1,1),j,     k)
                  if(lperiod.and.i.eq.1) then
                    depthi(0,1,k) =
     &              depthi(0,1,k-1) + dp(idm,     j,     k)
                  endif
                  depthi(1,0,k) = 
     &            depthi(1,0,k-1) + dp(    i, max(j-1,1),k)
                  dukm1 = duk
                  duk   = min( depthu(i,j),
     &                         0.5*(depthi(0,1,k) + depthi(1,1,k)) )
                  dvkm1 = dvk
                  dvk   = min( depthv(i,j),
     &                         0.5*(depthi(1,0,k) + depthi(1,1,k)) )
                  dpu(i,j,k) = max( 0.0, duk-dukm1 )
                  dpv(i,j,k) = max( 0.0, dvk-dvkm1 )
                enddo !k
              enddo !i
            enddo !j
c
          call floats_advect_back(fwbw)


c
          enddo !jar
c
        endif !nintrp
c
      enddo  ! iar
c
      call floats_restart
c
      stop
c
      end
