! ************************************************************************ 
! 
!      PROGRAM GW3D: TRANSIENT, NON HOMOGENEOUS  
!			       MAINFRAME F90 VERSION 
! 
!  This is a Finite-Difference program for simulating 3D-Groundwater 
!   flow in variably-saturated porous media. 
! 
!        WRITTEN BY G.A. OLYPHANT  
!        Modified for Huntertown full domain, detailed Lagro, variable surface flux - July 2011
!        Options available for submodel domains, and option added to save output at each time step
! 
! ************************************************************************* 
! *************************************************************************
	module data_store
	save

      real,parameter :: RHO0=1.0 
      real,parameter :: G=976.0 
      real,parameter :: MU=0.01787 
      real,parameter :: GAMMA=G*RHO0 
      real,parameter :: BETA=4.4*10.0**(-10) 

      integer,parameter :: NIDC=8
      integer,parameter :: NROWIN=394
      integer,parameter :: NCOLIN=730
!     integer,parameter :: NLAYIN=138
      integer,parameter :: NLAYIN=69
!     NTIMEIN=30  ! number of time steps
      integer,parameter :: NTIMEIN=30  ! number of time steps
      integer,parameter :: NUMUNITS=33
	  
      integer,parameter :: SUBROW=394
      integer,parameter :: SUBCOL=730
      integer,parameter :: STARTROW=1
      integer,parameter :: ENDROW=394+2
      integer,parameter :: STARTCOL=1
      integer,parameter :: ENDCOL=730+2
		  	  	   
      integer,parameter :: NROW=NROWIN+2 
      integer,parameter :: NCOLUMN=NCOLIN+2 
      integer,parameter :: NLAYER=NLAYIN+2 

      integer,parameter :: ri = NROWIN+2		!row index
      integer,parameter :: ci = NCOLIN+2		!col index
      integer,parameter :: li = NLAYIN+2		!layer index

      integer,parameter :: itmax=7999
      integer,parameter :: start_itime=1        !mid-start
	  
      real,parameter :: EPSILON=0.1
!     real,parameter :: EPSILON=100   !just for quick testing
      real,parameter :: DXIN=50.0
      real,parameter :: DYIN=50.0 
      real,parameter :: DZL=0.5
!     real,parameter :: DZU=1
!     real,parameter :: Z0=210
      real,parameter :: Z0=208
      real,parameter :: H0=228-10           !! lowest expected saturation elevation
	  
!	   
!  CONVERT GRID SPACING AND Z0 FROM M TO CM 
! 
      real,parameter :: DX=DXIN*100.0 
      real,parameter :: DY=DYIN*100.0 

!     CONVERT COMPRESSIBILITY OF WATER TO UNIT OF CM^-1 
!     RECOMMENDED IN FREEZE, 1971; BOTTOM OF P 349 
!     ALPHA AND BETA ASUMED TO BE PA-1 
! 
      real,parameter :: GAMMASI=9.8*10.0**3 
      real,parameter :: BETAP=BETA*GAMMASI/100.0 

      integer,parameter :: IFIXP1=71
      integer,parameter :: JFIXP1=339
      integer,parameter :: IFIXP2=124
      integer,parameter :: JFIXP2=338
      integer,parameter :: IFIXP3=252
      integer,parameter :: JFIXP3=400
      integer,parameter :: IFIXP4=103
      integer,parameter :: JFIXP4=253
      integer,parameter :: IFIXP5=270
      integer,parameter :: JFIXP5=233
      integer,parameter :: IFIXP6=146
      integer,parameter :: JFIXP6=234
      integer,parameter :: IFIXP7=50
      integer,parameter :: JFIXP7=680
      integer,parameter :: IFIXP8=94
      integer,parameter :: JFIXP8=85
      integer,parameter :: IFIXP9=202
      integer,parameter :: JFIXP9=470
      integer,parameter :: IFIXP10=20
      integer,parameter :: JFIXP10=300
	  	  
!!	  integer,dimension(610,461) :: KTOP,KBOT
!!    real,dimension(610,461) :: ELEV

      integer,allocatable,dimension(:,:) :: ktop,kbot,isoils
      real,allocatable,dimension(:,:) :: elev
      
!!      integer,dimension(610,461,123) :: IDCELL
!!      real,dimension(610,461,123) :: QX,QY,QZ,PSI,H, PSIPRED,HIN,PSIIN,PSIOLD,PSIOLD2,PERM1,PERM2,PERM3,PERM4,PERM5,PERM6,COND

      integer,allocatable,dimension(:,:,:) :: IDCELL
      real,allocatable,dimension(:,:,:) ::                                    &
           QX,QY,QZ,PSI,H,PSIPRED,HIN,PSIIN,PSIOLD,PSIOLD2,PERM1,PERM2,PERM3,PERM4,PERM5,PERM6,COND
      real,allocatable,dimension(:,:) :: evin
      
      real,dimension(nlayer) :: DELTAZ,Z,IXE,IXW,IYN,IYS
      real,dimension(nidc) :: PRES,PA,ALPHA,POROS,PM,PN,HV 

    end module data_store

! ******************************************************************** 
! 
!          MAIN PROGRAM 
! 
    PROGRAM GW3D 
! 
      use data_store
      real,DIMENSION(ri,ci,ntimein):: izt
      REAL IZB
      allocate (psiin(ri,ci,li))
      allocate (psi(ri,ci,li))
      allocate (psiold(ri,ci,li))
      allocate (psiold2(ri,ci,li))
      allocate (psipred(ri,ci,li))
      allocate (h(ri,ci,li))
      allocate (qx(ri,ci,li),qy(ri,ci,li),qz(ri,ci,li))
      allocate (evin(ri,ci))

! 
!      COMMON/G1/ KTOP(610,461),KBOT(610,461),IDCELL(610,461,123),ELEV(610,461)  
!      COMMON/G2/ QX(610,461,123),QY(610,461,123),QZ(610,461,123),PSI(610,461,123),H(610,461,123) 
!      COMMON/G3/ PSIPRED(610,461,123),PSIOLD(610,461,123),PSIOLD2(610,461,123),IXE(123),IXW(123),IYN(123),IYS(123),PERM1(610,461,123),PERM2(610,461,123),PERM3(610,461,123),PERM4(610,461,123),PERM5(610,461,123),PERM6(610,461,123) 
!      COMMON/G4/ HIN(610,461,123),PSIIN(610,461,123) 
!	   COMMON/G5/ DELTAZ(123),Z(123) 
!      COMMON/G7/ COND(610,461,123) 
!
!APPEND
      OPEN(UNIT=50,FILE='PSI_PROFILES1_71_339.DAT',STATUS='UNKNOWN',position='append')
      OPEN(UNIT=52,FILE='PSI_PROFILES2_124_338.DAT',STATUS='UNKNOWN',position='append')
      OPEN(UNIT=53,FILE='PSI_PROFILES3_252_400.DAT',STATUS='UNKNOWN',position='append')
      OPEN(UNIT=54,FILE='PSI_PROFILES4_103_253.DAT',STATUS='UNKNOWN',position='append')
      OPEN(UNIT=55,FILE='PSI_PROFILES5_270_233.DAT',STATUS='UNKNOWN',position='append')
      OPEN(UNIT=70,FILE='PSI_PROFILES6_146_234.DAT',STATUS='UNKNOWN',position='append')
      OPEN(UNIT=71,FILE='PSI_PROFILES7_50_680.DAT',STATUS='UNKNOWN',position='append')
      OPEN(UNIT=72,FILE='PSI_PROFILES8_94_85.DAT',STATUS='UNKNOWN',position='append')
      OPEN(UNIT=73,FILE='PSI_PROFILES9_202_470.DAT',STATUS='UNKNOWN',position='append')
      OPEN(UNIT=74,FILE='PSI_PROFILES10_20_300.DAT',STATUS='UNKNOWN',position='append')
! 
!  INITIALIZE ARRAYS AND SET CONSTANTS 
! 
! 
      PSIIN=0.0 
      PSIOLD=0.0 
      PSIOLD2=0.0
      PSIPRED=0.0 

      QX=0.0 
      QY=0.0 
      QZ=0.0 
! 
      IXE=0.0 
      IXW=0.0 
      IYN=0.0 
      IYS=0.0 
      IZB=0.0
      evin=0.0
      izt=0.0
! 
!!      DO 30 ITIME=1,365
!!      IZT(ITIME)=0.0 
!!   30 CONTINUE 

      IZT=0.0 
	  
!
!   SPECIFY RUN PARAMETERS
!     INCODE(1)= STEADY-STATE HEADS; INCODE(2)=STEADY-STATE WATER TABLE
      INCODE=1
!     IHDM = SURFACE FLUX DATA: 1=HOURLY 2=DAILY
      IHDM=2
!     DTHR= CALCULATION TIME-STEP (in hours)
!     DTHR=0.5
!     DTHR=6
      DTHR=4
!
      CALL READIN(INCODE)
!
!
!    TIME-DEPENDENT INPUT
!
!
      SSOPT=0.0
      KKOUNT=1
!
!     OPEN(UNIT=51,FILE='..\input\INF_CM_Jan2011.DAT',STATUS='OLD')
!     OPEN(UNIT=51,FILE='..\..\input\INF_CM_Mar2011_shortrun.DAT',STATUS='OLD')  
      OPEN(UNIT=51,FILE='..\..\input\#SFLUX_Huntertown_30days_2003.txt',STATUS='OLD')  !!! block of distributed surface flux, day 1 at top
!	  READ(51,*)(IZT(ITIME),ITIME=1,NTIMEIN) 

      do itime=1,NTIMEIN
         READ(51,*)((IZT(I,J,ITIME),J=2,NCOLUMN-1),I=2,NROW-1)   !!unformatted read - kinda dangerous.
!        do i=2,nrow-1
!           do j=2,ncolumn-1
!              if(izt(i,j,itime).gt.0.0)then
!                 izt(i,j,itime)=2.55    !homogenous rain input
!              end if
!           end do
!        end do
      end do
      
                 
!   WRITE(6,2)(ITIME,IZT(ITIME),ITIME=1,NTIMEIN)   !!!
2     FORMAT(I3,F7.1)
!
!  CONVERT THE SURFACE FLUX VALUES TO CM/S
!
!	DO 100 ITIME=start_itime,NTIMEIN   !!!!!!!!!!!!!!!mid-start
!      IF(IHDM.EQ.1)THEN                  !HOURLY
!      IZT(ITIME)=IZT(ITIME)/3600.0
!      ELSE
!      IZT(ITIME)=IZT(ITIME)/(24.0*3600.0)   !DAILY
!     ENDIF
!	  if(izt(itime).lt.0)izt(itime)=0.65*izt(itime)
!  100 CONTINUE

      IF(IHDM.EQ.1)THEN       !HOURLY
         IZT=IZT/3600.0
      ELSE
         IZT=IZT/(24.0*3600.0)   !DAILY
      ENDIF

!
      IF(IHDM.EQ.1)THEN
         NDT=1.0/DTHR
         TIME=0.0
      ELSE
         NDT=24.0/DTHR
         TIME=0.0
      ENDIF
      DTSEC=DTHR*3600.0
      ITIME=start_itime
!     ITIME=0       !mid-start
      if(itime.eq.0)then
         WRITE(50,3)ITIME,(PSIIN(IFIXP1,JFIXP1,K),K=NLAYER-1,2,-1)
         WRITE(52,3)ITIME,(PSIIN(IFIXP2,JFIXP2,K),K=NLAYER-1,2,-1)
         WRITE(53,3)ITIME,(PSIIN(IFIXP3,JFIXP3,K),K=NLAYER-1,2,-1)
         WRITE(54,3)ITIME,(PSIIN(IFIXP4,JFIXP4,K),K=NLAYER-1,2,-1)
         WRITE(55,3)ITIME,(PSIIN(IFIXP5,JFIXP5,K),K=NLAYER-1,2,-1)
         WRITE(70,3)ITIME,(PSIIN(IFIXP6,JFIXP6,K),K=NLAYER-1,2,-1)
         WRITE(71,3)ITIME,(PSIIN(IFIXP7,JFIXP7,K),K=NLAYER-1,2,-1)
         WRITE(72,3)ITIME,(PSIIN(IFIXP8,JFIXP8,K),K=NLAYER-1,2,-1)
         WRITE(73,3)ITIME,(PSIIN(IFIXP9,JFIXP9,K),K=NLAYER-1,2,-1)
         WRITE(74,3)ITIME,(PSIIN(IFIXP10,JFIXP10,K),K=NLAYER-1,2,-1)

!     (50) PSI_Profiles.dat
!!    3 FORMAT(I3,63F7.1)         
3        FORMAT(I3,<li>F7.1)  
      end if
!
!  INITIAL CONDITION FOR TRANSIENT SIMULATIONS: 
!        PSI=PSIOLD=PSIOLD2=PSIPRED=PSIIN
!
      PSI=PSIIN
      PSIPRED=PSIIN
      PSIOLD=PSIIN
      PSIOLD2=PSIIN
      H=HIN
!
!  TIME LOOP FOR TRANSIENT SIMULATIONS
!
      KOUNT=0
      WRITE(6,*)'                           '
      WRITE(6,*)'  '
      WRITE(6,*)' <<<< ITERATING .....PLEASE WAIT....>>>>> '
!
      DO ITIME=start_itime,NTIMEIN     !mid-start
         KOUNT=KOUNT+1
         II=194			!! test prism
         JJ=447			!! test prism
         DO K=NLAYER-1,2,-1
            WRITE(6,*)K,IDCELL(II,JJ,K),COND(II,JJ,K),PSI(II,JJ,K)
         end do
         WRITE(6,4)ITIME
4        FORMAT(/,2X,'.......now in time step',i4,'......')
!
!  INNER TIME-LOOP PERTAINING TO CALCULATION TIME-STEP
!
!    go to 99    !!uncomment to check prisms before beginning


         DO IDT=1,NDT
            DELTAT9=DTSEC
            DELTAT1=DELTAT9
            TIME=TIME+DTHR
            EVIN(:,:)=IZT(:,:,itime)  
!
!  CALL ITERATION SUBROUTINE
!
            CALL ITERATE(IDT,DELTAT9,IZB,SSOPT,ITIME,DELTAT1,ICONVIT)
!
         ENDDO
!
         IF(KOUNT.EQ.KKOUNT)THEN
!
!  WRITE SELECTED RESULTS TO OUTPUT FILE
!
            WRITE(50,3)ITIME,(PSI(IFIXP1,JFIXP1,K),K=NLAYER-1,2,-1)
            WRITE(52,3)ITIME,(PSI(IFIXP2,JFIXP2,K),K=NLAYER-1,2,-1)
            WRITE(53,3)ITIME,(PSI(IFIXP3,JFIXP3,K),K=NLAYER-1,2,-1)
            WRITE(54,3)ITIME,(PSI(IFIXP4,JFIXP4,K),K=NLAYER-1,2,-1)
            WRITE(55,3)ITIME,(PSI(IFIXP5,JFIXP5,K),K=NLAYER-1,2,-1)
            WRITE(70,3)ITIME,(PSI(IFIXP6,JFIXP6,K),K=NLAYER-1,2,-1)
            WRITE(71,3)ITIME,(PSI(IFIXP7,JFIXP7,K),K=NLAYER-1,2,-1)
            WRITE(72,3)ITIME,(PSI(IFIXP8,JFIXP8,K),K=NLAYER-1,2,-1)
            WRITE(73,3)ITIME,(PSI(IFIXP9,JFIXP9,K),K=NLAYER-1,2,-1)
            WRITE(74,3)ITIME,(PSI(IFIXP10,JFIXP10,K),K=NLAYER-1,2,-1)
!     (50) PSI_PROFILES.DAT
 	  
            call midwrite(itime)
 	
            KOUNT=0
         ENDIF
!
!
      ENDDO
!
!
!   PRINT SELECTED OUTPUT TO FILES
!
!
      CALL PPRINT(SSOPT,ITIME)
!
! 
! 
!   SCREEN PRINT OPTION 
! 
99    WRITE (6,*) 'SPECIFY COORDINATES OF VERTICAL PRISM TO BE PRINTED' 
      READ (5,*) II,JJ 
      IF(II.EQ.999)GO TO 2000 
      IF(II.GT.endROW.OR.JJ.GT.endCOL.OR.II.LT.startrow.OR.JJ.LT.startcol)GO TO 99 
! 
	  WRITE (6,86)  
   86 FORMAT(//,10X,'*******************',//) 
      WRITE(6,*)'       TOP       BOTTOM' 
	  WRITE(6,*)KTOP(II,JJ),KBOT(II,JJ) 
	  WRITE(6,86) 
 	  WRITE (6,*)'    K  IDC    COND         H       PSI         Z'      
      DO 400 K=NLAYER-1,2,-1 
      WRITE (6,87) K,IDCELL(II,JJ,K),COND(II,JJ,K),(H(II,JJ,K)/100.0),PSI(II,JJ,K),(Z(K)/100.0),PSIOLD(II,JJ,K) 
  400 CONTINUE 
   87 FORMAT(2I5,2X,E8.2,4F10.2) 
      GO TO 99
! 

 2000 STOP 
      END 
! 
! 
! 
! ********************************************************************** 
! 
!
      SUBROUTINE READIN(INCODE) 
! 
      use data_store
 !    REAL MU 
!	  DIMENSION ELEVM(610,461),HFIXED(610,461),IDSURF(610,461),IDFIX(610,461),IDCR(610,461),ELEVCR(610,461) 
	  real,allocatable,dimension(:,:) :: ELEVM,ELEVCR,HFIXED
	  integer,allocatable,dimension(:,:) :: IDSURF,IDFIX  !IDCR
	  character*30 header
	  integer val,inodata
	  	  
	  allocate (idcell(ri,ci,li))
	  allocate (hin(ri,ci,li))
	  allocate (cond(ri,ci,li))
!	  allocate (ELEVM(ri,ci),HFIXED(ri,ci),IDSURF(ri,ci),IDFIX(ri,ci),IDCR(ri,ci),ELEVCR(ri,ci))
	  allocate (ELEVM(ri,ci),HFIXED(ri,ci),IDSURF(ri,ci),IDFIX(ri,ci),ELEVCR(ri,ci))
	  allocate (kbot(ri,ci),ktop(ri,ci),elev(ri,ci))
	  allocate (isoils(ri,ci))
	  
! 
!!      COMMON/G1/ KTOP(610,461),KBOT(610,461),IDCELL(610,461,123)  
!!      COMMON/G4/ HIN(610,461,123),PSIIN(610,461,123) 
!!      COMMON/G5/ DELTAZ(123),Z(123) 
!!      COMMON/G6/ PRES(10),PA(10),ALPHA(10),POROS(10),PM(10),PN(10),K0(10),HV(10) 
!!      COMMON/G7/ COND(610,461,123) 
! 
! 
    OPEN(UNIT=28,FILE='..\..\input\ht_strmel2k_2009may5.txt',STATUS='old') 
!   OPEN(UNIT=28,FILE='..\..\input\ht_net2k_2009may5.txt',STATUS='old') 
!   OPEN(UNIT=29,FILE='..\input\ht_net2k_2009may5.txt',STATUS='old') 
    OPEN(UNIT=29,FILE='..\..\input\ht_domain_id.txt',STATUS='old')       ! 0 = outside subdomain, 1 = inside subdomain, 3 = stream cell
	OPEN(UNIT=31,FILE='IDMATRIX.TXT',STATUS='UNKNOWN') 
    OPEN(UNIT=32,FILE='ELEVDIF.DAT',STATUS='UNKNOWN') 
	OPEN(UNIT=35,FILE='..\..\input\Ht6_conductivitiesLagroDetail_2011Aug16.txt',STATUS='old') 
	OPEN(UNIT=36,FILE='IDGROUND.DAT',STATUS='UNKNOWN') 
    OPEN(UNIT=37,FILE='ELEVATIONS_M.DAT',STATUS='UNKNOWN') 
	OPEN(UNIT=39,FILE='ID_CHECK.DAT',STATUS='UNKNOWN') 
	OPEN(UNIT=135,FILE='..\..\input\ht_soilhydgrp.asc',STATUS='UNKNOWN')
 !
! 
!  SET INPUT FORMATS 
! 
! 
    1 FORMAT(4I4,2F5.0,I5,2F6.0,2F4.0) 
    2 FORMAT(<ncolin>I3) 
    3 FORMAT(<ncolin>F7.1) 
    4 FORMAT(<ncolin>F7.2) 
    5 FORMAT(<ncolin>f5.2) 
    6 FORMAT(<ncolin>I6) 
!    7 FORMAT(<ci>F7.2) 
!    8 FORMAT(<ci>I2) 
    7 FORMAT(<ncolin>F7.2) 
    8 FORMAT(<ncolin>I3) 
    9 FORMAT(<ncolin>F3.0) 
    10 FORMAT(<ncolin>I2) 
! 
!   READ THE NUMBERS OF COLUMNS, ROWS, LAYERS, AND MEDIA TYPES; DATUM, 
!        ITERATION PARAMETERS, AND GRID SPACINGS 
! 
!!      READ(30,*) NROWIN,NCOLIN,NLAYIN,NIDC,Z0,EPSILON,itmax,DX,DY,DZL,DZU 
!!      WRITE(6,*)NROWIN,NCOLIN,NLAYIN,NIDC,Z0,EPSILON,ITMAX,DX,DY,DZL,DZU 
!!      NCOLUMN=NCOLIN+2 
!!      NROW=NROWIN+2 
!!      NLAYER=NLAYIN+2 
!!	  H0=180.0 
! 
!  INITIALIZE ARRAYS 
! 
	HIN=0.0
	KBOT=0 
	KTOP=0 
!	IDCR=0 
	ELEVCR=0.0 
	IDFIX=0
	COND=0.0
	IDCELL=0 
	isoils=0

!     WRITE(6,*)'ARRAYS INITIALIZED' 
! 
!  CONVERT GRID SPACING AND Z0 FROM M TO CM 
! 
      DO 130 K=1,NLAYER 
!!	  IF(K.le.11)THEN 
	  DELTAZ(K)=DZL*100.0 
!!	  ELSE if(K.ge.12.and.K.le.84)then
!!	  DELTAZ(K)=DZU*100.0 
!!	  else 
!!	  DELTAZ(K)=DZL*100.0
!!	  ENDIF
	   
  130 CONTINUE 
! 
!  READ NODE IDENTIFICATION CODES 
! 
    DO 140 K=2,NLAYER-1	 
    READ(35,2)((IDCELL(I,J,K),J=2,NCOLUMN-1),I=2,NROW-1) 
	WRITE(31,2)((IDCELL(I,J,K),J=2,ncolumn-1),I=2,nrow-1)  
  140 CONTINUE  

 
! 
!  READ DOMAIN/MASK ARRAY THAT DEFINES THE SUBMODEL DOMAIN AND STREAM LOCATIONS
!  0 = outside subdomain, 1 = inside subdomain, 3 = stream cells (also inside subdomain) 
!  SET ALL BUT ONE LAYER OF UNIT 2 (Trafalgar) TO 0    !!!!!!!!!!!!!!! Ht CUSTOMIZATION
!
     READ(29,10)((IDFIX(I,J),J=2,NCOLUMN-1),I=2,NROW-1)   
          
     close(unit=35,status='keep') 
               close(unit=29,status='keep') 
               
     read(135,*)header,val,header,val,header,val,header,val,header,val,header,inodata
	 read(135,*)((isoils(i,j),J=2,NCOLUMN-1),I=2,NROW-1) 
	 close(unit=135,status='keep')

! for first two layers of Huntertown GFM, substitute soils data (hydrologic group codes, 1-5, ABCD and water/wetland)	 
	DO I=2,NROW-1 
	DO J=2,NCOLUMN-1 
    if(ktop(i,j).gt.0)then
    do k=ktop(i,j)-1,ktop(i,j)
	if(IDCELL(I,J,K).eq.6)IDCELL(I,J,K)=isoils(i,j)
	end do
	end if
	end do
	end do               
                  
! 
!  MASK OUTSIDE AND BELOW SUBDOMAIN
! 
	DO 150 I=2,NROW-1 
	DO 150 J=2,NCOLUMN-1 
	    IF(IDFIX(I,J).EQ.0)THEN 
	    DO 155 K=2,NLAYER-1  
	    IDCELL(I,J,K)=0 
	    155 CONTINUE
	    end if
    
150 CONTINUE 

    DO 151 I=2,NROW-1 
	DO 151 J=2,NCOLUMN-1 
	DO 151 K=2,NLAYER-1  
	IF(IDCELL(I,J,K).EQ.NUMUNITS.AND.IDCELL(I,J,K+1).EQ.NUMUNITS)THEN 
	    IDCELL(I,J,K)=0 
	ENDIF 
	151 CONTINUE

! 
!  FIND TOP AND BOTTOM OF EACH VERTICAL PRISM 
! 
	DO 160 I=startrow,endrow
    DO 160 J=startcol,endcol 
!   DO 170 K=1,NLAYER 
! 
!  OPTION TO SET BOUNDARIES AS FIXED HEADS 
!      IF(I.EQ.2.OR.J.EQ.2.AND.IDCELL(I,J,K).EQ.2)IDCELL(I,J,K)= 
!     &-IDCELL(I,J,K) 
!      IF(I.EQ.NROW-1.OR.J.EQ.NCOLUMN-1.AND.IDCELL(I,J,K).EQ.2) 
!     &IDCELL(I,J,K)=-IDCELL(I,J,K) 
!  170 CONTINUE 
    DO 180 K=2,NLAYER-1 
	IF(IDCELL(I,J,K).NE.0.AND.IDCELL(I,J,K-1).EQ.0)THEN 
	KBOT(I,J)=K 
	ENDIF 
    IF(IDCELL(I,J,K).NE.0.AND.IDCELL(I,J,K+1).EQ.0)THEN 
	KTOP(I,J)=K 
	ENDIF 
  180 CONTINUE 
  160 CONTINUE 
! 
!  SPECIFY HYDRAULIC PROPERTIES OF POROUS MEDIA: INITIALIZE 1ST 
! 

	PRES=0.0 
	PA=0.0 
	ALPHA=0.0 
	POROS=0.0 
	PM=0.0 
	PN=0.0 
	HV=1.0

      POROS(1)=0.427 
      POROS(2)=0.375
      POROS(3)=0.375
	  POROS(4)=0.484
	  POROS(5)=0.391
	  POROS(6)=0.375
	  POROS(7)=0.457
	  POROS(8)=0.457
	  PRES(1)=0.1
      PRES(2)=0.053
      PRES(3)=0.053
	  PRES(4)=0.09
	  PRES(5)=0.049
	  PRES(6)=0.053
	  PRES(7)=0.1
	  PRES(8)=0.1
      PA(1)=0.011 
      PA(2)=0.035 
      PA(3)=0.035
	  PA(4)=0.0076
	  PA(5)=0.032
	  PA(6)=0.035 
	  PA(7)=0.011
	  PA(8)=0.011
	  PN(1)=1.27
      PN(2)=3.18 
	  PN(3)=3.18
	  PN(4)=1.53
	  PN(5)=1.76
	  PN(6)=3.18
	  PN(7)=1.27
	  PN(8)=1.27
	  ALPHA(1)=1.0E-07 
	  ALPHA(2)=1.0E-09
      ALPHA(3)=1.0E-09
	  ALPHA(4)=1.0E-08 
	  ALPHA(5)=1.0E-08
	  ALPHA(6)=1.0E-09 
	  ALPHA(7)=1.0E-07
	  ALPHA(8)=1.0E-07
	  HV(1)=1.0
	  HV(2)=3.0
	  HV(3)=10.0
	  HV(4)=100.0
	  HV(5)=3.0
	  HV(6)=10.0
	  HV(7)=0.001		!fractured till   !!  10-5 vertical to 10-8 horizontal
	  HV(8)=1.0
!	  HV(7)=1.0 		!regular till
! 
!  CONVERT ALPHA(IDC) TO UNIT OF CM-1 (SEE P. 349 OF FREEZE, 1971) 
!               CALCULATE M FROM N USING STANDARD ASSUMPTION 
      DO 195 IP=1,NIDC 
      ALPHA(IP)=ALPHA(IP)*GAMMASI/100.0 !Gammasi=9800 (9.8m/s2 * 1000 kgr/m3) 
      PM(IP)=1.0-(1.0/PN(IP)) 
  195 CONTINUE   
! 
! 
!  READ INPUT HYDRAULIC HEAD VALUES (TRANSIENT SIMULATIONS ONLY) 
! 
      IF(INCODE.EQ.1)THEN 
 	  OPEN(UNIT=38, FILE='..\..\input\SSOUT_H_24aug2011.DAT', STATUS='OLD') 
! 	  OPEN(UNIT=38, FILE='..\..\input\OUTHEADS_MID.DAT',STATUS='OLD')   !mid-start
      DO 200 K=2,NLAYER-1 
      READ(38,4)((HIN(I,J,K),J=2,NCOLUMN-1),I=2,NROW-1) 
  200 CONTINUE 
!    4 FORMAT(254F7.0) 
      ENDIF 
      
      close(unit=38,status='keep')  
! 
!  INITIALIZE HEADS, IN M, FOR STEADY STATE SIMULATIONS 
! 
      IF(INCODE.EQ.0)THEN 
      DO 210 I=1,NROW 
	  DO 210 J=1,NCOLUMN 
	  DO 210 K=1,NLAYER 
	  HIN(I,J,K)=H0 
  210 CONTINUE 
      ENDIF 
      
!  CONVERT INPUT HEAD VALUES FROM M (OR FT) TO CM 
! 
      DO 220 I=startrow,endrow
      DO 220 J=startcol,endcol 
      DO 220 K=1,NLAYER  
!      DO 220 I=1,NROW 
!      DO 220 J=1,NCOLUMN 
!      DO 220 K=1,NLAYER 
      IF(IDCELL(I,J,K).EQ.0)THEN 
	  HIN(I,J,K)=0.0 
	  ELSE 
      HIN(I,J,K)=HIN(I,J,K)*100.0 
	  ENDIF 
      IF(IDCELL(I,J,K).NE.0.AND.HIN(I,J,K).EQ.0.0)THEN 
	  HIN(I,J,K)=H0*100.0 
	  WRITE(6,*)'INPUT PROBLEM AT CELL ',I,J,K,IDCELL(I,J,K) 
	  ENDIF 
  220 CONTINUE 
! 
! 
!  COMPUTE ELEVATION HEAD VALUES AND DETERMINE INITIAL PRESSURE HEADS 
!           
      Z(1)=Z0*100.0 
	  ZSUM=Z(1) 
      DO 300 K=2,NLAYER 
      Z(K)=ZSUM+((DELTAZ(K-1)/2.0)+(DELTAZ(K)/2.0)) 
	  ZSUM=Z(K) 
  300 CONTINUE 
!
	  ELEV=0.0 
	  ELEVM=0.0 
      DO 310 I=startrow,endrow 
      DO 310 J=startcol,endcol 
      DO 320 K=2,NLAYER 
	  IF(IDCELL(I,J,K).NE.0)THEN 
	  PSIIN(I,J,K)=HIN(I,J,K)-Z(K) 
      ENDIF 
  320 CONTINUE 
  310 CONTINUE 
!
      DO 330 I=startrow,endrow
	  DO 330 J=startcol,endcol
	  IF(KTOP(I,J).GT.0)THEN
      KT=KTOP(I,J)
      ELEV(I,J)=Z(KT) 
	  ELEVM(I,J)=ELEV(I,J)/100.0
	  ENDIF
  330 CONTINUE 
!      write(37,77)((elevm(i,j),j=2,ncolumn-1),i=2,nrow-1)
!   77 format(<ncolin>f4.0)
!      stop
! 
!    SET HEADS OF LAKE BOUNDARY CELLS = 175.5 M 
! 
!      HLAKE=175.5*100.0 
!	  DO 340 I=1,NROW 
!      DO 340 J=1,NCOLUMN 
!	  IF(KTOP(I,J).NE.0)THEN 
!	  IF(IDFIX(I,J).EQ.0.AND.IDFIX(I,J-1).EQ.2)THEN 
!	  KT=KTOP(I,J) 
!	  KB=KBOT(I,J) 
!	  DO 345 K=KB,51 
!	  HIN(I,J,K)=HLAKE 
!      PSIIN(I,J,K)=HIN(I,J,K)-Z(K) 
!	  IDCELL(I,J,K)=-1 
!  345 CONTINUE 
!      ENDIF 
!	  ENDIF 
!  340 CONTINUE 
! 
!  READ ARRAYS THAT DEFINES CREEK LOCATIONS AND ELEVATIONS 
!       SET FIXED HEADS=ELEV AT THOSE LOCATIONS  
! 
!     READ(29,*)((IDCR(I,J),J=2,NCOLUMN-1),I=2,NROW-1) 
	  READ(28,*)((ELEVCR(I,J),J=2,NCOLUMN-1),I=2,NROW-1)    !unformatted read, variable precision
! 
      DO 350 I=startrow,endrow
	  DO 350 J=startcol,endcol 
!	  IF(IDCR(I,J).NE.0)THEN 
      IF(IDFIX(I,J).EQ.3)THEN 
!      write(6,*)'creek cell ',i,j 
	  IF(ELEVCR(I,J).EQ.0.0)WRITE(6,*)'CREEK PROBLEM',I,J
	  KT=KTOP(I,J) 
	  KB=KBOT(I,J)
!	  IF(KB.LT.KT-2)THEN
!	  KT2=KT-2
!	  ELSE
!	  KT2=KB
!	  ENDIF
!	  IF(KT.LT.25)KT2=KT
	    IF(KB.EQ.0.OR.KT.EQ.0)THEN
	    WRITE(6,*)'NO KT OR KB',I,J
	    GO TO 350
	    ENDIF 
	  DO 360 K=KB,KT 
	  HIN(I,J,K)=ELEVCR(I,J)*100.0 
	  PSIIN(I,J,K)=HIN(I,J,K)-Z(K) 
!	  IDCELL(I,J,K)=-IDCELL(I,J,K) 
  360 CONTINUE
      IDCELL(I,J,KT)=-IDCELL(I,J,KT)
	  IDCELL(I,J,KT-1)=-IDCELL(I,J,KT-1)
	  ELEVM(I,J)=ELEVCR(I,J) 
      ENDIF 
  350 CONTINUE 
  
      close(unit=28,status='keep')  
! 
! 
!  READ HYDRAULIC CONDUCTIVITY/UNIT CODES  
! 
!      DO 400 K=2,NLAYER-1 
!	  READ (35,9) ((COND10(I,J,K),J=2,NCOLUMN-1),I=2,NROW-1) 
!  400 CONTINUE 
! 
! 
!   CALCULATE HYDRAULIC CONDUCTIVITIES BASED ON CONDUCTIVITY CODES AND 
!  RE-DEFINE MATERIAL (ID) CODES BASED ON THE INPUT CONDUCTIVITY CODES 
! 
	  DO I=startrow,endrow
	  DO J=startcol,endcol
	  DO K=1,NLAYER 
      IF(IDCELL(I,J,K).NE.0)THEN 
 
!     IF (COND10(I,J,K).NE.0.0) THEN 
	  KOND=ABS(IDCELL(I,J,K)) 
	  select case (kond) 
 
	  case(33,27,25,23,22,20,19,15,13) 
	  IF(KOND.EQ.33)COND(I,J,K)=7.5*10.0**(-8)
	  IF(KOND.EQ.27)COND(I,J,K)=1.0*10.0**(-8)
	  IF(KOND.EQ.25)COND(I,J,K)=1.0*10.0**(-6)
	  IF(kond.eq.23.or.kond.eq.22.or.kond.eq.20.or.kond.eq.19.or.kond.eq.15.or.kond.eq.13)&
	  &COND(I,J,K)=1.0*10.0**(-6)
	  if(IDCELL(I,J,K).LT.0)then 
	  IDCELL(I,J,K)=-1 
	  else 
      IDCELL(I,J,K)=1 
	  end if 
 	   
	  case(14,7,1)
	  IF(KOND.EQ.1)COND(I,J,K)=3.5*10.0**(-3)
	  IF(KOND.EQ.7.or.kond.eq.14)COND(I,J,K)=5.0*10.0**(-3) 
	  if(IDCELL(I,J,K).LT.0)then 
	  IDCELL(I,J,K)=-2 
	  else 
	  IDCELL(I,J,K)=2 
	  end if 
 	   
	  case(32,31,28,26,17,11) 
	  IF(KOND.EQ.32)COND(I,J,K)=4.79*10.0**(-2) 
	  IF(KOND.EQ.31)COND(I,J,K)=3.27*10.0**(-2) 
	  IF(KOND.EQ.28)COND(I,J,K)=5.2*10.0**(-2) 
	  IF(KOND.EQ.26.or.kond.eq.17.or.kond.eq.11)COND(I,J,K)=4.94*10.0**(-2) 
	  if(IDCELL(I,J,K).LT.0)then 
	  IDCELL(I,J,K)=-3 
	  else 
	  IDCELL(I,J,K)=3 
	  end if 
	   
	  case(30,18) 
	  IF(KOND.EQ.30.or.kond.eq.18)COND(I,J,K)=2.88*10.0**(-3) 
	  if(IDCELL(I,J,K).LT.0)then 
	  IDCELL(I,J,K)=-4 
	  else 
	  IDCELL(I,J,K)=4
	  end if 
	  	   
	  case(24,21,16,12,9) 
	  IF(KOND.EQ.24.or.kond.eq.21.or.kond.eq.16.or.kond.eq.12.or.kond.eq.9)COND(I,J,K)=1.08*10.0**(-2)
 	  if(IDCELL(I,J,K).LT.0)then 
	  IDCELL(I,J,K)=-5 
	  else 
	  IDCELL(I,J,K)=5 
	  end if 
 
	  case(29)  
      COND(I,J,K)=5.65*10.0**(-2)
	  if(IDCELL(I,J,K).LT.0)then 
	  IDCELL(I,J,K)=-6 
	  else 
      IDCELL(I,J,K)=6
	  end if 

	  case(10,8,6,4,3,2)  
	  IF(KOND.EQ.2)COND(I,J,K)=3.675*10.0**(-4)
	  IF(KOND.EQ.3.or.kond.eq.5)COND(I,J,K)=5.74*10.0**(-5)
	  IF(KOND.EQ.4)COND(I,J,K)=1.0*10.0**(-6)
	  IF(KOND.EQ.6.or.KOND.EQ.8.or.KOND.EQ.10)COND(I,J,K)=7.5*10.0**(-8)
	  if(IDCELL(I,J,K).LT.0.and.k.ge.ktop(i,j)-6)then 
	  IDCELL(I,J,K)=-7 
	  else if(IDCELL(I,J,K).ge.0.and.k.ge.ktop(i,j)-6)then
      IDCELL(I,J,K)=7
	  else if(IDCELL(I,J,K).LT.0.and.k.lt.ktop(i,j)-6)then 
	  IDCELL(I,J,K)=-8
	  else if(IDCELL(I,J,K).ge.0.and.k.lt.ktop(i,j)-6)then
      IDCELL(I,J,K)=8
	  end if 

      case default  
      COND(I,J,K)=0.0 
      IDCELL(I,J,K)=0 
      IF(KOND.GT.0.AND.COND(I,J,K).EQ.0.0)WRITE(6,*)I,J,K,KOND,COND(I,J,K)
	  end select 
  
	  end if
	  end do 
	  end do 
	  end do 
 
!     DO 650 I=1,NROW 
!	  DO 650 J=1,NCOLUMN 
!	  DO 650 K=1,NLAYER 
!	  IF(IDCELL(I,J,K).NE.0.AND.COND(I,J,K).EQ.0.0)WRITE(6,*)I,J,K,COND10(I,J,K) 
!  650 CONTINUE 
! 
! 
!  CHECK FOR PROBLEMS WITH ELEVATION ASSIGNMENTS
! 
      DO 410 I=startrow,endrow 
	  DO 410 J=startcol,endcol 
	  IF(ELEVCR(I,J).GT.0.0)THEN
	  ELEVDIF=ELEVCR(I,J)-ELEVM(I,J)
	  AELEVD=ABS(ELEVDIF)
	  IF(AELEVD.GT.2.0)WRITE(32,66)I,J,KTOP(I,J),ELEVCR(I,J),ELEVM(I,J),ELEVDIF 
   66 FORMAT(3I5,3F6.1)
!      IF(ELEVDIF.LT.-2.0.AND.KTOP(I,J).GT.25)THEN
!      KTOLD=KTOP(I,J)
!	  KTNEW=KTOLD-1
!	  IDCELL(I,J,KTOLD)=0 
!      ELEV(I,J)=Z(KTNEW) 
!	  ELEVM(I,J)=ELEV(I,J)/100.0 
!      KTOP(I,J)=KTNEW 
!      ENDIF 
	  ENDIF
  410 CONTINUE
!      STOP
! 
! 
!  PRINT OUTPUT FILES FOR CHECKING AND ADJUSTING 
! 
	  DO 460 K=1,NLAYER 
	  WRITE(39,8)((IDCELL(I,J,K),J=2,NCOLUMN-1),I=2,NROW-1) 
  460 CONTINUE 
! 
      DO 470 I=1,NROW 
	  DO 470 J=1,NCOLUMN 
	  IDSURF(I,J)=0 
	  IF(KTOP(I,J).NE.0)THEN 
	  DO 480 K=2,NLAYER-1 
	  IF(K.EQ.KTOP(I,J))THEN 
	  IDSURF(I,J)=IDCELL(I,J,K) 
	  ENDIF 
  480 CONTINUE 
      ENDIF 
  470 CONTINUE 
! 
      WRITE(36,8)((IDSURF(I,J),J=2,NCOLUMN-1),I=2,NROW-1) 
! 
      DO 490 I=1,NROW 
	  DO 490 J=1,NCOLUMN 
	  IF(KTOP(I,J).GT.0.AND.ELEVM(I,J).EQ.0.0)THEN
	  KT=KTOP(I,J)
	  ELEVM(I,J)=Z(KT)/100.0
	  ENDIF
	  IF(ELEVM(I,J).LT.Z0.AND.IDSURF(I,J).NE.0)WRITE(6,*)'ELEV PROBLEM',I,J,ELEVM(I,J),Z0,KTOP(I,J)
	  IF(IDSURF(I,J).LT.0)THEN 
	  ELEVM(I,J)=-ELEVM(I,J) 
	  ENDIF 
  490 CONTINUE 
! 
      WRITE(37,77)((ELEVM(I,J),J=2,NCOLUMN-1),I=2,NROW-1) 
   77 FORMAT(<ncolin>F7.1)    

!!	  deallocate(ELEVM,HFIXED,IDSURF,IDFIX,IDCR,ELEVCR)
      deallocate(ELEVM,HFIXED,IDSURF,IDFIX,ELEVCR)
! 
!      STOP 
! 
	  RETURN 
      END 
! 
! 
! 
! ********************************************************************** 
! 
! 
      SUBROUTINE ITERATE(IDT,DELTAT9,IZB,SSOPT,ITIME,DELTAT1,ICONVITY) 
! 
      use data_store
      DOUBLE PRECISION A,BT1,BT2,BT3,BT4,BT5,BT6,B1,B,C 
       
      INTEGER FAIL,TEXTRA 
! 
      REAL IZB,KBSAT
! 
      DIMENSION PERM(7),BPSI(7),IDCODE(7),KBSAT(7) 
      
	  real,allocatable,dimension(:,:,:) :: PSIOLDI,PSIPOLD
!	  real,allocatable,dimension(:,:) :: EVIN
	  
	  allocate (psioldi(ri,ci,li),psipold(ri,ci,li)) 
	  allocate (perm1(ri,ci,li),perm2(ri,ci,li),perm3(ri,ci,li),perm4(ri,ci,li),perm5(ri,ci,li),perm6(ri,ci,li))
!	  allocate (evin(ri,ci))
! 

!!      COMMON/G1/ KTOP(610,461),KBOT(610,461),IDCELL(610,461,123),ELEV(610,461)  
!!      COMMON/G2/ QX(610,461,123),QY(610,461,123),QZ(610,461,123),PSI(610,461,123),H(610,461,123) 
!!      COMMON/G3/ PSIPRED(610,461,123),PSIOLD(610,461,123),PSIOLD2(610,461,123),IXE(123),IXW(123),IYN(123),IYS(123),PERM1(610,461,123),PERM2(610,461,123),PERM3(610,461,123),PERM4(610,461,123),PERM5(610,461,123),PERM6(610,461,123) 
!!      COMMON/G5/ DELTAZ(123),Z(123) 
!!      COMMON/G7/ COND(610,461,123) 

!	  psiold2=0.0
!	  psiold2=psiin

	  PERM1=0.0 
	  PERM2=0.0 
	  PERM3=0.0 
	  PERM4=0.0 
	  PERM5=0.0 
 	  PERM6=0.0 
 	  
 	  idcode=0
 
	  OPEN(UNIT=60,FILE='iter_log.dat',status='unknown',position='append') 
!	WRITE(6,*)'IN ITERATE' 
! 
      KKNT=0 
      ITER=0 
	  EPSI=EPSILON 
      DT=DELTAT9 
      LAMBDA=0.0 
      IF(SSOPT.EQ.1.0)THEN 
	  OMEGA=1.88
      ELSE 
	  OMEGA=0.8
!	  OMEGA=1.2
      ENDIF 
! 
!  COMPUTE VARIABLE TIME-STEP COEFFICIENT USED IN CALCULATION 
!          OF PSIPRED:  FREEZE(71, EQ 12) 
! 
      T=DT/(DELTAT1*2.0) 
! 
!  INITIALIZE ITERATION BY SETTING PSIPOLD AND PSIOLDI 
!     (SEE  EQN 13 IN FREEZE, 1971) EQUAL TO ZERO 
! 
!  NOTES: PSIOLD = PSI IN PREVIOUS TIME STEP = PSIIN AT T=1 FOR TRANSIENT SIMULATIONS  
!            BUT = 0 FOR STEADY-STATE SIMULATIONS 
!         PSIOLDI= PSI IN PREVIOUS ITERATION = 0 IN 1ST ITERATION 
!         PSIPRED= PREDICTED PSI VALUE; SEE EQN'S BELOW 
!   
!!       DO 100 I=2,NROW-1 
!!       DO 100 J=2,NCOLUMN-1 
!!       DO 100 K=2,NLAYER-1 
!!       PSIPOLD(I,J,K)=0.0 
!!       PSIOLDI(I,J,K)=0.0   
!!  100 CONTINUE 

		PSIPOLD=0.0 
	    PSIOLDI=0.0   
! 
!  MAIN ITERATION LOOP 
! 
!      DO 1000 ITER=1,ITMAX 
!	  IF(KICKOUT.EQ.1)EXIT 
 1000 DIFMAX=0.0 
	  KSEEP=0 
	  ITER=ITER+1 
! 
!  OUTER DO-LOOPS FOR CYCLING THROUGH ALL PRISMS IN HORIZONTAL MESH 
! 
      DO 200 I=startrow,endrow 
      DO 210 J=startcol,endcol 
! 
!  INNER DO-LOOP CYCLES THROUGH ALL NODES OF VERTICAL PRISM 
! 
      DO 220 K=2,NLAYER-1 
! 
!  CHECK CELL STATUS AND SKIP PAST CONSTANT HEAD NODES 
! 
!     IF(IDCELL(I,J,K).GT.0.AND.IDCELL(I,J,K).LT.96)THEN 
      IF(IDCELL(I,J,K).LE.0)GO TO 220 
	  KKNT=KKNT+1 
! 
!  COMPUTE PSIPRED VALUES BASED ON EQ'S 11 - 13; FREEZE(1971) 
! 
!     EQUATION 11 
! 
	  IF(ITER .EQ. 1 .AND. ITIME.EQ.1)THEN 
      PSIPRED(I,J,K)=PSI(I,J,K) 
      ENDIF 
! 
!     EQUATION 12 
! 
      IF(ITER.EQ.1.AND.ITIME.GT.1)THEN 
      PSIPRED(I,J,K)=(T+1.0)*PSIOLD(I,J,K)-T*PSIOLD2(I,J,K) 
      ENDIF 
! 
!     EQUATION 13 
! 
	  IF(ITER.GT.1)THEN 
      PSIPRED(I,J,K)=PSIPOLD(I,J,K)+LAMBDA*(PSI(I,J,K)-PSIPOLD(I,J,K)) 
      ENDIF 
! 
!  COMPUTE PREDICTED PSI VALUES AT EACH BLOCK BOUNDARY 
! 
	  BPSI(7)=PSIPRED(I,J,K) 
      BPSI(1)=(BPSI(7)+PSIPRED(I,J+1,K))/2.0 
      BPSI(2)=(BPSI(7)+PSIPRED(I,J-1,K))/2.0 
      BPSI(3)=(BPSI(7)+PSIPRED(I+1,J,K))/2.0 
      BPSI(4)=(BPSI(7)+PSIPRED(I-1,J,K))/2.0 
      BPSI(5)=(BPSI(7)+PSIPRED(I,J,K+1))/2.0 
      BPSI(6)=(BPSI(7)+PSIPRED(I,J,K-1))/2.0 
! 
!  SET APPROPRIATE BOUNDARY VALUES EQUAL TO ZERO 
! 
    IF(IDCELL(I,J+1,K).EQ.0)BPSI(1)=0.0 
    IF(IDCELL(I,J-1,K).EQ.0)BPSI(2)=0.0 
    IF(IDCELL(I+1,J,K).EQ.0)BPSI(3)=0.0 
    IF(IDCELL(I-1,J,K).EQ.0)BPSI(4)=0.0 
    IF(IDCELL(I,J,K+1).EQ.0)BPSI(5)=0.0 
    IF(IDCELL(I,J,K-1).EQ.0)BPSI(6)=0.0 
! 
! 
!  CALCULATE CELL BOUNDARY VALUES OF PSI AND CONDUCTIVITY 
!                   THEN 
!  CALL SUBROUTINE TO ESTIMATE FUNCTIONAL PARAMETERS FROM 
!       VAN GUNUCHTEN CHARACTERISTIC EQUATIONS 
! 
! 
	IDCODE(1)=ABS(IDCELL(I,J+1,K)) 
	IF (IDCELL(I,J+1,K).NE.0) THEN 
    KBSAT(1)=SQRT(COND(I,J,K)*COND(I,J+1,K)) 
	ELSE 
	KBSAT(1)=0.0 
	ENDIF 
! 
	IDCODE(2)=ABS(IDCELL(I,J-1,K)) 
	IF (IDCELL(I,J-1,K).NE.0) THEN 
    KBSAT(2)=SQRT(COND(I,J,K)*COND(I,J-1,K)) 
	ELSE 
 	KBSAT(2)=0.0 
	ENDIF 
! 
	IDCODE(3)=ABS(IDCELL(I+1,J,K)) 
    IF (IDCELL(I+1,J,K).NE.0) THEN 
    KBSAT(3)=SQRT(COND(I,J,K)*COND(I+1,J,K)) 
	ELSE 
	KBSAT(3)=0.0 
	ENDIF 
! 
	IDCODE(4)=ABS(IDCELL(I-1,J,K)) 
	IF (IDCELL(I-1,J,K).NE.0) THEN 
    KBSAT(4)=SQRT(COND(I,J,K)*COND(I-1,J,K)) 
	ELSE 
	KBSAT(4)=0.0 
	ENDIF 
! 
    IDCODE(5)=ABS(IDCELL(I,J,K+1)) 
	IF (IDCELL(I,J,K+1).NE.0) THEN 
    KBSAT(5)=SQRT(COND(I,J,K)*COND(I,J,K+1)) 
	ELSE 
	KBSAT(5)=0.0 
	ENDIF 
! 
    IDCODE(6)=ABS(IDCELL(I,J,K-1)) 
	IF (IDCELL(I,J,K-1).NE.0) THEN 
    KBSAT(6)=SQRT(COND(I,J,K)*COND(I,J,K-1)) 
	ELSE 
	KBSAT(6)=0.0 
	ENDIF 
! 
    IDCODE(7)=ABS(IDCELL(I,J,K)) 
	KBSAT(7)=COND(I,J,K) 
! 
! 
!    IF(ITER.EQ.2)THEN
!	WRITE(6,*)I,J,K,(BPSI(II),II=1,7)
!	ENDIF
    CALL CHARACT(IDCODE,BPSI,PERM,THETA7,POROS7,CPSI7,ALPHAP,KBSAT) 
!	SAVTHET(I,J,K)=THETA7 
! 
!  ASSIGN PERM VALUES FOR FLOW COMPONENT CALCULATIONS 
! 
	PERM1(I,J,K)=PERM(1) 
	PERM2(I,J,K)=PERM(2) 
	PERM3(I,J,K)=PERM(3) 
	PERM4(I,J,K)=PERM(4) 
	PERM5(I,J,K)=PERM(5) 
	PERM6(I,J,K)=PERM(6) 
! 
!  UPDATE FOR NEXT ITERATION BY SETTING PSIOLDI=PSI 
! 
      PSIOLDI(I,J,K)=PSI(I,J,K) 
! 
! 
!  COMPUTE COEFFICIENTS TO BE UTILIZED IN GAUSS-SIEDEL TYPE (LSOR) 
!      ITERATIVE SOLUTION OF EQUATION 9 IN FREEZE, 1971 
! 
!      EQUATION 9 IN COEFFICIENT FORM IS: 
!       -A*PSI(I,J,K+1)+B*PSI(I,J,K)-C*PSI(I,J,K-1)=D 
! 
!    WHERE  B=B1+(G/MU)*(BT1+BT2+ . . . BT6) 
!    AND    D=DT1-DT2+ . . . -DT6+DT7 
! 
! 
!  COMPUTE A 
! 
      A=G*(RHO0**2)*PERM(5)/(MU*DELTAZ(K)*(DELTAZ(K)+DELTAZ(K+1))) 
! 
!  COMPUTE COMPONENTS OF B 
! 
	  B1=((RHO0*THETA7/POROS7)*(ALPHAP+BETAP*POROS7)+RHO0*CPSI7)/DT 
!  
      BT1=(RHO0**2)*PERM(1)/(DX*(DX+DX)) 
      BT2=(RHO0**2)*PERM(2)/(DX*(DX+DX)) 
      BT3=(RHO0**2)*PERM(3)/(DY*(DY+DY)) 
      BT4=(RHO0**2)*PERM(4)/(DY*(DY+DY)) 
      BT5=(RHO0**2)*PERM(5)/(DELTAZ(K)*(DELTAZ(K)+DELTAZ(K+1))) 
      BT6=(RHO0**2)*PERM(6)/(DELTAZ(K)*(DELTAZ(K)+DELTAZ(K-1))) 
! 
! 
!  COMPUTE C, ZERO OUT IF ON LOWER BOUNDARY 
! 
      C=G*(RHO0**2)*PERM(6)/(MU*DELTAZ(K)*(DELTAZ(K)+DELTAZ(K-1))) 
! 
! 
!  COMPUTE COMPONENTS OF D 
! 
! 
  	  DT1=(G*(RHO0**2)*PERM(1)/(MU*DX*(DX+DX)))*(PSI(I,J+1,K)+PSIOLD(I,J+1,K)-PSIOLD(I,J,K)) 
      DT2=(G*(RHO0**2)*PERM(2)/(MU*DX*(DX+DX)))*(PSIOLD(I,J,K)-PSI(I,J-1,K)-PSIOLD(I,J-1,K)) 
      DT3=(G*(RHO0**2)*PERM(3)/(MU*DY*(DY+DY)))*(PSI(I+1,J,K)+PSIOLD(I+1,J,K)-PSIOLD(I,J,K)) 
      DT4=(G*(RHO0**2)*PERM(4)/(MU*DY*(DY+DY)))*(PSIOLD(I,J,K)-PSI(I-1,J,K)-PSIOLD(I-1,J,K)) 
      DT5=(G*(RHO0**2)*PERM(5)/(MU*DELTAZ(K)))*(1.0+(1.0/(DELTAZ(K)+DELTAZ(K+1)))*(PSIOLD(I,J,K+1)-PSIOLD(I,J,K))) 
      DT6=(G*(RHO0**2)*PERM(6)/(MU*DELTAZ(K)))*(1.0+(1.0/(DELTAZ(K)+DELTAZ(K-1)))*(PSIOLD(I,J,K)-PSIOLD(I,J,K-1))) 
      DT7=((RHO0*THETA7/POROS7)*(ALPHAP+POROS7*BETAP)+(RHO0*CPSI7))*(PSIOLD(I,J,K)/DT) 
! 
! 
!  IMPOSE BOUNDARY CONDITIONS BY REPLACING APPROPRIATE SUBTERMS 
!          WITH SPECIFIED FLUX VALUES 
! 
!     TOP BOUNDARY 
!		  
      IF(IDCELL(I,J,K+1).EQ.0)THEN 
      A=0.0 
      BT5=0.0 
	  IF (EVIN(i,j).GT.COND(I,J,K)) THEN 
 	  SFLUX=COND(I,J,K) 
	  ELSE 
	  SFLUX=EVIN(i,j) 
	  ENDIF 
      DT5=RHO0*SFLUX/DELTAZ(K) 
      ENDIF 
! 
!     BOTTOM BOUNDARY 
! 
      IF(IDCELL(I,J,K-1).EQ.0)THEN 
      C=0.0 
      BT6=0.0 
      DT6=RHO0*IZB/DELTAZ(K) 
      ENDIF 
! 
!     NORTH (UPPER GRID) BOUNDARY 
! 
      IF(IDCELL(I-1,J,K).EQ.0)THEN 
      BT4=0.0 
      DT4=RHO0*IYN(K)/DY 
      ENDIF 
! 
!     SOUTH (LOWER GRID) BOUNDARY 
! 
      IF(IDCELL(I+1,J,K).EQ.0)THEN 
      BT3=0.0 
      DT3=RHO0*IYS(K)/DY 
      ENDIF 
! 
!     WEST (LEFT GRID) BOUNDARY 
! 
      IF(IDCELL(I,J-1,K).EQ.0)THEN 
      BT2=0.0 
      DT2=RHO0*IXW(K)/DX 
      ENDIF 
! 
!     EAST (RIGHT GRID) BOUNDARY 
! 
      IF(IDCELL(I,J+1,K).EQ.0)THEN 
      BT1=0.0 
      DT1=RHO0*IXE(K)/DX 
      ENDIF 
! 
!  MAKE APPROPRIATE CORRECTIONS FOR STEADY-STATE SIMULATION 
!         AND CALCULATE B 
! 
      B2=(G/MU)*(BT1+BT2+BT3+BT4+BT5+BT6) 
! 
      IF(SSOPT.EQ.1.0) THEN 
      A=2.0*A 
      B1=0.0 
      B2=2.0*B2 
      C=2.0*C 
      DT1=2.0*DT1 
      DT2=2.0*DT2 
      DT3=2.0*DT3 
      DT4=2.0*DT4 
!	  DT5=0.0 
!	  DT6=0.0 
      DT7=0.0 
      ENDIF 
! 
      B=B1+B2 
      D=DT1-DT2+DT3-DT4+DT5-DT6+DT7 
! 
! 
!  SOLVE FOR NEW PSI VALUE BASED ON EQUATION 10 IN FREEZE, 1971 
! 
      PSI(I,J,K)=(D+A*PSI(I,J,K+1)+C*PSI(I,J,K-1))/B 
! 
!  OVERRELAXATION OR UNDERRELAXATION ACCORDING TO VALUE OF OMEGA 
! 
      PSI(I,J,K)=OMEGA*PSI(I,J,K)+(1.0-OMEGA)*PSIOLDI(I,J,K) 
	  if(psi(i,j,k).lt.-1000.0)psi(i,j,k)=-1000
!
      IF(K.NE.KTOP(I,J))GO TO 155 
	  PSICHK=PSI(I,J,K)+Z(K) 
	  IF(PSICHK.GT.Z(K))THEN 
	  PSI(I,J,K)=1.0
	  ENDIF  
! 
!  CHECK FOR DEVELOPMENT OF SEEPAGE FACES  
! 
  155 IF(PSI(I,J,K).LT.0.0.OR.K.EQ.KTOP(I,J))GO TO 225 
!     CHECK EAST AND WAST SIDE FOR ZERO (EMPTY CELL) 
      IF(IDCELL(I,J+1,K).EQ.0.OR.IDCELL(I,J-1,K).EQ.0)GO TO 221 
      GO TO 222 
!     VERIFY EMPTY EAST OR WEST CELL IS AIR AND NOT BOUNDARY 
  221 IF(KBOT(I,J+1).LT.K.AND.KTOP(I,J+1).LT.K.AND.KTOP(I,J+1).GT.0)& 
     &GO TO 224 
      IF(KBOT(I,J-1).LT.K.AND.KTOP(I,J-1).LT.K.AND.KTOP(I,J-1).GT.0)& 
     &GO TO 224 
!     CHECK NORTH AND SOUTH SIDE FOR ZERO (EMPTY CELL) 
  222 IF(IDCELL(I+1,J,K).EQ.0.OR.IDCELL(I-1,J,K).EQ.0)GO TO 223 
      GO TO 225 
!     VERIFY EMPTY NORTH OR SOUTH CELL IS AIR AND NOT BOUNDARY 
  223 IF(KBOT(I+1,J).LT.K.AND.KTOP(I+1,J).LT.K.AND.KTOP(I+1,J).GT.0)& 
     &GO TO 224  
	  IF(KBOT(I-1,J).LT.K.AND.KTOP(I-1,J).LT.K.AND.KTOP(I-1,J).GT.0)& 
     &GO TO 224 
      GO TO 225 
  224 PSI(I,J,K)=0.1 
	  KSEEP=KSEEP+1 
! 
! 
!  UPDATE ITERATION DIFFERENCE AND SAVE PSIPRED FOR NEXT ITERATION 
! 
 225 DIFF=ABS(PSI(I,J,K)-PSIOLDI(I,J,K)) 
	 IF(DIFF.GT.DIFMAX) THEN 
	 DIFMAX=DIFF 
	 IBAD=I 
 	 JBAD=J 
	 KBAD=K 
	 ENDIF 
	 PSIPOLD(I,J,K)=PSIPRED(I,J,K) 
! 
!  END OF MAIN ITERATION LOOP 
! 
!      ENDIF    
!      WRITE(6,*)PSI(I,J,K)   
  220 CONTINUE 
  210 CONTINUE 
  200 CONTINUE 
! 
! 
!  CHECK FOR CONVERGENCE THEN ITERATE OR CYCLE AS APPROPRIATE 
! 
!        
!	IF(SSOPT.EQ.1.0)THEN
      XIZT=evin(194,447)*3600.0*24.0   !sflux for test prism
	  WRITE(6,1)ITIME,IDT,ITER,DIFMAX,XIZT,KSEEP 
	  WRITE(60,1)ITIME,IDT,ITER,DIFMAX,XIZT,KSEEP 
!	  WRITE(6,*) IBAD,JBAD,KBAD,PSI(IBAD,JBAD,KBAD) 
!    &,PSIOLDI(IBAD,JBAD,KBAD) 
!     ENDIF 
    1 FORMAT(1X,'T = ',I3,'(',I2,')',2X,'ITERATION #',I4,2X,'DIFMAX =',F10.4,2X,'IZT =',F5.2,2X,'# OF SEEPS =',I5) 
!     IF(DIFMAX.LE.EPSI.OR.ITER.EQ.ITMAX-1)THEN 
!	  KICKOUT=1 
!	  ENDIF 
! 1000 continue		!Greg - check this continue location		 
! 
      IF(DIFMAX.GT.EPSI.AND.ITER.LT.ITMAX)GO TO 1000 
      ICONVIT=ITER 
      WRITE(6,2)ITIME,IDT,ICONVIT,DIFMAX,KSEEP 
	  WRITE(60,2)ITIME,IDT,ICONVIT,DIFMAX,KSEEP 
    2 FORMAT(1X,'TIME-STEP # ',I3,'(',I2,')',4X,'CONVERGENCE ACHIEVED IN ',I5,' ITERATIONS',5X,'DIFMAX = ',F10.2,4X,'# OF SEEPS =',I5) 
 
! 
! 
!  RECOMPUTE HYDRAULIC HEAD VALUES THEN UPDATE BY SETTING PSIOLD2=PSIOLD 
!                              AND PSIOLD=PSI 
! 
      DO 300 I=startrow,endrow 
      DO 300 J=startcol,endcol 
      DO 300 K=2,NLAYER-1 
      IF(IDCELL(I,J,K) .NE. 0)THEN 
      H(I,J,K)=Z(K)+PSI(I,J,K) 
      PSIOLD2(I,J,K)=PSIOLD(I,J,K) 
      PSIOLD(I,J,K)=PSI(I,J,K) 
      ENDIF 
  300 CONTINUE 
! 
!   CALCULATE FLOW COMPONENTS IN CM/DAY 
! 
      DO 400 K=NLAYER-1,2,-1 
      DO 410 I=endrow,startrow,-1 
      DO 420 J=endcol,startcol,-1 
! 
!   QZ FOR SURFACE CELL 
	  if (idcell(i,j,k).ne.0.and.idcell(i,j,k+1).eq.0) then 
!     write(6,*)i,j,k
	  QZ(I,J,K)=-EVIN(i,j)*86400.0 
	  QX(i,j,k)=0.0 
	  QY(i,j,k)=0.0 
	  GO TO 420 
	  ENDIF 
! 
!  	QZ FOR INTERNAL CELL 
      IF(IDCELL(I,J,K).NE.0.AND.IDCELL(I,J,K-1).NE.0)THEN 
      QZ(I,J,K)=((-PERM5(I,J,K)*GAMMA/MU)*((H(I,J,K)-H(I,J,K-1))/DELTAZ(K)))*86400.0 
      ENDIF 
! 
!   QZ FOR BOTTOM CELL	 
      IF(IDCELL(I,J,K).NE.0.AND.IDCELL(I,J,K-1).EQ.0)THEN 
      QZ(I,J,K)=((-PERM5(I,J,K)*GAMMA/MU)*((H(I,J,K+1)-H(I,J,K))/DELTAZ(K)))*86400.0 
      ENDIF 
! 
!   QX AND QY FOR INTERNAL CELL	OR EAST (RIGHT) BOUNDARY OR NORTH (TOP) BOUNDARY 
      IF(IDCELL(I,J,K).NE.0.AND.IDCELL(I-1,J,K).GT.0.AND.IDCELL(I,J-1,K).GT.0)THEN	 
      QX(I,J,K)=((-PERM1(I,J,K)*GAMMA/MU)*(H(I,J,K)-H(I,J-1,K))/DX)*86400.0 
      QY(I,J,K)=((PERM4(I,J,K)*GAMMA/MU)*((H(I,J,K)-H(I-1,J,K))/DY))*86400.0 
	  ENDIF 
! 
!   QX FOR CELL ADJACENT TO FIXED HEAD 
	  IF (IDCELL(I,J,K).NE.0.AND.IDCELL(I,J+1,K).LT.0)THEN 
      QX(I,J,K)=((-PERM1(I,J,K)*GAMMA/MU)*(H(I,J+1,K)-H(I,J,K))/DX)*86400.0 									  
      ENDIF 
! 
!   QY FOR CELL ADJACENT TO FIXED HEAD 
	  IF (IDCELL(I,J,K).NE.0.AND.IDCELL(I+1,J,K).LT.0)THEN  
      QY(I,J,K)=((PERM4(I,J,K)*GAMMA/MU)*((H(I+1,J,K)-H(I,J,K))/DY))*86400.0 
	  ENDIF 
! 
!     QX FOR WEST (LEFT) BOUNDARY 
      IF(IDCELL(I,J,K).NE.0.AND.IDCELL(I,J-1,K).EQ.0)THEN 
      QX(I,J,K)=((-PERM1(I,J,K)*GAMMA/MU)*(H(I,J+1,K)-H(I,J,K))/DX)*86400.0 
      ENDIF 
! 
!     QY FOR SOUTH (LOWER) BOUNDARY  
      IF (IDCELL(I,J,K).NE.0.AND.IDCELL(I-1,J,K).EQ.0)THEN 
	  QY(I,J,K)=((PERM4(I,J,K)*GAMMA/MU)*((H(I+1,J,K)-H(I,J,K))/DY))*86400.0 
      ENDIF 
! 
  420 CONTINUE 
  410 CONTINUE    	 
  400 CONTINUE   
!      ENDIF 
! 
	deallocate (perm1,perm2,perm3,perm4,perm5,perm6,PSIOLDI,PSIPOLD)

!  RETURN TO MAIN PROGRAM 
! 
      RETURN 
      END 
! 
! 
! ********************************************************************* 
! 
! 
	SUBROUTINE CHARACT(IDCODE,BPSI,PERM,THETA7,POROS7,CPSI7,ALPHAP,KBSAT) 
! 
      use data_store
      REAL KBSAT
! 
      DIMENSION BPSI(7),IDCODE(7),PERM(7),KBSAT(7) 
! 
!!      COMMON/G6/ PRES(10),PA(10),ALPHA(10),POROS(10),PM(10),PN(10),K0(10),HV(10) 
! 
! CALCULATE CONDUCTIVITIES ON GRID-CELL BOUNDARIES 
! 
      IDC7=IDCODE(7) 
! 
      DO 100 I=1,7 
      IF(IDCODE(I).EQ.0)THEN 
	  IDCODE(I)=IDC7 
	  ENDIF 

      IDC=IDCODE(I) 
      
!     IF(IDC.GE.10)IDC=IDC/10 
!!    IF(IDC.GE.nidc)IDC=IDC/10
! 
!     NOTE CORRECTION MUST BE MADE FOR WELL ID CODE 
!       IF (IDC .EQ. 9) IDC=7 
! 
! SET SATURATED-ZONE VALUES AND SKIP CALCULATIONS 
! 
	  IF(BPSI(I).GE.0.0)THEN	 
      CONDUCT=KBSAT(I) 
      COND7=KBSAT(7) 
      SE=1.0 
      SE7=1.0 
      GO TO 101 
      ELSE 
! 
! VAN GENUCHTEN (1980) EQUATIONS 
! 
	  SE1=1.0+((PA(IDC)*ABS(BPSI(I)))**PN(IDC)) 
	  SE17=1.0+((PA(IDC7)*ABS(BPSI(7)))**PN(IDC7)) 
	  SE=(1.0/SE1)**PM(IDC) 
      SE7=(1.0/SE17)**PM(IDC7) 
! 
      TP1=1.0-SE**(1.0/PM(IDC)) 
      TP17=1.0-SE7**(1.0/PM(IDC7)) 
      TP2=1.0-TP1**PM(IDC) 
      TP27=1.0-TP17**PM(IDC7) 
      TP3=TP2*TP2 
      TP37=TP27*TP27 
      CONDUCT=KBSAT(I)*SQRT(SE)*TP3 
      COND7=KBSAT(7)*SQRT(SE7)*TP37 
      ENDIF 
! 
  101 PERM(I)=CONDUCT*MU/GAMMA 
      PERM7=COND7*MU/GAMMA 
      PERM(I)=(PERM(I)+PERM7)/2.0
	  IF(I.EQ.5.OR.I.EQ.6.and.BPSI(7).GE.0.0)THEN
	  PERM(I)=PERM(I)/HV(IDC7)
      ENDIF 
  100 CONTINUE 
! 
!     CALCULATE SPECIFIC MOISTURE CAPACITY 
! 
      IF(BPSI(7).GE.0.0)THEN 
      CPSI7=0.0 
      THETA7=POROS(IDC7) 
      ALPHAP=ALPHA(IDC7) 
      ELSE 
! 
	  C1=(POROS(IDC7)-PRES(IDC7))*PM(IDC7)*PN(IDC7) 
      C2=(PA(IDC7)**PN(IDC7))*(ABS(BPSI(7))**(PN(IDC7)-1.0)) 
      CPSI7=C1*C2*((1.0/SE17)**(PM(IDC7)-1.0))/(SE17*SE17) 
      THETA7=PRES(IDC7)+(POROS(IDC7)-PRES(IDC7))*SE7 
      ALPHAP=0.0 
	  ENDIF 
! 
      POROS7=POROS(IDC7) 
! 
      RETURN 
      END 
!	     
! 
! ******************************************************************** 
! 
! 
      SUBROUTINE PPRINT(SSOPT,ITIME) 
! 
       use data_store

!!     COMMON/G1/ KTOP(610,461),KBOT(610,461),IDCELL(610,461,123),ELEV(610,461)  
!!     COMMON/G2/ QX(610,461,123),QY(610,461,123),QZ(610,461,123),PSI(610,461,123),H(610,461,123) 
!!	   COMMON/G5/ DELTAZ(123),Z(123) 
! 
!	  DIMENSION WT(610,461),LSEEP(610,461) 
	  real,allocatable,dimension(:,:) :: wt
	  integer,allocatable,dimension(:,:) :: lseep
         
	  allocate (wt(ri,ci),lseep(ri,ci))
	  
! 
	  OPEN(UNIT=40,FILE='OUTHEADS_H.DAT',STATUS='UNKNOWN') 
	  OPEN(UNIT=41,FILE='OUTPSI_H.DAT',STATUS='UNKNOWN') 
	  OPEN(UNIT=42,FILE='WT_H_Sept2011.TXT',STATUS='unknown')
!	  OPEN(UNIT=42,FILE='WT_H_july.TXT',STATUS='UNKNOWN')	   
      OPEN(UNIT=44,FILE='SSOUT_H.DAT',STATUS='unknown') 
!	  OPEN(UNIT=45,FILE='FLUXOUT_H.CRD',STATUS='UNKNOWN') 
!	  OPEN(UNIT=46,FILE='FLUXOUT_H.TXT ',STATUS='UNKNOWN') 
	  OPEN(UNIT=47,FILE='VECTORS_H.TXT',STATUS='UNKNOWN') 
	  OPEN(UNIT=48,FILE='SEEPS.DAT',STATUS='UNKNOWN')  
! 
    1 FORMAT(<ncolin>F7.2) 
    2 FORMAT (3F18.9,E13.7) 
    3 FORMAT (3I5,3F18.6) 
    4 FORMAT(<ci>F7.2) 
	5 FORMAT(<ncolin>F8.3)
	6 FORMAT(<ncolin>I2)
! 
      WRITE (6,*) "PPRINT ITIME", ITIME 
      
      IF(SSOPT.EQ.1.0)THEN		 
      DO 123 K=2,NLAYER-1 
      WRITE(44,1)(((H(I,J,K)/100.0),J=2,NCOLUMN-1),I=2,NROW-1) 
  123 CONTINUE  
      DO 111 K=2,NLAYER-1 
      WRITE(41,1)(((PSI(I,J,K)/100.0),J=2,NCOLUMN-1),I=2,NROW-1) 
  111 CONTINUE 
      ELSE 
      DO 115 K=2,NLAYER-1 
      WRITE(40,1)(((H(I,J,K)/100.0),J=2,NCOLUMN-1),I=2,NROW-1) 
  115 CONTINUE  
      DO 116 K=2,NLAYER-1 
      WRITE(41,5)(((PSI(I,J,K)/10.0),J=2,NCOLUMN-1),I=2,NROW-1) 
  116 CONTINUE 
      ENDIF 
      
! 
!  FIND THE WATER TABLE AT EACH PRISM AND PRINT TO FILES 
! 
	  WT=0.0 
	  LSEEP=0
	  rows: do I=startrow,endrow 
	  cols: do J=startcol,endcol
	  IF(KTOP(I,J).NE.0)THEN 
	  KT=KTOP(I,J) 
	  DO 121 K=1,NLAYER-1 
	  IF(PSI(I,J,K).EQ.0.1)THEN
	  LSEEP(I,J)=1
	  ENDIF 
	  IF(PSI(I,J,K).GE.0.0.AND.PSI(I,J,K+1).LT.0.0)THEN 
      WT(I,J)=H(I,J,K)/100.0 
      cycle cols
	  ENDIF 
  121 CONTINUE 
      IF(WT(I,J).EQ.0.0)THEN 
	  WT(I,J)=H(I,J,KT)/100.0 
	  ENDIF 
	  ENDIF 
      end do cols
      end do rows
      WRITE(42,1)((WT(I,J),J=2,NCOLUMN-1),I=2,NROW-1) 
	  WRITE(48,6)((LSEEP(I,J),J=2,NCOLUMN-1),I=2,NROW-1) 
!       
  200 CONTINUE 
! 
	  DO 210 I=2,NROW-1 
	  DO 220 J=2,NCOLUMN-1 
      DO 230 K=2,NLAYER-1 
	  IF(IDCELL(I,J,K).EQ.0)GO TO 230 
      IF(IDCELL(I+1,J,K).NE.0.AND.IDCELL(I-1,J,K).NE.0.AND.IDCELL(I,J+1,K).NE.0.AND.IDCELL(I,J-1,K).NE.0)THEN	 
	  WRITE(47,3) I,J,K,QX(I,J,K),QY(I,J,K),QZ(I,J,K) 
	  end if 
  230 CONTINUE 
  220 CONTINUE 
  210 CONTINUE	 	
  
  deallocate (wt,lseep,qx,qy,qz) 
!   
! 
	  RETURN 
      END 
! 
! ! 
! ******************************************************************** 
! 
! 
      SUBROUTINE MIDWRITE(itime)
! 
       use data_store

!!     COMMON/G1/ KTOP(610,461),KBOT(610,461),IDCELL(610,461,123),ELEV(610,461)  
!!     COMMON/G2/ QX(610,461,123),QY(610,461,123),QZ(610,461,123),PSI(610,461,123),H(610,461,123) 
!!	   COMMON/G5/ DELTAZ(123),Z(123) 
! 
!	  DIMENSION WT(610,461),LSEEP(610,461) 
	  real,allocatable,dimension(:,:) :: wt
      character*15 watertable
      character*3 count
      character*4 extension
         
	  allocate (wt(ri,ci))
	  
! 
	  OPEN(UNIT=140,FILE='OUTHEADS_MID.DAT',STATUS='REPLACE') 
!	  OPEN(UNIT=42,FILE='WT_H.TXT',STATUS='UNKNOWN')	   
! 
    1 FORMAT(<ncolin>F7.2) 
    2 FORMAT (3F18.9,E13.7) 
    3 FORMAT (3I5,3F18.6) 
    4 FORMAT(<ci>F7.2) 
	5 FORMAT(<ncolin>F8.3)
	6 FORMAT(<ncolin>I2)
! 
      WRITE (6,*) "PPRINT ITIME", ITIME 
      
      extension = ".txt"
     
      
      if(itime.lt.10)then
      write(count, '(I1)')itime
      else if(itime.ge.10.and.itime.lt.100)then
      write(count, '(I2)')itime
      else
      write(count, '(I3)')itime
      end if
      
      count = trim(count)
      
      watertable = 'wt_h_'//count//extension
      watertable = trim(watertable)
      OPEN(UNIT=42,FILE=watertable,STATUS='UNKNOWN')	
! 
!  FIND THE WATER TABLE AT EACH PRISM AND PRINT TO FILES 
! 
	  WT=0.0 
	  midrows: do I=startrow,endrow 
	  midcols: do J=startcol,endcol
	  IF(KTOP(I,J).NE.0)THEN 
	  KT=KTOP(I,J) 
	  DO 121 K=1,NLAYER-1 
	  IF(PSI(I,J,K).EQ.0.1)THEN
	  ENDIF 
	  IF(PSI(I,J,K).GE.0.0.AND.PSI(I,J,K+1).LT.0.0)THEN 
      WT(I,J)=H(I,J,K)/100.0 
      cycle midcols
	  ENDIF 
  121 CONTINUE 
      IF(WT(I,J).EQ.0.0)THEN 
	  WT(I,J)=H(I,J,KT)/100.0 
	  ENDIF 
	  ENDIF 
      end do midcols
      end do midrows
      WRITE(42,1)((WT(I,J),J=2,NCOLUMN-1),I=2,NROW-1) 
!       
  
      DO K=2,NLAYER-1 
      WRITE(140,1)(((H(I,J,K)/100.0),J=2,NCOLUMN-1),I=2,NROW-1) 
      end do
     
      CLOSE(UNIT=42,status='keep')   
      close(unit=140,status='keep')
! 
  
  deallocate (wt) 
	  
	  RETURN 
      END 
! 
