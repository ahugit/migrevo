MODULE pNelder_Mead
!ahu 010417: TOOK THIS FROM COHABITC/CLUSTER 040615. ONLY CHANGE FOR NOW WILL BE THE INCLUDE MPIF.H AND UNCOMMENTING OUT USE GLOBAL. 
!AND ALSO SOME CHANGED AT THE BEGINNING ABOUT HOW REALRANK IS OBTAINED (WAS EQUATING TO MYID BEFORE WHICH REQUIRED THE USE OF GLOBAL. 
!BUT NOW UNCOMMENTED OUT THE CALL MPI AT THE BEGINNING AND I OBTAIN REALRANK THAT WAY). 
use nrtype, only: dp
!use global  !, only: myid
!use mpi

IMPLICIT NONE
INCLUDE 'mpif.h'
!ahu f14 putting this in the subroutine so that it does not clash with nrtype INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
PRIVATE
PUBLIC :: pminim

CONTAINS


SUBROUTINE pminim(p, step, nop, func, maxfn, iprint, stopcr, nloop, iquad,  &
                 simp, var, functn, functn2, ifault, myrank,nprocs,T0, Tstep, Tfreq, randseed)
!     A PROGRAM FOR FUNCTION MINIMIZATION USING THE SIMPLEX METHOD.

!     FOR DETAILS, SEE NELDER & MEAD, THE COMPUTER JOURNAL, JANUARY 1965

!     PROGRAMMED BY D.E.SHAW,
!     CSIRO, DIVISION OF MATHEMATICS & STATISTICS
!     P.O. BOX 218, LINDFIELD, N.S.W. 2070

!     WITH AMENDMENTS BY R.W.M.WEDDERBURN
!     ROTHAMSTED EXPERIMENTAL STATION
!     HARPENDEN, HERTFORDSHIRE, ENGLAND

!     Further amended by Alan Miller
!     CSIRO Division of Mathematical & Information Sciences
!     Private Bag 10, CLAYTON, VIC. 3169

!     Fortran 90 conversion by Alan Miller, June 1995
!     Alan.Miller @ vic.cmis.csiro.au
!     Latest revision - 5 December 1999

!     Parallelization and Simulated Annealing added by Steven Laufer, March 2012
!     Parallelization follows algorithm of Lee and Wiswall, Comput Econ (2007) 30:171-187
!     Simulated Annealing follows Numerical Recipes (10.9)
!     Latest revision - 3 May 2012


!     ARGUMENTS:-
!     P()     = INPUT, STARTING VALUES OF PARAMETERS
!               OUTPUT, FINAL VALUES OF PARAMETERS
!     STEP()  = INPUT, INITIAL STEP SIZES
!     NOP     = INPUT, NO. OF PARAMETERS, INCL. ANY TO BE HELD FIXED
!     FUNC    = OUTPUT, THE FUNCTION VALUE CORRESPONDING TO THE FINAL
!                 PARAMETER VALUES.
!     maxfn     = INPUT, THE MAXIMUM NO. OF FUNCTION EVALUATIONS ALLOWED.
!               Say, 20 times the number of parameters, NOP.
!               If negative, the absolute value is the maximum number of function
!                   evaluations per processor. 
!     IPRINT  = INPUT, PRINT CONTROL PARAMETER
!                 < 0 NO PRINTING
!                 = 0 PRINTING OF PARAMETER VALUES AND THE FUNCTION
!                     VALUE AFTER INITIAL EVIDENCE OF CONVERGENCE.
!                 > 0 AS FOR IPRINT = 0 PLUS PROGRESS REPORTS AFTER EVERY
!                     IPRINT EVALUATIONS, PLUS PRINTING FOR THE INITIAL SIMPLEX.
!     STOPCR  = INPUT, STOPPING CRITERION.
!               The criterion is applied to the standard deviation of
!               the values of FUNC at the points of the simplex.
!     NLOOP   = INPUT, THE STOPPING RULE IS APPLIED AFTER EVERY NLOOP
!               FUNCTION EVALUATIONS.   Normally NLOOP should be slightly
!               greater than NOP, say NLOOP = 2*NOP.
!     IQUAD   = INPUT, = 1 IF FITTING OF A QUADRATIC SURFACE IS REQUIRED
!                      = 0 IF NOT
!               N.B. The fitting of a quadratic surface is strongly
!               recommended, provided that the fitted function is
!               continuous in the vicinity of the minimum.   It is often
!               a good indicator of whether a premature termination of
!               the search has occurred.
!     SIMP    = INPUT, CRITERION FOR EXPANDING THE SIMPLEX TO OVERCOME
!               ROUNDING ERRORS BEFORE FITTING THE QUADRATIC SURFACE.
!               The simplex is expanded so that the function values at
!               the points of the simplex exceed those at the supposed
!               minimum by at least an amount SIMP.
!     VAR()   = OUTPUT, CONTAINS THE DIAGONAL ELEMENTS OF THE INVERSE OF
!               THE INFORMATION MATRIX.
!     FUNCTN  = INPUT, NAME OF THE USER'S SUBROUTINE - ARGUMENTS (P,FUNC)
!               WHICH RETURNS THE FUNCTION VALUE FOR A GIVEN SET OF
!               PARAMETER VALUES IN ARRAY P.
!     FUNCTN2 = INPUT, FUNCTION TO BE CALLED BY MASTER AFTER EACH
!               NLOOP EVALUATIONS ON BEST POINT IN CURRENT SIMPLX,
!               IF BEST POINT OF SIMPLEX IS P WITH VALUE H, CALLS FUNCTN2(P,H)
!               e.g. Check if have point better than previous best, if so
!               calculate some numbers of interest and write them to a file
!****    FUNCTN, FUNCTN2 MUST BE DECLARED EXTERNAL IN THE CALLING PROGRAM.
!     IFAULT  = OUTPUT, = 0 FOR SUCCESSFUL TERMINATION
!                 = 1 IF MAXIMUM NO. OF FUNCTION EVALUATIONS EXCEEDED
!                 = 2 IF INFORMATION MATRIX IS NOT +VE SEMI-DEFINITE
!                 = 3 IF NOP < 1
!                 = 4 IF NLOOP < 1
! NOTE DEFINITIONS OF MYRANK AND NPROCS. COUNTER-INTUITIVE.
!     MYRANK   = WHICH GROUP THIS PROCESSOR IS PART OF
!     NPROCS   = NUMBER OF GROUPS
!     T0       = STARTING TEMPERATURE FOR ANNEALING
!     Tstep    = FRACTION BY WHICH TEMPERATURE IS DECREASED AT EACH STEP
!     Tfreq    = FREQUENCY (NUMBER OF FUNCTION EVALUATIONS) AT WHICH
!                TEMPERATURE IS DECREASED
!     randseed = RANDOM SEED FOR THERMAL FLUCTUATIONS
!     N.B. P, STEP AND VAR (IF IQUAD = 1) MUST HAVE DIMENSION AT LEAST NOP
!          IN THE CALLING PROGRAM.

!*****************************************************************************

INTEGER, INTENT(IN)        :: nop, maxfn, iprint, nloop, iquad
INTEGER, INTENT(OUT)       :: ifault
REAL(8), INTENT(IN)      :: stopcr, simp
REAL(8), INTENT(IN OUT)  :: p(:), step(:)
REAL(8), INTENT(OUT)     :: var(:), func
INTEGER, INTENT(IN)      :: myrank,nprocs
EXTERNAL functn
EXTERNAL functn2
REAL(8), INTENT(IN)      :: T0,Tstep
INTEGER, INTENT(IN)      :: Tfreq,randseed
!INTERFACE
!  SUBROUTINE functn(p, func)
!    IMPLICIT NONE
!    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
!    REAL, INTENT(IN)  :: p(:)
!    REAL, INTENT(OUT) :: func
!  END SUBROUTINE functn
!  SUBROUTINE functn2(p, func)
!    IMPLICIT NONE
!    REAL, INTENT(IN)  :: p(:)
!    REAL, INTENT(IN) :: func
!  END SUBROUTINE functn2
!END INTERFACE


!     Local variables
REAL(8)   :: g(nop+1,nop), h(nop+1), pbar(nop), pstar(nop), pstst(nop), &
               aval(nop), pmin(nop), temp(nop), bmat(nop*(nop+1)/2),  &
               vc(nop*(nop+1)/2), ymin, rmax, hstst, a0, hmin, test, hmean, &
               hstd, hstar, hmax, savemn, savehstd

REAL(8), PARAMETER :: zero = 0._dp, half = 0.5_dp, one = 1._dp, two = 2._dp
INTEGER     :: i, i1, i2, iflag, ii, ij, imax, imin, irank, irow, j, j1, jj, &
               k, l, loop, nap, neval, nmore, np1, nullty

!     A = REFLECTION COEFFICIENT, B = CONTRACTION COEFFICIENT, AND
!     C = EXPANSION COEFFICIENT.

REAL(8), PARAMETER :: a = 1._dp, b = 0.5_dp, c = 2._dp

!     SET LOUT = LOGICAL UNIT NO. FOR OUTPUT

INTEGER, PARAMETER :: lout = 6

! Delcaraltions related to parallelization
INTEGER :: nrounds,nround,nevalp,myimax
REAL(8) :: htilde
REAL(8), DIMENSION(nop) :: ptilde
INTEGER, PARAMETER :: maxprocs=200  !ahu s15 changed from 48
INTEGER, DIMENSION(maxprocs) :: mycase
REAL(8), DIMENSION(maxprocs) :: myval
REAL(8), DIMENSION(maxprocs,nop) :: mypoint
INTEGER :: npoints
REAL(8), DIMENSION(nop*(nop-1)/2,nop) :: points
REAL(8), DIMENSION(nop*(nop-1)/2) :: values
INTEGER, DIMENSION(1) :: best1
! additional variables for simulated annealing: current temperature plus thermal fluctuations
REAL(8)   :: SAtemp, htherm(nop+1), thermsimp(nop+1),thermstar(nop),thermstst(nop),SArandom(2*ABS(maxfn),3*nop)
INTEGER   :: tstepnext,r,SAnrand

integer :: realrank,mpierr,mpistat(MPI_STATUS_SIZE)

real(dp) :: timing(2)

CALL MPI_Comm_rank(MPI_COMM_WORLD,realrank,mpierr)
!CALL MPI_Comm_rank(comm,i,mpierr)
!print*, "Here is ME AT COMM ", myrank,i !,myrank_world
!realrank=myid
if (realrank==0) then 
	OPEN(UNIT=lout,FILE='opt.txt',STATUS='replace')
!	write(lout,'("realrank,mygroup,myrank,mygrank,iprint:",5I4)') realrank,mygroup,myrank,mygrank,iprint
    write(lout,'("nop,maxfn,iprint,nloop,iquad: ",5I4)') nop, maxfn, iprint, nloop, iquad
    write(lout,'("stopcr,simp: ",2g14.6)') stopcr,simp
    write(lout,'("T0,Tstep: ",2g14.6)') T0,Tstep
    write(lout,'("Tfreq: ",I4)') Tfreq
!	write(*,'("realrank,mygroup,myrank,mygrank,iprint:",5I4)') realrank,mygroup,myrank,mygrank,iprint
    write(*,'("nop,maxfn,iprint,nloop,iquad: ",5I4)') nop, maxfn, iprint, nloop, iquad
    write(*,'("stopcr,simp: ",2g14.6)') stopcr,simp
    write(*,'("T0,Tstep: ",2g14.6)') T0,Tstep
    write(*,'("Tfreq: ",I4)') Tfreq
end if 
!write(*,'("realrank,mygroup,myrank,mygrank,iprint:",5I4)') realrank,mygroup,myrank,mygrank,iprint



!     CHECK INPUT ARGUMENTS

ifault = 0
IF (nop <= 0) ifault = 3
IF (nloop <= 0) ifault = 4
!IF ((realrank<nprocs).AND.(realrank/=myrank)) ifault = 5
IF (ifault /= 0) RETURN

!     SET NAP = NO. OF PARAMETERS TO BE VARIED, I.E. WITH STEP /= 0

nap = COUNT(step /= zero)
neval = 0
nevalp= 0
loop = 0
iflag = 0

!     IF NAP = 0 EVALUATE FUNCTION AT THE STARTING POINT AND RETURN

IF (nap <= 0) THEN
  CALL functn(p,func)
  RETURN
END IF


! IF NAP<=NPROCS, RETURN
IF (nap <=nprocs) THEN
  CALL functn(p,func)
  RETURN
ENDIF

! If printing, print degree of parallelization
IF ((iprint>=0).AND.(realrank==0))  WRITE(lout,4900) nap,nprocs
!     IF PROGRESS REPORTS HAVE BEEN REQUESTED, PRINT HEADING
IF ((iprint > 0).AND.(realrank==0)) WRITE (lout,5000) iprint

!     SET UP THE INITIAL SIMPLEX
		  !IF (iprint > 0) THEN
		!	WRITE (lout,'("Set up initial simplex ")') 
		  !END IF
20 g(1,:) = p
irow = 2
DO i = 1, nop
  IF (step(i) /= zero) THEN
    g(irow,:) = p
    g(irow,i) = p(i) + step(i)
    irow = irow + 1
  END IF
END DO

np1 = nap + 1

!SUBROUTINE function_distribute(npoints,minrank,maxrank,points,vals,nevalps,func)
!print*, "realrank,nprocs,myrank ", realrank,nprocs,myrank 	!ahu april13
CALL function_distribute(np1,0,nprocs-1,g(1:np1,:),h(1:np1),nevalp,functn,myrank) 
neval=np1
IF ((iprint > 0).AND.(realrank==0)) THEN
  DO i=1,np1
    WRITE (lout,5100) i, h(i), g(i,:)
  ENDDO
  write(lout,*) "finished setting up initial simplex"
  write(lout,*) "count number of function evaluations which here is neval=np1", neval,np1  
END IF

! Initialize simulated annealing
SAtemp=T0 ! set initial temperature
tstepnext=neval+Tfreq ! next time will lower temperature
SAnrand=1 ! which set of random numbers up to
! set seed of RANDOM_NUMBER function and draw random numbers
CALL RANDOM_SEED (SIZE=r)
CALL RANDOM_SEED (PUT = randseed*(/(k,k=1,r)/))
CALL RANDOM_NUMBER(SArandom)

!     START OF MAIN CYCLE.

!     FIND MAX. & MIN. VALUES FOR CURRENT SIMPLEX (HMAX & HMIN).

Main_loop: DO
  loop = loop + 1
	!if (realrank==0.or.realrank==1) PRINT*, "HERE I AM starting loop", MYRANK,MYRANK_WORLD,loop
    IF ((iprint > 0).AND.(realrank==0)) THEN
        WRITE(LOUT,*) 
        WRITE(LOUT,'("beginning of Main_loop, just augmented loop: loop=loop+1")') 
        write(LOUT,'("T0,tstep,tfreq (these do not change): ",2g14.6)') T0,tstep,tfreq
        write(LOUT,'("current SAtemp,tstepnext are: ",2g14.6)') SAtemp,tstepnext
        write(LOUT,'("loop is: ",I8)') loop
        write(LOUT,'("neval is: ",I8)') neval
    END IF 

  IF ((iprint > 0).AND.(realrank==0)) THEN
    WRITE(LOUT,*) 
    WRITE(LOUT,*) 
    WRITE(LOUT,'("Should I lower temperature?")') 
    WRITE(LOUT,'("check: neval>tstepnext?")') 
    write(LOUT,'("---if yes: lower temp and get the new temp as SAtemp=SAtemp*Tstep and new tstepnext as tstepnext=tstepnext+Tfreq : ")')
    write(LOUT,'("---if no:  leave SAtemp and tstepnext as they are  lower temp and get the new temp as SAtemp=SAtemp*Tstep: ")') 
    WRITE(LOUT,*) 
    WRITE(LOUT,*) 
  END IF 

  If (neval>tstepnext) THEN ! lower temperature
      SAtemp=SAtemp*Tstep
      tstepnext=tstepnext+Tfreq ! next time will lower temperature
      IF ((iprint > 0).AND.(realrank==0)) THEN
        WRITE(LOUT,*) 
        WRITE(LOUT,*) 
        WRITE(LOUT,'("neval>tstepnex so lower away!t")') 
        write(LOUT,'("---> lower temp and get the new temp as SAtemp=SAtemp*Tstep and new tstepnext as tstepnext=tstepnext+Tfreq : ")') 
        write(LOUT,'("     new SAtemp,tstepnext are : ",2g14.6)') SAtemp,tstepnext
        WRITE(LOUT,*) 
      END IF 
        !if (realrank==0) write(*,'("lowering temperature now ",2I4,g14.6,I4)') realrank,loop,SAtemp,tstepnext
  ENDIF

  IF ((iprint > 0).AND.(realrank==0)) THEN
    WRITE(LOUT,*) "After checking wtr I should lower temperature i.e. out of the neval>tstepnext if statement "
    write(LOUT,'("next time, will lower temprature when neval>tstepnext: ",I8)') tstepnext          
    write(LOUT,'("neval,tstepnext: ",2I8)') neval,tstepnext              
  END IF 

  IF ((iprint > 0).AND.(realrank==0)) THEN
    WRITE(LOUT,*) 
    write(LOUT,'("compute and add thermal fluctuations: thermsimp,thermstar,thermstst ")')       
    WRITE(LOUT,*) 
  END IF 

  ! POSITIVE THERMAL FLUCTUATIONS
  thermsimp=-1*SAtemp*LOG(SArandom(SAnrand,1:nop+1))
  thermstar(1:nprocs)=-1*SAtemp*LOG(SArandom(SAnrand,nop+2:nop+nprocs+1))
  thermstst(1:nprocs)=-1*SAtemp*LOG(SArandom(SAnrand,2*nop+2:2*nop+nprocs+1))
    IF ((iprint > 0).AND.(realrank==0)) THEN
        write(lout,*) 
        WRITE(LOUT,'(tr3,"i",tr6,"SArandom",tr3,"logSArandom",tr4,"SAtemp*log",tr5,"thermsimp",tr6,"h",tr7,"htherm=h+thermsimp")')
        do i=1,np1
            WRITE(LOUT,'(I4,6g14.6)') i,SArandom(SAnrand,i),LOG(SArandom(SAnrand,i)),SAtemp*LOG(SArandom(SAnrand,i)),thermsimp(i),h(i),h(i)+thermsimp(i)
        end do 
        write(*,'("loop,thermsimp(1:2) ",I8,2g14.6)') loop,thermsimp(1:2)
    END IF 
  SAnrand=SAnrand+1
  IF (SAnrand>2*ABS(maxfn)) THEN ! used up all random numbers- don't think this should happen
    SAnrand=1
  ENDIF
  !thermsimp=0.0_dp !ahu 112213
  !thermstar=0.0_dp
  !thermstst=0.0_dp
  
! SORT SIMPLEX SO VALUES OF h (with thermal fluctuations)  ARE INCREASING IN INDEX. SORT g AND h
  IF ((iprint > 0).AND.(realrank==0)) THEN
    write (lout,*)
    WRITE (lout,*) "sort simplex so values of h (with thermal fluctuations) are increasing in index. sort g and h acordingly"
    write (lout,*)
  END IF 
    
  htherm(1:np1)=h(1:np1)+thermsimp(1:np1)
  CALL QsortC3(htherm(1:np1),g(1:np1,:))
  htherm(1:np1)=h(1:np1)+thermsimp(1:np1)
  CALL QsortC2(htherm(1:np1),h(1:np1))

    !IF ((iprint > 0).AND.(realrank==0)) THEN
    !    write(lout,*) 
    !    write(lout,*) 
    !    WRITE(LOUT,'(tr3,"i",tr13,"h",tr8,"htherm")')
    !    do i=1,np1
    !        WRITE(LOUT,'(I4,2g14.6)') i,h(i),htherm(i)
    !    end do 
    !END IF 
  
  
  imax=np1
  hmax=h(np1)
  imin=1
  hmin=h(1)

      !FIND THE CENTROID OF THE VERTICES OF BEST np1 MINUS nprocs POINTS
  IF ((iprint > 0).AND.(realrank==0)) THEN
    write (lout,*)
    WRITE (lout,*) "calculate pbar [the centroid fo the vertices of best np1 MINUS nprocs points]"
    WRITE (lout,*) "do i=1,np1-nprocs .... pbar = pbar + g(i,:)"
    write (lout,*)
  END IF 
    
  pbar = zero
  DO i = 1, np1-nprocs
    pbar = pbar + g(i,:)
  END DO
  pbar = pbar / (np1-nprocs)
  !IF ((iprint > 0).AND.(realrank==0)) THEN
	!WRITE (lout,*) 
	!WRITE (lout,*) 
	!WRITE (lout,*) 
	!WRITE (lout,'("Position of centroid PBAR (midway between all points other than imax) is calculated. ")') 
	!WRITE (lout,*) 
	!write(lout,'(TR2,"maxloc",TR2,"minloc",TR4,"maxval",TR4,"minval",TR2,"thermsimp")')
	!WRITE (lout,'(2(4x,I4),3F10.2)') maxloc(htherm(1:np1)),minloc(htherm(1:np1)),maxval(htherm(1:np1)),minval(htherm(1:np1)),thermsimp(1)
	!WRITE (lout,*) 
	!WRITE (lout,'("Here is G(IMAX,:) ")') 
	!WRITE (lout,5100) neval, hmax,g(imax,:)
	!WRITE (lout,*) 
	!WRITE (lout,'("Here is the centroid PBAR ")') 
	!WRITE (lout,5100) neval, 0.0_dp, pbar
	!WRITE (lout,*) 
  !END IF

! ASSIGN EACH PROCESSOR ONE OF THE NPROCS WORST POINTS
  myimax=np1-myrank
  IF ((iprint > 0).AND.(realrank==0)) THEN
    write (lout,*)
    WRITE (lout,*) "assign each processor one of the nprocs worst points"
    WRITE (lout,*) "myimax=np1-myrank"
    WRITE (lout,*) "for example: if numprocs is 2 so that myrank is 0 or 1, we will have myimax=np1 for myrank=0 and myimax=np1-1 for myrank=1. "
    WRITE (lout,*) "then reflect the worst point [g(myimax,:) ] through pbar and call this pstar and get obj value at pstar and call it hstar "
    write (lout,*)
  END IF 
!For example: if numprocs is 2 so that myrank is 0 or 1, we will 
!have myimax=np1 for myrank=0 and myimax=np1-1 for myrank=1. 
!Thes are the worst and the second worst points in the current simplex values with thermal fluctuations
!   A reflection of the worst-response point W is performed and the response RR of the reflected point R is evaluated. 
!	Reflect the worst (i.e.G(MYIMAX,:)) point through PBAR
!	Call this reflection point PSTAR and evaluate the obj function at PSTAR
!	HSTAR = function value at PSTAR
 pstar = a * (pbar - g(myimax,:)) + pbar
 timing(1)=MPI_WTIME() !ahu 0317
 CALL functn(pstar,hstar)
 timing(2)=MPI_WTIME() !ahu 0317
 if (realrank==0) WRITE (lout,'("Just calling func ",2I4,F14.2)') realrank,myrank,timing(2)-timing(1) !ahu 0317
 IF ((iprint > 0).AND.(realrank==0)) THEN
	WRITE (lout,*) 
	WRITE (lout,*) 
	WRITE (lout,'("reflection --> reflecting MY worst point i.e. g( myimax,:) through PBAR to get the reflection point PSTAR. ")') 
	WRITE (lout,'("(R) PSTAR = a * (PBAR - g(imax,:)) + PBAR ")') 
	WRITE (lout,'("(RR) HSTAR=functn(PSTAR,HSTAR) ")') 
	WRITE (lout,*) 
	WRITE (lout,'("now each processor follows different instructions based on their hstar  ")') 
	WRITE (lout,'("and they return a point mypoint with value myval in appropriate index corresponding to which proccessor. ")') 
	WRITE (lout,'("record what happened in appropriate engry of mycase ")') 
	WRITE (lout,*) 
END IF

	! NOW EACH PROCESSOR DOES THE FOLLOWING SEPARATELY:
	! Each processor follows instrcutions based on hstar, returns a point MYPOINT with value MYVAL
	! in approproate index corresponding to which processor. Record what happened in entry of MYCASE

	! case 1: hstar < hmin 
	! reflection did well i.e. hstar is beter than the best point hmin 
	! ---> expand further i.e. expand pbar in the direction of pstar
	! ---> g and h at expansion pt: PSTST,HSTST
	IF (hstar-thermstar(myrank+1) < htherm(1)) THEN ! case 1 - return the better of the expansion point and the reflection point
		mycase(myrank+1)=1
		pstst = c * (pstar - pbar) + pstar
		CALL functn(pstst,hstst)
		IF (hstst<hstar) THEN ! return expansion point (pstst)
			mypoint(myrank+1,:)=pstst
			myval(myrank+1)=hstst
		ELSE ! return reflection point (pstar)
			mypoint(myrank+1,:) = pstar
			myval(myrank+1) = hstar
		END IF
		
	! case 2: hstar > hmin  and hstar<htherm(np1-1) or (np1-2) ... (np1-myrank-1) .. etc. 	
	! reflection did ok i.e. hstar not beter than the best point, BUT better than the second (or third etc. depending on myrank) worst point 
	! ---> do nothing and just return the original reflection point pstar.  
	! ---> Note that this is the only case in which the processor does not do an extra function call
	ELSEIF (hstar-thermstar(myrank+1)<htherm(np1-myrank-1)) THEN ! case 2 - return the relction point
		mycase(myrank+1)=2
		mypoint(myrank+1,:)=pstar
		myval(myrank+1) = hstar

	! case 3 or 4: hstar > hmin AND hstar > hterm(..) 
	! reflection did not do any good i.e. hstar not better than the best AND not better than any of the worst points
	! ---> contract by getting the convex combo of pstar and pbar
	! ---> contraction point: PSTST
	! ---> fnctn value at contraction point: HSTST
	ELSE 
		! calculate contraction point, pstst
		pstst= b * pstar + (one-b) * pbar
		CALL functn(pstst,hstst)
		! compare the contraction point to the better of pstar and the original point of the simplex
		IF (hstar-thermstar(myrank+1)<htherm(myimax)) THEN ! pstar is better of the two
			IF (hstst<hstar) THEN ! case 3 - return the contraction point
			mycase(myrank+1)=3
			mypoint(myrank+1,:)=pstst
			myval(myrank+1) = hstst
			ELSE ! case 4 - return the better of pstar and the original point, which is here pstar
			mycase(myrank+1)=4
			mypoint(myrank+1,:)=pstar
			myval(myrank+1) = hstar
			ENDIF
		ELSE ! ther original point of the simplex is better than pstar
			IF (hstst-thermstst(myrank+1)<htherm(myimax)) THEN ! case 3 - return the contraction point
			mycase(myrank+1)=3
			mypoint(myrank+1,:)=pstst
			myval(myrank+1) = hstst
			ELSE ! case 4 - return the better of pstar and the original point of the simplex,here g(myimax,:)
			mycase(myrank+1)=4
			mypoint(myrank+1,:)=g(myimax,:)
			myval(myrank+1) = h(myimax)
			ENDIF
		ENDIF 
	ENDIF

! SHARE OUTCOMES AND THEN FIGURE OUT WHAT THE NEW SIMPLEX WILL LOOK LIKE  
! First share outcomes MYCASE, MYPOINT and MYVAL, count number of function evaluations, print as appropriate
  IF ((iprint > 0).AND.(realrank==0)) THEN
  write(lout,*) "share outcomes and figure out what the new simplex will look like"
  write(lout,*) "First share outcomes MYCASE, MYPOINT and MYVAL"
  write(lout,*) "right after mycase,mypoint,myval assignments: count number of function evaluations, print as appropriate"
  write(lout,*) "nrounds is no of evals per processor (it can be >1 depending on mycase"
  write(lout,*) "nevalp=nevalp+nrounds"
  END IF 
  nrounds=1 ! the number of evals per processor
  DO ii=0,nprocs-1
    CALL MPI_BCAST(mycase(ii+1),1,MPI_INTEGER,ii,MPI_COMM_WORLD,mpierr)
    CALL MPI_BCAST(mypoint(ii+1,:),nop,MPI_REAL8,ii,MPI_COMM_WORLD,mpierr)
    CALL MPI_BCAST(myval(ii+1),1,MPI_REAL8,ii,MPI_COMM_WORLD,mpierr)
    IF (mycase(ii+1) == 2) THEN ! only did one function evaluation at point MYPOINT(ii+1)
      neval=neval+1
      IF ((iprint > 0).AND.(realrank==0)) THEN
        IF (MOD(neval,iprint) == 0) WRITE (lout,5100) neval, myval(ii+1), mypoint(ii+1,:)
      END IF
    ELSE ! two function evaluations
      IF ((ii>0) .AND. (iprint>0)) THEN
         ! Send pstar, hstar, pstst,hstst to master for printing. If ii=0, this is master, have values already.
         IF (realrank==ii) THEN
           CALL MPI_SEND(pstar,nop, MPI_REAL8,0,1,MPI_COMM_WORLD,mpierr)
           CALL MPI_SEND(hstar,1, MPI_REAL8,0,2, MPI_COMM_WORLD,mpierr)
           CALL MPI_SEND(pstst,nop, MPI_REAL8,0,3, MPI_COMM_WORLD,mpierr)
           CALL MPI_SEND(hstst,1, MPI_REAL8,0,4, MPI_COMM_WORLD,mpierr)
         ENDIF
         IF (realrank==0) THEN
           CALL MPI_RECV(pstar,nop, MPI_REAL8,ii,1,MPI_COMM_WORLD, mpistat,mpierr)
           CALL MPI_RECV(hstar,1, MPI_REAL8,ii,2,MPI_COMM_WORLD, mpistat,mpierr)
           CALL MPI_RECV(pstst,nop, MPI_REAL8,ii,3,MPI_COMM_WORLD, mpistat,mpierr)
           CALL MPI_RECV(hstst,1, MPI_REAL8,ii,4,MPI_COMM_WORLD, mpistat,mpierr)
         ENDIF
      ENDIF
      neval=neval+1 
      IF ((iprint > 0).AND.(realrank==0)) THEN
        IF (MOD(neval,iprint) == 0) WRITE (lout,5100) neval, hstar, pstar
      ENDIF
      neval=neval+1 
      IF ((iprint > 0).AND.(realrank==0)) THEN
        IF (MOD(neval,iprint) == 0) WRITE (lout,5100) neval, hstst, pstst
      ENDIF 
      nrounds=2 ! at least one processor did two evaluations rather than one, update nevalp by 2
    ENDIF
  ENDDO
  nevalp=nevalp+nrounds
  
	
	if (iprint>0.and.realrank==0) then
    write(lout,*) "do 0 to nprocs and writes down what their mycase situation is (after broadcast)"
    write(lout,*) "For case=1 where the reflection PSTAR gives a HSTAR lower than current min, we do  expansion to calculate a new PSTST through expanding and then compare that PSTST to PSTAR"
    write(lout,*) "For case=2 where the reflection PSTAR gives a HSTAR higher than current min BUT lower than htherm(myimax-1), we do nothing (no new PSTST) "
    write(lout,*) "For case=3, where the reflection PSTAR gives a HSTAR higher than both the current min AND htherm(myimax-1)"
    write(lout,*) "............we do contraction to calculate a new PSTST by convx combo of pstar and pbar "
    write(lout,*) "............then we compare PSTST to EITHER pstar OR the worst point of the original simplex (depending on which is better):  "
    write(lout,*)
    write(lout,*)
		do ii=0,nprocs
        write(lout,'(TR2,"myrank",TR2,"mycase",TR5,"myval")')
        write(lout,'(2(4x,I4),F10.2)') myrank,mycase(ii+1),myval(ii+1)  
      if (mycase(ii+1)==1) then
        write(lout,*) "mycase(ii+1)=1: hstar-thermstar(myrank+1) < htherm(1) where htherm(1) is the minimum"
        write(lout,*) "nonparal: hstar < hmin. reflection did well i.e. hstar is lower than the current min (htherm(1))"
        write(lout,*) " ---> expand further i.e. expand pbar in the direction of pstar --> call it pstst"
        write(lout,*) " ---> call functn to get h at expansion pt pstst and call it hstst: PSTST,HSTST"
        write(lout,*) " -------------> if hstst<hstar then return expansion point (pstst) i.e. mypoint(myrank+1,:)=pstst,myval(myrank+1)=hstst"
        write(lout,*) " -------------> if hstst>=hstar then return reflection point (pstar) i.e. mypoint(myrank+1,:)=pstar,myval(myrank+1)=hstar"  
        write(lout,*)   
      else if (mycase(ii+1)==2) then
        write(lout,*) "mycase(ii+1)=2: hstar-thermstar(myrank+1) > htherm(1) BUT  hstar-thermstar(myrank+1)<htherm(np1-myrank-1) i.e hstar<htherm(np1-1) or (np1-2) ... (np1-myrank-1) .. etc"
        write(lout,*) "mycase(ii+1)=2: hstar > hmi but lower than htherm(myimax-1) "
				write(lout,*) "nonparal: hstar is not <hmin. but lower than some point, replace p(imax) with pstar and hmax with hstar"
				write(lout,*) "nonparal: note that in nonparallel version it's different, it doesn't sort the objvalue array and it compares hstar to aNY point other than imax"
				write(lout,*) " ---> do nothing and just return the reflection point pstar i.e. mycase(..)=2,mypoint(myrank+1,:)=pstar,myval(myrank+1)=hstar."
				write(lout,*) " ---> Note that this is the only case in which the processor does not do an extra function call"
			else if (mycase(ii+1)>2) then
        write(lout,*) "mycase(ii+1)=3 or 4: hstar-thermstar(myrank+1) > htherm(1) AND  hstar-thermstar(myrank+1)>htherm(np1-myrank-1)"
				write(lout,*) "nonparal: hstar > all func values except possibly hmax.  "
				write(lout,*) "nonparal: if hstar<=hmax, replace p(imax) with pstar/ hmax with hstar test whether it is < func value at some point other than p(imax)"
        write(lout,*) "reflection did not do any good i.e. hstar not better than the best AND not better than MY worst point (after imax)"
				write(lout,*) "---> contract by getting the convex combo of pstar and pbar: PSTST= b * pstar + (one-b) * pbar"
				write(lout,*) "---> compare the contraction point to EITHER pstar OR the worst point of the original simplex (depending on which is better):  "
				write(lout,*) "---> IF (hstar-thermstar(myrank+1)<htherm(myimax)) THEN i.e. (pstar is better of the two)"
        write(lout,*) "----------------> IF (hstst<hstar) THEN ! case 3 - return the contraction point"
        write(lout,*) "----------------------------> mycase(myrank+1)=3, mypoint..pstst, myval..hstst"
        write(lout,*) "----------------> ELSE - return the better of pstar and the original point, which is here pstar"
        write(lout,*) "----------------------------> mycase(myrank+1)=4, mypoint..pstar, myval..hstar"
				write(lout,*) "---> ELSE  ther original point of the simplex is better than pstar "
        write(lout,*) "----------------> IF (hstst-thermstst(myrank+1)<htherm(myimax)) THEN ! case 3 - return the contraction point"
        write(lout,*) "----------------------------> mycase(myrank+1)=3, mypoint..pstst, myval..hstst"
        write(lout,*) "---------------->ELSE ! case 4 - return the better of pstar and the original point of the simplex,here g(myimax,:)"
        write(lout,*) "----------------------------> mycase(myrank+1)=4, mypoint..g(myimax,:), myval..h(myimax)"
      end if   
			write(lout,*) 
			write(lout,*) 
		end do   
	end if 
	
	if (iprint>0.and.realrank==0) then
    write(lout,*) 
    write(lout,*) 
		write(lout,*) "CHECK cases"
		write(lout,*) "-----IF at least one proces did not have case 4, then replace worst points and restart loop"   
		write(lout,*) "-----IF all processes had case 4, then shrink"
    write(LOUT,'("current SAtemp,tstepnext are: ",2g14.6)') SAtemp,tstepnext
    write(LOUT,'("loop is: ",I8)') loop
    write(LOUT,'("neval is: ",I8)') neval
    write(lout,*) "............"
    write(lout,*) "............"
	end if 

  ! if at least one processor did not have case 4, replace the worst nprocs points of simplex
  if (MINVAL(mycase(1:nprocs))<4) then ! use points and values from each processor
    if (iprint>0.and.realrank==0) then
      write(lout,*) "............"
      write(lout,*) "............"
      write(lout,*) "At least one process did not have case 4 ---> replace the worst nprocs points of simplex with mypoint(1:nprocs) and restart loop"
      write(lout,*) "-----> g(np1-nprocs+1:np1,:)=mypoint(1:nprocs,:)"   
      write(lout,*) "-----> h(np1-nprocs+1:np1)=myval(1:nprocs)"
      write(lout,*) "-----> go to 250 (which starts main cycle if loop<nloop)"
      write(LOUT,'("        loop,nloop are: ",2I8)') loop,nloop
      write(lout,*) "right before GO TO 250"
      write(lout,*) 
      write(lout,*) 
    end if 
    g(np1-nprocs+1:np1,:)=mypoint(1:nprocs,:)
    h(np1-nprocs+1:np1)=myval(1:nprocs)
    GO TO 250
  endif

!  Otherwise, all the processors had case 4
!	---> shrink the simplex by replacing each point other than the current minimum 
!	     by a point mid-way between its current position and the minimum"
	if (iprint>0.and.realrank==0) then
    write(lout,*) "............"
    write(lout,*) "............"
		write(lout,*) "All the processors had case 4"
		write(lout,*) "---> shrink the simplex by replacing each point other than the current minimum"
		write(lout,*) "     by a point mid-way between its current position and the minimum"
    write(lout,*) 
    write(lout,*) 
  end if 
		
  DO i = 2, np1
      DO j = 1, nop
        IF (step(j) /= zero) g(i,j) = (g(i,j) + g(1,j)) * half
      END DO
  ENDDO

 if (iprint>0.and.realrank==0) then
  write(lout,*) 
  write(lout,*) 
	write(lout,*) "after shrink"
  write(lout,'("before func_dist: nevalp,nrounds (this is a value left from share outcome loop so 1 or 2):",2I8)') nevalp,nrounds
	write(lout,*) "note that before this the last time nrounds was assigned a value was in the share outcome loop (1 or 2)"
  write(lout,*) "............"
  write(lout,*) "............"
end if 
  !SUBROUTINE function_distribute(npoints,minrank,maxrank,points,vals,nevalps,func)
  CALL function_distribute(nap,0,nprocs-1,g(2:np1,:),h(2:np1),nrounds,functn,myrank)
  ! increase count of number of evals per processor
  nevalp=nevalp+nrounds
  if (iprint>0.and.realrank==0) then
    write(lout,*) 
    write(lout,*) 
    write(lout,*) "after func_dist"
  	write(lout,*) "increase count of number of evals per procesor (nevalp=nevalp+nrounds) "
    write(lout,'("after func_dist: nevalp,nrounds:",2I8)') nevalp,nrounds
    write(lout,*) "note that before this the last time nrounds was assigned a value was in the share outcome loop (1 or 2)"
    write(lout,*) 
    write(lout,*) 
  end if 


if (iprint>0.and.realrank==0) then
	write(lout,*) 
	write(lout,*) 
  write(lout,*) "starting do loop for updating neval=neval+1 within do i=1,np1 which has some printing (now commented out)"
end if 
  ! count evals and print 
  DO i=2,np1
      neval = neval + 1
      IF ((iprint > 0).AND.(realrank==0)) THEN
        !IF (MOD(neval,iprint) == 0) WRITE (lout,5100) neval, h(i), g(i,:)
      END IF
  END DO
 if (iprint>0.and.realrank==0) then
	write(lout,*) 
  write(lout,'("after updating neval in do loop, neval:",I8)') neval 
 	write(lout,*) 
	write(lout,*) 
 end if 	 

!     IF LOOP = NLOOP TEST FOR CONVERGENCE, OTHERWISE REPEAT MAIN CYCLE.
 if (iprint>0.and.realrank==0) then
	write(lout,*) 
	write(lout,*) 
  write(lout,'("loop,nloop:",2I8)') loop,nloop
  write(lout,*) "CHECK loop<?nloop"
  write(lout,'("-----IF loop<nloop, repeat main cycle")')
  write(lout,'("-----IF loop=nloop, calc hmean&stdev, write best of h so far and test for convergence (there are other conditions so see code) ")') 
  write(lout,'("--------------------------------------------------IF not convgd cond then repeat main cycle")')
 end if 
 if (iprint>0.and.realrank==0) then
  write(lout,*) 
  write(lout,*) "right before: 250 IF (loop < nloop) CYCLE Main_loop "
	write(lout,*) 
end if 
  250 IF (loop < nloop) CYCLE Main_loop

  IF (iprint>0.and.realrank==0) THEN	
    write(lout,*) 
    write(lout,*) "right after: 250 IF (loop < nloop) CYCLE Main_loop "
  	write(lout,'("loop,nloop: ",2I8)') loop,nloop
    write(lout,'("loop was equal to nloop so did not start main cycle again, now will do the following: ")') 
  	write(lout,*) 
    write(lout,*) 
  ENDIF

  IF (iprint>0.and.realrank==0) THEN	
    write(lout,*) 
    write(lout,*) 
    write(lout,'("calculate mean & stdev of func values for current simplex (hmean,hstd)")') 
    write(lout,*) "write hstd"
  ENDIF

!     CALCULATE MEAN & STANDARD DEVIATION OF FUNCTION VALUES FOR THE
!     CURRENT SIMPLEX.
  hmean = SUM( h(1:np1) ) / np1
  hstd = SUM( (h(1:np1) - hmean) ** 2 )
  hstd = SQRT(hstd / np1)

	IF ((iprint > 0).AND.(realrank==0)) THEN
		WRITE (lout,5300) hstd
	END IF
  IF (iprint>0.and.realrank==0) THEN	
    write(lout,*) 
    write(lout,'("hmean,hstd : ",2g14.6)') hmean,hstd
    write(lout,*) 
  ENDIF
  IF (iprint>0.and.realrank==0) THEN	
    write(lout,*) 
    write(lout,'("writing best so far by calling functn2 ")') 
    write(lout,*) "but just writing params not moments. see the change in functn2"
    best1=MINLOC(h(1:np1))
    CALL functn2(g(best1(1),:),h(best1(1)),neval,hmean,hstd)
    write(lout,'("hmean,hstd : ",2g14.6)') hmean,hstd
    write(lout,'("h(1),h(best1),h(np1) : ",3g14.6)') h(1),h(best1(1)), h(np1)
    write(lout,'("htherm(1),htherm(best1),htherm(np1) : ",3g14.6)') htherm(1),htherm(best1(1)), htherm(np1)
  ENDIF

  IF (iprint>0.and.realrank==0) THEN	
    write(lout,*) 
    write(lout,*) 
    write(lout,*) "CHECK hstd>?stopcr"
    write(lout,'("---------> IF hstd>stopcr (and maxfn,neval.nevalp conditions as well) then set iflag and loop to zero and go to the start of the main cycle again ")') 
    write(lout,'("---------> IF hstd<=stopcr then find the centroid of the current simplex and get the function value there ")') 
    write(lout,*) "............"
    write(lout,*) "............"
  END IF
!     IF THE RMS > STOPCR, SET IFLAG & LOOP TO ZERO AND GO TO THE
!     START OF THE MAIN CYCLE AGAIN.
  IF (hstd > stopcr .AND. (((maxfn>=0).AND.(neval <= maxfn)).OR.((maxfn<0).AND.(nevalp <= -1*maxfn))  )) THEN
    iflag = 0
    loop = 0
    IF (iprint>0.and.realrank==0) THEN	
      write(lout,*) "............"
      write(lout,*) "............"
      write(lout,*) " (hstd > stopcr .AND. (((maxfn>=0).AND.(neval <= maxfn)).OR.((maxfn<0).AND.(nevalp <= -1*maxfn))  )) " 
      write(lout,*) " so set iflag=0 and loop=0 and start the main cycle again " 
      write(lout,*) " right before CYCLE MAIN_LOOP " 
      write(lout,*) 
      write(lout,*)   
    ENDIF
    CYCLE Main_loop
  END IF

  IF (iprint>0.and.realrank==0) THEN	
    write(lout,*) "............"
    write(lout,*) "............"
    write(lout,*) " .NOT. (hstd > stopcr .AND. (((maxfn>=0).AND.(neval <= maxfn)).OR.((maxfn<0).AND.(nevalp <= -1*maxfn))  )) " 
    write(lout,*) " so find the centroid of current simplex and get teh function value there " 
    write(lout,*) " call functn(p,func) "
    write(lout,*) " neval = neval + 1"
    write(lout,*) " nevalp = nevalp + 1"  
    write(lout,*) " IF ((iprint > 0).AND.(realrank==0)) THEN "
    write(lout,*) "    IF (MOD(neval,iprint) == 0) WRITE (lout,5100) neval, func, p "
    write(lout,*) " END IF  "
    write(lout,*) " CHECK wtr the no. of func values allowed, maxfn, has been overrun. if so, exit with ifault= 1 " 
    write(lout,*) " IF maxfn>0 check based on neval. if maxfn<0 check based on nevalp. "
    write(lout,*) " i.e.: ---->IF (((maxfn>=0).AND.(neval > maxfn)).OR.((maxfn<0).AND.(nevalp > -1*maxfn))) THEN"
    write(lout,*) "       ---->IF the above if statement true then something ... RETURN "
    write(lout,*) "       ---->IF the above if statement not true then you get out of that if statement ... CONVGENCE CRTI SATISFIED - RETURN "
    write(lout,*) "............"
    write(lout,*) "............"
  ENDIF
!     FIND THE CENTROID OF THE CURRENT SIMPLEX AND THE FUNCTION VALUE THERE.
  DO i = 1, nop
    IF (step(i) /= zero) THEN
      p(i) = SUM( g(1:np1,i) ) / np1
    END IF
  END DO
  ! Note: all processors will do this eval
  CALL functn(p,func)
  neval = neval + 1
  nevalp= nevalp+1
  IF ((iprint > 0).AND.(realrank==0)) THEN
    IF (MOD(neval,iprint) == 0) WRITE (lout,5100) neval, func, p
  END IF
 
!     TEST WHETHER THE NO. OF FUNCTION VALUES ALLOWED, maxfn, HAS BEEN
!     OVERRUN; IF SO, EXIT WITH IFAULT = 1.
!  If MAXFN>0, check based on neval. If MAXFN<0, check based on nevalp.
  IF (((maxfn>=0).AND.(neval > maxfn)).OR.((maxfn<0).AND.(nevalp > -1*maxfn))) THEN
    ifault = 1
    IF (iprint < 0) RETURN
    IF ((iprint > 0).AND.(realrank==0)) THEN
      IF (maxfn>=0) THEN
        WRITE (lout,5200) maxfn
      ELSE
        WRITE (lout,5201) -1*maxfn
      ENDIF
      WRITE (lout,5300) hstd
      WRITE (lout,5400) p
      WRITE (lout,5500) func
    ENDIF
    RETURN
  END IF


!     CONVERGENCE CRITERION SATISFIED.
!     IF IFLAG = 0, SET IFLAG & SAVE HMEAN.
!     IF IFLAG = 1 & CHANGE IN HMEAN <= STOPCR THEN SEARCH IS COMPLETE.

  IF ((iprint >= 0).AND.(realrank==0)) THEN
    WRITE (lout,5600)
    WRITE (lout,5400) p
    WRITE (lout,5500) func
  END IF

  IF (iflag == 0 .OR. ABS(savemn-hmean) >= stopcr) THEN
    iflag = 1
    savemn = hmean
    loop = 0
  ELSE
    EXIT Main_loop
  END IF

END DO Main_loop

IF((iprint >= 0).AND.(realrank==0)) THEN
  WRITE (lout,5700) neval
  IF (maxfn<0) THEN
    WRITE (lout,5701) nevalp
  ENDIF
  WRITE (lout,5800) p
  WRITE (lout,5900) func
  !CLOSE(LOUT) ! ahu f14 ahu 092822 if you close this then you won't have the below (quadratic surface, second derivatives etc.)
END IF
IF (iquad <= 0) RETURN

!------------------------------------------------------------------

!     QUADRATIC SURFACE FITTING

IF ((iprint >= 0).AND.(realrank==0)) WRITE (lout,6000)

!     EXPAND THE FINAL SIMPLEX, IF NECESSARY, TO OVERCOME ROUNDING
!     ERRORS.
!     NOTE: This step does not take advantage of prallelization
hmin = func
nmore = 0
DO i = 1, np1
  DO
    test = ABS(h(i)-func)
    IF (test < simp) THEN
      DO j = 1, nop
        IF (step(j) /= zero) g(i,j) = (g(i,j)-p(j)) + g(i,j)
        pstst(j) = g(i,j)
      END DO
      CALL functn(pstst,h(i))
      nmore = nmore + 1
      neval = neval + 1
      nevalp = nevalp + 1
      IF (h(i) >= hmin) CYCLE
      hmin = h(i)
      IF ((iprint >= 0).AND.(realrank==0)) WRITE (lout,5100) neval, hmin, pstst
    ELSE
      EXIT
    END IF
  END DO
END DO

!     FUNCTION VALUES ARE CALCULATED AT AN ADDITIONAL NAP POINTS.

DO i = 1, nap
  i1 = i + 1
  points(i,:) = (g(1,:) + g(i1,:)) * half
ENDDO
!SUBROUTINE function_distribute(npoints,minrank,maxrank,points,vals,nevalps,func)
CALL function_distribute(nap,0,nprocs-1,points(1:nap,:),aval(1:nap),nrounds,functn,myrank)
nmore = nmore + nap
neval = neval + nap
nevalp = nevalp + nrounds

!     THE MATRIX OF ESTIMATED SECOND DERIVATIVES IS CALCULATED AND ITS
!     LOWER TRIANGLE STORED IN BMAT.

a0 = h(1)
! Get points to evaluate at
ii=1
DO i = 1, nap
  i1 = i - 1
  i2 = i + 1
  DO j = 1, i1
    j1 = j + 1
    points(ii,:) = (g(i2,:) + g(j1,:)) * half
    ii=ii+1
  ENDDO
ENDDO
npoints=nap*(nap-1)/2
! Evaluate points
!SUBROUTINE function_distribute(npoints,nprocs,points,vals,nevalps,func)
CALL function_distribute(npoints,0,nprocs-1,points(1:npoints,:),values(1:npoints),nrounds,functn,myrank)
! Enter values into bmat
ii=1
DO i = 1, nap
  i1 = i - 1
  i2 = i + 1
  DO j = 1, i1
    j1 = j + 1
    l = i * (i-1) / 2 + j
    bmat(l) = two * (values(ii) + a0 - aval(i) - aval(j))
    ii=ii+1
  END DO
END DO
! Update counts of function evaluations
nmore = nmore + nap*(nap-1)/2
neval = neval + nap*(nap-1)/2
nevalp = nevalp + nrounds
! ORIGINAL NON-PARALLEL VERSION
!DO i = 1, nap
!  i1 = i - 1
!  i2 = i + 1
!  DO j = 1, i1
!    j1 = j + 1
!    pstst = (g(i2,:) + g(j1,:)) * half
!    CALL functn(pstst,hstst)
!    nmore = nmore + 1
!    neval = neval + 1
!    nevalp = nevalp + 1
!    l = i * (i-1) / 2 + j
!    bmat(l) = two * (hstst + a0 - aval(i) - aval(j))
!  END DO
!END DO

l = 0
DO i = 1, nap
  i1 = i + 1
  l = l + i
  bmat(l) = two * (h(i1) + a0 - two*aval(i))
END DO

!     THE VECTOR OF ESTIMATED FIRST DERIVATIVES IS CALCULATED AND
!     STORED IN AVAL.

DO i = 1, nap
  i1 = i + 1
  aval(i) = two * aval(i) - (h(i1) + 3._dp*a0) * half
END DO

!     THE MATRIX Q OF NELDER & MEAD IS CALCULATED AND STORED IN G.

pmin = g(1,:)
DO i = 1, nap
  i1 = i + 1
  g(i1,:) = g(i1,:) - g(1,:)
END DO

DO i = 1, nap
  i1 = i + 1
  g(i,:) = g(i1,:)
END DO

!     INVERT BMAT

CALL syminv(bmat, nap, bmat, temp, nullty, ifault, rmax)
IF (ifault == 0) THEN
  irank = nap - nullty
ELSE                                 ! BMAT not +ve definite
                                     ! Resume search for the minimum
  IF ((iprint >= 0).AND.(realrank==0)) WRITE (lout,6100)
  ifault = 2
  IF (((maxfn>=0).AND.(neval > maxfn)).OR.((maxfn<0).AND.(nevalp > -1*maxfn))) RETURN
  IF (realrank==0) WRITE (lout,6200)
  step = half * step
  GO TO 20
END IF

!     BMAT*A/2 IS CALCULATED AND STORED IN H.

DO i = 1, nap
  h(i) = zero
  DO j = 1, nap
    IF (j <= i) THEN
      l = i * (i-1) / 2 + j
    ELSE
      l = j * (j-1) / 2 + i
    END IF
    h(i) = h(i) + bmat(l) * aval(j)
  END DO
END DO

!     FIND THE POSITION, PMIN, & VALUE, YMIN, OF THE MINIMUM OF THE
!     QUADRATIC.

ymin = DOT_PRODUCT( h(1:nap), aval(1:nap) )
ymin = a0 - ymin
DO i = 1, nop
  pstst(i) = DOT_PRODUCT( h(1:nap), g(1:nap,i) )
END DO
pmin = pmin - pstst
IF ((iprint >= 0).AND.(realrank==0)) THEN
  WRITE (lout,6300) ymin, pmin
  WRITE (lout,6400)
END IF

!     Q*BMAT*Q'/2 IS CALCULATED & ITS LOWER TRIANGLE STORED IN VC

DO i = 1, nop
  DO j = 1, nap
    h(j) = zero
    DO k = 1, nap
      IF (k <= j) THEN
        l = j * (j-1) / 2 + k
      ELSE
        l = k * (k-1) / 2 + j
      END IF
      h(j) = h(j) + bmat(l) * g(k,i) * half
    END DO
  END DO

  DO j = i, nop
    l = j * (j-1) / 2 + i
    vc(l) = DOT_PRODUCT( h(1:nap), g(1:nap,j) )
  END DO
END DO

!     THE DIAGONAL ELEMENTS OF VC ARE COPIED INTO VAR.

j = 0
DO i = 1, nop
  j = j + i
  var(i) = vc(j)
END DO
IF (iprint < 0) RETURN
IF (realrank==0) THEN
  WRITE (lout,6500) irank
  CALL print_tri_matrix(vc, nop, lout)
  WRITE (lout,6600)
ENDIF
CALL syminv(vc, nap, bmat, temp, nullty, ifault, rmax)

!     BMAT NOW CONTAINS THE INFORMATION MATRIX

IF (realrank==0) THEN
  WRITE (lout,6700)
  CALL print_tri_matrix(bmat, nop, lout)
ENDIF

ii = 0
ij = 0
DO i = 1, nop
  ii = ii + i
  IF (vc(ii) > zero) THEN
    vc(ii) = one / SQRT(vc(ii))
  ELSE
    vc(ii) = zero
  END IF
  jj = 0
  DO j = 1, i - 1
    jj = jj + j
    ij = ij + 1
    vc(ij) = vc(ij) * vc(ii) * vc(jj)
  END DO
  ij = ij + 1
END DO

IF (realrank==0) WRITE (lout,6800)
ii = 0
DO i = 1, nop
  ii = ii + i
  IF (vc(ii) /= zero) vc(ii) = one
END DO
IF (realrank==0) CALL print_tri_matrix(vc, nop, lout)

!     Exit, on successful termination.

IF (realrank==0) WRITE (lout,6900) nmore
RETURN

4900 FORMAT (' Estimating ',i4,' parameters on ',i3,' communicators')
5000 FORMAT (' Progress Report every',i4,' function evaluations'/  &
             ' EVAL.   FUNC.VALUE.          PARAMETER VALUES')
5100 FORMAT (/' ', i4, '  ', g12.5, '  ', 5G11.4, 3(/t22, 5G11.4/)/)
5200 FORMAT (' No. of function evaluations > ',i5)
5201 FORMAT (' No. of function evaluations per processor > ',i5)
5300 FORMAT (' RMS of function values of last simplex =', g14.6)
5400 FORMAT (' Centroid of last simplex =',4(/' ', 6G13.5/))
5500 FORMAT (' Function value at centroid =', g14.6)
5600 FORMAT (/' EVIDENCE OF CONVERGENCE'/)
5700 FORMAT (/' Minimum found after',i5,' function evaluations'/)
5701 FORMAT (/' Minimum found after',i5,' function evaluations per processor'/)
5800 FORMAT (' Minimum at',4(/' ', 6G13.6/))
5900 FORMAT (' Function value at minimum =', g14.6)
6000 FORMAT (/' Fitting quadratic surface about supposed minimum'/)
6100 FORMAT (/' MATRIX OF ESTIMATED SECOND DERIVATIVES NOT +VE DEFN.'/ &
             ' MINIMUM PROBABLY NOT FOUND'/)
6200 FORMAT (/t11, 'Search restarting'/)
6300 FORMAT (' Minimum of quadratic surface =',g14.6,' at',4(/' ', 6G13.5/))
6400 FORMAT (' IF THIS DIFFERS BY MUCH FROM THE MINIMUM ESTIMATED ',      &
             'FROM THE MINIMIZATION,'/' THE MINIMUM MAY BE FALSE &/OR THE ',  &
             'INFORMATION MATRIX MAY BE INACCURATE')
6500 FORMAT (' Rank of information matrix =',i3/ &
             ' Inverse of information matrix:-')
6600 FORMAT (/' If the function minimized was -LOG(LIKELIHOOD),'/         &
             ' this is the covariance matrix of the parameters.'/         &
             ' If the function was a sum of squares of residuals,'/       &
             ' this matrix must be multiplied by twice the estimated ',   &
             'residual variance'/' to obtain the covariance matrix.'/)
6700 FORMAT (' INFORMATION MATRIX:-')
6800 FORMAT (/' CORRELATION MATRIX:-'/)
6900 FORMAT (/' A further',i4,' function evaluations have been used'/)

call MPI_Finalize(mpierr)

END SUBROUTINE pminim




SUBROUTINE syminv(a, n, c, w, nullty, ifault, rmax)

!     ALGORITHM AS7, APPLIED STATISTICS, VOL.17, 1968.

!     ARGUMENTS:-
!     A()    = INPUT, THE SYMMETRIC MATRIX TO BE INVERTED, STORED IN
!                LOWER TRIANGULAR FORM
!     N      = INPUT, ORDER OF THE MATRIX
!     C()    = OUTPUT, THE INVERSE OF A (A GENERALIZED INVERSE IF C IS
!                SINGULAR), ALSO STORED IN LOWER TRIANGULAR FORM.
!                C AND A MAY OCCUPY THE SAME LOCATIONS.
!     W()    = WORKSPACE, DIMENSION AT LEAST N.
!     NULLTY = OUTPUT, THE RANK DEFICIENCY OF A.
!     IFAULT = OUTPUT, ERROR INDICATOR
!                 = 1 IF N < 1
!                 = 2 IF A IS NOT +VE SEMI-DEFINITE
!                 = 0 OTHERWISE
!     RMAX   = OUTPUT, APPROXIMATE BOUND ON THE ACCURACY OF THE DIAGONAL
!                ELEMENTS OF C.  E.G. IF RMAX = 1.E-04 THEN THE DIAGONAL
!                ELEMENTS OF C WILL BE ACCURATE TO ABOUT 4 DEC. DIGITS.

!     LATEST REVISION - 1 April 1985

!***************************************************************************

REAL(8), INTENT(IN OUT) :: a(:), c(:), w(:)
INTEGER, INTENT(IN)       :: n
INTEGER, INTENT(OUT)      :: nullty, ifault
REAL(8), INTENT(OUT)    :: rmax

REAL(8), PARAMETER :: zero = 0._dp, one = 1._dp
INTEGER              :: i, icol, irow, j, jcol, k, l, mdiag, ndiag, nn, nrow
REAL(8)            :: x

nrow = n
ifault = 1
IF (nrow > 0) THEN
  ifault = 0

!     CHOLESKY FACTORIZATION OF A, RESULT IN C

  CALL chola(a, nrow, c, nullty, ifault, rmax, w)
  IF (ifault == 0) THEN

!     INVERT C & FORM THE PRODUCT (CINV)'*CINV, WHERE CINV IS THE INVERSE
!     OF C, ROW BY ROW STARTING WITH THE LAST ROW.
!     IROW = THE ROW NUMBER, NDIAG = LOCATION OF LAST ELEMENT IN THE ROW.

    nn = nrow * (nrow+1) / 2
    irow = nrow
    ndiag = nn
    10 IF (c(ndiag) /= zero) THEN
      l = ndiag
      DO i = irow, nrow
        w(i) = c(l)
        l = l + i
      END DO
      icol = nrow
      jcol = nn
      mdiag = nn

      30 l = jcol
      x = zero
      IF (icol == irow) x = one / w(irow)
      k = nrow
      40 IF (k /= irow) THEN
        x = x - w(k) * c(l)
        k = k - 1
        l = l - 1
        IF (l > mdiag) l = l - k + 1
        GO TO 40
      END IF

      c(l) = x / w(irow)
      IF (icol == irow) GO TO 60
      mdiag = mdiag - icol
      icol = icol - 1
      jcol = jcol - 1
      GO TO 30
    END IF ! (c(ndiag) /= zero)

    l = ndiag
    DO j = irow, nrow
      c(l) = zero
      l = l + j
    END DO

    60 ndiag = ndiag - irow
    irow = irow - 1
    IF (irow /= 0) GO TO 10
  END IF
END IF
RETURN

END SUBROUTINE syminv



SUBROUTINE chola(a, n, u, nullty, ifault, rmax, r)

!     ALGORITHM AS6, APPLIED STATISTICS, VOL.17, 1968, WITH
!     MODIFICATIONS BY A.J.MILLER

!     ARGUMENTS:-
!     A()    = INPUT, A +VE DEFINITE MATRIX STORED IN LOWER-TRIANGULAR
!                FORM.
!     N      = INPUT, THE ORDER OF A
!     U()    = OUTPUT, A LOWER TRIANGULAR MATRIX SUCH THAT U*U' = A.
!                A & U MAY OCCUPY THE SAME LOCATIONS.
!     NULLTY = OUTPUT, THE RANK DEFICIENCY OF A.
!     IFAULT = OUTPUT, ERROR INDICATOR
!                 = 1 IF N < 1
!                 = 2 IF A IS NOT +VE SEMI-DEFINITE
!                 = 0 OTHERWISE
!     RMAX   = OUTPUT, AN ESTIMATE OF THE RELATIVE ACCURACY OF THE
!                DIAGONAL ELEMENTS OF U.
!     R()    = OUTPUT, ARRAY CONTAINING BOUNDS ON THE RELATIVE ACCURACY
!                OF EACH DIAGONAL ELEMENT OF U.

!     LATEST REVISION - 1 April 1985

!***************************************************************************

REAL(8), INTENT(IN)   :: a(:)
INTEGER, INTENT(IN)     :: n
REAL(8), INTENT(OUT)  :: u(:)
INTEGER, INTENT(OUT)    :: nullty, ifault
REAL(8), INTENT(OUT)  :: rmax, r(:)

!     ETA SHOULD BE SET EQUAL TO THE SMALLEST +VE VALUE SUCH THAT
!     1._dp + ETA IS CALCULATED AS BEING GREATER THAN 1._dp IN THE ACCURACY
!     BEING USED.

REAL(8), PARAMETER :: eta = EPSILON(1.0_dp), zero = 0._dp
INTEGER              :: i, icol, irow, j, k, l, m
REAL(8)            :: rsq, w

ifault = 1
IF (n > 0) THEN
  ifault = 2
  nullty = 0
  rmax = eta
  r(1) = eta
  j = 1
  k = 0

!     FACTORIZE COLUMN BY COLUMN, ICOL = COLUMN NO.

  DO  icol = 1, n
    l = 0

!     IROW = ROW NUMBER WITHIN COLUMN ICOL

    DO  irow = 1, icol
      k = k + 1
      w = a(k)
      IF (irow == icol) rsq = (w*eta) ** 2
      m = j
      DO  i = 1, irow
        l = l + 1
        IF (i == irow) EXIT
        w = w - u(l) * u(m)
        IF (irow == icol) rsq = rsq + (u(l)**2*r(i)) ** 2
        m = m + 1
      END DO

      IF (irow == icol) EXIT
      IF (u(l) /= zero) THEN
        u(k) = w / u(l)
      ELSE
        u(k) = zero
        IF (ABS(w) > ABS(rmax*a(k))) GO TO 60
      END IF
    END DO

!     END OF ROW, ESTIMATE RELATIVE ACCURACY OF DIAGONAL ELEMENT.

    rsq = SQRT(rsq)
    IF (ABS(w) > 5.*rsq) THEN
      IF (w < zero) RETURN
      u(k) = SQRT(w)
      r(i) = rsq / w
      IF (r(i) > rmax) rmax = r(i)
    ELSE
      u(k) = zero
      nullty = nullty + 1
    END IF

    j = j + icol
  END DO
  ifault = zero
END IF
60 RETURN

END SUBROUTINE chola



SUBROUTINE print_tri_matrix(a, n, lout)

INTEGER, INTENT(IN)    :: n, lout
REAL(8), INTENT(IN)  :: a(:)

!     Local variables
INTEGER  :: i, ii, i1, i2, l

l = 1
DO l = 1, n, 6
  ii = l * (l-1) / 2
  DO i = l, n
    i1 = ii + l
    ii = ii + i
    i2 = MIN(ii,i1+5)
    WRITE (lout,'(1X,6G13.5)') a(i1:i2)
  END DO
END DO
RETURN

END SUBROUTINE print_tri_matrix

! QsortC2
!  IF have and arrays A and B
! with each element of B corresponding to an entry of A,
! sort both A and the B by the entries of A.
! Like row sort.

recursive subroutine QsortC2(A,B)
  real(8), intent(in out), dimension(:) :: A
  real(8), intent(in out), dimension(:) :: B
  integer :: iq

  if(size(A) > 1) then
     call Partition2(A,B,iq)
     call QsortC2(A(:iq-1),B(:iq-1))
     call QsortC2(A(iq:),B(iq:))
  endif
end subroutine QsortC2

subroutine Partition2(A,B,marker)
  real(8), intent(in out), dimension(:) :: A
  real(8), intent(in out), dimension(:) :: B
  integer, intent(out) :: marker
  integer :: i, j
  real(8) :: temp
  real(8) :: temp1
  real(8) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j), B(i) and B(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        temp1 = B(i)
        B(i) = B(j)
        B(j) = temp1
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition2

! QsortC3
!  IF have and array A and 2D array B
! with each row of B corresponding to an entry of A,
! sort both A and the rows of B by the entries of A.
! Like row sort.

recursive subroutine QsortC3(A,B)
  real(8), intent(in out), dimension(:) :: A
  real(8), intent(in out), dimension(:,:) :: B
  integer :: iq

  if(size(A) > 1) then
     call Partition3(A,B,iq)
     call QsortC3(A(:iq-1),B(:iq-1,:))
     call QsortC3(A(iq:),B(iq:,:))
  endif
end subroutine QsortC3

subroutine Partition3(A,B,marker)
  real(8), intent(in out), dimension(:) :: A
  real(8), intent(in out), dimension(:,:) :: B
  integer, intent(out) :: marker
  integer :: i, j
  real(8) :: temp
  real(8), dimension(size(B,2)) :: temp1
  real(8) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j), B(i,:) and B(j,:)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        temp1 = B(i,:)
        B(i,:) = B(j,:)
        B(j,:) = temp1
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition3

! function_distribute: Evaluate a function at an array of points, dividing
!     the function calls among the available processors
SUBROUTINE function_distribute(npoints,minrank,maxrank,points,vals,nevalsp,functn,myrank)
    INTEGER, INTENT(IN) :: npoints ! number of points to be evaluated
    INTEGER, INTENT(IN) :: minrank,maxrank ! range of processes/groups to distribute over
    REAL(8), DIMENSION(:,:), INTENT(IN) ::  points ! npoints x nop points to evaluate at
    REAL(8), DIMENSION(:), INTENT(OUT) :: vals ! function values at each point
    INTEGER, INTENT(OUT) :: nevalsp ! number of evaluations per processor
    EXTERNAL functn ! function to be evaluated
    INTEGER, INTENT(IN) :: myrank ! group/communicator of current processor

    INTEGER :: nrounds,nround,mpierr,i,nprocs,IERR

    real(dp) :: timing(2)
    
    !CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank,mpierr)
	!CALL MPI_Comm_rank(comm,i,mpierr)
	!print*, "Here is ", myrank,i !,myrank_world
    !CALL MPI_Comm_size(MPI_COMM_WORLD,nprocs,mpierr)
    !print*, "myrank,myrank_world ",myrank,myrank_world
    
    nprocs=maxrank-minrank+1
    nrounds=npoints/nprocs ! number of full rounds
    DO nround=1,nrounds
      DO i=1,nprocs
        IF (i==myrank-minrank+1) THEN
            timing(1)=MPI_WTIME()
		CALL functn(points(i+nprocs*(nround-1),:),vals(i+nprocs*(nround-1)))	
            timing(2)=MPI_WTIME()		
            !write(*,'("After functn ",2I4,F14.4,F10.2)') myrank,i+nprocs*(nround-1),vals(i+nprocs*(nround-1)),timing(2)-timing(1)
        ENDIF
      ENDDO
      DO i=1,nprocs
		!PRINT*, "Bfre bcast ",MYRANK,i+nprocs*(nround-1),vals(i+nprocs*(nround-1))
        CALL MPI_BCAST(vals(i+nprocs*(nround-1)),1,MPI_REAL8,i-1,MPI_COMM_WORLD,mpierr) 		
		!IF (MYRANK_WORLD==I) PRINT*, MYRANK,i+nprocs*(nround-1)
		!PRINT*, "After bcast ",MYRANK,i+nprocs*(nround-1),vals(i+nprocs*(nround-1))
	  ENDDO
    ENDDO
	  		!print*, myrank,myrank_world,vals(1:15)

			
			
    DO i=1,npoints-nprocs*nrounds
      IF (i==myrank-minrank+1) THEN
	  !print*, "second if ",myrank,myrank_world,  i+nprocs*(nround-1)
        CALL functn(points(i+nprocs*nrounds,:),vals(i+nprocs*nrounds))
      ENDIF
    ENDDO
    DO i=1,npoints-nrounds*nprocs
	CALL MPI_BCAST(vals(i+nrounds*nprocs),1,MPI_REAL8,i-1,MPI_COMM_WORLD,mpierr)
    ENDDO
    
	!if (myrank==1) then 
	!	print*, "myrank and vals ",myrank,
	
	!	print*, vals
	!end if 
    IF (npoints==nprocs+nrounds) THEN
      nevalsp=nrounds
    ELSE
      nevalsp=nrounds+1
    ENDIF

END SUBROUTINE function_distribute

END MODULE pNelder_Mead

