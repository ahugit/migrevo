!what are all these compile options?
!#mpif90 -O3 -fPIC -unroll -ip -axavx -xsse4.2 -openmp -vec-report -par-report -openmp-report -o train.x program.f90
!mpif90 -O3 -fPIC -unroll -ip -axavx -xsse4.2 -openmp -vec -par -openmp -o train.x program.f90
!#mpif90 train.x program.f90
	
program main 
	use params !, only: myid,iter,parname,npars,stepmin
	use objf
	use pnelder_mead  
    use alib, only: FINDinv !for FINDinv

	!use alib, only: logitinv
	implicit none	
	!include 'mpif.h'
	real(8), dimension(npars) :: pars1,pars2,pars3,parsfuck,pars0
	integer, parameter :: nj=10,nj2=2*nj
	integer :: i,j,k,pp,ierror,k1,k2
	! declarations for simulated annealing: at temperature t, a new point that is worse by an amount will be accepted with probability exp(-q/t)
	real(8) :: pars(npars),incr,realpars(npars),realpars1(npars),val
	real(sp) :: q1val(nj2),step(nj2),lb,ub
    integer :: maxfn,iprint,nloop,iquad,ifault
    real(8) :: qval,stopcr,simp,var(npars)
	!ahu 121818: the below are no longer declared as parameters. their value is set below in the program. 
    real(8) :: tstart !=15. 			! starting temperature (set to zero to turn off simulated annealing)
	real(8) :: tstep !=0.8			! fraction temperature reduced at each step
	integer :: tfreq !=30			! number of function calls between temperature reductions
	integer :: saseed !=1			! seed for random numbers
	!ahu 121818: the above are no longer declared as parameters. their value is set below in the program. 
    
    ! for the new sa routine that does grouping
	integer, dimension(0:ninp-1) :: myhomes != (/ 1,2,3,4,5,6,7,8,9 /) !,1,2,3,4,5,6,7,8,9 /)              ! types !mytyps = (/1,1,2,2/)
	integer, dimension(0:ninp-1) :: mytyps  != 1 !(/ 1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2 /)              ! types 
	integer, dimension(0:ninp-1) :: mycos != 1
	!integer, dimension(0:ninp-1) :: myindexs != (/ 1,2,3,4,5,6,7,8,9 /)  !,10,11,12,13,14,15,16,17,18 /)    ! types    !mytyps = (/1,1,2,2/)
	! keep track of time 
	real(4):: begintime,runtime
	real(8) :: temp
	integer :: numgroup,numworld,iam,mygrank0 !mygrank,ngroups,mygroup,nprocs_world,iam
	integer :: mpi_group_world,mygsize
	integer :: mpierr,mpistat(mpi_status_size)
    real(4), allocatable, dimension(:,:) :: mytime 
    !for standard errors:
    real(8) :: D(nmom,npars),WD(nmom,npars),QWD(nmom,npars)
    real(8), dimension(npars,npars) :: DpWDinv,SEmat,parvo1,realparvo1
    real(8) :: dtheta(npars),worksp(npars)
    integer :: activepari(npars)
    real(dp) :: QQ(nmom),stderrs(npars),stepstderr(npars),val1
    integer(i4b) :: kactive,eflag,im,kk,ka
    integer :: nactive,itermin1 !to take into account the fact that iter is augmented by 1 at the end of objf so need to undo that
    
    !integer, parameter :: lenrunid=5
    !character(LEN=lenrunid) :: runid='test' ! SEY OM OPTIONS.TXT

    !OPEN(5,FILE='options.txt')
    !READ(5,*) runid
    !READ(5,*) estimate
    !READ(5,*) stepsize0
    !READ(5,*) fixtaua
    !READ(5,*) fixtauth
    !READ(5,*) min_wage
    !READ(5,*) bootrun
    !CLOSE(5)
    

	call mpi_init(mpierr)
	call mpi_comm_rank(mpi_comm_world,iam,mpierr)
	call mpi_comm_size(mpi_comm_world,numworld,mpierr)
	print*, "numworld,Iam ", numworld,iam
    mysay=iam
	conditional_moments=.true.	;     terminalval=.TRUE.	
    if (iam==0) print*, "Here is numworld", numworld

    allocate(mytime(numworld,2))
   
	if (groups) then 	
		if (mod(numworld,ninp)>0) then ; print*, "numworld needs to be a multiple of ninp " ; end if 

		numgroup=int(numworld/ninp)		!number of groups  36/18=2 
		mygroup=mod(iam,numgroup) ! which group each process belongs to: myid mygroup: 0 0 / 1 1 / 2 0 / 3 1 / 4 0 / 5 1 / ..... / 34 0 / 35 1 
		mygrank=int(iam/numgroup) !id within group: myid mygroup: 0 0 / 1 0 / 2 1 / 3 1 / 4 2 / 5 2 / 6 3 / 7 3 / 8 4 / 9 4 / ..../ 34 17 / 35 17 
        !ahu 0317: BEWARE OF DOING THINGS LIKE BELOW. IF OYU SET THE MYGROUP THE WAY IT IS DONE BELOW THEN THE FUNCTION_DISTRIBUTE 
        !WILL NOT WORK PROPERLY. WHY? BEACAUSE MPI_BCAST DOES NOT DO WHAT IT SHOULD DO IN THIS CASE. 
        !BECAUSE THEN FOR I=1,NPROCS WHERE NPROCS IS THE NUMBER OF GROUPS, FOR I=1, THE ROOT PROCESS IS 0 AND FOR I=2, THE ROOT PROCESS IN BCAST IS 1. 
        !IF YOU DO BELLOW, WHEN BCAST SAYS SEND VALUE FROM ROOT PROCESS 1 TO EVERYONE, IT DOES NOT TAKE ROOT PROCESS 1 TO BE THE PROCEESSSES IN GROUP 2 
        !BUT IT JUST TAKES IT TO BE THE PROCESSOR WITH IAM=1. SO THEN WHEN YOU'RE GOING THROUGH I=1,NPROCS IN THE BCAST PART OF FUNCTION_DISTRIBUTE
        !YOU ARE NOT ACTUALLY GOING THROUGH GROUPS. 
        !BUT IF YOU DO THE ABOVE SETTING OF MYGROUP, THEN THE TWO ARE EQUIVALENT SINCE THE IAM'S ALSO CORRESPOND TO GROUP NUMBERS. (JUST DUE TO THE WAY MYGROUP IS SET)
        !numgroup=int(numworld/ninp)		!number of groups  36/18=2
        !mygroup=int(iam/ninp) ! which group each process belongs to: myid mygroup: 0 0 / 1 0 / 2 0 / 3 0 / 4 0 /.../17 0 / 18 1 /...../35 1 
		!mygrank=mod(iam,ninp) !id within group: myid mygroup: 0 0 / 1 1 / 2 2 / 3 3 / 4 4 / 5 5 / 6 6 / 7 7 / 8 8 ..../ 17 17 / 18 0 / 19 1 / ..../35 17 
       
        call mpi_comm_split(mpi_comm_world,mygroup,mygrank,comm,mpierr)
        !ahu 0317 call mpi_comm_split(mpi_comm_world,mygroup,iam,comm,mpierr) !ahu 0317
        !print*, "after split ", iam
        !call MPI_Comm_rank(comm, mygrank,mpierr) !ahu 0317
        !print*, "after rank ", iam
        !call MPI_Comm_size(comm, mygsize,mpierr) !ahu 0317
        !print*, "after size ", iam
        !if (mygrank0.ne.mygrank) then 
        !    print*, 'The two ways of calculating mygrank are not equal'
        !    stop
        !end if 
       ! write(*,'("iam,numgroup,mygroup,mygrank,mygsize:")') !ahu 0317
       ! write(*,'(5i4)') iam,numgroup,mygroup,mygrank,mygsize !ahu 0317
        
        pp=0
        do i=1,ncop
            do j=1,ntypp
                do k=1,nhomep
                    mycos(pp)=i
                    mytyps(pp)=j
                    myhomes(pp)=k
                    !myindexs(pp)=pp
                    pp=pp+1
                end do 
            end do 
        end do 
        if ( (pp-1).ne.(ncop*ntypp*nhomep-1)) then
            print*, 'something wrong with pp',pp,ncop*ntypp*nhomep
            stop
        end if
        if (iam==0) then 
            write(*,'("mycos,mytyps,myhomes,myindexs")') 
            if (ntypp==1) then 
                write(*,'(9I4)') mycos
                write(*,'(9I4)') mytyps
                write(*,'(9I4)') myhomes
                !write(*,'(9I4)') myindexs
            else if (ntypp==2) then
                write(*,'(18I4)') mycos
                write(*,'(18I4)') mytyps
                write(*,'(18I4)') myhomes
                !write(*,'(18I4)') myindexs
            end if
        end if
		
		!myindex = myindexs(mygrank)
		myco   = mycos(mygrank)
		mytyp  = mytyps(mygrank)
		myhome = myhomes(mygrank)
	else 
		numgroup=0 ; mygroup=0 ; mygrank=0 ; myco=0 ; mytyp=0 ; myhome=0
    end if
    print*, iam,mycos(iam),mytyps(iam),myhomes(iam)
	begintime=secnds(0.0)
    iter=1
	iwritegen=0
	if (iam==0) iwritegen=1
	if (iwritegen==1.and.chkobj) then						! ahu april13: to check how much objval changes with each stepsize and make sure that the objval changes by similar amounts when changing each parameter. otherwise the simplex does weird things. and also to check how much the pminim routine changes objval with each iteration during estimation but the latter is not as important. 
		open(unit=63, file='chkobj.txt',status='replace')	
	end if 
	if (iam==0.and.chkstep) then 
		open(unit=64, file='chkstep.txt',status='replace')
		open(unit=65, file='step.txt',status='replace')
	end if 
	if (iam==0.and.comparepars) then 
		open(unit=64, file='chkstep.txt',status='replace')
	end if     
	!ahu030622 skriv=.false. ; if (iam==24.and.iter==1) skriv=.false.  !ahu 0327 iam=9
	skriv=.false. ; if (iam==5) skriv=.false.  !and iter=1 !ahu030622
    if (.not.groups) skriv=.false. !ahu 0327
    iwritemom=0 ; if (iam==0) iwritemom=1 
	if (skriv) then 	
		open(unit=12, file='chkdat.txt',status='replace') !written in call getmom/yaz_getmom	
		open(unit=40,file='chkq.txt',status='replace') 
		open(unit=50,file='chkpars.txt',status='replace')    !written in objf/yaz0
		open(unit=100,file='chksolgetdec_s.txt',status='replace')
        !written in sol/couples/getdec/yaz_getdec (callfrom==40), sol/marmkt/getdec/yaz_getdec (callfrom=50)
		open(unit=200,file='chksolgetdec_c.txt',status='replace') 
		open(unit=201,file='chk2b.txt',status='replace')		
        !400 is written in simulate, simulate/getdec/yaz_getdec (callfrom=80), yaz_decs,yaz_sim,yaz_simmatch,yaz_simdecmar 
        open(unit=400, file='chksim.txt',status='replace')   
		open(unit=500, file='chksimpath.txt',status='replace')	 !written by simulate/yaz_simpath
	end if 
    writestderr=.FALSE.
    if (iam==0.and.getstderr) writestderr=.TRUE.    
    if (writestderr) OPEN(13,FILE='stderrors.txt',STATUS='REPLACE')  
    policytax=0  
	!call cohab_stderr(parvector,stepsize,stderrs)
    nonneg=.TRUE.
    !onthejobsearch=.TRUE.
    nonlabinc=0.0_dp
    numperdat=6812 !6921  !6825 !7067 !7106 !7209 !7538
    allocate(init(numperdat)) !dat is deallocated at the end of the numit if statement but init is not, because I need it for the simulations in all calls! 
    numperobsdat=105599 !105978  !102138 !108022 !108252 !114537
    numpersim=numperdat*nsimeach
    datafilename='familymig2022-cendiv.txt'
    !taxset=1
    allocate(nummove_save(numpersim),nummovemar_save(numpersim),nummovesin_save(numpersim) )      
    nummove_save=-99 ; nummovemar_save=-99 ; nummovesin_save=-99
    policytax=0 ; policyubenefit=0
    !terminalval=.TRUE.    
    !taxset=0
    !numperdat=5233
    !numperobsdat=84507
    !numpersim=numperdat*nsimeach
    !datafilename='familymigpsid.txt'
    METHODWAGEDISCRETIZE=0
    TAXRATES=1
    JOINTDIST=.TRUE.
    POLICYTAX=0
    JOINTMARSHOCK=.FALSE.
    UPFRONTTRANSFERS=.FALSE.
    nonlabinc=0.0_dp
    JOINTMARSHOCK=.TRUE.
    TAXRATES=1
    JOINTDIST=.TRUE.
    POLICYTAX=6         !MEANS NO TAXES
    JOINTMARSHOCK=.FALSE.
    UPFRONTTRANSFERS=.FALSE.
    nonlabinc=0.0_dp
    JOINTMARSHOCK=.TRUE.
    POLICYUBENEFIT=0
    NOMISSINGOBS=.FALSE.
    TERMINALVAL=.FALSE. !SET NTERMVAL IF TERMINALVAL=TRUE


    if (icheck_eqvmvf) then
		pars(3)=pars(2)
		pars(8)=pars(7)
		pars(10)=pars(9)
		pars(34:45)=pars(22:33)
		pars(47)=pars(46)
		pars(59:60)=pars(57:58)
		pars(53:56)=pars(49:52)
	end if 	
	if (icheck_eqvcvs) then 
		pars(1)=-huge(0.0)       ! sig_mar
		pars(6)=-huge(0.0)       ! divpenalty
		pars(npars)=0.0_dp          ! mu_mar
	end if 
    

    !open(unit=2,file='o062323_2bpobj.txt',status='old',action='read') ; read(2,*) pars2	; close(2)
    pars=pars2
    pars(75)=pars0(75)
    pars(81)=pars0(81)
    pars(87)=pars0(87)
    pars(93)=pars0(93)
    pars(22:23)=pars0(22:23)
    pars(24)=pars0(24)
    pars(70)=0.01507_dp  !     0.01467       0.01554      0.0005
    pars(71)=-0.00018_dp 
    pars(15)=pars(15)+1000.0_dp 
    pars(11)=pars(11)-0.6_dp 
    pars(12)=pars(12)-0.5_dp                                                     
    pars(27)=pars(27)-1200.0_dp 
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 

    pars(28)=-22045.0_dp 
    pars(74)=-14045.0_dp 
    pars(7)=pars(7)-0.5_dp
    pars(10)=pars(10)+0.5_dp 
    pars(90)=8.75_dp
    pars(91)=8.59_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 
    pars0=pars

    !open(unit=2,file='o062623_bpobj.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(70:71)=pars0(70:71)
    pars(70)= 0.01507    
    pars(71)=-0.00018_dp
    pars(80)=-27325.0_dp    !cstedf
    pars(3)=pars0(3)
    pars(11)=pars(11)-0.25_dp
    pars(22:23)=pars(22:23)-4000.0_dp
    pars(74)=pars(74)-3000.0_dp
    UPFRONTTRANSFERS=.FALSE.
    policyubenefit=0
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 
    !pars1=pars

    open(unit=2,file='paramsBEST070323.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    UPFRONTTRANSFERS=.FALSE.
    policyubenefit=0
    PARS(25)=0.0_DP !CSTU(MLES)
    PARS(65)=0.0_DP !CSTU(FEMALES)
    pars(93)=pars(93)-4500.0_dp 
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 

    open(unit=2,file='paramsBEST070723.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    UPFRONTTRANSFERS=.FALSE.
    policyubenefit=0
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 
    pars1=pars

    pars(29)=-16000.0_dp  !pars(29)-4000.0_dp
    pars(80)=-16000.0_dp  !pars(80)-4000.0_dp
    pars(37)=3658.0_dp
    pars(15)=pars(15)+2000.0_dp 
    pars(80)=pars(80)-2000.0_dp 
    pars(12)=pars(12)+0.3_dp 
    pars(15)=pars(15)+1000.0_dp 
    pars(80)=pars(80)-2000.0_dp
    pars(80)=pars(80)-2000.0_dp
    pars(15)=pars(15)-1000.0_dp 
    pars(29)=pars(29)-4000.0_dp
    pars(15)=pars(15)-1000.0_dp 
    pars(15)=pars(15)-1000.0_dp 
    pars(29)=pars(29)-2000.0_dp
    pars(14)=pars(14)-1000.0_dp 
    pars(12)=pars(12)+0.05_dp
    pars(29)=-22472.0_dp
    pars(80)=-22153.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 


    open(unit=2,file='paramsBEST071023.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    UPFRONTTRANSFERS=.FALSE.
    policyubenefit=0
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 
    pars1=pars

    pars(14)=6085.0_dp
    pars(28)=-22372.0_dp 
    pars(74)=-15845.0_dp 
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 

    pars(37)=4658.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 

    pars(37)=5658.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 

    !72 73  78 79  84 85   90 91 

    !*************** 071213: GOING BACK TO 071013 MOM1 WHIHCH WAS SOME KIND OF GOING BACK IN AND OF ITSELF. SEE COMMITS. 
    open(unit=2,file='paramsBEST070323.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    UPFRONTTRANSFERS=.FALSE.
    policyubenefit=0
    PARS(25)=0.0_DP !CSTU(MLES)
    PARS(65)=0.0_DP !CSTU(FEMALES)
    pars(93)=pars(93)-4500.0_dp 
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 


    open(unit=2,file='paramsBEST070723.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    UPFRONTTRANSFERS=.FALSE.
    policyubenefit=0
    call getpars(pars,realpars)
    call objfunc(pars,qval) ; realpars=realpartemp 



    pars1=pars
    if(policyruns) then
       !UPFRONTTRANSFERS=.FALSE.
       !policyubenefit=0
       !pars(28)=-999999.0_dp
       !pars(29)=-999999.0_dp
       !pars(74)=-999999.0_dp
       !pars(80)=-999999.0_dp
       !call getpars(pars,realpars)
       !call objfunc(pars,qval) ; realpars=realpartemp     
       !pars=pars1



        pars=pars1
        UPFRONTTRANSFERS=.FALSE.
        policyubenefit=2
        call getpars(pars,realpars)
        call objfunc(pars,qval) ; realpars=realpartemp 

        pars=pars1
        UPFRONTTRANSFERS=.TRUE.
        policyubenefit=0
        call getpars(pars,realpars)
        call objfunc(pars,qval) ; realpars=realpartemp 

        !pars=pars1
        !UPFRONTTRANSFERS=.FALSE.
        !nonlabinc=3000.0_dp
        !policyubenefit=0
        !call getpars(pars,realpars)
        !call objfunc(pars,qval) ; realpars=realpartemp   
!
        !pars=pars1
        !UPFRONTTRANSFERS=.FALSE.
        !nonlabinc=6000.0_dp
        !policyubenefit=0
        !call getpars(pars,realpars)
        !call objfunc(pars,qval) ; realpars=realpartemp   
        !
!
        !pars=pars1
        !UPFRONTTRANSFERS=.FALSE.
        !nonlabinc=10000.0_dp
        !policyubenefit=0
        !call getpars(pars,realpars)
        !call objfunc(pars,qval) ; realpars=realpartemp   
!
        !pars=pars1
        !UPFRONTTRANSFERS=.FALSE.
        !nonlabinc=35000.0_dp
        !policyubenefit=0
        !call getpars(pars,realpars)
        !call objfunc(pars,qval) ; realpars=realpartemp   

     
    end if 
    !extra policies that are not being used anymore: 
    !pars=pars1
    !UPFRONTTRANSFERS=.FALSE.
    !nonlabinc=50000.0_dp
    !policyubenefit=0
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 


    !pars=pars1
    !UPFRONTTRANSFERS=.TRUE.
    !policyubenefit=0
    !pars(28)=-999999.0_dp
    !pars(29)=-999999.0_dp
    !pars(74)=-999999.0_dp
    !pars(80)=-999999.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 


!
!
    !pars=pars1
    !UPFRONTTRANSFERS=.FALSE.
    !policyubenefit=3
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 
!
!
    !pars=pars1
    !UPFRONTTRANSFERS=.FALSE.
    !policyubenefit=4
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp 


    if(policyruns_mumar) then
        pars=pars1
        UPFRONTTRANSFERS=.FALSE.
        policyubenefit=0
        PARS(75)=pars(87) !2500.0_dp 
        PARS(81)=PARS(75)
        PARS(87)=PARS(75)
        PARS(93)=PARS(75)
        call getpars(pars,realpars)
        call objfunc(pars,qval) ; realpars=realpartemp 

        pars=pars1
        PARS(75)=pars(81) !4500.0_dp 
        PARS(81)=PARS(75)
        PARS(87)=PARS(75)
        PARS(93)=PARS(75)
        call getpars(pars,realpars)
        call objfunc(pars,qval) ; realpars=realpartemp 

        pars=pars1
        PARS(75)=4500.0_dp 
        PARS(81)=PARS(75)
        PARS(87)=PARS(75)
        PARS(93)=PARS(75)
        call getpars(pars,realpars)
        call objfunc(pars,qval) ; realpars=realpartemp 

        pars=pars1
        !PARS(75)=8500.0_dp 
        PARS(81)=PARS(75)
        PARS(87)=PARS(75)
        PARS(93)=PARS(75)
        call getpars(pars,realpars)
        call objfunc(pars,qval) ; realpars=realpartemp 


        pars=pars1
        PARS(75)=19500.0_dp 
        PARS(81)=PARS(75)
        PARS(87)=PARS(75)
        PARS(93)=PARS(75)
        call getpars(pars,realpars)
        call objfunc(pars,qval) ; realpars=realpartemp 


        PARS(75)=32500.0_dp 
        PARS(81)=PARS(75)
        PARS(87)=PARS(75)
        PARS(93)=PARS(75)
        call getpars(pars,realpars)
        call objfunc(pars,qval) ; realpars=realpartemp 
    end if


    !pars(76:77)=0.0_dp 
    !pars(82:83)=0.0_dp
    !pars(88:89)=0.0_dp
    !pars(96:101)=0.0_dp

    stepmin=0.0_dp
    stepmin=stepos

    !stepmin=0.0_dp
    !stepmin(28:29)=2500.0_dp 
    !stepmin(74)=2500.0_dp 
    !stepmin(80)=2500.0_dp 
    !stepmin(25)=2000.0_dp 
    !stepmin(65)=2000.0_dp 
    !stepmin(72:73)=0.12_dp
    !stepmin(78:79)=0.12_dp
    !stepmin(84:85)=0.12_dp
    !stepmin(90:91)=0.12_dp
    !stepmin(19)=0.9_dp
    !stepmin(20)=1000.0_dp
    !stepmin(21)=-1.5_dp
    !stepmin(14:15)=-750.0_dp
    !stepmin(9:12)=0.6_dp

    !stepmin(17:18)=0.3_dp
    !stepmin(51:52)=0.3_dp
    !stepmin(63:64)=0.3_dp
    !!stepmin(66:67)=0.3_dp
    !stepmin(70)=0.0004
    !stepmin(71)=-0.00003
    !stepmin(94:95)=stepmin(70:71)
    !stepmin(72:73)=0.1_dp
    !stepmin(78:79)=0.1_dp
    !stepmin(84:85)=0.1_dp
    !stepmin(90:91)=0.1_dp
    !stepmin(16)=stepos(16)

    if (optimize) then 
		! call simplex: the call to minim (the simplex routine) is set up to work with starting vector parvector
		!if (iam==0) then ; iprint=1 ; else ; iprint=-1 ; end if 
		!open(unit=64,file='step11_10.txt',status='old',action='read') ; read(64,*) stepmin	; close(64) 			
		!open(unit=64,file='step111416.txt',status='old',action='read') ; read(64,*) stepmin	; close(64) 			
        !open(unit=64,file='stepp.txt',status='old',action='read') ; read(64,*) stepmin	; close(64) 
        !stepmin=0.5_dp

        !stepos(1:21)=0.0_dp
        !stepos(24:69)=0.0_dp 
        !optimization parameters:
        !stepmin=stepos 
        !stepmin(33:49)=0.0_dp
        !stepmin(68:npars)=0.0_dp
        !stepmin(1:12)=0.0_dp
        !stepmin=0.0_dp ; stepmin=0.0_dp
        !stepmin(14)=stepos(14)
        !stepmin(23)=stepos(23)
        !stepmin(80)=stepos(80)
        !stepmin(33:41)=stepos(33:41)
        !stepmin(68:69)=stepos(68:69)
        maxfn=8000
        iprint=20
        !iprint=-1
        stopcr=0.0001_dp*qval
        nloop=npars+2
        iquad=0 !ag092522 agsept2022: This was 0 before but I don't think that's right so I'm trying this one. actually I don't think it matters (it only matters at end after convgence)
        simp=0.0_dp
        !sim annealing parameters: 
        tstart=0.15_dp*qval  !T0: starting temp (can set to zero to turn off sim annealing)
        tstep=0.5_dp        !fraction temp reduced at each tfreq
        tfreq=COUNT(stepmin /= zero) !number of function calls between temp reductions
        saseed=1            !seed for random numbers
        !call minim(p,    step,    nop, func,  maxfn,  iprint, stopcr, nloop, iquad,  simp, var, functn, ifault)
        !call minim(pars, stepmin, npars, qval, maxfn, iprint, stopcr, nloop, iquad,  simp, var, objfunc, writebest, ifault)        
    	if (groups) then 
			call pminim(pars, stepmin, npars, qval, maxfn, iprint, stopcr, nloop, iquad, simp, var, objfunc, writebest,ifault,mygroup,numgroup,tstart, tstep	, tfreq, saseed)
		else 
			call pminim(pars, stepmin, npars, qval, maxfn, iprint, stopcr, nloop, iquad, simp, var, objfunc, writebest,ifault,iam,numworld,tstart, tstep	, tfreq, saseed)
		end if 

   
		if (iam==0) then ; print*, "out of minim now and here is ifault ", ifault ; end if 	
        if (iam==0) then  
            close(6538)
        end if 
    else 
        if (getstderr) then 
            D=-99. ; WD=-99. ; QWD=-99. ; DpWDinv=-99. ; SEmat=-99.
            dtheta=-99. ; pars1=-99. ; realpars=-99. ; realpars1=-99.
            parvo1=-99. ; realparvo1=-99.
            call getpars(pars,realpars)
            call objfunc(pars,val) ; realpars=realpartemp
            if (iam==0) print*, "Just ran iter", iter
            stepstderr=stepos !0.0_dp
            !Got the following: alphaed(m,ed)   not identified 12.84093 2.00000  29     0.11954     0.11954
            !Got the following: alphakid(m)     not identified  0.02182 2.00000  31     0.11954     0.11954
            !so ignore those for now since I'm identifying the problematic parameters anyway
            !stepstderr(1:26)=  stepos(1:26)     !ok  1:26 33:49
            !stepstderr(33:49)=  stepos(33:49)   !ok 1:26 33:49
            !stepstderr(50:63)=  stepos(50:63)   !ok 50:63 and also together with 1:26 and 33:49
            !stepstderr(64:69)=  stepos(64:69)   !ok 
            !stepstderr(70:npars)=  stepos(70:npars)   !ok
            !stepstderr=zero
            !stepstderr(14:15)=stepos(14:15)
            !stepstderr(24)=stepos(24)
            !stepstderr(65)=stepos(65)
            !stepstderr(66:67)=stepos(66:67)
            !stepstderr(75)=stepos(75)
            !stepstderr(81)=stepos(81)
            !stepstderr(87)=stepos(87)
            !stepstderr(93)=stepos(93)
            !stepstderr=0.0_dp 
            !stepstderr(1:69)=stepos(1:69)
            !stepstderr(1:19)=0.0_dp
            !stepstderr(70:npars)=0.0_dp
            nactive = COUNT(abs(stepstderr) > 0)
            if (iam==0) print*, "Here is nactive", nactive
            !if (writestderr) OPEN(13,FILE='stderrors'//runid//'.txt',STATUS='REPLACE')    
            !do im=1,nmom
            !    if (cntdat_save(im,1)>0) then
            !        QQ(im)=vardat_save(im,1)/real(cntdat_save(im,1))   !numperdat
            !    else 
            !        QQ(im)=0.0_dp
            !    end if
            !end do 
            !q1val(j)=val
            QQ=vardat_save(:,1)/numperdat
            kactive=0 ! count of active parameters
write(13,'(2x,"iter",tr2,"itm1",15x,"parname", tr2,"kk", tr1,  tr1,"ka",tr1,    tr3,"pars(kk)",  tr2,"pars1(kk)",     tr2,"realp(kk)",   tr1,"realp1(kk)",  tr5,"step",   tr3,"dtheta",   tr6,"val",   tr5,"val1",9x)') 
            DO kk=1,npars
                IF (abs(stepstderr(kk)) > 0.0_dp)  THEN
                    pars1=pars
                    pars1(kk)=pars1(kk)+stepstderr(kk)
                    if (iam==0) print*, "About to do iter for ", kk, parname(kk)
                    call getpars(pars1,realpars1)
                    call objfunc(pars1,val1) ; realpars1=realpartemp ; itermin1=iter-1
                    if (iam==0) print*, "Just ran iter", iter   , val,val1    
                    dtheta(kk)=realpars1(kk)-realpars(kk)
                    realparvo1(:,kk)=realpars1(:); parvo1(:,kk)=pars1(:)
IF (writestderr.and.(MAXVAL(ABS(QQ*momwgt*msm_wgt*(momsim_save(:,itermin1)-momsim_save(:,1))))==0) )  THEN
    if (mysay==0) print*,  "not ident ", iter,    itermin1,  parname(kk),  kk,  kactive,     pars(kk),       pars1(kk),          realpars(kk),realpars1(kk),stepstderr(kk),dtheta(kk),val,val1
    WRITE(13,'(2x,I4,     2x,I4,      1A22,        I4,    I4,        F11.2,           F11.2,               F11.2,         F11.2,    4f9.2,                                1A9)') &
    &            iter,    itermin1,  parname(kk),  kk,  kactive,     pars(kk),       pars1(kk),          realpars(kk),realpars1(kk),stepstderr(kk),dtheta(kk),val,val1,' notident'
    STOP
ELSE
    kactive=kactive+1 ! number of parameters actually iterating on.   
    activepari(kactive)=kk ! indexes of active parameters
    D(:,kactive)=(momsim_save(:,itermin1)-momsim_save(:,1))/dtheta(kk) ! derivative of moments wrt paramter
    WD(:,kactive)=momwgt*msm_wgt*D(:,kactive)
    QWD(:,kactive)=QQ*WD(:,kactive)
    if (writestderr) then
    WRITE(13,'(2x,I4,    2x,I4,      1A22,         I4,    I4,        F11.2,           F11.2,               F11.2,         F11.2,    4f9.2,                                9x)') &
     &           iter,    itermin1,  parname(kk),  kk,  kactive,     pars(kk),       pars1(kk),          realpars(kk),realpars1(kk),stepstderr(kk),dtheta(kk),val,val1
        !DO im=1,nmom
        !    WRITE(13,'(I4,11F12.5)') im,momsim_save(im,itermin1),momsim_save(im,1),momsim_save(im,itermin1)-momsim_save(im,1),D(im,1:8)
        !ENDDO        
    end if 
ENDIF
                ENDIF
            ENDDO
            if (kactive.ne.nactive) then ; print*, "Something wrong with kactive",kactive,nactive ; stop ; end if
            ! kactive now holds total number of active parameters
            !Q_MSM=SUM(weights*MSM_weights*(datamoments-simmoments)**2)
            !FINDInv(matrix, inverse, n, errorflag)
            CALL FINDinv(MATMUL(TRANSPOSE(WD(:,1:kactive)),D(:,1:kactive)),DpWDinv(1:kactive,1:kactive),kactive,eflag) ! get (D'WD)^-1
            ! Standard errors are sqrt of diagional elements of matrix ((D'WD)^-1 * D'WQWD * (D'WD)^-1)
            SEmat(1:kactive,1:kactive)=MATMUL(MATMUL(DpWDinv(1:kactive,1:kactive),TRANSPOSE(WD(:,1:kactive))),MATMUL(QWD(:,1:kactive),DpWDinv(1:kactive,1:kactive)))
            stderrs=-1. ! default value to indicate that the paramter was not estimated
            DO ka=1,kactive
                stderrs(activepari(ka))=SQRT(SEmat(ka,ka))
            ENDDO
            ! TESTING
            if (writestderr) then
                WRITE(13,*) 'standard errors:'
                DO ka=1,kactive
                    WRITE(13,'(1A15,6F12.5,1A5,1F12.5,1A1)') parname(activepari(ka)),pars(activepari(ka)),realpars(activepari(ka)),&
                    & parvo1(activepari(ka),activepari(ka)),realparvo1(activepari(ka),activepari(ka)),dtheta(activepari(ka)),realpars(activepari(ka)),'    (',stderrs(activepari(ka)),')'
                ENDDO    
                WRITE(13,*) 'D'
                DO im=1,nmom
                    WRITE(13,'(I4,8F12.5)') im,D(im,1:8)
                ENDDO
                WRITE(13,*)
                WRITE(13,*) 'QWD,WD,momwgt,msm_wgt,momwgt*msm_wgt,cntdat,numperdat,vardatinv' 
                DO im=1,nmom
                    !WRITE(*,'(2F12.5)') WD(im,1:2)
                    !WRITE(*,'(4F12.5)') WD(im,1),weights(im), MSM_weights(im),D(im,1)
                    WRITE(13,'(5F14.5,2I8,F14.5)') QWD(im,1),WD(im,1),momwgt(im),msm_wgt(im),momwgt(im)*msm_wgt(im),cntdat_save(im,1),numperdat,vardat_save(im,1)**(-1) ! weighting vector (main diagonal of diaginal weighting matrix)
                ENDDO
                WRITE(13,*)
                WRITE(13,*) 'DpWDinv(kk,kk)'
                DO ka=1,kactive
                    WRITE(13,'(1F20.5)') DpWDinv(ka,ka)
                ENDDO
                WRITE(13,*)
                WRITE(13,*) 'SEmat'
                DO ka=1,kactive
                    WRITE(13,'(1F20.5)') SEmat(ka,ka)
                ENDDO
                WRITE(13,*)
                CLOSE(13)
            end if
        end if 
    
        !call random_seed (size=p)
	    !p=12
	    !call random_seed (put = (/(k,k=1,p)/))
	    !if (skriv) then
	    !    print*, 'Here is p',p
	    !    print*, 'Here is k',(/(k,k=1,p)/)
     !			call random_number(rand)
	!		print*, 'rand: ', rand
    ! 			call random_number(rand)
!			print*, 'rand: ', rand
 !       end if
	!else     
	!    allocate(newseed(p))
	!    newseed=10
	!    call random_seed( put=newseed(1:p)  )
	!    deallocate(newseed)
	!end if 
!		call objfunc(pars,q) ; realpars=realpartemp
 
        
        !open(unit=2,file='bestpar110516.txt',status='old',action='read') ; read(2,*) pars2	; close(2) 		
        !pars(4)=-13.7901
		!call objfunc(pars,q) ; realpars=realpartemp
        if (comparepars) then
            if (iam==0) then 
                write(64,*)
                write(64,*)
                write(64,'(tr3,"i",tr3,"j",tr2,"iter",2x,tr20,tr13,"q",tr9,"q1(j)",tr6,"pars",tr5,"pars1",5x,tr6,"realpars",5x,tr5,"realpars1",5x,tr3,"step(j)",5x,tr2,"abs(q1(j)-q)")')
            end if
            do i=1,npars
                !if (abs(stepmin(i)) >0.0_dp) 
                    pars1=pars
                    pars1(i)=pars2(i)
                    call objfunc(pars1,val) ; realpars1=realpartemp
                    j=1 ; q1val(j)=val !j is just a nuissance index here which is kept just to be consistent with the chkstep output
                    if (iam==0) then 
                    	write(64,'(2i4,i6,". ",1a20,2f14.2,2F10.2,2(5x,f14.2),(5x,f10.2),(5x,f14.2) )') i,j,iter-1,parname(i),qval,q1val(j),pars(i),pars1(i),realpars(i),realpars1(i),step(j),abs(q1val(j)-qval)
                    end if 				
                !end if 
            end do 
        end if 
        
		if (chkstep) then !don't forget to set stepmin equal to stepos (i.e. assign it some value! with getpars!)
            lb=0.4_sp ; ub=1.25_sp ; incr=0.2_dp
            do i=1,npars
                !if (  parname(i)=='cst(1)'.or.parname(i)=='cst(2)'.or.parname(i)=='ecst' ) then
                !    stepmin(i)=0.2_dp
                !else 
                !    stepmin(i)=0.0_dp
                !end if 
                if ( abs(stepmin(i)) > 0.0_dp ) then	
                !if (  parname(i)=='uloc'.or.  & 
                 ! & parname(i)=='alf13'.or.parname(i)=='alf22'.or.i==47.or.i==48 .or. parname(i)=='mu_o') then
                    !if (i==4) then !i==40.or.i==41.or.i==42.or.i==50.or.i==53.or.i==54.or.i==56.or.i==57.or.i==59.or.i==60.or.i==62.or.i==63.or.i>=65) then
					if (iam==0) then 
						write(64,*)
						write(64,*)
						write(64,'(tr3,"i",tr3,"j",tr2,"iter",2x,tr20,tr13,"q",tr9,"q1(j)",tr6,"pars",tr5,"pars1",5x,tr6,"realpars",5x,tr5,"realpars1",5x,tr3,"step(j)",5x,tr2,"abs(q1(j)-q)")')
					end if 
					step=pen 
					q1val=pen 
					k=0
					incr=0.2_dp
					jloop: do j=1,nj2
						if (parname(i)=='cst(1)'.or.parname(i)=='cst(2)') then
                            step(j)=(j-nj)*0.2_dp
						else if (parname(i)=='alf13') then
                            step(j)=(j-nj)*0.2_dp*abs(pars(i))
						else if (parname(i)=='alf22') then
                            step(j)=(j-nj)*0.2_dp*abs(pars(i))
                        !else if (i>=50.and.i<=67) then
						!	step(j)=(j-nj)*0.25_dp   
                        else 
                            step(j)=(j-nj)*incr*abs(pars(i))
						end if 						
						pars1=pars						
						pars1(i)=pars(i)+step(j)
						call objfunc(pars1,val) ; realpars1=realpartemp
						q1val(j)=val
						!if ( abs(q1(j)-q) < incr) incr=0.1_dp
						if (iam==0) then 
							write(64,'(2i4,i6,". ",1a20,2f14.2,2F10.2,2(5x,f14.6),(5x,f10.2),(5x,f14.2) )') i,j,iter-1,parname(i),qval,q1val(j),pars(i),pars1(i),realpars(i),realpars1(i),step(j),abs(q1val(j)-qval)
						end if 												
					end do jloop
					call sort2(q1val,step)
					j=0
					do while (lb+j*0.1_sp<=0.95_sp)
						k=locate(q1val,(lb+j*0.1_sp)*real(qval)) 
						if (iam==0) write(64,'("first try ",i6,3f14.4,i6)') j,q1val(1),q1val(nj2),(lb+j*0.1_sp)*real(qval),k
						if (k>0.and.k<nj2) exit 
						j=j+1
					end do
					j=0
					if (k==0) then 
						do while (ub-j*0.1_sp>=1.05_sp)
							k=locate(q1val,(ub-j*0.1_sp)*real(qval)) 
							if (iam==0) write(64,'("second try ",i6,3f14.4,i6)') j,q1val(1),q1val(nj2),(ub-j*0.1_sp)*real(qval),k
							if (k>0.and.k<nj2) exit 
							j=j+1
						end do 
					end if 					
					if (iam==0) write(64,'("and now it is done ",i6)') k
					if (k>0) then ; stepmin(i)=step(k) ; else ; stepmin(i)=-99.0 ; end if 
					if (iam==0) then 
						!write(64,*) "here is the best bumps for the initial simplex"
						write(64,*)
						write(64,'("lb,ub,q1min,q1max",4x,4f14.4)') lb*qval,ub*qval,q1val(1),q1val(nj2)
						write(64,*)
						if (k>0) then 
							write(64,'(2i6,4f14.4,4(5x,f14.4) )') i,k,qval,q1val(k),pars(i),pars1(i),realpars(i),realpars1(i),step(k),abs(q1val(k)-qval)
						end if 
						write(64,*) 
					end if 												

				end if 	
				if (iam==0) write(65,'(f14.4)') stepmin(i)
			end do 
		end if 
	end if 





	runtime=secnds(begintime)
	if (iam==0) then ; write(*,'("run time: ", f14.2)') runtime ; end if 
	!ag090122 agsept2022 if (groups) then 
	!ag090122 agsept2022 	call mpi_finalize(ierror)   
	!ag090122 agsept2022 end if 
    !ag090122 agsept2022 commenting out the aove groups if statment because 
    !I get mpi error when I run with groups FALSE saying some process couldnt finish because mpi finalze wasn't called etc.
    !since I invoke mpi regardless of groups, I have to call mpi finalize regardless of groups
    deallocate(init)
    deallocate(nummove_save,nummovemar_save,nummovesin_save)

    
    !*************************    
    !mytime(iam+1,1)=secnds(0.0)
    !mytime(iam+1,2)=secnds(mytime(iam+1,1))        
    !if (iam==0) print*, 'Here is qval: ', qval
    !if (iam==0) print*, 'iam,qval,mytime ', iam,qval,mytime(iam+1,2)
    deallocate(mytime)
    call mpi_finalize(ierror)   
end program main 
