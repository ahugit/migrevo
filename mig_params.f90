	module params
	!insert space and keep tabs option in visual studio
    ! testing
	use nrtype 
	use nr 
	use alib, only: one,logit,logitinv,min2pls,min2plsinv,lin2ndim,ndim2lin,sqrtpi,multinom,condmom
	implicit none 
    include 'mpif.h'
    integer :: mysay,mygroup
	integer :: iter,comm,iwritegen 
    real(dp), parameter ::  delta=0.93_dp  !set in main =0.96_dp
    !real(dp) ::  delta   !=0.93_dp  !set in main =0.96_dp
	!real(dp) :: one=1.0_dp
	!integer(i4b), parameter :: rp = kind(1.0d0)			! kind(1.0) !!!
    integer(i4b) :: taxrates !assigned in main. determines whether taxes should be read from marginal or from average. 
    integer(i4b) :: policytax,policyubenefit !set in main to determine whether running policy experiment and which one
    integer(i4b) :: ntermval !set in main for now
    logical :: nonneg,terminalval,JOINTMARSHOCK,JOINTDIST,UPFRONTTRANSFERS,NOMISSINGOBS
    real(dp), dimension(2) :: nonlabinc !=(/ 0.0_dp,0.0_dp /) !(/ 300.0_dp,1100.0_dp /) !ahu summer18 051418: changing it back to parameter and changing dimension to 2 (not educ and educ) !ahu summer18 042318 changing this so it is set at main again
    integer(i4b), parameter :: nl=9,ndecile=15
	integer(i4b), parameter :: np=3 !ag090122 agsept2022 changed nz frmo 1 to 5 !ahu 121818 changed from 3 to 6 !ahu 0327 changed np from 5 to 2
	integer(i4b), parameter :: ncs = 3 !nl+2 !3 !nl+2 !3 !ahu october2022 HUGE MAJOR CHANGE HUGE MAJOR CHANGE nl+2
	integer(i4b), parameter :: nc  = 9 !nl+9 !9 !nl+12 !9   !nl+8 !ahu october2022 HUGE MAJOR CHANGE HUGE MAJOR CHANGE
    integer(i4b), parameter :: nepsmove=3,nz=3 !ag090122 agsept2022 changed nepsmove frmo 2 to 5 !ahumarch1022 changed nepsmove to 2 from 13
    integer(i4b) :: METHODWAGEDISCRETIZE !all the bad ideas
    integer(i4b), parameter :: numit    = 3
	integer(i4b), parameter :: nmom     = 20200 !1500 !4350 !25200 !50550 !900 !16000 !900 !16000 !900 !15600 ! 750 !2525 !2400 !815 !2950 !760 !1400 !600 !1270 !100 !6500 !605 !305 1800 !ahu summer18 050418: changed from 4200 to 498
    integer(i4b), parameter :: DETAILEDOUTPUT=10 !When this is 1, then getmom writes detailed output (all moments by type etc.), when it's large it writes normal  
    logical, parameter :: groups=.true.,onlysingles=.FALSE.,onlymales=.false.,onlyfem=.false.,choicesetlarge=.FALSE.
    logical, parameter :: optimize=.FALSE.,chkstep=.false.,chkobj=.true.,condmomcompare=.false.,comparepars=.false.
    logical, parameter :: typemoments=.FALSE.,momdisplay=.FALSE., savetime=.FALSE.,savetime2=.FALSE.,policyruns=.TRUE.,policyruns_mumar=.FALSE.
    logical, parameter :: getstderr=.FALSE.
    logical, parameter :: onthejobsearch=.TRUE.
	logical :: writestderr !set in main
    real(dp), parameter :: eps = 1.0d-6,zero=0.0_dp,epstest=2.0_dp					! could do tiny(.) but that gives a number that is way too small and therefore not appropriate for those places where we check the inequalities or equalities	
	real(dp), parameter :: eps2= 1.0d-6
	integer(i4b), parameter :: nhome=1,nhomep=nl
	logical :: conditional_moments		! can set this in main
	logical :: skriv,yaz,insol,yazmax 
    logical, parameter :: runbig=.false.  
	character(len=1), parameter :: runid='r'		! string inserted into output filenames to identify which run !ahu 062413 set this in main instead 
	integer(i4b), parameter :: nco=1,ncop=1
	integer(i4b), parameter :: ntyp=1,ntypp=4   ! types !ahu030622 changed ntypp to 1 (was 4)
	integer(i4b), parameter :: nin  = nco * ntyp * nhome
	integer(i4b), parameter :: ninp = ncop * ntypp * nhomep
    integer(i4b) :: nindex !determined in objf according to groups, it's either nin or ninp
	integer(i4b) :: iwritemom,myhome,mytyp,myco,myindex,mygrank
	logical, parameter :: indsimplexwrite=.false.		! in optimiation with parallel simplex, re-solve model and write moments to file whenever check for,find new best point 
	logical, parameter :: write_welfare=.false.			! write things needed for welfare calculations
	logical, parameter :: grid_usage=.false.			! keep track of how much of each grid is being used in simulations
	logical, parameter :: icheck_eqvmvf=.false.,icheck_eqvcvs=.false.,icheck_probs=.false.
	integer(i4b), parameter :: npars    = 101
    character(len=15), dimension(npars) :: parname ! names of each parameter   !ahu 121118 now declkare as global it here instead of getsteps
    character(len=25), dimension(nl) :: locname
    real(dp), dimension(npars) :: stepmin,realpartemp,parsforcheck,stepos !ahu 121118
	!character(len=15), dimension(npars) :: parname 
	integer(i4b), parameter :: mna=18,mxa=50   !,mxai=50		!ahu 070713 40	!if you ever change mina from 16 that might be problematic, because of psiddata%home and simdata%home definitions. look in read_data and read simdata for this
    integer(i4b), parameter :: MNAD=MNA-1,MXAD=MXA-1            !ahu jan19 010219
    integer(i4b), parameter :: nh=2,nexp=2,nsimeach=10,neduc=2,nkid=2 !kid is 1 if no kid,2 if yes kid !ahu 0327 changed nsimeach from 10 to 5
	integer(i4b), parameter :: nqs = (np+2) *  nl
	integer(i4b), parameter :: nq  = (np+2) * (np+2) * nl
	integer(i4b), parameter :: np1=np+1, np2=np+2   !w=np1 is getting laid off in the shock space q and w=np2 is nothing happening. In the state space q, np1 is unemployment and np2 is not in the state space q (I try to check that it is not there, at various points in code)
	integer(i4b), parameter :: nxs = neduc * nexp * nkid
	integer(i4b), parameter :: nx  = nxs * nxs
	integer(i4b) :: numperdat,numperobsdat,numpersim !previously ndata,ndataobs,nsim ALL set in main now 
	!integer(i4b), parameter :: ndataobs = 84507  set in main now 
	!integer(i4b), parameter :: nsim     = ndata*nsimeach  set in main now
    integer(i4b), dimension(nmom) :: calcvar,calcorr
    integer(i4b), dimension(nmom,11) :: momwhich !for momdisplay purposes
	integer(i4b), parameter :: maxrellength=10
	integer(i4b), parameter :: namelen=55					!if you change this, don't forget to also change a100 in writemoments	
	integer(i4b), parameter :: ma=1,fe=2
    INTEGER(I4B), PARAMETER :: NOCOLLEGE=1,COLLEGE=2
	integer(i4b), parameter, dimension(2) :: agestart=(/ mna,20 /)		!changed this from 18,22 !chanage this back ahu 070312 (/18,22/) !starting age for simulations for each education level
	real(dp), parameter :: mult1=10000.0_dp !ahu jan19 012519
    real(dp), parameter :: multmar=50000.0_dp,multsigo=300000.0_dp,multdiv=100000.0_dp,multcst=30000.0_dp !ahu jan19 012019  !ahu030622 VERY IMPORTANT CHANGE MULTMAR
    real(dp), parameter :: maxhgrid=8.0_dp 
	real(dp), parameter :: tottime=16.0_dp
	real(dp), parameter :: hhours_conv=250.0_dp					! multuiply hours per day by this to get hours per year for hours worked
	real(dp), parameter :: maxh=hhours_conv*maxhgrid				! hhours_conv*hmgrid(ntwork) ! truncation point of labor market hours in data !ahu 071712
	real(dp), parameter, dimension(nh) :: hgrid=(/ 0.0_dp,maxhgrid /) 
	real(dp), parameter, dimension(nh) :: hgrid_ann=hgrid*hhours_conv
	real(dp), parameter :: d1=1.0_dp						! note that now doing calculations as single real(sp) not double real (i modified module1 to allow this)
	real(dp), parameter :: h_parttime=1000.0_dp					! number of hours to be considred part time !ahu 071712
	real(dp), parameter :: h_fulltime=1000.0_dp					! ahu 071712 changed from 1000 too 2000 !ahu 062812 changed this to 1000. 2000. ! number of hours to be considred full time 
        !ahu 021817: note that the employment is decided according to whether hours is more than h_fulltime which is 1000 but annual wages is calculated using wge*h_wmult where h_wmult is 2000. 
        !This is how you did it in the original version. I do the same thing so that the wage numbers are consistent with the previous version. See page 15 last paragraph in the original text.	
    real(dp), parameter :: h_wmult=2000.0_dp                    !what you multiply hourly wages with in order to turn them into annual wages
    real(dp), parameter :: hbar=h_parttime
	real(dp), parameter :: minw=1.0_dp					! lower truncation point of male log wage
	real(dp), parameter :: maxw=150.0_dp                ! upper truncation point of male log wage
	real(dp), parameter :: pen=-99999999.0_dp
	integer(i4b), parameter :: ipen=-99999	, MALES=1, FEMALES=2
    !real(dp), parameter :: wwage=60.0_dp,wwvar=10000.0_dp,wdifww=1.0_dp,whour=10.0_dp,wcorr=1000.0_dp,wwdecile=30.0_dp,wmovenum=0.0_dp , wmovemar=10.0_dp, wmovesin=10.0_dp , wmoveall=10.0_dp
    !real(dp), parameter :: wrel=1.0_dp,whome=10.0_dp,wkid=1.0_dp,wprop=10.0_dp,wswitch=1.0_dp,wwagejoint=50.0_dp 	
    real(dp), parameter :: wwage=30.0_dp,wwvar=10000.0_dp,wdifww=1.0_dp,whour=20.0_dp,wcorr=1000.0_dp,wwdecile=30.0_dp,wmovenum=0.0_dp , wmovemar=10.0_dp, wmovesin=10.0_dp , wmoveall=10.0_dp
    real(dp), parameter :: wrel=1.0_dp,whome=10.0_dp,wkid=1.0_dp,wprop=10.0_dp,wswitch=1.0_dp,wwagejoint=50.0_dp 		    
    
    !ahu 121918 changed wmove to 10 from 1 and changed wmovemar from 10 to 100		! weights for moments for married couples. set in objfunc.
    real(dp), parameter :: wdivstay=500.0_dp
    !real(dp), parameter :: wwage=100.0_dp,wwvar=1500.0_dp,wdifww=3.0_dp,whour=10.0_dp,wcorr=100.0_dp,wwdecile=100.0_dp,wmovenum=1.0_dp
    !real(dp), parameter :: wrel=10.0_dp,wmove=30.0_dp,whome=10.0_dp,wkid=1.0_dp,wprop=10.0_dp,wswitch=3.0_dp		!ahu 121918 changed wmove to 10 from 1 and changed wmovemar from 10 to 100		! weights for moments for married couples. set in objfunc.
    !real(dp), parameter :: wwage=1.0_dp,wwvar=1.0_dp,wdifww=0.01_dp,whour=1.0_dp,wcorr=1.0_dp,wwdecile=1.0_dp
    !real(dp), parameter :: wrel=1.0_dp,wmove=10.0_dp,whome=10.0_dp,wkid=1.0_dp,wprop=1.0_dp		!ahu 121918 changed wmove 
    !real(dp), parameter :: wtrans=100.0_dp,wwaged=10.0_dp,wdifww=100.0_dp,wrel=1.0_dp,wmove=10.0_dp,whour=1.0_dp,wwvar=100.0_dp
!    real(dp), parameter :: wwage=1.0_dp,wkid=1.0_dp,wmovemar=1.0_dp,wmovesin=1.0_dp,wwagebymove=1.0_dp		!ahu 121918 changed wmove to 10 from 1 and changed wmovemar from 10 to 100		! weights for moments for married couples. set in objfunc.
    character(len=30) :: datafilename  != 'familymigpsid.txt' ! data filename set in main now
	!character(len=23), parameter :: initcondfile= 'familymiginit111113.txt' ! filename of initial conditions
	character(len=23) :: momentfile='mom.txt'
    character(len=23) :: momentonlyfile='momonly.txt'
    integer(i4b) :: mominfo(0:5,nmom)
	!model parameters 
	!real(sp):: umdist_condf(np,np),ufdist_condm(np,np)
	real(dp), parameter :: alf=0.5_dp   !ahumarch1022 delta=0.96_dp changing to 0 to figure out the mumardecrease problem
	real(dp), parameter :: mu_wge(2)=0.0_dp
	real(dp) :: sig_wge(2,MALES:FEMALES),mu_mar(ntypp),ro,mu_o , sigo_m,sigo_f,sig_mar  !,inc_coef !kennan's income coef
	real(dp) :: uhomet(2),alphaed(2,neduc),alphakid(2),alphab(2) !ahu october2022: changing alphakid so that it doesn't have that obsolete dimension anymore (ahu jan2023: but really this should be alphakid(2) NOT alphakid(nkid) )
	real(dp) :: uloc(nl),cstned(MALES:FEMALES),csted(MALES:FEMALES),CSTU(MALES:FEMALES),cstadjacent,mu_errormeasure,sig_errormeasure,diveduc(2)
	real(dp) :: alf10(nl),alf11,alf12,alf1t(ntypp),alf13,alf14          ! types
	real(dp) :: alf20(nl),alf21,alf22,alf2t(ntypp),alf23,alf24           ! types
	real(dp) :: ptype,omega(2),ptypehs(ntypp,MALES:FEMALES),ptypecol(ntypp,MALES:FEMALES) ! types
	real(dp) :: psio(12),psil(2),psihup_e   !,psihdown_u !getting rid of psih2,3,4 so that now psih is just a scalar
    real(dp) :: pkidmar=0.07_dp,pkidsingle=0.03_dp
	real(dp) :: popsize(nl),heat(nl),ubenefitbyloc(nl), meetsame(2), alfmeet(2),corral
	integer(i4b) :: distance(nl,nl)
	real(dp) :: best
	type :: initcond
        integer(i4b) :: id
		integer(i4b) :: co				
		integer(i4b) :: sexr			
		integer(i4b) :: hme
		integer(i4b) :: endage,minage		
		integer(i4b) :: edr		
        integer(i4b) :: typ	
    end type
	type, extends(initcond) :: statevar
        integer(i4b) :: expr
        integer(i4b) :: kidr      
        integer(i4b) :: hhr		! annual hours worked by r (turned into discrete in read_data)
		real(dp) :: logwr,wr	! wr is annual income and logwr is log of annual income. hourly wage (wr_perhour or wsp_perhour) is read from the data (see familymig_2.do to see how it's turned into hourly) and that hourly wage is turned into annual by multiplying it by h_wmult	  
        integer(i4b) :: l		! location		
        integer(i4b) :: rel		! relationship status. -1: not observed, 0: single, 1: married, 2: cohabiting
		integer(i4b) :: rellen	! length of current relationship starting from 1 in first period
        integer(i4b) :: edsp
        integer(i4b) :: expsp
        integer(i4b) :: kidsp
        integer(i4b) :: hhsp	! annual hours worked by spouse (turned into discrete in read_data)
		real(dp) :: logwsp,wsp  ! see explanation for logwr,wr.  
        integer(i4b) :: lsp
        integer(i4b) :: nomiss
        integer(i4b) :: nn,mm,r,job,jobswitch
        real(dp) :: emax,emaxsp,consr,consp,incsum,consum !adding this to the statevar type but that's obviously not a statevar. just for convenience.
	end type
	type :: shock
		real(dp) :: meet 
		real(dp) :: marie
		real(dp) :: meetq 
		real(dp) :: meetx 
		real(dp) :: q
		real(dp) :: x
        real(dp) :: move
        real(dp) :: typ
	end type	
	type(statevar) :: ones
    type(initcond) :: ones_init
    type :: taxo
        real(sp) :: pwages
        real(sp) :: swages 
        real(dp) :: statesin
        real(dp) :: fedsin
        real(dp) :: statemar
        real(dp) :: fedmar
    end type 
    type :: prod
        real(dp) :: state
        real(dp) :: fed
    end type 
    type :: joint
        real(dp) :: hub
        real(dp) :: wfe
    end type
    !The below types shockdraws and decisions are defined last minute just for output purposes
	type :: writedraws
		integer(i4b) :: qshock
        integer(i4b) :: wmshock,wfshock,lshock
		integer(i4b) :: iephub,iepwfe,iepsingle
        integer(i4b) :: zshock
	end type	
	type :: writedecisions
        integer(i4b) :: q0dec,qdec,qmdec,qfdec
		integer(i4b) :: rel0dec,reldec
        integer(i4b) :: wm0dec,wf0dec
		integer(i4b) :: wmdec,wfdec
        integer(i4b) :: l0dec,ldec,lmdec,lfdec
	end type	
	type :: mdec
        integer(i4b) :: qmdec,wmdec,lmdec 
	end type	
	type :: fdec
        integer(i4b) :: qfdec,wfdec,lfdec 
	end type	
	type :: writeefficiency
        !real(dp) :: ve1,ve2,ve3,ve4,ve5,sur
        real(dp) :: tra1,tra2
        !logical :: haveenough,haveenoughforNB
        integer(i4b) :: maxA,maxB
        !integer(i4b) :: maxwmA,maxwfA,maxlmA,maxlfA
        !integer(i4b) :: maxwmB,maxwfB,maxlmB,maxlfB
	end type	
    type(writedraws) :: ones_draws
    type(writedecisions) :: ones_decisions
    type(mdec) :: ones_mdecisions
    type(fdec) :: ones_fdecisions
    type(writeefficiency) :: ones_efficiency
    integer(i4b), parameter :: numbin=31
    integer(i4b), parameter :: numtaxes=nl*numbin*numbin 
    !integer(i4b) :: taxset !set in main
    type(taxo), dimension(numbin,numbin,nl) :: tax !pwages,swages.myreg
    real(4) :: pbracket(numbin),sbracket(numbin) !because locate() needs it to be real(4). otherwise I have to define a separate variable. 
    type(prod), dimension(numbin,nl) :: bracketprod_s,bracketprodsum_s
    type(prod), dimension(numbin,numbin,nl) :: bracketprod_h,bracketprod_w,bracketprodsum_h,bracketprodsum_w
	type(initcond), dimension(:), allocatable :: init !array for initial conditions (size just num of persons in actual data)
	integer(i4b), dimension(:), allocatable :: nummove_save,nummovemar_save,nummovesin_save   
	type(writedraws), dimension(:,:), allocatable :: draws
	type(writedecisions), dimension(:,:), allocatable :: decisions
	type(mdec), dimension(:,:), allocatable :: mdecisions
	type(fdec), dimension(:,:), allocatable :: fdecisions
    type(writeefficiency) :: myefficiency
	type(writeefficiency), dimension(:,:), allocatable :: myefo
contains

	! get parameters from transformed values. input is free
	! to be any real(sp) resulting parameters are appropriately constrained
	subroutine getpars(par,realpar)
	real(dp), dimension(npars), intent(in)  :: par ! vector of parameters
	real(dp), dimension(npars), intent(out) :: realpar ! vector of parameters
	integer(i4b) :: g,i,j,ed,indust1(ntypp,MALES:FEMALES),indust2(ntypp,MALES:FEMALES),co,typ,homeloc
    integer(i4b) :: dw0,de,dsex
    real(dp) :: temprob(nl),junk,junkpar
    junk=999.0_dp
    stepos=0.0_dp
    indust1=0
    indust2=0
    temprob=0.0_dp 
	realpar=pen 
    parname=''
    j=1
    !ahu jan19 012819: not iterating on ed offers anymore. replacing them with curloc and ofloc offers instead 
    !note that realpar's for psio parameters are reassigned at the end of this file just for visual purpoes, to write those in writemoments.
    !but the actual values that are used are assigned to psio right here and those are the ones that are used in fnprof.
	realpar(j)=par(j)               ; parname(j)='emp,cur,m' ; stepos(j)=0.5_dp*par(j)  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(1)=realpar(j)	            ; j=j+1
	
    realpar(j)=par(j)               ; parname(j)='emp,cur,m' ; stepos(j)=0.5_dp*par(j)  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(2)=realpar(j)            ; j=j+1
	
    realpar(j)=par(j)               ; parname(j)='emp,of,m' ; stepos(j)=0.5_dp*par(j)  ; if (onlyfem) stepos(j)=0.0_dp!psio is for the offer function
	psio(3)=realpar(j)	            ; j=j+1
	
    realpar(j)=0.0_dp               ; parname(j)='emp,of,m' ; stepos(j)=0.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function !NO MORE OFLOC LAYOFF NONSENSE
	psio(4)=realpar(j)	            ; j=j+1
	
    realpar(j)=par(j)               ; parname(j)='emp,cur,f' ; stepos(j)=0.5_dp*par(j) ; if (onlymales) stepos(j)=0.0_dp !psio is for the offer function
	psio(5)=realpar(j)	            ; j=j+1
	
    realpar(j)=par(j)               ; parname(j)='emp,cur,f' ; stepos(j)=0.5_dp*par(j) ; if (onlymales) stepos(j)=0.0_dp  !psio is for the offer function
	psio(6)=realpar(j)	            ; j=j+1
	
    realpar(j)=par(j)               ; parname(j)='emp,of,f' ; stepos(j)=0.5_dp*par(j) ; if (onlymales) stepos(j)=0.0_dp  !psio is for the offer function
	psio(7)=realpar(j)	            ; j=j+1
	
    realpar(j)=0.0_dp               ; parname(j)='emp,of,f' ; stepos(j)=0.0_dp ; if (onlymales) stepos(j)=0.0_dp  !psio is for the offer function !NO MORE OFLOC LAYOFF NONSENSE
	psio(8)=realpar(j)	            ; j=j+1
	
    realpar(j)=par(j)               ; parname(j)='u,cur,m' ; stepos(j)=0.5_dp*par(j)  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(9)=realpar(j)	            ; j=j+1
	
    realpar(j)=par(j)               ; parname(j)='u,of,m' ; stepos(j)=0.5_dp*par(j)  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(10)=realpar(j)	            ; j=j+1
	
    realpar(j)=par(j)               ; parname(j)='u,cur,f' ; stepos(j)=0.5_dp*par(j)  ; if (onlymales) stepos(j)=0.0_dp !psio is for the offer function
	psio(11)=realpar(j)	            ; j=j+1
	
    realpar(j)=par(j)               ; parname(j)='u,of,f' ; stepos(j)=0.5_dp*par(j)  ; if (onlymales) stepos(j)=0.0_dp !psio is for the offer function
	psio(12)=realpar(j)	            ; j=j+1
    !print*, 'Here is psio12',j-1
    !note that realpar's for psio parameters are reassigned at the end of this file just for visual purpoes, to write those in writemoments.
    !but the actual values that are used are assigned to psio right here and those are the ones that are used in fnprof.

	realpar(j)=par(j)               ; parname(j)='psil(1)' ; stepos(j)=0.5_dp !0.5_dp 
	psil(1)=realpar(j)	            ; j=j+1
	realpar(j)= par(j)             ; parname(j)='uhomet men'	; stepos(j)=0.5_dp*par(j) 
	uhomet(1)=realpar(j)            ; j=j+1
	realpar(j)= par(j)            ; parname(j)='uhomet fem'	; stepos(j)=0.5_dp*par(j) !2.0_dp*(1.0_dp/(1.0_dp+exp(-par(j))))-1.0_dp 
	uhomet(2)=realpar(j)            ; j=j+1
    !if (iwritegen==1) print*, "Here is ro", ro, par(j-1)
    !note that realpar's for psih parameters are reassigned at the end of this file just for visual purpoes, to write those in writemoments.
    !but the actual values that are used are assigned to psih right here and those are the ones that are used in fnprhc.
	realpar(j)=par(j)             ; parname(j)='p(ex=2|ex=1),e' ; stepos(j)=0.5_dp !this is psih, the only one that governs fnprhc.  
	psihup_e=realpar(j)	            ; j=j+1

    realpar(j)=logit(par(j))                    ; parname(j)='sig_wge M'	; stepos(j)=0.5_dp	  ; if (onlyfem) stepos(j)=0.0_dp  !66:67
	sig_wge(1,MALES)=realpar(j)             ; j=j+1
    realpar(j)=logit(par(j))                    ; parname(j)='sig_wge M'	; stepos(j)=0.5_dp	  ; if (onlyfem) stepos(j)=0.0_dp  !66:67
	sig_wge(2,MALES)=realpar(j)             ; j=j+1

    realpar(j)=logit(par(j))              ; parname(j)='corral' ; stepos(j)=0.5_dp ; if (onlysingles) stepos(j)=0.0_dp  
	corral=realpar(j)	            ; j=j+1
    !note that realpar's for psih parameters are reassigned at the end of this file just for visual purpoes, to write those in writemoments.
    !but the actual values that are used are assigned to psih right here and those are the ones that are used in fnprhc.
    !ahu jan23: adding psih down 
    
    realpar(j)=(par(j) )       ; parname(j)='cstadjacent' ; stepos(j)=0.5_dp*par(j)  
	cstadjacent=realpar(j)                 ; j=j+1
    realpar(j)=par(j)              ; parname(j)='meetsame 1' ; stepos(j)=0.5_dp*par(j)  ; if (onlysingles) stepos(j)=0.0_dp  
	meetsame(1)=realpar(j)	            ; j=j+1


    !realpar(j) = mult1c * logit(par(j))              ; parname(j)='uhome(1)' ; stepos(j)=0.2_dp	!mult3*logit(par(2:3)) !22:23
	!uhome(1)=realpar(j)                             ; j=j+1
    !realpar(j) = mult1c * logit(par(j))              ; parname(j)='uhome(2)' ; stepos(j)=0.2_dp	 ; if (onlymales) stepos(j)=0.0_dp !mult3*logit(par(2:3)) !22:23
	!uhome(2)=realpar(j)                             ; j=j+1

    realpar(j) = par(j)             ; parname(j)='diveduc NED ' ; stepos(j)=-0.5_dp*PAR(J)	 ; if (onlysingles) stepos(j)=0.0_dp  
    DIVEDUC(1)=realpar(j)                             ; j=j+1
    realpar(j) = par(j)              ; parname(j)='diveduc ED' ; stepos(j)=-0.5_dp*PAR(J)	 ; if (onlysingles) stepos(j)=0.0_dp  
	DIVEDUC(2)=realpar(j)                    ; j=j+1

    realpar(j)=par(j)            ; parname(j)='sig_mar'	; stepos(j)=-0.5_dp*PAR(J)  ; if (onlysingles.or.nz==1) stepos(j)=0.0_dp  !*par(j) !24 !-1.0_dp*mult1c * logit(par(j)) !ahu 112718 changing to only minus from: mult1 * min2pls(par(j))     ! types
    sig_mar=realpar(j)                                     ; j=j+1               ! types
    
    !sigom and sigof: 68:69
    realpar(j)=par(j)                               ; parname(j)='CSTU M'	; stepos(j)=0.0_dp*PAR(J) ; if (onlysingles.or.nepsmove==1) stepos(j)=0.0_dp ; if (onlyfem) stepos(j)=0.0_dp
    !print*, "Here it is sigom", j,par(j),realpar(j)
    CSTU(MALES)=realpar(j)                                ; j=j+1
    
    
    realpar(j) = par(j)         ; parname(j)='alphab m'	; stepos(j)=2.0_dp*PAR(J)  !; if (onlysingles) stepos(j)=0.0_dp !26 !ahu 112718 changing to only minus from: mult1 * min2pls(par(6))                         !ahu summer18 050418: changed from 1000 to 10,000 (mult to mult1)
	alphab(MALES)=realpar(j)                               ; j=j+1
    !print*, 'Here is divpenalty',j-1,divpenalty 

    !realpar(j:j+1) = mult1 * logit(par(j:j+1))          ; parname(j)='alphaed(m,ned)' ; parname(j+1)='alphaed(f,ned)'    !27:28   !ahu jan19 012719 changing it yet again back to logit because there is not that much of different in objval between alpha=0 and alpha=-49000    !ahu jan19 012019 changing it back to min2pls  ! noed !ahu 112718 changing to only plus from: mult1*min2pls(par(7:8))   !mult1 * logit(par(7))	
    !realpar(j:j+1) = exp(par(j:j+1))          ; parname(j)='alphaed(m,ned)' ; parname(j+1)='alphaed(f,ned)'   
    !realpar(j:j+1) = exp(par(j:j+1))          ; parname(j)='alphaed(m,ed)' ; parname(j+1)='alphaed(f,ed)'    
    
    realpar(j) = par(j)             ; parname(j)='alphab f'     ;   stepos(j)=-0.5_dp*par(j)          
    alphab(FEMALES)=realpar(j)               ; j=j+1 
    realpar(j) = PAR(J)             ; parname(j)='CSTNED M'     ;   stepos(j)=0.5_dp*par(j)          
    CSTNED(MALES)=realpar(j)               ; j=j+1 
    realpar(j) = PAR(J)             ; parname(j)='CSTNED F'     ;   stepos(j)=0.5_dp*par(j)          
    CSTNED(FEMALES)=realpar(j)               ; j=j+1 
    realpar(j) = par(j)             ; parname(j)='alfmeet 1 '     ;   stepos(j)=0.5_dp*par(j)     ; if (onlysingles) stepos(j)=0.0_dp       
    alfmeet(1)=realpar(j)               ; j=j+1 
    
    !realpar(j:j+1)=mult1 * logit(par(j:j+1))            ; parname(j)='alphakid(m)' ; parname(j+1)='alphakid(f)'          !31:32           !ahu 112718 changing to only plus from: mult1 * min2pls(par(j:j+1))	 !mult1 * logit(par(9:10))	
    realpar(j)=par(j)            ; parname(j)='alphakid(m)'     ;  stepos(j)=2.0_dp*par(j)    ; if (onlyfem) stepos(j)=0.0_dp    !31:32          	    
    realpar(j+1)=par(j+1)   ; parname(j+1)='alphakid(f)'   ; stepos(j+1)=-0.5_dp*par(j)   ; if (onlymales) stepos(j:j+1)=0.0_dp        !31:32         
    alphakid(:)=realpar(j:j+1)                        ; j=j+2         
    !print*, 'Here is uloc',j
	
    !uloc: 33-41
    !do ed=1,2
        do i=1,nl
            if (i==3) then
			    realpar(j) = 0.0_dp  ; stepos(j)=0.0_dp
			    uloc(i)=0.0_dp
		    else 
			    realpar(j) = par(j) ; stepos(j)=0.5_dp*PAR(J)    !mult1 * min2pls( par(j) )
			    uloc(i)=realpar(j)
		    end if 
            parname(j)='uloc' 
            j=j+1
        end do
	!end do
    !print*, 'Here is alf10',j
	!wage 42: 65
    do i=1,nl
        if (i==7) then
            realpar(j)=0.0_dp ; stepos(j)=0.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !1.5_dp*min2pls(par(j))+8.5_dp 
            alf10(i)=realpar(j)
        else 
            realpar(j)=par(j) ; stepos(j)=0.15_dp  ; if (onlyfem) stepos(j)=0.0_dp !1.5_dp*min2pls(par(j))+8.5_dp 
            alf10(i)=realpar(j)
        end if     
        parname(j)='alf10'  
        j=j+1
    end do 
    !print*, 'Here is alf11 etc',j
    realpar(j)=logit(par(j))                        ; parname(j)='alf11' ; stepos(j)=0.5_dp  ; if (onlyfem) stepos(j)=0.0_dp
	alf11=realpar(j)                                ; j=j+1
    !print*, 'Here is alf12',j	
    realpar(j)=logit(par(j))                 ; parname(j)='alf12' ; stepos(j)=0.5_dp  ; if (onlyfem) stepos(j)=0.0_dp
    alf12=realpar(j)                                ; j=j+1
    !print*, 'Here is alf13',j	
    realpar(j)=par(j)                               ; parname(j)='psil(2)' ; stepos(j)=0.5_dp   ; if (onlyfem) stepos(j)=0.0_dp  !-1.0_dp*logit(par(j)) 
    psil(2)=realpar(j)	                            ; j=j+1
    !print*, 'Here is alf20',j 54:62
    do i=1,nl
        if (i==7) then
            realpar(j)=0.0_dp ; stepos(j)=0.0_dp  ; if (onlymales) stepos(j)=0.0_dp !1.5_dp*min2pls(par(j))+8.5_dp 
            alf20(i)=realpar(j)
        else 
            realpar(j)=par(j) ; stepos(j)=0.15_dp  ; if (onlymales) stepos(j)=0.0_dp !1.5_dp*min2pls(par(j))+8.5_dp 
            alf20(i)=realpar(j)
        end if     
        parname(j)='alf20'  
		j=j+1
	end do 
    !print*, 'Here is alf21 etc',j
	realpar(j)=logit(par(j))                        ; parname(j)='alf21' ; stepos(j)=0.5_dp  ; if (onlymales) stepos(j)=0.0_dp 
	alf21=realpar(j)                                ; j=j+1
    !print*, 'Here is alf22',j	    
    realpar(j)=logit(par(j))                 ; parname(j)='alf22' ; stepos(j)=0.5_dp  ; if (onlymales) stepos(j)=0.0_dp 
	alf22=realpar(j)                                ; j=j+1
    !print*, 'Here is alf23',j	
    realpar(j)=PAR(J)                            ; parname(j)='CSTU F' ; stepos(j)=0.5_dp*par(j) ; if (onlysingles) stepos(j)=0.0_dp  
	CSTU(FEMALES)=realpar(j)	                            ; j=j+1

    realpar(j)=logit(par(j))                    ; parname(j)='sig_wge F'	; stepos(j)=0.5_dp	  ; if (onlyfem) stepos(j)=0.0_dp  !66:67
	sig_wge(1,FEMALES)=realpar(j)             ; j=j+1
    realpar(j)=logit(par(j))                    ; parname(j)='sig_wge F'	; stepos(j)=0.5_dp	  ; if (onlyfem) stepos(j)=0.0_dp  !66:67
	sig_wge(2,FEMALES)=realpar(j)             ; j=j+1

    !sigom and sigof: 68:69
    realpar(j)=multsigo * logit(par(j))                               ; parname(j)='sigo_m/f'	; stepos(j)=0.5_dp ; if (nepsmove==1) stepos(j)=0.0_dp ; if (onlyfem) stepos(j)=0.0_dp
    !print*, "Here it is sigom", j,par(j),realpar(j)
    sigo_m=realpar(j)                                ; j=j+1
    sigo_f=sigo_m

    realpar(j)=par(j)                               ; parname(j)='alfmeet 2'	; stepos(j)=0.5_dp*par(j) ; if (onlysingles) stepos(j)=0.0_dp  
    alfmeet(2)=realpar(j)                                ; j=j+1


    do i=1,ntypp !The below are parameters 70 to 93. So 70-75 is type1, 76-81 is type2, 82-87 is type3, 88-93 is type4.
        if (i==1) then 
            realpar(j)=par(j)                           ; parname(j)='alf13' ; stepos(j)=0.005_dp 
            alf13=realpar(j)                            ;  j=j+1
            realpar(j)=par(j)                           ; parname(j)='alf14'  ; stepos(j)=-0.00005_dp
            alf14=realpar(j)                            ;  j=j+1

            realpar(j)=par(j)                           ; parname(j)='alf1t'     ; stepos(j)=0.2_dp  ; if (onlyfem) stepos(j)=0.0_dp
            alf1t(i)=realpar(j)                         ; j=j+1
            realpar(j)=par(j)                          ; parname(j)='alf2t'     ; stepos(j)=0.2_dp  ; if (onlymales) stepos(j)=0.0_dp
	        alf2t(i)=realpar(j)                         ; j=j+1
            !realpar(j)= -1.0_dp*mult1c * logit(par(j))   ; parname(j)='cst'       ; stepos(j)=0.5_dp
            !cst(i)=realpar(j)                           ; j=j+1 
            realpar(j)=par(j)                          ; parname(j)='csted M'       ; stepos(j)=0.5_dp*PAR(J) !par(j) !not iterating on this anymore. see notes. under cost vs. sigo. they are just not sep ident I think. 
            cstED(MALES)=realpar(j)                           ; j=j+1 
            !ahu082822 august2022 print*, 'mumar(1)',j,par(j),multmar, min2pls(par(j)),multmar*min2pls(par(j))
            realpar(j)=par(j)          ; parname(j)='mu_mar'     ; stepos(j)=-0.5_dp*PAR(J)    ; if (onlysingles) stepos(j)=0.0_dp 	    
            mu_mar(i)=realpar(j)                        ; j=j+1      
        else
            realpar(j)=par(j)                            ; parname(j)='ptypehs' ; stepos(j)=0.5_dp
            ptypehs(i,MALES)=exp(realpar(j))                  ; indust1(i,MALES)=j ; j=j+1
            realpar(j)=par(j)                            ; parname(j)='ptypecol'  ; stepos(j)=0.5_dp
            ptypecol(i,MALES)=exp(realpar(j))                 ; indust2(i,MALES)=j ; j=j+1
            realpar(j)=par(j)                           ; parname(j)='alf1t'     ; stepos(j)=0.2_dp  ; if (onlyfem) stepos(j)=0.0_dp
            alf1t(i)=realpar(j)                         ; j=j+1
            realpar(j)=par(j)                           ; parname(j)='alf2t'     ; stepos(j)=0.2_dp  ; if (onlymales) stepos(j)=0.0_dp
	        alf2t(i)=realpar(j)                         ; j=j+1
            !if (i<4) then
            !    realpar(j)= par(j)                         ; parname(j)='cst'       ; stepos(j)=1.0_dp*par(j)
            !    cst(i)=realpar(j)                           ; j=j+1 
            !else if (i==4) then
            if (i==2) then 
                realpar(j)= par(j)                         ; parname(j)='csted F'       ; stepos(j)=0.5_dp*PAR(J)
                csted(FEMALES)=realpar(j)                           ; j=j+1 
            else if (i==3) then 
                realpar(j)=-1.0_dp*0.05_dp*logit(par(j))               ; parname(j)='mu_errormeasure' ; stepos(j)=2.0_dp !this is just for visuals. not a parameter anymore. 
                mu_errormeasure=realpar(j)	            ; j=j+1	
            else if (i==4) then 
                realpar(j)=(-2.0*mu_errormeasure)**0.5               ; parname(j)='sig_errormeasure' ; stepos(j)=0.0_dp !based on the flinn italy paper
                sig_errormeasure=realpar(j)	            ; j=j+1                
            end if
            realpar(j)=par(j)                        ; parname(j)='mu_mar'     ; stepos(j)=-0.5_dp*PAR(J)   ; if (onlysingles) stepos(j)=0.0_dp 	    
            mu_mar(i)=realpar(j)                        ; j=j+1      
        end if 
    	!print*, 'Here is cost',cst(i),par(j-1),realpar(j-1)
        !ahu 122818 changed mult1 to multmar 
    end do 
    !ptypehs(1,FEMALES)=0.0_dp
    !ptypecol(1,FEMALES)=0.0_dp
    !70:71 76:77 82:83 88:89 and now 94,95   94:95=70:71  96:97=76:77  98:99=82:83    100:101=88:89
    realpar(j)=par(j)                           ; parname(j)='alf23' ; stepos(j)=0.005_dp 
    alf23=realpar(j)                            ;  j=j+1
    realpar(j)=par(j)                           ; parname(j)='alf24'  ; stepos(j)=-0.00005_dp
    alf24=realpar(j)                            ;  j=j+1
    do i=2,ntypp
        realpar(j)=par(j)                            ; parname(j)='ptypehs F' ; stepos(j)=0.5_dp
        ptypehs(i,FEMALES)=exp(realpar(j))                  ; indust1(i,FEMALES)=j ; j=j+1
        realpar(j)=par(j)                            ; parname(j)='ptypecol F'  ; stepos(j)=0.5_dp
        ptypecol(i,FEMALES)=exp(realpar(j))                 ; indust2(i,FEMALES)=j ; j=j+1
    end do 
    !print*, "pars70,71 and 94:95", pars(70:71),pars(94:95)
    !print*, "pars76,77 and 96:97", pars(76:77),pars(96:97)
    !print*, "pars82,83 and 98:99", pars(76:77),pars(96:97)
    !print*, "pars88,89 and 100:101", pars(76:77),pars(96:97)


    ptypehs(1,MALES)=exp(0.0_dp)                
    ptypecol(1,MALES)=exp(0.0_dp)               
    ptypehs(1,FEMALES)=exp(0.0_dp)              
    ptypecol(1,FEMALES)=exp(0.0_dp)             

    ptypehs(:,MALES)=ptypehs(:,MALES)/sum(ptypehs(:,MALES))
    ptypecol(:,MALES)=ptypecol(:,MALES)/sum(ptypecol(:,MALES))
    do i=2,ntypp
        realpar(indust1(i,MALES))=ptypehs(i,MALES)
        realpar(indust2(i,MALES))=ptypecol(i,MALES)
    end do 
    ptypehs(:,FEMALES)=ptypehs(:,FEMALES)/sum(ptypehs(:,FEMALES))
    ptypecol(:,FEMALES)=ptypecol(:,FEMALES)/sum(ptypecol(:,FEMALES))
    do i=2,ntypp
        realpar(indust1(i,FEMALES))=ptypehs(i,FEMALES)
        realpar(indust2(i,FEMALES))=ptypecol(i,FEMALES)
    end do 

    !print*, "here it is PTYPE FOR FEMALES",  ptypehs(:,FEMALES)

    !ptypehs(:)=ptypehs(:)/sum(ptypehs)
    !ptypecol(:)=ptypecol(:)/sum(ptypecol)
    !do i=1,ntypp
    !    realpar(indust1(i))=ptypehs(i)
    !    realpar(indust2(i))=ptypecol(i)
    !end do 

    !******************** ahu october2022 **********************************************
    ! writing of realpar for psio parameters in getpars by calling fnprof 
    ! writing of realpar for psih parameters in getpars by calling fnprof 
    ! note that the below is just for visual purposes. realpar is not used in the solution or simulation. the actual parameters psio are used. 
    ! and here I am not reassigning values of psio. just reassignign values of realpar so I can write it in writemoments. 
    ! see writemoments file for more detailed writing of the fnprof parameters. 
    temprob(1:3)=fnprof(np,5,1) !emp curloc m
    realpar(1:2)=temprob(1:2) 
    temprob(1:3)=fnprof(np,10,1) !emp ofloc m
    realpar(3:4)=temprob(1:2) 
    temprob(1:3)=fnprof(np,5,2) !emp curloc f 
    realpar(5:6)=temprob(1:2) 
    temprob(1:3)=fnprof(np,10,2) !emp ofloc f
    realpar(7:8)=temprob(1:2) 
    temprob(1:3)=fnprof(np1,5,1) !unemp curloc m
    realpar(9)=temprob(1) 
    temprob(1:3)=fnprof(np1,10,1) !unemp ofloc m
    realpar(10)=temprob(1) 
    temprob(1:3)=fnprof(np1,5,2) !unemp curloc f 
    realpar(11)=temprob(1) 
    temprob(1:3)=fnprof(np1,10,2) !unemp ofloc f
    realpar(12)=temprob(1) 
    !if (nexp>2) then ; print*, "it is only ok to write this way if nexp is 2! so beware" ; stop ; end if
    temprob(1:nexp)=fnprhc(1,np) !when experience is 1 and when working
    !if (mysay==0) then
    !    print*, "fnprhc(1,np)", fnprhc(1,np)
    !    print*, "fnprhc(2,np)", fnprhc(2,np)
    !    print*, "fnprhc(3,np)", fnprhc(3,np)
    !    print*, "fnprhc(1,np1)", fnprhc(1,np1)
    !    print*, "fnprhc(2,np1)", fnprhc(2,np1)
    !    print*, "fnprhc(3,np1)", fnprhc(3,np1)
    !end if  
    !realpar(16)=temprob(2) !this is prob of moving to experience=2 when your experience is 1. !
    !realpar(17)=temprob(1) !this is prob of moving to experience=1 when your experience is 1. 
    !temprob(1:nexp)=fnprhc(nexp,np1) !when experience is 2 and when not working
    !realpar(18)=temprob(2) !this is prob of moving to experience=2 when your experience is 2.  
    !realpar(19)=temprob(1) !this is prob of moving to experience=1 when your experience is 2.
    !note that we are not writing the fnprhc(.,np1) because that is just prob of staying where you are is 1. Check this in the writing of fnprhc in writemoments. 
   
    !******************** ahu october2022 **********************************************
    temprob(1:nl)=fnprloc(1,1) !if origin location is loc1 and homeloc is 2 (cuz index is 2), what is the probability of drawing location 1 
    !realpar(13)=temprob(1)
    if (mysay==0.and.skriv) THEN
        call index2cotyphome(1,co,typ,homeloc)
        print*, "here is co,typ,homeloc",co,typ,homeloc
        print*, "fnprloc(1,1)"
        print*, temprob(:)
    end if 
    temprob(1:nl)=fnprloc(1,5) !if origin location is loc1 and homeloc is 2(cuz index is 2), what is the probability of drawing location 2 if index is 1 (so you are at your homeloc)
    !realpar(53)=temprob(2)
    if (mysay==0.and.skriv) THEN
        call index2cotyphome(5,co,typ,homeloc)
        print*, "here is co,typ,homeloc",co,typ,homeloc
        print*, "fnprloc(1,5)"
        print*, temprob(:)
        call index2cotyphome(9,co,typ,homeloc)
        print*, "here is co,typ,homeloc",co,typ,homeloc
        print*, "fnprloc(1,9)"
        print*, fnprloc(1,9)
    end if 
    !******************** ahu october2022 **********************************************
    !Assign Cendiv names 
    locname(1)='New England        '
    locname(2)='Mid Atlantic       '
    locname(3)='East North Central '
    locname(4)='West North Central '
    locname(5)='South Atlantic     '
    locname(6)='East South Central '
    locname(7)='West South Central '
    locname(8)='Mountain           '
    locname(9)='Pacific            '
    !******************** ahu october2022 **********************************************

    !alphaed(:,2)=alphaed(:,1)
    mu_o=0.0_dp
    ro=0.0_dp !no longer a parameter 

    !***********************
    !ahu 041118 del and remove later:
    !alphaed(2,:)=alphaed(1,:)
    !psio(5:8)=psio(1:4)
    !psio(11:12)=psio(9:10)
    !alphakid(2,:)=alphakid(1,:)
    
    !alf20=alf10
    !alf21=alf11
    !alf22=alf12
    !alf23=alf13
    !uhome(2)=uhome(1)
    !***********************
    
    !psio(1:4)=psio(5:8)
    !psio(9:10)=psio(11:12)
    
    !pkid=0.0_dp
    !alphakid=0.0_dp
    !kcst=0.0_dp 
    !ro=0.0_dp
    !scst=0.0_dp
    
    !alf11=0.0_dp
    !alf21=0.0_dp
    !psio(3:4)=psio(1:2)
    !psio(7:8)=psio(5:6)
    !psio(10)=psio(9)
    !psio(12)=psio(11)
    !uloc(:,2)=uloc(:,1)
    !alphaed(:,2)=alphaed(:,1)
    !alf1t(2)=0.0_dp
    !alf2t(2)=0.0_dp
    !ptypecol=ptypehs
    
    
    !cst=0.0_dp
    !kcst=0.0_dp
    
    !ro=0.0_dp !0.98_dp
    !alpha=0.0_dp
    !kcst=0.0_dp
    !pkid=0.0_dp
    !alf1t(2)=alf1t(1)
    !alf2t(2)=alf2t(1)
    !cst(1)=cst(2)
    !uhome=0.0_dp

    !uhome=0.0_dp
    !cst=-150000.0_dp
    !kcst=-150000.0_dp
    !divpenalty=0.0_dp
    !pkid=0.0_dp
    !kcst=0.0_dp
    !alpha(:,2)=alpha(:,1)
    
    !sig_o=sig_mar    
    !alpha(1,:)=alpha(2,:) 
    	!if ((.not.optimize).and.(.not.chkstep) ) print*, "sig_mar,mu_mar,npars;,; ", sig_mar,mu_mar,j
    
	!print*, logitinv(alf11),logitinv(alf12),logitinv(alf13)
	!print*, logitinv(alf21),logitinv(alf22),logitinv(alf23)

	!if (j/=npars) then ; print*, "something wrong in getpar! ",j,npars ; stop ; end if
    

    
        realpartemp=realpar
	end subroutine getpars

	subroutine getsteps(par,step)
	real(dp), dimension(:), intent(in) :: par 
	!character(len=15), dimension(:), intent(out) :: name ! names of each parameter
	real(dp), dimension(:), intent(out) :: step 
	integer(i4b) :: i,j
    step=0.0_dp
	end subroutine getsteps

	subroutine getdistpop
	integer(i4b) :: i,j
    real(dp) :: temp(nl) 
	distance=0
	popsize=0.0_dp 
!    !***************************************************************************************
!    !*THE BELOW IS OBSOLETE NOW AS IT WAS DEFINED FOR MYDIV'S
!    !*REDO DISTANCE BACK TO WHAT IT WAS IN CENDIV DAYS. AS WELL AS U BENEFITS 
!    !ahu 030717: redid this adjacence thing. see map and appendix of the draft. 
!	distance(1,2)=1 !1 is CT/NY and 2 is WV/OH/PA 
!	
!    distance(2,1)=1 !2 is WV/OH/PA, 1 is tje CT/NY crowd, 7 is the GA/NC crowd, 4 is MO/IL/IN/KY, 5 is NM/CO, 8 is the NV crowd
!	distance(2,4)=1 
!	distance(2,5)=1
!	distance(2,7)=1
!
!	distance(3,2)=1 !3 is Minn, 2 is WV/OH/PA, 4 is MO/IL/IN/KY, 5 is NM/CO..., 8 is the NV crowd
!	distance(3,4)=1 
!	distance(3,5)=1 
!	distance(3,8)=1 
!	
!    distance(4,2)=1 !4 is MO/IL/IN/KY, 2 is WV/OH/PA, 3 is Minn, 5 is NM/CO, 6 is TX/FL, 7 is GA/NC etc.
!	distance(4,3)=1 
!	distance(4,5)=1 
!	distance(4,6)=1 
!	distance(4,7)=1 
!
!	distance(5,8)=1 !5 is NM/CO.., 8 is NV crowd, 3 is Minn, 4 is MO/IL/IN/KY, 6 is TX/FL
!	distance(5,3)=1 
!	distance(5,4)=1
!	distance(5,6)=1
!
!	distance(6,7)=1 !6 is TX/FL, 7 is GA/NC etc. 5 is NM/CO.. , 4 is MO/IL/IN/KY,
!	distance(6,5)=1 
!	distance(6,4)=1 
!
!	distance(7,2)=1 !7 is GA/NC/SC/.../MD/DC, 6 is the TX/FL crowd, 4 is MO/IL/IN/KY, 2 is WV/OH/PA
!	distance(7,4)=1 
!	distance(7,6)=1 
!
!	distance(8,3)=1 !8 is NV crowd, 5 is the NM/CO crowd, 3 is the Minn crowd
!	distance(8,5)=1 !8 is NV crowd, 5 is the NM/CO crowd, 3 is the Minn crowd
!	distance(8,9)=1 !8 is NV crowd, 5 is the NM/CO crowd, 3 is the Minn crowd
!	
!    distance(9,8)=1 !9 CA and 8 is the NV crowd
!    
!    do i=1,nl 
!		do j=1,nl 
!			if (j==i) then 
!				if (distance(j,i)==1) then 		! ahu 061513: for some locaitons, you did have that their distance(i,i) was 1 so corrected this
!                    print*, 'distance(i,i) should not be 1 because that is not adjacence'
!                    stop
!                end if 
!			end if 
!		end do 
!	end do 
!
!    heat(1)=48.
!    heat(2)=51.
!    heat(3)=44.
!    heat(4)=53.
!    heat(5)=55.
!    heat(6)=66.
!    heat(7)=59.
!    heat(8)=47.
!    heat(9)=57.
!
!!
!!    |         Summary of heatdiv2
!!    mydiv |        Mean   Std. dev.       Freq.
!!------------+------------------------------------
!!        1 |   48.288521           0           8
!!        2 |   50.976101           0           4
!!        3 |   44.609676           0           4
!!        4 |    53.72319           0           4
!!        5 |   55.905964           0           7
!!        6 |   66.642204           0           6
!!        7 |   59.903889           0           6
!!        8 |   47.206505           0           8
!!        9 |   57.954315           0           4
!!------------+------------------------------------
!!    Total |   53.796699   6.7434864          51
!	popsize(1)=0.9112_dp	!new england
!	popsize(2)=2.695_dp		!middle atlantic 
!	popsize(3)=2.9598_dp	!east north central
!	popsize(4)=1.2468_dp	!west north central
!	popsize(5)=2.5736_dp	!south atlantic
!	popsize(6)=1.0089_dp	!east south central
!	popsize(7)=1.6121_dp	!west south central
!	popsize(8)=0.7661_dp	!mountain
!	popsize(9)=2.2757_dp	!pacific

!    !1985 Unemployment benefits by location (see state-mydiv.do for more information) 
!    ubenefitbyloc(1)=4983.
!    ubenefitbyloc(2)=5194.
!    ubenefitbyloc(3)=4578.
!    ubenefitbyloc(4)=3452.
!    ubenefitbyloc(5)=4241.
!    ubenefitbyloc(6)=4229.
!    ubenefitbyloc(7)=3864.
!    ubenefitbyloc(8)=4902.
!    ubenefitbyloc(9)=4436.

!        |      Summary of ubenefit_bydiv
!        mydiv |        Mean   Std. dev.       Freq.
!    ------------+------------------------------------
!            1 |   4983.8359           0           6
!            2 |   5194.4204           0           5
!            3 |   4578.1226           0           5
!            4 |   3452.7996           0           4
!            5 |   4241.0591           0           8
!            6 |   4229.2441           0           6
!            7 |   3864.4907           0           6
!            8 |   4902.4458           0           7
!            9 |   4436.0693           0           4
!    ------------+------------------------------------
!        Total |   4453.5139    503.4298          51

!!*************************************************************************************************************************************************

!*************************************************************************************************************************************************
    !DEFINE ADJACENCY FOR CENDIV'S
    distance(1,2)=1 
	distance(2,1)=1 
	distance(2,5)=1 
	distance(2,3)=1 
	distance(3,2)=1 
	distance(3,4)=1 
	distance(3,6)=1 
	distance(3,5)=1 
	distance(4,3)=1 
	distance(4,8)=1 
	distance(4,7)=1 
	distance(4,6)=1 
	distance(5,8)=1 
	distance(5,3)=1 
	distance(5,6)=1 
	distance(6,5)=1 
	distance(6,7)=1 
	distance(6,4)=1 
	distance(6,3)=1 
	distance(7,6)=1 
	distance(7,8)=1 
	distance(7,4)=1 
	distance(8,7)=1 
	distance(8,9)=1 
	distance(8,4)=1 
	distance(9,8)=1 
    do i=1,nl 
		do j=1,nl 
			if (j==i) then 
				if (distance(j,i)==1) then 		! ahu 061513: for some locaitons, you did have that their distance(i,i) was 1 so corrected this
                    print*, 'distance(i,i) should not be 1 because that is not adjacence'
                    stop
                end if 
			end if 
		end do 
	end do 

	popsize(1)=0.9112_dp	!new england
	popsize(2)=2.695_dp		!middle atlantic 
	popsize(3)=2.9598_dp	!east north central
	popsize(4)=1.2468_dp	!west north central
	popsize(5)=2.5736_dp	!south atlantic
	popsize(6)=1.0089_dp	!east south central
	popsize(7)=1.6121_dp	!west south central
	popsize(8)=0.7661_dp	!mountain
	popsize(9)=2.2757_dp	!pacific
	
    !1985 Unemployment benefits by location (see state-mydiv.do for more information) 
    ubenefitbyloc(1)=4983.
    ubenefitbyloc(2)=5164.
    ubenefitbyloc(3)=4126.
    ubenefitbyloc(4)=3997.
    ubenefitbyloc(5)=3958.
    ubenefitbyloc(6)=3221.
    ubenefitbyloc(7)=4914.
    ubenefitbyloc(8)=4212.
    ubenefitbyloc(9)=4586.

!    tab location, su(ubenefit_bydiv)
!    |      Summary of ubenefit_bydiv
!location |        Mean   Std. dev.       Freq.
!------------+------------------------------------
!  1 |   4983.8359           0           6
!  2 |   5164.1675           0           3
!  3 |   4126.6826           0           5
!  4 |   3997.0962           0           7
!  5 |   3958.4944           0           9
!  6 |   3221.4905           0           4
!  7 |   4914.3491           0           4
!  8 |   4212.5503           0           8
!  9 |   4586.2427           0           5
!------------+------------------------------------
!Total |   4290.3927    519.9386          51

    temp=-99999.0_dp
    temp(:)=ubenefitbyloc(:)
    if (policyubenefit==1) then
        ubenefitbyloc=temp*1.2_dp
    else if (policyubenefit==2) then
        ubenefitbyloc=temp*1.4_dp
    else if (policyubenefit==3) then
        ubenefitbyloc(:)=temp(6)
    else if (policyubenefit==4) then
        ubenefitbyloc(:)=temp(2)
    end if
    

	end subroutine getdistpop

	subroutine getones
		ones%co=-99
		ones%sexr=-99
		ones%hme=-99
		ones%endage=-99
		ones%edr=-99
		ones%expr=-99
        ones%kidr=-99
		ones%hhr=-99
        ones%logwr=-999.0_dp
        ones%wr=-9999.0_dp
        ones%wsp=-9999.0_dp
        ones%l=-99
        ones%rel=-99
        ones%rellen=-99
        ones%edsp=-99
        ones%expsp=-99
        ones%kidsp=-99
        ones%hhsp=-99
        ones%logwsp=-999.0_dp
		ones%lsp=-99        
		ones%nomiss=-99
        ones%nn=-99
        ones%mm=-99
        ones%r=-99
        ones%job=-99
        ones%jobswitch=-99 
        ones%emax=-9999.0_dp
        ones%emaxsp=-9999.0_dp
        ones%consr=-9099.0_dp
        ones%consp=-9099.0_dp
        ones%incsum=-9099.0_dp
        ones%consum=-9099.0_dp

        ones_init%id=-99
        ones_init%co=-99
        ones_init%sexr=-99
        ones_init%hme=-99
        ones_init%endage=-99
        ones_init%edr=-99
        ones_init%typ=-99

        ones_draws%qshock=-99
        ones_draws%wmshock=-99
        ones_draws%wfshock=-99
        ones_draws%lshock=-99
        ones_draws%iephub=-99
        ones_draws%iepwfe=-99
        ones_draws%iepsingle=-99
        ones_draws%zshock=-99

        ones_decisions%rel0dec=-99
        ones_decisions%q0dec=-99
        ones_decisions%wm0dec=-99
        ones_decisions%wf0dec=-99
        ones_decisions%l0dec=-99
        ones_decisions%reldec=-99
        ones_decisions%qdec=-99
        ones_decisions%qmdec=-99
        ones_decisions%qfdec=-99
        ones_decisions%wmdec=-99
        ones_decisions%wfdec=-99
        ones_decisions%ldec=-99
        ones_decisions%lmdec=-99
        ones_decisions%lfdec=-99


        !ones_efficiency%ve1=-9999.0_dp
        !ones_efficiency%ve2=-9999.0_dp
        !ones_efficiency%ve3=-9999.0_dp
        !ones_efficiency%ve4=-9999.0_dp
        !ones_efficiency%ve5=-9999.0_dp
        !ones_efficiency%sur=-9999.0_dp
        ones_efficiency%tra1=-9999.0_dp
        ones_efficiency%tra2=-9999.0_dp
        !ones_efficiency%haveenough=.FALSE.
        !ones_efficiency%haveenoughforNB=.FALSE.
        ones_efficiency%maxA=-99
        ones_efficiency%maxB=-99

        ones_mdecisions%qmdec=-99
        ones_mdecisions%wmdec=-99
        ones_mdecisions%lmdec=-99
        ones_fdecisions%qfdec=-99
        ones_fdecisions%wfdec=-99
        ones_fdecisions%lfdec=-99


        end subroutine getones	

        function fnwge(dg,dtyp,dl,dw,de,dr,dia)		!de is educ here but in fnprof it's no longer educ				
            integer(i4b), intent(in) :: dg,dtyp,dl,de,dr,dia						! gender,typ,location,education,experience
            real(dp), intent(in) :: dw								! wage draw
            real(dp) :: fnwge
            if (dg==1) then 
                fnwge=exp(alf1t(dtyp)+alf10(dl)+alf11*one(de==2) + alf12*(dr-1) + dw + alf13*dia + alf14 * (dia**2) ) 
            else if (dg==2) then 
                fnwge=exp(alf2t(dtyp)+alf20(dl)+alf21*one(de==2) + alf22*(dr-1) + dw  + alf23*dia + alf24 * (dia**2) ) 
            end if 
            end function fnwge
        
            function fnprof(dw0,de,dsex) !ahu october2022: note that de was ed before but now it's wtr it's curloc or ofloc (takes on values 5 or 10)
            integer(i4b), intent(in) :: dw0,de,dsex
            real(dp), dimension(3) :: fnprof
            if (onthejobsearch) then 
                fnprof=0.0_dp
                if ( dw0 <= np ) then  
                    if ( de==5 .and. dsex==1 ) then 
                        fnprof(1:2)=exp(psio(1:2)) !exp( psio(1) + psio(2) * abs(dsex==1) + psio(3) * abs(de==2) )	! offer 	 
                    else if ( de==10 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(3)) 
                        fnprof(2)=0.0_dp !psio4 nuissance parameter
                    else if ( de==5 .and. dsex==2 ) then 
                        fnprof(1:2)=exp(psio(5:6)) 
                    else if ( de==10 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(7)) 
                        fnprof(2)=0.0_dp !psio8 nuissance parameter
                    end if 
                    fnprof(3)=exp( 0.0_dp )		! nothing happens												
                else if (dw0 == np1) then 
                    if ( de==5 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(9)) !exp( psio(1) + psio(2) * abs(dsex==1) + psio(3) * abs(de==2) )	! offer 	 
                    else if ( de==10 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(10) )
                    else if ( de==5 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(11)) 
                    else if ( de==10 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(12)) 
                    end if 
                    fnprof(2)=0.0_dp		! 0 since you can't get laid off if you don't have a job! 
                    fnprof(3)=exp(0.0_dp)		! nothing happens												
                else  
                    print*, "in fnprof: dw0 > np1 which doesnt' make sense as that's a state variable " , dw0,de,dsex
                    stop
                end if 
                fnprof(1:3)=fnprof/sum(fnprof)
            else 
                fnprof=0.0_dp
                if ( dw0 <= np ) then  
                    if ( de==5 .and. dsex==1 ) then 
                        fnprof(1:2)=0.0_dp !exp(psio(1:2)) !exp( psio(1) + psio(2) * abs(dsex==1) + psio(3) * abs(de==2) )	! offer 	 
                    else if ( de==10 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(3)) 
                        fnprof(2)=0.0_dp !psio4 nuissance parameter
                    else if ( de==5 .and. dsex==2 ) then 
                        fnprof(1:2)=0.0_dp !exp(psio(5:6)) 
                    else if ( de==10 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(7)) 
                        fnprof(2)=0.0_dp !psio8 nuissance parameter
                    end if 
                    fnprof(3)=exp( 0.0_dp )		! nothing happens												
                else if (dw0 == np1) then 
                    if ( de==5 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(9)) !exp( psio(1) + psio(2) * abs(dsex==1) + psio(3) * abs(de==2) )	! offer 	 
                    else if ( de==10 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(10) )
                    else if ( de==5 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(11)) 
                    else if ( de==10 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(12)) 
                    end if 
                    fnprof(2)=0.0_dp		! 0 since you can't get laid off if you don't have a job! 
                    fnprof(3)=exp(0.0_dp)		! nothing happens												
                else  
                    print*, "in fnprof: dw0 > np1 which doesnt' make sense as that's a state variable " , dw0,de,dsex
                    stop
                end if 
                fnprof(1:3)=fnprof/sum(fnprof)
            end if 
            if (skriv) then 
                if ( abs(  sum(fnprof) - 1.0_dp  ) > eps) then ; print*, "error in getfnprof : offer does not add up " , sum(fnprof) ; stop ; end if 
            end if 
            end function fnprof
        
        
            
            
            function fnprloc(orig,index)
            integer(i4b), intent(in) :: orig,index		! origin location and index to indicate what the home location is
            real(dp), dimension(nl) :: fnprloc
            integer(i4b) :: j,co,typ,homeloc
            real(dp), dimension(nl) :: temp
            fnprloc=0.0_dp
            temp=0.0_dp
            call index2cotyphome(index,co,typ,homeloc)
            !ahu 030717 fnprloc(orig)=logit(psil(1)) !ahu 030717: putting psil(1) inside the exp instead because otherwise very high prloc from origin and low from others and 
                                                     !            we get very low moving rates especially for married people. 
            !do j=1,nl	
                !ahu 030717 if (j /= orig) then 
                !ahu 030717 	fnprloc(j)= exp( psil(2) * distance(j,orig) + psil(3) * popsize(j) ) 
                !ahu 030717 	sum_sans_orig = sum_sans_orig + fnprloc(j) 
                !ahu 030717 end if 
                !if (orig.eq.homeloc) then
                !    temp(j)= exp( psil(1) * one(j==orig) )   ! + psil(2) * distance(j,orig) )    !+ psil(3) * popsize(j) ) 
                !else 
            !        temp(j)= exp( psil(1) * one(j==orig) + psil(2) * one(j==homeloc) )   ! + psil(2) * distance(j,orig) )    !+ psil(3) * popsize(j) ) 
            !end do 
!            if (orig==homeloc) then
!                do j=1,nl	
!                        temp(j)= exp( psil(1) * one(j==orig) )   ! + psil(2) * distance(j,orig) )    !+ psil(3) * popsize(j) ) 
!                end do 
!            else 
!            end if
            do j=1,nl	
                temp(j)= exp(-2.0_dp + psil(1) * one(j==orig) + psil(2) * one(j==homeloc) )   ! + psil(2) * distance(j,orig) )    !+ psil(3) * popsize(j) ) 
            end do 
            do j=1,nl	
                !ahu 030717 if (j /= orig) then 
                !ahu 030717 	fnprloc(j)=(1.0_dp-fnprloc(orig) ) * fnprloc(j)/sum_sans_orig
                !ahu 030717 end if 
                fnprloc(j)=temp(j)/sum(temp(:))
            end do 
            if ( abs(  sum(fnprloc) - 1.0_dp  ) > eps) then ; print*, "error in getfnprloc : offer does not add up " , sum(fnprloc) ; stop ; end if 
            end function fnprloc
            
            function fnprhc(dr,dw)
                integer(i4b), intent(in) :: dr,dw		! experience and employment: w<=np work, w==np1 not work,  w=np2 nothing/can't be a state variable here so if you get this, there's something wrong
                real(dp), dimension(nexp) :: fnprhc
                integer(i4b) :: j
                if (skriv) then 
                    if ( dw > np1 ) then ; print*, "in fnprof: dw0 > np1 which doesnt' make sense as that's a state variable " ; stop ; end if 
                end if 
                !ahu october2022: note that when j=nexp, there will be no j such that j-dr=+1, so fnprhc(nexp) will be 1/0+1+1 = 1 and all other fnprhc(j)'s are 0. 
                fnprhc=0.0_dp
                do j=1,nexp	
                    if ( dw <= np ) then 
                        if ( j-dr == +1 ) then  
                            fnprhc(j)=exp(psihup_e)      !ahu jan19 011719 changing to logit
                        else if (j==dr) then 
                            fnprhc(j)=exp(0.0_dp)  
                        else if ( j-dr == -1 ) then  !ahu jan19 011519 getting rid of probdown
                            fnprhc(j)=0.0_dp !exp(psihdown_e)
                        else 
                            fnprhc(j)=0.0_dp
                        end if 
                    else if ( dw == np1 ) then !ahu october2022 no exp increase or decrease if unemp. so j such that j=dr is 1/0+1+0 =1 and all other fnprhc(j)'s are 0. 
                        if ( j-dr == +1 ) then  
                            fnprhc(j)= 0.0_dp    !ahu jan19 011719 changing to logit
                        else if (j==dr) then 
                            fnprhc(j)=exp(0.0_dp)      !exp(0.0_dp)   !ahu jan19 011719 changing to logit
                        else if ( j-dr == -1 ) then  !ahu jan19 011519 getting rid of probdown
                            fnprhc(j)=0.0_dp !exp(psihdown_u)
                        else 
                            fnprhc(j) = 0.0_dp
                        end if 
                    end if 
                end do 	
                fnprhc(:)=fnprhc(:)/sum(fnprhc)
                !print*, dr,fnprhc(:)
                
                if (skriv) then 
                    if ( abs(sum(fnprhc(:))-1.0_dp) > eps ) then ; print*, " error in fnprhc: prhc does not add up " , dw , sum(fnprhc(:)) ; stop ; end if 
                end if 
                end function fnprhc

 
                function fnmove(empo,kid,educo,trueindex,gender) 
                    integer(i4b), intent(in) :: empo,kid,educo,trueindex,gender
                    real(dp) :: fnmove
                    integer(i4b) :: c,t,h
                    call index2cotyphome(trueindex,c,t,h)			
                    !if (relstat==1) then !are you married 
                    !    fnmove = cstmar  !+ csted*one(educo==2)  !+ cstmar_u * one(empo==np1)  !+ cstkid * one(kid>1) !+ kcst * one(empo==np1)  !kid 1 is no kid, kid 2 is yes kid
                    !else if (relstat==0) then !or are you single
                    !    fnmove = cstsin  !+ csted*one(educo==2)  !+ cstsin_u * one(empo==np1)   !+ cstkid * one(kid>1)  !+ kcst * one(empo==np1)  !kid 1 is no kid, kid 2 is yes kid
                    !else 
                    !    print*, "There is something wrong with relstat in fnmove!"
                    !    stop 
                    !end if 
                    if (educo==1) then !are you married 
                        fnmove = cstned(gender)  !+ csted*one(educo==2)  !+ cstmar_u * one(empo==np1)  !+ cstkid * one(kid>1) !+ kcst * one(empo==np1)  !kid 1 is no kid, kid 2 is yes kid
                    else if (educo==2) then !or are you single
                        fnmove = csted(gender)  !+ csted*one(educo==2)  !+ cstsin_u * one(empo==np1)   !+ cstkid * one(kid>1)  !+ kcst * one(empo==np1)  !kid 1 is no kid, kid 2 is yes kid
                    else 
                        print*, "There is something wrong with relstat in fnmove!"
                        stop 
                    end if 

                    if (empo==np1) then !are you married 
                        fnmove = fnmove+CSTU(gender)  !+ csted*one(educo==2)  !+ cstmar_u * one(empo==np1)  !+ cstkid * one(kid>1) !+ kcst * one(empo==np1)  !kid 1 is no kid, kid 2 is yes kid
                    end if 


                        !ahu october2022: no idea why this was kcst * one(empo==np1) 
                    !fnmove = fnmove / div
                end function fnmove
      
                !function fndiv(educo,sexo)
                !    integer(i4b), intent(in) :: educo,sexo
                !    real(dp) :: fndiv
                !        fndiv= 
                !end function
                !function uhomeia(age)
                !    integer(i4b), intent(in) :: age
                !    real(dp) :: uhomeia
                !    uhomeia= agehome*(age-18)
                !end function
            !function fnprkid(kid0)
            !integer(i4b), intent(in) :: kid0
            !real(dp), dimension(0:maxkid) :: fnprkid
            !integer(i4b) :: j
            !fnprkid=0.0_dp
            !do j=kid0,maxkid
            !	fnprkid(j)=exp(  pkid * (j-kid0) )
            !	if ( abs(j-kid0) > 1 ) then 
            !		fnprkid(j)=0.0_dp	!can only move one step up or down
            !	end if 
            !end do 						
            !fnprkid(0:maxkid)=fnprkid(0:maxkid)/sum(fnprkid(0:maxkid))
            !if (skriv) then 
            !	if ( abs(sum(fnprkid(0:maxkid))-1.0_dp)>eps ) then ; print*, "error in fnprkid: prkid does not add up " , kid0 , sum(fnprkid(0:maxkid)) ; stop ; end if 
            !end if 
            !end function fnprkid


                subroutine q2wloc(dq,dw,dl)
                    ! extract indeces w,l from q 
                    integer(i4b), intent(in) :: dq		
                    integer(i4b), intent(out) :: dw,dl
                    integer(i4b), dimension(2) :: indeces	
                        indeces=lin2ndim( (/ np2 , nl /) , dq )
                        dw=indeces(1)
                        dl=indeces(2)
                        if (skriv) then  
                            if ( dq > nqs ) then ; print*, "q2wl: q > nqs", dq, nqs,indeces ; stop ; end if  
                            if ( dw > np2 ) then ; print*, "q2wl: w > np2" ; stop ; end if  
                            if ( dl > nl  ) then ; print*, "q2wl: l > nl" ; stop ; end if  
                        end if 
                    end subroutine
                    subroutine wloc2q(dq,dw,dl)
                    ! construct combined q from w,l
                    integer(i4b), intent(out) :: dq		
                    integer(i4b), intent(in) :: dw,dl		
                        dq = ndim2lin( (/ np2 , nl /),(/ dw,dl /) )
                        if (skriv) then 		
                            if ( dq > nqs ) then ; print*, "wl2q: q > nqs" ; stop ; end if  
                            if ( dw > np2 ) then ; print*, "wl2q: w > np2" ; stop ; end if  
                            if ( dl > nl  ) then ; print*, "wl2q: l > nl" ; stop ; end if  
                        end if 
                    end subroutine
                
                    subroutine x2edexpkid(dx,de,dr, dkid)
                    ! extract indeces educ,experience from x
                    integer(i4b), intent(in) :: dx		
                    integer(i4b), intent(out) :: de,dr,dkid
                    integer(i4b), dimension(3) :: indeces	
                        indeces=lin2ndim( (/ neduc, nexp, nkid /) , dx )
                        de=indeces(1)
                        dr=indeces(2)
                        dkid=indeces(3)
                    end subroutine
                    subroutine edexpkid2x(dx,de,dr,dkid)
                    !construct combined x from educ,experience
                    integer(i4b), intent(out) :: dx		
                    integer(i4b), intent(in) :: de,dr,dkid
                        dx=ndim2lin( (/ neduc, nexp, nkid /),(/ de,dr,dkid /) )
                    end subroutine
                
                    subroutine index2cotyphome(index,co,typ,home)
                    ! extract indeces cohort,type,educ,homeloc from combined index
                    integer(i4b), intent(in) :: index		
                    integer(i4b), intent(out) :: co,typ,home	
                    integer(i4b), dimension(3) :: indeces	
                    !if (groups) then 
                        !indeces=lin2ndim((/nco,ntyp,nl/),index)
                        !print*, 'this should not be called if groups!'
                        !stop
                        !co=myco !indeces(1)
                        !typ=mytyp !indeces(2)
                        !home=myhome !indeces(3)
                    !else 
                        indeces=lin2ndim((/ncop,ntypp,nhomep/),index)
                        co=indeces(1)
                        typ=indeces(2)
                        home=indeces(3)
                    !end if 
                    end subroutine index2cotyphome
                    subroutine cotyphome2index(index,co,typ,home)
                    !construct combined index from co,typ,home
                    integer(i4b), intent(out) :: index		! combined index
                    integer(i4b), intent(in) :: co,typ,home 
                    !if (groups) then 
                    !	index=1 !ndim2lin((/nco,ntyp,nl/),(/co,typ,home/))
                    !else 
                        index=ndim2lin((/ncop,ntypp,nhomep/),(/co,typ,home/))
                    !end if 
                    end subroutine cotyphome2index
                 
end module params
