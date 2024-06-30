!headloc(ihead)=im; headstr(ihead)='labor market hours by gender/rel/ia';ihead=ihead+1
!ahu 061211: have to control for ia here because the two brs have different ia compositions
!ahu 061211: br 2 has no hours/kids/cohmar simultaneously in the biannual years so if you condition on all that you will just get something until they are ia 28 or something (depending on what the br grouping is)
!ahu 061211: and so if we don't control for ia, it looks as if br 2 females who are cohabiting have decreased their hours of work. but this is just a composition effect.
!ahu 061211: excluding ia 20 because, something looks weird. br 2 works too few hours at ia 20 (for females,coh,nokid). so then when i include them, it looks as if br 2 coh females with no kids work less in the later br. 

module mom
	use params
	use share
	use myaz 
	use sol, only: getdec !ag090122 agsept2022
	implicit none
contains

FUNCTION random(iseed)
! When first call, iseed must be a large positive integer.
! iseed will be changed when exit and be used for next calling.
! The range of the generated random number is between 1 and -1
!
implicit none
integer, intent(inout) :: iseed
real :: random
!
iseed = mod(8121*iseed+28411, 134456) ! 0 =< iseed < 134456
random = real(iseed)/134456. ! 0 < random < 1
!
end FUNCTION random


subroutine random_stduniform(u)
    implicit none
    real,intent(out) :: u
    real :: r
    call random_number(r)
    u = 1 - r
end subroutine random_stduniform
        

subroutine random_stdnormal(x)
    implicit none
    real,intent(out) :: x
    real,parameter :: pi=3.14159265
    real :: u1,u2
    call random_stduniform(u1)
    call random_stduniform(u2)
    x = sqrt(-2*log(u1))*cos(2*pi*u2)
end subroutine random_stdnormal
    

subroutine read_taxes
    integer(i4b) :: pp,ss,mydiv,kk
    real(sp) :: pwages,swages !becuause of lcoate() using pwages and swages as input and locate needs them to be real(sp)
    real(dp) :: statesin,fedsin,statemar,fedmar 
    tax%pwages=-999
    tax%swages=-999
    tax%statesin=-999
    tax%fedsin=-999
    tax%statemar=-999
    tax%fedmar=-999         ! initialize
    if (taxrates==1) then !taxrates 1 means "read marginal tax rates"
    	open(unit=68856,file='taxes1983-marginal-bydiv.txt')
    else if (taxrates==2) then !taxrates 2 means "read average tax rates"
    	open(unit=68856,file='taxes1983-avg-bydiv.txt')
    end if 

    do kk=1,numtaxes !for singles, this is 9 regions times 31 brackets = 279 and for mar it is 9 regions times 31 times 31 brackets which is 
		read(68856,*) mydiv,pp,pwages,ss,swages,statesin,fedsin,statemar,fedmar 
        !Note that I'm using statesin,fedsin,staemar,fedmar vaariables for both marginal adn average here
        !but in stata the marginal rates are defined as srate, frate whereas averages are (user generated) statesin,fedsin,statemar,fedmar etc.
        !don't get confused by this 
        !I'm only doing this so I don't have to redefine variables here but just go with it
		pp=pp+1 !changing the first index to 1 instead of 0 because otherwise the values between 0 and 1000 do not get assigned anything
        ss=ss+1 
        tax(pp,ss,mydiv)%pwages=pwages
		tax(pp,ss,mydiv)%swages=swages
		tax(pp,ss,mydiv)%statesin=statesin
		tax(pp,ss,mydiv)%fedsin=fedsin
		tax(pp,ss,mydiv)%statemar=statemar
		tax(pp,ss,mydiv)%fedmar=fedmar
    end do
    pbracket(1:numbin)=tax(1:numbin,numbin,nl)%pwages
    sbracket(1:numbin)=tax(numbin,1:numbin,nl)%swages
    close(68856)

    !CALCULATE THE PRODUCT BETWEEN THE MARGINAL RATE FOR EACH BRACKET AND THE BRACKET AMOUNT AND THEN SUM THAT UP OVER UP UNTIL THAT BRACKET
    bracketprod_s%state=0. ; bracketprod_s%fed=0. 
    bracketprodsum_s%state=0. ; bracketprodsum_s%fed=0. 
    bracketprod_h%state=0. ; bracketprod_h%fed=0. 
    bracketprodsum_h%state=0. ; bracketprodsum_h%fed=0. 
    bracketprod_w%state=0. ; bracketprod_w%fed=0.
    bracketprodsum_w%state=0. ; bracketprodsum_w%fed=0.
    do mydiv=1,nl
        !SINGLES STATE TAXES
        ss=numbin
        do pp=2,numbin
            bracketprod_s(pp,mydiv)%state = tax(pp,ss,mydiv)%statesin * ( pbracket(pp) - pbracket(pp-1) )
        end do 
        do pp=1,numbin
            bracketprodsum_s(pp,mydiv)%state = sum(bracketprod_s(1:pp,mydiv)%state)
        end do 

        !SINGLES FED TAXES 
        ss=numbin
        do pp=2,numbin
            bracketprod_s(pp,mydiv)%fed = tax(pp,ss,mydiv)%fedsin * ( pbracket(pp) - pbracket(pp-1) )
        end do 
        do pp=1,numbin
            bracketprodsum_s(pp,mydiv)%fed = sum(bracketprod_s(1:pp,mydiv)%fed)
        end do 


        !MARRIED STATE TAXES 
        do pp=2,numbin
            do ss=2,numbin
                bracketprod_h(pp,ss,mydiv)%state  = tax(pp,ss,mydiv)%statemar * ( pbracket(pp) - pbracket(pp-1) )  
                bracketprod_w(pp,ss,mydiv)%state  = tax(pp,ss,mydiv)%statemar * ( sbracket(ss) - sbracket(ss-1) ) 
            end do 
        end do 
        do pp=1,numbin
            do ss=1,numbin
                bracketprodsum_h(pp,ss,mydiv)%state = sum(bracketprod_h(1:pp,ss,mydiv)%state  )
                bracketprodsum_w(pp,ss,mydiv)%state = sum(bracketprod_w(pp,1:ss,mydiv)%state  )
            end do 
        end do 
        
        !MARRIED FED TAXES
        do pp=2,numbin
            do ss=2,numbin
                bracketprod_h(pp,ss,mydiv)%fed  = tax(pp,ss,mydiv)%fedmar * ( pbracket(pp) - pbracket(pp-1) )  
                bracketprod_w(pp,ss,mydiv)%fed  = tax(pp,ss,mydiv)%fedmar * ( sbracket(ss) - sbracket(ss-1) ) 
            end do 
        end do 
        do pp=1,numbin
            do ss=1,numbin
                bracketprodsum_h(pp,ss,mydiv)%fed = sum(bracketprod_h(1:pp,ss,mydiv)%fed  )
                bracketprodsum_w(pp,ss,mydiv)%fed = sum(bracketprod_w(pp,1:ss,mydiv)%fed  )
            end do 
        end do 
        end do !location
        

    if (mysay==0.and.skriv) then 
        mydiv=1
        open(unit=68859,file='checkmarginals_sin.txt')
        ss=numbin
        do pp=1,numbin
            write(68859,'(I4,F9.1,F9.4,3F14.2)') pp,pbracket(pp),tax(pp,ss,mydiv)%statesin,& 
                & bracketprod_s(pp,mydiv)%state,bracketprodsum_s(pp,mydiv)%state,pbracket(pp)*tax(pp,ss,mydiv)%statesin
        end do                                                             
        close(68859)

        open(unit=68859,file='checkmarginals_mar.txt')
        do pp=1,numbin
            do ss=1,numbin
                write(68859,'(2I4,2F9.1,2F9.4,6F14.2)') pp,ss,pbracket(pp),sbracket(ss),tax(pp,ss,mydiv)%statesin,tax(pp,ss,mydiv)%statemar,& 
                                                            & bracketprod_s(pp,mydiv)%state,bracketprodsum_s(pp,mydiv)%state,& 
                                                            & bracketprod_h(pp,ss,mydiv)%state,bracketprodsum_h(pp,ss,mydiv)%state,& 
                                                            & bracketprod_w(pp,ss,mydiv)%state,bracketprodsum_w(pp,ss,mydiv)%state
            end do 
        end do                                                             
        close(68859)
    end if !mysay for printing


end subroutine read_taxes


	! read psid data, save into global structure psiddata, calculate moments
	! notes on data
	!   record total number of observations in nrealdatobs
	!	need ages of observations to be strictly increasing
	!	ids should go from 1 to nrealdat
	!	data entry should look like: idnum, age sex rel kids edm edf incm incf hhm hhf ddm ddf rellen with <0 being blank
	!	what can be missing? anything.either entire observation is missing (rel<0)
	!		or observation is there but no record of housework (if year \=81) (housework hours <0)
	!		if not working, wage will be missing
	! ahu 071712: you can see in familymig.do that everyone starts in the data at age 16 (i drop the other people, but that's not many people)
	! this is because i wanted their location at age 16, because that's my definition of home location
	! keep this in mind, if this ever changes 
	! only want their information after they have completed school (whether high school or college)	!ahu 071712 
	! in the simulation this is automatically taken care of by the fact that the simulation starts at age age0sim(edsim(r))
	subroutine read_actualdata(dat,nper,nperobs) !nper is for num of persons /      nperobs is for num of person-periods
    integer(i4b), intent(in) :: nper,nperobs !in calling program, nper should be set to the sample size of actual data and ndatobs is the total size of all person-periods in actual data
	!type(initcond), dimension(nper), intent(out) :: init !now declared in params, allocated in main and assigned values by calling this from objf
	type(statevar), dimension(mnad:mxa,nper), intent(out) :: dat	!data set. first entry is ia index, second observation number
	integer(i4b) :: kk,id,age,cohort,sexr,rel,kid,edr,edsp,hhr,hhsp,rellen,loc,homeloc,minage,endage,jobswitch,nomiss,ierr,checkminage(nper)
	real(dp) :: wr_perhour,wsp_perhour
	dat=ones 	            ! initialize
    init=ones_init      ! initialize
    checkminage=1
	open(unit=77,file=datafilename)
	do kk=1,nperobs
		read(77,*) id, age, cohort, sexr, rel, kid, edr, edsp, wr_perhour, wsp_perhour, hhr, hhsp, rellen,loc,minage,endage,homeloc,jobswitch
		if (rel/=0.and.rel/=1.and.rel/=-1) then ; print*, "data has other relationship states! ",rel ; stop ; end if !ahu 102112: count cohabiting people as married
		!if (age==mna) hme=loc      !ahu 022517: changing home back to state grew up. 
                                    !because otherwise I have to drop those who I don't observe starting from age 18 and 
                                    !then get a different cohort composition than the original sample and the avg wage per loc's end up being much smaller 
                                    !than the original numbers. 
                                    !see more details on the reason in the explanations for ahu 022517 in familymig_2.do
		dat(age,id)%id=id           !initcond		
        dat(age,id)%co=cohort       !initcond

        dat(age,id)%sexr=sexr       !initcond
		dat(age,id)%hme=homeloc     !initcond
		dat(age,id)%endage=endage   !initcond     
		dat(age,id)%edr=edr         !initcond
        dat(age,id)%minage=minage !initcond 
		dat(age,id)%expr=-99        ! (there is no exp variable in the actual data so this is always -99)
		if (kid==0) then
            dat(age,id)%kidr=1 
        else if (kid==1) then
            dat(age,id)%kidr=2 
        else if (kid==-1) then 
            dat(age,id)%kidr=-99
        else 
            print*, 'something wrong with kid in data'
            stop
        end if 
        !kid  !min(kid,maxkid)      !initcond (it should be just 0 in the beginning but might not be in the actual data so just read it from the data here as initcond)
		!ahu 030217 if (age==mna) init(id)=dat(age,id)%initcond 
        !if (id==6740) then 
        !    print*, "I am 6740"
        !    print*, nper,nperobs
        !    print*, id, age, cohort, sexr, rel, kid, edr, edsp, wr_perhour, wsp_perhour, hhr, hhsp, rellen,loc,minage,endage,homeloc
        !end if
        if (age==minage) then
            checkminage(id)=0
            init(id)=dat(age,id)%initcond
        end if 
        
        !ahu jan 19 010219
        !initial conditions are the state variables of age 17,which are the state variables that agents use inorder to make decisions at the beginning of age 18
        !in order to have moving rates in the data at age 17 (because we do have sim rates at age 17 in simulation since we record initconditions at age 17 and then their decisions become the age 18 variables)
        !if (age==mna) then 
        !    dat(age-1,id)%initcond=init(id)
        !    dat(age-1,id)%rel=0
        !    dat(age-1,id)%l=loc
        !    dat(age-1,id)%hr=0
        !end if 

		dat(age,id)%l=loc
		call get_dathrwge(hhr,wr_perhour,dat(age,id)%hhr,dat(age,id)%wr,dat(age,id)%logwr)        
        dat(age,id)%jobswitch=jobswitch
		dat(age,id)%rel=rel
		if (rel==1) then 
			dat(age,id)%rellen=rellen
		else if (rel==0) then 
			dat(age,id)%rellen=-99 !99
		end if 
		
		if (rel==1) then 
			call get_dathrwge(hhsp,wsp_perhour,dat(age,id)%hhsp,dat(age,id)%wsp,dat(age,id)%logwsp)
		else if (rel==0) then 
			dat(age,id)%hhsp=-99 !99
			dat(age,id)%wsp=-99.0_dp !99.0_dp
            dat(age,id)%logwsp=-99.0_dp
		end if 

		dat(age,id)%edsp=edsp    !for some reason, don't read edsp from the data yet. just set it to -99 and do read it later. ahu 021617: now I read it!  
		dat(age,id)%expsp=-99   !no experience variable in the data 
		dat(age,id)%kidsp=-99    !kidsp is just a simulation concept
				
		call getmiss(dat(age,id),nomiss)
		dat(age,id)%nomiss=nomiss
		
		!some checks:
		if (edr.ne.1.and.edr.ne.2) then
			print*, 'There is something wrong with edr in read_actual data'
			!edr=1
            !dat(age,id)%edr=1
            !stop 
		end if 
		if (dat(age,id)%co<0.or.dat(age,id)%sexr<0) then ; print*, "cohort and sex are negative!", age,dat(age,id)%co,cohort,dat(age,id)%sexr,sexr ; stop ; end if 

        !why is this age<agestart-1 and not age<agestart
        !because for ed, agestart is 22 in sim the initial conditions are assigned to age 21 and then 
        !the move rates at age 21 for ed people in sim turn out to be very large because they move immediately from their starting 
        !home location to the best uloc location. In the data, this moment is not there because of the setting of dat=ones 
        !for ages age<agestart. In order to let the estimation to its thing and be able to compare that large sim moving rate 
        !at age 21 to the data, I am setting the dat to ones for ages 0-20 rather than 0-21. doing this is also more consistent because 
        !now (i.e. when age<agestart case) the simulation has age 21 for the ed people but not the data. 
        !I am not doing this for noed though since for noed people there is no age read from the data that is less than mna
        !check this
        if (age<agestart(edr)-1) then 
			dat(age,id)=ones
		end if
        if (age<mna-1) then !17 is ok now because I now include the age17 observations in psidanalysis_2.do. this happened in feb 2023. age17 is like their period -1. 
            print*, 'There is something wrong with age in read_actual data'
            stop
        end if 
    enddo 
    if (sum(checkminage)>0) then !this would mean that we don't get age=minage for some people which is not right
        print*, 'something is wrong with checkminage',checkminage
        stop
    end if
    if (minval(init(:)%id)<0.or.minval(init(:)%co)<0.or.minval(init(:)%sexr)<0.or.minval(init(:)%hme)<0.or.minval(init(:)%endage)<0.or.minval(init(:)%edr)<0 ) then 
        print*, 'something is wrong with init everything' !because this would mean things did not get set at minage. either minage was not around (which can't be right) or something else. 
        stop
    end if
	close(77)
	end subroutine read_actualdata
	
	subroutine simulate(sim,nper1,nper2) !nperdat is for num of persons in actual data, npersim is for num of persons in sim data (which should be nperdat*nsimeach)
    integer(i4b), intent(in) :: nper1,nper2 !in the calling program, ndat should be set to the # of people for actual data and is read here to get the dimensions of initial conditions
	!type(initcond), dimension(nper1), intent(in) :: init
	type(statevar), dimension(mnad:mxa,nper2), intent(out) :: sim
	type(shock), allocatable, dimension(:,:) :: epsim
	integer(i4b) :: q0,x0,q,x,z,dec(3)  
	integer(i4b) :: relnext,qnext,xnext
	integer(i4b) :: i0,n0,i,n
	integer(i4b) :: rel0,index,trueindex,ia,r,g,endage,qmdec,qfdec
	integer(i4b) :: nn,mm,typsim,hh(2),l(2),nomiss
	integer(i4b) :: p,k,imax				! for random_seed
	real(dp) :: rand,vmax(2),logw(2),inc(2),wage(2),dum(2) !dum is just a dummy variable for yaz_checknb
	logical  :: welldef,newrel0,meet,checkgroups,defj(nc)
	integer(i4b) :: qmatch,xmatch,dd(12),iepjoint,iephub,iepwfe,iepsingle,indeces(2),ww(2),ed(2),expe(2),kid(2)
	real(dp) :: vec(5),mcost(2),surplus,surplusmax,surplusj(nc),val(2),vcheck(2),transfers(2),cons0(2),consnext(2),transo(2)
	integer(i4b) :: l0,j,jmax,qmax,relmax,de(1),job0,jobnext
	!integer, allocatable :: newseed(:)
	!integer ::     seedo
	integer(i4b), dimension(12) :: seed_init !98765
    integer(i4b) :: callfrom !ag090122 agsept2022
    real(dp) :: valso(2)
    real(4), allocatable, dimension(:,:) :: errormeasure1,errormeasure2
    real(4) :: rand_measure
    
     !seedo=98765
	!seedo=seed_init    
	allocate(epsim(mna:mxa,nper2),errormeasure1(mna:mxa,nper2),errormeasure2(mna:mxa,nper2))
	
	!if (iter==1) then
	    call random_seed (size=p)
	    !p=12
	    call random_seed (put = (/(k,k=1,p)/))
       !seed_init=3
       !call random_seed (put = seed_init )
	    if (skriv) then
	        print*, 'Here is p',p
	        print*, 'Here is k',(/(k,k=1,p)/)
	    end if
	!else     
	!    allocate(newseed(p))
	!    newseed=10
	!    call random_seed( put=newseed(1:p)  )
	!    deallocate(newseed)
	!end if 

    !call random_seed()
    !call random_seed(size=seed_size) 
    !print*, "seed size for random seed is:", seed_size
    !allocate(seedonzo(seed_size))
    !call random_seed(get=seedonzo)
    !print*, "system generated random seed is:"
    !print*, seedonzo
    !seedonzo=300000
    !call random_seed(put=seedonzo)
    !print*, "my generated random seed is:"
    !print*, seedonzo
    !call random_number(rhino)
    !print*, "rhino:", rhino 
    !deallocate(seedonzo)
    
    !idumo=-1
    !rhino=ran(idumo) 
    !print*,'here is rhino today', iam,rhino
    !rhino=ran(idumo) 
    !print*,'here is rhino tomorrow', iam,rhino
   
	draws=ones_draws
    decisions=ones_decisions
    myefo=ones_efficiency
	sim=ones 	! initialize to smallest valid value
	r=0
	do nn=1,nper1			! person being simulated this is the same number of people (ndat) as actual data !ahu 031911 big change	
		do mm=1,nsimeach	! number of simulations for this person
			r=r+1
			do ia=mna,mxa
				call random_number(rand)
				epsim(ia,r)%q=rand !random(seedo)
                call random_number(rand)
				epsim(ia,r)%x=rand !random(seedo) 
				call random_number(rand)
				epsim(ia,r)%marie=rand 
				call random_number(rand)
				epsim(ia,r)%meet=rand !random(seedo) 
				call random_number(rand)
				epsim(ia,r)%meetq=rand !random(seedo)
				call random_number(rand)
				epsim(ia,r)%meetx=rand !random(seedo) 
				call random_number(rand)
				epsim(ia,r)%move=rand !random(seedo) 
				!call minus1(sim(ia,r))		!initialize to smallest valid value
                call random_number(rand)
                epsim(ia,r)%typ=rand !random(seedo)
			end do 
		end do 
	end do 
    if (nper2/=nper1*nsimeach) then 
    print*, "problem with npersim"
    stop 
    end if 
	!print*, 'epsim',epsim(18,5)%q,epsim(18,5)%x
	!call random_number(rand) 
    !print*, 'rand',rand
	!if (skriv) print*, 'epsim',epsim(18,5)%q,epsim(18,5)%x
	!if (skriv) call random_number(rand) 
    !if (skriv) print*, 'rand',rand
    !if (skriv) print*, 'random(seedo)',random(seedo),random(seedo)
	!if (skriv) print*, 'random(seedo)',random(seedo),random(seedo)

    !if (skriv) print*, 'random(seedo)',random(seedo),random(seedo)
	!if (skriv) print*, 'random(seedo)',random(seedo),random(seedo)
	
    !measurement error
    if (mysay==0.and.skriv.and.iter==1) then 
        open(unit=1000, file='chkrando.txt',status='replace')
    end if
    r=0 ; errormeasure1=-999.0_dp ; errormeasure2=-999.0_dp
    do nn=1,nper1			! person being simulated this is the same number of people (ndat) as actual data !ahu 031911 big change	
        do mm=1,nsimeach	! number of simulations for this person
            r=r+1
            do ia=mna,mxa
                call  random_stdnormal(rand_measure)
                errormeasure1(ia,r)=rand_measure * sig_errormeasure + mu_errormeasure
                call  random_stdnormal(rand_measure)
                errormeasure2(ia,r)=rand_measure * sig_errormeasure + mu_errormeasure
                if (mysay==0.and.skriv) write(1000,*) errormeasure1(ia,r),errormeasure2(ia,r)
            end do 
        end do 
    end do 
    if (mysay==0.and.skriv) close(1000)    
    if (nper2/=nper1*nsimeach) then 
        print*, "problem with npersim"
        stop 
    end if 

    checkgroups=.true.
	r=0								! overall count of total number of simulations
    typsim=-1
	iddat: do nn=1,nper1    			! person being simulated. this is the same # of people as in actual data !ahu 031911 big change	
		idsim: do mm=1,nsimeach		! number of simulations for this person
			r=r+1
					!sim(mna,r)%initcond=init(nn)   !co, sex, hme, endage   !ahu0115del

			!if (skriv.and.mod(r,1000)==0.0_dp) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !mod in order to avoid huge output
            !if (  skriv  .and.  (  (nn==348.and.mm==6).or.(nn==482.and.mm==7).or.(nn==482.and.mm==8).or.(nn==504.and.mm==6).or.(nn==559.and.mm==8)     )    ) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !mod in order to avoid huge output
            !if (  skriv  .and.  (  (nn==1.and.mm==2).or.(nn==1.and.mm==10)    )    ) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !mod in order to avoid huge output
            if (  skriv  .and.  (  (nn==4.and.mm==1).or.(nn==6.and.mm==8).or.(nn==7.and.mm==7)    )    ) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !mod in order to avoid huge output

        
		    IF (init(nn)%edr==1) THEN ! low ed, pick type based on ptypehs(co)
			    typsim=multinom(ptypehs(:,init(nn)%sexr),epsim(MNA,r)%typ) 
		    ELSE IF (init(nn)%edr==2) THEN  ! high ed type, pick type based on ptypecol(co)  
                typsim=multinom(ptypecol(:,init(nn)%sexr),epsim(MNA,r)%typ) 
            ELSE 
                print*, 'something is wrong with type',init(nn)%edr,typsim,ptypehs,ptypecol,epsim(MNA,r)%typ
                stop
		    ENDIF
            if (ntypp==1.and.typsim.ne.1) then
                print*, 'something wrong with ntypp and type',ntypp,typsim
                stop
            end if 
            !typ=1
            	                    
            if (init(nn)%edr.ne.1.and.init(nn)%edr.ne.2) then ; print*, 'There is something wrong with init edr' ; stop ; end if  
			call cotyphome2index(trueindex,init(nn)%co,typsim,init(nn)%hme)
            !index is set to 1 if groups, because each processor has its own value function and so that they don't need to be declared with dimension ninp
            if (groups) then 
                index=1
            else 
                index=trueindex
            end if
			!if (skriv.and.trueindex==1) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !ahu 0327 trueindex 2
			!if (skriv) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !ahu 040918 del and remove later. the above yaz was the original yaz
            if (typsim<0) then 
                print*, 'typ did not get assigned'
                stop
            end if 
			ind: if ( (.not.groups) .or.  (myindex==trueindex) ) then !(myco==init(nn)%co.and.mytyp==typsim.and. myhome==init(nn)%hme)  ) then	  
                !if (checkgroups) then 
                !    write(*,'("here I am",7I4)') mysay,init(nn)%co,typsim,init(nn)%hme,myco,mytyp,myhome
                !    checkgroups=.false.
                !end if
                !write(*,'("Here I am", 6I4)') index,trueindex,myindex,init(nn)%co,typsim,init(nn)%hme

				if (init(nn)%edr.ne.1.and.init(nn)%edr.ne.2) then ; print*, 'There is something wrong with init edr' ; stop ; end if  
				rel0=0 
				newrel0=.false.
                q0 = wl2q(np1,init(nn)%hme)  ; job0=np1 !; jobswitch=-999 !everyone starts out unemployed
				x0 = erkid2x( init(nn)%edr, 1, 1) !ed,exp which is 1 in the beginning and kidr which is 1 in the beginning (for kid, 1 means no kid) 
				g=init(nn)%sexr
				endage=min(init(nn)%endage,mxa)				
                age: do ia=agestart(init(nn)%edr),endage  !mna,endage
        			!if (skriv.and.trueindex==1.and.ia==18) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !ahu 0327 trueindex 2 !ahu 040918 del and remove later
                    !if (skriv.and.ia<=20) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !ahu 0327 trueindex 2 !ahu 040918 del and remove later
					call qx2hrwge(g,rel0,q0,x0,trueindex,hh,l,wage,logw,ia-1)
					if (rel0==0) call x2edexpkid(x0,ed(g),expe(g),kid(g))
					if (rel0==1) then
						ed(:)=xx2e(:,x0)
						expe(:)=xx2r(:,x0)
						kid(:)=xx2kid(:,x0)
					end if
                    if (ed(g).ne.init(nn)%edr) then ; print*, 'something wrong with ed and init%ed' ; stop ; end if
					if (q0<=0) print*, "q0 is no good ",r,ia-1,init(nn)%hme,wl2q(np1,init(nn)%hme)                    
					!if (yaz) then ; write(400,'(4/,2I8)') r,ia ; write(400,'("State Variables:")') ; call yaz_sim(g,rel0,q0,x0) ; end if 
                    if (yaz) then ; write(400,*) epsim(ia,r)%q,epsim(ia,r)%x,logw(g) ; end if !ahu 012019 del
					!!ag 110416: changed to have sim(ia-1,r) at the beginning of the sim loop rather than sim(ia,r) in order to not have 0 emp at age 18
					init(nn)%typ=typsim ; sim(ia-1,r)%initcond=init(nn)   !co, sex, ed, type, hme, endage
                    !sim(ia-1,r)%typ=typsim   !co, sex, hme, endage
                    if (sim(ia-1,r)%edr.ne.ed(g)) then ; print*, 'something wrong with ed' ; stop ; end if
					sim(ia-1,r)%expr=expe(g)
					sim(ia-1,r)%kidr=kid(g)                    
					sim(ia-1,r)%hhr = hh(g)
                    sim(ia-1,r)%wr = wage(g)                          !VERY IMPORTANT: WAGE DURING STAY PERIODS DOES FALL EVEN WITHOUT MEASUREMENT ERROR BECAUSE OF THE AGESQ TERM (WHICH MOSTLY APPLIES TO MEN BECAUSE THEIR AGESQ IS HIGHER) 
                    sim(ia-1,r)%logwr = logw(g) + errormeasure1(ia,r) !VERY IMPORTANT TO NOTE THAT ERRORMEASURE IS ONLY BEING ADDED TO LOG WAGES WHICH ARE THE ONLY WAGES USED IN COMPARING TO THE DATA
                                                                      !THE OTHER WAGE VARIABLE, WHICH IS WR, DOES NOT INCLUDE ERRORMEASURE. IT IS THE TRUE SIMULATED WAGE AND IS ONLY USED IN UNDERSTANDING MODEL PURPOSES. 
                                                                      !WR COMPARISON BETWEEN SIMULATED DATA AND ACTUAL DATA DO NOT CONTRIBUTE TO THE OBJECTIVE CRITERION 
                    sim(ia-1,r)%l = l(g)
					sim(ia-1,r)%rel=rel0
                    sim(ia-1,r)%job=job0 !to get job switch
                    !sim(ia-1,r)%jobswitch=jobswitch                                
                    if (rel0>0) then                      ! r is male/sp is female or r is female/sp is male
						i = one(g==1) * 2 + one(g==2) * 1 ! r is male/sp is female or r is female/sp is male
						if (newrel0 .or. ia==agestart(init(nn)%edr) ) then	!agestart(simdat(mna,nn)%edr)   ) then		!if (newrelsimsimdat(ia,r) .or.(ia==minage)) then correct this in cohabitation also. cohabitation correction. 
							sim(ia-1,r)%rellen=1
						else
							sim(ia-1,r)%rellen=sim(ia-2,r)%rellen+1
						endif                        
						sim(ia-1,r)%edsp=ed(i)
						sim(ia-1,r)%expsp=expe(i)
						sim(ia-1,r)%kidsp=kid(i)                                            
                        sim(ia-1,r)%hhsp   = hh(i)
                        sim(ia-1,r)%wsp = wage(i)
                        sim(ia-1,r)%logwsp = logw(i)  + errormeasure2(ia,r)
                        sim(ia-1,r)%lsp = l(i)  !just for checking purposes
                        if (g==1) then ; sim(ia-1,r)%emax=emaxm_c(x0,q0,ia) ; sim(ia-1,r)%consr=cons0(1) ; end if
                        if (g==2) then ; sim(ia-1,r)%emax=emaxf_c(x0,q0,ia) ; sim(ia-1,r)%consr=cons0(2) ; end if
                        if (g==1) then ; sim(ia-1,r)%emaxsp=emaxf_c(x0,q0,ia) ; sim(ia-1,r)%consp=cons0(2) ; end if
                        if (g==2) then ; sim(ia-1,r)%emaxsp=emaxm_c(x0,q0,ia) ; sim(ia-1,r)%consp=cons0(1) ; end if
                        if (l(g).ne.l(i)) then ; print*, 'lg and li not equal' ; stop ; end if 
                        if ( sim(ia-1,r)%hhr==1 .and. sim(ia-1,r)%hhsp==1) then 
                            sim(ia-1,r)%incsum = sim(ia-1,r)%wr + sim(ia-1,r)%wsp  + ubenefit_c(q0)
                            if (ubenefit_c(q0)>0) then ; print*, "Ubenefit should be 0 but it is not",ubenefit_c(q0) ; stop ; end if 
                        else if ( sim(ia-1,r)%hhr==1 .and. sim(ia-1,r)%hhsp==0) then 
                            sim(ia-1,r)%incsum = sim(ia-1,r)%wr + ubenefit_c(q0)
                        else if ( sim(ia-1,r)%hhr==0 .and. sim(ia-1,r)%hhsp==1) then 
                            sim(ia-1,r)%incsum = sim(ia-1,r)%wsp  + ubenefit_c(q0)
                        else if ( sim(ia-1,r)%hhr==0 .and. sim(ia-1,r)%hhsp==0) then 
                            sim(ia-1,r)%incsum =  ubenefit_c(q0) 
                        end if 
                        sim(ia-1,r)%consum=sim(ia-1,r)%consr + sim(ia-1,r)%consp
					else if (rel0==0) then 
						sim(ia-1,r)%rellen=-99		
						sim(ia-1,r)%edsp=-99
						sim(ia-1,r)%expsp=-99
						sim(ia-1,r)%kidsp=-99                                            
                        sim(ia-1,r)%hhsp   = -99	
                        sim(ia-1,r)%wsp = -9999.0_dp
						sim(ia-1,r)%logwsp = -9999.0_dp
						sim(ia-1,r)%emaxsp = -9999.0_dp ; sim(ia-1,r)%consp=-9999.0_dp
						sim(ia-1,r)%lsp = -99  !just for checking purposes
                        if (g==1) sim(ia-1,r)%emax=emaxm_s(x0,q0,ia)
                        if (g==2) sim(ia-1,r)%emax=emaxf_s(x0,q0,ia)
                        sim(ia-1,r)%incsum=-9999.0_dp
                        sim(ia-1,r)%consum=-9999.0_dp
                        if (sim(ia-1,r)%hhr==1) then ; sim(ia-1,r)%consr=wage(g) ; else if (sim(ia-1,r)%hhr==0) then  ; sim(ia-1,r)%consr=ubenefit_s(q0) ; else ; sim(ia-1,r)%consr=-9999.0_dp ; end if
					endif
					sim(ia-1,r)%nn = nn
					sim(ia-1,r)%mm = mm
					sim(ia-1,r)%r = r
					call getmiss(sim(ia-1,r),nomiss)
					sim(ia-1,r)%nomiss=nomiss
					if (yaz) call yaz_simpath(ia-1,nn,mm,r,sim(ia-1,r))
					!if (ia==47 .and. sim(ia-1,r)%hhr ==1) then 
                    !    print*, 'I found it!',ia-1,sim(ia-1,r)%co,sim(ia-1,r)%sexr,sim(ia-1,r)%rel
                    !    stop
                    !end if 
					if (rel0==0) then 
                        decisions(ia-1,r)%rel0dec=0
                        decisions(ia-1,r)%q0dec=q0 
                        if (g==1) decisions(ia-1,r)%wm0dec=q2w(q0)  
                        if (g==2) decisions(ia-1,r)%wf0dec=q2w(q0) 
                        decisions(ia-1,r)%l0dec=q2l(q0)                                 
						q = multinom( ppsq(:,q0,trueindex,g) , epsim(ia,r)%q )     
						x = multinom( ppsx(:,q0,x0)   , epsim(ia,r)%x ) 
						meet=( epsim(ia,r)%meet<=pmeet(ia) )
						z	= multinom( ppmarie(:)			, epsim(ia,r)%marie) !MARSHOCK Z HAS NOT BEEN AROUND FOR A WHILE (NZ HAS BEEN JUST 1)
						iepsingle = multinom( ppmovesingle(:)   , epsim(ia,r)%move )  
                        draws(ia-1,r)%qshock=q ; draws(ia-1,r)%lshock=q2l(q) ; draws(ia-1,r)%iepsingle=iepsingle
                        if (g==1) draws(ia-1,r)%wmshock=q2w(q)
                        if (g==2) draws(ia-1,r)%wfshock=q2w(q)
                        if (g==1) dec(1) = decm_premarmkt(iepsingle,x,q,q0,ia,index) 
						if (g==2) dec(1) = decf_premarmkt(iepsingle,x,q,q0,ia,index) 
						dec(2)=x 
						relnext=0 ; qnext=dec(1) ; xnext=dec(2) ; jobnext=q2w(qnext)   ! next period's variabeles are these unless marriage market happens:
						qmatch	= multinom( ppmeetq(:, dec(1) )	, epsim(ia,r)%meetq) 
						!xmatch	= multinom( ppmeetx(:, dec(2) )	, epsim(ia,r)%meetx) 
                        xmatch	= multinom( probmeetx(:, dec(2),TYPSIM,g )	, epsim(ia,r)%meetx) 
						if (yaz) then 
                            write(400,'("person number ",3I4)') nn,mm,iter
                            write(400,'("age           ",I4)') ia
                            write(400,'("in simulate, rel0=0,r,nn,mm,ia,iter: ",3I8,2I4)') r,nn,mm,ia,iter
                            write(400,'("Trueindex:",I4)') trueindex
						    write(400,'("Draws epsimq,epsimx:",2F10.2)')  epsim(ia,r)%q, epsim(ia,r)%x  !, epsim(ia,r)%marie
						    write(400,'("Draws q,x:",2I10)')  q,x
						    write(400,'("Draws:")') ; call yaz_sim(g,rel0,q,x)
						    write(400,'("Single Dec Before Mar Mkt:")') ; call yaz_sim(g,rel0,qnext,xnext)
						    write(400,'("Match:")') ;  call yaz_simmatch(meet,qmatch,xmatch)
                        end if                         
						if (meet) then 	
							if (g==1) then ; q = q2qq(dec(1),qmatch) ; x = x2xx(dec(2),xmatch) ; end if 
							if (g==2) then ; q = q2qq(qmatch,dec(1)) ; x = x2xx(xmatch,dec(2)) ; end if 
							relnext=dec_mar(z,x,q,ia,index)
							if (relnext==1) then 
								qnext=q ; xnext=x ; jobnext=qq2w(g,qnext)
							else if (relnext==0) then
								qnext=dec(1) ; xnext=dec(2) ; jobnext=q2w(qnext)
							end if 
						    if (yaz) then 
                                write(400,'("here is decmar:",5I8)') q,x,ia,index,dec_mar(z,x,q,ia,index)
							    write(400,'("decision at mar mkt:")') ;  call yaz_simdecmar(relnext)
						    end if      !!!                   
                        end if !meet      
                            decisions(ia-1,r)%reldec=relnext
                            if (relnext==1) then 
                                qmdec=qq2q(1,qnext)       ; qfdec=qq2q(2,qnext)
                                decisions(ia-1,r)%qmdec=qmdec       ; decisions(ia-1,r)%qfdec=qfdec 
                                decisions(ia-1,r)%wmdec=qq2w(1,qnext)  ; decisions(ia-1,r)%wfdec=qq2w(2,qnext) 
                                decisions(ia-1,r)%ldec=qq2l(1,qnext)         
                            else if (relnext==0) then 
                                if (g==1) then 
                                    qmdec=qnext     ; qfdec=-88
                                    decisions(ia-1,r)%qmdec=qmdec       ; decisions(ia-1,r)%qfdec=qfdec 
                                    decisions(ia-1,r)%wmdec=q2w(qmdec)  ; decisions(ia-1,r)%wfdec=-88
                                    decisions(ia-1,r)%lmdec=q2l(qmdec)  ; decisions(ia-1,r)%lfdec=-88        
                                else if (g==2) then 
                                    qmdec=-88     ; qfdec=qnext
                                    decisions(ia-1,r)%qmdec=qmdec       ; decisions(ia-1,r)%qfdec=qfdec 
                                    decisions(ia-1,r)%wmdec=-88        ; decisions(ia-1,r)%wfdec=q2w(qfdec)
                                    decisions(ia-1,r)%lmdec=-88        ; decisions(ia-1,r)%lfdec=q2l(qfdec)  
                                end if 
                            end if 
                        else if (rel0==1) then 
                        decisions(ia-1,r)%rel0dec=1 
                        decisions(ia-1,r)%q0dec=q0 
                        decisions(ia-1,r)%wm0dec=qq2w(1,q0)  ; decisions(ia-1,r)%wf0dec=qq2w(2,q0) 
                        decisions(ia-1,r)%l0dec=qq2l(1,q0)                                 
						q	= multinom( ppcq(:,q0,trueindex)	, epsim(ia,r)%q) 
						x	= multinom( ppcx(:,q0,x0)	, epsim(ia,r)%x) 	
						z	= multinom( ppmarie(:)		, epsim(ia,r)%marie) 
						iepjoint = multinom( ppmovejoint(:)   , epsim(ia,r)%move ) 
                        indeces=lin2ndim( (/ nepsmove , nepsmove /) , iepjoint )
                        iephub=indeces(1) ; iepwfe=indeces(2) 
                        draws(ia-1,r)%qshock=q ; draws(ia-1,r)%wmshock=qq2w(1,q) ; draws(ia-1,r)%wfshock=qq2w(2,q) ; draws(ia-1,r)%lshock=qq2l(1,q)  ; draws(ia-1,r)%iephub=iephub ; draws(ia-1,r)%iepwfe=iepwfe ; draws(ia-1,r)%zshock=z 
                        if (yaz) then 
                            write(400,*) 
                            write(400,*) 
                            write(400,'("person number ",3I4)') nn,mm,iter
                            write(400,'("age           ",I4)') ia
                            write(400,'("in simulate, rel0=1,r,nn,mm,ia,iter: ",3I8,2I4)') r,nn,mm,ia,iter
                            write(400,'("Trueindex:",I4)') trueindex
						    write(400,'("Draws epsimq,epsimx,empsimmarie,epsimmove:",4F10.2)')  epsim(ia,r)%q, epsim(ia,r)%x, epsim(ia,r)%marie, epsim(ia,r)%move
						    write(400,'("Draws q,x,z,iepjoint:",4I10)')  q,x,z,iepjoint
						    write(400,'("Draws:")') ; call yaz_sim(g,rel0,q,x)
						end if                         
                        valso=pen !ag090122 agsept2022
						!callfrom is for telling getdec and therefore yaz_getdec where we are when yaz_getdec is called from within getdec
                        !(if skriv and yaz true) when callfrom=80 then within getdec, I call yaz_getdec to write in 400 the altspecific value functions and bargaining stuff
                        !in getdec when callfrom is 40 or 80, then dd(8) altj and dd(9) altq are filled in within the choice loop 
                        callfrom=80 !ag090122 agsept2022
                        dd = (/ia,trueindex,q,x,z,q0,callfrom,-1,-1,-1,iepjoint,g /) 	! (/ ia,index,q,x,z,q0,callfrom,altj,altq,relmax,iepsmove,gender /)  	  
                        myefficiency=ones_efficiency
                        call getdec(dd,vmax,valso,transo)  
                        myefo(ia-1,r)=myefficiency                        
                        relnext=dd(10)
						if (yaz) then ; call yaz_decision(dd,vmax) ; end if	    !callrom tells yaz_decision where we are    
						if (relnext==1) then 
							qnext=dd(9) ; xnext=dd(4) ;  jobnext=qq2w(g,qnext) ; consnext(1:2)=transo(1:2)
                            qmdec=qq2q(1,qnext)       ; qfdec=qq2q(1,qnext)
                            decisions(ia-1,r)%qmdec=qmdec       ; decisions(ia-1,r)%qfdec=qfdec 
                            decisions(ia-1,r)%wmdec=q2w(qmdec)  ; decisions(ia-1,r)%wfdec=q2w(qfdec)
                            decisions(ia-1,r)%lmdec=q2l(qmdec)  ; decisions(ia-1,r)%lfdec=q2l(qfdec)        
                            decisions(ia-1,r)%reldec=1    
                            decisions(ia-1,r)%qdec=qnext 
                            decisions(ia-1,r)%wmdec=qq2w(1,qnext)  ; decisions(ia-1,r)%wfdec=qq2w(2,qnext) 
                            decisions(ia-1,r)%ldec=qq2l(1,qnext)     

                            !To see whether they are tied stayers or tied movers or neither 
                            i = qq2q(1,q) ; n=xx2x(1,x) ; i0 = qq2q(1,q0)	! translate couple q and x and q0 into singles
                            qmdec=decm_postdiv(iephub,n,i,i0,ia,index) 
                            i = qq2q(2,q) ; n=xx2x(2,x) ; i0 = qq2q(2,q0)	! translate couple q and x and q0 into singles
                            qfdec=decf_postdiv(iepwfe,n,i,i0,ia,index) 
                            mdecisions(ia-1,r)%qmdec=qmdec       ; fdecisions(ia-1,r)%qfdec=qfdec 
                            mdecisions(ia-1,r)%wmdec=q2w(qmdec)  ; fdecisions(ia-1,r)%wfdec=q2w(qfdec)
                            mdecisions(ia-1,r)%lmdec=q2l(qmdec)  ; fdecisions(ia-1,r)%lfdec=q2l(qfdec)        

						else if (relnext==0) then 
							i = qq2q(g,q) ; n=xx2x(g,x) ; i0 = qq2q(g,q0)	! translate couple q and x and q0 into singles
							if (g==1) dec(1) = decm_postdiv(iephub,n,i,i0,ia,index) 
							if (g==2) dec(1) = decf_postdiv(iepwfe,n,i,i0,ia,index) 
							dec(2)=n
							qnext=dec(1) ; xnext=dec(2)    ;   jobnext=q2w(qnext)   ; consnext(1:2)=transo(1:2)   
                            !the below is just for writing decisions
                            i = qq2q(1,q) ; n=xx2x(1,x) ; i0 = qq2q(1,q0)	! translate couple q and x and q0 into singles
                            qmdec=decm_postdiv(iephub,n,i,i0,ia,index) 
                            i = qq2q(2,q) ; n=xx2x(2,x) ; i0 = qq2q(2,q0)	! translate couple q and x and q0 into singles
                            qfdec=decf_postdiv(iepwfe,n,i,i0,ia,index) 
                            decisions(ia-1,r)%qmdec=qmdec       ; decisions(ia-1,r)%qfdec=qfdec 
                            decisions(ia-1,r)%wmdec=q2w(qmdec)  ; decisions(ia-1,r)%wfdec=q2w(qfdec)
                            decisions(ia-1,r)%lmdec=q2l(qmdec)  ; decisions(ia-1,r)%lfdec=q2l(qfdec)        
                            decisions(ia-1,r)%reldec=0    
                            !the above is just for writing decisions
                            if (consnext(1)>pen .or. consnext(2)>pen )  then ; print*, "consnext is not pen even though they hate each other now" ; end if 
						end if
					end if ! rel
                    !sim(ia-1,r)%draw=draw(ia-1,r)  ; sim(ia-1,r)%decisions=decisions(ia-1,r)
					if (yaz) then 
						write(400,'("Next: ")') ; call yaz_sim(g,relnext,qnext,xnext) 
					end if 
                    !jobswitch=-999
                    if ( job0>=1 .and. job0<=np .and. jobnext>=1 .and. jobnext<=np ) then 
                        if (job0==jobnext) then  
                            sim(ia-1,r)%jobswitch=0
                        else if (job0/=jobnext) then 
                            sim(ia-1,r)%jobswitch=1 
                        end if
                    end if 
					newrel0 = (rel0==0.and.relnext==1) 
					rel0 = relnext 
					q0 = qnext
					x0 = xnext
                    job0=jobnext ; cons0=consnext
					relnext=-99 ; qnext=-99 ; xnext=-99 ; q=-99 ; x=-99 ; z=-99 ; iepjoint=-99 ; iepsingle=-99 ; iephub=-99 ; iepwfe=-99 ; qmatch=-99 ; xmatch=-99 ; dec=-99 ; meet=.FALSE.  ; jobnext=-999 ; consnext=pen
					if (init(nn)%edr==2.and.ia<agestart(init(nn)%edr)) sim(ia,r)=ones
                    if (  (ia<init(nn)%minage.or.ia>init(nn)%endage) .and.NOMISSINGOBS) sim(ia,r)=ones
				end do	age
			end if ind
		end do idsim
	end do iddat
	deallocate(epsim)		
                    !if (writeallsim) then 
                    !    CALL MPI_SEND(sim(ia-1,r)%sex    ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%typ    ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%hme    ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%l      ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%job    ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%rel    ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%edr    ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%edsp   ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%hhr    ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%hhsp   ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%kidr   ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%kidsp  ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%expr   ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%expsp  ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%wr     ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%wsp    ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%consr  ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%consp  ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%incsum ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%emax   ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_SEND(sim(ia-1,r)%emaxsp ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
!
                    !    CALL MPI_RECV(simall_int(1,ia-1,r)  ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(2,ia-1,r)  ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(3,ia-1,r)  ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(4,ia-1,r)  ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(5,ia-1,r)  ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(6,ia-1,r)  ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(7,ia-1,r)  ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(8,ia-1,r)  ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(9,ia-1,r)  ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(10,ia-1,r) ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(11,ia-1,r) ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(12,ia-1,r) ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(13,ia-1,r) ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_int(14,ia-1,r) ,1, MPI_INT,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_rea(15,ia-1,r) ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_rea(16,ia-1,r) ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_rea(17,ia-1,r) ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_rea(18,ia-1,r) ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_rea(19,ia-1,r) ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_rea(20,ia-1,r) ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !    CALL MPI_RECV(simall_rea(21,ia-1,r) ,1, MPI_REAL8,0,mysay, MPI_COMM_WORLD,mpistat,mpierr)
                    !end i
	end subroutine simulate

	subroutine qx2hrwge(g,rel,q,x,trueindex,hh,l,wage,logw,age)
	integer(i4b), intent(in) :: g,rel,q,x,trueindex,age
	integer(i4b), dimension(2) :: ed,expe
	integer(i4b), intent(out) :: hh(2),l(2)
	real(dp), intent(out) :: wage(2),logw(2) 	
	integer(i4b) :: i,w(2)
	hh=-99
    wage=-99.0_dp
	logw=-99.0_dp
	l=-99
    if (rel == 0) then			
		w(g)=q2w(q) 
		l(g)=q2l(q) 
		ed(g)=x2e(x) 
		expe(g)=x2r(x) 
		if ( w(g)<=np ) then 
			hh(g)=1
            wage(g) = ws(g,x,q,age,trueindex)   !fnwge( g,typ,l(g) ,wg(w(g),g) ,ed(g) ,expe(g))
			logw(g) = log( wage(g) )        !log ( fnwge( g,typ,l(g) ,wg(w(g),g) ,ed(g) ,expe(g))  )			
		else
			hh(g)=0
			wage(g) = -99.0_dp
            logw(g) = -99.0_dp
		endif
		if ( w(g)<=0.or.w(g)>np1 )  then ; print*, "error in read_dat: w<=0 or w>np1 !!" ; stop ; end if 
	else if (rel > 0) then 		
		w(:)=qq2w(:,q)
		l(:)=qq2l(:,q)
		ed(:)=xx2e(:,x)
		expe(:)=xx2r(:,x)
		do i=1,2
			if ( w(i)<=np ) then 
				hh(i)=1
				wage(i) = wc(i,x,q,age,trueindex)   !fnwge(i,typ,l(i) ,wg(w(i),i) ,ed(i) ,expe(i)) 		
                logw(i) = log( wage(i) )        !log ( fnwge(i,typ,l(i) ,wg(w(i),i) ,ed(i) ,expe(i)) )			
			else
				hh(i)=0
                wage(i) = -99.0_dp
				logw(i) = -99.0_dp
			endif
			if ( w(i)<=0.or.w(i)>np1 )  then ; print*, "error in read_dat: w<=0 or w>np1 !!" ; stop ; end if 
		end do 
	end if 
	end subroutine qx2hrwge
				
	subroutine getmiss(dat,nomiss)
	type(statevar), intent(in) :: dat ! data set. first entry is ia index, second observation number
	integer(i4b),intent(out) :: nomiss 
	nomiss=1
	if (dat%co<0) nomiss=0
	if (dat%sexr<0) nomiss=0
	if (dat%hme<0) nomiss=0		
	if (dat%endage<0) nomiss=0		
	if (dat%edr<0) nomiss=0
	!if (dat%expr<0) nomiss=0   !expr is alwyays missing in the data since there is no experience variable
	if (dat%rel==1.and.dat%kidr<0) nomiss=0
	if (dat%l<0) nomiss=0
	if (dat%hhr<0) nomiss=0
	if (dat%logwr<0.0_dp) nomiss=0
	if (dat%rel<0.0_dp) nomiss=0
	if (dat%rel==1.and.dat%rellen<0) nomiss=0
	if (dat%rel==1.and.dat%hhsp<0) nomiss=0    
	if (dat%rel==1.and.dat%logwsp<0.0_dp) nomiss=0
	!if (dat%edsp<0) nomiss=0   !edsp is always missing in the data because for now it is not read from the data
	!if (dat%expsp<0) nomiss=0  !expsp is alwyays missing in the data since there is no experience variable
	if (dat%kidsp<0) nomiss=0   !kidsp is always missing in the data since kidsp is just a simulation construct
	
	end subroutine getmiss 

	subroutine get_dathrwge(hr,inc_perhour,dhr,dw,dlogw)
	integer(i4b), intent(in) :: hr
	real(dp), intent(in) :: inc_perhour
	integer(i4b), intent(out) :: dhr
	real(dp), intent(out) :: dw,dlogw	
	!if (hr>0.and.wge>0) then 
        !note that the employment is decided according to whether hours is more than h_fulltime which is 1000 but annual wages is calculated using wge*h_wmult where h_wmult is 2000. 
        !This is how you did it in the original version. I do the same thing so that the wage numbers are consistent with the previous version. See page 15 last paragraph in the original text.
	!	dhr=one(hr>=h_fulltime)
	!	dwge=wge*h_wmult !ahu 021817: 1) got rid of the minmax thing 2) changed h_fulltime to h_wmult in the way annual wages is calculated      !min(max(minw,wge),maxw)*h_fulltime
	!	dwge=log(dwge)
	!	dwge=one(hr>=h_fulltime)*dwge
	!else if (hr==0.and.wge==0) then 
	!	dhr=hr 
	!	dwge=wge
	!else if (hr==-1.and.wge==-1) then 
	!	dhr=hr 
	!	dwge=wge
	!else 
	!	print*, "error in get_dathrwge ", hr,wge
	!	print*, "by construction, it should be that it's both 0, both -1 or both positive "
	!	stop
	!end if	
    
	if (hr>=0) then 
        !note that the employment is decided according to whether hours is more than h_fulltime which is 1000 but annual wages is calculated using wge*h_wmult where h_wmult is 2000. 
        !This is how you did it in the original version. I do the same thing so that the wage numbers are consistent with the previous version. See page 15 last paragraph in the original text.
		dhr=one(hr>=h_fulltime)
		if (inc_perhour>=minw.and.inc_perhour<=maxw) then
            dw=inc_perhour*h_wmult    
            dlogw=log(inc_perhour*h_wmult) !ahu 021817: 1) got rid of the minmax thing 2) changed h_fulltime to h_wmult in the way annual wages is calculated      !min(max(minw,wge),maxw)*h_fulltime            
        else 
            dw=-99.0_dp
            dlogw=-99.0_dp
        end if
    else if (hr<0) then 
        dhr=-99
        dw=-99.0_dp
        dlogw=-99.0_dp
	else 
		print*, "error in get_dathrwge ", hr,inc_perhour
		print*, "by construction, it should be that it's both 0, both -1 or both positive "
		stop
	end if	
    
    end subroutine get_dathrwge	
	
	subroutine get_mom(dat,nper,mom,cnt,var,name,headstr,headloc,weights)
	integer(i4b), intent(in) :: nper ! number of observations in dat array. nper is just the sample size for actual data and is (sample size)*nsimeach for sim data. 
	type(statevar), dimension(MNAD:MXA,nper), intent(in) :: dat ! data set. first entry is ia index, second is the observation number
	real(dp), dimension(nmom), intent(out) :: mom	 ! conditional mom
	integer(i4b), dimension(nmom), intent(out) :: cnt	 ! number of observations contributing to each moment
	real(dp), dimension(nmom), intent(out) :: var		 ! unconditional variance of each conditional moemnt
	character(len=namelen), dimension(nmom), intent(out) :: name ! names of mom
	character(len=500),dimension(nmom), intent(out) :: headstr ! headers for different sections of mom
	integer(i4b), dimension(nmom), intent(out) :: headloc	 ! location of headers for different sections of mom
	real(dp), dimension(nmom), intent(out) :: weights		 ! some real assoicated to each moment, maybe used as weights
	integer(i4b) :: ia,co,im,ihead,g,i,j,jj,ii,ddd,edo,ier,iesp,iloc,jloc,typosto,kk,mybuff,keeptrack
    integer(i4b) :: minsex,maxsex,maxrelo,INCR,LA,UA,agerange,checkco
    integer(i4b) :: maxwmA,maxwfA,maxlmA,maxlfA,maxwmB,maxwfB,maxlmB,maxlfB
	logical,dimension(MNAD:MXA,nper) :: cosexedtyp,cosex,cosexrel,corel,coho
    logical, dimension(nper) :: cosexedtypid
    type(initcond), dimension(nper) :: inito
	integer(i4b), dimension(MNAD:MXA,nper) :: norelchg,move,emph,empw,edh,edw,dur,dee,deue,moveadjacent
	real(dp), dimension(MNAD:MXA,nper) :: logwh,logww,wdif
    integer(i4b), dimension(nper) :: nummove,nummovemar,nummovesin    !nummove_ma,nummove_si !integer(i4b), dimension(MNAD:MXA,ndat) :: nummov,nummove_mar,nummove_sin
    real(dp) :: decilegrid(ndecile),wmove,tempsum(MNA:MXA),lifetimeinc(nper)
    INTEGER(I4B),DIMENSION(MNAD:MXA,nper) :: kidtrans,homemove,moverank,relid,nummove_relid 
    INTEGER(I4B), DIMENSION(0:40) :: NUMMOVETEMP !can be less than 40 as that is just for relid. so at most there would be like 30 relid's or something. 
    INTEGER(I4B),DIMENSION(nper) :: movesum
    REAL(dp),DIMENSION(MNAD:MXA,nper) :: mean4h,deltawage4,deltawage
    INTEGER(I4B), dimension(MNAD:MXA,nper) :: iacat,empr,nummovemar2,nummovesin2,jointemp
    INTEGER(I4B), dimension(5,MNAD:MXA,nper) :: etr
    LOGICAL,DIMENSION(MNAD:MXA,nper) :: obs4h

        !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!
        !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!
        !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!

    !initiate
	mom=-99.0_dp
	cnt=9999
	var=-99.0_dp
	headloc=-99 
	weights=0.0_dp 
    maxsex=-1
    maxrelo=-1
    coho=.FALSE.
    cosex=.FALSE.   ; cosexedtyp=.FALSE. ; cosexedtypid=.FALSE. ; coho=.FALSE.  ; inito=ones_init 
    cosexrel=.FALSE.
    corel=.FALSE.
    norelchg=-99
    move=-99 ; moveadjacent=-99
    emph=-99 
    empw=-99 
    edh=-99 
    edw=-99 
    dur=-99
    dee=-99
    deue=-99
    logwh=-99.0_dp
    logww=-99.0_dp ; wdif=-99.0_dp
    nummove=-99 
    jointemp=-99
    calcvar=0   !declared globally
    calcorr=0   !declared globally
    mominfo=-1  !declared globally
    momwhich=-1 !declared globally
	im=1
	ihead=1
    !decilegrid=0.0_dp
    !dur=0
    kidtrans=-99
    homemove=-99 
    wdif=-99.0_dp 

    headloc(ihead)=im
	if (skriv) call yaz_getmom(dat,nper) 
	
    do ddd=1,ndecile
        decilegrid(ddd)=8.0_dp+0.2_dp*(ddd-1)
	end do
    !print*, decilegrid
    
	WHERE ((dat(MNAD:MXAD,:)%rel>=0) .AND. dat(MNA:MXA,:)%rel>=0  )
		norelchg(MNAD:MXAD,:)=one(  dat(MNAD:MXAD,:)%rel==dat(MNA:MXA,:)%rel   ) 
	ENDWHERE
	WHERE ((dat(MNAD:MXAD,:)%l>0) .AND. (dat(MNA:MXA,:)%l>0))
		move(MNAD:MXAD,:)=one(dat(MNAD:MXAD,:)%l/=dat(MNA:MXA,:)%l)
	ENDWHERE
    do ia=MNAD,MXAD
        do i=1,nper
            if (   (dat(ia,i)%l>0) .AND. (dat(ia+1,i)%l>0)  ) then 
                iloc=dat(ia,i)%l ; jloc=dat(ia+1,i)%l 
                moveadjacent(ia,i)=one( iloc/=jloc .and. distance(iloc,jloc)==1)
            end if
        end do 
    end do 

	WHERE ((dat(MNAD:MXA,:)%rel==1) .AND.  (dat(MNAD:MXA,:)%sexr==1)  )
		emph(MNAD:MXA,:)=dat(MNAD:MXA,:)%hhr 
		empw(MNAD:MXA,:)=dat(MNAD:MXA,:)%hhsp
		edh(MNAD:MXA,:)=dat(MNAD:MXA,:)%edr 
		edw(MNAD:MXA,:)=dat(MNAD:MXA,:)%edsp
		logwh(MNAD:MXA,:)=dat(MNAD:MXA,:)%logwr 
		logww(MNAD:MXA,:)=dat(MNAD:MXA,:)%logwsp
    ENDWHERE
	WHERE ((dat(MNAD:MXA,:)%rel==1) .AND.  (dat(MNAD:MXA,:)%sexr==2)  )
		emph(MNAD:MXA,:)=dat(MNAD:MXA,:)%hhsp
		empw(MNAD:MXA,:)=dat(MNAD:MXA,:)%hhr
		edh(MNAD:MXA,:)=dat(MNAD:MXA,:)%edsp 
		edw(MNAD:MXA,:)=dat(MNAD:MXA,:)%edr
		logwh(MNAD:MXA,:)=dat(MNAD:MXA,:)%logwsp 
		logww(MNAD:MXA,:)=dat(MNAD:MXA,:)%logwr
    ENDWHERE

    !dee for wdif conditionoing on ee transitions 
    WHERE (dat(MNAD:MXAD,:)%hhr==1 .AND. dat(MNA:MXA,:)%hhr==1 .AND. dat(MNAD:MXAD,:)%logwr>=0 .AND. dat(MNA:MXA,:)%logwr>=0  ) 
        dee(MNAD:MXAD,:)=1
        wdif(MNAD:MXAD,:) =  dat(MNA:MXA,:)%logwr-dat(MNAD:MXAD,:)%logwr
    ENDWHERE
    
    !deue for wdif conditioning on eue transitions
    WHERE ( dat(MNA-1:MXA-2,:)%hhr==1 .AND. dat(MNA:MXA-1,:)%hhr==0  .AND. dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA-1:MXA-2,:)%logwr>=0   .AND. dat(MNA+1:MXA,:)%logwr>=0  ) 
        deue(MNA-1:MXA-2,:)=1
    ENDWHERE
 
	WHERE ((dat(MNAD:MXA,:)%hhr==0)   )
		dur(MNAD:MXA,:)=1
    ENDWHERE

    !kidtrans=-1.
    WHERE (  (dat(MNAD:MXAD,:)%rel>=0) .AND. (dat(MNAD:MXAD,:)%kidr==1) .AND. (dat(MNA:MXA,:)%kidr>=1)  )
    kidtrans(MNAD:MXAD,:)=one(dat(MNA:MXA,:)%kidr==2) !kid=1 is nokid and kid=2 is yeskid
    ENDWHERE
    WHERE (move(MNAD:MXAD,:)==1)
    homemove(MNAD:MXAD,:)=one(dat(MNA:MXA,:)%hme==dat(MNA:MXA,:)%l)
    ENDWHERE

    movesum(:)=sum(move(MNAD:MXAD,:),1,move(MNAD:MXAD,:)>=0)
    moverank=-1
    do ia=MNAD,MXAD
    moverank(ia,:)=sum(move(MNAD:ia,:),1,move(MNAD:ia,:)>=0)   
    end do 

    !ddd=1
    !do i=5,10
    !    print*, 'Here it is ', d1*one( (dat(25,i)%logwr>decilegrid(ddd) .and. dat(25,i)%logwr<decilegrid(ddd+1) ) )
    !end do     
    lifetimeinc(:)=0.0_dp ; tempsum=0.0_dp
    !if (momdisplay) then
        do i=1,nper
            do ia=MNA,MXA
                if (        (dat(ia,i)%hhr==1)      ) then
                    tempsum(ia)=(dat(ia,i)%wr) * ( delta**(ia-MNA) )
                else 
                    tempsum(ia)=0.0_dp
                end if
            end do 
            lifetimeinc(i)=sum(tempsum(MNA:MXA))
        end do 
    !end if 

    do i=1,nper
        do ia=MNAD,MXAD
	        if (        (dat(ia,i)%hhr==0)    .AND. (dat(ia+1,i)%hhr==0)   ) then
		        dur(ia+1,i)=dur(ia,i)+1
            end if 
            if ( dat(ia,i)%sexr==1) then 
                if (        (dat(ia,i)%hhr==1)    .AND. (dat(ia,i)%hhsp==1)   ) then
                    jointemp(ia,i)=1
                else if (        (dat(ia,i)%hhr==1)    .AND. (dat(ia,i)%hhsp==0)   ) then
                    jointemp(ia,i)=2                    
                else if (        (dat(ia,i)%hhr==0)    .AND. (dat(ia,i)%hhsp==1)   ) then
                    jointemp(ia,i)=3
                else if (        (dat(ia,i)%hhr==0)    .AND. (dat(ia,i)%hhsp==0)   ) then
                    jointemp(ia,i)=4
                end if 
            else if ( dat(ia,i)%sexr==2) then 
                if (        (dat(ia,i)%hhr==1)    .AND. (dat(ia,i)%hhsp==1)   ) then
                    jointemp(ia,i)=1
                else if (        (dat(ia,i)%hhr==0)    .AND. (dat(ia,i)%hhsp==1)   ) then
                    jointemp(ia,i)=2                    
                else if (        (dat(ia,i)%hhr==1)    .AND. (dat(ia,i)%hhsp==0)   ) then
                    jointemp(ia,i)=3
                else if (        (dat(ia,i)%hhr==0)    .AND. (dat(ia,i)%hhsp==0)   ) then
                    jointemp(ia,i)=4
                end if 
            end if
        end do 
    end do 

    relid=-100 ; nummove_relid=-999 ; nummovetemp=0
    do i=1,nper
        !Just generating cohort and sex of someone because don't have initcond in this routine
        !cohogen(i)=  inito(i)%co    !maxval(dat(MNAD:MXA,i)%co)
        !sexgen(i)=   inito(i)%sex   !maxval(dat(MNAD:MXA,i)%sexr)	                
        inito(i)%co = maxval(dat(MNAD:MXA,i)%co)
        inito(i)%sexr = maxval(dat(MNAD:MXA,i)%sexr)
        inito(i)%edr = maxval(dat(MNAD:MXA,i)%edr)
        inito(i)%typ = maxval(dat(MNAD:MXA,i)%typ)
        !Total number of moves overall
        !nummove(i)=sum(move(MNA:MXAI,i),  MASK=( move(MNA:MXAI,i)>0  .and. norelchg(MNA:MXAI,i)==1  )   )
        nummove(i)=sum(move(MNAD:MXAD,i),  MASK=( move(MNAD:MXAD,i)>=0  .and. norelchg(MNAD:MXAD,i)==1  )   )
        nummovemar(i)=sum(move(MNAD:MXAD,i),  MASK=( move(MNAD:MXAD,i)>=0  .and. norelchg(MNAD:MXAD,i)==1 .and. dat(MNAD:MXAD,i)%rel==1 )   )
        nummovesin(i)=sum(move(MNAD:MXAD,i),  MASK=( move(MNAD:MXAD,i)>=0  .and. norelchg(MNAD:MXAD,i)==1 .and. dat(MNAD:MXAD,i)%rel==0 )   )
        keeptrack=0 !reset relid tracker
        do ia=MNA,MXA
            if (dat(ia,i)%rellen==1.and.dat(ia,i)%rel==1) then 
                keeptrack=keeptrack+1
                relid(ia,i)=keeptrack
            else if (dat(ia,i)%rellen>1.and.dat(ia,i)%rel==1) then  
                relid(ia,i)=keeptrack
            !else if (dat(ia,i)%rel==0) then
            !    relid(ia,i)=-100
            end if
        end do 
        nummovetemp=0
        do ia=MNA,MXA
            if (dat(ia,i)%rel==1.and.move(ia,i)==1.and.norelchg(ia,i)==1.and.dat(ia,i)%rellen>=1) then 
                nummovetemp(relid(ia,i))=nummovetemp(relid(ia,i))+1
            end if
        end do 
        do ia=MNA,MXA
            if (relid(ia,i)>0) then 
                nummove_relid(ia,i)=nummovetemp(relid(ia,i))
            end if
        end do 
    end do 

    if (nper==numpersim.and.iter==1) then 
        nummove_save(:)=nummove(:)
        nummovemar_save(:)=nummovemar(:)
        nummovesin_save(:)=nummovesin(:)
    end if 
    if (nper==numpersim.and.iter>1) then 
        nummove(:)=nummove_save(:)
        nummovemar(:)=nummovemar_save(:)
        nummovesin(:)=nummovesin_save(:)
    end if 

    do ia=MNAD,MXA
        nummovemar2(ia,:)=nummovemar(:)
        nummovesin2(ia,:)=nummovesin(:)
    end do 

    ddd=1
	if (skriv.and.nper==numpersim.and. iter==1) then
	open(unit=94632, file='dursim1.txt',status='replace')
    do j=1,nper
        do ia=mnad,mxa
            if (dat(ia,j)%hme==1) then 
                write(94632,'(6i6,F10.2,3I6)' ) j,dat(ia,j)%nn,dat(ia,j)%mm,ia,dat(ia,j)%L,dat(ia,j)%hhr,dat(ia,j)%logwr,dat(ia,j)%EXPR,dat(ia,j)%HME,dat(ia,j)%sexr   !, d1*one( (dat(ia,j)%logwr>decilegrid(ddd) .and. dat(ia,j)%logwr<decilegrid(ddd+1) ) ), decilegrid(ddd),decilegrid(ddd+1)   !,dat(ia,j)%nn,dat(ia,j)%mm,dat(ia,j)%r
	        end if 
        end do 
	end do
    close(94632)
    end if 

    ddd=1
	if (skriv.and.nper==numperdat .and. iter==2) then
	open(unit=94632, file='dursim2.txt',status='replace')
    do j=1,nper
        do ia=mnad,mxa
            if (dat(ia,j)%hme==1) then 
                write(94632,'(6i6,F10.2,3I6)' ) j,dat(ia,j)%nn,dat(ia,j)%mm,ia,dat(ia,j)%L,dat(ia,j)%hhr,dat(ia,j)%logwr,dat(ia,j)%EXPR,dat(ia,j)%HME,dat(ia,j)%sexr   !, d1*one( (dat(ia,j)%logwr>decilegrid(ddd) .and. dat(ia,j)%logwr<decilegrid(ddd+1) ) ), decilegrid(ddd),decilegrid(ddd+1)   !,dat(ia,j)%nn,dat(ia,j)%mm,dat(ia,j)%r
	        end if 
        end do 
	end do
    close(94632)
    end if 

    if ( nper==numpersim.and. iter==1.AND.(.NOT.SAVETIME).and.DETAILEDOUTPUT==10) then !.AND.DETAILEDOUTPUT==1) then 
        mybuff=94630 + mysay
        if (mysay==0)  open(unit=mybuff, file='checkdat0.txt',status='replace')
        if (mysay==1)  open(unit=mybuff, file='checkdat1.txt',status='replace')
        if (mysay==2)  open(unit=mybuff, file='checkdat2.txt',status='replace')
        if (mysay==3)  open(unit=mybuff, file='checkdat3.txt',status='replace')
        if (mysay==4)  open(unit=mybuff, file='checkdat4.txt',status='replace')
        if (mysay==5)  open(unit=mybuff, file='checkdat5.txt',status='replace')
        if (mysay==6)  open(unit=mybuff, file='checkdat6.txt',status='replace')
        if (mysay==7)  open(unit=mybuff, file='checkdat7.txt',status='replace')
        if (mysay==8)  open(unit=mybuff, file='checkdat8.txt',status='replace')
        if (mysay==9)  open(unit=mybuff, file='checkdat9.txt',status='replace')
        if (mysay==10) open(unit=mybuff, file='checkdat10.txt',status='replace')
        if (mysay==11) open(unit=mybuff, file='checkdat11.txt',status='replace')
        if (mysay==12) open(unit=mybuff, file='checkdat12.txt',status='replace')
        if (mysay==13) open(unit=mybuff, file='checkdat13.txt',status='replace')
        if (mysay==14) open(unit=mybuff, file='checkdat14.txt',status='replace')
        if (mysay==15) open(unit=mybuff, file='checkdat15.txt',status='replace')
        if (mysay==16) open(unit=mybuff, file='checkdat16.txt',status='replace')
        if (mysay==17) open(unit=mybuff, file='checkdat17.txt',status='replace')
        if (mysay==18) open(unit=mybuff, file='checkdat18.txt',status='replace')
        if (mysay==19) open(unit=mybuff, file='checkdat19.txt',status='replace')
        if (mysay==20) open(unit=mybuff, file='checkdat20.txt',status='replace')
        if (mysay==21) open(unit=mybuff, file='checkdat21.txt',status='replace')
        if (mysay==22) open(unit=mybuff, file='checkdat22.txt',status='replace')
        if (mysay==23) open(unit=mybuff, file='checkdat23.txt',status='replace')
        if (mysay==24) open(unit=mybuff, file='checkdat24.txt',status='replace')
        if (mysay==25) open(unit=mybuff, file='checkdat25.txt',status='replace')
        if (mysay==26) open(unit=mybuff, file='checkdat26.txt',status='replace')
        if (mysay==27) open(unit=mybuff, file='checkdat27.txt',status='replace')
        if (mysay==28) open(unit=mybuff, file='checkdat28.txt',status='replace')
        if (mysay==29) open(unit=mybuff, file='checkdat29.txt',status='replace')
        if (mysay==30) open(unit=mybuff, file='checkdat30.txt',status='replace')
        if (mysay==31) open(unit=mybuff, file='checkdat31.txt',status='replace')
        if (mysay==32) open(unit=mybuff, file='checkdat32.txt',status='replace')
        if (mysay==33) open(unit=mybuff, file='checkdat33.txt',status='replace')
        if (mysay==34) open(unit=mybuff, file='checkdat34.txt',status='replace')
        if (mysay==35) open(unit=mybuff, file='checkdat35.txt',status='replace')
        !write(mybuff,'(tr5,"j",tr4,"nn",tr4,"mm",tr4,"ia",tr3,"hhr",tr2,"expr",tr3,"job",tr3,"swi",tr5,"logwr",tr6,"wdif",tr1,"hme",tr1,"sex",tr1,"edr",tr3,"l",tr10,"emax",tr8,"emaxsp",tr9,"consr",tr9,"consp",tr8,"consum",tr8,"incsum",tr12,"wr",tr11,"wsp" ,tr1,"lifeinc",tr1,"adj",tr1,"rel",tr1,"ubenefit",tr1,"hsp",tr1,"nurelid",tr4,"logwsp",tr1,"kid",tr1,"typ",tr2,"it",tr3,"norel",tr4,"move",tr5,"num",tr2,"nummar",tr2,"numsin",tr4,"mina",tr4,"maxa",tr4,"edsp",tr3,"relid",tr2,"rellen")') 
    end if 

    if ( nper==numpersim .AND.(.NOT.SAVETIME).and.DETAILEDOUTPUT==10) then
        mybuff=94630 + mysay
        do j=1,nper
            do ia=mnad,mxa
                if (dat(ia,j)%hme>0) then 
                    maxwmA=-99 ; maxwfA=-99 ; maxlmA=-99  ; maxlfA=-99
                    maxwmB=-99 ; maxwfB=-99 ; maxlmB=-99  ; maxlfB=-99
                    ddd=myefo(ia,j)%maxA
                    if (ddd>0) then ; maxwmA=qq2w(1,ddd) ; maxwfA=qq2w(2,ddd) ; maxlmA=qq2l(1,ddd)  ; maxlfA=qq2l(2,ddd)   ; end if
                    ddd=myefo(ia,j)%maxB
                    if (ddd>0) then ; maxwmB=qq2w(1,ddd) ; maxwfB=qq2w(2,ddd) ; maxlmB=qq2l(1,ddd)  ; maxlfB=qq2l(2,ddd)   ; end if                        
                    !write(mybuff,'(8i6,2F10.1,4I4,8F14.0,i8,i4,3I4,2F9.1,3I4,8I8)' ) j,dat(ia,j)%nn,dat(ia,j)%mm,ia,dat(ia,j)%hhr,dat(ia,j)%EXPR,dat(ia,j)%job,dat(ia,j)%jobswitch,dat(ia,j)%logwr,wdif(ia,j),dat(ia,j)%HME,dat(ia,j)%sexr,dat(ia,j)%edr,dat(ia,j)%L,dat(ia,j)%emax,dat(ia,j)%emaxsp,dat(ia,j)%consr,dat(ia,j)%consp,dat(ia,j)%consum,dat(ia,j)%incsum,dat(ia,j)%wr,dat(ia,j)%wsp,j,ia,dat(ia,j)%rel,dat(ia,j)%hhr,dat(ia,j)%hhsp,dat(ia,j)%logwr,dat(ia,j)%logwsp,dat(ia,j)%sexr,& 
                    !& dat(ia,j)%typ,iter ,norelchg(ia,j),move(ia,j),nummove(j),nummovemar(j),nummovesin(j),dat(ia,j)%minage,dat(ia,j)%endage,dat(ia,j)%edsp
                    write(mybuff,'(8I6,2F10.1,4I4,9F14.0,2I4,F9.1,I4,I8,F9.1,3I4,10I8,I8,7I4,4I8,10I4,2I8,2F10.1,8I4,I8,2I4,I8,2I4)' ) j,dat(ia,j)%nn,dat(ia,j)%mm,ia,dat(ia,j)%hhr,dat(ia,j)%EXPR,dat(ia,j)%job,dat(ia,j)%jobswitch &
                    &,dat(ia,j)%logwr,wdif(ia,j)&
                    &,dat(ia,j)%HME,dat(ia,j)%sexr,dat(ia,j)%edr,dat(ia,j)%L&
                    &,dat(ia,j)%emax,dat(ia,j)%emaxsp,dat(ia,j)%consr,dat(ia,j)%consp,dat(ia,j)%consum,dat(ia,j)%incsum,dat(ia,j)%wr,dat(ia,j)%wsp,lifetimeinc(j)&
                    &,moveadjacent(ia,j),dat(ia,j)%rel &
                    &,ubenefitbyloc(dat(ia,j)%l),dat(ia,j)%hhsp,nummove_relid(ia,j)&
                    &,dat(ia,j)%logwsp&
                    &,dat(ia,j)%kidr,dat(ia,j)%typ,iter &
                    &,norelchg(ia,j),move(ia,j),nummove(j),nummovemar(j),nummovesin(j),dat(ia,j)%minage,dat(ia,j)%endage,dat(ia,j)%edsp,relid(ia,j),dat(ia,j)%rellen &
                    &,draws(ia,j)%qshock &
                    &,draws(ia,j)%wmshock,draws(ia,j)%wfshock,draws(ia,j)%lshock,draws(ia,j)%iephub,draws(ia,j)%iepwfe,draws(ia,j)%iepsingle,draws(ia,j)%zshock &
                    &,decisions(ia,j)%q0dec,decisions(ia,j)%qdec,decisions(ia,j)%qmdec,decisions(ia,j)%qfdec &
                    &,decisions(ia,j)%rel0dec,decisions(ia,j)%reldec &
                    &,decisions(ia,j)%wm0dec,decisions(ia,j)%wf0dec &
                    &,decisions(ia,j)%wmdec,decisions(ia,j)%wfdec &
                    &,decisions(ia,j)%l0dec,decisions(ia,j)%ldec,decisions(ia,j)%lmdec,decisions(ia,j)%lfdec  &
                    &,myefo(ia,j)%maxA,myefo(ia,j)%maxB &
                    &,myefo(ia,j)%tra1,myefo(ia,j)%tra2 &                     
                    &,maxwmA,maxwfA,maxlmA,maxlfA &
                    &,maxwmB,maxwfB,maxlmB,maxlfB &
                    &,mdecisions(ia,j)%qmdec,mdecisions(ia,j)%wmdec,mdecisions(ia,j)%lmdec &
                    &,fdecisions(ia,j)%qfdec,fdecisions(ia,j)%wfdec,fdecisions(ia,j)%lfdec 
                    !&,myefo(1:nc,ia,j)%sur,myefo(1:nc,ia,j)%haveenough,myefo(1:nc,ia,j)%haveenoughforNB,myefo(1,ia,j)%maxA,myefo(1,ia,j)%maxB,myefo(1:2,ia,j)%tra  
                    !&,myefo(1:nc,ia,j)%ve1,myefo(1:nc,ia,j)%ve2,myefo(1:nc,ia,j)%ve3,myefo(1:nc,ia,j)%ve4,myefo(1:nc,ia,j)%ve5 & 
                    !, d1*one( (dat(ia,j)%logwr>decilegrid(ddd) .and. dat(ia,j)%logwr<decilegrid(ddd+1) ) ), decilegrid(ddd),decilegrid(ddd+1)   !,dat(ia,j)%nn,dat(ia,j)%mm,dat(ia,j)%r
                end if 
            end do 
        end do
    end if 

    if ( nper==numpersim .and. iter>=numit .AND.(.NOT.SAVETIME).and.DETAILEDOUTPUT==10) then
        mybuff=94630 + mysay
        close(mybuff)
    end if 


!		if (dat(ia,j)%expr>-99.and.dat(ia,j)%expr<0) then 
	!			print*, 'dat exp is negative',ia,j,dat(ia,j)%expr,dat(ia,j)%hme,dat(ia,j)%nn,dat(ia,j)%mm,dat(ia,j)%r
	!			stop
	!		end if
    
	!if ( iter==1) then
	!do ia=mna,mxa
	!	do j=1,ndat
            !The below check does give you an error. Despite the fact that everyone starts from age 18 and there is no missing cohort variable in the data (all is 1920-49).
            !Why?
            !This is because, for those people who are ED, the agestart is 22. So in read_actualdata, after reading these people's initcond's from their age 18 (and for all of them I have all those initconds), 
            !I set all the variables from 18 to 21 to -99, including cohort (after reading cohort info into init). Hence for such people, you
            !will have a cohort variable value of -99 at age 18-21. But starting from 22, you will have the cohort value. 
            !For some, there are some ages for which there is no observation in the PSID. For example, for idnum who is ED, the person is observed from age 18-20
            !but then there is nothing until 24. So for this person, the cohort values (as well as other variable values) would be all filled in with -99 until 24. 
            !SO instead of conditioning the nummove moments on dat(mna,:)%co=co, I condition them on coho(mna,:).			
            !(I generate coho with		coho(mna:mxai,:)= (  maxval(dat(mna:mxai,:)%co)==co  )			
            !This solves the problem. 
            !if (dat(ia,j)%co<1.or.dat(ia,j)%co>nco ) then
                !print*, 'something wrong with dat%co A',ia,j,dat(ia,j)%co,ndat
                !stop
            !end if
            !if ( maxval(dat(:,j)%co) /= dat(ia,j)%co  ) then
                !print*, 'something wrong with dat%co B',ia,j,dat(ia,j)%co,ndat
                !stop
            !end if
    !    end do 
	!end do
	!end if 
    
    !if (onlysingles) then ! if doing only singles, condition all mom on singles 
	!	nomiss(mna:mxai-1,:) = (  nomiss(mna:mxai-1,:) .and. dat(mna:mxai-1,:)%rel==0 .and. dat(mna+1:mxai,:)%rel==0 ) 
	!end if 

    if (onlysingles) then 
        maxrelo=0
    else 
        maxrelo=1
    end if 
    if (onlymales.and.(.not.onlyfem)) then 
        minsex=1 ; maxsex=1
    else if ( (.not.onlymales).and.onlyfem) then 
        minsex=2 ; maxsex=2
    else if ( (.not.onlymales).and.(.not.onlyfem)) then 
        minsex=1 ; maxsex=2
    else 
        print*, 'Something is wrong in the neighborhood'
        stop
    end if 
    
    !if (onlymales) then 
    !    maxsex=1
    !else 
    !    maxsex=2
    !end if 
    
    !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!
    !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!
    !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!

    !ag092922 sept2022: note that I got rid of all the if stateements following weights(im) i.e. if(onlysingles and j==1) weights(im)=0.0_dp
    !got rid of them because I don't need them if I have maxrel0 since maxrel0 is set to 1 if onlysingles in the above if statement before the co do loop starts. 
    !this if (onlysingles and j==1) was leading to bugs i.e. different processors would have different momwgt (all else same tho at least)
    !because I wasn't assigning j in the loops for just co or cosex (j is for rel) and when I moved around the moments do loops I forgot that j there 
    !so it was being assigned something different for each processor and they were each gibing different momwgt because of that. 
    !but just get rid of onlysingles if statements after weights(im) since you don't need that anyway (due tot he presence of maxrel)
	cohort: do co=1,nco
    
        headloc(ihead)=im
        headstr(ihead)='everyone rates '
        cosexedtypid(:)            = (inito(:)%co==co  )			
        cosexedtyp(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co  )		


        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummove(:)>=0 ),d1*one(nummove(:)==0),mom,cnt,var)
        write(name(im),'("nummove=0 ",tr10)')	
        weights(im)=wmovenum 
        im=im+1
        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummove(:)>=0 ),d1*one(nummove(:)==1),mom,cnt,var)
        write(name(im),'("nummove=1 ",tr10)')	
        weights(im)=wmovenum 
        im=im+1
        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummove(:)>=0 ),d1*one(nummove(:)==2),mom,cnt,var)
        write(name(im),'("nummove=2 ",tr10)')	
        weights(im)=wmovenum
        im=im+1
        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummove(:)>=0 ),d1*one(nummove(:)>=3),mom,cnt,var)
        write(name(im),'("nummove>=3 ",tr10)')	
        weights(im)=wmovenum 
        im=im+1

        do g=minsex,maxsex
            CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%sexr==g.AND. nummove(:)>=0 ),d1*one(nummove(:)==0),mom,cnt,var)
            write(name(im),'("nummove=0 by sex ",I4)') g	
            weights(im)=wmovenum 
            im=im+1
            CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%sexr==g.AND. nummove(:)>=0 ),d1*one(nummove(:)==1),mom,cnt,var)
            write(name(im),'("nummove=1 by sex ",I4)') g	
            weights(im)=wmovenum 
            im=im+1
            CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%sexr==g.AND. nummove(:)>=0 ),d1*one(nummove(:)==2),mom,cnt,var)
            write(name(im),'("nummove=2 by sex ",I4)') g	
            weights(im)=wmovenum
            im=im+1
            CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%sexr==g .AND. nummove(:)>=0 ),d1*one(nummove(:)>=3),mom,cnt,var)
            write(name(im),'("nummove>=3 by sex ",I4)') g	
            weights(im)=wmovenum 
            im=im+1
        end do 

        do ier=1,2
            CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%edr==ier.AND. nummove(:)>=0 ),d1*one(nummove(:)==0),mom,cnt,var)
            write(name(im),'("nummove=0 by ed ",I4)') ier
            weights(im)=wmovenum 
            im=im+1
            CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%edr==ier.AND. nummove(:)>=0 ),d1*one(nummove(:)==1),mom,cnt,var)
            write(name(im),'("nummove=1 by ed ",I4)') ier	
            weights(im)=wmovenum 
            im=im+1
            CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%edr==ier.AND. nummove(:)>=0 ),d1*one(nummove(:)==2),mom,cnt,var)
            write(name(im),'("nummove=2 by ed ",I4)') ier	
            weights(im)=wmovenum
            im=im+1
            CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%edr==ier.AND. nummove(:)>=0 ),d1*one(nummove(:)>=3),mom,cnt,var)
            write(name(im),'("nummove>=3 by ed ",I4)') ier	
            weights(im)=wmovenum 
            im=im+1
        end do 

        do g=minsex,maxsex
            do ier=1,2
                CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%sexr==g.AND. inito(:)%edr==ier.AND. nummove(:)>=0 ),d1*one(nummove(:)==0),mom,cnt,var)
                write(name(im),'("nummove=0 by sex and ed ",2I4)') g,ier
                weights(im)=wmovenum 
                im=im+1
                CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%sexr==g.AND. inito(:)%edr==ier.AND. nummove(:)>=0 ),d1*one(nummove(:)==1),mom,cnt,var)
                write(name(im),'("nummove=1 by sex and ed ",2I4)') g,ier	
                weights(im)=wmovenum 
                im=im+1
                CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%sexr==g.AND. inito(:)%edr==ier.AND. nummove(:)>=0 ),d1*one(nummove(:)==2),mom,cnt,var)
                write(name(im),'("nummove=2 by sex and ed ",2I4)') g,ier	
                weights(im)=wmovenum
                im=im+1
                CALL condmom(im,( COSEXEDTYPID(:) .AND. inito(:)%sexr==g.AND. inito(:)%edr==ier.AND. nummove(:)>=0 ),d1*one(nummove(:)>=3),mom,cnt,var)
                write(name(im),'("nummove>=3 by sex and ed ",2I4)') g,ier	
                weights(im)=wmovenum 
                im=im+1
            end do 
        end do 

        
        CALL condmom(im,( COSEXEDTYP(:,:) .AND. nummove_relid(:,:)>=0 .AND. dat(:,:)%rel==1 .AND. dat(:,:)%rellen==1 ),d1*one(nummove_relid(:,:)==0),mom,cnt,var)
        write(name(im),'("nummove relid=0 ",tr10)')	 ; weights(im)=wmovenum 
        im=im+1
        CALL condmom(im,( COSEXEDTYP(:,:) .AND. nummove_relid(:,:)>=0 .AND. dat(:,:)%rel==1 .AND. dat(:,:)%rellen==1 ),d1*one(nummove_relid(:,:)==1),mom,cnt,var)
        write(name(im),'("nummove relid=1 ",tr10)')	 ; weights(im)=wmovenum 
        im=im+1
        CALL condmom(im,( COSEXEDTYP(:,:) .AND. nummove_relid(:,:)>=0 .AND. dat(:,:)%rel==1 .AND. dat(:,:)%rellen==1 ),d1*one(nummove_relid(:,:)>1),mom,cnt,var)
        write(name(im),'("nummove relid>1 ",tr10)')	 ; weights(im)=wmovenum 
        im=im+1
        do i=1,ntypp
            CALL condmom(im,( COSEXEDTYP(:,:) .AND. nummove_relid(:,:)>=0 .AND. dat(:,:)%rel==1 .AND. dat(:,:)%rellen==1 .AND. dat(:,:)%typ==i ),d1*one(nummove_relid(:,:)==0),mom,cnt,var)
            write(name(im),'("nummove relid=0 typ ",i4)') i 	 ; weights(im)=wmovenum 
            im=im+1
            CALL condmom(im,( COSEXEDTYP(:,:) .AND. nummove_relid(:,:)>=0 .AND. dat(:,:)%rel==1 .AND. dat(:,:)%rellen==1 .AND. dat(:,:)%typ==i),d1*one(nummove_relid(:,:)==1),mom,cnt,var)
            write(name(im),'("nummove relid=1 typ",i4)') i	 ; weights(im)=wmovenum 
            im=im+1
            CALL condmom(im,( COSEXEDTYP(:,:) .AND. nummove_relid(:,:)>=0 .AND. dat(:,:)%rel==1 .AND. dat(:,:)%rellen==1 .AND. dat(:,:)%typ==i),d1*one(nummove_relid(:,:)>1),mom,cnt,var)
            write(name(im),'("nummove relid>1 typ",i4)') i	  ; weights(im)=wmovenum 
            im=im+1
        end do 


        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummovemar(:)>=0 ),d1*one(nummovemar(:)==0),mom,cnt,var)
        write(name(im),'("nummovemar=0 ",tr10)')	
        weights(im)=wmovenum 
        im=im+1
        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummovemar(:)>=0 ),d1*one(nummovemar(:)==1),mom,cnt,var)
        write(name(im),'("nummovemar=1 ",tr10)')	
        weights(im)=wmovenum 
        im=im+1
        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummovemar(:)>=0 ),d1*one(nummovemar(:)==2),mom,cnt,var)
        write(name(im),'("nummovemar=2 ",tr10)')	
        weights(im)=wmovenum
        im=im+1
        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummovemar(:)>=0 ),d1*one(nummovemar(:)>=3),mom,cnt,var)
        write(name(im),'("nummovemar>=3 ",tr10)')	
        weights(im)=wmovenum 
        im=im+1
        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummovesin(:)>=0 ),d1*one(nummovesin(:)==0),mom,cnt,var)
        write(name(im),'("nummovesin=0 ",tr10)')	
        weights(im)=wmovenum 
        im=im+1
        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummovesin(:)>=0 ),d1*one(nummovesin(:)==1),mom,cnt,var)
        write(name(im),'("nummovesin=1 ",tr10)')	
        weights(im)=wmovenum 
        im=im+1
        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummovesin(:)>=0 ),d1*one(nummovesin(:)==2),mom,cnt,var)
        write(name(im),'("nummovesin=2 ",tr10)')	
        weights(im)=wmovenum
        im=im+1
        CALL condmom(im,( COSEXEDTYPID(:) .AND. nummovesin(:)>=0 ),d1*one(nummovesin(:)>=3),mom,cnt,var)
        write(name(im),'("nummovesin>=3 ",tr10)')	
        weights(im)=wmovenum 
        im=im+1


        if (momdisplay) then 
            do ia=MNAD,MNA+1
                call condmom(im,( COSEXEDTYP(ia,:)  ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                write(name(im),'("emax ",tr5,i4)') ia  
                weights(im)=0.0_dp
                im=im+1
            end do
            do ia=19,22
                call condmom(im,( COSEXEDTYP(ia,:) ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                write(name(im),'("emax ",tr5,i4)') ia  
                weights(im)=0.0_dp
                im=im+1
            end do
            do ia=MNAD,MXA,10
                call condmom(im,( COSEXEDTYP(ia,:) ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                write(name(im),'("emax ",tr5,i4)') ia  
                weights(im)=0.0_dp
                im=im+1
            end do
        end if 
        CALL condmom(im,(   COSEXEDTYP(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*dat(mna:mxa,:)%logwr,mom,cnt,var)
        WRITE(name(im),'("wage all ",tr1)') 
        weights(im)=wwage 
        calcvar(im)=1
        im=im+1
        CALL condmom(im,(   COSEXEDTYP(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*  (dat(mna:mxa,:)%logwr**2),mom,cnt,var)
        WRITE(name(im),'("wvar all ",tr1)') 
        weights(im)=wwvar 
        calcvar(im)=5
        im=im+1
        CALL condmom(im,(   COSEXEDTYP(mna:mxa,:) ) ,d1*one( dat(mna:mxa,:)%hhr==1 ) ,mom,cnt,var)
        WRITE(name(im),'("emp all ",tr1)') 
        weights(im)=whour
        im=im+1

        do ia=MNAD,21,1
            CALL condmom(im,( COSEXEDTYP(ia,:) .AND. move(ia,:)>=0 ),d1*move(ia,:),mom,cnt,var)
            WRITE(name(im),'("move by age",tr3,I4)') ia
            weights(im)=wmove
            im=im+1
        end do
        do ia=22,50,10
            CALL condmom(im,( COSEXEDTYP(ia,:) .AND. move(ia,:)>=0 ),d1*move(ia,:),mom,cnt,var)
            WRITE(name(im),'("move by age",tr3,I4)') ia
            weights(im)=wmove
            im=im+1
        end do     
        call condmom(im,( COSEXEDTYP(MNAD:MXAD,:) .AND. move(MNAD:MXAD,:)==1 ),   d1*one( dat(MNA:MXA,:)%l==dat(MNA:MXA,:)%hme  ),mom,cnt,var)		
        write(name(im),'("move home ",tr3)')  !added this ahu 121718
        weights(im)=whome
        im=im+1 
        !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!
        !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!
        !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!
        !home dist ned      18     0.9780        0.9216        0.4401        0.0000          32670      1250       
        !home dist ed      18      0.9688        0.8084        1.6618        0.0000          37280      1748       
        !. count if agec==18 & edr==1          1,250
        !tab athome if agec==18 & edr==1
        !     athome |      Freq.     Percent        Cum.
        !------------+-----------------------------------
        !          0 |         84        6.87        6.87
        !          1 |      1,138       93.13      100.00
        !*------------+-----------------------------------
        !      Total |      1,222      100.00
        
        !count if agec==18 & edr==2
        !  1,748
        !tab athome if agec==18 & edr==2
        !     athome |      Freq.     Percent        Cum.
        !------------+-----------------------------------
        !          0 |        107        7.13        7.13
        !          1 |      1,393       92.87      100.00
        !------------+-----------------------------------
        !      Total |      1,500      100.00
        !This discrepancy between moments and data arose because I was not conditinoning on loc>0 in these homedist moments. And there ARE loc missing cases even though I dropped everstatemis cases because these loc missing are relhdgen<0 cases. I did 
        !had NOT dropped those in order to still beable to keep mover out nonresponse observations and I think that's appropriate. 
        !So there is nothing wrong with the data, but what is wrong here is that I Was not conditinoing properly.
        !And this discrepancy did make a big difference. For example, in model moments the home dist for ed people at age 18 was 80 percent. whereas in the data it is 90 something percent. 
        !Becaise the conditioning cell in the model moments case was bigger (and also included loc missing situations), the proportion of cases with loc=home at 18 seemed to be lower in th esimulation. 
        !*****************************************************************************************                    
        do ia=agestart(NOCOLLEGE)-1,19
            call condmom(im,( COSEXEDTYP(ia,:) .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%l>0 .AND. dat(ia,:)%hme>0 ),   d1*one( dat(ia,:)%l==dat(ia,:)%hme  ),mom,cnt,var)		
            write(name(im),'("home dist ned ",tr3,i4)') ia !added this ahu 121718
            weights(im)=whome
            im=im+1 
        end do
        do ia=agestart(COLLEGE)-1,24
            call condmom(im,( COSEXEDTYP(ia,:)  .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%l>0 .AND. dat(ia,:)%hme>0 ),   d1*one( dat(ia,:)%l==dat(ia,:)%hme  ),mom,cnt,var)		
            write(name(im),'("home dist ed ",tr3,i4)') ia !added this ahu 121718
            weights(im)=whome
            im=im+1 
        end do
        !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!
        !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!
        !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!


        LA=MNA ; UA=MXAD
        do jloc=1,nl         
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA+1:UA+1,:)%l==jloc .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%l==dat(LA+1:UA+1,:)%hme  ),mom,cnt,var)		
            write(name(im),'("%hme-mvs to",tr3,i4)') jloc
            weights(im)=whome
            im=im+1 
        end do   

        do g=1,2
            LA=MNA ; UA=MXAD
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. dat(LA:UA,:)%sexr==g .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
            write(name(im),'("wdif<0 | stay by sex ",i4)')  g
            weights(im)=wdifww  
            calcvar(im)=0
            im=im+1 
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. dat(LA:UA,:)%sexr==g .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
            write(name(im),'("wdif<0 | move by sex ",i4)')  g
            weights(im)=wdifww  
            calcvar(im)=0
            !momwhich(im,1)=16 ; momwhich(im,2)=-5 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
            im=im+1 
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. dat(LA:UA,:)%sexr==g .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%l/=dat(LA+1:UA+1,:)%hme  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),    d1*one(wdif(LA:UA,:)<0)  ,mom,cnt,var)		
            write(name(im),'("wdif<0 | hmemve=0 by sex ",i4)')  g
            weights(im)=wdifww  
            calcvar(im)=0
            im=im+1 
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. dat(LA:UA,:)%sexr==g .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%l==dat(LA+1:UA+1,:)%hme  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),    d1*one(wdif(LA:UA,:)<0)  ,mom,cnt,var)		
            write(name(im),'("wdif<0 | hmemve=1 by sex ",i4)') g  
            weights(im)=wdifww  
            calcvar(im)=0
            im=im+1 
        end do 



        do jloc=1,NL
            CALL condmom(im,( COSEXEDTYP(LA:UA,:).AND.(move(LA:UA,:)==1).AND.(dat(LA+1:UA+1,:)%l>=0).AND.(norelchg(LA:UA,:)==1) .and. dat(LA+1:UA+1,:)%l==jloc .AND. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),d1*one( (dat(LA+1:UA+1,:)%logwr-dat(LA:UA,:)%logwr<0) ),mom,cnt,var)
            WRITE(name(im),'("prop-of-moves-to , wdif move<0",I4)') jloc
            weights(im)=wdifww
            im=im+1 
        end do  
        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			

        do jloc=1,NL
            CALL condmom(im,( COSEXEDTYP(LA:UA,:).AND.(move(LA:UA,:)==1).AND.(norelchg(LA:UA,:)==1) .and. (dat(LA+1:UA+1,:)%l==jloc) .and. dat(LA:UA,:)%sexr==MALES .and. dat(LA:UA,:)%rel==0),d1*one( (dat(LA+1:UA+1,:)%hhr==0) ),mom,cnt,var)
            WRITE(name(im),'("u at dest men,sin",I4)') jloc
            weights(im)=0.
            im=im+1 
        end do  
        do jloc=1,NL
            CALL condmom(im,( COSEXEDTYP(LA:UA,:).AND.(move(LA:UA,:)==1).AND.(norelchg(LA:UA,:)==1) .and. (dat(LA+1:UA+1,:)%l==jloc) .and. dat(LA:UA,:)%sexr==MALES .and. dat(LA:UA,:)%rel==1),d1*one( (dat(LA+1:UA+1,:)%hhr==0) ),mom,cnt,var)
            WRITE(name(im),'("u at dest men,mar",I4)') jloc
            weights(im)=0.
            im=im+1 
        end do  
        do jloc=1,NL
            CALL condmom(im,( COSEXEDTYP(LA:UA,:).AND.(move(LA:UA,:)==1).AND.(norelchg(LA:UA,:)==1) .and. (dat(LA+1:UA+1,:)%l==jloc) .and. dat(LA:UA,:)%sexr==FEMALES .and. dat(LA:UA,:)%rel==0),d1*one( (dat(LA+1:UA+1,:)%hhr==0) ),mom,cnt,var)
            WRITE(name(im),'("u at dest fem,sin",I4)') jloc
            weights(im)=0.
            im=im+1 
        end do  
        do jloc=1,NL
            CALL condmom(im,( COSEXEDTYP(LA:UA,:).AND.(move(LA:UA,:)==1).AND.(norelchg(LA:UA,:)==1) .and. (dat(LA+1:UA+1,:)%l==jloc) .and. dat(LA:UA,:)%sexr==FEMALES .and. dat(LA:UA,:)%rel==1),d1*one( (dat(LA+1:UA+1,:)%hhr==0) ),mom,cnt,var)
            WRITE(name(im),'("u at dest fem,mar",I4)') jloc
            weights(im)=0.
            im=im+1 
        end do  

        IF (.NOT.SAVETIME) THEN
            !this was previously prop of moves to and from (conditionaoning on all the moves, what is the proportion of those moves that are to and from this location
            !but that one is very much influenced by the proportion of people in each location. so instead I use the below moments. Feb2023
            LA=MNA ; UA=MXAD
            do iloc=1,NL
                CALL condmom(im,( COSEXEDTYP(LA:UA,:) .AND.(dat(LA:UA,:)%l>=0).AND.(norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l==iloc)  .AND. (move(LA:UA,:)>=0) ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                WRITE(name(im),'("moves-from ",I4)') iloc
                weights(im)=wmove
                im=im+1 
            end do  
            LA=MNA ; UA=MXAD        
            do jloc=1,NL
                CALL condmom(im,( COSEXEDTYP(LA:UA,:) .AND.(dat(LA:UA,:)%l>=0).AND.(norelchg(LA:UA,:)==1.AND.dat(LA+1:UA+1,:)%l==jloc)  .AND. (move(LA:UA,:)>=0) ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                WRITE(name(im),'("moves-to  ",I4)') jloc
                weights(im)=wmove
                im=im+1 
            end do  

            LA=MNAD ; UA=MXAD
            do iloc=1,nl
                    call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA:UA,:)%l==iloc .AND. move(LA:UA,:)==1 ),   d1*one( moveadjacent(LA:UA,:)==1  ),mom,cnt,var)		
                    write(name(im),'("%adjacent mves from ",i4)') iloc
                    weights(im)=wmove
                    im=im+1       
                    call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA+1:UA+1,:)%l==iloc .AND. move(LA:UA,:)==1 ),   d1*one( moveadjacent(LA:UA,:)==1  ),mom,cnt,var)		
                    write(name(im),'("%adjacent mves to ",i4)') iloc
                    weights(im)=wmove
                    im=im+1       
            end do 

            !LA=MNAD ; UA=MXAD
            !do iloc=1,nl
            !    do jloc=1,nl         
            !    if (iloc/=jloc) then 
            !        call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA:UA,:)%l==iloc .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%l==jloc  ),mom,cnt,var)		
            !        write(name(im),'("%i to j",tr3,2i4)') iloc,jloc
            !        weights(im)=whome
            !        im=im+1       
            !        !call condmom(im,( coho(LA:UA,:) .AND. dat(LA:UA,:)%l==jloc .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%l==iloc  ),mom,cnt,var)		
            !        !write(name(im),'("%j to i",tr3,2i4)') jloc,iloc
            !        !weights(im)=whome
            !        !im=im+1 
            !    end if 
            !    end do     
            !end do 

            !headloc(ihead)=im
            !headstr(ihead)='everyone after the first move what proportion is return to home? (all)'
            !ihead=ihead+1
            LA=MNAD ; UA=MXAD     
            CALL condmom(im,( COSEXEDTYP(LA:UA,:).AND.  (moverank(LA:UA,:)>1).AND.(move(LA:UA,:)==1)  ),d1*one(homemove(LA:UA,:)==0),mom,cnt,var)
            WRITE(name(im),'("non-home ")') 
            weights(im)=whome
            im=im+1 
            CALL condmom(im,( COSEXEDTYP(LA:UA,:).AND.  (moverank(LA:UA,:)>1).AND.(move(LA:UA,:)==1)  ),d1*one(homemove(LA:UA,:)==1),mom,cnt,var)
            WRITE(name(im),'("home ")') 
            weights(im)=whome
            im=im+1 
            CALL condmom(im,(COSEXEDTYP(MNA:MXAD,:).AND.  (moverank(MNA:MXAD,:)>1).AND.(move(MNA:MXAD,:)==1)  ),d1*one(homemove(MNA:MXAD,:)==1),mom,cnt,var)
            WRITE(name(im),'("home     ")') 
            weights(im)=whome
            im=im+1 
        END IF !savetime
        !*****************************************************************************************

        !*************************************************************************************
        !BEGINNING OF MOMENTS THAT PERTAIN TO JUST MARRIAGE (UNIVERSE: EVERYONE OR MARRIED)
        if (.not.onlysingles) then
            !headloc(ihead)=im
            !headstr(ihead)='everyone mar and getmar and getdiv rates'
            !ihead=ihead+1
            !call condmom(im,( COSEXEDTYP(MNA:MXA,:) ), d1*one(dat(MNA:MXA,:)%rel==-1),mom,cnt,var)
            !write(name(im),'("margen=-1!!! ",tr10)')	
            !weights(im)=0.0_dp !; if (onlysingles) weights(im)=0.0_dp   !forget about that if statement usually. just putting there to get im to go up by 1 for the next loop
            !im=im+1

            call condmom(im,( COSEXEDTYP(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 ), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
            write(name(im),'("mar ",tr10)')	
            weights(im)=wrel !; if (onlysingles) weights(im)=0.0_dp   !forget about that if statement usually. just putting there to get im to go up by 1 for the next loop
            im=im+1
            call condmom(im,( COSEXEDTYP(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%edr==1), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
            write(name(im),'("mar ned ",tr10)')	
            weights(im)=wrel !; if (onlysingles) weights(im)=0.0_dp   !forget about that if statement usually. just putting there to get im to go up by 1 for the next loop
            im=im+1
            call condmom(im,( COSEXEDTYP(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%edr==2), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
            write(name(im),'("mar ed ",tr10)')	
            weights(im)=wrel !; if (onlysingles) weights(im)=0.0_dp   !forget about that if statement usually. just putting there to get im to go up by 1 for the next loop
            im=im+1

            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==0 .AND. dat(MNA+1:MXA,:)%rel>=0  ), d1*one(dat(MNA+1:MXA,:)%rel==1),mom,cnt,var)
            write(name(im),'("getmar all ")')
            weights(im)=wrel
            im=im+1
            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0  .AND. move(MNA:MXAD,:)==1), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move all ")')
            weights(im)=wrel
            im=im+1
            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0  .AND. move(MNA:MXAD,:)==0), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv stay all ")')
            weights(im)=wrel
            im=im+1        

            do ia=mna,mxad,10 !ahu030622  changed from maxai-1 to mxad
                call condmom(im,( COSEXEDTYP(ia,:) .AND. dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0  ), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
                write(name(im),'("getmarbyia ",i4)') ia
                weights(im)=wrel
                im=im+1
            end do      

            do ia=mna,mxad,10 !ahu030622 changed from maxai-1 to mxad    
                call condmom(im,( COSEXEDTYP(ia,:) .AND. dat(ia,:)%rel==1 .AND. dat(ia+1,:)%rel>=0 .AND. move(ia,:)==0), d1*one(dat(ia+1,:)%rel==0),mom,cnt,var)
                write(name(im),'("getdiv stay by ia ",i4)') ia
                weights(im)=wrel
                im=im+1
            end do 

            do ia=mna,mxad,10 !ahu030622 changed from maxai-1 to mxad    
                call condmom(im,( COSEXEDTYP(ia,:) .AND. dat(ia,:)%rel==1 .AND. dat(ia+1,:)%rel>=0 .AND. move(ia,:)==1), d1*one(dat(ia+1,:)%rel==0),mom,cnt,var)
                write(name(im),'("getdiv move by ia ",i4)') ia
                weights(im)=wrel
                im=im+1
            end do 

            ia=MNAD
            call condmom(im,( COSEXEDTYP(ia,:) .AND. dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0   .AND. dat(ia,:)%edr==1 ), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
            write(name(im),'("getmarbyia,ned ",i4)') ia
            weights(im)=wrel
            im=im+1
            !do ia=mna,mxad,10 !ahu030622  changed from maxai-1 to mxad
            !    call condmom(im,( COSEXEDTYP(ia,:) .AND. dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0   .AND. dat(ia,:)%edr==1 ), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
            !    write(name(im),'("getmarbyia,ned ",i4)') ia
            !    weights(im)=wrel
            !    im=im+1
            !end do      
            ia=agestart(College)-1
            call condmom(im,( COSEXEDTYP(ia,:) .AND. dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0   .AND. dat(ia,:)%edr==2 ), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
            write(name(im),'("getmarbyia, ed ",i4)') ia
            weights(im)=wrel
            im=im+1
            !do ia=mna,30,10 !ahu030622  changed from maxai-1 to mxad
            !    call condmom(im,( COSEXEDTYP(ia,:) .AND. dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0   .AND. dat(ia,:)%kidr==2 ), d1*one(dat(ia+1,:)!%rel==1),mom,cnt,var)
            !    write(name(im),'("getmarbyia, kid ",i4)') ia
            !    weights(im)=wrel
            !    im=im+1
            !end do

            do ia=mna,mxad,6 !ahu030622 changed from maxai-1 to mxad    
                call condmom(im,( COSEXEDTYP(ia,:) .AND. dat(ia,:)%rel==1 .AND. dat(ia+1,:)%rel>=0 .AND. move(ia,:)==1 .AND. dat(ia,:)%edr==1), d1*one(dat(ia+1,:)%rel==0),mom,cnt,var)
                write(name(im),'("getdiv move by ia,ned ",i4)') ia
                weights(im)=wrel
                im=im+1
            end do 
            do ia=mna,mxad,6 !ahu030622 changed from maxai-1 to mxad    
                call condmom(im,( COSEXEDTYP(ia,:) .AND. dat(ia,:)%rel==1 .AND. dat(ia+1,:)%rel>=0 .AND. move(ia,:)==1 .AND. dat(ia,:)%edr==2), d1*one(dat(ia+1,:)%rel==0),mom,cnt,var)
                write(name(im),'("getdiv move by ia,ed ",i4)') ia
                weights(im)=wrel
                im=im+1
            end do 
            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%edr==1  .and. dat(MNA:MXAD,:)%edsp==1 ), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv stay ned,ned ")') 
            weights(im)=wrel
            im=im+1
            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%edr==1  .and. dat(MNA:MXAD,:)%edsp==2 ), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv stay ned,ed ")') 
            weights(im)=wrel
            im=im+1

            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%edr==2  .and. dat(MNA:MXAD,:)%edsp==1 ), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv stay ed,ned ")') 
            weights(im)=wrel
            im=im+1

            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%edr==2  .and. dat(MNA:MXAD,:)%edsp==2 ), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv stay ed,ed ")') 
            weights(im)=wrel
            im=im+1


            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%edr==1  .and. dat(MNA:MXAD,:)%edsp==1 ), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move ned,ned ")') 
            weights(im)=wrel
            im=im+1

            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%edr==1  .and. dat(MNA:MXAD,:)%edsp==2 ), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move ned,ed ")') 
            weights(im)=wrel
            im=im+1

            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%edr==2  .and. dat(MNA:MXAD,:)%edsp==1 ), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move ed,ned ")') 
            weights(im)=wrel
            im=im+1

            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%edr==2  .and. dat(MNA:MXAD,:)%edsp==2 ), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move ed,ed ")') 
            weights(im)=wrel
            im=im+1

            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==1), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move all ")') 
            weights(im)=wrel
            im=im+1
            call condmom(im,( COSEXEDTYP(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==0), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv stay all ")') 
            weights(im)=wdivstay
            im=im+1            
            LA=18 ; UA=49
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==1), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move all   18:49 ")') 
            weights(im)=wrel
            im=im+1
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==0), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv stay all   18:49 ")') 
            weights(im)=wdivstay
            im=im+1            

            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            LA=18 ; UA=28 
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==1), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move all 18:28 ")') 
            weights(im)=wrel
            im=im+1
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==0), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv stay all  18:28 ")') 
            weights(im)=wdivstay
            im=im+1            
            LA=29 ; UA=34
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==1), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move all  29:34 ")') 
            weights(im)=wrel
            im=im+1
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==0), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv stay all  29:34 ")') 
            weights(im)=wdivstay
            im=im+1            
            LA=29 ; UA=40
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==1), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move all   29:40 ")') 
            weights(im)=wrel
            im=im+1
            call condmom(im,( COSEXEDTYP(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==0), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv stay all   29:40 ")') 
            weights(im)=wdivstay
            im=im+1            

            !if (j==1) then !only do kid by age for married. for singles, it looks weird in the data (.30, .8, .16,.19) 
            !    do ia=MNA,MXAD,8
            !        CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%kidr>=1 ),d1*one( dat(ia,:)%kidr==2 ),mom,cnt,var)
            !        WRITE(name(im),'("kid  by age",tr3,I4)') ia
            !        weights(im)=whour 
            !        im=im+1
            !    end do            
            !end if

            !CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%rel==1).and.(norelchg(MNA:MXAD,:)==1).AND.(kidtrans(MNA:MXAD,:)>=0)),d1*kidtrans(MNA:MXAD,:),mom,cnt,var)
            !WRITE(name(im),'("kidtrans-married ")') 
            !weights(im)=wkid
            !im=im+1
            !CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%rel==1).and.(norelchg(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%kidr>=1)),d1*dat(MNA:MXAD,:)%kidr,mom,cnt,var)
            !WRITE(name(im),'("kids-married     ")') 
            !weights(im)= wkid 
            !im=im+1
            !CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%rel==0).and.(norelchg(MNA:MXAD,:)==1).AND.(kidtrans(MNA:MXAD,:)>=0)),d1*kidtrans(MNA:MXAD,:),mom,cnt,var)
            !WRITE(name(im),'("kidtrans-single ")') 
            !weights(im)=wkid
            !im=im+1
            !CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%rel==0).and.(norelchg(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%kidr>=1)),d1*dat(MNA:MXAD,:)%kidr,mom,cnt,var)
            !WRITE(name(im),'("kids-single     ")') 
            !weights(im)= wkid 
            !im=im+1

        end if !not onlysingles
        !END OF MOMENTS THAT PERTAIN TO JUST MARRIAGE
        !*************************************************************************************



        
        !*******************************************************************************************************************************************************************************************************************************
        !START OF LOOP BY SEX AND REL
        do g=minsex,maxsex   
            do j=maxrelo,0,-1
                do typosto=0,ntypp,  DETAILEDOUTPUT
                !!!!!!if ( onlysingles ) then !  (.not.onlysingles).or.(onlysingles.and.j==0) ) then 
                !!!!!!    cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g )			                    
                !!!!!!else 
                !!!!!!    cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 .and. dat(MNAD:MXAD,:)%typ==typosto )			
                !!!!!!end if             
                if (j==1) wmove=wmovemar
                if (j==0) wmove=wmovesin
                if (typosto==0.and.(.not.onlysingles) ) then 
                    cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                else if (typosto==0.and.onlysingles) then 
                    cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g   )			
                else 
                    cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 .and. dat(MNAD:MXAD,:)%typ==typosto )			
                end if 

                headloc(ihead)=im
                if (g==1.and.j==0.and.typosto==0) headstr(ihead)='single men move rates all types'
                if (g==1.and.j==1.and.typosto==0) headstr(ihead)='married men move rates all types'
                if (g==2.and.j==0.and.typosto==0) headstr(ihead)='single fem move rates all types'
                if (g==2.and.j==1.and.typosto==0) headstr(ihead)='married fem move rates all types'
                if (g==1.and.j==0.and.typosto==1) headstr(ihead)='single men move rates type 1'
                if (g==1.and.j==1.and.typosto==1) headstr(ihead)='married men move rates type 1'
                if (g==2.and.j==0.and.typosto==1) headstr(ihead)='single fem move rates type 1'
                if (g==2.and.j==1.and.typosto==1) headstr(ihead)='married fem move rates type 1'
                if (g==1.and.j==0.and.typosto==2) headstr(ihead)='single men move rates type 2'
                if (g==1.and.j==1.and.typosto==2) headstr(ihead)='married men move rates type 2'
                if (g==2.and.j==0.and.typosto==2) headstr(ihead)='single fem move rates type 2'
                if (g==2.and.j==1.and.typosto==2) headstr(ihead)='married fem move rates type 2'
                if (g==1.and.j==0.and.typosto==3) headstr(ihead)='single men move rates type 3'
                if (g==1.and.j==1.and.typosto==3) headstr(ihead)='married men move rates type 3'
                if (g==2.and.j==0.and.typosto==3) headstr(ihead)='single fem move rates type 3'
                if (g==2.and.j==1.and.typosto==3) headstr(ihead)='married fem move rates type 3'
                if (g==1.and.j==0.and.typosto==4) headstr(ihead)='single men move rates type 4'
                if (g==1.and.j==1.and.typosto==4) headstr(ihead)='married men move rates type 4'
                if (g==2.and.j==0.and.typosto==4) headstr(ihead)='single fem move rates type 4'
                if (g==2.and.j==1.and.typosto==4) headstr(ihead)='married fem move rates type 4'

                ihead=ihead+1

                LA=MNAD ; UA=MXAD
                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==0  .AND. move(LA:UA,:)>=0 ),   d1* move(LA:UA,:) ,mom,cnt,var)	
                write(name(im),'("move | u MNAD ",tr5)')  
                weights(im)=wmove
                im=im+1
                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1  .AND. move(LA:UA,:)>=0 ),   d1* move(LA:UA,:) ,mom,cnt,var)	
                write(name(im),'("move | e MNAD ",tr5)')  
                weights(im)=wmove
                im=im+1  
                LA=MNA ; UA=MXAD
                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==0  .AND. move(LA:UA,:)>=0 ),   d1* move(LA:UA,:) ,mom,cnt,var)	
                write(name(im),'("move | u MNA ",tr5)')  
                weights(im)=wmove
                im=im+1
                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1  .AND. move(LA:UA,:)>=0 ),   d1* move(LA:UA,:) ,mom,cnt,var)	
                write(name(im),'("move | e MNA ",tr5)')  
                weights(im)=wmove
                im=im+1  

        

                !do ia=MNAD,20,1
                !    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==0 ),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                !    write(name(im),'("e | u stay ia ",tr3,I4)') ia  
                !    weights(im)=whour
                !    im=im+1 
                !    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==1 ),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                !    write(name(im),'("e | u move ia",tr3,I4)') ia  
                !    weights(im)=whour
                !    im=im+1 
                !    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==1 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==0 ),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                !    write(name(im),'("e | e stay ia",tr3,I4)')  ia
                !    weights(im)=whour
                !    im=im+1 
                !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%hhr==1 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==1 ),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                !    write(name(im),'("e | e move ia",tr3,I4)') ia  
                !    weights(im)=whour
                !    im=im+1 
                !end do
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | u stay",tr3)')  
                weights(im)=whour
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=0  ; momwhich(im,11)=6
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | u move",tr3)')  
                weights(im)=whour
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=6
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | e stay",tr3)')  
                weights(im)=whour
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=0  ; momwhich(im,11)=6
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:)  .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | e move",tr3)')  
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6
                weights(im)=whour
                im=im+1 


                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1  .AND. moveadjacent(MNA:MXAD,:)==0 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | u move NONADJA ",tr3)')  
                weights(im)=whour
                momwhich(im,1)=584 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=99 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=15  ; momwhich(im,11)=99
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1  .AND. moveadjacent(MNA:MXAD,:)==1 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | u move ADJA ",tr3)')  
                weights(im)=whour
                momwhich(im,1)=584 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=99 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=18  ; momwhich(im,11)=99
                im=im+1 

                call condmom(im,( cosexrel(MNA:MXAD,:)  .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. moveadjacent(MNA:MXAD,:)==0),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | e move NONADJA ",tr3)')  
                weights(im)=whour
                momwhich(im,1)=584 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=99 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=15  ; momwhich(im,11)=99
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:)  .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. moveadjacent(MNA:MXAD,:)==1),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | e move ADJA ",tr3)')  
                weights(im)=whour
                momwhich(im,1)=584 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=99 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=18  ; momwhich(im,11)=99
                im=im+1 


                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid                
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1  .AND. dat(MNA+1:MXA,:)%l/=dat(MNA+1:MXA,:)%hme ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | u move NONHOME ",tr3)')  
                weights(im)=whour
                momwhich(im,1)=585 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=99 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=15  ; momwhich(im,11)=99
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1  .AND. dat(MNA+1:MXA,:)%l==dat(MNA+1:MXA,:)%hme ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | u move HOME ",tr3)')  
                weights(im)=whour
                momwhich(im,1)=585 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=99 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=18  ; momwhich(im,11)=99
                im=im+1 

                call condmom(im,( cosexrel(MNA:MXAD,:)  .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. dat(MNA+1:MXA,:)%l/=dat(MNA+1:MXA,:)%hme),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | e move NONHOME ",tr3)')  
                weights(im)=whour
                momwhich(im,1)=585 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=99 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=15  ; momwhich(im,11)=99
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:)  .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. dat(MNA+1:MXA,:)%l==dat(MNA+1:MXA,:)%hme),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | e move HOME ",tr3)')  
                weights(im)=whour
                momwhich(im,1)=585 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=99 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=18  ; momwhich(im,11)=99
                im=im+1 


                LA=MNA ; UA=MXAD
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                write(name(im),'("wdif<0 | move ",2i4)') ier,iesp  
                weights(im)=wdifww  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA:UA,:)%edr==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                write(name(im),'("wdif<0 | move NED ")') 
                weights(im)=wdifww  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA:UA,:)%edr==2 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                write(name(im),'("wdif<0 | move ED ")') 
                weights(im)=wdifww  
                im=im+1 

                if (j==1) then 
                    LA=MNA ; UA=MXAD
                    do iesp=1,2
                        do ier=1,2    
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA:UA,:)%edr==ier .and. dat(LA:UA,:)%edsp==iesp .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                            write(name(im),'("wdif | move edr,edsp ",2i4)') ier,iesp  
                            weights(im)=wdifww  
                            im=im+1 
                        end do 
                    end do 

                    LA=MNA ; UA=MXAD
                    do iesp=0,1
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .AND. dat(LA:UA,:)%HHSP==iesp .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                            write(name(im),'("wdif | move SPOUSE EMP ",i4)') iesp  
                            weights(im)=wdifww  
                            im=im+1 
                    end do 




                    LA=MNA ; UA=MXAD
                    do iesp=1,2
                        do ier=1,2    
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA:UA,:)%edr==ier .and. dat(LA:UA,:)%edsp==iesp .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                            write(name(im),'("wdif<0 | move edr,edsp ",2i4)') ier,iesp  
                            weights(im)=wdifww  
                            im=im+1 
                        end do 
                    end do 

                    LA=MNA ; UA=MXAD
                    do iesp=0,1
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .AND. dat(LA:UA,:)%HHSP==iesp .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                            write(name(im),'("wdif<0 | move SPOUSE EMP ",i4)') iesp  
                            weights(im)=wdifww  
                            im=im+1 
                    end do 
                end if !j=1
    
                if (j==1 ) then !.and. momdisplay) then 
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND. moveadjacent(MNA:MXAD,:)==0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                    write(name(im),'("ee|ee move NONADJA")') 
                    weights(im)=whour
                    im=im+1 
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND. moveadjacent(MNA:MXAD,:)==1  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                    write(name(im),'("ee|ee move ADJA")') 
                    weights(im)=whour
                    im=im+1 


                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND. moveadjacent(MNA:MXAD,:)==0   ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                    write(name(im),'("eu|ee move NONADJA")') 
                    weights(im)=whour
                    im=im+1 
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND. moveadjacent(MNA:MXAD,:)==1   ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                    write(name(im),'("eu|ee move ADJA")') 
                    weights(im)=whour
                    im=im+1 


                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND. moveadjacent(MNA:MXAD,:)==0   ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                    write(name(im),'("ue|ee move NONADJA")') 
                    weights(im)=whour
                    im=im+1 
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND. moveadjacent(MNA:MXAD,:)==1   ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                    write(name(im),'("ue|ee move ADJA")') 
                    weights(im)=whour
                    im=im+1 
                    
                    
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND. moveadjacent(MNA:MXAD,:)==0   ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                    write(name(im),'("uu|ee move NONADJA")')
                    weights(im)=whour
                    im=im+1 

                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND. moveadjacent(MNA:MXAD,:)==1   ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                    write(name(im),'("uu|ee move ADJA")')
                    weights(im)=whour
                    im=im+1 



                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND.  dat(MNA+1:MXA,:)%l/=dat(MNA+1:MXA,:)%hme   ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                    write(name(im),'("ee|ee move NONHOME")') 
                    weights(im)=whour
                    im=im+1 
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND.  dat(MNA+1:MXA,:)%l==dat(MNA+1:MXA,:)%hme   ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                    write(name(im),'("ee|ee move HOME")') 
                    weights(im)=whour
                    im=im+1 


                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND.  dat(MNA+1:MXA,:)%l/=dat(MNA+1:MXA,:)%hme    ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                    write(name(im),'("eu|ee move NONHOME")') 
                    weights(im)=whour
                    im=im+1 
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND.  dat(MNA+1:MXA,:)%l==dat(MNA+1:MXA,:)%hme    ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                    write(name(im),'("eu|ee move HOME")') 
                    weights(im)=whour
                    im=im+1 


                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND.  dat(MNA+1:MXA,:)%l/=dat(MNA+1:MXA,:)%hme    ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                    write(name(im),'("ue|ee move NONHOME")') 
                    weights(im)=whour
                    im=im+1 
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND.  dat(MNA+1:MXA,:)%l==dat(MNA+1:MXA,:)%hme    ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                    write(name(im),'("ue|ee move HOME")') 
                    weights(im)=whour
                    im=im+1 
                    
                    
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND. dat(MNA+1:MXA,:)%l/=dat(MNA+1:MXA,:)%hme   ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                    write(name(im),'("uu|ee move NONHOME")')
                    weights(im)=whour
                    im=im+1 

                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1 .AND.  dat(MNA+1:MXA,:)%l==dat(MNA+1:MXA,:)%hme    ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                    write(name(im),'("uu|ee move HOME")')
                    weights(im)=whour
                    im=im+1 

                end if !REL AND MOMDISPLAY

                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%kidr==1),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                !write(name(im),'("e | e stay no kid",tr3)')  
                !weights(im)=0.0_dp
                !im=im+1 
                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%kidr==2),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                !write(name(im),'("e | e stay   kid",tr3)')  
                !weights(im)=0.0_dp
                !im=im+1 

                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%kidr==1),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | u stay no kid",tr3)')  
                weights(im)=whour
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=0  ; momwhich(im,11)=1
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%kidr==1),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | u move no kid",tr3)')  
                weights(im)=whour
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=1
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%kidr==1),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | e stay no kid",tr3)')  
                weights(im)=whour
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=0  ; momwhich(im,11)=1
                im=im+1 

                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%kidr==1),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | e move no kid",tr3)')  
                weights(im)=whour
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=1
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%kidr==2),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | u stay   kid",tr3)')  
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=0  ; momwhich(im,11)=2
                weights(im)=whour
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%kidr==2),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | u move   kid",tr3)')  
                weights(im)=whour
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=2
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%kidr==2),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | e stay   kid",tr3)')  
                weights(im)=whour
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=0  ; momwhich(im,11)=2
                im=im+1 
                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%kidr==2),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e | e move   kid",tr3)')  
                weights(im)=whour
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=2
                im=im+1 



                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do edo=1,2
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%edr==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                    write(name(im),'("e | u stay by ed",tr3,i4)')  edo
                        momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=edo ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=0  ; momwhich(im,11)=6
                    weights(im)=whour
                    im=im+1
                end do 
                do edo=1,2
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%edr==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                    write(name(im),'("e | u move by ed",tr3,i4)')  edo
                    weights(im)=whour
                        momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=edo ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=6
                    im=im+1 
                end do 
                do edo=1,2
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0  .and. dat(MNA:MXAD,:)%edr==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                    write(name(im),'("e | e stay by ed",tr3,i4)') edo  
                    weights(im)=whour
                        momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=edo ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=0  ; momwhich(im,11)=6
                    im=im+1 
                end do 
                do edo=1,2
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%edr==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                    write(name(im),'("e | e move by ed",tr3,i4)') edo  
                    weights(im)=whour
                    momwhich(im,1)=500 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=edo ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6
                    im=im+1 
                end do 




                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                if (j==1) then 
                    LA=MNA ; UA=MXAD
                    do iesp=1,2
                        do ier=1,2
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | u move jointed typ ",2i4)')  ier,iesp
                            weights(im)=0.0_dp
                            momwhich(im,1)=518 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=6
                            im=im+1 
                        end do 
                    end do

                    LA=MNA ; UA=MXAD
                    do iesp=1,2
                        do ier=1,2                            
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | e move jointed typ ",2i4)')  ier,iesp
                            weights(im)=0.0_dp
                            momwhich(im,1)=518 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6
                            im=im+1 
                        end do 
                    end do
                end if





                call condmom(im,( cosexrel(MNAD:MXAD,:) .AND. dat(MNAD:MXAD,:)%edr==1  .AND. move(MNAD:MXAD,:)>=0 ),   d1* move(MNAD:MXAD,:) ,mom,cnt,var)	
                write(name(im),'("move | ned ",tr1)')  
                weights(im)=wmove
                momwhich(im,1)=270 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=1 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                im=im+1
                call condmom(im,( cosexrel(MNAD:MXAD,:) .AND. dat(MNAD:MXAD,:)%edr==2  .AND. move(MNAD:MXAD,:)>=0 ),   d1* move(MNAD:MXAD,:) ,mom,cnt,var)	
                write(name(im),'("move | ed ",tr1)')  
                weights(im)=wmove
                momwhich(im,1)=270 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=2 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                im=im+1
                call condmom(im,( cosexrel(MNAD:MXAD,:)  .AND. move(MNAD:MXAD,:)>=0 ),   d1* move(MNAD:MXAD,:) ,mom,cnt,var)	
                write(name(im),'("move  ",tr1)')  
                weights(im)=wmove
                momwhich(im,1)=270 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                im=im+1


                call condmom(im,( cosexrel(MNAD:MXAD,:) .AND. dat(MNAD:MXAD,:)%kidr==1  .AND. move(MNAD:MXAD,:)>=0 ),   d1* move(MNAD:MXAD,:) ,mom,cnt,var)	
                write(name(im),'("move | nokid ",tr1)')  
                weights(im)=wmove
                im=im+1

                call condmom(im,( cosexrel(MNAD:MXAD,:) .AND. dat(MNAD:MXAD,:)%kidr==2  .AND. move(MNAD:MXAD,:)>=0 ),   d1* move(MNAD:MXAD,:) ,mom,cnt,var)	
                write(name(im),'("move | kid ",tr1)')  
                weights(im)=wmove
                im=im+1


                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                !LA=MNA ; UA=MXAD
                !call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0 .AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme  ),   d1*one( dat(LA+1:UA+1,:)%l/=dat(LA+1:UA+1,:)%hme .and. dat(LA+1:UA+1,:)%l==dat(LA:UA,:)%l  ),mom,cnt,var)		
                !write(name(im),'("AT home - stay    ",2i4)') LA,UA !added this ahu 121718
                !weights(im)=whome
                !im=im+1 
                !call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0 .AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme  ),   d1*one( dat(LA+1:UA+1,:)%l==dat(LA+1:UA+1,:)%hme .and. dat(LA+1:UA+1,:)%l/=dat(LA:UA,:)%l ),mom,cnt,var)		
                !write(name(im),'("AT home - move ",2i4)') LA,UA  !added this ahu 121718
                !weights(im)=whome
                !im=im+1 
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                LA=MNA ; UA=MXAD
                call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0 .AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme  ),   d1*one( dat(LA+1:UA+1,:)%l/=dat(LA+1:UA+1,:)%hme .and. dat(LA+1:UA+1,:)%l==dat(LA:UA,:)%l  ),mom,cnt,var)		
                write(name(im),'("NOT home - stay    ",2i4)') LA,UA !added this ahu 121718
                weights(im)=whome
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0 .AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme  ),   d1*one( dat(LA+1:UA+1,:)%l==dat(LA+1:UA+1,:)%hme .and. dat(LA+1:UA+1,:)%l/=dat(LA:UA,:)%l ),mom,cnt,var)		
                write(name(im),'("NOT home - go home ",2i4)') LA,UA  !added this ahu 121718
                weights(im)=whome
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0 .AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme  ),   d1*one( dat(LA+1:UA+1,:)%l/=dat(LA+1:UA+1,:)%hme .and. dat(LA+1:UA+1,:)%l/=dat(LA:UA,:)%l  ),mom,cnt,var)		
                write(name(im),'("NOT home - go else ",2i4)') LA,UA !added this ahu 121718                
                weights(im)=whome
                im=im+1 


                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. move(MNA:MXAD,:)==1 ),   d1*one( dat(MNA+1:MXA,:)%l==dat(MNA+1:MXA,:)%hme  ),mom,cnt,var)		
                write(name(im),'("move home ",tr3)')  !added this ahu 121718
                weights(im)=whome
                im=im+1 




                !****************************************************************************************************************************************************************************************************************
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                LA=MNAD ; UA=MXAD
                CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme  .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                WRITE(name(im),'("move 17-49 - at home ALL ED ")') 
                weights(im)=wprop
                momwhich(im,1)=280 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                im=im+1 
                LA=MNAD ; UA=25
                CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme  .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                WRITE(name(im),'("move 17-25 - at home ALL ED ")') 
                weights(im)=wprop
                momwhich(im,1)=280 ; momwhich(im,2)=1725 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                im=im+1 
                LA=26 ; UA=34
                CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme  .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                WRITE(name(im),'("move 26-34 - at home ALL ED ")') 
                momwhich(im,1)=280 ; momwhich(im,2)=2634 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                weights(im)=wprop
                im=im+1
                do edo=1,2
                    LA=MNAD ; UA=MXAD
                    CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme .and.  dat(LA:UA,:)%edr==edo .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                    WRITE(name(im),'("move 17-49 - at home by ED ",i4)') edo
                    weights(im)=wprop
                    momwhich(im,1)=280 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=edo ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    im=im+1 
                    LA=MNAD ; UA=25
                    CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme   .and.  dat(LA:UA,:)%edr==edo .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                    WRITE(name(im),'("move 17-25 - at home by ED ",i4)') edo
                    weights(im)=wprop
                    momwhich(im,1)=280 ; momwhich(im,2)=1725 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=edo ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    im=im+1 
                    LA=26 ; UA=34
                    CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme   .and.  dat(LA:UA,:)%edr==edo .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                    WRITE(name(im),'("move 26-34 - at home by ED ",i4)') edo
                    momwhich(im,1)=280 ; momwhich(im,2)=2634 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=edo ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=wprop
                    im=im+1
                end do 
                !****************************************************************************************************************************************************************************************************************
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                LA=MNAD ; UA=MXAD
                CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme  .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                WRITE(name(im),'("move 17-49 - not home ALL ED ")') 
                weights(im)=wprop
                momwhich(im,1)=281 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                im=im+1 
                LA=MNAD ; UA=25
                CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme  .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                WRITE(name(im),'("move 17-25 - not home ALL ED ")') 
                weights(im)=wprop
                momwhich(im,1)=281 ; momwhich(im,2)=1725 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                im=im+1 
                LA=26 ; UA=34
                CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme  .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                WRITE(name(im),'("move 26-34 - not home ALL ED ")') 
                momwhich(im,1)=281 ; momwhich(im,2)=2634 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                weights(im)=wprop
                im=im+1
                do edo=1,2
                    LA=MNAD ; UA=MXAD
                    CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme .and.  dat(LA:UA,:)%edr==edo .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                    WRITE(name(im),'("move 17-49 - not home by ED ",i4)') edo
                    weights(im)=wprop
                    momwhich(im,1)=281 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=edo ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    im=im+1 
                    LA=MNAD ; UA=25
                    CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme   .and.  dat(LA:UA,:)%edr==edo .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                    WRITE(name(im),'("move 17-25 - not home by ED ",i4)') edo
                    weights(im)=wprop
                    momwhich(im,1)=281 ; momwhich(im,2)=1725 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=edo ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    im=im+1 
                    LA=26 ; UA=34
                    CALL condmom(im,( COSEXREL(LA:UA,:) .AND.dat(LA:UA,:)%l>=0.AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme   .and.  dat(LA:UA,:)%edr==edo .AND. move(LA:UA,:)>=0 ),d1*one(move(LA:UA,:)==1),mom,cnt,var)
                    WRITE(name(im),'("move 26-34 - not home by ED ",i4)') edo
                    momwhich(im,1)=281 ; momwhich(im,2)=2634 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=edo ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=wprop
                    im=im+1
                end do 
                !****************************************************************************************************************************************************************************************************************


                !ahu jan19 010219: so I am still writing the age 17 (mad) mar rate. age 17 mar rate does not exit in the data since we don't record age 17 rel there. 
                !so the weight on this moment is 0 but then it's wmovebyrel
                !do ia=MNAD,21,1
                !    CALL condmom(im,( cosexrel(ia,:) .AND. move(ia,:)>=0 ),d1*move(ia,:),mom,cnt,var)
                !    WRITE(name(im),'("move by age",tr3,I4)') ia
                !    weights(im)=wmove
                !    im=im+1
                !end do
                !do ia=22,50,4
                !    CALL condmom(im,( cosexrel(ia,:) .AND. move(ia,:)>=0 ),d1*move(ia,:),mom,cnt,var)
                !    WRITE(name(im),'("move by age",tr3,I4)') ia
                !    weights(im)=wmove
                !    im=im+1
                !end do     




                    !IF (.NOT.SAVETIME) THEN
                    !        call condmom(im,( cosexrel(MNA:MXAD,:)  .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .AND. dat(MNA:MXAD,:)%jobswitch>=0),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA:MXAD,:)%jobswitch==0 ),mom,cnt,var)		
                    !        write(name(im),'("e | e stay,jobswitch=0",tr3)')  
                    !        weights(im)=whour
                    !        im=im+1 
                    !        do ia=MNAD,22,1
                    !            call condmom(im,( cosexrel(ia,:)   .AND. dee(ia,:)==1  .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==0 .AND. dat(ia,:)%jobswitch>=0 ),   d1*one( dat(ia,:)%jobswitch==1 ),mom,cnt,var)		
                    !            write(name(im),'("e | e stay,jobswitch=1",i4)') ia 
                    !            weights(im)=whour
                    !            im=im+1 
                    !        end do     
                    !        do ia=23,MXAD,6
                    !            call condmom(im,( cosexrel(ia,:)   .AND. dee(ia,:)==1  .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==0 .AND. dat(ia,:)%jobswitch>=0 ),   d1*one( dat(ia,:)%jobswitch==1 ),mom,cnt,var)		
                    !            write(name(im),'("e | e stay,jobswitch=1",i4)') ia 
                    !            weights(im)=whour
                    !            im=im+1 
                    !        end do     
                    !        do ia=MNAD,22,1
                    !            call condmom(im,( cosexrel(ia,:)  .AND. dee(ia,:)==1  .AND. move(ia,:)==0 .and. dat(ia,:)%jobswitch==0 .and. dat(ia+1,:)%logwr>=0 .AND. dat(ia,:)%logwr>=0 ),   d1*wdif(ia,:) ,mom,cnt,var)		
                    !            write(name(im),'("wdif | stay,jobwitch==0 ",I4)') ia  
                    !            weights(im)=wdifww  
                    !            calcvar(im)=0 !1
                    !            im=im+1 
                    !        end do     
                    !        do ia=23,MXAD,6
                    !            call condmom(im,( cosexrel(ia,:)  .AND. dee(ia,:)==1  .AND. move(ia,:)==0 .and. dat(ia,:)%jobswitch==0  .and. dat(ia+1,:)%logwr>=0 .AND. dat(ia,:)%logwr>=0 ),   d1*wdif(ia,:) ,mom,cnt,var)		
                    !            write(name(im),'("wdif | stay,jobwitch==0 ",I4)') ia  
                    !            weights(im)=wdifww  
                    !            calcvar(im)=0 !1
                    !            im=im+1 
                    !        end do     
                    !        do ia=MNAD,22,1
                    !            call condmom(im,( cosexrel(ia,:)   .AND. dee(ia,:)==1  .AND. move(ia,:)==0 .and. dat(ia,:)%jobswitch==1  .and. dat(ia+1,:)%logwr>=0 .AND. dat(ia,:)%logwr>=0 ),   d1*wdif(ia,:), mom,cnt,var)		
                    !            write(name(im),'("wdif | stay,jobwitch==1 ",I4)') ia  
                    !            weights(im)=wswitch  
                    !            calcvar(im)=0 !1
                    !            im=im+1 
                    !        end do     
                    !        do ia=23,MXAD,6
                    !            call condmom(im,( cosexrel(ia,:)   .AND. dee(ia,:)==1  .AND. move(ia,:)==0 .and. dat(ia,:)%jobswitch==1  .and. dat(ia+1,:)%logwr>=0 .AND. dat(ia,:)%logwr>=0  ),   d1*wdif(ia,:), mom,cnt,var)		
                    !            write(name(im),'("wdif | stay,jobwitch==1 ",I4)') ia  
                    !            weights(im)=wswitch  
                    !            calcvar(im)=0 !1
                    !            im=im+1 
                    !        end do     
                    !        !do ia=MNAD,MXAD,2
                    !        !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%typ==i .AND. dee(ia,:)==1  .AND. move(ia,:)==0 .and. dat(ia,:)%jobswitch==0 ),   d1*one( dat(ia+1,:)%logwr-dat(ia,:)%logwr<0 ),mom,cnt,var)		
                    !        !    write(name(im),'("wdif<0 | stay,jobwitch==0 ",I4)') ia  
                    !        !    weights(im)=wdifww  
                    !        !    calcvar(im)=0 !1
                    !        !    im=im+1 
                    !        !end do     
    !
    !
                    !        !do ia=agestart(NOCOLLEGE)-1,27,5
                    !        !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr==1 .AND. dat(ia+1,:)%logwr>=0  .AND. dat(ia,:)%edr==1),   d1*dat(ia+1,:)%logwr ,mom,cnt,var)		
                    !        !    write(name(im),'("w|u by age ned ",i4)') ia
                    !        !    weights(im)=whour
                    !        !    im=im+1 
                    !        !end do 
                    !        !do ia=agestart(COLLEGE)-1,27,5
                    !        !    call condmom(im,( cosexrel(ia,:)   .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr==1 .AND. dat(ia+1,:)%logwr>=0  .AND. dat(ia,:)%edr==2),   d1*dat(ia+1,:)%logwr ,mom,cnt,var)		
                    !        !    write(name(im),'("w|u by age ed ",i4)') ia
                    !        !    weights(im)=whour
                    !        !    im=im+1 
                    !        !end do 
                    !    END IF !SAVETIME



                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%kidr==2  .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
                !write(name(im),'("move | kid ",tr3)')  
                !weights(im)=wmove
                !im=im+1


                if (j==1) then


                    do iesp=1,2 !educ spouse
                        do ier=1,2 !educ r
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND.norelchg(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)>=0  .and. dat(MNA:MXAD,:)%edr==ier  .and. dat(MNA:MXAD,:)%edsp==iesp .AND. move(MNA:MXAD,:)==1),   d1*one( moveadjacent(MNA:MXAD,:)==1  ),mom,cnt,var)		
                                write(name(im),'("move adjacent (cond.) by jointed",2i4)') ier,iesp 
                                weights(im)=wmove
                                im=im+1 
                        end do 
                    end do 
             

                    do iesp=1,2 !educ spouse
                        do ier=1,2 !educ r
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND.norelchg(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)>=0  .and. dat(MNA:MXAD,:)%edr==ier  .and. dat(MNA:MXAD,:)%edsp==iesp ),   d1*one( moveadjacent(MNA:MXAD,:)==1  ),mom,cnt,var)		
                                write(name(im),'("move adjacent (uncond) by jointed",2i4)') ier,iesp 
                                weights(im)=wmove
                                im=im+1 
                        end do 
                    end do 


                    !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                    do edo=1,2
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%edr==1 .and. dat(MNA:MXAD,:)%edsp==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | u stay NED by sped",tr3,i4)')  edo
                        weights(im)=whour
                        momwhich(im,1)=350 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=1 ; momwhich(im,7)=edo ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        im=im+1
                    end do 
                    do edo=1,2
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%edr==2 .and. dat(MNA:MXAD,:)%edsp==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | u stay  ED by sped",tr3,i4)')  edo
                        weights(im)=whour
                        momwhich(im,1)=350 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=2 ; momwhich(im,7)=edo ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        im=im+1
                    end do 
                    do edo=1,2
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%edr==1 .and. dat(MNA:MXAD,:)%edsp==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | u move NED by sped",tr3,i4)')  edo
                        weights(im)=whour
                        momwhich(im,1)=350 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=1 ; momwhich(im,7)=edo ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        im=im+1
                    end do 
                    do edo=1,2
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%edr==2 .and. dat(MNA:MXAD,:)%edsp==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | u move  ED by sped",tr3,i4)')  edo
                        weights(im)=whour
                        momwhich(im,1)=350 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=2 ; momwhich(im,7)=edo ;  momwhich(im,8)=99 ; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        im=im+1
                    end do 
                    do edo=1,2
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%edr==1 .and. dat(MNA:MXAD,:)%edsp==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | e stay NED by sped",tr3,i4)')  edo
                        weights(im)=whour
                        momwhich(im,1)=350 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=1 ; momwhich(im,7)=edo ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        im=im+1
                    end do 
                    do edo=1,2
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .and. dat(MNA:MXAD,:)%edr==2 .and. dat(MNA:MXAD,:)%edsp==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | e stay  ED by sped",tr3,i4)')  edo
                        weights(im)=whour
                        momwhich(im,1)=350 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=2 ; momwhich(im,7)=edo ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        im=im+1
                    end do 
                    do edo=1,2
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%edr==1 .and. dat(MNA:MXAD,:)%edsp==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | e move NED by sped",tr3,i4)')  edo
                        weights(im)=whour
                        momwhich(im,1)=350 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=1 ; momwhich(im,7)=edo ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        im=im+1
                    end do 
                    do edo=1,2
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .and. dat(MNA:MXAD,:)%edr==2 .and. dat(MNA:MXAD,:)%edsp==edo),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | e move  ED by sped",tr3,i4)')  edo
                        weights(im)=whour
                        momwhich(im,1)=350 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=2 ; momwhich(im,7)=edo ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        im=im+1
                    end do 
                end if !j=1

                    !call condmom(im,( COSEXREL(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==1), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
                    !write(name(im),'("getdiv move by sex ")') 
                    !weights(im)=wrel
                    !im=im+1

                    !call condmom(im,( COSEXREL(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1 .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==0), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,!var)
                    !write(name(im),'("getdiv stay by sex ")') 
                    !weights(im)=wrel
                    !im=im+1            
                    !Moved these here conditioning on sex because note that the respondents can be both male as well as female! 
                    !headloc(ihead)=im; headstr(ihead)='joint emp transition rates ';ihead=ihead+1
                    IF (j==1) then !only for married
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==0 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                        write(name(im),'("ee|ee stay")')  
                        weights(im)=whour  
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=0  ; momwhich(im,11)=6  
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                        write(name(im),'("eu|ee stay")') 
                        weights(im)=whour 
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=2 ; momwhich(im,10)=0  ; momwhich(im,11)=6  
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==0 ),   d1*one( dat(MNA+1:MXA,:)%hhr==0..AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                        write(name(im),'("ue|ee stay")') 
                        weights(im)=whour 
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=3 ; momwhich(im,10)=0  ; momwhich(im,11)=6  
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==0 ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                        write(name(im),'("uu|ee stay")')  
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=4 ; momwhich(im,10)=0  ; momwhich(im,11)=6  
                        im=im+1 

                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1   ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                        write(name(im),'("ee|ee move")') 
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1    ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                        write(name(im),'("eu|ee move")') 
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=2 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1    ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                        write(name(im),'("ue|ee move")') 
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=3 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1    ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                        write(name(im),'("uu|ee move")')
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=4 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                        im=im+1 

        
                        !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid

                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%kidr==1 ) ,   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                        write(name(im),'("ee|ee move NO KID")')
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=1  
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%kidr==2 ) ,   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                        write(name(im),'("ee|ee move    KID")')
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=2
                        im=im+1 


                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%kidr==1  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                        write(name(im),'("eu|ee move NO KID")') 
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=2 ; momwhich(im,10)=1  ; momwhich(im,11)=1
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%kidr==2  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                        write(name(im),'("eu|ee move    KID")') 
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=2 ; momwhich(im,10)=1  ; momwhich(im,11)=2
                        im=im+1 
                        
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%kidr==2  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                        write(name(im),'("ue|ee move  NO KID")') 
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=3 ; momwhich(im,10)=1  ; momwhich(im,11)=1
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%kidr==2  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                        write(name(im),'("ue|ee move    KID")') 
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=3 ; momwhich(im,10)=1  ; momwhich(im,11)=2
                        im=im+1 

                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%kidr==2  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                        write(name(im),'("uu|ee move  NO KID")') 
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=4 ; momwhich(im,10)=1  ; momwhich(im,11)=1
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%kidr==2  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                        write(name(im),'("uu|ee move    KID")') 
                        weights(im)=whour
                        momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=4 ; momwhich(im,10)=1  ; momwhich(im,11)=2
                        im=im+1 
                        do iesp=1,2 !educ spouse
                        do ier=1,2 !educ r
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%edr==ier  .and. dat(MNA:MXAD,:)%edsp==iesp),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                                write(name(im),'("ee|ee move",2i4)') ier,iesp 
                                weights(im)=whour
                                momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp ;  momwhich(im,8)=99 ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=99  
                                im=im+1 
                            end do 
                        end do 

                        do iesp=1,2 !educ spouse
                        do ier=1,2 !educ r
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%edr==ier  .and. dat(MNA:MXAD,:)%edsp==iesp),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                                write(name(im),'("eu|ee move",2i4)') ier,iesp  
                                weights(im)=whour
                                momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp ;  momwhich(im,8)=99 ; momwhich(im,9)=2 ; momwhich(im,10)=1  ; momwhich(im,11)=99  
                                im=im+1 
                            end do 
                        end do 

                        do iesp=1,2 !educ spouse
                        do ier=1,2 !educ r
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%edr==ier  .and. dat(MNA:MXAD,:)%edsp==iesp),   d1*one( dat(MNA+1:MXA,:)%hhr==0..AND. dat(MNA+1:MXA,:)%hhsp==1  ),mom,cnt,var)		
                                write(name(im),'("ue|ee move",2i4)') ier,iesp  
                                weights(im)=whour
                                momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp ;  momwhich(im,8)=99 ; momwhich(im,9)=3 ; momwhich(im,10)=1  ; momwhich(im,11)=99  
                                im=im+1 
                            end do 
                        end do 


                        do iesp=1,2 !educ spouse
                            do ier=1,2 !educ r
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. move(MNA:MXAD,:)==1  .and. dat(MNA:MXAD,:)%edr==ier  .and. dat(MNA:MXAD,:)%edsp==iesp),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhsp==0  ),mom,cnt,var)		
                                write(name(im),'("uu|ee move",2i4)') ier,iesp 
                                weights(im)=whour
                                momwhich(im,1)=400 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp ;  momwhich(im,8)=99 ; momwhich(im,9)=4 ; momwhich(im,10)=1  ; momwhich(im,11)=99  
                                im=im+1 
                            end do !iesp
                        end do !ier

                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0 ), d1*  (dat(MNA:MXAD,:)%logwr)  ,mom,cnt,var)		
                        write(name(im),'("lwager")')  
                        weights(im)=wwage
                        calcorr(im)=1
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0 ), d1*  (dat(MNA:MXAD,:)%logwr**2)  ,mom,cnt,var)		
                        write(name(im),'("lwagersq")')  
                        weights(im)=wwage
                        calcorr(im)=5
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0 ), d1*  (dat(MNA:MXAD,:)%logwsp)  ,mom,cnt,var)		
                        write(name(im),'("lwagesp")')  
                        weights(im)=wwage
                        calcorr(im)=5
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0 ), d1*  (dat(MNA:MXAD,:)%logwsp**2)  ,mom,cnt,var)		
                        write(name(im),'("lwagespsq")')  
                        weights(im)=wwage
                        calcorr(im)=5
                        im=im+1 
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0 ), d1*  (dat(MNA:MXAD,:)%logwr*dat(MNA:MXAD,:)%logwsp)  ,mom,cnt,var)		
                        write(name(im),'("LW*LWSP")')  
                        weights(im)=wcorr
                        calcorr(im)=5
                        im=im+1
                        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0 ), d1*  (dat(MNA:MXAD,:)%logwr*dat(MNA:MXAD,:)%logwsp)  ,mom,cnt,var)		
                        write(name(im),'("LW*LWSP")')  
                        weights(im)=wcorr
                        calcorr(im)=5
                        im=im+1


                        do iesp=0,1 !emp spouse
                            do ier=0,1 !emp r
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND.dat(MNA:MXAD,:)%l>=0.AND.norelchg(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)>=0  .and. dat(MNA:MXAD,:)%hhr==ier  .and. dat(MNA:MXAD,:)%hhsp==iesp),   d1*one(move(MNA:MXAD,:)==1),mom,cnt,var)		
                                write(name(im),'("move by jointemp",2i4)') ier,iesp 
                                weights(im)=wmove
                                im=im+1 
                            end do 
                        end do                     
                        do iesp=1,2 !educ spouse
                            do ier=1,2 !educ r
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND.dat(MNA:MXAD,:)%l>=0.AND.norelchg(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)>=0  .and. dat(MNA:MXAD,:)%edr==ier  .and. dat(MNA:MXAD,:)%edsp==iesp),   d1*one(move(MNA:MXAD,:)==1),mom,cnt,var)		
                                write(name(im),'("move by jointed",2i4)') ier,iesp 
                                weights(im)=wmove
                                im=im+1 
                            end do 
                        end do 
                        do iesp=1,2 !educ spouse
                            do ier=1,2 !educ r
                                call condmom(im,( cosexrel(MNA:MXAD,:) ), d1*  one(dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp)  ,mom,cnt,var)		
                                write(name(im),'("COUPLES JOINTED ",2i4)') ier,iesp
                                weights(im)=wrel*10.0_dp
                                im=im+1 
                            end do 
                        end do 


                        do iesp=1,2 !educ spouse
                        do ier=1,2 !educ r
                                i=1
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0 .AND. dat(MNA:MXAD,:)%rellen==i .AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwr)  ,mom,cnt,var)		
                                write(name(im),'("lwager   ",tr9)') 
                                weights(im)=wwagejoint
                                calcorr(im)=1
                                im=im+1 
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0 .AND. dat(MNA:MXAD,:)%rellen==i .AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwr**2)  ,mom,cnt,var)		
                                write(name(im),'("lwagersq ",tr9)') 
                                weights(im)=wwagejoint
                                calcorr(im)=5
                                im=im+1 
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0  .AND. dat(MNA:MXAD,:)%rellen==i.AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwsp)  ,mom,cnt,var)		
                                write(name(im),'("lwagesp  ",tr9)')
                                weights(im)=wwagejoint
                                calcorr(im)=5
                                im=im+1 
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0  .AND. dat(MNA:MXAD,:)%rellen==i.AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwsp**2)  ,mom,cnt,var)		
                                write(name(im),'("lwagespsq",tr9)') 
                                weights(im)=wwagejoint
                                calcorr(im)=5
                                im=im+1 
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0  .AND. dat(MNA:MXAD,:)%rellen==i .AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwr*dat(MNA:MXAD,:)%logwsp)  ,mom,cnt,var)		
                                write(name(im),'("CORRELATION LW*LWSP ")')
                                weights(im)=wcorr
                                calcorr(im)=5
                                im=im+1
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0  .AND. dat(MNA:MXAD,:)%rellen==i .AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwr*dat(MNA:MXAD,:)%logwsp)  ,mom,cnt,var)		
                                write(name(im),'("CORRELATION LW*LWSP ",3i4)') ier,iesp,i
                                weights(im)=wcorr
                                calcorr(im)=5
                                im=im+1
                                do ddd=1,ntypp
                                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*one(dat(MNA:MXAD,:)%typ==ddd),mom,cnt,var)
                                write(name(im),'("TYP | MAR JOINTED   ",3i4)') ier,iesp,ddd
                                weights(im)=0.0_dp
                                im=im+1
                                end do 
                            end do !ier
                        end do !iesp
                        !IF (.NOT.SAVETIME) THEN
                        !    do iesp=1,2 !educ spouse
                        !    do ier=1,2 !educ r
                        !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0  .AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwr)  ,mom,cnt,var)		
                        !        write(name(im),'("lwager ALL jed ",tr9,2i4)') ier,iesp
                        !        weights(im)=wwage
                        !        calcorr(im)=1
                        !        im=im+1 
                        !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0  .AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwr**2)  ,mom,cnt,var)		
                        !        write(name(im),'("lwager ALL jed ",tr9,2i4)') ier,iesp
                        !        weights(im)=wwage
                        !        calcorr(im)=5
                        !        im=im+1 
                        !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0 .AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwsp)  ,mom,cnt,var)		
                        !        write(name(im),'("lwagesp ALL jed ",tr9,2i4)') ier,iesp
                        !        weights(im)=wwage
                        !        calcorr(im)=5
                        !        im=im+1 
                        !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0 .AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwsp**2)  ,mom,cnt,var)		
                        !        write(name(im),'("lwagesp ALL jed ",tr9,2i4)') ier,iesp
                        !        weights(im)=wwage
                        !        calcorr(im)=5
                        !        im=im+1 
                        !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0  .AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwr*dat(MNA:MXAD,:)%logwsp)  ,mom,cnt,var)		
                        !        write(name(im),'("LW*LWSP ALL JED ",tr9,2i4)') ier,iesp
                        !        weights(im)=wcorr
                        !        calcorr(im)=5
                        !        im=im+1
                        !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0  .AND. dat(MNA:MXAD,:)%edr==ier .AND. dat(MNA:MXAD,:)%edsp==iesp), d1*  (dat(MNA:MXAD,:)%logwr*dat(MNA:MXAD,:)%logwsp)  ,mom,cnt,var)		
                        !        write(name(im),'("LW*LWSP ALL JED ",tr9,2i4)') ier,iesp
                        !        weights(im)=wcorr
                        !        calcorr(im)=5
                        !        im=im+1
                        !    end do 
                        !end do 
                        !END IF 
                    end if !j rel
                            ! IF (.NOT.SAVETIME) THEN
                                    !do i=1,20 !rellen
                                    !    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0 .AND. dat(MNA:MXAD,:)%rellen==i ), d1*  (dat(MNA:MXAD,:)%logwr)  ,mom,cnt,var)		
                                    !    write(name(im),'("lwager by len ",tr9,i4)') i	 
                                    !    weights(im)=0.0_dp
                                    !    calcorr(im)=1
                                    !    im=im+1 
                                    !    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0  .AND. dat(MNA:MXAD,:)%rellen==i), d1*  (dat(MNA:MXAD,:)%logwsp)  ,mom,cnt,var)		
                                    !    write(name(im),'("lwagesp by len ",tr9,i4)') i	
                                    !    weights(im)=0.0_dp
                                    !    calcorr(im)=5
                                    !    im=im+1 
                                    !    call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%hhsp==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%logwsp>=0  .AND. dat(MNA:MXAD,:)%rellen==i), d1*  (dat(MNA:MXAD,:)%logwr*dat(MNA:MXAD,:)%logwsp)  ,mom,cnt,var)		
                                    !    write(name(im),'("LW*LWSP by len ",tr9,i4)') i	
                                    !    weights(im)=wcorr
                                    !    calcorr(im)=5
                                    !    im=im+1
                                    !end do
                            !  END IF !SAVETIME
                !end if !REL=1


                headloc(ihead)=im
                if (g==1.and.j==0.and.typosto==0) headstr(ihead)='single men wage all types'
                if (g==1.and.j==1.and.typosto==0) headstr(ihead)='married men wage all types'
                if (g==2.and.j==0.and.typosto==0) headstr(ihead)='single fem wage all types'
                if (g==2.and.j==1.and.typosto==0) headstr(ihead)='married fem wage all types'
                if (g==1.and.j==0.and.typosto==1) headstr(ihead)='single men wage type 1'
                if (g==1.and.j==1.and.typosto==1) headstr(ihead)='married men wage type 1'
                if (g==2.and.j==0.and.typosto==1) headstr(ihead)='single fem wage type 1'
                if (g==2.and.j==1.and.typosto==1) headstr(ihead)='married fem wage type 1'
                if (g==1.and.j==0.and.typosto==2) headstr(ihead)='single men wage type 2'
                if (g==1.and.j==1.and.typosto==2) headstr(ihead)='married men wage type 2'
                if (g==2.and.j==0.and.typosto==2) headstr(ihead)='single fem wage type 2'
                if (g==2.and.j==1.and.typosto==2) headstr(ihead)='married fem wage type 2'
                if (g==1.and.j==0.and.typosto==3) headstr(ihead)='single men wage type 3'
                if (g==1.and.j==1.and.typosto==3) headstr(ihead)='married men wage type 3'
                if (g==2.and.j==0.and.typosto==3) headstr(ihead)='single fem wage type 3'
                if (g==2.and.j==1.and.typosto==3) headstr(ihead)='married fem wage type 3'
                if (g==1.and.j==0.and.typosto==4) headstr(ihead)='single men wage type 4'
                if (g==1.and.j==1.and.typosto==4) headstr(ihead)='married men wage type 4'
                if (g==2.and.j==0.and.typosto==4) headstr(ihead)='single fem wage type 4'
                if (g==2.and.j==1.and.typosto==4) headstr(ihead)='married fem wage type 4'
                ihead=ihead+1


                do ia=agestart(NOCOLLEGE),mxad-1,2 !DETAILEDOUTPUT
                    CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%logwr>=0 ) ,d1*dat(ia,:)%logwr,mom,cnt,var)
                    WRITE(name(im),'("lnw  ned by age",tr1,I2)') ia
                    weights(im)=wwage 
                    calcvar(im)=0
                    im=im+1
                    !CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%logwr>=0 ) ,d1*(dat(ia,:)%logwr**2),mom,cnt,var)
                    !WRITE(name(im),'("w2  ned by age",tr1,I2)') ia
                    !weights(im)=wwvar
                    !calcvar(im)=5
                    !im=im+1
                end do  
                do ia=agestart(COLLEGE),mxad-1,2 !DETAILEDOUTPUT
                    CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%logwr>=0 ) ,d1*dat(ia,:)%logwr,mom,cnt,var)
                    WRITE(name(im),'("lnw  ed by age",tr1,I2)') ia
                    weights(im)=wwage 
                    calcvar(im)=0
                    im=im+1
                    !CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%logwr>=0 ) ,d1*(dat(ia,:)%logwr**2),mom,cnt,var)
                    !WRITE(name(im),'("w2  ed by age",tr1,I2)') ia
                    !weights(im)=wwvar
                    !calcvar(im)=5 
                    !im=im+1
                end do  
                do ia=agestart(NOCOLLEGE),mxad-1,DETAILEDOUTPUT
                    CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%logwr>=0 ) ,d1*dat(ia,:)%logwr,mom,cnt,var)
                    WRITE(name(im),'("lnw  by age",tr1,I2)') ia
                    weights(im)=wwage 
                    calcvar(im)=0
                    im=im+1
                    !CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%logwr>=0 ) ,d1*(dat(ia,:)%logwr**2),mom,cnt,var)
                    !WRITE(name(im),'("w2  ed by age",tr1,I2)') ia
                    !weights(im)=wwvar
                    !calcvar(im)=5 
                    !im=im+1
                end do  




                !************************************************************************************
                IF (momdisplay) THEN  
                        !LA=MNA ; UA=MXAD
                        !CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1  .AND. dat(LA:UA,:)%wr>=0 ) ,d1*dat(LA:UA,:)%wr,mom,cnt,var)
                        !WRITE(name(im),'("w  ALL AGES",tr1,I2)') ia
                        !weights(im)=wwage 
                        !momwhich(im,1)=10 ; momwhich(im,2)=-5 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7:11)=99 
                        !im=im+1
                        !CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==1 .AND. dat(LA:UA,:)%wr>=0 ) ,d1*dat(LA:UA,:)%wr,mom,cnt,var)
                        !WRITE(name(im),'("w  ned ALL AGES",tr1,I2)') ia
                        !weights(im)=wwage 
                        !momwhich(im,1)=10 ; momwhich(im,2)=-5 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=1 ; momwhich(im,7:11)=99 
                        !im=im+1
                        !CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==2 .AND. dat(LA:UA,:)%wr>=0 ) ,d1*dat(LA:UA,:)%wr,mom,cnt,var)
                        !WRITE(name(im),'("w  ed ALL AGES",tr1,I2)') ia
                        !weights(im)=wwage 
                        !momwhich(im,1)=10 ; momwhich(im,2)=-5 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=2 ; momwhich(im,7:11)=99 
                        !im=im+1
        !
        !
                        !LA=20 ; UA=25
                        !CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1  .AND. dat(LA:UA,:)%logwr>=0 ) ,d1*dat(LA:UA,:)%logwr,mom,cnt,var)
                        !WRITE(name(im),'("lnw  LA UA ",2I4)') LA,UA
                        !weights(im)=wwage 
                        !momwhich(im,1)=10 ; momwhich(im,2)=-4 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7:11)=99 
                        !im=im+1
                        !CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==1 .AND. dat(LA:UA,:)%logwr>=0 ) ,d1*dat(LA:UA,:)%logwr,mom,cnt,var)
                        !WRITE(name(im),'("lnw NED LA UA ",2I4)') LA,UA
                        !weights(im)=wwage 
                        !momwhich(im,1)=10 ; momwhich(im,2)=-4 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=1 ; momwhich(im,7:11)=99 
                        !im=im+1
                        !CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==2 .AND. dat(LA:UA,:)%logwr>=0 ) ,d1*dat(LA:UA,:)%logwr,mom,cnt,var)
                        !WRITE(name(im),'("lnw ED LA UA ",2I4)') LA,UA
                        !weights(im)=wwage 
                        !momwhich(im,1)=10 ; momwhich(im,2)=-4 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=2 ; momwhich(im,7:11)=99 
                        !im=im+1
        !
        !
                        !LA=26 ; UA=34
                        !CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1  .AND. dat(LA:UA,:)%logwr>=0 ) ,d1*dat(LA:UA,:)%logwr,mom,cnt,var)
                        !WRITE(name(im),'("lnw  LA UA ",2I4)') LA,UA
                        !weights(im)=wwage 
                        !momwhich(im,1)=10 ; momwhich(im,2)=-3 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7:11)=99 
                        !im=im+1
                        !CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==1 .AND. dat(LA:UA,:)%logwr>=0 ) ,d1*dat(LA:UA,:)%logwr,mom,cnt,var)
                        !WRITE(name(im),'("lnw NED LA UA ",2I4)') LA,UA
                        !weights(im)=wwage 
                        !momwhich(im,1)=10 ; momwhich(im,2)=-3 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=1 ; momwhich(im,7:11)=99 
                        !im=im+1
                        !CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==2 .AND. dat(LA:UA,:)%logwr>=0 ) ,d1*dat(LA:UA,:)%logwr,mom,cnt,var)
                        !WRITE(name(im),'("lnw ED LA UA ",2I4)') LA,UA
                        !weights(im)=wwage 
                        !momwhich(im,1)=10 ; momwhich(im,2)=-3 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=2 ; momwhich(im,7:11)=99 
                        !im=im+1
!
!
                        !!1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        !!1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        do ier=1,2
                            do ia=MNAD,MXAD
                                CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%wr>=0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                                WRITE(name(im),'("w by ed by ia ",2i4)') ier,ia
                                weights(im)=wwage 
                                momwhich(im,1)=10 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99  
                                im=im+1
                                CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%wr>=0 .AND. nummovemar2(ia,:)==0) ,d1*dat(ia,:)%wr,mom,cnt,var)
                                WRITE(name(im),'("w by ed mvemar=0 by ia ",2i4)') ier,ia
                                weights(im)=wwage 
                                momwhich(im,1)=10 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99  
                                im=im+1
                                CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%wr>=0 .AND. nummovemar2(ia,:)>=1) ,d1*dat(ia,:)%wr,mom,cnt,var)
                                WRITE(name(im),'("w by ed mvemar>=1 by ia ",2i4)') ier,ia
                                weights(im)=wwage 
                                momwhich(im,1)=10 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99  
                                im=im+1
                            end do
                        end do 
                        do ier=1,2
                            do ia=MNAD,MXAD
                                CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%wr>=0 .AND. nummovesin2(ia,:)==0) ,d1*dat(ia,:)%wr,mom,cnt,var)
                                WRITE(name(im),'("w by ed mvesin=0 by ia ",2i4)') ier,ia
                                weights(im)=wwage 
                                momwhich(im,1)=10 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=3 ; momwhich(im,11)=99  
                                im=im+1
                                CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%wr>=0 .AND. nummovesin2(ia,:)>=1) ,d1*dat(ia,:)%wr,mom,cnt,var)
                                WRITE(name(im),'("w by ed mvesin>=1 by ia ",2i4)') ier,ia
                                weights(im)=wwage 
                                momwhich(im,1)=10 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=5 ; momwhich(im,11)=99  
                                im=im+1
                            end do 
                        end do
                        !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        LA=MNAD ; UA=49
                        do ier=1,2
                            CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%wr>=0 ) ,d1*dat(LA:UA,:)%wr,mom,cnt,var)
                            WRITE(name(im),'("w by ed 1749 ",I2)') ier
                            weights(im)=wwage 
                            momwhich(im,1)=10 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99  
                            im=im+1
                            CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%wr>=0 .AND. nummovemar2(LA:UA,:)==0) ,d1*dat(LA:UA,:)%wr,mom,cnt,var)
                            WRITE(name(im),'("w by ed mvemar=0 1749 ",I2)') ier
                            weights(im)=wwage 
                            momwhich(im,1)=10 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99  
                            im=im+1
                            CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%wr>=0 .AND. nummovemar2(LA:UA,:)>=1) ,d1*dat(LA:UA,:)%wr,mom,cnt,var)
                            WRITE(name(im),'("w by ed mvemar>=1 1749 ",I2)') ier
                            weights(im)=wwage 
                            momwhich(im,1)=10 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99  
                            im=im+1
                        end do
                        LA=MNAD ; UA=49
                        do ier=1,2
                            CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%wr>=0 .AND. nummovesin2(LA:UA,:)==0) ,d1*dat(LA:UA,:)%wr,mom,cnt,var)
                            WRITE(name(im),'("w by ed mvesin=0 1749 ",I2)') ier
                            weights(im)=wwage 
                            momwhich(im,1)=10 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=3 ; momwhich(im,11)=99  
                            im=im+1
                            CALL condmom(im,( cosexrel(LA:UA,:) .AND.   dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%wr>=0 .AND. nummovesin2(LA:UA,:)>=1) ,d1*dat(LA:UA,:)%wr,mom,cnt,var)
                            WRITE(name(im),'("w by ed mvesin>=1 1749 ",I2)') ier
                            weights(im)=wwage 
                            momwhich(im,1)=10 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=5 ; momwhich(im,11)=99  
                            im=im+1
                        end do
                        

 
                        IF (J==1) THEN !ONLY FOR MARRIED
                            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                            do ier=1,2
                                do ia=MNAD,34,1
                                    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==1 .AND. dat(ia,:)%hhsp==1 .AND. dat(ia,:)%wr>0 .AND. dat(ia,:)%wsp>0 .AND. dat(ia,:)%edr==ier ), d1*  (dat(ia,:)%wr/dat(ia,:)%wsp)  ,mom,cnt,var)		
                                    write(name(im),'("wr/wsp edr age ",2i4)') ier,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=30 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99 
                                    im=im+1
                                    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==1 .AND. dat(ia,:)%hhsp==1 .AND. dat(ia,:)%wr>0 .AND. dat(ia,:)%wsp>0 .AND. dat(ia,:)%edr==ier .AND. nummovemar2(ia,:)==0), d1*  (dat(ia,:)%wr/dat(ia,:)%wsp)  ,mom,cnt,var)		
                                    write(name(im),'("wr/wsp mvemar=0 edr age ",2i4)') ier,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=30 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                                    im=im+1
                                    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==1 .AND. dat(ia,:)%hhsp==1 .AND. dat(ia,:)%wr>0 .AND. dat(ia,:)%wsp>0 .AND. dat(ia,:)%edr==ier .AND. nummovemar2(ia,:)>=1), d1*  (dat(ia,:)%wr/dat(ia,:)%wsp)  ,mom,cnt,var)		
                                    write(name(im),'("wr/wsp mvemar=1 edr age ",2i4)') ier,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=30 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                                    im=im+1
                                end do 
                            end do
                            LA=MNAD ; UA=49
                            do ier=1,2
                                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA:UA,:)%wr>0 .AND. dat(LA:UA,:)%wsp>0 .AND. dat(LA:UA,:)%edr==ier ), d1*  (dat(LA:UA,:)%wr/dat(LA:UA,:)%wsp)  ,mom,cnt,var)		
                                write(name(im),'("wr/wsp edr 17-49 ",i4)') ier
                                weights(im)=wwage
                                momwhich(im,1)=31 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99 
                                im=im+1
                                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA:UA,:)%wr>0 .AND. dat(LA:UA,:)%wsp>0 .AND. dat(LA:UA,:)%edr==ier .AND. nummovemar2(LA:UA,:)==0), d1*  (dat(LA:UA,:)%wr/dat(LA:UA,:)%wsp)  ,mom,cnt,var)		
                                write(name(im),'("wr/wsp mvemar=0 edr 17-49 ",i4)') ier
                                weights(im)=wwage
                                momwhich(im,1)=31 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                                im=im+1
                                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA:UA,:)%wr>0 .AND. dat(LA:UA,:)%wsp>0 .AND. dat(LA:UA,:)%edr==ier .AND. nummovemar2(LA:UA,:)>=1), d1*  (dat(LA:UA,:)%wr/dat(LA:UA,:)%wsp)  ,mom,cnt,var)		
                                write(name(im),'("wr/wsp mvemar>=1 edr 17-49 ",i4)') ier
                                weights(im)=wwage
                                momwhich(im,1)=31 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                                im=im+1
                            end do 
                        END IF !j=1
                            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        do iesp=1,2 !educ spouse
                            do ier=1,2 !educ r
                                do i=1,12,1
                                    LA=MNAD ; UA=34
                                    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA:UA,:)%wr>0 .AND. dat(LA:UA,:)%wsp>0  .AND. dat(LA:UA,:)%rellen==i.AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%edsp==iesp), d1*  (dat(LA:UA,:)%wsp/dat(LA:UA,:)%wr)  ,mom,cnt,var)		
                                    write(name(im),'("wagesp/wr  ",5i4)') ier,iesp,i,LA,UA
                                    weights(im)=wwage
                                    momwhich(im,1)=32 ; momwhich(im,2)=1734 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=i ; momwhich(im,9:11)=99 
                                    im=im+1 
                                    LA=MNAD ; UA=49
                                    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA:UA,:)%wr>0 .AND. dat(LA:UA,:)%wsp>0  .AND. dat(LA:UA,:)%rellen==i.AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%edsp==iesp), d1*  (dat(LA:UA,:)%wsp/dat(LA:UA,:)%wr)  ,mom,cnt,var)		
                                    write(name(im),'("wagesp/wr  ",5i4)') ier,iesp,i,LA,UA
                                    weights(im)=wwage
                                    momwhich(im,1)=32 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=i ; momwhich(im,9:11)=99 
                                    im=im+1 
                                end do 
                            end do 
                        end do 
                        !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        do iesp=1,2 !educ spouse
                            do ier=1,2 !educ r
                                do ia=MNAD,34,1
                                    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==1 .AND. dat(ia,:)%hhsp==1 .AND. dat(ia,:)%wr>0 .AND. dat(ia,:)%wsp>0 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp), d1*  (dat(ia,:)%wsp/dat(ia,:)%wr)  ,mom,cnt,var)		
                                    write(name(im),'("wagesp/wr jted age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=33 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99 
                                    im=im+1
                                    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==1 .AND. dat(ia,:)%hhsp==1 .AND. dat(ia,:)%wr>0 .AND. dat(ia,:)%wsp>0 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovemar2(ia,:)==0), d1*  (dat(ia,:)%wsp/dat(ia,:)%wr)  ,mom,cnt,var)		
                                    write(name(im),'("wagesp/wr mvemar=0 jted age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=33 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                                    im=im+1
                                    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==1 .AND. dat(ia,:)%hhsp==1 .AND. dat(ia,:)%wr>0 .AND. dat(ia,:)%wsp>0 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovemar2(ia,:)>=1), d1*  (dat(ia,:)%wsp/dat(ia,:)%wr)  ,mom,cnt,var)		
                                    write(name(im),'("wagesp/wr mvemar=1 jted age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=33 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                                    im=im+1
                                end do 
                            end do 
                        end do 
                        !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        LA=MNAD ; UA=49
                        do iesp=1,2 !educ spouse
                            do ier=1,2 !educ r
                                    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA:UA,:)%wr>0 .AND. dat(LA:UA,:)%wsp>0 .AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%edsp==iesp), d1*  (dat(LA:UA,:)%wsp/dat(LA:UA,:)%wr)  ,mom,cnt,var)		
                                    write(name(im),'("wagesp/wr jted age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=33 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99 
                                    im=im+1
                                    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA:UA,:)%wr>0 .AND. dat(LA:UA,:)%wsp>0 .AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%edsp==iesp.AND. nummovemar2(LA:UA,:)==0), d1*  (dat(LA:UA,:)%wsp/dat(LA:UA,:)%wr)  ,mom,cnt,var)		
                                    write(name(im),'("wagesp/wr mvemar=0 jted age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=33 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                                    im=im+1
                                    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA:UA,:)%wr>0 .AND. dat(LA:UA,:)%wsp>0 .AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%edsp==iesp.AND. nummovemar2(LA:UA,:)>=1), d1*  (dat(LA:UA,:)%wsp/dat(LA:UA,:)%wr)  ,mom,cnt,var)		
                                    write(name(im),'("wagesp/wr mvemar=1 jted age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=33 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                                    im=im+1
                            end do 
                        end do 
                END IF !momdisplay
                !************************************************************************************

                CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==1 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*dat(mna:mxa,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("w|noed",tr1)') 
                weights(im)=wwage 
                calcvar(im)=1
                im=im+1
                CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==1 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*  (dat(mna:mxa,:)%logwr**2),mom,cnt,var)
                WRITE(name(im),'("wvar|noed",tr1)') 
                weights(im)=wwvar           
                calcvar(im)=5
                im=im+1
                !CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==1 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*  (dat(mna:mxa,:)%logwr-mom(im-2))**2,mom,cnt,var)
                !WRITE(name(im),'("wrng|noed",tr1)') 
                !weights(im)=0.0_dp
                !im=im+1            
                CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==2 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*dat(mna:mxa,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("w|  ed",tr1)') 
                weights(im)=wwage 
                calcvar(im)=1
                im=im+1
                CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==2 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*  (dat(mna:mxa,:)%logwr**2),mom,cnt,var)
                WRITE(name(im),'("wvar|ed",tr1)') 
                weights(im)=wwvar 
                calcvar(im)=5
                im=im+1
                LA=18 ; UA=23
                do ddd=1,ndecile-1
                    CALL condmom(im,(   cosexrel(LA:UA,:) .AND.  dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==1 .AND. dat(LA:UA,:)%logwr>=0 ) ,d1*one( (dat(LA:UA,:)%logwr>decilegrid(ddd) .and. dat(LA:UA,:)%logwr<decilegrid(ddd+1) ) ),mom,cnt,var)
                    WRITE(name(im),'("wdecile|ned",tr1,i4)') ddd
                    weights(im)=wwdecile !; if (onlysingles.and.j==1) weights(im)=0.0_dp !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess 
                    im=im+1
                end do
                LA=20 ; UA=25
                do ddd=1,ndecile-1
                    CALL condmom(im,(   cosexrel(LA:UA,:) .AND.  dat(LA:UA,:)%hhr==1 .AND. dat(LA:UA,:)%edr==2 .AND. dat(LA:UA,:)%logwr>=0 ) ,d1*one( (dat(LA:UA,:)%logwr>decilegrid(ddd) .and. dat(LA:UA,:)%logwr<decilegrid(ddd+1) ) ),mom,cnt,var)
                    WRITE(name(im),'("wdecile| ed",tr1,i4)') ddd
                    weights(im)=wwdecile !; if (onlysingles.and.j==1) weights(im)=0.0_dp !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess 
                    im=im+1
                end do


                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid

                !if (momdisplay) then 
                    do ia=agestart(NOCOLLEGE),MXAD-5,4
                        call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%hhr>=0  .AND. dat(ia,:)%edr==2),   d1*one( dat(ia,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("emp by age EDNED ",i4)') ia
                        momwhich(im,1)=888 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=whour
                        im=im+1 
                    end do 
                    do ia=agestart(NOCOLLEGE),MXAD-5,4
                        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr>=0  .AND. dat(ia,:)%edr==1),   d1*one( dat(ia,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("emp by age ned ",i4)') ia
                        momwhich(im,1)=888 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=1 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=whour
                        im=im+1 
                    end do 
                    do ia=agestart(COLLEGE),MXAD-5,4
                        call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%hhr>=0  .AND. dat(ia,:)%edr==2),   d1*one( dat(ia,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("emp by age ed ",i4)') ia
                        momwhich(im,1)=888 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=2 ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=whour
                        im=im+1 
                    end do 
                !end if !momdisplay





                headloc(ihead)=im
                if (g==1.and.j==0.and.typosto==0) headstr(ihead)='single men  wdif all types'
                if (g==1.and.j==1.and.typosto==0) headstr(ihead)='married men wdif  all types'
                if (g==2.and.j==0.and.typosto==0) headstr(ihead)='single fem  wdif all types'
                if (g==2.and.j==1.and.typosto==0) headstr(ihead)='married fem wdif  all types'
                if (g==1.and.j==0.and.typosto==1) headstr(ihead)='single men  wdif type 1'
                if (g==1.and.j==1.and.typosto==1) headstr(ihead)='married men wdif  type 1'
                if (g==2.and.j==0.and.typosto==1) headstr(ihead)='single fem  wdif type 1'
                if (g==2.and.j==1.and.typosto==1) headstr(ihead)='married fem wdif  type 1'
                if (g==1.and.j==0.and.typosto==2) headstr(ihead)='single men  wdif type 2'
                if (g==1.and.j==1.and.typosto==2) headstr(ihead)='married men wdif  type 2'
                if (g==2.and.j==0.and.typosto==2) headstr(ihead)='single fem  wdif type 2'
                if (g==2.and.j==1.and.typosto==2) headstr(ihead)='married fem wdif  type 2'
                if (g==1.and.j==0.and.typosto==3) headstr(ihead)='single men  wdif type 3'
                if (g==1.and.j==1.and.typosto==3) headstr(ihead)='married men wdif  type 3'
                if (g==2.and.j==0.and.typosto==3) headstr(ihead)='single fem  wdif type 3'
                if (g==2.and.j==1.and.typosto==3) headstr(ihead)='married fem wdif  type 3'
                if (g==1.and.j==0.and.typosto==4) headstr(ihead)='single men  wdif type 4'
                if (g==1.and.j==1.and.typosto==4) headstr(ihead)='married men wdif  type 4'
                if (g==2.and.j==0.and.typosto==4) headstr(ihead)='single fem  wdif type 4'
                if (g==2.and.j==1.and.typosto==4) headstr(ihead)='married fem wdif  type 4'
                ihead=ihead+1


                if (.not.savetime) then
                    LA=MNA ; UA=MXAD
                    call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*wdif(LA:UA,:) ,mom,cnt,var)		
                    write(name(im),'("wdif | stay ",tr2)')  
                    weights(im)=wdifww  
                    calcvar(im)=0 !1
                    im=im+1 

                    LA=18 ; UA=18
                    call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                    write(name(im),'("wdif | move ",tr2,2I4)')  LA,UA
                    weights(im)=wdifww 
                    calcvar(im)=0 !dont forget to set these to 0 if there is no wdif2 following this moment 
                    im=im+1  
                    LA=18 ; UA=20
                    call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                    write(name(im),'("wdif | move ",tr2,2I4)')  LA,UA
                    weights(im)=wdifww 
                    calcvar(im)=0 !dont forget to set these to 0 if there is no wdif2 following this moment 
                    im=im+1  
                    LA=20 ; UA=22
                    call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                    write(name(im),'("wdif | move ",tr2,2I4)')  LA,UA
                    weights(im)=wdifww 
                    calcvar(im)=0 !dont forget to set these to 0 if there is no wdif2 following this moment 
                    im=im+1  
                    LA=23 ; UA=25
                    call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                    write(name(im),'("wdif | move ",tr2,2I4)')   LA,UA
                    weights(im)=wdifww 
                    calcvar(im)=0 !dont forget to set these to 0 if there is no wdif2 following this moment 
                    im=im+1  
                    LA=26 ; UA=28
                    call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                    write(name(im),'("wdif | move ",tr2,2I4)')   LA,UA
                    weights(im)=wdifww 
                    calcvar(im)=0 !dont forget to set these to 0 if there is no wdif2 following this moment 
                    im=im+1  
                end if 

                LA=MNA ; UA=MXAD
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                write(name(im),'("wdif<0 | stay ",tr2)')  
                weights(im)=wdifww  
                calcvar(im)=0
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                write(name(im),'("wdif<0 | move ",tr2)')  
                weights(im)=wdifww  
                calcvar(im)=0
                !momwhich(im,1)=16 ; momwhich(im,2)=-5 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%l/=dat(LA+1:UA+1,:)%hme  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),    d1*one(wdif(LA:UA,:)<0)  ,mom,cnt,var)		
                write(name(im),'("wdif<0 | hmemve=0 ",tr2)')  
                weights(im)=wdifww  
                calcvar(im)=0
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%l==dat(LA+1:UA+1,:)%hme  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),    d1*one(wdif(LA:UA,:)<0)  ,mom,cnt,var)		
                write(name(im),'("wdif<0 | hmemve=1 ",tr2)')  
                weights(im)=wdifww  
                calcvar(im)=0
                im=im+1 




                LA=MNAD ; UA=MXAD
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%l/=dat(LA+1:UA+1,:)%hme  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                write(name(im),'("wdif | hmemve=0 ",tr2)')  
                weights(im)=wdifww
                calcvar(im)=1 !dont forget to set these to 0 if there is no wdif2 following this moment 
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%l/=dat(LA+1:UA+1,:)%hme  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),   d1*( wdif(LA:UA,:)**2 ),mom,cnt,var)		
                write(name(im),'("wdif2 | hmemve=0 ",tr2)')  
                weights(im)=wdifww
                calcvar(im)=5 !dont forget to set these to 0 if there is no wdif2 following this moment 
                im=im+1
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%l==dat(LA+1:UA+1,:)%hme  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*wdif(LA:UA,:) ,mom,cnt,var)		
                write(name(im),'("wdif | hmemve=1 ",tr2)')  
                weights(im)=wdifww
                calcvar(im)=1 !dont forget to set these to 0 if there is no wdif2 following this moment 
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%l==dat(LA+1:UA+1,:)%hme  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),    d1*( wdif(LA:UA,:)**2 ),mom,cnt,var)	
                write(name(im),'("wdif2 | hmemve=1 ",tr2)')  
                weights(im)=wdifww
                calcvar(im)=5 !dont forget to set these to 0 if there is no wdif2 following this moment 
                im=im+1
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                write(name(im),'("wdif | move ",tr2)')  
                weights(im)=wdifww 
                calcvar(im)=1 !dont forget to set these to 0 if there is no wdif2 following this moment 
                im=im+1  
                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1 .AND. move(LA:UA,:)==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*( wdif(LA:UA,:)**2 ),mom,cnt,var)	
                write(name(im),'("wdif2 | move ",tr2)')  
                weights(im)=wdifww
                calcvar(im)=5
                im=im+1


                IF (.NOT.SAVETIME2) THEN
                        LA=MNA-1 ; UA=MXA-2
                        call condmom(im,( cosexrel(LA:UA,:) .AND. deue(LA:UA,:)==1   .AND. dat(LA:UA,:)%l/=dat(LA+1:UA+1,:)%l .AND. dat(LA+1:UA+1,:)%l==dat(LA+2:UA+2,:)%l .AND. dat(LA:UA,:)%l>0  .and. dat(LA+2:UA+2,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*( dat(LA+2:UA+2,:)%logwr-dat(LA:UA,:)%logwr ),mom,cnt,var)		
                        write(name(im),'("wdif | eue,m A ",tr2)')  
                        weights(im)=0.0_dp  
                        calcvar(im)=0
                        im=im+1 
                        !call condmom(im,( cosexrel(LA:UA,:) .AND. deue(LA:UA,:)==1   .AND. dat(LA:UA,:)%l/=dat(LA+1:UA+1,:)%l .AND. dat(LA+1:UA+1,:)%l==dat(LA+2:UA+2,:)%l .AND. dat(LA:UA,:)%l>0 .and. dat(LA+2:UA+2,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*( dat(LA+2:UA+2,:)%logwr-dat(LA:UA,:)%logwr )**2,mom,cnt,var)		
                        !write(name(im),'("wdif2 | eue,m A ",tr2)')  
                        !weights(im)=0.0_dp
                        !calcvar(im)=5
                        !im=im+1 
                        !call condmom(im,( cosexrel(LA:UA,:) .AND. deue(LA:UA,:)==1   .AND. dat(LA:UA,:)%l==dat(LA+1:UA+1,:)%l .AND. dat(LA+1:UA+1,:)%l/=dat(LA+2:UA+2,:)%l .AND. dat(LA:UA,:)%l>0 .and. dat(LA+2:UA+2,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*( dat(LA+2:UA+2,:)%logwr-dat(LA:UA,:)%logwr ),mom,cnt,var)		
                        !write(name(im),'("wdif | eue,m B ",tr2)')  
                        !weights(im)=0.0_dp  
                        !calcvar(im)=1
                        !im=im+1 
                        !call condmom(im,( cosexrel(LA:UA,:) .AND. deue(LA:UA,:)==1   .AND. dat(LA:UA,:)%l==dat(LA+1:UA+1,:)%l .AND. dat(LA+1:UA+1,:)%l/=dat(LA+2:UA+2,:)%l .AND. dat(LA:UA,:)%l>0 .and. dat(LA+2:UA+2,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*( dat(LA+2:UA+2,:)%logwr-dat(LA:UA,:)%logwr )**2,mom,cnt,var)		
                        !write(name(im),'("wdif2 | eue,m B ",tr2)')  
                        !weights(im)=0.0_dp
                        !calcvar(im)=5
                        !im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. deue(LA:UA,:)==1   .AND. dat(LA:UA,:)%l==dat(LA+1:UA+1,:)%l .AND. dat(LA+1:UA+1,:)%l==dat(LA+2:UA+2,:)%l .AND. dat(LA:UA,:)%l>0 .and. dat(LA+2:UA+2,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*( dat(LA+2:UA+2,:)%logwr-dat(LA:UA,:)%logwr ),mom,cnt,var)		
                        write(name(im),'("wdif | eue,s ",tr2)')  
                        weights(im)=0.0_dp  
                        calcvar(im)=0
                        im=im+1 
                        !call condmom(im,( cosexrel(LA:UA,:) .AND. deue(LA:UA,:)==1   .AND. dat(LA:UA,:)%l==dat(LA+1:UA+1,:)%l .AND. dat(LA+1:UA+1,:)%l==dat(LA+2:UA+2,:)%l .AND. dat(LA:UA,:)%l>0 .and. dat(LA+2:UA+2,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*( dat(LA+2:UA+2,:)%logwr-dat(LA:UA,:)%logwr )**2,mom,cnt,var)		
                        !write(name(im),'("wdif2 | eue,s ",tr2)')  
                        !weights(im)=0.0_dp  
                        !calcvar(im)=5
                        !im=im+1 



                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        write(name(im),'("wdif | move ",tr2)')  
                        weights(im)=wdifww 
                        calcvar(im)=0 !dont forget to set these to 0 if there is no wdif2 following this moment 
                        im=im+1  
                        !call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1 .AND. move(LA:UA,:)==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),     d1*( wdif(LA:UA,:)**2 ),mom,cnt,var)		
                        !write(name(im),'("wdif2 | move ",tr2)')  
                        !weights(im)=wdifww
                        !calcvar(im)=5
                        !im=im+1   
                        !call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        !write(name(im),'("wdif | stay ",tr2)')  
                        !weights(im)=wdifww  
                        !calcvar(im)=0 !1
                        !im=im+1 
                        !call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==0  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        !write(name(im),'("wdif | stay,jobwitch==0 ",tr2)')  
                        !weights(im)=wdifww  
                        !calcvar(im)=0 !1
                        !im=im+1 
                        !call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        !write(name(im),'("wdif | stay,jobwitch==1 ",tr2)')  
                        !weights(im)=wdifww  
                        !calcvar(im)=0 !1
                        !im=im+1 

                        !call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1 .AND. move(LA:UA,:)==0 ),   d1*( dat(LA+1:UA+1,:)%logwr-dat(LA:UA,:)%logwr )**2,mom,cnt,var)		
                        !write(name(im),'("wdif2 | stay ",tr2)')  
                        !weights(im)=wdifww
                        !calcvar(im)=5
                        !im=im+1             


                        !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        LA=MNA ; UA=MXAD
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0  .and. dat(LA:UA,:)%edr==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        write(name(im),'("wdif | stay ned ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=0 !1
                        !momwhich(im,1)=15 ; momwhich(im,2)=-5 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=1 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0  .and. dat(LA:UA,:)%edr==2 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        write(name(im),'("wdif | stay ed ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=0 !1
                        !momwhich(im,1)=15 ; momwhich(im,2)=-5 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=2 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA:UA,:)%edr==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        write(name(im),'("wdif | move ned ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=0 !1
                        !momwhich(im,1)=15 ; momwhich(im,2)=-5 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=1 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA:UA,:)%edr==2 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        write(name(im),'("wdif | move ed ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=0 !1
                        !momwhich(im,1)=15 ; momwhich(im,2)=-5 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=2 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        im=im+1 

                        LA=MNAD ; UA=30
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        write(name(im),'("wdif 17:30 ned | move ",tr2)')  
                        weights(im)=wdifww 
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==2 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        write(name(im),'("wdif 17:30 ed | move ",tr2)')  
                        weights(im)=wdifww 
                        im=im+1 

                        if (j==1) then
                            call condmom(im,( cosexrel(LA:UA,:).AND. dee(LA:UA,:)==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA+1:UA+1,:)%hhsp==1 .AND. move(LA:UA,:)==1   .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .and. dat(LA+1:UA+1,:)%logwsp>=0 .AND. dat(LA:UA,:)%logwsp>=0  ),    d1*one( wdif(LA:UA,:)>0 .AND. (dat(LA+1:UA+1,:)%logwsp - dat(LA:UA,:)%logwsp>0) )  ,mom,cnt,var)		
                            write(name(im),'("wrdif>0,wspdif>0 | move ")') 
                            weights(im)=wdifww  
                            calcvar(im)=0 !1
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA+1:UA+1,:)%hhsp==1 .AND. move(LA:UA,:)==1   .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .and. dat(LA+1:UA+1,:)%logwsp>=0 .AND. dat(LA:UA,:)%logwsp>=0  ),     d1*one( wdif(LA:UA,:)>0 .AND. (dat(LA+1:UA+1,:)%logwsp - dat(LA:UA,:)%logwsp<0) )  ,mom,cnt,var)		
                            write(name(im),'("wrdif>0,wspdif<0 | move ")') 
                            weights(im)=wdifww  
                            calcvar(im)=0 !1
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA+1:UA+1,:)%hhsp==1 .AND. move(LA:UA,:)==1   .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .and. dat(LA+1:UA+1,:)%logwsp>=0 .AND. dat(LA:UA,:)%logwsp>=0  ),     d1*one( wdif(LA:UA,:)<0 .AND. (dat(LA+1:UA+1,:)%logwsp - dat(LA:UA,:)%logwsp>0) ) ,mom,cnt,var)		
                            write(name(im),'("wrdif<0,wspdif>0 | move ")') 
                            weights(im)=wdifww  
                            calcvar(im)=0 !1
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1 .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA+1:UA+1,:)%hhsp==1 .AND. move(LA:UA,:)==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .and. dat(LA+1:UA+1,:)%logwsp>=0 .AND. dat(LA:UA,:)%logwsp>=0  ),     d1*one( wdif(LA:UA,:)<0 .AND. (dat(LA+1:UA+1,:)%logwsp - dat(LA:UA,:)%logwsp<0) ) ,mom,cnt,var)		
                            write(name(im),'("wrdif<0,wspdif<0 | move ")') 
                            weights(im)=wdifww  
                            calcvar(im)=0 !1
                            im=im+1 


                            do ier=1,2
                                do iesp=1,2
                                    call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA:UA,:)%edr==ier  .and. dat(LA:UA,:)%edsp==iesp  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                                    write(name(im),'("wdif | move ier,iesp ",2i4)') ier,iesp 
                                    weights(im)=wdifww  
                                    calcvar(im)=0 !1
                                    !momwhich(im,1)=18 ; momwhich(im,2)=-5 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                                    im=im+1 
                                end do 
                            end do 
                            do ier=1,2
                                do iesp=1,2
                                    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhsp==1 .AND. dat(LA+1:UA+1,:)%hhsp==1 .AND. move(LA:UA,:)==1  .and. dat(LA:UA,:)%edr==ier  .and. dat(LA:UA,:)%edsp==iesp  .and. dat(LA+1:UA+1,:)%logwsp>=0 .AND. dat(LA:UA,:)%logwsp>=0  ),    d1*(dat(LA+1:UA+1,:)%logwsp - dat(LA:UA,:)%logwsp) ,mom,cnt,var)		
                                    write(name(im),'("wsp dif | move ier,iesp ",2i4)') ier,iesp 
                                    weights(im)=wdifww  
                                    calcvar(im)=0 !1
                                    im=im+1 
                                end do 
                            end do 

                           do iesp=0,1
                                do ier=0,1
                                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND.dat(MNA:MXAD,:)%l>=0.AND.norelchg(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)>=0  .and. dat(MNA:MXAD,:)%HHR==ier  .and. dat(MNA:MXAD,:)%HHSP==iesp),   d1*one(move(MNA:MXAD,:)==1),mom,cnt,var)		
                                    write(name(im),'("move by EMP,EMPSP",2i4)') ier,iesp 
                                    weights(im)=wmove
                                    im=im+1 
                                end do 
                            end do 
                            LA=20 ; UA=30
                            do iesp=1,2
                                do ier=1,2    
                                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA:UA,:)%edr==ier .and. dat(LA:UA,:)%edsp==iesp .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                                write(name(im),'("wdif | move edr,edsp 20-30 ",2i4)') ier,iesp  
                                weights(im)=wdifww  
                                im=im+1 
                                end do 
                            end do 
                            do iesp=0,1 !spouse emp 
                                call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .AND. dat(LA:UA,:)%HHSP==iesp .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),    d1*wdif(LA:UA,:) ,mom,cnt,var)		
                                write(name(im),'("wdif | move SPOUSE EMP 20-30 ",i4)') iesp  
                                weights(im)=wdifww  
                                im=im+1 
                            end do 
                        end if !j rel 

                        
                        LA=MNA ; UA=MXAD
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        write(name(im),'("wdif | stay ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=0
                        im=im+1 
                        LA=MNA ; UA=MXAD
                        !call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*(wdif(LA:UA,:)**2) ,mom,cnt,var)		
                        !write(name(im),'("wdif2 | stay ",tr2)')  
                        !weights(im)=wdifww  
                        !calcvar(im)=5
                        !im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)<0),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | wdif<0, stay ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=0
                        !momwhich(im,1)=16 ; momwhich(im,2)=-5 ; momwhich(im,3)=g ; momwhich(im,4)=J ; momwhich(im,5)=typosto ; momwhich(im,6)=6 ; momwhich(im,7)=6 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        im=im+1 
                        !call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)<0),   d1*(wdif(LA:UA,:)**2) ,mom,cnt,var)		
                        !write(name(im),'("wdif2 | wdif<0, stay ",tr2)')  
                        !weights(im)=wdifww  
                        !calcvar(im)=5
                        !im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)>0 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | wdif>0, stay ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=0
                        im=im+1 
                        !call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  .AND. wdif(LA:UA,:)>0),   d1*(wdif(LA:UA,:)**2) ,mom,cnt,var)		
                        !write(name(im),'("wdif2 | wdif>0, stay ",tr2)')  
                        !weights(im)=wdifww  
                        !calcvar(im)=5
                        !im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)<0),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | wdif<0, move ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=0
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)>0),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | wdif>0, move ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=0
                        im=im+1 



                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==0  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        write(name(im),'("wdif 0 | stay,jobwitch==0 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=1
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==0  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*(wdif(LA:UA,:)**2) ,mom,cnt,var)		
                        write(name(im),'("wdif2 0 | stay,jobwitch==0 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=5
                        im=im+1 

                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==0  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)<0 ),   d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        write(name(im),'("wdif 0 | stay,jobwitch==0,wdif<0 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=1
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==0  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)<0 ),   d1*(wdif(LA:UA,:)**2) ,mom,cnt,var)		
                        write(name(im),'("wdif2 0 | stay,jobwitch==0,wdif<0 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=5
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==0  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)>0 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif 0 | stay,jobwitch==0,wdif>0 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=1
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==0  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)>0 ),   d1*(wdif(LA:UA,:)**2) ,mom,cnt,var)		
                        write(name(im),'("wdif2 0 | stay,jobwitch==0,wdif>0 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=5
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==0  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                        write(name(im),'("wdif<0 | stay,jobwitch==0 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=0
                        im=im+1 

                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*wdif(LA:UA,:) ,mom,cnt,var)		
                        write(name(im),'("wdif 1 | stay,jobwitch==1 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=1
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*(wdif(LA:UA,:)**2) ,mom,cnt,var)		
                        write(name(im),'("wdif2 1 | stay,jobwitch==1 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=5
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)<0 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif 1 | stay,jobwitch==1,wdif<0 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=1
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)<0 ),   d1*(wdif(LA:UA,:)**2) ,mom,cnt,var)		
                        write(name(im),'("wdif2 1 | stay,jobwitch==1,wdif<0 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=5
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)>0 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif 1 | stay,jobwitch==1,wdif>0 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=1
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. wdif(LA:UA,:)>0 ),   d1*(wdif(LA:UA,:)**2) ,mom,cnt,var)		
                        write(name(im),'("wdif2 1 | stay,jobwitch==1,wdif>0 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=5
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA:UA,:)%jobswitch==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0  ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                        write(name(im),'("wdif<0 | stay,jobwitch==1 ",tr2)')  
                        weights(im)=wdifww  
                        calcvar(im)=0
                        im=im+1 
                END IF !SAVETIME



               if (momdisplay) then 
                    headloc(ihead)=im
                    if (g==1.and.j==0.and.typosto==0) headstr(ihead)='single men  emax all types'
                    if (g==1.and.j==1.and.typosto==0) headstr(ihead)='married men emax  all types'
                    if (g==2.and.j==0.and.typosto==0) headstr(ihead)='single fem  emax all types'
                    if (g==2.and.j==1.and.typosto==0) headstr(ihead)='married fem emax  all types'
                    if (g==1.and.j==0.and.typosto==1) headstr(ihead)='single men  emax type 1'
                    if (g==1.and.j==1.and.typosto==1) headstr(ihead)='married men emax  type 1'
                    if (g==2.and.j==0.and.typosto==1) headstr(ihead)='single fem  emax type 1'
                    if (g==2.and.j==1.and.typosto==1) headstr(ihead)='married fem emax  type 1'
                    if (g==1.and.j==0.and.typosto==2) headstr(ihead)='single men  emax type 2'
                    if (g==1.and.j==1.and.typosto==2) headstr(ihead)='married men emax  type 2'
                    if (g==2.and.j==0.and.typosto==2) headstr(ihead)='single fem  emax type 2'
                    if (g==2.and.j==1.and.typosto==2) headstr(ihead)='married fem emax  type 2'
                    if (g==1.and.j==0.and.typosto==3) headstr(ihead)='single men  emax type 3'
                    if (g==1.and.j==1.and.typosto==3) headstr(ihead)='married men emax  type 3'
                    if (g==2.and.j==0.and.typosto==3) headstr(ihead)='single fem  emax type 3'
                    if (g==2.and.j==1.and.typosto==3) headstr(ihead)='married fem emax  type 3'
                    if (g==1.and.j==0.and.typosto==4) headstr(ihead)='single men  emax type 4'
                    if (g==1.and.j==1.and.typosto==4) headstr(ihead)='married men emax  type 4'
                    if (g==2.and.j==0.and.typosto==4) headstr(ihead)='single fem  emax type 4'
                    if (g==2.and.j==1.and.typosto==4) headstr(ihead)='married fem emax  type 4'
                    ihead=ihead+1
                    do ia=MNAD,25
                        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==1 ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax ned ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=1  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==2 ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax ed ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=2  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:) ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax  ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==1 .AND. nummovemar2(ia,:)==0 ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax mvemar=0 ned ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=1  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==2  .AND. nummovemar2(ia,:)==0),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax mvemar=0 ed ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=2  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:)  .AND.nummovemar2(ia,:)==0),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax  mvemar=0 ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==1 .AND. nummovemar2(ia,:)>=1 ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax mvemar=1 ned ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=1  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==2  .AND. nummovemar2(ia,:)>=1 ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax mvemar=1 ed ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=2  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:)  .AND. nummovemar2(ia,:)>=1 ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax  mvemar=1 ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==1 .AND. nummovesin2(ia,:)==0 ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax mvesin=0 ned ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=1  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=3 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==2  .AND. nummovesin2(ia,:)==0),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax mvesin=0 ed ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=2  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=3 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:)  .AND. nummovesin2(ia,:)==0),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax  mvesin=0 ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=3 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==1 .AND. nummovesin2(ia,:)>=1 ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax mvesin=1 ned ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=1  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=5 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==2  .AND. nummovesin2(ia,:)>=1 ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax mvesin=1 ed ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=2  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=5 ; momwhich(im,11)=99 
                        im=im+1
                        call condmom(im,( cosexrel(ia,:)  .AND. nummovesin2(ia,:)>=1 ),   d1* dat(ia,:)%emax ,mom,cnt,var)	
                        write(name(im),'("emax  mvesin=1 ",tr5,i4)') ia  
                        weights(im)=0.0_dp
                        momwhich(im,1)=80 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=6  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=5 ; momwhich(im,11)=99 
                        im=im+1
                    end do


!                    
!                    !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                    if (j==1) then
                        do iesp=1,2 !educ spouse
                            do ier=1,2 !educ r
                                do ia=MNAD,25
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp ), d1*  (dat(ia,:)%emax)  ,mom,cnt,var)		
                                    write(name(im),'("emaxr by jted, age   ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=40 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99         
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp), d1*  (dat(ia,:)%emaxsp)  ,mom,cnt,var)		
                                    write(name(im),'("emaxsp by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=41 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp   ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99   
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp), d1* (dat(ia,:)%emaxsp/(dat(ia,:)%emax+dat(ia,:)%emaxsp) )  ,mom,cnt,var)		
                                    write(name(im),'("emaxsp/emaxsum by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=42 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp   ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99   
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp .AND. nummovemar2(ia,:)==0), d1*  (dat(ia,:)%emax)  ,mom,cnt,var)		
                                    write(name(im),'("emaxr mvemar=0 by jted, age   ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=40 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99         
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovemar2(ia,:)==0), d1*  (dat(ia,:)%emaxsp)  ,mom,cnt,var)		
                                    write(name(im),'("emaxsp mvemar=0 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=41 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp   ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99   
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovemar2(ia,:)==0), d1* (dat(ia,:)%emaxsp/(dat(ia,:)%emax+dat(ia,:)%emaxsp) )  ,mom,cnt,var)		
                                    write(name(im),'("emaxsp/emaxsum mvemar=0 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=42 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp   ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99   
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp .AND. nummovemar2(ia,:)>=1), d1*  (dat(ia,:)%emax)  ,mom,cnt,var)		
                                    write(name(im),'("emaxr mvemar=1 by jted, age   ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=40 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp   ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99           
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovemar2(ia,:)>=1), d1*  (dat(ia,:)%emaxsp)  ,mom,cnt,var)		
                                    write(name(im),'("emaxsp mvemar=1 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=41 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp    ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99    
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovemar2(ia,:)>=1), d1* (dat(ia,:)%emaxsp/(dat(ia,:)%emax+dat(ia,:)%emaxsp) )  ,mom,cnt,var)		
                                    write(name(im),'("emaxsp/emaxsum mvemar=1 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=42 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp    ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99   
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp .AND. nummovesin2(ia,:)==0), d1*  (dat(ia,:)%emax)  ,mom,cnt,var)		
                                    write(name(im),'("emaxr mvesin=0 by jted, age   ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=40 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp    ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=3 ; momwhich(im,11)=99            
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovesin2(ia,:)==0), d1*  (dat(ia,:)%emaxsp)  ,mom,cnt,var)		
                                    write(name(im),'("emaxsp mvesin=0 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=41 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp    ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=3 ; momwhich(im,11)=99   
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovesin2(ia,:)==0), d1* (dat(ia,:)%emaxsp/(dat(ia,:)%emax+dat(ia,:)%emaxsp) )  ,mom,cnt,var)		
                                    write(name(im),'("emaxsp/emaxsum mvesin=0 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=42 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp    ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=3 ; momwhich(im,11)=99    
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp .AND. nummovesin2(ia,:)>=1), d1*  (dat(ia,:)%emax)  ,mom,cnt,var)		
                                    write(name(im),'("emaxr mvesin=1 by jted, age   ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=40 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp    ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=5 ; momwhich(im,11)=99            
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovesin2(ia,:)>=1), d1*  (dat(ia,:)%emaxsp)  ,mom,cnt,var)		
                                    write(name(im),'("emaxsp mvesin=1 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=41 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp    ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=5 ; momwhich(im,11)=99   
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%emax>0 .AND. dat(ia,:)%emaxsp>0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovesin2(ia,:)>=1), d1* (dat(ia,:)%emaxsp/(dat(ia,:)%emax+dat(ia,:)%emaxsp) )  ,mom,cnt,var)		
                                    write(name(im),'("emaxsp/emaxsum mvesin=1 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=42 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp    ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=5 ; momwhich(im,11)=99   
                                    im=im+1 
                                end do !ia
                            end do !ier
                        end do !iesp


         !             !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        do iesp=1,2 !educ spouse
                            do ier=1,2 !educ r
                                do ia=MNAD,30
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp), d1*  (dat(ia,:)%consr)  ,mom,cnt,var)		
                                    write(name(im),'("consr by jted,age   ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=70 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99              
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp), d1*  (dat(ia,:)%consp)  ,mom,cnt,var)		
                                    write(name(im),'("consp by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=71 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99   
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp), d1*  (dat(ia,:)%consr/dat(ia,:)%incsum)  ,mom,cnt,var)		
                                    write(name(im),'("consr/incsum by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=72 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp   ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99   
                                    im=im+1 

                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp .AND. nummovemar2(ia,:)==0), d1*  (dat(ia,:)%consr)  ,mom,cnt,var)		
                                    write(name(im),'("consr mvemar=0 by jted,age   ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=70 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99              
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovemar2(ia,:)==0), d1*  (dat(ia,:)%consp)  ,mom,cnt,var)		
                                    write(name(im),'("consp mvemar=0 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=71 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99   
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovemar2(ia,:)==0), d1*  (dat(ia,:)%consr/dat(ia,:)%incsum)  ,mom,cnt,var)		
                                    write(name(im),'("consr/incsum mvemar=0 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=72 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp   ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99   
                                    im=im+1 

                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0 .AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp .AND. nummovemar2(ia,:)>=1), d1*  (dat(ia,:)%consr)  ,mom,cnt,var)		
                                    write(name(im),'("consr mvemar>=1 by jted,age   ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=70 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99              
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovemar2(ia,:)>=1), d1*  (dat(ia,:)%consp)  ,mom,cnt,var)		
                                    write(name(im),'("consp mvemar>=1 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=71 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99   
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0.AND. dat(ia,:)%edr==ier .AND. dat(ia,:)%edsp==iesp.AND. nummovemar2(ia,:)>=1), d1*  (dat(ia,:)%consr/dat(ia,:)%incsum)  ,mom,cnt,var)		
                                    write(name(im),'("consr/incsum mvemar>=1 by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=72 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp   ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99   
                                    im=im+1 
                                end do !ia
                            end do !ier
                        end do !iesp
         !             !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        LA=MNAD ; UA=49
                        do iesp=1,2 !educ spouse
                            do ier=1,2 !educ r
                                    call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%consp>=0 .AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%edsp==iesp), d1*  (dat(LA:UA,:)%consr)  ,mom,cnt,var)		
                                    write(name(im),'("consr by jted,age   ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=70 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99              
                                    im=im+1 
                                    call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%consp>=0.AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%edsp==iesp), d1*  (dat(LA:UA,:)%consp)  ,mom,cnt,var)		
                                    write(name(im),'("consp by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=71 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99   
                                    im=im+1 
                                    call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%consp>=0.AND. dat(LA:UA,:)%edr==ier .AND. dat(LA:UA,:)%edsp==iesp), d1*  (dat(LA:UA,:)%consr/dat(LA:UA,:)%incsum)  ,mom,cnt,var)		
                                    write(name(im),'("consr/incsum by jted,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=72 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp   ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99   
                                    im=im+1
                            end do 
                        end do 

                        !BY JOINT EMP STATUS
                        !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        do iesp=0,1 !working or not 
                            do ier=0,1 !working or not
                                do ia=MNAD,30
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0 .AND. dat(ia,:)%hhr==ier .AND. dat(ia,:)%hhsp==iesp), d1*  (dat(ia,:)%consr)  ,mom,cnt,var)		
                                    write(name(im),'("consr by jtemp,age   ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=90 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99 ; momwhich(im,9:11)=99         
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0.AND. dat(ia,:)%hhr==ier .AND. dat(ia,:)%hhsp==iesp), d1*  (dat(ia,:)%consp)  ,mom,cnt,var)		
                                    write(name(im),'("consp by jtemp,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=91 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99 ; momwhich(im,9:11)=99 
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0.AND. dat(ia,:)%hhr==ier .AND. dat(ia,:)%hhsp==iesp), d1*  (dat(ia,:)%consr/dat(ia,:)%incsum)  ,mom,cnt,var)		
                                    write(name(im),'("consr/incsum by jtemp,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=92 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99 ; momwhich(im,9:11)=99 
                                    im=im+1 
                                    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%consp>=0.AND. dat(ia,:)%hhr==ier .AND. dat(ia,:)%hhsp==iesp), d1*  ((dat(ia,:)%consr+dat(ia,:)%consp)/dat(ia,:)%incsum)  ,mom,cnt,var)		
                                    write(name(im),'("(consr+consp)/incsum by jtemp,age ",3i4)') ier,iesp,ia
                                    weights(im)=wwage
                                    momwhich(im,1)=93 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99 ; momwhich(im,9:11)=99 
                                    im=im+1 
                                end do !ia
                            end do !ier
                        end do !iesp
                        !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        do iesp=0,1
                        do ier=0,1 !working or not
                                LA=MNAD ; UA=49
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%consp>=0 .AND. dat(LA:UA,:)%hhr==ier .AND. dat(LA:UA,:)%hhsp==iesp .AND. move(LA:UA,:)==0), d1*  (dat(LA:UA,:)%consr)  ,mom,cnt,var)		
                                write(name(im),'("consr by jtemp mvemar=0 17-49   ",2i4)') ier,iesp
                                weights(im)=wwage
                                momwhich(im,1)=90 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99      
                                im=im+1 
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%consp>=0.AND. dat(LA:UA,:)%hhr==ier .AND. dat(LA:UA,:)%hhsp==iesp .AND. move(LA:UA,:)==0), d1*  (dat(LA:UA,:)%consp)  ,mom,cnt,var)		
                                write(name(im),'("consp by jtemp mvemar=0 17-49 ",2i4)') ier,iesp
                                weights(im)=wwage
                                momwhich(im,1)=91 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                                im=im+1 
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%consp>=0.AND. dat(LA:UA,:)%hhr==ier .AND. dat(LA:UA,:)%hhsp==iesp .AND. move(LA:UA,:)==0), d1*  (dat(LA:UA,:)%consr/dat(LA:UA,:)%incsum)  ,mom,cnt,var)		
                                write(name(im),'("consr/incsum by jtemp mvemar=0 17-49 ",2i4)') ier,iesp
                                weights(im)=wwage
                                momwhich(im,1)=92 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                                im=im+1 
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%consp>=0.AND. dat(LA:UA,:)%hhr==ier .AND. dat(LA:UA,:)%hhsp==iesp .AND. move(LA:UA,:)==0), d1*  ((dat(LA:UA,:)%consr+dat(LA:UA,:)%consp)/dat(LA:UA,:)%incsum)  ,mom,cnt,var)		
                                write(name(im),'("(consr+consp)/incsum by mvemar=0 jtemp 17-49 ",2i4)') ier,iesp
                                weights(im)=wwage
                                momwhich(im,1)=93 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                                im=im+1 


                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%consp>=0 .AND. dat(LA:UA,:)%hhr==ier .AND. dat(LA:UA,:)%hhsp==iesp .AND. move(LA:UA,:)==1), d1*  (dat(LA:UA,:)%consr)  ,mom,cnt,var)		
                                write(name(im),'("consr by jtemp mvemar>=1 17-49   ",2i4)') ier,iesp
                                weights(im)=wwage
                                momwhich(im,1)=90 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99      
                                im=im+1 
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%consp>=0.AND. dat(LA:UA,:)%hhr==ier .AND. dat(LA:UA,:)%hhsp==iesp .AND. move(LA:UA,:)==1), d1*  (dat(LA:UA,:)%consp)  ,mom,cnt,var)		
                                write(name(im),'("consp by jtemp mvemar>=1 17-49 ",2i4)') ier,iesp
                                weights(im)=wwage
                                momwhich(im,1)=91 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                                im=im+1 
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%consp>=0.AND. dat(LA:UA,:)%hhr==ier .AND. dat(LA:UA,:)%hhsp==iesp .AND. move(LA:UA,:)==1), d1*  (dat(LA:UA,:)%consr/dat(LA:UA,:)%incsum)  ,mom,cnt,var)		
                                write(name(im),'("consr/incsum by jtemp mvemar>=1 17-49 ",2i4)') ier,iesp
                                weights(im)=wwage
                                momwhich(im,1)=92 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                                im=im+1 
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%consp>=0.AND. dat(LA:UA,:)%hhr==ier .AND. dat(LA:UA,:)%hhsp==iesp .AND. move(LA:UA,:)==1), d1*  ((dat(LA:UA,:)%consr+dat(LA:UA,:)%consp)/dat(LA:UA,:)%incsum)  ,mom,cnt,var)		
                                write(name(im),'("(consr+consp)/incsum by mvemar>=1 jtemp 17-49 ",2i4)') ier,iesp
                                weights(im)=wwage
                                momwhich(im,1)=93 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=iesp  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                                im=im+1 
                        end do !iesp
                        end do !ier
                    end if !rel


                        !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                    do ier=1,2 !educ r
                        do ia=MNAD,40
                            call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0  .AND. dat(ia,:)%edr==ier ), d1*  (dat(ia,:)%consr)  ,mom,cnt,var)		
                            write(name(im),'("consr by ed,age   ",2i4)') ier,ia
                            weights(im)=wwage
                            momwhich(im,1)=73 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99         
                            im=im+1 
                            call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%edr==ier  ), d1*  (dat(ia,:)%consp)  ,mom,cnt,var)		
                            write(name(im),'("consp  by ed,age ",2i4)') ier,ia
                            weights(im)=wwage
                            momwhich(im,1)=74 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99 
                            im=im+1 
                            call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%edr==ier  ), d1*  (dat(ia,:)%consr/dat(ia,:)%incsum)  ,mom,cnt,var)		
                            write(name(im),'("consr/incsum  by ed,age ",2i4)') ier,ia
                            weights(im)=wwage
                            momwhich(im,1)=75 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99 
                            im=im+1 

                            call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0  .AND. dat(ia,:)%edr==ier .AND. nummovemar(:)==0), d1*  (dat(ia,:)%consr)  ,mom,cnt,var)		
                            write(name(im),'("consr mvemar=0 by ed,age   ",2i4)') ier,ia
                            weights(im)=wwage
                            momwhich(im,1)=73 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99         
                            im=im+1 
                            call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%edr==ier  .AND. nummovemar(:)==0), d1*  (dat(ia,:)%consp)  ,mom,cnt,var)		
                            write(name(im),'("consp mvemar=0  by ed,age ",2i4)') ier,ia
                            weights(im)=wwage
                            momwhich(im,1)=74 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                            im=im+1 
                            call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%edr==ier  .AND. nummovemar(:)==0), d1*  (dat(ia,:)%consr/dat(ia,:)%incsum)  ,mom,cnt,var)		
                            write(name(im),'("consr/incsum mvemar=0  by ed,age ",2i4)') ier,ia
                            weights(im)=wwage
                            momwhich(im,1)=75 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                            im=im+1 

                            call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0  .AND. dat(ia,:)%edr==ier .AND. nummovemar(:)>=1), d1*  (dat(ia,:)%consr)  ,mom,cnt,var)		
                            write(name(im),'("consr mvemar=1 by ed,age   ",2i4)') ier,ia
                            weights(im)=wwage
                            momwhich(im,1)=73 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99         
                            im=im+1 
                            call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%edr==ier  .AND. nummovemar(:)>=1), d1*  (dat(ia,:)%consp)  ,mom,cnt,var)		
                            write(name(im),'("consp mvemar=1  by ed,age ",2i4)') ier,ia
                            weights(im)=wwage
                            momwhich(im,1)=74 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                            im=im+1 
                            call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%consr>=0 .AND. dat(ia,:)%edr==ier  .AND. nummovemar(:)>=1), d1*  (dat(ia,:)%consr/dat(ia,:)%incsum)  ,mom,cnt,var)		
                            write(name(im),'("consr/incsum mvemar=1  by ed,age ",2i4)') ier,ia
                            weights(im)=wwage
                            momwhich(im,1)=75 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                            im=im+1 
                        end do !ia
                    end do !ier
                    !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                    LA=MNAD ; UA=49
                    do ier=1,2 !educ r
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%edr==ier ), d1*  (dat(LA:UA,:)%consr)  ,mom,cnt,var)		
                            write(name(im),'("consr by ed 17-49   ",i4)') ier
                            weights(im)=wwage
                            momwhich(im,1)=73 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99         
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%edr==ier  ), d1*  (dat(LA:UA,:)%consp)  ,mom,cnt,var)		
                            write(name(im),'("consp  by ed 17-49   ",i4)') ier
                            weights(im)=wwage
                            momwhich(im,1)=74 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99 
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%edr==ier  ), d1*  (dat(LA:UA,:)%consr/dat(LA:UA,:)%incsum)  ,mom,cnt,var)		
                            write(name(im),'("consr/incsum  by ed 17-49   ",i4)') ier
                            weights(im)=wwage
                            momwhich(im,1)=75 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=6 ; momwhich(im,11)=99 
                            im=im+1 


                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0  .AND. dat(LA:UA,:)%edr==ier .AND. nummovemar2(LA:UA,:)==0), d1*  (dat(LA:UA,:)%consr)  ,mom,cnt,var)		
                            write(name(im),'("consr mvemar=0 by ed 17-49   ",i4)') ier
                            weights(im)=wwage
                            momwhich(im,1)=73 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99         
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%edr==ier  .AND. nummovemar2(LA:UA,:)==0), d1*  (dat(LA:UA,:)%consp)  ,mom,cnt,var)		
                            write(name(im),'("consp mvemar=0  by ed 17-49   ",i4)') ier
                            weights(im)=wwage
                            momwhich(im,1)=74 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%edr==ier  .AND. nummovemar2(LA:UA,:)==0), d1*  (dat(LA:UA,:)%consr/dat(LA:UA,:)%incsum)  ,mom,cnt,var)		
                            write(name(im),'("consr/incsum mvemar=0  by ed 17-49   ",i4)') ier
                            weights(im)=wwage
                            momwhich(im,1)=75 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=0 ; momwhich(im,11)=99 
                            im=im+1 

                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0  .AND. dat(LA:UA,:)%edr==ier .AND. nummovemar2(LA:UA,:)>=1), d1*  (dat(LA:UA,:)%consr)  ,mom,cnt,var)		
                            write(name(im),'("consr mvemar=1 by ed 17-49   ",i4)') ier
                            weights(im)=wwage
                            momwhich(im,1)=73 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99         
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%edr==ier  .AND. nummovemar2(LA:UA,:)>=1), d1*  (dat(LA:UA,:)%consp)  ,mom,cnt,var)		
                            write(name(im),'("consp mvemar=1  by ed 17-49   ",i4)') ier
                            weights(im)=wwage
                            momwhich(im,1)=74 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%consr>=0 .AND. dat(LA:UA,:)%edr==ier  .AND. nummovemar2(LA:UA,:)>=1), d1*  (dat(LA:UA,:)%consr/dat(LA:UA,:)%incsum)  ,mom,cnt,var)		
                            write(name(im),'("consr/incsum mvemar=1  by ed 17-49   ",i4)') ier
                            weights(im)=wwage
                            momwhich(im,1)=75 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j ; momwhich(im,5)=typosto ; momwhich(im,6)=ier  ; momwhich(im,7)=99  ;  momwhich(im,8)=99  ;  momwhich(im,9)=99 ;  momwhich(im,10)=1 ; momwhich(im,11)=99 
                            im=im+1 
                    end do !ier

                end if !momdisplay




                end do !typosto
            end do !j rel
        end do !g sex
        !END OF MOMENTS BY SEX AND REL 
        !*******************************************************************************************************************************************************************************************************************************



        !*******************************************************************************************************************************************************************************************************************************
        !START OF THE "BY GENDER ONLY" LOOP (NO CONDITIONING ON TYPE EVEN WHEN DETAILED OUTPUT)
        do g=1,2 
            cosex(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co  .and. dat(MNAD:MXA,:)%sexr==g  )			
            wmove=wmoveall
            headloc(ihead)=im
            if (co==1.and.g==1) headstr(ihead)='all men wage by loc '
            if (co==1.and.g==2) headstr(ihead)='all fem wage by loc '
            ihead=ihead+1

            do i=1,nl
                ia=20 ; INCR=5
                call condmom(im,( cosex(ia:ia+INCR,:) .AND. dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%l==i .AND. dat(ia:ia+INCR,:)%logwr>=0 ),   d1*dat(ia:ia+INCR,:)%logwr ,mom,cnt,var)		
                write(name(im),'("w|loc",tr9,2i4)') ia,i			
                weights(im)=wwage
                im=im+1 
                ia=26 ; INCR=8
                call condmom(im,( cosex(ia:ia+INCR,:) .AND. dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%l==i .AND. dat(ia:ia+INCR,:)%logwr>=0 ),   d1*dat(ia:ia+INCR,:)%logwr ,mom,cnt,var)		
                write(name(im),'("w|loc",tr9,2i4)') ia,i			
                weights(im)=wwage 
                im=im+1 
            end do 

            if (.NOT.SAVETIME2) THEN 
                ia=20 ; INCR=5
                CALL condmom(im,(   cosex(ia:ia+INCR,:) .AND.  dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%edr==1 .AND. dat(ia:ia+INCR,:)%logwr>=0 ) ,d1*dat(ia:ia+INCR,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("wned|ia ",i4)') ia
                calcvar(im)=0
                weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                im=im+1
                ia=26 ; INCR=8
                CALL condmom(im,(   cosex(ia:ia+INCR,:) .AND.  dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%edr==1 .AND. dat(ia:ia+INCR,:)%logwr>=0 ) ,d1*dat(ia:ia+INCR,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("wned|ia ",i4)') ia
                calcvar(im)=0
                weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                im=im+1
                ia=35 ; INCR=15
                CALL condmom(im,(   cosex(ia:ia+INCR,:) .AND.  dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%edr==1 .AND. dat(ia:ia+INCR,:)%logwr>=0 ) ,d1*dat(ia:ia+INCR,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("wned|ia ",i4)') ia
                calcvar(im)=0
                weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                im=im+1

                ia=20 ; INCR=5
                CALL condmom(im,(   cosex(ia:ia+INCR,:) .AND.  dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%edr==2 .AND. dat(ia:ia+INCR,:)%logwr>=0 ) ,d1*dat(ia:ia+INCR,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("wed|ia ",i4)') ia
                calcvar(im)=0
                weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                im=im+1
                ia=26 ; INCR=8
                CALL condmom(im,(   cosex(ia:ia+INCR,:) .AND.  dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%edr==2 .AND. dat(ia:ia+INCR,:)%logwr>=0 ) ,d1*dat(ia:ia+INCR,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("wed|ia ",i4)') ia
                calcvar(im)=0
                weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                im=im+1
                ia=35 ; INCR=15
                CALL condmom(im,(   cosex(ia:ia+INCR,:) .AND.  dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%edr==2 .AND. dat(ia:ia+INCR,:)%logwr>=0 ) ,d1*dat(ia:ia+INCR,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("wed|ia ",i4)') ia
                calcvar(im)=0
                weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                im=im+1


                !do i=1,nl
                !    call condmom(im,( cosex(MNA:MXA,:) .AND. dat(MNA:MXA,:)%hhr==1 .AND. dat(MNA:MXA,:)%l==i .AND. dat(MNA:MXA,:)%logwr>=0 ),   d1*dat(MNA:MXA,:)%logwr ,mom,cnt,var)		
                !    write(name(im),'("w|loc",tr9,i4)') i			
                !    weights(im)=0.0_dp 
                !    im=im+1 
                !end do 


                do i=1,nl
                    ia=20 ; INCR=5
                    CALL condmom(im,(   cosex(ia:ia+INCR,:) .AND.  dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%edr==1 .AND. dat(ia:ia+INCR,:)%logwr>=0 .AND. dat(ia:ia+INCR,:)%l==i) ,d1*dat(ia:ia+INCR,:)%logwr,mom,cnt,var)
                    WRITE(name(im),'("wned|ia by loc",2i4)') ia,i
                    calcvar(im)=0
                    weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                    im=im+1
                    ia=26 ; INCR=8
                    CALL condmom(im,(   cosex(ia:ia+INCR,:) .AND.  dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%edr==1 .AND. dat(ia:ia+INCR,:)%logwr>=0 .AND. dat(ia:ia+INCR,:)%l==i) ,d1*dat(ia:ia+INCR,:)%logwr,mom,cnt,var)
                    WRITE(name(im),'("wned|ia by loc",2i4)') ia,i
                    calcvar(im)=0
                    weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                    im=im+1
                end do 

                do i=1,nl
                    ia=20 ; INCR=5
                    CALL condmom(im,(   cosex(ia:ia+INCR,:) .AND.  dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%edr==2 .AND. dat(ia:ia+INCR,:)%logwr>=0 .AND. dat(ia:ia+INCR,:)%l==i) ,d1*dat(ia:ia+INCR,:)%logwr,mom,cnt,var)
                    WRITE(name(im),'("w ed|ia by loc",2i4)') ia,i
                    calcvar(im)=0
                    weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                    im=im+1
                    ia=26 ; INCR=8
                    CALL condmom(im,(   cosex(ia:ia+INCR,:) .AND.  dat(ia:ia+INCR,:)%hhr==1 .AND. dat(ia:ia+INCR,:)%edr==2 .AND. dat(ia:ia+INCR,:)%logwr>=0 .AND. dat(ia:ia+INCR,:)%l==i) ,d1*dat(ia:ia+INCR,:)%logwr,mom,cnt,var)
                    WRITE(name(im),'("w ed|ia by loc",2i4)') ia,i
                    calcvar(im)=0
                    weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                    im=im+1
                end do
            END IF !SAVETIME
        end do !GENDER
        !BY GENDER ONLY
        !**************************************************************************************************************************************






        !*****************************************************************************************
        !MOMENTS BY TYPE
        !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!         !BE CAREFUL WITH THE CONDITIONING STATEMENTS!!!!!
        if (typemoments) then
            coho(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co  )	

            headloc(ihead)=im
            headstr(ihead)='mar by type'
            ihead=ihead+1
            call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 ), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
            write(name(im),'("mar ",tr10)')	
            weights(im)=0. !; if (onlysingles) weights(im)=0.0_dp   !forget about that if statement usually. just putting there to get im to go up by 1 for the next loop
            im=im+1
            call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%edr==1), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
            write(name(im),'("mar ned ",tr10)')	
            weights(im)=0. !; if (onlysingles) weights(im)=0.0_dp   !forget about that if statement usually. just putting there to get im to go up by 1 for the next loop
            im=im+1
            call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%edr==2), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
            write(name(im),'("mar ed ",tr10)')	
            weights(im)=0. !; if (onlysingles) weights(im)=0.0_dp   !forget about that if statement usually. just putting there to get im to go up by 1 for the next loop
            im=im+1



            call condmom(im,( coho(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==1), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move by LA:UA " )') 
            weights(im)=0.0_dp
            im=im+1
            do i=1,ntypp            
                call condmom(im,( coho(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0  .AND. dat(LA:UA,:)%typ==i .AND. move(LA:UA,:)==1), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
                write(name(im),'("getdiv move by LA:UA,typ ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do            

            call condmom(im,( coho(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==0), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv stay LA:UA " )') 
            weights(im)=0.0_dp
            im=im+1
            do i=1,ntypp            
                call condmom(im,( coho(LA:UA,:) .AND. dat(LA:UA,:)%rel==1 .AND. dat(LA+1:UA+1,:)%rel>=0  .AND. dat(LA:UA,:)%typ==i .AND. move(LA:UA,:)==0), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
                write(name(im),'("getdiv stay LA:UA,typ ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do            





            do g=minsex,maxsex
                do i=1,ntypp
                    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%typ==i.and.dat(MNA:MXA,:)%sexr==g), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
                    write(name(im),'("mar by typ edned sex ",2i4)') i,g
                    weights(im)=0.0_dp
                    im=im+1
                end do         
                do i=1,ntypp
                    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%typ==i .AND. dat(MNA:MXA,:)%edr==1.and.dat(MNA:MXA,:)%sexr==g), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
                    write(name(im),'("mar by typ ned sex ",2i4)') i,g
                    weights(im)=0.0_dp
                    im=im+1
                end do         
                do i=1,ntypp
                    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%typ==i .AND. dat(MNA:MXA,:)%edr==2.and.dat(MNA:MXA,:)%sexr==g), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
                    write(name(im),'("mar by typ ed sex ",2i4)') i,g
                    weights(im)=0.0_dp
                    im=im+1
                end do                     
            end do
            do i=1,ntypp
                call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%typ==i), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
                write(name(im),'("mar by typ edned BTH sex ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do         
            do i=1,ntypp
                call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%typ==i .AND. dat(MNA:MXA,:)%edr==1), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
                write(name(im),'("mar by typ ned BTH sex ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do         
            do i=1,ntypp
                call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%typ==i .AND. dat(MNA:MXA,:)%edr==2), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
                write(name(im),'("mar by typ ed BTH sex ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do                                 



            do g=minsex,maxsex
                do i=1,ntypp
                    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==1 .and. dat(MNA:MXA,:)%sexr==g), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | mar edned sex ",2i4)') i,g
                    weights(im)=0.0_dp
                    im=im+1
                end do 
                do i=1,ntypp
                    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==1 .AND. dat(MNA:MXA,:)%edr==1.and. dat(MNA:MXA,:)%sexr==g), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | mar ned sex ",2i4)') i,g
                    weights(im)=0.0_dp
                    im=im+1
                end do 
                do i=1,ntypp
                    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==1 .AND. dat(MNA:MXA,:)%edr==2.and. dat(MNA:MXA,:)%sexr==g), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | mar ed sex ",2i4)') i,g
                    weights(im)=0.0_dp
                    im=im+1
                end do 
            end do 
            do i=1,ntypp
                call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==1), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                write(name(im),'("typ | mar edned BTH sex ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do 
            do i=1,ntypp
                call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==1 .AND. dat(MNA:MXA,:)%edr==1), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                write(name(im),'("typ | mar ned BTH sex ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do 
            do i=1,ntypp
                call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==1 .AND. dat(MNA:MXA,:)%edr==2), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                write(name(im),'("typ | mar ed BTH sex ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do 
            



            

            do g=minsex,maxsex
                do i=1,ntypp
                    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==0 .and. dat(MNA:MXA,:)%sexr==g), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | sin edned sex ",2i4)') i,g
                    weights(im)=0.0_dp
                    im=im+1
                end do 
                do i=1,ntypp
                    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==0 .AND. dat(MNA:MXA,:)%edr==1.and. dat(MNA:MXA,:)%sexr==g), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | sin ned sex ",2i4)') i,g
                    weights(im)=0.0_dp
                    im=im+1
                end do 
                do i=1,ntypp
                    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==0 .AND. dat(MNA:MXA,:)%edr==2.and. dat(MNA:MXA,:)%sexr==g), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | sin ed sex ",2i4)') i,g
                    weights(im)=0.0_dp
                    im=im+1
                end do 
            end do 
            do i=1,ntypp
                call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==0), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                write(name(im),'("typ | sin edned BTH sex ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do 
            do i=1,ntypp
                call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==0 .AND. dat(MNA:MXA,:)%edr==1), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                write(name(im),'("typ | sin ned BTH sex ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do 
            do i=1,ntypp
                call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==0 .AND. dat(MNA:MXA,:)%edr==2), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                write(name(im),'("typ | sin ed BTH sex ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do 
            



            !***********************************************************************************************************************************************************************************************************************************
            if (momdisplay) then 
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do i=1,ntypp
                    do ia=MNA,MXA,10
                    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel>=0 .AND. dat(ia,:)%rel==0 .AND. dat(ia,:)%edr==1), d1*one(dat(ia,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | sin ia ned ",2i4)') i,ia
                    weights(im)=0.0_dp
                    momwhich(im,1)=900 ; momwhich(im,2)=ia ; momwhich(im,3)=99 ; momwhich(im,4)=99 ; momwhich(im,5)=i ; momwhich(im,6)=1 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99  
                    im=im+1
                    end do 
                end do 

                do i=1,ntypp
                    do ia=MNA,MXA,10
                    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel>=0 .AND. dat(ia,:)%rel==0 .AND. dat(ia,:)%edr==2), d1*one(dat(ia,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | sin ia ed ",2i4)') i,ia
                    weights(im)=0.0_dp
                    momwhich(im,1)=900 ; momwhich(im,2)=ia ; momwhich(im,3)=99 ; momwhich(im,4)=99 ; momwhich(im,5)=i ; momwhich(im,6)=2 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99  
                    im=im+1
                    end do 
                end do 


                do i=1,ntypp
                    do ia=MNA,MXA,10
                    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel>=0 .AND. dat(ia,:)%rel==0 ), d1*one(dat(ia,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | sin ia EDNED ",2i4)') i,ia
                    weights(im)=0.0_dp
                    momwhich(im,1)=900 ; momwhich(im,2)=ia ; momwhich(im,3)=99 ; momwhich(im,4)=99 ; momwhich(im,5)=i ; momwhich(im,6)=6 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99  
                    im=im+1
                    end do 
                end do 


                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do i=1,ntypp
                    do ia=MNA,MXA,10
                    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel>=0 .AND. dat(ia,:)%rel==1 .AND. dat(ia,:)%edr==1), d1*one(dat(ia,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | mar ia ned ",2i4)') i,ia
                    weights(im)=0.0_dp
                    momwhich(im,1)=901 ; momwhich(im,2)=ia ; momwhich(im,3)=99 ; momwhich(im,4)=99 ; momwhich(im,5)=i ; momwhich(im,6)=1 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99  
                    im=im+1
                end do 
                end do 

                do i=1,ntypp
                    do ia=MNA,MXA,10
                    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel>=0 .AND. dat(ia,:)%rel==1 .AND. dat(ia,:)%edr==2), d1*one(dat(ia,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | mar ia ed ",2i4)') i,ia
                    weights(im)=0.0_dp
                    momwhich(im,1)=901 ; momwhich(im,2)=ia ; momwhich(im,3)=99 ; momwhich(im,4)=99 ; momwhich(im,5)=i ; momwhich(im,6)=2 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99  
                    im=im+1
                    end do 
                end do 


                do i=1,ntypp
                    do ia=MNA,MXA,10
                    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel>=0 .AND. dat(ia,:)%rel==1 ), d1*one(dat(ia,:)%typ==i),mom,cnt,var)
                    write(name(im),'("typ | mar ia EDNED ",2i4)') i,ia
                    weights(im)=0.0_dp
                    momwhich(im,1)=901 ; momwhich(im,2)=ia ; momwhich(im,3)=99 ; momwhich(im,4)=99 ; momwhich(im,5)=i ; momwhich(im,6)=6 ; momwhich(im,7)=99 ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99  
                    im=im+1
                    end do 
                end do 


                !headloc(ihead)=im
                !headstr(ihead)='getmar (ned) by typ and getdiv by type'
                !ihead=ihead+1
                !do ia=mna,mxad,10 !ahu030622  changed from maxai-1 to mxad
                !    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0 .AND. dat(ia,:)%edr==1), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
                !    write(name(im),'("getmar (ned) by ia ",i4)') ia
                !    weights(im)=0.0_dp
                !    im=im+1
                !end do 
                !do ia=mna,mxad,10
                !    do i=1,ntypp
                !        call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0  .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%typ==i), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
                !        write(name(im),'("getmar (ned) by ia,typ ",2i4)') ia,i
                !        weights(im)=0.0_dp
                !        im=im+1
                !    end do      
                !end do
                !
                !headloc(ihead)=im
                !headstr(ihead)='getmar (ed) by typ and getdiv by type'
                !ihead=ihead+1
                !do ia=mna,mxad,10 !ahu030622  changed from maxai-1 to mxad
                !    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0 .AND. dat(ia,:)%edr==2), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
                !    write(name(im),'("getmar (ed) by ia ",i4)') ia
                !    weights(im)=0.0_dp
                !    im=im+1
                !end do 
                !do ia=mna,mxad,10
                !    do i=1,ntypp
                !        call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0  .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%typ==i), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
                !        write(name(im),'("getmar (ed) by ia,typ ",2i4)') ia,i
                !        weights(im)=0.0_dp
                !        im=im+1
                !    end do      
                !end do 
                !do ia=mna,mxad,10 !ahu030622 changed from maxai-1 to mxad    
                !    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==1 .AND. dat(ia+1,:)%rel>=0 ), d1*one(dat(ia+1,:)%rel==0),mom,cnt,var)
                !    write(name(im),'("getdiv by ia ",i4)') ia
                !    weights(im)=0.0_dp
                !    im=im+1
                !    do i=1,ntypp
                !        call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==1 .AND. dat(ia+1,:)%rel>=0  .AND. dat(ia,:)%typ==i), d1*one(dat(ia+1,:)%rel==0),mom,cnt,var)
                !        write(name(im),'("getdiv by ia,typ ",2i4)') ia,i
                !        weights(im)=0.0_dp
                !        im=im+1
                !    end do            
                !end do 
            end if !momdisplay
            !***********************************************************************************************************************************************************************************************************************************


            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            !do g=minsex,maxsex
            !    LA=MNAD ; UA=MXAD
            !    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==0 )			                    
            !    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. dat(LA:UA,:)%edr==1 ), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
            !    write(name(im),'("getmar NED SEX ",i4)') g
            !    momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
            !    weights(im)=0.0_dp
            !    im=im+1
            !    do i=1,ntypp            
            !        call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. dat(LA:UA,:)%edr==1 .AND. dat(LA:UA,:)%typ==i), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
            !        write(name(im),'("getmar NED SEX TYP ",2i4)') g,i
            !        momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
            !        weights(im)=0.0_dp
            !        im=im+1
            !    end do            
            !end do 
!
            !do g=minsex,maxsex
            !    LA=MNAD ; UA=MXAD
            !    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==0 )			                    
            !    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. dat(LA:UA,:)%edr==2 ), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
            !    write(name(im),'("getmar ED SEX ",i4)') g
            !    momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
            !    weights(im)=0.0_dp
            !    im=im+1
            !    do i=1,ntypp            
            !        call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. dat(LA:UA,:)%edr==2 .AND. dat(LA:UA,:)%typ==i), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
            !        write(name(im),'("getmar ED SEX TYP ",2i4)') g,i
            !        momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
            !        weights(im)=0.0_dp
            !        im=im+1
            !    end do            
            !end do         
!
            !do g=minsex,maxsex
            !    LA=MNAD ; UA=MXAD
            !    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==0 )			                    
            !    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0  ), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
            !    write(name(im),'("getmar EDNED SEX ",i4)') g
            !    momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
            !    weights(im)=0.0_dp
            !    im=im+1
            !    do i=1,ntypp            
            !        call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. dat(LA:UA,:)%typ==i), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
            !        write(name(im),'("getmar EDNED SEX TYP ",2i4)') g,i
            !        momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
            !        weights(im)=0.0_dp
            !        im=im+1
            !    end do            
            !end do         



            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            !LA=MNAD ; UA=MXAD
            !cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==0 )			                    
            !call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. dat(LA:UA,:)%edr==1 ), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
            !write(name(im),'("getmar NED ALL SEX ")') 
            !momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=6 ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
            !weights(im)=0.0_dp
            !im=im+1
            !do i=1,ntypp            
            !    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. dat(LA:UA,:)%edr==1 .AND. dat(LA:UA,:)%typ==i), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
            !    write(name(im),'("getmar NED ALL SEX TYP ",i4)') i
            !    momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=6 ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
            !    weights(im)=0.0_dp
            !    im=im+1
            !end do            
!
            !LA=MNAD ; UA=MXAD
            !cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==0 )			                    
            !call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. dat(LA:UA,:)%edr==2 ), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
            !write(name(im),'("getmar ED ALL SEX ")') 
            !momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=6 ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
            !weights(im)=0.0_dp
            !im=im+1
            !do i=1,ntypp            
            !    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. dat(LA:UA,:)%edr==2 .AND. dat(LA:UA,:)%typ==i), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
            !    write(name(im),'("getmar ED ALL SEX TYP ",i4)') i
            !    momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=6 ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
            !    weights(im)=0.0_dp
            !    im=im+1
            !end do            
   
            LA=MNAD ; UA=MXAD
            cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==0 )			                    
            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0  ), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
            write(name(im),'("getmar EDNED ALL SEX ")') 
            momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=6 ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
            weights(im)=0.0_dp
            im=im+1
            do i=1,ntypp            
                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. dat(LA:UA,:)%typ==i), d1*one(dat(LA+1:UA+1,:)%rel==1),mom,cnt,var)   
                write(name(im),'("getmar EDNED ALL SEX TYP ",i4)') i
                momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=6 ; momwhich(im,4)=0  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                weights(im)=0.0_dp
                im=im+1
            end do            
            



            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            do g=minsex,maxsex
                LA=MNAD ; UA=MXAD
                cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==1 )			                    
                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==0), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
                write(name(im),'("getdiv stay SEX ",i4)') g
                momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=1  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=99; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                weights(im)=0.0_dp
                im=im+1
                do i=1,ntypp            
                    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==0 .AND. dat(LA:UA,:)%typ==i ), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
                    write(name(im),'("getdiv stay SEX TYP ",2i4)') g,i
                    momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=1  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=99; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                end do            
            end do 

            do g=minsex,maxsex
                LA=MNAD ; UA=MXAD
                cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==1 )			                    
                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==1), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
                write(name(im),'("getdiv move SEX ",i4)') g
                momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=1  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=99; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                weights(im)=0.0_dp
                im=im+1
                do i=1,ntypp            
                    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%typ==i ), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
                    write(name(im),'("getdiv move SEX TYP ",2i4)') g,i
                    momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=1  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=99; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                end do            
            end do 




            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            LA=MNAD ; UA=MXAD
            cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co   .and. dat(LA:UA,:)%rel==1 )			                    
            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==0), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv STAY ALL SEX ")') 
            momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=6 ; momwhich(im,4)=1  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=99; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
            weights(im)=0.0_dp
            im=im+1
            do i=1,ntypp            
                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==0 .AND. dat(LA:UA,:)%typ==i ), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
                write(name(im),'("getdiv STAY ALL SEX TYP ",i4)') i
                momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=6 ; momwhich(im,4)=1  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=99; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                weights(im)=0.0_dp
                im=im+1
            end do            

            LA=MNAD ; UA=MXAD
            cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co   .and. dat(LA:UA,:)%rel==1 )			                    
            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==1), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("getdiv move ALL SEX ")') 
            momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=6 ; momwhich(im,4)=1  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=99; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
            weights(im)=0.0_dp
            im=im+1
            do i=1,ntypp            
                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA+1:UA+1,:)%rel>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%typ==i ), d1*one(dat(LA+1:UA+1,:)%rel==0),mom,cnt,var)
                write(name(im),'("getdiv move ALL SEX TYP ",i4)') i
                momwhich(im,1)=700 ; momwhich(im,2)=99 ; momwhich(im,3)=6 ; momwhich(im,4)=1  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=99; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                weights(im)=0.0_dp
                im=im+1
            end do            




            do ia=MNAD,40,5
                call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel>=0  .AND. dat(ia,:)%edr==1), d1*one(dat(ia,:)%rel==1),mom,cnt,var)
                write(name(im),'("mar by age ned",i4)') ia
                weights(im)=0.0_dp
                im=im+1
            end do         
            do ia=MNAD,40,5
                call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel>=0  .AND. dat(ia,:)%edr==2), d1*one(dat(ia,:)%rel==1),mom,cnt,var)
                write(name(im),'("mar by age ed",i4)') ia
                weights(im)=0.0_dp
                im=im+1
            end do         

            !***********************************************************************************************************************************************************************************************************************************
            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            if (momdisplay) then 
                g=MALES
                do j=maxrelo,0,-1
                    headloc(ihead)=im
                    if (g==1.and.j==0) headstr(ihead)='single men move by age'
                    if (g==1.and.j==1) headstr(ihead)='married men move by age'
                    ihead=ihead+1                       
                    do ia=MNA,40,4
                        cosexrel(ia,:)= (dat(ia,:)%co==co .and. dat(ia,:)%sexr==g  .and. dat(ia,:)%rel==j .AND.norelchg(ia,:)==1)			                    		                    
                        call condmom(im,( cosexrel(ia,:) .AND.  move(ia,:)>=0 ),   d1* move(ia,:) ,mom,cnt,var)	
                            write(name(im),'("move M by rel,age  ",2I4)') j,ia
                        momwhich(im,1)=590 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                    end do 
                    do i=1,ntypp
                        do ia=MNA,40,4
                            cosexrel(ia,:)= (dat(ia,:)%co==co .and. dat(ia,:)%sexr==g  .and. dat(ia,:)%rel==j .AND.norelchg(ia,:)==1)			                    		                    
                            call condmom(im,( cosexrel(ia,:) .AND.  move(ia,:)>=0  .AND. dat(ia,:)%typ==i ),   d1* move(ia,:) ,mom,cnt,var)	
                            write(name(im),'("move M rel,typ,age ",3I4)') j,i,ia
                            momwhich(im,1)=590 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                            weights(im)=0.0_dp 
                            im=im+1
                        end do !age 
                    end do  !type
                end do !j rel


                g=MALES
                do j=maxrelo,0,-1
                    headloc(ihead)=im
                    if (g==1.and.j==0) headstr(ihead)='single men NED move by age'
                    if (g==1.and.j==1) headstr(ihead)='married men NED move by age'
                    ihead=ihead+1                       
                    do ia=MNA,40,4
                        cosexrel(ia,:)= (dat(ia,:)%co==co .and. dat(ia,:)%sexr==g  .and. dat(ia,:)%edr==1  .and. dat(ia,:)%rel==j .AND.norelchg(ia,:)==1)			                    		                    
                        call condmom(im,( cosexrel(ia,:) .AND.  move(ia,:)>=0 ),   d1* move(ia,:) ,mom,cnt,var)	
                            write(name(im),'("move M NED by rel,age  ",2I4)') j,ia
                        momwhich(im,1)=590 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                    end do 
                    do i=1,ntypp
                        do ia=MNA,40,2    
                            cosexrel(ia,:)= (dat(ia,:)%co==co .and. dat(ia,:)%sexr==g .and. dat(ia,:)%edr==1  .and. dat(ia,:)%rel==j .AND.norelchg(ia,:)==1)			                    		                    
                            call condmom(im,( cosexrel(ia,:) .AND.  move(ia,:)>=0  .AND. dat(ia,:)%typ==i ),   d1* move(ia,:) ,mom,cnt,var)	
                            write(name(im),'("move M NED rel,typ,age ",3I4)') j,i,ia
                            momwhich(im,1)=590 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                            weights(im)=0.0_dp 
                            im=im+1
                        end do !age 
                    end do  !type
                end do !j rel

                g=MALES
                do j=maxrelo,0,-1
                    headloc(ihead)=im
                    if (g==1.and.j==0) headstr(ihead)='single men ED move by age'
                    if (g==1.and.j==1) headstr(ihead)='married men ED move by age'
                    ihead=ihead+1                       
                    do ia=MNA,40,4
                        cosexrel(ia,:)= (dat(ia,:)%co==co .and. dat(ia,:)%sexr==g  .and. dat(ia,:)%edr==2  .and. dat(ia,:)%rel==j .AND.norelchg(ia,:)==1)			                    		                    
                        call condmom(im,( cosexrel(ia,:) .AND.  move(ia,:)>=0 ),   d1* move(ia,:) ,mom,cnt,var)	
                        write(name(im),'("move M ED by rel,age ",2I4)') j,ia
                        momwhich(im,1)=590 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                    end do 
                    do i=1,ntypp
                        do ia=MNA,40,2    
                            cosexrel(ia,:)= (dat(ia,:)%co==co .and. dat(ia,:)%sexr==g .and. dat(ia,:)%edr==2  .and. dat(ia,:)%rel==j .AND.norelchg(ia,:)==1)			                    		                    
                            call condmom(im,( cosexrel(ia,:) .AND.  move(ia,:)>=0  .AND. dat(ia,:)%typ==i ),   d1* move(ia,:) ,mom,cnt,var)	
                            write(name(im),'("move M ED rel,typ,age ",3I4)') j,i,ia
                            momwhich(im,1)=590 ; momwhich(im,2)=ia ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                            weights(im)=0.0_dp 
                            im=im+1
                        end do !age 
                    end do  !type
                end do !j rel


                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                g=MALES
                do j=maxrelo,0,-1
                    headloc(ihead)=im
                    if (g==1.and.j==0) headstr(ihead)='single men emptrans by age'
                    if (g==1.and.j==1) headstr(ihead)='married men emptrans by age'
                    ihead=ihead+1                       
                    do ia=MNA,40,2
                        cosexrel(ia,:)= (dat(ia,:)%co==co .and. dat(ia,:)%sexr==g  .and. dat(ia,:)%rel==j .and. norelchg(ia,:)==1 )			
                        call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==1 ),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | u move M by rel,age ",2i4)') j,ia 
                        weights(im)=0.0_dp
                        im=im+1 
                    end do 
                    do i=1,ntypp
                        do ia=MNA,40,2
                            cosexrel(ia,:)= (dat(ia,:)%co==co .and. dat(ia,:)%sexr==g  .and. dat(ia,:)%rel==j .and. norelchg(ia,:)==1 )			
                            call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%typ==i .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==1 ),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | u move M rel,type,age ",3i4)') j,i,ia
                            weights(im)=0.0_dp
                            im=im+1
                        end do  
                    end do 
                end do 

            end if !momdisplay
            !***********************************************************************************************************************************************************************************************************************************




            !****************************************************************************************************************************************************************
            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE BY SEX,REL,TYPE, 17-49'
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1)			                    			                    
                    call condmom(im,( cosexrel(LA:UA,:) .AND.  move(LA:UA,:)>=0 ),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move by sex,rel     17-49 ",2I4)') g,j
                    momwhich(im,1)=600 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND.  move(LA:UA,:)>=0  .AND. dat(LA:UA,:)%typ==i ),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move by sex,rel,typ 17-49 ",3I4)') g,j,i
                        momwhich(im,1)=600 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99               ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !j rel
            end do !g sex 

            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE BY SEX,REL,ED,TYPE, 17-49'
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                do edo=1,2
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%edr==edo )			                    
                    call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 ),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move BY sex,rel,ed     17-49 ",3I4)') g,j,edo
                    momwhich(im,1)=600 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                      ; momwhich(im,5)=6 ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND. dat(LA:UA,:)%typ==i),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move BY sex,rel,ed,typ 17-49 ",4I4)') g,j,edo,i
                        momwhich(im,1)=600 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !edo
                end do !j rel
            end do !g sex 
            !****************************************************************************************************************************************************************
            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE HOME BY SEX,REL,TYP, 17-49'
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1)
                    call condmom(im,( cosexrel(LA:UA,:) .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%l==dat(LA+1:UA+1,:)%hme  ),mom,cnt,var)		
                        write(name(im),'("movehme by sex,rel    17-49 ",2I4)') g,j
                    momwhich(im,1)=605 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%typ==i ),   d1*one( dat(LA+1:UA+1,:)%l==dat(LA+1:UA+1,:)%hme  ),mom,cnt,var)		
                        write(name(im),'("movehme by sex,rel,typ 17-49 ",3I4)') g,j,i
                        momwhich(im,1)=605 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99               ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !j rel
            end do !g sex 
            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE HOME BY SEX,ED,REL,TYP, 17-49'
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                do edo=1,2
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1.AND.dat(LA:UA,:)%edr==edo )
                    call condmom(im,( cosexrel(LA:UA,:) .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%l==dat(LA+1:UA+1,:)%hme  ),mom,cnt,var)		
                        write(name(im),'("movehme by sex,rel,ed     17-49 ",3I4)') g,j,edo
                    momwhich(im,1)=605 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%typ==i ),   d1*one( dat(LA+1:UA+1,:)%l==dat(LA+1:UA+1,:)%hme  ),mom,cnt,var)		
                        write(name(im),'("movehme by sex,rel,ed,typ 17-49 ",4I4)') g,j,edo,i
                        momwhich(im,1)=605 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99               ; momwhich(im,5)=i ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !edo
                end do !j rel
            end do !g sex 


            !****************************************************************************************************************************************************************
            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE AT HOME -BY SEX,REL,TYPE 17-49'
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1)			                    
                    call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-at hme-by sex,rel     17-49 ",2I4)') g,j
                    momwhich(im,1)=610 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                      ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme .AND. dat(LA:UA,:)%typ==i),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-at hme-by sex,rel,typ 17-49 ",3I4)') g,j,i
                        momwhich(im,1)=610 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !j rel
            end do !g sex 

            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE AT HOME -BY SEX,REL,ED,TYPE 17-49'
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                do edo=1,2
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1 .and. dat(LA:UA,:)%edr==edo)			                    
                    call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-at hme-by sex,rel,ed     17-49 ",3I4)') g,j,edo
                    momwhich(im,1)=610 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                      ; momwhich(im,5)=6 ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme .AND. dat(LA:UA,:)%typ==i),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-at hme-by sex,rel,ed,typ 17-49 ",4I4)') g,j,edo,i
                        momwhich(im,1)=610 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !edo
                end do !j rel
            end do !g sex 
            !****************************************************************************************************************************************************************
            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE NOT HOME -BY SEX,REL,TYPE 17-49'
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1)			                    
                    call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-not hme-by sex,rel    17-49 ",2I4)') g,j
                    momwhich(im,1)=612 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                      ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme .AND. dat(LA:UA,:)%typ==i),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-not hme-sex,rel,typ   17-49 ",3I4)') g,j,i
                        momwhich(im,1)=612 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !j rel
            end do !g sex 

            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE NOT HOME -BY SEX,REL,ED,TYPE 17-49'
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                do edo=1,2
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co .and. dat(LA:UA,:)%sexr==g  .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1 .and. dat(LA:UA,:)%edr==edo)			                    
                    call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-nt hme-by sex,rel,ed     17-49 ",3I4)') g,j,edo
                    momwhich(im,1)=612 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                      ; momwhich(im,5)=6 ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme .AND. dat(LA:UA,:)%typ==i),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-nt hme-by sex,rel,ed,typ 17-49 ",4I4)') g,j,edo,i
                        momwhich(im,1)=612 ; momwhich(im,2)=1749 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !edo
                end do !j rel
            end do !g sex 
            !****************************************************************************************************************************************************************



            !****************************************************************************************************************************************************************
            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE AT HOME -BTH SEX,REL,TYPE 17-49'
            ihead=ihead+1
                do j=maxrelo,0,-1
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co   .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1)			                    
                    call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-at hme-BTH sex,rel     17-49 ",I4)') j
                    momwhich(im,1)=610 ; momwhich(im,2)=1749 ; momwhich(im,3)=6 ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                      ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme .AND. dat(LA:UA,:)%typ==i),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-at hme-BTH sex,rel,typ 17-49 ",2I4)') j,i
                        momwhich(im,1)=610 ; momwhich(im,2)=1749 ; momwhich(im,3)=6 ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !j rel

            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE AT HOME -BTH SEX,REL,ED,TYPE 17-49'
            ihead=ihead+1
                do j=maxrelo,0,-1
                do edo=1,2
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co   .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1 .and. dat(LA:UA,:)%edr==edo)			                    
                    call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-at hme-BTH sex,rel,ed     17-49 ",2I4)') j,edo
                    momwhich(im,1)=610 ; momwhich(im,2)=1749 ; momwhich(im,3)=6 ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                      ; momwhich(im,5)=6 ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l==dat(LA:UA,:)%hme .AND. dat(LA:UA,:)%typ==i),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-at hme-BTH sex,rel,ed,typ 17-49 ",3I4)') j,edo,i
                        momwhich(im,1)=610 ; momwhich(im,2)=1749 ; momwhich(im,3)=6 ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !edo
                end do !j rel
            !****************************************************************************************************************************************************************
            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE NOT HOME -BTH SEX,REL,TYPE 17-49'
            ihead=ihead+1
                do j=maxrelo,0,-1
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co   .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1)			                    
                    call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-not hme-BTHsex,rel    17-49 ",I4)') j
                    momwhich(im,1)=612 ; momwhich(im,2)=1749 ; momwhich(im,3)=6 ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                      ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme .AND. dat(LA:UA,:)%typ==i),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-not hme-BTH sex,rel,typ   17-49 ",2I4)') j,i
                        momwhich(im,1)=612 ; momwhich(im,2)=1749 ; momwhich(im,3)=6 ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !j rel
            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            headloc(ihead)=im
            headstr(ihead)='MOVE NOT HOME -BTH SEX,REL,ED,TYPE 17-49'
            ihead=ihead+1
                do j=maxrelo,0,-1
                do edo=1,2
                    LA=MNAD ; UA=MXAD
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co   .and. dat(LA:UA,:)%rel==j .AND.norelchg(LA:UA,:)==1 .and. dat(LA:UA,:)%edr==edo)			                    
                    call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-nt hme-BTH sex,rel,ed     17-49 ",2I4)') j,edo
                    momwhich(im,1)=612 ; momwhich(im,2)=1749 ; momwhich(im,3)=6 ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                      ; momwhich(im,5)=6 ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        call condmom(im,( cosexrel(LA:UA,:) .AND.dat(LA:UA,:)%l>=0..AND.move(LA:UA,:)>=0 .AND.dat(LA:UA,:)%l/=dat(LA:UA,:)%hme .AND. dat(LA:UA,:)%typ==i),   d1* move(LA:UA,:) ,mom,cnt,var)	
                        write(name(im),'("move-nt hme-BTH sex,rel,ed,typ 17-49 ",3I4)') j,edo,i
                        momwhich(im,1)=612 ; momwhich(im,2)=1749 ; momwhich(im,3)=6 ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=edo; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do !edo
                end do !j rel
            !****************************************************************************************************************************************************************








            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            do g=minsex,maxsex
                do edo=1,2
                    cosexedtypid(:)            = (inito(:)%co==co .and. inito(:)%sexr==g  .and. inito(:)%edr==edo  )						
                    headloc(ihead)=im
                    if (co==1.and.g==1.and.edo==1) headstr(ihead)='LIFETIME INC men NED '
                    if (co==1.and.g==2.and.edo==1) headstr(ihead)='LIFETIME INC fem NED '
                    if (co==1.and.g==1.and.edo==2) headstr(ihead)='LIFETIME INC men ED '
                    if (co==1.and.g==2.and.edo==2) headstr(ihead)='LIFETIME INC fem ED '
                    ihead=ihead+1
                    !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                    CALL condmom(im,( COSEXEDTYPID(:) .AND. lifetimeinc(:)>=0 ),d1*lifetimeinc(:),mom,cnt,var)
                    write(name(im),'("LIFETIME INC ",tr10)')	
                    momwhich(im,1)=55 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=99 ; momwhich(im,5)=6 ; momwhich(im,6)=edo ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                    weights(im)=0.0_dp 
                    im=im+1
                    do i=1,ntypp
                        CALL condmom(im,( COSEXEDTYPID(:) .AND. lifetimeinc(:)>=0 .and. inito(:)%typ==i 	),d1*lifetimeinc(:),mom,cnt,var)
                        write(name(im),'("LIFETIME INC by type ",i4)')	i 
                        momwhich(im,1)=55 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=99 ; momwhich(im,5)=i ; momwhich(im,6)=edo ; momwhich(im,7)=99  ;  momwhich(im,8)=99 ; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                        weights(im)=0.0_dp 
                        im=im+1
                    end do 
                end do 
            end do 




                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                            headloc(ihead)=im
                            if (g==1.and.j==0) headstr(ihead)='single men wage by type '
                            if (g==1.and.j==1) headstr(ihead)='married men wage by type '
                            if (g==2.and.j==0) headstr(ihead)='single fem wage by type '
                            if (g==2.and.j==1) headstr(ihead)='married fem wage by type '
                            ihead=ihead+1     
                            call condmom(im,( cosexrel(LA:UA,:)  .and. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%hhr==1),    d1*dat(LA:UA,:)%logwr ,mom,cnt,var)		
                            write(name(im),'("wages ALL types")') 
                            weights(im)=0.
                            im=im+1 
                            do i=1,ntypp
                                LA=MNA ; UA=MXAD
                                call condmom(im,( cosexrel(LA:UA,:).AND. dat(LA:UA,:)%typ==i  .and. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%hhr==1),    d1*dat(LA:UA,:)%logwr ,mom,cnt,var)		
                                write(name(im),'("wages ", i4)') i  
                                weights(im)=0.
                                im=im+1 
                            end do 
                        end do 
                end do 

                


                LA=MNA ; UA=MXAD
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men wdif stay NED by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif stay NED by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem wdif stay NED by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif stay NED by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==1  ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | stay NED ALL types ",tr2)')  
                        momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==1 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                            write(name(im),'("wdif | stay NED ",tr2,i4)')  i
                            momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men wdif stay ED by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif stay ED by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem wdif stay ED by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif stay ED by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==2  ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | stay ED ALL types ",tr2)')  
                        momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==2 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                            write(name(im),'("wdif | stay ED ",tr2,i4)')  i
                            momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men wdif stay by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif stay by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem wdif stay by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif stay by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | stay ALL types ",tr2)')  
                        momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                            write(name(im),'("wdif | stay ",tr2,i4)')  i
                            momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         



                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men wdif  move NED by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif move NED by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem wdif  move NED by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif move NED by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==1 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | move NED ALL types ",tr2)')  
                        momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==1 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                            write(name(im),'("wdif | move NED ",tr2,i4)')  i
                            momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         

                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        LA=MNAD ; UA=MXAD
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men wdif  move ED by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif move ED by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem wdif  move ED by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif move ED by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==2 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | move ED ALL types ",tr2)')  
                        momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==2 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                            write(name(im),'("wdif | move ED ",tr2,i4)')  i
                            momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         
    
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        LA=MNAD ; UA=MXAD
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men wdif  move by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif move by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem wdif  move by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif move by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | move ALL types ",tr2)')  
                        momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                            write(name(im),'("wdif | move ",tr2,i4)')  i
                            momwhich(im,1)=160 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         
    
               
    





                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        LA=MNAD ; UA=MXAD
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men wdif<0 stay NED by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif<0 stay NED by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem wdif<0 stay NED by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif<0 stay NED by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==1 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                        write(name(im),'("wdif<0 | stay NED ALL types ",tr2)')  
                        momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==1 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                            write(name(im),'("wdif<0 | stay NED ",tr2,i4)')  i
                            momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         

                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        LA=MNAD ; UA=MXAD
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men wdif<0 stay ED by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif<0 stay ED by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem wdif<0 stay ED by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif<0 stay ED by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==2 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                        write(name(im),'("wdif<0 | stay ED ALL types ",tr2)')  
                        momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==2 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                            write(name(im),'("wdif<0 | stay ED ",tr2,i4)')  i
                            momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         
    
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        LA=MNAD ; UA=MXAD
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men wdif<0 stay by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif<0 stay by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem wdif<0 stay by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif<0 stay by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                        write(name(im),'("wdif<0 | stay ALL types ",tr2)')  
                        momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                            write(name(im),'("wdif<0 | stay ",tr2,i4)')  i
                            momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=0  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         













   
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        LA=MNAD ; UA=MXAD
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men  wdif<0  move NED by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif<0 move NED by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem  wdif<0  move NED by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif<0 move NED by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==1 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                        write(name(im),'("wdif<0 | move NED ALL types ",tr2)')  
                        momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==1 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                            write(name(im),'("wdif<0 | move NED ",tr2,i4)')  i
                            momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         

                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        LA=MNAD ; UA=MXAD
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men  wdif<0  move ED by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif<0 move ED by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem  wdif<0  move ED by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif<0 move ED by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==2 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                        write(name(im),'("wdif<0 | move ED ALL types ",tr2)')  
                        momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==2 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                            write(name(im),'("wdif<0 | move ED ",tr2,i4)')  i
                            momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         
    
                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        LA=MNAD ; UA=MXAD
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men  wdif<0  move by type '
                        if (g==1.and.j==1) headstr(ihead)='married men wdif<0 move by type '
                        if (g==2.and.j==0) headstr(ihead)='single fem  wdif<0  move by type '
                        if (g==2.and.j==1) headstr(ihead)='married fem wdif<0 move by type '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                        write(name(im),'("wdif<0 | move ALL types ",tr2)')  
                        momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                            write(name(im),'("wdif<0 | move ",tr2,i4)')  i
                            momwhich(im,1)=165 ; momwhich(im,2)=99 ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         


                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do iesp=1,2
                    do ier=1,2
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==1  .and. dat(MNAD:MXAD,:)%rel==1 .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (ier==1.and.iesp==1) headstr(ihead)='wdif M by jointed NED NED '
                        if (ier==2.and.iesp==1) headstr(ihead)='wdif M by jointed  ED NED '
                        if (ier==1.and.iesp==2) headstr(ihead)='wdif M by jointed NED  ED '
                        if (ier==2.and.iesp==2) headstr(ihead)='wdif M by jointed  ED  ED '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | move by jointed ",tr2)')  
                        momwhich(im,1)=167 ; momwhich(im,2)=99 ; momwhich(im,3)=1 ; momwhich(im,4)=1  ; momwhich(im,7)=iesp ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=ier; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp  ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                            write(name(im),'("wdif | move by jointed typ ",tr2,i4)')  i
                            momwhich(im,1)=167 ; momwhich(im,2)=99 ; momwhich(im,3)=1 ; momwhich(im,4)=1  ; momwhich(im,7)=iesp ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=ier; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         


                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do iesp=1,2
                    do ier=1,2
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==1  .and. dat(MNAD:MXAD,:)%rel==1 .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (ier==1.and.iesp==1) headstr(ihead)='wdif<0 M by jointed NED NED '
                        if (ier==2.and.iesp==1) headstr(ihead)='wdif<0 M by jointed  ED NED '
                        if (ier==1.and.iesp==2) headstr(ihead)='wdif<0 M by jointed NED  ED '
                        if (ier==2.and.iesp==2) headstr(ihead)='wdif<0 M by jointed  ED  ED '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                        write(name(im),'("wdif<0 | move by jointed ",tr2)')  
                        momwhich(im,1)=169 ; momwhich(im,2)=99 ; momwhich(im,3)=1 ; momwhich(im,4)=1  ; momwhich(im,7)=iesp ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=ier; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp  ),    d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                            write(name(im),'("wdif<0 | move by jointed typ ",tr2,i4)')  i
                            momwhich(im,1)=169 ; momwhich(im,2)=99 ; momwhich(im,3)=1 ; momwhich(im,4)=1  ; momwhich(im,7)=iesp ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=ier; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do      





                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do ier=1,2
                    do iesp=1,2
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==2  .and. dat(MNAD:MXAD,:)%rel==1 .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (ier==1.and.iesp==1) headstr(ihead)='wdif F by jointed NED NED '
                        if (ier==1.and.iesp==2) headstr(ihead)='wdif F by jointed  ED NED '
                        if (ier==2.and.iesp==1) headstr(ihead)='wdif F by jointed NED  ED '
                        if (ier==2.and.iesp==2) headstr(ihead)='wdif F by jointed  ED  ED '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                        write(name(im),'("wdif | move by jointed ",tr2)')  
                        momwhich(im,1)=167 ; momwhich(im,2)=99 ; momwhich(im,3)=2 ; momwhich(im,4)=1  ; momwhich(im,7)=iesp ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=ier; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp  ),   d1*(wdif(LA:UA,:)) ,mom,cnt,var)		
                            write(name(im),'("wdif | move by jointed typ ",tr2,i4)')  i
                            momwhich(im,1)=167 ; momwhich(im,2)=99 ; momwhich(im,3)=2 ; momwhich(im,4)=1  ; momwhich(im,7)=iesp ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=ier; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do         


                !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do ier=1,2
                    do iesp=1,2
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==2  .and. dat(MNAD:MXAD,:)%rel==1 .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (ier==1.and.iesp==1) headstr(ihead)='wdif<0 F by jointed NED NED '
                        if (ier==1.and.iesp==2) headstr(ihead)='wdif<0 F by jointed  ED NED '
                        if (ier==2.and.iesp==1) headstr(ihead)='wdif<0 F by jointed NED  ED '
                        if (ier==2.and.iesp==2) headstr(ihead)='wdif<0 F by jointed  ED  ED '
                        ihead=ihead+1     
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp ),   d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                        write(name(im),'("wdif<0 | move by jointed ",tr2)')  
                        momwhich(im,1)=169 ; momwhich(im,2)=99 ; momwhich(im,3)=2 ; momwhich(im,4)=1  ; momwhich(im,7)=iesp ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=ier; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                        weights(im)=0.0_dp
                        im=im+1
                        do i=1,ntypp            
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp  ),    d1*one(wdif(LA:UA,:)<0) ,mom,cnt,var)		
                            write(name(im),'("wdif<0 | move by jointed typ ",tr2,i4)')  i
                            momwhich(im,1)=169 ; momwhich(im,2)=99 ; momwhich(im,3)=2 ; momwhich(im,4)=1  ; momwhich(im,7)=iesp ;  momwhich(im,8)=99                   ; momwhich(im,5)=i ; momwhich(im,6)=ier; momwhich(im,9)=99 ; momwhich(im,10)=1  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1
                        end do           
                    end do  
                end do      




            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        do ddd=1,2
                            LA=MNAD ; UA=MXAD
                            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                            headloc(ihead)=im
                            if (g==1.and.j==0) headstr(ihead)='single men emp overall by type 18-49 or 20-34'
                            if (g==1.and.j==1) headstr(ihead)='married men emp overall by type 18-49 or 20-34'
                            if (g==2.and.j==0) headstr(ihead)='single fem emp overall by type 18-49 or 20-34'
                            if (g==2.and.j==1) headstr(ihead)='married fem emp overall by type 18-49 or 20-34'
                            ihead=ihead+1                       
                            if (ddd==1) agerange=1849
                            if (ddd==2) agerange=2034 
                            if (ddd==1) then ;  LA=18 ; UA=49 ; end if 
                            if (ddd==2) then ;  LA=20 ; UA=34 ; end if 
    
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr>=0 ),   d1*one( dat(LA:UA,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("emp all ALL ")') 
                            momwhich(im,1)=3480 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1 
                            do i=1,ntypp
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr>=0 .AND. dat(LA:UA,:)%typ==i),   d1*one( dat(LA:UA,:)%hhr==1 ),mom,cnt,var)		
                                write(name(im),'("emp all ALL by type ",i4)') i  
                                momwhich(im,1)=3480 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99               ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                                weights(im)=0.0_dp
                                im=im+1 
                            end do 
                        end do 
                    end do 
                end do 

    
    



            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        do ddd=1,2
                            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 .and. dat(MNAD:MXAD,:)%EDR==1 )			
                            headloc(ihead)=im
                            if (g==1.and.j==0) headstr(ihead)='single men NED emp overall by type 18-49 or 20-34'
                            if (g==1.and.j==1) headstr(ihead)='married men NED emp overall by type 18-49 or 20-34'
                            if (g==2.and.j==0) headstr(ihead)='single fem NED emp overall by type 18-49 or 20-34'
                            if (g==2.and.j==1) headstr(ihead)='married fem NED emp overall by type 18-49 or 20-34'
                            ihead=ihead+1                       
                            if (ddd==1) agerange=1849
                            if (ddd==2) agerange=2034 
                            if (ddd==1) then ;  LA=18 ; UA=49 ; end if 
                            if (ddd==2) then ;  LA=20 ; UA=34 ; end if 
    
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr>=0 ),   d1*one( dat(LA:UA,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("emp NED all ALL ")') 
                            momwhich(im,1)=3480 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1 
                            do i=1,ntypp
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr>=0 .AND. dat(LA:UA,:)%typ==i),   d1*one( dat(LA:UA,:)%hhr==1 ),mom,cnt,var)		
                                write(name(im),'("emp NED all ALL by type ",i4)') i  
                                momwhich(im,1)=3480 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99               ; momwhich(im,5)=i ; momwhich(im,6)=1; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                                weights(im)=0.0_dp
                                im=im+1 
                            end do 
                        end do 
                    end do 
                end do 


            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                do g=minsex,maxsex
                    do j=maxrelo,0,-1
                        do ddd=1,2
                            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 .and. dat(MNAD:MXAD,:)%EDR==2 )			
                            headloc(ihead)=im
                            if (g==1.and.j==0) headstr(ihead)='single men ED emp overall by type 18-49 or 20-34'
                            if (g==1.and.j==1) headstr(ihead)='married men ED emp overall by type 18-49 or 20-34'
                            if (g==2.and.j==0) headstr(ihead)='single fem ED emp overall by type 18-49 or 20-34'
                            if (g==2.and.j==1) headstr(ihead)='married fem ED emp overall by type 18-49 or 20-34'
                            ihead=ihead+1                       
                            if (ddd==1) agerange=1849
                            if (ddd==2) agerange=2034 
                            if (ddd==1) then ;  LA=18 ; UA=49 ; end if 
                            if (ddd==2) then ;  LA=20 ; UA=34 ; end if 
    
                            call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr>=0 ),   d1*one( dat(LA:UA,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("emp ED all ALL ")') 
                            momwhich(im,1)=3480 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                            weights(im)=0.0_dp
                            im=im+1 
                            do i=1,ntypp
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr>=0 .AND. dat(LA:UA,:)%typ==i),   d1*one( dat(LA:UA,:)%hhr==1 ),mom,cnt,var)		
                                write(name(im),'("emp ED all ALL by type ",i4)') i  
                                momwhich(im,1)=3480 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99               ; momwhich(im,5)=i ; momwhich(im,6)=2; momwhich(im,9)=99 ; momwhich(im,10)=99  ; momwhich(im,11)=99
                                weights(im)=0.0_dp
                                im=im+1 
                            end do 
                        end do 
                    end do 
                end do 





                LA=MNA ; UA=MXAD
                headloc(ihead)=im
                headstr(ihead)='.|EE STAY ALL (MEN added 091023)'
                ihead=ihead+1                       
                cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 )			
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==1 ),mom,cnt,var)		
                write(name(im),'("ee|ee move ALL ")')   
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==2  ),mom,cnt,var)		
                write(name(im),'("eu|ee move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==3  ),mom,cnt,var)		
                write(name(im),'("ue|ee move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==4 ),mom,cnt,var)		
                write(name(im),'("uu|ee move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                
                
                headloc(ihead)=im
                headstr(ihead)='.|EU STAY ALL (MEN added 091023)'
                ihead=ihead+1                       
                cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 )			
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==1  ),mom,cnt,var)		
                write(name(im),'("ee|eu move ALL ")')   
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==2  ),mom,cnt,var)		
                write(name(im),'("eu|eu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==3  ),mom,cnt,var)		
                write(name(im),'("ue|eu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==4  ),mom,cnt,var)		
                write(name(im),'("uu|eu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                
                headloc(ihead)=im
                headstr(ihead)='.|UE STAY ALL (MEN added 091023)'
                ihead=ihead+1                       
                cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 )			
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==1  ),mom,cnt,var)		
                write(name(im),'("ee|ue move ALL ")')   
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==2  ),mom,cnt,var)		
                write(name(im),'("eu|ue move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==3  ),mom,cnt,var)		
                write(name(im),'("ue|ue move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==4 ),mom,cnt,var)		
                write(name(im),'("uu|ue move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                
                
                headloc(ihead)=im
                headstr(ihead)='.|UU STAY ALL (MEN added 091023)'
                ihead=ihead+1                       
                cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 )			
                call condmom(im,( cosexrel(LA:UA,:) .AND.  jointemp(LA:UA,:)==4  .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==1 ),mom,cnt,var)		
                write(name(im),'("ee|uu move ALL ")')   
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND.  jointemp(LA:UA,:)==4  .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==2 ),mom,cnt,var)		
                write(name(im),'("eu|uu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND.  jointemp(LA:UA,:)==4  .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==3 ),mom,cnt,var)		
                write(name(im),'("ue|uu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND.  jointemp(LA:UA,:)==4 .AND. move(LA:UA,:)==0 ),   d1*one( jointemp(LA+1:UA+1,:)==4 ),mom,cnt,var)		
                write(name(im),'("uu|uu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 



                

                LA=MNA ; UA=MXAD
                headloc(ihead)=im
                headstr(ihead)='.|EE MOVE ALL (MEN added 091023)'
                ihead=ihead+1                       
                cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 )			
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==1 ),mom,cnt,var)		
                write(name(im),'("ee|ee move ALL ")')   
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==2  ),mom,cnt,var)		
                write(name(im),'("eu|ee move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==3  ),mom,cnt,var)		
                write(name(im),'("ue|ee move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==4 ),mom,cnt,var)		
                write(name(im),'("uu|ee move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                
                
                headloc(ihead)=im
                headstr(ihead)='.|EU MOVE ALL (MEN added 091023)'
                ihead=ihead+1                       
                cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 )			
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==1  ),mom,cnt,var)		
                write(name(im),'("ee|eu move ALL ")')   
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==2  ),mom,cnt,var)		
                write(name(im),'("eu|eu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==3  ),mom,cnt,var)		
                write(name(im),'("ue|eu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==4  ),mom,cnt,var)		
                write(name(im),'("uu|eu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                
                headloc(ihead)=im
                headstr(ihead)='.|UE MOVE ALL (MEN added 091023)'
                ihead=ihead+1                       
                cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 )			
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==1  ),mom,cnt,var)		
                write(name(im),'("ee|ue move ALL ")')   
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==2  ),mom,cnt,var)		
                write(name(im),'("eu|ue move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==3  ),mom,cnt,var)		
                write(name(im),'("ue|ue move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==4 ),mom,cnt,var)		
                write(name(im),'("uu|ue move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                
                
                headloc(ihead)=im
                headstr(ihead)='.|UU MOVE ALL (MEN added 091023)'
                ihead=ihead+1                       
                cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 )			
                call condmom(im,( cosexrel(LA:UA,:) .AND.  jointemp(LA:UA,:)==4  .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==1 ),mom,cnt,var)		
                write(name(im),'("ee|uu move ALL ")')   
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND.  jointemp(LA:UA,:)==4  .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==2 ),mom,cnt,var)		
                write(name(im),'("eu|uu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND.  jointemp(LA:UA,:)==4  .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==3 ),mom,cnt,var)		
                write(name(im),'("ue|uu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                call condmom(im,( cosexrel(LA:UA,:) .AND.  jointemp(LA:UA,:)==4 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==4 ),mom,cnt,var)		
                write(name(im),'("uu|uu move ALL ")')  
                weights(im)=0.0_dp  
                im=im+1 
                

                

                
                
                
                
                
                
                







                LA=MNA ; UA=MXAD
                headloc(ihead)=im
                headstr(ihead)='.|EE MOVE BY TYPE (MEN added 091023)'
                ihead=ihead+1                       
                do i=1,ntypp
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 .AND. dat(LA:UA,:)%typ==i )			
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==1 ),mom,cnt,var)		
                    write(name(im),'("ee|ee move type ",i4)') i  
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==2 ),mom,cnt,var)		
                    write(name(im),'("eu|ee move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1  .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==3 ),mom,cnt,var)		
                    write(name(im),'("ue|ee move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==1 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==4 ),mom,cnt,var)		
                    write(name(im),'("uu|ee move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                end do

                headloc(ihead)=im
                headstr(ihead)='.|EU MOVE BY TYPE (MEN added 091023)'
                ihead=ihead+1                       
                do i=1,ntypp
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 .AND. dat(LA:UA,:)%typ==i )			
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==1 ),mom,cnt,var)		
                    write(name(im),'("ee|eu move type ",i4)') i  
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==2 ),mom,cnt,var)		
                    write(name(im),'("eu|eu move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2  .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==3 ),mom,cnt,var)		
                    write(name(im),'("ue|eu move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==2  .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==4 ),mom,cnt,var)		
                    write(name(im),'("uu|eu move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                end do

                headloc(ihead)=im
                headstr(ihead)='.|UE MOVE BY TYPE (MEN added 091023)'
                ihead=ihead+1                       
                do i=1,ntypp
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 .AND. dat(LA:UA,:)%typ==i )			
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==1 ),mom,cnt,var)		
                    write(name(im),'("ee|ue move type ",i4)') i  
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==2 ),mom,cnt,var)		
                    write(name(im),'("eu|ue move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==3 ),mom,cnt,var)		
                    write(name(im),'("ue|ue move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==3 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==4 ),mom,cnt,var)		
                    write(name(im),'("uu|ue move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                end do


                headloc(ihead)=im
                headstr(ihead)='.|UU MOVE BY TYPE (MEN added 091023)'
                ihead=ihead+1                       
                do i=1,ntypp
                    cosexrel(LA:UA,:)= (dat(LA:UA,:)%co==co  .and. dat(LA:UA,:)%rel==1 .and. dat(LA:UA,:)%sexr==1 .and. norelchg(LA:UA,:)==1 .AND. dat(LA:UA,:)%typ==i )			
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==4 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==1 ),mom,cnt,var)		
                    write(name(im),'("ee|uu move type ",i4)') i  
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==4 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==2 ),mom,cnt,var)		
                    write(name(im),'("eu|uu move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==4 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==3 ),mom,cnt,var)		
                    write(name(im),'("ue|uu move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                    call condmom(im,( cosexrel(LA:UA,:) .AND. jointemp(LA:UA,:)==4 .AND. move(LA:UA,:)==1 ),   d1*one( jointemp(LA+1:UA+1,:)==4),mom,cnt,var)		
                    write(name(im),'("uu|uu move type ",i4)') i 
                    weights(im)=0.0_dp  
                    im=im+1 
                end do

            !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                    do ddd=1,2
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                        headloc(ihead)=im
                        if (g==1.and.j==0) headstr(ihead)='single men emptrans by type 18-49 or 20-34'
                        if (g==1.and.j==1) headstr(ihead)='married men emptrans by type 18-49 or 20-34'
                        if (g==2.and.j==0) headstr(ihead)='single fem emptrans by type 18-49 or 20-34'
                        if (g==2.and.j==1) headstr(ihead)='married fem emptrans by type 18-49 or 20-34'
                        ihead=ihead+1                       
                        if (ddd==1) agerange=1849
                        if (ddd==2) agerange=2034 
                        if (ddd==1) then ;  LA=18 ; UA=49 ; end if 
                        if (ddd==2) then ;  LA=20 ; UA=34 ; end if 

                        call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | u move ALL ",i6)') agerange  
                        momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=6 
                        weights(im)=whour
                        im=im+1 
                        do i=1,ntypp
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | u move type ",i4)') i  
                            momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99               ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=6
                            weights(im)=whour
                            im=im+1 
                        end do 



                        call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | e move ALL ",i6)') agerange   
                        weights(im)=whour
                        momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6
                        im=im+1 
                        do i=1,ntypp
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | e move type ",i4)') i  
                            momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99               ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6
                            weights(im)=whour
                            im=im+1 
                        end do 




                        call condmom(im,( cosexrel(LA:UA,:).AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | u move NO KID ",i6)') agerange   
                        momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=1  
                        weights(im)=whour
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:).AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | u move    KID ",i6)') agerange   
                        momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=2  
                        weights(im)=whour
                        im=im+1 
                        do i=1,ntypp
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | u move NO KID ",i4)') i
                            momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99              ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=1  
                            weights(im)=whour
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | u move    KID ",i4)') i  
                            momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99              ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=2  
                            weights(im)=whour
                            im=im+1 
                        end do 
                    



                        
                        call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | u move NO ED ",i6)') agerange  
                        weights(im)=whour
                        momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=1; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=6 
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | u move    ED ",i6)') agerange 
                        momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=2; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                        weights(im)=whour
                        im=im+1 
                        do i=1,ntypp
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | u move NO ED ",i4)')  i
                            momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99               ; momwhich(im,5)=i ; momwhich(im,6)=1; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                            weights(im)=whour
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | u move    ED ",i4)')  i
                            momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99               ; momwhich(im,5)=i ; momwhich(im,6)=2; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                            weights(im)=whour
                            im=im+1 
                        end do 
                        

                        call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | e move NO KID ",i6)') agerange  
                        momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=1  
                        weights(im)=whour
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | e move    KID ",i6)') agerange 
                        momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=6; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=2  
                        weights(im)=whour
                        im=im+1 
                        do i=1,ntypp
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | e move NO KID ",i4)')  i
                            momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99               ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=1  
                            weights(im)=whour
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | e move    KID ",i4)')  i
                            momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99              ; momwhich(im,5)=i ; momwhich(im,6)=6; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=2  
                            weights(im)=whour
                            im=im+1 
                        end do 
                        

                        call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | e move NO ED ",i6)') agerange   
                        momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=1; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6 
                        weights(im)=whour
                        im=im+1 
                        call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e | e move    ED ",i6)') agerange   
                        momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99                   ; momwhich(im,5)=6 ; momwhich(im,6)=2; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6 
                        weights(im)=whour
                        im=im+1 
                        do i=1,ntypp
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | e move NO ED ",i4)')  i
                            momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99              ; momwhich(im,5)=i ; momwhich(im,6)=1; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                            weights(im)=whour
                            im=im+1 
                            call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e | e move    ED ",i4)')  i
                            momwhich(im,1)=3500 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ; momwhich(im,7)=99 ;  momwhich(im,8)=99              ; momwhich(im,5)=i ; momwhich(im,6)=2; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                            weights(im)=whour
                            im=im+1 
                        end do 



                        !1 is which moment, 2 is age, 3 gender, 4 rel, 5 type, 6 edr, 7 edsp, 8 rellen, 9 empstat, 10 staymove, 11 is kid
                        do iesp=1,2 !educ spouse
                            do ier=1,2 !educ r
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                                write(name(im),'("e | u move jointed ",3i4)')  ier,iesp,agerange
                                weights(im)=0.0_dp
                                momwhich(im,1)=3518 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ;  momwhich(im,8)=99              ; momwhich(im,5)=6 ; momwhich(im,6)=ier ; momwhich(im,7)=iesp ; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                                im=im+1 
                                do i=1,ntypp
                                    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                                    write(name(im),'("e | u move jointed typ ",3i4)')  ier,iesp,i
                                    momwhich(im,1)=3518 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ;  momwhich(im,8)=99              ; momwhich(im,5)=i ; momwhich(im,6)=ier ; momwhich(im,7)=iesp ; momwhich(im,9)=0 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                                    weights(im)=0.0_dp
                                    im=im+1 
                                end do !typ 
                            end do !iesp
                        end do !ier
    
                        do iesp=1,2 !educ spouse
                            do ier=1,2 !educ r
                                call condmom(im,( cosexrel(LA:UA,:)  .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                                write(name(im),'("e | e move jointed ",3i4)')  ier,iesp,agerange
                                weights(im)=0.0_dp
                                momwhich(im,1)=3518 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ;  momwhich(im,8)=99              ; momwhich(im,5)=6 ; momwhich(im,6)=ier ; momwhich(im,7)=iesp ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                                im=im+1 
                                do i=1,ntypp
                                    call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==ier.AND. dat(LA:UA,:)%edsp==iesp),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
                                    write(name(im),'("e | e move jointed typ ",3i4)')  ier,iesp,i
                                    weights(im)=0.0_dp
                                    momwhich(im,1)=3518 ; momwhich(im,2)=agerange ; momwhich(im,3)=g ; momwhich(im,4)=j  ;  momwhich(im,8)=99              ; momwhich(im,5)=i ; momwhich(im,6)=ier ; momwhich(im,7)=iesp ; momwhich(im,9)=1 ; momwhich(im,10)=1  ; momwhich(im,11)=6  
                                    im=im+1 
                                end do 
                            end do 
                        end do
        


                    end do  !ddd
                    end do  !j rel
                end do      ! g gender
                
    
                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 .AND. dat(MNA:MXAD,:)%jobswitch>=0),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA:MXAD,:)%jobswitch==0 ),mom,cnt,var)		
                !write(name(im),'("e | e stay,jobswitch=0",tr3)')  
                !weights(im)=whour
                !im=im+1 
                !do ia=MNAD,MXAD,1
                !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%typ==i .AND. dee(ia,:)==1  .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==0 .AND. dat(ia,:)%jobswitch>=0 ),   d1*one( dat(ia,:)%jobswitch==1 ),mom,cnt,var)		
                !    write(name(im),'("e | e stay,jobswitch=1",i4)') ia 
                !    weights(im)=whour
                !    im=im+1 
                !end do     
                !do ia=MNAD,MXAD,1
                !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%typ==i .AND. dee(ia,:)==1  .AND. move(ia,:)==0 .and. dat(ia,:)%jobswitch==0  .and. dat(ia+1,:)%logwr>=0 .AND. dat(ia,:)%logwr>=0),   d1*wdif(ia,:) ,mom,cnt,var)		
                !    write(name(im),'("wdif | stay,jobwitch==0 ",I4)') ia  
                !    weights(im)=wdifww  
                !    calcvar(im)=0 !1
                !    im=im+1 
                !end do     
                !do ia=MNAD,MXAD,1
                !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%typ==i .AND. dee(ia,:)==1  .AND. move(ia,:)==0 .and. dat(ia,:)%jobswitch==1 .and. dat(ia+1,:)%logwr>=0 .AND. dat(ia,:)%logwr>=0),   d1*wdif(ia,:), mom,cnt,var)		
                !    write(name(im),'("wdif | stay,jobwitch==1 ",I4)') ia  
                !    weights(im)=wdifww  
                !    calcvar(im)=0 !1
                !    im=im+1 
                !end do     
                !do ia=MNAD,MXAD,2
                !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%typ==i .AND. dee(ia,:)==1  .AND. move(ia,:)==0 .and. dat(ia,:)%jobswitch==0 ),   d1*one( dat(ia+1,:)%logwr-dat(ia,:)%logwr<0 ),mom,cnt,var)		
                !    write(name(im),'("wdif<0 | stay,jobwitch==0 ",I4)') ia  
                !    weights(im)=wdifww  
                !    calcvar(im)=0 !1
                !    im=im+1 
                !end do     


                !do ia=agestart(NOCOLLEGE)-1,19,1
                !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%typ==i .AND. dat(ia,:)%hhr==0 .AND. dat(ia,:)%hhr>=0  .AND. dat(ia,:)%edr==1),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                !    write(name(im),'("e|u by age ned ",i4)') ia
                !    weights(im)=whour
                !    im=im+1 
                !end do 
    
                !do ia=agestart(COLLEGE)-1,22,1
                !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%typ==i .AND. dat(ia,:)%hhr==0 .AND. dat(ia,:)%hhr>=0  .AND. dat(ia,:)%edr==2),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                !    write(name(im),'("e|u by age ed ",i4)') ia
                !    weights(im)=whour
                !    im=im+1 
                !end do 
                !do ia=agestart(NOCOLLEGE)-1,27,5
                !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%typ==i .AND. dat(ia,:)%hhr>=0  .AND. dat(ia,:)%edr==1),   d1*one( dat(ia,:)%hhr==1 ),mom,cnt,var)		
                !    write(name(im),'("emp by age ned ",i4)') ia
                !    weights(im)=whour
                !    im=im+1 
                !end do 
                !do ia=agestart(COLLEGE)-1,27,5
                !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%typ==i  .AND. dat(ia,:)%hhr>=0  .AND. dat(ia,:)%edr==2),   d1*one( dat(ia,:)%hhr==1 ),mom,cnt,var)		
                !    write(name(im),'("emp by age ed ",i4)') ia
                !    weights(im)=whour
                !    im=im+1 
                !end do 
    
                !do ia=agestart(NOCOLLEGE)-1,27,5
                !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%typ==i .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr==1 .AND. dat(ia+1,:)%logwr>=0  .AND. dat(ia,:)%edr==1),   d1*dat(ia+1,:)%logwr ,mom,cnt,var)		
                !    write(name(im),'("w|u by age ned ",i4)') ia
                !    weights(im)=whour
                !    im=im+1 
                !end do 
                !do ia=agestart(COLLEGE)-1,27,5
                !    call condmom(im,( cosexrel(ia,:)  .AND. dat(ia,:)%typ==i .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr==1 .AND. dat(ia+1,:)%logwr>=0  .AND. dat(ia,:)%edr==2),   d1*dat(ia+1,:)%logwr ,mom,cnt,var)		
                !    write(name(im),'("w|u by age ed ",i4)') ia
                !    weights(im)=whour
                !    im=im+1 
                !end do 
            !********************************************************************************************************************************************************************************************


        end if !typemoments
        
    
    




        !*******************************************************************************************************************************************************************************************************************************
        !START OF LOOP BY SEX AND REL FOR COMMENTED OUT MOMENTS
        !do g=minsex,maxsex   
        !    do j=maxrelo,0,-1
        !        do typosto=0,ntypp,  DETAILEDOUTPUT
        !        !!!!!!if ( onlysingles ) then !  (.not.onlysingles).or.(onlysingles.and.j==0) ) then 
        !        !!!!!!    cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g )			                    
        !        !!!!!!else 
        !        !!!!!!    cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 .and. dat(MNAD:MXAD,:)%typ==typosto )			
        !        !!!!!!end if             
        !        if (j==1) wmove=wmovemar
        !        if (j==0) wmove=wmovesin
        !        if (typosto==0.and.(.not.onlysingles) ) then 
        !            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
        !        else if (typosto==0.and.onlysingles) then 
        !            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g   )			
        !        else 
        !            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 .and. dat(MNAD:MXAD,:)%typ==typosto )			
        !        end if 

                !do ia=mnad,mxad,1
                !    call condmom(im,( cosexrel(ia,:) .AND. dee(ia,:)==1  .AND. move(ia,:)==1 ),   d1*( dat(ia+1,:)%logwr-dat(ia,:)%logwr ),mom,cnt,var)		
                !    write(name(im),'("wdif | move ia ",tr2,i6)')  ia
                !    weights(im)=wdifww
                !    im=im+1 
                !end do   

                !----> do ia=mnad,35,15
                !---->     call condmom(im,( cosexrel(ia,:) .AND. dee(ia,:)==1  .AND. move(ia,:)==1 ),   d1*( dat(ia+1,:)%logwr-dat(ia,:)%logwr ),mom,cnt,var)		
                !---->     write(name(im),'("wdif | move ia ",tr2,i6)')  ia
                !---->     weights(im)=wdifww
                !---->     calcvar(im)=1
                !---->     im=im+1 
                    !call condmom(im,( cosexrel(ia,:) .AND. dee(ia,:)==1  .AND. move(ia,:)==1 ),   d1*( (dat(ia+1,:)%logwr-dat(ia,:)%logwr)**2 ),mom,cnt,var)		
                    !write(name(im),'("wdif2 | move ia ",tr2,i6)')  ia
                    !weights(im)=wdifww
                    !calcvar(im)=5
                    !im=im+1 
                !----> end do   

                
                !----> ia=agestart(NOCOLLEGE)-1
                !----> do iloc=1,nl
                !---->     call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%l==dat(ia,:)%hme ),   d1*one( dat(ia,:)%l==iloc  ),mom,cnt,var)		
                !---->     write(name(im),'("home dist by loc at 17, ned ",tr3,i4)') iloc !added this ahu 121718
                !---->     weights(im)=whome
                !---->     im=im+1 
                !----> end do 
                !----> ia=agestart(COLLEGE)-1
                !----> do iloc=1,nl
                !---->     call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%l==dat(ia,:)%hme ),   d1*one( dat(ia,:)%l==iloc  ),mom,cnt,var)		
                !---->     write(name(im),'("home dist by loc at 21, ed ",tr3,i4)') iloc !added this ahu 121718
                !---->     weights(im)=whome
                !---->     im=im+1 
                !----> end do 


                !----> LA=MNAD ; UA=MXAD
                !do jloc=1,NL
                !    CALL condmom(im,( cosexrel(LA:UA,:).AND.(move(LA:UA,:)==1).AND.(dat(LA+1:UA+1,:)%l>=0).AND.(norelchg(LA:UA,:)==1) .and. (dat(LA:UA,:)%edr==1) .and. dat(LA+1:UA+1,:)%l==jloc ),d1*one( (dat(LA+1:UA+1,:)%logwr-dat(LA:UA,:)%logwr<0) ),mom,cnt,var)
                !    WRITE(name(im),'("prop-of-moves-to ned , wdif move<0",I4)') jloc
                !    weights(im)=wprop
                !    im=im+1 
                !end do  
                !do jloc=1,NL
                !    CALL condmom(im,( cosexrel(LA:UA,:).AND.(move(LA:UA,:)==1).AND.(dat(LA+1:UA+1,:)%l>=0).AND.(norelchg(LA:UA,:)==1) .and. (dat(LA:UA,:)%edr==2) .and. dat(LA+1:UA+1,:)%l==jloc ),d1*one( (dat(LA+1:UA+1,:)%logwr-dat(LA:UA,:)%logwr<0) ),mom,cnt,var)
                !    WRITE(name(im),'("prop-of-moves-to ed , wdif move<0",I4)') jloc
                !    weights(im)=wprop
                !    im=im+1 
                !end do  

                !do ia=MNAD,22,1
                !    call condmom(im,( cosexrel(ia,:) .AND. move(ia,:)==1 .AND. dat(ia,:)%edr==1),   d1*one( dat(ia+1,:)%l==dat(ia+1,:)%hme  ),mom,cnt,var)		
                !    write(name(im),'("move home ned by ia ",tr3,i4)') ia !added this ahu 121718
                !    weights(im)=whome
                !    im=im+1             
                !    call condmom(im,( cosexrel(ia,:) .AND. move(ia,:)==1 .AND. dat(ia,:)%edr==2),   d1*one( dat(ia+1,:)%l==dat(ia+1,:)%hme  ),mom,cnt,var)		
                !    write(name(im),'("move home ed by ia ",tr3,i4)') ia !added this ahu 121718
                !    weights(im)=whome
                !    im=im+1 
                !end do

                !do ia=23,40,6
                !    call condmom(im,( cosexrel(ia,:) .AND. move(ia,:)==1 .AND. dat(ia,:)%edr==1),   d1*one( dat(ia+1,:)%l==dat(ia+1,:)%hme  ),mom,cnt,var)		
                !    write(name(im),'("move home ned by ia ",tr3,i4)') ia !added this ahu 121718
                !    weights(im)=whome
                !    im=im+1 
                !    call condmom(im,( cosexrel(ia,:) .AND. move(ia,:)==1 .AND. dat(ia,:)%edr==2),   d1*one( dat(ia+1,:)%l==dat(ia+1,:)%hme  ),mom,cnt,var)		
                !    write(name(im),'("move home ed by ia ",tr3,i4)') ia !added this ahu 121718
                !    weights(im)=whome
                !    im=im+1 
                !end do


                !call condmom(im,( cosexrel(MNA:MXAD-1,:) .AND. deue(MNA:MXAD-1,:)==1   .AND. dat(MNA+1:MXAD,:)%l==dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l/=dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l>0 ),   d1*( dat(MNA+2:MXA,:)%logwr-dat(MNA:MXAD-1,:)%logwr ),mom,cnt,var)		
                !write(name(im),'("wdif | eue,m ",tr2)')  
                !weights(im)=0.0_dp  
                !calcvar(im)=0 !1
                !im=im+1 
                !call condmom(im,( cosexrel(MNA:MXAD-1,:) .AND. deue(MNA:MXAD-1,:)==1   .AND. dat(MNA+1:MXAD,:)%l==dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l/=dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l>0 ),   d1*( dat(MNA+2:MXA,:)%logwr-dat(MNA:MXAD-1,:)%logwr )**2,mom,cnt,var)		
                !write(name(im),'("wdif2 | eue,m ",tr2)')  
                !weights(im)=0.0_dp
                !calcvar(im)=5
                !im=im+1 
                !call condmom(im,( cosexrel(MNA:MXAD-1,:) .AND. deue(MNA:MXAD-1,:)==1   .AND. dat(MNA+1:MXAD,:)%l==dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l/=dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l>0 ),   d1*( dat(MNA+2:MXA,:)%logwr-dat(MNA:MXAD-1,:)%logwr -mom(im-2))**2,mom,cnt,var)		
                !write(name(im),'("wdif2p | eue,m ",tr2)')  
                !weights(im)=100.0_dp 
                !im=im+1 
                
                !do ia=mna,mxad,10
                !    call condmom(im,( cosexrel(ia,:) .AND. dee(ia,:)==1  .AND. move(ia,:)==0 ),   d1*( dat(ia+1,:)%logwr-dat(ia,:)%logwr ),mom,cnt,var)		
                !    write(name(im),'("wdif | stay ia ",tr2,i6)')  ia
                !    weights(im)=0.0_dp
                !    im=im+1 
                !end do   
                !do ia=mna,mxad,10
                !    call condmom(im,( cosexrel(ia,:) .AND. dee(ia,:)==1  .AND. move(ia,:)==1 ),   d1*( dat(ia+1,:)%logwr-dat(ia,:)%logwr ),mom,cnt,var)		
                !    write(name(im),'("wdif | move ia ",tr2,i6)')  ia
                !    weights(im)=0.0_dp
                !    im=im+1 
                !end do   
                !do ia=mna,mxad,5
                !    call condmom(im,( cosexrel(ia,:) .AND. dee(ia,:)==1  .AND. move(ia,:)==0 ),   d1*( dat(ia+1,:)%logwr-dat(ia,:)%logwr ),mom,cnt,var)		
                !    write(name(im),'("wdif | stay ia ",tr2,i6)')  ia
                !    weights(im)=wdifww
                !    calcvar(im)=0
                !    im=im+1 
                    !call condmom(im,( cosexrel(ia,:) .AND. dee(ia,:)==1  .AND. move(ia,:)==0 ),   d1*( dat(ia+1,:)%logwr-dat(ia,:)%logwr )**2,mom,cnt,var)		
                    !write(name(im),'("wdif2 | stay ia ",tr2,i6)')  ia
                    !weights(im)=wdifww
                    !calcvar(im)=5
                    !im=im+1 
                !end do   
                !do ia=mna,mxad,10
                !    call condmom(im,( cosexrel(ia,:) .AND. dee(ia,:)==1  .AND. move(ia,:)==0 .AND. ( dat(ia+1,:)%logwr-dat(ia,:)%logwr>0 ) ),   d1*( dat(ia+1,:)%logwr-dat(ia,:)%logwr ),mom,cnt,var)		
                !    write(name(im),'("wdif | stay ia,wdif>0 ",tr2,i6)')  ia
                !    weights(im)=wdifww
                !    calcvar(im)=0
                !    im=im+1 
                    !call condmom(im,( cosexrel(ia,:) .AND. dee(ia,:)==1  .AND. move(ia,:)==0 .AND. ( dat(ia+1,:)%logwr-dat(ia,:)%logwr>0 ) ),   d1*( dat(ia+1,:)%logwr-dat(ia,:)%logwr )**2,mom,cnt,var)		
                    !write(name(im),'("wdif2 | stay ia,wdif>0 ",tr2,i6)')  ia
                    !weights(im)=wdifww
                    !calcvar(im)=5
                    !im=im+1 
                !end do   
                
                
                

                !----> headloc(ihead)=im
                !----> if (g==1.and.j==0) headstr(ihead)='single men emptrans by ed'
                !----> if (g==1.and.j==1) headstr(ihead)='married men emptrans by ed'
                !----> if (g==2.and.j==0) headstr(ihead)='single fem emptrans by ed'
                !----> if (g==2.and.j==1) headstr(ihead)='married fem emptrans by ed'
                !----> ihead=ihead+1
                

                !do ia=MNAD,MXAD,4
                !    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                !    write(name(im),'("e|u by ia ",i4)') ia
                !    weights(im)=whour
                !    im=im+1 
                !end do 

                !do ia=MNAD,MXAD
                !    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .and. move(ia,:)==0),   d1*one( dat(ia+1,:)%hhr==1  ),mom,cnt,var)		
                !    write(name(im),'("e|u stay by ia ",i4)') ia
                !    weights(im)=whour
                !    im=im+1 
                !end do 

                !do ia=MNAD,MXAD
                !    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .and. move(ia,:)==1),   d1*one( dat(ia+1,:)%hhr==1  ),mom,cnt,var)		
                !    write(name(im),'("e|u move by ia ",i4)') ia
                !    weights(im)=whour
                !    im=im+1 
                !end do 

                !do ddd=1,ntypp
                !    do jj=1,3
                !        call condmom(im,( cosexrel(MNAD:MXAD,:) .AND. dat(MNAD:MXAD,:)%hhr==0 .AND. dat(MNA:MXA,:)%hhr>=0  .AND. dur(MNAD:MXAD,:)==jj .AND. dat(MNAD:MXAD,:)%edr==1 .and. dat(MNAD:MXAD,:)%typ==ddd),   d1*one( dat(MNA:MXA,:)%hhr==1 ),mom,cnt,var)		
                !        write(name(im),'("e|u by dur,typ ned ",2i4)') jj,ddd
                !        weights(im)=whour
                !        im=im+1 
                !    end do 
                !end do
                !do ddd=1,ntypp
                !    do jj=1,3
                !        call condmom(im,( cosexrel(MNAD:MXAD,:) .AND. dat(MNAD:MXAD,:)%hhr==0 .AND. dat(MNA:MXA,:)%hhr>=0  .AND. dur(MNAD:MXAD,:)==jj .AND. dat(MNAD:MXAD,:)%edr==2 .and.  dat(MNAD:MXAD,:)%typ==ddd),   d1*one( dat(MNA:MXA,:)%hhr==1 ),mom,cnt,var)		
                !        write(name(im),'("e|u by dur,typ ed ",2i4)') jj,ddd
                !        weights(im)=whour
                !        im=im+1 
                !    end do 
                !end do 


                !do jj=1,5,2
                !    call condmom(im,( cosex(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%logwr>=0  .AND. dur(MNA:MXAD,:)==jj ),   d1*dat(MNA+1:MXA,:)%logwr ,mom,cnt,var)		
                !    write(name(im),'("w|u by dur ",2i4)') g,jj
                !    weights(im)=wwage
                !    im=im+1 
                !end do 
                !cell sample size too small 
            !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1  .and.  dat(MNA:MXAD,:)%edr==1 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            ! write(name(im),'("e | u move ned",tr3)')  
                !weights(im)=wtrans 
                !im=im+1 
                !cell sample size too small 
                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1  .and.  dat(MNA:MXAD,:)%edr==2 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                !write(name(im),'("e | u move ed",tr3)')  
                !weights(im)=wtrans 
                !im=im+1 
                !cell sample size too small 
                !do ia=mna, 38,10
                !    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==1 .and.  dat(ia,:)%edr==1),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                !    write(name(im),'("e | u move by ia (ned)",i4)') ia  
                !    weights(im)=0.0_dp
                !    im=im+1 
                !end do 
                !cell sample size too small 
                !do ia=mna, 38,10
                !    call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==1 .and.  dat(ia,:)%edr==2),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                !    write(name(im),'("e | u move by ia (ed)",i4)') ia
                !    weights(im)=0.0_dp
                !    im=im+1 
                !end do 
                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. move(MNA:MXAD,:)==0),mom,cnt,var)		
                !write(name(im),'("e-stay | e ",tr3)')  
                !weights(im)=wtrans 
                !im=im+1 
                
                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. move(MNA:MXAD,:)==0),mom,cnt,var)		
                !write(name(im),'("u-stay | e ",tr3)')  
                !weights(im)=wtrans  
                !im=im+1 
                
                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. move(MNA:MXAD,:)==1),mom,cnt,var)		
                !write(name(im),'("e-move | e ",tr3)')  
                !weights(im)=wtrans  
                !im=im+1 
                
                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. move(MNA:MXAD,:)==1),mom,cnt,var)		
                !write(name(im),'("u-move | e ",tr3)')  
                !weights(im)=wtrans  
                !im=im+1 

                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. move(MNA:MXAD,:)==0),mom,cnt,var)		
                !write(name(im),'("e-stay | u ",tr3)')  
                !weights(im)=wtrans
                !im=im+1 
                
                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. move(MNA:MXAD,:)==0),mom,cnt,var)		
                !write(name(im),'("u-stay | u ",tr3)')  
                !weights(im)=wtrans  
                !im=im+1 
                
                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. move(MNA:MXAD,:)==1),mom,cnt,var)		
                !write(name(im),'("e-move | u ",tr3)')  
                !weights(im)=wtrans  
                !im=im+1 
                
                !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. move(MNA:MXAD,:)==1),mom,cnt,var)		
                !write(name(im),'("u-move | u ",tr3)')  
                !weights(im)=wtrans  
                !im=im+1 

                
                !do i=1,nl
                !    CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==1 .AND. dat(mna:mxa,:)%logwr>=0 .AND. dat(mna:mxa,:)%l==i) ,d1*dat(mna:mxa,:)%logwr,mom,cnt,var)
                !    WRITE(name(im),'("w|noed l",tr1,i4)') i
                !    weights(im)=wwage    
                !    calcvar(im)=1
                !    im=im+1
                !    CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==1 .AND. dat(mna:mxa,:)%logwr>=0  .AND. dat(mna:mxa,:)%l==i) ,d1*  (dat(mna:mxa,:)%logwr**2),mom,cnt,var)
                !    WRITE(name(im),'("wvar|noed l",tr1,i4)') i
                !    weights(im)=wwage         
                !    calcvar(im)=5
                !    im=im+1
                !end do 
                

                !do i=1,nl
                !    CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==2 .AND. dat(mna:mxa,:)%logwr>=0 .AND. dat(mna:mxa,:)%l==i) ,d1*dat(mna:mxa,:)%logwr,mom,cnt,var)
                !    WRITE(name(im),'("w|ed l",tr1,i4)') i
                !    weights(im)=wwage           
                !    calcvar(im)=1
                !    im=im+1
                !    CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==2 .AND. dat(mna:mxa,:)%logwr>=0  .AND. dat(mna:mxa,:)%l==i) ,d1*  (dat(mna:mxa,:)%logwr**2),mom,cnt,var)
                !    WRITE(name(im),'("wvar|ed l",tr1,i4)') i
                !    weights(im)=wwage     
                !    calcvar(im)=5
                !    im=im+1
                !end do 
                !do i=1,nl
                !    CALL condmom(im,(   cosexrel(18:21,:) .AND.  dat(18:21,:)%hhr==1 .AND. dat(18:21,:)%edr==1 .AND. dat(18:21,:)%logwr>=0  .AND. dat(18:21,:)%l==i) ,d1*  (dat(18:21,:)%logwr),mom,cnt,var)
                !    WRITE(name(im),'("w|1821 l",tr1,i4)') i
                !    weights(im)=wwage        
                !    im=im+1
                !    !CALL condmom(im,(   cosexrel(18:21,:) .AND.  dat(18:21,:)%hhr==1 .AND. dat(18:21,:)%edr==1 .AND. dat(18:21,:)%logwr>=0  .AND. dat(18:21,:)%l==i) ,d1*  (dat(18:21,:)%logwr**2-mom(im-1)**2),mom,cnt,var)
                !    !WRITE(name(im),'("w2|1821 l",tr1,i4)') i
                !    !weights(im)=wwage             
                !    !im=im+1
                !    !CALL condmom(im,(   cosexrel(18:21,:) .AND.  dat(18:21,:)%hhr==1 .AND. dat(18:21,:)%edr==1 .AND. dat(18:21,:)%logwr>=0  .AND. dat(18:21,:)%l==i) ,d1*  (dat(18:21,:)%logwr**3-mom(im-1)**3),mom,cnt,var)
                !    !WRITE(name(im),'("w3|1821 l",tr1,i4)') i
                !    !weights(im)=wwage     
                !    !im=im+1
                !end do 
            

                    !headloc(ihead)=im
                    !if (g==1.and.j==0) headstr(ihead)='single men employment overall '
                    !if (g==1.and.j==1) headstr(ihead)='married men employment overall'
                    !if (g==2.and.j==0) headstr(ihead)='single fem employment overall '
                    !if (g==2.and.j==1) headstr(ihead)='married fem employment overall '
                    !ihead=ihead+1

                    !call condmom(im,( cosexrel(mna:mxa,:) .AND. dat(mna:mxa,:)%hhr>=0 ), d1*one( dat(mna:mxa,:)%hhr==1 ),mom,cnt,var)		
                    !write(name(im),'("e ",tr13)')			
                    !weights(im)=whour 
                    !im=im+1 
                    
                    !call condmom(im,( cosexrel(mna:mxa,:) .AND. dat(mna:mxa,:)%hhr>=0 .AND. dat(mna:mxa,:)%edr==1 ),   d1*one( dat(mna:mxa,:)%hhr==1 ),mom,cnt,var)	
                    !write(name(im),'("e | noed",tr6)')  
                    !weights(im)=whour  
                    !im=im+1 
                    
                    !call condmom(im,( cosexrel(mna:mxa,:) .AND. dat(mna:mxa,:)%hhr>=0 .AND. dat(mna:mxa,:)%edr==2 ),   d1*one( dat(mna:mxa,:)%hhr==1 ),mom,cnt,var)		
                    !write(name(im),'("e |   ed",tr6)')  
                    !weights(im)=whour 
                    !im=im+1 

                    !Note about conditioning on age:
                    !The max endage in data is 47. 
                    !The way sim is done, for those whose endage is 47, the last age where variables get recorded is ia-1=46. 
                    !This is because dat(ia-1,.) is recorded for each ia. So for the last age 47, the variables for 46 gets written
                    !But then there is nothing after age 46, despite the fact that we do have people whose endage is 46 (namely 47).      
                    !If one of the conditioning statements is norelchg, then there is nothing after age 45. 
                    !Since any age 46 no relchg would need a age 47 rel. 
                    !This is why emp and wage by age cond on cosexrel only goes until 45. 
                    !do ia=MNAD,23 
                    !    CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr>=0 .and. dat(ia,:)%edr==1 ) ,d1*one(dat(ia,:)%hhr==1),mom,cnt,var)
                    !    WRITE(name(im),'("emp by age ned",tr4,I2)') ia
                    !    weights(im)=whour 
                    !    if (ia<mna) weights(im)=0.0_dp
                    !    im=im+1
                    !end do            
                    !do ia=24,MXA,10 
                    !    CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr>=0 .and. dat(ia,:)%edr==1 ) ,d1*one(dat(ia,:)%hhr==1),mom,cnt,var)
                    !    WRITE(name(im),'("emp by age ned ",tr4,I2)') ia
                    !    weights(im)=whour 
                    !    if (ia<mna) weights(im)=0.0_dp
                    !    im=im+1
                    !end do            
                    !do ia=MNAD,23 
                    !    CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr>=0 .and. dat(ia,:)%edr==2 ) ,d1*one(dat(ia,:)%hhr==1),mom,cnt,var)
                    !    WRITE(name(im),'("emp by age ed",tr4,I2)') ia
                    !    weights(im)=whour 
                    !    if (ia<mna) weights(im)=0.0_dp
                    !    im=im+1
                    !end do            
                    !do ia=24,MXA,10 
                    !    CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr>=0 .and. dat(ia,:)%edr==2 ) ,d1*one(dat(ia,:)%hhr==1),mom,cnt,var)
                    !    WRITE(name(im),'("emp by age ed ",tr4,I2)') ia
                    !    weights(im)=whour 
                    !    if (ia<mna) weights(im)=0.0_dp
                    !    im=im+1
                    !end do         


                    !do ia=24,MXA,10 
                    !    CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr>=0 ) ,d1*one(dat(ia,:)%hhr==1),mom,cnt,var)
                    !    WRITE(name(im),'("emp by age",tr4,I2)') ia
                    !    weights(im)=whour 
                    !    if (ia<mna) weights(im)=0.0_dp
                    !    im=im+1
                    !end do       
                    !call condmom(im,( cosexrel(MNA:MXA,:) .AND. dat(MNA:MXA,:)%kidr==1 .AND. dat(MNA:MXA,:)%hhr>=0 ),   d1*one( dat(MNA:MXA,:)%hhr==1 ),mom,cnt,var)		
                    !write(name(im),'("e | nokid",tr5)') 
                    !weights(im)=whour 
                    !im=im+1 
                    
                    !call condmom(im,( cosexrel(MNA:MXA,:) .AND. dat(MNA:MXA,:)%kidr==2 .AND. dat(MNA:MXA,:)%hhr>=0 ),   d1*one( dat(MNA:MXA,:)%hhr==1 ),mom,cnt,var)		
                    !write(name(im),'("e |   kid",tr5)') 
                    !weights(im)=whour 
                    !im=im+1 


        !*******************************************************************************************************************************************************************************************************************************
        !START EVERYONE MOMENTS
        !headloc(ihead)=im
        !headstr(ihead)='everyone: i to j with wdif'
        !ihead=ihead+1
        !LA=MNAD ; UA=MXAD
        !do iloc=1,nl
        !     do jloc=1,nl         
        !     call condmom(im,( coho(LA:UA,:) .AND. dat(LA:UA,:)%l==iloc .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%l==jloc  ),mom,cnt,var)		
        !     write(name(im),'("%i to j",tr3,2i4)') iloc,jloc
        !     weights(im)=whome
        !     im=im+1       
        !     call condmom(im,( coho(LA:UA,:) .AND. dat(LA:UA,:)%l==jloc .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%l==iloc  ),mom,cnt,var)		
        !     write(name(im),'("%j to i",tr3,2i4)') jloc,iloc
        !     weights(im)=whome
        !     im=im+1 
        !     call condmom(im,( coho(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%l==iloc .and. dat(LA+1:UA+1,:)%l==jloc),    d1*wdif(LA:UA,:), mom,cnt,var)		
        !     write(name(im),'("wdif move i to j ",2i4)') iloc,jloc
        !     weights(im)=wdifww 
        !     im=im+1 
        !     call condmom(im,( coho(LA:UA,:) .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%l==jloc .and. dat(LA+1:UA+1,:)%l==iloc),     d1*wdif(LA:UA,:), mom,cnt,var)		
        !     write(name(im),'("wdif move j to i ",2i4)') jloc,iloc
        !     weights(im)=wdifww 
        !     im=im+1 
        !     end do      
        ! end do     

        !loop by only sex
        !do g=minsex,maxsex
        !    cosex(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co .and. dat(MNAD:MXA,:)%sexr==g  )


            !CALL condmom(im,(   cosex(18:19,:) .AND.  dat(18:19,:)%hhr==1 .AND. dat(18:19,:)%edr==1 .AND. dat(18:19,:)%logwr>=0) ,d1*dat(18:19,:)%logwr,mom,cnt,var)
            !WRITE(name(im),'("wnned|1819 ")') 
            !weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp   !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess    
            !calcvar(im)=1
            !im=im+1
            !CALL condmom(im,(   cosex(18:19,:) .AND.  dat(18:19,:)%hhr==1 .AND. dat(18:19,:)%edr==1 .AND. dat(18:19,:)%logwr>=0) ,d1*  (dat(18:19,:)%logwr**2),mom,cnt,var)
            !WRITE(name(im),'("wvarned|1819 ")') 
            !weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
            !calcvar(im)=5
            !im=im+1            
            !CALL condmom(im,(   cosex(20:25,:) .AND.  dat(20:25,:)%hhr==1 .AND. dat(20:25,:)%edr==1 .AND. dat(20:25,:)%logwr>=0) ,d1*dat(20:25,:)%logwr,mom,cnt,var)
            !WRITE(name(im),'("wnned|2025 ")') 
            !weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp   !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess    
            !calcvar(im)=1
            !im=im+1
            !CALL condmom(im,(   cosex(20:25,:) .AND.  dat(20:25,:)%hhr==1 .AND. dat(20:25,:)%edr==1 .AND. dat(20:25,:)%logwr>=0) ,d1*  (dat(20:25,:)%logwr**2),mom,cnt,var)
            !WRITE(name(im),'("wvarned|2025 ")') 
            !weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
            !calcvar(im)=5
            !im=im+1            
            !do i=1,nl
                !CALL condmom(im,(   cosex(18,:) .AND.  dat(18,:)%hhr==1 .AND. dat(18,:)%edr==1 .AND. dat(18,:)%logwr>=0 .AND. dat(18,:)%l==i) ,d1*dat(18,:)%logwr,mom,cnt,var)
                !WRITE(name(im),'("wned|18 by loc small!",i4)') i
                !weights(im)=0.0_dp !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                !calcvar(im)=1
                !im=im+1
                !CALL condmom(im,(   cosex(18,:) .AND.  dat(18,:)%hhr==1 .AND. dat(18,:)%edr==1 .AND. dat(18,:)%logwr>=0  .AND. dat(18,:)%l==i) ,d1*  (dat(18,:)%logwr**2),mom,cnt,var)
                !WRITE(name(im),'("wvarned|18 by loc small!",i4)') i
                !weights(im)=0.0_dp !; if (onlysingles.and.j==1) weights(im)=0.0_dp  !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess           
                !calcvar(im)=5
                !im=im+1
                !CALL condmom(im,(   cosex(18:19,:) .AND.  dat(18:19,:)%hhr==1 .AND. dat(18:19,:)%edr==1 .AND. dat(18:19,:)%logwr>=0 .AND. dat(18:19,:)%l==i) ,d1*dat(18:19,:)%logwr,mom,cnt,var)
                !WRITE(name(im),'("wned|18:19 by loc",i4)') i
                !weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                !calcvar(im)=1
                !im=im+1
                !CALL condmom(im,(   cosex(18:19,:) .AND.  dat(18:19,:)%hhr==1 .AND. dat(18:19,:)%edr==1 .AND. dat(18:19,:)%logwr>=0  .AND. dat(18:19,:)%l==i) ,d1*  (dat(18:19,:)%logwr**2),mom,cnt,var)
                !WRITE(name(im),'("wvarned|18:19 by loc",i4)') i
                !weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp  !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess           
                !calcvar(im)=5
                !im=im+1
                !CALL condmom(im,(   cosex(20:25,:) .AND.  dat(20:25,:)%hhr==1 .AND. dat(20:25,:)%edr==1 .AND. dat(20:25,:)%logwr>=0 .AND. dat(20:25,:)%l==i) ,d1*dat(20:25,:)%logwr,mom,cnt,var)
                !WRITE(name(im),'("wned|20:25 by loc",i4)') i
                !weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                !calcvar(im)=1
                !im=im+1
                !CALL condmom(im,(   cosex(20:25,:) .AND.  dat(20:25,:)%hhr==1 .AND. dat(20:25,:)%edr==1 .AND. dat(20:25,:)%logwr>=0  .AND. dat(20:25,:)%l==i) ,d1*  (dat(20:25,:)%logwr**2),mom,cnt,var)
                !WRITE(name(im),'("wvarned|20:25 by loc",i4)') i
                !weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp  !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess           
                !calcvar(im)=5
                !im=im+1
            !end do 
            !ia=18
            !do i=1,nl
            !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%logwr>=0 .AND. dat(ia,:)%l==i) ,d1*dat(ia,:)%logwr,mom,cnt,var)
            !    WRITE(name(im),'("w|noed l 18",tr1,i4)') i
            !    weights(im)=0.0_dp ; if (onlysingles.and.j==1) weights(im)=0.0_dp            
            !    im=im+1
            !end do 
            
            !ia=18
            !do i=1,nl
            !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%logwr>=0 .AND. dat(ia,:)%l==i) ,d1*dat(ia,:)%logwr,mom,cnt,var)
            !    WRITE(name(im),'("w|ed l  18",tr1,i4)') i
            !    weights(im)=0.0_dp ; if (onlysingles.and.j==1) weights(im)=0.0_dp            
            !    im=im+1
            !end do 
        !end do !g for sex
        !******************************************


        
        !headloc(ihead)=im; headstr(ihead)=' ';ihead=ihead+1		
        !do i=1,nl
        !    call condmom(im,( coho(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%l==i  .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)		
        !    write(name(im),'("mvfr | loc",tr4,i4)') i
        !    weights(im)=wmove 
        !    im=im+1 
        !end do

        !do i=1,nl
        !    call condmom(im,( coho(MNA:MXAD,:) .AND. dat(MNA+1:MXA,:)%l==i .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)		
        !    write(name(im),'("mvto | loc",tr4,i4)') i
        !    weights(im)=wmove
        !    im=im+1 
        !end do 

        !do g=minsex,maxsex
        !    cosex(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co .and. dat(MNAD:MXA,:)%sexr==g  )
        !    ia=MNA
        !    do i=1,nl
        !        call condmom(im,( cosex(ia,:) .AND. dat(ia,:)%hhr==1 .AND. dat(ia,:)%l==i .AND. dat(ia,:)%logwr>=0 ),   d1*dat(ia,:)%logwr ,mom,cnt,var)		
        !        write(name(im),'("w|loc a=18 dur,sex",2i4)') i,g	
        !        weights(im)=0.0_dp
        !        im=im+1 
        !    end do      
        !    ia=MNA
        !    do i=1,nl,8
        !        call condmom(im,( cosex(ia:ia+3,:) .AND. dat(ia:ia+3,:)%hhr==1 .AND. dat(ia:ia+3,:)%l==i .AND. dat(ia:ia+3,:)%logwr>=0 ),   d1*dat(ia:ia+3,:)%logwr ,mom,cnt,var)		
        !        write(name(im),'("w|loc a=18:21 dur,sex",2i4)') i,g
        !        weights(im)=wwage
        !        im=im+1 
        !    end do 
        !end do !g sex
        !END OF EVERYONE MOMENTS

        !*******************************************************************************************************************************************************************************************************************************



    
    ENDDO cohort
    
    
    
    
    end subroutine get_mom
    end module mom
	






        !*******************************************************************************************************************************************************************************************************************************
        !START OF THE GENDER,ED,TYPE LOOP
        !do g=0,2 !, DETAILEDOUTPUT
            !do edo=0,2 !, DETAILEDOUTPUT
                !do typosto=0,ntypp !, DETAILEDOUTPUT
                    !IF (  (g==0.and.edo==0.and.typosto==0) .or. (g==0.and.edo==0.and.typosto>0)   .or.(g>0.and.edo==0.and.typosto==0) .or. (g>0.and.edo>0.and.typosto==0) ) THEN  !.or.(g>0.and.edo==0.and.typosto>0)
                    !    checkco=0	            
                    !if (co==1.and.g==0.and.edo==0.and.typosto==0) then                  !ALL SEXES, ALL ED, ALL TYPES
                    !    cosexedtypid(:)            = (inito(:)%co==co  )			
                    !    cosexedtyp(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co  )		
                    !    checkco=1
                    !else if (co==1.and.g==0.and.edo==0.and.typosto>0) then              ! ALL SEXES, ALL ED, BUT BY TYPE
                    !    cosexedtypid(:)            = (inito(:)%co==co .and. inito(:)%typ==typosto )			
                    !    cosexedtyp(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co  .and. dat(MNAD:MXA,:)%typ==typosto )	
                    !    checkco=1    
                    !else if (co==1.and.g>0.and.edo==0.and.typosto==0) then              !BY SEX, ALL ED, ALL TYPE
                    !    cosexedtypid(:)            = (inito(:)%co==co .and. inito(:)%sexr==g   )			
                    !    cosexedtyp(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co  .and. dat(MNAD:MXA,:)%sexr==g   )			
                    !    checkco=1
                    !else  if (co==1.and. g>0.and.edo>0.and.typosto==0) then             !BY SEX, BY ED, ALL TYPES
                    !    cosexedtypid(:)            = (inito(:)%co==co .and. inito(:)%sexr==g  .and. inito(:)%edr==edo  )			
                    !    cosexedtyp(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co  .and. dat(MNAD:MXA,:)%sexr==g  .and. dat(MNAD:MXA,:)%edr==edo  )			
                    !    checkco=1	
                    !end if 
                    !if (checkco==0) then ; print*, "STOP COSEX NOT RIGHT " ; stop ; end if 
                    !wmove=wmoveall
                    !headloc(ihead)=im
                    !if (co==1.and.g==0.and.edo==0.and.typosto==0) headstr(ihead)='everyone rates '
                    !if (co==1.and.g==0.and.edo==0.and.typosto==1) headstr(ihead)='everyone type=1 '
                    !if (co==1.and.g==0.and.edo==0.and.typosto==2) headstr(ihead)='everyone type=2 '
                    !if (co==1.and.g==0.and.edo==0.and.typosto==3) headstr(ihead)='everyone type=3 '
                    !if (co==1.and.g==0.and.edo==0.and.typosto==4) headstr(ihead)='everyone type=4 '
                    !if (co==1.and.g==1.and.edo==0.and.typosto==0) headstr(ihead)='all men '
                    !if (co==1.and.g==2.and.edo==0.and.typosto==0) headstr(ihead)='all fem '
                    !if (co==1.and.g==1.and.edo==1.and.typosto==0) headstr(ihead)='all men ned '
                    !if (co==1.and.g==2.and.edo==1.and.typosto==0) headstr(ihead)='all fem ned '
                    !if (co==1.and.g==1.and.edo==2.and.typosto==0) headstr(ihead)='all men ed '
                    !if (co==1.and.g==2.and.edo==2.and.typosto==0) headstr(ihead)='all fem ed '
                    !if (co==1.and.g==1.and.edo==0.and.typosto==1) headstr(ihead)='all men type 1 '
                    !if (co==1.and.g==2.and.edo==0.and.typosto==1) headstr(ihead)='all fem type 1 '
                    !if (co==1.and.g==1.and.edo==0.and.typosto==2) headstr(ihead)='all men type 2 '
                    !if (co==1.and.g==2.and.edo==0.and.typosto==2) headstr(ihead)='all fem type 2 '
                    !if (co==1.and.g==1.and.edo==0.and.typosto==3) headstr(ihead)='all men type 3 '
                    !if (co==1.and.g==2.and.edo==0.and.typosto==3) headstr(ihead)='all fem type 3 '
                    !if (co==1.and.g==1.and.edo==0.and.typosto==4) headstr(ihead)='all men type 4 '
                    !if (co==1.and.g==2.and.edo==0.and.typosto==4) headstr(ihead)='all fem type 4 '
                    !ihead=ihead+1








            !do ia=mna,mxad,10
            !    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==1 .AND. dat(ia+1,:)%rel>=0 .AND. move(ia,:)==1), d1*one(dat(ia+1,:)%rel==0),mom,cnt,var)
            !    write(name(im),'("getdiv move by ia ",i4)') ia
            !    weights(im)=0.0_dp
            !    im=im+1
            !    do i=1,ntypp            
            !        call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==1 .AND. dat(ia+1,:)%rel>=0  .AND. dat(ia,:)%typ==i .AND. move(ia,:)==1), d1*one(dat(ia+1,:)%rel==0),mom,cnt,var)
            !        write(name(im),'("getdiv move by ia,typ ",2i4)') ia,i
            !        weights(im)=0.0_dp
            !        im=im+1
            !    end do            
            !end do 
!
            !do ia=mna,mxad,10
            !    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==1 .AND. dat(ia+1,:)%rel>=0 .AND. move(ia,:)==0), d1*one(dat(ia+1,:)%rel==0),mom,cnt,var)
            !    write(name(im),'("getdiv stay by ia ",i4)') ia
            !    weights(im)=0.0_dp
            !    im=im+1
            !    do i=1,ntypp            
            !        call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==1 .AND. dat(ia+1,:)%rel>=0  .AND. dat(ia,:)%typ==i .AND. move(ia,:)==0), d1*one(dat(ia+1,:)%rel==0),mom,cnt,var)
            !        write(name(im),'("getdiv stay by ia,typ ",2i4)') ia,i
            !        weights(im)=0.0_dp
            !        im=im+1
            !    end do            
            !end do 




            !headloc(ihead)=im
            !headstr(ihead)='wdif by typ'
            !ihead=ihead+1
            !do g=minsex,maxsex
            !    do j=maxrelo,0,-1
            !        if ( onlysingles ) then !  (.not.onlysingles).or.(onlysingles.and.j==0) ) then 
            !            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g )			                    
            !        else 
            !            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
            !        end if 

            !        do i=1,ntypp
            !            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1  .AND. move(MNA:MXAD,:)==1 .and. dat(MNA+1:MXA,:)%l/=dat(MNA+1:MXA,:)%hme  .AND. dat(MNA:MXAD,:)%typ==i ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr ),mom,cnt,var)		
            !            write(name(im),'("wdif | hmemve0 sex,rel,typ ",3i4)') g,j,i
            !            weights(im)=wdifww
            !            calcvar(im)=0 !dont forget to set these to 0 if there is no wdif2 following this moment 
            !            im=im+1 
            !            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1  .AND. move(MNA:MXAD,:)==1 .and. dat(MNA+1:MXA,:)%l==dat(MNA+1:MXA,:)%hme  .AND. dat(MNA:MXAD,:)%typ==i ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr ),mom,cnt,var)		
            !            write(name(im),'("wdif | hmemve1 sex,rel,typ ",3i4)') g,j,i 
            !            weights(im)=wdifww
            !            calcvar(im)=0 !dont forget to set these to 0 if there is no wdif2 following this moment 
            !            im=im+1 
            !            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1  .AND. move(MNA:MXAD,:)==1  .AND. dat(MNA:MXAD,:)%typ==i ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr ),mom,cnt,var)		
            !            write(name(im),'("wdif | move sex,rel,typ    ",3i4)') g,j,i
            !            weights(im)=wdifww 
            !            calcvar(im)=0 !dont forget to set these to 0 if there is no wdif2 following this moment 
            !            im=im+1             
            !        end do !typ
            !    end do !j rel
            !end do !g sex 

            !headloc(ihead)=im
            !headstr(ihead)='wage by type'
            !ihead=ihead+1
            !do g=minsex,maxsex
            !    do j=maxrelo,0,-1
            !        !if ( onlysingles ) then !  (.not.onlysingles).or.(onlysingles.and.j==0) ) then 
            !        !    cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g )			                    
            !        !else 
            !            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
            !        !end if 
    
            !        CALL condmom(im,(   cosexrel(MNA:MXAD,:) .AND.  dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%edr==1 .AND. dat(MNA:MXAD,:)%logwr>=0  ) ,d1*dat(MNA:MXAD,:)%logwr,mom,cnt,var)
            !        WRITE(name(im),'("w|ned by sex,rel",2I4)') g,j
            !        weights(im)=0.0_dp
            !        im=im+1                        
            !        do i=1,ntypp
            !            CALL condmom(im,(   cosexrel(MNA:MXAD,:) .AND.  dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%edr==1 .AND. dat(MNA:MXAD,:)%logwr>=0  .AND. dat(MNA:MXAD,:)%typ==i ) ,d1*dat(MNA:MXAD,:)%logwr,mom,cnt,var)
            !            WRITE(name(im),'("w|ned by sex,rel,typ",3I4)') g,j,i
            !            weights(im)=0.0_dp
            !            im=im+1
            !        end do

            !        CALL condmom(im,(   cosexrel(MNA:MXAD,:) .AND.  dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%edr==2 .AND. dat(MNA:MXAD,:)%logwr>=0  ) ,d1*dat(MNA:MXAD,:)%logwr,mom,cnt,var)
            !        WRITE(name(im),'("w|ed by sex,rel ",2I4)') g,j
            !        weights(im)=0.0_dp
            !        im=im+1
            !        do i=1,ntypp
            !            CALL condmom(im,(   cosexrel(MNA:MXAD,:) .AND.  dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%edr==2 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%typ==i ) ,d1*dat(MNA:MXAD,:)%logwr,mom,cnt,var)
            !            WRITE(name(im),'("w|ed by sex,rel,typ",3I4)') g,j,i
            !            weights(im)=0.0_dp
            !            im=im+1
            !        end do !typ
                
            !    end do !j rel
            !end do !g sex


            !do g=minsex,maxsex   
            !    do j=maxrelo,0,-1
            !        if ( onlysingles ) then !  (.not.onlysingles).or.(onlysingles.and.j==0) ) then 
            !            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g )			                    
            !        else 
            !            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
            !        end if 
            !        headloc(ihead)=im
            !        if (g==1.and.j==0) headstr(ihead)='single men eumv by ia'
            !        if (g==1.and.j==1) headstr(ihead)='married men eumv by ia'
            !        if (g==2.and.j==0) headstr(ihead)='single fem eumv by ia'
            !        if (g==2.and.j==1) headstr(ihead)='married fem eumv by ia'
            !        ihead=ihead+1
            !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | u move",tr3)')  
            !        weights(im)=0.0_dp
            !        im=im+1 
           ! 
           !         do ia=MNA,MXAD
            !        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==1 ),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | u move by ia",I4)') ia  
            !        weights(im)=0.0_dp
            !        im=im+1 
            !        end do 
            !    end do 
            !end do 

            !headloc(ihead)=im
            !headstr(ihead)='e|u dur by type '
            !ihead=ihead+1
            !do g=minsex,maxsex
            !    do j=maxrelo,0,-1
            !        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
            !        do jj=1,7,3
            !            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  .AND. dur(MNA:MXAD,:)==jj),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !            write(name(im),'("e|u by sex,rel,dur",3i4)') g,j,jj
            !            weights(im)=whour
            !            im=im+1 
            !            do i=1,ntypp
            !                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  .AND. dur(MNA:MXAD,:)==jj  .AND. dat(MNA:MXAD,:)%typ==i),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e|u by sex,rel,dur,typ ",4i4)') g,j,jj,i
            !                weights(im)=0.0_dp 
            !                im=im+1 
            !            end do 
            !        end do 
            !    end do !rel
            !end do !sex
            

            !headloc(ihead)=im
            !headstr(ihead)='w|u dur by type '
            !ihead=ihead+1
            !do g=minsex,maxsex
            !    do j=maxrelo,0,-1
            !        !if ( onlysingles ) then !  (.not.onlysingles).or.(onlysingles.and.j==0) ) then 
            !        !    cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g )			                    
            !        !else 
            !            cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
            !        !end if 

            !        do jj=1,5,4
            !            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%logwr>=0  .AND. dur(MNA:MXAD,:)==jj ),   d1*dat(MNA+1:MXA,:)%logwr ,mom,cnt,var)		
            !            write(name(im),'("w|u by sex,rel,dur ",3I4)') g,j,jj
            !            weights(im)=wwage
            !            im=im+1 
            !            do i=1,ntypp
            !                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%logwr>=0  .AND. dur(MNA:MXAD,:)==jj   .AND. dat(MNA:MXAD,:)%typ==i ),   d1*dat(MNA+1:MXA,:)%logwr ,mom,cnt,var)		
            !                write(name(im),'("w|u by sex,rel,dur,typ ",4i4)') g,j,jj,i
            !                weights(im)=0.0_dp
            !                im=im+1 
            !            end do   !typ i              
            !        end do !dur j 
            !    end do !rel j 
            !end do !sex g
            
            
            !********************************************************************************************************************************************************************************************
            !do g=minsex,maxsex
            !do j=maxrelo,0,-1
            !do i=1,ntypp
            !        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
            !        headloc(ihead)=im
            !        if (i==1) then 
            !            if (g==1.and.j==0) headstr(ihead)='single men type 1'
            !            if (g==1.and.j==1) headstr(ihead)='married men type 1'
            !            if (g==2.and.j==0) headstr(ihead)='single fem type 1'
            !            if (g==2.and.j==1) headstr(ihead)='married fem type 1'
            !        else if (i==2) then 
            !            if (g==1.and.j==0) headstr(ihead)='single men type 2'
            !            if (g==1.and.j==1) headstr(ihead)='married men type 2'
            !            if (g==2.and.j==0) headstr(ihead)='single fem type 2'
            !            if (g==2.and.j==1) headstr(ihead)='married fem type 2'
            !        else if (i==3) then 
            !            if (g==1.and.j==0) headstr(ihead)='single men type 3'
            !            if (g==1.and.j==1) headstr(ihead)='married men type 3'
            !            if (g==2.and.j==0) headstr(ihead)='single fem type 3'
            !            if (g==2.and.j==1) headstr(ihead)='married fem type 3'
            !        else if (i==4) then 
            !            if (g==1.and.j==0) headstr(ihead)='single men type 4'
            !            if (g==1.and.j==1) headstr(ihead)='married men type 4'
            !            if (g==2.and.j==0) headstr(ihead)='single fem type 4'
            !            if (g==2.and.j==1) headstr(ihead)='married fem type 4'
            !        end if 
            !        ihead=ihead+1
            !        !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        !write(name(im),'("e | u stay",tr3))') 
            !        !im=im+1 
            !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | u move",tr3)')  
            !        weights(im)=whour
            !        im=im+1 
            !        !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        !write(name(im),'("e | e stay",tr3)')  
            !        !weights(im)=whour
            !        !im=im+1 
            !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | e move",tr3)')  
            !        weights(im)=whour
            !        im=im+1 
!
            !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. dat(MNA:MXAD,:)%kidr==1),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | u move NO KID",tr3)')  
            !        weights(im)=whour
            !        im=im+1 
            !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. dat(MNA:MXAD,:)%kidr==2),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | u move    KID",tr3)')  
            !        weights(im)=whour
            !        im=im+1 
!
            !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. dat(MNA:MXAD,:)%edr==1),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | u move NO ED",tr3)')  
            !        weights(im)=whour
            !        im=im+1 
            !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. dat(MNA:MXAD,:)%edr==2),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | u move    ED",tr3)')  
            !        weights(im)=whour
            !        im=im+1 
!
!
            !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. dat(MNA:MXAD,:)%kidr==1),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | e move NO KID",tr3)')  
            !        weights(im)=whour
            !        im=im+1 
            !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. dat(MNA:MXAD,:)%kidr==2),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | e move    KID",tr3)')  
            !        weights(im)=whour
            !        im=im+1 
            !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. dat(MNA:MXAD,:)%edr==1),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | e move NO ED",tr3)')  
            !        weights(im)=whour
            !        im=im+1 
            !        call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%typ==i .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 .AND. dat(MNA:MXAD,:)%edr==2),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            !        write(name(im),'("e | e move    ED",tr3)')  
            !        weights(im)=whour
            !        im=im+1 
            !    end do 
            !end do
            !end do 
            


            !headloc(ihead)=im
            !g=2 ; headstr(ihead)='fem emptrans by kid,ed,type and also age'
            !ihead=ihead+1   
            !DO KK=1,2
            !    DO EDO=1,2
            !        do i=1,ntypp
            !            do ia=MNA,35
            !                cosexrel(ia,:)= (dat(ia,:)%co==co .and. dat(ia,:)%sexr==g  .and. dat(ia,:)%rel==1 .and. norelchg(ia,:)==1 )			
            !                call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%typ==i .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==0  .AND. dat(ia,:)%kidr==kk.AND. dat(ia,:)%edr==edo),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e|ustay MAR KID,ED,TY ",4i4)') KK,EDO,I,ia
            !                weights(im)=whour
            !                im=im+1 
            !                cosexrel(ia,:)= (dat(ia,:)%co==co .and. dat(ia,:)%sexr==g  .and. dat(ia,:)%rel==0 .and. norelchg(ia,:)==1 )			
            !                call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%typ==i .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==0  .AND. dat(ia,:)%kidr==kk.AND. dat(ia,:)%edr==edo),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e|ustay SIN KID,ED,TY ",4i4)') KK,EDO,I,ia
            !                weights(im)=whour
            !                im=im+1 
            !            end do
            !        end do 
            !    end do 
            !end do 
            !               

            !do g=minsex,maxsex
            !    do j=maxrelo,0,-1
            !        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
            !            headloc(ihead)=im
            !            if (g==1.and.j==0) headstr(ihead)='single men emptrans by type all ages'
            !            if (g==1.and.j==1) headstr(ihead)='married men emptrans by type all ages'
            !            if (g==2.and.j==0) headstr(ihead)='single fem emptrans by type all ages'
            !            if (g==2.and.j==1) headstr(ihead)='married fem emptrans by type all ages'
            !            ihead=ihead+1               
            !            LA=MNA ; UA=MXAD
            !            do i=1,ntypp
            !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==0 ),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e | u stay type ALL AGES ",i4)') i  
            !                weights(im)=whour
            !                im=im+1 
            !            end do 
            !            do i=1,ntypp
            !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e | u move type ALL AGES ",i4)') i  
            !                weights(im)=whour
            !                im=im+1 
            !            end do 
            !            do i=1,ntypp
            !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 ),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e | e move type ",i4)') i  
            !                weights(im)=whour
            !                im=im+1 
            !            end do 
            !            do i=1,ntypp
            !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e | u move NO KID ",i4)') i  
            !                weights(im)=whour
            !                im=im+1 
            !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e | u move    KID ",i4)') i  
            !                weights(im)=whour
            !                im=im+1 
            !            end do 
            !            
            !            do i=1,ntypp
            !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e | u move NO ED ",i4)')  i
            !                weights(im)=whour
            !                im=im+1 
            !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==0 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e | u move    ED ",i4)')  i
            !                weights(im)=whour
            !                im=im+1 
            !            end do 
            !            
            !            do i=1,ntypp
            !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e | e move NO KID ",i4)')  i
            !                weights(im)=whour
            !                im=im+1 
            !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%kidr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e | e move    KID ",i4)')  i
            !                weights(im)=whour
            !                im=im+1 
            !            end do 
            !            
            !            
            !            do i=1,ntypp
            !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==1),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e | e move NO ED ",i4)')  i
            !                weights(im)=whour
            !                im=im+1 
            !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dat(LA:UA,:)%hhr==1 .AND. dat(LA+1:UA+1,:)%hhr>=0 .AND. move(LA:UA,:)==1 .AND. dat(LA:UA,:)%edr==2),   d1*one( dat(LA+1:UA+1,:)%hhr==1 ),mom,cnt,var)		
            !                write(name(im),'("e | e move    ED ",i4)')  i
            !                weights(im)=whour
            !                im=im+1 
            !            end do 
            !
            !        end do 
            !    end do 
                



                !do g=minsex,maxsex
                !    do j=maxrelo,0,-1
                !        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                !            headloc(ihead)=im
                !            if (g==1.and.j==0) headstr(ihead)='single men wdif by type '
                !            if (g==1.and.j==1) headstr(ihead)='married men wdif by type '
                !            if (g==2.and.j==0) headstr(ihead)='single fem wdif by type '
                !            if (g==2.and.j==1) headstr(ihead)='married fem wdif by type '
                !            ihead=ihead+1     
                !            do i=1,ntypp
                !                LA=MNA ; UA=MXAD
                !                call condmom(im,( cosexrel(LA:UA,:).AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==0 .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*wdif(LA:UA,:) ,mom,cnt,var)		
                !                write(name(im),'("wdif | stay ",tr2)')  
                !                weights(im)=wdifww  
                !                calcvar(im)=0 !1
                !                im=im+1 
                !                call condmom(im,( cosexrel(LA:UA,:) .AND. dat(LA:UA,:)%typ==i .AND. dee(LA:UA,:)==1  .AND. move(LA:UA,:)==1  .and. dat(LA+1:UA+1,:)%logwr>=0 .AND. dat(LA:UA,:)%logwr>=0),   d1*wdif(LA:UA,:) ,mom,cnt,var)		
                !                write(name(im),'("wdif | move ",tr2)')  
                !                weights(im)=wdifww 
                !                calcvar(im)=0 !dont forget to set these to 0 if there is no wdif2 following this moment 
                !                im=im+1  
                !            end do 
                !        end do 
                !end do 

