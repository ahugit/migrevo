
module objf
	use params
	use sol, only: solve
	use mom
	implicit none 
	!include 'mpif.h'
	real(dp) :: q_save(numit),qcont_save(nmom,numit)
	real(dp), dimension(nmom) :: msm_wgt,momwgt
	real(dp), dimension(npars,numit) :: par_save,realpar_save
	real(dp), dimension(nmom,numit) :: momdat_save,momsim_save,vardat_save
	integer(i4b), dimension(nmom,numit) :: cntdat_save,cntsim_save
	character(len=namelen), dimension(nmom) :: name
	character(len=500), dimension(nmom) :: header
	integer(i4b), dimension(nmom) :: headerloc
	real(dp), dimension(nepsmove,numit) :: moveshocksave,totcostsave
	real(dp), dimension(2,numit) :: parcostsave
contains
	subroutine objfunc(parvec,objval)
	! puts everything together: takes the parameter vector, computes msm objective function at that value
	real(8), dimension(npars), intent(in) :: parvec	
	real(8), intent(out) :: objval					
	real(8), dimension(npars) :: realparvec
	type(statevar), dimension(:,:), allocatable :: dat
	real(dp), dimension(nmom), save :: momdat,vardat    !deljan03
	integer(i4b), dimension(nmom), save :: cntdat       !deljan03
	real(dp), dimension(nmom) :: momsim,varsim,mymom,myvar,qcont
	integer(i4b), dimension(nmom) :: cntsim,mycnt
	integer(i4b), parameter :: datcountbar=10 !ahu jan19 010219 changing from 20 to 10 in order for emp at age 18 to count (thereis only 12 in that cell) because it identifies offer rates
	real(dp) :: time(10),mutemp(6)
    real(4) :: timing(11)
	integer(i4b) :: i,t,mycommrank,mpierr   
    !print*, "iter,mysay,iwritegen",iter,mysay,iwritegen !ag090522 agsept2022
    !initiate
	objval=-99.0_dp
    realparvec=-99.0_dp
    momsim=-99.0_dp
    varsim=-99.0_dp
    mymom=-99.0_Dp
    myvar=-99.0_dp
    msm_wgt=-99.0_dp 
    momwgt=-99.0_dp 
    qcont=-99.0_dp    
    cntsim=9999
    mycnt=9999
    
	time(1)=secnds(0.0)
	timing(1)=secnds(0.0)
    if (groups) then 
        nindex=nin
    else 
        nindex=ninp
    end if 
    allocate( decm_postdiv(nepsmove,nxs,nqs,nqs,mna:mxa,nindex),decf_postdiv(nepsmove,nxs,nqs,nqs,mna:mxa,nindex))
    allocate(decm_premarmkt(nepsmove,nxs,nqs,nqs,mna:mxa,nindex),decf_premarmkt(nepsmove,nxs,nqs,nqs,mna:mxa,nindex))
	allocate(vm_postdiv(nepsmove,nxs,nqs,nqs,mna:mxa,nindex),vf_postdiv(nepsmove,nxs,nqs,nqs,mna:mxa,nindex))
    allocate(dec_mar(nz,nx,nq,mna:mxa,nindex))
	allocate(vm0_c(nx,nq,mna:mxa,nindex),vf0_c(nx,nq,mna:mxa,nindex))
    allocate(vm0ctemp(nq,nx),vf0ctemp(nq,nx))

	call getpars(parvec,realparvec) 
    parsforcheck=parvec !jusr for checking where it says emax is negative. can get rid of later.
    !if (iter>1) skriv=.false.
    !if (mysay==0.and.iter==3) then 
    !    skriv=.true.
    !    print*, 'here mysay,iter ',mysay, iter
    !else 
    !    skriv=.false.
    !end if
    !if (iter==1.and.iwritegen==1) then  !ahu 121118
    !    open(unit=98799, file='parameters.txt',status='replace')
    !end if 
	!if (iwritegen==1) then ; write(98799,'("iter ",i6)') iter ; end if  !ahu 121118
	!if (iwritegen==1) then !ahu 121118
    !    do i=1,npars
    !        write(98799,*) parvec(i)
    !    end do 
    !    write(98799,*)
    !end if
    !if (iter==1.and.iwritegen==1) then  !ahu 121118
    !    open(unit=987991, file='paros.txt',status='replace')
    !    do i=1,npars
    !        write(987991,'(1a15,f14.4)') parname(i),realparvec(i)
	!    end do 
    !    close(987991)
    !end if
    
	allocate(decisions(mnad:mxa,numpersim),draws(mnad:mxa,numpersim))
	allocate(myefo(mnad:mxa,numpersim))
	allocate(mdecisions(mnad:mxa,numpersim),fdecisions(mnad:mxa,numpersim))
	
	if (iter==1) then 
		best=huge(1.0_dp)
		momdat_save=0.0_dp
		vardat_save=0.0_dp
		cntdat_save=0
		momsim_save=0.0_dp
		cntsim_save=0
		qcont_save=0.0_dp
		q_save=0.0_dp
		par_save=0.0_dp
		realpar_save=0.0_dp
		parcostsave=0.0_dp
		moveshocksave=0.0_dp
		totcostsave=0.0_dp
	end if 	
	if (iter<=numit) then
        !initiate
        momdat=-99.0_dp ; vardat=-99.0_dp ; cntdat=9999
		call getones                !ones
		call getdistpop	            !distance and popsize	
		call getq2q					! trans btw couple and single state space
		call getx2x					! trans btw couple and single state space
		call getch_single			! choice set cond on state variable and shock for singles
		call getch_couple			! choice set cond on state variable and shock for singles
		allocate(dat(mnad:mxa,numperdat)) !numperdat is the number of persons in actual data) 
		! read taxes
		call read_taxes
		 !init is allocated above !numperdat: num of persons in actual data, numperobsdat: num of person-periods in actual data
		call read_actualdata(dat,numperdat,numperobsdat)
		!  calculate moments from data, store in datamoments
		!  momentname,momentheaders,headerlocs,weights will be filled in when calculate moments from simulations
		!  so that they're not all floating around as globals. loaded into temporary variables here
		call get_mom(dat,numperdat,momdat,cntdat,vardat,name,header,headerloc,momwgt)
        !print*, name(70),momdat(70),cntdat(70)
        do i=1,nmom
            if (calcvar(i)==1) then 
                if ( cntdat(i) > 0) then
                    mutemp(1:2)=momdat(i:(i+1))    !E(X) and E(Xsq)  						!REMEMBER THAT THESE SHOULD BE DIVIDED BY CNTSIM IF YOU WERE DOING THIS CALCULATION FOR SIM MOMENTS
                    momdat(i)= mutemp(1)    	
                    momdat(i+1)= mutemp(2)  - mutemp(1)**2
                end if 
            else if (calcorr(i)==1) then 
                if ( cntdat(i) > 0) then
                    mutemp(1:6)=momdat(i:(i+5))   !( E(X) E(Xsq) E(Y) E(Ysq) E(XY) E(XY)) 	!REMEMBER THAT THESE SHOULD BE DIVIDED BY CNTSIM IF YOU WERE DOING THIS CALCULATION FOR SIM MOMENTS
					call getcorrelation(mutemp)
                    momdat(i:(i+5))=mutemp(1:6)   !( E(X) E(Xsq) E(Y) E(Ysq) E(XY) E(XY)) 	!REMEMBER THAT THESE SHOULD BE DIVIDED BY CNTSIM IF YOU WERE DOING THIS CALCULATION FOR SIM MOMENTS
                end if 
            end if
        end do 
		deallocate(dat)
	end if 
	call getgauss
	call getppsq 
	call getppcq
	call getppsx 
	call getppcx
	call getppmeet    
    
	if (skriv) call yaz0	
    call solve		
    timing(2)=secnds(timing(1))
    timing(3)=secnds(0.0)
    allocate(dat(mnad:mxa,numpersim)) !numpersim=numperdat*nsimeach here and its num of persons in simulation !ag 110416: this is allocated mna-1 rather than mna because simulation is changed to have sim(ia-1,r) at the beginning of the sim loop rather than sim(ia,r) in order to not have 0 emp at age 18
	call simulate(dat,numperdat,numpersim)
    timing(4)=secnds(timing(3))
    timing(5)=secnds(0.0)    
    if (groups) then 
		call get_mom(dat,numpersim,mymom,mycnt,myvar,name,header,headerloc,momwgt)
        !deljan03
        !if ( mysay==0.and.iter==3) then 
        !    do i=1,nmom
        !        write(*,'("i,mymom,momdat,mycnt,q ",I4,2F10.2,I8,F10.2)') i,mymom(i),momdat(i),mycnt(i),q
        !    end do
        !end if        
        !deljan03
	else
		call get_mom(dat,numpersim,momsim,cntsim,varsim,name,header,headerloc,momwgt)
	end if 
	deallocate(dat)
    timing(6)=secnds(timing(5))
    timing(7)=secnds(0.0)    
	if (groups) then 
		call mpi_comm_rank(comm,mycommrank,mpierr)
		call mpi_allreduce(mycnt,cntsim,nmom,mpi_integer,mpi_sum,comm,mpierr)
		call mpi_allreduce(mymom*mycnt,momsim,nmom,mpi_real8,mpi_sum,comm,mpierr)
        !deljan03
        !if ( mysay==0.and.(iter==3.or.iter==2)) then         
        !    do i=1,nmom
        !        write(*,'("mysay,iter,i,momsim,momdat,cntsim,q ",2I4,I8,2F10.2,I8,F10.2)') mysay,iter,i,momsim(i),momdat(i),cntsim(i),q
        !    end do
        !end if 
        !deljan03        
        !print*, 'mygroup,mysay,mycnt,cntsim',mygroup,mysay,mycnt(112),cntsim(112)
	end if 

    !allocate(Iamtrying(mna:mxa))
    !Iamtrying=2
    !Iamtrying(mxa)=1
    !print*, 'Here it is',Iamtrying(1),Iamtrying(mxa)
    !deallocate(Iamtrying)
    
	do i=1,nmom
		!if (vardatamom(im)>0d0.or.countdatamom(im)>=datacountbar) !this "or" does not make any sense! !cohabitation correction
		if (vardat(i)*cntdat(i)>0.0_dp .and. cntdat(i)>15 ) then
			msm_wgt(i)= vardat(i)**(-1) !( real(cntdat(i))/real(numperdat) ) * vardat(i)**(-1) !1.0_dp ahu 041219    !AHU JAN19 012919
		else 
			msm_wgt(i)=0.0_dp
			!else if (cntdat(i)==0) then
		!	msm_wgt(i)=0.0_dp
		!else
		!	msm_wgt(i)=1.0_dp
		end if  
		!if (vardat(i)>0.0_dp .and. cntdat(i)>=datcountbar ) then
		!	msm_wgt(i)=vardat(i)**(-1) !1.0_dp ahu 041219    !AHU JAN19 012919
		!else 
		!	msm_wgt(i)=0.0_dp
		!end if 
		!the version in train_stderr:
		!if (vardat_save(im,1)*(cntdat_save(im,1))>0.) then !these are declared in the objf module and assigned values in objf
		!    msm_wgt(im)=( cntdat_save(im,1) ) * vardat_save(im,1)**(-1)
		!elseif (cntdat_save(im,1)==0) then
		!    msm_wgt(im)=0.
		!else
		!    msm_wgt(im)=1.
		!end if
		!AG090122 AGSEPT2022
		!IF WANT TO COMPARE OBJVAL BETWEEN RUNS GROUPS=TRUE AND GROUPS=FALSE THEN 
		!NEED TO RUN GROUPS=TRUE WITH THE BELOW AND THEN GROUPS=FALSE WITH THE BELOW 
		!INSTEAD OF THE LATTER IF STATEMENT (RIGHT AFTER IT)
		!BECAUSE OTHERWISE WHEN GROUPS=FALSE, IT DOES NOT CALCULATE THE MOMENTS 
		!WITH CALCVAR=1 AND CALCORR=1 THE WAY THEY ARE CALCULATED WHEN GROUPS=TRUE 
		!(I.E. THE WAY THEY ARE CALCULATED WITHIN THE IF STATEMENT WITH MUTEMP'S)
		!BUT IF COMPARING IS NOT THE PURPOSE THEN WHEN GROUPS=TRUE, 
		!THE CALCVAR AND CALCORR MOMENTS NEED TO BE CALCULATED WIHT THOSE MUTEMP'S 
		!I AM NTO JUST GETTING RID OF THESE MOMENTS BECAUSE I NEED THEM (OR THINK SO)
		!if (groups) then 
        !        if ( cntsim(i) > 0) then
		!		    momsim(i)=momsim(i)/cntsim(i)		
        !        end if 
		!end if 
		if (groups) then 
            if (calcvar(i)==0 .and. calcorr(i)==0 ) then
                if ( cntsim(i) > 0) then
				    momsim(i)=momsim(i)/cntsim(i)		!else ; simom(i)=0.0_dp  
                end if 
            else if (calcvar(i)==1) then 
                if ( cntsim(i) > 0) then
                    mutemp(1:2)=momsim(i:(i+1))/cntsim(i:(i+1))  !E(X) and E(Xsq)  		
                    momsim(i)= mutemp(1)   
                    momsim(i+1)= mutemp(2)  - mutemp(1)**2
					!ag090122 agsept2022 : 
					!IF THIS IS NOT COMMENTED OUT:
					!GROUPS TRUE AND FALSE They give the same moments so they should give same objval. 
					!But they don't give same objval mainly because of this calculation above. 
					!When groups false, the momsims dont need to be dividded by cntsims and all the momsims 
					!are already calculated the way they should be by getmom (when groups false).
					!So when groups false, we do not get in these calculations within the if (groups) statement above. 
					!In this if statement, some momsims are calculated with these mutemp's (when calcvar=1 and calcorr=1) (when groups true) 
					!so there is no way to this mutemp calculation when groups true. (for calcvar=1 and calcorr=1 cases)
					!so then instead of trying to figure out how to do the calcvar and calcorr calcuations for groups true
					!i just comment them out and compare objval to groups false. 
                end if 
            else if (calcorr(i)==1) then 
                if ( cntsim(i) > 0) then
                    mutemp(1:6)=momsim(i:(i+5)) /cntsim(i:(i+5))   !( E(X) E(Xsq) E(Y) E(Ysq) E(XY) E(XY)) 	!REMEMBER THAT THESE SHOULD BE DIVIDED BY CNTSIM IF YOU WERE DOING THIS CALCULATION FOR SIM MOMENTS
					call getcorrelation(mutemp)
                    momsim(i:(i+5))=mutemp(1:6)   !( E(X) E(Xsq) E(Y) E(Ysq) E(XY) E(XY)) 	!REMEMBER THAT THESE SHOULD BE DIVIDED BY CNTSIM IF YOU WERE DOING THIS CALCULATION FOR SIM MOMENTS
                end if 
            end if
		end if 
	end do
	qcont=momwgt*msm_wgt*(momdat-momsim)**2
	objval=sum(qcont) 
	timing(8)=secnds(timing(7))
	timing(9)=secnds(0.0)
    !if (iter==1) print*, 'my name is ',mysay,' and iwritegen is ',iwritegen
	!if (iwritegen==1) then ; write(*,'("iter,obj,time: ",i6,f20.2,4f14.2)') iter,objval,timing(2),timing(4),timing(6),timing(8) ; end if  
	!ahu 0317 write(*,'("iter,obj: ",3i6,f20.2,3f14.2)') mygroup,mysay,iter,q,timing(2),timing(4),timing(6)  
    
	! save the moments and objective function values from the first iteration, for comparison to the later ones: 
	t=min(numit,iter)       ! ahu 030517: will always have numit at most 2 (i.e.1 or 2). if you want more comparisons (i.e. numit>2) you have to change this. 
                            ! if numit is 1 and this is coded as t=min(2,iter), then we get out of bounds run time error for momdat_save(.,t) etc.
                            ! so be careful. but setting t=min(numit,iter) should solve that problem and also never having t=iter case as in the if statement below.
    !if (optimize .or. chkstep) then 
	!	t=min(numit,iter)
	!else 
	!	t=iter
	!end if 
	!if (t>numit) then ; print*, "error: iter>numit!!!" ; stop ; end if 
	momdat_save(:,t)=momdat
	vardat_save(:,t)=vardat
	cntdat_save(:,t)=cntdat
	momsim_save(:,t)=momsim
	cntsim_save(:,t)=cntsim
	qcont_save(:,t)=qcont
	q_save(t)=objval
	par_save(:,t)=parvec 
	realpar_save(:,t)=realparvec
	if (iwritegen==1.and.chkobj) then
		i=maxloc(abs(par_save(:,t)-par_save(:,1)),1)
	    write(63,'(i4,". ",1a20,2f20.2,4f14.4)') i,parname(i),q_save(1),q_save(t),par_save(i,1),par_save(i,t),realpar_save(i,1),realpar_save(i,t)
    end if 
	if ((.not.optimize).or.(optimize.and.objval<best)) then	
		best=objval
		if (iwritegen==1 ) then 
			call writemoments(objval) 
			open(unit=749,file='bpobj.txt',status='replace')
			do i=1,npars 
				write(749,*) parvec(i) 
			end do
			close(749)
		end if 
	end if 
    deallocate(decisions,draws,myefo,mdecisions,fdecisions)
    deallocate(decm_postdiv,decf_postdiv)
    deallocate(decm_premarmkt,decf_premarmkt)
    deallocate(vm_postdiv,vf_postdiv)
    deallocate(dec_mar)
    deallocate(vm0_c,vf0_c)
	deallocate(vm0ctemp,vf0ctemp)
	timing(10)=secnds(timing(9))
	timing(11)=secnds(timing(1))
	!if (iwritegen==1) then ; write(*,'(7f14.2)') timing(2),timing(4),timing(6),timing(8),timing(10),timing(2)+timing(4)+timing(6)+timing(8)+timing(10),timing(11) ; end if  
	if (iwritegen==1) then ; write(*,'("iter,obj,time: ",i6,7f14.2)') iter,objval,timing(2),timing(4),timing(6),timing(8),timing(10),timing(11) ; end if  
	iter=iter+1	
	end subroutine objfunc

	! subroutine to be called by parallel simplex optimization whenever it's time to check if have new best point
	! note that called by master. master will have current vale of 'best' from the last time this was called.
	subroutine writebest(parvector,hminvalue,nevalno,hmeanvalue,hstdev)
		real(dp), dimension(npars), intent(in) ::parvector ! transformed vector of parameters
		integer, intent(in) :: nevalno
		real(dp), intent(in) :: hminvalue(1),hmeanvalue,hstdev 
		integer :: i
		!naming these nelder because writebest is called from nelder_mead and these report the hmin of simplex so far
		!potentially might differ from the bestpar written by a call from within objfunc (depends on temperature of annealing)
		open(unit=61, file='bvnel.txt',status='replace')
			write(61,*) 'neval,hmin,hmean,hstd so far from pnmead are:'
			write(61,*) nevalno,hminvalue(1),hmeanvalue,hstdev
			write(61,*) 
			write(61,*) "parameters that correspond to hminvalue so far are:"
			do i=1,npars ; write(61,*) parvector(i) ; end do
		close(61)

		open(unit=66,file='bpnel.txt',status='replace')
		do i=1,npars ; write(66,*) parvector(i) ; end do
		close(66)
	end subroutine

	subroutine writemoments(objval)
	real(8), intent(in) :: objval
	integer(i4b) :: i,t,ihead,j,k,trueindex,iepjoint,iephub,iepwfe,indeces(2),co,typ,homeloc,loca,wdraw,typo
	real(dp) :: temprob(nl)
	real(8) :: wel,fel
	character(len=15) :: parcostname(2),moveshockname,totcostname
	if (momdisplay) then
		open(unit=6091115, file='momdisplay.txt',status='replace')
		do t=1,numit
			do i=1,nmom
			if (momwhich(i,1)>=10.and.momwhich(i,1)<=9900) then
				write(6091115,'(i4,i8,2F20.4,11i6)')	t,i,momsim_save(i,t),momdat_save(i,t),momwhich(i,1:11)
			end if 
			end do
		end do 
		close(6091115)
	end if 
	open(unit=60, file=momentfile,status='replace')
    !open(unit=61 change this 61 to another number since bestval is also 61 maybe among other things, file=momentonlyfile,status='replace')
	!parcostsave(1,iter)=sigo_m
	!parcostsave(2,iter)=cst(4)
	!parcostname(1)='sigo_m'
	!parcostname(2)='cst(4)'
	!moveshocksave(1:nepsmove,iter)=moveshock_m(1:nepsmove)
	!moveshockname='moveshock_m'
	!totcostsave(1:nepsmove,iter)=moveshock_m(1:nepsmove)+cst(4)
	!totcostname='total cost'
	do i=1,npars
        write(60,'(1a15,12f14.5)') parname(i),realpar_save(i,1:numit) 
	end do 
	!do i=1,npars
    !    write(60,'(1a15,9f10.2)') parname(i),realpar_save(i,10:18) 
	!end do 
	!do i=1,npars
    !    write(60,'(1a15,9f10.2)') parname(i),realpar_save(i,19:numit) 
	!end do 
!	do i=1,npars
!        write(60,'(1a15,18f10.2)') parname(i),realpar_save(i,37:54) 
!	end do 
!	do i=1,npars
!        write(60,'(1a15,18f10.2)') parname(i),realpar_save(i,55:numit) 
!	end do


	write(60,*)
	write(60,*) 
    write(60,'("TYPE 1  ",2F9.2)')    alf1t(1),alf2t(1)
	write(60,'("....HS  ",2F9.2)')	ptypeHS(1,MALES) ,   ptypeHS(1,FEMALES)
	write(60,'("....COL ",2F9.2)')	ptypeCOL(1,MALES) , ptypeCOL(1,FEMALES)
	write(60,*)
    write(60,'("TYPE 2  ",2F9.2)')     alf1t(2),alf2t(2)
	write(60,'("....HS  ",2F9.2)')	 ptypeHS(2,MALES) ,   ptypeHS(2,FEMALES)
	write(60,'("....COL ",2F9.2)')	ptypeCOL(2,MALES) , ptypeCOL(2,FEMALES)
	write(60,*)
    write(60,'("TYPE 3  ",2F9.2)')     alf1t(3),alf2t(3)
	write(60,'("....HS  ",2F9.2)')	 ptypeHS(3,MALES) ,   ptypeHS(3,FEMALES)
	write(60,'("....COL ",2F9.2)')	ptypeCOL(3,MALES) , ptypeCOL(3,FEMALES)
	write(60,*)
	write(60,'("TYPE 4  ",2F9.2)')     alf1t(4),alf2t(4)
	write(60,'("....HS  ",2F9.2)')	 ptypeHS(4,MALES) ,   ptypeHS(4,FEMALES)
	write(60,'("....COL ",2F9.2)')	ptypeCOL(4,MALES) , ptypeCOL(4,FEMALES)
	write(60,*)
	write(60,*)

	write(60,'("fnprloc(1,1):")') 
    write(60,'("origin=1,homeloc=1    ",2x,9F9.2)') fnprloc(1,1)  !prob(nthing) is the residual so only psio(1:2) governs this
	write(60,'("fnprloc(1,5):")') 
    write(60,'("origin=1,homeloc=2    ",2x,9F9.2)') fnprloc(1,5)  !prob(nthing) is the residual so only psio(1:2) governs this


	!write(60,'(1a15,2f9.1)') parcostname(1),parcostsave(1,1:numit) 
	!write(60,'(1a15,2f9.1)') parcostname(2),parcostsave(2,1:numit) 
	!do i=1,nepsmove
	!	write(60,'(1a15,2f9.1)') moveshockname,moveshocksave(i,1:numit) 
    !end do
	!do i=1,nepsmove
	!	write(60,'(1a15,2f9.1)') totcostname,totcostsave(i,1:numit) 
    !end do
	!write(60,'(1a15,2f9.1)') parname(16),realpar_save(16,1:numit) 
	!write(60,'(1a15,2f9.1)') parname(52),realpar_save(52,1:numit) 
	!write(60,'(1a15,2f9.1)') parname(66),realpar_save(66,1:numit) 
	write(60,*)
    !write(60,'(tr2,"np",tr1,"np1",tr1,"np2",tr2,"nl",tr1,"neduc",tr2,"nexp ",tr2,"nkid",tr5,"nqs",tr6,"nq",tr6,"nx",tr5,"nxs",tr2,"nepsmv")') !ahumarch1122
	!write(60,'(4i4,3(2x,i4),4i8,2i4)') np,np1,np2,nl,neduc,nexp,nkid,nqs,nq,nx,nxs,nepsmove
	!write(60,'(tr2,"nz",tr2,"nh",tr1,"ncs",tr2,"nc",tr5,"ndata",tr3,"nsimeach",tr6,"nsim",tr2,"ndataobs",tr6,"nmom")') 
	!write(60,'(4i4,5i10)') nz,nh,ncs,nc,ndata,nsimeach,nsim,ndataobs,nmom
	write(60,'("wage grid wg() first row ned/ed men, third row ned/ed fem, fifth row wgt:")')    !     tr6,"m(1)",tr6,"m(2)",tr6,"m(3)",tr6,"m(4)",tr6,"m(5)",tr4,"h(1)",tr4,"h(2)")') 
	write(60,*) wg(:,1,MALES)   !       ,mg(:),hgrid(:)
    write(60,*) wg(:,2,MALES)
	write(60,*) wg(:,1,FEMALES)   !       ,mg(:),hgrid(:)
    write(60,*) wg(:,2,FEMALES)
	write(60,*) wgt(:)    
	write(60,*)
    write(60,'("moveshock grid and its wgts (ppmovesingle and ppmovejoint):")')    
	do iepjoint=1,nepsmove*nepsmove
		indeces=lin2ndim( (/ nepsmove , nepsmove /) , iepjoint )
		iephub=indeces(1) ; iepwfe=indeces(2)
		write(60,'(3I4,4F12.4)') iephub,iepwfe,iepjoint,moveshockmar(iephub,1),moveshockmar(iepwfe,2),moveshockjoint(iephub,iepwfe)%hub,moveshockjoint(iephub,iepwfe)%wfe                  
		write(60,'(3I4,3F12.4)') iephub,iepwfe,iepjoint,ppmovesingle(iephub),ppmovesingle(iepwfe),ppmovejoint(iepjoint)               
	end do 
    !write(60,'("ppso:")')    
    !write(60,*)
    write(60,'("mar grid:")')    !this writing doesn't really work because the second index of mg is trueindex not type
	do i=1,nz
		write(60,'(2F10.2)') marshock(i),ppmarie(i)         
	end do
	!write(60,*) mg(:,2)          
	!write(60,*) mg(:,3)          
	!write(60,*) mg(:,4)  
    !do trueindex=1,ninp
    !    call index2cotyphome(trueindex,i,j,k)
	!    write(60,*) "trueindex,co,typ,home",trueindex,i,j,k
    !    write(60,*) mg(:,trueindex) 
    !end do 
    write(60,*)
    !fnprof(dw0,de,dsex) !(empstat,cur/ofloc,sex)
    !psio(1:2) fnprof(np,5,1)  emp cur m , if m and working and draw curloc: get offer, get laid off, nothing happen
    !psio(3:4) fnprof(np,10,1) emp of m    if m and working and draw ofloc:  get offer, get laid off, nothing happen 
    !psio(5:6) fnprof(np,5,2)   emp cur f if f and working and draw curloc: get offer, get laid off, nothing happen
    !psio(7,8) fnprof(np,10,2) emp of f   if f and working and draw curloc: get offer, get laid off, nothing happen
    !psio(9) fnprof(np1,5,1)  u cur m if m and unemp and draw curloc: get offer,  0 , nothing happen
    !psio(10) fnprof(np1,10,1)  u of m if m and unemp and draw ofloc: get offer,  0 , nothing happen
    !psio(11) fnprof(np1,5,2)  u cur f if f and unemp and draw curloc: get offer,  0 , nothing happen
    !psio(12) fnprof(np1,10,2)  u of f if f and unemp and draw ofloc: get offer,  0 , nothing happen
    write(60,'("offer probabilities (fnprof, governed by psio):")') 
	write(60,'("employed men:")') 
    write(60,'(10x,2x,         4x,"offer",4x,"ldoff",3x,"nthing")') 
    write(60,'("curloc    ",2x,3F9.2)') fnprof(np,5,1)  !prob(nthing) is the residual so only psio(1:2) governs this
    write(60,'("psio(1:2) ",2x,2F9.2)') psio(1:2)
    write(60,'("ofloc     ",2x,3F9.2)') fnprof(np,10,1)  !prob(nthing) is the residual so only psio(3:4) governs this
    write(60,'("psio(3:4) ",2x,2F9.2)') psio(3:4)
	write(60,'("employed fem:")') 
    write(60,'(10x,2x,         4x,"offer",4x,"ldoff",3x,"nthing")') 
    write(60,'("curloc    ",2x,3F9.2)') fnprof(np,5,2)  !prob(nthing) is the residual so only psio(5:6) governs this
    write(60,'("psio(5:6) ",2x,2F9.2)') psio(5:6)
    write(60,'("ofloc     ",2x,3F9.2)') fnprof(np,19,2)  !prob(nthing) is the residual so only psio(7:8) governs this
    write(60,'("psio(7:8) ",2x,2F9.2)') psio(7:8)
    write(60,*) 
    !write(60,*) 
	write(60,'("unemp men:")') 
    write(60,'(10x,2x,         4x,"offer",4x,"ldoff",3x,"nthing")') 
    write(60,'("curloc    ",2x,3F9.2)') fnprof(np1,5,1)  !prob(nthing) is the residual so only psio(9) governs this
    write(60,'("psio(9) ",2x,F9.2)') psio(9)
    write(60,'("ofloc     ",2x,3F9.2)') fnprof(np1,10,1)  !prob(nthing) is the residual so only psio(11) governs this
    write(60,'("psio(10) ",2x,F9.2)') psio(10)
    write(60,'("unemp fem:")') 
    write(60,'(10x,2x,         4x,"offer",4x,"ldoff",3x,"nthing")') 
    write(60,'("curloc    ",2x,3F9.2)') fnprof(np1,5,2)  !prob(nthing) is the residual so only psio(12) governs this
    write(60,'("psio(11) ",2x,F9.2)') psio(11)
    write(60,'("ofloc     ",2x,3F9.2)') fnprof(np1,10,2)  !prob(nthing) is the residual so only psio(13) governs this
    write(60,'("psio(12) ",2x,F9.2)') psio(12)
	write(60,*) 
    !fnprhc(dr,dw) !(curexp,empstat) !fnprhc only depends on psih(1)   
    write(60,'("experience transitions (fnprhc, governed only by psih(1):")') 
    do i=1,nexp
    write(60,'("curexp: ",i4,"empstat: np "," fnprhc(curexp,empstat=np) is: ")') i
    write(60,'(2F9.2)') fnprhc(i,np)
    write(60,'("curexp: ",i4,"empstat: np1 "," fnprhc(curexp,empstat=np1) is: ")') i
    write(60,'(2F9.2)') fnprhc(i,np1)
    end do
    !******************** ahu october2022 **********************************************
    !temprob(1:nl)=fnprloc(1,1) !if origin location is loc1 and homeloc is 2 (cuz index is 2), what is the probability of drawing location 1 
    !realpar(13)=temprob(1)
    !if (mysay==0.and.skriv) THEN
        call index2cotyphome(1,co,typ,homeloc)
        write(60,*) "index,co,typ,homeloc",1,co,typ,homeloc 
        write(60,*) "fnprloc(1,1)",fnprloc(1,1)
    !end if 
    !temprob(1:nl)=fnprloc(1,5) !if origin location is loc1 and homeloc is 2(cuz index is 2), what is the probability of drawing location 2 if index is 1 (so you are at your homeloc)
    !realpar(53)=temprob(2)
    !if (mysay==0.and.skriv) THEN
        !call index2cotyphome(5,co,typ,homeloc)
        !print*, "here is co,typ,homeloc",co,typ,homeloc
        !print*, "fnprloc(1,5)"
        !print*, temprob(:)
        call index2cotyphome(5,co,typ,homeloc)
        write(60,*) "index,co,typ,homeloc",5,co,typ,homeloc 
        write(60,*) "fnprloc(1,5)",fnprloc(1,5)
        call index2cotyphome(9,co,typ,homeloc)
        write(60,*) "index,co,typ,homeloc",9,co,typ,homeloc 
        write(60,*) "fnprloc(1,9)",fnprloc(1,9)
    !end if 
    !******************** ahu october2022 **********************************************
		!fnwge(dg,dtyp,dl,dw,de,dr,dia)	
    	write(60,*)
	    write(60,'("loc,wdraw,wages1-4men,wages1-4fem NED")')   
		do loca=1,nl
		DO WDRAW=1,NP
			wel=wg(wdraw,1,MALES)
			fel=wg(wdraw,1,FEMALES)
			write(60,'(2i4,8F10.0)') loca,wdraw,fnwge(MALES,1,loca,wel,1,1,25), fnwge(MALES,2,loca,wel,1,1,25),fnwge(MALES,3,loca,wel,1,1,25),fnwge(MALES,4,loca,wel,1,1,25)	, fnwge(FEMALES,1,loca,fel,1,1,25), fnwge(FEMALES,2,loca,fel,1,1,25),fnwge(FEMALES,3,loca,fel,1,1,25),fnwge(FEMALES,4,loca,fel,1,1,25)							
		END DO 
		end do 
	    write(60,'("loc,wdraw,wages1-4men,wages1-4fem ED")')   
		do loca=1,nl
		DO WDRAW=1,NP
			wel=wg(wdraw,2,MALES)
			fel=wg(wdraw,2,FEMALES)
			write(60,'(2i4,8F10.0)') loca,wdraw,fnwge(MALES,1,loca,wel,2,1,25), fnwge(MALES,2,loca,wel,2,1,25),fnwge(MALES,3,loca,wel,2,1,25),fnwge(MALES,4,loca,wel,2,1,25)	, fnwge(FEMALES,1,loca,fel,2,1,25), fnwge(FEMALES,2,loca,fel,2,1,25),fnwge(FEMALES,3,loca,fel,2,1,25),fnwge(FEMALES,4,loca,fel,2,1,25)							
		END DO 
		end do 

    	write(60,*)
    	write(60,*)
    	write(60,*)
		write(60,'("typo,wdraw,wages loc1-9  MEN NED ")')   
		do typo=1,ntypp
			DO WDRAW=1,NP
				wel=wg(wdraw,1,MALES)
				fel=wg(wdraw,1,FEMALES)
				write(60,'(2i4,9F10.0)') typo,wdraw,fnwge(MALES,typo,1,wel,1,1,22),fnwge(MALES,typo,2,wel,1,1,22),fnwge(MALES,typo,3,wel,1,1,22),fnwge(MALES,typo,4,wel,1,1,22),fnwge(MALES,typo,5,wel,1,1,22),fnwge(MALES,typo,6,wel,1,1,22),fnwge(MALES,typo,7,wel,1,1,22),fnwge(MALES,typo,8,wel,1,1,22),fnwge(MALES,typo,9,wel,1,1,22)
			END DO 
		end do 
		write(60,'("typo,wdraw,wages loc1-9  MEN ED ")')   
		do typo=1,ntypp
			DO WDRAW=1,NP
				wel=wg(wdraw,1,MALES)
				fel=wg(wdraw,1,FEMALES)
				write(60,'(2i4,9F10.0)') typo,wdraw,fnwge(MALES,typo,1,wel,2,1,22),fnwge(MALES,typo,2,wel,2,1,22),fnwge(MALES,typo,3,wel,2,1,22),fnwge(MALES,typo,4,wel,2,1,22),fnwge(MALES,typo,5,wel,2,1,22),fnwge(MALES,typo,6,wel,2,1,22),fnwge(MALES,typo,7,wel,2,1,22),fnwge(MALES,typo,8,wel,2,1,22),fnwge(MALES,typo,9,wel,2,1,22)
			END DO 
		end do 

		
		write(60,'("typo,wdraw,wages loc1-9  FEM NED ")')   
		do typo=1,ntypp
			DO WDRAW=1,NP
				wel=wg(wdraw,1,MALES)
				fel=wg(wdraw,1,FEMALES)
				write(60,'(2i4,9F10.0)') typo,wdraw,fnwge(FEMALES,typo,1,fel,1,1,22),fnwge(FEMALES,typo,2,fel,1,1,22),fnwge(FEMALES,typo,3,fel,1,1,22),fnwge(FEMALES,typo,4,fel,1,1,22),fnwge(FEMALES,typo,5,fel,1,1,22),fnwge(FEMALES,typo,6,fel,1,1,22),fnwge(FEMALES,typo,7,fel,1,1,22),fnwge(FEMALES,typo,8,fel,1,1,22),fnwge(FEMALES,typo,9,fel,1,1,22)
			END DO 
		end do 
		write(60,'("typo,wdraw,wages loc1-9  FEM ED ")')   
		do typo=1,ntypp
			DO WDRAW=1,NP
				wel=wg(wdraw,1,MALES)
				fel=wg(wdraw,1,FEMALES)
				write(60,'(2i4,9F10.0)') typo,wdraw,fnwge(FEMALES,typo,1,fel,2,1,22),fnwge(FEMALES,typo,2,fel,2,1,22),fnwge(FEMALES,typo,3,fel,2,1,22),fnwge(FEMALES,typo,4,fel,2,1,22),fnwge(FEMALES,typo,5,fel,2,1,22),fnwge(FEMALES,typo,6,fel,2,1,22),fnwge(FEMALES,typo,7,fel,2,1,22),fnwge(FEMALES,typo,8,fel,2,1,22),fnwge(FEMALES,typo,9,fel,2,1,22)
			END DO 
		end do 
    	write(60,*)
		!write(60,'("hgrid:")')    
		!write(60,*) hgrid(:)          
		!write(60,*)
		!write(60,'("wgts:  ",1("move",tr4,"hour",tr4,"wage",tr4),tr1,"rel",tr5,"kid")') 
		!write(60,'(3x,f8.2,5f8.2)') wmove,whour,wwage,wrel,wkid
		!write(60,*) 
		write(60,'(tr2,"groups",tr3,"nhome",tr2,"nhomep",tr2,"onlysingles",tr2,"optimize",tr3,"chkstep",tr5,"skriv",tr1,"numit")') 
		write(60,'(2x,L6,2(3x,I4),7x,L6,3(4x,L6),I6)') groups,nhome,nhomep,onlysingles,optimize,chkstep,skriv,numit
		write(60,'(tr2,"nonneg",tr2,tr10,"eps2",tr11,"eps",tr5,"nonlabinc")') 
		write(60,*) nonneg,eps2,eps,nonlabinc
		write(60,*)
		write(60,'("objective function:")') 
		write(60,*) q_save(:)
		write(60,*) 

    write(60, '(35x,tr7,"sim",tr7,"dat",tr7,"obj",tr7,"dif",5x,tr2,"countsim",tr2,"countdat",tr4,"vardat" )' ) 
    ihead=1   
    do i=1,nmom
		if (headerloc(ihead)==i) then
			write(60,*)
            write(60,*)
			write(60,'(a120)') header(ihead)
			ihead=ihead+1
		end if
		do t=1,numit
            if (condmomcompare) then
                if (t<=2) then
    			    !write(60,'(1i5,". ",1a23,4f10.4,5x,2i10,f10.4)')	i,name(i),momsim_save(i,t),momdat_save(i,t),qcont_save(i,t),qcont_save(i,t)-qcont_save(i,1),cntsim_save(i,t),cntdat_save(i,t),vardat_save(i,t)
                else 
                    !write(60,'(1i5,". ",1a23,4f10.4,5x,2i10,f10.4)')	i,name(i),momsim_save(i,t),momdat_save(i,t),qcont_save(i,t),qcont_save(i,t)-qcont_save(i,3),cntsim_save(i,t),cntdat_save(i,t),vardat_save(i,t)
                end if 
            else 
                write(60,'(2i5,". ",1a55,2F14.4,2f14.4,5x,2i10,f14.4)')	i,t,name(i),momsim_save(i,t),momdat_save(i,t),qcont_save(i,t),qcont_save(i,t)-qcont_save(i,1),cntsim_save(i,t),cntdat_save(i,t),vardat_save(i,t)
                !write(61,'(8i5,2f20.4)')	i,t,mominfo(0:5,i),momsim_save(i,t),momdat_save(i,t)
            end if 
        end do
		!write(60,*)
	end do
	close(60)


	end subroutine writemoments



	subroutine getcorrelation(dumvec)
		real(dp), intent(inout) :: dumvec(6)
		real(dp) :: mux,sdx,muy,sdy,xy,out(6)
			mux = dumvec(1)
			sdx = (dumvec(2) - dumvec(1)**2 ) ** (0.5_dp)    
			muy = dumvec(3)
			sdy = (dumvec(4) - dumvec(3)**2 ) ** (0.5_dp)    
			xy = dumvec(5)
			out(1)=mux
			out(2)=sdx 
			out(3)=muy 
			out(4)=sdy
			out(5)=(xy - mux*muy  )
			out(6)=(xy - mux*muy  ) 	/ (sdx*sdy)
			dumvec(1:6)=out(1:6)
	end subroutine
	
end module objf

!msm_weights(im)=(countdatamom(im)**2)*vardatamom(im)**(-1) ! weighting vector (main diagonal of diaginal weighting matrix)
!making a slight modification to the above. multiplying by (n1/n)**2 rather than n1**2 (where n1 is the countdatamom of first moment here). 
!otherwise the weights blow up, as countdatamom is usually in the thousands.   
!ahu 073112 msm_weights(im)=(real(countdatamom(im))/real(ndata))*vardatamom(im)**(-1) ! weighting vector (main diagonal of diaginal weighting matrix)
!ahu 081912 msm_weights(im)=(1d0/real(ndata))*vardatamom(im)**(-1) !ahu 073112
!ahu 061813 msm_wgt(i)=1.0_dp !(real(countdatmom(i))/real(ndata))*vardatmom(i)**(-1) 
!ahu 061713 msm_wgt(i)= (real(countdatmom(i))/real(ndataobs))  *  ( vardatmom(i)**(-1) ) !ahu 061413. now using ndataobs instead of ndaata. need to check this stuff. still! 

!	if (iwritegen==1 .and. optimize .and. (obj<best.or.chkobj==2)) then
!		write(61,*) 'found a better point', obj
!		do i=1,npars ; write(61,'(f15.10)') parvec(i) ; end do
!		open(unit=66,file='parbest.txt',status='replace')
!		do i=1,npars ; write(66,'(f15.10)') parvec(i) ; end do
!		close(66)
!		best=obj	
!	end if 

		!if (datmom(i)<0.01_dp) then 
		!		datmom(i)=1000.0_dp *datmom(i)
		!		simom(i)=1000.0_dp *simom(i)
		!else if (datmom(i)>=0.01_dp.and.datmom(i)<=0.1_dp) then 
		!		datmom(i)=100.0_dp *datmom(i)
		!		simom(i)=100.0_dp *simom(i)
		!else if (datmom(i)>0.1_dp.and.datmom(i)<=1.0_dp) then 
		!		datmom(i)=10.0_dp *datmom(i)
		!		simom(i)=10.0_dp *simom(i)
		!else if (datmom(i)>4.0_dp) then 
		!		datmom(i)=0.1_dp *datmom(i)
		!		simom(i)=0.1_dp *simom(i)
		!end if 
