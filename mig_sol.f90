! information about the value functions
! get vmsep/vfsep=values of being single for males
! this is the value of remaining single, either because you started the period single and you reject the partner you meet or because you leave a relationship at
! the start of the period. not the value of starting period single, as it does not include the value of potentially entering into a relationship		  			
!get each spouse's consumption transfers for each j alternative using foc to the nb problem to get:
!share(2,cho): spouses' respective consumption shares of total income (pi) for each alternative cho
!these share solutions are conditional on the nb problem being well-defined. otherwise they don't mean anything
!i check for whether the nb problem is well-defined in the getnb subroutine later because i wait for marie
! checknb: check if the nb problem is well-defined meaning that there exists some allocation in the utility possibility set that is mutually beneficial
	!testing(1,1)=1 ; testing(1,2)=2 ; testing(1,3)=3 ; testing(2,1)=4 ; testing(2,2)=5 ; testing(2,3)=6
	!print*, "maxloc(testing,1),maxloc(testing,2)",maxloc(testing,1),maxloc(testing,2),"maxval",maxval(testing(2,:),1)   !,maxval(testing,2)
	!print*, nqs*nqs,nq*nq,count(pps(:,:,1)),count(ppc(:,:))
module sol
	use params
	use share
	use myaz
	implicit none
	real(sp) :: begintime,time(5)   
    !real(dp) :: vsingtest(2) !ag090822 agsept2022
contains		
	subroutine solve
	real(dp) :: vmax(2),val(2),vsum(2),valso(2),probmux(nx,nx,nq),transo(2)
    real(dp), dimension(nxs,nqs) :: vm0_s,vf0_s,vm_marmkt,vf_marmkt,prob_s
	real(dp), dimension(nx,nq) :: vm_c,vf_c,prob
	real(dp), dimension(nepsmove,nxs,nqs,nqs) :: vm_premarmkt,vf_premarmkt
	integer(i4b) :: ia,q0,q,x0,x,z,index,trueindex,dd(12),ed(2),qm,xm,qf,xf,callfrom,nnn,i,n,i0,iepjoint,iephub,iepwfe,indeces(2),iepsingle
    integer(i4b) :: whenmakingdec,trueco,truetyp,truehome
    !print*, 'Here is iter',iter
	begintime=secnds(0.0)
    yaz=.false.
    if (groups) then 
        call cotyphome2index(myindex,myco,mytyp,myhome)
    else 
        myindex=0
    end if 
    do x0=1,nx
    do q0=1,nq
        !do q=1,nq
        !probmuq(x0,q,q0)= ppcq(q,q0,x0)
        !end do 
        do x=1,nx
        probmux(x0,x,q0)=ppcx(x,q0,x0) !ppcq(q,q0,x0)*ppcx(x,q0,x0)
        end do 
    end do 
end do 


    call get_util_w
    if (onlysingles) then ; utilc=pen ; end if 
	insol=.true.
	!prob=0.0_dp
	decm_postdiv = ipen	; decf_postdiv = ipen  ; decm_premarmkt = ipen	; decf_premarmkt = ipen ; dec_mar=ipen ! dum is a dummy for yaz_checknb
    vm_postdiv   = pen ; vf_postdiv  = pen ; vm0_c = pen ; vf0_c = pen
    do trueindex=1,ninp 
        call index2cotyphome(trueindex,trueco,truetyp,truehome)
        if (groups) then 
            index=1
        else 
            index=trueindex
        end if
        ind: if ( (.not.groups) .or.  (myindex==trueindex) ) then !(myco==init(nn)%co.and.mytyp==typsim.and. myhome==init(nn)%hme)  ) then	  
		time(1)=secnds(0.0)
		emaxm_s=0.0_dp	; emaxf_s=0.0_dp ; emaxm_c=0.0_dp	; emaxf_c = 0.0_dp
		do ia=mxa,mna,-1
			!print*, "mysay,ia,trueindex ",iter,mysay,ia,trueindex !ahu 030622
			vm_premarmkt  = pen ; vf_premarmkt = pen	
			time(1)=secnds(0.0)
			if (ia==mxa) then 
				vm0_s               = utils(1,:,:,ia,trueindex) !+ wsnet(1,:,:,ia,trueindex)		! v0: value function without movecost 
				vf0_s               = utils(2,:,:,ia,trueindex) !+ wsnet(2,:,:,ia,trueindex)		! v0: value function without movecost
				vm0_c(:,:,:,ia,index) = utilc(1,:,:,ia,trueindex)	! v0: value function without movecost, umar, consumption
				vf0_c(:,:,:,ia,index) = utilc(2,:,:,ia,trueindex)	! v0: value function without movecost, umar, consumption
            else 
				vm0_s               = utils(1,:,:,ia,trueindex) + delta * emaxm_s(:,:,ia+1)  !+ wsnet(1,:,:,ia,trueindex)			! v0: value function without movecost 
				vf0_s               = utils(2,:,:,ia,trueindex) + delta * emaxf_s(:,:,ia+1)  !+ wsnet(2,:,:,ia,trueindex)			! v0: value function without movecost 
				vm0_c(:,:,:,ia,index) = utilc(1,:,:,ia,trueindex)	+ delta * emaxm_c(:,:,:,ia+1)	! v0: value function without movecost, umar, consumption
				vf0_c(:,:,:,ia,index) = utilc(2,:,:,ia,trueindex)	+ delta * emaxf_c(:,:,:,ia+1)	! v0: value function without movecost, umar, consumption
			end if 		
            
			
			!if (skriv.and.trueindex==1.and.(ia<=19.or.ia==45)) then ; yaz=.true. ; else ; yaz=.false. ; end if !ahu 0327 trueindex==2
            dd = -1 ; dd(1:2) = (/ia,trueindex/) 
            !vm0_s and vf0_s is the value of just being a location/emp BEFORE any shocks AS a single person. 
            !getdec_s: takes v0_s and decides the optimal loc/emp AFTER shocks 
            !getdec_s: inuput is v0_s and output is vm_postdiv and vf_postdiv (the value AFTER shocks) if this is POSTDIV (i.e. whenmakingdec=1)
            !if postdiv is TRUE: then getdec_s calculates optimal decision and value at optimal loc/emp for a single person who JUST GOT DIVORCED
            !the only difference from the other case is that the movecost and moveshock are different 
            whenmakingdec=1 !1 means this is postdiv
            dd(7)=1 ; call getdec_s(dd,vm0_s(:,:),decm_postdiv(:,:,:,:,ia,index),vm_postdiv(:,:,:,:,ia,index),whenmakingdec) ; yaz=.false.
			dd(7)=2 ; call getdec_s(dd,vf0_s(:,:),decf_postdiv(:,:,:,:,ia,index),vf_postdiv(:,:,:,:,ia,index),whenmakingdec)	
            !dd(7)=1 ; call getdec_s(dd,vm0_s,decm0_s(:,:,:,:,ia,index),vm(:,:,:,:,ia,index),valdiv) ; yaz=.false.
			!dd(7)=2 ; call getdec_s(dd,vf0_s,decf0_s(:,:,:,:,ia,index),vf(:,:,:,:,ia,index),valdiv)	
			if (onlysingles) then 	
                whenmakingdec=2 !2 means this is pre marriage market 
                dd(7)=1 ; call getdec_s(dd,vm0_s(:,:),decm_premarmkt(:,:,:,:,ia,index),vm_premarmkt(:,:,:,:),whenmakingdec) ; yaz=.false.
                dd(7)=2 ; call getdec_s(dd,vf0_s(:,:),decf_premarmkt(:,:,:,:,ia,index),vf_premarmkt(:,:,:,:),whenmakingdec)	
            else 
                !for loop opt 
                do q=1,nq
                    do x=1,nx
                        vm0ctemp(q,x)=vm0_c(x,q,ia,index) 
                        vf0ctemp(q,x)=vf0_c(x,q,ia,index) 
                        wmctemp(q,x,ia,trueindex)=wcnet(1,x,q,ia,trueindex)
                        wfctemp(q,x,ia,trueindex)=wcnet(2,x,q,ia,trueindex)
                    end do 
                end do 

				do q0=1,nq	
					vm_c=pen ; vf_c=pen
					prob=0.0_dp
					do q=1,nq
                    do x=1,nx	
						
                            !ahumarch1122 
                            !if (skriv.and.(ia==50.or.ia==48.or.ia==18).and.q0==4.and.q==92.and.x==19.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if !ahu030822
                            if (skriv.and.(ia==50.or.ia==48.or.ia==18).and.(q0>=4.and.q0<=8).and.(q>=92.and.q<=108).and.x==19.and.trueindex==1) then ; yaz=.false. ; else ; yaz=.false. ; end if 
                            !ahu030822 Make sure you don't have iep in the if statements here

							if ( ppc(q,q0) ) then
                                valso=pen
								val=0.0_dp
                                marshocks: do z=1,nz
                                    moveshocks: do iepjoint=1,nepsmove*nepsmove
                                        if (   (JOINTMARSHOCK .and.ppmovejoint(iepjoint)>0.)   .OR.   (.NOT.JOINTMARSHOCK) ) then 
                                            !callfrom is for telling getdec and therefore yaz_getdec where we are when yaz_getdec is called from within getdec
                                            !(if skriv and yaz true) when callfrom=40, then within getdec I call yaz_getdec to write in 200 the altspecific value functions and bargaining stuff
                                            !in getdec when callfrom is 40 or 80, then dd(8) altj and dd(9) altq are filled in within the choice loop
                                            callfrom=40 
                                            dd = (/ia,trueindex,q,x,z,q0,callfrom,-1,-1,-1,iepjoint,-1 /) 	! (/ ia,index,q,x,iepjoint (used to be z),q0,callfrom,altj/jmax,altq/qmax,relmax,iepsingle,gender /)  	
                                            call getdec(dd,vmax,valso,transo)	!valso is set to pen. getdec here gets the decision of married household. outside option is vm_postdiv and vf_postdiv. note that valso is supposed to tbe the outside option input but I set that value inside getdec instead (at least when calling it from here. other times i.e. when callfrom is 50, I set valso to the outside option value. for no reason. ! dd is inout and dd(8:10) tells you jmax,qmax,relmax 
                                            if (yaz) then ; call yaz_decision(dd,vmax) ; end if	!callfrom tells yaz_decision where we are 									
                                            val=val+ppmarie(z)*ppmovejoint(iepjoint)*vmax !*mgt(z)
                                        end if 
                                    end do moveshocks
                                end do marshocks
								vm_c(x,q)=val(1) ; vf_c(x,q)=val(2)
							end if !pc0(q0)=.true. i.e. w(1) or w(2) are not np2 and pc
						end do !x
					end do !q
					!if (icheck_eqvmvf) then
					!	call check_eqvmvf(dd,vm_c,vf_c) 
					!end if 


                        emaxm_c(:,q0,ia)=0.0_dp
						emaxf_c(:,q0,ia)=0.0_dp
                        if ( maxval(qq2w(:,q0)) <= np1 ) then !state variable part of the q space i.e. w <= np1
						    !prob=matmul( reshape(ppcq(:,q0,x0),(/nq,1/)) , reshape(ppcx(:,q0,x0),(/1,nx/)) )
						    !emaxm_c(x0,q0,ia)=sum(prob*vm_c)
						    !emaxf_c(x0,q0,ia)=sum(prob*vf_c)                    
                            do q=1,nq
                            do x=1,nx
                                if (terminalval.and.ia==MXA) then 
                                    valso=0.0_dp
                                    do nnn=0,ntermval
                                        valso(1) = valso(1)  + (delta**nnn)*vm_c(x,q)
                                        valso(2) = valso(2)  + (delta**nnn)*vf_c(x,q)
                                    end do
                                else 
                                    valso(1)=vm_c(x,q)
                                    valso(2)=vf_c(x,q)
                                end if                     
                            do x0=1,nx
                                emaxm_c(x0,q0,ia)=emaxm_c(x0,q0,ia)+ppcq(q,q0,trueindex)*probmux(x0,x,q0)*valso(1)
                                emaxf_c(x0,q0,ia)=emaxf_c(x0,q0,ia)+ppcq(q,q0,trueindex)*probmux(x0,x,q0)*valso(2)
                            end do 
                            end do 
                            end do	 				
                        end if !state spce check 
                end do !q0
				!if (Trueindex==1) print*, "I even came here and nothing happened"

                !timeline: when single, you first decide on loc and work. and then when the day comes to an end, you go to the marriage market in whatever location you chose for that period.
				vm_marmkt=0.0_dp ; vf_marmkt=0.0_dp 
				do qm=1,nqs	
                do xm=1,nxs					
                    ed(1)=x2e(xm)   !ahu summer18 051418: adding ed dimension to nonlabinc
						!no need for this ahu040518 vsingle=0.0_dp
						do qf=1,nqs
                        do xf=1,nxs			
                            ed(2)=x2e(xf)   !ahu summer18 051418: adding ed dimension to nonlabinc
							
                                !i = dd(8) ; n=dd(9) 
                                valso(1) = vm0_s(xm,qm) !outside option when making decision in marriage market
                                !i = dd(10) ; n=dd(11) 
                                valso(2) = vf0_s(xf,qf) !outsdie option when making decisoin in marriage market
                                !vsingtest(1:2)=valso(1:2) !ag090822
                    
                                if (skriv.and.(ia==18).and.(qm==18).and.xm==1) then ; yaz=.true. ; else ; yaz=.false. ; end if     !ahu 0327 trueindex==2								
                                !no need for this ahu040518 sex: do g=1,2
                                !no need for this ahu040518 if (g==1) then ; qm=qr ; qf=qp ; xm=xr ; xf=xp ; end if ! determine the match and the joint q from the combination
                                !no need for this ahu040518 if (g==2) then ; qm=qp ; qf=qr ; xm=xp ; xf=xr ; end if ! determine the match and the joint q from the combination
                                if (ppmeetq(qm,qf)>0.0_dp ) then !.and. ppmeetx(xm,xf)>0.0_dp) then 
                                    q0=-1	; q = q2qq(qm,qf)  ; x = x2xx(xm,xf) 
                                
                                    if (q2w(qm)>np1.or.q2w(qf)>np1.or.q==0.or.x==0) stop   !if (ps0(qm) .and. ps0(qf) .and. q > 0 .and. x >0) then ! all w should be <= np1	and q,x should exist
                                    do z=1,nz
                                        callfrom=50
                                        dd = (/ ia,trueindex,q,x,z,q0,callfrom,-1,q,-1,-1,-1 /) ! (/ ia,index,q,x,iepjoint (used to be z),q0,callfrom,altj,altq,rel,iepsingle,gender /)  	                                        
                                        !callfrom is for telling getdec and therefore yaz_getdec where we are when yaz_getdec is called from within getdec
                                        !(if skriv and yaz true) when callfrom=50, then within getdec I call yaz_getdec to write in 200 the altspecific value functions and bargaining stuff
                                        !in getdec when callfrom is 40 or 80, then dd(8) altj and dd(9) altq are filled in within the choice loop
                                        !in getdec when callfrom is 50, then dd(8) altj is -1 and dd(9) altq is just q (same as dd(3) ) 
                                                            ! g and q0 is just -1 here (/ ia,index,q,x,iepjoint (used to be z),q0,gender,j,q,rel,iepsingle /) 
                                                            ! when in the marriage market, the index that corresponds 
                                                            ! to altq (i.e. i in yaz for example or dim 9 in dd) is just set to q since there is no altq in marriage market
                                        call getdec(dd,vmax,valso,transo) !decision whether to get married. with out option just v0_s (value of being at a loc/emp BEFORE any shocks)
                                        dec_mar(z,x,q,ia,index)=dd(10) 
                                        val=vmax(1:2)         
                                        if (yaz) then ; call yaz_decision(dd,vmax) ; end if	!callfrom tells yaz_decision where we are 				

                                        !if (  icheck_eqvmvf.and.qr==qp.and.xr==xp.and.  (abs(vec(1)-vec(2)) > eps .or. abs(vsum(1)-vsum(2)) > eps) ) then 
                                        !	print*, "vm not eq to vf!", ia,qr,xr,vec(1:4),vsum
                                        !	stop 
                                        !end if 
                                        !no need for this ahu040518 vsingle(g,qp,xp,z)=val(g)	
                                        !vm_marmkt(xm,qm)=vm_marmkt(xm,qm) + ppmarie(z) * ppmeetq(qm,qf) * ppmeetx(xm,xf) * val(1)
                                        !vf_marmkt(xf,qf)=vf_marmkt(xf,qf) + ppmarie(z) * ppmeetq(qm,qf) * ppmeetx(xm,xf) * val(2)
                                        vm_marmkt(xm,qm)=vm_marmkt(xm,qm) + ppmarie(z) * ppmeetq(qm,qf) * probmeetx(XF,XM,TRUETYP,MALES) * val(1)
                                        vf_marmkt(xf,qf)=vf_marmkt(xf,qf) + ppmarie(z) * ppmeetq(qm,qf) * probmeetx(XM,XF,TRUETYP,FEMALES) * val(2)

                                    end do 
                                end if 
							end do 
						end do 							
						!moving this after the loop ahu040518 vm_s(xr,qr) = pmeet*vm_s(xr,qr) + (1.0_dp-pmeet)*vm0_s(xr,qr)
						!moving this after the loop ahu040518 vf_s(xr,qr) = pmeet*vf_s(xr,qr) + (1.0_dp-pmeet)*vf0_s(xr,qr)
						!if (  icheck_eqvmvf.and.abs(vm_s(xr,qr)-vf_s(xr,qr)) > eps ) then 
						!	print*, "vm not equal to vf!", ia,qr,xr,vm_s(xr,qr),vf_s(xr,qr)
						!	stop 
						!end if 
					end do 
				end do
				vm_marmkt(:,:) = pmeet(ia)*vm_marmkt(:,:) + (1.0_dp-pmeet(ia))*vm0_s(:,:) !ahu 040518
				vf_marmkt(:,:) = pmeet(ia)*vf_marmkt(:,:) + (1.0_dp-pmeet(ia))*vf0_s(:,:) !ahu 040518
                !vm_s=vm0_s
                !vf_s=vf0_s
                !dec_mar=0
				yaz=.false.
				if (skriv.and.(ia==mxa.or.ia==29).and.trueindex==1) yaz=.false.
				dd = -1 ; dd(1:2) = (/ia,trueindex/) 							! (/ ia,index,q,x,iepjoint (used to be z),q0,g,j,altq,rel,iepsingle /)
                whenmakingdec=2 !2 indicates we are looking at decisions BEFORE the marriage market
                !getdec_s: now takes vm_s and vf_s as input which are vm0_s and vf0_s AFTER taking into account the marriage market
				dd(7)=1 ; call getdec_s(dd,vm_marmkt,decm_premarmkt(:,:,:,:,ia,index),vm_premarmkt(:,:,:,:) ,whenmakingdec    ) 
				dd(7)=2 ; call getdec_s(dd,vf_marmkt,decf_premarmkt(:,:,:,:,ia,index),vf_premarmkt(:,:,:,:) ,whenmakingdec    ) 			
			end if 	! only singles or no 

			do q0=1,nqs
                if (q2w(q0)<=np1) then !state variable part of the q space i.e. w <= np1
                do x0=1,nxs			
						do q=1,nqs
                            do x=1,nxs 
                                do iepsingle=1,nepsmove
                                    !prob_s=matmul( reshape(ppsq(:,q0,x0,1),(/nqs,1/)) , reshape(ppsx(:,q0,x0),(/1,nxs/)) )
						            !emaxm_s(q0,x0,ia)	= sum( prob_s * vmr(:,:,q0) )		
                                    !prob_s=matmul( reshape(ppsq(:,q0,x0,2),(/nqs,1/)) , reshape(ppsx(:,q0,x0),(/1,nxs/)) )
						            !emaxf_s(q0,x0,ia)	= sum( probf_s * vfr(:,:,:,q0) )
                                    if (terminalval.and.ia==MXA) then 
                                        valso=0.0_dp
                                        do nnn=0,ntermval
                                            valso(1) = valso(1)  + (delta**nnn)*vm_premarmkt(iepsingle,x,q,q0) 
                                            valso(2) = valso(2)  + (delta**nnn)*vf_premarmkt(iepsingle,x,q,q0) 
                                        end do
                                    else 
                                        valso(1)=vm_premarmkt(iepsingle,x,q,q0) 
                                        valso(2)=vf_premarmkt(iepsingle,x,q,q0) 
                                    end if                     
                                    emaxm_s(x0,q0,ia)	= emaxm_s(x0,q0,ia) + ppmovesingle(iepsingle) * ppsq(q,q0,trueindex,1) * ppsx(x,q0,x0) * valso(1) 
                                    emaxf_s(x0,q0,ia)	= emaxf_s(x0,q0,ia) + ppmovesingle(iepsingle) * ppsq(q,q0,trueindex,2) * ppsx(x,q0,x0) * valso(2) 	
                                end do 
                            end do 
                        end do
				end do 
                    end if !state variable part of the q space i.e. w <= np1
			end do 

			yaz=.false.			
			!if (skriv) call yaz1(index,ia) !ahu 0327. this was trueindex and corrected it to index.  
		end do ! age 
        end if ind
	end do !index
	!if ( (.not.chkstep).and.(.not.optimize) ) print*, "Finished solve and mysay is: ", mysay
	if (onlysingles) then
		dec_mar=0
	end if 
	insol=.false.
	end subroutine solve





    !ag090122 agsept2022 Changed getdec_c to getdec in order to include mar market decisinos in getdec as well
    !                     so that everything is all in one same place and if decisino protocol changes are made they have to be made to only one place
	!                     OLDER VERSIONS OF GETDEC_C CAN BE FOUND IN MIGSOLGETDEC.F90 FILE
    One of the main goals in this paper is to show that gender wage differences do not exist in a vacuum. In other words, what they imply about what
     happens within households is an important aspect of the problem. As I highlight in the introduction of the first draft, the welfare implications ofgende rwage differences 

    subroutine getdecR(dd,vmax,vsing,transfers) 
        integer(i4b), intent(inout) :: dd(:)
        real(dp), intent(out) :: vmax(2),transfers(2) 
        real(dp), intent(in) :: vsing(2) 
        real(dp) :: vecj(5,nc),vec(5),surplusj(nc),surplus,vdif(2),mcost(2),asum,yazvec(5),vmaxcheck(2)
        integer(i4b) :: ia,index,q,x,z,q0,g,jmax,qmax,relmax,i,i0,n,trueindex,j,de(1),ed(2),locch,loc0,callfrom
        logical :: haveenough(nc),haveenoughtomakeindiff,haveenoughtomakeindiff_alt,intsol(3),pc(2),pc_alt(2),haveenoughforNB(nc),haveenoughNB
        integer(i4b) :: trueco,truetyp,truehome,iepjoint,iephub,iepwfe,indeces(2),tempmax
        callfrom=dd(7)
        ia=dd(1) 
        trueindex=dd(2) 
        if (groups) then 
            index=1
        else 
            index=trueindex
        end if
        call index2cotyphome(trueindex,trueco,truetyp,truehome)
        q=dd(3) 
        x=dd(4) ; ed(:)=xx2e(:,x)   
        z=dd(5)
        q0=dd(6)  
        vec=pen				! initialize
        jmax=-1				! initialize
        qmax=-1				! initialize
        relmax=-1           ! initialize
        vmax=pen			! initialize
        transfers=pen 
        if (Callfrom==40.or.callfrom==80) then !calling from sol/couples/getdec (40) OR simul/couples/getdec (80) 

            vec=pen
            loc0=qq2l(1,q0)     !ahumarch2022 ahu032022
            do locch=1,nl ; locdif(locch)=one( locch /= loc0 ) ; end do !locations array
            iepjoint=dd(11) !note that dd(5) used to be z and iepsmove used to be dd(11) before
            indeces=lin2ndim( (/ nepsmove , nepsmove /) , iepjoint )
            iephub=indeces(1) ; iepwfe=indeces(2)  
            array_m(:,loc0) = utilmar(MALES,x,z,trueindex,ia)   + locdif(:) * (movecost(xx2x(1,x),qq2q(1,q0),trueindex,MALES)    + moveshockmar(iephub,1)    + distance(:,loc0)*cstadjacent ) 
            array_f(:,loc0) = utilmar(FEMALES,x,z,trueindex,ia) + locdif(:) * (movecost(xx2x(2,x),qq2q(2,q0),trueindex,FEMALES)  + moveshockmar(iepwfe,2)    + distance(:,loc0)*cstadjacent )
            
            i = qq2q(1,q) ; n=xx2x(1,x) ; i0 = qq2q(1,q0) 
            vec(1) = vm_postdiv(iephub,n,i,i0,ia,index) + diveduc(x2e(n)) !+ divpenalty * one( q2l(qchoicesingle) /= loc0 )
            i = qq2q(2,q) ; n=xx2x(2,x) ; i0 = qq2q(2,q0) 	
            vec(2) = vf_postdiv(iepwfe,n,i,i0,ia,index) + diveduc(x2e(n)) !+ divpenalty * one( q2l(qchoicesingle) /= loc0 )
            
            vecj=pen ; surplusj=pen ; haveenough=.FALSE. ; yazvec=pen 
            do ll=1,nl
                muii=mui0
                GO TO 358

                
                
                
                
                
                
                
                
                cc=0
                358 choice1: do j=1,nc	
                    i = ch(j,q,q0)	!alternative q
                    if (i>0 ) then		
                        locch=qq2l(1,i)     !ahumarch2122 ahu032122
                        if (ll==locch) then
                            cc=cc+1  
                        if (callfrom==40) then !calling from sol
                            vmchoice(cc) = vm0ctemp(i,muii,x)       + array_m(locch,loc0)
                            vfchoice(cc) = vf0ctemp(i,muii,x)       + array_f(locch,loc0)
                        else if (callfrom==80) then !calling from simulation
                            vmchoice(cc) = vm0_c(x,i,muii,ia,index) + array_m(locch,loc0)
                            vfchoice(cc) = vf0_c(x,i,muii,ia,index) + array_f(locch,loc0)
                        end if 
                    end if      
                end do choice1


                pmax(1)=MAXLOC(muvec(muii)*vmchoice(:)+(1-muvec(muii))*vfchoice(:))
                qmax(ll)=pmax(1)
                vmax(1,ll)=vmchoice(pmax(1))
                vmax(2,ll)=vfchoice(pmax(1))
                pc(1:2)	= ( vmax(1:2,ll) - vec(1:2) + eps >= 0.0_dp )
            
                if (.not.pc(1) .or. .not.pc(2) ) then 
                muii=mui0
                DO WHILE ((vmax(1,ll)<vec(1)) .AND. (muii<=Nmu))
                    muii=muii+1
                    cc=0
                    choice2: do j=1,nc	
                    i = ch(j,q,q0)	!alternative q
                    if (i>0 ) then		
                        locch=qq2l(1,i)     !ahumarch2122 ahu032122
                        if (ll==locch) then
                            cc=cc+1  
                        if (callfrom==40) then !calling from sol
                            vmchoice(cc) = vm0ctemp(i,muii,x)       + array_m(locch,loc0)
                            vfchoice(cc) = vf0ctemp(i,muii,x)       + array_f(locch,loc0)
                        else if (callfrom==80) then !calling from simulation
                            vmchoice(cc) = vm0_c(x,i,muii,ia,index) + array_m(locch,loc0)
                            vfchoice(cc) = vf0_c(x,i,muii,ia,index) + array_f(locch,loc0)
                        end if 
                    end if      
                    end do choice2
                    pmax(1)=MAXLOC(muvec(muii)*vmchoice(:)+(1-muvec(muii))*vfchoice(:))
                    qmax(ll)=pmax(1)
                    vmax(1,ll)=vmchoice(pmax(1))
                    vmax(2,ll)=vfchoice(pmax(1))
                    pc(1:2)	= ( vmax(1:2,ll) - vec(1:2) + eps >= 0.0_dp )
                END DO !while vmmaxloc is still less than outside value
                end do !muivec


				if locch=1,nl
                    call sort1(vmuim)
                    k=locate(vmuim(:,locch),vec(1)) 
                    if (iam==0) write(64,'("first try ",i6,3f14.4,i6)') j,q1val(1),q1val(nj2),(lb+j*0.1_sp)*real(qval),k
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


                end do !mui 
            end if !participation constraint check                   
            
            
            
            
            if ( .not.pc(1)   ) THEN
                ! only male's PC binds, increase mu until doesn't
                DO WHILE ((vmax(1,muii)<vec(1)) .AND. (muii<=Nmu))
                    muii=muii+1
                ENDDO
                IF ((muii>Nmu) .OR. (vmax(2,muii)<vec(2))  ) THEN
                    ! can't find mu that satisfies both PC -> separate
                    Vcm(wagexj,thetai,nkidj,prevkidsi,mui)=vmout+cohabpenalty(myco)+kidpenalty*nkidj
                    Vcf(wagexj,thetai,nkidj,prevkidsi,mui)=vfout+cohabpenalty(myco)+kidpenalty*nkidj
                    IF (age<=maxagesim) THEN
                        muc(mui,thetai,prevkidsi,nkidj,wagexj,age)=0
                    ENDIF
                    if (indwriteoutput.eq.1) then
                        WRITE(2,'(16I4)') 3,mui,thetai,nkidm,nkidf,nkidj,com,wagem,typm,dlawm,cof,wagef,typf,dlawf,age,0
                    end if 



            IF (.NOT.UPFRONTTRANSFERS) THEN
                de=maxloc(surplusj,MASK=(haveenough.AND.haveenoughforNB)) 
            ELSE 
                de=maxloc(surplusj,MASK=(surplusj>=0.0_dp)) 
            END IF 
            if (de(1)>0) then 
                jmax=de(1) ; qmax=ch(jmax,q,q0) ; relmax=1
                surplus=surplusj(jmax)                
                vdif(1:2)=vecj(1:2,jmax)    
                vec(3:5)=vecj(3:5,jmax)
                transfers(1) = alf * surplus -  vdif(1)                                                           
                transfers(2) = (1.0_dp-alf) * (vec(5) + vdif(1) ) - alf*vdif(2)  
                if (    abs(transfers(1)+transfers(2) - vec(5) ) >  0.00001_dp  ) then 
                    print*, "Stop it right now",transfers(1:2),vec(5)
                end if 
                vmaxcheck(1:2)=vec(3:4)+transfers(1:2)
                vmax(1:2)=vec(1:2)+alf*surplus
                if (    abs(vmaxcheck(1)-vmax(1))>0.00001_dp   .or.  abs(vmaxcheck(2)-vmax(2))>0.00001_dp ) then 
                    print*, "Stop it right now"
                end if 
                !if (skriv) then 
                    !myefficiency(:)%ve1=vecj(1,:)    
                    !myefficiency(:)%ve2=vecj(2,:)    
                    !myefficiency(:)%ve3=vecj(3,:)    
                    !myefficiency(:)%ve4=vecj(4,:)    
                    !myefficiency(:)%ve5=vecj(5,:)    
                    !myefficiency(:)%sur=surplusj(:) 
                    !myefficiency(:)%haveenough=haveenough(:)
                    !myefficiency(:)%haveenoughforNB=haveenoughforNB(:)
                    de=maxloc(surplusj,MASK=(haveenough.AND.haveenoughforNB))
                    tempmax=de(1)
                    if (de(1)>0) then 
                        myefficiency%maxA=ch(tempmax,q,q0) 
                    else 
                        myefficiency%maxA=0
                    end if
                    de=maxloc(surplusj,MASK=(surplusj>=0.0_dp)) 
                    tempmax=de(1)
                    if (de(1)>0) then 
                        myefficiency%maxB=ch(tempmax,q,q0) 
                    else 
                        myefficiency%maxB=0
                    end if
                    myefficiency%tra1=transfers(1) 
                    myefficiency%tra2=transfers(2) 
                !end if
                !SAVETIME IF (.NOT.UPFRONTTRANSFERS) THEN
                !SAVETIME     if (minval(transfers)<-0.00001 .or. abs(vmaxcheck(1)-vmax(1))>0.00001_dp .or.  abs(vmaxcheck(2)-vmax(2))>0.00001_dp ) then 
                !SAVETIME         print*, "Seriously though", transfers,vmaxcheck,vmax 
                !SAVETIME         stop
                !SAVETIME     end if 
                !SAVETIME     if ( vec(5) + eps < abs(vdif(1)-vdif(2)) )  then 
                !SAVETIME         print*, "why here should not be here like at all ", callfrom,vec(5),vdif(1)-vdif(2)
                !SAVETIME         stop
                !SAVETIME     end if
                !SAVETIME END IF !NO UPFRONTTRANSFERS.
                !SAVETIME if (yaz) then ; write(200,*)  ; write(200,*) ; write(400,*) ; write(400,*) ; end if
                !SAVETIME if (yaz.and.callfrom==40) write(200,*) "in sol/getdec couples de(1)>0: jmax,qmax,relmax,surplusJ(jmax),vecj(1:5,jmax) NOW GO TO 2013"
                !SAVETIME if (yaz.and.callfrom==40) write(200,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplusj(jmax),vecj(1:5,jmax)
                !SAVETIME if (yaz.and.callfrom==80) write(400,*) "in sim/getdec couples de(1)>0: jmax,qmax,relmax,surplusJ(jmax),vecj(1:5,jmax) NOW GO TO 2013"
                !SAVETIME if (yaz.and.callfrom==80) write(400,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplusj(jmax),vecj(1:5,jmax)
                GO TO 2017
            else 
                jmax=-1 ; qmax=-1  ; relmax=0 ; vmax(1:2)=vec(1:2) ; transfers=pen
                !SAVETIME if (yaz) then ; write(200,*)  ; write(200,*) ; write(400,*) ; write(400,*) ; end if
                !SAVETIME if (yaz.and.callfrom==40) write(200,*) "in sol/getdec couples de(1)<=0: jmax,qmax,relmax,de(1) NOW GO TO 2017"
                !SAVETIME if (yaz.and.callfrom==40) write(200,*) jmax,qmax,relmax,de(1)
                !SAVETIME if (yaz.and.callfrom==80) write(400,*) "in sim/getdec couples de(1)<=0: jmax,qmax,relmax,de(1) NOW GO TO 2017"
                !SAVETIME if (yaz.and.callfrom==80) write(400,*) jmax,qmax,relmax,de(1)
                GO TO 2017
            end if 
        end if 

        if (callfrom==50) then !calling from sol/marmkt (this getdec instance is only called in sol (not sim because I save all decisions to dec array)
            if (dd(11).ne.-1) then ; print*, "iepsingle is not -1 in call to getdec and callfrom is from sol/marmkt!!" ; stop ; end if 
            vec(1:2)=vsing(1:2) !ag090822 agsept2022 vsing(1:2) 
            vec(3) = vm0ctemp(q,x) + utilmar(MALES,x,z,trueindex,ia) 
            vec(4) = vf0ctemp(q,x) + utilmar(FEMALES,x,z,trueindex,ia) 
            vec(5) = wmctemp(q,x,ia,trueindex) +  wfctemp(q,x,ia,trueindex) + nonlabinc(ed(1)) + nonlabinc(ed(2))    + ubenefit_c(q)                                                
            vdif(1)=vec(3)-vec(1)   !ahumarch1522 ahu031522 adding cornersol
            vdif(2)=vec(4)-vec(2)   !ahumarch1522 ahu031522 adding cornersol
            !ahumarch2122 ahu032122 replacing this with vdif for saving time surplus = vec(5) + vec(3) - vec(1) + vec(4) - vec(2)  
            surplus = vec(5) + sum(vdif(1:2))  !ahumarch2122 ahu032122 replacing with vdif for saving time
            pc(1:2)	= ( vdif + eps >= 0.0_dp )	!pc(1:2)    = ( vec(3:4) - vec(1:2) >= 0.0_dp )						        
            !ahu032122 ahumarch2122 commenting out to save time pc_alt(1:2)=( vdif >= 0.0_dp )	
            asum = sum(  one(.not. pc)  *   abs( vdif )   ) 
            haveenoughtomakeindiff=(  vec(5) + eps - asum  >= 0.0_dp  ) ; haveenoughNB= (vec(5) +eps >= abs(vdif(1)-vdif(2)) )
            IF (.NOT.UPFRONTTRANSFERS) THEN
                de(1)=one( surplus+eps > 0.0_dp .and. haveenoughtomakeindiff .and. haveenoughNB)             
            ELSE 
                de(1)=one( surplus+eps > 0.0_dp )             
            END IF 
            !SAVETIME if (yaz) then 
                !callfrom=50 is sol/marmkt/getdec and before calling getdec, dd(9) is just set to q i.e. equal to dd(3) and dd dd(8)=-1 ; dd(9)=-1 !to tell yaz_getdec which alternative is being evaluated (where j is the choice and i is the corresponding q)
                !but when calling from singles mar market there is no j or i so set equal to -1 
                !SAVETIME call yaz_getdec(dd,vec(1:5),surplus,pc(1:2),asum,haveenoughtomakeindiff) !called from sol/getdec/marmkt 
            !SAVETIME end if 
            if (de(1)>0) then 
                jmax=-1 ; qmax=-1 ; relmax=1
                transfers(1) = alf * surplus -  vdif(1)                                                           
                transfers(2) = (1.0_dp-alf) * (vec(5) + vdif(1) ) - alf*vdif(2)  
                vmaxcheck(1:2)=vec(3:4)+transfers(1:2)
                vmax(1:2)=vec(1:2)+alf*surplus
                if (    abs(vmaxcheck(1)-vmax(1))>0.00001_dp   .or.  abs(vmaxcheck(2)-vmax(2))>0.00001_dp ) then 
                    print*, "Stop it right now"
                end if 
                if (    abs(transfers(1)+transfers(2) - vec(5) ) >  0.00001_dp  ) then 
                    print*, "Stop it right now this is sol mar market ",transfers(1:2),vec(5)
                end if 
                !SAVETIME IF (.NOT.UPFRONTTRANSFERS) THEN
                !SAVETIME     if (minval(transfers)<-0.00001 .or. abs(vmaxcheck(1)-vmax(1))>0.00001_dp .or.  abs(vmaxcheck(2)-vmax(2))>0.00001_dp ) then 
                !SAVETIME         print*, "Seriously though", transfers,vmaxcheck,vmax 
                !SAVETIME         stop
                !SAVETIME     end if 
                !SAVETIME     if ( vec(5) + eps < abs(vdif(1)-vdif(2)) )  then 
                !SAVETIME         print*, "why here should not be here like at all ", callfrom
                !SAVETIME         stop
                !SAVETIME     end if
                !SAVETIME END IF !UPFRONT transfers   
                !SAVETIME if (yaz) write(200,*) "in sol/getdec marmkt de(1)>0: jmax,qmax,relmax,surplus,vec(1:5) NOW GO TO 2013"
                !SAVETIME if (yaz) write(200,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplus,vec(1:5) 
                GO TO 2017 !have to calculate vmax and transfers before going to 2017 
            else 
                jmax=-1 ; qmax=-1  ; relmax=0 ; vmax(1:2)=vec(1:2) ; transfers=pen 
                !SAVETIME if (yaz) write(200,*) "in sol/getdec marmkt de(1)<=0: jmax,qmax,relmax,surplus,vec(1:5) NOW GO TO 2017"
                !SAVETIME if (yaz) write(200,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplus,vec(1:5)                 
                GO TO 2017 ! go directly to 2017
            end if 
        end if 

        !CORNER SOLUTIONS NOW OBSOLETE
        !else if ( vec(5) <= vdif(1)-vdif(2) + eps  )  then 
        !    transfers(1)=0.0_dp                         !ahumarch1522 adding cornersol
        !    transfers(2)=vec(5)                         !ahumarch1522 adding cornersol
        !    vmax(1:2)=vec(3:4)+inc_coef*transfers(1:2)
        !    GO TO 2017
        !else if ( vec(5) <= vdif(2)-vdif(1) + eps  )   then 
        !    transfers(1)=vec(5)                         !ahumarch1522 adding cornersol
        !    transfers(2)=0.0_dp                         !ahumarch1522 adding cornersol
        !    vmax(1:2)=vec(3:4)+inc_coef*transfers(1:2)
        !    GO TO 2017
        !!!!!!!!else if (relmax==0) then 
        !!!!!!!    vmax(1:2)=vec(1:2)
        !!!!!!!!    GO TO 2017
        2017 dd(8)=jmax
        dd(9)=qmax
        dd(10)=relmax 
        if (yaz.and.callfrom==40) write(200,*) "in sol/getdec couples end: jmax,qmax,relmax,surplus),vec(1:5),transfers"
        if (yaz.and.callfrom==40) write(200,'(I4,I8,I4,8F9.2)') jmax,qmax,relmax,surplus,vec(1:5),transfers 
        if (yaz.and.callfrom==50) write(200,*) "in sol/getdec mar mkt end: jmax,qmax,relmax,surplus),vec(1:5),transfers"
        if (yaz.and.callfrom==50) write(200,'(I4,I8,I4,8F9.2)') jmax,qmax,relmax,surplus,vec(1:5),transfers 
        if (yaz.and.callfrom==80) write(400,*) "in sim/getdec couples end: jmax,qmax,relmax,surplus),vec(1:5),transfers"
        if (yaz.and.callfrom==80) write(400,'(I4,I8,I4,8F9.2)') jmax,qmax,relmax,surplus,vec(1:5),transfers  
    end subroutine getdecR    



    !ag090122 agsept2022 Changed getdec_c to getdec in order to include mar market decisinos in getdec as well
    !                     so that everything is all in one same place and if decisino protocol changes are made they have to be made to only one place
	!                     OLDER VERSIONS OF GETDEC_C CAN BE FOUND IN MIGSOLGETDEC.F90 FILE
    subroutine getdec(dd,vmax,vsing,transfers) 
        integer(i4b), intent(inout) :: dd(:)
        real(dp), intent(out) :: vmax(2),transfers(2) 
        real(dp), intent(in) :: vsing(2) 
        real(dp) :: vecj(5,nc),vec(5),surplusj(nc),surplus,vdif(2),mcost(2),asum,yazvec(5),vmaxcheck(2)
        integer(i4b) :: ia,index,q,x,z,q0,g,jmax,qmax,relmax,i,i0,n,trueindex,j,de(1),ed(2),locch,loc0,callfrom
        logical :: haveenough(nc),haveenoughtomakeindiff,haveenoughtomakeindiff_alt,intsol(3),pc(2),pc_alt(2),haveenoughforNB(nc),haveenoughNB
        integer(i4b) :: trueco,truetyp,truehome,iepjoint,iephub,iepwfe,indeces(2),tempmax
        callfrom=dd(7)
        ia=dd(1) 
        trueindex=dd(2) 
        if (groups) then 
            index=1
        else 
            index=trueindex
        end if
        call index2cotyphome(trueindex,trueco,truetyp,truehome)
        q=dd(3) 
        x=dd(4) ; ed(:)=xx2e(:,x)   
        z=dd(5)
        q0=dd(6)  
        vec=pen				! initialize
        jmax=-1				! initialize
        qmax=-1				! initialize
        relmax=-1            ! initialize
        vmax=pen			! initialize
        transfers=pen 
        if (Callfrom==40.or.callfrom==80) then !calling from sol/couples/getdec (40) OR simul/couples/getdec (80) 
            vec=pen
            loc0=qq2l(1,q0)     !ahumarch2022 ahu032022
            iepjoint=dd(11) !note that dd(5) used to be z and iepsmove used to be dd(11) before
            indeces=lin2ndim( (/ nepsmove , nepsmove /) , iepjoint )
            iephub=indeces(1) ; iepwfe=indeces(2)  
            i = qq2q(1,q) ; n=xx2x(1,x) ; i0 = qq2q(1,q0) 
            vec(1) = vm_postdiv(iephub,n,i,i0,ia,index) + diveduc(x2e(n)) !+ divpenalty * one( q2l(qchoicesingle) /= loc0 )
            i = qq2q(2,q) ; n=xx2x(2,x) ; i0 = qq2q(2,q0) 	
            vec(2) = vf_postdiv(iepwfe,n,i,i0,ia,index) + diveduc(x2e(n)) !+ divpenalty * one( q2l(qchoicesingle) /= loc0 )
            vecj=pen ; surplusj=pen ; haveenough=.FALSE. ; yazvec=pen
            choice: do j=1,nc	
                i = ch(j,q,q0)	!alternative q
                if (i>0 ) then		
                    locch=qq2l(1,i)     !ahumarch2122 ahu032122
                    if (callfrom==40) then !calling from sol
                    vecj(3,j) = vm0ctemp(i,x)       + utilmar(MALES,x,z,trueindex,ia)       + one( locch /= loc0 ) * (movecost(xx2x(1,x),qq2q(1,q0),trueindex,MALES)    + moveshockmar(iephub,1)    + distance(locch,loc0)*cstadjacent )  
                    vecj(4,j) = vf0ctemp(i,x)       + utilmar(FEMALES,x,z,trueindex,ia)     + one( locch /= loc0 ) * (movecost(xx2x(2,x),qq2q(2,q0),trueindex,FEMALES)  + moveshockmar(iepwfe,2)    + distance(locch,loc0)*cstadjacent )
                    else if (callfrom==80) then !calling from simulation
                    vecj(3,j) = vm0_c(x,i,ia,index)  + utilmar(MALES,x,z,trueindex,ia)      + one( locch /= loc0 ) * (movecost(xx2x(1,x),qq2q(1,q0),trueindex,MALES)    + moveshockmar(iephub,1)    + distance(locch,loc0)*cstadjacent )  
                    vecj(4,j) = vf0_c(x,i,ia,index)  + utilmar(FEMALES,x,z,trueindex,ia)    + one( locch /= loc0 ) * (movecost(xx2x(2,x),qq2q(2,q0),trueindex,FEMALES)  + moveshockmar(iepwfe,2)    + distance(locch,loc0)*cstadjacent )
                    end if 
                    !vec(5) = wc(1,x,i,trueindex) + wc(2,x,i,trueindex) + one( qq2w(1,q0)<=np )* ubc(1,x,i,trueindex) + one( qq2w(2,q0)<=np )* ubc(2,i,x,trueindex) + nonlabinc + nonlabinc 	 !ahu summer18 050318: added the ubc
                    vecj(5,j) = wmctemp(i,x,ia,trueindex) + wfctemp(i,x,ia,trueindex) + nonlabinc(ed(1)) + nonlabinc(ed(2))  + ubenefit_c(i)     
                    yazvec(1:2)=vec(1:2) ; yazvec(3:5)=vecj(3:5,j) !just for yaz
                    vecj(1,j)=vecj(3,j)-vec(1)   !ahumarch1522 ahu031522 adding cornersol
                    vecj(2,j)=vecj(4,j)-vec(2)   !ahumarch1522 ahu031522 adding cornersol
                    !ahumarch2122 ahu032122 replacing this with vdif for saving time surplusj(j) = vec(5) + vec(3) - vec(1) + vec(4) - vec(2)  
                    surplusj(j) = vecj(5,j) + sum(vecj(1:2,j))  !ahumarch2122 ahu032122 replacing with vdif for saving time    
                    pc(1:2)	= ( vecj(1:2,j) + eps >= 0.0_dp )
                    asum = sum(  one(.not. pc)  *   abs( vecj(1:2,j) )   ) 
                    haveenough(j)=(  vecj(5,j) + eps - asum  >= 0.0_dp  ) ; haveenoughforNB(j)= (  vecj(5,j) + eps  >=  abs(vecj(1,j) -  vecj(2,j)) ) 
                    if (yaz) then 
                        dd(8)=j ; dd(9)=i !to tell yaz_getdec which alternative is being evaluated (where j is the choice and i is the corresponding q)
                        call yaz_getdec(dd,yazvec(1:5),surplusj(j),pc(1:2),asum,haveenough(j))
                    end if
                end if      
            end do choice    
            IF (.NOT.UPFRONTTRANSFERS) THEN
                de=maxloc(surplusj,MASK=(haveenough.AND.haveenoughforNB)) 
            ELSE 
                de=maxloc(surplusj,MASK=(surplusj>=0.0_dp)) 
            END IF 
            if (de(1)>0) then 
                jmax=de(1) ; qmax=ch(jmax,q,q0) ; relmax=1
                surplus=surplusj(jmax)                
                vdif(1:2)=vecj(1:2,jmax)    
                vec(3:5)=vecj(3:5,jmax)
                transfers(1) = alf * surplus -  vdif(1)                                                           
                transfers(2) = (1.0_dp-alf) * (vec(5) + vdif(1) ) - alf*vdif(2)  
                if (    abs(transfers(1)+transfers(2) - vec(5) ) >  0.00001_dp  ) then 
                    print*, "Stop it right now",transfers(1:2),vec(5)
                end if 
                vmaxcheck(1:2)=vec(3:4)+transfers(1:2)
                vmax(1:2)=vec(1:2)+alf*surplus
                if (    abs(vmaxcheck(1)-vmax(1))>0.00001_dp   .or.  abs(vmaxcheck(2)-vmax(2))>0.00001_dp ) then 
                    print*, "Stop it right now"
                end if 
                !if (skriv) then 
                    !myefficiency(:)%ve1=vecj(1,:)    
                    !myefficiency(:)%ve2=vecj(2,:)    
                    !myefficiency(:)%ve3=vecj(3,:)    
                    !myefficiency(:)%ve4=vecj(4,:)    
                    !myefficiency(:)%ve5=vecj(5,:)    
                    !myefficiency(:)%sur=surplusj(:) 
                    !myefficiency(:)%haveenough=haveenough(:)
                    !myefficiency(:)%haveenoughforNB=haveenoughforNB(:)
                    de=maxloc(surplusj,MASK=(haveenough.AND.haveenoughforNB))
                    tempmax=de(1)
                    if (de(1)>0) then 
                        myefficiency%maxA=ch(tempmax,q,q0) 
                    else 
                        myefficiency%maxA=0
                    end if
                    de=maxloc(surplusj,MASK=(surplusj>=0.0_dp)) 
                    tempmax=de(1)
                    if (de(1)>0) then 
                        myefficiency%maxB=ch(tempmax,q,q0) 
                    else 
                        myefficiency%maxB=0
                    end if
                    myefficiency%tra1=transfers(1) 
                    myefficiency%tra2=transfers(2) 
                !end if
                !SAVETIME IF (.NOT.UPFRONTTRANSFERS) THEN
                !SAVETIME     if (minval(transfers)<-0.00001 .or. abs(vmaxcheck(1)-vmax(1))>0.00001_dp .or.  abs(vmaxcheck(2)-vmax(2))>0.00001_dp ) then 
                !SAVETIME         print*, "Seriously though", transfers,vmaxcheck,vmax 
                !SAVETIME         stop
                !SAVETIME     end if 
                !SAVETIME     if ( vec(5) + eps < abs(vdif(1)-vdif(2)) )  then 
                !SAVETIME         print*, "why here should not be here like at all ", callfrom,vec(5),vdif(1)-vdif(2)
                !SAVETIME         stop
                !SAVETIME     end if
                !SAVETIME END IF !NO UPFRONTTRANSFERS.
                !SAVETIME if (yaz) then ; write(200,*)  ; write(200,*) ; write(400,*) ; write(400,*) ; end if
                !SAVETIME if (yaz.and.callfrom==40) write(200,*) "in sol/getdec couples de(1)>0: jmax,qmax,relmax,surplusJ(jmax),vecj(1:5,jmax) NOW GO TO 2013"
                !SAVETIME if (yaz.and.callfrom==40) write(200,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplusj(jmax),vecj(1:5,jmax)
                !SAVETIME if (yaz.and.callfrom==80) write(400,*) "in sim/getdec couples de(1)>0: jmax,qmax,relmax,surplusJ(jmax),vecj(1:5,jmax) NOW GO TO 2013"
                !SAVETIME if (yaz.and.callfrom==80) write(400,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplusj(jmax),vecj(1:5,jmax)
                GO TO 2017
            else 
                jmax=-1 ; qmax=-1  ; relmax=0 ; vmax(1:2)=vec(1:2) ; transfers=pen
                !SAVETIME if (yaz) then ; write(200,*)  ; write(200,*) ; write(400,*) ; write(400,*) ; end if
                !SAVETIME if (yaz.and.callfrom==40) write(200,*) "in sol/getdec couples de(1)<=0: jmax,qmax,relmax,de(1) NOW GO TO 2017"
                !SAVETIME if (yaz.and.callfrom==40) write(200,*) jmax,qmax,relmax,de(1)
                !SAVETIME if (yaz.and.callfrom==80) write(400,*) "in sim/getdec couples de(1)<=0: jmax,qmax,relmax,de(1) NOW GO TO 2017"
                !SAVETIME if (yaz.and.callfrom==80) write(400,*) jmax,qmax,relmax,de(1)
                GO TO 2017
            end if 
        end if 

        if (callfrom==50) then !calling from sol/marmkt (this getdec instance is only called in sol (not sim because I save all decisions to dec array)
            if (dd(11).ne.-1) then ; print*, "iepsingle is not -1 in call to getdec and callfrom is from sol/marmkt!!" ; stop ; end if 
            vec(1:2)=vsing(1:2) !ag090822 agsept2022 vsing(1:2) 
            vec(3) = vm0ctemp(q,x) + utilmar(MALES,x,z,trueindex,ia) 
            vec(4) = vf0ctemp(q,x) + utilmar(FEMALES,x,z,trueindex,ia) 
            vec(5) = wmctemp(q,x,ia,trueindex) +  wfctemp(q,x,ia,trueindex) + nonlabinc(ed(1)) + nonlabinc(ed(2))    + ubenefit_c(q)                                                
            vdif(1)=vec(3)-vec(1)   !ahumarch1522 ahu031522 adding cornersol
            vdif(2)=vec(4)-vec(2)   !ahumarch1522 ahu031522 adding cornersol
            !ahumarch2122 ahu032122 replacing this with vdif for saving time surplus = vec(5) + vec(3) - vec(1) + vec(4) - vec(2)  
            surplus = vec(5) + sum(vdif(1:2))  !ahumarch2122 ahu032122 replacing with vdif for saving time
            pc(1:2)	= ( vdif + eps >= 0.0_dp )	!pc(1:2)    = ( vec(3:4) - vec(1:2) >= 0.0_dp )						        
            !ahu032122 ahumarch2122 commenting out to save time pc_alt(1:2)=( vdif >= 0.0_dp )	
            asum = sum(  one(.not. pc)  *   abs( vdif )   ) 
            haveenoughtomakeindiff=(  vec(5) + eps - asum  >= 0.0_dp  ) ; haveenoughNB= (vec(5) +eps >= abs(vdif(1)-vdif(2)) )
            IF (.NOT.UPFRONTTRANSFERS) THEN
                de(1)=one( surplus+eps > 0.0_dp .and. haveenoughtomakeindiff .and. haveenoughNB)             
            ELSE 
                de(1)=one( surplus+eps > 0.0_dp )             
            END IF 
            !SAVETIME if (yaz) then 
                !callfrom=50 is sol/marmkt/getdec and before calling getdec, dd(9) is just set to q i.e. equal to dd(3) and dd dd(8)=-1 ; dd(9)=-1 !to tell yaz_getdec which alternative is being evaluated (where j is the choice and i is the corresponding q)
                !but when calling from singles mar market there is no j or i so set equal to -1 
                !SAVETIME call yaz_getdec(dd,vec(1:5),surplus,pc(1:2),asum,haveenoughtomakeindiff) !called from sol/getdec/marmkt 
            !SAVETIME end if 
            if (de(1)>0) then 
                jmax=-1 ; qmax=-1 ; relmax=1
                transfers(1) = alf * surplus -  vdif(1)                                                           
                transfers(2) = (1.0_dp-alf) * (vec(5) + vdif(1) ) - alf*vdif(2)  
                vmaxcheck(1:2)=vec(3:4)+transfers(1:2)
                vmax(1:2)=vec(1:2)+alf*surplus
                if (    abs(vmaxcheck(1)-vmax(1))>0.00001_dp   .or.  abs(vmaxcheck(2)-vmax(2))>0.00001_dp ) then 
                    print*, "Stop it right now"
                end if 
                if (    abs(transfers(1)+transfers(2) - vec(5) ) >  0.00001_dp  ) then 
                    print*, "Stop it right now this is sol mar market ",transfers(1:2),vec(5)
                end if 
                !SAVETIME IF (.NOT.UPFRONTTRANSFERS) THEN
                !SAVETIME     if (minval(transfers)<-0.00001 .or. abs(vmaxcheck(1)-vmax(1))>0.00001_dp .or.  abs(vmaxcheck(2)-vmax(2))>0.00001_dp ) then 
                !SAVETIME         print*, "Seriously though", transfers,vmaxcheck,vmax 
                !SAVETIME         stop
                !SAVETIME     end if 
                !SAVETIME     if ( vec(5) + eps < abs(vdif(1)-vdif(2)) )  then 
                !SAVETIME         print*, "why here should not be here like at all ", callfrom
                !SAVETIME         stop
                !SAVETIME     end if
                !SAVETIME END IF !UPFRONT transfers   
                !SAVETIME if (yaz) write(200,*) "in sol/getdec marmkt de(1)>0: jmax,qmax,relmax,surplus,vec(1:5) NOW GO TO 2013"
                !SAVETIME if (yaz) write(200,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplus,vec(1:5) 
                GO TO 2017 !have to calculate vmax and transfers before going to 2017 
            else 
                jmax=-1 ; qmax=-1  ; relmax=0 ; vmax(1:2)=vec(1:2) ; transfers=pen 
                !SAVETIME if (yaz) write(200,*) "in sol/getdec marmkt de(1)<=0: jmax,qmax,relmax,surplus,vec(1:5) NOW GO TO 2017"
                !SAVETIME if (yaz) write(200,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplus,vec(1:5)                 
                GO TO 2017 ! go directly to 2017
            end if 
        end if 

        !CORNER SOLUTIONS NOW OBSOLETE
        !else if ( vec(5) <= vdif(1)-vdif(2) + eps  )  then 
        !    transfers(1)=0.0_dp                         !ahumarch1522 adding cornersol
        !    transfers(2)=vec(5)                         !ahumarch1522 adding cornersol
        !    vmax(1:2)=vec(3:4)+inc_coef*transfers(1:2)
        !    GO TO 2017
        !else if ( vec(5) <= vdif(2)-vdif(1) + eps  )   then 
        !    transfers(1)=vec(5)                         !ahumarch1522 adding cornersol
        !    transfers(2)=0.0_dp                         !ahumarch1522 adding cornersol
        !    vmax(1:2)=vec(3:4)+inc_coef*transfers(1:2)
        !    GO TO 2017
        !!!!!!!!else if (relmax==0) then 
        !!!!!!!    vmax(1:2)=vec(1:2)
        !!!!!!!!    GO TO 2017
        2017 dd(8)=jmax
        dd(9)=qmax
        dd(10)=relmax 
        if (yaz.and.callfrom==40) write(200,*) "in sol/getdec couples end: jmax,qmax,relmax,surplus),vec(1:5),transfers"
        if (yaz.and.callfrom==40) write(200,'(I4,I8,I4,8F9.2)') jmax,qmax,relmax,surplus,vec(1:5),transfers 
        if (yaz.and.callfrom==50) write(200,*) "in sol/getdec mar mkt end: jmax,qmax,relmax,surplus),vec(1:5),transfers"
        if (yaz.and.callfrom==50) write(200,'(I4,I8,I4,8F9.2)') jmax,qmax,relmax,surplus,vec(1:5),transfers 
        if (yaz.and.callfrom==80) write(400,*) "in sim/getdec couples end: jmax,qmax,relmax,surplus),vec(1:5),transfers"
        if (yaz.and.callfrom==80) write(400,'(I4,I8,I4,8F9.2)') jmax,qmax,relmax,surplus,vec(1:5),transfers  
    end subroutine getdec

    subroutine getdec_s(dd,vs,qs,vs1,whenmakingdec)
	integer(i4b), intent(in) :: dd(:), whenmakingdec
	real(dp), dimension(nxs,nqs), intent(in) :: vs
	integer(i4b), dimension(nepsmove,nxs,nqs,nqs), intent(out) :: qs
	real(dp), dimension(nepsmove,nxs,nqs,nqs), intent(out) :: vs1
	integer(i4b) :: q0,q,l,l0,w,x,i,j,iepsingle,sex,index,trueindex,jstar(1),age,ss(size(dd,1)),ia
    real(dp) :: mcost(nxs,nqs),moveshock(nepsmove),vcho(ncs),addon
    integer(i4b) :: trueco,truetyp,truehome

    if (dd(5).ne.-1) then ; print*, "z is not -1 in call to getdec_s!" ; stop ; end if 
    if (dd(11).ne.-1) then ; print*, "iepsingle is not -1 in call to getdec_s!" ; stop ; end if 
    trueindex=dd(2)
    sex=dd(7) ; ia=dd(1)
	qs = -1 ; vs1 = 0.0_dp
    if (groups) then 
        index=1
    else 
        index=trueindex
    end if
    call index2cotyphome(trueindex,trueco,truetyp,truehome)
    if (whenmakingdec==1 ) then !1 indicates teh decision is being made AFTER DIVORCE
        mcost(:,:)=movecost(:,:,trueindex,SEX) !+divpenmove !+ alfdivpen*ia !move cost as a married person right after divorce
        moveshock(:)=moveshockmar(:,sex)
    else if (whenmakingdec==2 ) then !2 indicates teh decision is being made BEFORE MARRIAGE MARKET
        mcost(:,:)=movecost(:,:,trueindex,SEX) !move cost as a single person
        moveshock(:)=moveshocksin(:,sex)
    else 
        print*, 'decs something wrong'
        stop
    end if 
	do q0=1,nqs
        l0=q2l(q0)
		if ( q2w(q0)<=np1 ) then 
            ofloc: do l=1,nl											
            wagedraw: do w=1,np2
            q = wl2q(w,l)
            do x=1,nxs
                !mcost=movecost_s(x,q0,trueindex)    !movecost of a single person
                !if (whenmakingdec==1 ) then !1 indicates teh decision is being made AFTER DIVORCE
                !    addon= divpenbyed(x,sex) 
                !else if (whenmakingdec==2 ) then !2 indicates teh decision is being made BEFORE MARRIAGE MARKET
                !    addon=0.
                !end if 
            
                do iepsingle=1,nepsmove !moveshocks

                    vcho=pen
                    choice: do j=1,ncs	
                        i = chs(j,q,q0)	!alternative q
                        if (i>0 ) then		
                            vcho(j) = vs(x,i) + one(q2l(i)/=l0) * (mcost(x,q0)  + moveshock(iepsingle)  + distance(q2l(i),l0)*cstadjacent) !+ addon
                        end if 
                    end do choice
                    
                                                
                    vs1(iepsingle,x,q,q0) =	maxval(vcho)
                    jstar(1)=maxloc(vcho,1)
                    qs(iepsingle,x,q,q0)  =  chs(jstar(1),q,q0)
                    if ( vs1(iepsingle,x,q,q0) < pen+1.0_dp) then 
                        print*, 'There is something wrong in decs new new new',vs1(iepsingle,x,q,q0),dd(1),moveshock(iepsingle),mcost(x,q0)
                        print*, 'vcho',vcho
                        stop
                    end if 
                    
                    if (yaz .and.q2w(q0)==1.and.l==8.and.x==1) then
                        ss = dd  !(/ ia,index,q,x,iepjoint (used to be z),q0,gender,j,altq,iepsingle /)  													
                        ss(3:6)= (/q,x,-1,q0/)
                        ss(11)=iepsingle
                        !!!!!!call yaz_decs(ss,vcho)   !look here
                    end if 

            end do !moveshocks		
            end do !x
            end do wagedraw
            end do ofloc             
		    end if		! w0<=np1 i.e. ps0(q0) is true
	end do			! q0
	end subroutine getdec_s

    
	subroutine get_util_w
	integer(i4b) :: q,x,w(2),l(2),kid(2),ed(2),expe(2),trueindex,trueco,truetyp,truehome,g,j,k
    real(dp) :: epsw(2) !RHO(2,2),CD(2,2),
    real(sp) :: wsgross,wcgross(2) !pwages(1:numbin),swages(1:numbin)
    integer(i4b) :: bracket,bracketprev,bracketnext,sbrack,pbrack,z,loc,age
    real(dp) :: statetaxtot,fedtaxtot,statetaxtotc(2),fedtaxtotc(2)

    utils=pen
    utilc=pen
    ws=pen
    wc=pen   
    
    !if (mysay==1) then 
    !    open(unit=68857,file='checktax1.txt')
    !    open(unit=68858,file='checktax2.txt')
    !end if 
    do age=MNA,MXA
    do g=MALES,FEMALES
    do trueindex=1,ninp
        do z=1,nz 
            do x=1,nx
                utilmar(g,x,z,trueindex,age)= mg(trueindex)*one(xx2kid(1,x)==2)  + marshock(z) !+ umared*one(xx2e(1,x)==2) !+ umarkid*one(xx2kid(1,x)==2)  + umared*one(maxval(xx2e(:,x))==2)  !+ umaria * age
            !    if (xx2kid(1,x).ne.xx2kid(2,x)) then ; print*, 'sth wrong with kids' ; stop ; end if
            end do 
        end do 
    end do 
    end do 
    end do 

    ubenefit_s=0.0_dp
    do q=1,nqs
        if (q2w(q) <= np ) then
            ubenefit_s(q) = 0.0_dp 
        else if (q2w(q) == np1 ) then
            ubenefit_s(q)= ubenefitbyloc(q2l(q))
        else 
            ubenefit_s(q)= 0.0_dp
        end if 
    end do 
    ubenefit_c=0.0_dp
    do q=1,nq
        w(:)=qq2w(:,q) ; loc=qq2l(1,q)
        if ( w(1) <= np .and. w(2)<=np  ) then
            ubenefit_c(q) = 0.0_dp 
        else if ( w(1) == np1 .and. w(2)<=np  ) then
            ubenefit_c(q)= ubenefitbyloc(loc)
        else if ( w(1) <= np .and. w(2) == np1  ) then
            ubenefit_c(q)= ubenefitbyloc(loc) 
        else if ( w(1) == np1 .and. w(2) == np1  ) then
            ubenefit_c(q)= ubenefitbyloc(loc) + ubenefitbyloc(loc) 
        else 
            ubenefit_c(q)= 0.0_dp
        end if 
    end do 

    do trueindex=1,ninp    
        do age=MNA,MXA
        call index2cotyphome(trueindex,trueco,truetyp,truehome)
		qs: do q=1,nqs 
            xs: do x=1,nxs    
            !call x2edexpkid(x,ed,exp,kid)    
            do g=1,2
                movecost(x,q,trueindex,g)=fnmove( q2w(q),x2kid(x),x2e(x),trueindex,g) !0 to indicate singles 
                w(g) = q2w(q)						! wage 
			    l(g) = q2l(q)						! location
                !pwages(1:numbin)=tax(1:numbin,numbin,l(g))%pwages
			    if ( w(g) <= np ) then	
                    epsw(g)=wg(w(g),x2e(x),g) !sig_wge(g)*wg(w(g),g) !FEB7 2023: BIG SIGLOC(1:NL) CHANGE TO WAGE PROCESS AND LATER SIGWGE CHANGED TO SIGWGE BY ED. NO MORE LOC. 
				    ws(g,x,q,age,trueindex)	= fnwge(g,truetyp, l(g),epsw(g), x2e(x), x2r(x),age) 
                    wsgross=ws(g,x,q,age,trueindex)
                    bracket=locate(  pbracket(1:numbin)  ,  wsgross ) 
                    if (bracket<1.or.bracket>numbin) then               
                        print*, "There is something wrong with locate bracket",bracket
                        print*, "pwages(1:numbin)",pbracket(1:numbin)
                        print*, "wsgross         ",wsgross
                        stop
                    end if 
                    !if (mysay==1) then 
                        !if (bracket<numbin) then ; bracketnext=bracket+1 ; else ; bracketnext=numbin ; end if 
                        !if (bracket>1) then ; bracketprev=bracket-1 ; else ; bracketprev=1 ; end if 
                        !write(68857,'(tr3,"wsgross",tr1,"bracket",tr4,"prev",tr4,"here",tr4,"next")')
                        !write(68857,'(f10.1,i4,3f10.1)') wsgross,bracket,pwages(bracketprev),pwages(bracket),pwages(bracketnext)
                        !write(68857,*) 
                        !write(68858,'(tr2,"truind",tr7,"q",tr7,"x",tr1,"sex",tr1,"typ",tr3,"l",tr2,"ed",tr1,"exp",tr3,"wsgross",tr5,"wsnet",tr2,"statetax",tr4,"fedtax")')
                        !write(68858,'(3i8,5i4,4f10.2)') trueindex,q,x,g,truetyp,l(g),x2e(x),x2r(x),wsgross,wsnet(g,x,q,trueindex),staterate,fedrate
                        !write(68858,*) 
                    !end if 
                else if ( w(g) == np1 ) then 
                    ws(g,x,q,age,trueindex)	 = 0.0_dp
                    wsgross=ws(g,x,q,age,trueindex)
                    !wsnet(g,x,q,trueindex)=0.0_dp
                    bracket=1
                end if 
                
                call get_taxliability_s(l(g),bracket,wsgross,statetaxtot,fedtaxtot)
                wsnet(g,x,q,age,trueindex)=wsgross-(statetaxtot+fedtaxtot) 
                !if (mysay==0) write(*,'(3F12.1)') wsgross,statetaxtot,fedtaxtot

                !if (taxset==0) then
                !    wsnet(g,x,q,trueindex)=wsgross !ws(g,x,q,trueindex)
                !end if 

                
                if ( w(g) <= np ) then						
				    utils(g,x,q,trueindex)	= uhomet(g) * one(l(g)==truehome) + uloc(l(g))  + nonlabinc(x2e(x)) !+ ulocheat*heat(l(g))
			    else if ( w(g) == np1 ) then 
				    !utils(g,x,q,trueindex)	= uhomet(g) * one(l(g)==truehome) + uloc(l(g))  + (alphaed(g,x2e(x))  + alphakid(g)) * one(x2kid(x)>1)  + nonlabinc(x2e(x)) + ubenefit_s(q)!+ ulocheat*heat(l(g))
                    utils(g,x,q,trueindex)	= uhomet(g) * one(l(g)==truehome) + uloc(l(g))  + alphakid(g) * one(x2kid(x)>1)  + nonlabinc(x2e(x)) + ubenefit_s(q)  + alphab(g) !+ ulocheat*heat(l(g))
                end if
                
                !ahu october2022: 
                !alphakid(g): kid(g)=1 is no kid, and kid(g)=2 and above is yes kid.   this alphakid used to have two dimensions for no reason. 
                !ahu later in january 2023: what is this kid(g) thing? and why is alphakid declared as alphakid(nkid) and yet it's really just alphakid(sex) i.e. alphakid(g)
                !it doesnt' affect anythign I suppose since in all the code it is always alphakid(g) even though it's declared as alphakid(nkid)
                !since nkid is 2 and number of sex is 2 as well, this doesn't really affect anythign but still it is not right. 
                !nkid is never above 2 though (i set it that way in data too for numkids above 2 is just always 2 in data). 
                !alphaed(g,ed(g)):   ed(g)=1 is noed and ed(g)=2 is yes ed (where  neduc is 2)    
            end do !gender
		end do xs
	end do qs
    end do !age
    !close(68857)
    !close(68858)

    do age=MNA,MXA
    qc: do q=1,nq
	xc: do x=1,nx
        ed(:)=xx2e(:,x)    
        expe(:)=xx2r(:,x)    
        kid(:)=xx2kid(:,x)           
        w(:) = qq2w(:,q)						! wage 
        l(:) = qq2l(:,q)						! location
        if (l(1).ne.l(2)) then ; print*, 'lm not equal to lf' ; stop ; end if 

        !******************************
        !ahu summer18 050318 
        !if (w(1)<=np) then 
        !    ubc(1,q,x,trueindex)	= 0.0_dp           
        !else if (w(1)==np1) then 
        !    epsw(1)=0.0_dp
        !    ubc(1,q,x,trueindex)	= replacement_rate*fnwge(1,truetyp, l(1),epsw(1), ed(1), expe(1) ) 
        !end if 
        !if (w(2)<=np) then 
        !    ubc(2,q,x,trueindex)	= 0.0_dp           
        !else if (w(2)==np1) then 
        !    epsw(2)=0.0_dp
        !    ubc(2,q,x,trueindex)	= replacement_rate*fnwge(2,truetyp, l(2),epsw(2), ed(2), expe(2) ) 
        !end if 
        !ahu summer18 050318 
        !******************************

        !pwages(1:numbin)=tax(1:numbin,numbin,l(1))%pwages
        !swages(1:numbin)=tax(numbin,1:numbin,l(2))%swages
        if ( w(1) <= np .and. w(2) <= np ) then		
            epsw(1)=wg(w(1),ed(1),1) !CD(1,1)*wg(w(1),1) !FEB7 2023: BIG SIGLOC(1:NL) CHANGE TO WAGE PROCESS
            epsw(2)=wg(w(2),ed(2),2) !CD(2,1)*wg(w(1),1) + CD(2,2)*wg(w(2),2) !FEB7 2023: BIG SIGLOC(1:NL) CHANGE TO WAGE PROCESS
            wc(1,x,q,age,trueindex)	= fnwge(1,truetyp, l(1),epsw(1), ed(1), expe(1) ,age ) 
            wc(2,x,q,age,trueindex)	= fnwge(2,truetyp, l(2),epsw(2), ed(2), expe(2) ,age )  
            wcgross(1:2)= wc(1:2,x,q,age,trueindex)
            pbrack=locate(  pbracket(1:numbin)  ,  wcgross(1) ) 
            sbrack=locate(  sbracket(1:numbin)  ,  wcgross(2) ) 
        else if ( w(1) <= np .and. w(2) == np1 ) then		
            epsw(1)=wg(w(1),ed(1),1) !sig_wge(1)*wg(w(1),1) !FEB7 2023: BIG SIGLOC(1:NL) CHANGE TO WAGE PROCESS
            wc(1,x,q,age,trueindex)	= fnwge(1,truetyp, l(1),epsw(1), ed(1), expe(1) ,age) 
            wc(2,x,q,age,trueindex)	= 0.0_dp
            wcgross(1:2)= wc(1:2,x,q,age,trueindex)
            pbrack=locate(  pbracket(1:numbin)  ,  wcgross(1) ) 
            sbrack=1
        else if ( w(1) == np1 .and. w(2) <= np ) then		
            epsw(2)=wg(w(2),ed(2),2) !sig_wge(2)*wg(w(2),2) !FEB7 2023: BIG SIGLOC(1:NL) CHANGE TO WAGE PROCESS
            wc(1,x,q,age,trueindex)	= 0.0_dp
            wc(2,x,q,age,trueindex)	= fnwge(2,truetyp, l(2),epsw(2), ed(2), expe(2) ,age)   
            wcgross(1:2)= wc(1:2,x,q,age,trueindex)
            pbrack=1
            sbrack=locate(  sbracket(1:numbin)  ,  wcgross(2) ) 
        else if ( w(1) == np1 .and. w(2) == np1 ) then		
            wc(1,x,q,age,trueindex)	= 0.0_dp
            wc(2,x,q,age,trueindex)	= 0.0_dp           
            wcgross(1:2)= 0.0_dp
            pbrack=1
            sbrack=1
        end if 

    
        call get_taxliability_c(l(1),pbrack,sbrack,wcgross(1:2),statetaxtotc,fedtaxtotc)
        wcnet(1:2,x,q,age,trueindex)=wcgross(1:2)- ( sum(statetaxtotc(1:2)) + sum(fedtaxtotc(1:2)) ) 

        !if (taxset==0) then
        !    wcnet(1:2,x,q,trueindex)	= wc(1:2,x,q,trueindex)
        !end if 

            
        do g=1,2
            if ( w(g) <= np ) then						
                utilc(g,x,q,trueindex)	= uhomet(g) * one(l(g)==truehome)  + uloc(l(g)) !+ ulocheat*heat(l(g))
            else if ( w(g) == np1 ) then 
                !utilc(g,x,q,trueindex)	= uhomet(g) * one(l(g)==truehome)  + uloc(l(g)) + (alphaed(g,ed(g) ) + alphakid(g) ) * one(kid(g)>1)  !+ ulocheat*heat(l(g))
                utilc(g,x,q,trueindex)	= uhomet(g) * one(l(g)==truehome)  + uloc(l(g)) + alphakid(g) * one(kid(g)>1) + alphab(g) !+ ulocheat*heat(l(g))
                !ahu october2022: 
                !alphakid(g): kid(g)=1 is no kid, and kid(g)=2 and above is yes kid.   this alphakid used to have two dimensions for no reason. 
                !nkid is never above 2 though (i set it that way in data too for numkids above 2 is just always 2 in data). 
                !alphaed(g,ed(g)):   ed(g)=1 is noed and ed(g)=2 is yes ed (where  neduc is 2)
            end if
        end do   
        
    end do xc
	end do qc
    end do !age

    end do !trueindex

    !Male and female wage shocks for each period and from each location
!    do it=1,nt
!	    do iloc=1,nloc
!		    do rn=1,nn
!			    up1(it,rn,iloc)=CD(1,1)*rv1(it,rn,iloc)
!			    up2(it,rn,iloc)=CD(2,1)*rv1(it,rn,iloc)+CD(2,2)*rv2(it,rn,iloc)
!		    end do 
!	    end do
!    end do 

	end subroutine get_util_w
    


	subroutine get_util_wR
        integer(i4b) :: q,x,w(2),l(2),kid(2),ed(2),expe(2),trueindex,trueco,truetyp,truehome,g,j,k
        real(dp) :: epsw(2) !RHO(2,2),CD(2,2),
        real(sp) :: wsgross,wcgross(2) !pwages(1:numbin),swages(1:numbin)
        integer(i4b) :: bracket,bracketprev,bracketnext,sbrack,pbrack,z,loc,age
        real(dp) :: statetaxtot,fedtaxtot,statetaxtotc(2),fedtaxtotc(2)
    
        utils=pen
        utilc=pen
        ws=pen
        wc=pen   
        
        !if (mysay==1) then 
        !    open(unit=68857,file='checktax1.txt')
        !    open(unit=68858,file='checktax2.txt')
        !end if 
        do age=MNA,MXA
        do g=MALES,FEMALES
        do trueindex=1,ninp
            do z=1,nz 
                do x=1,nx
                    utilmar(g,x,z,trueindex,age)= mg(trueindex)*one(xx2kid(1,x)==2)  + marshock(z) !+ umared*one(xx2e(1,x)==2) !+ umarkid*one(xx2kid(1,x)==2)  + umared*one(maxval(xx2e(:,x))==2)  !+ umaria * age
                !    if (xx2kid(1,x).ne.xx2kid(2,x)) then ; print*, 'sth wrong with kids' ; stop ; end if
                end do 
            end do 
        end do 
        end do 
        end do 
    
        ubenefit_s=0.0_dp
        do q=1,nqs
            if (q2w(q) <= np ) then
                ubenefit_s(q) = 0.0_dp 
            else if (q2w(q) == np1 ) then
                ubenefit_s(q)= ubenefitbyloc(q2l(q))
            else 
                ubenefit_s(q)= 0.0_dp
            end if 
        end do 
        ubenefit_c=0.0_dp
        do q=1,nq
            w(:)=qq2w(:,q) ; loc=qq2l(1,q)
            if ( w(1) <= np .and. w(2)<=np  ) then
                ubenefit_c(q) = 0.0_dp 
            else if ( w(1) == np1 .and. w(2)<=np  ) then
                ubenefit_c(q)= ubenefitbyloc(loc)
            else if ( w(1) <= np .and. w(2) == np1  ) then
                ubenefit_c(q)= ubenefitbyloc(loc) 
            else if ( w(1) == np1 .and. w(2) == np1  ) then
                ubenefit_c(q)= ubenefitbyloc(loc) + ubenefitbyloc(loc) 
            else 
                ubenefit_c(q)= 0.0_dp
            end if 
        end do 
    
        do trueindex=1,ninp    
            do age=MNA,MXA
            call index2cotyphome(trueindex,trueco,truetyp,truehome)
            qs: do q=1,nqs 
                xs: do x=1,nxs    
                !call x2edexpkid(x,ed,exp,kid)    
                do g=1,2
                    movecost(x,q,trueindex,g)=fnmove( q2w(q),x2kid(x),x2e(x),trueindex,g) !0 to indicate singles 
                    w(g) = q2w(q)						! wage 
                    l(g) = q2l(q)						! location
                    !pwages(1:numbin)=tax(1:numbin,numbin,l(g))%pwages
                    if ( w(g) <= np ) then	
                        epsw(g)=wg(w(g),x2e(x),g) !sig_wge(g)*wg(w(g),g) !FEB7 2023: BIG SIGLOC(1:NL) CHANGE TO WAGE PROCESS AND LATER SIGWGE CHANGED TO SIGWGE BY ED. NO MORE LOC. 
                        ws(g,x,q,age,trueindex)	= fnwge(g,truetyp, l(g),epsw(g), x2e(x), x2r(x),age) 
                        wsgross=ws(g,x,q,age,trueindex)
                        bracket=locate(  pbracket(1:numbin)  ,  wsgross ) 
                        if (bracket<1.or.bracket>numbin) then               
                            print*, "There is something wrong with locate bracket",bracket
                            print*, "pwages(1:numbin)",pbracket(1:numbin)
                            print*, "wsgross         ",wsgross
                            stop
                        end if 
                        !if (mysay==1) then 
                            !if (bracket<numbin) then ; bracketnext=bracket+1 ; else ; bracketnext=numbin ; end if 
                            !if (bracket>1) then ; bracketprev=bracket-1 ; else ; bracketprev=1 ; end if 
                            !write(68857,'(tr3,"wsgross",tr1,"bracket",tr4,"prev",tr4,"here",tr4,"next")')
                            !write(68857,'(f10.1,i4,3f10.1)') wsgross,bracket,pwages(bracketprev),pwages(bracket),pwages(bracketnext)
                            !write(68857,*) 
                            !write(68858,'(tr2,"truind",tr7,"q",tr7,"x",tr1,"sex",tr1,"typ",tr3,"l",tr2,"ed",tr1,"exp",tr3,"wsgross",tr5,"wsnet",tr2,"statetax",tr4,"fedtax")')
                            !write(68858,'(3i8,5i4,4f10.2)') trueindex,q,x,g,truetyp,l(g),x2e(x),x2r(x),wsgross,wsnet(g,x,q,trueindex),staterate,fedrate
                            !write(68858,*) 
                        !end if 
                    else if ( w(g) == np1 ) then 
                        ws(g,x,q,age,trueindex)	 = 0.0_dp
                        wsgross=ws(g,x,q,age,trueindex)
                        !wsnet(g,x,q,trueindex)=0.0_dp
                        bracket=1
                    end if 
                    
                    call get_taxliability_s(l(g),bracket,wsgross,statetaxtot,fedtaxtot)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     wsnet(g,x,q,age,trueindex)=wsgross-(statetaxtot+fedtaxtot) 
                    !if (mysay==0) write(*,'(3F12.1)') wsgross,statetaxtot,fedtaxtot
    
                    !if (taxset==0) then
                    !    wsnet(g,x,q,trueindex)=wsgross !ws(g,x,q,trueindex)
                    !end if 
    
                    
                    if ( w(g) <= np ) then						
                        utils(g,x,q,age,trueindex)	= uhomet(g) * one(l(g)==truehome) + uloc(l(g)) + lama(g)*ln( ws(g,x,q,age,trueindex) + nonlabinc(x2e(x)) ) + (1.0_dp-lama(g)) * ln(leishours(FULLTIMEWORK))  !+ nonlabinc(x2e(x)) !+ ulocheat*heat(l(g))
                    else if ( w(g) == np1 ) then 
                        !utils(g,x,q,trueindex)	= uhomet(g) * one(l(g)==truehome) + uloc(l(g))  + (alphaed(g,x2e(x))  + alphakid(g)) * one(x2kid(x)>1)  + nonlabinc(x2e(x)) + ubenefit_s(q)!+ ulocheat*heat(l(g))
                        utils(g,x,q,age,trueindex)	= uhomet(g) * one(l(g)==truehome) + uloc(l(g))  + alphakid(g) * one(x2kid(x)>1)  + ubenefit_s(q)  + alphab(g) & !+ nonlabinc(x2e(x))  !+ ulocheat*heat(l(g))
                        & + lama(g)*ln( nonlabinc(x2e(x)) ) + (1.0_dp-lama(g)) * ln(leishours(NOWORK)) 
                    end if
                    
                    !ahu october2022: 
                    !alphakid(g): kid(g)=1 is no kid, and kid(g)=2 and above is yes kid.   this alphakid used to have two dimensions for no reason. 
                    !ahu later in january 2023: what is this kid(g) thing? and why is alphakid declared as alphakid(nkid) and yet it's really just alphakid(sex) i.e. alphakid(g)
                    !it doesnt' affect anythign I suppose since in all the code it is always alphakid(g) even though it's declared as alphakid(nkid)
                    !since nkid is 2 and number of sex is 2 as well, this doesn't really affect anythign but still it is not right. 
                    !nkid is never above 2 though (i set it that way in data too for numkids above 2 is just always 2 in data). 
                    !alphaed(g,ed(g)):   ed(g)=1 is noed and ed(g)=2 is yes ed (where  neduc is 2)    
                end do !gender
            end do xs
        end do qs
        end do !age
        !close(68857)
        !close(68858)
    
        do age=MNA,MXA
        qc: do q=1,nq
        xc: do x=1,nx
            ed(:)=xx2e(:,x)    
            expe(:)=xx2r(:,x)    
            kid(:)=xx2kid(:,x)           
            w(:) = qq2w(:,q)						! wage 
            l(:) = qq2l(:,q)						! location
            if (l(1).ne.l(2)) then ; print*, 'lm not equal to lf' ; stop ; end if 
    
            !******************************
            !ahu summer18 050318 
            !if (w(1)<=np) then 
            !    ubc(1,q,x,trueindex)	= 0.0_dp           
            !else if (w(1)==np1) then 
            !    epsw(1)=0.0_dp
            !    ubc(1,q,x,trueindex)	= replacement_rate*fnwge(1,truetyp, l(1),epsw(1), ed(1), expe(1) ) 
            !end if 
            !if (w(2)<=np) then 
            !    ubc(2,q,x,trueindex)	= 0.0_dp           
            !else if (w(2)==np1) then 
            !    epsw(2)=0.0_dp
            !    ubc(2,q,x,trueindex)	= replacement_rate*fnwge(2,truetyp, l(2),epsw(2), ed(2), expe(2) ) 
            !end if 
            !ahu summer18 050318 
            !******************************
    
            !pwages(1:numbin)=tax(1:numbin,numbin,l(1))%pwages
            !swages(1:numbin)=tax(numbin,1:numbin,l(2))%swages
            if ( w(1) <= np .and. w(2) <= np ) then		
                epsw(1)=wg(w(1),ed(1),1) !CD(1,1)*wg(w(1),1) !FEB7 2023: BIG SIGLOC(1:NL) CHANGE TO WAGE PROCESS
                epsw(2)=wg(w(2),ed(2),2) !CD(2,1)*wg(w(1),1) + CD(2,2)*wg(w(2),2) !FEB7 2023: BIG SIGLOC(1:NL) CHANGE TO WAGE PROCESS
                wc(1,x,q,age,trueindex)	= fnwge(1,truetyp, l(1),epsw(1), ed(1), expe(1) ,age ) 
                wc(2,x,q,age,trueindex)	= fnwge(2,truetyp, l(2),epsw(2), ed(2), expe(2) ,age )  
                wcgross(1:2)= wc(1:2,x,q,age,trueindex)
                pbrack=locate(  pbracket(1:numbin)  ,  wcgross(1) ) 
                sbrack=locate(  sbracket(1:numbin)  ,  wcgross(2) ) 
            else if ( w(1) <= np .and. w(2) == np1 ) then		
                epsw(1)=wg(w(1),ed(1),1) !sig_wge(1)*wg(w(1),1) !FEB7 2023: BIG SIGLOC(1:NL) CHANGE TO WAGE PROCESS
                wc(1,x,q,age,trueindex)	= fnwge(1,truetyp, l(1),epsw(1), ed(1), expe(1) ,age) 
                wc(2,x,q,age,trueindex)	= 0.0_dp
                wcgross(1:2)= wc(1:2,x,q,age,trueindex)
                pbrack=locate(  pbracket(1:numbin)  ,  wcgross(1) ) 
                sbrack=1
            else if ( w(1) == np1 .and. w(2) <= np ) then		
                epsw(2)=wg(w(2),ed(2),2) !sig_wge(2)*wg(w(2),2) !FEB7 2023: BIG SIGLOC(1:NL) CHANGE TO WAGE PROCESS
                wc(1,x,q,age,trueindex)	= 0.0_dp
                wc(2,x,q,age,trueindex)	= fnwge(2,truetyp, l(2),epsw(2), ed(2), expe(2) ,age)   
                wcgross(1:2)= wc(1:2,x,q,age,trueindex)
                pbrack=1
                sbrack=locate(  sbracket(1:numbin)  ,  wcgross(2) ) 
            else if ( w(1) == np1 .and. w(2) == np1 ) then		
                wc(1,x,q,age,trueindex)	= 0.0_dp
                wc(2,x,q,age,trueindex)	= 0.0_dp           
                wcgross(1:2)= 0.0_dp
                pbrack=1
                sbrack=1
            end if 
    
        
            call get_taxliability_c(l(1),pbrack,sbrack,wcgross(1:2),statetaxtotc,fedtaxtotc)
            wcnet(1:2,x,q,age,trueindex)=wcgross(1:2)- ( sum(statetaxtotc(1:2)) + sum(fedtaxtotc(1:2)) ) 
    
            !if (taxset==0) then
            !    wcnet(1:2,x,q,trueindex)	= wc(1:2,x,q,trueindex)
            !end if 
    
                
            do g=1,2
                if ( w(g) <= np ) then						
                    utilc(g,x,q,age,trueindex)	= uhomet(g) * one(l(g)==truehome)  + uloc(l(g)) + lama(g)*ln( wc(1,x,q,age,trueindex) + wc(2,x,q,age,trueindex) + nonlabinc(x2e(x)) ) + (1.0_dp-lama(g)) * ln(leishours(FULLTIMEWORK))  !+ ulocheat*heat(l(g))
                else if ( w(g) == np1 ) then 
                    !utilc(g,x,q,trueindex)	= uhomet(g) * one(l(g)==truehome)  + uloc(l(g)) + (alphaed(g,ed(g) ) + alphakid(g) ) * one(kid(g)>1)  !+ ulocheat*heat(l(g))
                    utilc(g,x,q,age,trueindex)	= uhomet(g) * one(l(g)==truehome)  + uloc(l(g)) + alphakid(g) * one(kid(g)>1) + alphab(g) + lama(g)*ln( wc(1,x,q,age,trueindex) + wc(2,x,q,age,trueindex) + nonlabinc(x2e(x)) ) + (1.0_dp-lama(g)) * ln(leishours(NOWORK)) !+ ulocheat*heat(l(g))
                    !ahu october2022: 
                    !alphakid(g): kid(g)=1 is no kid, and kid(g)=2 and above is yes kid.   this alphakid used to have two dimensions for no reason. 
                    !nkid is never above 2 though (i set it that way in data too for numkids above 2 is just always 2 in data). 
                    !alphaed(g,ed(g)):   ed(g)=1 is noed and ed(g)=2 is yes ed (where  neduc is 2)
                end if
            end do   
            
        end do xc
        end do qc
        end do !age
    
        end do !trueindex
    
        !Male and female wage shocks for each period and from each location
    !    do it=1,nt
    !	    do iloc=1,nloc
    !		    do rn=1,nn
    !			    up1(it,rn,iloc)=CD(1,1)*rv1(it,rn,iloc)
    !			    up2(it,rn,iloc)=CD(2,1)*rv1(it,rn,iloc)+CD(2,2)*rv2(it,rn,iloc)
    !		    end do 
    !	    end do
    !    end do 
    
        end subroutine get_util_wR



    !subroutine get_taxliability_s(sex,lastbracketrate,gross,lastbracketprod)
    !    integer(i4b), intent(in) :: sex
    !    real(dp), intent(in) :: lastbracketrate,gross
    !    real(dp), intent(out) :: lastbracketprod
    !    if (sex==1) then
    !        lastbracketprod = lastbracketrate * [gross - pwages(pbrack)]    
    !    else if (sex==2) then 
    !        lastbracketprod = lastbracketrate * [gross - pwages(pbrack)]    
    !    end if 
    !end subroutine get_tax_liability_s    

    !subroutine get_taxliability_c(lastbracketrate,gross,lastbracketprod)
    !    real(dp), intent(in), dimension(2) :: lastbracketrate,gross
    !    real(dp), intent(out), dimension(2) :: lastbracketprod
    !    lastbracketprod(1) = lastbracketrate(1) * [gross(1) - pwages(pbrack)]  !husband       
    !    lastbracketprod(2) = lastbracketrate(2) * [gross(2) - swages(sbrack)]  !wife
    !end subroutine get_taxliability_c


    subroutine get_taxliability_s(loc,pp,gross,statetaxtotal,fedtaxtotal)
        integer(i4b), intent(in) :: loc,pp
        real(sp), intent(in) :: gross  !has to be sp because in getutil it is used as input to locate() which needs it sp
        real(dp), intent(out) :: statetaxtotal,fedtaxtotal
        integer(i4b) :: ppmin1,ss
        real(dp) :: lastbracketrate,lastbracketprod
        if (pp==1) ppmin1=1 
        ss=numbin
        !since pwages and swages are the same, the gender should not make a difference here. so just use pwages.
        if (policytax<4) then !tax baseline and the policies that don't affect singles
        
            lastbracketrate=tax(pp,ss,loc)%statesin
            lastbracketprod = lastbracketrate * (gross - pbracket(pp))  
            statetaxtotal = bracketprodsum_s(ppmin1,loc)%state + lastbracketprod
    
            lastbracketrate=tax(pp,ss,loc)%fedsin
            lastbracketprod = lastbracketrate * (gross - pbracket(pp))  
            fedtaxtotal = bracketprodsum_s(ppmin1,loc)%fed + lastbracketprod 

        else if (policytax==4) then !sets state taxes equal to 0
    
            statetaxtotal = 0.0_dp
        
            lastbracketrate=tax(pp,ss,loc)%fedsin
            lastbracketprod = lastbracketrate * (gross - pbracket(pp))
            fedtaxtotal = bracketprodsum_s(ppmin1,loc)%fed + lastbracketprod 
    
        else if (policytax==5) then !sets fed taxes equal to 0
        
            lastbracketrate=tax(pp,ss,loc)%statesin
            lastbracketprod = lastbracketrate * (gross - pbracket(pp))
            statetaxtotal = bracketprodsum_s(ppmin1,loc)%state + lastbracketprod
            
            fedtaxtotal = 0.0_dp
         
        else if (policytax==6) then !sets BOTH state and fed taxes equal to 0

            statetaxtotal = 0.0_dp
            fedtaxtotal = 0.0_dp
         
        end if
        !loc,pbrack,pwages(pbrack),gross(1),sbrack,swages(sbrack),gross(2),tax(pbrack,sbrack,loc)%statemar,taxsum(pbrack,sbrack,loc)%statemar
        !wcnet(1:2,x,q,trueindex)	= (1.0_dp - (staterate+fedrate) )*wcgross(1:2)    
    end subroutine get_taxliability_s





    subroutine get_taxliability_c(loc,pp,ss,gross,statetaxtotal,fedtaxtotal)
        integer(i4b), intent(in) :: loc,pp,ss
        real(sp), dimension(2), intent(in) :: gross(2) !has to be sp because in getutil it is used as input to locate() which needs it sp
        real(dp), dimension(2), intent(out) :: statetaxtotal,fedtaxtotal
        integer(i4b) :: ppmin1,ssmin1 
        real(dp) :: lastbracketrate,lastbracketprod(2)    
        
        if (pp==1) ppmin1=1 
        if (ss==1) ssmin1=1 
        
        if (policytax==0) then !tax baseline

            lastbracketrate=tax(pp,ss,loc)%statemar
            lastbracketprod(1) = lastbracketrate * (gross(1) - pbracket(pp))  !husband       
            lastbracketprod(2) = lastbracketrate * (gross(2) - sbracket(ss))  !wife
            statetaxtotal(1) = bracketprodsum_h(ppmin1,ss,loc)%state + lastbracketprod(1) 
            statetaxtotal(2) = bracketprodsum_w(pp,ssmin1,loc)%state + lastbracketprod(2)
        
            lastbracketrate=tax(pp,ss,loc)%fedmar
            lastbracketprod(1) = lastbracketrate * (gross(1) - pbracket(pp))  !husband       
            lastbracketprod(2) = lastbracketrate * (gross(2) - sbracket(ss))  !wife
            fedtaxtotal(1) = bracketprodsum_h(ppmin1,ss,loc)%fed + lastbracketprod(1) 
            fedtaxtotal(2) = bracketprodsum_w(pp,ssmin1,loc)%fed + lastbracketprod(2)
        
        else if (policytax==1) then !sets state tax rate for mar equal to state tax rate for singles
        
            lastbracketrate=tax(pp,ss,loc)%statesin
            lastbracketprod(1) = lastbracketrate * (gross(1) - pbracket(pp))  !husband       
            lastbracketprod(2) = lastbracketrate * (gross(2) - sbracket(ss))  !wife
            statetaxtotal(1) = bracketprodsum_s(ppmin1,loc)%state + lastbracketprod(1) 
            statetaxtotal(2) = bracketprodsum_s(ssmin1,loc)%state + lastbracketprod(2)
        
            lastbracketrate=tax(pp,ss,loc)%fedmar
            lastbracketprod(1) = lastbracketrate * (gross(1) - pbracket(pp))  !husband       
            lastbracketprod(2) = lastbracketrate * (gross(2) - sbracket(ss))  !wife
            fedtaxtotal(1) = bracketprodsum_h(ppmin1,ss,loc)%fed + lastbracketprod(1) 
            fedtaxtotal(2) = bracketprodsum_w(pp,ssmin1,loc)%fed + lastbracketprod(2)
        
        else if (policytax==2) then !sets fed tax rate for mar equal to fed tax rate for singles
        
            lastbracketrate=tax(pp,ss,loc)%statemar
            lastbracketprod(1) = lastbracketrate * (gross(1) - pbracket(pp))  !husband       
            lastbracketprod(2) = lastbracketrate * (gross(2) - sbracket(ss))  !wife
            statetaxtotal(1) = bracketprodsum_h(ppmin1,ss,loc)%state + lastbracketprod(1) 
            statetaxtotal(2) = bracketprodsum_w(pp,ssmin1,loc)%state + lastbracketprod(2)
        
            lastbracketrate=tax(pp,ss,loc)%fedsin
            lastbracketprod(1) = lastbracketrate * (gross(1) - pbracket(pp))  !husband       
            lastbracketprod(2) = lastbracketrate * (gross(2) - sbracket(ss))  !wife
            fedtaxtotal(1) = bracketprodsum_s(ppmin1,loc)%fed + lastbracketprod(1) 
            fedtaxtotal(2) = bracketprodsum_s(ssmin1,loc)%fed + lastbracketprod(2)
        
        else if (policytax==3) then !sets fed AND state tax rate for mar equal to fed AND state tax rate for singles
        
            lastbracketrate=tax(pp,ss,loc)%statesin
            lastbracketprod(1) = lastbracketrate * (gross(1) - pbracket(pp))  !husband       
            lastbracketprod(2) = lastbracketrate * (gross(2) - sbracket(ss))  !wife
            statetaxtotal(1) = bracketprodsum_s(ppmin1,loc)%state + lastbracketprod(1) 
            statetaxtotal(2) = bracketprodsum_s(ssmin1,loc)%state + lastbracketprod(2)
        
            lastbracketrate=tax(pp,ss,loc)%fedsin
            lastbracketprod(1) = lastbracketrate * (gross(1) - pbracket(pp))  !husband       
            lastbracketprod(2) = lastbracketrate * (gross(2) - sbracket(ss))  !wife
            fedtaxtotal(1) = bracketprodsum_s(ppmin1,loc)%fed + lastbracketprod(1) 
            fedtaxtotal(2) = bracketprodsum_s(ssmin1,loc)%fed + lastbracketprod(2)
        
        else if (policytax==4) then !sets state taxes equal to 0
        
            statetaxtotal = 0.0_dp
        
            lastbracketrate=tax(pp,ss,loc)%fedmar
            lastbracketprod(1) = lastbracketrate * (gross(1) - pbracket(pp))  !husband       
            lastbracketprod(2) = lastbracketrate * (gross(2) - sbracket(ss))  !wife
            fedtaxtotal(1) = bracketprodsum_h(ppmin1,ss,loc)%fed + lastbracketprod(1) 
            fedtaxtotal(2) = bracketprodsum_w(pp,ssmin1,loc)%fed + lastbracketprod(2)
        
        else if (policytax==5) then !sets fed taxes equal to 0 
        
            lastbracketrate=tax(pp,ss,loc)%statemar
            lastbracketprod(1) = lastbracketrate * (gross(1) - pbracket(pp))  !husband       
            lastbracketprod(2) = lastbracketrate * (gross(2) - sbracket(ss))  !wife
            statetaxtotal(1) = bracketprodsum_h(ppmin1,ss,loc)%state + lastbracketprod(1) 
            statetaxtotal(2) = bracketprodsum_w(pp,ssmin1,loc)%state + lastbracketprod(2)
        
            fedtaxtotal = 0.0_dp
        
        else if (policytax==6) then !sets BOTH state and fed taxes equal to 0
        
            statetaxtotal = 0.0_dp
            fedtaxtotal = 0.0_dp
        
        end if
        
        !loc,pbrack,pbracket(pbrack),gross(1),sbrack,sbracket(sbrack),gross(2),tax(pbrack,sbrack,loc)%statemar,taxsum(pbrack,sbrack,loc)%statemar
        !wcnet(1:2,x,q,trueindex)	= (1.0_dp - (staterate+fedrate) )*wcgross(1:2)    
        end subroutine get_taxliability_c
        
end module sol





hold off
 for i=1:9

mu=[(i*(1/10))  (1-i*(1/10)) ] 
mum(i)=mu(1)

lambm=0.5;
lambf=0.3;
w=[5 5];
T=15;
hh=optimvar('hh',1,2);
obj=-( (lambm*log(w(1)*hh(1)+w(2)*hh(2)) + (1-lambm)*log( T-hh(1) )   )^mu(1)) * ( (lambf*log(w(1)*hh(1)+w(2)*hh(2)) + (1-lambf)*log( T-hh(2) )   )^mu(2));
prob = optimproblem('Objective',obj);
show(prob);
hh0.hh=[1 1];
[sol,fval,exitflag,output,xsol] = solve(prob,hh0)
hh=sol.hh
vmaxm(i)= lambm*log(w(1)*hh(1)+w(2)*hh(2)) + (1-lambm)*log( T-hh(1) );   
vmaxf(i)= lambf*log(w(1)*hh(1)+w(2)*hh(2)) + (1-lambf)*log( T-hh(2) );

 end 



 mum1=mum
 vmaxm1=vmaxm
 vmaxf1=vmaxf
 plot(vmaxm,vmaxf,'--gs')
 hold on

 for i=1:9

mu=[(i*(1/10))  (1-i*(1/10)) ] 

lambm=0.5;
lambf=0.3;
w=[8 3];
T=15;
hh=optimvar('hh',1,2);
obj=-( (lambm*log(w(1)*hh(1)+w(2)*hh(2)) + (1-lambm)*log( T-hh(1) )   )^mu(1)) * ( (lambf*log(w(1)*hh(1)+w(2)*hh(2)) + (1-lambf)*log( T-hh(2) )   )^mu(2));
prob = optimproblem('Objective',obj);
show(prob);
hh0.hh=[1 1];
[sol,fval,exitflag,output,xsol] = solve(prob,hh0)
hh=sol.hh
vmaxm(i)= lambm*log(w(1)*hh(1)+w(2)*hh(2)) + (1-lambm)*log( T-hh(1) )   
vmaxf(i)= lambf*log(w(1)*hh(1)+w(2)*hh(2)) + (1-lambf)*log( T-hh(2) )

 end 


 mum1
 vmaxm1
 vmaxf1

 mum
 vmaxm
 vmaxf
 plot(vmaxm,vmaxf,'--gs')
 hold on






lambm=0.3;
lambf=0.3;
w=[10 5];
T=15;
hhm=optimvar('hhm');
obj=- ( lambm*log(w(1)*hhm) + (1-lambm)*log( T-hhm ) );
prob = optimproblem('Objective',obj);
hh0.hhm=1;
show(prob)
[sol,fval,exitflag,output,xsol] = solve(prob,hh0)
hhm=sol.hhm
vmaxm_s=lambm*log(w(1)*hhm) + (1-lambm)*log( T-hhm )






if (m(1)>=m(2) & m(2)>=m(3) )
        if      ( f(1)>=f(2) & f(2)>=f(3) ) %1 ABS ABS  1
            decx(i)=1                       %   BB AA SS       
            decy(i)=1                       %   BB AA SS       
        elseif ( f(1)>=f(3) & f(3)>=f(2))  %2  ABS ASB2      
            decx(i)=1                       %   AA  BS S       
            decy(i)=1                       %   AA  BB S       
        elseif ( f(2)>=f(1) & f(1)>=f(3) ) %3  ABS BAS3      
            decx(i)=1                       %   AA BB SS       
            decy(i)=2                       %   AA BB SS                                        
        end
elseif ( m(1)>=m(3) & m(3)>=m(2) ) 
        if      ( f(1)>=f(2) & f(2)>=f(3) ) %1 ASB ABS  1
            decx(i)=1                       %   BB AA SS       
            decy(i)=1                       %   BS AA SS       
        elseif ( f(1)>=f(3) & f(3)>=f(2))  %2  ASB ASB2      
            decx(i)=1                       %   AA  BS S       
            decy(i)=1                       %   AA  BS S       
        elseif ( f(2)>=f(1) & f(1)>=f(3) ) %3  ASB BAS3      
            decx(i)=1                       %   AA BB SS       
            decy(i)=1                       %   AA BS SS              
        end
elseif ( m(2)>=m(1) & m(1)>=m(3) ) 
        if      ( f(1)>=f(2) & f(2)>=f(3) ) %1 BAS ABS  1
            decx(i)=2                       %   BB AA SS       
            decy(i)=1                       %   BB AA SS       
        elseif ( f(1)>=f(3) & f(3)>=f(2))  %2  BAS ASB2      
            decx(i)=1                       %   AA  BS S       
            decy(i)=0                       %   AA  BB S       
        elseif ( f(2)>=f(1) & f(1)>=f(3) ) %3  BAS BAS3      
            decx(i)=2                       %   AA BB SS       rrrdtx
            decy(i)=2                       %   AA BB SS     
        end
elseif ( m(3)>max(m(1),m(2)) | f(3)>max(f(1),f(2)) )      
    decx(i)=0
    decy(i)=0
end
            
            decx(i)=0                       %   AS BS SS       
            decy(i)=0                       %   SS                              
        elseif ( f(3)>=f(2) & f(2)>=f(1) ) %5  ABS SBA5      
            decx(i)=0                       %   AS BS SS       
            decy(i)=0                       %   SS                              
        elseif ( f(3)>=f(1) & f(1)>=f(2) ) %4  ASB SAB4      
            decx(i)=0                       %   AS BS SS       
            decy(i)=0                       %   SS AA BS                             
        elseif ( f(3)>=f(2) & f(2)>=f(1) ) %5  ASB SBA5      
            decx(i)=0                       %   AS BS SS       
            decy(i)=0                       %   SS BS AA                                           

        elseif ( f(3)>=f(1) & f(1)>=f(2) ) %4  BAS SAB4      
            decx(i)=0                       %   AS BS SS       
            decy(i)=0                       %   SS                              
        elseif ( f(3)>=f(2) & f(2)>=f(1) ) %5  BAS SBA5      
            decx(i)=0                       %   AS BS SS       
            decy(i)=0                       %   SS                                            



SAB ABS  1
SS
SS

SAB ASB  2
SS
SS

SAB BAS 3
SS
SS

SAB SAB   4
SS
SS
SAB SBA  5
SS
SS




