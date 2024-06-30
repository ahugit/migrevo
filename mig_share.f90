module share
	use params
	!use invo !probabilities
	implicit none
	integer(i4b) :: qstart(2,nl),xstart(2,neduc)
	integer(i4b) :: qq2q(2,nq),xx2x(2,nx),q2qq(nqs,nqs),x2xx(nxs,nxs)
	integer(i4b) :: qq2w(2,nq),qq2l(2,nq),xx2e(2,nx),xx2r(2,nx),xx2kid(2,nx)
	integer(i4b) :: q2w(nqs),q2l(nqs),x2e(nxs),x2r(nxs) ,x2kid(nxs)
	integer(i4b) :: wl2q(np2,nl),erkid2x(neduc,nexp,nkid) 
	integer(i4b) :: chs(ncs,nqs,nqs),ch(nc,nq,nq)
	real(dp) :: ppsq(nqs,nqs,ninp,2) , ppcq(nq,nq,ninp)  !, ptemp(nqs,nqs,nxs,2)
	real(dp) :: ppmq_check(nqs,nqs,nxs),ppfq_check(nqs,nqs,nxs)
	real(dp) :: ppsx  (nxs,nqs,nxs) , ppcx(nx,nq,nx) 
	real(dp) :: ppmeetq(nqs,nqs),empdec(mna:mxa,2),pmeet(MNA:MXA)
	real(dp) :: probmeetX(nxs,nxs,ntypp,MALES:FEMALES) 
	real(dp), dimension(2,nxs,nqs,MNA:MXA,ninp) :: ws     !ahu summer18 050318: added ubs
	real(dp), dimension(2,nx,nq,MNA:MXA,ninp) :: wc       !ahu summer18 050318: added ubc
	real(dp), dimension(2,nxs,nqs,MNA:MXA,ninp) :: wsnet	!after-tax income 
	real(dp), dimension(2,nx,nq,MNA:MXA,ninp) :: wcnet 	    !after-tax income 
	real(dp), dimension(2,nxs,nqs,ninp) :: utils 
	real(dp), dimension(2,nx,nq,ninp) :: utilc
    real(dp), dimension(MALES:FEMALES,nx,nz,ninp,MNA:MXA) :: utilmar
	real(dp), dimension(nxs,nqs,ninp,MALES:FEMALES) :: movecost
	!real(dp), dimension(nx,nq,ninp) :: movecostmar_m,movecostmar_f
	real(dp), dimension(nxs,nqs,mna:mxa) :: emaxm_s,emaxf_s
	real(dp), dimension(nx,nq,mna:mxa)   :: emaxm_c,emaxf_c
	!integer(i4b), dimension(nepsmove,nxs,nqs,nqs,mna:mxa,nin) :: decm0_s,decf0_s,decm_s,decf_s
	!real(dp), dimension(nepsmove,nxs,nqs,nqs,mna:mxa,nin) :: vm,vf
    !integer(i4b), dimension(nz,nx,nq,mna:mxa,nin) :: dec_mar   
	!real(dp), dimension(nx,nq,mna:mxa,nin) :: vm0_c,vf0_c	
    !the below is allocatable because nindex (the last dimension) is determined according to groups
    integer(i4b), allocatable, dimension(:,:,:,:,:,:) :: decm_postdiv,decf_postdiv
    integer(i4b), allocatable, dimension(:,:,:,:,:,:) :: decm_premarmkt,decf_premarmkt
	real(dp), allocatable, dimension(:,:,:,:,:,:) :: vm_postdiv,vf_postdiv !need these because simulation calls getdec which needs these as outside options
    integer(i4b), allocatable, dimension(:,:,:,:,:) :: dec_mar   
	real(dp), allocatable, dimension(:,:,:,:) :: vm0_c,vf0_c
	real(dp), allocatable, dimension(:,:) :: vf0ctemp,vm0ctemp	
real(dp), dimension(nq,nx,MNA:MXA,ninp) :: wmctemp,wfctemp
	logical, parameter :: writecho=.FALSE. 
    logical :: ppc(nq,nq),pps(nqs,nqs,2)
    real(dp) :: ppmovesingle(nepsmove),moveshocksin(nepsmove,2),moveshockmar(nepsmove,2)
	real(dp) :: ppmovejoint(nepsmove*nepsmove)
	real(dp) :: marshock(nz),ppmarie(nz)
	type(joint) :: moveshockjoint(nepsmove,nepsmove)
	real(dp) :: wg(np,neduc,MALES:FEMALES),wgt(np),mg(ninp),wgtmen(np),wgtfem(np),wgtjoint(np,np) !these used to be declared in params !FEB7 2023 BIG WG CHANGE FROM WG(NP,2) TO WG(NP,NL)
	real(dp) :: ubenefit_s(nqs),ubenefit_c(nq)
contains 

	subroutine getq2q
	integer(i4b) :: q,qq,i,j,w(2),l(2),z,check
	qq2q=0 ; q2qq=0 ; qq=0 ; qq2w=0 ; qq2l=0 ; wl2q=0
	do z=1,3
	qm: do i=1,nqs
		call q2wloc(i, w(1) , l(1) )			! put all this into arrays below so that you don't have to keep calling this
		q2w(i) = w(1)					! these are all the same for males and females
		q2l(i) = l(1)					! these are all the same for males and females
		wl2q( w(1), l(1) ) = i				! these are all the same for males and females
        !the below is just for checking to make sure wl2q and wloc2q are the same (one is a function and the oth
        call wloc2q(check,w(1),l(1))
        if (check .ne. wl2q( w(1), l(1) ) ) then ; print*, 'something wrong in wl2q' ; stop ; end if
        qf: do j=1,nqs
			call q2wloc(j, w(2) , l(2) )		! do not need to do the above again since q2w,q2l and wl2q are the same for males and females
			if (l(1) == l(2) ) then 
				if ( z==1.and.w(1)==np1.and.w(2)==np1 ) then  
					qq=qq+1
					q2qq(i,j)	= qq 
					qq2q(:,qq)	= (/ i , j /) 
					qq2w(:,qq)	= w(:)
					qq2l(:,qq)	= l(:)
					if (skriv) then 
						write(40,10) 
						write(40,20) i,w(1),l(1),j,w(2),l(2),qq
					end if 
				else if (z==2.and.(w(1)<=np.or.w(2)<=np)) then 
					qq=qq+1
					q2qq(i,j)	= qq 
					qq2q(:,qq)	= (/ i , j /) 
					qq2w(:,qq)	= w(:)
					qq2l(:,qq)	= l(:)
					if ( skriv) then 
						write(40,10) 
						write(40,20) i,w(1),l(1),j,w(2),l(2),qq
					end if 
				else if (z==3.and.(w(1)>=np1.and.w(2)>=np1).and.(w(1)==np2.or.w(2)==np2) ) then 
					qq=qq+1
					q2qq(i,j)	= qq 
					qq2q(:,qq)	= (/ i , j /) 
					qq2w(:,qq)	= w(:)
					qq2l(:,qq)	= l(:)						
					if (skriv) then 
						write(40,10) 
						write(40,20) i,w(1),l(1),j,w(2),l(2),qq
					end if 					
				end if				 
			end if 
		end do qf
	end do qm
	end do 
	10 format (/ 1x,tr5,'qm',tr2,'w1',tr2,'l1',tr6,'qf',tr2,'w2',tr2,'l2',tr6,'qq' /) 
	20 format (i8,2i4,i8,2i4,i8)
	if (skriv) then 
		if ( qq /= nq ) then ; print*, "qq not equal to nq! ",qq,nq ; stop ; end if 
		if ( qq /= np2 * np2 * nl ) then ; print*, "qq not equal to npr*npr*nl! ",qq,np2*np2*nl ; stop ; end if 
	!	print*, "done with getq2q "
	end if
	if (skriv) then 
		do q=1,nqs
			i=q2w(q)
			j=q2l(q)
			if (wl2q(i,j)/=q) then 
				print*, "in getq2q: wl2q does not work ",q,wl2q(i,j),i,j
				stop 
			end if 
		end do 
		do qq=1,nq
			w=qq2w(:,qq)
			l=qq2l(:,qq)
			i=qq2q(1,qq)
			j=qq2q(2,qq)			
			if (wl2q(w(1),l(1))/=i .or. wl2q(w(2),l(2))/=j) then 
				print*, "in getq2q: qq2w,qq2l,wl2q does not work "  !,q,wl2q(i,j),i,j
				stop 
			end if 
			if ( q2qq(i,j) /= qq) then 
				print*, "in getq2q: q2qq does not work "   !,q,wl2q(i,j),i,j
				stop 
			end if 
		end do 
	end if 
	end subroutine getq2q

	subroutine getx2x
	integer(i4b) :: x,xx,i,j,e(2),r(2),kid(2)		!e is education and r is experience 
	xx2x=0 ; x2xx=0 ; xx=0 ; xx2e=0 ; xx2r=0 
	xm: do i=1,nxs
		call x2edexpkid(i, e(1) , r(1) , kid(1) )
		x2e(i) = e(1) 
		x2r(i) = r(1)
        x2kid(i) = kid(1)
		erkid2x( e(1), r(1) , kid(1) ) = i 
		xf: do j=1,nxs
			call x2edexpkid(j, e(2) , r(2) , kid(2) )
			xx=xx+1
			x2xx(i,j)	= xx 
			xx2x(:,xx)	= (/ i , j /)
			xx2e(:,xx)	= e(:)		! education
			xx2r(:,xx)	= r(:)		! experience
			xx2kid(:,xx)= kid(:)    ! kid          
		end do xf
	end do xm
	if (skriv) then 
		if ( xx /= nx ) then ; print*, "xx not equal to nx! ",xx,nx ; stop ; end if 
		if ( xx /= neduc*nexp*nkid*neduc*nexp*nkid ) then ; print*, "xx not equal to neduc*nexp*nkid*neduc*nexp*nkid! " ; stop ; end if 
	!	print*, "done with getx2x "
	end if 
	end subroutine getx2x

	subroutine getppsq
	! calculates transition probs between q0 and q. getgauss needs to be called before this as you need wgt
	integer(i4b) :: g,i,j,r,w,l,w0,l0,index
	real(dp):: prof(3),prloc(nl),dum(np2,nl)
	ppsq=0.0_dp 
	sex: do g=1,2					! sex
		!x0: do n=1,nxs				! x0 !this has not been valid for a while. in other words, since prof doesn't chge with ed anymore, we haven't needed this for a while
		!	e=x2e(n)			! e is education 
		do index=1,ninp 		!to incorporate home location in prloc 
		q0: do i=1,nqs			! q0
				if ( q2w(i)<=np1 ) then	!state variable part of the q space i.e. w /= np2
					w0=q2w(i) 
					l0=q2l(i)
					!prof=abs(fnprof(w0,e,1)-fnprof(w0,e,2)) 
					!if ( icheck_eqvmvf.and.maxval(prof)>eps ) then 
					!	print*, "not even ",w,e,prof
					!	!stop
					!end if 
					prloc(:) = fnprloc(l0,index) 				 
					q: do j=1,nqs
						w=q2w(j) 
						l=q2l(j)	
                        !ahu jan19 012819: separating the offer probs by whether it's curloc or ofloc. no more ed. 
                        if (q2l(j)==q2l(i) ) then   !if curloc
    					    prof(:)  = fnprof (w0,5,g)
                        else                        !if ofloc
    					    prof(:)  = fnprof (w0,10,g)
                        end if                         

                        !ahu jan19 011219: no on the job seearch since not identified
                        !if (onthejobsearch ) then 
                        !proftemp=prof
                        !else
                        !    if (q2w(i)<=np ) then           !if working
                        !        if (q2l(j)==q2l(i) ) then   !and the offer is from curloc
                        !            proftemp(1)=0.0_dp ; proftemp(2)=prof(2) ; proftemp(3)=prof(1)+prof(3)  !then there is no job offer so that p(offer)=0 and p(layoff) is just the same i.e. prof(layoff) and p(nothing) is prof(offer)+prof(nothing)
                        !        else                        !and the offer is from ofloc
                        !            proftemp(:)=prof(:)     !then there is potentially a job offer
                        !        end if 
                        !    else                            !if not working, no such complications. can get job offer regardless of curloc or ofloc
                        !        proftemp(:)=prof(:)
                        !    end if
                        !end if !onthejobsearch
                        !ahu jan19 011219: no on the job seearch since not identified
                        
						if	( w <= np ) then 
							if (g==1) then
								dum(w,l) = prloc(l) * prof(1) * wgtmen(w)  !ahu jan19 011219:
							else 
								dum(w,l) = prloc(l) * prof(1) * wgtfem(w)  
							end if 
						else if ( w == np1 ) then 
							dum(w,l) = prloc(l) * prof(2)           !ahu jan19 011219:
						else if ( w == np2 ) then 
							dum(w,l) = prloc(l) * prof(3)           !ahu jan19 011219:
						end if 				
                        ppsq(j,i,index,g) = dum(w,l) !ahu040917 prob changes 
					end do q 
					!ahu040917 prob changes ppsq(:,i,n,g) = reshape( dum, (/ np2 * nl /) ) 
					if ( q2w(i) == np2 ) then ; print*, "error in ppsq: q2w does not agree with ps0 " ; stop ; end if  
                    if ( abs( sum(ppsq(:,i,index,g))-1.0_dp   ) > eps   ) then ; print*, 'ppsq does not sum to 1',sum(ppsq(:,i,index,g)) ; stop ; end if 
				else 
					ppsq(:,i,index,g) = pen 
				end if 									
			end do	q0
		end do !index

		pps(:,:,g)=.false.			
		do i=1,nqs
			do j=1,nqs				
				if (  q2w(i)<=np1 .and. maxval(ppsq(j,i,:,g)) > 0.0_dp ) then 
					pps(j,i,g) = .true.
				end if 
			end do 
		end do 
	end do sex 

	if (icheck_eqvmvf) then 
	do index=1,ninp
		do i=1,nqs
			do j=1,nqs
				if ( abs(ppsq(j,i,index,1)-ppsq(j,i,index,2))>eps ) then 
					print*, "m f par same but ppsq not equal! ", index,i,j,ppsq(j,i,index,:)
					stop
				end if 
			end do 
		end do 
	end do 
	end if 
	!if (skriv) print*, "done with getppsq "
	end subroutine getppsq

	subroutine getppcq
	integer(i4b) :: i,j,is(2),js(2),index
    integer(i4b) :: w0(2),l0(2),w(2),l(2)
    real(dp) :: profm(3),proff(3),prloc(nl)
    
    !ahu 040917 prob changes: 
	ppcq=0.0_dp 	
	!x0: do n=1,nx		!this has not been valid for a while. in other words, since prof doesn't chge with ed anymore, we haven't needed this for a while
        !ns(:)=xx2x(:,n)
        !e0(1)=x2e( ns(1) )
        !e0(2)=x2e( ns(2) )
	do index=1,ninp	!for changing prloc to incorporate home location
		q0: do i=1,nq	
        
			if ( maxval(qq2w(:,i)) <= np1 ) then !state variable part of the q space i.e. w <= np1

                is(:)=qq2q(:,i)   
			    w0(1)=q2w( is(1) )
			    w0(2)=q2w( is(2) )       
			    l0(1)=q2l( is(1) )
			    l0(2)=q2l( is(2) )       
                !profm(:)  = fnprof (w0(1),e0(1),1)
                !proff(:)  = fnprof (w0(2),e0(2),2)
                prloc(:) = fnprloc(l0(1),index) 				   
                if (l0(1).ne.l0(2)) then ; print*, 'l0 not equal' ; stop ; end if 
                q: do j=1,nq
                    js(:)=qq2q(:,j)   
			        w(1)=q2w( js(1) )
			        w(2)=q2w( js(2) )       
			        l(1)=q2l( js(1) )
			        l(2)=q2l( js(2) )
                    
                    !ahu jan19 012819: separating the offer probs by whether it's curloc or ofloc. no more ed. 
                    if (q2l(js(1))==q2l(is(1)) ) then   !if curloc
    					profm(:)  = fnprof (w0(1),5,1)
    					proff(:)  = fnprof (w0(2),5,2)
                    else                                !if ofloc
    					profm(:)  = fnprof (w0(1),10,1)
    					proff(:)  = fnprof (w0(2),10,2)
                    end if                         

                    !ahu jan19 011219: no on the job seearch since not identified
                    !if (onthejobsearch ) then
                    !profmtemp=profm
                    !profftemp=proff
                    !else
                    !    if ( w0(1) <=np ) then           !if working
                    !        if (l(1)==l0(1) ) then       !and the offer is from curloc
                    !            profmtemp(1)=0.0_dp ; profmtemp(2)=profm(2) ; profmtemp(3)=profm(1)+profm(3)  !then there is no job offer so that p(offer)=0 and p(layoff) is just the same i.e. prof(layoff) and p(nothing) is prof(offer)+prof(nothing)
                    !        else                        !and the offer is from ofloc
                    !            profmtemp(:)=profm(:)     !then there is potentially a job offer
                    !        end if 
                    !    else                            !if not working, no such complications. can get job offer regardless of curloc or ofloc
                    !        profmtemp(:)=profm(:)
                    !    end if
                    !    if ( w0(2) <=np ) then           !if working
                    !        if (l(2)==l0(2) ) then       !and the offer is from curloc
                    !            profftemp(1)=0.0_dp ; profftemp(2)=proff(2) ; profftemp(3)=proff(1)+proff(3)  !then there is no job offer so that p(offer)=0 and p(layoff) is just the same i.e. prof(layoff) and p(nothing) is prof(offer)+prof(nothing)
                    !        else                        !and the offer is from ofloc
                    !            profftemp(:)=proff(:)     !then there is potentially a job offer
                    !        end if 
                    !    else                            !if not working, no such complications. can get job offer regardless of curloc or ofloc
                    !        profftemp(:)=proff(:)
                    !    end if
                    !end if !onthejobsearch
                    !ahu jan19 011219: no on the job seearch since not identified

                    if (l(1).ne.l(2)) then ; print*, 'l not equal' ; stop ; end if 
                    
					IF (JOINTDIST) THEN

						if	( w(1) <= np .and. w(2)<=np ) then !both get offer
							ppcq(j,i,index) = prloc(l(1)) * profm(1) * wgtjoint(w(1),w(2)) * proff(1)
						else if	( w(1) <= np .and. w(2)==np1 ) then !he gets offer, she gets laid off
							ppcq(j,i,index) = prloc(l(1)) * profm(1) * wgtmen(w(1)) * proff(2)
						else if	( w(1) <= np .and. w(2)==np2 ) then !he gets offer, nothing happens for her
							ppcq(j,i,index) = prloc(l(1)) * profm(1) * wgtmen(w(1)) * proff(3)
							
						else if	( w(1) == np1 .and. w(2)<=np ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(2) * proff(1) * wgtfem(w(2))
						else if ( w(1) == np1 .and. w(2)==np1 ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(2) * proff(2) 
						else if ( w(1) == np1 .and. w(2)==np2 ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(2) * proff(3) 
						
						else if	( w(1) == np2 .and. w(2)<=np ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(3) * proff(1) * wgtfem(w(2))
						else if ( w(1) == np2 .and. w(2)==np1 ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(3) * proff(2) 
						else if ( w(1) == np2 .and. w(2)==np2 ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(3) * proff(3) 
						end if 		
					ELSE 
						if	( w(1) <= np .and. w(2)<=np ) then !both get offer
							ppcq(j,i,index) = prloc(l(1)) * profm(1) * wgtmen(w(1)) * proff(1) * wgtfem(w(2))
						else if	( w(1) <= np .and. w(2)==np1 ) then !he gets offer, she gets laid off
							ppcq(j,i,index) = prloc(l(1)) * profm(1) * wgtmen(w(1)) * proff(2)
						else if	( w(1) <= np .and. w(2)==np2 ) then !he gets offer, nothing happens for her
							ppcq(j,i,index) = prloc(l(1)) * profm(1) * wgtmen(w(1)) * proff(3)
							
						else if	( w(1) == np1 .and. w(2)<=np ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(2) * proff(1) * wgtfem(w(2))
						else if ( w(1) == np1 .and. w(2)==np1 ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(2) * proff(2) 
						else if ( w(1) == np1 .and. w(2)==np2 ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(2) * proff(3) 
						
						else if	( w(1) == np2 .and. w(2)<=np ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(3) * proff(1) * wgtfem(w(2))
						else if ( w(1) == np2 .and. w(2)==np1 ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(3) * proff(2) 
						else if ( w(1) == np2 .and. w(2)==np2 ) then 
							ppcq(j,i,index) = prloc(l(1)) * profm(3) * proff(3) 
						end if 		
                    END IF

                end do q          
				if (maxval(qq2w(:,i)) == np2) then ; print*, "error in ppcq: qq2w does not agree with pc0 " ; stop ; end if  	
                if ( abs( sum(ppcq(:,i,index))-1.0_dp   ) > eps   ) then ; print*, 'ppcq does not sum to 1',sum(ppcq(:,i,index)) ; stop ; end if 
			else 
				ppcq(:,i,index)=pen 
			end if 
                     
		end do q0
	end do !index

	!ahu 040917 prob changes ppcq=0.0_dp		
	!ahu 040917 prob changes x0: do n=1,nx				
    !ahu 040917 prob changes ns=xx2x(:,n)
		!ahu 040917 prob changes q0: do i=1,nq	
			!ahu 040917 prob changes if ( maxval(qq2w(:,i)) <= np1 ) then !state variable part of the q space i.e. w <= np1
				!ahu 040917 prob changes is=qq2q(:,i)
				!ahu 040917 prob changes q: do j=1,nq
					!ahu 040917 prob changes js=qq2q(:,j)	
					!ahu 040917 prob changes ppcq(j,i,n) = ppsq(js(1),is(1),ns(1),1) * ppsq(js(2),is(2),ns(2),2) 	! normalize: i.e. taking into account that they only receive offer from same loc (by modeling restriction, not because its prob is 0) 					
                                                                                            ! the q for couples just does not include points where the l is different for the partners. so this is how this restriction is imposed. 
				!ahu 040917 prob changes end do q			 
				!ahu 040917 prob changes ppcq(:,i,n)=ppcq(:,i,n)/sum( ppcq(:,i,n) )
				!ahu 040917 prob changes if (maxval(qq2w(:,i)) == np2) then ; print*, "error in ppcq: qq2w does not agree with pc0 " ; stop ; end if  			
			!ahu 040917 prob changes else 
				!ahu 040917 prob changes ppcq(:,i,n)=pen 
			!ahu 040917 prob changes end if 
		!ahu 040917 prob changes end do q0
	!ahu 040917 prob changes end do x0
	!if ( q2l(qm) == q2l(qf) ) then 
	!	tmpsum( q2qq(qm,qf) ) = tmpsum( q2qq(qm,qf) ) + ppsq(qm,is(1),ns(1),1) * ppsq(qf,is(2),ns(2),2) 						
	!end if 
	!if ( maxval( qq2w(:,j) ) >= np1 ) then		! they can't both get offers
	! by assumption, if they do get offers at all, both spouses' offer location is the same
	! incorporate this as a conditioning statement when calculating the joint offer/loc probabilities of castanza
	! this is more transparent than the ptemp and prloc calculation that we were doing in the previous versions  i.e. ppcq(j,i,n) = ptemp(js(1),is(1),ns(1),1) * ptemp(js(2),is(2),ns(2),2) * prloc(loc)
	! the above boils down to calculating the prob of the conditioning statement and dividing by it to get the conditional probability	
	!defining ppc below just so that we don't need to use x0 in the big loop when trying to save time by skipping the zero prob q points
	ppc=.false.
	do i=1,nq
		do j=1,nq
			if ( maxval(qq2w(:,i)) <= np1 .and.  maxval(ppcq(j,i,:)) > 0.0_dp  ) then 
				ppc(j,i) = .true. 
			end if 
		end do 
	end do 	
	!if (skriv) print*, "done with getppcq "
	end subroutine getppcq

	subroutine getppsx
	real(dp), dimension(nexp) :: prhc
	integer(i4b) :: n,i,j,e0,e1,r0,r1,w0,l,kid0,kid1
    real(dp) :: probkid
	prhc=0.0_dp
	ppsx=0.0_dp 	
	x0: do n=1,nxs	
		e0=x2e(n)
		r0=x2r(n)
        kid0=x2kid(n)
		q0: do i=1,nqs				
            if ( q2w(i)<=np1 ) then	!state variable part of the q space i.e. w /= np2

                w0=q2w(i)
                prhc(1:nexp) = fnprhc( r0 , w0 ) 
                
			    x: do j=1,nxs	
				    e1=x2e(j)
				    r1=x2r(j)
                    kid1=x2kid(j)
                    probkid=-99.0_dp
                    if (kid0==1.and.kid1==1) then
                        probkid=1-pkidsingle !if single, can't have kids jan2023
                    else if (kid0==1.and.kid1==2) then
                        probkid=pkidsingle 
                    else if (kid0==2.and.kid1==1) then
                        probkid=0.0_dp
                    else if (kid0==2.and.kid1==2) then
                        probkid=1.0_dp
                    end if          
                    if (probkid<0) then ; print*, 'probkid is negative in getppsx' ; stop ; end if 
				    !ppsx(j,i,n) = prhc(r1) * one(e0==e1) * one(kid0==kid1)	! trans prob is 0 if educ is not equal. it is also 0 if kid state is not equal. 
			        ppsx(j,i,n) = prhc(r1) * one(e0==e1) * probkid	! trans prob is 0 if educ is not equal. it is also 0 if kid state is not equal. 
                end do x 
                
                if ( abs( sum(ppsx(:,i,n))-1.0_dp   ) > eps   ) then ; print*, 'ppsx does not sum to 1',sum(ppsx(:,i,n)) ; stop ; end if 
                
			else 
				ppsx(:,i,n)=pen 
			end if 
                
        end do q0
	end do 	x0
	!if (skriv) print*, "done with getppsx "
	end subroutine getppsx
    
	subroutine getppcx
	integer(i4b) :: n,i,j,ns(2),is(2),js(2)
    integer(i4b) :: e0(2),r0(2),kid0(2),w0(2)
    integer(i4b) :: e1(2),r1(2),kid1(2)
    real(dp) :: probkid,prhc_m(nexp),prhc_f(nexp)
	ppcx=0.0_dp 	
	x0: do n=1,nx				
        ns(:)=xx2x(:,n)
        e0(1)=x2e( ns(1) )
		r0(1)=x2r( ns(1) )
        kid0(1)=x2kid( ns(1) )
        e0(2)=x2e( ns(2) )
		r0(2)=x2r( ns(2) )
        kid0(2)=x2kid( ns(2) )        
		q0: do i=1,nq	
			if ( maxval(qq2w(:,i)) <= np1 ) then !state variable part of the q space i.e. w <= np1

                is(:)=qq2q(:,i)   
			    w0(1)=q2w( is(1) )
			    w0(2)=q2w( is(2) )       
 
                prhc_m(1:nexp) = fnprhc( r0(1) , w0(1) ) 
                prhc_f(1:nexp) = fnprhc( r0(2) , w0(2) ) 
            
			    x: do j=1,nx
				    js(:)=xx2x(:,j)
				    e1(1)=x2e( js(1) )
				    r1(1)=x2r( js(1) )
                    kid1(1)=x2kid( js(1) )
				    e1(2)=x2e( js(2) )
				    r1(2)=x2r( js(2) )
                    kid1(2)=x2kid( js(2) )
                    probkid=-99.0_dp
                    if ( kid0(1) == 1 .and. kid0(2)==1 ) then 
                        if (kid1(1)==1 .and. kid1(2)==1 ) then 
                                probkid=1.0_dp-pkidmar
                        else if (kid1(1)==2 .and. kid1(2)==1 ) then 
                                probkid=0.0_dp
                        else if (kid1(1)==1 .and. kid1(2)==2 ) then 
                                probkid=0.0_dp
                        else if (kid1(1)==2 .and. kid1(2)==2 ) then 
                                probkid=pkidmar
                        end if 
                    else if ( kid0(1) == 2 .and. kid0(2)==1 ) then 
                        !IMPORTANT! ASSIGNING PROBKID HERE BUT THIS REALLY SHOULD NOT BE A REACHABLE STATE SPACE POINT SINCE I AM SETTING PPMEET(XM,XF) IN SUCH A WAY THAT IF YOU HAVE KID, THE PROB OF MEETING SOMEONE WITHOUT KID IS 0. 
                        if (kid1(1)==1 .and. kid1(2)==1 ) then 
                                probkid=0.0_dp
                        else if (kid1(1)==2 .and. kid1(2)==1 ) then 
                                probkid=1.0_dp-pkidmar
                        else if (kid1(1)==1 .and. kid1(2)==2 ) then 
                                probkid=0.0_dp
                        else if (kid1(1)==2 .and. kid1(2)==2 ) then 
                                probkid=pkidmar
                        end if 
                    else if ( kid0(1) == 1 .and. kid0(2)==2 ) then 
                        !IMPORTANT! ASSIGNING PROBKID HERE BUT THIS REALLY SHOULD NOT BE A REACHABLE STATE SPACE POINT SINCE I AM SETTING PPMEET(XM,XF) IN SUCH A WAY THAT IF YOU HAVE KID, THE PROB OF MEETING SOMEONE WITHOUT KID IS 0. 
						if (kid1(1)==1 .and. kid1(2)==1 ) then 
                                probkid=0.0_dp
                        else if (kid1(1)==2 .and. kid1(2)==1 ) then 
                                probkid=0.0_dp
                        else if (kid1(1)==1 .and. kid1(2)==2 ) then 
                                probkid=1.0_dp-pkidmar
                        else if (kid1(1)==2 .and. kid1(2)==2 ) then 
                                probkid=pkidmar
                        end if 
                    else if ( kid0(1) == 2 .and. kid0(2)==2 ) then 
                        if (kid1(1)==1 .and. kid1(2)==1 ) then 
                                probkid=0.0_dp
                        else if (kid1(1)==2 .and. kid1(2)==1 ) then 
                                probkid=0.0_dp
                        else if (kid1(1)==1 .and. kid1(2)==2 ) then 
                                probkid=0.0_dp
                        else if (kid1(1)==2 .and. kid1(2)==2 ) then 
                                probkid=1.0_dp
                        end if 
                    end if 
                    if (probkid<0) then ; print*, 'something wrong with probkid',probkid,kid0,kid1 ; stop ; end if 
                    !ppcx(j,i,n) = ppsx( js(1), is(1), ns(1) ) * ppsx( js(2), is(2), ns(2) )
                    ppcx(j,i,n) = prhc_m( r1(1) ) * prhc_f( r1(2) ) * one( e0(1) == e1(1) ) * one( e0(2) == e1(2) ) * probkid  
			    end do x
                
                if ( abs( sum(ppcx(:,i,n))-1.0_dp   ) > eps   ) then ; print*, 'ppcx does not sum to 1',sum(ppcx(:,i,n)) ; stop ; end if 
			else 
				ppcx(:,i,n)=pen 
			end if 

		end do q0
	end do x0
	!if (skriv) print*, "done with getppcx "
	end subroutine getppcx

	subroutine getch_single_bignc
		! choice vectors that indicate whether an alternative is feasible given the status quo and given the offer situation
		! w=np+1 --> this is layoff if you're working and otherwise it's nothing. and choosing to go anywhere unemp is always an option in any case
		! w=np+2 --> this is just the "nothing happens" case and choosing that doesn't make sense 
		integer(i4b) :: i,j,w0,w,c,l0,l,dum(2,nl)
		chs=0 ; c=0
		do i=1,nl							! choice index
			c=c+1
			w=np1							! w=np+1 denotes unemployment 
			j=wl2q(w,i)						! get what choosing unemployment (w=np+1) at location l corresponds to in terms of the composite index q
			chs(c,:,:)=J !THIS IS A HUGE CHANGE ahu october2022					! moving anywhere unemployed is always an option					
			dum(1,c)=w
			dum(2,c)=i
		end do 	
		if (skriv) then 		
			if ( c /= nl ) then ; print*, "c is not equal to nl! ",c,nl ; stop ; end if 		
		end if 
		!c=1 !HUGE MAJOR CHANGE AHU OCTOBER2022 ahu october2022
		q0: do i=1,nqs
			w0 = q2w(i)
			l0 = q2l(i)
			if (w0<=np1) then 
				q: do j=1,nqs					
					!call q2wloc(j,w,l) 
					!chs(c,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option
					w = q2w(j)
					l = q2l(j)
					
					if (w0<=np) then        !currently employed 
						if (l==l0) then         !curloc draw
							!chs(c,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option
							if (w<=np) then         !gets wage offer
								chs(c+1,j,i) = i    !status quo
								chs(c+2,j,i) = j    !accept offer
							else if (w==np1) then   !gets laid off
								chs(:,j,i) = 0      !when gets laid off, nothing is feasible other than unemployment at l0
								chs(l0,j,i) = wl2q(np1,l0)	     !when gets laid off, nothing is feasible other than unemployment at l0
								chs(c+1,j,i) = 0    !when gets laid off, status quo not feasible
								chs(c+2,j,i) = 0    !when gets laid off, no offer to accept
							else if (w==np2) then   !nothing happens
								chs(l0,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option !setting this with chs(c,j,i)=wl2q(np1,l0) instead
								chs(c+1,j,i) = i    !status quo
								chs(c+2,j,i) = 0    !when nothing happens, no offer to accept
							end if
						else if (l/=l0) then         !ofloc draw
							!chs(c,j,i) = wl2q(np1,l)	     !moving to ofloc unemp is an option 
							if (w<=np) then         !gets wage offer
								!chs(l0,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option !setting this with chs(c,j,i)=wl2q(np1,l0) instead
								chs(c+1,j,i) = i    !status quo
								chs(c+2,j,i) = j    !accept offer
							else if (w==np1) then   !gets laid off !BUT NOTE THT THIS DOES NOT HAPPEN WHEN THE DRAW IS OFLOC
								chs(c+1,j,i) = 0    !when gets laid off, status quo not feasible
								chs(c+2,j,i) = 0    !when gets laid off, no offer to accept
							else if (w==np2) then   !nothing happens
								chs(c+1,j,i) = i    !status quo
								chs(c+2,j,i) = 0    !when nothing happens, no offer to accept
							end if
						end if 
					else if (w0==np1) then      !currently unemployed 
						if (l==l0) then         !curloc draw
							!chs(c,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option
							if (w<=np) then         !gets wage offer
								chs(c+1,j,i) = i    !status quo
								chs(c+2,j,i) = j    !accept offer
							else if (w==np1) then   !gets laid off BUT NOTE THAT THIS DOES NOT HAPPEN WHEN UNEMPLOYED
								chs(c+1,j,i) = 0    !when gets laid off, status quo not feasible
								chs(c+2,j,i) = 0    !when gets laid off, no offer to accept
							else if (w==np2) then   !nothing happens
								chs(c+1,j,i) = i    !status quo
								chs(c+2,j,i) = 0    !when nothing happens, no offer to accept
							end if
						else if (l/=l0) then         !ofloc draw
							!chs(c,j,i) = wl2q(np1,l)	     !moving to ofloc unemp is an option 
							if (w<=np) then         !gets wage offer
								chs(c+1,j,i) = i    !status quo
								chs(c+2,j,i) = j    !accept offer
							else if (w==np1) then   !gets laid off  BUT NOTE THAT THIS DOES NOT HAPPEN WHEN UNEMPLOYED
								chs(c+1,j,i) = 0    !when gets laid off, status quo not feasible
								chs(c+2,j,i) = 0    !when gets laid off, no offer to accept
							else if (w==np2) then   !nothing happens
								chs(c+1,j,i) = i    !status quo
								chs(c+2,j,i) = 0    !when nothing happens, no offer to accept
							end if
						end if 
					end if     
							
	
					!if ( w /= np1 .and. w0 /= np1 ) then	! w=np+1 denotes getting laid off (in the case where q denotes shocks) 
										! and w0=np+1 denotes the state variable w0 indicating that you start that period as unemployed 
					!	chs(c+1,j,i) = i		! c+1 is nl+1 which denotes the status quo option. can choose the status quo only if you don't get laid off. also don't need to allow this choice if the status quo is unemployment 	
					!end if 
					!if ( w <= np ) then			! offer is only the cases where w<=np. in other offer situation where w>np, it means either that you get laid off or nothing happens, so those things do not really give you any options. 
					!	chs(c+2,j,i) = j		! c+2 is nl+2 which denotes the option of taking the offer q
					!end if 		
				end do q
			else 
				chs(:,:,i)=0
			end if 
		end do q0
		!if (skriv) print*, "done with getch_single "
	end subroutine 

	subroutine getch_single
	! choice vectors that indicate whether an alternative is feasible given the status quo and given the offer situation
	! w=np+1 --> this is layoff if you're working and otherwise it's nothing. and choosing to go anywhere unemp is always an option in any case
	! w=np+2 --> this is just the "nothing happens" case and choosing that doesn't make sense 
	integer(i4b) :: i,j,w0,w,c,l0,l,dum(2,nl),qcho
	if (mysay==0.and.writecho) then
		open(unit=89898688, file='chos.txt',status='replace')
	end if 
	chs=0 ; c=0
	!do i=1,nl							! choice index
	!	c=c+1
	!	w=np1							! w=np+1 denotes unemployment 
	!	j=wl2q(w,i)						! get what choosing unemployment (w=np+1) at location l corresponds to in terms of the composite index q
	!	chs(c,:,:)=J !THIS IS A HUGE CHANGE ahu october2022					! moving anywhere unemployed is always an option					
	!	dum(1,c)=w
	!	dum(2,c)=i
	!end do 	
	!if (skriv) then 		
	!	if ( c /= nl ) then ; print*, "c is not equal to nl! ",c,nl ; stop ; end if 		
	!end if 
	c=1 !HUGE MAJOR CHANGE AHU OCTOBER2022 ahu october2022
	q0: do i=1,nqs
		w0 = q2w(i)
		l0 = q2l(i)
		if (w0<=np1) then 
			q: do j=1,nqs					
				!call q2wloc(j,w,l) 
				chs(1,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option
				w = q2w(j)
				l = q2l(j)
                
                if (w0<=np) then        !currently employed 
                    if (l==l0) then         !curloc draw
                        if (w<=np) then         !gets wage offer
							chs(1,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option
                            chs(2,j,i) = i    !status quo
                            chs(3,j,i) = j    !accept offer
                        else if (w==np1) then   !gets laid off
                            !chs(:,j,i) = 0      !when gets laid off, nothing is feasible other than unemployment at l0
                            !chs(l0,j,i) = wl2q(np1,l0)	     !when gets laid off, nothing is feasible other than unemployment at l0
							chs(1,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option
                            chs(2,j,i) = 0    !when gets laid off, status quo not feasible
                            chs(3,j,i) = 0    !when gets laid off, no offer to accept
                        else if (w==np2) then   !nothing happens
                            !chs(l0,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option !setting this with chs(c,j,i)=wl2q(np1,l0) instead
							chs(1,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option
							chs(2,j,i) = i    !status quo
                            chs(3,j,i) = 0    !when nothing happens, no offer to accept
                        end if
                    else if (l/=l0) then         !ofloc draw
                        if (w<=np) then         !gets wage offer
                            !chs(l0,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option !setting this with chs(c,j,i)=wl2q(np1,l0) instead
							chs(1,j,i) = wl2q(np1,l)	     !moving to ofloc unemp is an option 
                            chs(2,j,i) = i    !status quo
                            chs(3,j,i) = j    !accept offer
                        else if (w==np1) then   !gets laid off !BUT NOTE THT THIS DOES NOT HAPPEN WHEN THE DRAW IS OFLOC
							chs(1,j,i) = wl2q(np1,l)	     !moving to ofloc unemp is an option 
							chs(2,j,i) = 0    !when gets laid off, status quo not feasible
                            chs(3,j,i) = 0    !when gets laid off, no offer to accept
                        else if (w==np2) then   !nothing happens
							chs(1,j,i) = wl2q(np1,l)	     !moving to ofloc unemp is an option 
							chs(2,j,i) = i    !status quo
                            chs(3,j,i) = 0    !when nothing happens, no offer to accept
                        end if
                    end if 
                else if (w0==np1) then      !currently unemployed 
                    if (l==l0) then         !curloc draw
                        if (w<=np) then         !gets wage offer
							chs(1,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option
							chs(2,j,i) = i    !status quo
                            chs(3,j,i) = j    !accept offer
                        else if (w==np1) then   !gets laid off BUT NOTE THAT THIS DOES NOT HAPPEN WHEN UNEMPLOYED
							chs(1,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option
                            chs(2,j,i) = 0    !when gets laid off, status quo not feasible
                            chs(3,j,i) = 0    !when gets laid off, no offer to accept
                        else if (w==np2) then   !nothing happens
							chs(1,j,i) = wl2q(np1,l0)	     !unemployment at l0 always an option
                            chs(2,j,i) = i    !status quo
                            chs(3,j,i) = 0    !when nothing happens, no offer to accept
                        end if
                    else if (l/=l0) then         !ofloc draw
                        if (w<=np) then         !gets wage offer
							chs(1,j,i) = wl2q(np1,l)	     !moving to ofloc unemp is an option 
							chs(2,j,i) = i    !status quo
                            chs(3,j,i) = j    !accept offer
                        else if (w==np1) then   !gets laid off  BUT NOTE THAT THIS DOES NOT HAPPEN WHEN UNEMPLOYED
							chs(1,j,i) = wl2q(np1,l)	     !moving to ofloc unemp is an option 
                            chs(2,j,i) = 0    !when gets laid off, status quo not feasible
                            chs(3,j,i) = 0    !when gets laid off, no offer to accept
                        else if (w==np2) then   !nothing happens
							chs(1,j,i) = wl2q(np1,l)	     !moving to ofloc unemp is an option 
                            chs(2,j,i) = i    !status quo
                            chs(3,j,i) = 0    !when nothing happens, no offer to accept
                        end if
                    end if 
                end if     
                        
				if (mysay==0.and.writecho) then
					write(89898688,*), "**********************************************"
					write(89898688,'(1x,tr3,"c",tr2,"w",tr3,"l")')
					write(89898688,'(1x,3i4)') c,q2w(i),q2l(i)
					write(89898688,'(1x,3i4)') c,q2w(j),q2l(j)
					qcho=chs(2,j,i)
					if (qcho==0) then 
						write(89898688,'(1x,3i4)') 2,0,0					
					else if ( qcho>0 ) then 
						write(89898688,'(1x,3i4)') 2,q2w(qcho),q2l(qcho)
					end if
					qcho=chs(3,j,i)
					if (qcho==0) then 
						write(89898688,'(1x,3i4)') 3,0,0					
					else if ( qcho>0 ) then 
						write(89898688,'(1x,3i4)') 3,q2w(qcho),q2l(qcho)
					end if
					write(89898688,*), "**********************************************"
				end if 


				!if ( w /= np1 .and. w0 /= np1 ) then	! w=np+1 denotes getting laid off (in the case where q denotes shocks) 
									! and w0=np+1 denotes the state variable w0 indicating that you start that period as unemployed 
				!	chs(c+1,j,i) = i		! c+1 is nl+1 which denotes the status quo option. can choose the status quo only if you don't get laid off. also don't need to allow this choice if the status quo is unemployment 	
				!end if 
				!if ( w <= np ) then			! offer is only the cases where w<=np. in other offer situation where w>np, it means either that you get laid off or nothing happens, so those things do not really give you any options. 
				!	chs(c+2,j,i) = j		! c+2 is nl+2 which denotes the option of taking the offer q
				!end if 		
			end do q
		else 
			chs(:,:,i)=0
		end if 
	end do q0
	!if (skriv) print*, "done with getch_single "
	if (mysay==0.and.writecho) then
		close(89898688)
	end if 
    end subroutine 


	subroutine getch_couple
	! choice vectors that indicate whether an alternative is feasible given the status quo and given the offer situation
	integer(i4b) :: i,j,a,b,cc,qc(2),is(2),js(2),l(2), d(nc),k
	IF (CHOICESETLARGE) THEN
		if (mysay==0.and.writecho) then
			open(unit=89898689, file='choc.txt',status='replace')
		end if 
			ch=0
			qq0: do i=1,nq
					if (  qq2w(1,i)<=np1 .or.  qq2w(2,i)<=np1 )  then
						qq: do j=1,nq			
							cc=0						! counter for number of choices for couple
							is(:)=qq2q(:,i)				! transform qq0 into single's q0
							js(:)=qq2q(:,j)				! transform qq  into single's q

							do a=1,nl				! choice of husband
								do b=1,nl			! choice of wife
									qc(1)=chs(a,js(1),is(1))	! transform husband's choice into q
									qc(2)=chs(b,js(2),is(2))	! transform wife's choice into q
									if ( minval(qc) > 0 ) then	! the option has to be feasible for both of them which is denoted by qc(1) and qc(2) both being larger than 0.
										l(1) = q2l( qc(1) )	! get the location that corresponds to husband choice qc(1) here 
										l(2) = q2l( qc(2) )	! get the location that corresponds to wife whoice qc(2) here 
										if ( l(1) == l(2) ) then !a choice is feasible for the couple only if the respective locations of that choice is equal. the others are not options.
											cc=cc+1	
											if (mysay==0.and.writecho) then
												write(89898689,*), "**********************************************"
												write(89898689,'(1x,"a",tr3,"b",tr3,"cc",tr2,"wh",tr2,"wf",tr2,"lh",tr2,"lf")')
												write(89898689,'(1x,7i4)') a,b,cc,q2w(is(1)),q2w(is(2)),q2l(is(1)),q2l(is(2))
												write(89898689,'(1x,7i4)') a,b,cc,q2w(js(1)),q2w(js(2)),q2l(js(1)),q2l(js(2))
												write(89898689,'(1x,7i4)') a,b,cc,q2w(qc(1)),q2w(qc(2)),q2l( qc(1) ),q2l( qc(2) )
												write(89898689,*), "**********************************************"
											end if 
											ch(cc,j,i)=q2qq( qc(1) , qc(2) )	!transform single's choice qm and qf into couple's choice qq  
										end if 
									end if 
								end do 
							end do
							!if (mysay==0) print*, "Here is cc after uu", cc

							do a=nl+1,max(ncs,nl+1)						! choice of husband
								do b=1,nl				! choice of wife
									qc(1)=chs(a,js(1),is(1))	! transform husband's choice into q
									qc(2)=chs(b,js(2),is(2))	! transform wife's choice into q
									if ( minval(qc) > 0 ) then	! the option has to be feasible for both of them which is denoted by qc(1) and qc(2) both being larger than 0.
										l(1) = q2l( qc(1) )	! get the location that corresponds to husband choice qc(1) here 
										l(2) = q2l( qc(2) )	! get the location that corresponds to wife whoice qc(2) here 
										if ( l(1) == l(2) ) then !a choice is feasible for the couple only if the respective locations of that choice is equal. the others are not options.
											cc=cc+1	
											if (mysay==0.and.writecho) then
												write(89898689,*), "**********************************************"
												write(89898689,'(1x,"a",tr3,"b",tr3,"cc",tr2,"wh",tr2,"wf",tr2,"lh",tr2,"lf")')
												write(89898689,'(1x,7i4)') a,b,cc,q2w(is(1)),q2w(is(2)),q2l(is(1)),q2l(is(2))
												write(89898689,'(1x,7i4)') a,b,cc,q2w(js(1)),q2w(js(2)),q2l(js(1)),q2l(js(2))
												write(89898689,'(1x,7i4)') a,b,cc,q2w(qc(1)),q2w(qc(2)),q2l( qc(1) ),q2l( qc(2) )
												write(89898689,*), "**********************************************"
											end if 
											ch(cc,j,i)=q2qq( qc(1) , qc(2) )	!transform single's choice qm and qf into couple's choice qq  
										end if 
									end if 
								end do !b
							end do !a
							!if (mysay==0) print*, "Here is cc after eu", cc


							do a=1,nl						! choice of husband
								do b=nl+1,max(ncs,nl+1)				! choice of wife
									qc(1)=chs(a,js(1),is(1))	! transform husband's choice into q
									qc(2)=chs(b,js(2),is(2))	! transform wife's choice into q
									if ( minval(qc) > 0 ) then	! the option has to be feasible for both of them which is denoted by qc(1) and qc(2) both being larger than 0.
										l(1) = q2l( qc(1) )	! get the location that corresponds to husband choice qc(1) here 
										l(2) = q2l( qc(2) )	! get the location that corresponds to wife whoice qc(2) here 
										if ( l(1) == l(2) ) then !a choice is feasible for the couple only if the respective locations of that choice is equal. the others are not options.
											cc=cc+1	
											if (mysay==0.and.writecho) then
												write(89898689,*), "**********************************************"
												write(89898689,'(1x,"a",tr3,"b",tr3,"cc",tr2,"wh",tr2,"wf",tr2,"lh",tr2,"lf")')
												write(89898689,'(1x,7i4)') a,b,cc,q2w(is(1)),q2w(is(2)),q2l(is(1)),q2l(is(2))
												write(89898689,'(1x,7i4)') a,b,cc,q2w(js(1)),q2w(js(2)),q2l(js(1)),q2l(js(2))
												write(89898689,'(1x,7i4)') a,b,cc,q2w(qc(1)),q2w(qc(2)),q2l( qc(1) ),q2l( qc(2) )
												write(89898689,*), "**********************************************"
											end if 
											ch(cc,j,i)=q2qq( qc(1) , qc(2) )	!transform single's choice qm and qf into couple's choice qq  
										end if 
									end if 
								end do !b
							end do !a
							!if (mysay==0) print*, "Here is cc after ue", cc

							do a=nl+1,max(ncs,nl+1)							! choice of husband
								do b=nl+1,max(ncs,nl+1)					! choice of wife
									qc(1)=chs(a,js(1),is(1))	! transform husband's choice into q
									qc(2)=chs(b,js(2),is(2))	! transform wife's choice into q
									if ( minval(qc) > 0 ) then	! the option has to be feasible for both of them which is denoted by qc(1) and qc(2) both being larger than 0.
										l(1) = q2l( qc(1) )	! get the location that corresponds to husband choice qc(1) here 
										l(2) = q2l( qc(2) )	! get the location that corresponds to wife whoice qc(2) here 
										if ( l(1) == l(2) ) then !a choice is feasible for the couple only if the respective locations of that choice is equal. the others are not options.
											cc=cc+1	
											if (mysay==0.and.writecho) then
												write(89898689,*), "**********************************************"
												write(89898689,'(1x,"a",tr3,"b",tr3,"cc",tr2,"wh",tr2,"wf",tr2,"lh",tr2,"lf")')
												write(89898689,'(1x,7i4)') a,b,cc,q2w(is(1)),q2w(is(2)),q2l(is(1)),q2l(is(2))
												write(89898689,'(1x,7i4)') a,b,cc,q2w(js(1)),q2w(js(2)),q2l(js(1)),q2l(js(2))
												write(89898689,'(1x,7i4)') a,b,cc,q2w(qc(1)),q2w(qc(2)),q2l( qc(1) ),q2l( qc(2) )
												write(89898689,*), "**********************************************"
											end if 
											ch(cc,j,i)=q2qq( qc(1) , qc(2) )	!transform single's choice qm and qf into couple's choice qq  
										end if 
									end if 
								end do !b
							end do !a
							!if (mysay==0) print*, "Here is cc after qq", cc

						end  do qq
					else 
						ch(:,:,i)=0
					end if 

			end do 	qq0
			if (skriv) then 
				if (cc > nc ) then ; print*, " cc > nc " , cc , nc ; end if		!theoretically, the max number of options they have is nl+8 (nl is the uu for each l option, and then it's 3*3-1 because it's combinatorial of --- --- where each have 3 things they can be (u,q0,q). so the total combo comes to 3*3. but then you subtract 1 because uu is already accounted for in the first nl options
			!	print*, "done with getch_couple "
			end if 

			if (mysay==0.and.writecho) then
				close(89898689)
			end if 
			
	ELSE if (.NOT.CHOICESETLARGE) then 

		!stay both unemp, stay emp,u, stay u,emp, stay u, u     move both unemp, move emp u, move u emp, move uu 
		if (mysay==0.and.writecho) then
			open(unit=89898689, file='choc.txt',status='replace')
		end if 
		ch=0
		do i=1,nq
				if (  qq2w(1,i)<=np1 .or.  qq2w(2,i)<=np1 )  then
					do j=1,nq			
						cc=0						! counter for number of choices for couple
						is(:)=qq2q(:,i)				! transform qq0 into single's q0
						js(:)=qq2q(:,j)				! transform qq  into single's q  
						do a=1,ncs						! choice of husband
							do b=1,ncs				! choice of wife
								qc(1)=chs(a,js(1),is(1))	! transform husband's choice into q
								qc(2)=chs(b,js(2),is(2))	! transform wife's choice into q
								if ( minval(qc) > 0 ) then	! the option has to be feasible for both of them which is denoted by qc(1) and qc(2) both being larger than 0.
									l(1) = q2l( qc(1) )	! get the location that corresponds to husband choice qc(1) here 
									l(2) = q2l( qc(2) )	! get the location that corresponds to wife whoice qc(2) here 
									if ( l(1) == l(2) ) then !a choice is feasible for the couple only if the respective locations of that choice is equal. the others are not options.
										cc=cc+1	
										!if (mysay==0.and.writecho) then
										!	write(89898689,*), "**********************************************"
										!	write(89898689,'(1x,tr3,"a",tr3,"b",tr2,"cc",tr2,"wh",tr2,"wf",tr2,"lh",tr2,"lf")')
										!	write(89898689,'(1x,7i4)') a,b,cc,q2w(is(1)),q2w(is(2)),q2l(is(1)),q2l(is(2))
										!	write(89898689,'(1x,7i4)') a,b,cc,q2w(js(1)),q2w(js(2)),q2l(js(1)),q2l(js(2))
										!	write(89898689,'(1x,7i4)') a,b,cc,q2w(qc(1)),q2w(qc(2)),q2l( qc(1) ),q2l( qc(2) )
										!	write(89898689,*), "**********************************************"
										!end if 
										ch(cc,j,i)=q2qq( qc(1) , qc(2) )	!transform single's choice qm and qf into couple's choice qq  
									end if 
								end if 
							end do !b
						end do !a
						!if (mysay==0) print*, "Here is cc after qq", cc

					end  do !qq
				else 
					ch(:,:,i)=0
				end if 

		end do 	!qq0
		if (skriv) then 
			if (cc > nc ) then ; print*, " cc > nc " , cc , nc ; end if		!theoretically, the max number of options they have is nl+8 (nl is the uu for each l option, and then it's 3*3-1 because it's combinatorial of --- --- where each have 3 things they can be (u,q0,q). so the total combo comes to 3*3. but then you subtract 1 because uu is already accounted for in the first nl options
		!	print*, "done with getch_couple "
		end if 

		if (mysay==0.and.writecho) then
			close(89898689)
		end if 


	END IF !CHOICESETLARGE DETERMINATION
	

	end subroutine getch_couple

	subroutine getgauss
	! gauss quadrature hermite weights and abscissas for wage draws and marriage utility shocks
	integer(i4b) :: g,i,j,k,trueindex,edo
	real(sp) :: abs1(np),wgt1(np)
	real(sp) :: abs2(nz),wgt2(nz)
	real(sp) :: abs3(nepsmove),wgt3(nepsmove)
    real(sp) :: dum1(np),dum2(np)
	real(sp) :: dum3(nz),dum4(nz)
	real(sp) :: dum5(nepsmove),dum6(nepsmove)
    real(dp) :: RHO(2,2),CD(2,2),epsw(2)
	integer(i4b) :: indeces(2),loca !for iepjoint indices
	real(8) :: array5(5),array9(9),array3(3)
	real(8), dimension(np) :: wgtA,wgtB
	real(8), dimension(np,np) :: wgtBgivenA

    !**********************************************************************
    !*Call Cholesky routine that does the following:                      *
    !*Construct the covariance matrix btw male and female wage draws         *
    !*Take Cholesky Decomposition of the covariance matrix                *
    !**********************************************************************
    !Construct the correlation matrix 
    !And then using the corr matrix, construct the cov matrix
    !sigma(1)=sig_wge(1)
    !sigma(2)=sig_wge(2)
    RHO(1,1)=1.0_dp		
    RHO(2,2)=1.0_dp	
    RHO(2,1)=ro		
    RHO(1,2)=RHO(2,1)
    !if (myid.eq.0) then
    !	print*, "Here are the sigmas: "
    !	print*, sig_w1,sig_w2,ro1
    !	print*, "Here is the correlation matrix: "
    !	print*, RHO(1,1),RHO(1,2)
    !	print*, RHO(2,1),RHO(2,2)
    !end if
    !Now turn correlation matrix into varcov matrix
    do j=1,2
	    do k=1,2
		    !RHO(j,k)=RHO(j,k)*sig_wge(j)*sig_wge(k)
	    end do 
    end do 
    if (skriv) then
    	print*, "Here is the covariance matrix: "
    	print*, RHO(1,1),RHO(1,2)
    	print*, RHO(2,1),RHO(2,2)
    end if
    CALL cholesky_m(RHO,2,CD)
    if (skriv) then
    	print*, "Here is the cholesky decomposition: "
    	print*, CD(1,1),CD(1,2)
    	print*, CD(2,1),CD(2,2)
    end if

	call gauher(abs1 , wgt1 )				
	do i=1,np
		dum1(i)=abs1(np-i+1)
		dum2(i)=wgt1(np-i+1)
	end do 
	abs1	=	dum1
	wgt1	=	dum2
	wgt	=	wgt1 / sqrtpi						! weights for wage distributions
	!sex: do g=1,2
		!ahu 0327 wg(:,g) = sqrt(2.0_dp) * sig_wge(g) * abs1(:) + mu_wge(g)	! abscissas for wage draw distribution (g for males and females)
        !wg(:,g) = sqrt(2.0_dp) * abs1(:) 	!ahu 0327  abscissas for wage draw distribution (g for males and females)
    !wg(:,1)=CD(1,1)*sqrt(2.0_dp) * abs1(:)
    !wg(:,2)=CD(2,1)*sqrt(2.0_dp) * abs1(:) + CD(2,2)*sqrt(2.0_dp) * abs1(:)
	!end do sex 
	do edo=1,neduc !sig wge is different by education
	    wg(:,edo,MALES)=sig_wge(edo,MALES)*sqrt(2.0_dp) * abs1(:)
	    wg(:,edo,FEMALES)=sig_wge(edo,FEMALES)*sqrt(2.0_dp) * abs1(:)
	end do 

	wgtA=wgt
	wgtBgivenA=0.0_dp !initiate
	do i=1,np
		do j=1,np
			if (j==i) then 
				wgtBgivenA(j,i)=corral
			end if 
		end do 
	end do 

	do i=1,np
		do j=1,np
			if (j/=i) then 
				wgtBgivenA(j,i)= (1.0_dp - wgtBgivenA(i,i)) / (np-1)
			end if 
		end do 
	end do 

	do i=1,np 
		do j=1,np
			wgtjoint(i,j)= wgtBgivenA(j,i) * wgtA(i)
		end do 
	end do

	do j=1,np
		wgtB(j) = sum(wgtjoint(:,j))
	end do 

	wgtmen=wgtA
	wgtfem=wgtB

	if (abs( sum(wgtmen) - 1.0_dp ) > 1.0d-6  ) then 
		print*, "sum of wgtmen is not 1! " ,sum(wgtmen),wgtmen
		stop 
	end if 
	if (abs( sum(wgtfem) - 1.0_dp ) > 1.0d-6  ) then 
		print*, "sum of wgtfem is not 1! " ,sum(wgtfem),wgtfem
		stop 
	end if 
	if (abs( sum(wgtjoint) - 1.0_dp ) > 1.0d-4  ) then 
		print*, "sum of wgtjoint is not 1! " ,sum(wgtjoint),wgtjoint
		stop 
	end if 

	if (mysay==0.and.skriv) then 
		print*, "WAGE PROBS FROM DISCRETIZATION "
		do i=1,np
			do j=1,np
		print*, i,j,wgtmen(i),wgtfem(j),wgtjoint(i,j)
			end do 
		end do 
	end if 

	if (METHODWAGEDISCRETIZE==1) then
		array3= (/ -0.99_dp,  0.0_dp,  0.99_dp /)
		do edo=1,neduc
			!wg(:,edo,MALES) = array3 * sig_wge(edo,MALES)
			!wg(:,edo,FEMALES) = array3 * sig_wge(edo,FEMALES)
		end do 
		wgtmen=1.0_dp/np
		wgtfem=1.0_dp/np
		if (corral==0.2) then
			wgtjoint(1,1)=.37*wgtmen(1)
			wgtjoint(1,2)=.34*wgtmen(1)
			wgtjoint(1,3)=.29*wgtmen(1)
			wgtjoint(2,1)=.33*wgtmen(2) 
			wgtjoint(2,2)=.34*wgtmen(2)
			wgtjoint(2,3)=.33*wgtmen(2)
			wgtjoint(3,1)=.30*wgtmen(3)
			wgtjoint(3,2)=.33*wgtmen(3)
			wgtjoint(3,3)=.37*wgtmen(3)
		else if (corral==0.5) then
			wgtjoint(1,1)=.55*wgtmen(1)
			wgtjoint(1,2)=.32*wgtmen(1)
			wgtjoint(1,3)=.13*wgtmen(1)
			wgtjoint(2,1)=.31 *wgtmen(2)
			wgtjoint(2,2)=.37*wgtmen(2)
			wgtjoint(2,3)=.32*wgtmen(2)
			wgtjoint(3,1)=.14*wgtmen(3) 
			wgtjoint(3,2)=.31*wgtmen(3)
			wgtjoint(3,3)=.55*wgtmen(3)
		else if (corral==0.8) then
			wgtjoint(1,1)=.72*wgtmen(1)
			wgtjoint(1,2)=.25*wgtmen(1)
			wgtjoint(1,3)=.03*wgtmen(1)
			wgtjoint(2,1)=.25 *wgtmen(2)
			wgtjoint(2,2)=.50*wgtmen(2)
			wgtjoint(2,3)=.25*wgtmen(2)
			wgtjoint(3,1)=.03*wgtmen(3)
			wgtjoint(3,2)=.25*wgtmen(3)
			wgtjoint(3,3)=.72*wgtmen(3)
		else if (corral==0.9) then
			wgtjoint(1,1)=.80*wgtmen(1)
			wgtjoint(1,2)=.19*wgtmen(1)
			wgtjoint(1,3)=.01*wgtmen(1)
			wgtjoint(2,1)=.19*wgtmen(2) 
			wgtjoint(2,2)=.62*wgtmen(2)
			wgtjoint(2,3)=.19*wgtmen(2)
			wgtjoint(3,1)=.01*wgtmen(3) 
			wgtjoint(3,2)=.19*wgtmen(3)
			wgtjoint(3,3)=.80*wgtmen(3)
		else if (corral==1) then
			wgtjoint(1,1)=1.0*wgtmen(1)
			wgtjoint(1,2)=0.*wgtmen(1)
			wgtjoint(1,3)=0.*wgtmen(1)
			wgtjoint(2,1)=0. *wgtmen(2)
			wgtjoint(2,2)=1.0 *wgtmen(2)
			wgtjoint(2,3)=0.*wgtmen(2)
			wgtjoint(3,1)=0.*wgtmen(3)
			wgtjoint(3,2)=0.*wgtmen(3)
			wgtjoint(3,3)=1.0*wgtmen(3)
		end if 
		!print*, "SHOULD NOT BE HERE "
		!STOP

		if (abs( sum(wgtmen) - 1.0_dp ) > 1.0d-6  ) then 
			print*, "sum of wgtmen is not 1! " ,sum(wgtmen),wgtmen
			stop 
		end if 
		if (abs( sum(wgtfem) - 1.0_dp ) > 1.0d-6  ) then 
			print*, "sum of wgtfem is not 1! " ,sum(wgtfem),wgtfem
			stop 
		end if 
		if (abs( sum(wgtjoint) - 1.0_dp ) > 1.0d-3  ) then 
			print*, "sum of wgtjoint is not 1! " ,sum(wgtjoint),wgtjoint
			stop 
		end if 
	
		if (mysay==0.and.skriv) then 
			print*, "WAGE PROBS FROM DISCRETIZATION "
			do i=1,np
				do j=1,np
			print*, i,j,wgtmen(i),wgtfem(j),wgtjoint(i,j)
				end do 
			end do 
		end if 	
	end if 

	
	!wg(1,1)=-0.54_dp
	!wg(2,1)=-0.33_dp
	!wg(3,1)=-0.20_dp
	!wg(4,1)=-0.11_dp
	!wg(5,1)=0.0_dp
	!wg(6,1)=0.11_dp
	!wg(7,1)=0.20_dp
	!wg(8,1)=0.33_dp
	!wg(9,1)=0.54_dp

	!wg(1,2)=-0.77_dp
	!wg(2,2)=-0.47_dp
	!wg(3,2)=-0.29_dp
	!wg(4,2)=-0.15_dp
	!wg(5,2)=0.0_dp
	!wg(6,2)=0.14_dp
	!wg(7,2)=0.29_dp
	!wg(8,2)=0.47_dp
	!wg(9,2)=0.78_dp

	!wg(1,1)=-0.44_dp
	!wg(2,1)=-0.18_dp
	!wg(3,1)=0.0_dp
	!wg(4,1)=0.18_dp
	!wg(5,1)=0.44_dp

	!wg(1,2)=-0.54_dp
	!wg(2,2)=-0.22_dp
	!wg(3,2)=0.0_dp
	!wg(4,2)=0.22_dp
	!wg(5,2)=0.54_dp

	!wg(1,1)=  (-1.3)*sig_wge(1)
	!wg(2,1)=  (-1.0)*sig_wge(1)
    !wg(3,1)=(-2.0_dp/3.0)*sig_wge(1)
    !wg(4,1)=(-1.0_dp/3.0)*sig_wge(1)
    !wg(5,1)=0.0_dp
	!wg(6,1)=(1.0_dp/3.0)*sig_wge(1)
    !wg(7,1)=(2.0_dp/3.0)*sig_wge(1)
    !wg(8,1)=(1.0)*sig_wge(1)
    !wg(9,1)=(1.3)*sig_wge(1)


	!wg(1,2)=  (-1.3)*sig_wge(2)
	!wg(2,2)=  (-1.0)*sig_wge(2)
    !wg(3,2)=(-2.0_dp/3.0)*sig_wge(2)
    !wg(4,2)=(-1.0_dp/3.0)*sig_wge(2)
    !wg(5,2)=0.0_dp
	!wg(6,2)=(1.0_dp/3.0)*sig_wge(2)
    !wg(7,2)=(2.0_dp/3.0)*sig_wge(2)
    !wg(8,2)=(1.0)*sig_wge(2)
    !wg(9,2)=(1.3)*sig_wge(2)
	!wgt=1.0_dp/np

	if (np==5 ) then 
		!*********************************
		!for np=5
		!got the inver normal cdf values manually from: https://burymathstutor.co.uk/inverseNormal.html
		!Follwing the method of Kennan-Discrete-Approx and its notation
		!support points x_i have to be such that F(x_i)=(2i-1)/2n
		!So for example for n=5: 
		!i=1 -----> F(x_1)=1/10
		!i=2 -----> F(x_2)=3/10
		!i=3 -----> F(x_3)=5/10
		!i=4 -----> F(x_4)=7/10
		!i=5 -----> F(x_5)=9/10
		!where all the support points have equal weights
		!for standard normal, some quantiles are the following. To get for sigma, just multiple by sigma.
		!-1.28    P=0.1   stdev =1 
		!-0.52    P=0.3   stdev =1 
		!0     P = 0.5  stdev =1
		!0.5244  P=0.7   stdev=1
		!1.28     P=0.9
		!array5 = (/ -1.28_dp, -0.52_dp , 0.0_dp, 0.52_dp , 1.28_dp /)
		!do edo=1,neduc
	!		wg(:,edo) = array5 * sig_wge(edo)
!		end do 
	else if (np==9) then 
		!*********************************
		!for np=9
		!got the inverse normal cdf values manually from: https://burymathstutor.co.uk/inverseNormal.html
		!	-1.5937      1/18 = 0.0555  this is F(x_1) for np=9
		!	-0.96745    2*2-1 / 18 = 1/6 = 0.16666
		!  -0.58948    2*3-1 / 18 = 5/18 = 0.27777
		!	-0.28224   2*4-1 / 18 = 7/18 = 0.38888
		!	2*5-1 / 18 = 1/2
		!0.28221   2*6-1 / 18 = 11/18 = 0.611111
		!0.58945     2*7-1 / 18 = 13/18 = 0.722222
		! 0.96741   2*8-1 / 18 = 15/18 = 0.833333
		!1.5932     2*9-1 / 18 = 17/18 = 0.94444
		!array9= (/ -1.5937_dp, -0.96745_dp , -0.58948_dp, -0.28224_dp, 0.0_dp , 0.28221_dp, 0.58945_dp, 0.96741_dp, 1.5932_dp /)
		!do edo=1,neduc
			!wg(:,edo) = array9 * sig_wge(edo)
		!end do 
	else if (np==3) then 
		!array3= (/ -0.99_dp,  0.0_dp,  0.99_dp /)
		!do edo=1,neduc
		!	wg(:,edo,MALES) = array3 * sig_wge(edo,MALES)
		!	wg(:,edo,FEMALES) = array3 * sig_wge(edo,FEMALES)
		!end do 
		!wgt=1.0_dp/np
	end if
		
		!if (np==5) then 
		!	do edo=1,2
		!		wg(1,edo)=  (-1.28)*sig_wge(edo)
		!		wg(2,edo)=  (-0.52)*sig_wge(edo)
		!		wg(3,edo)= 0.0_dp
		!		wg(4,edo)= 0.52*sig_wge(edo)
		!		wg(5,edo)= 1.28*sig_wge(edo)
		!	end do 
		!else if (np==9) then 
		!end if
	!	wgt=1.0_dp/np !weights of the support points are equal. as in kennan-discrete-approx.

	!do edo=1,2
	!	do i=1,np
			!cumprob=(2.0*i-1.)/(2.0*np)
			!if (mysay==0) print*, "Here is cumprob", i,cumprob
			!call inverse_gaussian_cdf(cumprob,meanpdf,sigmapdf, xstar_i )
			!call cdfnorminv( 1, cumprob, xstar_i )
			!if (mysay==0) print*, "Here is cumprob", i,cumprob,xstar_i
	!		wg(i,edo)=xstar_i
	!	end do
	!end do 
	!wgt=1.0_dp/np !weights of the support points are equal. as in kennan-discrete-approx.

	call gauher( abs2 , wgt2)				
	do i=1,nz
		dum3(i)=abs2(nz-i+1)
		dum4(i)=wgt2(nz-i+1)
    end do 
    abs2=dum3
	wgt2=dum4
	ppmarie(:) = wgt2 / sqrtpi							! weights for mar distribution
	!call rand2num(mu_mar,sig_mar,abs2,mg)
    do i=1,nz
	    marshock(i) = sqrt(2.0_dp) * sig_mar * abs2(i)			! abscissas for marrige utility shock distribution
	    !mg(:,trueindex) = sqrt(2.0_dp) * sig_mar * dble(abs2(:)) + mu_mar(j)		! abscissas for marrige utility shock distribution
    end do 
	!mumar by type (Doesn't really need to be calculated here)
    do trueindex=1,ninp
        call index2cotyphome(trueindex,i,j,k)
        mg(trueindex) = mu_mar(j)	
    end do 
    
	call gauher( abs3 , wgt3)				
	do i=1,nepsmove
		dum5(i)=abs3(nepsmove-i+1)
		dum6(i)=wgt3(nepsmove-i+1)
    end do 
	!if (skriv.and.(.not.chkstep).and.(.not.optimize)) print*, "abs2 ", abs2
    !if (skriv.and.(.not.chkstep).and.(.not.optimize)) print*, "dum3 ", dum3
    abs3=dum5
	wgt3=dum6
	ppmovesingle = wgt3 / sqrtpi							! weights for moveshock distribution
	!call rand2num(mu_mar,sig_mar,abs2,mg)
	moveshockmar(:,1) = sqrt(2.0_dp) * sigo_m * abs3(:) + mu_o			! abscissas for marrige utility shock distribution
	moveshockmar(:,2) = sqrt(2.0_dp) * sigo_f * abs3(:) + mu_o			! abscissas for marrige utility shock distribution

	moveshocksin(:,1)=moveshockmar(:,1)      != sqrt(2.0_dp) * sigo_sin * abs3(:) + mu_o			! abscissas for marrige utility shock distribution
	moveshocksin(:,2)=moveshockmar(:,2)
	!if ((.not.chkstep).and.(.not.optimize) ) print*, "mg :", mg(:)
	!moveshock(:) = sqrt(2.0_dp) * sig_o * dble(abs3(:)) + mu_o			! abscissas for marrige utility shock distribution
    if (nepsmove==1) then 
        moveshockmar=0.0_dp
        moveshocksin=0.0_dp
    end if
    do i=1,nepsmove
		do j=1,nepsmove
			moveshockjoint(i,j)%hub=moveshockmar(i,1)
			moveshockjoint(i,j)%wfe=moveshockmar(j,2)
			k=ndim2lin( (/ nepsmove , nepsmove /),(/ i,j /) )
			ppmovejoint(k)=ppmovesingle(i)*ppmovesingle(j)
			indeces=lin2ndim( (/ nepsmove , nepsmove /) , k )
			if (i.ne.indeces(1)) then ; print*, "sth wrong with indeces ",i,j,k ; stop ; end if 
			if (j.ne.indeces(2)) then ; print*, "sth wrong with indeces ",i,j,k ; stop ; end if 
			if (JOINTMARSHOCK) then 
				if (i==j) then
					ppmovejoint(k)=ppmovesingle(i)
				else 
					ppmovejoint(k)=0.0_dp
				end if 
			end if 
		end do 
	end do
    if (skriv) print*, 'here is ppmovesingle',ppmovesingle,sum(ppmovesingle)    
    if (skriv) print*, 'here is ppmovejoint',ppmovejoint,sum(ppmovejoint)    
	if (abs( sum(ppmovejoint) - 1.0_dp ) > 1.0d-6  ) then 
		print*, "sum of ppmovejoint is not 1! " ,sum(ppmovejoint),ppmovejoint
		stop 
	end if 
    end subroutine getgauss
	
	subroutine getppmeet
	integer(i4b) :: i,j,w(2),age,notg,g,typo
	real(dp), dimension(ntypp,MALES:FEMALES) :: ptypeHS_cond_type,ptypeCOL_cond_type
	real(dp) :: edmeet,kidmeet,expmeet,ptype_uncond,temp(nxs),checkprob(1:2,NXS,ntypp,MALES:FEMALES)

	ppmeetq=0.0_dp
	do j=1,nqs	
	!	w(1) = q2w( i ) 
		do i=1,nqs		
			if ( q2w(j) <= np1 .and. q2w(i) <= np1 .and. q2qq(i,j) > 0 ) then 
					ppmeetq(i,j)=exp( meetsame(1)*one(q2w(j)==q2w(i) ) )  
			end if 
		end do 
		ppmeetq(:,j)=ppmeetq(:,j)/sum(ppmeetq(:,j))
	end do 
	!ppmeetx=0.0_dp
	!do j=1,nxs	
	!	do i=1,nxs		
	!		if ( x2xx(i,j) > 0 ) then
	!			ppmeetx(i,j)=exp(  meetsame(1)*one(x2r(i).eq.x2r(j)) + meetsame(2)*one(x2e(i).eq.x2e(j))   )
	!		end if
	!		if (x2kid(i).ne.x2kid(j)) then !can only meet someone with the same number of kid as you 
	!			ppmeetx(i,j)=0.0			
	!		end if 
	!	end do 
	!	ppmeetx(:,j)=ppmeetx(:,j)/sum(ppmeetx(:,j))
	!end do



	!FIRST get the probability of COLLEGE or HS conditional on TYPE using Bayes' Rule
	do g=MALES,FEMALES
		do typo=1,ntypp
		ptype_uncond=ptypehs(typo,g)*0.5_dp + ptypecol(typo,g)*0.5_dp
		ptypeHS_cond_type(typo,g)=  (  ptypeHS(typo,g)*0.5_dp   ) / ptype_uncond 
		ptypeCOL_cond_type(typo,g)=  ( ptypeCOL(typo,g)*0.5_dp  ) / ptype_uncond 
		end do
	end do
	
	!Next, use those probs to calculated pmeetx
	probmeetx=0.0_dp
	do g=MALES,FEMALES
		do typo=1,ntypp
	
			do j=1,nxs		!xs j
				temp=0.0_dp ; edmeet=0.0_dp ; kidmeet=0.0_dp ; expmeet=0.0_dp
				do i=1,nxs	!xs i	
					if ( x2xx(i,j) > 0 ) then				
						if (g==1) notg=2
						if (g==2) notg=1
						if ( x2e(i).eq.1 ) then 
							edmeet=ptypeHS_cond_type(TYPO,NOTG)
						else if (x2e(i).eq.2) then
							edmeet=ptypeCOL_cond_type(TYPO,NOTG)
						end if
						if (x2kid(i).ne.x2kid(j)) then !can only meet someone with the same number of kid as you 
							kidmeet=0.0_dp
						else 
							kidmeet=1.0_dp			
						end if
						
						expmeet=1.0_dp
			
						temp(i)=kidmeet*expmeet*edmeet	
					end if
				end do !xs i  
				probmeetx(:,j,TYPO,g)=temp(:)/sum(temp(:))
				probmeetx(:,j,TYPO,g)=temp(:)/sum(temp(:))
			end do 		!xs j

			!checkprob=0.0_dp
			do j=1,nxs
				if ( abs( sum(probmeetx(:,j,TYPO,g))-1.0_dp) > eps ) then ; print*, "probxmeet doesn't add up ", sum(probmeetx(:,j,TYPO,g)) ; stop ; end if 	
				!do i=1,nxs
				!	checkprob(x2e(i),j,TYPO,g)=checkprob(x2e(i),j,TYPO,g)+probmeetx(i,j,TYPO,g)
				!end do 
			end do 
	
		end do 	!typo 
	end do 		!g MALES AND FEMALES





	IF (mysay==0 .AND. ITER==1) THEN
		open(unit=94639, file='checkprob.txt',status='replace')
		DO TYPO=1,NTYPP
			write(94639,*) "MALES / TYPO= ", TYPO
			DO j=1,nxs
				do i=1,nxs
				write(94639,'(tr3,"j",tr3,"i",tr1,"edj",tr1,"edi",tr3,"prob",tr3,tr2,"xi",tr2,"xj",tr2,"kidj",tr2,"kidi")') 
				write(94639,'(4I4,F10.2,2I4,2I6)' ) j,i,x2e(j),x2e(i),probmeetx(i,j,TYPO,MALES),x2r(j),x2r(i),x2kid(j),x2kid(i)
				end do
			END DO 
		END DO
		write(94639,*) 
		write(94639,*) 
		write(94639,*) 
		write(94639,*) 


			DO TYPO=1,NTYPP
				write(94639,'(I4,2F10.2)' ) TYPO,ptypeHS_cond_type(TYPO,MALES),ptypeCOL_cond_type(TYPO,MALES)
			END DO 
			write(94639,*) 
			write(94639,*) 
			DO TYPO=1,NTYPP
				write(94639,'(I4,2F10.2)' ) TYPO,ptypeHS_cond_type(TYPO,FEMALES),ptypeCOL_cond_type(TYPO,FEMALES)
			END DO 
	
				do j=1,nxs
					write(94639,*) 
					write(94639,*) "MALES. THIS IS ANOTHER J ", J
					DO TYPO=1,NTYPP
						write(94639,*) "MALES / TYPO=", TYPO
						write(94639,'(2F10.2)' ) checkprob(1:2,J,TYPO,MALES)   
					END DO
				end do 
				write(94639,*) 
				do j=1,nxs
					write(94639,*) 
					write(94639,*) "FEMALES. THIS IS ANOTHER J ", J
					DO TYPO=1,NTYPP
						write(94639,*) "FEMALES/ TYPO=", TYPO
						write(94639,'(2F10.2)' ) checkprob(1:2,J,TYPO,FEMALES)   
					END DO 
				end do
				!write(94639,'(tr5,"TYPO",tr4,"HIM",tr4,"HER",tr4,"PROB")') 
		CLOSE(94639)
	END IF 






	pmeet=0.0_dp 
	do age=mna,mxa
		pmeet(age)=logit(alfmeet(1)-alfmeet(2)*age)
	end do








	!		w(2) = q2w( j ) 
	!		if ( w(1) <= np1 .and. w(2) <= np1 ) then 
	!			if ( q2qq(i,j) > 0 ) then 
	!				    if (w(1) <= np .and. w(2) <= np ) then
	!					    udif = aw( w(1) , 1 ) - aw( w(2) , 2 )
	!				    else 
	!					    udif = 0.0_sp
	!				    end if 
	!				    hdif=abs ( ( w(1) <= np .and. w(2) == np+1 ) .or.  ( w(1) == np+1 .and. w(2) <= np )  )
	!				    ppmeetq(1,j,i)=max(0.0_dp,1.0_dp-omega(1)*udif**2-omega(2)*hdif )
	!				    ppmeetq(2,i,j)=max(0.0_dp,1.0_dp-omega(1)*udif**2-omega(2)*hdif )
	!			end if 
	!		 else 
	!			ppmeetq(:,j,i)=0.0_dp
	!		end if         
	!	end do qf
	! end do qm
	!do i=1,nqs
	! if ( q2w(i)>=np1 ) then
	! ppmeetq(1,:,i)=ppmeetq(1,:,i)/sum(ppmeetq(1,:,i))
	! ppmeetq(2,:,i)=ppmeetq(2,:,i)/sum(ppmeetq(2,:,i))
	!end if 
	! if (icheck==1) then 
	!	if ( q2w(i)>=np1 .and. abs( sum(ppmeetq(1,:,i))-1d0)>eps ) then ; print*, "ppmeetq for m doesn't add up ", sum(ppmeetq(1,:,i)) ; stop ; end if 
	!	if ( q2w(i)>=np1 .and. abs( sum(ppmeetq(2,:,i))-1d0)>eps ) then ; print*, "ppmeetq for f doesn't add up ", sum(ppmeetq(2,:,i)) ; stop ; end if 
	!	if ( q2w(i)>np1 .and. sum(ppmeetq(:,:,i)) /= 0.0_dp ) then ; print*, " w>np1 and ppmeetq is not 0 ", sum(ppmeetq(:,:,i)) ; stop ; end if 
	! end if 
	! end do 	
	!ppmeetx=0.0_dp
	!xm: do i=1,nxs			
	!	xf: do j=1,nxs		
	!		if ( x2xx(i,j) > 0 ) then 
	!			ppmeetx(i,j)=exp( 0.0_dp )
	!		end if 
	!	end do xf
	!end do xm
	!do i=1,nxs 
	!	ppmeetx(:,i)=ppmeetx(:,i)/sum(ppmeetx(:,i))
	!	if (icheck==1) then 		
	!		if ( abs( sum(ppmeetx(:,i))-1.0_dp) > eps ) then ; print*, "ppmeetx doesn't add up ", sum(ppmeetx(:,i)) ; stop ; end if 
	!	end if 
	!end do 		
	!if (skriv) print*, "done with ppmeet "
	end subroutine getppmeet

	subroutine checkgauss(npoint,lb,ub,cdf,expval)
	integer(i4b), intent(in) :: npoint
	real(sp) :: x0(npoint),w0(npoint)
	real(dp) :: x(npoint),w(npoint),dum
	real(dp), intent(in) :: lb,ub
	real(dp), intent(out) :: cdf(2),expval
	real(dp) :: sigma=1.0_dp,mu=0.0_dp
	integer(i4b) :: i
	dum=(6.28318_dp)**(-0.5_dp)*exp(-(0.6_dp)**2/2.0_dp)  / 0.72_dp			! dum=sqrt(2.0*pi)*exp(-(0.6)**2.0/2.0)  / 0.72
	call gauher( x0, w0)
	w	= w0 / sqrtpi		
	x	= sqrt(2.0_dp) * sigma * x0 + mu
	cdf=0.0_dp
	expval=0.0_dp
	do i=1,np
		if ( x(i) < lb ) then ; cdf(1) = cdf(1) + w(i) ; end if 
		if ( x(i) < ub ) then ; cdf(2) = cdf(2) + w(i) ; end if 
		if ( x(i)>lb .and. x(i)<ub ) then 
			expval = expval + w(i) * x(i)
		end if 
	end do
	expval=expval / ( cdf(2)-cdf(1) )  					!  / ( sqrtpi * (cdf(2)-cdf(1)) ) no need to divide by sqrtpi anymore since already do that with weights so it's already incorporated in both cdf and expval
	!write(*,'(/1x,a,f12.6,a,f12.6)') 'check value:',expval,'  should be:', dum
	end subroutine checkgauss

SUBROUTINE cholesky_m (a,n,p)
IMPLICIT NONE
      integer,intent(in)::n
      real(8),dimension(n,n),intent(in)::a
      real(8),dimension(n,n),intent(out)::p
      real(8),dimension(n,n)::aa
	real(8) sum
      integer i,j,k
	sum = 0.0d-0
      aa(1:n,1:n)=a(1:n,1:n)
      do 13 i = 1,n
         do 12 j = i,n
         sum=aa(i,j)
	   do 11 k = i-1,1,-1
	      sum = sum - aa(i,k)*aa(j,k)
11	   continue 
         if (i.eq.j) then
           if (sum.le.0.0d-0) then 
            print*, 'choldc failed'
            stop
           end if
           p(i,i) = dsqrt(sum)
           else
             aa(j,i) = sum/p(i,i)
	       p(j,i) = aa(j,i)
         end if
12       continue
13    continue
      return
END SUBROUTINE cholesky_m


	
end module share

