

        subroutine checkdecmar(vec,vdif,transfers,intsol)
            real(dp), intent(in) :: vec(5),vdif(2),transfers(2)
            logical :: intsol(3)
            !ahumarch1522 Now check if the interior condition is really indeed satisfied but get rid of this later
              if (vec(5)-abs(vdif(1)-vdif(2)) < 0.0_dp ) then !ahumarch1522 if the transfers statement is correct, then this should be correct too!
              !   print*, "There is something wrong!", ia,q0,q,x,vec(5),abs(vdif(1)-vdif(2))  !ahumarch1522
                  print*, "intsol1",intsol(1)
                  print*, "vec(5)",vec(5)
                  print*, "Vdif(1)",Vdif(1)
                  print*, "Vdif(2)",Vdif(2)
                  print*, "Vdif1-Vdif2",Vdif(1)-Vdif(2)
                  print*, "vec5 >= abs(Vdif1-Vdif2)?",vec(5),abs(Vdif(1)-Vdif(2))
                  print*, "minval(transfers) >= 0?",minval(transfers),0.0_dp                         
                  print*, "maxval(transfers) <= vec5?",maxval(transfers),vec(5)                     
                  stop                                                    !ahumarch1522
              end if                                                      !ahumarch1522 
              if ( minval(transfers) + eps2 < 0.0_dp .or. maxval(transfers) > vec(5)+eps2 ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
            !      print*, "Int cond fine but transfers not fine!", ia,q0,q,x,z,transfers(1:2),vec(5)  !ahumarch1522
                  print*, "intsol1",intsol(1)
                  print*, "vec(5)",vec(5)
                  print*, "Vdif(1)",Vdif(1)
                  print*, "Vdif(2)",Vdif(2)
                  print*, "Vdif1-Vdif2",Vdif(1)-Vdif(2)
                  print*, "vec5 >= abs(Vdif1-Vdif2)?",vec(5),abs(Vdif(1)-Vdif(2))
                   print*, "minval(transfers) >= 0?",minval(transfers),0.0_dp                         
                   print*, "maxval(transfers) <= vec5?",maxval(transfers),vec(5)                     
                   stop                                                                    !ahumarch1522                            
               end if 
               if ( minval(transfers) < 0.0_dp .or. maxval(transfers) > vec(5) ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
            !       print*, "Int cond fine but transfers not fine 2!", ia,q0,q,x,z,transfers(1:2),vec(5)  !ahumarch1522
                   print*, "intsol1",intsol(1)
                   print*, "vec(5)",vec(5)
                   print*, "Vdif(1)",Vdif(1)
                   print*, "Vdif(2)",Vdif(2)
                   print*, "Vdif1-Vdif2",Vdif(1)-Vdif(2)
                   print*, "vec5 >= abs(Vdif1-Vdif2)?",vec(5),abs(Vdif(1)-Vdif(2))
                   print*, "minval(transfers) >= 0?",minval(transfers),0.0_dp                         
                   print*, "maxval(transfers) <= vec5?",maxval(transfers),vec(5)   
                   !!transfers(2) = (1.0_dp-alf) * (vec(5) + vec(3) - vec(1) ) - alf*( vec(4)-vec(2) )   
                   !print*, "1-alf",1.0_dp-alf
                   !print*, "(vec5+vdif(1))",  (vec(5) + vec(3) - vec(1) )
                   !print*, "(1-alf)(vec5+vdif(1))",  (1.0_dp-alf) * (vec(5) + vec(3) - vec(1) )
                   !print*, "alf(vdif(2))", alf*( vec(4)-vec(2) )   
                   !print*, "transfers2", (1.0_dp-alf) * (vec(5) + vec(3) - vec(1) ) - alf*( vec(4)-vec(2) )   
                   !print*, "transfers2 alt",  (vec(5) + vec(3) - vec(1) ) - alf*( vec(5) + vec(3) - vec(1) + vec(4)-vec(2)  )   
                   !print*, "test1",  ( vec(5) + vec(3) - vec(1) + vec(4)-vec(2)  )   
                   !print*, "test2",  alf*( vec(5) + vec(3) - vec(1) + vec(4)-vec(2)  )   
                    !write(*,'(4F30.20)'),  (vec(5) + vec(3) - vec(1) ),vec(5),vec(3),vec(1)
                    ! testing(1)=vec(5)+vdif(1)
                    !write(*,'(F30.20)'), testing(1)
                    !testing(2)=alf*surplus
                    !write(*,'(F30.20)'), testing(2)
                    ! testing(3)=testing(1)-testing(2)
                    !write(*,'(F30.20)'), testing(3)                             
                !intsol1 T
                ! vec(5)   5327.9026198192287
                ! Vdif(1)   8749.1068802403315
                ! Vdif(2)   3421.2042604211028
                ! Vdif1-Vdif2   5327.9026198192287
                ! vec5 >= abs(Vdif1-Vdif2)?   5327.9026198192287        5327.9026198192287
                ! minval(transfers) >= 0?   0.0000000000000000        0.0000000000000000
                ! maxval(transfers) <= vec5?   5327.9026198192296        5327.9026198192287
                ! 1-alf  0.50000000000000000
                ! (vec5+vdif(1))   14077.009500059561
                ! (1-alf)(vec5+vdif(1))   7038.5047500297806
                ! alf(vdif(2))   1710.6021302105514
                ! transfers2   5327.9026198192296
                ! transfers2 alt   5327.9026198192296
                ! test1   17498.213760480663
                ! test2   8749.1068802403315
                !"out.txt" 133L, 10998B                        
                   stop                                                                    !ahumarch1522                            
            end if 

            !**********************************************************************************
            !THE BELOW ARE IF STATEMENTS FOR CHECKING AT THE BEGINNING OF SOL ROUTINE (AT THE THE BEGINNING OF COUPLES DO LOOP)
            !if (skriv.and.ia==48.and.q0==4.and.q==92.and.x==19.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if !ahu030822
            !if (skriv.and.ia==18.and.q0==4.and.q==92.and.x==19.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if !ahu030822              
            !if (skriv.and.(ia==48).and.(q0==4).and.q<=92.and.x==19.and.iepsmove==10.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if  !ahu030622 !ahu030822
            !if (skriv.and.(ia==48).and.(q0==4).and.(q>=190.and.q<=201).and.x==1.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if  !ahu030622 !ahu030822
                
            !if (skriv.and.(ia<=18.or.ia==45.or.ia==49.or.ia==50).and.(q0<=15).and.(q>=190.and.q<=201).and.x==1.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if  !ahu030622 !ahu030822
                !if (skriv.and.(ia==18.or.ia==49.or.ia==50).and.x==1.and.trueindex==1.and.q0<=50) then ; yazmax=.true. ; else ; yazmax=.false. ; end if
            
            !if (skriv.and.ia==48.and.q0==4.and.q==92.and.iepsmove==10.and.x==19.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if !ahu030822            
            !if (skriv.and.(ia==50).and.(q0<=15).and.(q>=190.and.q<=201).and.x==1.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if  !ahumarch1122

            !ahumarch1122 if ( skriv.and.q0==10.and. q2w( qq2q(1,q) )==np2 .and. q2w( qq2q(2,q) )==np2 .and. q2l( qq2q(1,q) )==1 .and. q2l( qq2q(2,q) )==1 ) then  !ahu030622 
                !ahumarch1122 print*, 'Is this it', q,q2w( qq2q(1,q) ),q2w( qq2q(2,q) ), q2l( qq2q(1,q) ) !ahu030622 
            !ahumarch1122 end if  !ahu030622 
            !if (skriv.and.ia==48.and.trueindex==1.and.x==19.and.q0==4.and.(q==92) ) then ; yazmax=.false. ; else ; yazmax=.false. ; end if
            !if (skriv.and.(ia==18.or.ia==50).and.(q>33.and.q<=37).and.(q0>=32.and.q0<=35).and.x==2) then ; yaz=.true. ; else ; yaz=.false. ; end if  
            !if (skriv.and.(q0>=32.and.q0<=35).and.x==2.and.(q>=33.and.q<=37).and.(ia==mxa.or.ia==47)) then ; yaz=.false. ; else ; yaz=.false. ; end if  
            !if (skriv.and.(ia>=40).and.(q0==18).and.q==24.and.x==1.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if  !ahu 0327 trueindex==2    
           !**********************************************************************************


        end subroutine 	


        !!!!!if (relmax==1) then
            !ahu032122 ahumarch2122 commenting out to save time haveenoughtomakeindiff_alt=(  vec(5) - asum  >= 0.0_dp  )
            !ii is the max of surplusj. 
                !check if that ii option is welldef and interior/corner
                !welldef ----> you have found your best ii and set jmax equal to that
                !not welldef ----> set surplusj(ii)=pen and update nn and go on to the next best ii 
                !    
            !In orer for NB to be well defined, we need: 1) surplus>=0 2) there exists some transfer such that transfer1>=0 and transfer2>=0
            !Condition 2 follows from the requirement that the utility transfer has to only come from current wsum and not the V's. In other words, 
            !w has to be such that it can cover any utility transfers that are needed o have the PC's hold. 
            !See notes in black notebook about why Condition 2 translates into that asum condition. 
            !Once we establish that NB is well defined, then we can go ahead and see if the optimum is interior or corner. 
            !intsol(1)=( vec(5) + eps >= abs(vdif(1)-vdif(2)) )  !interior
            !intsol(2)=( vec(5) <= vdif(1)-vdif(2) + eps  )      !corner where c1=0 and c2=wsum !ahumarch2122 ahu032122 moving eps to the other side
            !intsol(3)=( vec(5) <= vdif(2)-vdif(1) + eps  )      !corner where c1=wsum and c2=0 !ahumarch2122 ahu032122 moving eps to the other side
        !!!!!    if ( intsol(1) ) then     !IF welldef, then check whether the solution is interior or corner  
        !!!!!        transfers(1) = alf * surplus -  vdif(1)                     !ahumarch2122 ahu032122 replacing with vdif to save time                                           
        !!!!!        transfers(2) = (1.0_dp-alf) * (vec(5) + vdif(1) ) - alf*vdif(2) !ahumarch2122 ahu032122 replacing with vdif to save time   
                !vmarioj(1:2,j)=vec(1:2)+alf*surplusj(j)
                !vmario(1:2)=vec(3:4)+transfers(1:2)
                !if (vec(5)-abs(vdif(1)-vdif(2)) < 0.0_dp ) then !ahumarch1522 if the transfers statement is correct, then this should be correct too!
                !    print*, "There is something wrong!", vec(5),abs(vdif(1)-vdif(2))  !ahumarch1522
                !    stop                                                    !ahumarch1522
                !end if                                                      !ahumarch1522 
                !if ( minval(transfers) + eps2 < 0.0_dp .or. maxval(transfers) > vec(5)+eps2 ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
                !    print*, "Int cond fine but transfers not fine!", transfers(1:2),vec(5)  !ahumarch1522
                !    stop                                                                    !ahumarch1522                            
                !end if 
                !if ( minval(transfers) < 0.0_dp .or. maxval(transfers) > vec(5) ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
                !    print*, "Int cond fine but transfers not fine 2!", transfers(1:2),vec(5)  !ahumarch1522
                !    stop                                                                    !ahumarch1522                            
                !end if 
        !!!!!    else if ( intsol(2) ) then                      !ahumarch1522 adding cornersol
        !!!!!        transfers(1)=0.0_dp                         !ahumarch1522 adding cornersol
        !!!!!        transfers(2)=vec(5)                         !ahumarch1522 adding cornersol
                !vmario(1:2)=vec(3:4)+transfers(1:2)
                !vec(1:2)=vec(3:4)+transfers(1:2)
                !if (vec(5) >= vdif(1)-vdif(2) ) then        !ahumarch1522 adding cornersol
                !    print*, "There is somethingwrong in corner1!", vec(5),vdif(1)-vdif(2)  !ahumarch1522 adding cornersol
                !    stop
                !end if                                       !ahumarch1522 adding cornersol
        !!!!!    else if ( intsol(3) ) then                      !ahumarch1522 adding cornersol
        !!!!!        transfers(1)=vec(5)                         !ahumarch1522 adding cornersol
        !!!!!        transfers(2)=0.0_dp                         !ahumarch1522 adding cornersol
                !vmario(1:2)=vec(3:4)+transfers(1:2)
                !if (vec(5) >= vdif(2)-vdif(1) ) then        !ahumarch1522 adding cornersol
                !    print*, "There is somethingwrong in corner2!", vec(5),vdif(1)-vdif(2)  !ahumarch1522 adding cornersol
                !    stop
                !end if                                       !ahumarch1522 adding cornersol
        !!!!!    else 
        !!!!!        print*, "I should not end up here because we already established that we haveenoughtomakeindiff"
        !!!!!        stop
        !!!!!    end if
        !!!!!    vmax(1:2)=vec(3:4)+transfers(1:2)
        !!!!!    if (yazmax) then !ahu030822
        !!!!!        write(201,'(2x,tr1,"age",tr2,"q0",tr3,"q",tr2,"qt",tr2,"ie",TR3,"X",tr1,"rel")' ) !ahu030822
                !write(201,'(2x,7i4)')	ia,q0,q,qmax,iepsmove,X,relmax    ! altq is just the q that altrnative j corresponds to       !ahu030822
        !!!!!    end if !ahu030822
        !!!!!else  if (relmax==0) then !Second condition for welldef NB not well defined either because don't have enough to make indiff
        !!!!!    vmax(1:2)=vec(1:2)   
        !!!!!else 
        !!!!!    print*, "I should not end up here because we already established that we haveenoughtomakeindiff"
        !!!!!    stop   
        !!!!!end if 
    



! Look at dropbox/familymig/v130814/Source2 and Source11 and temp_010413 for the old versions of this
	subroutine checknb( dd, vec , def ,vsum)   
	integer, intent(in) :: dd(:) 
	real(dp), intent(in) :: vec(5)        
	logical, intent(out) :: def
	real(dp), intent(out) :: vsum(2)
	real(dp) :: surplus,transfers(2)
	vsum = pen ; transfers=pen 
	def  = .FALSE.
    surplus = vec(5) + vec(3) - vec(1) + vec(4) - vec(2)  
	if ( surplus + eps >= 0.0_dp )  then
		transfers(1) = alf * surplus - ( vec(3)-vec(1) )                                                    ! wage
		transfers(2) = (1.-alf) * (vec(5) + vec(3) - vec(1) ) - alf*( vec(4)-vec(2) )                       ! firm payoff
        !if (transfers(1) < 0.0_dp) cornersol=(/ 0._dp,vec(5) /) 
		!if (transfers(2) < 0.0_dp) cornersol=(/ vec(5),0._dp /) 
        !if not corner if statement. but not sure if this is right. look at the other tries below. 
        !because the below checks for corner solution o rinterior solution (checking transfers>0 or <w afterwards is the same as 
        !by the way checking the condition for checking FOC at corners 0 and w. see notes p.3 . 
        !but before this you need to figure otu whether there are 
        !any feasible utility allocations through doing that checking: 
        !asum = sum(  abs(.not. pc)  *   abs( vdif )   )
	    !def  = ( w + eps - asum  >= 0.0_dp )
    
		if ( minval(transfers) + eps >= 0.0_dp .and. maxval(transfers) < vec(5) ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
            def=.TRUE.
		    vsum(1) = transfers(1) + vec(3)      ! vsum is the value function after all transfers are made. i.e. vsum = wage + beta*EV_worker         ! vec(3) is just beta*EV i.e. what the worker has before transfers are made
		    vsum(2) = transfers(2) + vec(4)      ! vsum = firmpayoff + beta*EV_firm     ! " 
            if ( vsum(1) + eps < vec(1) .or. vsum(2) + eps < vec(2) ) then
                print*, "Transfers positive and less than wages but someone is getting a bad deal "
		        write(*,'(2x,2(tr6,"vbar"),2(tr8,"vc"),tr8,"wc",tr6,"surp")' )    
		        write(*,'(2x,6f10.2)') vec,surplus
		        write(*,'(2x,2(tr5,"trans"),2(tr6,"vsum") )' )    
		        write(*,'(2x,4f10.2)') transfers,vsum
                write(*,*) vsum(1)-vec(1),vsum(2)-vec(2)
            end if 
        end if
	end if
	end subroutine checknb

	!if ( pc(1) .and. pc(2) ) then	
	!	def  = .true.				
	!else 
	!	asum = sum(  abs(.not. pc)  *   abs( vdif )   )
	!	def  = ( w + eps - asum  >= 0.0_dp )
	!end if 
	!if (def) then 
	!	vsum	= vec(1:2) + 0.5_dp * ( w + eps + sum(vdif) )
	!	if ( vsum(1) < vec(1) .or. vsum(2) < vec(2) ) then 
	!		print*, "error in checknb " ,vsum(1)-vec(1),vsum(2)-vec(2),( w - asum  >= 0.0_dp ),w,asum,pc,w + sum(vdif) 
	!		print*, "vdif(2),wf_s,wf_s+vdif(2),wf_s+vdif(2) ", vdif(2)  !,wagetemp(2),wagetemp(2) + vdif(2),wagetemp(2) + vec(4) - vec(2) 
	!		stop 
	!	end if
	!end if 
	!end subroutine checknb

	!vsum = penaltyval 
	!transfers = penaltyval
	!def  = .FALSE.
	!criter=.FALSE.
	!surplus = vec(5) + vec(3) - vec(1) + vec(4) - vec(2)  
	! CASE 1: NO CONSTRAINT ON TRANSFERS. I.E. WAGES AND FIRM SHARES CAN BE NEGATIVE. 
	!         When there are no limits on what can be transferred, the only criterion for NB to be well-defined is simply that the match surplus is positive.  
	!if (nonneg_constraint==0) then                      ! no nonneg constraints for wages nor firm payoff  
	!	def = ( surplus + eps >= 0.0 )                  ! When there are no limits on what can be transferred btw parties, the only criterion for NB to be well defined is that surplus is positive. 
		
	! CASE 2: NONNEGATIVITY CONSTRAINTS ON TRANSFERS: 
	!           So that transfers can only come from current period resources. Cannot borrow from the continuation values to pay the other party today. 
	!           This limits the set of of feasible allocations and therefore the number of times NB is welldefined is less in this case. Because now any divsion of total surplus can only be implemented by using current period resources, which might not be enough sometimes. 
	!else if (nonneg_constraint==2) then                 ! see notes for how the nonneg constraints mean that the match surplus (outputnet) needs to satisfy the following three criteria
	!	criter(1) = ( vec(5) + eps >= vec(1) - vec(3) ) 
	!	criter(2) = ( vec(5) + eps >= vec(2) - vec(4) ) 
	!	criter(3) = ( vec(5) + eps >= vec(1) - vec(3) + vec(2) - vec(4) ) 
	!	def = ( criter(1) .and. criter(2) .and. criter(3) ) 
	!end if 
	

	! Check interior opt conditions 
	! If NB is well defined, then take FOC to get optimal transfers (wage and firm payoffs)
	!if (def) then 
	!	transfers(1) = alpha * surplus - ( vec(3)-vec(1) )                                                      ! wage
	!	transfers(2) = (1.-alpha) * (vec(5) + vec(3) - vec(1) ) - alpha*( vec(4)-vec(2) )                       ! firm payoff
	!	if (transfers(1) < 0.) transfers=(/ 0.,vec(5) /) 
	!	if (transfers(2) < 0.) transfers=(/ vec(5),0. /) 
	!	vsum(1) = transfers(1) + vec(3)      ! vsum is the value function after all transfers are made. i.e. vsum = wage + beta*EV_worker         ! vec(3) is just beta*EV i.e. what the worker has before transfers are made
	!	vsum(2) = transfers(2) + vec(4)      ! vsum = firmpayoff + beta*EV_firm     ! " 
	!	if (writeval) then ! if want to check whether things look right: check whether vsum's and the below temp's are equal. they should be since they are the same way of calculating value functions after the transfers.        
	!   end if 
	!end if 
