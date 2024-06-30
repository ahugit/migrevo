! train_stderr.f90
! Steven Laufer, 11/12/15
! Calculate standard errors

MODULE train_se

USE train_header
USE train_solver
USE train_simulate
USE train_moments

IMPLICIT NONE

CONTAINS

SUBROUTINE untranspars(parvector,rpv)
    REAL, DIMENSION(npars), INTENT(IN) :: parvector
    REAL, DIMENSION(npars), INTENT(OUT) :: rpv ! real parvector

    rpv(1)=logit(parvector(1))
    rpv(2)=2*parvector(2)
    rpv(3)=EXP(parvector(3))
    rpv(4)=EXP(parvector(4))
    rpv(5)=EXP(parvector(5))
    IF (NE==1) THEN
        rpv(6)=EXP(parvector(6))
        rpv(7)=EXP(parvector(7))
    ELSE
        rpv(6)=EXP(parvector(6))
        rpv(7)=EXP(parvector(7))
    ENDIF
    rpv(8)=logit(parvector(8))
    rpv(9)=logit(parvector(9))
    IF (NE==1) THEN
        rpv(10)=-1.+2*logit(parvector(10))
        rpv(11)=-1.+2*logit(parvector(11))
    ELSE
        rpv(10)=EXP(parvector(10))
        rpv(11)=EXP(parvector(11))
    ENDIF
    rpv(12)=EXP(parvector(12))
    rpv(13)=EXP(parvector(13))
    rpv(14)=EXP(parvector(14))
    rpv(15)=parvector(15)
    rpv(16)=EXP(parvector(16))
    rpv(17)=parvector(17)
    rpv(18)=parvector(18)
    rpv(19)=parvector(19)
    rpv(20:19+ntauobspar)=parvector(20:19+ntauobspar)
    IF (NE>1) THEN
        rpv(21+ntauobspar)=EXP(parvector(21+ntauobspar))
        rpv(22+ntauobspar)=EXP(parvector(22+ntauobspar))
    ENDIF
    rpv(20+ntauobspar)=EXP(parvector(20+ntauobspar))

END SUBROUTINE

! Calculate standard errors. (Only master will have computed values.)
SUBROUTINE train_stderr(parvector,sestepsize,stderrs)
    REAL, DIMENSION(npars), INTENT(IN) :: parvector ! estimated parameter vector at which to calculate standard errors
    REAL, DIMENSION(npars), INTENT(IN) :: sestepsize ! Stepsize to use for numerical derivatives. Zero if parameter wasn't estimated
    REAL, DIMENSION(npars), INTENT(OUT) :: stderrs ! calculated standard errors
   
    REAL, DIMENSION(Nmoments) :: simmoments,MSM_weights,varsimmom,weights, simmoments0,QQ,vardatamom,datamoments
    INTEGER, DIMENSION(Nmoments) :: countsimmom,activepari,countdatamom
    CHARACTER(LEN=40), DIMENSION(Nmoments) :: momentname
    REAL, DIMENSION(npars,npars) :: DpWDinv ! (D'WD)^-1
    REAL, DIMENSION(Nmoments,npars) :: D,WD,QWD
    INTEGER :: im, kactive, kk,eflag
    REAL, DIMENSION(npars) :: parvecbump,realparvecbump,realparvector,dtheta
    REAL, DIMENSION(npars,npars) :: SEmat
    TYPE(train_obs),DIMENSION(maxt,nsimtot) :: simdat
    TYPE(train_obs),DIMENSION(maxt,Ndata) :: realdat
    CHARACTER(LEN=22), DIMENSION(npars) :: parnames

        CALL get_data_outcomes(realdat)
        CALL wage_percentiles(realdat,wageps)
        CALL train_mom(realdat,Ndata,datamoments,countdatamom,vardatamom,momentname,weights)
        !QQ=(datamoments**2*(countdatamom/Ndata-(countdatamom/Ndata)**2)+vardatamom*countdatamom/Ndata)/Ndata
        QQ=vardatamom/Ndata
        DO im=1,Nmoments
            IF (vardatamom(im)*(countdatamom(im))>0.) THEN
                MSM_weights(im)=(countdatamom(im))*vardatamom(im)**(-1)
            ELSEIF (countdatamom(im)==0) THEN
                MSM_weights(im)=0.
            ELSE
                MSM_weights(im)=1.
            ENDIF
        ENDDO
        OPEN(13,FILE='stderrors'//runid//'.txt',STATUS='REPLACE')
    ! solve and get moments at estimated parameter (stored in simmoments0)
    CALL getpars(parvector)
    CALL train_setup
    CALL train_solve
    CALL train_sim
    CALL get_sim_outcomes(simdat)
    CALL nonmissing_wages(simdat)
    CALL train_mom(simdat,Nsimtot,simmoments0,countsimmom,varsimmom,momentname,weights)
    CALL train_diffmom(momentname,simmoments0,datamoments,countdatamom,vardatamom,weights,MSM_weights)
    CALL untranspars(parvector,realparvector)
    kactive=0 ! count of active parameters
    DO kk=1,npars
        IF (sestepsize(kk) .NE. 0.0)  THEN
            parvecbump=parvector
            parvecbump(kk)=parvector(kk)+sestepsize(kk)
            CALL untranspars(parvecbump,realparvecbump)
            dtheta(kk)=realparvecbump(kk)-realparvector(kk)
            CALL getpars(parvecbump)
            CALL train_setup
            CALL train_solve
            CALL train_sim
            CALL get_sim_outcomes(simdat)
            CALL nonmissing_wages(simdat)
            CALL train_mom(simdat,Nsimtot,simmoments,countsimmom,varsimmom,momentname,weights)
            CALL train_diffmom(momentname,simmoments,datamoments,countdatamom,vardatamom,weights,MSM_weights)
            IF (MAXVAL(ABS(QQ*weights*MSM_weights*(simmoments-simmoments0)))==0) THEN
                WRITE(13,'(1A22,1A16,1F8.5)') parnames(kk), ' not identified ', dtheta(kk)
            ELSE
                kactive=kactive+1 ! number of parameters actually iterating on.
                activepari(kactive)=kk ! indexes of active parameters
                D(:,kactive)=(simmoments-simmoments0)/dtheta(kk) ! derivative of moments wrt paramter
                WD(:,kactive)=weights*MSM_weights*D(:,kactive)
                QWD(:,kactive)=QQ*WD(:,kactive)
            ENDIF
        ENDIF
    ENDDO
    ! kactive now holds total number of active parameters
    !Q_MSM=SUM(weights*MSM_weights*(datamoments-simmoments)**2)
    !FINDInv(matrix, inverse, n, errorflag)
    CALL writepars(parnames)
        CALL FINDinv(MATMUL(TRANSPOSE(WD(:,1:kactive)),D(:,1:kactive)),DpWDinv(1:kactive,1:kactive),kactive,eflag) ! get (D'WD)^-1
        ! Standard errors are sqrt of diagional elements of matrix ((D'WD)^-1 * D'WQWD * (D'WD)^-1)
        SEmat(1:kactive,1:kactive)=MATMUL(MATMUL(DpWDinv(1:kactive,1:kactive),TRANSPOSE(WD(:,1:kactive))),MATMUL(QWD(:,1:kactive),DpWDinv(1:kactive,1:kactive)))
        stderrs=-1. ! default value to indicate that the paramter was not estimated
        DO kk=1,kactive
            stderrs(activepari(kk))=SQRT(SEmat(kk,kk))
            WRITE(13,'(1A22,1F12.6,1A5,1F6.3,1A1)') parnames(activepari(kk)),stderrs(activepari(kk)),'    (',dtheta(activepari(kk)),')'
        ENDDO
        CLOSE(13)
        ! TESTING
        WRITE(*,*) 'D'
        DO im=1,nmoments
            WRITE(*,'(2F12.5)') D(im,1:2)
        ENDDO
        WRITE(*,*)
        WRITE(*,*) 'WD' 
        DO im=1,nmoments
            !WRITE(*,'(2F12.5)') WD(im,1:2)
            !WRITE(*,'(4F12.5)') WD(im,1),weights(im), MSM_weights(im),D(im,1)
            WRITE(*,'(5F14.5)') MSM_weights(im),real(countdatamom(im)),real(Ndata),vardatamom(im),vardatamom(im)**(-1) ! weighting vector (main diagonal of diaginal weighting matrix)
        ENDDO
        WRITE(*,*)
        WRITE(*,*) 'QWD'
        DO im=1,nmoments
            WRITE(*,'(1F12.5)') QWD(im,1)
        ENDDO
        WRITE(*,*)
        WRITE(*,*) 'DpWDinv'
        DO im=1,1
            WRITE(*,'(1F12.5)') DpWDinv(im,1)
        ENDDO
        WRITE(*,*)
        WRITE(*,*) 'SEmat'
        DO im=1,1
            WRITE(*,'(1F12.5)') SEmat(im,1:1)
        ENDDO
        WRITE(*,*)

END SUBROUTINE

END MODULE
