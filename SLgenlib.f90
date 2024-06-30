


MODULE SLgenlib
! Steven Laufer, Updated 6/10/10
! Library of general use functions and subroutines
IMPLICIT NONE
integer :: dint1=1
CONTAINS

    ! linear inperpolation
	REAL(8) FUNCTION interplin(x,y,xi)

		REAL(8), INTENT(IN), DIMENSION(:) :: x
		REAL(8), INTENT(IN), DIMENSION(:) :: y
		REAL(8), INTENT(IN) :: xi

		INTEGER :: ii
		REAL(8) :: z

		IF (SIZE(x)/=SIZE(Y)) THEN
			!error
		END IF
		ii=2
		DO
			IF ((xi<=x(ii)) .OR. (ii>=SIZE(x))) EXIT
			ii=ii+1
		END DO
		ii=ii-1
		z=(xi-x(ii))/(x(ii+1)-x(ii))
		interplin=(1-z)*y(ii)+z*y(ii+1)
	
	END FUNCTION

	! one-dimensional nearest neighbor interpretation
	REAL(8) FUNCTION interpnn(x,y,xi)

		REAL(8), INTENT(IN), DIMENSION(:) :: x
		REAL(8), INTENT(IN), DIMENSION(:) :: y
		REAL(8), INTENT(IN) :: xi

		INTEGER :: ii
		REAL(8) :: z

		IF (SIZE(x)/=SIZE(Y)) THEN
			!error
		END IF
		ii=2
		DO
			IF ((xi<=x(ii)) .OR. (ii>=SIZE(x))) EXIT
			ii=ii+1
		END DO
		ii=ii-1
		z=(xi-x(ii))/(x(ii+1)-x(ii))
		interpnn=y(ii+dint1*(z>0.5)) !ag 120422: will not compile will just the parenthesis so added dint1 which is decalred at head of module
	
	END FUNCTION

    ! linear interpolation acting on vector
	REAL(8) FUNCTION interplinv(x,y,xi,nxi)

		REAL(8), INTENT(IN), DIMENSION(:) :: x
		REAL(8), INTENT(IN), DIMENSION(:) :: y
		REAL(8), INTENT(IN), DIMENSION(:) :: xi
		INTEGER, INTENT(IN) :: nxi
		DIMENSION interplinv(nxi)
		!REAL(8), INTENT(OUT), DIMENSION(nxi) :: interplinv

		INTEGER :: n
		INTEGER :: nn
	
		INTEGER, DIMENSION(nxi) :: in
		REAL(8), DIMENSION(nxi) :: z

		in=2
		nn=2
		DO	
			IF (nn>=SIZE(x)) EXIT
			WHERE(xi>x(in))
				in=in+1
			END WHERE
		
			nn=nn+1
		END DO
		in=in-1
	
		z=(xi-x(in))/(x(in+1)-x(in))
	
		interplinv=(1-z)*y(in)+z*y(in+1)
	
	END FUNCTION
					

	SUBROUTINE checklogic

		!IF (.TRUE. /= 1) THEN !ag 120422: will not compile will just the parenthesis so uncommenting this out because I can't
		!	WRITE(*,*) 'WARNING: Logical true does not equal one.'
		!END IF
		
	END SUBROUTINE

    ! kroniker product
	REAL(8) FUNCTION kron(A,B)
		REAL(8), INTENT(IN), DIMENSION(:,:) :: A
		REAL(8), INTENT(IN), DIMENSION(:,:) :: B
		DIMENSION kron(SIZE(A,1)*SIZE(B,1),SIZE(A,2)*SIZE(B,2))

		INTEGER :: n1,n2,sizea1,sizea2,sizeb1,sizeb2

		sizea1=SIZE(A,1)
		sizeb1=SIZE(B,1)
		sizea2=SIZE(A,2)
		sizeb2=SIZE(B,2)

		DO n1=1,sizea1
			DO n2=1,sizea2
				kron(sizeb1*(n1-1)+1:sizeb1*n1,sizeb2*(n2-1)+1:sizeb2*n2)=A(n1,n2)*B
			END DO
		END DO		

	END FUNCTION

	INTEGER FUNCTION firsttrue(A)
		! Finds index of first true element of A. Returns UBOUND(A)+1 if no true elements.
		LOGICAL, INTENT(IN), DIMENSION(:) :: A
		
		IF (ALL(.NOT. A)) THEN ! no true elements
			firsttrue=UBOUND(A,1)+1
		ELSE ! at least one true element
			firsttrue=LBOUND(A,1)
			DO
				IF (A(firsttrue)) EXIT
				firsttrue=firsttrue+1
			END DO
		END IF

	END FUNCTION

REAL(8) FUNCTION normcdf ( x )
! Downloaded 3/19/08 from http://orion.math.iastate.edu/burkardt/f_src/prob/prob.html


!*******************************************************************************
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!
!  Reference: 
!
!    A G Adams,
!    Areas Under the Normal Curve,
!    Algorithm 39, 
!    Computer j., 
!    Volume 12, pages 197-198, 1969.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, REAL(8) X, the argument of the CDF.
!
!    Output, REAL(8) CDF, the value of the CDF.
!
  implicit none
!
  REAL(8), parameter :: a1 = 0.398942280444E+00
  REAL(8), parameter :: a2 = 0.399903438504E+00
  REAL(8), parameter :: a3 = 5.75885480458E+00
  REAL(8), parameter :: a4 = 29.8213557808E+00
  REAL(8), parameter :: a5 = 2.62433121679E+00
  REAL(8), parameter :: a6 = 48.6959930692E+00
  REAL(8), parameter :: a7 = 5.92885724438E+00
  REAL(8), parameter :: b0 = 0.398942280385E+00
  REAL(8), parameter :: b1 = 3.8052E-08
  REAL(8), parameter :: b2 = 1.00000615302E+00
  REAL(8), parameter :: b3 = 3.98064794E-04
  REAL(8), parameter :: b4 = 1.98615381364E+00
  REAL(8), parameter :: b5 = 0.151679116635E+00
  REAL(8), parameter :: b6 = 5.29330324926E+00
  REAL(8), parameter :: b7 = 4.8385912808E+00
  REAL(8), parameter :: b8 = 15.1508972451E+00
  REAL(8), parameter :: b9 = 0.742380924027E+00
  REAL(8), parameter :: b10 = 30.789933034E+00
  REAL(8), parameter :: b11 = 3.99019417011E+00
  REAL(8) q
  REAL(8) x
  REAL(8) y
!
!  |X| <= 1.28.
!
  if ( abs ( x ) <= 1.28 ) then

    y = 0.5E+00 * x**2

    q = 0.5E+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7E+00 ) then

    y = 0.5E+00 * x**2

    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 &
      + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 &
      - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0E+00

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0E+00 ) then
    normcdf = q
  else
    normcdf = 1.0E+00 - q
  end if

  return
end function

	! 2 dim linear inerpolation
	REAL(8) FUNCTION interplin2(x1,x2,V,y1,y2)
	
		REAL(8), INTENT(IN), DIMENSION(:) :: x1,x2
		REAL(8), INTENT(IN), DIMENSION(:,:) :: V
		REAL(8), INTENT(IN) :: y1, y2
	
		INTEGER, DIMENSION(2) :: sizes ! sizes of each dim of V
		INTEGER :: ii, i1,i2 ! determine cell of array
		REAL(8) :: z1,z2
		INTEGER :: a1,a2

		sizes(1)=SIZE(x1)
		sizes(2)=SIZE(x2)

		ii=2
		DO
			IF (x1(ii)>=y1 .OR. ii>=sizes(1)) EXIT
			ii=ii+1
		END DO
		i1=ii-1
	
		ii=2
		DO
			IF (x2(ii)>=y2 .OR. ii>=sizes(2)) EXIT
			ii=ii+1
		END DO
		i2=ii-1

		z1=1.*(y1-x1(i1))/(x1(i1+1)-x1(i1))
		z2=1.*(y2-x2(i2))/(x2(i2+1)-x2(i2))
		
	
		interplin2=0

		DO a1=0,1
			DO a2=0,1
				interplin2=interplin2+(V(i1+a1,i2+a2))*(a1*z1+(1-a1)*(1-z1))*(a2*z2+(1-a2)*(1-z2))
			END DO
		END DO

	END FUNCTION

	! 2 dim nearest neighbor inerpolation
	REAL(8) FUNCTION interpnn2(x1,x2,V,y1,y2)
	
		REAL(8), INTENT(IN), DIMENSION(:) :: x1,x2
		REAL(8), INTENT(IN), DIMENSION(:,:) :: V
		REAL(8), INTENT(IN) :: y1, y2
	
		INTEGER, DIMENSION(2) :: sizes ! sizes of each dim of V
		INTEGER :: ii, i1, i2 ! determine cell of array
		REAL(8) :: z1,z2
		INTEGER :: a1,a2

		sizes(1)=SIZE(x1)
		sizes(2)=SIZE(x2)

		ii=2
		DO
			IF (x1(ii)>=y1 .OR. ii>=sizes(1)) EXIT
			ii=ii+1
		END DO
		i1=ii-1
	
		ii=2
		DO
			IF (x2(ii)>=y2 .OR. ii>=sizes(2)) EXIT
			ii=ii+1
		END DO
		i2=ii-1

		z1=1.*(y1-x1(i1))/(x1(i1+1)-x1(i1))
		z2=1.*(y2-x2(i2))/(x2(i2+1)-x2(i2))
		a1=dint1*(z1>.5) !ag 120422: will not compile will just the parenthesis so added dint1 which is decalred at head of module
		a2=dint1*(z2>.5) !ag 120422: will not compile will just the parenthesis so added dint1 which is decalred at head of module

		interpnn2=V(i1+a1,i2+a2)

	END FUNCTION



! compute covariance matrix of MLE estimator by outer product of gradients
SUBROUTINE MLEcov_opg(func,xval,step,cov,ndata)
	
	REAL(8), DIMENSION(:), INTENT(IN) :: xval ! value to evaluete at (ndim)
	REAL(8), DIMENSION(:), INTENT(IN) :: step ! step sizes for numerical derivatices (ndim)
	REAL(8), DIMENSION(:,:),INTENT(OUT) :: cov ! output cov matrix (nap x nap)
	EXTERNAL func ! function (xval, (Ndata x 1) vector of likelihoods)
	INTEGER, INTENT(IN) :: ndata ! number of observations

	REAL(8), DIMENSION(:,:), POINTER :: grad,ee
	REAL(8), DIMENSION(ndata) :: f0 ! function value at xval
	INTEGER :: ndim,nap,ii,astat,ncol

	ndim=SIZE(xval)
	nap = COUNT(step /= 0.)

	ALLOCATE (ee(1:ndim,1:ndim),STAT=astat)
	ee= 0
	DO ii=1,ndim
		ee(ii,ii)=step(ii)
	ENDDO

	IF ((SIZE(cov,1) /=nap) .OR. (SIZE(cov,2) /=nap)) THEN
		WRITE(*,*) 'MLEcov_opg: covariance matrix is wrong size'	
	ENDIF

	ALLOCATE (grad(1:ndata,1:nap),STAT=astat)
	CALL func(xval,f0)
	ncol=1
	DO ii=1,ndim
		IF (step(ii) /= 0.) THEN
			CALL func(xval+ee(ii,:),grad(:,ncol))
			grad(:,ncol)=(grad(:,ncol)-f0)/step(ii)
			ncol=ncol+1
		ENDIF
	ENDDO
	CALL FINDInv(MATMUL(TRANSPOSE(grad),grad)/ndata, cov, nap, astat)
	!cov=MATMUL(TRANSPOSE(grad),grad)/ndata

END SUBROUTINE

SUBROUTINE hessian(func,xval,step,H)
	
	REAL(8), DIMENSION(:), INTENT(IN) :: xval ! value to evaluete at (ndim)
	REAL(8), DIMENSION(:), INTENT(IN) :: step ! step sizes for numerical derivatices (ndim)
	REAL(8), DIMENSION(:,:),INTENT(OUT) :: H ! output cov matrix (nap x nap)
	REAL(8), EXTERNAL :: func ! function: function with vector input the length of xval

	REAL(8), DIMENSION(:,:), POINTER :: ee
	INTEGER :: ndim,nap,ii,jj,astat,nrow,ncol

	ndim=SIZE(xval)
	ALLOCATE (ee(1:ndim,1:ndim),STAT=astat)
	nap = COUNT(step /= 0.)

	IF ((SIZE(H,1) /=nap) .OR. (SIZE(H,2) /=nap)) THEN
		WRITE(*,*) 'hessian: output matrix is wrong size'	
	ENDIF
	 
	ee= 0
	DO ii=1,ndim
		ee(ii,ii)=step(ii)
	ENDDO
	
	nrow=1
	DO ii=1,ndim
		IF (step(ii) /= 0.) THEN
			ncol=1
			DO jj=1,ii
				IF (step(jj) /= 0.) THEN
					H(nrow,ncol)=(func(xval+ee(ii,:)+ee(jj,:))-func(xval+ee(ii,:)-ee(jj,:)) &
						-func(xval-ee(ii,:)+ee(jj,:))+func(xval-ee(ii,:)-ee(jj,:)))/(4*step(ii)*step(jj))
					H(ncol,nrow)=H(nrow,ncol)
					ncol=ncol+1
				ENDIF
			ENDDO
			nrow=nrow+1
		ENDIF
	ENDDO

	DEALLOCATE(ee,STAT=astat)


END SUBROUTINE

! pdf of standard normal
REAL(8) FUNCTION normpdf(x)

	REAL(8), INTENT(IN) :: x

	normpdf=(6.28318)**(-.5)*EXP(-x**2/2)
END FUNCTION


!Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm           
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
!Downloaded 7/23/08 from: http://www.dreamincode.net/code/snippet1272.htm

SUBROUTINE FINDInv(matrix, inverse, n, errorflag)
	IMPLICIT NONE
	!Declarations
	INTEGER, INTENT(IN) :: n
	INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
	REAL(8), INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
	REAL(8), INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
	
	LOGICAL :: FLAG = .TRUE.
	INTEGER :: i, j, k, l
	REAL(8) :: m
	REAL(8), DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
	! n=1 case added by SL, 7/23/08
	IF (n==1) THEN
		IF (matrix(1,1)==0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ELSE
			inverse(1,1) = 1./matrix(1,1)
			errorflag = 0
			return
		ENDIF
	ENDIF

	!Augment input matrix with an identity matrix
	DO i = 1, n
		DO j = 1, 2*n
			IF (j <= n ) THEN
				augmatrix(i,j) = matrix(i,j)
			ELSE IF ((i+n) == j) THEN
				augmatrix(i,j) = 1
			Else
				augmatrix(i,j) = 0
			ENDIF
		END DO
	END DO
	
	!Reduce augmented matrix to upper traingular form
	DO k =1, n-1
		IF (augmatrix(k,k) == 0) THEN
			FLAG = .FALSE.
			DO i = k+1, n
				IF (augmatrix(i,k) /= 0) THEN
					DO j = 1,2*n
						augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
					END DO
					FLAG = .TRUE.
					EXIT
				ENDIF
				IF (FLAG .EQV. .FALSE.) THEN
					PRINT*, "Matrix is non - invertible"
					inverse = 0
					errorflag = -1
					return
				ENDIF
			END DO
		ENDIF
		DO j = k+1, n			
			m = augmatrix(j,k)/augmatrix(k,k)
			DO i = k, 2*n
				augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
			END DO
		END DO
	END DO
	
	!Test for invertibility
	DO i = 1, n
		IF (augmatrix(i,i) == 0) THEN
			PRINT*, "Matrix is non - invertible"
			inverse = 0
			errorflag = -1
			return
		ENDIF
	END DO
	
	!Make diagonal elements as 1
	DO i = 1 , n
		m = augmatrix(i,i)
		DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / m)
		END DO
	END DO
	
	!Reduced right side half of augmented matrix to identity matrix
	DO k = n-1, 1, -1
		DO i =1, k
		m = augmatrix(i,k+1)
			DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
			END DO
		END DO
	END DO				
	
	!store answer
	DO i =1, n
		DO j = 1, n
			inverse(i,j) = augmatrix(i,j+n)
		END DO
	END DO
	errorflag = 0
END SUBROUTINE FINDinv

	! Find the standard deviation of a vector
	REAL(8) FUNCTION stdvec(x)

		REAL(8), INTENT(IN), DIMENSION(:) :: x
		
		REAL(8) :: meanvec
		INTEGER :: nn

		nn=SIZE(x)

		stdvec=SQRT(SUM((x-SUM(x)/nn)**2)/nn)

	END FUNCTION

    ! convert from n-dimensional index to linear one
    INTEGER FUNCTION ndim2lin(dims,nindex)
        INTEGER, DIMENSION(:), INTENT(IN) :: dims !like outpit of size(A)- vector
        INTEGER, DIMENSION(:), INTENT(IN) :: nindex !n-dimensional index- row vector
        ! indexout is linear index
        
        INTEGER, DIMENSION(SIZE(dims)-1) :: cumprod ! running cumulative product of dims
        INTEGER :: ii,ndims
        
        ndims=SIZE(dims)
        cumprod(1)=dims(1)
        DO ii=2,ndims-1
            cumprod(ii)=dims(ii)*cumprod(ii-1)
        ENDDO
        ndim2lin=nindex(1)+SUM(cumprod*(nindex(2:ndims)-1))
    END FUNCTION
    
    ! convert from linear index to n-dimensional one
    INTEGER FUNCTION lin2ndim(dims,linindex)
        INTEGER, DIMENSION(:), INTENT(IN) :: dims ! like outpit of size(A)
        INTEGER, INTENT(IN) :: linindex ! linear index
        
        DIMENSION lin2ndim(SIZE(dims)) 
        INTEGER :: nn, prod, linindexcopy,ndims       
        
        ndims=SIZE(dims)
        linindexcopy=linindex
        DO nn=NDims,2,-1
            !WRITE(*,*) linindexcopy
            prod=PRODUCT(dims(1:nn-1))
            lin2ndim(nn)=CEILING(REAL(linindexcopy)/prod)
            linindexcopy=MOD(linindexcopy,prod)
            IF (linindexcopy==0) THEN
                linindexcopy=prod
            ENDIF
        ENDDO
        lin2ndim(1)=linindexcopy
        
    END FUNCTION

    ! draw from a multinomal distribution.
    ! Inputs are the vector of probabilities for each outcome and a U[0,1] random draw
    ! Output is the index of the outcome
    INTEGER FUNCTION multinom(probvec,randdraw)
        REAL(8), DIMENSION(:), INTENT(IN) :: probvec ! probabilities of each state
        REAL(8), INTENT(IN) :: randdraw ! U[0,1]

        INTEGER :: nn
        REAL(8) :: cutoff
        REAL(8), DIMENSION(SIZE(probvec)) :: probvecnorm

        probvecnorm=probvec/SUM(probvec) ! normalize if probabilities don't sum to one
        
        ! cutoff is upper bound of interval that if randdraw falls in that interval, get current outcome
        cutoff=probvecnorm(1)
        nn=1
        DO WHILE (randdraw>cutoff)
            cutoff=cutoff+probvecnorm(nn+1)
            nn=nn+1
        ENDDO
        multinom=nn
        
    END FUNCTION
        
    ! gaussbin- takes in paramters of continuous normal distribution (mu, sigma) and vector
    ! of allowed values (values) and calculates fraction of outcomes closest to
    ! each element in vector, i.e. a discrete approximation to that normal distribution
    SUBROUTINE gaussbin(mu,sigma,values0,probvec)
        REAL(8), INTENT(IN):: mu ! mean of normal distribution
        REAL(8), INTENT(IN):: sigma ! standard deviation of normal distribution
        REAL(8), DIMENSION(:), INTENT(IN) :: values0 ! allowed discrete values
        REAL(8), DIMENSION(SIZE(values0)),INTENT(OUT) :: probvec ! probability weight at each allowed value
        REAL(8), DIMENSION(SIZE(values0)) :: values
		INTEGER :: nn,kk

		values=(values0-mu)/sigma
		
        
        nn=SIZE(values)
	IF (nn>1) THEN
        	probvec(1)=normcdf((values(1)+values(2))/2)
        	probvec(nn)=1-normcdf((values(nn)+values(nn-1))/2)
        	DO kk=2,nn-1,1
            		probvec(kk)=normcdf((values(kk)+values(kk+1))/2)-normcdf((values(kk)+values(kk-1))/2)
        	ENDDO
	ELSE
		probvec(1)=1.
	ENDIF

    END SUBROUTINE
            
    ! Use  algorithm from Tauchen(1986) to construct discrete approximation to an AR(1)
    SUBROUTINE tauchenAR1(rho,sigma,nstates,maxsd,statevec,PP)
        REAL(8), INTENT(IN) :: rho ! autocorrelation of AR(1)
        REAL(8), INTENT(IN) :: sigma ! SD of innovations of AR(1)
        INTEGER, INTENT(IN) :: nstates ! number of states in Markov proces
        REAL(8), INTENT(IN) :: maxsd !number of standard deviations from mean for first, last states 
        REAL(8), DIMENSION(nstates),INTENT(OUT) :: statevec ! Markov state vector
        REAL(8), DIMENSION(nstates,nstates),INTENT(OUT) :: PP ! Markov transition matrix
    
        INTEGER :: ii
        REAL(8) :: step
    
        statevec(1)=-1.0*maxsd*(sigma/SQRT(1-rho**2))
        step=2*maxsd*(sigma/SQRT(1-rho**2))/(nstates-1)
        DO ii=2,nstates
            statevec(ii)=statevec(ii-1)+step
        ENDDO
        DO ii=1,nstates
            CALL gaussbin(rho*statevec(ii),sigma,statevec,PP(ii,:))
        ENDDO

    END SUBROUTINE


END MODULE SLgenlib
		 
