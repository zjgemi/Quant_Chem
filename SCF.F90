PROGRAM MAIN
! FORTRAN 90 code for Hatree-Fock method
! Framework written by Xinzijian Liu
! Gaussian integration written by Ning Zhang
! DIIS written by Yuhang Yao
USE gauss_integration_s
IMPLICIT NONE
INTEGER, PARAMETER :: double = 8
REAL(KIND=double) :: autoang = 0.52917720859d0
INTEGER, ALLOCATABLE :: nutyp(:), nbas(:)
REAL(KIND=double), ALLOCATABLE :: cent(:,:), alp(:)
REAL(KIND=double), ALLOCATABLE :: S(:,:), T(:,:)
REAL(KIND=double), ALLOCATABLE :: V1e(:,:), F1e(:,:)
REAL(KIND=double), ALLOCATABLE :: V2e(:,:,:,:), F2e(:,:,:), Fnow(:,:)  ! modified by DIIS to save F of each step
REAL(KIND=double), ALLOCATABLE :: E(:), C(:,:)
REAL(KIND=double), ALLOCATABLE :: cor(:,:)
REAL(KIND=double), ALLOCATABLE :: Dnow(:,:), Dold(:,:)
REAL(KIND=double), ALLOCATABLE :: X(:,:), work(:)
REAL(KIND=double), ALLOCATABLE :: orth(:,:), eigval(:)
REAL(KIND=double), ALLOCATABLE :: ErrF(:,:,:)   ! DIIS ErrF 
REAL(KIND=double) :: Enow, Eold, eps1, eps2, Te, Vne, Vee, Vnn
INTEGER :: i, j, k, l, kmax, ierr
INTEGER :: natom, charge, mult, nbastot, nele, nocc
LOGICAL :: fexist, fdiis

!STO-3G
REAL(KIND=double) :: a1s(3), d1s(3), a2sp(3), d2s(3), d2p(3)
REAL(KIND=double) :: z1s(10), z2sp(10)

NAMELIST/ctr/ natom, charge, mult, kmax, eps1, eps2, fdiis

CALL INIT
! Calculate integrals
CALL init_param
DO i = 1, nbastot
   DO j = 1, i
      S(i,j) = overlap_s(alp(i),cent(:,i),alp(j),cent(:,j))
      S(j,i) = S(i,j)

      T(i,j) = kinetic_s(alp(i),cent(:,i),alp(j),cent(:,j))
      T(j,i) = T(i,j)

      V1e(i,j) = 0d0
      DO k = 1, natom
         V1e(i,j) = V1e(i,j) - nutyp(k) * e_attract_s(alp(i),cent(:,i),alp(j),cent(:,j),cor(:,k))
      END DO
      V1e(j,i) = V1e(i,j)
   END DO
END DO
DO i = 1, nbastot
   DO j = 1, i
      DO k = 1, i
         DO l = 1, k
            V2e(i,j,k,l) = e_repulsion_s(alp(i),cent(:,i),alp(j),cent(:,j),alp(k),cent(:,k),alp(l),cent(:,l))
            V2e(j,i,k,l) = V2e(i,j,k,l)
            V2e(i,j,l,k) = V2e(i,j,k,l)
            V2e(j,i,l,k) = V2e(i,j,k,l)
            V2e(k,l,i,j) = V2e(i,j,k,l)
            V2e(l,k,i,j) = V2e(i,j,k,l)
            V2e(k,l,j,i) = V2e(i,j,k,l)
            V2e(l,k,j,i) = V2e(i,j,k,l)
         END DO
      END DO
   END DO
END DO
F1e = T + V1e

! Calculate transformation matrix X = S^(-1/2)
orth = S
CALL dsyev('V','L',nbastot,orth,nbastot,eigval,work,nbastot*10,ierr)
DO i = 1, nbastot
   IF(eigval(i) .eq. 0d0) THEN
      WRITE(*,*) 'No.', i, 'eigenvalue of S is 0.'
      STOP
   END IF
   X(:,i) = orth(:,i) * 1d0/sqrt(eigval(i))
END DO
X = matmul(X,transpose(orth))

k = 1
DO WHILE(.TRUE.)

   DO i = 1, nbastot
      DO j = 1, nbastot
         F2e(i,j,k) = F1e(i,j) + sum(Dnow(:,:)*(V2e(i,j,:,:)*2-V2e(i,:,j,:))) ! modified by DIIS to save F of each step
      END DO
   END DO

   Fnow = F2e(:,:,k)
   IF(fdiis) THEN
      ErrF(:,:,k) = MATMUL(MATMUL(F2e(:,:,k),Dnow),S) - MATMUL(MATMUL(S,Dnow),F2e(:,:,k)) ! calculate the ErrF of No. k step for DIIS
      CALL DIIS(F2e,Fnow)
   END IF
   CALL EIGEN(Fnow,E,C,Dnow,Enow)

   IF(sum((Dnow-Dold)**2) <= eps1) THEN
      WRITE(*,*) 'Density converges.'
      EXIT
   ELSE IF(abs(Enow-Eold) <= eps2) THEN
      WRITE(*,*) 'Energy converges.'
      EXIT
   ELSE IF (k > kmax) THEN
      WRITE(*,*) 'Iteration exceeds max number.'
      EXIT
   END IF
   Dold = Dnow
   Eold = Enow
   k = k + 1
END DO

WRITE(*,*) 'HF Energy:', Enow
Te = trace(matmul(Dnow,T))*2
Vne = trace(matmul(Dnow,V1e))*2
Vee = Enow - Te - Vne
Vnn = 0d0
DO i = 1, natom
   DO j = 1, i - 1
      Vnn = Vnn + nutyp(i)*nutyp(j)/sqrt(sum((cor(:,i)-cor(:,j))**2))
   END DO
END DO
WRITE(*,*) 'Te:', Te
WRITE(*,*) 'Vne:', Vne
WRITE(*,*) 'Vee:', Vee
WRITE(*,*) 'Vnn:', Vnn
WRITE(*,*) 'Etot:', Enow + Vnn
WRITE(*,*) 'Job finishes.'

DEALLOCATE(S,T,V1e,V2e,F1e,F2e,E,C,Dnow,Dold,ErrF)

CONTAINS


SUBROUTINE INIT
IMPLICIT NONE

! Read input
INQUIRE(FILE='ctrfile.dat',EXIST=fexist)
IF(fexist /= .TRUE.) THEN
   WRITE(*,*) "Error : ctrfile not found."
   STOP
ENDIF
OPEN(10,FILE='ctrfile.dat')
READ(10,ctr)
CLOSE(10)

ALLOCATE(nutyp(natom))
ALLOCATE(cor(3,natom))
ALLOCATE(nbas(natom))

! Read configuration
INQUIRE(FILE='config.dat',EXIST=fexist)
IF(fexist /= .TRUE.) THEN
   WRITE(*,*) "Error : config not found."
   STOP
ENDIF
OPEN(10,FILE='config.dat')
DO i = 1, natom
   READ(10,*) nutyp(i)
   READ(10,*) cor(1:3,i)
END DO
CLOSE(10)
cor = cor / autoang

nele = sum(nutyp) - charge
IF(nele < 0) THEN
   WRITE(*,*) 'Wrong electron number.'
   STOP
END IF

IF((mod(nele,2) .eq. mod(mult,2)) .or. mult <= 0&
   &.or. mult >= nele + 1) THEN
   WRITE(*,*) 'Wrong spin multiplicity.'
   STOP
END IF
nocc = (nele - (mult - 1)) / 2 + mult - 1

DO i = 1, natom
   IF(nutyp(i)>=1 .and. nutyp(i)<=2) THEN
      nbas(i) = 3
   ELSE IF(nutyp(i)>=3 .and. nutyp(i) <=10) THEN
      nbas(i) = 6
   ELSE
      WRITE(*,*) 'Wrong nuclear type.'
      STOP
   END IF
END DO
nbastot = sum(nbas)

ALLOCATE(cent(3,nbastot), alp(nbastot))
ALLOCATE(S(nbastot,nbastot), T(nbastot,nbastot))
ALLOCATE(V1e(nbastot,nbastot), V2e(nbastot,nbastot,nbastot,nbastot))
ALLOCATE(F1e(nbastot,nbastot), F2e(nbastot,nbastot,kmax))  ! modified by DIIS to save F of each step
ALLOCATE(Fnow(nbastot,nbastot))  ! modified by DIIS to save F of each step
ALLOCATE(orth(nbastot,nbastot), eigval(nbastot))
ALLOCATE(work(nbastot*10), X(nbastot,nbastot))
ALLOCATE(E(nbastot), C(nbastot,nbastot))
ALLOCATE(Dnow(nbastot,nbastot), Dold(nbastot,nbastot))
ALLOCATE(ErrF(nbastot,nbastot,kmax))  ! Allocate ErrF for DIIS

! Initialize basis set
WRITE(*,*) 'Basis: STO-3G'
INQUIRE(FILE='param.dat',EXIST=fexist)
IF(fexist /= .TRUE.) THEN
   WRITE(*,*) "Error : param not found."
   STOP
ENDIF
OPEN(10,FILE='param.dat') ! P184 of Szabo's book
READ(10,*) a1s(:)
READ(10,*) d1s(:)
READ(10,*) a2sp(:)
READ(10,*) d2s(:)
READ(10,*) d2p(:)
DO i = 1, 10
   READ(10,*) z1s(i), z2sp(i)
END DO
k = 0
DO i = 1, natom
   IF(nutyp(i)>=1 .and. nutyp(i)<=2) THEN
      DO j = 1, 3
         cent(:,k+j) = cor(:,i)
      END DO
      alp(k+1:k+3) = a1s(:)
      k = k + 3
   ELSE IF(nutyp(i)>=3 .and. nutyp(i) <=10) THEN
      DO j = 1, 6
         cent(:,k+j) = cor(:,i)
      END DO
      alp(k+1:k+3) = a1s(:)
      alp(k+4:k+6) = a2sp(:)
      k = k + 6
   END IF
END DO
WRITE(*,*) 'Job starts.'

! Initial guess
Dnow = 0d0
DO i = 1, nocc
   Dnow = 1d0
END DO
Dold = Dnow

END SUBROUTINE INIT


SUBROUTINE DIIS(F2e,Fnow)
IMPLICIT NONE
REAL(KIND=double), INTENT(INOUT) :: F2e(:,:,:)
REAL(KIND=double), INTENT(OUT) :: Fnow(:,:)
REAL(KIND=double), ALLOCATABLE :: Bnow(:,:),Bequl(:),Cdiis(:), work2(:)
INTEGER, ALLOCATABLE :: ipiv(:)
ALLOCATE(Bnow(k+1,k+1),Bequl(k+1),Cdiis(k+1),work2(10*k+10),ipiv(k+1))
Bequl = 0d0
Bequl(k+1) = -1d0
DO i = 1, k
   DO j = 1, k
      Bnow(i,j)=trace(matmul(ErrF(:,:,i),ErrF(:,:,j)))
   END DO
END DO
Bnow(:,k+1) = -1d0
Bnow(k+1,:) = -1d0
Bnow(k+1,k+1) = 0d0
Cdiis = Bequl
CALL dsysv('U',k+1,1,Bnow,k+1,ipiv,Cdiis,k+1,work2,10*k+10,ierr)
Fnow = 0d0
DO i = 1, k
   Fnow = Fnow + Cdiis(i)*F2e(:,:,i)
END DO
DEALLOCATE(Bnow,Bequl,Cdiis,work2,ipiv)

RETURN
END SUBROUTINE DIIS


SUBROUTINE EIGEN(Fnow,E,C,Dnow,Enow)
IMPLICIT NONE
REAL(KIND=DOUBLE), INTENT(IN) :: Fnow(:,:)
REAL(KIND=DOUBLE), INTENT(OUT) :: E(:), C(:,:), Dnow(:,:), Enow
C = matmul(X,Fnow(:,:))
C = matmul(C,X)
CALL dsyev('V','L',nbastot,C,nbastot,E,work,nbastot*10,ierr)
C = matmul(X,C)
Dnow = matmul(C(:,1:nocc),transpose(C(:,1:nocc)))
Enow = trace(matmul(Dnow,F1e+F2e(:,:,k)))
WRITE(*,'(A,I5,A,F15.7)') ' Step', k, '   Energy', Enow
END SUBROUTINE EIGEN


FUNCTION trace(mat)
IMPLICIT NONE
REAL(KIND=double) :: mat(:,:), trace
INTEGER :: ndim, i
ndim = size(mat,1)
IF(size(mat,2) /= ndim) THEN
   WRITE(*,*) 'Not a square matrix.'
END IF
trace = 0d0
DO i = 1, ndim
   trace = trace + mat(i,i)
END DO
RETURN
END FUNCTION trace


END PROGRAM
