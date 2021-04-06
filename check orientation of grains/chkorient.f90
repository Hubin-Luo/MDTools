	PROGRAM main
!=====================================================================
! Try to find out the orientation of a grain with non-cubic structure
! A uniaxial direction should be given.
! The centers of grains are given by the standard output of atomsk
!                                                              H.B.Luo
!=====================================================================
	IMPLICIT NONE
!
	INTEGER :: i,j,k,m, natom, ntype, iTypnum, nv
	REAL(8),ALLOCATABLE :: posseq(:,:)
	REAL(8) :: grainCenter(3),  celldim(3,2)
	REAL(8) :: dist, pairdist, temp(6), pairvec(3), orientvec(3)
	REAL(8) :: axisdist, matchthd, hfcell(3)
	REAL(8) :: pi, theta, anglethd, acute, obtuse
        CHARACTER(20) :: fnm
!
        pi=ACOS(-1.d0)
        READ(*,*) fnm
        READ(*,*) (grainCenter(i),i=1,3)
	READ(*,*) dist, iTypnum
        READ(*,*) axisdist, matchthd, anglethd
	acute=pi*anglethd; obtuse=pi*(1.d0-anglethd)
!
	OPEN(8,file=TRIM(fnm),status='old')
	READ(8,*); READ(8,*)
	READ(8,*) natom
	READ(8,*) ntype
	ALLOCATE(posseq(natom,6))
	READ(8,*)
	DO i=1,3
	   READ(8,*) (celldim(i,j),j=1,2)
           hfcell(i)=(celldim(i,2)-celldim(i,1))/2
	END DO
	READ(8,*); READ(8,*); READ(8,*)
	DO i=1, natom
	   READ(8,*) (posseq(i,j),j=1,5)
	END DO
	CLOSE(8)
!
!shift the grain center to the origin and calculate all modulus of position vecs
!
        DO i=1,natom
           posseq(i,3:5)=posseq(i,3:5)-grainCenter(:)
           DO j=1,3
              IF(posseq(i,j+2)>hfcell(j)) posseq(i,j+2)=posseq(i,j+2)-2.d0*hfcell(j)
              IF(posseq(i,j+2)<-hfcell(j))posseq(i,j+2)=posseq(i,j+2)+2.d0*hfcell(j)
           END DO
           posseq(i,6)=SQRT(DOT_PRODUCT(posseq(i,3:5),posseq(i,3:5)))
        END DO
!
!sort the positions up
!
        CALL quicksort(posseq, 1, natom, natom, 6, 6)
!
!find out all pairs of atoms of given type within dist and calculate the distance
!
        orientvec=0.d0; nv=0
        DO i=1,natom
           IF(posseq(i,6)>dist) THEN
              EXIT
           ELSE IF(NINT(posseq(i,2))/=iTypnum) THEN
              CYCLE
           ELSE
              DO j=i+1, natom
                 IF(posseq(j,6)>dist) THEN
                    EXIT
                 ELSE IF(NINT(posseq(j,2))/=iTypnum) THEN
                    CYCLE
                 ELSE
                    pairvec=posseq(j,3:5)-posseq(i,3:5)
                    pairdist=SQRT(DOT_PRODUCT(pairvec,pairvec))
                    IF(ABS(pairdist-axisdist)<matchthd) THEN
                       IF(nv/=0) THEN
                          theta=ACOS(DOT_PRODUCT(pairvec,orientvec)/pairdist&
                               &/SQRT(DOT_PRODUCT(orientvec,orientvec)))
                       ELSE
                          theta=0.d0
                       END IF
                       IF(theta<acute .OR. theta>obtuse) THEN
                          IF(DOT_PRODUCT(pairvec,orientvec)<0.d0) THEN
                             orientvec(:)=orientvec(:)-pairvec(:)
                          ELSE
                             orientvec(:)=orientvec(:)+pairvec(:)
                          END IF
                          nv=nv+1
                       END IF
                    END IF
                 END IF
              END DO
           END IF
        END DO
        orientvec(:)=orientvec(:)/nv
!
!print the orientation of the axis
!
        WRITE(*,'(1X,"The grain orientation is: ",3F7.2)') (orientvec(i),i=1,3)
!
        DEALLOCATE(posseq)
!
        END PROGRAM

!===============================================================================
        RECURSIVE SUBROUTINE quicksort(a, first, last, nlow, ncolum, icolum)
        IMPLICIT NONE
        INTEGER :: first, last, nlow, ncolum, icolum
        INTEGER :: i, j
        REAL(8) :: x
        REAL(8) :: a(nlow,ncolum)
        REAL(8) :: t(ncolum)

        x = a( (first+last) / 2, icolum )
        i = first
        j = last
        DO
           DO WHILE (a(i,icolum) < x)
              i=i+1
           END DO
           DO WHILE (x < a(j,icolum))
              j=j-1
           END DO
           IF (i >= j) EXIT
           t(:) = a(i,:);  a(i,:) = a(j,:);  a(j,:) = t(:)
           i=i+1
           j=j-1
       END DO
       IF (first < i-1) CALL quicksort(a, first, i-1, nlow, ncolum, icolum)
       IF (j+1 < last)  CALL quicksort(a, j+1, last, nlow, ncolum, icolum)
       END SUBROUTINE quicksort
