        MODULE DumpInfo
!===========================================================
! read the whole dump file. Only orthorgonal box is supported
! The atom are sorted using IDs from low to high in the arrays 
!                                                  H.B.LUO
!===========================================================
        IMPLICIT NONE
        TYPE :: frames
           CHARACTER(80) :: timestring
           CHARACTER(80) :: numstring
           CHARACTER(80) :: boxstring
           CHARACTER(80) :: atomstring
           CHARACTER(30),ALLOCATABLE :: field(:)    ! fields: use like IF(TRIM(field(i))=="id")
           INTEGER :: timestep
           INTEGER :: natom
           REAL(8) :: box(3,2)
           REAL(8),ALLOCATABLE :: atomInfo(:,:)
        END TYPE
        TYPE(frames),ALLOCATABLE :: frame(:)
        INTEGER :: nframe
!
        CONTAINS
!
        SUBROUTINE readDump(dumpfile)
        IMPLICIT NONE
        CHARACTER(50) :: dumpfile, string, ctemp
        INTEGER :: i,j,ierr,nf,iframe
        
        OPEN(8,file=TRIM(dumpfile),status='old')
!
        nframe=0
        DO WHILE(.TRUE.)
           READ(8,'(A50)',iostat=ierr) string
           IF(string(7:14)=="TIMESTEP") nframe=nframe+1
           IF(ierr/=0) EXIT
        END DO
        IF(ALLOCATED(frame)) DEALLOCATE(frame)
        ALLOCATE(frame(nframe))
        REWIND(8)
!
        iframe=1
        DO WHILE(iframe<=nframe)
           READ(8,'(A50)') frame(iframe)%timestring
           READ(8,*) frame(iframe)%timestep
           READ(8,'(A50)') frame(iframe)%numstring
           READ(8,*) frame(iframe)%natom
           READ(8,'(A50)') frame(iframe)%boxstring
           DO i=1,3
              READ(8,*) (frame(iframe)%box(i,j),j=1,2)
           END DO
           READ(8,'(A50)') frame(iframe)%atomstring
           nf=getFieldNum(frame(iframe)%atomstring)-2
           ALLOCATE(frame(iframe)%field(nf))
           READ(frame(iframe)%atomstring,*) &
               &ctemp,ctemp,(frame(iframe)%field(i),i=1,nf) 
!
           ALLOCATE(frame(iframe)%atomInfo(frame(iframe)%natom,nf))
           DO i=1,frame(iframe)%natom
              READ(8,*) (frame(iframe)%atomInfo(i,j),j=1,nf)
           END DO
           CALL quicksort(frame(iframe)%atomInfo(:,:),1,frame(iframe)%natom,&
                         &frame(iframe)%natom,nf,1)
           iframe=iframe+1
        END DO
        CLOSE(8)
!
        END SUBROUTINE
!
        INTEGER FUNCTION getFieldNum(string)
!==================================================
! get the field number in a string with space(s) as
! the delimiter                            H.B.LUO
!==================================================
        IMPLICIT NONE
        CHARACTER(*), INTENT(IN) :: string
        INTEGER :: i
        getFieldNum=0
        IF(LEN_TRIM(ADJUSTL(string)) > 0 .AND.&
          &LEN_TRIM(ADJUSTL(string)) <= 2) THEN
           getFieldNum=getFieldNum+1
        ELSE
           DO i=2, LEN_TRIM(ADJUSTL(string))
              IF(string(i:i)==" " .AND. &
                &string(i-1:i-1)/=" ") THEN
                 getFieldNum=getFieldNum+1
              END IF
           END DO
           getFieldNum=getFieldNum+1
        END IF
        END FUNCTION
!
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
!
        END MODULE
