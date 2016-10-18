! ##############################################################################################

! @@@@   MODULES   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ##############################################################################################

module Global
implicit none

integer :: nInd, nSnp, nIndLH, maxSibSize, MaxMismatch, maxOppHom, &
nRounds, nC(2), nAgeClasses, nPairs, BY1, maxAgePO, mxA, XP
integer, allocatable, dimension(:) :: Sex, BY, PairType, nFS
integer, allocatable, dimension(:,:) :: Genos, AgeDiff, Parent, OppHomM, &
    nS, PairID, FSID
integer, allocatable, dimension(:,:,:) :: SibID, GpID
double precision :: thLR, thLRrel, Er, OcA(3,3), AKA2P(3,3,3), OKA2P(3,3,3)
double precision, allocatable, dimension(:) ::  Lind, PairDLLR, AF
double precision, allocatable, dimension(:,:) :: AHWE, OHWE, LLR_O, LindX, &
    LR_parent, AgePriorM, CLL
double precision, allocatable, dimension(:,:,:) :: AKAP, OKAP, OKOP, LR_GP, &
  LindG, PHS, PFS, DumBY
double precision, allocatable, dimension(:,:,:,:) :: DumP
double precision, allocatable, dimension(:,:,:,:,:) :: XPr
character(len=2) :: DumPrefix(2)
character(len=200) :: ParentageFileName, PedigreeFileName, AgePriorFileName
character(len=30), allocatable, dimension(:) :: Id, NameLH
character(len=30), allocatable, dimension(:,:) :: ParentName

  contains
pure function MaxLL(V)
double precision, intent(IN) :: V(:)
double precision :: MaxLL
if (ANY(V < 0)) then
    MaxLL = MAXVAL(V, mask = (V<0 .and. V>-HUGE(0.0D0)), DIM=1)
else
    MaxLL = MINVAL(V, DIM=1)  ! 777: can't do; 888: already is; 999: not calc'd
endif
end function MaxLL

end module Global

! ##############################################################################################

! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing
! Made F conformant by Walt Brainerd

! Adapted by J Huisman (j.huisman@ed.ac.uk) to output rank, in order to enable 
! sorting of parallel vectors, and changed to decreasing rather than increasing order

module qsort_c_module
implicit none
public :: QsortC
private :: Partition

 contains
recursive subroutine QsortC(A, Rank)
  double precision, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: Rank
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq, Rank)
     call QsortC(A(:iq-1), Rank(:iq-1))
     call QsortC(A(iq:), Rank(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker, Rank)
  double precision, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: Rank
  integer, intent(out) :: marker
  integer :: i, j, TmpI
  double precision :: temp
  double precision :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1
  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp 
        
        TmpI = Rank(i) 
        Rank(i) = Rank(j)
        Rank(j) = TmpI
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module


! ##############################################################################################

! @@@@   PROGRAMS   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ##############################################################################################


subroutine duplicates(nDupGenoID, nDupLhID, nDupGenos, nSexless)
!program Duplicates
use Global
implicit none

integer, intent(OUT) :: nDupGenoID, nDupLhID, nDupGenos, nSexless   
integer :: i, j, l, Match, CountMismatch
integer, allocatable, dimension(:,:) :: dupGenoIDs, dupLhIDs, DupGenos
integer, allocatable, dimension(:) :: Sexless, nMisMatch

 call ReadData
 
nDupGenoID = 0
nDupLhID = 0
nDupGenos = 0
nSexless = 0
allocate(DupGenoIDs(nInd,2))
allocate(DupLhIDs(nIndLH,2))
allocate(DupGenos(nInd,2))
allocate(Sexless(nInd))
allocate(nMismatch(nInd))

! subroutine intpr (label, nchar, data, ndata)
! call intpr ("nInd Geno: ", 11, nInd, 1)
! call intpr ("nInd LH: ", 9, nIndLH, 1)
! any duplicate IDs?
do i=1,nInd-1
    do j=i+1, nInd
        if (Id(i) == Id(j)) then
            nDupGenoID = nDupGenoID + 1
            dupGenoIDs(nDupGenoID,1) = i
            dupGenoIDs(nDupGenoID,2) = j
        endif
    enddo
enddo
 
do i=1,nIndLH-1
    do j=i+1, nIndLH
        if (NameLH(i) == NameLH(j)) then
            nDupLhID = nDupLhID + 1
            DupLhIDs(nDupLhID,1) = i
            DupLhIDs(nDupLhID,2) = j
        endif
    enddo
enddo

! identical genotypes?
do i=1,nInd-1
    do j=i+1, nInd
        Match=1
        CountMismatch=0
        do l=1, nSnp
            if (Genos(l,i)==-9 .or. Genos(l,j)==-9) cycle
            if (Genos(l,i) /= Genos(l,j)) then
                CountMismatch=CountMismatch+1
                if (CountMismatch > MaxMismatch) then
                    Match=0
                    exit
                endif
            endif
        enddo
        if (Match==1) then
            nDupGenos = nDupGenos + 1
            DupGenos(nDupGenos,1) = i
            DupGenos(nDupGenos,2) = j
            nMisMatch(nDupGenos) = CountMismatch
        endif
    enddo
enddo

! check if all genotyped individuals have a gender in the lifehistory data
do i=1,nInd
    if (Sex(i)/=1 .and. Sex(i)/=2) then
        nSexless = nSexless + 1
        Sexless(nSexless) = i
    endif
enddo

! call intpr ( "dup: ",5, (/ nDupGenoID, nDupLhID, nDupGenos, nSexless /), 4)
! ##########################

open (unit=201,file="DuplicatesFound.txt",status="replace")
write (201, '(a15, a10, a30, a10, a30, a10)') "Type", "Row1", "ID1", "Row2", "ID2", "nDiffer"
if (nDupGenoID>0) then
    do i=1,nDupGenoID
        write (201,'(a15, i10, " ", a30, i10, " ",a30)') "GenoID", dupGenoIDs(i,1), Id(dupGenoIDs(i,1)), &
        dupGenoIDs(i,2), Id(dupGenoIDs(i,2))
    enddo   
endif
if (nDupLhID>0) then
    do i=1,nDupLhID
        write (201,'(a15, i10," ", a30, i10, " ",a30)') "LifehistID", DupLhIDs(i,1), NameLH(DupLhIDs(i,1)), &
        DupLhIDs(i,2), NameLH(DupLhIDs(i,2))
    enddo   
endif
if (nDupGenos>0) then
    do i=1,nDupGenos
        write (201,'(a15, i10, " ",a30, i10, " ",a30, i10)') "Genotype", dupGenos(i,1), Id(dupGenos(i,1)), &
        dupGenos(i,2), Id(dupGenos(i,2)), nMismatch(i)  
    enddo   
endif
if (nSexless>0) then
    do i=1, nSexless
        write (201,'(a15, i10," ", a30, i10, " ",a30)') "NoSex", Sexless(i), Id(Sexless(i)), 0, "NA"
    enddo
endif
 close (201)

if (allocated (DupGenoIDs)) deallocate(DupGenoIDs)
if (allocated (DupLhIDs))  deallocate(DupLhIDs)
if (allocated (DupGenos))  deallocate(DupGenos)
if (allocated (Sexless))   deallocate(Sexless)
if (allocated (nMismatch)) deallocate(nMismatch)

 call DeAllocAll

end subroutine duplicates
!end program Duplicates

! ##############################################################################################

subroutine parents(nParents, nAmbiguous)
!program AssignParents
use qsort_c_module
use Global
implicit none

integer, intent(INOUT) :: nParents(2), nAmbiguous
integer :: i, j, k, Round, isP(2), topX, maybe
integer, allocatable, dimension(:) :: BYRank
integer, allocatable, dimension(:,:) :: OppHomDF
double precision, allocatable, dimension(:) :: SortBY
double precision :: LLtmp(7,3), dLL, TotLL(42)
 character(len=2) :: RelName(8)
               
 call ReadData
 call PrecalcProbs

maxOppHom = MaxMismatch - FLOOR(-nSNP * Er)   ! round up to nearest integer  
allocate(OppHomM(nInd, nInd))
OppHomM = 999 
allocate(LLR_O(nInd, nInd))
LLR_O = 999
allocate(Parent(nInd,2)) 
Parent = 0
nAmbiguous = 0
nC = 0
!============================

call UpdateAllProbs
call intpr ( "Parentage ... ", -1, 0, 0)
call dblepr("Initial total LL : ", -1, SUM(Lind), 1) 

call CalcOppHom   ! also checks no. SNPs typed in both
 
do i=1, nInd-1
    do j=i+1,nInd 
        if (OppHomM(i,j) > maxOppHom) cycle    
        if (AgeDiff(i,j) /= 999) then
            if (AgeDiff(i,j) > 0) then  ! j older than i
                call CalcPO(i, j, LLR_O(i,j))  ! LLR PO/U
            else if (AgeDiff(i,j) < 0) then
                call CalcPO(j, i, LLR_O(j,i))
            endif
        else
            call CalcPO(i, j, LLR_O(i,j))
            call CalcPO(j, i, LLR_O(j,i))
        endif  
    enddo
enddo

! get birthyear ranking (increasing)
allocate(SortBY(nInd))
allocate(BYRank(nInd))
SortBY = REAL(BY, 8)
WHERE (SortBY < 0) SortBY = HUGE(0.0D0) 
BYRank = (/ (i, i=1, nInd, 1) /)
call QsortC(SortBy, BYRank)
 
! loop through, assign parents where possible. start with oldest.
allocate(LR_parent(nInd,3))
LR_parent = 999

TotLL = 0
TotLL(1) = SUM(Lind)
do Round=1,42  
    call Parentage(BYrank)
     
    call UpdateAllProbs
    
    ! assign sex if assigned as parent >= 2 times
    do i=1,nInd
        if (Sex(i)==3) then
            isP = 0
            do k=1,2
                do j=1,nInd
                    if (Parent(j,k) == i) then
                        isP(k) = isP(k) + 1
                    endif
                enddo
            enddo
            if (isP(1)>0 .and. isP(2)>0) then
                call rwarn("Assigned as both dam & sire: "//ID(i))
            else
                do k=1,2
                    if (isP(k)>1) then
                        Sex(i) = k
                    endif
                enddo
            endif
        endif
    enddo
    call UpdateAllProbs
    TotLL(Round + 1) = SUM(Lind)
    if (TotLL(Round + 1) - TotLL(Round) < thLR) exit  ! run til convergence
enddo

 call UpdateAllProbs
 call dblepr("Post-parentage total LL : ", -1, SUM(Lind), 1) 
 call CalcParentLLR

allocate(OppHomDF(nInd,2))
OppHomDF = -9
nParents = 0
do i=1,nInd
    do k=1,2
        if (Parent(i,k)>0) then
            ParentName(i,k) = Id(Parent(i,k))
            nParents(k) = nParents(k) + 1
            OppHomDF(i,k) = OppHomM(i, Parent(i,k))  
        endif
    enddo
enddo
    
open (unit=201,file=trim(ParentageFileName), status="unknown")  ! "Parents_assigned.txt"
write (201, '(3a30, 5a12, 3a6)') "ID", "Dam", "Sire", "LLR_dam", "LLR_sire", "LLR_pair", & 
"OH_Mother", "OH_Father", "RowO", "RowD", "RowS"   
do i=1,nInd
    write (201,'(3a30, 3f10.2, 2i12, 3i6)') Id(i), ParentName(i,1:2), LR_parent(i,1:3), &
    OppHomDF(i,1:2), i, Parent(i,1:2)
enddo   
 close (201)
! last 3 columns (9-11) read in for sibship assignment - update subroutine ReadParents when changing number of columns!!

! output non-assigned likely parents-offspring pairs. (unknown sex or birth year etc.) 
RelName = (/ "PO", "FS", "HS", "GP", "FA", "HA", "U ", "XX" /) 
open (unit=201,file="Unassigned_relatives_par.txt",status="unknown")
write (201, '(2a30, 2a5, 4a10, a5)') "ID1", "ID2", "Sex1", "Sex2", "AgeDif", "BestRel", "LLR_R_U", "LLR_R1_R2", "OH"
do i=1,nInd-1  ! TODO: in separate subroutine
    do j=i+1,nInd
        if (OppHomM(i,j) > maxOppHom .or. OppHomM(i,j)<0) cycle
        if (LLR_O(i,j)==999 .or. LLR_O(i,j)< 0) cycle
        if (ANY(Parent(i,:)==j) .or. ANY(Parent(j,:)==i)) cycle
        if (Parent(i,1)/=0 .and. Parent(i,2)/=0 .and. Parent(i,1)==Parent(j,1) .and. &
        Parent(i,2)==Parent(j,2)) cycle  ! full sibs.  
        if (Sex(j)/=3) then
            call CalcPair(i, j, Sex(j), .FALSE., LLtmp(:,1), 1)
        else
            call CalcPair(i, j, 1, .FALSE., LLtmp(:,1), 1)
        endif
        if (Sex(i)/=3) then
            call CalcPair(j, i, Sex(i), .FALSE., LLtmp(:,2), 1)
        else
            call CalcPair(j, i, 1, .FALSE., LLtmp(:,2), 1)
        endif
        do k=1,7
            LLtmp(k,3) = MaxLL(LLtmp(k,1:2)) 
        enddo
        call BestRel(LLtmp(:,3), 1, topX, dLL)
        maybe = 1
        if (topX>3) then
            do k=1,2
                if (Parent(i,k)>0) then
                    if (ANY(Parent(Parent(i,k),:)==j))  maybe = 0  ! GP
                    if (Parent(i,k) == Parent(j,k))  maybe = 0
                endif
                if (Parent(j,k)>0) then
                    if (ANY(Parent(Parent(j,k),:)==i))  maybe = 0
                endif
            enddo
        endif
        if (maybe == 1) then
            nAmbiguous = nAmbiguous + 1
            write (201,'(2a30, 2i5, i10, a10, 2f10.2, i5)') Id(i), ID(j), Sex(i), Sex(j), &
              AgeDiff(i,j), RelName(TopX), LLR_O(i,j), dLL, OppHomM(i,j)
        endif
    enddo
enddo   
 close (201)

deallocate(BYRank)
deallocate(OppHomDF)
deallocate(SortBY)
 call DeAllocAll
 
end subroutine parents
!end program AssignParents

! ##############################################################################################

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ##############################################################################################

subroutine sibships
!program Sibships
use Global
implicit none

integer :: Round, i, j, s, k, n, RX, LastR, topX
 character(len=1) :: RoundC
 character(len=2) :: RelName(8)
 character(len=4) :: DumTmp
 character(len=30), allocatable, dimension(:,:) :: DumName
 character(len=30), allocatable, dimension(:,:,:) :: GpName
double precision :: TotLL(42), dLL, tmpLL(7), LRS

RX = 1  ! no. of initial rounds, pairs-cluster-merge only
    
 call ReadData
 call ReadParents
 call PrecalcProbs
 
XP = 5
allocate(PairID(XP*nInd, 2)) 
allocate(PairDLLR(XP*nInd))
allocate(PairType(XP*nInd))  ! maternal (1), paternal (2) or unknown (3)
allocate(CLL(nInd/2,2))
 CLL = 999
allocate(nS(nInd/2,2))
nS = 0
allocate(SibID(maxSibSize, nInd/2, 2))
SibID = 0
nC = 0
allocate(GpName(2,nInd/2,2))
allocate(DumName(nInd/2, 2))
allocate(LR_parent(nInd,3))
LR_parent = 999
allocate(LR_GP(3, nInd/2,2))
LR_GP = 999
allocate(DumBY(nAgeClasses +maxAgePO, nInd/2, 2)) 
!=======================

!open (unit=222,file="Sequoia.log",status="unknown")
!write(222, *) "== Sibship clustering == "

! find current FS (based on real parents)
do i=1,nInd-1
    do j=i+1,nInd
        if (Parent(i,1)==Parent(j,1) .and. Parent(i,1)/=0 .and. &
            Parent(i,2)==Parent(j,2) .and. Parent(i,2)/=0) then
            call MakeFS(i, j)
        endif
    enddo
enddo

call UpdateAllProbs
! call dblepr("Sibships - Initial Total LL : ", -1, SUM(Lind), 1)
call intpr ( "Sibships ... ", -1, 0, 0)
! write(222, '("Total LL IN: ", f12.1)') SUM(Lind) 
 
LastR = 0
do Round=1, Nrounds
    call UpdateAllProbs
    TotLL(Round) = SUM(Lind)
    write(RoundC, '(i1)') Round  ! for filenames
    if(Round > 1) then
        if (ABS(TotLL(Round-1) - TotLL(Round)) < thLR) then
            LastR = 1
        endif
    endif
    if (Round==1 .and. nAgeClasses>1) then
        call FindPairs(.FALSE.)   ! do not use age priors yet to diff. HS/GP/FA
    else
        call FindPairs(.TRUE.)
    endif
!    call intpr ( "Found pairs: ",-1, nPairs, 1)
!    write(222, '("Found Pairs: ", i5)') nPairs
!    open(unit=405, file="PairwiseLLR_"//RoundC//".txt", status="replace")
!    write(405, '(2a15, 2a10, a15)') "ID1", "ID2", "PatMat", "dLLR"
!    do i=1,nPairs
!        write(405, '(2a15, i10, f15.4)') ID(PairID(i,1:2)), PairType(i),  PairDLLR(i)
!    enddo
!    close(405)
   
    if (Round==1) then
        call Clustering(-1)
    else
        call Clustering(LastR)
    endif
    call UpdateAllProbs
!    write(222, '("Found clusters: ", 2i5)') nC
!   call intpr ( "Clusters: ",-1, nC, 2)
    
    call Merging
    call UpdateAllProbs
    
!    write(222, *) "expand clusters ..."
    if (Round > RX .or. Round==Nrounds .or. LastR==1) then
        call GrowClusters
        call UpdateAllProbs
    endif
   
 !   open (unit=401,file="Sibships_"//RoundC//".txt",status="replace")
 !   write (401,'(2a10, a5)')  "PatMat", "Sibship", "ID"
 !   do k=1,2
 !       do s=1,nC(k)
 !           do n=1, nS(s,k)
 !               write (401,'(2i10, " ", a30)') k, s, Id(SibID(n,s,k))
 !           enddo
 !       enddo
 !   enddo
 !   close (401)
 
    if (nAgeClasses > 1 .and. (Round > RX .or. Round==Nrounds .or. LastR==1)) then 
!        write(222, *) "add parents ..."
        call SibParent  ! replace dummy parents by indivs
        call UpdateAllProbs
        
 !       write(222, *) "parental pairs ..."
        call MoreParent  !  assign additional parents to singletons (e.g. unknown BY/sex)
        call UpdateAllProbs      
        

        if (Round > RX+1) then
  !     write(222, *) "add grandparents ..."
            call GGpairs
            call UpdateAllProbs
        endif
    endif
    
!    write(222, *) "add grandparents ..."
    call SibGrandparents
    call UpdateAllProbs
    
    !=======================
    call dblepr("Round "//RoundC//", Total LL : ", -1, SUM(Lind), 1) 
!    write(222, '("Round ", i5, ", Total LL: ", f12.1)') Round, SUM(Lind) 
    
    ! give names to dummy parents
    do k=1,2
        do s=1, nC(k)
            write(DumTmp, '(i4.4)') s
            DumName(s,k) = trim(DumPrefix(k))//DumTmp
        enddo
    enddo    
    
    ParentName = "NA"
    do i=1,nInd
        do k=1,2
            if (Parent(i,k)>0) then
                ParentName(i,k) = Id(Parent(i,k))
            else if (Parent(i,k)<0) then
                ParentName(i,k) = DumName(-Parent(i,k),k)
            endif
        enddo
    enddo
    
    GpName = "NA"
    do k=1,2
        do s=1,nC(k)
            do n=1,2
                if (GpID(n,s,k)>0) then
                    GpName(n,s,k) = Id(GpID(n,s,k))
                else if (GpID(n,s,k)<0) then
                    GpName(n,s,k) = DumName(-GpID(n,s,k), n)
                endif
            enddo
        enddo
    enddo
        
    if (Round == nRounds .or. LastR==1) then
        call UpdateAllProbs
        call CalcParentLLR    
!        close(222)
        open (unit=407,file=trim(PedigreeFileName),status="unknown")
        write (407,'(3a30, 3a10)')  "ID", "Dam", "Sire", "LLR_dam", "LLR_sire", "LLR_pair"
        do i=1,nInd
            write (407,'(3a30, 3f10.2)') Id(i), ParentName(i,1:2), LR_parent(i,:)
        enddo
        do k=1,2
            do s=1,nC(k)
                write (407,'(3a30, 3f10.2)') DumName(s,k), GpName(:, s,k), LR_GP(:, s, k)
            enddo
        enddo
        close (407)
        
        ! output non-assigned likely 1st/2nd degree relatives pairs. 
        RelName = (/ "PO", "FS", "HS", "GP", "FA", "HA", "U ", "XX" /) 
        open (unit=201, file="Unassigned_relatives_sib.txt", status="unknown")
        write (201, '(2a15, 2a5, 4a10)') "ID1", "ID2", "Sex1", "Sex2", "AgeD", "BestRel", "LLR_R_U", "LLR_R1_R2"
        do i=1,  nInd-1
            do j=i+1,nInd
                do k=1,2
                    if (k==2 .and. topX==4 .and. Parent(i,1)==0 .and. Parent(i,2)==0) cycle
                    if (Parent(i,k)/=0 .and. Parent(j,k)/=0) cycle  
                    ! assume assigned 1st/2nd degree relatives are not futher related.
                    if (ANY(Parent(i,:)==j) .or. ANY(Parent(j,:)==i))  cycle  ! PO
                    if (Parent(i,k)/=0 .and. Parent(i,k)==Parent(j,k)) cycle
                    if (Parent(i,k) < 0) then
                        if (ANY(GpID(:,-Parent(i,k),k)==j))  cycle
                    endif
                    if (Parent(j,k) < 0) then
                        if (ANY(GpID(:,-Parent(j,k),k)==i))  cycle  
                    endif
                    if (k==2 .and. Parent(i,1)==0 .and. Parent(j,1)==0) cycle ! already done
                    call PairQHS(i, j, LRS)  ! quick check
                    if (LRS < thLRrel) cycle           
                    if (AgeDiff(i,j)==999 .or. AgeDiff(i,j)>=0) then
                        call CalcPair(i, j, k, .FALSE., tmpLL, 4)  ! 4: gp -> return all
                    else
                        call CalcPair(j, i, k, .FALSE., tmpLL, 4)
                    endif
                    call BestRel(tmpLL, 3, topX, dLL)
                    if (topX>4) cycle
                    write (201,'(2a30, 2i5, i10, a10, 2f10.2)') Id(i), ID(j), Sex(i), Sex(j), &
                      AgeDiff(i,j), RelName(TopX), LRS, dLL  
                enddo
            enddo
        enddo 
        close (201) 
        
!        call WriteDummies  ! slow! 
        
        call intpr ( "Finished. ", -1, 0, 0)  
        exit
!    else
!        open (unit=407,file="Pedigree_"//RoundC//".txt",status="replace")
!        write (407,'(3a20)')  "ID", "Dam", "Sire"
!        do i=1,nInd
!            write (407,'(3a20,3f15.4)') Id(i), ParentName(i,1:2)
!        enddo
!        do k=1,2
!            do s=1,nC(k)
!                write (407,'(3a20,3f15.4)') DumName(s,k), GpName(:, s,k)
!            enddo
!        enddo
!        close (407)       
    endif
enddo

deallocate(DumName)
deallocate(GpName)
 call DeAllocAll
    
end subroutine Sibships
!end program Sibships

! ##############################################################################################

! @@@@   SUBROUTINES   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ##############################################################################################

subroutine CalcOppHom  ! nInd x nInd matrix with no. opp. hom. loci
use Global
implicit none

integer :: i, j, Lboth

OppHomM = -999
do i=1, nInd-1
    do j=i+1,nInd
        call CalcOH(i, j, OppHomM(i,j))
        OppHomM(j,i) = OppHomM(i,j)
        if (OppHomM(i,j) <= maxOppHom) then
            Lboth = COUNT(Genos(:,i)/=-9 .and. Genos(:,j)/=-9)    
            if (Lboth < nSnp/4.0) then   ! more than 3/4th of markers missing
                OppHomM(i,j) = -Lboth 
                OppHomM(j,i) = -Lboth
            endif
        endif
    enddo
enddo

end subroutine CalcOppHom

! ##############################################################################################

subroutine CalcOH(A,B,OH)
use Global
implicit none

integer, intent(IN) :: A, B
integer, intent(OUT) :: OH
integer :: l

OH = 0
do l=1,nSnp
    if ((Genos(l,A)==1).and.(Genos(l,B)==3)) then
        OH = OH+1
        if (OH > maxOppHom) exit
    endif                       
    if ((Genos(l,A)==3).and.(Genos(l,B)==1)) then
        OH = OH+1
        if (OH > maxOppHom) exit
    endif                       
enddo

end subroutine CalcOH
! ##############################################################################################

subroutine CalcPO(A,B, LLR)  ! LLR of A as offspring from B, vs A as random sample from pop
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LLR
integer :: l
double precision :: PrL(nSnp)

LLR = 999
PrL = 0

do l=1,nSnp
    if (Genos(l,A)/=-9 .and. Genos(l,B)/=-9) then
        PrL(l) = LOG10(OKOP(Genos(l,A), Genos(l,B), l))
    else if (Genos(l,A)/=-9) then
        PrL(l) = LOG10(OHWE(Genos(l,A),l))
    endif
enddo

LLR = SUM(PrL) - Lind(A)

end subroutine CalcPO

! ########################################################################################

subroutine Parentage(BYrank)
use Global
implicit none

integer, intent(IN) :: BYrank(nInd)  ! TODO: make global
integer :: i, j, x, y, k, CandPar(50, 2), nCP(2), u, v, AncJ(2,mxA), offspr
logical :: skip

do x=1, nInd
!    if (MOD(x,500)==0) then
!          call intpr ( " ",2, x, 1)
!   endif
    i = BYRank(x)
    nCP = 0
    CandPar = 0
    do y=1,nInd 
        j = BYRank(y)
        if (i==j) cycle
        if (ANY(Parent(j,:)==i) .or. ANY(Parent(i,:)==j)) cycle
        if (OppHomM(i,j) > maxOppHom .or. OppHomM(i,j)<0) cycle 
        if (LLR_O(i,j) < -thLR) cycle  
        if (AgeDiff(i,j) == 999) then
            if (ALL(Parent(i,:)==0) .and. ALL(Parent(j,:)==0)) cycle
            call GetAncest(j,1,AncJ)
            if (ANY(AncJ == i)) cycle
            skip = .FALSE.
            if (BY(i)/=-999 .and. BY(j)==-999) then
                do k=1,2
                    do u=2, mxA
                        if (AncJ(k,u)>0) then
                            if (AgeDiff(i, AncJ(k,u)) <= 0)  skip = .TRUE.
                        endif
                    enddo
                enddo
            endif
            if (skip) cycle
            if (BY(i)==-999 .and. BY(j)/=-999 .and. Sex(i)/=3) then
                if (ANY(Parent(:, Sex(i))==i)) then
                    offspr = MINLOC(Parent(:, Sex(i)) - i, DIM=1)  ! note: 1st offspring only
                    if (AgeDiff(offspr, j) <= 0)  cycle
                endif
            endif
        else 
            if (AgeDiff(i,j) <= 0)  cycle
            if (Sex(j)==3 .and. Parent(i,1)==0 .and. Parent(i,2)==0) cycle 
        endif    
        do k=1,2
            if (Sex(j) /= 3 .and. Sex(j)/= k) cycle
            if (nCP(k)==50) cycle
            nCP(k) = nCP(k) + 1
            CandPar(nCP(k), k) = j
        enddo
    enddo
    
    if (ALL(nCP>0) .and. ANY(nCP>1)) then   ! test all combo's
        do u = 1, nCP(1)
            do v = 1, nCP(2)
                if (CandPar(u, 1)>0 .and. CandPar(u, 1)==CandPar(v, 2))  cycle  ! when unknown sex.
                if (Parent(i, 1)/= CandPar(u,1))  call calcPOZ(i, CandPar(u, 1), .FALSE.) 
                if (Parent(i, 2)/= CandPar(v, 2))  call calcPOZ(i, CandPar(v, 2), .FALSE.)  
            enddo
        enddo
    else if (ANY(nCP>1)) then  ! & nCP(3-k)=0
        do k=1,2
            if (nCP(k)>1) then
                do v = 1, nCP(k)
                    call calcPOZ(i, CandPar(v, k), .FALSE.)
                enddo
            endif
        enddo
    else if (CandPar(1, 1)>0 .and. CandPar(1, 1)==CandPar(1, 2)) then
        cycle
    else if (ALL(nCP==1)) then
        do k=1,2
            call calcPOZ(i, CandPar(1, k), .FALSE.)
        enddo
    else
        do k=1,2
            if (nCP(k)>0) then
                call CalcPOZ(i, CandPar(1, k), .FALSE.)
            endif
        enddo
    endif
enddo

end subroutine Parentage

! ########################################################################################
subroutine CalcPOZ(A, B, UseAge)  ! replace a current parent of A by B? k=sex(B)
use Global
implicit none 

integer, intent(IN) :: A, B
logical, intent(IN) :: UseAge
integer :: m, CurPar(2), TopX, k, x, CY(3), kY(3), y
double precision :: LLA(2,7,7), TopLL, LLcp(3), LLBA(7), TopBA, dLL, LLtmp(2), ALR(3)
 
 CurPar = Parent(A,:)

if (AgeDiff(A,B)==999 .and. (ALL(Parent(A,:)==0) .and. ALL(Parent(B,:)==0))) then
    return  ! can't tell if A or B is parent
endif
if (Sex(B)/=3) then
    k = Sex(B) 
else if (ALL(Parent(A,:)==0) .or. ALL(Parent(A,:)/=0)) then
    return  ! can't tell if B is father or mother
else 
    do m=1,2
        if (Parent(A,m)==0) then
            k = m  ! try B as opposite-sex parent
        endif
    enddo   
endif
LLtmp = 999
if (Sex(B)==3 .and. CurPar(3-k)<0) then
    call AddParent(B, -CurPar(3-k),3-k, LLtmp(1))      
    call CalcU(B,3-k, CurPar(3-k),3-k, LLtmp(2))
    if ((LLtmp(1) - LLtmp(2)) > -thLRrel) then
        return  ! B likely to replace CurPar(3-k) instead
    endif
endif
if (AgeDiff(A,B)==999) then  
    do m=1,2
        if ((Sex(A)==m .or. Sex(A)==3) .and. Parent(B,m) < 0) then
            call AddParent(A, -Parent(B,m),m, LLtmp(1)) 
            call CalcU(A,m, Parent(B,m),m, LLtmp(2))
            if (LLtmp(1)<0 .and. (LLtmp(1) - LLtmp(2) > -thLRrel)) then 
                return    ! A likely to replace Parent(B,m) instead
            endif
        endif
    enddo
endif

if (Parent(A,3-k) < 0) call CalcCLL(-Parent(A,3-k), 3-k)
 call CalcLind(A)
 call CalcLind(B)

TopX = 0
LLA = 999
LLBA = 999
if (Parent(A,1)==0 .and. Parent(A,2)==0) then   
    call CalcPair(A, B, k, UseAge, LLA(1,:,7), 1)   
    call BestRel(LLA(1,:,7), 1, TopX, dLL)
    TopLL = MaxLL(LLA(1,:,7))
    if (AgeDiff(A,B)==999) then
        if (Sex(A)/=3) then
            m = Sex(A)
        else if (Parent(B,2)==0) then
            m = 2
        else
            m = 1
        endif
        if (Parent(B,m)==0) then
            call CalcPair(B, A, m, UseAge, LLBA, 1)
            if (Parent(B,3-m) < 0) then  ! include changes in its CLL
                call CalcU(B,m, Parent(B,3-m),3-m, LLtmp(1))
                Parent(B,m) = A
                call CalcU(B,m, Parent(B,3-m),3-m, LLtmp(2))
                Parent(B,m) = 0
                LLBA(1) = LLBA(1) + (LLtmp(2) - LLtmp(1))
            endif
            TopBA = MaxLL(LLBA)  
            if (ABS(TopBA - TopLL) < thLRrel .or. (TopBA > TopLL))  return
        endif
    endif
    if (TopX==1) then
        Parent(A,k) = B
    endif
else 
    LLcp = 0  ! need LL over all (2-4) indiv involved
    TopLL = 999
    call CalcU(CurPar(1), 1, CurPar(2), 2, LLcp(3))
    do m=1,2
        call CalcU(CurPar(3-m), 3-m, B, k, LLcp(m))
        if (curPar(3-m) < 0) then
            LLCP(m) = LLCP(m) - Lind(A)
            LLCP(3) = LLCP(3) - Lind(A)  ! never called when both parents <0
        endif
    enddo
    ! age based prior prob.
    ALR = 0
    if (UseAge .and. BY(A)/=-999) then
        CY = (/ CurPar(1), CurPar(2), B /)
        kY = (/ 1, 2, k /)
        do y=1,3
            if (CY(y) > 0) then
                if (BY(CY(y))/=-999) then
                    ALR(y) = LOG10(AgePriorM(ABS(AgeDiff(A, CY(y)))+1, 6+k))   
                endif
            else if (CY(y) < 0) then
                if (y<3)  call RemoveSib(A, -CY(y), kY(y))
                call CalcAgeLR(A,0, CY(y),kY(y), 0,1, ALR(y))
                if (y<3)  call DoAdd(A, -CY(y), kY(y))
            endif
        enddo
    endif
       
    do m=1,2
        if (CurPar(m)==0) cycle  ! LLA(m,:,:) empty if only 1 CurPar
        call CalcPair(A, B, k, UseAge, LLA(m,:,1), 1)   ! CurPar(m)=Par + B_7        
        Parent(A,m) = 0 
        if (CurPar(m) < 0)  call RemoveSib(A, -CurPar(m), m)
        
        call CalcPair(A, B, k, UseAge, LLA(m,:,7), 1)   ! B_7
        if (CurPar(m)>0) then
            call CalcPair(A, CurPar(m), m, UseAge, LLA(m,7,:), 1)  ! CurPar(m)_7  
            Parent(A,k) = B
            call CalcLind(A)
            call CalcPair(A, CurPar(m), m, UseAge, LLA(m,1,:), 1)  ! B=Par, CurPar(m)_7
        else if (CurPar(m)<0) then
            call checkAdd(A, -CurPar(m), m, LLA(m,7,:), 4)  ! CurPar(m)_7
            LLA(m,7,2) = 333   ! FS does not count here.
            call ReOrderAdd(LLA(m,7,:))
            Parent(A,k) = B
            if (Parent(A,k)<0)  call RemoveSib(A, -CurPar(k), k)
            call checkAdd(A, -CurPar(m), m, LLA(m,1,:), 4)  ! B=Par, CurPar(m)_7 
            call ReOrderAdd(LLA(m,1,:))         
        endif
        if (Parent(B, m)==CurPar(m) .and.  CurPar(m)/= 0) then  ! HS implies curPar = Par
            LLA(m,3,7) = 888
            if (Parent(B, 3-m)==CurPar(3-m) .and.  CurPar(3-m)/= 0) then
                LLA(m,2,7) = 888
            endif
        endif
        
        Parent(A,:) = CurPar  ! restore
        if (CurPar(m) < 0) call DoAdd(A, -CurPar(m), m)
        if (m/=k .and. CurPar(k)<0)  call DoAdd(A, -CurPar(k), k)
        call CalcLind(A)
        WHERE (LLA(m,2:6,1)<0) LLA(m,2:6,1) = LLA(m,2:6,1) + LLcp(3) + ALR(1) + ALR(2)
        WHERE (LLA(m,2:6,7)<0) LLA(m,2:6,7) = LLA(m,2:6,7) + LLcp(3) + ALR(3-m)
        WHERE (LLA(m,1,:)<0) LLA(m,1,:) = LLA(m,1,:) + LLcp(m) + ALR(3)
        if (m==k) LLA(m,1,:) = LLA(m,1,:) + ALR(3-m)
        WHERE (LLA(m,7,:)<0) LLA(m,7,:) = LLA(m,7,:) + LLcp(m) + ALR(3-m) 
    enddo  
    
    TopLL = MaxLL(RESHAPE(LLA(:,:,:), (/2*7*7/))) ! MAXLOC doesn't cope well with ties

    if (AgeDiff(A,B)==999) then
        if (Sex(A)/=3) then
            call CalcPair(B, A, Sex(A), .FALSE., LLBA, 1)
        else
            if (Parent(B,1)/=0 .and. Parent(B,2)/=0) return
            do m=1,2
                if (Parent(B,m)==0) then
                    call CalcPair(B, A, m, .FALSE., LLBA, 1)  ! try A as opposite-sex parent
                endif
            enddo   
        endif
        WHERE (LLBA<0) LLBA = LLBA + LLcp(3)
        TopBA = MaxLL(LLBA)
        if (ABS(TopBA - TopLL) < thLRrel .or. (TopBA > TopLL)) then
            return
        endif
    endif
    
    if (LLA(3-k,1,1)==TopLL .or. ANY(LLA(k,1,:) == TopLL)) then  ! B + CurPar(3-k) 
        Parent(A, k) = B
        Parent(A, 3-k) = CurPar(3-k)
    else if ((TopLL - MaxLL(RESHAPE(LLA(:,2:7,1), (/2*6/)))) < 0.01) then  ! keep CurPar (both)
        Parent(A,:) = CurPar
    else if (ANY(LLA(:,1,:)==TopLL)) then  ! only B
        call BestRel(LLA(k,:,1),1,topX, TopBA)
        if ((MaxLL(LLA(k,1,:)) - TopBA) > thLRrel) then
            Parent(A,k) = B
        else
            Parent(A,k) = 0  ! unclear.
        endif
        Parent(A,3-k) = 0
    else if (ANY(LLA(k, 2:7, 2:7) == TopLL)) then  ! keep CurPar(3-k)
        Parent(A, k) = 0
    else if (ANY(LLA(3-k, 2:7, 2:7) == TopLL)) then  ! keep CurPar(k)
        Parent(A, k) = CurPar(k)
        Parent(A, 3-k) = 0
    else
        Parent(A,:) = CurPar
    endif
endif 

do m = 1, 2
    if (CurPar(m)<0 .and. Parent(A,m)==0) then   ! remove A from sibship.
        call RemoveSib(A, -CurPar(m), m)
        if (nS(-CurPar(m), m)==1 .and. ALL(GpID(:, -CurPar(m), m)==0)) then  ! single sib left; remove sibship
            y = SibID(1,-CurPar(m), m)
            call RemoveSib(y, -CurPar(m), m)
            call DoMerge(0, -CurPar(m), m)
        endif
    endif
enddo
if (Parent(A,3-k) < 0) then
    if (parent(A,k)/=0 .and. curPar(k)/=parent(A,k)) then  ! check for FS
        do y=1, nS(-Parent(A,3-k),3-k)  
            x = SibID(y,-Parent(A,3-k),3-k)
            if (x==A) cycle
            if (Parent(x, k)==Parent(A, k) .and. nFS(x)/=0) then
                call MakeFS(A, x)
            endif
        enddo
    endif
    call CalcCLL(-Parent(A,3-k), 3-k)
endif
call CalcLind(A)

end subroutine CalcPOZ

! ##############################################################################################

subroutine ReOrderAdd(LL)  ! reorder output from CheckAdd for compatibility with CalcPair (for POZ)
use Global
implicit none

double precision, intent(INOUT) :: LL(7)
double precision :: LLtmp(7)

LLtmp = 999
if (LL(4) < 0 .and. LL(4) - MaxLL(LL(2:3)) > -thLRrel) then  ! more likely to be parent of dummy than offspring / unclear
    LLtmp(1) = 222
else
    LLtmp(1) = MaxLL(LL(2:3))
endif
LLtmp(2:3) = LL(5:6)
LLtmp(7) = LL(7) 

LL = LLtmp

end subroutine ReOrderAdd

! ##############################################################################################

subroutine FindPairs(UseAgePrior)
use Global
use qsort_c_module
implicit none

logical, intent(IN) :: UseAgePrior
integer :: k, i, j, top, PairTypeTmp(XP*nInd), PairIDtmp(XP*nInd,2)
double precision :: dLL, PairLLRtmp(XP*nInd), tmpLL(7), LRS(2)
integer, allocatable, dimension(:) :: Rank
double precision, allocatable, dimension(:) :: SortBy

nPairs = 0
PairID = -9
PairDLLR = 999
PairType = 0

do i=1,  nInd-1
!    if (MODULO(i,500)==0) then 
!        call intpr ( " ",1, i, 1)
!    endif
    if (Parent(i,1)/=0 .and. Parent(i,2)/=0) cycle
    do j=i+1,nInd
        do k=1,2
            LRS = 0
            if (Parent(i,k)/=0 .or. Parent(j,k)/=0) cycle
            if (k==2 .and. Parent(i,1)==0 .and. Parent(j,1)==0) cycle ! already done
            call PairQHS(i, j, LRS(1))  ! quick check
            if (LRS(1) < -thLR) cycle  ! true FS more likely to be HS than U
            if (ALL(Parent(i,:)==0) .and. ALL(Parent(j,:)==0)) then
                call PairQFS(i, j, LRS(2)) 
                if (LRS(2) < -thLR) cycle   ! can't distinguish HS from FA and GP
            endif
            if (AgeDiff(i,j)==999 .or. AgeDiff(i,j)>=0) then
                call CalcPair(i, j, k, UseAgePrior, tmpLL, 3)
            else
                call CalcPair(j, i, k, UseAgePrior, tmpLL, 3)
            endif
            call BestRel(tmpLL, 3, top, dLL)
            if (top==2 .or. top==3) then
                if (nPairs >= XP*nInd) cycle  ! do in next round
                nPairs = nPairs+1
                PairID(nPairs, :) = (/ i, j /)
                PairDLLR(nPairs) = dLL
                if (k==1 .and. Parent(i,2)==0 .and. Parent(j,2)==0) then
                    pairType(nPairs) = 3  
                else
                    PairType(nPairs) = k
                endif
            endif
        enddo
    enddo
enddo
 
 ! sort by decreasing dLL
PairIDtmp = 0
PairLLRtmp = 0
allocate(Rank(nPairs))
allocate(SortBy(nPairs))
Rank = (/ (i, i=1, nPairs, 1) /)
SortBy = PairDLLR(1:nPairs)
 
 call QsortC(SortBy, Rank(1:nPairs))
do i=1,nPairs
    PairTypeTmp(i) = PairType(Rank(nPairs-i+1))  ! decreasing order
    PairIDtmp(i,1:2) = PairID(Rank(nPairs-i+1), 1:2)  
    PairLLRtmp(i) = PairDLLR(Rank(nPairs-i+1)) 
enddo 

PairType = PairTypeTmp
PairID = PairIDtmp 
PairDLLR = PairLLRtmp
deallocate(Rank)
deallocate(SortBy)

end subroutine FindPairs

! ##############################################################################################

subroutine PairQHS(A, B, LR)  ! quick check, not conditioning on parents.
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LR
integer :: l
double precision :: PrL(nSnp)

PrL = 0
do l=1,nSnp
    if (Genos(l,A)==-9 .or. Genos(l,B)==-9) cycle
    PrL(l) = LOG10(PHS(Genos(l,A), Genos(l,B), l))   ! note: >0 for FS pairs
enddo
LR = SUM(PrL)

end subroutine PairQHS

! ##############################################################################################

subroutine PairQFS(A, B, LR)  ! quick check, not conditioning on parents.
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LR
integer :: l
double precision :: PrL(nSnp)

PrL = 0
do l=1,nSnp
    if (Genos(l,A)==-9 .or. Genos(l,B)==-9) cycle
    PrL(l) = LOG10(PFS(Genos(l,A), Genos(l,B), l))  
enddo
LR = SUM(PrL)

end subroutine PairQFS

! ##############################################################################################

subroutine CalcPair(A, B, k, InclAge, LL, focal)  ! joined LL A,B under each hypothesis
use Global
implicit none
! to consider: self, full 1st cousins, GGP  (+ double 1st cousins?)

integer, intent(IN) :: A,B,k, focal
logical, intent(IN) :: InclAge  ! include age prior y/n (always checks if age compatible)
double precision, intent(OUT) :: LL(7)  ! PO,FS,HS,GG,FAU,HAU,U
integer :: x, cgp
double precision :: LLg(7), LLtmpA(2,3), LLtmpGGP, LLCC, LRS, ALR, LLX(4), LLZ(2), LLC(7)

LLg = 999
LL = 999
LLtmpGGP = 999
LRS = 999

if (AgeDiff(A, B)/=999) then
    if (AgePriorM(ABS(AgeDiff(A, B))+1, k) == 0.0) then
        LLg(2) = 777  
        LLg(3) = 777  ! not sibs
    endif
    if (AgePriorM(ABS(AgeDiff(A, B))+1, 3-k) == 0.0) then
        LLg(2) = 777  
    endif
    if (AgePriorM(ABS(AgeDiff(A, B))+1,6)==0.0) then   
        LLg(5)=777
        LLg(6)=777
    endif
    if (AgeDiff(A,B) <= 0) then  ! B younger than A
        LLg(1) = 777
        LLg(4) = 777
    else if (Sex(B) /= 3 .and. Sex(B)/=k) then
        LLg(1) = 777
    else if (AgePriorM(AgeDiff(A,B)+1, k+2) == 0.0 .and. &
        AgePriorM(AgeDiff(A,B)+1, 5) == 0.0) then
        LLg(4) = 777 ! not GP
    endif
endif
if (LLg(focal) == 777 .and. focal/=1) then
    LL = LLg
    return
endif

 call CalcU(A,k,B,k, LLg(7))
 
! PO?
if (LLg(1)==999) then
    call PairPO(A, B, k, LLg(1))
endif
if (focal/=1 .and. Sex(B)/=k .and. Sex(B)/=3 .and. AgeDiff(A,B)>0) then
    call PairPO(A, B, Sex(B), LLg(1))   
endif

! FS? 
if (LLg(2)==999) then
    call PairFullSib(A, B, LLg(2))  ! includes check for parent mismatch
endif

! HS?
if (LLg(3)==999) then
    if (Parent(A, 3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
        LLg(3) = 777
    else
        call PairHalfSib(A, B, k, LLg(3))  ! includes check for parent mismatch
    endif
endif

! GP?
if (LLg(4)==999) then
    call PairGP(A, B, k, focal, LLg(4))
    call PairGGP(A, B, k, LLtmpGGP)
    if (focal==3 .and. MaxLL(LLg(2:3))>LLg(4) .and. Sex(B)/=3) then  ! and AgeDiff==999
        cgp = 0
        if (Parent(A,3-k)>0)  cgp = Parent(Parent(A,3-k),Sex(B))
        if (Parent(A,3-k)<0)  cgp = GpID(Sex(B), -Parent(A,3-k), 3-k)
        if (cgp == 0) then
            call PairGP(A, B, 3-k, focal, LLCC)
            LLg(4) = MaxLL((/ LLg(4), LLCC /))    ! Note: wrong ageprior will be used.
        endif
    endif
endif

! FA/HA?
if (LLg(5)==999) then  ! FA & HA have same ageprior.
    LLtmpA = 999
    do x=1,3  ! mat, pat, FS
        call PairUA(A, B, k, x, LLtmpA(1,x))
        call PairUA(B, A, k, x, LLtmpA(2,x))
    enddo
    LLg(5) = MaxLL(LLtmpA(:,3))
    LLg(6) = MaxLL(RESHAPE(LLtmpA(:,1:2), (/2*2/) ))
endif

LLCC = 999
 call PairCC(A, B, k, LLCC) 

LL = LLg
if (AgeDiff(A,B)/=999 .and. InclAge) then
    if (LLg(1) < 0) then
        LL(1) = LLg(1) + LOG10(AgePriorM(ABS(AgeDiff(A, B))+1, 6+k))
    endif
    if (LLg(2) < 0) then
        LL(2) = LLg(2) + LOG10(AgePriorM(ABS(AgeDiff(A, B))+1, 1)) &
                + LOG10(AgePriorM(ABS(AgeDiff(A, B))+1, 2)) 
    endif
    if (LLg(3) < 0) then
        LL(3) = LLg(3) + LOG10(AgePriorM(ABS(AgeDiff(A, B))+1, k))
    endif
    do x=5,6  ! same age prior for full & half aunts/uncles
        if (LLg(x) < 0) then
            LL(x) = LLg(x) + LOG10(AgePriorM(ABS(AgeDiff(A, B))+1, 6)) 
        endif
    enddo
endif

LL(6) = MaxLL( (/LL(6), LLtmpGGP, LLCC/) )  ! use most likely 3rd degree relative

if (LLg(4) < 0 .and. AgeDiff(A,B)/=999) then
        x = 0
        if (k==1) then
            if ((AgeDiff(A,B) > 0 .and. Sex(B)==1) .or. &
            (AgeDiff(A,B) < 0 .and. Sex(A)==1)) then
                x = 3  ! mat. grandmother
            else 
                x = 5  ! mat. grandfather
            endif
        else if (k==2) then
            if ((AgeDiff(A,B) > 0 .and. Sex(B)==2) .or. &
            (AgeDiff(A,B) < 0 .and. Sex(A)==2)) then
                x = 4  ! pat. grandfather
            else 
                x = 5  ! pat. grandmother
            endif
        endif
        ALR = AgePriorM(ABS(AgeDiff(A, B))+1, x)
        if (ALR /= 0) then
        if (InclAge) then
            LL(4) = LLg(4) + LOG10(ALR)
        else
            LL(4) = LLg(4)
        endif
    else
        LL(4) = 777
    endif
endif 

LLX = 999
if (LL(2)<0 .and. (focal==2 .or. focal==3) .and. (LL(2) - MaxLL(LL)) > -thLRrel) then
    if (ALL(Parent(A,:)==0) .and. ALL(Parent(B,:)==0)) then  ! consider if HS + AU (also catches inbred A or B)
        call PairHSHA(A, B, LLX(1))
        call PairHSHA(B, A, LLX(2))
        if (LL(2) - MaxLL(LLX(1:2)) < thLRrel) then  ! conservative.
            LL(2) = 222
        endif
    endif
    if (LL(2)/=222 .and. (Parent(A,3-k)==0 .or. Parent(B,3-k)==0)) then
        call PairHSGP(A, B,k, LLX(3))
        call PairHSGP(B, A,k, LLX(4))
        if (LL(2) - MaxLL(LLX(3:4)) < thLRrel) then  
            LL(2) = 222
        endif
    endif
    if (LL(2)/=222) then ! check if inbred FS (& boost LL)
        call PairFSHA(A, B, k, LLZ(1))
        call PairFSHA(A, B, 3-k, LLZ(2))
        if (MaxLL(LLZ) > LL(2)) then
            LL(2) = MaxLL(LLZ)
        endif
    endif
    if (Parent(A,3-k) < 0 .and. Parent(B,3-k)==0) then  ! check if add instead
        call CheckAdd(B, -Parent(A,3-k), 3-k, LLC, 3)
        if ((LLC(2) - MaxLL(LLC)) < thLRrel) then  ! not FS; conservative
            LL(2) = 222
        endif
    else if (Parent(A,3-k)==0 .and. Parent(B,3-k)<0) then
        call CheckAdd(A, -Parent(B,3-k), 3-k, LLC, 3)
        if ((LLC(2) - MaxLL(LLC)) < thLRrel) then  ! not FS; conservative
            LL(2) = 222
        endif
    endif
endif
end subroutine CalcPair

! ##############################################################################################

subroutine PairPO(A, B, k, LL)
use Global
implicit none

integer, intent(IN) :: A, B, k
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp,3), PrX(3,3,3), PrPA(3), PrB(3), PrPB(3), LLtmp(3), PrPAB(3)

LL = 999
if(Parent(A,k)>0) then  ! allow dummy to be replaced (need for AddFS)
    if (Parent(A,k)==B) then
        LL = 888
    else
        LL = 777
    endif
endif

if (Sex(B)/=3 .and. Sex(B)/=k) then
    LL = 777
endif
if (LL/=999) return

PrL = 0
do l=1,nSnp
    if (Genos(l,A)==-9) cycle
    call ParProb(l, Parent(A,3-k), 3-k, A,0, PrPA)
    call ParProb(l, B, k, 0,0, PrB)
    if (Parent(A,3-k)==Parent(B,3-k)) then  ! incl. ==0
        call ParProb(l, Parent(A,3-k), 3-k, A, B, PrPAB)
        call ParProb(l, Parent(B,k), k, B, 0, PrPB)
    endif
    
    do x=1,3
        do y=1,3
            PrX(x,y,1) = OKA2P(Genos(l, A), x, y) * PrB(x) * PrPA(y)
            if (Parent(A,3-k)==0) then  ! consider close inbreeding
                PrX(x,y,2) = OKA2P(Genos(l, A), x, y) * PrB(x) * AKAP(y,x,l)
            endif
            if (Parent(A,3-k)==Parent(B,3-k)) then  ! both 0 or both /=0
                PrX(x,y,3) = OKA2P(Genos(l, A), x, y) * PrPAB(y) * &
                  SUM(AKA2P(x,y,:) * PrPB)
                if (Genos(l,B)/=-9) then
                    PrX(x,y,3) = PrX(x,y,3) * PrB(x)  ! PrB = OcA
                endif
            endif
        enddo
    enddo
    PrL(l,1) = LOG10(SUM(PrX(:,:,1)))
    if (Parent(A,3-k)==0)  PrL(l,2) = LOG10(SUM(PrX(:,:,2)))
    if (Parent(A,3-k)==Parent(B,3-k))  PrL(l,3) = LOG10(SUM(PrX(:,:,3)))
enddo

LLtmp = SUM(PrL, DIM=1)
LLtmp(1:2) = LLtmp(1:2) + Lind(B)
if (Parent(A,3-k)==0) then
    LL = MaxLL(LLtmp)
else if (Parent(A,3-k)==Parent(B,3-k)) then
    LL = LLtmp(3)
else
    LL = LLtmp(1)
endif

end subroutine PairPO

! ##############################################################################################

subroutine PairFullSib(A, B, LL)
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LL
integer :: x, y, l,k, Par(2), AncA(2,mxA), AncB(2,mxA)
double precision :: PrL(nSnp), PrXY(3,3), Px(3,2), LUX(2), LLtmp

LL = 999
Par = 0  ! joined parents of A & B
if (Parent(A,1)==Parent(B,1) .and. Parent(A,1)/=0 .and. &
    Parent(A,2)==Parent(B,2) .and. Parent(A,2)/=0) then ! already FS
    LL = 888
else 
    do k=1,2
        if (Parent(A,k) == B .or. Parent(B,k) == A) then
            LL = 777
            exit
        else if (Parent(A,k)/=Parent(B,k) .and. ((Parent(A,k)>0 .and. Parent(B,k)>0) .or. &
          (Parent(A,k)<0 .and. Parent(B,k)<0))) then
            LL = 777
            exit
        else if (Parent(A,k)/=0) then
            Par(k) = Parent(A,k)
        else
            Par(k) = Parent(B,k)
        endif        
    enddo
endif  
if (LL/=999) return

 call GetAncest(A,1,AncA)
 call GetAncest(B,1,AncB)
if (ANY(AncA == B) .or. ANY(AncB == A)) then
    LL = 777
    return
endif
   
PrL = 0 
if (Par(1) < 0 .and. Par(2)<0) then  ! call AddFS
    call CalcU(A, 1, B, 1, LUX(1))
    do k=1,2 
        if (Parent(A,k)==Par(k) .and. Parent(B,k)==0) then
            call CalcU(B, k, Par(k), k, LUX(2))
            call addFS(B, -Par(k), k, 0, k, LLtmp) 
            LL = LLtmp - LUX(2) + LUX(1)
            exit
        else if (Parent(B,k)==Par(k) .and. Parent(A,k)==0) then
            call CalcU(A, k, Par(k), k, LUX(2))
            call addFS(A, -Par(k), k, 0, k, LLtmp) 
            LL = LLtmp - LUX(2) + LUX(1)
            exit
        endif
    enddo     
else  
    do l=1, nSnp
        if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
        do k=1,2
            if (Parent(A,k)==Parent(B,k)) then  ! doesn't matter if Par(k)==0
                call ParProb(l, Par(k), k, A, B, Px(:,k))
            else if (Parent(A,k)==Par(k)) then
                call ParProb(l, Par(k), k, A, 0, Px(:,k))
            else if (Parent(B,k)==Par(k)) then
                call ParProb(l, Par(k), k, B, 0, Px(:,k))
            else
                call ParProb(l, Par(k), k, 0, 0, Px(:,k))
            endif       
        enddo 
    
        do x=1,3
            do y=1,3
                PrXY(x,y) = Px(x,1) * Px(y,2)
                if(Genos(l,A)/=-9) then
                    PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)
                endif
                if(Genos(l,B)/=-9) then
                    PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,B), x, y)
                endif
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXY))
    enddo
    LL = SUM(PrL) 
endif    

end subroutine PairFullSib

! ##############################################################################################

subroutine PairHalfSib(A, B, k, LL)
use Global
implicit none

integer, intent(IN) :: A, B, k
double precision, intent(OUT) :: LL
integer :: x,y, l, Par, Inbr, AB(2), i
double precision :: PrL(nSnp), PrX(3), PrPx(3,2), PrXY(3,3)

LL = 999
Par = 0  ! parent K
if (Parent(A,k)/=0) then
    if (Parent(A,k)/=Parent(B,k) .and. Parent(B,k)/=0) then
        LL = 777 ! mismatch
    else if (Parent(A,k)==Parent(B,k)) then
        LL = 888 ! already captured under H0
    else
        Par = Parent(A,k)
        if (Par>0) then
            if (AgeDiff(B, Par) <= 0) then  ! Par(k) younger than B
                LL = 777
            endif
        endif
    endif
else if (Parent(B,k)/=0) then
    Par = Parent(B,k)
    if (Par>0) then
        if (AgeDiff(A, Par) <= 0) then  ! Par(k) younger than A
            LL = 777
        endif
    endif
endif
if (LL/=999) return
 
AB = (/ A, B /)
Inbr = 0
if (Parent(A,3-k)==B)  Inbr = 1
if (Parent(B,3-k)==A)  Inbr = 2
PrL = 0
do l=1,nSnp
   if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    if (Par==Parent(A,k) .and. Par/=0) then
        call ParProb(l, Par, k, A, 0, PrX)
    else if (Par==Parent(B,k) .and. Par/=0) then
        call ParProb(l, Par, k, B, 0, PrX)
    else
        call ParProb(l, Par, k, 0, 0, PrX)    
    endif
    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPx(:,1))
    call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPx(:,2))
    if (Inbr==0) then
        do x=1,3
            do i=1,2
                if (Genos(l,AB(i))/=-9) then
                    PrX(x) = PrX(x) * SUM(OKA2P(Genos(l,AB(i)), x, :) * PrPx(:,i))
                endif
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrX))
    else 
        do x=1,3
            do y=1,3
                if (Genos(l,AB(Inbr))/=-9) then
                    PrXY(x,y) = PrX(x) * OKA2P(Genos(l,AB(Inbr)), x, y) * &
                      SUM(AKA2P(y, x, :) * PrPx(:,3-Inbr))
                endif
                if (Genos(l,AB(3-Inbr))/=-9) then
                    PrXY(x,y) = PrXY(x,y) * OcA(Genos(l, AB(3-Inbr)), y)
                endif
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXY))
    endif
enddo
LL = SUM(PrL)

end subroutine PairHalfSib

! ##############################################################################################

subroutine pairHSHA(A, B, LL)   ! HS via k, & parent A is HS of B via 3-k
use Global
implicit none

integer, intent(IN) :: A,B
double precision, intent(OUT) :: LL
integer :: l, x, y, z
double precision :: PrL(nSnp), PrXYZ(3,3,3)

if (ANY(Parent(A,:)/=0) .or. ANY(Parent(B,:)/=0)) then
    LL = 777
    return
endif   ! else not necessary.

PrL = 0
do l=1, nSnp
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    do x=1,3
        do y=1,3    
            do z=1,3
                PrXYZ(x,y,z) = AHWE(y,l) * AHWE(z,l) * AKAP(x, y, l)
                if (Genos(l,B)/=-9) then
                    PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,B), y, z)
                endif
                if (Genos(l,A)/=-9) then
                    PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,A), x, z)
                endif
            enddo
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine pairHSHA

! ##############################################################################################

subroutine pairFSHA(A, B, k, LL)   ! inbred FS: parent k offspring of parent 3-k
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp), PrXY(3,3), PrY(3)

if (Parent(A,k)/=0 .or. Parent(B,k)/=0) then
    LL = 444   ! TODO (prob. not necessary)
    return
endif  

PrL = 0
do l=1, nSnp
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    call ParProb(l, Parent(A,3-k), 3-k, -1,0, PrY)  ! TODO: remainder of sibship
    do x=1,3
        do y=1,3    
            PrXY(x,y) = PrY(y) * AKAP(x, y, l)
            if (Genos(l,B)/=-9) then
                PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,B), x, y)
            endif
            if (Genos(l,A)/=-9) then
                PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)
            endif
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine pairFSHA

! ##############################################################################################

subroutine pairHSGP(A, B,k, LL)   ! HS via k, B is GP of A via 3-k
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y, z
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrPB(3)

if (Parent(A,3-k)/=0) then
    LL = 444
    return
endif 

PrL = 0
do l=1, nSnp
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPB)
    do x=1,3  ! parent 3-k of A, offspring of B
        do y=1,3  ! shared parent k 
            do z=1,3  ! B
                PrXYZ(x,y,z) = AKAP(x, z, l) * AHWE(y,l) * SUM(AKA2P(z,y,:) * PrPB)
                if (Genos(l,A)/=-9) then
                    PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,A), x, y)
                endif
                if (Genos(l,B)/=-9) then
                    PrXYZ(x,y,z) = PrXYZ(x,y,z) * OcA(Genos(l,B), z)
                endif
            enddo
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine pairHSGP

! ##############################################################################################

subroutine PairGP(A, B, k, focal, LL)  ! calculates LL that B is maternal(k=1) or paternal(k=2) grandparent of A
use Global
implicit none

integer, intent(IN) :: A,B,K, focal
double precision, intent(OUT) :: LL
integer :: l, x, y, curGP(2), m, z, AncB(2, mxA), v
double precision :: PrL(nSnp,3), PrPA(3,2), PrG(3), LLtmp(2), PrXZ(3,3,3,3), PrB(3), &
  PrGx(3), PrPB(3), PrV(3)

LL = 999
 curGP = 0  
if (Parent(A,k)>0) then
    curGP = Parent(Parent(A,k),:) ! current grandparents of A (shortcut)
else if (Parent(A,k)<0) then
    curGP = GpID(:, -Parent(A,k), k)
endif

if (AgeDiff(A, B) <= 0) then  ! B younger than A
    LL = 777
else if (Sex(B)/=3) then
    if (curGP(Sex(B))/=0) then
        if (curGP(Sex(B))/=B) then
            LL = 777  ! conflict
        endif
    endif
    m = Sex(B)
else 
    if (curGP(1)/=0 .and. curGP(2)/=0) then
        do y=1,2
            if (curGP(y) == B) then
                LL = 888
                exit
            else
                LL = 777
            endif
        enddo
    endif
    if (curGP(1)==0) then
        m = 1  ! doesn't really matter.
    else
        m = 2 
    endif
endif

LLtmp = 999
if (Parent(A,k)>0 .and. LL==999) then
    if (AgeDiff(Parent(A,k), B) <= 0) then  ! B younger than Parent(A,k)
        LL = 777 
    else
        if (Sex(B)/=3) then
            call PairPO(Parent(A,k), B, Sex(B), LLtmp(1))
        else
            call PairPO(Parent(A,k), B, 1, LLtmp(1))
        endif
        if (LLtmp(1) > 0) then    ! impossible
            LL = 777
        else 
            call CalcU(Parent(A,k), k, B,k, LLtmp(2))
            if (LLtmp(1) - LLtmp(2) < thLRrel) then
                LL = 777
            endif
        endif
    endif
endif
if (LL/=999) return

 call GetAncest(B, k, AncB)
if (ANY(AncB == A)) then
    LL = 777
endif
if (LL/=999) return

PrL = 0
do l=1,nSnp
    if (Genos(l,A)==-9) cycle
    call ParProb(l, curGP(3-m), 3-m, 0,0, PrG)
    call ParProb(l, Parent(A,k), k, A, -4, PrPA(:, k))  ! excl. both GPs & A
    if (Parent(A,3-k)==Parent(B,3-k)) then
        call ParProb(l, Parent(A,3-k), 3-k, A, B, PrPA(:, 3-k))
    else
        call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA(:, 3-k))
    endif
    call ParProb(l, B, 0, 0, 0, PrB)
    call ParProb(l, Parent(B,k), k, B, 0, PrPB)
    if (Parent(A,3-k) < 0) then
        call ParProb(l, GpID(3-m, -Parent(A,3-k),3-k), 3-m, 0, 0, PrGx)   
    endif

    PrXZ = 1
    do x=1,3  ! PA(k)
        if (Parent(A,k)/=0) then
            PrXZ(x,:,:,:) = PrPA(x, k)
        endif      
        do y=1,3  ! PA(3-k)
            do z=1,3  !  PrG(3-m)
                PrXZ(x,y,z,:) = PrXZ(x,y,z,:) * OKA2P(Genos(l,A), x, y) * PrG(z)
                if (.not. (Parent(B,3-k)==Parent(A,3-k) .and. Parent(A,3-k)/=0)) then
                    PrXZ(x,y,z,1) = PrXZ(x,y,z,1) * PrPA(y,3-k) * SUM(AKA2P(x, :, z) * PrB)  !non-inbred
                endif
                if (Parent(A,3-k)==0) then
                    PrXZ(x,y,z,2) =  PrXZ(x,y,z,2) * SUM(AKA2P(x, :, z) * AKAP(y,:,l) * PrB)  !inbreeding loop
                endif
                if (Parent(B,3-k)==0 .or. Parent(B,3-k)==Parent(A,3-k)) then
                    do v=1,3
                        PrV(v) = AKA2P(x, v, z) * SUM(AKA2P(v,y,:) * PrPB) * PrPA(y,3-k)
                        if (Genos(l,B)/=-9) then
                            PrV(v) = PrV(v) * OcA(Genos(l,B), v)
                        endif
                    enddo
                    PrXZ(x,y,z,3) =  PrXZ(x,y,z,3) * SUM(PrV)
                endif
            enddo
        enddo
    enddo
    do z=1,3  ! inbred/non-inbred
        PrL(l,z) = LOG10(SUM(PrXZ(:,:,:,z)))
    enddo
enddo

if (ANY(SUM(PrL, DIM=1) < 0)) then
    if (focal==4) then
        LL = MaxLL(SUM(PrL, DIM=1)) + Lind(B)
    else if (Parent(B,3-k)==Parent(A,3-k)) then
        LL = MaxLL((/ SUM(PrL(:,1)), SUM(PrL(:,3)) /)) + Lind(B) 
    else if (SUM(PrL(:,1)) < 0) then
        LL = SUM(PrL(:,1)) + Lind(B)
    else
        LL = 777
    endif
else
    LL = 777
endif
if (LL < 0 .and. Parent(A,k)>0) then
    LL = LL - Lind(Parent(A,k))
endif
end subroutine PairGP

! ##############################################################################################

subroutine PairGGP(A, B, k, LL)   
! calculates LL that B is maternal(k=1) or paternal(k=2) great-grandparent of A
! todo?: merge with GGP for clusters
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y,z,w, AncB(2, mxA), m, MaybeF
double precision :: PrL(nSnp,2), PrXY(3,3), PrXZ(3,3,3,3, 2), LLtmp(3), PrPA(3,2), &
  PrB(3), PrG(3), PrPAX(3)  

LL = 999
if (AgeDiff(A, B) <= 0) then  ! B younger than A
    LL = 777
else if (B==Parent(A,k)) then
    LL = 777
else
    call GetAncest(B, k, AncB)
    if (ANY(AncB == A)) then
        LL = 777
    else if (Parent(A,3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
        LL = 444   ! not impossible, just unlikely. TODO: implement below.
    endif
endif
if (LL==777) return

if (Parent(A,k)>0) then 
    if (ANY(Parent(Parent(A,k), :)/=0)) then
        LL = 444    ! should be picked up elsewere
    else
        call PairGP(Parent(A,k), B, k, 4, LLtmp(1))
        if (LLtmp(1) > 0) then    
            LL = LLtmp(1)
        endif
    endif
else if (Parent(A,k)<0) then
    if (ANY(GpID(:,-Parent(A,k),k)/=0)) LL = 444
endif
if (LL/=999) return

MaybeF = 0
if (Parent(A,3-k) < 0 .and. Parent(A,k)==0) then
    if (ANY(GpID(:, -Parent(A,3-k),3-k)==0)) then
        MaybeF = 1  ! maybe inbreeding loop
    endif
endif

PrL = 0    
do l=1,nSnp
    PrXY = 0
    PrXZ = 0
    if (Genos(l,A)==-9) cycle
    call ParProb(l, B, 0, 0, 0, PrB)
    if (Parent(A,k)==0) then
        PrPA(:,k) = 1
    else
        call ParProb(l, Parent(A,k), k, A, 0, PrPA(:,k))
    endif
    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA(:,3-k))
    if (MaybeF==1) then  ! at least 1 GP ==0
        call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX)  ! exclude GP contribution
        if (ALL(GpID(:, -Parent(A,3-k),3-k)==0)) then
            PrG = AHWE(:,l)
        else
            do m=1,2
                if (GpID(m, -Parent(A,3-k),3-k)/=0) then   ! either one ==0 
                    call ParProb(l, GpID(m, -Parent(A,3-k),3-k), m, 0, 0, PrG)    
                endif
            enddo
        endif
    endif
    
    do x=1,3  
        do y=1,3 
            PrXY(x,y) = SUM(OKA2P(Genos(l,A),x,:) * PrPA(:,3-k)) * PrPA(x,k) * &
              AKAP(x,y,l) * SUM(AKAP(y, :, l) * PrB)
            do z=1,3  !consider double GGP (2x k, or k & 3-k)
                do w=1,3
!                    PrXZ(x,y,z,w,1) = OKA2P(Genos(l,A), x,z) * PrPA(x,k) * PrPA(z,3-k) * &
 !                     AKA2P(x,y,w) * SUM(AKAP(y,:,l) * AKAP(w,:,l) * PrB)  ! Parent(A,k) inbred - indistinguishable from GP
                    if (Parent(A,3-k)==0) then  
                        PrXZ(x,y,z,w,2) = OKA2P(Genos(l,A), x,z) * PrPA(x,k) * &
                          AKAP(x,y,l) * AKAP(z,w,l) * SUM(AKAP(y,:,l) * AKAP(w,:,l) * PrB)
                    else if (MaybeF==1) then
                        PrXZ(x,y,z,w,2) = OKA2P(Genos(l,A), x,z) * PrPAX(z) * SUM(AKA2P(z,w,:) * PrG) * &
                          AKAP(x,y,l) * SUM(AKAP(y,:,l) * AKAP(w,:,l) * PrB)
                    endif
                enddo
            enddo
        enddo
    enddo
    PrL(l,1) = LOG10(SUM(PrXY))
    if (Parent(A,3-k)==0 .or. MaybeF==1) then
        PrL(l,2) = LOG10(SUM(PrXZ(:,:,:,:,2)))
    endif
enddo

if (Parent(A,3-k)==0 .or. MaybeF==1) then
    LL = MaxLL(SUM(PrL, DIM=1)) + Lind(B)
else
    LL = SUM(PrL(:,1)) + Lind(B)
endif

end subroutine PairGGP

! ##############################################################################################

subroutine PairUA(A, B, kA, kB, LL)
! B half sib or full sib (kB=3) of parent kA of A?
use Global
implicit none

integer, intent(IN) :: A,B,kA, kB  ! kB=3 : full sibs
double precision, intent(OUT) :: LL
integer :: l, x, g, y, z, GG(2), GA(2), PB(2), PA, i, nA, r, u,j,e,Ei,m, &
  AncA(2,mxA), AncG(2, 2,mxA), AA(maxSibSize), cat(maxSibSize), &
  doneB(maxSibSize), BB(maxSibSize), nB, catG, GGP, catB(maxSibSize), &
  nBx(2), BBx(maxSibSize, 2), Bj
double precision :: PrL(nSnp), PrG(3,2), PrXYZ(3,3,3), PrPA(3), PrA(3), &
  PrPB(3), PrGA(3), PrAB(3,3,3,2), PrE(3), PrH(3), PrGG(3)

if (A>0) then  
    nA = 1
    AA(1) = A
    PA = Parent(A,kA)
    if (PA<0) then
        LL = 444  ! possible but not implemented; TODO
        return
    else if (PA>0) then
        GA = Parent(PA, :)
    else
        GA = 0
    endif
else
    nA = nS(-A, kA)
    AA(1:nA) = SibID(1:nA, -A, kA)
    PA = A
    GA = GpID(:, -A, kA)
endif

if (B > 0) then
    nB = 1
    BB(1) = B
    PB = Parent(B,:)
    do m=1,2
        if (kB/=3 .and. m/=kB)  cycle
        if (Parent(B, m) >=0) then
            nBx(m) = 1
            BBx(1, m) = B
        else 
            nBx(m) = nS(-Parent(B, m), m)
            BBx(1:nBx(m), m) = SibID(1:nBx(m), -Parent(B, m), m)
        endif
    enddo
else if (B < 0) then
    nB = nS(-B, kB)
    BB(1:nB) = SibID(1:nB, -B, kB)
    PB(kB) = B
    PB(3-kB) = 0
endif

! prep
LL = 999
if (kB < 3) then
    if (B < 0 .and. GA(kB) == B) then
        LL = 888
    else if (B > 0 .and. GA(kB) == PB(kB) .and. PB(kB)/=0) then
        LL = 888
    else if (GA(3-kB)==PB(3-kB) .and. GA(3-kB)/=0) then  ! B>0; would become FA not HA
        LL = 777
    endif
else if (GA(1)==PB(1) .and. GA(2)==PB(2) .and. GA(1)/=0 .and. GA(2)/=0) then  ! kB==3
    LL = 888
endif
if (LL /= 999) return

GG = 0  ! parent of B, GP of A
AncG = 0
 call GetAncest(A, kA, AncA)
do x=1,2
    if (x/=kB .and. kB/=3) cycle
    if (GA(x)==0) then
        GG(x) = PB(x)
    else if (GA(x)/=PB(x) .and. PB(x)/=0) then
        LL = 777
    else
        GG(x) = GA(x)
    endif
    if (ANY(AA(1:nA)==GG(x))) then
        LL = 777
    endif
    call GetAncest(GG(x), x, AncG(x, :, :))
enddo

if (A > 0) then
    if (ANY(AncG == A)) then
        LL = 777
    endif
else if (A < 0) then
    if (ANY(AncG(:, kA, 2:mxA) == A)) then
        LL = 777
    endif
endif
if (B > 0) then
    if (ANY(AncG == B)) then
        LL = 777
    endif
else if (B < 0) then
    if (ANY(AncG(:, kB, 3:mxA) == B)) then
        LL = 777
    endif
endif
if (kB<3) then
    if (B<0) then
        if (ANY(AncA(kB,3:mxA)==B))  LL = 444  ! B is GGP; possible but not (yet) implemented
    else if (B>0) then
        if (ANY(AncA(:,3:mxA)==B))  LL = 444
    endif
endif
if (LL /= 999) return

do x=2,mxA
    do y=1,2
        do g=1,2
            if (AncG(g,y,x) > 0) then
                if (A > 0) then
                    if (AgeDiff(A, AncG(g,y,x)) < 0) then
                        LL = 777  ! A older than putative ancestor
                    endif 
                else if (A<0) then
                    if (ANY(AgeDiff(SibID(1:nS(-A,kA),-A,kA), AncG(g,y,x)) < 0)) then
                        LL = 777  ! A older than putative ancestor
                    endif
                endif
                if (x==2) cycle  ! parent of A doesn't need to be older than B
                if (B > 0) then
                    if (AgeDiff(B, AncG(g,y,x)) < 0) then
                        LL = 777  
                    endif 
                else if (B<0) then
                    if (ANY(AgeDiff(SibID(1:nS(-B,kB),-B,kB), AncG(g,y,x)) < 0)) then
                        LL = 777 
                    endif
                endif
            endif
        enddo
    enddo
enddo
if (LL /= 999) return

!==============================================

if (A>0 .and.  B>0) then  ! quicker.
    if (ALL(Parent(A,:)>=0) .and. ALL(Parent(B,:)>=0)) then
        PrL = 0
        do l=1, nSnp
            if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
            call ParProb(l, Parent(A,3-kA), 3-kA, A, 0, PrPA)
            if (kB == 3) then
                do g=1,2
                    call ParProb(l, GG(g), g, 0, 0, PrG(:,g))  ! >=0
                enddo        
            else
                call ParProb(l, GG(kB), kB, 0, 0, PrG(:,kB))  ! >=0
                call ParProb(l, GA(3-kB), 3-kB, PA, 0, PrGA)
                call ParProb(l, PB(3-kB), 3-kB, B, 0, PrPB)
            endif
            
            PrXYZ = 0
            do z=1,3
                do y=1,3
                    do x=1,3
                        if (kB == 3) then
                            PrXYZ(x,y,z) = AKA2P(x,y,z) * PrG(y, 1) * PrG(z, 2)
                        else
                            PrXYZ(x,y,z) = AKA2P(x,y,z) * PrG(y, kB) * PrGA(z)   
                        endif
                        if (Genos(l,A)/=-9) then
                            PrXYZ(x,y,z) = PrXYZ(x,y,z) * SUM(OKA2P(Genos(l,A), x, :) * PrPA)  
                        endif
                        if (Genos(l,B)/=-9) then
                            if (kB==3) then
                                PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,B), y, z)   
                            else
                                PrXYZ(x,y,z) = PrXYZ(x,y,z) * SUM(OKA2P(Genos(l,B), y, :) * PrPB)
                            endif
                        endif
                        if (Parent(A,kA)>0) then
                            if (Genos(l,Parent(A,kA))/=-9) then
                                PrXYZ(x,y,z) = PrXYZ(x,y,z) * OcA(Genos(l,Parent(A,kA)), x)
                            endif
                        endif
                    enddo
                enddo
            enddo
            PrL(l) = LOG10(SUM(PrXYZ))
        enddo
        LL = SUM(PrL)

        return
    endif
endif

!==============================================

cat=0
catG = 0
catB = 0
GGP = 0
do i = 1, nA
    if (kA/=kB .and. GG(3-kA)<0 .and. Parent(AA(i), 3-kA) == GG(3-kA)) then  !incl. kB=3
        cat(i) = 1  
    else if (kA==kB .and. Parent(AA(i), 3-kA)<0) then
        if (Parent(AA(i), 3-kA)==GA(3-kA)) then
            cat(i) = 2
        else
            do j=1, nB
                if (AA(i) == BB(j)) cycle
                if (Parent(AA(i), 3-kA) == Parent(BB(j), 3-kA)) then
                    cat(i) = 3
                endif
            enddo
        endif
    endif
    if (Parent(AA(i),3-kA)/=0) then
        do g=1,2
            if (kB/=g .and. kB/=3) cycle
            if (GG(g) > 0) then
                if (Parent(AA(i),3-kA) == Parent(GG(g), 3-kA)) then
                    cat(i) = 5  ! TODO? 4+g when kB==3
                    catG = 2
                endif
            else if (GG(g) < 0) then
                if (Parent(AA(i),3-kA) == GpID(3-kA, -GG(g),g)) then
                    cat(i) = 5  ! TODO? 4+g when kB==3
                    catG = 2
                endif
            endif
        enddo
    endif
enddo
if (kB/=3) then
    do j=1, nB
        if (Parent(BB(j),3-kB)==0) cycle
        do g=1,2
            if (GG(g) > 0) then
                if (Parent(BB(j),3-kB) == Parent(GG(g), 3-kB)) then
                    catB(j) = 5  
                    catG = 3
                endif
            else if (GG(g) < 0) then
                if (Parent(BB(j),3-kB) == GpID(3-kB, -GG(g),g)) then
                    catB(j) = 5  
                    catG = 3
                endif
            endif
        enddo
    enddo
endif

if (kA/=kB .and. kB/=3) then
    do j=1,nB
        if (Parent(BB(j), 3-kB) == PA .and. A>0 .and. PA<0) then
            cat(nA+1) = 4  ! only possible if A>0 (else cat(i)=1)
        endif
    enddo
endif
if (kB/=3) then
    if (GG(kB) > 0) then
        if (Parent(GG(kB), 3-kB) == GA(3-kB) .and. GA(3-kB)/=0) then
            catG = 1
            GGP = Parent(GG(kB), kB)
        endif
    else if (GG(kB) < 0) then
        if (GpID(3-kB, -GG(kB), kB) == GA(3-kB) .and. GA(3-kB)/=0) then
            catG = 1
            GGP = GpID(kB, -GG(kB), kB)
        endif
    endif
endif

PrL = 0
DoneB = 0
do l=1,nSnp
    if (A>0 .and. B>0) then
        if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    endif
    do g=1,2
        if (g/=kB .and. kB/=3) cycle
        if (catG==0) then
            if (ALL(cat==0) .and. ALL(GG >=0)) then
                call ParProb(l, GG(g), g, B, 0, PrG(:,g))  ! B>0
            else
                call ParProb(l, GG(g), g, -1, 0, PrG(:,g))
            endif
        else
            if (GG(g) > 0) then  ! TODO: catG/=0, kB==3
                call ParProb(l, GG(g), g, 0, 0, PrG(:,g)) 
                if (catG==1) then
                    call ParProb(l, GGP, kB, GG(g), 0, PrGG)
                else if (catG==2) then
                    call ParProb(l, GGP, 3-kA, GG(g), 0, PrGG)
                else if (catG==3) then
                    call ParProb(l, GGP, 3-kB, GG(g), 0, PrGG)
                endif
            else if (GG(g) < 0) then
                PrG(:,g) = 1
                if (catG==1) then
                    call ParProb(l, GGP, kB, 0, 0, PrGG)
                else if (catG==2) then
                    call ParProb(l, GGP, 3-kA, 0, 0, PrGG)
                else if (catG==3) then
                    call ParProb(l, GGP, 3-kB, 0, 0, PrGG)
                endif
            endif
        endif
    enddo
        if (kB/=3) then  ! TODO: if(ANY(cat==2))
        if (ANY(cat==2)) then
            call ParProb(l, GA(3-kB), 3-kB, -1, 0, PrGA)
        else if (PA>0) then
            call ParProb(l, GA(3-kB), 3-kB, PA, 0, PrGA)
        else if (catG==1 .and. GG(kB)>0) then
            call ParProb(l, GA(3-kB), 3-kB, GG(kB), 0, PrGA)
        else
            call ParProb(l, GA(3-kB), 3-kB, 0, 0, PrGA)
        endif
        if (B>0) then
            call ParProb(l, PB(3-kB), 3-kB, B, 0, PrPB)  
        endif
    endif

    ! === 
    
    if (ALL(cat==0) .and. ALL(GG >=0) .and. catG==0 .and. ALL(catB==0)) then
        if (A < 0) then
            PrA = Xpr(1,:,l, -A,kA)
        else if (A>0) then
            call ParProb(l, Parent(A,3-kA), 3-kA, A, 0, PrPA)
            if (Genos(l,A)/=-9) then
                do x=1,3 
                    PrA(x) = SUM(OKA2P(Genos(l,A), x, :) * PrPA)
                enddo
            else
                PrA = 1 ! ? AHWE(:,l)  
            endif
            if (Parent(A,kA)>0) then
                if (Genos(l,Parent(A,kA))/=-9) then
                    PrA = PrA * OcA(Genos(l,Parent(A,kA)), :)  ! prob. x given obs.
                endif
            else if (Parent(A,kA)<0) then
                do x=1,3
                    PrH(x) = Xpr(1,x,l, -Parent(A,kA),kA) / SUM(OKA2P(Genos(l,A), x, :) * PrPA)
                enddo
                PrH = PrH / SUM(PrH)  ! total contribution from other sibs summing to 1
                do x=1,3 
                    PrA(x) = SUM(OKA2P(Genos(l,A),x,:) * PrPA * PrH(x))
                enddo
            endif
        endif
    
        do x=1,3  ! PA, kA
            do y=1,3  ! PrG, kB
                do z=1,3  ! PrGA, 3-kB / PrG, 3-kB
                    if (kB==3 .and. B>0) then  ! SA/PA FS of B; call AddFS for SB/=0
                        PrXYZ(x,y,z) = PrA(x) * AKA2P(x,y,z) * PrG(y,3-kA) * PrG(z,kA)    
                        if (Genos(l,B)/=-9) then
                            PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,B), y, z)
                        endif
                    else 
                        PrXYZ(x,y,z) = PrA(x) * AKA2P(x,y,z) * PrGA(z)
                        if (B>0) then                 
                            if (Genos(l,B)/=-9) then
                                PrXYZ(x,y,z) = PrXYZ(x,y,z) * PrG(y, kB) * &
                                  SUM(OKA2P(Genos(l,B),y,:) * PrPB)
                            else
                                PrXYZ(x,y,z) = PrXYZ(x,y,z) * PrG(y, kB)
                            endif
                        else if (B<0) then
                            PrXYZ(x,y,z) = PrXYZ(x,y,z) * XPr(3,y,l, -B, kB)
                        endif
                    endif     
                enddo
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXYZ))   
    else 
        PrAB = 0
        if (PA>0) then
            call ParProb(l, PA, kA, 0, 0, PrPA)
        else
            PrPA = 1
        endif
        
        do y=1,3  ! PrG, kB
            do x=1,3  ! PA, kA
                do z=1,3
                    if (kB==3) then
                        PrAB(x,y,z,:) = PrPA(x) * AKA2P(x,y,z) * PrG(z, kA) * PrG(y, 3-kA)
                    else if (catG==1) then
                        PrAB(x,y,z,:) = PrPA(x) * AKA2P(x,y,z) * PrGA(z) * &
                          SUM(AKA2P(y,z,:) * PrGG) * PrG(y, kB) 
                    else
                        PrAB(x,y,z,:) = PrPA(x) * AKA2P(x,y,z) * PrGA(z) * PrG(y, kB) 
                    endif
                enddo
            enddo 
                
            do x=1,3
                doneB = 0
                do r=1, nA
                    if (A<0 .and. NFS(AA(r))==0) cycle 
                    if (cat(r)==0 .or. cat(r)>2) then
                        call ParProb(l, Parent(AA(r), 3-kA), 3-kA, -1, 0, PrE)
                    else
                        PrE = 1
                    endif
                    if (Parent(AA(r), 3-kA) < 0 .and. Parent(AA(r), 3-kA) /= GG(3-kA)) then   ! not: Parent(AA(r), 3-kA)== PB(kB)
                        do e=1,3
                            do g=1, nS(-Parent(AA(r), 3-kA), 3-kA)
                                Ei = SibID(g, -Parent(AA(r), 3-kA), 3-kA)
                                if (nFS(Ei) == 0) cycle 
                                if (Parent(Ei, kA) == PA .and. PA/=0) cycle
                                if (kB<3) then
                                    if (Parent(Ei, kB)== PB(kB) .and. PB(kB)/=0) cycle
                                endif
!                                if (ANY(GG == Ei)) cycle
                                call ParProb(l, Parent(Ei, kA), kA, Ei, -1, PrH)  ! Assume is not one of GPs
                                do i=1, nFS(Ei)
                                    if (A>0 .and. FSID(i, Ei)==A) cycle
                                    if (B>0 .and. FSID(i, Ei)==B) cycle
                                    if (Genos(l,FSID(i, Ei))==-9) cycle
                                    PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                                enddo
                                PrE(e) = PrE(e) * SUM(PrH)
                            enddo
                        enddo
                        if (A>0) then
                            do i=1, MAX(nFS(AA(r)),1)
                                if (FSID(i, AA(r))==A) cycle
                                if (ANY(GG == FSID(i, AA(r)))) cycle
                                if (Genos(l,FSID(i, AA(r)))==-9) cycle
                                PrE = PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)
                            enddo
                        endif
                        if (cat(r)==3 .and. B>0) then  ! kA==kB, Parent(AA(r), 3-kA)/=0   .or. cat(nA+1)==4
                            do j=1,nB
                                if (Parent(BB(j), 3-kB) /= Parent(AA(r), 3-kB)) cycle
                                do i=1, MAX(nFS(BB(j)),1)
                                    if (B<0 .and. nFS(BB(j))==0) cycle
                                    if (B>0 .and. FSID(1, BB(j))==B) cycle
                                    if (Genos(l,FSID(i, BB(j)))==-9) cycle
                                    PrE = PrE * OKA2P(Genos(l,FSID(i,BB(j))), y, :)
                                enddo
                            enddo
                        endif 
                    endif 
                    if (cat(r)==5) then
                        do e=1,3
                            PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)
                        enddo
                    endif
                    if (cat(r)==0 .or. cat(r)==3 .or. cat(r)==5) then  ! incl. cat(nA+1)==4
                        if (SUM(PrE)<3)  PrAB(x,y,:,1) = PrAB(x,y,:,1) * SUM(PrE)
                    else if (cat(r)==1) then  ! Parent(AA(i), 3-kA) == GG(3-kA), kA/=kB
                        PrAB(x,y,:,1) = PrAB(x,y,:,1) * PrE(y)
                    else if (cat(r)==2) then  ! Parent(AA(i), 3-kA)==GA(3-kA), kA==kB
                        PrAB(x,y,:,1) = PrAB(x,y,:,1) * PrE
                    endif          
                    do i=1, MAX(nFS(AA(r)),1)
                        if (A>0 .and. FSID(i, AA(r))/=A) cycle
                        if (Genos(l,FSID(i, AA(r)))==-9) cycle
                        PrE = PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :) 
                    enddo
                    if (cat(r)==3) then  ! kA==kB, Parent(AA(r), 3-kA)/=0   .or. cat(nA+1)==4
                        do j=1,nB
                            if (Parent(BB(j), 3-kB) /= Parent(AA(r), 3-kB)) cycle
                            if (B<0 .and. nFS(BB(j))==0) cycle
                            do i=1, MAX(nFS(BB(j)),1)
                                if (B>0 .and. FSID(i, BB(j))/=B) cycle
                                if (Genos(l,FSID(i, BB(j)))==-9) cycle
                                PrE = PrE * OKA2P(Genos(l,FSID(i,BB(j))), y, :)
                            enddo
                            DoneB(j) = 1
                        enddo
                    endif
                    if (cat(r)==0 .or. cat(r)==3 .or. cat(r)==5) then  ! incl. cat(nA+1)==4
                        if (SUM(PrE)<3)  PrAB(x,y,:,2) = PrAB(x,y,:,2) * SUM(PrE)
                    else if (cat(r)==1) then
                        PrAB(x,y,:,2) = PrAB(x,y,:,2) * PrE(y)
                    else if (cat(r)==2) then
                        PrAB(x,y,:,2) = PrAB(x,y,:,2) * PrE
                    endif
                enddo  ! r
            enddo  ! x
            
            if (B<0) then            
                do j=1,nB
                    if (nFS(BB(j))==0) cycle
                    if (DoneB(j)==1) cycle
                    if (kA/=kB .and. A<0 .and. Parent(BB(j), 3-kB)==PA) cycle  ! TODO double check A<0
                    call ParProb(l, Parent(BB(j), 3-kB), 3-kB, -1, 0, PrE)
                    if (Parent(BB(j), 3-kB) < 0) then  
                        do e=1,3
                            do g=1, nS(-Parent(BB(j), 3-kB), 3-kB)
                                Ei = SibID(g, -Parent(BB(j), 3-kB), 3-kB)
                                if (nFS(Ei) == 0) cycle
                                if (Parent(Ei, kB) == PB(kB) .and. PB(kB)/=0) cycle  
                                if (Parent(Ei, kA)== PA .and. PA/=0) cycle
                                call ParProb(l, Parent(Ei, kB), kB, Ei, -1, PrH) 
                                do i=1, nFS(Ei)
                                    if (A>0 .and. FSID(i, Ei)==A) cycle
                                    if (Genos(l,FSID(i, Ei))/=-9) then
                                        PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                                    endif
                                enddo
                                PrE(e) = PrE(e) * SUM(PrH) 
                            enddo
                        enddo                   
                    endif
                    if (catB(j)==5) then
                        do e=1,3
                            PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)
                        enddo
                    endif
                    PrAB(:,y,:,1) = PrAB(:,y,:,1) * SUM(PrE)
                    do u=1, nFS(BB(j))
                        if (Genos(l,FSID(u, BB(j)))==-9) cycle
                        PrE = PrE * OKA2P(Genos(l,FSID(u,BB(j))), y, :)
                    enddo
                    PrAB(:,y,:,2) = PrAB(:,y,:,2) * SUM(PrE) 
                enddo   ! j
            else if (B>0) then 
                do z=1,3  
                    do m=1,2
                        if (kB/=3 .and. m/=kB)  cycle
!                        if (Parent(B, m) < 0) then
                        do j=1, nBx(m)
                            Bj = BBx(j, m)
                            if (nFS(Bj)==0 .and. Bj/=B) cycle
                            if (kB==3 .and. Parent(Bj, 3-m)==GG(3-m) .and. GG(3-m)/=0) then  ! FS of B
                                if (.not. (Parent(B,1)<0 .and. Parent(B,2)<0 .and. m==2)) then
                                    do u=1, nFS(Bj)
                                        if (FSID(u,Bj)==B) cycle
                                        if (ANY(AA(1:nA)==FSID(u,Bj))) cycle
                                        if (Genos(l,FSID(u, Bj))==-9) cycle
                                        PrAB(:,y,z,:) = PrAB(:,y,z,:) * OKA2P(Genos(l,FSID(u,Bj)), y, z) 
                                    enddo
                                endif
                            else 
                                call ParProb(l, Parent(Bj, 3-m), 3-m, -1, 0, PrE)
                                if (Parent(Bj, 3-m) < 0) then
                                    do e=1,3
                                        do g=1, nS(-Parent(Bj, 3-m), 3-m)
                                            Ei = SibID(g, -Parent(Bj, 3-m), 3-m)
                                            if (nFS(Ei) == 0) cycle
                                            if (kB<3 .and. Parent(Ei,m)==GG(m) .and. GG(m)/=0) cycle 
                                            if (ANY(AA(1:nA)==Ei)) cycle
                                            call ParProb(l, Parent(Ei, m), m, Ei, -1, PrH)  
                                            do i=1, nFS(Ei)
                                                if (ANY(AA(1:nA)==FSID(i, Ei))) cycle
                                                if (FSID(i, Ei)==B) cycle
                                                if (FSID(i, Ei)==PA) cycle
                                                if (Genos(l,FSID(i, Ei))/=-9) then
                                                    PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                                                endif
                                            enddo
                                            PrE(e) = PrE(e) * SUM(PrH)  
                                        enddo
                                    enddo 
                                endif
                                do u=1, MAX(nFS(Bj),1)
                                    if (FSID(u,Bj)==B) cycle
                                    if (FSID(u, Bj)==PA) cycle
                                    if (ANY(AA(1:nA)==FSID(u,Bj))) cycle
                                    if (Genos(l,FSID(u, Bj))==-9) cycle
                                    if (kB==3 .and. m==kA) then
                                        PrE = PrE * OKA2P(Genos(l,FSID(u,Bj)), z, :) 
                                    else
                                        PrE = PrE * OKA2P(Genos(l,FSID(u,Bj)), y, :) 
                                    endif
                                enddo
                                if (catB(j)==5) then
                                    do e=1,3
                                        PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)
                                    enddo
                                endif
                                if (kB<3) then
                                    PrAB(:,y,z,1) = PrAB(:,y,z,1) * SUM(PrE)
                                    do u=1, MAX(nFS(Bj),1)
!                                    if (A==-22 .and. kA==2 .and. B==4429 .and. kB==2)  cycle !! TEMP !!
                                        if (FSID(u,Bj)==B .and. Genos(l,B)/=-9) then
                                            PrE = PrE * OKA2P(Genos(l,B), y, :) 
                                        endif
                                    enddo
                                    PrAB(:,y,z,2) = PrAB(:,y,z,2) * SUM(PrE) 
                                else
                                    PrAB(:,y,z,:) = PrAB(:,y,z,:) * SUM(PrE)
                                endif
                            endif
                        enddo  ! j_m
                    enddo  ! m
                    if (kB==3 .and. Genos(l,B)/=-9) then
                        PrAB(:,y,z,2) = PrAB(:,y,z,2) * OKA2P(Genos(l,B), y, z)
                    endif
                enddo  ! z
            endif  ! kB==3
        enddo  ! y
        PrL(l) = LOG10(SUM(PrAB(:,:,:,2))) - LOG10(SUM(PrAB(:,:,:,1)))
    endif
enddo

LL = SUM(PrL)

end subroutine PairUA 

! ##############################################################################################

subroutine pairCC(A,B,k, LL)  ! 1st cousins
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y, u, v
double precision :: PrL(nSnp), PrXY(3,3), PrUV, PrPA(3), PrPB(3), PrC(3,3)

LL = 999
if (Parent(A,k)/=0 .and. Parent(B,k)/=0) then
    if (Parent(A,k)==Parent(B,k)) then
        LL = 777
    endif
endif

if (Parent(A,k)/=0 .or. Parent(B,k)/=0) then  ! TODO: this is temp.
    LL = 444
endif
if (Parent(A,3-k) == Parent(B,3-k) .and. Parent(A, 3-k)/=0) then
    LL = 444
endif
if (LL/=999) return

PrL = 0
do l=1, nSnp
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    call ParProb(l, Parent(A,3-k), 3-k, A,0, PrPA)
    call ParProb(l, Parent(B,3-k), 3-k, B,0, PrPB)
    
    do u=1,3  ! GG1
        do v=1,3  ! GG2
            PrUV = AHWE(u,l) * AHWE(v,l)
            do x=1,3  !PA
                do y=1,3    !PB
                    PrXY(x,y) = AKA2P(x,u,v) * AKA2P(y,u,v) * PrUV
                    if (Genos(l,A)/=-9) then
                        PrXY(x,y) = PrXY(x,y) * SUM(OKA2P(Genos(l,A), x, :) * PrPA)
                    endif
                    if (Genos(l,B)/=-9) then
                        PrXY(x,y) = PrXY(x,y) * SUM(OKA2P(Genos(l,B), y, :) * PrPB)
                    endif
                enddo
            enddo
            PrC(u,v) = SUM(PrXY)
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrC))
enddo

LL = SUM(PrL)

end subroutine pairCC

! ##############################################################################################

subroutine Clustering(LR)
use Global
implicit none

integer, intent(IN) :: LR  ! -1: first round, +1: last round
integer :: k, x, n, m, ij(2), sx(2), topX, PK, u
double precision :: LLtmp(7), dLL, LLx(7, 2)

do x=1, nPairs
    LLtmp = 999
    ij = PairID(x,:)
    PK = PairType(x)
    do k=1,2
        sx(1) = -Parent(ij(1),k)  ! interested in sibships, which have neg. numbers
        sx(2) = -Parent(ij(2),k)
        if (sx(1)==0 .and. sx(2)==0) then
            if (AgeDiff(ij(1),ij(2))==999) then
                call CalcPair(ij(1), ij(2), k, LR>0, LLx(:,1), 3)  ! HS/FS symmetrical, rest not.
                call CalcPair(ij(2), ij(1), k, LR>0, LLx(:,2), 3)
                do u=1,7
                    LLtmp(u) = MaxLL(LLx(u,:))
                enddo
            else if(AgeDiff(ij(1),ij(2))>=0) then
                call CalcPair(ij(1), ij(2), k, LR>0, LLtmp, 3)  ! only rely on ageprior in last round
            else
                call CalcPair(ij(2), ij(1), k, LR>0, LLtmp, 3)
            endif
            call BestRel(LLtmp, 2, topX, dLL)  
            if (PK==3 .and. LR<=0 .and. (topX==3 .or. (topX==2 .and. dLL<2*thLRrel))) cycle  ! unclear if MHS or PHS
            call BestRel(LLtmp, 3, topX, dLL)
            if (topX==2 .or. (topX==3 .and. (dLL > 2*thLRrel .or. LR>=0))) then
                nC(k) = nC(k)+1  ! new sibship (pair)
                if (Parent(ij(1), 3-k)/=0 .and. Parent(ij(1), 3-k)==Parent(ij(2), 3-k)) then
                    call MakeFS(ij(1), ij(2))
                endif                
                nS(nC(k), k) = 2
                SibID(1:2, nC(k), k) = ij
                Parent(ij(1), k) = -nC(k)
                Parent(ij(2), k) = -nC(k)
                if (Parent(ij(1), 3-k) /=0 .and. Parent(ij(1), 3-k)==Parent(ij(2), 3-k)) then
                    call MakeFS(ij(1), ij(2))
                endif
                call CalcCLL(nC(k), k)
                do n=1, 2
                    u = SibID(n, nC(k), k)
                    if (Parent(u,3-k) < 0) then
                        call CalcCLL(-Parent(u,3-k), 3-k)
                    endif                
                    call CalcLind(SibID(n, nC(k), k))
                enddo
                call CalcCLL(nC(k), k)
            endif
        else if (sx(1)>0 .and. sx(2)>0 .and. sx(1) /= sx(2)) then   ! else do nothing
            call CheckMerge(sx(1), sx(2), k,k, LLtmp, 1)
            call BestRel(LLtmp, 3,topX, dLL)              
            if (topX==1 .and. dLL > thLRrel * MIN(nS(sx(1),k), nS(sx(2),k))) then
                call DoMerge(sx(1), sx(2), k)
            endif
        else
            do m=1,2
                if (sx(m)>0 .and. sx(3-m)==0) then
                    call CheckAdd(ij(3-m), sx(m), k, LLtmp, 3)
                    call BestRel(LLtmp, 3,topX, dLL)
                    if (topX==2 .or. topX==3) then
                       call DoAdd(ij(3-m), sx(m), k)
                    endif
                endif
            enddo
        endif
    enddo
enddo

end subroutine Clustering

! ##############################################################################################

subroutine Merging  ! check if any of the existing clusters can be merged
use Global
implicit none

integer :: k, s, r, n, topX, xr
double precision :: LLtmp(7,2), LLm(7), dLL

do k=1,2
    do s=1,nC(k)-1
        r = s
        do xr=s+1, nC(k)
            r = r + 1
            if (r > nC(k)) exit   ! possible due to merged sibships
            topX = 0
            LLtmp = 999
            call CheckMerge(s, r, k, k, LLtmp(:,1), 1)           
            do n = 1,7
                LLm(n) = MaxLL(LLtmp(n,:))
            enddo
            call BestRel(LLm, 1,topX, dLL)
            if (topX==1 .and. dLL > thLRrel * MIN(nS(s,k), nS(r,k))) then
                call DoMerge(s, r, k)
                r = r-1  ! otherwise a cluster is skipped
            endif
        enddo
    enddo
enddo

end subroutine Merging

! ##############################################################################################

subroutine GrowClusters
! for each individual, check if they can be added to any sibship clusters
use Global
implicit none

integer :: k, s, i, nMaybe, ClM(100), topX
double precision :: LLtmp(7), dLL, dLLM(100)
    
do i=1, nInd
    do k=1,2
        if (Parent(i,k)/=0) cycle  
        nMaybe = 0
        Clm = 0
        dLLM = 999
        do s=1,nC(k)  ! sibships
            if (GpID(1,s,k)==i .or. GpID(2,s,k)==i) cycle          
            LLtmp = 999
            call CheckAdd(i, s, k, LLtmp, 3) 
            call BestRel(LLtmp, 3,topX, dLL)
            if (topX==2 .or. topX==3) then 
                nMaybe = nMaybe+1
                Clm(nMaybe) = s 
                dLLM(nMaybe) = dLL
            endif
        enddo
        if (nMaybe==1) then
            call DoAdd(i, Clm(1), k)
        endif
    enddo
enddo
 
end subroutine GrowClusters

! ##############################################################################################

subroutine SibParent  ! for each sibship, check if a parent can be assigned
use Global
implicit none

integer :: k, s, xs, i, n,maybe, topX, CurNumC(2), &
  OH, Par, MaybeOpp, CurPar(2)
double precision :: LLtmp(7), dLL

CurNumC = nC
maxOppHom = MaxMismatch - FLOOR(-nSNP * Er)   ! round up to nearest integer  
do k=1,2
    s = 0
    do xs=1, CurNumC(k)
        s = s+1
        if (s > nC(k)) exit   ! possible due to dissolved sibships after parent assigned

        do i=1,nInd
            if (Sex(i)/=k .and. Sex(i)/=3) cycle
            if (Parent(i,k)==-s) cycle
            if (GpID(1,s,k)==i .or. GpID(2,s,k)==i) cycle
            maybe=1
            CurPar = 0
            do n=1,nS(s,k)
                if (AgeDiff(i,SibID(n,s,k)) > 0 .and. AgeDiff(i,SibID(n,s,k))/= 999) then  ! candidate i younger than sib j 
                    maybe = 0
                    exit
                endif
                call CalcOH(i, SibID(n,s,k), OH)
                if (OH > maxOppHom) then
                    maybe = 0
                    exit
                endif   
            enddo
            if (maybe==0) cycle
            do n=1,2
                if (GpID(n,s,k) <= 0) cycle
                if (AgeDiff(i,GpID(n,s,k)) < 0 .and. AgeDiff(i,GpID(n,s,k))/= 999) then  ! candidate i younger than sib j 
                    maybe = 0
                    exit
                endif
                call CalcOH(i, GpID(n,s,k), OH)
                if (OH > maxOppHom) then
                    maybe = 0
                    exit
                endif   
            enddo
            if (maybe==0) cycle
            if (BY(i)==-999) then  ! check if any in S accidentally assigned as parent of i
                do n=1,2
                    if (ANY(SibID(:,s,k)==Parent(i,n))) then
                        CurPar = Parent(i,:)
                        Parent(i,n) = 0
                    endif
                enddo
            endif
            LLtmp = 999
            call CheckAdd(i, s, k, LLtmp, 1)
            call BestRel(LLtmp, 1, topX, dLL)
            if (topX/=1) then
                maybe = 0
                if (ANY(CurPar/=0))  Parent(i,:) = CurPar
                cycle
            endif
            if (Sex(i)==3) then  ! check if parent of opposite sex instead
                MaybeOpp = 1
                call getFSpar(s, k, Par)
                if (Par > 0)  MaybeOpp = 0
                if (Par==0 .and. ANY(Parent(SibID(1:nS(s,k), s,k),3-k)>0))  MaybeOpp = 0
                if (Par/=0 .and. Parent(i, 3-k) == Par)  MaybeOpp = 0  ! already are HS via opp parent   
                if (MaybeOpp == 1) then
                    if (Par < 0) then  ! may have more/fewer sibs
                        call CheckAdd(i, -Par, 3-k, LLtmp, 1)
                        call BestRel(LLtmp, 1, topX, dLL)
                        if (topX==1)  maybe = 0
                    else if (Par == 0) then
                        call PairPO(SibID(1, s, k), i, 3-k, LLtmp(1))   ! prob. not necessary todo with all sibs
                        call CalcU(SibID(1, s, k), k, i, 3-k, LLtmp(2))
                        if (LLtmp(1)<0 .and. LLtmp(1) - LLtmp(2) > thLRrel)  maybe = 0
                    endif
                endif
            endif
            if (maybe==0) cycle
            do n=1,nS(s,k)  ! assign parent i to all sibs in S  (note: no change in FS groups)
                Parent(SibID(n, s, k), k) = i
            enddo   
            if (Sex(i)==3)  Sex(i) = k
            call DoMerge(0, s, k)  !removes cluster s, shifts all subsequent ones, fixes GPs            
            s = s-1  ! otherwise a cluster is skipped
            exit
        enddo  ! i
    enddo ! s
enddo ! k
    
end subroutine SibParent

! ##############################################################################################

subroutine MoreParent  ! for each individual, check if a parent can be assigned now.
use Global
implicit none

integer :: i, j, l, OH, Lboth, AncJ(2,mxA), k, u
double precision :: LRP

do i=1, nInd
    if (ALL(Parent(i,:)/=0)) cycle 
    do j=1, nInd  ! candidate parent.
        if (i==j) cycle
        if (ANY(Parent(j,:)==i) .or. ANY(Parent(i,:)==j)) cycle
        if (AgeDiff(i,j) <= 0)  cycle  ! note: unknown = 999 > 0
        if (AgeDiff(i,j) == 999) then
            if (ALL(Parent(i,:)==0) .and. ALL(Parent(j,:)==0)) cycle
            if (BY(i)/=-999 .and. BY(j)==-999) then
                call GetAncest(j,1,AncJ)
                if (ANY(AncJ == i)) cycle
                do k=1,2
                    do u=2, mxA
                        if (AncJ(k,u)>0) then
                            if (AgeDiff(i, AncJ(k,u)) <= 0)  cycle
                        endif
                    enddo
                enddo
            endif
        endif
        if (Sex(j) == 3) then
            if (Parent(i,1)==0 .and. Parent(i,2)==0) cycle 
        else
            if (Parent(i, Sex(j)) /= 0) cycle 
        endif
        OH = 0
        do l=1,nSnp
            if ((Genos(l,i)==1).and.(Genos(l,j)==3)) then
                OH = OH+1
                if (OH > maxOppHom) exit
            endif                       
            if ((Genos(l,i)==3).and.(Genos(l,j)==1)) then
                OH = OH+1
                if (OH > maxOppHom) exit
            endif                       
        enddo
        if (OH > maxOppHom) cycle  
        if (OH <= maxOppHom) then
            Lboth = COUNT(Genos(:,i)/=-9 .and. Genos(:,j)/=-9)    
            if (Lboth < nSnp/4.0)  cycle   ! more than 3/4th of markers missing
        endif               
        call CalcPO(i, j, LRP)
        if (LRP < -thLR) cycle
        call CalcPOZ(i,j, .TRUE.)  ! assigns parent as side effect; use age prior
    enddo
enddo

end subroutine MoreParent

! ##############################################################################################

subroutine SibGrandparents 
! for each sibship, find the most likely male and/or female grandparent
use Global
implicit none

integer :: k, s, i,j, r, m, par, xs, candGP(20, 2), nCG(2), curGP(2), u, v, AncR(2,mxA), curNG
double precision :: LRG, curCL, ALRtmp, LL(7)

do k=1,2
    s = 0
    do xs=1, nC(k) 
        s = s+1
        if (s > nC(k)) exit
        nCG = 0
        CandGP = 0
        curNG = 0
        CurGP = GpID(:,s,k)
        curCL = CLL(s,k)
        do m=1,2
            if (GpID(m,s,k)/=0) then
                nCG(m) = 1
                CandGP(1,m) = GpID(m,s,k)
                curNG = curNG + 1
            endif
        enddo
        call CalcCLL(s,k)
        if (ALL(Parent(SibID(1:nS(s,k), s, k), 3-k) < 0)) then
            call getFSpar(s, k, par)
            if (par < 0) then
                if (nS(-par, 3-k) == nS(s,k))  cycle    ! all sibs are FS, cannot tell if mat or pat GP
            endif          
        endif
        do i=1,nInd
            if (Parent(i,k)==-s) cycle
            if (GpID(1,s,k)==i .or. GpID(2,s,k)==i) cycle
            if (Sex(i)==3 .and. ALL(GpID(:,s,k)==0)) cycle ! cannot tell if GM or GF
            call CalcAgeLR(-s,k, i,0, 0,1, ALRtmp)
            if (ALRtmp == 777) cycle
            if(GPID(1,s,k)/=0 .and. GpID(1,s,k)==Parent(i,1) .and. &  ! FS of dummy
                GPID(2,s,k)/=0 .and. GpID(2,s,k)==Parent(i,2)) cycle
            if (nS(s,k)>1) then
                call QGP(i, 0, s, k, LRG)  ! quick check
                if (LRG < -thLR)  cycle    !  *nS(s,k)
            else if (ns(s,k)==1) then
                call PairQHS(i, SibID(1,s,k), LRG)
                if (LRG < -thLR)  cycle
            endif
            call GetAncest(i,1,AncR)
            if (ANY(AncR(k,:) == -s)) cycle
            call CheckAdd(i, s, k, LL, 4)
            if (LL(4)>222 .or. LL(4) - LL(7) < thLRrel) cycle  !NOT bestrel: keep when LL(4)==LL(5)
            if (Sex(i)/=3) then   ! 222 = rely on agePrior only
                if (ncG(Sex(i))<20) then  ! arbitrary threshold to limit comp. time
                    nCG(Sex(i)) = nCG(Sex(i)) + 1
                    CandGP(nCG(Sex(i)), Sex(i)) = i    
                endif
            else
                do m=1,2
                    if (ncG(m)<20) then
                        ncG(m) = nCG(m) + 1
                        CandGP(nCG(m), m) = i
                    endif
                enddo
            endif
        enddo
        
        do m=1,2
            do r=1, nC(m) 
                if (m==k .and. s==r) cycle
                if (GPID(m,s,k) == -r) cycle   ! already assigned.
                if (GPID(k,r,m) == -s) cycle 
                if (ANY(SibID(1:nS(s,k),s,k) == GpID(1,r,m)) .or. &
                  ANY(SibID(1:nS(s,k),s,k) == GpID(2,r,m))) cycle
                if (nS(r,m)==1 .and. ANY(SibID(1:nS(s,k),s,k) == SibID(1,r,m)))  cycle
                if (ALL(GpID(:,s,k)==0) .and. COUNT(nFS(SibID(1:nS(r,m),r,m))>0)==1) cycle  ! FS family, can't tell if F or M is GP
                call QGP(-s, k, r, m,  LRG)  ! TODO: age checks? (CalcAgeLRCAU) 
                if (LRG < -thLR) cycle  ! conservative.  *nS(s,k)
                call GetAncest(-r,m,AncR)
                if (ANY(AncR(k,:)==-s)) cycle
                do j=1, nS(s,k)
                    if (ANY(AncR == SibID(j,s,k))) cycle ! covers 1 extra generation
                enddo  
                call CalcAgeLR(-s,k, -r,m, 0,1, ALRtmp)
                if (ALRtmp == 777) cycle
                call CheckMerge(s, r, k, m, LL, 4)
                if (LL(4)>222 .or. LL(4) - LL(7) < thLRrel) cycle  ! 222 = rely on agePrior only
                if (ncG(m) < 20) then
                    nCG(m) = nCG(m) + 1
                    CandGP(nCG(m), m) = -r
                endif                
            enddo
        enddo
        if (SUM(nCG) == curNG) cycle  ! no new candidate GPs
        
        if (ALL(nCG>0) .and. ANY(nCG>1)) then   ! test all combo's
            do u = 1, nCG(1)
                do v = 1, nCG(2)
                    if (CandGP(u, 1)>0 .and. CandGP(u, 1)==CandGP(v, 2))  cycle  ! when unknown sex.
                    call CalcGPz(s, k, CandGP(u, 1), 1) ! GpID(1,s,k) = CandGP(u, 1)
                    call CalcGPz(s, k, CandGP(v, 2), 2)   ! removes unlikely GPs
                enddo
            enddo
        else if (ANY(nCG>1)) then  ! & nCG(3-m)=0
            do m=1,2
                if (nCG(m)>1) then
                    do v = 1, nCG(m)
                        call CalcGPz(s, k, CandGP(v, m), m)
                    enddo
                endif
            enddo
        else if (CandGP(1, 1)>0 .and. CandGP(1, 1)==CandGP(1, 2)) then
            cycle
        else if (ALL(nCG==1)) then
            do m=1,2
                call CalcGPz(s, k, CandGP(1, m), m)
            enddo
        else
            do m=1,2
                if (nCG(m)>0) then
                    call CalcGPz(s, k, CandGP(1, m), m)
                endif
            enddo
        endif
!        if (CLL(s,k) - curCL < -thLRrel) then  ! made worse somehow.
 !           GpID(:,s,k) = CurGP
 !           call CalcCLL(s,k)
 !       endif
        if (nS(s,k)==1 .and. ALL(GpID(:, s,k)==0)) then  ! remove cluster s
            v = SibID(1,s,k)
            Parent(v,k) = 0
            call DoMerge(0, s, k) 
            call CalcLind(v)
        endif 
    enddo  ! s
enddo  ! k
                        
end subroutine SibGrandparents

! ##############################################################################################

subroutine CalcGPz(SA, kA, B, kB) ! B or SB grandparent of sibship SA?
use Global
implicit none

integer, intent(IN) :: SA, kA, B, kB
integer :: m,n, CurGP(2), topX, x, mid(5), y, CY(4), kY(4), v
double precision :: LLA(2,7,7), dLL, TopLL, LLcp(3,2), LLINOUT(2), LLU(3), LLtmp(3), ALR(3)
logical :: ConPar(4,4)

curGP = GpID(:,SA,kA)
dLL = 999
TopLL = 999
LLA = 999
call CalcCLL(SA,kA)
if (B < 0) then
    n = kB
else if (B>0) then
    if (Sex(B)/=3) then
        n = Sex(B)
    else if (GpID(1,SA,kA)==0) then
        n = 1
    else !if (GpID(2,SA,kA)==0) then
        n = 2  ! TODO: check in both configs (as in POZ)
    endif
endif

if (curGP(1)==0 .and. curGP(2)==0) then
    if (B > 0) then
        call checkAdd(B, SA, kA, LLA(1,:,7), 4)
    else if (B < 0) then
        call checkMerge(SA, -B, kA, kB, LLA(1,:,7), 4) 
    endif
    call BestRel(LLA(1,:,7), 4, topX, dLL)
    if (topX==4) then   !   .and. dLL>thLRrel*nS(SA,kA)
        GpID(n,SA,kA) = B
    endif       
else
    ConPar = .FALSE.
    GpID(:,SA,kA) = 0
    call CalcCLL(SA,kA)
    if (ANY(curGP<0) .or. B<0) then
        do m=1,2
            if (curGP(m)==0) cycle
            call Connected(CurGP(m), m, -SA, kA, ConPar(4,m))
            call Connected(CurGP(m), m, B, n, ConPar(3,m)) 
            call Connected(curGP(3-m),3-m,CurGP(m), m, ConPar(2,1))
        enddo
        call Connected(B,n, -SA, kA, ConPar(4,3))
    endif
    
    LLU = 0 
    do m=1,2
        if (curGP(m)==0) cycle
        call CalcU(curGP(m),m, -SA,kA, LLtmp(1))
        LLU(m) = LLtmp(1) - CLL(SA,kA)
    enddo  
    call CalcU(B,n, -SA,kA, LLtmp(2)) 
    LLU(3) = LLtmp(2) - CLL(SA,kA)
    do m=1,2
        LLcp(m,m) = LLU(3-m) + LLU(3) 
    enddo
    LLcp(3,:) = LLU(1) + LLU(2)  ! B  
    
    CY = (/ CurGP(1), CurGP(2), B, -SA /)
    kY = (/ 1, 2, n, kA /)
    if (ANY(ConPar)) then 
        do m=1,2        
            do y=1,3  ! focal: CurGP, B 
                if (y/=m .and. y/=3) cycle           
                if (y==1) then
                    call CalcU(CY(2),kY(2), CY(3),kY(3), LLCP(1,m))
                else if (y==2) then
                    call CalcU(CY(1),kY(1), CY(3),kY(3), LLCP(2,m))    
                else if (y==3) then
                    call CalcU(CY(1),kY(1), CY(2),kY(2), LLCP(3,m))
                endif
                do x=1,3
                    if (ConPar(4,x) .and. x/=y) then ! focal relationship: SA, CY(y)
                        call CalcU(CY(x),kY(x), CY(4),kY(4), LLtmp(1))
                        call CalcU(CY(x),kY(x), 0,0, LLtmp(2))
                        call CalcU(CY(4),kY(4), 0,0, LLtmp(3))
                        LLCP(y,m) = LLCP(y,m) + (LLtmp(1)-LLtmp(2)-LLtmp(3))
                    endif
                    do v=1,2
                        if (ConPar(x,v) .and. (x==y .or. v==y)) then
                            call CalcU(CY(x),kY(x), curGP(v),v, LLtmp(1))
                            call CalcU(CY(x),kY(x), 0,0, LLtmp(2))
                            call CalcU(curGP(v),v, 0,0, LLtmp(3))
                            LLCP(y,m) = LLCP(y,m) + (LLtmp(1)-LLtmp(2)-LLtmp(3))
                        endif
                    enddo
                enddo
            enddo
        enddo
    endif
    ! age based prior prob.
    ALR = 0
    do y=1,3
        call CalcAgeLR(-SA,kA, CY(y),kY(y), 0,1, ALR(y))
    enddo 
    
    mid = (/1,2,3,5,6/)
    GpID(:,SA,kA) = CurGP 
    call CalcCLL(SA,kA)
    dLL = 0
    LLtmp = 999
    do m=1,2 ! sex currently assigned GP
        if (CurGP(m)==0) cycle    
        
        if (B > 0) then
            call checkAdd(B, SA, kA, LLA(m,:,4), 4)   ! CurGP(m)=GP + A_7 
            if (Sex(B)==3)  dLL = LLA(m,4,4)
        else if (B < 0) then
            call checkMerge(SA, -B, kA, kB, LLA(m,:,4), 4) 
        endif 
        GpID(m,SA,kA) = 0
        call CalcCLL(SA,kA)
        if (B > 0) then
            call checkAdd(B, SA, kA, LLA(m,:,7), 4)   ! A_7
            if (curGP(m)/=0 .and. Parent(B,m)==curGP(m)) then ! HS w. sib = curGP being GP
                LLA(m, 6, 7) = 888
            endif
            if (curGP(m)>0) then
                if (Parent(CurGP(m),n)==B)  LLA(m, 6, 7) = 888   
            else if (curGP(m) < 0) then
                if (GpID(n,-CurGP(m),m)==B)  LLA(m, 6, 7) = 888
            endif
        else if (B < 0) then
            call checkMerge(SA, -B, kA, kB, LLA(m,:,7), 4) 
            if (curGP(m)/=0 .and. (GpID(m,-B,kB)==curGP(m) .or. &  ! parentHFS w. sib = curGP being GP
              ANY(SibID(1:nS(-B,kB),-B,kB)==curGP(m)))) then  ! SB being GGP = curGP being GP
                LLA(m, 6, 7) = 888   ! TODO: consider other 3rd degree rel only  / but LL ought to be identical to LLU
            endif
        endif
        if (curGP(m) > 0) then
            call checkAdd(CurGP(m), SA,kA, LLA(m,7,:), 4)  ! CurGP(m)_7
        else if (curGP(m) < 0) then
            call checkMerge(SA, -CurGP(m), kA, m, LLA(m,7,:), 4)
            if (m/=n .and. ANY(Parent(SibID(1:nS(-curGP(m),m),-curGP(m),m),3-m)==B)) then
                call PairUA(-SA, CurGP(m), kA, m, LLA(m,7,4))   ! not counting SA FS with a sib
            endif
        endif
        
        if (LLA(m,4,4)<0 .or. LLA(m,4,7)<0) then  ! Else not possible
            GpID(n,SA,kA) = B
            call CalcCLL(SA,kA)
            if (curGP(m) > 0) then
                call checkAdd(CurGP(m),SA,kA, LLA(m,4,:), 4)  ! B=GP + CurGP(m)_7
            else if (curGP(m) < 0) then
                call checkMerge(SA, -CurGP(m), kA, m, LLA(m,4,:), 4)
            endif
        endif
        WHERE (LLA(m,mid,4)<0) LLA(m,mid,4) = LLA(m,mid,4) + LLcp(3,m) + ALR(m)
        WHERE (LLA(m,mid,7)<0) LLA(m,mid,7) = LLA(m,mid,7) + LLcp(3,m)
        WHERE (LLA(m,4,:)<0) LLA(m,4,:) = LLA(m,4,:) + LLcp(m,m) + ALR(3)
        WHERE (LLA(m,7,:)<0) LLA(m,7,:) = LLA(m,7,:) + LLcp(m,m)
        
        if (n==m) then
            if (B>0) then
                if (Sex(B) == 3) then
                    LLA(m,4,4) = dLL  ! from checkAdd(B, SA,..)
                else
                    LLA(m,4,4) = 777
                endif
            else   
                LLA(m,4,4) = 777 ! cannot have 2 same-sex GPs
            endif
        endif
        
        GpID(:,SA,kA) = CurGP  ! restore
        call CalcCLL(SA,kA)  
    enddo    
    TopLL = MaxLL(RESHAPE(LLA(:,:,:), (/2*7*7/)))
    
    if (LLA(3-n,4,4)==TopLL .or. ANY(LLA(n,4,:)==TopLL)) then
        GpID(n,SA,kA) = B
        GPID(3-n, SA, kA) = CurGP(3-n)
    else if ((TopLL - MaxLL(RESHAPE(LLA(:,:,4), (/2*7/)))) < thLRrel) then ! keep CurGP (both)
        GpID(:,SA,kA) = curGP
    else if (ANY(LLA(:,4,:) == TopLL)) then  ! only B
        GpID(n,   SA, kA) = B
        GpID(3-n, SA, kA) = 0
    else if (ANY(LLA(n, :, :) == TopLL)) then  ! keep CurGP(3-n)
        GpID(n,   SA, kA) = 0
        GpID(3-n, SA, kA) = curGP(3-n)
    else if (ANY(LLA(3-n, :, :) == TopLL)) then  ! keep CurGP(n)
        GpID(n,   SA, kA) = curGP(n)
        GpID(3-n, SA, kA) = 0
    endif
endif

call CalcCLL(SA, kA)
do n=1, nS(SA,kA)
    if (Parent(SibID(n,sA,kA),3-kA) < 0) then
        call CalcCLL(-Parent(SibID(n,sA,kA),3-kA), 3-kA)
    endif         
    call calcLind(SibID(n, SA, kA))
enddo
call CalcCLL(SA, kA)
 
end subroutine CalcGPz

! ##############################################################################################

subroutine GGpairs  ! find & assign grandparents of singletons
use Global
implicit none

integer :: i, j, k, x, TopX
double precision :: LLtmp(7), dLL, LRS

do i=1, nInd
    if (Parent(i,1)==0 .and. Parent(i,2)==0) cycle
    if (Parent(i,1)/=0 .and. Parent(i,2)/=0) cycle
    if (BY(i) < 0) cycle
    do k=1,2
        if (Parent(i,k)/=0) cycle
        do j=1, nInd
            if (AgeDiff(i,j) <= 0)  cycle
            if (AgeDiff(i,j) == 999) cycle
            if (Sex(j)==3) cycle
            x = 0
            if (k==1) then
                if (Sex(j)==1) then
                    x = 3 ! mat. grandmother
                else 
                    x = 5  ! mat. grandfather
                endif
            else if (k==2) then
                if (Sex(j)==2) then
                    x = 4  ! pat. grandfather
                else 
                    x = 5  ! pat. grandmother
                endif
            endif
            if (LOG10(AgePriorM(AgeDiff(i, j)+1, x)) < 0.1) cycle  ! unlikely
            call PairQHS(i, j, LRS)  ! quick check
            if (LRS < -thLR) cycle
            
            call CalcPair(i, j, k, .FALSE., LLtmp, 4)
            call BestRel(LLtmp, 4, topX, dLL)
            if (topX==4 .and. dLL > 2*thLRrel .and. (LLtmp(4)-LLtmp(7) > 2*thLR)) then  
                nC(k) = nC(k) + 1
                nS(nC(k),k) = 1
                SibID(1, nC(k), k) = i
                Parent(i, k) = -nC(k)
                GpID(Sex(j), nc(k), k) = j
                call CalcCLL(nC(k), k)
                call CalcLind(i)
                if (Parent(i, 3-k) < 0) then
                    call CalcCLL(-Parent(i, 3-k), 3-k)
                endif
                exit  ! consider alternative GP's with GPZ
            endif 
        enddo
    enddo
enddo

end subroutine GGpairs

! ##############################################################################################

subroutine Qadd(A, SB, kB, LR)
use Global
implicit none

integer, intent(IN) :: A, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x
double precision :: PrL(nSnp), PrX(3)

PrL = 0
do l=1,nSnp
    if (Genos(l,A)==-9) cycle
    do x=1,3
        PrX(x) = OKAP(Genos(l,A), x, l) * DumP(x,l,SB,kB) / AHWE(x,l)
    enddo   ! simple LL identical for HS and GP
    PrL(l) = LOG10(SUM(PrX))
enddo
LR = SUM(PrL)

end subroutine Qadd

! ##############################################################################################

subroutine QGP(A, kA, SB, kB, LR)  ! A indiv or dummy
use Global
implicit none

integer, intent(IN) :: A, kA, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x, y
double precision :: PrL(nSnp), PrX(3), PrXY(3,3)

PrL = 0
do l=1,nSnp
    do x=1,3
        PrX(x) = XPr(1,x,l,SB,kB) * AHWE(x,l)
        do y=1,3
            if (A>0) then
                PrXY(x,y) = XPr(1,x,l,SB,kB) * AKAP(x,y,l) * LindG(y,l,A)
            else if (A<0) then
                PrXY(x,y) = XPr(1,x,l,SB,kB) * AKAP(x,y,l) * DumP(y,l,-A,kA)
            endif
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY)) - LOG10(SUM(PrX))
enddo
LR = SUM(PrL)

end subroutine QGP

! ##############################################################################################

subroutine CheckAdd(A, SB, k, LL, focal)
use Global
implicit none

integer, intent(IN) :: A, SB, k, focal
double precision, intent(OUT) :: LL(7)
double precision :: LLg(7),  LLtmp(2,2), ALR(7), LLz(5), LRHS, ALRx(2,2), LLM(3)
integer :: x, y, Par, MaybeOpp, i, ParTmp(2), npt
 
LL = 999
LLg = 999
LLz = 999
LRHS = 999
ALR = 999

! ensure LL up to date
 call CalcCLL(SB, k)
 call CalcLind(A)
 if (Parent(A,3-k)<0) then
    call CalcCLL(-Parent(A,3-k), 3-k)
endif 
! quick check
call Qadd(A, SB, k, LRHS)  ! 2nd degree relatives vs unrelated
 if (LRHS < -thLR .and. (focal/=4 .and. focal/=7)) return   ! need all LL's for GPZ
if (focal==1) then
    if (Sex(A)/=3 .and. Sex(A)/=k) then
        LLg(1) = 777
        return
    endif
endif

call CalcU(A,k, -SB, k, LLg(7))   ! unrelated
LL(7) = LLg(7)

 !=======
if (LLg(1)/=777) then
    call AddParent(A, SB, k, LLg(1))  ! A parent of SB
endif
LL(1) = LLg(1)
if (focal==1 .and. (LL(1) > 0 .or. LL(1) - LL(7) < thLRrel)) return
 call AddFS(A, SB, k,0,k, LLg(2))
 call AddSib(A, SB, k, LLg(3))
call CalcAgeLR(A,0, -SB,k, 0,1, ALR(3))
ALR(2) = ALR(3)
do x=2,3
    if (LLg(x) < 0 .and. ALR(x) /= 777) then
        LL(x) = LLg(x) + ALR(x)
    else
        LL(x) = 777
    endif
enddo
if (focal==3 .and. ((LL(2) > 0 .and. LL(3)>0) .or. &
  ((LL(2) - LL(7) < thLRrel) .and. (LL(3) - LL(7) < thLRrel)))) return
call AddGP(A, SB, k, LLg(4))
call CalcAgeLR(-SB,k, A,0, 0,1, ALR(4))
if (LLg(4) < 0) then 
    if (ALR(4) /= 777) then
        LL(4) = LLg(4) + ALR(4)
    else
        LL(4) = 777
    endif
else
    LL(4) = LLg(4)
endif

! FAU
 call pairUA(-SB, A, k, 3, LLg(5))    ! A FS with SB?
if (Parent(A, 3-k)==0) then  
    call getFSpar(SB, k, Par)
    LLtmp = 999
    if (Par /= 0) then
        call pairUA(A, -SB, k, k, LLtmp(1,1))
        LLg(5) = MaxLL((/ LLg(5), LLtmp(1,1) /))
    endif
endif

! HAU
LLtmp = 999
do x=1,2
    call pairUA(A, -SB, x, k, LLtmp(1,x))  
    if (Parent(A,x)>0) then
        call CalcAgeLR(Parent(A,x),0, -SB,k, 0,1, ALRx(1,x))
    else
        call CalcAgeLR(A,0, -SB,k, x,4, ALRx(1,x))
    endif
    call pairUA(-SB, A, k, x, LLtmp(2,x))
    call CalcAgeLR(-SB,k, A,0, x,3, ALRx(2,x))
    do y=1,2
        if (LLtmp(y,x) < 0 .and. ALRx(y,x) /= 777) then
            LLtmp(y,x) = LLtmp(y,x) + ALRx(y,x)
        else
            LLtmp(y,x) = 777
        endif
    enddo
enddo

if (LLg(5)<0 .and. ALRx(2,1)/=777 .and. ALRx(2,2)/=777) then
    LL(5) = LLg(5) + ALRx(2,1) + ALRx(2,2)  ! A FS of SB
else
    LL(5) = 777
endif
  LL(6) = MaxLL(RESHAPE(LLtmp, (/2*2/) ))

if ((LL(focal)<0 .and. LL(focal)>LL(7)) .or. focal==4) then
    call AddGGP(A, SB, k, LLz(1))
    if (LLz(1) < 0 .and. ALR(4)/=777 .and. LLg(4)<0) then 
        call CalcAgeLR(-SB,k, A,0, 1,4, ALR(5))
        call CalcAgeLR(-SB,k, A,0, 2,4, ALR(6))
        if (ALR(5)/=777 .or. ALR(6)/=777) then
            LLz(1) = LLz(1) + MAXVAL(ALR(5:6), MASK=ALR(5:6)/=777)
        else
            LLz(1) = 777
        endif
    endif
    call ParentHFS(A, 0,1, SB, k,3, LLz(2))
    call ParentHFS(A, 0,2, SB, k,3, LLz(3))
    do x=1,2   ! as checkmerge: full great-uncle  (2x 1/4)
        if (GpID(x,SB,k) <0) then  ! else cond. indep.
            call PairUA(GpID(x,SB,k), A, x, 3, LLz(3+x))
            if (LLz(3+x) < 0) then
                LLz(3+x) = LLz(3+x) - CLL(-GpID(x,SB,k), x) + CLL(SB,k)  ! approx.
            endif
        endif
    enddo
    LL(6) = MaxLL((/LL(6), LLz/))
endif
            
LLM = 999    ! TODO: this piece slows it down quite a bit; do more checks
if ((MaxLL(LL)==LL(3) .or. MaxLL(LL)==LL(2)) .and. (focal==2 .or. focal==3) .and. &
  Parent(A,3-k)==0) then  ! check if not HS via opp parent only
    MaybeOpp = 1
    call getFSpar(SB, k, Par)
    if (Par > 0)  MaybeOpp = 0
    if (Par/=0 .and. Parent(A, 3-k) == Par)  MaybeOpp = 0  ! already are HS via opp parent  
    if (Par==0) then
        if (ANY(Parent(SibID(1:nS(SB,k), SB,k),3-k)>0)) then
            MaybeOpp = 0
        else   ! check if opp. parent possibly to be merged
            npt = 0
            ParTmp = 0
            do i=1, nS(SB,k)
                if (Parent(SibID(i,SB,k), 3-k)<0 .and. &
                  .not. ANY(ParTmp == Parent(SibID(i,SB,k), 3-k))) then
                    npt = npt + 1
                    if (npt > 2) then
                        MaybeOpp = 0
                        exit
                    else
                        ParTmp(npt) = Parent(SibID(i,SB,k), 3-k)
                    endif
                endif
            enddo
            if (MaybeOpp == 1 .and. npt==2) then
                call CalcU(ParTmp(1), 3-k, ParTmp(2), 3-k, LLM(1))
                call MergeSibs(-ParTmp(1), -ParTmp(2), 3-k, LLM(2))
                if (LLM(2) - LLM(1) < -thLRrel*nS(SB,k))  MaybeOpp = 0
            endif
        endif
    endif    
    if (MaybeOpp == 1 .and. Par < 0) then
        call Qadd(A, -Par, 3-k, LLM(1))  ! 2nd degree relatives vs unrelated    
        if (LLM(1) < -thLR)  MaybeOpp = 0
    endif
    if (MaybeOpp == 1) then
        LLM = 999
        if (Par < 0) then  ! may have more/fewer sibs
            call AddFS(A, -Par, 3-k,0,3-k, LLM(1))
            call AddSib(A, -Par, 3-k, LLM(2))
            call CalcU(A, 3-k, Par, 3-k, LLM(3))
        else if (Par == 0  .and. nS(SB,k)>0) then
            call PairFullSib(A, SibID(1, SB, k), LLM(1))  ! prob. not necessary todo with all sibs
            call PairHalfSib(A, SibID(1, SB, k), 3-k, LLM(2))
            call CalcU(A, 3-k, SibID(1, SB, k), 3-k, LLM(3))
        endif
        if (LLM(1) < 0 .and. LLM(1) - LLM(2) > thLRrel) then
            LL(2) = LL(2)   ! proceed.
        else if (LLM(2) < 0 .and. LLM(2) - LLM(3) > thLRrel) then
            LL(2:3) = 222  ! more likely to be added to opp. parent, or unclear
        endif
    endif
endif

if (focal==4 .and. MaxLL(LL)==LL(4) .and. &
  (ABS(LLg(3)-LLg(4))<.1 .or. ABS(LLg(5)-LLg(4))<.1)) then
    LL(4) = 222 ! relies on age prior only. 
endif   ! TODO: use age prior if last round. 
end subroutine CheckAdd

! ##############################################################################################

subroutine Qmerge(SA, SB, k, LR)
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LR
integer :: l, x, y
double precision :: PrL(nSnp), PrX(3), PrXY(3,3)

PrL = 0
do l=1,nSnp
    do x=1,3
        PrX(x) = XPr(1,x,l,SA,k)*XPr(1,x,l,SB,k)* AHWE(x,l)
        do y=1,3
            PrXY(x,y) = XPr(1,x,l,SA,k)*XPr(1,y,l,SB,k)* AHWE(x,l) * AHWE(y,l)
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrX)) - LOG10(SUM(PrXY))
enddo

LR = SUM(PrL)

end subroutine Qmerge

! ##############################################################################################

subroutine CheckMerge(SA, SB, kA, kB, LL, focal) 
use Global
implicit none

integer, intent(IN) :: SA, SB, kA, kB, focal
double precision, intent(OUT) :: LL(7)
double precision :: LLtmp(2), ALR(7), LLx(6), LLz(2), LRHS, LLM(6), LUX, LLHA(3), &
  dLH(nS(SB,kB)), ALRx(6)
integer :: i, j, x, Par(2), SX(2), kx, m, MaybeOpp(2)
logical :: ShareOpp

LL = 999
ALR = 999
ShareOpp = .FALSE.
if (kA /= kB .and. focal==1) then
    LL(1) = 777
    return
endif
do i=1, nS(SA, kA)
    if (SibID(i, SA, kA)==GpID(1,SB,kB) .or. SibID(i, SA, kA)==GpID(2,SB,kB)) then
        LL(1) = 777
        exit
    endif
    do j=1, nS(SB, kB)
        if (SibID(j, SB, kB)==GpID(1,SA,kA) .or. SibID(j, SB, kB)==GpID(2,SA,kA)) then
            LL(1) = 777
            exit
        endif
        if (AgeDiff(SibID(i,SA,kA), SibID(j,SB,kB))==999) cycle
        if (AgePriorM(ABS(AgeDiff(SibID(i,SA,kA), SibID(j,SB,kB)))+1,kA) == 0.0) then   
            LL(1) = 777
            exit
        endif
        if (LL(1)==777) exit
        if (kA==kB) then
            if (Parent(SibID(i, SA, kA), 3-kA) == Parent(SibID(j, SB, kB), 3-kB) .and. &
              Parent(SibID(i, SA, kA), 3-kA) /=0) then
                ShareOpp = .TRUE.
            endif
        endif
    enddo
enddo 
if (LL(1) == 777 .and. focal==1) return

if (.not. ShareOpp) then
    call Qmerge(SA, SB, kB,  LRHS)
    if (LRHS < -thLR*MIN(nS(SA,kA), nS(SB,kB))) then
        LL(1) = 777
    endif
    if (LL(1) == 777 .and. focal==1) return
endif

call CalcU(-SA,kA, -SB,kB, LL(7))

if (LL(1)/=777 .and. kA==kB) then
    call MergeSibs(SA, SB, kA, LL(1))   ! SB parent of A's
    if (focal==1 .and. (LL(1) > 0 .or. LL(1) - LL(7) < thLRrel)) return
    call CalcALRmerge(SA, SB, kA, ALR(1))
    if (ALR(1) /=777 .and. LL(1) < 0) then  
        LL(1) = LL(1) + ALR(1)
    else
        LL(1) = 777
    endif
else
    LL(1) = 777
endif
if (focal==1 .and. (LL(1)==777 .or. LL(1) - LL(7) < thLRrel)) return

 call addFS(0, SA, kA, SB, kB, LL(2))  ! SB FS with an A
 call PairUA(-SB, -SA, kB, kA, LL(3))  ! SB HS with an A
if (LL(2) < 0 .or. LL(3) < 0) then
    call CalcAgeLR(-SB,kB, -SA,kA, 0,1, ALR(2))
endif
if (ALR(2) /= 777) then
    do x=2,3
        if (LL(x) < 0) then
            LL(x) = LL(x) + ALR(2)
        endif
    enddo
else
    LL(2:3) = 777
endif
 
LLtmp = 999
 call addFS(0, SB, kB, SA, kA, LLtmp(1))   ! SB GP of A's
 call PairUA(-SA, -SB, kA, kB, LLtmp(2))  ! SB GP of A's
 LL(4) = MaxLL(LLtmp)
call CalcAgeLR(-SA,kA, -SB,kB, 0,1, ALR(4))
if (ALR(4)/=777) then
    if (LL(4) < 0) then
        LL(4) = LL(4) + ALR(4)
    endif
else
    LL(4) = 777
endif
  
 call ParentHFS(0, SA, kA, SB, kB,3, LL(5))  ! SB FA of A's 
! TODO: PairUA for FS clusters

LLx = 999
do x=1,2
    call ParentHFS(0, SA, kA, SB, kB, x, LLx(x))  !SB HA of A's
enddo
if (nAgeClasses > 2) then
    call dummyGP(SA, SB, kA, kB, LLx(3))  ! SB GGP of A's
    call dummyGP(SB, SA, kB, kA, LLx(4))  ! SA GGP of B's
endif
do x=1,4
    if (LLX(x)<0) then
        if (x==1 .or. x==2) then
            call CalcAgeLR(-SA,kA, -SB,kB, x,3, ALRx(x))
        else if (x==3) then
            call CalcAgeLR(-SA,kA, -SB,kB, 3-kB,4, ALRx(x))  ! TODO? mat & pat gp
        else if (x==4) then
            call CalcAgeLR(-SB,kB, -SA,kA, 3-kA,4, ALRx(x))
        endif
        if (ALRx(x)/=777) then
            LLX(x) = LLX(x) + ALRx(x)
        else
            LLX(x) = 777
        endif
    endif
enddo
LLz = 999
do x=1,2
    if (GpID(x, SA, kA) > 0) then   ! TODO: more general
        if (Parent(GpID(x, SA, kA), kB)==-SB) then
            LLz(x) = 777
        else
            call PairUA(-SB, GpID(x, SA, kA), kB, 3, LLz(x))
        endif
        if (LLz(x) < 0) then
            LLx(4+x) = LLz(x) + CLL(SA, kA) - Lind(GpID(x, SA, kA))  ! approx.
        endif
    endif
enddo
LL(6) = MaxLL(LLx)  ! most likely 3rd degree relative

if (LL(4)<0 .and. LLtmp(1)<0 .and. LLtmp(1)>LLtmp(2) .and. &
  LLtmp(1) > MaxLL((/LL(1:3), LL(5:7)/))) then  !LLtmp1: SA FS with Bi. check if not HS via opp par
    dLH = 999
    do j=1, nS(SB, kB)
        if (GpID(3-kB,SA,kB)==0) then
            shareOpp = .TRUE.
            par(1) = Parent(SibID(j, SB, kB), 3-kB)
        else if (Parent(SibID(j, SB, kB), 3-kB)==GpID(3-kB,SA,kB) .or. &
          Parent(SibID(j, SB, kB), 3-kB)==0) then
            shareOpp = .TRUE.
            par(1) = GpID(3-kB,SA,kB)
        else
            shareOpp = .FALSE.
        endif
        if (shareOpp) then
            call PairUA(-SA, SibID(j, SB, kB), kA, 3, LLHA(1))
            call PairUA(-SA, SibID(j, SB, kB), kA, 3-kB, LLHA(2))
            call CalcU(-SA, kA, SibID(j, SB, kB), kB, LLHA(3))
            if (LLHA(1)<0)  dLH(j) = LLHA(1) - MaxLL(LLHA(2:3))
        endif
    enddo
    if (MAXVAL(dLH, MASK=dLH<999) < thLRrel) then
        LL(4) = LLtmp(2) + ALR(4)
    endif
endif
LLM = 999
Par = 0
if (kA == kB .and. focal==1 .and. ABS(MaxLL(LL) - LL(1))<thLRrel) then    
! check if they're all FS, in which case merging might also/instead go via 3-k
    SX = (/SA, SB/)
    kx = kA  !=kB
    MaybeOpp = 0
    do m=1,2
        call getFSpar(SX(m), kx, Par(m))
        if (Par(m)>0) cycle
        if (Par(m)==0 .and. ANY(Parent(SibID(1:nS(SX(m),kx),SX(m),kx),3-kx)>0)) cycle
        MaybeOpp(m) = 1
    enddo
    if (MaybeOpp(1)==1 .and. MaybeOpp(2)==1) then   
        if (Par(1)==Par(2) .and. Par(1)/=0) then ! .and. nS(-Par(1),3-kx)==nS(SA,kA)+nS(SB,kB)) then
            MaybeOpp = 0
        else if (Par(1)<0 .and. Par(2)<0) then
            call CalcALRmerge(-Par(1), -Par(2), 3-kA, ALR(7))
            if (ALR(7)==777) MaybeOpp = 0
        endif
    endif
    if (MaybeOpp(1)==1 .and. MaybeOpp(2)==1) then
        call FSMerge(SA,SB,kA, LLM)  ! LLM(1): not, 2:k, 3:3-k, 4:FS
        LLM(2) = MaxLL((/LLM(2), LL(1)-ALR(1)/))  ! merge via k
        LLM(1) = MaxLL((/LLM(1), LL(7)/))  ! do not merge
        if (par(1)<0 .and. par(2)<0) then
            call MergeSibs(-par(1), -par(2), 3-kA, LLM(3))  ! may include indivs not in SA and SB
            call CalcU(Par(1), 3-kA, Par(2), 3-kA, LLM(1))
        endif
        if (MaxLL(LLM)==LLM(4) .and. LLM(4)-LLM(2) > thLRrel*MIN(nS(SA,kA),nS(SB,kB))) then  
            LL(1) = LLM(4) + ALR(1)  ! FS merge most likely - go ahead.
        else !if (LLM(3)<0 .and. LLM(3)-LLM(1) > thLRrel) then
 !           if (LLM(3)-LLM(4) > thLRrel) then
                LL(1) = 222 ! more likely that opp. parent need to be merged, or unclear
 !           else if (Par(1) < 0 .and. Par(2)<0) then
 !!               if (nS(-Par(1),3-kx)+ns(-Par(2),3-kx) > nS(SA,kA)+nS(SB,kB)) then
 !                   LL(1) = 222   ! LLM(3) and (4) not comparable
 !               endif
 !           endif
        endif
    endif
    ! consider merge par, & SA,SB are PO (pairUA). identical LL if complete overlap
    if (LL(1) /= 222 .and.  par(1)<0 .and. par(2)<0) then
        if (nS(SA,kA) + nS(SB,kB) == nS(-Par(1),3-kx)+ns(-Par(2),3-kx)) then  
            if (GpID(kA,SB,kB)==0) then
                GpID(kA,SB,kB) = -SA   ! temp assignment
                call calcCLL(SB, kB)
                call MergeSibs(-par(1), -par(2), 3-kx, LLM(5)) 
                GpID(kA,SB,kB) = 0
                call calcCLL(SB, kB)
                if (LLM(5) - LLM(4) > thLRrel) then 
                    LL(1) = 222
                endif
            endif
            if (GpID(kB,SA,kA)==0) then
                GpID(kB,SA,kA) = -SB   ! temp assignment
                call calcCLL(SA, kA)
                call MergeSibs(-par(1), -par(2), 3-kx, LLM(6)) 
                GpID(kB,SA,kA) = 0
                call calcCLL(SA, kA)
                if (LLM(6) - LLM(4) > thLRrel) then 
                    LL(1) = 222
                endif
            endif
        endif
    endif
endif

end subroutine CheckMerge

! ##############################################################################################

subroutine getFSpar(SA, kA, par)  ! all individuals in SA are FS to eachother
use Global
implicit none

integer, intent(IN) :: SA,  kA
integer, intent(OUT) :: Par
integer :: i, j

Par = 0
do i=1, nS(SA,kA)
    if (Parent(SibID(i,SA,kA), 3-kA)/=0) then
        Par = Parent(SibID(i,SA,kA), 3-kA)
        do j= i+1, nS(SA, kA)
            if (Parent(SibID(j,SA,kA), 3-kA) /= Par .and. &
              Parent(SibID(j,SA,kA), 3-kA)/=0) then
                Par = 0
                return
            endif
        enddo
    endif
enddo

end subroutine getFSpar

! ##############################################################################################

subroutine FSmerge(SA,SB,k, LL)  ! calc LL if SA and SB merged via both pat & mat
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LL(4)
integer :: l, x, y, i, u,v, G(2),z
double precision :: PrL(nSnp,4), PrXY(3,3), PrUV(3,3), PrXV(3,3,3,3,4), &
    PrG(3,2), PrX(3), PrTmp(3)

! TODO: currently assumes no grandparents of sibship 3-k, no close inbreeding
! check done in CheckMerge

LL = 999
do i=1,2
    if (GpID(i,SA,k)/=0) then
        if(GpID(i,SA,k)/=GpID(i,SB,k) .and. GpID(i,SB,k)/=0) then
            G(i) = 0  ! shouldn't happen
        else
            G(i) = GpID(i,SA,k)
        endif
    else
        G(i) = GpID(i,SB,k)
    endif
enddo

PrL = 0
do l=1,nSnp 
    do i=1,2
        call ParProb(l, G(i), k, 0,0, PrG(:,i))
    enddo
    do x=1,3
        do z=1,3
            PrTmp(z) = SUM(AKA2P(x,:,z) * PrG(:,1) * PrG(z,2))
        enddo
        PrX(x) = SUM(PrTmp)
    enddo

    do x=1,3  ! P1
        do y=1,3  ! P2
            PrXY(x,y) = 1  ! XPr(2,x,l, sA,k) * AHWE(y,l)
            do i=1,nS(SA,k)
                if (Genos(l, SibID(i,SA,k))/=-9) then
                    PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l, SibID(i,SA,k)), x, y)
                endif
            enddo
        enddo
    enddo
    do u=1,3
        do v=1,3
            PrUV(u,v) = 1  ! XPr(2,u,l, sB,k) * AHWE(v,l)
            do i=1,nS(SB,k)
                if (Genos(l, SibID(i,SB,k))/=-9) then
                    PrUV(u,v) = PrUV(u,v) * OKA2P(Genos(l, SibID(i,SB,k)), u, v)
                endif
            enddo
        enddo
    enddo

    PrXV = 0
    do x=1,3 
        do y=1,3
            do u=1,3
                do v=1,3
                    PrXV(x,y,u,v,1) = PrXY(x,y) * XPr(2,x,l, sA,k) * AHWE(y,l) * &
                      PrUV(u,v) * XPr(2,u,l, sB,k) * AHWE(v,l)
                    PrXV(x,y,x,v,2) = PrXY(x,y) * PrX(x) * AHWE(y,l) * &
                      PrUV(x,v) * AHWE(v,l)
                enddo
                PrXV(x,y,u,y,3) = PrXY(x,y) * XPr(2,x,l, sA,k) * AHWE(y,l) * &
                  PrUV(u,y) * XPr(2,u,l, sB,k)
            enddo
            PrXV(x,y,x,y,4) = PrXY(x,y) * PrX(x) * AHWE(y,l) * PrUV(x,y)
        enddo
    enddo
    do x=1,4
        PrL(l,x) = LOG10(SUM(PrXV(:,:,:,:,x)))
    enddo
enddo
LL = SUM(PrL,DIM=1)

end subroutine FSmerge

! ##############################################################################################

subroutine MakeFS(A, B)
use Global
implicit none

integer, intent(IN) :: A,B
integer :: x, i, j, Ai, Bj

Ai = 0
Bj = 0
if (nFS(A)>0) then
    Ai = A
else if (ANY(FSID(1:nFS(B),B)==A)) then
    return ! already are FS.
else
    do i=1, nInd
        if (Parent(i,1)==Parent(A,1) .and. Parent(i,2)==Parent(A,2) .and. nFS(i)>0) then
            Ai = i
            exit
        endif
    enddo
endif
if (nFS(B)>0) then
    Bj = B
else if (ANY(FSID(1:nFS(Ai),Ai)==B)) then
    return 
else
    do i=1, nInd
        if (Parent(i,1)==Parent(B,1) .and. Parent(i,2)==Parent(B,2) .and. nFS(i)>0) then
            Bj = i
            exit
        endif
    enddo
endif

i = MIN(Ai,Bj)
j = MAX(Ai,Bj)
do x=1, nFS(j)   
    FSID(nFS(i)+x, i) = FSID(x, j)
enddo
nFS(i) = nFS(i) + nFS(j)
FSID(:,j) = 0
FSID(1,j) = j
nFS(j) = 0

end subroutine MakeFS

! ##############################################################################################

subroutine DoAdd(A, SB, k)
use Global
implicit none

integer, intent(IN) :: A, SB, k
integer :: i, n

if (nS(SB,k) +1 >= maxSibSize) then
    call rexit("Reached maxSibSize ")
endif

do n=1, nS(SB,k)  ! check for FS
    i = SibID(n,SB,k)
    if (nFS(i)==0) cycle
    if (Parent(A, 3-k)/=0 .and. Parent(A, 3-k)==Parent(i, 3-k)) then
        call MakeFS(A, i)
    endif
enddo

SibID(nS(SB,k)+1, SB, k) = A  ! add A to sibship
Parent(A, k) = -SB
nS(SB,k) = nS(SB,k) + 1
 call calcCLL(SB,k)

do n=1,nS(SB,k)   ! update LL of connected sibships
    i = SibID(n,SB,k)
    if (Parent(i,3-k) < 0) then
        call CalcCLL(-Parent(i,3-k), 3-k)
    endif                    
    call CalcLind(i)
enddo
 call calcCLL(SB,k)
 call CalcLind(A)
 
end subroutine DoAdd

! ##############################################################################################

subroutine DoMerge(SA, SB, k)
use Global
implicit none

integer, intent(IN) :: SA, SB, k
integer :: i, j, n, m, x

if (SA/=0) then
    if (nS(SA,k) + nS(SB,k) >= maxSibSize) then
        call rexit("reached maxSibSize")
    endif

    do n=1,nS(SA,k)    ! check for FS
        i = SibID(n,SA,k)
        if (nFS(i)==0) cycle
        do m=1, nS(SB,k)
            j = SibID(m,SB,k)
            if (nFS(j)==0) cycle
            if (Parent(i, 3-k)/=0 .and. Parent(i, 3-k)==Parent(j, 3-k)) then
                call MakeFS(i, j)
            endif
        enddo
    enddo
    do m=1,nS(SB,k)  ! add sibship SB to SA
        SibID(nS(SA,k)+m, SA, k) = SibID(m, SB, k)
        Parent(SibID(m, SB, k), k) = -SA
        do i=1,2
            if (GpID(i, SA, k)==0 .and. GpID(i, SB, k)/=0) then
                GpID(i, SA, k) = GpID(i, SB, k)  ! checked for mismatches earlier
            endif
        enddo
    enddo
    nS(SA,k) = nS(SA,k) + nS(SB,k)
    
    call calcCLL(SA,k)
    do n=1,nS(SA,k) 
        i = SibID(n,SA,k)
        if (Parent(i,3-k) < 0) then
            call CalcCLL(-Parent(i,3-k), 3-k)
        endif                    
        call CalcLind(i)
    enddo
endif
 
do x=SB, nC(k)-1  !remove cluster SB, shift all subsequent ones
    SibID(:, x, k) = SibID(:, x+1, k)
    nS(x, k) = nS(x+1, k)
    GpID(:, x,k) = GpID(:, x+1,k)
    do n=1, nS(x,k)
        Parent(SibID(n,x,k),k) = -x  !Parent(SibID(n,x,k),k)+1 ! shift towards zero.
    enddo
    CLL(x,k) = CLL(x+1, k)
    XPr(:,:,:,x,k) = XPr(:,:,:,x+1,k)
    DumP(:,:,x,k) = DumP(:,:,x+1,k)
    DumBY(:,x,k) = DumBY(:,x+1,k)
enddo
SibID(:,nC(k),k) = 0
GpID(:,nC(k),k) = 0
nS(nC(k), k) = 0

do x=SB, nC(k)  !fix GPs
    do m=1,2
        do n=1, nC(m)
            if (GpID(k, n, m) == -x)  GpID(k, n, m) = -x+1 
        enddo
    enddo
enddo
nC(k) = nC(k) -1

end subroutine DoMerge

! ##############################################################################################

subroutine CalcU(A, kAIN, B, kBIN, LL)  ! A, SB, k, SA, LL
! TODO: gives incorrect result for A>0, kA==kB, B<0, parent(A,3-kA)=GB
! TODO: integrate out per opposite parent, as in CalcCLL
use Global
implicit none

integer, intent(IN) :: A, kAIN, B, kBIN
double precision, intent(OUT) :: LL
integer :: m, n, cat, par(2), Ai, Bj, SA, SB, kA, kB
logical :: swap, con

LL = 999
if (A>0) then
    call CalcLind(A)
else if (A<0) then
    kA = kAIN
    call CalcCLL(-A, kA)
endif
if (B>0) then
    call CalcLind(B)
else if (B<0) then
    kB = kBIN
    call CalcCLL(-B, kB)
endif
!==================================

if (A==0) then
    if (B==0) then
        LL = 0
    else if (B>0) then
        LL = Lind(B)
    else if (B<0) then
        LL = CLL(-B, kB)
    endif
    return
else if (B==0) then
    if (A>0) then
        LL = Lind(A)
    else if (A<0) then
        LL = CLL(-A,kA)
    endif
    return
else if (A>0 .and. B<0) then
    if (Parent(A,kB)==B) then
        LL = CLL(-B,kB)
        return
    else if (ANY(GpID(:,-B,kB) == A)) then
        LL = CLL(-B,kB) + Lind(A)  ! CLL already conditional on A
        return
    else if (ALL(Parent(A,:)>=0)) then  
        LL = Lind(A) + CLL(-B, kB)
        return
    else
        call Connected(A,1,-B,kB, con)
        if (.not. con) then
            LL = Lind(A) + CLL(-B, kB)
            return
        endif
    endif
else if (B>0 .and. A<0) then
    if (Parent(B,kA)==A) then
        LL = CLL(-A, kA)
        return
    else if (ANY(GpID(:,-A,kA) == B)) then
        LL = CLL(-A,kA) + Lind(B)  
        return
    else if (ALL(Parent(B,:)>=0)) then 
        LL = CLL(-A, kA) + Lind(B) 
        return
    else
        call Connected(B,1,-A,kA, con)
        if (.not. con) then
            LL = CLL(-A, kA) + Lind(B) 
            return
        endif
    endif
endif

!==================================
! determine relationship between focal individuals/clusters

Ai = 0
Bj = 0
SA = 0
SB = 0
cat = 0

if (A>0 .and. B>0) then  ! == pairs ==
    do m=1,2
        if (Parent(A,m)/=0 .and. Parent(A,m)==Parent(B,m)) then
           par(m) = Parent(A,m)
        else
            par(m) = 0  ! unknown or unequal
        endif
    enddo

    if (par(1)/=0 .and. par(2)/=0) then
        cat = 2  ! FS
    else if (par(1)/=0 .or. par(2)/=0) then
        cat = 3  ! HS
    else 
        do m=1,2
            if (parent(A,m) < 0) then
                do n=1,2
                    if (GpID(n, -parent(A,m), m) == B) then
                        cat = 0 !4  ! already conditioned on.
                    else
                        if (GpID(n, -parent(A,m), m) == Parent(B, n) .and. Parent(B, n)<0) then
                            cat = 5
                        endif
                    endif
                enddo
            else if (parent(B,m) < 0) then
                if (ANY(GpID(:, -parent(B,m), m) == A)) then
                    cat = 0 !4
                    swap = .TRUE.
                else
                    do n=1,2
                        if (GpID(n, -parent(B,m), m) == Parent(A, n) .and. Parent(A, n)<0) then
                            cat = 5
                            swap = .TRUE.
                        endif
                    enddo
                endif
            endif
        enddo
    endif

    if (cat==0 .or. cat==5) then  ! TODO? cat=5
        LL = Lind(A) + Lind(B)
        return
    else if (cat==2 .and. par(1)<0 .and. par(2)<0) then
        Ai = A
        Bj = B
        SA = -par(1)
        kA = 1
        SB = -par(2)
        kB = 2
        cat=0
    else
        call Upair(A, B, cat, LL)
        return
    endif

else if (A>0 .and. B<0) then
    SB = -B
    Ai = A
    if (SA==0) then
        if (Parent(A,3-kB) < 0) then
            SA = -Parent(A,3-kB)
            kA = 3-kB    
        else if (Parent(A,kB) < 0) then
            SA = -Parent(A,kB)
            kA = kB
        endif ! else: Lind + CLL (earlier)
    endif  
else if (B>0 .and. A<0) then
    SA = -A
    Bj = B
    if (SB==0) then
        if (Parent(B,3-kA) < 0) then
            SB = -Parent(B,3-kA)
            kB = 3-kA
        else if (Parent(B,kA) < 0) then
            SB = -Parent(B,kA)
            kB = kA 
        endif
    endif

else if (A<0 .and. B<0) then
    SA = -A  
    SB = -B
endif

cat = 0
swap = .FALSE.
if (GpID(kB, SA, kA) == -SB) then
    cat = 1  ! PO
else if (GpID(kA, SB, kB) == -SA) then
    cat = 1
    swap = .TRUE.
else 
    do m=1,2
        if (GpID(m, SA, kA)==GpID(m, SB, kB) .and. GpID(m, SA, kA)/=0) then
            if (GpID(3-m, SA, kA)==GpID(3-m, SB, kB) .and. GpID(3-m, SA, kA)/=0) then  
                cat = 2  ! FS
            else
                cat = 3  ! HS
            endif
        else 
            if (GpID(m, SA, kA)<0) then
                if (GpID(kB, -GpID(m, SA, kA), m) == -SB) then
                    cat = 4  ! GP
                endif
            endif
            if (GpID(m, SB, kB)<0) then
                if (GpID(kA, -GpID(m, SB, kB),m) == -SA) then
                    cat = 4
                    swap = .TRUE.
                endif
            endif
        endif
    enddo  ! FA between SA, SB not currently considered.
endif

if (.not. swap) then
    call UClust(-SA, -SB, kA, kB, cat, Ai, Bj, LL)
else
    call UClust(-SB, -SA, kB, kA, cat, Bj, Ai, LL)
endif

end subroutine CalcU

! ##############################################################################################

subroutine Upair(A, B, cat, LL)
use Global
implicit none

integer, intent(IN) :: A, B, cat
double precision, intent(OUT) :: LL
integer :: m, l, n, x, y, par(2)
double precision :: PrL(nSnp), PrP(3,2), PrPA(3), PrPB(3), PrXY(3,3)

LL = 999
call CalcLind(A)
call CalcLind(B)

do m=1,2
    if (Parent(A,m)/=0 .and. Parent(A,m)==Parent(B,m)) then
       par(m) = Parent(A,m)
    else
        par(m) = 0  ! unknown or unequal
    endif
enddo
          
PrL = 0
do l=1, nSnp
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) then
        cycle
    else if (Genos(l,A)==-9) then
        PrL(l) = LindX(l,B)
        cycle
    else if (Genos(l,B)==-9) then
        PrL(l) = LindX(l,A)
        cycle
    endif
    
    if (cat==2) then
        do m=1,2
            call ParProb(l, Par(m), m, A, B, PrP(:,m))
        enddo   
        do x=1,3
            do y=1,3
                PrXY(x,y) = OKA2P(Genos(l,A),x,y) * OKA2P(Genos(l,B),x,y) * &
                  PrP(x,1) * PrP(y,2)
            enddo
        enddo
    else if (cat==3) then  ! HS
        do m=1,2
            if (Par(m)==0) cycle
            call ParProb(l, Par(m), m, A, B, PrP(:,m))
            call ParProb(l, Parent(A, 3-m), 3-m, A, 0, PrPA)
            call ParProb(l, Parent(B, 3-m), 3-m, B, 0, PrPB)
            do x=1,3  ! shared parent
                do y=1,3  ! parent A
                    PrXY(x,y) = OKA2P(Genos(l,A),x,y) * PrP(x,m) * PrPA(y) * &
                       SUM(OKA2P(Genos(l,B),x,:) * PrPB)
                enddo
            enddo
        enddo
    else if (cat==4) then
        do m=1,2
            if (Parent(A,m)<0) then
                do n=1,2
                    if (GpID(n, -parent(A,m), m) == B) then
                        call ParProb(l, parent(A,m), m, A, -4, PrP(:,m))  ! excl A, excl GPs
                        call ParProb(l, parent(A,3-m), 3-m, A, 0, PrPA)
                        call ParProb(l, GpID(3-n, -parent(A,m), m), 3-n, 0, 0, PrPB)
                        call ParProb(l, B, n, 0, 0, PrP(:,3-m))
                        do x=1,3  ! in-between parent
                            do y=1,3  ! other parent of A
                                PrXY(x,y) = SUM(OKA2P(Genos(l,A),x,:) * PrPA) * PrP(x,m) * &
                                   SUM(AKA2P(x, y,:) * PrP(y,3-m) * PrPB)
                            enddo
                        enddo
                    endif
                enddo
            endif
        enddo
 !   else if (cat==5) then
        
    endif
    PrL(l) = LOG10(SUM(PrXY))
enddo

LL = SUM(PrL)

end subroutine Upair

! ##############################################################################################

subroutine UClust(A, B, kA, kB, cat, Ai, Bj, LL)
use Global
implicit none

integer, intent(IN) :: A, B, kA, kB, cat, Ai, Bj
double precision, intent(OUT) :: LL
integer, allocatable, dimension(:) :: AA, BB
integer :: nA, nB, l,x,y, v, i, j, z,m, f, e, u, &
  DoneA(maxSibSize), catA(maxSibSize), catB(maxSibSize), GA(2), GB(2), g, Ei
double precision :: PrL(nSnp,2), PrGA(3,2), PrGB(3,2), PrGGP(3), &
  PrUZ(3,3, 3,3,3,3,2), PrE(3), PrH(3)

call CalcCLL(-A, kA)
call CalcCLL(-B, kB)
LL = 999

nA = nS(-A, kA)
allocate(AA(nA))
AA = SibID(1:nA, -A, kA)
GA = GpID(:, -A, kA)

nB = nS(-B, kB)
allocate(BB(nB))
BB = SibID(1:nB, -B, kB)
GB = GpID(:, -B, kB)

catA = 0
catB = 0
do i = 1, nA
    do j = 1, nB
        if (kA /= kB) then
            if (Parent(AA(i), kB) == B) then
                catA(i) = 1
            endif
            if (Parent(BB(j), kA) == A) then
                catB(j) = 1
            endif
        else if (kA == kB) then
            if (Parent(AA(i), 3-kA) == Parent(BB(j), 3-kB) .and. Parent(BB(j), 3-kB)<0) then  ! TODO: /=0    
                catA(i) = 7
                catB(j) = 7
            endif
        endif
        do m=1,2
            if (GA(m) < 0) then
                if (GA(m) == Parent(AA(i), m) .and. m/=kA) then
                    catA(i) = 2
                endif
                if (GA(m) == Parent(BB(j), m) .and. m/=kB) then
                    catB(j) = 2
                endif
            endif
            if (GB(m) < 0) then
                if (GB(m) == Parent(AA(i), m) .and. m/=kA) then
                    catA(i) = 3
                endif
                if (GB(m) == Parent(BB(j), m) .and. m/=kB) then
                    catB(j) = 3  
                endif
            endif
        enddo
        if (Parent(AA(i), 3-kA) < 0) then
            if (GpID(kB, -Parent(AA(i), 3-kA), 3-kA)==B) then
                catA(i) = 6
            endif
        endif
        if (Parent(BB(j), 3-kB) < 0) then
            if (GpID(kA, -Parent(BB(j),3-kB), 3-kB)==A) then
                catB(j) = 6
            endif 
        endif
    enddo
enddo

!==================================
if (cat==0 .and. ALL(catA==0) .and. ALL(CatB==0)) then
    if (A<0 .and.B>0) then
        LL = CLL(-A,kA) + Lind(B)
    else if (A>0 .and. B<0) then
        LL = Lind(A) + CLL(-B,kB)
    else if (A<0 .and. B<0) then
        LL = CLL(-A,kA) + CLL(-B,kB)
    endif
    return
endif

!==================================

PrL = 0
do l=1, nSnp
    PrUZ = 0
    
    do m=1, 2
        if ((ANY(catA==2) .and. m/=kA) .or. (ANY(catB==2) .and. m/=kB))  then
            call ParProb(l, GA(m), m, -1, 0, PrGA(:, m))
        else
            call ParProb(l, GA(m), m, 0, 0, PrGA(:, m))
        endif
        if ((ANY(catA==3) .and. m/=kA) .or. (ANY(catB==3) .and. m/=kB))  then
            call ParProb(l, GB(m), m, -1, 0, PrGB(:, m))
        else
            call ParProb(l, GB(m), m, 0, 0, PrGB(:, m))
        endif
    enddo
    if (cat==4) then
        do m=1,2
            if (GA(m)<0) then
                if (GpID(kB, -GA(m), m) == B) then
                    call ParProb(l, GpID(3-kB, -GA(m), m), 3-kB, 0, 0, PrGGP) 
                endif
            endif
        enddo
    endif
    
    ! == grandparents ==
    do x=1,3 
        do y=1,3
            do u=1,3  ! GP A, kB
                do z=1,3  ! GP A, 3-kB
                    do v=1,3  ! GP B, kB
                        if (cat == 1) then
                            if (kA==kB .and. GA(3-kB)==GB(3-kB) .and. GA(3-kB)/=0) then  ! inbreeding loop
                                PrUZ(x,y,y,z,v,z,1) = AKA2P(x,y,z) * AKA2P(y,v,z) * PrGA(z,3-kB) * PrGB(v,kB)
                            else
                                PrUZ(x,y,y,z,v,:,1) = AKA2P(x,y,z) * AKA2P(y,v,:) * &
                                  PrGA(z,3-kB) * PrGB(v,kB) * PrGB(:,3-kB)
                            endif
                        else if (cat==2) then
                            PrUZ(x,y,u,z,u,z,1) = AKA2P(x,u,z) * AKA2P(y,u,z) * PrGA(u,kB) * PrGA(z,3-kB)
                        else if (cat==3) then
                            do m=1,2
                                if (GA(m)/=0 .and. GA(m) == GB(m)) then
                                    if (m==kB) then
                                        PrUZ(x,y,u,z,u,:,1) = AKA2P(x,u,z) * AKA2P(y,u,:) * &
                                           PrGA(u,m) * PrGA(z,3-m) * PrGB(:,3-m)   
                                    else
                                        PrUZ(x,y,u,z,v,z,1) = AKA2P(x,u,z) * AKA2P(y,v,z) * &
                                           PrGA(u,3-m) * PrGA(z,m) * PrGB(v,3-m)
                                    endif
                                endif
                            enddo
                        else if (cat==4) then 
                            do m=1,2
                                if (GA(m)<0) then
                                    if (GpID(kB, -GA(m), m) == B) then
                                        if (m==kB) then
                                            PrUZ(x,y,u,z,v,:,1) = AKA2P(x,u,z) * SUM(AKA2P(u,y,:) * PrGGP) * &
                                              PrGA(z,3-kB) * AKA2P(y,v,:) * PrGB(v,kB) * PrGB(:, 3-kB) 
                                        else
                                            PrUZ(x,y,u,z,v,:,1) = AKA2P(x,u,z) * SUM(AKA2P(z,y,:) * PrGGP) * &
                                              PrGA(u,kB) * AKA2P(y,v,:) * PrGB(v,kB) * PrGB(:, 3-kB)
                                        endif
                                    endif
                                endif
                            enddo
                        else
                            PrUZ(x,y,u,z,v,:,1) = AKA2P(x,u,z) * AKA2P(y,v,:) * &
                              PrGA(u,kB) * PrGA(z,3-kB) * PrGB(v,kB) * PrGB(:, 3-kB)
                        endif
                        PrUZ(x,y,u,z,v,:,2) = PrUZ(x,y,u,z,v,:,1)
                    enddo
                enddo
            enddo
        enddo
    enddo
   
    ! == siblings ==   
    if (ALL(catA==0) .and. ALL(catB==0) .and. Ai==0 .and. Bj==0) then
        do x=1,3 
            do y=1,3
                PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * XPr(1,x,l,-A,kA) * XPr(1,y,l,-B,kB)
            enddo  ! TODO: needs special for cat<4 ?
        enddo
    
    else

    do x=1,3  ! SA
        doneA = 0
        do y=1,3  ! SB
            do j=1, nB
                if (nFS(BB(j))==0) cycle
                if (catB(j)==1 .or. catB(j)==2 .or. catB(j)==3) then
                    PrE = 1
                else if (catB(j)==6) then   ! GpID(kA, -Parent(BB(j),3-kB), 3-kB)==PA
                    call ParProb(l, GpID(3-kA, -Parent(BB(j), 3-kB), 3-kB), 3-kA, 0, 0, PrH)
                    do e=1,3
                        PrE(e) = SUM(AKA2P(e,x,:) * PrH) 
                    enddo
                else
                    call ParProb(l, Parent(BB(j),3-kB), 3-kB, -1, 0, PrE)
                endif
                
                if (Parent(BB(j),3-kB) < 0) then 
                    do e=1,3
                        do g=1, nS(-Parent(BB(j),3-kB), 3-kB)
                            Ei = SibID(g, -Parent(BB(j),3-kB), 3-kB)
                            if (nFS(Ei) == 0) cycle
                            if (Parent(Ei, kB) == B) cycle  
                            if (Parent(Ei, kA) == A) cycle 
!                            if (CatB(j)==2 .or. CatB(j)==3) cycle  ! TODO: or cat=4 special case
                            call ParProb(l, Parent(Ei, kB), kB, Ei, -1, PrH) 
                            do i=1, nFS(Ei)
                                if (Genos(l,FSID(i, Ei))==-9) cycle
                                PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                            enddo
                            PrE(e) = PrE(e) * SUM(PrH)
                        enddo
                    enddo
                endif
                
                if (Bj/=0 .or. catB(j)==7) then
                    do f=1, nFS(BB(j))
                        if (Bj==0 .or. FSID(f, BB(j))==Bj) cycle
                        if (Parent(FSID(f, BB(j)),kA)==A .and. (Ai==0 .or. FSID(f, BB(j))==Ai)) cycle
                        if (Genos(l,FSID(f, BB(j)))==-9) cycle
                        PrE = PrE * OKA2P(Genos(l,FSID(f,BB(j))), y, :)
                    enddo
                endif
                
                if (catB(j)==7 .and. Ai/=0) then 
                    do i=1,nA
                        if (Parent(AA(i), kB) == B) cycle
                        do f=1, nFS(AA(i))
                            if (FSID(f, AA(i))==Ai) cycle
                            if (Genos(l,FSID(f, AA(i)))==-9) cycle
                            PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
                        enddo
                    enddo
                endif
                
                if (catB(j)==1) then  ! Parent(BB(j), 3-kB)==PA
                    PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * PrE(x)
                else if (CatB(j)==2) then
                    do z=1,3
                        PrUZ(x,y,:,z,:,:,1) = PrUZ(x,y,:,z,:,:,1) * PrE(z)   
                    enddo
                else if (CatB(j)==3) then
                    do z=1,3
                        PrUZ(x,y,:,:,:,z,1) = PrUZ(x,y,:,:,:,z,1) * PrE(z)   
                    enddo                       
                else 
                    PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * SUM(PrE)
                endif
                
                do f=1, nFS(BB(j)) ! includes some AA if cat=1 
                    if (Bj==0 .or. FSID(f, BB(j))==Bj .or. &
                      (Parent(BB(j),kA)==A .and. (Ai==0 .or. BB(j)==Ai))) then
                        if (Genos(l,FSID(f, BB(j)))==-9) cycle
                        PrE = PrE * OKA2P(Genos(l,FSID(f, BB(j))), y, :)
                    endif
                enddo

                if (catB(j)==7) then  !  Parent(BB(j), 3-kB)/=0, shared opposite parent
                    do i=1,nA
                        if (Parent(AA(i), 3-kB) /= Parent(BB(j), 3-kB)) cycle
                        if (Parent(AA(i), kB) == B) cycle
                        if (Ai==0 .and. nFS(AA(i))==0) cycle
                        do f=1, nFS(AA(i))
                            if (FSID(f, AA(i))/=Ai) cycle
                            if (Genos(l,FSID(f, AA(i)))==-9) cycle
                            PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
                            DoneA(i) = 1
                        enddo
                    enddo
                endif
                
                if (catB(j)==1) then  ! Parent(BB(j), 3-kB)==PA
                    PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * PrE(x)
                else if (CatB(j)==2) then
                    do z=1,3
                        PrUZ(x,y,:,z,:,:,2) = PrUZ(x,y,:,z,:,:,2) * PrE(z)   
                    enddo
                else if (CatB(j)==3) then
                    do z=1,3
                        PrUZ(x,y,:,:,:,z,2) = PrUZ(x,y,:,:,:,z,2) * PrE(z)   
                    enddo                       
                else 
                    PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * SUM(PrE)
                endif
                
            enddo  ! B_j
        
            do i=1, nA
                if (DoneA(i)==1) cycle
                if (nFS(AA(i))==0) cycle
                if (Parent(AA(i),kB)==B) cycle
                if (catA(i)>1 .and. catA(i)<4) then  ! catA==1 already done
                    PrE = 1
                else if (catA(i)==6) then ! Parent(AA(i), 3-kA) <0
                    call ParProb(l, GpID(3-kB, -Parent(AA(i), 3-kA), 3-kA), 3-kB, 0, 0, PrH)
                    do e=1,3
                        PrE(e) = SUM(AKA2P(e,y,:) * PrH) 
                    enddo
                else    
                    call ParProb(l, Parent(AA(i), 3-kA), 3-kA, -1, 0, PrE)   
                endif
                
                if (Parent(AA(i), 3-kA) < 0) then 
                    do e=1,3
                        do g=1, nS(-Parent(AA(i), 3-kA), 3-kA)
                            Ei = SibID(g, -Parent(AA(i), 3-kA), 3-kA)
                            if (nFS(Ei) == 0) cycle
                            if (Parent(Ei, kA) == A) cycle 
                            if (Parent(Ei, kB) == B) cycle
!                            if (CatA(i)==2 .or. CatA(i)==3) cycle  ! TODO: or cat=4 special case
                            call ParProb(l, Parent(Ei, kA), kA, Ei, -1, PrH)  ! Assume is not one of GPs
                            do j=1, nFS(Ei)                  
                                if (Genos(l,FSID(j, Ei))==-9) cycle
                                PrH = PrH * OKA2P(Genos(l,FSID(j,Ei)), :, e)
                            enddo
                            PrE(e) = PrE(e) * SUM(PrH)
                        enddo
                    enddo
                endif
                
                if (Ai/=0) then
                    do f=1, nFS(AA(i))
                        if (FSID(f, AA(i))==Ai) cycle
                        if (Genos(l,FSID(f, AA(i)))==-9) cycle
                        PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
                    enddo
                endif
                
                if (catA(i)==2) then
                    do z=1,3
                        if (kA==kB) then
                            PrUZ(x,y,:,z,:,:,1) = PrUZ(x,y,:,z,:,:,1) * PrE(z)
                        else
                            PrUZ(x,y,z,:,:,:,1) = PrUZ(x,y,z,:,:,:,1) * PrE(z)         
                        endif
                    enddo
                else if (catA(i)==3) then  !  .or. (catA(nA+1)==4 .and. Parent(AA(i), 3-kA) == GB(3-kA) .and. GB(3-kA)/=0)
                    do z=1,3
                        if (kA==kB) then
                            PrUZ(x,y,:,:,:,z,1) = PrUZ(x,y,:,:,:,z,1) * PrE(z)
                        else
                            PrUZ(x,y,:,:,z,:,1) = PrUZ(x,y,:,:,z,:,1) * PrE(z)
                        endif
                    enddo
                else
                    PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * SUM(PrE)    
                endif
                
                do f=1, nFS(AA(i)) 
                    if (Ai/=0 .and. FSID(f, AA(i))/=Ai) cycle 
                    if (Genos(l,FSID(f, AA(i)))==-9) cycle
                    PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
                    DoneA(i)=2  ! for printing only 
                enddo
                
                if (catA(i)==2) then
                    do z=1,3
                        if (kA==kB) then
                            PrUZ(x,y,:,z,:,:,2) = PrUZ(x,y,:,z,:,:,2) * PrE(z)
                        else
                            PrUZ(x,y,z,:,:,:,2) = PrUZ(x,y,z,:,:,:,2) * PrE(z)         
                        endif
                    enddo
                else if (catA(i)==3) then  !  .or. (catA(nA+1)==4 .and. Parent(AA(i), 3-kA) == GB(3-kA) .and. GB(3-kA)/=0)
                    do z=1,3
                        if (kA==kB) then
                            PrUZ(x,y,:,:,:,z,2) = PrUZ(x,y,:,:,:,z,2) * PrE(z)
                        else
                            PrUZ(x,y,:,:,z,:,2) = PrUZ(x,y,:,:,z,:,2) * PrE(z)
                        endif
                    enddo
                else
                    PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * SUM(PrE)    
                endif
            enddo  ! i
        enddo  ! x
    enddo  ! y
    endif
    do f=1,2
        PrL(l,f) = LOG10(SUM(PrUZ(:,:,:,:,:,:,f)))
    enddo
enddo

LL = SUM(PrL(:,2)) - SUM(PrL(:,1))

deallocate(AA)
deallocate(BB)

end subroutine UClust

! ##############################################################################################

subroutine AddSib(A, SB, k, LL)  ! TODO?: don't assume 3-k's remain fixed
use Global
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l, x, Bj, m, f, AncB(2,mxA), Inbr
double precision :: PrL(nSnp), PrX(3), PrY(3)

LL = 999
if (Parent(A,k)==-SB) then
    LL = 888
else if (Parent(A,k)/=0) then
    LL = 777
endif
if (LL/=999) return
do f=1, nS(SB,k)
    Bj = SibID(f, SB, k)
    if (Parent(A, 3-k) /= 0) then
        if (Parent(Bj, 3-k) == Parent(A, 3-k)) then
            LL = 777  ! use addFS() instead
        endif
    endif
    if (AgeDiff(A,Bj)==999) cycle
    if (AgePriorM(ABS(AgeDiff(A, Bj))+1,k)==0.0) then   
        LL=777
    endif 
enddo
if (LL/=999) return

 call GetAncest(-SB, k, AncB)
if (ANY(AncB == A)) then  ! A>0
    LL = 777
else if (BY(A)/=-999) then
    do x=3, mxA
        do m=1,2
            if (AncB(m, x) > 0) then
                if (AgeDiff(A, AncB(m, x)) <=0) then  ! A older than putative ancestor
                    LL = 777
                endif
            endif
        enddo
    enddo
endif
if (LL == 777) return

Inbr = 0
if (Parent(A,3-k) < 0) then
    if (Parent(A,3-k) == GpID(3-k, SB, k)) then
        Inbr = 1  ! inbreeding loop created
    endif
endif
do f=1, nS(SB,k)
    Bj = SibID(f, SB, k)
    if (Parent(A,3-k) == Bj)  Inbr = 1
    if (Parent(Bj,3-k) == A)  Inbr = 1
enddo

if (Inbr == 0) then
    PrL = 0
    do l=1,nSnp
        call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrY)
        do x=1,3
            if (nS(SB,k)>0) then
                PrX(x) = XPr(3,x,l, SB,k)
            else
                PrX(x) = XPr(2,x,l, SB,k)
            endif
            if (Genos(l,A) /= -9) then
                PrX(x) = PrX(x) * SUM(OKA2P(Genos(l,A), x, :) * PrY)
            endif
        enddo
        PrL(l) = LOG10(SUM(PrX))   
    enddo
    LL = SUM(PrL)
else
    call DoAdd(A, SB, k)
    call CalcU(A, k, -SB, k, LL)
    call RemoveSib(A, SB, k)
endif

end subroutine AddSib

! ##############################################################################################

subroutine MergeSibs(SA, SB, k, LL)  
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LL
integer :: l,x,y, r,v, Bj, Ai, i, G(2), m, Ei, e,f, &
  AncA(2,mxA), AncB(2,mxA), cat(2), catB(ns(SB,k)), catA(ns(SA,k))
double precision :: PrL(nSnp), PrG(3,2), PrXY(3,3), PrE(3), PrH(3)

LL = 999
G = 0  
do m=1,2  
    if (GpID(m,SA,k) /= 0) then
        if (GpID(m,SB,k) /= 0 .and. GpID(m,SA,k)/=GpID(m,SB,k)) then
            LL = 777  ! incompatible grandparents
        else
            G(m) = GpID(m,SA,k)  ! including if GP(B) is dummy
        endif
    else if (GpID(m,SA,k) == 0) then
        G(m) = GpID(m,SB,k)
    endif
enddo
if (GpID(k, SA,k)==-SB .or. GpID(k, SB, k)==-SA) then
    LL = 777
endif
if (LL==777) return
 call GetAncest(-SA, k, AncA)
 call GetAncest(-SB, k, AncB)
if (ANY(AncA(k, 3:mxA) == -SB)) LL = 777
if (ANY(AncB(k, 3:mxA) == -SA)) LL = 777
if (LL==777) return 

cat = 0
if (G(3-k)/=0) then
    do r = 1, nS(SA, k)
        Ai = SibID(r, SA, k) 
        if (nFS(Ai)==0) cycle
        if (Parent(Ai,3-k)==G(3-k)) then
            cat(1) = Ai
        endif
    enddo
    do r = 1, nS(SB, k)
        Bj = SibID(r, SB, k) 
        if (nFS(Bj)==0) cycle
        if (Parent(Bj,3-k)==G(3-k)) then
            cat(2) = Bj
        endif
    enddo
endif 

! TODO: check if sibship parent is inbred

catB = 0
catA = 0
do r = 1, nS(SB, k)
    Bj = SibID(r, SB, k) 
    if (nFS(Bj)==0 .or. Parent(Bj,3-k)==0) cycle
    do v=1,nS(SA,k)
        Ai = SibID(v, SA, k)   
        if (nFS(Ai)==0) cycle
        if (Parent(Ai,3-k) == Parent(Bj,3-k)) then
            catB(r) = v
            catA(v) = 1
        endif
    enddo
enddo

PrL = 0
do l=1,nSnp
    do m=1,2
        if (m/=k .and. cat(1)/=0 .or. cat(2)/=0) then
            call ParProb(l, G(m), m, -1, 0, PrG(:,m))
        else
            call ParProb(l, G(m), m, 0, 0, PrG(:,m))
        endif
    enddo
    do x=1,3
        do y=1,3
            PrXY(x,y) = SUM(AKA2P(x, y, :) * PrG(y,3-k) * PrG(:,k)) 
        enddo
    enddo
    do x=1,3
        do r=1, nS(SB,k)  ! TODO?: collapse into ID vector AB
            Bj = SibID(r, SB, k) 
            if (NFS(Bj) == 0) cycle
            if (cat(2)==Bj) then
                PrE = 1
            else
                call ParProb(l, Parent(Bj,3-k), 3-k, -1, 0, PrE)
            endif
            if (Parent(Bj,3-k) < 0) then 
                do e=1,3
                    do f=1, nS(-Parent(Bj,3-k), 3-k)
                        Ei = SibID(f, -Parent(Bj,3-k), 3-k)
                        if (nFS(Ei) == 0) cycle
                        if (Parent(Ei, k)==-SB .or. Parent(Ei,k)==-SA) cycle  
                        call ParProb(l, Parent(Ei, k), k, Ei, -1, PrH) 
                        do i=1, nFS(Ei)
                            if (Genos(l,FSID(i, Ei))==-9) cycle
                            PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                        enddo
                        PrE(e) = PrE(e) * SUM(PrH)
                    enddo
                enddo
                if (SUM(PrE)<3) then
                    PrE = PrE / SUM(PrE)
                endif
            endif
            do f=1, nFS(Bj)
                if (Genos(l,FSID(f, Bj))==-9) cycle
                PrE = PrE * OKA2P(Genos(l,FSID(f, Bj)), x, :)
            enddo
            if (catB(r)/=0) then
                Ai = SibID(catB(r), SA, k)
                do f=1, nFS(Ai) 
                    if (Genos(l,FSID(f, Ai))==-9) cycle
                    PrE = PrE * OKA2P(Genos(l,FSID(f, Ai)), x, :)
                enddo
            endif
            if (cat(2)==Bj) then
                PrXY(x,:) = PrXY(x,:) * PrE
            else
                PrXY(x,:) = PrXY(x,:) * SUM(PrE)
            endif
        enddo
       
        do r=1, nS(SA,k)
            Ai = SibID(r, SA, k) 
            if (NFS(Ai) == 0) cycle
            if (catA(r)==1) cycle
            if (cat(1)==Ai) then
                PrE = 1
            else
                call ParProb(l, Parent(Ai,3-k), 3-k, -1, 0, PrE)
            endif
            if (Parent(Ai,3-k) < 0) then 
                do e=1,3
                    do f=1, nS(-Parent(Ai,3-k), 3-k)
                        Ei = SibID(f, -Parent(Ai,3-k), 3-k)
                        if (nFS(Ei) == 0) cycle
                        if (Parent(Ei, k)==-SB .or. Parent(Ei,k)==-SA) cycle  
                        call ParProb(l, Parent(Ei, k), k, Ei, -1, PrH) 
                        do i=1, nFS(Ei)
                            if (Genos(l,FSID(i, Ei))==-9) cycle
                            PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                        enddo
                        PrE(e) = PrE(e) * SUM(PrH)
                    enddo
                enddo
                if (SUM(PrE)<3) then
                    PrE = PrE / SUM(PrE)
                endif
            endif
            do f=1, nFS(Ai)
                if (Genos(l,FSID(f, Ai))==-9) cycle
                PrE = PrE * OKA2P(Genos(l,FSID(f, Ai)), x, :)
            enddo
            if (cat(1)==Ai) then
                PrXY(x,:) = PrXY(x,:) * PrE
            else
                PrXY(x,:) = PrXY(x,:) * SUM(PrE)
            endif
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY)) 
enddo
LL = SUM(PrL)

end subroutine MergeSibs

! ##############################################################################################

subroutine AddFS(A, SB, kB, SA, kA, LL)  ! A/SA FS with any B?
use Global
implicit none

integer, intent(IN) :: A, SB, kB, SA, kA
double precision, intent(OUT) :: LL
integer :: l, x, y, Par(nS(SB,kB)), i, Bj, Ei, f, g, MaybeFS(nS(SB,kB)), z, &
 PA, AncA(2,mxA), AncB(2,mxA), AncPF(2,mxA), m, h
double precision :: PrL(nSnp, nS(SB,kB),2), PrY(3,2), PrX(3,2), PrZ(3), &
  dLL(nS(SB,kB)), LLtmp(2), LLUX

PrL = 0
LL = 999

Par = 0  ! shared parent 3-kB  (cand. parent(kB) == SB)
MaybeFS = 1

if (nS(SB,kB)==0) then
    LL = 777
    return   ! nobody to be FS with
endif
call GetAncest(-SB, kB , AncB)
PA = 0
if (A /= 0) then
    PA = Parent(A, 3-kB)
    call GetAncest(A, kA, AncA)
else if (SA /= 0) then   ! TODO: does it matter if kA=kB?
    PA = GpID(3-kB, SA, kA)
     call GetAncest(-SA, kA, AncA)
endif

if (A/=0) then
    if (Parent(A,kB)/=0 .and. Parent(A,kB)/=-SB) then
        LL = 777
    else if (ANY(AncB == A)) then 
        LL = 777
    else
        do x=2, mxA
            do m=1,2
                if (AncB(m, x) > 0) then
                    if (AgeDiff(A, AncB(m, x)) <=0) then  ! A older than putative ancestor
                        LL = 777
                    endif
                endif
            enddo
        enddo
    endif
else if (SA/=0) then
if (GpID(kB, SA, kA)/=0 .and. GpID(kB, SA, kA)/=-SB) then
        LL = 777
    else if (ANY(AncB(kA, 2:mxA) == -SA)) then
        LL = 777
    else
        do x=2, mxA
            do m=1,2
                if (AncB(m, x) > 0) then
                    if (ANY(AgeDiff(SibID(1:nS(SA,kA),SA,kA), AncB(m, x)) <=0)) then  ! A older than putative ancestor
                        LL = 777
                    endif
                endif
            enddo
        enddo
    endif
endif
if (LL /= 999) return

if (ANY(AncA(kB, 3:mxA) == -SB)) then
    LL = 444   ! TODO: check
    return
endif

do f=1, nS(SB,kB)
    if (NFS(SibID(f, SB, kB))==0) then
        MaybeFS(f) = -1
        cycle
    endif   
    do i=1,nFS(SibID(f, SB, kB))
        Bj = FSID(i, SibID(f, SB, kB))
        if (A == Bj) then
            LL = 888
        else if (A >0) then
            if (Parent(A,3-kB) == Bj) then  
               MaybeFS(f) = 0          ! possible, but unlikely    
            else if (Parent(Bj, 3-kB) == A) then
                MaybeFS(f) = 0
            else if (AgeDiff(A,Bj)/=999) then
                if (AgePriorM(ABS(AgeDiff(A, Bj))+1, kB)==0.0) then   
                    LL=777
                else if (AgePriorM(ABS(AgeDiff(A, Bj))+1,3-kB)==0.0) then   
                    MaybeFS(f) = 0
                endif
            endif
        else if (SA/=0 .and. kA/=kB) then
            if (Parent(Bj, 3-kB) == -SA) then
                MaybeFS(f) = 0  ! cannot be FS with own parent
                LL = 444  ! TODO: implement. 
                cycle
            endif
        endif
        if (Bj == PA .or. (A/=0 .and. A == Parent(Bj, 3-kB))) then
            MaybeFS(f) = 0
            cycle
        endif
        if (PA>0) then
            if (Parent(PA,1)==Bj .or. Parent(PA,2)==Bj) then
                MaybeFS(f) = 0
                cycle
            endif
        endif
        
        Par(f) = Parent(Bj, 3-kB)
        if (PA/=0 .and. PA/=Par(f) .and. Par(f)/=0) then
            MaybeFS(f) = 0
        else if (Par(f)==0) then
            Par(f) = PA
        endif
    enddo
    if (Par(f)<0 .and. SA/=0 .and. kA==kB) then
        do i=1, nS(SA,kA)
            if (Parent(SibID(i,SA,kA), 3-kA) == Par(f)) then
                LL = 444   ! TODO: implement (P-O inbreeding)
            endif
        enddo
    endif
enddo
if (LL /= 999) return
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
    LL = 777
    return
endif

do f=1, nS(SB,kB)
    if (nFS(SibID(f, SB, kB))==0) cycle
    if (MaybeFS(f)/=1 .or. Par(f)==0 .or. Par(f)==PA) cycle
    call getAncest(Par(f), 3-kB, AncPF)
    if (Par(f)>0) then
        if (A/=0 .and. ANY(AncPF(:,2:mxA)==A)) MaybeFS(f) = 0
        if (SA/=0 .and. ANY(AncPF(kA,2:mxA)==-SA)) MaybeFS(f) = 0
    else if (Par(f) < 0) then
        if (A/=0 .and. ANY(AncPF(:,3:mxA)==A)) MaybeFS(f) = 0
        if (SA/=0 .and. ANY(AncPF(kA,3:mxA)==-SA)) MaybeFS(f) = 0
        if (SA/=0 .and. GpID(kA,SB,kB)==0 .and. GpID(kA,-Par(f),3-kB)==0) then
            call PairUA(SibID(f,SB,kB),-SA,kB, kA, LLtmp(1))
            if (LLtmp(1) < 0) then ! can't tel if FS or double GP
                MaybeFS(f) = 0
            endif 
        endif
    endif
enddo

if (A/=0 .and. nAgeClasses>1) then  ! A may be not-yet-assigned Parent of Sib f
    do f=1, nS(SB,kB)
        if (MaybeFS(f)/=1 .or. Par(f)/=0 .or. Parent(SibID(f, SB, kB), 3-kB)/=0) cycle  
        if (AgeDiff(SibID(f, SB, kB), A) <= 0) cycle
        call CalcU(-SB, kB, A, kB, LLtmp(1))
        Parent(SibID(f, SB, kB), 3-kB) = A
        call CalcU(-SB, kB, A, kB, LLtmp(2))
        Parent(SibID(f, SB, kB), 3-kB) = 0
        call CalcCLL(SB, kB)
        if (LLtmp(1) - LLtmp(2) < thLRrel) then
            MaybeFS(f) = 0
        endif
    enddo
endif
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
    LL = 777
    return
endif
dLL = 999
do l=1,nSnp
    do f=1, nS(SB,kB)
        if (MaybeFS(f) /= 1) cycle
        do x=1,3
            PrX(x,:) = XPr(2,x,l, SB, kB)
            do g=1,nS(SB,kB)
                Bj = SibID(g, SB, kB)
                if (NFS(Bj) == 0) cycle
                if (g==f) then
                    call ParProb(l, Par(g), 3-kB, -1,0, PrY(:,1))
                else
                    call ParProb(l, Parent(Bj, 3-kB), 3-kB, Bj,-1, PrY(:,1))
                endif
                PrY(:,2) = PrY(:,1)  ! 1: FS, 2: not FS (incl. HS via par(f))
                do y=1,3 
                    do i=1,nFS(Bj)
                        if (Genos(l,FSID(i, Bj))==-9) cycle
                        PrY(y,:) = PrY(y,:) * OKA2P(Genos(l,FSID(i,Bj)), x, y)
                    enddo
                    if (g==f) then
                        if (Par(g) < 0) then 
                            do h = 1, nS(-Par(g), 3-kB)
                                Ei = SibID(h, -Par(g), 3-kB)                               
                                if (Parent(Ei, kB) == -SB) cycle  
                                if (NFS(Ei) == 0) cycle  
                                call ParProb(l, Parent(Ei,kB), kB, Ei,-1, PrZ)
                                do z=1,3                 
                                    do i=1, nFS(Ei)
                                        if (FSID(i,Ei) == A) cycle
                                        if (Genos(l, FSID(i,Ei))/=-9) then
                                            PrZ(z) = PrZ(z) * OKA2P(Genos(l,FSID(i,Ei)),y,z)
                                        endif
                                    enddo
                                enddo
                                PrY(y,:) = PrY(y,:) * SUM(PrZ)
                            enddo
                        endif
                        if (A/=0) then
                            if (Genos(l, A)/=-9) then
                                PrY(y,1) = PrY(y,1) * OKA2P(Genos(l,A), x, y)
                                if (PA/=0) then
                                    PrY(y,2) = PrY(y,2) * OKAP(Genos(l,A), y, l)
                                else
                                    PrY(y,2) = PrY(y,2) * OHWE(Genos(l,A), l)
                                endif
                            endif
                        else if (SA/=0) then
                            PrY(y,1) = PrY(y,1) * SUM(XPr(1,:,l, SA,kA) * AKA2P(:, x, y))
                            if (PA/=0) then
                                PrY(y,2) = PrY(y,2) * SUM(XPr(1,:,l, SA,kA) * AKAP(:, y, l))
                            else
                                PrY(y,2) = PrY(y,2) * SUM(XPr(1,:,l, SA,kA) * AHWE(:,l))
                            endif
                        endif 
                    endif
                enddo  ! y
                do i=1,2
                    PrX(x,i) = PrX(x,i) * SUM(PrY(:,i))
                enddo
            enddo  ! g
        enddo  ! x
        PrL(l,f,:) = LOG10(SUM(PrX, DIM=1))
    enddo  ! f   
enddo

dLL = 777
do f = 1, nS(SB, kB)
    if (MaybeFS(f)/=1) cycle
    dLL(f) = SUM(PrL(:, f,1)) - SUM(PrL(:,f,2))
enddo

if (A/=0) then
    if (nS(SB,kB)==1 .and. ALL(GpID(:,SB,kB)==0)) then
        call CalcU(A,kA, SibID(1,SB,kB), kB, LLUX)
    else
        call CalcU(A,kA, -SB, kB, LLUX)
    endif
    LL = MAXVAL(dLL, MASK=dLL/=777) + LLUX   
else if (SA/=0) then
    call CalcU(-SA, kA, -SB, kB, LLUX)
    do f = 1, nS(SB, kB)
        if (MaybeFS(f)/=1) cycle
        if (Par(f)==0) then
            dLL(f) = SUM(PrL(:, f, 1)) - SUM(PrL(:, f, 2)) + LLUX
        else  ! consider changes in SA (e.g. inbreeding loops) 
            GpID(3-kB, SA, kA) = Par(f)  ! temporary assigned shared parent w. BB(f)
            call PairUA(-SA, -SB, kA, kB, dLL(f))  ! additionally parent via SB
        endif
    enddo
    GpID(3-kB, SA, kA) = PA
    LL = MaxLL(dLL) 
endif
    
end subroutine AddFS

! ##############################################################################################

subroutine AddParent(A, SB, k, LL)  ! is A parent of sibship SB?
use Global
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l,x,y,m, G(2)
double precision :: PrL(nSnp), PrX(3), PrXY(3,3), PrG(3, 2)

LL = 999
G = 0
if (ANY(AgeDiff(SibID(1:nS(SB,k),SB,k), A)<=0)) then  ! A too young to be P or GP of i
    LL = 777
    return
endif

do m=1,2
    if (Parent(A,m)/= 0) then   ! todo: allow for sibship/real parent
        if (GpID(m,SB,k)/= 0 .and. GpID(m,SB,k) /= Parent(A,m)) then
            LL = 777
            return
        else
            G(m) = Parent(A,m)
        endif
    else if(GpID(m,SB,k)/=0) then
        G(m) = GpID(m,SB,k)
    endif
enddo

PrL = 0
do l=1,nSnp
    if (Genos(l,A)==-9) then
        PrL(l) = LOG10(SUM(XPr(3,:,l, SB,k)))
    else
        call ParProb(l, A, k, 0,0, PrX)
        do m=1,2
            call ParProb(l, G(m), m, 0,0, PrG(:,m))    
        enddo
        do x=1,3
            do y=1,3
                PrXY(x,y) = XPr(1,x,l, SB,k) * PrX(x) * &
                  SUM(AKA2P(x, y, :) *  PrG(y,3-k) * PrG(:, k))
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXY))
    endif
enddo

LL = SUM(PrL)

end subroutine AddParent

! ##############################################################################################

subroutine AddGP(A, SB, k, LL)  ! add A as a grandparent to sibship SB
! TODO: check if sharing a parent 3-k
use Global
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l, x,y, m, i, cat, curGP
double precision :: PrL(nSnp), PrY(3), PrXY(3,3), LLtmp(3), PrA(3)

LL = 999
if (Sex(A)/=3) then
    m = Sex(A)
!    if (GpID(m,SB,k) /= 0) then  ! allow replacement  (else change CalcGPZ)
!        LL = 777
!    endif
else if (GpID(1,SB,k)==0) then
    m = 1
else if (GpID(2,SB,k)==0) then
    m = 2
else
    LL = 777
endif
if (LL==777) return

 cat = 0
if (GpID(3-k,SB,k) < 0) then
    if (Parent(A, 3-k)==GpID(3-k,SB,k)) then
        cat = 1
    else
        do i=1,nS(SB,k)
            if (Parent(SibID(i, SB, k), 3-k) == GpID(3-k,SB,k)) then 
                cat = 1
                exit
            endif
        enddo
    endif
endif

do i=1, nS(SB, k)
    if (AgeDiff(A, SibID(i,SB,k))==999) cycle
    if (AgeDiff(SibID(i,SB,k), A)<=0) then  ! A too young to be P or GP of i
        LL = 777
    else if (k==1 .and. m==1 .and. AgePriorM(AgeDiff(SibID(i,SB,k), A)+1, 3)==0.0) then
        LL = 777
    else if (k==2 .and. m==2 .and. AgePriorM(AgeDiff(SibID(i,SB,k), A)+1, 4)==0.0) then
        LL = 777 
    else if (k/=m .and. AgePriorM(AgeDiff(SibID(i,SB,k), A)+1, 5)==0.0) then
        LL = 777         
    endif
enddo
if (LL==777) return

if (cat==0) then
    PrL = 0
    do l=1,nSnp
        call ParProb(l, A, 0, 0, 0, PrA)
        call ParProb(l, GpID(3-m, SB, k), 3-m, 0, 0, PrY)
        do x=1,3
            do y=1,3
                PrXY(x,y) = XPr(1,x,l, SB,k) * SUM(AKA2P(x, :, y) * PrA * PrY(y))
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXY))
    enddo
    LL = SUM(PrL) + Lind(A)
else  ! inbreeding loop present / will be created
    if (GpID(3-m, SB, k) < 0) then
        call CalcU(-SB, k, GpID(3-m, SB, k), 3-m, LLtmp(1))
    endif
    curGP = GPID(m, SB, k)
    GpID(m, SB, k) = A
    call CalcU(-SB, k, A, 3-k, LL)
    if (GpID(3-m, SB, k) < 0) then
        call CalcU(-SB, k, GpID(3-m, SB, k), 3-m, LLtmp(2))
        LL = LL + (LLtmp(2) - LLtmp(1))
    endif
    GPID(m,SB,k) = CurGP
    call CalcCLL(SB, k)
endif
      
end subroutine AddGP

! ##############################################################################################

subroutine AddGGP(A, SB, k, LL)
use Global
implicit none
! A a great-grandparent of sibship SB? (only calculating over non-grandparent-assigned)

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l, x, m, y,z, AncG(2,mxA), i, v, catG, n,GG
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrZ(3), PrA(3), PrP(3), PrV(3)

LL = 999
AncG = 0
if (GpID(1, SB,k)/=0) then
    if (GpID(2, SB,k)/=0) then  ! should be assigned as parent-of-gp
        LL = 777   !(or 888)
    else
        m = 2
    endif
else
    m = 1  ! doesn't really matter (?); GpID(m, SB, k) == 0.
endif
if (LL==777) return

if (Sex(A)<3) then
    n = Sex(A)
else
    n = 1
endif

catG =0
GG = GpID(3-m, SB, k)
if (GG/=0) then
    if (GG==Parent(A,3-m)) then
        catG = 1
    endif
    call GetAncest(GG, 3-m, AncG)   
    if (ANY(AncG == A)) then
        if ((GG>0 .and. AncG(n,2)==A) .or. (GG<0 .and. AncG(n,3-m+2)==A)) then  
            catG = 2  ! already GGP via 3-m; check if double ggp
        else
            LL = 444  ! possible; not yet implemented
        endif
    else if ((Parent(A,1)/=0 .and. ANY(AncG(1, 2:4) == Parent(A,1))) .or. &
      (Parent(A,2)/=0 .and. ANY(AncG(2, 2:4) == Parent(A,2)))) then
        LL = 444   ! TODO: stricter implementation?
    endif
endif
if (GpID(3-k,SB,k) < 0) then
    do i=1,nS(SB,k)
        if (Parent(SibID(i, SB, k), 3-k) == GpID(3-k,SB,k)) then 
            LL = 444
            exit
        endif
    enddo
endif
if (LL==444) return ! TODO: implement.

if (catG/=0) then  ! age check
    if (GG>0) then
        if (AgeDiff(GG, A) >= 0 .and. catG==1 .and. AgeDiff(GG, A)/=999)  LL = 777  ! A older than GG
        if (AgeDiff(GG, A) <= 0 .and. catG==2)  LL = 777  ! GG older than A
    else if (GG<0) then
        do v=1, nS(-GG, 3-m)  ! TODO? AgePriorM; ancestors
            if (AgeDiff(SibID(v,-GG,3-m), A) >= 0 .and. catG==1 .and. &
              AgeDiff(SibID(v,-GG,3-m), A)/=999)  LL = 777
            if (AgeDiff(SibID(v,-GG,3-m), A) <= 0 .and. catG==2)  LL = 777
        enddo
    endif
endif
if (LL == 777) return

PrL = 0
do l=1,nSnp
    call ParProb(l, A, 0, 0, 0, PrA)
    if (catG==1) then
        call ParProb(l, GG, 3-m, A, 0, PrZ)
        call ParProb(l, Parent(A,m), m, -1, 0, PrP)
    else if (catG==2) then
        call ParProb(l, GG, 3-m, -4, 0, PrZ)
        if (GG > 0) then  ! TODO: use AncG
            call ParProb(l, Parent(GG,3-n), 3-n, GG, 0, PrP) 
        else if (GG < 0) then
            call ParProb(l, GpID(3-n, -GG,3-m), 3-n, 0, 0, PrP)
        else
            PrP = AHWE(:,l)
        endif
    else
        call ParProb(l, GG, 3-m, 0, 0, PrZ)
    endif
    do x=1,3  ! sibship parent
        do y=1,3
            do z=1,3
                PrXYZ(x,y,z) = XPr(1,x,l, SB,k) * AKA2P(x, y, z) * PrZ(z)  !note: PrZ catG specific
                if (catG==1) then
                    do v=1,3
                        PrV(v) = AKAP(y, v, l) * PrA(v) * SUM(AKA2P(v,z,:) * PrP)
                    enddo
                    PrXYZ(x,y,z) = PrXYZ(x,y,z) * SUM(PrV) 
                else if (catG==2) then
                    do v=1,3
                        PrV(v) = AKAP(y, v, l) * PrA(v) * SUM(AKA2P(z,v,:) * PrP)
                    enddo
                    PrXYZ(x,y,z) = PrXYZ(x,y,z) * SUM(PrV) 
                else
                    PrXYZ(x,y,z) = PrXYZ(x,y,z) * SUM(AKAP(y, :, l) * PrA)
                endif
            enddo
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYZ))           
enddo
if (catG==1) then
    LL = SUM(PrL)
else if (catG==2 .and. GG>0) then
    LL = SUM(PrL) + Lind(A) - Lind(GG)
else
    LL = SUM(PrL) + Lind(A)
endif
     
end subroutine AddGGP

! ##############################################################################################

subroutine ParentHFS(A, SA, kA, SB, kB, hf, LL)  ! parents of SA and SB HS/FS?
use Global
implicit none

integer, intent(IN) :: A, SA, kA, SB, kB, hf
double precision, intent(OUT) :: LL
integer :: m, G(2), l, x, y, u,v, AncA(2,mxA), AncB(2,mxA), i, j,z, r, Ei, GA, GB,e, &
    DoneA(MaxSibSize), Ai, Bj, nA, AA(maxSibSize), catA(maxSibSize), catB(nS(SB,kB)+1), catG, GGP
double precision :: PrG(3,2), PrL(nSnp), PrXV(3,3,3,3,3,2), PrPA(3, 2), LLm(2), &
    PrGA(3), PrGB(3), PrE(3), PrH(3), PrGG(3)

LLm = 999
G = 0  

 if (A/=0) then
    call GetAncest(A, kA, AncA)
    nA = 1
    AA(1) = A
else
    call GetAncest(-SA, kA, AncA)
    nA = nS(SA,kA)
    AA(1:nA) = SibID(1:nA, SA, kA)
endif
 call GetAncest(-SB, kB, AncB)
 
  G = 0 
do m=1,2
    if (m/=hf .and. hf/=3) cycle
    if (AncA(m, kA+2)/=0) then
        if (AncA(m, kA+2) == -SB) then
            LLm(m) = 777
        else if (AncB(m, kB+2)/=0 .and. AncA(m, kA+2)/=AncB(m, kB+2)) then
            LLm(m) = 777
        else if (AncB(m, kB+2)==0) then
            G(m) = AncA(m, kA+2)
        else if (AncB(m, kB+2)==AncA(m, kA+2)) then
            G(m) = AncA(m, kA+2)
            LLm(m) = 888  ! already are sibs
        else
            LLm(m) = 777
        endif
    else
        if (AncB(m,kB+2)/=0 .and. AncB(m,kB+2) == AncA(m, 2)) then
            LLm(m) = 777
        else 
            G(m) = AncB(m, kB+2)
        endif
    endif
    if (hf==3) then  ! FS
        if (ANY(AncA(kB, 3:mxA) == -SB)) then
            LLm = 777
        else if (A>0) then
            if (ANY(AncB == A)) then
                LLm = 777
            endif
        else if (SA/=0) then
            if (ANY(AncB(kA,3:mxA) == -SA)) then
                LLm = 777
            endif
        endif
    endif
    do x=2,mxA
        if (AncB(m,x) > 0) then
            if (A > 0) then
                if (AgeDiff(A, AncB(m,x)) < 0) then
                    LLm(m) = 777  ! A older than putative ancestor
                endif 
            else 
                if (ANY(AgeDiff(SibID(1:nS(SA,kA),SA,kA), AncB(m,x)) < 0)) then
                    LLm(m) = 777  ! Ai older than putative ancestor
                endif
            endif
        endif
        if (AncA(m,x) > 0) then
            if (ANY(AgeDiff(SibID(1:nS(SB,kB),SB,kB), AncA(m,x)) < 0)) then
                LLm(m) = 777 
            endif
        endif
    enddo
enddo

if (hf==3) then
    if (LLm(1)==777 .or. LLm(2)==777) then
        LL = 777
    else if (LLm(1)==888 .and. LLm(2)==888) then
        LL = 888  ! already are FS
    endif
else
    if (LLm(hf)==777) then 
        LL = 777
    else if (LLm(3-hf)==888) then
        LL = 777   ! already HS, would become FS
    else 
        GA = AncA(3-hf, kA+2)
        GB = AncB(3-hf, kB+2)
    endif
endif

if (ANY(AncA(kB, 5:mxA)==-SB)) then
    LL = 444  ! highly unlikely (but not strictly impossible: TODO)
else if (AncA(kA,2)/=0 .and. ANY(AncB(kA, 5:mxA) == AncA(kA,2))) then  ! PA
    LL = 777
endif
if (LL /=999) return  
   
 catA = 0  
 catB = 0
do i=1, nA
    if (kA/=kB) then
        if (Parent(AA(i), kB) == AncB(kB, 2) .and. AncB(kB, 2)<0) then
            catA(i) = 1
        endif
    else if (kA == kB .and. Parent(AA(i), 3-kA)<0) then  
        do j=1, nS(SB, kB)
            if (Parent(AA(i), 3-kA) == Parent(SibID(j,SB,kB), 3-kB)) then
                catA(i) = 2
                catB(j) = 2
            endif
        enddo
    endif
    if (Parent(AA(i), 3-kA) < 0) then
        if (G(3-kA) == Parent(AA(i), 3-kA)) then  ! incl. hf==3
            if (kA==kB) then
                catA(i) = 3  ! (u) 3-kA = 3-kB == hf 
            else if (kA/=kB) then
                catA(i) = 4  ! (z)
            endif 
        else if (kA==hf .and. GA == Parent(AA(i), 3-kA)) then
            catA(i) = 4  ! (z)
        else if (kA==hf .and. GB == Parent(AA(i), 3-kA)) then
            catA(i) = 5  ! (v)
        endif
    endif
enddo    

do i=1, nS(SB, kB)
    if (kA/=kB) then
        if (Parent(SibID(i,SB,kB), kA) == AncA(kA, 2) .and. AncA(kA,2)<0) then
            catB(i) = 1
        endif
    endif
    if (Parent(SibID(i,SB,kB), 3-kB) < 0) then
        if (G(3-kB) == Parent(SibID(i,SB,kB), 3-kB)) then
            catB(i) = 3  ! (u)  (for hf<3 .and. hf==3)
        else if (kB==hf .and. GA == Parent(SibID(i,SB,kB), 3-kB)) then
            catB(i) = 4  ! (z) (GA of type 3-kB if hf==kB) (hf==3: G_kB cannot be opp. parent)
        else if (kB==hf .and. GB == Parent(SibID(i,SB,kB), 3-kB)) then
            catB(i) = 5  ! (v)
        endif
    endif
enddo 

catG = 0
GGP = 0
if (hf/=3) then
    if (G(hf) > 0) then
        if (Parent(G(hf), 3-hf) == GA .and. GA/=0) then
            catG = 1
            GGP = Parent(G(hf), hf)
        else if (Parent(G(hf), 3-hf) == GB .and. GB/=0) then
            catG = 2
            GGP = Parent(G(hf), hf)
        endif
    else if (G(hf) < 0) then
        if (GpID(3-hf, -G(hf), hf) == GA .and. GA/=0) then
            catG = 1
            GGP = GpID(hf, -G(hf), hf)
        else if (GpID(3-hf, -G(hf), hf) == GB .and. GB/=0) then
            catG = 2
            GGP = GpID(hf, -G(hf), hf)
        endif
    endif
endif

PrL = 0
do l=1,nSnp
    do m=1,2
        if (m/=hf .and. hf/=3) cycle
        if (ANY(CatA==3) .or. ANY(CatB==3)) then
            call ParProb(l, G(m), m, -1,0, PrG(:,m)) 
        else if (catG/=0) then
            call ParProb(l, G(m), m, -4, 0, PrG(:,m))
            if (G(m) > 0) then
                call ParProb(l, GGP, hf, G(m), 0, PrGG) 
            else
                call ParProb(l, GGP, hf, 0, 0, PrGG)
            endif
        else
            call ParProb(l, G(m), m, 0,0, PrG(:,m)) 
        endif
    enddo
    if (hf < 3) then
        if (ANY(CatA==4) .or. ANY(CatB==4)) then
            call ParProb(l, GA, 3-hf, -1,0, PrGA)
        else
            call ParProb(l, GA, 3-hf, 0,0, PrGA)
        endif
        if (ANY(CatA==5) .or. ANY(CatB==5)) then
            call ParProb(l, GB, 3-hf, -1,0, PrGB)
        else
            call ParProb(l, GB, 3-hf, 0,0, PrGB)
        endif
    endif
    if (A>0) then
        if (Genos(l,A)==-9) then
            PrL(l) = LOG10(SUM(XPr(3,:,l, SB,kB)))
            cycle
        else  ! TODO: PrPA for Parent(A,kA) if /=0?
            do m=1,2
                call ParProb(l, Parent(A,m), m, A,0, PrPA(:,m))
            enddo
        endif
    endif
    
    PrXV = 0
    do x=1,3  ! SA/PA
        do y=1,3  ! SB
            do u=1,3  ! G_hf / G_3-kB (hf==3)
                do z=1,3  ! G_A (hf/=3) / G_kB (hf==3)
                    do v=1,3 ! G_B (hf/=3)
                        if (hf==3) then  ! 0 for z/=v
                            PrXV(x,y,u,z,z,:) = AKA2P(x,u,z) * AKA2P(y,u,z) * PrG(u,3-kB) * PrG(z, kB)
                        else
                            if (GA < 0 .and. GA == -SB) then
                                PrXV(x,y,u,y,v,:) = AKA2P(x,u,y) * AKA2P(y,u,v) * PrG(u,hf) * PrGB(v)
                            else if (GB < 0 .and. GB == -SA) then
                                PrXV(x,y,u,z,x,:) = AKA2P(x,u,z) * AKA2P(y,u,x) * PrG(u,hf) * PrGA(z)
                            else if (catG == 1) then
                                PrXV(x,y,u,z,v,:) = AKA2P(x,u,z) * AKA2P(y,u,v) * PrG(u,hf) * &
                                  SUM(AKA2P(u,z,:) * PrGG) * PrGA(z) * PrGB(v)
                            else if (catG == 2) then
                                PrXV(x,y,u,z,v,:) = AKA2P(x,u,z) * AKA2P(y,u,v) * PrG(u,hf) * &
                                  SUM(AKA2P(u,v,:) * PrGG) * PrGA(z) * PrGB(v)
                            else
                                PrXV(x,y,u,z,v,:) = AKA2P(x,u,z) * AKA2P(y,u,v) *PrG(u,hf) * PrGA(z) * PrGB(v)
                            endif
                        endif
                        if (A /=0) then
                            if (Parent(A, kA)/=0) then
                                PrXV(x,y,u,z,v,:) = PrXV(x,y,u,z,v,:) * PrPA(x, kA)
                            endif
                        endif
                    enddo
                enddo
            enddo
            
            DoneA = 0            
            if (ALL(catA==0) .and. ALL(catB==0)) then
                if (SA/=0) then
                    PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * XPr(1,x,l, SA,kA) * XPr(1,y,l, SB,kB)
                else if (A>0) then
                    PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * SUM(OKA2P(Genos(l,A), x, :) * PrPA(:,3-kA)) * &
                      XPr(1,y,l, SB,kB)
                endif            
            else
                do r=1, nS(SB,kB)
                    Bj = SibID(r, SB, kB) 
                    if (NFS(Bj) == 0) cycle 
                    if (catB(r)==0 .or. catB(r)==2) then
                        call ParProb(l, Parent(Bj,3-kB), 3-kB, -1, 0, PrE)
                    else
                        PrE = 1
                    endif                                           
                    if (Parent(Bj,3-kB) <0 .and. CatB(r)/=1) then
                        do e=1,3
                            do v=1, nS(-Parent(Bj,3-kB), 3-kB)
                                Ei = SibID(v, -Parent(Bj,3-kB), 3-kB)
                                if (nFS(Ei) == 0) cycle
                                if (Parent(Ei, kB) == -SB) cycle
                                if (Parent(Ei, kA) == Parent(AA(1),kA) .and. Parent(AA(1),kA)/=0) cycle
                                if (Ei==Parent(AA(1),kA)) cycle
                                call ParProb(l, Parent(Ei, kB), kB, Ei, -1, PrH) ! TODO: Parent(Ei, kB)==GPs
                                do i=1, nFS(Ei)
                                    if (Genos(l,FSID(i, Ei))==-9) cycle
                                    PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                                enddo
                                PrE(e) = PrE(e) * SUM(PrH)
                            enddo
                        enddo
                    endif
                    
                if (catB(r)==0 .or. catB(r)==2) then 
                        PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * SUM(PrE)
                    else if (catB(r)==1) then  ! Parent(Bj,3-kB) = PA, kA/=kB
                        PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * PrE(x)
                    else if (catB(r)==3) then  ! hf==kB, Parent(Bj,3-kB) = GA
                        do u=1,3
                            PrXV(x,y,u,:,:,1) = PrXV(x,y,u,:,:,1) * PrE(u)
                        enddo
                    else if (catB(r)==4) then  ! Parent(Bj,3-kB) = G
                        do z=1,3
                            PrXV(x,y,:,z,:,1) = PrXV(x,y,:,z,:,1) * PrE(z)
                        enddo
                    else if (catB(r)==5) then  ! hf==kB, Parent(Bj,3-kB) = GB
                        do v=1,3
                            PrXV(x,y,:,:,v,1) = PrXV(x,y,:,:,v,1) * PrE(v)
                        enddo
                    endif
         
                    do j=1, nFS(Bj) 
                        if (Genos(l,FSID(j, Bj))==-9) cycle
                        PrE =  PrE * OKA2P(Genos(l,FSID(j,Bj)), y, :)
                    enddo

                    if (catB(r)==2) then  ! kA==kB, share parent 3-kB
                        do v = 1, nA
                            Ai = AA(v)
                            if (SA/=0 .and. nFS(Ai) == 0) cycle
                            if (Parent(Ai, 3-kA)/=Parent(Bj,3-kB)) cycle
                            do i=1, nFS(Ai)  
                                if (A/=0 .and. FSID(i, Ai)/=A) cycle
                                if (Genos(l,FSID(i, Ai))==-9) cycle
                                PrE =  PrE * OKA2P(Genos(l,FSID(i,Ai)), x, :)
                            enddo
                            doneA(v) = 1
                        enddo
                    endif
                    
                    if (catB(r)==0 .or. catB(r)==2) then 
                        PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * SUM(PrE)
                    else if (catB(r)==1) then  ! Parent(Bj,3-kB) = PA, kA/=kB
                        PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * PrE(x)
                    else if (catB(r)==3) then  ! hf==kB, Parent(Bj,3-kB) = GA
                        do u=1,3
                            PrXV(x,y,u,:,:,2) = PrXV(x,y,u,:,:,2) * PrE(u)
                        enddo
                    else if (catB(r)==4) then  ! Parent(Bj,3-kB) = G
                        do z=1,3
                            PrXV(x,y,:,z,:,2) = PrXV(x,y,:,z,:,2) * PrE(z)
                        enddo
                    else if (catB(r)==5) then  ! hf==kB, Parent(Bj,3-kB) = GB
                        do v=1,3
                            PrXV(x,y,:,:,v,2) = PrXV(x,y,:,:,v,2) * PrE(v)
                        enddo
                    endif
                enddo
                
                do r = 1, nA
                    if (doneA(r)==1) cycle
                    if (SA/=0 .and. NFS(AA(r)) == 0) cycle
                    if (kA/=kB .and. Parent(AA(r),3-kA)==-SB) cycle  ! done
                    if (catA(r)==0) then
                        call ParProb(l, Parent(AA(r),3-kA), 3-kA, -1, 0, PrE)
                    else
                        PrE = 1
                    endif
                    if (Parent(AA(r), 3-kA) < 0 .and. (SA/=0 .or. ANY(FSID(1:nFS(AA(r)), AA(r))==A))) then
                        do e=1,3
                            do i=1, nS(-Parent(AA(r), 3-kA), 3-kA)
                                Ei = SibID(i, -Parent(AA(r), 3-kA), 3-kA)
                                if (nFS(Ei) == 0) cycle
                                if (Parent(Ei, kB) == -SB) cycle
                                if (SA/=0 .and. Parent(Ei, kA) == -SA) cycle
                                call ParProb(l, Parent(Ei, kA), kA, Ei, -1, PrH)  ! Assume is not one of GPs
                                do j=1, nFS(Ei)
                                    if (A/=0 .and. FSID(i, Ei)==A) cycle
                                    if (Genos(l,FSID(j, Ei))/=-9) then
                                        PrH = PrH * OKA2P(Genos(l,FSID(j,Ei)), :, e)
                                    endif
                                enddo
                                PrE(e) = PrE(e) * SUM(PrH)
                            enddo
                        enddo
                    endif
                    
                    if (catA(r)<3) then
                        PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * SUM(PrE)
                     else if (catA(r)==3) then 
                        do u=1,3
                            PrXV(x,y,u,:,:,1) = PrXV(x,y,u,:,:,1) * PrE(u)
                        enddo
                    else if (catA(r)==4) then 
                        do z=1,3
                            PrXV(x,y,:,z,:,1) = PrXV(x,y,:,z,:,1) * PrE(z)
                        enddo
                    else if (catA(r)==5) then 
                        do v=1,3
                            PrXV(x,y,:,:,v,1) = PrXV(x,y,:,:,v,1) * PrE(v)
                        enddo
                    endif
                    
                    do i=1, nFS(AA(r))
                        if (Genos(l,FSID(i, AA(r)))==-9) cycle
                        if (SA/=0 .or. FSID(i, AA(r))==A) then 
                            PrE =  PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)
                            doneA(r) = 2
                        else
                            PrE =  PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)
                        endif
                    enddo
                    
                    if (catA(r)<3) then
                        PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * SUM(PrE)
                     else if (catA(r)==3) then 
                        do u=1,3
                            PrXV(x,y,u,:,:,2) = PrXV(x,y,u,:,:,2) * PrE(u)
                        enddo
                    else if (catA(r)==4) then 
                        do z=1,3
                            PrXV(x,y,:,z,:,2) = PrXV(x,y,:,z,:,2) * PrE(z)
                        enddo
                    else if (catA(r)==5) then 
                        do v=1,3
                            PrXV(x,y,:,:,v,2) = PrXV(x,y,:,:,v,2) * PrE(v)
                        enddo
                    endif
                enddo
            endif
        enddo
    enddo
    PrL(l) = LOG10(SUM(PrXV(:,:,:,:,:,2))) - LOG10(SUM(PrXV(:,:,:,:,:,1)))
enddo

LL = SUM(PrL)

end subroutine ParentHFS

! ##############################################################################################

subroutine DummyGP(SA, SB, kA, kB, LL)  ! SB GP of SA? via observed or unobserved B
use Global
implicit none

integer, intent(IN) :: SA, SB, kA, kB
double precision, intent(OUT) :: LL
integer :: i, m, l, x, y, z, G(2), ggp(2), v, w, Bi, r, AncB(2,mxA), catB, catA, Ai
double precision :: LLGX(2), PrL(nSnp), PrZ(3), PrXYZ(3,3,3,3), PrG(3), &
    PrPG(3), PrW(3), LLtmp(maxSibSize, 2)

LL = 999
do i=1, nS(SB,kB)
    if (kA /= kB) then
        if (Parent(SibID(i,SB,kB), kA) == -SA) then
            LL = 444
        endif
    else if (kA == kB) then
        do r= 1, nS(SA, kA)
            if (Parent(SibID(i,SB,kB), 3-kB)==Parent(SibID(r,SA,kA), 3-kA) .and. &
              Parent(SibID(i,SB,kB), 3-kB)/=0) then
                LL = 444
            endif
        enddo
    endif
enddo
if (LL == 444)  return  ! TODO

 call GetAncest(-SB, kB, AncB)
if (ANY(AncB(kA,3:mxA) == -SA)) then
    LL = 777 
    return
else if (ANY(AncB(3-kA, 3:mxA) < 0)) then
    do r= 1, nS(SA, kA)
        if (ANY(AncB(3-kA, 3:mxA)/=0 .and. AncB(3-kA, 3:mxA) == Parent(SibID(r,SA,kA), 3-kA))) then
            LL = 444 ! not implemented
            return
        endif
    enddo
endif
G = GpID(:, SA, kA)

catA = 0
catB = 0
do r = 1, nS(sB,kB)   ! check if overlap
    Bi = SibID(r, sB, kB)
    if (NFS(Bi) == 0) cycle
    if (Parent(Bi, 3-kB) == G(3-kB) .and. G(3-kB)<0) then
        catB = Bi
    endif
enddo
do r = 1, nS(sA,kA)   ! check if inbreeding loop
    Ai = SibID(r,SA,kA)
    if (NFS(Ai)==0) cycle
    if (Parent(Ai, 3-kA) == G(3-kA) .and. G(3-kA)<0) then 
        catA = Ai
    endif
enddo           

LLGX = 999
LLtmp = 999
GGP = 0
do m=1,2
    if (m==kB .and. GpID(kB, SA, kA) == -SB) then
        LLGX(m) = 777
        cycle
    else if (G(m) > 0) then
        if (Parent(G(m), kB) /=0) then
            LLGX(m) = 777
        else 
            call AddSib(G(m), SB, kB, LLtmp(1,m))
            call AddFS(G(m), SB, kB,0,m, LLtmp(2,m))
            if (MaxLL(LLtmp(:,m)) < 0) then
                LLGX(m) = MaxLL(LLtmp(:,m)) + CLL(SA, kA)
                if (Parent(G(m), kB) /= -SB) then
                    LLGX(m) = LLGX(m) - Lind(G(m))
                endif
            else
                LLGX(m) = MaxLL(LLtmp(:,m))
            endif
        endif
        cycle
    else if (G(m) == 0) then
        do r=1, nS(sB,kB)
            Bi = SibID(r, sB, kB) 
            if (Sex(Bi)/=m .and. Sex(Bi)/=3) cycle
            call AddGP(Bi, SA, kA, LLtmp(r,m))
            if (LLtmp(r,m) < 0) then
                LLtmp(r,m) = LLtmp(r,m) - Lind(Bi) + CLL(SB,kB)
            endif
        enddo
        LLGX(m) = MaxLL(LLtmp(:,m))
    else if (G(m) < 0) then 
        if (GpID(kB, -G(m), m) /=0) then
            LL = 777
        else
            GGP(m) = GpID(3-kB, -G(m), m)
        endif
    endif
    
    PrL = 0
    do l=1,nSnp
        if (catB /= 0 .and. m==kB) then
            call ParProb(l, G(3-m), 3-m, catB, -1, PrZ)
        else if (catA/=0) then   ! TODO: catB/=0 .and. catA/=0
            call ParProb(l, G(3-m), 3-m, catA, -1, PrZ)
        else
            call ParProb(l, G(3-m), 3-m, 0, 0, PrZ)
        endif
        if (G(m) /= 0) then
            call ParProb(l, G(m), m, -4, 0, PrG)   !  G(m)'s offspring only, if <0
            PrG = PrG/SUM(PrG)
        else
            PrG = 1
        endif
        call ParProb(l, GGP(m), 3-kB, 0, 0, PrPG)

        PrXYZ = 0
        do x=1,3  ! SA
            do y=1,3  ! parent of SA, offspr of SB
                do v=1,3   ! SB 
                    do z=1,3  ! other parent of SA
                        PrXYZ(x,y,z,v) = AKA2P(x, y, z) * PrG(y) * PrZ(z) * &
                          SUM(AKA2P(y, v, :) * PrPG)
                    enddo
                    if (catA==0) then
                        PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * XPr(1,x,l, SA,kA)   
                    else 
                        do r=1, nS(sA,kA)
                            Ai = SibID(r, sA, kA)  
                            if (NFS(Ai) == 0) cycle
                            if (Bi == G(m)) cycle
                            if (catA==Ai) then
                                PrW = 1
                            else
                                call ParProb(l, Parent(Ai, 3-kA), 3-kA, Ai, -1, PrW)
                            endif
                            do w=1,3
                                do i=1, nFS(Ai)
                                    if (Genos(l,FSID(i, Ai))==-9) cycle
                                    PrW(w) = PrW(w) * OKA2P(Genos(l,FSID(i,Ai)), x, w)
                                enddo
                            enddo
                            if (catA==Ai) then
                                if (m==kA) then
                                    PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * PrW
                                else if (m/=kA) then
                                    PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * PrW(y)
                                endif
                            else
                                PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * SUM(PrW)
                            endif
                        enddo
                    endif
                    
                    if (catB==0) then
                        PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * XPr(3,v,l, SB,kB) 
                    else if (catB/=0) then
                        PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * XPr(2,v,l, SB,kB) 
                        do r=1, nS(sB,kB)
                            Bi = SibID(r, sB, kB)  
                            if (NFS(Bi) == 0) cycle
                            if (catB==Bi) then
                                PrW = 1
                            else
                                call ParProb(l, Parent(Bi, 3-kB), 3-kB, Bi, -1, PrW)
                            endif
                            do w=1,3
                                do i=1, nFS(Bi)
                                    if (Genos(l,FSID(i, Bi))==-9) cycle
                                    PrW(w) = PrW(w) * OKA2P(Genos(l,FSID(i,Bi)), v, w)
                                enddo
                            enddo
                            if (catB==Bi) then
                                PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * PrW
                            else
                                PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * SUM(PrW)
                            endif
                        enddo
                    endif  
                enddo
            enddo
        enddo
        PrL(l) = LOG10(SUM(PrXYZ))         
    enddo
    LLGX(m) = MaxLL((/ SUM(PrL),  LLGX(m)/))
enddo
LL = MaxLL(LLGX)

end subroutine DummyGP

! ##############################################################################################
!       Age priors  
! ##############################################################################################

subroutine CalcDumBY(s, k)   ! age prior for dummy parent, based on offspring & GP BY. 
! updates DumBY(:, s, k) as side effect
use Global
implicit none

integer, intent(IN) :: s, k
integer :: i, y, m, BYS(maxSibSize), BYGP(2)
double precision :: EstBY(nAgeClasses+ maxAgePO), z

BYS(1:nS(s,k)) =  BY(SibID(1:nS(s,k),s,k))
BYGP = -1
do m=1,2
    if (GpID(m,s,k)>0) then
        BYGP(m) = BY(GpID(m,s,k))
    endif  ! TODO: dummy GP
enddo
WHERE(BYS >= 0) BYS = BYS + maxAgePO    ! max age of parent;  shift to allow dummy BY pre-zero
WHERE(BYGP >= 0) BYGP = BYGP + maxAgePO

z = 0.0D0
EstBY = 0.0D0  ! lOG10(1)
do y=1,nAgeClasses + maxAgePO
    do i = 1, ns(s,k)
        if (BYS(i) < 0) cycle  ! unknown BY
        if (BYS(i) - y <= 0 .or. (BYS(i) - y +1) > nAgeClasses) then
            EstBY(y) = LOG10(z)  ! Sib i born prior to year y
            exit
        else
            EstBY(y) = EstBY(y) + LOG10(AgePriorM(BYS(i) - y +1, 6+k))
        endif
    enddo
    do m=1,2
        if (BYGP(m) < 0) cycle  ! no GP / dummy GP / unknown BY
        if (y - BYGP(m) <= 0) then
            EstBY(y) = LOG10(z) 
            exit
        else if (y - BYGP(m) < nAgeClasses) then
            EstBY(y) = EstBY(y) + LOG10(AgePriorM(y - BYGP(m) +1, 6+m))
!            else  ! TODO? or set to 0?  
        endif
    enddo
enddo
EstBY = 10**EstBY
if (SUM(EstBY)/=0.0) EstBY = EstBY / SUM(EstBY)  ! scale to sum to unity 
EstBY = LOG10(EstBY)
DumBY(:, s, k) = EstBY

end subroutine CalcDumBY

! ########################################################################################

subroutine CalcALRmerge(SA, SB, k, ALR)  ! change in ALR when merging
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: ALR
integer :: i, j

ALR = 0.0
do i = 1, nS(SA,k)
    if (BY(SibID(i,SA,k))<0) cycle
    do j=1, nS(SB, k)
        if (AgeDiff( SibID(i,SA,k), SibID(j,SB,k))/=999) then
            ALR = ALR + LOG10(AgePriorM(ABS(AgeDiff( SibID(i,SA,k), SibID(j,SB,k) ))+1, k)) 
        endif
    enddo
enddo

if (ALR < -HUGE(0.0D0)) then
    ALR = 777
endif

end subroutine CalcALRmerge

! ##############################################################################################

subroutine CalcAgeLR(A, kA, B, kB, m, focal, ALR)  ! m: mat/pat relatives
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB, m, focal
double precision, intent(OUT) :: ALR
integer :: x, y, z, gcat
double precision :: BYLR(nAgeClasses+ maxAgePO, 2), tmp, &
    ALRtmp(nAgeClasses+ maxAgePO, nAgeClasses+ maxAgePO)

! allocate(DumBY(nAgeClasses +maxAgePO, nInd/2, 2)) 
tmp = 0.0
BYLR = LOG10(tmp)  ! likelihood ratio to be born in year X 
ALRtmp = LOG10(tmp)
if (A < 0) then
    call CalcDumBY(-A, kA)
    BYLR(:, 1) = DumBY(:, -A, kA)
else if (A>0) then
    if (BY(A) < 0) then  ! unknown BY  - TODO: use parent & offspring age
        ALR = 777
        return   
    else
        BYLR(BY(A) + maxAgePO, 1) = LOG10(1.0)
    endif
endif
if (B < 0) then
    call CalcDumBY(-B, kB)
    BYLR(:, 2) = DumBY(:, -B, kB)
else if (B>0) then
    if (BY(B) < 0) then 
        ALR = 777
        return   
    else
        BYLR(BY(B) + maxAgePO, 2) = LOG10(1.0)
    endif
endif

if (focal==4) then
    if (m==1 .and. kB==1) then
        gcat = 3
    else if (m==2 .and. kB==2) then
        gcat = 4
    else
        gcat = 5
    endif
endif

do y=1,nAgeClasses + maxAgePO -1  ! B
    if (BYLR(y,2) < -HUGE(0.0D0)) cycle
    do x=1, nAgeClasses + maxAgePO  ! A 
        if (BYLR(x,1) < -HUGE(0.0D0)) cycle
        ALRtmp(x,y) = BYLR(x,1) + BYLR(y,2) 
        if (x > y .and. (x-y) < nAgeClasses) then
            if (focal == 1) then  ! B parent of A
                ALRtmp(x,y) = ALRtmp(x,y) + LOG10(AgePriorM(x - y +1, 6+kB))
            else if (focal == 4) then  ! B GP of A
                ALRtmp(x,y) = ALRtmp(x,y) + LOG10(AgePriorM(x - y +1, gcat))
            endif
        endif
        z = ABS(x - y + 1)
        if (z<1 .or. z > nAgeClasses) cycle
        if (focal == 2) then  ! FS
            ALRtmp(x,y) = ALRtmp(x,y) + LOG10(AgePriorM(z, 1)) + LOG10(AgePriorM(z, 2))
        else if (focal == 3) then  !HS
            ALRtmp(x,y) = ALRtmp(x,y) + LOG10(AgePriorM(z, m))
        else if (focal == 5 .or. focal==6) then  !FA / HA
            ALRtmp(x,y) = ALRtmp(x,y) + LOG10(AgePriorM(z, 6))
        endif  
    enddo
enddo
ALR = LOG10(SUM(10**ALRtmp))  ! sum rather than multiply across age differences
 
if (ALR < -HUGE(0.0D0)) then
    ALR = 777
endif
end subroutine CalcAgeLR

! ########################################################################################

! ########################################################################################

subroutine BestRel(LLIN, focal, X, dLL)
use Global
implicit none
! return which relationship is most likely, by threshold thLRrel
! assuming order PO,FS,HS,GG,FAU,HAU,U in LL vector 

double precision, intent(IN) :: LLIN(7)
integer, intent(IN) :: focal
integer, intent(OUT) :: X
double precision, intent(OUT) :: dLL   ! diff best vs next best
double precision :: LL(7)
integer :: i,j, maybe(6), Y(7)

X=0
dLL = 0
LL = LLIN

if (ALL(LL > 0)) then
    X = 8
    return
endif  
if (focal/=2 .and. LL(2) < 0 .and. LL(2)>=LL(3)) then  ! FS are also PHS/MHS, conditional on 1 parent
    LL(3) = 333
else if (focal==3 .and. LL(3)>LL(2) .and. LL(3)<0) then
    LL(2) = 333 ! want sib vs non-sib
endif

if ((LL(7) - MAXVAL(LL(1:6), MASK=LL(1:6)<0)) > thLRrel) then  
    X = 7  ! unrelated
else
    maybe = 1
    do i=1,6
        if (LL(i)>0) then
            maybe(i) = 0
        else
            do j=1,7
                if (i==j) cycle
                if (LL(j)>0) cycle
                if ((LL(i) - LL(j)) < thLRrel) then
                    maybe(i) = 0   ! i has no longer highest LL
                endif
            enddo
        endif
    enddo
    if (SUM(maybe)==0) then
        X = 8  ! unclear
    else if (SUM(maybe)==1) then
        X = MAXLOC(LL(1:6), MASK=maybe==1, DIM=1)
    endif       
endif

Y = (/(i, i=1,7, 1)/)
if (X/=8) then
    dLL = LL(X) - MAXVAL(LL, MASK=(LL<0 .and. Y/=X))
endif  
            
end subroutine BestRel

! ##############################################################################################

subroutine UpdateAllProbs
use Global
implicit none

integer :: i, k, s

do k=1,2
    do s=1,nC(k)
        call CalcCLL(s, k)
    enddo
enddo
do i=1,nInd  ! TODO? BY order
    call CalcLind(i)
enddo

end subroutine UpdateAllProbs

! ##############################################################################################

subroutine CalcLind(i)
use Global
implicit none

integer, intent(IN) :: i
integer :: l, x, y, k, z
double precision :: PrL(nSnp), Px(3,2), PrXYZ(3,3,3), PrX(3)

PrL = 0
do l=1,nSnp
    do k=1,2
        call ParProb(l, Parent(i,k), k, i,0, Px(:,k))
    enddo  ! TODO?: Joined parental probs?
    do x=1,3
        do y=1,3
            do z=1,3
                PrXYZ(x,y,z) = AKA2P(x, y, z) * Px(y,1) * Px(z,2)
            enddo
        enddo
        if (Genos(l,i)/=-9) then
            PrXYZ(x,:,:) = OcA(Genos(l,i), x) * PrXYZ(x,:,:)     
        endif
    enddo
    PrL(l) = LOG10(SUM(PrXYZ))
    do x=1,3
        PrX(x) = SUM(PrXYZ(x,:,:))    
    enddo
    LindG(:, l, i) = PrX / SUM(PrX)
enddo
                
Lind(i) = SUM(PrL)
LindX(:,i) = PrL

end subroutine CalcLind

! ###########################################################################

 subroutine CalcCLL(s, k)  ! TODO: parent(3-k) offspring of s, real or dummy
use Global
implicit none
! returns XPr: likelihood;  DumP: probability, scaled  (no age prior.),
! split into 1: sibs only 2: gp effect only, 3: all

integer, intent(IN) :: s, k  ! S: sibship number, k: maternal(1), paternal(2), or unknown(3)
integer :: l, x, i, Ei, r, y, z, g, Ri, v, cat, catG, Inbr(maxSibSize)
double precision :: PrL(nSnp), PrY(3,2), PrYp(3), PrGG(3,2), PrZ(3), PrXZ(3,3,2)

cat = 0
catG = 0
Inbr = 0
do r=1,nS(s,k)
    Ri = SibID(r, s, k)
    if (Parent(Ri, 3-k)/=0 .and. Parent(Ri, 3-k)==GpID(3-k,s,k)) then  
        cat = Ri
    endif
    if (nFS(Ri)==0) cycle
    do i=1, nFS(Ri)
        do v=1, nS(s,k)
            if (r==v) cycle
            if (Parent(SibID(v,s,k), 3-k) == FSID(i,Ri)) then
                Inbr(r) = -1
                Inbr(v) = r
            endif
        enddo
    enddo
enddo

if (ALL(GpID(:,s,k)<0) .and. cat==0) then  ! check if sibship parent is inbred
    if (GPID(1,s,k) == GPID(1, -GPID(2,s,k),2)) then
        catG = 2
    else if (GPID(2,s,k) == GPID(2, -GPID(1,s,k),1)) then
        catG = 1
    endif
endif

PrL = 0       
do l=1,nSnp
    do g=1,2   !grandparents
        if (g/=k .and. cat>0) then
            call ParProb(l, GpID(g,s,k), g, -1,0, PrGG(:,g))
        else if (catG==g) then
            PrGG(:,g) = XPr(3,:,l, -GpID(g,s,k),g)
        else
            call ParProb(l, GpID(g,s,k), g, 0,0, PrGG(:,g))
        endif
    enddo
    
    do x=1,3  ! genotype dummy parent
        do z=1,3
            if (catG==0) then
            PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(z,3-k) * PrGG(:,k))  ! GPs
            else if (catG==k) then
                PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(:,k))
            else if (catG==3-k) then
                PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(z,3-k))
            endif
        enddo
    enddo
    if (catG>0) then
        do r=1,2
            PrXZ(:,:,r) = PrXZ(:,:,r)/SUM(PrXZ(:,:,r)) 
        enddo
    endif
    
    do x=1,3
        XPr(2,x,l, s,k) = SUM(PrXZ(x,:,2))  ! GP 
        do r=1, nS(s,k)
            Ri = SibID(r, s, k)  ! array with IDs
            if (NFS(Ri) == 0) cycle  ! moved to its FS
            if (Inbr(r) > 0) cycle  ! inbred
            if (cat==Ri) then
                PrYp = 1
            else
                call ParProb(l, Parent(Ri, 3-k), 3-k, -1,0, PrYp) 
            endif 
            
            do y=1,3   ! parent 3-k           
                PrY(y,:) = PrYp(y)  ! dim2: 1:all, 2:non-sibs only
                do i=1, nFS(Ri)  ! default: nFS = 1
                    if (Inbr(r)==0) then
                        if (Genos(l,FSID(i, Ri))/=-9) then
                            PrY(y,1) = PrY(y,1) * OKA2P(Genos(l,FSID(i,Ri)), x, y)
                        endif
                    else
                        do z=1,3
                            PrZ(z) = AKA2P(z, x, y)
                            if (Genos(l,FSID(i, Ri))/=-9) then
                                PrZ(z) = PrZ(z) * OcA(Genos(l,FSID(i,Ri)), z)
                            endif
                            do v=1, nS(s,k)
                                if (Parent(SibID(v,s,k),3-k)==FSID(i,Ri)) then
                                    if (Genos(l,SibID(v,s,k))/=-9) then
                                        PrZ(z) = PrZ(z) * OKA2P(Genos(l,SibID(v,s,k)), x, z)
                                    endif
                                endif
                            enddo
                        enddo
                        PrY(y,1) = PrY(y,1) * SUM(PrZ)
                    endif
                enddo                
                
                if (Parent(Ri, 3-k) < 0) then
                    do v=1, nS(-Parent(Ri, 3-k), 3-k)
                        Ei = SibID(v, -Parent(Ri, 3-k), 3-k)  
                        if (NFS(Ei) == 0) cycle
                        if (Parent(Ei, k) == -s) cycle
                        call ParProb(l, Parent(Ei, k), k, Ei,-1, PrZ)        
                        do i=1, nFS(Ei) 
                            if (Genos(l,FSID(i, Ei))/=-9) then
                                PrZ = PrZ * OKA2P(Genos(l,FSID(i,Ei)), :, y)
                            endif
                        enddo
                        PrY(y,:) = PrY(y,:) * SUM(PrZ)  ! all; non-sibs only
                    enddo
                endif
            enddo 
            
            do z=1,2            
                if (cat==Ri) then
                    PrXZ(x,:,z) = PrXZ(x,:,z) * PrY(:,z)    ! all; non-sibs only
                else
                    PrXZ(x,:,z) = PrXZ(x,:,z) * SUM(PrY(:,z))
                endif
            enddo
        enddo ! r 
    enddo ! x
    do x=1,3  ! account for GP, dumm offspr & connected sibships
        XPr(3,x,l, s,k) = SUM(PrXZ(x,:,1))/ SUM(PrXZ(:,:,2))
        DumP(x,l, s,k) = XPr(3,x,l, s,k)/ SUM(XPr(3,:,l, s,k))
        XPr(1,x,l, s,k) = XPr(3,x,l, s,k) / XPr(2,x,l, s,k) 
    enddo 
    PrL(l) = LOG10(SUM(XPr(3,:,l, s,k))) 
enddo
CLL(s,k) = SUM(PrL)

end subroutine CalcCLL

! ##############################################################################################

subroutine ParProb(l, i, k, A, B, prob)  ! TODO: B<0 : no GP
use Global
implicit none

integer, intent(IN) :: l, i, k, A,B
double precision, intent(OUT) :: prob(3)
integer :: x,j, AB(2)
double precision :: PrP(3, 2), PrY(3)

 if (i == 0) then  ! no parent
    prob = AHWE(:, l)
else if (i > 0) then  ! real parent
    prob = LindG(:, l, i)
else if (i < 0) then  ! dummy parent
    if (A==0) then   ! probability
        prob = DumP(:,l, -i,k)    
    else if (A == -1) then  ! grandparent contribution only
        prob = XPr(2,:,l, -i, k) 
    else if (A==-4) then  ! offspring contribution only
        prob = XPr(1,:,l, -i, k)        
    else if (A>0) then   ! exclude indiv A from calc & standardise
        if ((Genos(l,A)==-9 .and. (nFS(A)<=1 .or. B>=0)) .or. Parent(A,k)/=i) then
            prob = DumP(:,l, -i,k)
        else
            AB = (/ A, B /)
            do j=1,2
                if (j==2 .and. B<=0) cycle
                if (Parent(AB(j), 3-k)==0) then  
                    PrP(:,j) = AHWE(:,l)
                else if (Parent(AB(j), 3-k)>0) then
                    if (Genos(l,Parent(AB(j), 3-k)) /= -9) then
                        PrP(:,j) = OcA(Genos(l, Parent(AB(j),3-k)), :)
                    else
                        PrP(:,j) = LindG(:, l, Parent(AB(j), 3-k))
                    endif
                else if (Parent(AB(j), 3-k)<0) then  
                    PrP(:,j) = DumP(:,l, -Parent(AB(j),3-k), 3-k)  ! TODO: recursive?
                endif
            enddo

            do x = 1, 3
                if (B>=0 .or. nFS(A)<=1) then
                    prob(x) = XPr(3,x,l, -i, k) / SUM(OKA2P(Genos(l,A), x, :) * PrP(:,1))    
                else if (B==-1) then  ! exclude all FS of A
                    PrY = PrP(:,1)
                    do j=1, nFS(A)
                        if (Genos(l,FSID(j, A))==-9) cycle
                        PrY = PrY * OKA2P(Genos(l,FSID(j, A)), x, :)
                    enddo
                    prob(x) = XPr(3,x,l, -i, k) / SUM(PrY)
                else if (B==-4) then ! exclude both GPs & A
                    prob(x) = XPr(1,x,l, -i,k)*AHWE(x,l) / SUM(OKA2P(Genos(l,A), x, :) * PrP(:,1))
                endif
                if (B>0) then
                    if (Genos(l,B)==-9) cycle
                    prob(x) = prob(x) / SUM(OKA2P(Genos(l,B), x, :) * PrP(:,2)) 
                endif
            enddo
            prob = prob/SUM(prob)
       endif
    endif
endif

end subroutine ParProb

! ##############################################################################################

subroutine Connected(A, kA, B, kB, Con)
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
logical, intent(OUT) :: Con
integer :: i, j, m, nA, nB, AA(maxSibsize), BB(maxSibsize), n

 Con = .FALSE.
if (A==0 .or. B==0) then
    Con = .FALSE.
    return
endif

if (A>0) then
    nA = 1
    AA(1) = A
else
    nA = nS(-A,kA)
    AA(1:nA) = SibID(1:nA, -A, kA)
endif

if (B>0) then
    nB = 1
    BB(1) = B
else
    nB = nS(-B,kB)
    BB(1:nB) = SibID(1:nB, -B, kB)
endif

do j=1, nB
    do i=1, nA
        do m=1,2  
            if (Parent(AA(i), m) < 0) then
                if (Parent(AA(i),m) == Parent(BB(j),m)) then
                    Con = .TRUE.
                    return
                else if(ANY(GpID(:,-Parent(AA(i), m),m) == BB(j))) then
                    Con = .TRUE.
                    return
                else if(ANY(GpID(:,-Parent(AA(i), m),m) < 0)) then
                    do n=1,2
                        if (GpID(n,-Parent(AA(i), m),m) == Parent(BB(j),n) .and. &
                          Parent(BB(j),n)<0) then
                            Con = .TRUE.
                            return
                        endif 
                    enddo
                endif
            endif
            if (Parent(BB(j),m)<0) then
                if (ANY(GpID(:,-Parent(BB(j),m),m) == AA(i))) then
                    Con = .TRUE.
                    return
                else if(ANY(GpID(:,-Parent(BB(j),m),m) < 0)) then
                    do n=1,2
                        if (GpID(n,-Parent(BB(j), m),m) == Parent(AA(i),n) .and. &
                          Parent(AA(i),n)<0) then
                            Con = .TRUE.
                            return
                        endif 
                    enddo
                endif
            endif
        enddo
    enddo
enddo

end subroutine Connected

! ##############################################################################################

subroutine GetAncest(A, k, Anc)
use Global
implicit none

integer, intent(IN) :: A, k
integer, intent(OUT) :: Anc(2, mxA)  ! 32 = 5 generations
integer :: m, j, i

if (A==0) return
Anc = 0
if (A > 0) then  
    if (Sex(A)/=3) then
        Anc(Sex(A),1) = A
    else
        Anc(1, 1) = A
    endif
    Anc(:, 2) = Parent(A, :)
else if (A < 0) then
    Anc(k, 2) = A
endif
if (ALL(Anc(:, 2)==0)) return
do j = 2, mxA/2  
    do m = 1, 2
        i = 2 * (j-1) + m
        if (Anc(m, j) > 0) then
            Anc(:, i) = Parent(Anc(m, j), :)
        else if (Anc(m, j) < 0) then
            Anc(:, i) = GpID(:, -Anc(m, j), m)  
        endif
    enddo
    if (j==2 .and. ALL(Anc(:, 3:4) == 0))  return
    if (j==4 .and. ALL(Anc(:, 5:8) == 0))  return
    if (j==8 .and. ALL(Anc(:, 9:16) == 0))  return
enddo

!if ((A>0 .and. ANY(Anc(:, 2:mxA)==A)) .or. (A<0 .and. ANY(Anc(k,3:mxA)==A))) then
!    call intpr ( "A1 ", -1, Anc(1,1:8), 4)
!    call intpr ( "A2 ", -1, Anc(2,1:8), 4)
!    call rexit("Individual is its own Ancestor ! ")
!endif

end subroutine GetAncest

! ##############################################################################################

subroutine CalcParentLLR  ! Calc parental LLR (vs next most likely relationship)
use Global
implicit none

integer :: i, k, s, CurPar(2), m, nonG(6), CurGP(2), g, curNFS
double precision :: LLA(7), LLtmp(2,2,2), LLCP(2,2)

LR_parent = 999 
do i=1, nInd
!    if (MODULO(i,500)==0) then 
!        call intpr(" ", 1, i, 1)
!    endif
    if (Parent(i,1)==0 .and. Parent(i,2)==0) cycle
    
    CurPar = Parent(i,:)
    curNFS = nFS(i)
    Parent(i,:) = 0
    do k=1,2  ! remove i from sibgroup
        if (CurPar(k)<0) then
            call RemoveSib(i, -CurPar(k), k)
        endif
    enddo
    nFS(i) = 1
    call CalcLind(i)
    
    LLtmp = 999
    do m=1,2  ! m=1: no opp. sex parent;  m=2: with opp. sex parent
        if (m==2 .and. (CurPar(1)==0 .or. CurPar(2)==0)) cycle
        do k=1,2  ! mother, father
            if (m==1 .and. CurPar(k) == 0) cycle
            if (m==2) then  ! temp. assign parent 3-k
                call CalcU(i, 3-k, CurPar(3-k), 3-k, LLCP(k,1))
                LLCP(k,1) = LLCP(k,1) - Lind(i)
                Parent(i, 3-k) = CurPar(3-k)
                if (CurPar(3-k)<0) then
                    call DoAdd(i, -CurPar(3-k), 3-k)
                endif
                call CalcLind(i)
                call CalcU(i, 3-k, CurPar(3-k), 3-k, LLCP(k,2))
                LLCP(k,2) = LLCP(k,2) - Lind(i)
            endif
            
            if (CurPar(k) > 0) then
                call calcPair(i, CurPar(k), k, .FALSE., LLA, 1)
                LLtmp(1,k,m) = LLA(1)
                LLtmp(2,k,m) = MaxLL(LLA(2:7))
            else if (CurPar(k) < 0) then
                call CheckAdd(i, -CurPar(k), k, LLA, 7)
                if (m==1) LLA(2) = 333   ! FS does not count here.
                LLtmp(1,k,m) =  MaxLL(LLA(2:3))
                LLtmp(2,k,m) =  MaxLL((/LLA(1), LLA(4:7)/))
            endif
            
            if (m==1) then
                LR_parent(i,k) = LLtmp(1,k,m) - LLtmp(2,k,m)
            else if (m==2) then  
                Parent(i, 3-k) = 0
                if (CurPar(3-k)<0) then
                    call RemoveSib(i, -CurPar(3-k), 3-k)
                endif
                call CalcLind(i)
            endif   
        enddo
    enddo
    if (CurPar(1)/=0 .and. CurPar(2)/=0) then
        LR_parent(i,3) = MINVAL(LLtmp(1,:,2) - MAX(LLtmp(1,:,1), LLtmp(2,:,2)))       
    endif
    
    Parent(i,:) = CurPar  ! restore
    if (ANY(CurPar < 0)) then
        do k=1,2
            if (CurPar(k)<0) then
                call DoAdd(i, -CurPar(k), k)
            endif
        enddo
    else
        nFS(i) = CurNFS
    endif
    call CalcLind(i)
enddo

!parents of dummies (Sibship GPs)
nonG = (/1,2,3,5,6,7/)
do k = 1,2
    do s=1, nC(k)
 !       if (MODULO(s,50)==0) then 
    !        call intpr(" ", 1, i, 1)
 !       endif
        CurGP = GpID(:, s, k)
        GpID(:, s, k) = 0
        call CalcCLL(s,k)
        
        LLtmp = 999
        do m=1,2
            if (m==2 .and. (CurGP(1)==0 .or. CurGP(2)==0)) cycle
            do g=1,2
                if (m==1 .and. CurGP(g) == 0) cycle
                if (m==2) then  ! temp. assign GP 3-g
                    GpID(3-g, s, k) = CurGP(3-g)
                    call CalcCLL(s,k)
                endif
                
                if (curGP(g) > 0) then
                    call checkAdd(CurGP(g),s,k, LLA, 7)  ! B=GP + CurGP(m)_7
                else if (curGP(g) < 0) then
                    call checkMerge(s, -CurGP(g), k, g, LLA, 4)  
                    if (m==1) then
                        call PairUA(-s, CurGP(g), k, g, LLA(4))  ! not counting SA FS with a B
                    endif
                endif
                LLtmp(1,g,m) = LLA(4)
                LLtmp(2,g,m) = MaxLL(LLA(nonG))                 
                if (m==1) then
                    LR_GP(g,s,k) = LLtmp(1,g,m) - LLtmp(2,g,m)
                else if (m==2) then  ! reset to 0
                    GpID(3-g, s, k) = 0
                    call CalcCLL(s,k)
                endif 
            enddo
        enddo
        if (CurGP(1)/=0 .and. CurGP(2)/=0) then
            LR_GP(3,s,k) = MINVAL(LLtmp(1,:,2) - MAX(LLtmp(1,:,1), LLtmp(2,:,2)))       
        endif       

        ! restore
        GpID(:,s,k) = CurGP
        call CalcCLL(s, k)
    enddo
enddo

end subroutine CalcParentLLR

! ##############################################################################################

subroutine RemoveSib(A, s, k)  ! removes individual A from sibship s
use Global
implicit none

integer, intent(IN) :: A, s, k
integer :: u, j, p, curFS, ox, o, h, v

do u=1,nS(s,k)   
    j = SibID(u,s,k)  ! 1st one in FS
    if (nFS(j)==0) cycle
    curFS = 0  ! Not sure what this was for ?
    do v=1, nFS(j)
        if (FSID(v, j) /= A) cycle   ! drop FS: move FS group to next lowest ID
        if (nFS(j)>1 .and. curFS==0) then
            p = MINVAL(FSID(1:nFS(j), j), MASK=(FSID(1:nFS(j), j)/=A))  ! note: p==j if A/=j
            curFS = p
            ox = 2  ! 1st one stays p
            do o=1, nFS(j)
                if (FSID(o,j)==A) cycle
                if (FSID(o,j)==p) cycle
                FSID(ox, p) = FSID(o, j)
                ox = ox+1
            enddo
            nFS(p) = nFS(j)-1
            nFS(A) = 1
            FSID(:,A) = 0
            FSID(1,A) = A
            exit
        endif
    enddo
enddo 

do u=1,nS(s,k)
    if (SibID(u,s,k)==A) then
        do h=u, nS(s, k)-1  ! drop HS
            SibID(h, s, k) = SibID(h+1, s, k)
        enddo
        SibID(nS(s,k), s, k) = 0
    endif
enddo
nS(s,k) = nS(s,k)-1

Parent(A, k) = 0
call calcCLL(s, k)
call CalcLind(A) 
do u=1,nS(s,k)   ! update LL of connected sibships
    j = SibID(u,s,k)
    if (Parent(j,3-k) < 0) then
        call CalcCLL(-Parent(j,3-k), 3-k)
    endif                    
    call CalcLind(j)
enddo
call calcCLL(s, k)
call CalcLind(A)
 
end subroutine RemoveSib

! ##############################################################################################

! @@@@   INPUT & PRECALC PROB.   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ##############################################################################################

subroutine ReadData
use Global
implicit none

integer :: i,j,k,l, minAgeRepro
integer, allocatable, dimension(:) :: SexTmp, ByTmp
integer, allocatable, dimension(:,:) :: GenosR
 character(len=2000) :: GenoFileName, LifehistFileName, dumC !, AgePriorFileName

! print *, "Reading specs ... "
open (unit=1,file="SequoiaSpecs.txt",status="old")
read (1,*) dumC, GenoFileName      ! 1
read (1,*) dumC, LifehistFileName
read (1,*) dumC, nInd              ! 3
read (1,*) dumC, nIndLH
read (1,*) dumC, nSnp              ! 5 
read (1,*) dumC, Er
read (1,*) dumC, MaxMismatch        ! 7
read (1,*) dumC, thLR           
read (1,*) dumC, thLRrel          ! 9    
read (1,*) dumC, minAgeRepro
read (1,*) dumC, ParentageFileName  ! 11
read (1,*) dumC, AgePriorFileName   ! 12
read (1,*) dumC, nAgeClasses
read (1,*) dumC, maxSibSize         ! 14
read (1,*) dumC, Nrounds
read (1,*) dumC, PedigreeFileName   ! 16
read (1,*) dumC, DumPrefix(1)    
read (1,*) dumC, DumPrefix(2)       ! 18 
 close (1)
 
 mxA = 32  ! max no. ancestors considered 
 
!=================
 
allocate(GenosR(nInd,nSnp))
allocate(Genos(nSnp, nInd))
Genos = -9
allocate(Id(nInd))
Id = "NA"

open (unit=101,file=trim(GenoFileName),status="old")
!read (101, *)  ! header
do i=1,nInd
!    read (101,*) dumI(1), Id(i), dumI(1:4), GenosR(i,:)
    read (101,*)  Id(i), GenosR(i,:)
    do l=1,nSnp
        if (GenosR(i,l)/=-9) then
            Genos(l,i) = GenosR(i,l)+1
        else
            Genos(l,i) = GenosR(i,l)
        endif
    enddo
enddo
 close (101)
 
!=================
! allele frequencies
allocate(AF(nSNP))
do l=1,nSnp
    AF(l)=float(SUM(GenosR(:,l), MASK=GenosR(:,l)/=-9))/(2*nInd)
enddo 

!=================
allocate(SexTmp(nIndLH))
allocate(NameLH(nIndLH))
allocate(ByTmp(nIndLH))

allocate(BY(nInd))
BY=-999
allocate(Sex(nInd))
Sex = 3
allocate(AgeDiff(nInd,nInd))
AgeDiff=999
 
open(unit=103, file=trim(LifehistFileName), status="old")
read(103, *)
do k=1, nIndLH
    read(103,*) NameLH(k), SexTmp(k), ByTmp(k)
enddo
 close(103)

 ! rearrange lifehistory info to same order as genotype file
do i=1,nInd
    do k=1,nIndLH
        if(Id(i)==NameLH(k)) then
            if (BYtmp(k)>0) then
                BY(i) = BYtmp(k)
            endif
            if (SexTmp(k)==1 .or. SexTmp(k)==2) then
                Sex(i)=SexTmp(k)
            endif
        endif
    enddo
enddo
BY1 = MINVAL(BY, MASK=BY>0)
WHERE (BY /= -999) BY = BY - BY1 + 1


!=================
do i=1, nInd
    do j=1, nInd
        if (BY(i)/=-999 .and. BY(j)/=-999) then
            AgeDiff(i,j) = BY(i) - BY(j)   ! if >0, then j older than i
        endif   
    enddo
enddo

!=================
allocate(AgePriorM(nAgeClasses, 8))
open(unit=102, file=trim(AgePriorFileName), status="old")
read(102,*)   ! header
do k=1,nAgeClasses
    read(102,*) AgePriorM(k,:)
enddo 
 close(102)
 
maxAgePO = 0  ! maximum PO age difference (needed for dummy parents)
do k = 1, nAgeClasses
    if (AgePriorM(k, 7)>0 .or. AgePriorM(k, 8)>0) then
        maxAgePO = k - 1  ! first k is agediff of 0
    endif
enddo
!=================
! allocate arrays
allocate(Lind(nInd))
Lind = 0
allocate(LindX(nSnp, nInd))
LindX = 0
allocate(FSID(MaxSibSize, nInd))
FSID(1, :) = (/ (i, i=1, nInd) /)
allocate(NFS(nInd))
NFS = 1
allocate(DumP(3,nSnp, nInd/2,2))
DumP = 0
allocate(XPr(3,3,nSNP, nInd/2,2))
XPr = 0  
allocate(GpID(2, nInd/2,2))
GpID = 0 
allocate(ParentName(nInd,2))
ParentName = "NA"

deallocate(SexTmp)
deallocate(ByTmp)
deallocate(GenosR)
 
end subroutine ReadData

! ##############################################################################################

subroutine ReadParents
use Global
implicit none

integer :: i,j, dumI(2)
real :: dumR(3)
 character(len=30) :: OffName(nInd), dumC(2)

allocate(Parent(nInd,2))
Parent = 0

open(unit=103, file=trim(ParentageFileName), status="old")
read(103,*)   ! header
do i=1,nInd
    read(103,*) OffName(i), dumC(1:2), dumR(1:3), dumI(1:2), j, Parent(i,1), Parent(i,2)
    if (OffName(i) /= Id(i) .or. j/=i) then 
        call rexit("Parentage file order differs from genotype file")
    endif
enddo
 close(103)

do i=1, nInd
    if (Sex(i)==3) then
        if (ANY(Parent(:,1) == i)) then
            Sex(i) = 1
        else if (ANY(Parent(:,2) == i)) then
            Sex(i) = 2
        endif
    endif
enddo
 
end subroutine ReadParents

! ##############################################################################################

subroutine WriteDummies
use Global
implicit none

integer :: s, k, n, y, EstDumBY(3,nInd/2,2), CI(2), mx
character(len=3) :: headTmp
character(len=4) :: DumNameTmp
character(len=5), allocatable, dimension(:) :: OffHeader
character(len=30) :: DumName(nInd/2,2), GpName(2,nInd/2,2)
double precision :: DBYP(nAgeClasses +maxAgePO), cumProp

EstDumBY = -9
do k=1,2
    do s=1,nC(k)
 !       if (MODULO(s,20)==0) then 
    !        call intpr(" ", 1, i, 1)
 !       endif
        call CalcDumBY(s,k)
        mx = MAXLOC(DumBY(:, s, k), DIM=1)
        DBYP = 10**DumBY(:, s, k) / SUM(10**DumBY(:, s, k))
        cumProp = DBYP(mx)
        CI = (/ mx, mx /)
        do y=1, nAgeClasses +maxAgePO 
            if (cumProp > 0.95) then
                EstDumBY(1,s,k) = mx  
                EstDumBY(2,s,k) = CI(1)
                EstDumBY(3,s,k) = CI(2)
                exit
            endif        
            if (CI(1) > 1 .and. CI(2) < nAgeClasses +maxAgePO) then
                if (DBYP(CI(1)-1) > DBYP(CI(2)+1)) then
                    CI(1) = CI(1)-1
                    cumProp = cumProp + DBYP(CI(1))
                else 
                    CI(2) = CI(2)+1
                    cumProp = cumProp + DBYP(CI(2))
                endif
            else if (CI(1) > 1) then
                CI(1) = CI(1)-1
                cumProp = cumProp + DBYP(CI(1))
            else if (CI(2) < nAgeClasses +maxAgePO) then
                CI(2) = CI(2)+1
                cumProp = cumProp + DBYP(CI(2))
            endif
        enddo
    enddo
enddo
EstDumBY = EstDumBY + BY1 -1 -maxAgePO
    
! TODO: this is a copy from code in Program Sibships
GpName = "NA"
do k=1,2
    do s=1,nC(k)
        write(DumNameTmp, '(i4.4)') s
        DumName(s,k) = trim(DumPrefix(k))//DumNameTmp
    enddo
enddo
do k=1,2
    do s=1,nC(k)
        do n=1,2
            if (GpID(n,s,k)>0) then
                GpName(n,s,k) = Id(GpID(n,s,k))
            else if (GpID(n,s,k)<0) then
                GpName(n,s,k) = DumName(-GpID(n,s,k),n) 
            endif
        enddo
    enddo
enddo

allocate(OffHeader(MAXVAL(ns)))    
do s=1, MAXVAL(ns)
    write(headTmp, '(i3.3)') s
    OffHeader(s) = "O_"//headTmp
enddo

open (unit=301,file="DummyParents.txt",status="replace")
    write (301,'(3a20, a4, 3a8, 500a10)')  "ID", "Dam", "Sire", "Sex", "est.BY", "min.BY", "max.BY", "NumOff", OffHeader
    do k=1,2
        do s=1,nC(k)
            write (301,'(3a20, i4, 4i8, " ", 500a10)') DumName(s,k), GpName(:, s,k), k, EstDumBY(:, s, k), &
              nS(s,k), ID(SibID(1:nS(s,k), s, k))
        enddo
    enddo
close (301)

deallocate(OffHeader)

end subroutine WriteDummies

! ##############################################################################################

subroutine PrecalcProbs
use Global
implicit none

integer :: h,i,j,k,l,m
double precision :: OjA(3,3,nSnp), Tmp1(3), Tmp2(3,3)

!###################
! arrays of inheritance prob. conditional on 0, 1 or both parental observed genotypes
! note: all based on observed genotypes, summed over all possible actual genotypes
allocate(AHWE(3,nSnp))
allocate(OHWE(3,nSnp))

! Prob. observed conditional on actual
OcA(1, 1:3) = (/ 1-Er, Er/2, 0.0D0 /)   ! obs=0
OcA(2, 1:3) = (/ Er, 1-Er, Er /)    ! obs=1
OcA(3, 1:3) = (/ 0.0D0, Er/2, 1-Er /)   ! obs=2

! probabilities actual genotypes under HWE
do l=1,nSnp
    AHWE(1,l)=(1 - AF(l))**2 
    AHWE(2,l)=2*AF(l)*(1-AF(l)) 
    AHWE(3,l)=AF(l)**2 
enddo

! joined probabilities actual & observed under HWE
do l=1,nSnp
    do i=1,3    ! obs
        do j=1,3    ! act
            OjA(i, j, l) = OcA(i,j) * AHWE(j, l)
        enddo
    enddo
enddo

! marginal prob. observed genotypes
do l=1,nSnp
    do i=1,3
        OHWE(i, l) = SUM(OjA(i, 1:3, l))
    enddo
enddo

! ########################
! inheritance conditional on 1 parent
allocate(AKAP(3,3,nSnp))
allocate(OKAP(3,3,nSnp))
allocate(OKOP(3,3,nSnp))

do l=1,nSnp
    AKAP(1, 1:3, l) = (/ 1-AF(l), (1-AF(l))/2, 0.0D0 /)
    AKAP(2, 1:3, l) = (/ AF(l), 0.5D0, 1-AF(l) /)
    AKAP(3, 1:3, l) = (/ 0D0, AF(l)/2, AF(l) /)
enddo

do l=1,nSnp
    do i=1,3  ! obs offspring
        do j=1,3    ! act parent
            Tmp1=0
            do k=1,3    ! act offspring
                Tmp1(k) = OcA(i,k) * AKAP(k,j,l)
            enddo
            OKAP(i,j,l) = SUM(Tmp1)
        enddo
    enddo
enddo

do l=1,nSnp
    do i=1,3  ! obs offspring
        do j=1,3    ! obs parent
            Tmp2=0
            do k=1,3    ! act offspring
                do m=1,3    ! act parent
                    Tmp2(k,m) = OcA(i,k) * OcA(j,m) * AKAP(k,m,l)
                enddo
            enddo
            OKOP(i,j,l) = SUM(Tmp2) 
        enddo
    enddo
enddo

! #########################
! inheritance conditional on both parents

AKA2P(1,1,:) = (/ 1.0, 0.5, 0.0 /)
AKA2P(1,2,:) = (/ 0.5, 0.25, 0.0 /)
AKA2P(1,3,:) = (/ 0.0, 0.0, 0.0 /)

AKA2P(2,1,:) = (/ 0.0, 0.5, 1.0 /)
AKA2P(2,2,:) = (/ 0.5, 0.5, 0.5 /)
AKA2P(2,3,:) = (/ 1.0, 0.5, 0.0 /)

AKA2P(3,1,:) = (/ 0.0, 0.0, 0.0 /)
AKA2P(3,2,:) = (/ 0.0, 0.25, 0.5 /)
AKA2P(3,3,:) = (/ 0.0, 0.5, 1.0 /)

do i=1,3  ! obs offspring
    do j=1,3    ! act parent 1
        do h=1,3    !act parent 2
            Tmp1=0
            do k=1,3    ! act offspring
                Tmp1(k) = OcA(i,k) * AKA2P(k,j,h) 
            enddo
            OKA2P(i,j,h) = SUM(Tmp1)
        enddo
    enddo
enddo

!=================
allocate(PHS(3,3,nSnp))
allocate(PFS(3,3,nSnp))
do l=1,nSnp
    do i=1,3  ! obs offspring 1
        do j=1,3    ! obs offspring 2
            Tmp1=0
            Tmp2=0
            do m=1,3    !act shared parent 
                Tmp1(m) = OKAP(i,m,l) * OKAP(j,m,l) * AHWE(m,l)
                do h=1,3
                    Tmp2(m,h) = OKA2P(i,m,h) * OKA2P(j,m,h) * AHWE(m,l) * AHWE(h,l)
                enddo
            enddo
            PHS(i,j,l) = SUM(Tmp1) / (AHWE(i,l) * AHWE(j,l))
            PFS(i,j,l) = SUM(Tmp2) / (AHWE(i,l) * AHWE(j,l))
        enddo
    enddo
enddo

allocate(LindG(3, nSnp, nInd))  ! used when missing genotype
do l=1,nSnp
    do i=1,3
        LindG(i,l,:) = AHWE(i,l)  
    enddo
enddo

deallocate(AF)

end subroutine PrecalcProbs

! ##############################################################################################

subroutine DeAllocAll
use Global
implicit none

! allocated in ReadData
if (allocated(Sex)) deallocate(Sex)
if (allocated(BY)) deallocate(BY)
if (allocated(PairType)) deallocate(PairType)
if (allocated(nFS)) deallocate(nFS)

if (allocated(Genos)) deallocate(Genos)
if (allocated(AgeDiff)) deallocate(AgeDiff)
if (allocated(Parent)) deallocate(Parent)
if (allocated(OppHomM)) deallocate(OppHomM)
if (allocated(nS)) deallocate(nS)
if (allocated(PairID)) deallocate(PairID)
if (allocated(FSID)) deallocate(FSID)

if (allocated(SibID)) deallocate(SibID)
if (allocated(GpID)) deallocate(GpID)

if (allocated(Lind)) deallocate(Lind)
if (allocated(PairDLLR)) deallocate(PairDLLR)
if (allocated(AF)) deallocate(AF)

if (allocated(AHWE)) deallocate(AHWE)
if (allocated(OHWE)) deallocate(OHWE) 
if (allocated(LLR_O)) deallocate(LLR_O)
if (allocated(LindX)) deallocate(LindX)
if (allocated(LR_parent)) deallocate(LR_parent)
if (allocated(AgePriorM)) deallocate(AgePriorM)
if (allocated(CLL)) deallocate(CLL)

if (allocated(AKAP)) deallocate(AKAP)
if (allocated(OKOP)) deallocate(OKOP)
if (allocated(OKAP)) deallocate(OKAP)
if (allocated(LR_GP)) deallocate(LR_GP)
if (allocated(LindG)) deallocate(LindG)
if (allocated(PHS)) deallocate(PHS)
if (allocated(PFS)) deallocate(PFS)
if (allocated(DumBY)) deallocate(DumBY)
if (allocated(DumP)) deallocate(DumP)
if (allocated(XPr)) deallocate(XPr)

if (allocated(Id)) deallocate(Id)
if (allocated(NameLH)) deallocate(NameLH)
if (allocated(ParentName)) deallocate(ParentName)

end subroutine DeAllocAll

! ##############################################################################################

! -9   NA
! 999  NA
! 888  Already assigned
! 777  impossible
! 444  not yet implemented (typically involves inbreeding)
! 222  as likely to go via opposite parent