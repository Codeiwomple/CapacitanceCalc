!Parallel code to implement the Jacobi iteration on a square grid with a boundary of constant value 0, 
!and areas of user defined constant value to represent parallel plates. Plate A has a constant value 
!of -1 and B+1. Code is only designed and tested to run on a square number of cores. Values outputted into a text file.

program parRewrite
use mpi
implicit none

!User defines variables
integer, parameter:: iterMax= 10000
real, parameter:: initValue=0.5
real, parameter:: tol=10.0E-5

integer, parameter:: dims=480

integer:: plateAco(4)=(/200, 280, 220, 230/) !Co- ords: plateTXMin, plateTXMax, plateTYMin, plateTYMax
integer:: plateBco(4)=(/200, 280, 250, 260/) !As above

!Program variables
integer:: points
real, dimension(:,:), allocatable:: procArray
integer:: gridDim
integer:: procRow, procCol
double precision:: tStart, tFinish, tElap
logical:: containsPlate=.FALSE.
common /global/ containsPlate
!counters
integer:: i
!Setup MPI environment
integer:: ierr, numProcs, ProcID
!Init MPI
call MPI_INIT (ierr)
!Get number of processes/ cores requested
call MPI_COMM_SIZE (MPI_COMM_WORLD, numProcs, ierr)
!Get rank of process
call MPI_COMM_RANK (MPI_COMM_WORLD, procID, ierr)

tStart=MPI_WTIME()
gridDim=sqrt(real(numProcs))
points=dims/gridDim
allocate(procArray(points+2,points+2))

!Determine which domain belongs to this process
do i=1, gridDim
    if((procID+1 .le. i*gridDim) .and. (procID+1 .ge. (i*gridDim-(gridDim-1)))) then
        procRow=i
    end if
end do
do i=1, gridDim
    if(mod(procID+1-i, gridDim) .eq. 0) then
        procCol=i
    end if
end do

call initArray(points, procArray, initValue, plateAco, plateBco, procRow, procCol, gridDim)
call applyAlgo(points, procArray, tol, iterMax, plateAco, plateBco, gridDim, numProcs, procRow, procCol, procId, dims, ierr)

deallocate(procArray)

tFinish= MPI_WTIME()
tElap=tFinish-tStart
if (procID.eq.0) then
open(1,file='finalValues.txt', ACCESS='APPEND')
write(1,*)'Time:', tElap
close(1)
end if 
call MPI_FINALIZE(ierr)
stop
end Program parRewrite

!***************************************************************************
!******************SUB ROUTINES*********************************************
!***************************************************************************

subroutine initArray(points, procArray, initValue, plateAco, plateBco, procRow, procCol, gridDim)
implicit none

!Parameters
integer:: points
real:: procArray(points+2, points+2)
real:: initValue
integer:: plateAco(4)
integer:: plateBco(4)
integer:: procRow
integer:: procCol
integer:: gridDim
!Variabes
integer:: globalRow, globalCol
logical:: containsPlate
common /global/ containsPlate
!counters
integer:: i, j

!Set all values to initValue 
do i=1, points+2
    do j=1, points+2
        procArray(i,j)=initValue
    end do
end do

!Set boundary conditions
if (procRow.eq.1) then
    !Set top values to 0 for ICs
    do i=2, points+1
        procArray(2,i)=0
    end do
end if

if (procRow.eq.gridDim) then
    !Set bottom values to 0 for ICs
    do i=2, points+1
        procArray(points+1,i)=0
    end do
end if

if (procCol.eq.1) then
    !Set left values to 0 for ICs
    do i=2, points+1
        procArray(i,2)=0
    end do
end if

if (procCol.eq.gridDim) then
    !Set right values to 0 for ICs
    do i=2, points+1
        procArray(i,points+1)=0
    end do
end if

!Set plate values 
!Loop 'internal' points
do i=0, points-1
    do j=0, points-1
        !Calculate global co- ordinate of current point 
        globalRow=procRow*points-(points-1)+i
        globalCol=procCol*points-(points-1)+j

        !If on plate A
        if ((globalCol.ge.plateAco(1).and.globalCol.le.plateAco(2)).AND.&
            & (globalRow.ge.plateAco(3).and.globalRow.le.plateAco(4))) then

            procArray(i+2,j+2)=-1
            containsPlate=.TRUE.
        end if
        !If on plate B
        if ((globalCol.ge.plateBco(1).and.globalCol.le.plateBco(2)).AND.&
            &(globalRow.ge.plateBco(3).and.globalRow.le.plateBco(4)))then

            procArray(i+2,j+2)=1
            containsPlate=.TRUE.
        end if
    end do
end do

return
end subroutine initArray

!***************************************************************************

subroutine applyAlgo(points, procArray, tol, iterMax, plateAco, plateBco, gridDim, numProcs, procRow, procCol, procID, dims, ierr)
use mpi
implicit none
!Parameters
integer:: points
real:: procArray(points+2, points+2)
real:: recvArray(points+2, points+2)
real:: tol
integer:: iterMax
integer:: plateAco(4)
integer:: plateBco(4)
integer:: gridDim
integer:: numProcs
integer:: procRow, procCol
integer:: procID
integer:: dims
integer:: ierr
!Variables 
real:: oldVal, newVal, diff, recvDiff
real:: maxDiff
real:: diffArray(numProcs)
real:: sendUp(points), sendDown(points), sendLeft(points), sendRight(points)
real:: recvAbove(points), recvBelow(points), recvLeft(points), recvRight(points)
integer:: startCol, endCol, startRow, endRow
integer:: globalRow, globalCol
real:: allProcArrays(numProcs,points+2,points+2)
real:: outArray(dims, dims)
logical:: containsPlate
common /global/ containsPlate
!MPI Variables 
integer:: status(MPI_STATUS_SIZE)
!Counters
integer:: i, j, k, iter
integer:: subDom, minorRow, minorCol, majorRow, majorCol
!Formats
1 format(480f7.3, 5x)!!Display array in matrix form

diff=tol+1
maxDiff=tol+1
iter=0
numIt: do while(iter.le.iterMax.AND. maxDiff.gt.tol)
maxDiff=tol
!Share data
!Set up data sharing arrays to be sent
do i=2, points+1
    sendDown(i-1)=procArray(points+1, i)
end do
do i=2, points+1
    sendUp(i-1)=procArray(2, i)
end do
do i=2, points+1
    sendLeft(i-1)=procArray(i, 2)
end do
do i=2, points+1
    sendRight(i-1)=procArray(i, points+1)
end do

!Send arrays (buffer, count, mpiType, dest, tag, comm, ierr)
call MPI_SEND (sendDown, points+2, MPI_REAL, modulo(procID+gridDim, numProcs), 0, MPI_COMM_WORLD, ierr)
call MPI_SEND (sendUp, points+2, MPI_REAL, modulo(procID-gridDim, numProcs), 1, MPI_COMM_WORLD, ierr)
call MPI_SEND (sendLeft, points+2, MPI_REAL, modulo(procID-1, numProcs), 2, MPI_COMM_WORLD, ierr)
call MPI_SEND (sendRight, points+2, MPI_REAL, modulo(procID+1, numProcs), 3, MPI_COMM_WORLD, ierr)

!Recieve arrays (buffer, count, mpiType, source, tag, comm, status, ierr)
call MPI_RECV (recvAbove, points+2, MPI_REAL, modulo(procID-gridDim, numProcs), 0, MPI_COMM_WORLD, status, ierr)
call MPI_RECV (recvBelow, points+2, MPI_REAL, modulo(procID+gridDim, numProcs), 1, MPI_COMM_WORLD, status, ierr)
call MPI_RECV (recvRight, points+2, MPI_REAL, modulo(procID+1, numProcs), 2, MPI_COMM_WORLD, status, ierr)
call MPI_RECV (recvLeft, points+2, MPI_REAL, modulo(procID-1, numProcs), 3, MPI_COMM_WORLD, status, ierr)

!Compile recieved values into the procedures array
!Add the top
do i=2, points+1
    procArray(1,i)=recvAbove(i-1)
end do
!Add the bottom
do i=2, points+1
    procArray(points+2,i)=recvBelow(i-1)
end do
!Add the left
do i=2, points+1
    procArray(i,1)=recvLeft(i-1)
end do
!Add the right
do i=2, points+1
    procArray(i,points+2)=recvRight(i-1)
end do

!Determine which points need to be looped 
!Init needed variables 
startCol=2
endCol=points+1
startRow=2
endRow=points+1
if (procCol.eq.1) then
    startCol=3
end if
if (procCol.eq.gridDim) then
    endCol=points
end if
if (procRow.eq.1) then
    startRow=3
end if
if (procRow.eq.gridDim) then
    endRow=points
end if

!Loop array and apply Jacobi iteration
If(containsPlate)then
    do i=startRow, endRow
        do j=startCol, endCol

        globalRow=procRow*points-(points-1)+i-2
        globalCol=procCol*points-(points-1)+j-2

plateCheck: if(.NOT.(((globalCol.ge.plateAco(1).and.globalCol.le.plateAco(2)).AND.&
                      & (globalRow.ge.plateAco(3).and.globalRow.le.plateAco(4))).OR.((globalCol.ge.plateBco(1)&
                                                                            & .and.globalCol.le.plateBco(2)).AND.(globalRow.ge.plateBco(3) .and.globalRow.le.&
                                                  & plateBco(4)))))then

    oldVal=procArray(i,j)
    newVal=0.25*(procArray(i+1,j)+procArray(i-1,j)+procArray(i,j+1)+procArray(i,j-1))
    procArray(i,j)=newVal
    diff= abs(newVal-oldVal)

    if(diff.ge.maxDiff) then
        maxDiff=diff
    end if
    end if plateCheck

    end do
    end do
else
    do i=startRow, endRow
        do j=startCol, endCol

            oldVal=procArray(i,j)
            newVal=0.25*(procArray(i+1,j)+procArray(i-1,j)+procArray(i,j+1)+procArray(i,j-1))
            procArray(i,j)=newVal
            diff= abs(newVal-oldVal)

            if(diff.ge.maxDiff) then
                maxDiff=diff
            end if

        end do
    end do
end if

!Send difference to root
if(ProcID.ne.0)then
    call MPI_SEND (maxDiff, 1, MPI_REAL, 0, procID+2000, MPI_COMM_WORLD, ierr)

else!Recieve differences
    diffArray(1)=maxDiff
    do k=2, numProcs 
        call MPI_RECV (recvDiff, 1, MPI_REAL, k-1, (k-1)+2000, MPI_COMM_WORLD, status, ierr)
        diffArray(k)=recvDiff
    end do

    !Set maxDiff to max value in array
    maxDiff=maxVal(diffArray)
end if
!Broadcast maxDiff to other procs
call MPI_BCAST(maxDiff, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
iter=iter+1
end do numIt

!Gather data
if(procID.ne.0)then !If not the root process, send data to it
    !(buffer, count, mpiType, dest, tag, comm, ierr)
    call MPI_SEND (procArray, (points+2)*(points+2), MPI_REAL, 0, procID+1000, MPI_COMM_WORLD, ierr)
else !Collect data in root process
    !Add array from proc 0
    do i=1, points+2
        do j=1, Points+2
            allProcArrays(1,i,j)=procArray(i,j)
        end do
    end do
    !Add the rest
    do k=2, numProcs
        call MPI_RECV (recvArray, (points+2)*(points+2), MPI_REAL, k-1, (k-1)+1000, MPI_COMM_WORLD, status, ierr)
        do i=1, points+2
            do j=1, points+2
            allProcArrays(k,i,j)=recvArray(i,j)
            end do
        end do
    end do
end if

!Output to file
if(procID.eq.0)then
    do majorRow=1,gridDim
    do majorCol=1,gridDim
        subDom=(majorRow-1)*gridDim+majorCol
        i=1
        do minorRow=majorRow*points-(points-1),majorRow*points
        i=i+1
        j=1
        do minorCol=majorCol*points-(points-1),majorCol*points 
        j=j+1
            outArray(minorRow,minorCol)=allProcArrays(subDom,i,j)
        end do
        end do
    end do
    end do

    open(1,file='finalValues.txt', ACCESS='APPEND')
    write(1,*)'***OUTPUT/ FINAL VALUES***'
    write(1,1)((outArray(i,j), j=1,dims), i=1,dims)
    write(1,*)
    write(1,*)'**************************************************************'
    write(1,*)
        if (iter.eq.iterMax)then
        write(1,*)"WARNING MAX ITERATIONS REACHED!"
        end if 
    write(1,*)'Cores:',numProcs, "dim:", dims, "tol:",tol, "iter", iter
    close(1)
end if

return
end subroutine applyAlgo
