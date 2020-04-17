!SOLVES LAPLACE EQUATION IN THE CONTEXT OF CAPACITANCE AND CREATES GRAPHICS
!OUTPUTS RESULTS AS TEXT FILES AND IMAGES
!NOTE- CERTAIN SUB-ROUTINES AND CALLS WILL ONLY WORK WITH SUITABLE PARAMETERS. NO DATA VALIDATION HAS BEEN CODED INTO THIS PROGRAM

program serialProg
implicit none 

    !User defined variables 
    integer, parameter:: iterMax=10000
    real, parameter:: initValue=0.5
    real, parameter:: tol=10.0E-5

    integer, parameter:: n=1440 !**IF N IS CHANGED, OUTPUT FORMATS MUST BE CHANGED ACCORDINGLY**!
    integer, parameter:: row= n
    integer, parameter:: col= n
    real:: array(row,col)

    real, parameter:: boundTL=0
    real, parameter:: boundTR=0
    real, parameter:: boundBL=0
    real, parameter:: boundBR=0

    logical, parameter:: plates=.TRUE.
    logical, parameter:: animate=.FALSE.
    real, parameter:: plateTVal=1
    real, parameter:: plateBVal=-1

    !Co- ords: plateTXMin, plateTXMax, plateTYMin, plateTYMax
        integer:: plateTco(4)=(/600, 840, 660, 690/)
    integer:: plateBco(4)=(/600, 840, 750, 780/) !As above
        integer, parameter:: shift=20 !The degree of the rightward shift of top plate
        integer, parameter:: numImages=30

        !**NOTE**Check the model is suitable for capacitance calculation i.e. plate dimensions, plate oriantation, plate charge etc
        logical, parameter:: capacitance=.FALSE.
        logical, parameter:: plateHoriz=.TRUE.!Horizontal orientation? 
        integer, parameter:: d=8 !This is the distance separating the plates
        real:: capac(30)
        common /global/ capac	

        logical, parameter:: contour=.TRUE.

        !Program variables
        real:: step
        integer:: i
        double precision:: tStart, tFinish
        1 format(1440f7.3, 5x)!!Display array in matrix form
        call cpu_time(tStart)
        if (animate) then
            step=nint((shift/numImages)+0.5)!0.5 ensures always rounded up
            do i=1, numImages

                !Call relivant routines 
                call initArray(row, col, initValue, array, boundTL, boundTR, boundBL,&
                                       & boundBR, plateTVal, plateBVal, plateTco, plateBco, plates)

                call applyAlgo(row, col, array, tol, iterMax, plateTco, plateBco, plates)

                if (capacitance) then
                call calcCapacitanceRatio(row, col, array, plateTco, plateTVal,plateBVal,&
                                                  & d, plateHoriz, i)
                end if 

                call createDiagram(row, col, array, contour)

                !Change array values accroding to shift 
                plateTco(1)=plateTco(1)+step
                plateTco(2)=plateTco(2)+step
                plateBco(1)=plateBco(1)-step
                plateBco(2)=plateBco(2)-step



            end do
        else 
            !Call relivant routines 
            call initArray(row, col, initValue, array, boundTL, boundTR, boundBL,&
                               & boundBR, plateTVal, plateBVal, plateTco, plateBco, plates)

            call applyAlgo(row, col, array, tol, iterMax, plateTco, plateBco, plates)

            !call createDiagram(row, col, array, contour)

            if (capacitance) then
                call calcCapacitanceRatio(row, col, array, plateTco, plateTVal,plateBVal,&
                                                  & d, plateHoriz)
            end if
        end if  
        !call plotCurve()
        call cpu_time(tFinish)

        open(1,file='capData.dat')
        do i=1, 30
            write(1,*)(i+7),capac(i)		
        end do
        close(1)

        write(*,*)"Time:", tFinish-tStart	
end program serialProg

!***************************************************************************
!******************SUB ROUTINES*********************************************
!***************************************************************************

subroutine initArray(row, col, initValue, array, boundTL, boundTR,&
                     & boundBL, boundBR, plateTVal, plateBVal, plateTco, plateBco, plates)

    !Variables
    integer:: row
    integer:: col
    real:: initValue
    real:: boundTL
    real:: boundTR
    real:: boundBL
    real:: boundBR
    real:: plateTVal
    real:: plateBVal
    real:: calcIncr
    real:: incr
    real:: array(row, col)
    integer:: plateTco(4)
    integer:: plateBco(4)
    logical:: plates
    !Counters
    integer:: i, j
    !Formats
    1 format(200f7.3, 5x)!!Display array in matrix form 

    !Set init values 
    outer1: do i=1, row
        inner1: do j=1, col

            array(i,j)=0

        end do inner1
    end do outer1

    !Calculate boundary incriments (Top row) and set
    incr=calcIncr(boundTR, boundTL, real(col-1))
    do i=2, col
        array(1, 1)=boundTL
        array(1, i)=array(1, i-1)+incr
    end do 
    !Calculate boundary incriments (Bottom row) and set
    incr=calcIncr(boundBR, boundBL, real(col-1))
    do i=2, col
        array(row, 1)=boundBL
        array(row, i)=array(row, i-1)+incr
    end do 
    !Calculate boundary incriments (Left col) and set
    incr=calcIncr(boundBL, boundTL, real(row-1))
    do i=2, row
        array(1, 1)=boundTL
        array(i, 1)=array(i-1, 1)+incr
    end do 
    !Calculate boundary incriments (Right col) and set
    incr=calcIncr(boundBR, boundTR, real(row-1))
    do i=2, row
        array(1, col)=boundTR
        array(i, col)=array(i-1, col)+incr
    end do 

    !Set internal values 
    outer: do i=2, row-1
        inner: do j=2, col-1

            array(i,j)=initValue

        end do inner
    end do outer		

    if (plates) then
        !Set plateT values 
        !Co- ords: plateTXMin, plateTXMax, plateTYMin, plateTYMax
        rowLoop2: do i=(plateTco(3)), plateTco(4)
            colLoop2: do j=plateTco(1), plateTco(2)

                array(i,j)=plateTval

            end do colLoop2
        end do rowLoop2

        !Set plateB values 
        rowLoop: do i=(plateBco(3)), plateBco(4)
            colLoop: do j=plateBco(1), plateBco(2)

                array(i,j)=plateBval

            end do colLoop
        end do rowLoop
    end if 

    open(1,file='output.txt')
    Write(1,*)'***INITIALISED VALUES***'
    write(1,1)((array(i,j), j=1,col), i=1,row)
    close(1)

return 
end subroutine initArray
!*************************************************************
subroutine applyAlgo(row, col, array, tol, iterMax, plateTco, plateBco, plates)

    !Variables
    integer:: row, col 
    integer:: iterMax
    integer:: x
    real:: array(row, col)
    real:: tol
    real:: oldVal, newVal
    real::diff
    integer:: plateTco(4)
    integer:: plateBco(4)
    logical:: withinTol
    logical:: plates
    !Counters
    integer:: i, j, k, iter
    !Formats
    1 format(200f7.3, 5x)!!Display array in matrix form

    !Initialise variables required for loop
    diff=tol+1
    iter=0
    if (plates) then
        !Loop internal points until change within tol
        numIt: do while(diff>=tol .AND. iter<iterMax)	
            arrayRow: do i=2, row-1
                arrayCol: do j=2, col-1

                    !If not on either plate 
                    !Co- ords: plateTXMin, plateTXMax, plateTYMin, plateTYMax
                    if(.NOT.(((j.ge.plateTco(1).and.j.le.plateTco(2)).AND.&
                                                  & (i.ge.plateTco(3).and.i.le.plateTco(4))).OR.((j.ge.plateBco(1)&
                                                                                        & .and.j.le.plateBco(2)).AND.(i.ge.plateBco(3) .and.i.le.&
                                                                      & plateBco(4)))))then 

                        oldVal= array(i,j)
                        newVal= 0.25*(array(i+1,j)+array(i-1,j)+array(i,j-1)+array(i,j+1))

                        if(abs(newVal-oldVal).ge.diff)then
                        diff= abs(newVal-oldVal)
                        end if 

                        array(i,j)=newVal
                    end if

                end do arrayCol
            end do arrayRow
            iter=iter+1
        end do numIt
    else
        numIt2: do while(diff>=tol .AND. iter<iterMax)	
            arrayRow2: do i=2, row-1
                arrayCol2: do j=2, col-1

                    oldVal= array(i,j)
                    newVal= 0.25*(array(i+1,j)+array(i-1,j)+array(i,j-1)+array(i,j+1))

                    if(abs(newVal-oldVal).ge.diff)then
                    diff= abs(newVal-oldVal)
                    end if 

                    array(i,j)=newVal

                end do arrayCol2
            end do arrayRow2
            iter=iter+1
        end do numIt2
    end if
    !Dsiplay message if max iteration was reached
    if (iter==iterMax) then
        write(*,*)'Max iteration was reached, results may not be accurate to defined tolerance!'
    end if

    open(1,file='output.txt', ACCESS='APPEND')
    write(1,*)
    write(1,*)'***OUTPUT/ FINAL VALUES***'
    write(1,1)((array(i,j), j=1,col), i=1,row)
    !write(1,*)
    !write(1,*)'*****************ITERATIONS:',iter
    close(1)

    write(*,*)"Iterations:", iter

return
end subroutine applyAlgo	
!*************************************************************
subroutine createDiagram(row, col, array, contour)
    use DISLIN
    implicit none

    !Variables 
    integer:: row, col 
    real:: array(row, col)
    real:: arrayT(col,row)
    real:: xArray(col)
    real:: yArray(row)
    real:: stepX, stepY, stepZ
    real:: ctrStep
    logical:: contour

    !counters
    integer:: i, j

    !Transform array for display and create x and y co- ordinate arrays
    do i=1, col
        xArray(i)=i-1
        do j=1, row
            yArray(j)=j-1
            arrayT(i,j)=array(row-(j-1),i)				
        end do
    end do 

    call metafl('PNG')
    call scrmod('REVERSE')
    call DISINI()

    call titlin('Colour plot of heat distribution within plate',-2)
    call name('Electric potential scale','Z')

    call autres(col,row)
    call axspos(300,1850)
    call ax3len(2200,1400,1400)
    call labdig(-2,'XYZ')

    stepX= (col/10)
    stepY= (row/10)
    stepZ= (maxval(arrayT)-minval(arrayT))/20

    call graf3(0.,real(col),0.,stepX,0.,real(row),0.,stepY,minval(arrayT),&
                   &maxval(arrayT),minval(arrayT),stepZ)!Coloured axis system
    call crvmat(arrayT,col,row,100,100)

    if (contour) then
        do i=2,10
            ctrStep=-1.+(i-1)*0.2
            if(ctrStep.ne.0)then
            call labels('FLOAT','CONTUR')
            call contur(xArray,col,yArray,row,arrayT,ctrStep)
            end if
        end do
    end if
    call DISFIN()

return 
end subroutine createDiagram	
!***********************************************************************************
subroutine calcCapacitanceRatio(row, col, array, plateTco, plateTVal,&
                                & plateBVal, d, plateHoriz, j)
implicit none

    !Variables
    real, parameter:: epsilon_0 = 8.854187817D-12
    integer:: plateTco(4)!Co- ords: plateTXMin, plateTXMax, plateTYMin, plateTYMax
    integer:: row, col 
    integer:: d
    real:: array(row, col)
    real:: plateTVal
    real:: plateBVal
    real:: potSumPlate=0
    real:: potSumOuter=0
    real:: voltage
    real:: capRatio
    real:: l
    logical:: plateHoriz
    real:: capac(30)
    common /global/ capac

    !Counters
    integer:: i, j

    potSumPlate=0
    potSumOuter=0
    !Loop the boundary of the plate
    top: do i=plateTco(1), plateTco(2)
        potSumPlate=potSumPlate+array(plateTco(3),i)
        potsumOuter=potsumOuter+array(plateTco(3)-1,i)
    end do top
    bottom: do i=plateTco(1), plateTco(2)
        potSumPlate=potSumPlate+array(plateTco(4),i)
        potsumOuter=potsumOuter+array(plateTco(4)+1,i)
    end do bottom
    left: do i=plateTco(3), plateTco(4)
        potSumPlate=potSumPlate+array(i,plateTco(1))
        potsumOuter=potsumOuter+array(i,plateTco(1)-1)
    end do left
    right: do i=plateTco(3), plateTco(4)
        potSumPlate=potSumPlate+array(i,plateTco(2))
        potsumOuter=potsumOuter+array(i,plateTco(2)+1)
    end do right

    !Subtract repeted values
    potSumPlate=potSumPlate-4*plateTval

    voltage= abs(plateTval-plateBval)

    !Calc length of plate
    if (plateHoriz) then 
        l=plateTco(2)-plateTco(1)
    else
        l=plateTco(4)-plateTco(3)
    end if 

    !Cal capacitance ratio
    capRatio=(d/(voltage*l))*abs(potSumOuter-potSumPlate)
    capac(j)=capRatio 

    open(3,file='capNumbers.txt', ACCESS='APPEND')
    write(3,*)
    write(3,*)'***CAPACITANCE***', j
    write(3,*)'Cap Ratio= ',capac(j)
    write(3,*)'d=',d,'pot outer=', potSumOuter, 'pot plate=', potSumPlate,&
        & 'diff=', abs(potSumOuter-potSumPlate), 'Voltage=', voltage, "l=", l
    close(3)

return
end subroutine calcCapacitanceRatio
subroutine plotCurve()
USE DISLIN
implicit none

    !Variables
    integer:: i, j, IC
    real:: capac(30)
    real:: yCo(30)
    real:: xCo(30)
    common /global/ capac

    do i=1, 30
        xCo(i)=i
        yCo(i)=capac(i)		
    end do

    CALL METAFL('PNG')
    CALL SCRMOD('REVERS')
    CALL DISINI()
    CALL PAGERA()
    CALL COMPLX()
    CALL AXSPOS(450,1800)
    CALL AXSLEN(2200,1200)

    CALL NAME('Separation','X')
    CALL NAME('Ratio','Y')

    call labdig(-2,'XY')
    CALL TICKS(10,'XY')

    CALL TITLIN('Capacitance ratio',1)
    CALL TITLIN('real/ infinite',3)

        IC=INTRGB(0.95,0.95,0.95)
    CALL AXSBGD(IC)

    CALL GRAF(1.,20.,40.,4.,0.,1.,0.,0.2)
        CALL SETRGB(0.7,0.7,0.7)
    CALL GRID(1,1)

    CALL COLOR('FORE')
    CALL TITLE()
        CALL COLOR('RED')
    CALL CURVE(xCo,yCo, 30)
    CALL DISFIN()


return
end subroutine plotCurve
!*************************************************************
!**********************FUNCTIONS******************************
!*************************************************************

function calcIncr(bound1, bound2, d) result(incr)

    real:: bound1, bound2, d
    real:: incr

    incr=(bound1-bound2)/d

return
end function calcIncr
