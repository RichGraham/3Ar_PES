!!CO2-Ar PES
!!Distances in Angstrom and energies in Hartrees
module PES_details
  double precision :: gpRMax = 9.0  
  double precision :: gpRMin = 2.5  
  double precision :: lCO = 1.1632
  double precision :: AngToBohr =  1.8897259885789
  interface PES_GP 
     function PES_GP(xStar) 
       implicit none 
       double precision:: PES_GP
       double precision, dimension(:) ::  xStar
     end function PES_GP
  end interface PES_GP
end module PES_details


module GP_variables
  double precision, allocatable :: alpha (:), lScale(:), xTraining(:,:), xTrainingPerm(:,:,:)
  double precision expVar,NuggVar, gpEmax
  integer :: nDim=3
  integer :: nTraining=161
  integer :: nPerms=6
end module GP_variables


! Test program
use PES_details
implicit none
 
integer k,i, choice

double precision e, PES
double precision, allocatable :: rab
double precision xStar(3)

!!allocate( xStar(nDim))


xStar(1:3)=(/0.3595111400000000068111206,0.2860995399999999855289445, &
     0.2801690199999999908442305/)

xStar(1:3)=(/0.3331502100000000021751134, 0.3260882999999999976026288, &
     0.2039624999999999910293980/)

xStar(1:3)=(/0.3992589799999999855550925, 0.3779625200000000240230236, &
     0.2954577199999999792545680/)

xStar(1:3)=(/0.3004692999999999947213780,0.2998323600000000199727879, &
     0.1819913700000000134071598/)

xStar(1:3)=(/3.4364261e-01 ,  1.1842219e-01 ,  8.8071872e-02/)

xStar(1:3)=(/3.4364261e-01 ,  2.3582520e-01,   1.3985175e-01/)

xStar(1:3)=(/3.4364261e-01  , 2.8349597e-01 ,  1.5534253e-01/)

!xStar(1:3)=(/0.2765329200000000153814028 ,0.2449013399999999951450746, &
!     0.1726287499999999974775733/)

call load_GP_Data
e=PES_GP( xStar)
write(6,*)e

!call EvsE
call fixedAngleSlice

end
!


subroutine EvsE()
  use PES_details
  use GP_variables
  implicit none
  double precision, allocatable:: rab(:), xStar(:)
  double precision dum(9), NonAdd, funcVal, PES, RMSE, mean
  integer i,j, count

  allocate (rab(nDim), xStar(nDim) )

  !====Read test data and compute error====
  open (unit = 7, file = "hiQ_test1_CCSDT_rInv-4000.lhc")
  open (unit=15, file="PES_Err.dat ", status='replace')

  
  RMSE=0
  mean=0
  count=0
  do i=1,4081 !!10078
        read(7,*) dum

     do j=1,3
        xStar(j)=dum(j)
     end do
     NonAdd=dum(9)

     rab(:) = 1/xStar(:)
     funcVal = PES(rab)
     if( dum(4)<0.005 .AND. dum(5)<0.005 .AND. dum(6)<0.005 &
        !.AND. Sqrt((NonAdd-funcVal)**2)/0.005*100> 0.03 &
        ) then
        write(15,*), i, rab(1), rab(2), rab(3), NonAdd, funcVal, Sqrt((NonAdd-funcVal)**2)
        RMSE = RMSE + Sqrt((NonAdd-funcVal)**2)
        mean = mean + Sqrt((NonAdd)**2)
        count = count +1
     endif
  enddo

  print *,RMSE/(1.0*count), mean/(1.0*count), mean/RMSE, count
     
end subroutine EvsE

subroutine fixedAngleSlice()
  use PES_details
  use GP_variables
  implicit none
  double precision, allocatable:: rab(:), xStar(:)
  double precision dum(8), NonAdd, funcVal, PES,r,e
  integer i,j,itot

  allocate (rab(nDim), xStar(nDim) )
  
  itot=500
  

  open (unit=15, file="PES_Out.dat ", status='replace')
  
  do i=0, itot

     ! specify centre-to-centre separation
     r = (  0.5 + 15.0*i/(1.0*itot) ) 

     rab(1)=2.91
     rab(2)=r-2.91
     rab(3)=r
     
     e=PES( rab)
     !e_GP = PES_GP( xStar)
     write(15,*) r , e , 3.000E-007
     
  enddo
  write(6,*)'Written to file: PES_Out.dat '
  close(15)

end subroutine fixedAngleSlice
  
  
  
subroutine computeDistances(r,beta1, rab)
  use PES_details
  implicit none
  double precision  rab(3), r, beta1
  integer ia, ib, ir,k

  rab(1)=r
  rab(2)= SQRT( (lCO*SIN(beta1))**2 +(lCO*COS(beta1)-r)**2  )
  rab(3)= SQRT( (lCO*SIN(beta1))**2 +(lCO*COS(beta1)+r)**2  )
   
end subroutine

double precision function asymp(rab)
  use PES_details
  implicit none
  double precision rab(3)
  double precision c1,c2,c3, cAr, o1Ar, o2Ar

  c1= (  rab(1)**2 + lCO**2 - rab(2)**2)/2.0/rab(1)/lCO
  
  c2= (  rab(2)**2 + lCO**2 - rab(1)**2)/2.0/rab(2)/lCO
  c3= (  rab(3)**2 + lCO**2 - rab(1)**2)/2.0/rab(3)/lCO


  cAr  = - ( 4.64 * (1 + 3*c1**2)  +  3.30 * (5 - 3*c1**2) ) / (rab(1)*AngToBohr)**6
  o1Ar = - ( 8.69 * (1 + 3*c2**2)  +  4.76 * (5 - 3*c2**2) ) / (rab(2)*AngToBohr)**6
  o2Ar = - ( 8.69 * (1 + 3*c3**2)  +  4.76 * (5 - 3*c3**2) ) / (rab(3)*AngToBohr)**6
  
  asymp= cAr + o1Ar + o2Ar
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Gaussian Process Code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pom


  
subroutine load_GP_Data
  use GP_variables
  use PES_details
  implicit none
  
  double precision, allocatable::  xStar(:)
  integer i,j,k
  double precision :: dum, expVar1, expVar2,expVar3
  character (len=90) :: filename
  integer, allocatable:: perm(:,:)

  allocate (alpha(nTraining), lScale(nDim), xTraining(nDim,nTraining),xTrainingPerm(nDim,nPerms,nTraining), xStar(nDim), &
       perm(nDim,nPerms))

  !====Load hyperparameters====
  write (filename, '( "TrainingData/HyperParams_Symm", I3.3, ".dat" )' )  nTraining
  open (unit = 7, file = filename)
  do i=1,nDim
     read (7,*) dum
     j=int(dum)
     !print *,i,j
     read (7,*) lScale(j)
  end do
  read (7,*) expVar
  read (7,*) NuggVar
  read (7,*) gpEMax
  close(7)
  
  
  do i=1,nDim
     !print *,i,lScale(i)
  end do
  !print *,"HyperParams",expVar,NuggVar, gpEmax

  
  !====Load alpha coefficients====
  write (filename, '( "TrainingData/alpha_Symm", I3.3, ".dat" )' )  nTraining
  open (unit = 7, file = filename)
  do i=1,nTraining
     read (7,*) alpha(i)
     !!print *,"alpha ",i, alpha(i)
  end do
  close(7)


  !====Load training data x values ====
  write (filename, '( "TrainingData/xTraining", I3.3, ".dat" )' )  nTraining
  open (unit = 7, file = filename)
    
  do i=1,nTraining
     read (7,*) xTraining(1,i), xTraining(2,i), xTraining(3,i)
     !print *,xTraining(1,i), xTraining(2,i), xTraining(3,i), xTraining(4,i), xTraining(5,i), &
     !     xTraining(6,i), xTraining(7,i), xTraining(8,i), xTraining(9,i)
  end do
  close(7)

  write (filename, '( "3Ar.sym" )' )
  open (unit = 7, file = filename)
  read(7,*) perm

  !do i=1,nPerms
   !  print *, perm(1,i), perm(2,i), perm(3,i), perm(4,i), perm(5,i), perm(6,i), perm(7,i),perm(8,i), perm(9,i)
  !end do
  
  !! Permute the training vectors
  
  do i=1,nDim
     do j=1,nPerms
        do k=1,nTraining
           xTrainingPerm(i,j,k)=xTraining(perm(i,j),k)
        end do
     end do
  end do

end subroutine load_GP_Data
  
function PES_GP(xStar)
  use GP_variables
  implicit none
  double precision, dimension(:) :: xStar
  double precision:: PES_GP
  integer i,j,k
  double precision kSqExpAllPerms, kSqExpJthPerm, kKernTotal


  !! Non-symmetric way
  !kKernTotal=0
  !do i=1,nTraining
  !   kSqExpAllPerms=1.0;
  !   do k=1,nDim
  !      kSqExpAllPerms=kSqExpAllPerms * &
  !           ( exp( - (xStar(k)-xTraining(k,i))**2 /2.0/lScale(k)**2) )
  !   end do
  !   kKernTotal = kKernTotal + alpha(i)*kSqExpAllPerms
  !end do

  !Symmetric way
  kKernTotal=0
  do i=1,nTraining
     kSqExpAllPerms=0
     do j=1,nPerms
        kSqExpJthPerm=1        
        do k=1,nDim
           kSqExpJthPerm  =  kSqExpJthPerm * &
                ( exp( - (xStar(k)-xTrainingPerm(k,j,i))**2 /2.0/lScale(k)**2) )
        end do !Dimensions (k)
        kSqExpAllPerms = kSqExpAllPerms + kSqExpJthPerm
     end do !Permuations (
     kKernTotal = kKernTotal + alpha(i) * kSqExpAllPerms
  end do !Training points (i)
  
  PES_GP=kKernTotal * expVar
end function PES_GP



function PES( rab )
  !! Takes in rab in Angstrom
  use PES_details
  use GP_variables
  implicit none
  double precision rab(3), xStar(3), asymp
  double precision  PES
  double precision repFactor

  repFactor=1.0
  
  !! For 3 body interactions any of the distances >RMax means set the non-additive part to zero
  if ( ANY( rab > gpRMax )  ) then
     PES=0.0        !!Use asymptotic function (which is zero for non-additive interactions)
     
  else if ( ANY( rab < gpRMin/repFactor)  ) then 
     PES=gpRMin !! Use repulsive approximation function
     
  else !! Use the Guassian Process function
     xStar(:) = 1/rab(:)
     PES = PES_GP( xStar)
  end if

  !PES=gpEmax/3.0* (1.0/rab(1)**12+1.0/rab(2)**12+1.0/rab(3)**12) *gpRMin **12

  
end function PES
