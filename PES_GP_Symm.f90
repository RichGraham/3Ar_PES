!!CO2-Ar PES
!!Distances in Angstrom and energies in Hartrees
module PES_details
  double precision :: gpRMax = 9.0  
  double precision :: gpRMin = 1.5  
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
  integer :: nTraining=148
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


xStar(1:3)=(/0.3914942399999999933513095,0.3269332899999999875184642, &
     0.2852068600000000064831340/)
xStar(1:3)=(/0.2420432499999999875761603, 0.2236785199999999917963578, &
     0.1429743199999999880844825/)
xStar(1:3)=(/0.3601643800000000061878325,0.3275462499999999832489550, &
     0.2331720899999999985219290/)
xStar(1:3)=(/0.3617689800000000177604420, 0.2993179400000000045523052,&
 0.2172357900000000119344890/)

!xStar(1:9)=(/ 0.24927635,  0.29087296,  0.3123213 ,  0.277394  ,  0.30196993, &
!        0.29306223,  0.28166435,  0.27929086,  0.25204695 /)

call load_GP_Data
e=PES_GP( xStar)
!e=PES( rab)
write(6,*)e

end
!

subroutine fixedAngleSlice()
  use PES_details
  implicit none
  double precision rab(9)
  integer i, itot
  double precision  r, beta1, e, e_GP, asymp, PES
    
  itot=500
  beta1 =  0   /180.0*3.14159265359

  open (unit=15, file="PES_Out.dat ", status='replace')
  
  do i=0, itot

     ! specify centre-to-centre separation
     r = (  0.5 + 15.0*i/(1.0*itot) ) 

     call computeDistances(r,beta1,rab)
     
     
     e=PES( rab)
     !e_GP = PES_GP( xStar)
     write(15,*) r , e 
     
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
     print *,i,j
     read (7,*) lScale(j)
  end do
  read (7,*) expVar
  read (7,*) NuggVar
  read (7,*) gpEMax
  close(7)
  
  
  do i=1,nDim
     print *,i,lScale(i)
  end do
  print *,"HyperParams",expVar,NuggVar, gpEmax

  
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
  double precision rab(9), xStar(9), asymp
  double precision  PES
  double precision repFactor

  repFactor=1.0
  
  if( rab(1) > gpRMax  .AND.  rab(2) > gpRMax .AND.  rab(3) > gpRMax &
       ) then !!Use asymptotic function
     PES = asymp(rab)
     
  else if (rab(1) < gpRMin/repFactor  .OR.  rab(2) < gpRMin/repFactor  .OR.  rab(3) < gpRMin/repFactor &
       ) then !! Use repulsive approximation function
     PES=gpEmax* (1.0/rab(1)**12+1.0/rab(2)**12+1.0/rab(3)**12) *gpRMin **12
     
  else !! Use the Guassian Process function
     xStar(:) = 1/rab(:)
     PES = PES_GP( xStar)
  end if

  !PES=gpEmax/3.0* (1.0/rab(1)**12+1.0/rab(2)**12+1.0/rab(3)**12) *gpRMin **12

  
end function PES
