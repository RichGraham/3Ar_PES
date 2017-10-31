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
  integer :: nDim=9
  integer :: nTraining=28
  integer :: nPerms=8
end module GP_variables


! Test program
use PES_details
implicit none
 
integer k,i, choice

double precision rab(9),e, PES, xStar(9)

xStar(1:9)=(/0.20971571,  0.22605512,  0.22877215,  0.20158676,  0.2298365 , &
     0.2503203 ,  0.18510657,  0.2182286 ,  0.25385964/)
xStar(1:9)=(/ 0.3018619 ,  0.33131717,  0.31699997,  0.30609471,  0.34708315, &
     0.34030303,  0.27652502,  0.3131383 ,  0.31583372 /)

xStar(1:9)=(/ 0.16918005,  0.20963065,  0.27441895,  0.15582496,  0.18869334, &
        0.23773219,  0.14123042,  0.16638506,  0.20074532 /)

xStar(1:9)=(/ 0.37918623,  0.43618246,  0.39957494,  0.32684715,  0.42808256, &
     0.48455317,  0.26286132,  0.34577986,  0.44897025 /)

xStar(1:9)=(/0.3791862300000000129962530, 0.4361824599999999940713735, &
     0.3995749399999999895705116, 0.3268471499999999752006374, 0.4280825600000000008549250,&
     0.4845531699999999775130277, 0.2628613200000000094114228, 0.3457798599999999944465401, &
     0.4489702500000000151558766/)

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
  !write (filename, '( "TrainingData/myTest.dat" )' )  
  !open (unit = 7, file = filename)
  !! Older symmetric way
  !Only need to read some as others are tied.
  !!read (7,*) lScale(1),lScale(2), lScale(5), expVar3, expVar2, expVar1,NuggVar, gpEmax
  
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
  
  !Copy over the tied values
  lScale(3) = lScale(1) ! 2=0
  lScale(9) = lScale(1) ! 8=0
  lScale(7) = lScale(1) ! 6=0
  
  lScale(4) = lScale(2) ! 3=1
  lScale(6) = lScale(2) ! 5=1
  lScale(8) = lScale(2) ! 7=1
  !expVar = expVar1 * expVar2 * expVar3
  print *,"HyperParams, lScale=",lScale(1), lScale(2),lScale(3),lScale(4), lScale(5), &
       lScale(6),lScale(7),lScale(8),lScale(9)
  
  print *,"HyperParams",expVar,NuggVar, gpEmax
  close(7)
  
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
     read (7,*) xTraining(1,i), xTraining(2,i), xTraining(3,i), xTraining(4,i), xTraining(5,i),&
          xTraining(6,i), xTraining(7,i), xTraining(8,i), xTraining(9,i)
     !print *,xTraining(1,i), xTraining(2,i), xTraining(3,i), xTraining(4,i), xTraining(5,i), &
     !     xTraining(6,i), xTraining(7,i), xTraining(8,i), xTraining(9,i)
  end do
  close(7)

  write (filename, '( "2CO2.sym" )' )
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
