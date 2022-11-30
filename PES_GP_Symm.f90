!!CO2-Ar PES
!!Distances in Angstrom and energies in Hartrees
module PES_details
  !!  double precision :: gpRMax = 8.5
  double precision :: gpRMax = 1000.0
  double precision :: gpRMin = 2.5  
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
  integer :: nPerms=6
  !character (len=20) :: label=''
  !integer :: nTraining=337
  
  character (len=20) :: label='lq_'
  integer :: nTraining=999
  !integer :: nTraining=350

end module GP_variables


! Test program
use PES_details
implicit none
 
integer k,i, choice

double precision e, PES
double precision, allocatable :: rab(:)
double precision xStar(3)

!!allocate( xStar(nDim))


xStar(1:3)=(/0.3595111400000000068111206,0.2860995399999999855289445, &
     0.2801690199999999908442305/)

xStar(1:3)=(/0.3436426100000000150025414,0.1184221899999999966235364, &
     0.0880718719999999954950454/)

xStar(1:3)=(/0.3436426100000000150025414,0.1244827099999999964197173, &
     0.0913805810000000023896050/)

!xStar(1:3)=(/0.3331502100000000021751134, 0.3260882999999999976026288, &
!     0.2039624999999999910293980/)

!xStar(1:3)=(/0.3992589799999999855550925, 0.3779625200000000240230236, &
!     0.2954577199999999792545680/)

!xStar(1:3)=(/0.3004692999999999947213780,0.2998323600000000199727879, &
!     0.1819913700000000134071598/)

!xStar(1:3)=(/3.4364261e-01 ,  1.1842219e-01 ,  8.8071872e-02/)

!xStar(1:3)=(/3.4364261e-01 ,  2.3582520e-01,   1.3985175e-01/)

!xStar(1:3)=(/3.4364261e-01  , 2.8349597e-01 ,  1.5534253e-01/)

!xStar(1:3)=(/0.2765329200000000153814028 ,0.2449013399999999951450746, &
!     0.1726287499999999974775733/)

xStar(1:3)=(/ 0.28306178,   0.23779215 ,  0.13116063 /)

call load_GP_Data
e=PES_GP( xStar)
write(6,*)e

!call EvsE
!call fixedAngleSlice

end
!


subroutine EvsE()
  use PES_details
  use GP_variables
  implicit none
  double precision, allocatable:: rab(:), xStar(:)
  double precision dum(9), NonAdd, funcVal, PES, RMSE, mean
  integer i,j, count, nPoints

  allocate (rab(nDim), xStar(nDim) )

  !====Read test data and compute error====
  !open (unit = 7, file = "hiQ_test1_CCSDT_rInv-4000.lhc")
  !nPoints=4081
  !open (unit = 7, file = "hiQ_scaleTh90_CCSDT_rInv-128.lhc")
  !nPoints=125
  open (unit = 7, file = "hiQ_scaleTh60_CCSDT_rInv-100.lhc")
  nPoints=99
  !open (unit = 7, file = "hiQ_scaleTh0_CCSDT_rInv-100.lhc")
  !nPoints=60
  open (unit=15, file="PES_Err.dat ", status='replace')

  
  RMSE=0
  mean=0
  count=0
  do i=1,nPoints !!10078
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
        write(15,*) 1.0/rab(1)/rab(2)/rab(3), rab(1), rab(2), rab(3), NonAdd, funcVal, Sqrt((NonAdd-funcVal)**2)
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
  double precision dum(8), NonAdd, funcVal, PES,r,e, cosTheta,r12,x
  integer i,j,itot
  allocate (rab(nDim), xStar(nDim) )
  itot=500
  
  cosTheta=0.5
  r12=5.0

  open (unit=15, file="PES_Out.dat ", status='replace')
  
  do i=0, itot

     ! specify centre-to-centre separation
     x = (  0.5 + 15.0*i/(1.0*itot) ) 

     rab(1)=r12
     rab(2)=Sqrt( r12**2/4.0  +  x**2  -  r12*x*cosTheta )
     rab(3)=Sqrt( r12**2/4.0  +  x**2  +  r12*x*cosTheta )
     
     e=PES( rab)

      !!rab(1)=r12
     !!rab(2)=Sqrt( r12**2/4.0  +  x**2  -  r12*x*cosTheta )
     !!rab(3)=Sqrt( r12**2/4.0  +  x**2  +  r12*x*cosTheta )
     !e_GP = PES_GP( xStar)
     write(15,*) rab(3) , e 
     
  enddo
  write(6,*)'Written to file: PES_Out.dat '
  close(15)

end subroutine fixedAngleSlice
  
  
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Gaussian Process Code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pom


  
subroutine load_GP_Data
  use GP_variables
  use PES_details
  implicit none
  
  double precision, allocatable::  xStar(:)
  integer i,j,k
  double precision :: dum, expVar1, expVar2,expVar3
  character (len=99) :: filename, subname
  integer, allocatable:: perm(:,:)

  allocate (alpha(nTraining), lScale(nDim), xTraining(nDim,nTraining),xTrainingPerm(nDim,nPerms,nTraining), xStar(nDim), &
       perm(nDim,nPerms))

  !====Load hyperparameters====
  subname = "TrainingData/"//trim(adjustl(label))//"HyperParams_Symm"
  write (filename,'(A90,I3.3,".dat")') trim(adjustl(subname)), nTraining
  open (unit = 7, file = adjustl(filename),Status='old')
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
  subname = "TrainingData/"//trim(adjustl(label))//"alpha_Symm"
  write (filename,'(A90,I3.3,".dat")') trim(adjustl(subname)), nTraining
  open (unit = 7, file = adjustl(filename),Status='old')
  do i=1,nTraining
     read (7,*) alpha(i)
     !!print *,"alpha ",i, alpha(i)
  end do
  close(7)


  !====Load training data x values ====
  subname = "TrainingData/"//trim(adjustl(label))//"xTraining"
  write (filename,'(A90,I3.3,".dat")') trim(adjustl(subname)), nTraining
  print *,adjustl(filename)
  open (unit = 7, file = adjustl(filename),Status='old')
    
  do i=1,nTraining
     read (7,*) xTraining(1,i), xTraining(2,i), xTraining(3,i)
     !print *,xTraining(1,i), xTraining(2,i), xTraining(3,i), xTraining(4,i), xTraining(5,i), &
     !     xTraining(6,i), xTraining(7,i), xTraining(8,i), xTraining(9,i)
  end do
  close(7)

  write (filename, '( "3Ar.sym" )' )
  open (unit = 7, file = filename,Status='old')
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
  double precision rab(3), xStar(3)
  double precision  PES
  double precision repFactor, oldr3, r13MIN,  r13new, r23new, r12new,temp, scale
  double precision cosTheta !![angle between Ar3 and centre of 1,2]
  double precision xnew !! x [distance between 3 and centre of 1,2]
  integer FLAG

  FLAG=0
  r13MIN=5.5
  repFactor=1.0


  r12new=rab(1)
  r13new=rab(2)
  r23new=rab(3)

  !!! Sort vectors
  if (r12new > r13new) then 
     temp=r12new
     r12new=r13new
     r13new=temp
  endif
  if (r13new > r23new) then
     temp=r13new
     r13new=r23new
     r23new=temp
  endif
  if (r12new > r13new) then 
     temp=r12new
     r12new=r13new
     r13new=temp
  endif

  !write(6,*) r12new, r13new,r23new, rab(1), rab(2), rab(3)

     
  if( r12new > r13MIN) then
     !!Scale down all distances so that r12=r13Min (5.5A)
     FLAG=1
     scale=r13Min/r12new
     !r12new=r13MIN
     r12new=r12new*scale
     r13new=r13new*scale
     r23new=r23new*scale
     !write(6,*)'Should be 5.5',r12new
  endif !!r12new > r13MIN
  
  if ( r13new > r13MIN ) then
     FLAG=1    
     !!====Slide along x until r13=r13MIN
     !! compute cosTheta (this is fixed throughout)
     cosTheta= (r23new**2-r13new**2)/r12new/Sqrt(2.0*(r13new**2+r23new**2)-r12new**2)
     !!cosTheta= (rab(3)**2-rab(2)**2)/rab(1)/Sqrt(2.0*(rab(2)**2+rab(3)**2)-rab(1)**2)

     !! change x so that r13=r13MIN
     xnew= 0.5 *( cosTheta*r12new + Sqrt( r12new**2*(-1.0 + cosTheta**2)  +  4.0 * r13MIN**2 ) )
     !xnew= 0.5 *( cosTheta*rab(1) + Sqrt( rab(1)**2*(-1.0 + cosTheta**2)  +  4.0 * r13MIN**2 ) )
     
     r13new = Sqrt(r12new**2/4.0 + xnew*xnew - r12new*xnew*cosTheta)
     r23new = Sqrt(r12new**2/4.0 + xnew*xnew + r12new*xnew*cosTheta)
     !r13new = Sqrt(rab(1)**2/4.0 + xnew*xnew - rab(1)*xnew*cosTheta)
     !r23new = Sqrt(rab(1)**2/4.0 + xnew*xnew + rab(1)*xnew*cosTheta)

     
     
  endif !!rab(2) > r13MIN

  if( r23new >  gpRMax) then
     FLAG=1
     
     !!Slide along x until r23=gpRMax
     !cosTheta= (r23new**2-r13new**2)/r12new/Sqrt(2.0*(r13new**2+r23new**2)-r12new**2)
     !xnew= 0.5 *( -cosTheta*r12new + Sqrt( r12new**2*(-1.0 + cosTheta**2)  +  4.0 * gpRMax**2 ) )
     !r13new = Sqrt(r12new**2/4.0 + xnew*xnew - r12new*xnew*cosTheta)
     !r23new = Sqrt(r12new**2/4.0 + xnew*xnew + r12new*xnew*cosTheta)

     !!Scale down all distances so that r23=gpRMax (8.5A)
     scale=gpRMax/r23new
     r12new=r12new*scale
     r13new=r13new*scale
     r23new=r23new*scale
     !!write(6,*) r12new,r13new, r23new
     !!write(6,*)'Should be 8.5',r23new
     
  endif !!rab(3)>gpRMax
  
  !!If flag has been set then compute with GP using scaled distances then use power law to extrapolate to required distances for rab
  if(FLAG>0) then
     xStar(1) = 1/r12new
     xStar(2)=1/r13new
     xStar(3)=1/r23new
     PES = PES_GP( xStar) * ( r12new/rab(1) * r13new/rab(2) * r23new/rab(3) )**3
     return
  endif

     
  if ( ANY( rab < gpRMin/repFactor)  ) then 
     PES=0.0 !! Use repulsive approximation function
     return
  endif
  
  !! Use the Guassian Process function
  xStar(:) = 1/rab(:)
  PES = PES_GP( xStar)
  

  !PES=gpEmax/3.0* (1.0/rab(1)**12+1.0/rab(2)**12+1.0/rab(3)**12) *gpRMin **12

  
end function PES
