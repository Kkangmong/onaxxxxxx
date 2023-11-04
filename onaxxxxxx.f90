!=========== here is the program======================
!=====================================================

Program Wave_2D_onadm
Implicit none
Include 'mpif.h'
!=====================================
!定义一个数据类型，专门读取segy格式的数据
!========================================================
	!/* TYPEDEFS */
	type segy 	!/* segy - trace identification header */
	sequence
	integer tracl 	!/* Trace sequence number within line		
	integer tracr 	!/* Trace sequence number within SEG Y file
	integer fldr 	!/* Original field record number		
	integer tracf 	!/* Trace number within original field record 	
	integer ep 		!/* energy source pointeger number
	integer cdp 	!/* Ensemble number (i.e. CDP, CMP, CRP,...) 		
	integer cdpt 	!/* trace number within the ensemble
		
	integer (kind=2) trid 	!/* trace identification code:		
	integer (kind=2) nvs 	!/* Number of vertically summed traces yielding		
	integer (kind=2) nhs 	!/* Number of horizontally summed traces yielding		
	integer (kind=2) duse 	!/* Data use:
	integer offset 	!/* Distance from the center of the source pointeger 
	integer gelev 	!/* Receiver group elevation from sea level
	integer selev 	!/* Surface elevation at source.
	integer sdepth 	!/* Source depth below surface (a positive number).
	integer gdel 	!/* Datum elevation at receiver group.
	integer sdel 	!/* Datum elevation at source.
	integer swdep 	!/* Water depth at source.
	integer gwdep 	!/* Water depth at receiver group.
	integer (kind=2) scalel 	!/* Scalar to be applied to the previous 7 entries
	integer (kind=2) scalco 	!/* Scalar to be applied to the next 4 entries
	integer  sx 	!/* Source coordinate - X 
	integer  sy 	!/* Source coordinate - Y 
	integer  gx 	!/* Group coordinate - X 
	integer  gy 	!/* Group coordinate - Y 

	 integer (kind=2) counit 	!/* Coordinate units: (for previous 4 entries and
	 integer (kind=2) wevel 	!/* Weathering velocity. 
	 integer (kind=2) swevel 	!/* Subweathering velocity. 
	 integer (kind=2) sut 	!/* Uphole time at source in milliseconds. 	
	 integer (kind=2) gut 	!/* Uphole time at receiver group in milliseconds. 
	 integer (kind=2) sstat 	!/* Source static correction in milliseconds. 	
	 integer (kind=2) gstat 	!/* Group static correction  in milliseconds.
	 integer (kind=2) tstat 	!/* Total static applied  in milliseconds.
	 integer (kind=2) laga 	!/* Lag time A, time in ms between end of 240-
	 integer (kind=2) lagb 	!/* lag time B, time in ms between the time break
	 integer (kind=2) delrt 	!/* delay recording time, time in ms between
	 integer (kind=2) muts 	!/* mute time--start 
	 integer (kind=2) mute 	!/* mute time--end 
     integer (kind=2) ns 	!/* number of samples in this trace 
	 integer (kind=2) dt 	!/* sample interval  in micro-seconds
	 integer (kind=2) gain 	!/* gain type of field instruments code:
     integer (kind=2) igc 	!/* instrument gain constant 
	 integer (kind=2) igi 	!/* instrument early or initial gain 
	 integer (kind=2) corr 	!/* correlated:
	 integer (kind=2) sfs 	!/* sweep frequency at start 
	 integer (kind=2) sfe 	!/* sweep frequency at end
	 integer (kind=2) slen 	!/* sweep length in ms 
	 integer (kind=2) styp 	!/* sweep type code:
	 integer (kind=2) stas 	!/* sweep trace length at start in ms
	 integer (kind=2) stae 	!/* sweep trace length at end in ms 
	 integer (kind=2) tatyp 	!/* taper type: 1=linear, 2=cos^2, 3=other 
	 integer (kind=2) afilf 	!/* alias filter frequency if used 
	 integer (kind=2) afils 	!/* alias filter slope
	 integer (kind=2) nofilf 	!/* notch filter frequency if used
	 integer (kind=2) nofils 	!/* notch filter slope
	 integer (kind=2) lcf 	!/* low cut frequency if used
	 integer (kind=2) hcf 	!/* high cut frequncy if used
	 integer (kind=2) lcs 	!/* low cut slope
	 integer (kind=2) hcs 	!/* high cut slope
     integer (kind=2) year 	!/* year data recorded
	 integer (kind=2) day 	!/* day of year
	 integer (kind=2) hour 	!/* hour of day (24 hour clock) 
	 integer (kind=2) minute 	!/* minute of hour
	 integer (kind=2) sec 	!/* second of minute
	 integer (kind=2) timbas 	!/* time basis code:
	 integer (kind=2) trwf 	!/* trace weighting factor, defined as 1/2^N
	 integer (kind=2) grnors 	!/* geophone group number of roll switch
	 integer (kind=2) grnofr 	!/* geophone group number of trace one within
	 integer (kind=2) grnlof 	!/* geophone group number of last trace within
	 integer (kind=2) gaps 	!/* gap size (total number of groups dropped)
	 integer (kind=2) otrav 	!/* overtravel taper code:

	real d1 	!/* sample spacing for non-seismic data
	real f1 	!/* first sample location for non-seismic data
	real d2 	!/* sample spacing between traces
	real f2 	!/* first trace location
	real ungpow 	!/* negative of power used for dynamic
	real unscale 	!/* reciprocal of scaling factor to normalize
	integer (kind=2) mark 	!/* mark selected traces
	integer (kind=2) mutb 	!/* mute time at bottom (start time)
	real dz 	!/* depth sampling interval in (m or ft)
	real fz 	!/* depth of first sample in (m or ft)
    integer (kind=2) n2 	!/* number of traces per cdp or per shot
    integer (kind=2) shortpad  !/* alignment padding
	integer ntr  	!/* number of traces
	integer (kind=2) unass(8) 	!/* unassigned
	real  trdata(1500) 
	end type segy 
!=================================================================================
!=================================================================================
!use segyhdr
type(segy):: tr(123540)
Integer NI_Global,NT_Global,NK_Global,NPz,NK_local,Kerrel,IOS,Record_length,ratio,M,M1
Real TT,ff0,Pi,dt,h1,h2,h3,h,z,temp1,temp2,cc,ftm,a
INTEGER counter,error,ntaper,ii,leftX
!速度横向为Ntrace，现在往前拓展20，变为1087
!这里为了第一炮能够跟其他炮规则相同，其实左边界的延拓是不符合事实的。这里要注意修改！！！

!Parameter(Ratio=1,M=5000,Ntrace=1067,ftm=0.3048,M1=M)
!Parameter(Ratio=1,M=5000,Ntrace=1107,ftm=0.3048,M1=M)
!Parameter(NI_Global=360*ratio,NT_Global=M,Nk_Global=402*ratio)
!原来是12+348，现在左边界往前拓展20层，顾变为360+20=380
!Parameter(NI_Global=410*ratio,NT_Global=M,Nk_Global=500*ratio)
!Parameter(NI_Global=360*ratio,NT_Global=M,Nk_Global=500*ratio)
Parameter(Ratio=1,M=5500,ftm=0.3048,M1=M)
Parameter(NI_Global=1067*ratio,NT_Global=M,Nk_Global=500*ratio)

!========================= ======================
!Parameter(NPz=20,ntaper=10) 
Parameter(NPz=100,ntaper=10)                                            
Parameter(NK_local=NK_Global/NPz)                          
Parameter(TT=1.1,ff0=20,Pi=3.1415926,Kerrel=4,leftX=0)
Parameter(dt=0.002/ratio,h1=0.075/ratio*ftm,h2=h1,h3=0.075/ratio*ftm,h=h1,z=h3,Record_length=TT/dt+1)   
!Mpi variable
INTEGER NProc                                  
INTEGER MyRank,MyTop,Myfloor                            
integer zglobal              
INTEGER Kstart,Kend,kstartf     
INTEGER Xtype,Xztype,XTtype
INTEGER status(Mpi_status_size),Ierr

Real ttime,min,max
Real ffs,ffc,fft,ff2t


INTEGER(4) startime,endtime,totletime,vibrate,transbegin,transend,transtime
!epicenter located (Fx,Fy,Fz);loop indics I,J,K
INTEGER FX,FZ,I,t,K,Norm,J,l

!array variable
DOUBLE PRECISION U(NI_Global,0:NK_local+1,NT_Global)
DOUBLE PRECISION Ux(NI_Global,0:NK_local+1,NT_Global)
DOUBLE PRECISION Uz(NI_Global,0:NK_local+1,NT_Global)
DOUBLE PRECISION w(NI_Global,NK_Global),aa(ntaper)
!DOUBLE PRECISION ps(NI_Global,0:NK_local+1,NT_Global)
!============反方向============
DOUBLE PRECISION U1(NI_Global,0:NK_local+1),U2(NI_Global,0:NK_local+1),U0(NI_Global,0:NK_local+1)
DOUBLE PRECISION U1x(NI_Global,0:NK_local+1),U2_x(NI_Global,0:NK_local+1),U0x(NI_Global,0:NK_local+1)
DOUBLE PRECISION U1z(NI_Global,0:NK_local+1),U2_z(NI_Global,0:NK_local+1),U0z(NI_Global,0:NK_local+1)

!save the 3 snapshot in xy-,xz-,yz-plants.
DOUBLE PRECISION c(NI_Global,NK_Global)
				  !,&&veltotal(Ntrace,NK_Global),img(Ntrace,NK_Global)
!				  ,&&dataa(NI_Global-6,M/2+1),XZU(NI_Global,NK_Global)
DOUBLE PRECISION dataa(348,1500)

!higher order derivative,scale
DOUBLE PRECISION U2x,U2z,Uxz,U3x,U3z,U2xz,Ux2z
DOUBLE PRECISION U4x,U4z,U3xz,Ux3z,U2x2z				
DOUBLE PRECISION U5x,U5z,U4xz,Ux4z,U3x2z,U2x3z
DATA U2x,U2z,Uxz,U3x,U3z,U2xz,Ux2z/7*0.0/
DATA U4x,U4z,U3xz,Ux3z,U2x2z/5*0.0/				

DOUBLE PRECISION P,Px,Pz,P2x,P2z,Pxz,P3x,P3z,P2xz,Px2z,P2t,Px2t,Pz2t !,Pt,Pxt,Pyt,Pzt
DATA P,Px,Pz,P2x,P2z,Pxz,P3x,P3z,P2xz,Px2z,P2t,Px2t,Pz2t/13*0.0/
!the receiver
DOUBLE PRECISION receiver1(Record_length)
Integer receiver1_processor

CALL Mpi_Init(Ierr)
CALL Mpi_Comm_size(Mpi_Comm_World,NProc,Ierr)
IF(Nproc.Ne.Npz) then
  Print *, 'error and mpi exit'
  CALL Mpi_Finalize(Ierr)
  Stop
EndIF

CALL Mpi_Comm_Rank(Mpi_Comm_World,MyRank,Ierr)

Kstart=0
kstartf=0
Kend=NK_local+1
MyTop=MyRank-1
Myfloor=MyRank+1

IF(Myrank.EQ.0) then 
  Mytop=Mpi_Proc_Null 
  Kstart=1
  kstartf=1
EndIF 

IF(Myrank.EQ.(NPz-1)) then 
  Myfloor=Mpi_Proc_Null
  Kend=NK_local 
EndIF 

!define the derivative dateype
CALL Mpi_Type_contiguous(NI_Global*1,MPI_DOUBLE_PRECISION,Xtype,Ierr)
CALL Mpi_Type_Commit(Xtype,Ierr)

Call MPI_Type_contiguous(NI_Global*NK_local,MPI_DOUBLE_PRECISION,xZtype,Ierr)
Call MPI_Type_Commit(XZtype,Ierr)

Call MPI_Type_contiguous(348*1500,MPI_DOUBLE_PRECISION,XTtype,Ierr)
Call MPI_Type_Commit(XTtype,Ierr)

!
vibrate=1
T=1
! ========read velocity==================
!========================================================================
!110 format(I4)
110 format(F8.2)
if(myrank.EQ.0)then
OPEN(FILE='str_zfrom2_1067_402.txt',UNIT=10)
  DO I=LeftX+1,NI_Global
      DO J=59,NK_Global-40
       READ(10,110,advance='yes')a
	 	  c(i,j)=a/1000*ftm
!         veltotal(i,j)=a
!		 img(i,j)=0
	  END DO
   END DO
CLOSE(10)

! DO I=1,leftX
 !     DO J=59,NK_Global-40      
!	 	 vel(i,j)=vel((LeftX+1),j)
!		 img(i,j)=0
!	  END DO
 ! END DO
DO I=1,NI_Global      
         c(i,1:58)=c(i,59)
!		 img(i,1:58)=0

		 c(i,NK_Global-59:NK_Global)=c(i,NK_Global-59)
!		 img(i,NK_Global-59:NK_Global)=0
 END DO
end if


!============================

!=========l循环开始===============
!do l=1,350
!do l=60,350
!实际上shot从i=13,15,17....,detx=75ft
!if((l.ne.543).or.(l.ne.545)) then
!===================================================
!=======give velocity================

 If (Myrank.EQ.0) Then
   Do I=1,NPz-1
     CALL MPi_Send(c(1,i*NK_local+1),1,xztype,i,vibrate+100,Mpi_Comm_World,Ierr)
   end do
  Else
        CALL Mpi_Recv(c(1,1),xztype,MPI_DOUBLE_PRECISION,0,&
             &vibrate+100,Mpi_Comm_World,Status,Ierr)
 EndIf
!===============================================================
!open(unit=20,file='sxnew.su',form='unformatted',&
 !     &access='direct',recl=1560*4)
!===============================================================

!==============================================================
!DO k=1,NK_local
!   DO I=1,NI_Global
!		 c(I,K)=4.
!	     zGlobal=k+Myrank*NK_local
!	  	 if(zGlobal>180) then
!		 c(I,K)=5.
!		 end if
!   endDO
!endDO

DO k=1,NK_local
   DO I=1,NI_Global
		 w(I,K)=0.
   endDO
endDO

!FX=NI_Global/2

do k=1,ntaper
aa(k)=0.5-0.5*cos(pi*(k-1)/ntaper)
enddo

transtime=0
startime=Mpi_Wtime()
!l就是炮数的变化，从1,700,2
do l=11,707,2
!l=NI_Global/2+50
!FX=7

!problem
Fx=12+l
FZ=60 
!Initial data for U and Ut.
T=1
DO k=Kstart,Kend
   DO I=1,NI_Global
         U(I,k,t)=0.
         Ux(I,k,t)=0.
         Uz(I,k,t)=0.
		  U(I,k,t+1)=0.
         Ux(I,k,t+1)=0.
         Uz(I,k,t+1)=0.
   endDO
endDO

 ttime=0.0
 vibrate=1
! receiver1_processor=(Fx+NI_Global/8)/NK_local
 !write(*,*) receiver1_processor
!==============================begin the recurrence===============================
DO T=2,M1-1
 ttime=ttime+dt
  !Vibrate=Vibrate+1
   ffs=sin(2.0*pi*ff0*ttime)*exp(-pi*pi*ff0*ff0*ttime*ttime/4.0)
   ffc=cos(2.0*pi*ff0*ttime)*exp(-pi*pi*ff0*ff0*ttime*ttime/4.0)
   fft=2.0*pi*ff0*ffc-pi**2*ff0**2*ttime*ffs/2.0 
   ff2t=(-9.0/2.0+pi**2*ff0**2*ttime**2/4.0)*pi**2*ff0**2*ffs-2.0*pi**3*ff0**3*ttime*ffc

!========================the true code of this program============================
!  ffs=-5.76*ff0**2*(1-16*(0.6*ff0*ttime-1)**2)*exp(-8*(0.6*ff0*ttime-1)**2)
!  fft=-3.456*ff0**3*(256*(0.6*ff0*ttime-1)**3-48*(0.6*ff0*ttime-1))&
!     &*exp(-8*(0.6*ff0*ttime-1)**2)
!  ff2t=-2.0736*ff0**4*(768*(0.6*ff0*ttime-1)**2-16*(0.6*ff0*ttime-1)-48)&
!     &*exp(-8*(0.6*ff0*ttime-1)**2)
!the inner grid

 DO k=Kstart+1,Kend-1
    DO I=2,NI_Global-1
            U2x=(2/(h1*h1))*(U(i+1,k,t)-2*U(i,k,t)+U(i-1,k,t))-(Ux(i+1,k,t)-Ux(i-1,k,t))/(2*h1) 
			U2z=(2/(h3*h3))*(U(i,k+1,t)-2*U(i,k,t)+U(i,k-1,t))-(Uz(i,k+1,t)-Uz(i,k-1,t))/(2*h3)   

            U3x=(15/(2*h1**3))*(U(i+1,k,t)-U(i-1,k,t))-(3/(2*h1**2))*(Ux(i+1,k,t)+8*Ux(i,k,t)+Ux(i-1,k,t)) 
 			U3z=(15/(2*h3**3))*(U(i,k+1,t)-U(i,k-1,t))-(3/(2*h3**2))*(Uz(i,k+1,t)+8*Uz(i,k,t)+Uz(i,k-1,t)) 
			
			U2xz=(1/(2*h1*h3))*(-Ux(i+1,k+1,t)-Ux(i-1,k-1,t)+Ux(i+1,k,t)+Ux(i-1,k,t)-2*Ux(i,k+1,t)+4*Ux(i,k,t)&
			&-2*Ux(i,k-1,t))+(1/(h1**2))*(Uz(i+1,k,t)-2*Uz(i,k,t)+Uz(i-1,k,t))+(1/(4*h1**2*h3))*(5*U(i+1,k+1,t)&
			&-5*U(i-1,k-1,t)+U(i+1,k-1,t)-U(i-1,k+1,t)-6*U(i+1,k,t)+6*U(i-1,k,t)-4*U(i,k+1,t)+4*U(i,k-1,t)) 
			Ux2z=(1/(2*h1*h3))*(-Uz(i+1,k+1,t)-Uz(i-1,k-1,t)+Uz(i,k+1,t)+Uz(i,k-1,t)-2*Uz(i+1,k,t)+4*Uz(i,k,t)&
			&-2*Uz(i-1,k,t))+(1/(h3**2))*(Ux(i,k+1,t)-2*Ux(i,k,t)+Ux(i,k-1,t))+(1/(4*h1*h3**2))*(5*U(i+1,k+1,t)&
			&-5*U(i-1,k-1,t)+U(i-1,k+1,t)-U(i+1,k-1,t)-6*U(i,k+1,t)+6*U(i,k-1,t)-4*U(i+1,k,t)+4*U(i-1,k,t)) 

            U4x=(6/(h1**3))*(Ux(i+1,k,t)-Ux(i-1,k,t))-(12/(h1**4))*(U(i+1,k,t)-2*U(i,k,t)+U(i-1,k,t)) 
            U4z=(6/(h3**3))*(Uz(i,k+1,t)-Uz(i,k-1,t))-(12/(h3**4))*(U(i,k+1,t)-2*U(i,k,t)+U(i,k-1,t)) 
 !second         
            U2x2z=(U(i+1,k+1,t)+U(i-1,k-1,t)+U(i+1,k-1,t)+U(i-1,k+1,t)+4*U(i,k,t)-2*U(i+1,k,t)-2*U(i-1,k,t)&
			&-2*U(i,k+1,t)-2*U(i,k-1,t))/(h1**2*h3**2)
 
            U5x=-90/(h1**5)*(U(i+1,k,t)-U(i-1,k,t))+30/(h1**4)*(Ux(i+1,k,t)+4*Ux(i,k,t)+Ux(i-1,k,t)) 
            U5z=-90/(h3**5)*(U(i,k+1,t)-U(i,k-1,t))+30/(h3**4)*(Uz(i,k+1,t)+4*Uz(i,k,t)+Uz(i,k-1,t)) 
              
            U4xz=6/(h1**3*h3)*(Ux(i+1,k+1,t)+Ux(i-1,k-1,t)+2*Ux(i,k+1,t)-4*Ux(i,k,t)+2*Ux(i,k-1,t)-&
			&Ux(i+1,k,t)-Ux(i-1,k,t))-3/(h1**4*h3)*(5*U(i+1,k+1,t)-5*U(i-1,k-1,t)+U(i+1,k-1,t)-&
			&U(i-1,k+1,t)-6*U(i+1,k,t)+6*U(i-1,k,t)-4*U(i,k+1,t)+4*U(i,k-1,t))   
			         
            Ux4z=6/(h1*h3**3)*(Uz(i+1,k+1,t)+Uz(i-1,k-1,t)+2*Uz(i+1,k,t)-4*Uz(i,k,t)+2*Uz(i-1,k,t)&
			&-Uz(i,k+1,t)-Uz(i,k-1,t))-3/(h1*h3**4)*(5*U(i+1,k+1,t)-5*U(i-1,k-1,t)+U(i-1,k+1,t)-&
			&U(i+1,k-1,t)-6*U(i,k+1,t)+6*U(i,k-1,t)-4*U(i+1,k,t)+4*U(i-1,k,t)) 
                    
            U3x2z=3/(h1**3*h3**2)*(U(i+1,k+1,t)-U(i-1,k-1,t)+U(i+1,k-1,t)-U(i-1,k+1,t)-2*U(i+1,k,t)&
			&+2*U(i-1,k,t))-6/(h1**2*h3**2)*(Ux(i,k+1,t)-2*Ux(i,k,t)+Ux(i,k-1,t)) 
            U2x3z=3/(h1**2*h3**3)*(U(i+1,k+1,t)-U(i-1,k-1,t)+U(i-1,k+1,t)-U(i+1,k-1,t)-2*U(i,k+1,t)&
			&+2*U(i,k-1,t))-6/(h1**2*h3**2)*(Uz(i+1,k,t)-2*Uz(i,k,t)+Uz(i-1,k,t))           
            		
! P及其各阶偏导数的计算（公式？）	
          cc=c(i,k)*c(i,k)   
        temp1=cc*cc	  
	    P=(U2x+U2z)*cc
        P2t=(U4x+U4z+2*U2x2z)*temp1  
        Px=(U3x+Ux2z)*cc

		Pz=(U2xz+U3z)*cc

		Px2t=(U5x+Ux4z+2*U3x2z)*temp1    
        Pz2t=(U4xz+U5z+2*U2x3z)*temp1    
       
!        ps(i,k,t)=p

	   zGlobal=k+Myrank*NK_local

	   IF(((i.EQ.(FX-2)).OR.(i.EQ.(FX+2))).AND.(zGlobal.GE.(Fz-2)).AND.(zGlobal.LE.(Fz+2))) THEN
		   P=P+ffs/(4)
		   P2t=P2t+ff2t/(4)
       END IF

		IF((i.GE.(FX-1)).AND.(i.LE.(FX+1)).AND.((zGlobal.EQ.(Fz-2)).OR.(zGlobal.EQ.(Fz+2)))) THEN
		   P=P+ffs/(4)
		   P2t=P2t+ff2t/(4)
        END IF

        IF((i.GE.(FX-1)).AND.(i.LE.(FX+1)).AND.(zGlobal.GE.(Fz-1)).AND.(zGlobal.LE.(Fz+1))) THEN
			IF((i.EQ.FX).AND.(zGlobal.EQ.Fz)) THEN
				P=P+ffs
				P2t=P2t+ff2t
		ELSE
				P=P+ffs/(2.0)
				P2t=P2t+ff2t/(2.0)
			END IF
    	end if 	
		temp2=dt*dt
		
		U(i,k,t+1)=2*U(i,k,t)-U(i,k,t-1)+temp2*(P+temp2*P2t/12.0)
       Ux(i,k,t+1)=2*Ux(i,k,t)-Ux(i,k,t-1)+temp2*(Px+temp2*Px2t/12.0)
       Uz(i,k,t+1)=2*Uz(i,k,t)-Uz(i,k,t-1)+temp2*(Pz+temp2*Pz2t/12.0)
	   
	 END DO
 END DO

!==========deal with the inner-boundary================================================
!=======top================= 
if(myrank.eq.0) then

DO I=1,NI_Global  !       		 
   U(I,1,t+1)=(2*U(I,1,t)-U(I,1,t-1)-U(I,3,t+1)+2*U(I,3,t)-U(I,3,t-1)+&
   &c(i,1)*dt*(U(I,3,t+1)-U(I,3,t-1)+U(I,1,t-1))/h3-2*(U(I,3,t)-&
   &2*U(I,2,t)+U(I,1,t))*((c(i,1)*dt/h3)**2))/(1+c(i,1)*dt/h3)
         		 
   Ux(I,1,t+1)=(2*Ux(I,1,t)-Ux(I,1,t-1)-Ux(I,3,t+1)+2*Ux(I,3,t)-Ux(I,3,t-1)+&
   &c(i,1)*dt*(Ux(I,3,t+1)-Ux(I,3,t-1)+Ux(I,1,t-1))/h3-2*(Ux(I,3,t)-&
   &2*Ux(I,2,t)+Ux(I,1,t))*((c(i,1)*dt/h3)**2))/(1+c(i,1)*dt/h3)
        		 
   Uz(I,1,t+1)=(2*Uz(I,1,t)-Uz(I,1,t-1)-Uz(I,3,t+1)+2*Uz(I,3,t)-Uz(I,3,t-1)+&
   &c(i,1)*dt*(Uz(I,3,t+1)-Uz(I,3,t-1)+Uz(I,1,t-1))/h3-2*(Uz(I,3,t)-&
   &2*Uz(I,2,t)+Uz(I,1,t))*((c(i,1)*dt/h3)**2))/(1+c(i,1)*dt/h3)  
END DO

endif
!
 DO k=1,NK_local      
   U(1,k,t+1)=(2*U(1,k,t)-U(1,k,t-1)-U(1+2,k,t+1)+2*U(1+2,k,t)-U(1+2,k,t-1)&
   &+c(1,k)*dt*(U(1+2,k,t+1)-U(1+2,k,t-1)+U(1,k,t-1))/h1-2*(U(1+2,k,t)-2*U&
   &(1+1,k,t)+U(1,k,t))*((c(1,k)*dt/h1)**2))/(1+c(1,k)*dt/h1)  
       
   Ux(1,k,t+1)=(2*Ux(1,k,t)-Ux(1,k,t-1)-Ux(1+2,k,t+1)+2*Ux(1+2,k,t)-Ux(1+2,k,t-1)&
   &+c(1,k)*dt*(Ux(1+2,k,t+1)-Ux(1+2,k,t-1)+Ux(1,k,t-1))/h1-2*(Ux(1+2,k,t)-2&
   &*Ux(1+1,k,t)+Ux(1,k,t))*((c(1,k)*dt/h1)**2))/(1+c(1,k)*dt/h1)  
         
   Uz(1,k,t+1)=(2*Uz(1,k,t)-Uz(1,k,t-1)-Uz(1+2,k,t+1)+2*Uz(1+2,k,t)-Uz(1+2,k,t-1)&
   &+c(1,k)*dt*(Uz(1+2,k,t+1)-Uz(1+2,k,t-1)+Uz(1,k,t-1))/h1-2*(Uz(1+2,k,t)-2&
   &*Uz(1+1,k,t)+Uz(1,k,t))*((c(1,k)*dt/h1)**2))/(1+c(1,k)*dt/h1)      
 END DO 


 DO k=1,NK_local      
   U(NI_Global,k,t+1)=(2*U(NI_Global,k,t)-U(NI_Global,k,t-1)-U(NI_Global-2,k,t+1)+&
   &2*U(NI_Global-2,k,t)-U(NI_Global-2,k,t-1)+c(NI_Global,k)*dt*(U(NI_Global-2&
   &,k,t+1)+U(NI_Global,k,t-1)-U(NI_Global-2,k,t-1))/h1-2*(U(NI_Global,k,t)-2&
   &*U(NI_Global-1,k,t)+U(NI_Global-2,k,t))&
   &*((c(NI_Global,k)*dt/h1)**2))/(1+c(NI_Global,k)*dt/h1)

   Ux(NI_Global,k,t+1)=(2*Ux(NI_Global,k,t)-Ux(NI_Global,k,t-1)-Ux(NI_Global&
   &-2,k,t+1)+2*Ux(NI_Global-2,k,t)-Ux(NI_Global-2,k,t-1)+c(NI_Global,k)*dt&
   &*(Ux(NI_Global-2,k,t+1)+Ux(NI_Global,k,t-1)-Ux(NI_Global-2,k,t-1))/h1-2*&
   &(Ux(NI_Global,k,t)-2*Ux(NI_Global-1,k,t)+Ux(NI_Global-2,k,t))*&
   &((c(NI_Global,k)*dt/h1)**2))/(1+c(NI_Global,k)*dt/h1)

   Uz(NI_Global,k,t+1)=(2*Uz(NI_Global,k,t)-Uz(NI_Global,k,t-1)-Uz(NI_Global&
   &-2,k,t+1)+2*Uz(NI_Global-2,k,t)-Uz(NI_Global-2,k,t-1)+c(NI_Global,k)*dt&
   &*(Uz(NI_Global-2,k,t+1)+Uz(NI_Global,k,t-1)-Uz(NI_Global-2,k,t-1))/h1-2&
   &*(Uz(NI_Global,k,t)-2*Uz(NI_Global-1,k,t)+Uz(NI_Global-2,k,t))*((c(NI_&
   &Global,k)*dt/h1)**2))/(1+c(NI_Global,k)*dt/h1)
 END DO 

IF(Myrank.EQ.(NPz-1)) then
 DO I=1,NI_Global 
 !change   		 
    U(I,NK_local,t+1)=(2*U(I,NK_local,t)-U(I,NK_local,t-1)-U(I,NK_local-2,t+1)+&
	&2*U(I,NK_local-2,t)-U(I,NK_local-2,t-1)+c(i,NK_local)*dt*(U(I,NK_local-2,t+1)&
	&+U(I,NK_local,t-1)-U(I,NK_local-2,t-1))/h3-2*(U(I,NK_local,t)-2*U(I,NK_local-1,t)&
	&+U(I,NK_local-2,t))*((c(I,NK_local)*dt/h3)**2))/(1+c(I,NK_local)*dt/h3)

   Ux(I,NK_local,t+1)=(2*Ux(I,NK_local,t)-Ux(I,NK_local,t-1)-Ux(I,NK_local-2,t+1)&
   &+2*Ux(I,NK_local-2,t)-Ux(I,NK_local-2,t-1)+c(I,NK_local)*dt*(Ux(I,NK_local-2,t+1)&
   &+Ux(I,NK_local,t-1)-Ux(I,NK_local-2,t-1))/h3-2*(Ux(I,NK_local,t)-2*Ux(I,NK_local-1,t)&
   &+Ux(I,NK_local-2,t))*((c(I,NK_local)*dt/h3)**2))/(1+c(I,NK_local)*dt/h3)

   Uz(I,NK_local,t+1)=(2*Uz(I,NK_local,t)-Uz(I,NK_local,t-1)-Uz(I,NK_local-2,t+1)&
   &+2*Uz(I,NK_local-2,t)-Uz(I,NK_local-2,t-1)+c(I,NK_local)*dt*(Uz(I,NK_local-2,t+1)&
   &+Uz(I,NK_local,t-1)-Uz(I,NK_local-2,t-1))/h3-2*(Uz(I,NK_local,t)-2*Uz(I,NK_local-1,t)&
   &+Uz(I,NK_local-2,t))*((c(I,NK_local)*dt/h3)**2))/(1+c(I,NK_local)*dt/h3)
 END DO
end if

transbegin=Mpi_Wtime()
!exchage date by using so called "odd-even method"
!==================================deal with U
 
 IF (mod(Myrank,2).EQ.0) Then
    CALL Mpi_Send(U(1,Kstart+1,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(U(1,Kend-1,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Recv(U(1,Kend,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(U(1,Kstart,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
 ELSE
    CALL Mpi_Recv(U(1,Kend,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(U(1,Kstart,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Send(U(1,Kstart+1,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(U(1,Kend-1,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
 ENDIF
!deal with Ux
 IF (mod(Myrank,2).EQ.0) Then
    CALL Mpi_Send(Ux(1,Kstart+1,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(Ux(1,Kend-1,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Recv(Ux(1,Kend,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(Ux(1,Kstart,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
 ELSE
    CALL Mpi_Recv(Ux(1,Kend,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(Ux(1,Kstart,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Send(Ux(1,Kstart+1,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(Ux(1,Kend-1,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
 ENDIF

!deal with Uz
 IF (mod(Myrank,2).EQ.0) Then
    CALL Mpi_Send(Uz(1,Kstart+1,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(Uz(1,Kend-1,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Recv(Uz(1,Kend,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(Uz(1,Kstart,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
 ELSE
    CALL Mpi_Recv(Uz(1,Kend,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(Uz(1,Kstart,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Send(Uz(1,Kstart+1,t+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(Uz(1,Kend-1,t+1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
 ENDIF


 transend=Mpi_Wtime()
 transtime=transtime+transend-transbegin

!=====deal with the nature- and artificial boundary,here is none====================

END DO
!===========反方向计算===========
if(myrank.EQ.0)then
open(unit=20,file='sxnew.su',form='unformatted',&
      &access='direct',recl=1560*4)
  do ii=1,348
    read(20,rec=ii+((l+1)/2-1)*348) tr(ii+((l+1)/2-1)*348)
     dataa(ii,1:1500)=tr(ii+((l+1)/2-1)*348)%trdata(1:1500)
 enddo
!数据文件关闭
close(20)
endif  

!注意看源在哪里就传到哪里
If((Myrank.EQ.0).or.(Myrank.EQ.(Fz/nk_local-1)))Then
!注意看源在哪里就传到哪里
  If (Myrank.EQ.0) Then
     CALL MPi_Send(dataa(1,1),1,XTtype,(Fz/nk_local-1),vibrate+100,Mpi_Comm_World,Ierr)
  Else 
        CALL Mpi_Recv(dataa(1,1),XTtype,MPI_DOUBLE_PRECISION,0,&
             &vibrate+100,Mpi_Comm_World,Status,Ierr)
  endif

EndIf


DO K=Kstart,Kend
   DO I=1,NI_Global
      U1(i,k)=0
	  U0(i,k)=0
	 U0x(i,k)=0
	 U1x(i,k)=0
	 U0z(i,k)=0
	 U1z(i,k)=0
   endDO
endDO 

  
DO T=M1-1,2,-1
!if(myrank.eq.0) then

!if(mod(t,2).eq.0) then 
!  DO I=7,NI_Global    
!     U1(i,2)=U1(i,2)+U(i,2,t)
!	  U1(i,10)=U1(i,10)+(dataa(i-6,(t+2)/2)+dataa(i-6,t/2))/2
!     U0(i,10)=U0(i,10)+dataa(i-6,(t+2)/2)	  
!   enddo
!	  else
!  DO I=7,NI_Global    
!     U1(i,10)=U1(i,10)+dataa(i-6,(t+1)/2)
!     U0(i,10)=U0(i,10)+(dataa(i-6,(t+1)/2)+dataa(i-6,(t+3)/2))/2
! 	 U0(i,2)=U0(i,2)+U(i,2,t+1)
!  endDO
!endif
!	 U0x(i,9)=Ux(i,9,t+1)
!	 U1x(i,9)=Ux(i,9,t)
!	 U0z(i,9)=Uz(i,9,t+1)
!	 U1z(i,9)=Uz(i,9,t)
!endif
  DO K=Kstartf+1,Kend-1
    DO I=2,NI_Global-1
!计算当前时刻当前点的各个差分近似
            U2x=(2/(h1*h1))*(U1(i+1,k)-2*U1(i,k)+U1(i-1,k))-(U1x(i+1,k)-U1x(i-1,k))/(2*h1) 
			U2z=(2/(h3*h3))*(U1(i,k+1)-2*U1(i,k)+U1(i,k-1))-(U1z(i,k+1)-U1z(i,k-1))/(2*h3)   

            U3x=(15/(2*h1**3))*(U1(i+1,k)-U1(i-1,k))-(3/(2*h1**2))*(U1x(i+1,k)+8*U1x(i,k)+U1x(i-1,k)) 
 			U3z=(15/(2*h3**3))*(U1(i,k+1)-U1(i,k-1))-(3/(2*h3**2))*(U1z(i,k+1)+8*U1z(i,k)+U1z(i,k-1)) 
			
			U2xz=(1/(2*h1*h3))*(-U1x(i+1,k+1)-U1x(i-1,k-1)+U1x(i+1,k)+U1x(i-1,k)-2*U1x(i,k+1)+4*U1x(i,k)&
			&-2*U1x(i,k-1))+(1/(h1**2))*(U1z(i+1,k)-2*U1z(i,k)+U1z(i-1,k))+(1/(4*h1**2*h3))*(5*U1(i+1,k+1)&
			&-5*U1(i-1,k-1)+U1(i+1,k-1)-U1(i-1,k+1)-6*U1(i+1,k)+6*U1(i-1,k)-4*U1(i,k+1)+4*U1(i,k-1)) 
			Ux2z=(1/(2*h1*h3))*(-U1z(i+1,k+1)-U1z(i-1,k-1)+U1z(i,k+1)+U1z(i,k-1)-2*U1z(i+1,k)+4*U1z(i,k)&
			&-2*U1z(i-1,k))+(1/(h3**2))*(U1x(i,k+1)-2*U1x(i,k)+U1x(i,k-1))+(1/(4*h1*h3**2))*(5*U1(i+1,k+1)&
			&-5*U1(i-1,k-1)+U1(i-1,k+1)-U1(i+1,k-1)-6*U1(i,k+1)+6*U1(i,k-1)-4*U1(i+1,k)+4*U1(i-1,k)) 

            U4x=(6/(h1**3))*(U1x(i+1,k)-U1x(i-1,k))-(12/(h1**4))*(U1(i+1,k)-2*U1(i,k)+U1(i-1,k)) 
            U4z=(6/(h3**3))*(U1z(i,k+1)-U1z(i,k-1))-(12/(h3**4))*(U1(i,k+1)-2*U1(i,k)+U1(i,k-1)) 
 !second         
           !U2x2z=(U1(i+1,k+1)+U1(i-1,k-1)+U1(i+1,k-1)+U1(i-1,k+1)+4*U1(i,k)-2*U1(i+1,k)-2*U1(i-1,k)-2*U1(i,k+1)-2*U1(i,k-1))/(h1**2*h3**2)
            U2x2z=(1/(h1*h1*h3*h3))*(2*(U1(i,k+1)+U1(i,k-1)+U1(i+1,k)-2*U1(i,k)+U1(i-1,k))-U1(i+1,k+1)&
			&-U1(i-1,k-1)-U1(i-1,k+1)-U1(i+1,k-1))+(1/(h1*h3*h3)/2)*(U1x(i+1,k+1)+U1x(i+1,k-1)&
			&-U1x(i-1,k-1)-U1x(i-1,k+1)-2*U1x(i+1,k)+2*U1x(i-1,k))          
		
            U5x=-90/(h1**5)*(U1(i+1,k)-U1(i-1,k))+30/(h1**4)*(U1x(i+1,k)+4*U1x(i,k)+U1x(i-1,k)) 
            U5z=-90/(h3**5)*(U1(i,k+1)-U1(i,k-1))+30/(h3**4)*(U1z(i,k+1)+4*U1z(i,k)+U1z(i,k-1)) 
              
           !U4xz=6/(h1**3*h3)*(U1x(i+1,k+1)+U1x(i-1,k-1)+2*U1x(i,k+1)-4*U1x(i,k)+2*U1x(i,k-1)-U1x(i+1,k)-U1x(i-1,k))-3/(h1**4*h3)*(5*U1(i+1,k+1)-5*U1(i-1,k-1)+U1(i+1,k-1)-U1(i-1,k+1)-6*U1(i+1,k)+6*U1(i-1,k)-4*U1(i,k+1)+4*U1(i,k-1))            
           U4xz=-6/(h1**4*h3)*(U1(i+1,k+1)-U1(i+1,k-1)-U1(i-1,k-1)+U1(i-1,k+1)-2*U1(i,k+1)+&
		   &2*U1(i,k-1))+3/(h1**3*h3)*(U1x(i+1,k+1)+U1x(i-1,k-1)-U1x(i+1,k-1)-U1x(i-1,k+1))


		   Ux4z=6/(h1*h3**3)*(U1z(i+1,k+1)+U1z(i-1,k-1)+2*U1z(i+1,k)-4*U1z(i,k)+2*U1z(i-1,k)-U1z(i,k+1)&
		   &-U1z(i,k-1))-3/(h1*h3**4)*(5*U1(i+1,k+1)-5*U1(i-1,k-1)+U1(i-1,k+1)-U1(i+1,k-1)&
		   &-6*U1(i,k+1)+6*U1(i,k-1)-4*U1(i+1,k)+4*U1(i-1,k)) 
           !Ux4z=-6/(h3**4*h1)*(U1(i+1,k+1)-U1(i-1,k-1)-U1(i-1,k+1)+U1(i+1,k-1)-2*U1(i+1,k)+2*U1(i-1,k))+3/(h3**3*h1)*(U1z(i+1,k+1)+U1z(i-1,k-1)-U1z(i-1,k+1)-U1z(i+1,k-1))       
            
		   U3x2z=3/(h1**3*h3**2)*(U1(i+1,k+1)-U1(i-1,k-1)+U1(i+1,k-1)-U1(i-1,k+1)-2*U1(i+1,k)+2*U1(i-1,k))&
		   &-6/(h1**2*h3**2)*(U1x(i,k+1)-2*U1x(i,k)+U1x(i,k-1)) 
           U2x3z=3/(h1**2*h3**3)*(U1(i+1,k+1)-U1(i-1,k-1)+U1(i-1,k+1)-U1(i+1,k-1)-2*U1(i,k+1)+2*U1(i,k-1))&
		   &-6/(h1**2*h3**2)*(U1z(i+1,k)-2*U1z(i,k)+U1z(i-1,k))           
                        		
		  !U3x2z=-3/(h1**3*h3**2)/2*(U1(i+1,k+1)-U1(i-1,k-1)+U1(i+1,k-1)-U1(i-1,k+1)-2*U1(i+1,k)+2*U1(i-1,k))+3/(h1**2*h3**2)/2*(U1x(i+1,k+1)+U1x(i-1,k-1)+U1x(i+1,k-1)+U1x(i-1,k+1)-2*U1x(i+1,k)-2*U1x(i-1,k)) 
          !U2x3z=-3/(h1**2*h3**3)/2*(U1(i+1,k+1)-U1(i-1,k-1)-U1(i+1,k-1)+U1(i-1,k+1)-2*U1(i,k+1)+2*U1(i,k-1))+3/(h1**2*h3**2)/2*(U1z(i+1,k+1)+U1z(i-1,k-1)+U1z(i+1,k-1)+U1z(i-1,k+1)-2*U1z(i,k+1)-2*U1z(i,k-1)) 

                       		
! P及其各阶偏导数的计算（公式？）
   
		cc=c(i,k)*c(i,k)  
        temp1=cc*cc	  
	    P=(U2x+U2z)*cc
        P2t=(U4x+U4z+2*U2x2z)*temp1  
        Px=(U3x+Ux2z)*cc

		Pz=(U2xz+U3z)*cc

		Px2t=(U5x+Ux4z+2*U3x2z)*temp1    
        Pz2t=(U4xz+U5z+2*U2x3z)*temp1    
      
     
        temp2=dt*dt
		

    	U2(i,k)=2*U1(i,k)-U0(i,k)+temp2*(P+temp2*P2t/12.0)
       u2_x(i,k)=2*U1x(i,k)-U0x(i,k)+temp2*(Px+temp2*Px2t/12.0)
       u2_z(i,k)=2*U1z(i,k)-U0z(i,k)+temp2*(Pz+temp2*Pz2t/12.0)
	   !if((k.EQ.2).and.(i.EQ.170)) then


	   zGlobal=k+Myrank*NK_local
	 ! if((zGlobal.EQ.2).and.(i.eq.250)) then
	! problem
	  if(zGlobal.EQ.FZ)then
	     if((i>(FX-1)).and.(i<(FX+348)))then

           if(mod(t,4).eq.1) then 
		   U2(i,k)=U2(i,k)+dataa(i-(FX-1),(t+3)/4) 
		   endif
	      if(mod(t,4).eq.2) then 
		    U2(i,k)=U2(i,k)+dataa(i-(FX-1),(t+7)/4)/4+3*dataa(i-(FX-1),(t+3)/4)/4
		  endif

         if(mod(t,4).eq.3) then
            U2(i,k)=dataa(i-(FX-1),(t+7)/4)/2+dataa(i-(FX-1),(t+3)/4)/2
		 end if

		 if(mod(t,4).eq.0) then
             U2(i,k)=3*dataa(i-(FX-1),(t+7)/4)/4+dataa(i-(FX-1),(t+3)/4)/4
		 end if
	               	     
!	    U2(i,k)=U2(i,k)+U(i,k,t)
!       u2_x(i,k)=u2_x(i,k)+Ux(i,k,t)
!       u2_z(i,k)=u2_z(i,k)+Uz(i,k,t)
       endif
	 end if
      W(i,k)=W(i,k)+U2(i,k)*U(i,k,t-1)
 

!===最原始成像条件(laplace 滤波后)====== 
     
!  W(i,k)=W(i,k)+ps(i,k,t)*U1(i,k)+2*cc*(ux(i,k,t)*u1x(i,k)+uz(i,k,t)*u1z(i,k))+u(i,k,t)*p		   	    
	 END DO
 END DO
!-------------------------------------------------------
!==========deal with the inner-boundary================================================

if(myrank.eq.0) then
DO I=1,NI_Global 

   u2(I,1)=(2*u1(I,1)-u0(I,1)-u2(I,3)+2*u1(I,3)-u0(I,3)+c(i,1)*dt*(u2(I,3)-u0(I,3)&
   &+u0(I,1))/h3-2*(u1(I,3)-2*u1(I,2)+u1(I,1))*((c(i,1)*dt/h3)**2))/(1+c(i,1)*dt/h3)

   u2_x(I,1)=(2*u1x(I,1)-u0x(I,1)-u2_x(I,3)+2*u1x(I,3)-u0x(I,3)+c(i,1)*dt*(u2_x(I,3)-u0x(I,3)&
   &+u0x(I,1))/h3-2*(u1x(I,3)-2*u1x(I,2)+u1x(I,1))*((c(i,1)*dt/h3)**2))/(1+c(i,1)*dt/h3)

   u2_z(I,1)=(2*u1z(I,1)-u0z(I,1)-u2_z(I,3)+2*u1z(I,3)-u0z(I,3)+c(i,1)*dt*(u2_z(I,3)-u0z(I,3)&
   &+u0z(I,1))/h3-2*(u1z(I,3)-2*u1z(I,2)+u1z(I,1))*((c(i,1)*dt/h3)**2))/(1+c(i,1)*dt/h3)

END DO
endif

DO K=1,NK_local 
   U2(1,K)=(2*U1(1,K)-U0(1,K)-U2(3,K)+2*U1(3,K)-U0(3,K)+c(1,K)*dt*(U2(3,K)-U0(3,K)&
   &+U0(1,K))/h-2*(U1(3,K)-2*U1(2,K)+U1(1,K))*((c(1,K)*dt/h)**2))/(1+c(1,K)*dt/h) 

   U2(NI_Global,K)=(2*U1(NI_Global,K)-U0(NI_Global,K)-U2(NI_Global-2,K)+&
   &2*U1(NI_Global-2,K)-U0(NI_Global-2,K)+c(NI_Global,K)*dt*(U2(NI_Global-2,K)&
   &+U0(NI_Global,K)-U0(NI_Global-2,K))/h-2*(U1(NI_Global,K)-2*U1(NI_Global-1,K)&
   &+U1(NI_Global-2,K))*((c(NI_Global,K)*dt/h)**2))/(1+c(NI_Global,K)*dt/h)
        
   W(1,K)=W(1,K)+U2(1,K)*U(1,K,T-1)	 
   W(NI_Global,K)=W(NI_Global,K)+U2(NI_Global,K)*U(NI_Global,K,T-1)	
 END DO

 
 IF(Myrank.EQ.(NPz-1)) then 
 DO I=1,NI_Global
   U2(I,NK_local)=(2*U1(I,NK_local)-U0(I,NK_local)-U2(I,NK_local-2)+&
   &2*U1(I,NK_local-2)-U0(I,NK_local-2)+c(I,NK_local)*dt*(U2(I,NK_local-2)&
   &+U0(I,NK_local)-U0(I,NK_local-2))/z-2*(U1(I,NK_local)-2*U1(I,NK_local-1)&
   &+U1(I,NK_local-2))*((c(I,NK_local)*dt/z)**2))/(1+c(I,NK_local)*dt/z)
   W(I,NK_local)=W(I,NK_local)+U2(I,NK_local)*U(I,NK_local,T-1)
 END DO
EndIF
 
DO K=1,NK_local       
   u2_x(1,K)=(2*U1x(1,K)-U0x(1,K)-u2_x(3,K)+2*U1x(3,K)-U0x(3,K)+c(1,K)&
   &*dt*(u2_x(3,K)-U0x(3,K)+U0x(1,K))/h-2*(U1x(3,K)-2*U1x(2,K)+U1x(1,K))&
   &*((c(1,K)*dt/h)**2))/(1+c(1,K)*dt/h) 

   u2_x(NI_Global,K)=(2*U1x(NI_Global,K)-U0x(NI_Global,K)-u2_x(NI_Global-2,K)&
   &+2*U1x(NI_Global-2,K)-U0x(NI_Global-2,K)+c(NI_Global,K)*dt*(u2_x(NI_Global&
   &-2,K)+U0x(NI_Global,K)-U0x(NI_Global-2,K))/h-2*(U1x(NI_Global,K)-2*U1x(NI_Global&
   &-1,K)+U1x(NI_Global-2,K))*((c(NI_Global,K)*dt/h)**2))/(1+c(NI_Global,K)*dt/h)
END DO

 IF(Myrank.EQ.(NPz-1)) then 
 DO I=1,NI_Global
  u2_x(I,NK_local)=(2*U1x(I,NK_local)-U0x(I,NK_local)-u2_x(I,NK_local-2)+&
  &2*U1x(I,NK_local-2)-U0x(I,NK_local-2)+c(I,NK_local)*dt*(u2_x(I,NK_local-2)&
  &+U0x(I,NK_local)-U0x(I,NK_local-2))/z-2*(U1x(I,NK_local)-2*U1x(I,NK_local-1)&
  &+U1x(I,NK_local-2))*((c(I,NK_local)*dt/z)**2))/(1+c(I,NK_local)*dt/z)
 END DO
 EndIF


DO K=1,NK_local       
   u2_z(1,K)=(2*U1z(1,K)-U0z(1,K)-u2_z(3,K)+2*U1z(3,K)-U0z(3,K)+c(1,K)&
   &*dt*(u2_z(3,K)-U0z(3,K)+U0z(1,K))/h-2*(U1z(3,K)-2*U1z(2,K)+U1z(1,K))&
   &*((c(1,K)*dt/h)**2))/(1+c(1,K)*dt/h) 

   u2_z(NI_Global,K)=(2*U1z(NI_Global,K)-U0z(NI_Global,K)-u2_z(NI_Global-2,K)&
   &+2*U1z(NI_Global-2,K)-U0z(NI_Global-2,K)+c(NI_Global,K)*dt*(u2_z(NI_Global&
   &-2,K)+U0z(NI_Global,K)-U0z(NI_Global-2,K))/h-2*(U1z(NI_Global,K)-2*U1z(NI_Global&
   &-1,K)+U1z(NI_Global-2,K))*((c(NI_Global,K)*dt/h)**2))/(1+c(NI_Global,K)*dt/h)
END DO

 IF(Myrank.EQ.(NPz-1)) then 
 DO I=1,NI_Global
  u2_z(I,NK_local)=(2*U1z(I,NK_local)-U0z(I,NK_local)-u2_z(I,NK_local-2)+&
  &2*U1z(I,NK_local-2)-U0z(I,NK_local-2)+c(I,NK_local)*dt*(u2_z(I,NK_local-2)&
  &+U0z(I,NK_local)-U0z(I,NK_local-2))/z-2*(U1z(I,NK_local)-2*U1z(I,NK_local-1)&
  &+U1z(I,NK_local-2))*((c(I,NK_local)*dt/z)**2))/(1+c(I,NK_local)*dt/z)
 END DO
 EndIF

transbegin=Mpi_Wtime()
!exchage date by using so called "odd-even method"
!==================================deal with U

 IF (mod(Myrank,2).EQ.0) Then
    CALL Mpi_Send(U2(1,Kstart+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(U2(1,Kend-1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Recv(U2(1,Kend),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(U2(1,Kstart),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
 ELSE
    CALL Mpi_Recv(U2(1,Kend),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(U2(1,Kstart),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Send(U2(1,Kstart+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(U2(1,Kend-1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
 ENDIF
!deal with Ux
 IF (mod(Myrank,2).EQ.0) Then
    CALL Mpi_Send(U2_x(1,Kstart+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(U2_x(1,Kend-1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Recv(U2_x(1,Kend),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(U2_x(1,Kstart),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
 ELSE
    CALL Mpi_Recv(U2_x(1,Kend),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(U2_x(1,Kstart),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Send(U2_x(1,Kstart+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(U2_x(1,Kend-1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
 ENDIF

!deal with Uz
 IF (mod(Myrank,2).EQ.0) Then
    CALL Mpi_Send(U2_z(1,Kstart+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(U2_z(1,Kend-1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Recv(U2_z(1,Kend),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(U2_z(1,Kstart),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
 ELSE
    CALL Mpi_Recv(U2_z(1,Kend),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Recv(U2_z(1,Kstart),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Status,Ierr)
    CALL Mpi_Send(U2_z(1,Kstart+1),1,Xtype,Mytop,vibrate+100,Mpi_Comm_World,Ierr)
    CALL Mpi_Send(U2_z(1,Kend-1),1,Xtype,Myfloor,vibrate+100,Mpi_Comm_World,Ierr)
 ENDIF

 transend=Mpi_Wtime()
 transtime=transtime+transend-transbegin

!==========taper==========

do k =Kstart,Kend
	do i=1,ntaper-1
	     u0(i,k)=u1(i,k)*aa(i)
         u1(i,k)=u2(i,k)*aa(i)
         u0x(i,k)=u1x(i,k)*aa(i)
         u1x(i,k)=u2_x(i,k)*aa(i)
         u0z(i,k)=u1z(i,k)*aa(i)
         u1z(i,k)=u2_z(i,k)*aa(i)

		 u0(NI_Global-i+1,k)=u1(NI_Global-i+1,k)*aa(i)
         u1(NI_Global-i+1,k)=u2(NI_Global-i+1,k)*aa(i)
         u0x(NI_Global-i+1,k)=u1x(NI_Global-i+1,k)*aa(i)
         u1x(NI_Global-i+1,k)=u2_x(NI_Global-i+1,k)*aa(i)
         u0z(NI_Global-i+1,k)=u1z(NI_Global-i+1,k)*aa(i)
         u1z(NI_Global-i+1,k)=u2_z(NI_Global-i+1,k)*aa(i)
	enddo
	do i = ntaper, NI_Global-ntaper+1 
		 u0(i,k)=u1(i,k)
         u1(i,k)=u2(i,k)
         u0x(i,k)=u1x(i,k)
         u1x(i,k)=u2_x(i,k)
         u0z(i,k)=u1z(i,k)
         u1z(i,k)=u2_z(i,k)
	enddo
enddo

if(myrank.EQ.0)then
 DO I=1,NI_Global
      DO k=1,ntaper-1
         u0(i,k)=u0(i,k)*aa(k)
         u1(i,k)=u1(i,k)*aa(k) 
		 u0x(i,k)=u0x(i,k)*aa(k)
         u1x(i,k)=u1x(i,k)*aa(k)  
		 u0z(i,k)=u0z(i,k)*aa(k)
         u1z(i,k)=u1z(i,k)*aa(k)   
	  END DO
   END DO
end if

!=======T循环的end=========
end do

!=======l循环的end=========
end do


endtime=Mpi_Wtime()
 totletime=endtime-startime


! input the results.
! time
 If (Myrank.eq.0) then
  OPEN(11,FILE='time.txt',STATUS='unknown')
    WRITE(11,*) totletime,'second is used to finish this program'
    WRITE(11,*) transtime,'second is used to transfer data form different node.'
  CLOSE(11)
 EndIf
  
! DO k=1,NK_local
!   DO I=1,NI_Global
!      xzu(i,k)=w(i,k)
!	 end do
!  end do

 
!u----------------------------
  If (Myrank.NE.0) Then
     CALL MPi_Send(w(1,1),1,xZtype,0,vibrate+100,Mpi_Comm_World,Ierr)
  Else
     Do I=1,NPz-1
        CALL Mpi_Recv(w(1,I*NK_local+1),NK_local*NI_Global,MPI_DOUBLE_PRECISION,I,&
            &vibrate+100,Mpi_Comm_World,Status,Ierr)
     EndDo
  EndIf

!==========write three component in the xy-, xz-, and yz-planes=========================

!case 1: input data as binary format which can be used to draw snapshot with eps format.
!  If (Myrank.EQ.0) Then
!write the snapshot in the xy-plane
!    Open (11,FILE='ONAD_2D_Uxz.dat',access='direct',recl=8*NI_Global*NK_Global,IOSTAT=ios)
!        Write (11, rec=1) ((XZU(i,k),k=1,NK_Global),i=1,NI_Global)
!     Close(11)
!  end if

!write the snapshot in the xz-plane
!u----------------------------
!  If (Myrank.NE.0) Then
!     CALL MPi_Send(c(1,1),1,xZtype,0,vibrate+100,Mpi_Comm_World,Ierr)
!  Else
!     Do I=1,NPz-1
!        CALL Mpi_Recv(XZU(1,I*NK_local+1),NK_local*NI_Global,MPI_DOUBLE_PRECISION,I,&
!            &vibrate+100,Mpi_Comm_World,Status,Ierr)
!     EndDo
!  EndIf

!==========write three component in the xy-, xz-, and yz-planes=========================

!write the snapshot in the xz-plane
 If (Myrank.EQ.0) Then
 OPEN(24,FILE='onalong.grd',STATUS='unknown')
 WRITE(24,240) 'DSAA'
 240 FORMAT(A4)
 WRITE(24,"(I3,1X,I4)") Nk_Global,NI_Global
 WRITE(24,"(I1,1X,I3)") 0,Nk_Global
 WRITE(24,"(I1,1X,I4)") 0,NI_Global
 max=w(1,1)
 min=max
 DO I=1,NI_Global
    DO K=1,Nk_Global
       IF(w(I,k).GT.MAX) THEN
	   max=w(I,k)
       END IF
       IF(w(I,k).LT.MIN) THEN
	   min=w(I,k)
       END IF
    END DO
 END DO
 WRITE(24,"(E15.7,1X,E15.7)") min,max
    DO I=1,NI_Global
	DO K=1,NK_Global
	 WRITE(24,*) w(I,k)
	END DO
    END DO
 CLOSE(24)
end if 

  CALL Mpi_Finalize(Ierr)
End
