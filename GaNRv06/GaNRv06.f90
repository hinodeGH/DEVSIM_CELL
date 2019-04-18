!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     バルクシリコン内の電子輸送のための
!     多粒子モンテカルロプログラム
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! rv00 2018.12.26 BULK_SI_rv0803の骨組み
! rv01 2018.12.27 骨組みにGaAsを組み込む
! rv01a 2019.01.17-- Γ谷のみのプログラムに単純化
!		2019.02.20-- phonon_count_ON_ratio: steady-state fraction
! GaNRv01 2019.02.27 change parameters from GaAs to GaN
! GaNRv02 2019.03.05 全体の1/10時間ごとに電子エネルギー分布を出力
! GaNRv03 2019.03.07 電子エネルギー分布を面積ノーマライズ(推移比較のため)
!
!===( 乱数発生関数 )=== 
! 2018.08.16 rnd() -> system random number:call random_number()
! 2018.08.23 rv04  -> Mersenne random number: call grnd()
! 2018.08.25 range mt19937 -> mt19937rv0a:[0 1)  -> (0 1]  for log(grnd)
! 2018.10.25 impliment function Ran1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
module M_random_number
!
!	mt19937rv0a 
	implicit none
	integer,private :: N, N1, M, MATA, UMASK, LMASK, TMASKB, TMASKC 
	parameter(& 
              & N = 624, & 
              & N1 = 625, & 
              & M = 397, & 
              & MATA = -1727483681, & 
              & UMASK = -2147483648, & 
              & LMASK = 2147483647, & 
              & TMASKB = -1658038656, & 
              & TMASKC = -272236544 & 
              & ) 
	integer,private :: mti = N1, mt(0:N-1), mag01(0:1) = (/0, MATA/) 
!
	contains
! 
subroutine sgrnd(seed) 
	integer,intent(in) :: seed 
! 
!		setting initial seeds to mt[N] using 
!		the generator Line 25 of Table 1 in 
!		[KNUTH 1981, The Art of Computer Programming 
!		Vol. 2 (2nd Ed.), pp102] 
! 
	mt(0) = iand(seed, -1) 
	do mti = 1, N - 1 
		mt(mti) = iand(69069 * mt(mti - 1), -1) 
	end do 
end subroutine sgrnd 
!
Double Precision function grnd()	! real(8) -> Double Precision *********
	integer :: y, kk 
!
	if(mti >= N) then 
!						generate N words at one time 
		if(mti == N + 1) then 
!						if sgrnd() has not been called, 
			call sgrnd(4357) 
!						a default initial seed is used 
		endif 
!
		do kk = 0, N - M - 1 
			y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK)) 
			mt(kk) = ieor(ieor(mt(kk + M), ishft(y, -1)), mag01(iand(y, 1))) 
		end do 
!
		do kk = N - M, N - 2 
			y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK)) 
			mt(kk) = ieor(ieor(mt(kk + (M - N)), ishft(y, -1)), mag01(iand(y, 1))) 
		end do 
!
		y = ior(iand(mt(N - 1), UMASK), iand(mt(0), LMASK)) 
		mt(N - 1) = ieor(ieor(mt(M - 1), ishft(y, -1)), mag01(iand(y, 1))) 
		mti = 0 
	endif 
!
	y = mt(mti) 
	mti = mti + 1 
	y = ieor(y, ishft(y, -11)) 
	y = ieor(y, iand(ishft(y, 7), TMASKB)) 
	y = ieor(y, iand(ishft(y, 15), TMASKC)) 
	y = ieor(y, ishft(y, -18)) 
!
	if(y < 0) then 
		grnd = (dble(y) + 2.d0 ** 32) / (2.d0 ** 32) 
	else 
		grnd = dble(y) / (2.d0 ** 32) 
	endif 
!
	grnd = 1.d0 - grnd  !  range [0 1)  -> (0 1]  for log(grnd)
!      
end function grnd 
!........................................................................................
! 
Double Precision function gasdev() !gaussian random number:Box Muller's method
	implicit none
	integer :: iset=0   ! 0 setting is only once at first function call
	Double Precision gset
	Double Precision  fac,rsq,v1,v2
	SAVE iset, gset
	integer cnt				!! debug
!
	if  (iset == 0) then
		v1=2.d0*grnd()-1.d0
		v2=2.d0*grnd()-1.d0
		rsq=v1*v1+v2*v2
		cnt=1				!! debug
		do while (rsq >= 1.0 .OR. rsq == 0.0)
			v1=2.d0*grnd()-1.d0
			v2=2.d0*grnd()-1.d0
			rsq=v1*v1+v2*v2
			cnt=cnt+1				!! debug
			if (cnt==10000) then				!! debug
				write(*,*) 'gasdev while cnt=', cnt				!! debug
				write(8,*) 'gasdev while cnt=', cnt				!! debug
			end if				!! debug
		end do
		fac=sqrt(-2.d0*log(rsq)/rsq)
		gset=v1*fac
		iset=1
		gasdev = v2*fac
	else
		iset=0
		gasdev = gset
	end if
end function gasdev
!========================================================================================
end module M_random_number
!****************************************************************************************
!						ran1 at Numerical Recipes in Fortran 2nd Ed. p.271	
module N_random_number
	implicit none
	integer,private :: IA,IM,IQ,IR,NTAB,NDIV
	Double Precision,private ::  AM,EPS,RNMX
	Parameter(IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=3.d-16,RNMX=1.-EPS)
!
!	Minimal random number generator of Park and Miller with Bays-Durham shuffle and added safegurds.
!	Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values).
!	Call with idum a negetive integer to initialize; thereafter, do not alter idum between successive deviates in a sequence.
!	RNMX should approximate the largest floating value that is less than 1.
!
	integer,private ::  j,k,iv(NTAB),iy,idum
	SAVE iv, iy, idum
	DATA iv /NTAB*0/, iy /0/
!
	contains
!
subroutine sgrnd(seedGen)
	integer seedGen
!
	idum=seedGen
	if (idum.LE.0. .OR. iy.EQ.0) then 		!	initialize.
		idum=max(-idum,1)					!	Be sure to prevent idum=0
		do j=NTAB+8, 1, -1					!	Load the shuffle table (after 8 warm-ups)
			k=idum/IQ
			idum=IA*(idum-k*IQ)-IR*k
			if (idum .LT. 0) idum=idum+IM
			if (j .LE. NTAB) iv(j)=idum
		end do
		iy=iv(1)
	end if
end subroutine sgrnd
!
Double Precision function grnd()
!
	k=idum/IQ								!	Start here when not initializing.
	idum=IA*(idum-k*IQ)-IR*k				!	Compute idum=mod(IA*idum,IM) without overflows by Schrange's method.
	if (idum .LT. 0) idum=idum+IM	
	j=1+iy/NDIV								!	Will be in the range 1:NTAB.
	iy=iv(j)								!	Output previously stored value and refill the shuffle table.
	iv(j)=idum
	grnd=min(AM*iy,RNMX)					!	Because users don't expect endpoint values.
	return
END function grnd
!........................................................................................
! 
Double Precision function gasdev() !gaussian random number:Box Muller's method
	implicit none
	integer :: iset=0   ! 0 setting is only once at first function call
	Double Precision gset
	Double Precision  fac,rsq,v1,v2
	SAVE iset, gset
	integer cnt				!! debug
!
	if  (iset == 0) then
		v1=2.d0*grnd()-1.d0
		v2=2.d0*grnd()-1.d0
		rsq=v1*v1+v2*v2
		cnt=1				!! debug
		do while (rsq >= 1.0 .OR. rsq == 0.0)
			v1=2.d0*grnd()-1.d0
			v2=2.d0*grnd()-1.d0
			rsq=v1*v1+v2*v2
			cnt=cnt+1				!! debug
			if (cnt==10000) then				!! debug
				write(*,*) 'gasdev while cnt=', cnt				!! debug
				write(8,*) 'gasdev while cnt=', cnt				!! debug
			end if				!! debug
		end do
		fac=sqrt(-2.d0*log(rsq)/rsq)
		gset=v1*fac
		iset=1
		gasdev = v2*fac
	else
		iset=0
		gasdev = gset
	end if
end function gasdev
!========================================================================================
end module N_random_number
!
!****************************************************************************************
module Tom2Pop_sub_program_variables    !----- 共通変数 -----
!		from Pop_sub_program_variables
	implicit none
!	Double Precision,parameter::EMIN = 0.d-10	! Pop 1.d-10 << 0.0005(=de)
	Double Precision,parameter::EMIN = 1.d-10	! 2019.03.13
!	Double Precision,parameter::EMAX = 2.0D0	! ~2eV
!	Double Precision,parameter::EMAX = 9.9D0	! 9.9eV
	Double Precision,parameter::EMAX = 4.0D0	! 4eV	GaN
	Double Precision,parameter::QMIN = 1.d0	! // nearly zero compared to QMAX
	Double Precision,parameter::QMAX = 1.256d+08		!  2*pi/aGaN(lattice constant, cm)
	Double Precision,parameter::Rws = 1.954d-8         !  aGaN*(3/(16*pi)^(1/3); aGaN=5.0d-8 cm
	Double Precision,parameter::m0 = 5.68562975d-16    ! eV*s^2/cm^2  !electron mass
!
	Double Precision,parameter::eps0 = 5.52634972d+05	!	// e/V/cm permittivity of free space
	Double Precision,parameter::epsGaN = 5.2500322d+06	!	// e/V/cm permittivity of GaN 9.5*eps0
	Double Precision,parameter::epsf =   2.9565971d+06	!	// e/V/cm permittivity of GaN 5.35*eps0
	Double Precision,parameter::QmidBZ = 0.5d0	!	middle of B.Z.
	Double Precision,parameter::QcntrBZ = 0.0d0	!	center of B.Z.
!
	integer,parameter:: LA=1
	integer,parameter:: TA=2
    integer,parameter:: LO=3
    integer,parameter:: TO=4
!
    Double Precision,parameter::pi = 3.14159265d0        !π
    Double Precision,parameter::rho = 3.8020d+12       ! density (eV*s^2/cm^5) (6.1 g/cm^3)
    Double Precision,parameter::kb = 8.61734318d-05    ! eV/K
    Double Precision,parameter::hbar = 6.58211915d-16  ! eV*s
!
    Double Precision,parameter::DLA = 10.0d0	!	8~12 eV
    Double Precision,parameter::DTA = 0.0d0
    Double Precision,parameter::Absorption = -1.d0
    Double Precision,parameter::EMS = 1.d0
    Double Precision,parameter::PQMAX = 8.26714d-08 ! hbar*QMAX(Wigner-Seitz cell radius)
    
	Double Precision,parameter::LTAO(4,3) = &
						&  reshape( (/	-8.65d-4,	-1.45d-3,	-1.88d-4,	1.98d-4, &
						&				4.91d5,		4.88d5,		-1.45d4,	2.52d4,   &
						&				0.0d0,		0.0d0,		1.38d14,	1.07d14 /),(/4,3/) )
!
	Double Precision,parameter::echarge = 1.d0		!  // electron charge 
!				#define echarge   -1.0             // electron charge 
!				#define ecoulom   -1.60217653d-19  // electron charge in Coulombs
!
!		from Tom_global_variables   !----- 共通変数 -----
	Double Precision,SAVE:: dt,Temp,fx
	Double Precision,SAVE:: de
	Double Precision,SAVE:: swk(10,4002)=0.d0
	Double Precision,SAVE:: mdos
	Double Precision,SAVE:: gm
	Double Precision,SAVE:: alpha
	Double Precision,SAVE:: qhbar
	Double Precision,SAVE:: kbTq
!!	Double Precision,SAVE:: ef		!	2019.03.28 check
	integer,SAVE:: inum,jtl
	integer,SAVE:: iemax
	integer, parameter::PX=5,PY=6,PZ=7,EE=8,TS=9,XXX=1,YYY=2
	integer, parameter::EElost=10
!
	Double Precision,SAVE:: Doping
	Double Precision,SAVE:: Egetlost_scat, Egetlost_drift		!	Energy conservation
!
	Double Precision,SAVE:: eee_convergence, epsLL		!	final_state_intra_ScatByLATA/final_state_intra_ScatByLOTO
	Double Precision, parameter::phonon_count_ON_ratio=0.5	!	This ratio of Last simulation time is supposed to be steady-state
!
end module Tom2Pop_sub_program_variables
!****************************************************************************************
module Pop_sub_programs
	use Tom2Pop_sub_program_variables
	use M_random_number
!	use N_random_number
	implicit none
	contains
!--------------- 以下に次のsubroutine/function ------------
!	trapz 		! 	数値積分の足し算
!	get_wq 		!	phonon q -> w
!	get_nq 		! 	phonon occupation
!	get_Mq 		! 	(3.8) 積分式の被積分関数
!	fabs_cost  	! 	(3.12) cos(φ)式の絶対値
!	find_root  	! 	|cos(φ)|-1=0の根を探索
!	acrate3    	!	E毎にPop(3.8) 積分
!	rateDOP		!	E毎のイオン化不純物散乱レート
!
!========================================================================================
Double Precision function trapz(xx, yy, n1, n2)
	implicit none
    Double Precision, intent(in)::xx(:), yy(:)
    integer, intent(in)::n1, n2
	Double Precision s
	integer i
!	
	s=0.
	do i=n1, n2-1
		s =s + 0.5*(yy(i+1)+yy(i))*(xx(i+1)-xx(i))
	end do
	trapz = s
end function trapz
!========================================================================================
Double Precision function get_wq(q, lt) ! -- ( phonon frequency ) --
	implicit none
	Double Precision q, qTemp
	integer, intent(in)::lt
	integer LTtemp
!	
	LTtemp=lt
	qTemp=q
!	
	if (qTemp > QMAX) then
		if (LTtemp == LA) LTtemp = LO			!	LA outside BZ extends into LO
		qTemp = 2.*QMAX - qTemp					!	and flip it back to first BZ
	else if (qTemp < 0.) then
		qTemp = -qTemp
	end if
!						here LTtemp = 1 (LA), 2 (TA), 3 (LO), 4 (TO) as in constants.h
	get_wq = ( LTAO(LTtemp,1)*qTemp*qTemp + LTAO(LTtemp,2)*qTemp + LTAO(LTtemp,3))
end function get_wq
!========================================================================================
Double Precision function get_nq(q, lt)       ! phonon occupation
	implicit none
	Double Precision, intent(in)::q
	integer, intent(in)::lt
	Double Precision xx
!	
	xx = hbar*get_wq(q,lt)/(kb*Temp)
!						/* if (xx < 3.5d0) return (1.d0/xx -0.5d0 +xx/12.d0 -xx*xx*xx/720.d0)
!						else return exp(-xx) */
	get_nq = (1.d0/(exp(xx)-1.d0))
end function get_nq
!========================================================================================
Double Precision function get_Mq(q, lt, aed)
	implicit none
	Double Precision, intent(inout)::q
	Double Precision, intent(in)::aed
	integer, intent(in)::lt
	Double Precision gq
!	
	if (q > QMAX) q = 2.d0*QMAX-q               !    // returns to the first BZ
!					#ifdef POP
		gq = 3.d0*(sin(q*Rws)-q*Rws*cos(q*Rws))/((q*Rws)**3.d0)
!					#else
!					gq = 1.d0             !  // overlap integral
!					#endif
	get_Mq = (q*q*gq*gq*(get_nq(q,lt)+0.5d0+0.5d0*aed)/get_wq(q,lt))
end function get_Mq
!========================================================================================
Double Precision function fabs_cost(q, k, E, lt, aed)
	implicit none
	Double Precision, intent(in)::q, k, E, aed
	integer, intent(in)::lt
	Double Precision wq
!					 |cos(t)| = |+/-q/(2k) + m*w(q)/(hbar*k*q)|   +em/-ab
	wq = get_wq(q,lt)
!					 *****  aed = -1. for absorption,  +1. for emission
	fabs_cost = abs( aed*q/(2.d0*k) + mdos*(1.d0 + alpha*(2.d0*E - aed*hbar*wq))*wq/(hbar*k*q) )
!
end function fabs_cost
!========================================================================================
Double Precision function find_root(q1, q2, k, E, lt, aed)
	implicit none
	Double Precision, intent(inout)::q1, q2
	Double Precision, intent(in):: k, E, aed
	integer, intent(in)::lt
	Double Precision f1, f2, fs, q, eps
	integer cnt				!! debug
!					 find root of |cos(t)| - 1 in interval q1..q2 */
	q=q2
	eps = (q2-q1)/40.d0               ! 20.
	f1 = fabs_cost(q1,k,E,lt,aed)
	f2 = fabs_cost(q2,k,E,lt,aed)
	fs = (f2-f1)/abs(f2-f1)
	cnt=1				!! debug
	do while (q2-q1 > eps) 
		q = 0.5d0*(q1+q2)
		if (fs*(fabs_cost(q,k,E,lt,aed)- 1.d0) > 0.) then
			q2 = q
		else
			q1 = q
		end if
		cnt=cnt+1				!! debug
		if (cnt==10000) then				!! debug
			write(*,*) 'find_root while cnt=', cnt				!! debug
			write(8,*) 'find_root while cnt=', cnt				!! debug
		end if				!! debug
	end do
	if (fs < 0.) then
		 find_root = q2 
	else 
		find_root =  q1
	end if
end function find_root
!========================================================================================
Double Precision function acrate3(E, lt, aed)
	implicit none
	Double Precision, intent(in)::E, aed
!					 *****  aed = -1. for absorption,  +1. for emission
	integer, intent(in)::lt
	Double Precision,allocatable ::q(:), integ(:)
	Double Precision q1, q2, dq, k, acost, C3, G, DP
	integer i, i1, i2, nloc, N, MULT
!                        /******  ******/
	N=200              ! 200 100
	MULT=20            ! 20  10
!      
	allocate(q(MULT*N))
	allocate(integ(MULT*N))
!						 /* k is the electron wavevector in Herring-Vogt space */
	k = sqrt(2.d0*mdos*E*(1.+alpha*E))/hbar
!
	q1 = QMIN
	q2 = QMAX
	i1 = 0
	i2 = N
	dq = (QMAX-QMIN)/real(N-1)
	q(1) = QMIN
!  
	do i=1, N
		if (i >= 2) then
			q(i) = q(i-1)+dq
!						/***  cos(t) = +/- q/(2k) + m*w(q)/(hbar*k*q)   +em/-ab  ***/
		end if
		acost = fabs_cost(q(i),k,E,lt,aed)
!						/* find range of q's [q1..q2] where the cosine is non-zero */
		if ((acost <= 1.) .and. (i1 == 0)) then
			i1 = i
			if (i >= 2) then
				q1 = find_root(q(i-1),q(i),k,E,lt,aed)
			end if
		else if ((acost >  1.) .and. (i2 == N) .and. (i1 == 1)) then
			i2 = i-1
			q2 = find_root(q(i-1),q(i),k,E,lt,aed)
		end if
	end do
!
	if (i1 == 0) then                !      /* in case no q-values were found */
		G = 0.
	else                             ! ⑥ begin /* found q-range to integrate over */
		if (i2 == N) then
			write(*,*) '*** Hit BZ End while integrating acoustic rates'
			write(*,*) 'at E=',E,'  Reduce EMAX!'
			write(*,*) ' or continue at your own risk (continuing works, sort of)'
!						exit(0)
			Return			! 2019.04.04
		end if
		nloc = MULT*int((q2-q1)/dq)
		dq = (q2-q1)/ real(nloc-1)
!	
		q(1) = q1
!						/* compute q-dependent part of the integrand */
		do i = 1, nloc
			if (i >= 2) then 
				q(i) = q(i-1)+dq
			end if
			integ(i) = q(i)*get_Mq(q(i),lt,aed)
		end do
!	
		if (lt == LA) then              !{  /* longitudinal */
			DP = DLA
		else                            !   /* transverse */
			DP = DTA
		end if
!
		C3 = DP*DP*mdos/(4.*pi*rho*hbar*hbar*k) 
		G = C3*trapz(q,integ,1,nloc) 
	end if                               ! ⑥ end
!  
	deallocate(q)
	deallocate(integ)
!
	acrate3 = G
end function acrate3
!========================================================================================
Double Precision function rateDOP(Ed)
	implicit none
	Double Precision Ed, beta,gamma, ks,v0,xx,Gx,bb
!	
	if (DOPING >= 1.d12) then
		if (Ed < DE) Ed = DE
		beta = sqrt(DOPING/(epsGaN*kb*Temp))
		gamma = Ed*(1.+alpha*Ed)
		ks = sqrt(2.d0*mdos*gamma)/hbar
		v0 = hbar*ks/(mdos*(1.d0+2.d0*alpha*Ed))
		xx = hbar*hbar*(2.d0*ks)*(2.d0*ks)/(8.d0*mdos*kb*Temp)
		Gx = (1.d0+0.165651d0*xx**2.+0.0557896d0*xx**4.d0)/ 	&
		&	(1.d0+0.829201d0*xx**2.d0+0.343064d0*xx**4.d0+0.111354d0*xx**6d0)
		bb = 4.d0*ks*ks/(beta*beta*Gx)
		rateDOP = DOPING/(2.d0*pi*hbar*hbar*epsGaN*epsGaN*v0*4.d0*ks*ks)
		rateDOP = 0.8d0*2.d0*rateDOP*(log(1.d0+bb) - bb/(1.d0+bb))
	else
		rateDOP = 0.d0
	end if
!
	return
end function rateDOP
!========================================================================================
end module Pop_sub_programs  
!========================================================================================
!****************************************************************************************   

!========================================================================================
module Tom_sub_programs
	use Tom2Pop_sub_program_variables
	use M_random_number
!	use N_random_number
	use Pop_sub_programs
	implicit none
	contains
!
!----- 以下に次のサブルーチン ------------
!     data_input:          総粒子数, 総時間ステップ, 時間増分, 格子温度, 電界
!     param:               物理定数, 材料定数および諸パラメータ;散乱レートの計算→表化
!     initia:              全粒子の初期条件
!     monte_carlo:         多粒子モンテカルロ計算：繰り返し
!       emc:               時間dtの間の多粒子モンテカルロ計算
!         drift:           ドリフト過程の計算
!         scat:            散乱過程の計算
!       out:               出力ルーチン
!========================================================================================
!===( Input Simulation condition )===
!
subroutine data_input(inum_dummy,SimTime)
	implicit none
	integer, intent(out):: inum_dummy
	Double Precision, intent(out):: SimTime

!---( 入力データ )---
!
!=== 総粒子数, 総時間ステップ, 時間増分, 格子温度, 電界
!
	read(5,*) inum     ! number of carriers (typ. 1E4~1E6)
	inum_dummy=inum    ! for allocation of Elec(inum,9)
	read(5,*) jtl      ! number of iteration (typ. 1000~10000)
	read(5,*) dt       ! time step [s] (typ. 5E-15 = 5fs)
	SimTime=jtl*dt		! Total Simulation Time [s] 
	read(5,*) Temp      ! temperature [K] (typ.300)
	read(5,*) fx       ! electric field [KV/cm] 
	read(5,*) Doping	!	Dopant Concentration [1/cm3]
	read(5,*) epsLL		!	convergence limit: 1%~1ppm?
	eee_convergence = epsLL
!
	write(*,*) inum/1D6,' M particles'
	write(*,*) jtl, ' iteration'
	write(*,*) dt, ' sec'
	write(*,*) Temp, ' K'
	write(*,*) fx, ' KV/cm'
	write(*,*) Doping, ' /cm3'
	write(*,*) epsLL, ' convergence limit'
!
	write(8,*) inum/1D6,' M particles'
	write(8,*) jtl, ' iteration'
	write(8,*) dt, ' sec'
	write(8,*) Temp, ' K'
	write(8,*) fx, ' KV/cm'
	write(8,*) Doping, ' /cm3'
	write(8,*) epsLL, ' convergence limit'
!	
	fx=fx*1000.d0   !  kV/cm -> V/cm
!
end subroutine data_input
!========================================================================================
!===( 散乱レート計算用関数 )===
   Double Precision function DOSe(Ed)
    implicit none
    Double Precision Ed,sgamma,dosconst
!    
    if (Ed > 0.0) then
	  sgamma = sqrt(Ed*(1.d0+alpha*Ed)) 
!			// includes 2x for spin, else it would be 4x in the denominator
	  dosconst = 2.d0*pi*pi*hbar*hbar*hbar
	  DOSe = ((2.d0*mdos)**1.5d0)*sgamma*(1.d0 + 2.d0*alpha*Ed)/dosconst
    else
      DOSe=0.d0
    end if
  end function DOSe
!========================================================================================
!===( 物理定数, 材料定数および諸パラメータ )===
!
subroutine param
	implicit none      
	Double Precision Egap
	Double Precision ei,ef,sei,sef
	Double Precision ep,hwLO,Nq,poe,poa
	Double Precision ph_qmax, ph_qmin
	integer ie, i
!
!!!			remind LT=   LA=1, TA=2, LO=3, TO=4
!		 aed = -1. for absorption,  +1. for emission
	ep = 1.0/(1.0/epsf-1.0/epsGaN)
!  /************************************************************************/
!  /**************** INTERVALLEY SCATTERING (LO+TO+LA+TA) ******************/
!  /************************************************************************/
!!	Double Precision scatOconst, Nq
!    
!---( 電子の縦有効質量, 横有効質量, 禁制帯幅, 誘電率 )---
!
!	Egap = 1.5216 - 8.871e-4*Temp*Temp/(Temp+572)
	Egap = 3.5
	write(*,*) 'Egap(GaN)=',Egap					! debug
!	write(8,*) 'Egap(GaN)=',Egap					! debug
!
!---( フォノン散乱の諸パラメータ )---
	kbTq = kb*Temp/echarge
	qhbar   = echarge/hbar     ! h->hbar
	write(*,*) hbar,QcntrBZ,QmidBZ, QMAX, LO			!	debug
!	hwLO = hbar*get_wq(QcntrBZ*QMAX,LO)			!	フォノンエネルギーの最大値はqがBZ原点
	hwLO = 0.092		!	92meV
	Nq = 1.d0/(exp(hwLO/(kb*Temp))-1.d0)
	write(*,*) 'hwLO=',hwLO,' hwLO/(kb*Temp)=',hwLO/(kb*Temp),' Nq=',Nq					! debug
	poe = hwLO*(Nq+1)/(8.0*pi*ep*hbar)
	poa = hwLO*Nq/(8.0*pi*ep*hbar)
!
!---( バンドの非等方性 )---
!
	mdos = 0.20*m0
!	
!
!---( バンドの非放物線性 )---
!
	alpha = (1.d0 - mdos/m0)**2/Egap
	write(*,*) 'alpha(GaN)=',alpha				! debug
!	write(8,*) 'alpha(GaN)=',alpha				! debug
!
!---( valley offset values )---
!
!---( 散乱レートの計算 )---
!
!!	de=0.0005d0   ! energy step  ! Tom(2meV)->Pop(0.5meV)
	de=0.001d0   ! energy step	GaN 1meV
!!!	de=0.002	!	2019.01.23
	iemax=MIN(4002,nint(EMAX/de)+2)	!  Tom(0-2eV;1000div)->Pop(0-1.5eV;3002div)->(0-2.0eV;4002div)
!
	do ie=1,iemax
!				ei=de*float(ie) ! energy table: de ~ de*iemax
		if (ie == 1) then
			ei = EMIN
		else
			ei=de*float(ie-1) ! energy table: EMIN=1E-10(~0) ~ de*(iemax-1)
		end if
		do i=1,10
			swk(i,ie) = 0.d0
		end do
!-----( 非有極性光学フォノン散乱 )--- GaN Non??
!-----( 有極性光学フォノン散乱 )---
!		--- phonon emission
		ef = ei - hwLO
		if (ef > 0.0) then
			sei = sqrt(ei)
			sef = sqrt(ef)
			ph_qmax = sei+sef
			ph_qmin = sei-sef
			swk(7,ie)= poe*sqrt(2.0*mdos)*sei*log(ph_qmax/ph_qmin)/(hbar*ei)
		else
			swk(7,ie)= 0.0
		end if
!		--- phonon absorption
		ef = ei + hwLO
		sei = sqrt(ei)
		sef = sqrt(ef)
		ph_qmax = sei+sef
		ph_qmin = -sei+sef
		swk(8,ie)= poa*sqrt(2.0*mdos)*sei*log(ph_qmax/ph_qmin)/(hbar*ei)
!
!=====[ 光学フォノン散乱 Pop }===== scattering type 6 - 9
! 			 // the value of the overlap integral for inter-valley scattering may
!  			 // be included in the coupling constants (Jacoboni 1983, p. 672)
!   
!=====[ 音響フォノンValley間散乱 Pop }===== scattering type 10 - 17
!
!=====[ 音響フォノン散乱 }=====  Pop scattering type #2 - 5
!
		swk(3,ie)= acrate3(ei, LA, EMS)         ! LA intra emission
		swk(4,ie)= acrate3(ei, LA, Absorption)  ! LA intra absorption
		swk(5,ie)= acrate3(ei, TA, EMS)         ! TA intra emission
		swk(6,ie)= acrate3(ei, TA, Absorption)  ! TA intra absorption
!
!=====[ イオン化不純物散乱  }===== Pop scattering type #18
		swk(2,ie)= rateDOP(ei)
!
	end do
!        
!----- debug scattering rate begin ------
!	do  ie=1,iemax
!		if (ie == 1) then
!			eee = EMIN
!		else
!			eee = de*float(ie-1)
!		end if
!		write (*,'(18E22.15)') eee,swk(2,ie),swk(3,ie),swk(4,ie),swk(5,ie),swk(6,ie),swk(7,ie),swk(8,ie)
!		write (8,'(18E22.15)') eee,swk(2,ie),swk(3,ie),swk(4,ie),swk(5,ie),swk(6,ie),swk(7,ie),swk(8,ie)
!	end do
!------ debug scattering rate end -------
!
!---( 散乱レートの和の計算 )---
!     
	gm=0.
	do  ie=1,iemax
		do i=2, 10
			swk(1,ie)=swk(1,ie)+swk(i,ie)
		end do  
		if (swk(1,ie) > gm) gm = swk(1,ie)
	end do
	do  ie=1,iemax
		swk(1,ie)=gm-swk(1,ie)
		do i=2, 10
			swk(i,ie)=swk(i-1,ie)+swk(i,ie)
		end do
		do i=1, 10
			swk(i,ie)=swk(i,ie)/gm
		end do
	end do     
!----- debug scattering rate begin ------
	write (*,*) 'total scattering rate=',gm
	write (8,*) 'total scattering rate=',gm
!	do  ie=1,iemax
!		if (ie == 1) then
!				eee = EMIN
!			else
!				eee = de*float(ie-1)
!		end if
!		write (*,'(19E22.15)') eee,swk(1,ie),swk(2,ie),swk(3,ie),swk(4,ie),swk(5,ie),swk(6,ie),swk(7,ie),swk(8,ie),swk(9,ie), &
!			&     swk(10,ie)
!		write (8,'(19E22.15)') eee,swk(1,ie),swk(2,ie),swk(3,ie),swk(4,ie),swk(5,ie),swk(6,ie),swk(7,ie),swk(8,ie),swk(9,ie), &
!			&     swk(10,ie)
!	end do
!------ debug scattering rate end -------
!
end subroutine param
!========================================================================================
!===( 全粒子の初期条件 )===
!
subroutine initia(t,Elec)
	implicit none
!      
	Double Precision, intent(out)::t
	Double Precision, intent(out)::Elec(:,:)
	Double Precision rn, ei, ak, cb, sb, fai, sf
	integer n
	Double Precision kx, ky, kz, gk
	integer MBdist  ! 1:initial=Maxwell-Boltzman distribution
	integer ReDo
!
	t=0.d0
	call sgrnd(-99)                    !乱数 seed=1 -> -99
	do  n=1,inum
	end do
	MBdist=1		!	Maxwell-Boltzman distribution
!!!	MBdist=0		!	approximate Tom(?) distribution
!
	do  n=1,inum
		ReDO=0	
!--------------------------------------------- begin isotropic Maxwell-Boltzman(Pop)
		if (MBdist == 1) then
			kx = sqrt(2.d0*mdos*echarge)/hbar*sqrt(kbTq/2.d0)*gasdev()    ! Local kx
			ky = sqrt(2.d0*mdos*echarge)/hbar*sqrt(kbTq/2.d0)*gasdev()    ! Local ky
			kz = sqrt(2.d0*mdos*echarge)/hbar*sqrt(kbTq/2.d0)*gasdev()    ! Local kz
!--------------------------------------------- begin isotropic initial p (Tom)
		else
			rn = grnd()
			ei  = -kbTq*log(rn)*1.5d0    !--- give E along exp{-E/(3/2kT)} ~39meV @300K
			ak  = sqrt(2.d0*mdos*echarge)/hbar*sqrt(ei*(1.d0+alpha*ei)) !--- E -> k
			rn = grnd()
			cb  = 1.d0-2.d0*rn ! cosθ [-1 1]
			sb  = sqrt(1.d0-cb*cb) ! sinθ
			rn = grnd()
			fai = 2.d0*pi*rn ! Φ
			sf  = sin(fai)
			kx = ak*cb*sf    ! Local kx
			ky = ak*sb*sf    ! Local ky
			kz = ak*cos(fai) ! Local kz
		end if
!--------------------------------------------- end isotropic initial k
		rn = grnd()
		Elec(n,TS) = -log(rn)/gm         ! mean free time along exp(-t/gm)
		Elec(n,XXX) = 0.d0 ! x
		Elec(n,YYY) = 0.d0 ! y
!
		Elec(n,PX)=hbar*kx
		Elec(n,PY)=hbar*ky
		Elec(n,PZ)=hbar*kz
!
		gk = Elec(n,PX)**2/(2.d0*mdos) + Elec(n,PY)**2/(2.d0*mdos) + Elec(n,PZ)**2/(2.d0*mdos)
!		Elec(n,EE)=(sqrt(1.d0+4.d0*alpha*gk)-1.d0)/(2.d0*alpha)  ! Tomizawa (1.7)
		Elec(n,EE)=2.d0*gk/(sqrt(1.d0+4.d0*alpha*gk)+1.d0)		!	2019.03.13
!
	end do		!	do  n=1,inum
end subroutine initia
!========================================================================================
!===( 多粒子モンテカルロ計算：繰り返し )===
!
subroutine monte_carlo(Elec,NELph0,NETph0,NBph,DEph)
	implicit none
	Double Precision, intent(inout)::Elec(:,:)
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph
	Double Precision, intent(in)::DEph
	Double Precision t
	integer jt
	integer phonon_count_ON
!
	phonon_count_ON=0
	do jt=0,jtl
		if (jt == Int(Int(jt/(jtl/10))*(jtl/10))) call energy_dist(Elec)	!2019.03.05 1/10期間毎に電子エネルギー分布を出力
		t=dt*float(jt)
		if (jt > phonon_count_ON_ratio*jtl) phonon_count_ON=1		!	count for last 50% iteration
		call emc(t,Elec,NELph0,NETph0,NBph,DEph,phonon_count_ON) 
		write(*,'(i10)',advance = 'no') jt   ! iteration no
		write(8,'(i10)',advance = 'no') jt   ! iteration no
		call out(t,Elec)
	end do
!
end subroutine monte_carlo
!========================================================================================
!===(時間dtの間の多粒子モンテカルロ計算 )===
!
	subroutine emc(t,Elec,NELph0,NETph0,NBph,DEph,phonon_count_ON)
	implicit none
	Double Precision, intent(in)::t
	Double Precision, intent(inout)::Elec(:,:)
	Double Precision kx,ky,kz,x,y,tss,eee
	Double Precision tdt, t1, tau, rn
	Double Precision EadjustDrift,EadjustScat
	integer n
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph,phonon_count_ON
	Double Precision, intent(in)::DEph
	Double Precision eeePre							!	Energy conservation
	Double Precision kxPre, kyPre, kzPre			!	debug
	Double Precision kxPreD, kyPreD, kzPreD			!	debug
	Double Precision gk, eeeCalc, eeePreCalc		!	debug
	Double Precision eeePreD						!	debug
!
	Egetlost_drift=0.d0				!	Energy conservation 
	Egetlost_scat=0.d0				!	Energy conservation
!
	tdt=t+dt
	do  n=1,inum
		kx = Elec(n,PX)/hbar      ! Global kx
		ky = Elec(n,PY)/hbar      ! Global ky
		kz = Elec(n,PZ)/hbar      ! Global kz
		eee = Elec(n,EE)
		tss = Elec(n,TS)
		x  = Elec(n,XXX)
		y  = Elec(n,YYY)
!
		t1=t
		do while (tss <= tdt)          ! mean free time > elapsed time -> drift
			tau=tss-t1                   ! time interval before scattering
			eeePreD=eee										!	Energy conservation
			kxPreD=kx			!	debug
			kyPreD=ky			!	debug
			kzPreD=kz			!	debug
!
			call drift(tau,kx,ky,kz,eee,x,y,EadjustDrift)
!
			Egetlost_drift=Egetlost_drift+(eee-eeePreD)		!	Energy conservation
!
			gk = hbar**2*(kx**2+ky**2+kz**2)/(2.d0*mdos)
!			eeePreCalc=(sqrt(1.d0+4.d0*alpha*gk)-1.d0)/(2.d0*alpha)		!	eee after scat (Tomizawa (1.7))
			eeePreCalc=2.d0*gk/(sqrt(1.d0+4.d0*alpha*gk)+1.d0)				!	桁落ち防止？2018.12.08
!
			eeePre=eee			!	Energy conservation
			kxPre=kx			!	debug
			kyPre=ky			!	debug
			kzPre=kz			!	debug
!
			call scat(kx,ky,kz,eee,NELph0,NETph0,NBph,DEph,EadjustScat,phonon_count_ON)
!
			if (EadjustScat>0.5) then		!	************ Debug!	********************
				gk = hbar**2*(kx**2+ky**2+kz**2)/2.d0*mdos
!				eeeCalc=(sqrt(1.d0+4.d0*alpha*gk)-1.d0)/(2.d0*alpha)		!	eee after scat (Tomizawa (1.7))
				eeeCalc=2.d0*gk/(sqrt(1.d0+4.d0*alpha*gk)+1.d0)				!	桁落ち防止？2018.12.08
				write(*,*) '***************** Scat ****************'		!	debug
			end if
!																!	debug
			Egetlost_scat=Egetlost_scat+(eee-eeePre)		!	Energy conservation
			t1=tss
			rn = grnd()
			tss=t1-log(rn)/gm         ! re-new mean free time
		end do ! ---- while end
!
		tau=tdt-t1
		eeePre=eee											!	Energy conservation
!
		call drift(tau,kx,ky,kz,eee,x,y,EadjustDrift)
!
		Egetlost_drift=Egetlost_drift+(eee-eeePre)			!	Energy conservation
!
		Elec(n,PX) = hbar*kx    ! Global px
		Elec(n,PY) = hbar*ky    ! Global py
		Elec(n,PZ) = hbar*kz    ! Global pz
		Elec(n,EE) = eee
		Elec(n,TS) = tss
		Elec(n,XXX) = x
		Elec(n,YYY) = y
		Elec(n,EElost) = EadjustDrift+EadjustScat
!
	end do
end subroutine emc
!========================================================================================
!===( ドリフト過程の計算 )===
!
subroutine drift(tau,kx,ky,kz,eee,x,y,Elost)  ! Global kx, ky, kz
	implicit none
	Double Precision, intent(inout)::  tau, kx,ky,kz,eee,x,y
	Double Precision dkx, cp, sq, gk
	Double Precision padjust, Elost		!	2019.03.27
!
	dkx=qhbar*fx*tau ! dk=q/h・F・τ
!
	cp  = tau*hbar/mdos					!	h・τ/m
!
	gk = eee*(1.d0+alpha*eee)	!	gk before drift
	sq = sqrt(1.d0+4.d0*alpha*gk)      !--- anisotropic & non-parabolic
!
	x = x+cp*(kx+0.5d0*dkx)/sq   !--- dx=h/m・τ(kx +0.5・dkx)/sq
	y = y
	kx=kx+dkx							! Global kx
!
!	gk = (hbar*kx)**2/(2.d0*mdos) + (hbar*ky)**2/(2.d0*mdos) + (hbar*kz)**2/(2.d0*mdos)

	gk = hbar**2*(kx**2+ky**2+kz**2)/(2.d0*mdos)
	eee=2.d0*gk/(sqrt(1.d0+4.d0*alpha*gk)+1.d0)				!	桁落ち防止？2018.12.08　値は非負
!
!---( ドリフト後のエネルギーのチェック )---
!
		Elost = 0.d0
!
		if (eee < EMIN) then	!	2019.03.13 activate this correction
			write(*,*) '****  Warning in Drift ***  Smaller than EMIN: eee=', eee
			padjust = sqrt((EMIN*(1.+alpha*EMIN))/(eee*(1.+alpha*eee)))
			kx = kx*padjust
			ky = ky*padjust
			kz = kz*padjust
			Elost = Elost + (eee - EMIN)     !				   // add up "gotton" energy
			eee=EMIN
		end if
!
		if (eee > EMAX) then	!
			write(*,*) '****  Warning in Drift ***  Larger than EMAX: eee=', eee
			padjust = sqrt((EMAX*(1.+alpha*EMAX))/(eee*(1.+alpha*eee)))
			kx = kx*padjust
			ky = ky*padjust
			kz = kz*padjust
			Elost = Elost + (eee - EMAX)     !				   // add up "lost" energy
			eee=EMAX
		end if
!   
end subroutine drift
!========================================================================================
!===( 散乱過程の計算 )===
!
subroutine scat(kx,ky,kz,eee,NELph0,NETph0,NBph,DEph,Elost,phonon_count_ON)   ! Global kx, ky, kz
	implicit none
	Double Precision, intent(inout)::  kx,ky,kz,eee
	Double Precision ki, gk, ei	
	Double Precision r1
	integer ie
	Double Precision qM, kspmax, Epmax
	Double Precision aed
	integer LT
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph,phonon_count_ON
	Double Precision, intent(in)::DEph
	Double Precision padjust, Elost
!
!---( 散乱前のエネルギーの計算 )---
!
	gk = eee*(1.d0+alpha*eee)
	ki = sqrt(2.d0*mdos*gk)/hbar   ! Local k
!			ei=(sqrt(1.d0+4.d0*alpha*gk)-1.d0)/(2.d0*alpha)  ! Tomizawa (1.7)
			ei=2.d0*gk/(sqrt(1.d0+4.d0*alpha*gk)+1.d0)				!	桁落ち防止？2018.12.08　値は非負
	ei = eee
	ie=int(ei/de)+1
	if (ie > iemax) ie=iemax   ! denomination to EMAX ?
	if (ie < 1) then			! 2019.03.12
		write (*,*) ' ****** ei=',ei,' ie=', ie,' gk=',gk,' de=',de
		write (8,*) ' ****** ei=',ei,' ie=', ie,' gk=',gk,' de=',de
		ie = 1						! denomination to 1
	end if
!
!---( 散乱機構の選択 )---
!
	r1 = grnd()
	if (r1 <= swk(1,ie)) then       ! self-scattering
!			no operation
	else if (r1 <= swk(2,ie)) then   ! ion scattering
		call final_state_impurity_scat(kx,ky,kz,eee)
	elseif (r1 <= swk(6,ie)) then   !  LA/TA intra emission/absorption
		if (r1 <= swk(3,ie)) then   !  LA intra emission
			LT=LA
			aed=+1.d0
			qM=2.d0*ki
		elseif (r1 <= swk(4,ie)) then   !  LA intra absorption
			LT=LA
			aed=-1.d0
			Epmax = eee + hbar*get_wq(QMAX,LT)					!	max energy of electron at final state
			kspmax = sqrt(2.*mdos*Epmax*(1.+alpha*Epmax))/hbar	!	max momentum of electron at final state
			qM = ki + kspmax									!	max momentum of phonon at final state
!		elseif (r1 <= swk(5,ie)) then   !  TA intra emission
!			LT=TA
!			aed=+1.d0  
!			qM=2.d0*ki            
!		else   !    r1 <= swk(6,ie)     TA intra absorption
!			LT=TA
!			aed=-1.d0
!			Epmax = eee + hbar*get_wq(QMAX,LT)
!			kspmax = sqrt(2.d0*mdos*Epmax*(1.d0+alpha*Epmax))/hbar
!			qM = ki + kspmax									!	max momentum of phonon at final state
		end if
!
		call final_state_intra_scatByLATA(kx,ky,kz,ki,eee,LT,aed,qM,NELph0,NETph0,NBph,DEph,phonon_count_ON)
!
	elseif (r1 <= swk(7,ie)) then       !  LO-intra-emission
		LT=LO
		aed=+1.d0
!
		call final_state_intra_scatByLOTO(kx,ky,kz,ki,eee,LT,aed,NELph0,NETph0,NBph,DEph,phonon_count_ON)
!
	elseif (r1 <= swk(8,ie)) then       !  LO-intra-absorption
		LT=LO
		aed=-1.d0
!
		call final_state_intra_scatByLOTO(kx,ky,kz,ki,eee,LT,aed,NELph0,NETph0,NBph,DEph,phonon_count_ON)
!
	end if
!
!---( 散乱後のエネルギーのチェック )---
		Elost = 0.d0
!
		if (eee < EMIN) then	!	2019.03.13 activate this correction
			write(*,*) '****  Warning in Scat ***  Smaller than EMIN: eee=', eee
			padjust = sqrt((EMIN*(1.+alpha*EMIN))/(eee*(1.+alpha*eee)))
			kx = kx*padjust
			ky = ky*padjust
			kz = kz*padjust
			Elost = Elost + (eee - EMIN)     !				   // add up "gotton" energy
			eee=EMIN
		end if
!
		if (eee > EMAX) then	!
			write(*,*) '****  Warning in Scat ***  Larger than EMAX: eee=', eee
			padjust = sqrt((EMAX*(1.+alpha*EMAX))/(eee*(1.+alpha*eee)))
			kx = kx*padjust
			ky = ky*padjust
			kz = kz*padjust
			Elost = Elost + (eee - EMAX)     !				   // add up "lost" energy
			eee=EMAX
		end if
!      
end subroutine scat
!========================================================================================
subroutine final_state_intra_scatByLOTO(kx,ky,kz,ki,eee,LT,aed,NELph0,NETph0,NBph,DEph,phonon_count_ON)
	implicit none
	Double Precision, intent(inout):: kx,ky,kz,ki,eee
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::LT,NBph, phonon_count_ON
	Double Precision, intent(in)::aed,DEph
	Double Precision ef,kf, ephon, qphon,qphon0, gk, ph_qmin,ph_qmax,r1
	integer cnt,REDO,ie,cntMAX
	Double Precision f,cb,sb,fai,cf,sf,skk,a11,a12,a13,a21,a22,a23,a32,a33,x1,x2,x3
	Double Precision eee00,kx00,ky00,kz00
!
!--- memorize initial value in case of scattering cancel due to the final-state-finding failure
	eee00 = eee
	kx00 = kx
	ky00 = ky
	kz00 = kz
!--- iteration to find final energy & momentum size to meet conservation rule
	cntMAX=1000
	cnt=1
	REDO=1
	r1 = grnd()		!	uniform random number
					!	scattering direction probability given => find p,e by iteration to satisfy conservation
	qphon0 = QmidBZ*QMAX		!	BZ原点と1stBZ端との中間値を初期値に
	do while (REDO > 0)
		ephon = hbar*get_wq(qphon0,LT)
		ef = eee - aed*ephon
!		if ((aed > 0.) .AND. (ef <= EMIN)) then
		if (ef <= EMIN) then		!	2019.03.14
!			write (*,*) '**** Warning No negative electron energy ex. due to emission in LOTO'
!		Cancel scattering & do nothing
			eee = eee00
			kx = kx00
			ky = ky00
			kz = kz00
			return
		end if
		gk = ef*(1.d0+alpha*ef)
		kf = sqrt(2.d0*mdos*gk)/hbar
		ph_qmax = ki + kf
		ph_qmin = Abs(ki - kf)
!!!		r1 = grnd()
		qphon=(ph_qmax**r1)*(ph_qmin**(1.-r1))	!	Transform uniform random number To
												!	-> qphon(phonon size~scattering direction~rate to final state)
		if (abs((qphon-qphon0)/(qphon+qphon0)) < epsLL) then	!	convergence check
			REDO=0
		end if
		if (cnt > cntMAX) then
			write(*,*) '******** Warning cnt > cntMAX in LOTO'
			REDO=0
		end if
		qphon0=qphon
		cnt=cnt+1
	end do
!		**  phonon counting(+emit & -absorb)
	if (phonon_count_ON == 1 .AND. cnt <= cntMAX) then
		ie=int(ephon/DEph)+1
		if (ie > NBph) then
			write(*,*) '** WARNING 8 ** too large ephon',ephon
			write(8,*) '** WARNING 8 ** too large ephon',ephon
		else if (LT==LA .OR. LT==LO) then
			NELph0(ie)=NELph0(ie)+aed*1.d0
		else if (LT==TA .OR. LT==TO) then
			NETph0(ie)=NETph0(ie)+aed*1.d0
		end if
	end if
!--- calculation of final energy & momentum vector
	kf  = (sqrt(2.0d0*mdos)/hbar)*sqrt(ef*(1.0d0+alpha*ef))
	f   = 2.0d0*ki*kf/(ki-kf)**2.0d0
	if (f <= 0.0d0) return
	cb  = (1.0d0+f-(1.0d0+2.0d0*f)**r1)/f		!	cosine of scattering angle; cb should be [-1 1];
	if (cb > 1.0d0) then
		write(*,*)'cb in intra_ScatByLOTO; cb, f, r1, kf, ki',cb,f,r1,kf,ki	!	2019.04.02
		write(8,*)'cb in intra_ScatByLOTO; cb, f, r1, kf, ki',cb,f,r1,kf,ki	!	2019.04.02
		cb=1.0d0
	end if
	if (cb < -1.0d0) then
		write(*,*)'cb in intra_ScatByLOTO; cb, f, r1, kf, ki',cb,f,r1,kf,ki	!	2019.04.02
		write(8,*)'cb in intra_ScatByLOTO; cb, f, r1, kf, ki',cb,f,r1,kf,ki	!	2019.04.02
		cb=-1.0d0
	end if
	sb  = sqrt(1.0-cb*cb)
	fai = 2.0*pi*grnd()
	cf  = cos(fai)
	sf  = sin(fai)
	skk = sqrt(kx*kx+ky*ky)
	a11 = ky/skk
	a12 = kx*kz/skk/ki
	a13 = kx/ki
	a21 =-kx/skk
	a22 = ky*kz/skk/ki
	a23 = ky/ki
	a32 =-skk/ki
	a33 = kz/ki
	x1  = kf*sb*cf
	x2  = kf*sb*sf
	x3  = kf*cb
	kx  = a11*x1+a12*x2+a13*x3
	ky  = a21*x1+a22*x2+a23*x3
	kz  =        a32*x2+a33*x3
	eee  = ef
!
	if(cnt > cntMAX) then	!	cancel scattering
		eee = eee00
		kx = kx00
		ky = ky00
		kz = kz00
	end if
!
	return
end subroutine final_state_intra_scatByLOTO
!========================================================================================
subroutine final_state_intra_scatByLATA(kx,ky,kz,ki,eee,LT,aed,qM,NELph0,NETph0,NBph,DEph,phonon_count_ON)
! common treatment for intra_ScatByLA&TA:
	implicit none
	Double Precision, intent(inout)::  kx,ky,kz,eee
	Double Precision ki
	integer LT
	Double Precision aed
	Double Precision C, qM, q0, f0, fx, ephon, acost
	Double Precision pprime, kprime
	Double Precision cbet, sbet, cgam, sgam, cgamp, cgampp, eta, ceta, seta, phi, sphi, cphi
	integer cnt
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph, phonon_count_ON
	Double Precision, intent(in)::DEph
	integer ie
	Double Precision gk,gk0,eeeCalc,kix,kiy,kiz,ramda,kx0,ky0,kz0,ramda0		!	trial
	integer REDO,REDO_Outer_Loop												!	trial
	Double Precision eee00,kx00,ky00,kz00
!
	eee00 = eee
	kx00 = kx
	ky00 = ky
	kz00 = kz
!
	C = qM*get_Mq(qM,LT,aed)*(1.d0+2.d0*alpha*(eee-aed*hbar*get_wq(qM,LT))) 
!
	q0 = QMIN + (qM-QMIN)*grnd()				!	set q0(up phonon momentum) b/w QMIN~qM
	ephon = hbar*get_wq(q0,LT)					!	phonon energy
	f0 = C*grnd()
	fx = q0*get_Mq(q0,LT,aed)*(1.d0+2.d0*alpha*(eee-aed*ephon))
	acost = fabs_cost(q0,ki,eee,LT,aed)
	cnt=1										!	debug
	do while ((f0 > fx) .OR. (acost > 1.d0))	!	M.C.procedure to choose q0 in accordance with the rate to final state
		q0 = QMIN + (qM-QMIN)*grnd()			!	setq0(up phonon momentum) b/w QMIN~qM
		ephon = hbar*get_wq(q0,LT)				!	phonon energy
		f0 = C*grnd()
		fx = q0*get_Mq(q0,LT,aed)*(1.d0+2.d0*alpha*(eee-aed*ephon))
		acost = fabs_cost(q0,ki,eee,LT,aed)
		cnt=cnt+1								! debug
		if (mod(cnt,100000)==0) then			! debug
			write(*,'(A,i10,2(A,E22.15),A,i10,A,E22.15)') 'cnt=',cnt,'eee=',eee,'acost=',acost,'LT=',LT,'aed=',aed    ! debug
!			write(8,'(A,i10,2(A,E22.15),A,i10,A,E22.15)') 'cnt=',cnt,'eee=',eee,'acost=',acost,'LT=',LT,'aed=',aed    ! debug
		end if
	end do
!               
!!	set the electron energy, consistent w/the momentum exchange above
	eee = eee - aed*ephon
	gk0 = eee*(1.d0+alpha*eee)
!	if ((aed > 0.) .AND. (ef <= EMIN)) then
	if (eee <= EMIN) then		!	2019.03.14
!		write (*,*) '**** Warning No negative electron energy ex. due to emission in LATA'
!	Cancel scattering & do nothing
		eee = eee00
		kx = kx00
		ky = ky00
		kz = kz00
		return
	end if
!
!        
!--- calculation of final state for LA/TA phonon emission/absorption
!	運動量空間で中心が(kix,kiy,kiz)、半径q0の球面と、
!	原点中心、半径が電子終状態エネルギー対応の回転楕円体表面点との、交線上の任意の点を繰り返しによって探す
!
	kix=kx			!	Global → Local
	kiy=ky			!	Global → Local
	kiz=kz			!	Global → Local
	pprime = sqrt(2.d0*mdos*eee*(1.d0+alpha*eee)) !
	kprime = pprime/hbar  ! k=p/hbar
	cbet = (ki*ki+kprime*kprime-q0*q0)/(2.d0*ki*kprime) ! low of cosine angle beta b/w ki and kprime shold be [-1 1]
!
	if (abs(cbet) > 1.0) then
		write(*,*) 'error abs(cbet) > 1.0'
		write(*,'(4(A,E22.15))') 'ki=',ki,' kprime=',kprime,' q0=',q0,' cbet=',cbet    !debug
!		write(8,'(4(A,E22.15))') 'ki=',ki,' kprime=',kprime,' q0=',q0,' cbet=',cbet    !debug
!		 	 exit(-1)
		stop
	end if
	sbet = sqrt(1.d0-cbet*cbet)		!	/*** gam = angle b/w ki and Z-axis ***/
	cgam = hbar*kz/(hbar*ki)
	if (abs(cgam) > 1.0) then
		write(*,*) 'error abs(cgam) > 1.0'
		write(*,'(3(A,E22.15))') 'cbet=',cbet,' sbet=',sbet,' cgam=',cgam    !debug
!		write(8,'(3(A,E22.15))') 'cbet=',cbet,' sbet=',sbet,' cgam=',cgam    !debug
!	 			 exit(-1)
		stop
	end if
	sgam = sqrt(1.d0-cgam*cgam)		!	/*** gamp = angle b/w ki and Y-axis ***/
	cgamp = ky/ki		!	/*** gampp = angle b/w ki and X-axis ***/
	cgampp = kx/ki		!	/*** angle b/w ki projection on X-Y plane and Y-axis ***/
	eta = atan(kx/ky)
	ceta = cos(eta)
	seta = sin(eta)
!
	REDO_Outer_Loop=100
	do while (REDO_Outer_Loop > 0)
		phi = 2.d0*pi*grnd()
		sphi = sin(phi) 
		cphi = cos(phi) 
!		--- first candidate k
		kz = (cbet*cgam   +sbet*cphi*sgam)*pprime/hbar						!	Local kz
		ky = (cbet*cgamp  -sbet*cphi*cgam*ceta -sbet*sphi*seta)*pprime/hbar	!	Local ky
		kx = (cbet*cgampp -sbet*cphi*cgam*seta +sbet*sphi*ceta)*pprime/hbar	!	Local kx
!
		cnt=1
		REDO=1
		gk=0.5d0*hbar*hbar*(kx**2+ky**2+kz**2)/mdos
!		eeeCalc=(SQRT(1.d0+4.d0*alpha*gk)-1.d0)/(2.d0*alpha)
		eeeCalc=2.d0*gk/(sqrt(1.d0+4.d0*alpha*gk)+1.d0)				!	桁落ち防止？2018.12.08
		if (ABS(eee-eeeCalc)/eee < eee_convergence) then
			REDO=0
			REDO_Outer_Loop=0
		end if
!
		do while(REDO > 0)
!			(kx,ky,kz)から最短距離の、(kix,kiy,kiz)を中心とし半径q0の球表面点→新たに(kx,ky,kz)
			ramda=q0/SQRT((kx - kix)**2+(ky - kiy)**2+(kz - kiz)**2)
			kx=kix + ramda*(kx - kix)
			ky=kiy + ramda*(ky - kiy)
			kz=kiz + ramda*(kz - kiz)
			kx0=kx							!	debug for print
			ky0=ky							!	debug for print
			kz0=kz							!	debug for print
			ramda0=ramda					!	debug for print
			gk=0.5d0*hbar*hbar*(kx**2+ky**2+kz**2)/mdos
!			eeeCalc=(SQRT(1.d0+4.d0*alpha*gk)-1.d0)/(2.d0*alpha)
			eeeCalc=2.d0*gk/(sqrt(1.d0+4.d0*alpha*gk)+1.d0)				!	桁落ち防止？2018.12.08
			if (ABS(eee-eeeCalc)/eee < eee_convergence) then
				REDO=0
				REDO_Outer_Loop=0
			end if
!
!			(kx,ky,kz)から～最短距離の、原点を中心とした回転楕円体表面点（球近似）→新たに(kx,ky,kz)
!			ramda=SQRT(gk0/gk)
			ramda=((kx**2+ky**2+kz**2)/mdos-2.*gk0/hbar**2)/(2.*(kx**2+ky**2+kz**2)/mdos**2)
			kx=kx/(1+ramda/mdos)
			ky=ky/(1+ramda/mdos)
			kz=kz/(1+ramda/mdos)
!
			if (cnt > 100) then
				REDO=0
				REDO_Outer_Loop=REDO_Outer_Loop-1
				if(REDO_Outer_Loop==0) then			!	Give Up ==>> No scattering
					REDO_Outer_Loop = -1
					write(*,'(A,I10)',advance = 'no') ' ** WARNING 2 (intra_scatByLATA) **  REDO_Outer_Loop= ',REDO_Outer_Loop	!	debug
					write(*,'(6(A,E22.15))') ' ',ramda,' ',kx/QMAX,' ',ky/QMAX,' ',kz/QMAX,' ',gk,' ',eee	!	debug
				end if
			end if
			cnt=cnt+1
!
		end do	!	while(REDO)
	end do	!	while(REDO_Outer_Loop)
!									!	Global k
!
	if(REDO_Outer_Loop == -1) then
		eee = eee00
		kx = kx00
		ky = ky00
		kz = kz00
	end if
!
!		**  phonon counting(+emit & -absorb)
	if (phonon_count_ON == 1 .AND. REDO_Outer_Loop >= 0) then
		ie=int(ephon/DEph)+1
		if (ie > NBph) then
			write(*,*) '** WARNING 1 ** too large ephon',ephon
!			write(8,*) '** WARNING 1 ** too large ephon',ephon
		else if (LT==LA .OR. LT==LO) then
			NELph0(ie)=NELph0(ie)+aed*1.d0
		else if (LT==TA .OR. LT==TO) then
			NETph0(ie)=NETph0(ie)+aed*1.d0
		end if
	end if
!
!--- end calculation of final state for LA/TA phonon emission/absorption (type=2-5)
end subroutine final_state_intra_ScatByLATA
!========================================================================================
!!subroutine final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph,phonon_count_ON) ! Global kx, ky, kz
!!end subroutine final_state_inter_scat_g
!========================================================================================
!!subroutine final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph,phonon_count_ON) ! Global kx, ky, kz
!!end subroutine final_state_inter_scat_f
!========================================================================================
!!subroutine change_equiv_valley(ivv)
!!end subroutine change_equiv_valley
!========================================================================================
!!subroutine change_noneq_valley(ivv)
!!end subroutine change_noneq_valley
!========================================================================================
subroutine final_state_impurity_scat(kx,ky,kz,eee)
	implicit none
	Double Precision, intent(inout)::  kx,ky,kz,eee
	Double Precision pprime,calp,salp,bet
!
	pprime = sqrt(2.d0*mdos*eee*(1.d0+alpha*eee)) !
!
	calp = 1.d0-2.d0*grnd()
	salp = sqrt(1.d0-calp*calp)
	bet = 2.d0*pi*grnd()
	kx = pprime*salp*cos(bet)/hbar
	ky = pprime*salp*sin(bet)/hbar
	kz = pprime*calp/hbar
!
end subroutine final_state_impurity_scat
!
!========================================================================================
!===( 出力ルーチン )===
!
subroutine out(t,Elec)
	implicit none
	Double Precision, intent(in)::t
	Double Precision, intent(in)::Elec(:,:)
	Double Precision ve, s, kx, gk
	Double Precision eee                       ! debug
	Double Precision eeeLost,eeeMin
	integer n
!
!---( 平均速度と距離とエネルギー保存の計算 )---
!
	ve = 0.d0
	eee = 0.d0   ! debug
	s = 0.d0
	eeeLost=0.d0
	eeeMin=999.d0
	do  n=1,inum
!!		ivv  = iv(n)
		gk = (Elec(n,PX)**2+Elec(n,PY)**2+Elec(n,PZ)**2)/(2.d0*mdos)
		kx = Elec(n,PX)/hbar   ! Local kx
		ve = ve + hbar*kx/mdos/sqrt(1.d0+4.d0*alpha*gk) 
		eee = eee + Elec(n,EE)
		if (Elec(n,EE) < eeeMin) eeeMin=Elec(n,EE)
		s  = s + Elec(n,XXX)
		eeeLost=eeeLost+Elec(n,EElost)
	end do
	ve = ve/float(inum)
	eee = eee/float(inum)
	s  = s /float(inum)
!
!---( 時間 vs 平均速度 vs 距離の出力 )---
	write(*,'(8(E22.15,A))') t,' ',s,' ',ve,' ',eee,' ',eeeLost,' ',eeeMin,' ',Egetlost_scat,' ',Egetlost_drift
	write(8,'(8(E22.15,A))') t,' ',s,' ',ve,' ',eee,' ',eeeLost,' ',eeeMin,' ',Egetlost_scat,' ',Egetlost_drift
end subroutine out
!  
!========================================================================================
!===( エネルギー分布計算 )===
!
subroutine energy_dist(Elec)
	implicit none
	Double Precision, intent(in)::Elec(:,:)                 !
	integer n
!
	Double Precision, allocatable ::ehist(:)
	integer:: nhist=301                !ヒストグラム分割数
	integer digitized_energy, iemax
	Double Precision eemax
!
	allocate(ehist(nhist))
!      
!---( 個々の粒子のエネルギー計算 )---
!
	eemax=0.d0                 ! do not use ::emax=0.d0
	do  n=1,inum
!
		if (Elec(n,EE) > eemax) then
			iemax = n              ! the same particle?
			eemax = Elec(n,EE)
		endif
	end do
!
	do  n=1,nhist
		ehist(n)=0.0
	end do
!
!			eemax = 0.01d0   ! Debug fine division  10meV/300div -> 0.03meV/div
	do  n=1,inum
		digitized_energy=int(real(nhist-1)*(Elec(n,EE)/eemax))+1
		if (digitized_energy <= nhist) then
			ehist(digitized_energy)=ehist(digitized_energy)+1.0
		end if
	end do   
!
!	--- histogram area normalization
	do  n=1,nhist
		ehist(n)=ehist(n)/eemax
	end do   
!
!---( エネルギーヒストグラムの出力 )---
!
	write(*,*) '* Emax=',eemax, '  particle#', iemax ! debug
	write(8,*) '* Emax=',eemax, '  particle#', iemax ! debug	!	blank is important for after excel treatment
!
	do n=1,nhist
		write(*,'(A,E22.15,A,E22.15)') ' * ',real(n)/real(nhist)*eemax,' ',ehist(n)
		write(8,'(A,E22.15,A,E22.15)') ' * ',real(n)/real(nhist)*eemax,' ',ehist(n)	!	blank is important for after excel treatment
	end do       
!		deallocate(eee)
	deallocate(ehist)
!
end subroutine energy_dist  
!
end module Tom_sub_programs  

!****************************************************************************************
program main
!
	use Tom_sub_programs
	use Pop_sub_programs
	implicit none
	Double Precision t ,SimTime
	Double Precision,allocatable ::Elec(:,:)
	integer inum_dummy, i
!!!	integer, parameter::NBph=280	!	2018.10.29 280->1000
!!!	integer, parameter::NBph=160	!	2019.02.20 280->160	Si->GaAs
	integer, parameter::NBph=200	!	2019.02.27 160->200	GaAs->GaN
!!!	Double Precision, parameter::DEph=42.d-3/float(NBph)	!	1.5e-4  ! 42meV/280div	Si
!!!	Double Precision, parameter::DEph=40.d-3/float(NBph)	!	2.5e-4  ! 40meV/160div	GaAs
	Double Precision, parameter::DEph=100.d-3/float(NBph)	!	5e-4  	! 100meV/200div	GaN
	Double Precision, allocatable ::NELph0(:), NETph0(:)
	Double Precision t1, t2
!
	call cpu_time( t1 )
!---( ファイルのオープン )---
!
	open(5,file='data2')
	open(8,file='outf')
!
!---( 計算条件：総粒子数, 総時間ステップ, 時間増分, 格子温度, 電界 )---
!
	call data_input(inum_dummy,SimTime)
	allocate(Elec(inum_dummy,10))
	allocate(NELph0(NBph), NETph0(NBph))
!
!---( 物理定数, 材料定数および諸パラメータ )---
!
	call param
!
!---( 初期条件の設定 )---
!
	call initia(t,Elec)  ! return with t=0.
!	call energy_dist(Elec)                    !  初期エネルギー分布出力
!		call out(t,Elec,iv)                         !  debug
	do i=1,NBph
		NELph0(i)=0.d0
		NETph0(i)=0.d0
	end do
!
!---( 多粒子モンテカルロ計算 )---
!
	call monte_carlo(Elec,NELph0,NETph0,NBph,DEph)
!
!---( 粒子の最終エネルギー分布 )---
	call energy_dist(Elec)
!
!---( フォノンの頻度分布 )---
	do i=1,NBph
		NELph0(i)=NELph0(i)*(1.0d17/inum_dummy)*(1.0/(1.0-phonon_count_ON_ratio)/SimTime)		!	phonon_count_ON time correction
		NETph0(i)=NETph0(i)*(1.0d17/inum_dummy)*(1.0/(1.0-phonon_count_ON_ratio)/SimTime)		!	phonon_count_ON time correction
		write (*,'(i10,3(A,E22.15))') i,' ',dble(i)*DEph,' ',NELph0(i),' ',NETph0(i)
		write (8,'(i10,3(A,E22.15))') i,' ',dble(i)*DEph,' ',NELph0(i),' ',NETph0(i)
	end do
!      
!---( ファイルのクローズ )---
!
	call cpu_time( t2 )
	write(*,*) 'cpu time:', t2-t1, 'seconds.'
	write(8,*) 'cpu time:', t2-t1, 'seconds.'
	deallocate(Elec)
	close(5)
	close(8)
	stop
end program main
