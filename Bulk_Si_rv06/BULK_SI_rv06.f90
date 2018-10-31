!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     バルクシリコン内の電子輸送のための
!     多粒子モンテカルロプログラム
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! rv01 2018.08.11 平均速度計算のバグ？をコメント化
! rv02 2018.08.17 rnd() -> random_number
! rv03 2018.08.22 c -> !
!                 do statementsを do~~end do型にし文番号削除
!                 散乱による谷遷移部をサブルーチンに
!                 散乱後の終状態計算をサブルーチンに
!                 文番号を完全に削除
!                 common変数のうち変化する変数は引数化
!                 最初に定義してその後は参照するだけのcommon変数はmoduleの共通変数化
!                 ファイルの前半にsubprogramを含んだmodule。最後にmain program
! rv03 2018.08.23
!                 繰り返し最終（～定常）状態の粒子エネルギー分布出力サブルーチン追加
! rv04 2018.08.23
!                 Mersenne Twister版の乱数　mt19937
! compile with option as >gfortran -fno-range-check Bulk_Si_rv04.f90
!
! rv05 2018.08.30
!                 Tom -> Pop (scattering type 5->18; include phonon dispersion)
!                 energy div 2->0.5meV; range 2->1.2eV
!
!      2018.08.31
!                 scattering type (1,2,3,4,5,self -> 6,7,8,9,5,1;2-4=0)
!      2018.09.20
!                 scat. type #5->#2,3,4,5(LA&TA x Emiss&Abs)...>smaller scat.rate/larger μ
!                 SI unit -> eV & cm unit(Pop)
!                 h -> hbar
!      2018.09.20(2)
!                 ee->eee,ts->tss,p(i,:) -> Elec(i,PX..);PX=5,PY=6,PZ=7,TS=9,XXX=1,YYY=2
!                 iv->ivv,ip(i)->iv(i); tentatively iv=1(+)&2(-)<-ip=1,iv=3(+)&4(-)<-ip=2,iv=5(+)&6(-)<-ip=3
!      2018.10.15
!					All Tom scattering processes were replaced by Pop processes (#1-17 except ion scat)
!	   2018.10.24
!					Elec(i,EElost=10)
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
	Double Precision,parameter::EMIN = 0.d-10	! Pop 1.d-10 << 0.0005(=de)
	Double Precision,parameter::EMAX = 1.5D0	! 1.5eV <2eV
	Double Precision,parameter::QMIN = 1.d0	! // nearly zero compared to QMAX
	Double Precision,parameter::QMAX = 1.1569113068d+08		!  2*pi/asi
!				#define asi    5.431d-08        // lattice constant, cm 
!!	Double Precision,parameter::QMAX = 1.15691130546861d+08		!	just Pop_C value
!    Double Precision,parameter::asi = 5.431d-08         ! // lattice constant, cm
	Double Precision,parameter::Rws = 2.1224d-8         !  asi*pow(3./(16.*pi), 1./3.)
!                                           Egap = 1.1756 - 8.8131e-05*T - 2.6814e-07*T*T
!!	Double Precision,parameter::mdos = 1.86314779783618D-16      ! pow(ml*mt*mt, 1./3.) eV*s^2/cm^2
!     ml = 0.9163*m0               // 5.2097d-16 eV*s^2/cm^2
!     mt = 0.1880*m0*1.1756/Egap   // 1.1220d-16 eV*s^2/cm^2
!				#define m0     5.68562975d-16      // eV*s^2/cm^2
!
	Double Precision,parameter::epsSi = 6.46582917d+06	!	// e/V/cm permittivity of Si 11.7*eps0
!				#define eps0   5.52634972d+05      // e/V/cm permittivity of free space
!
	integer,parameter:: LA=1
	integer,parameter:: TA=2
    integer,parameter:: LO=3
    integer,parameter:: TO=4
!
    Double Precision,parameter::pi = 3.14159265d0        !π
    Double Precision,parameter::rho = 1.4516d+12       ! density (eV*s^2/cm^5) (2.329 g/cm^3)
    Double Precision,parameter::kb = 8.61734318d-05    ! eV/K
    Double Precision,parameter::hbar = 6.58211915d-16  ! eV*s
!
    Double Precision,parameter::DLA = 6.39d0   ! debug 13.0= 6.39 x ~2
    Double Precision,parameter::DTA = 3.01d0   ! debug 6.0=3.01 x ~2
    Double Precision,parameter::Absorption = -1.d0
    Double Precision,parameter::EMS = 1.d0
!			  // f- and g-phonon norm, to be multiplied by QMAX
    Double Precision, parameter:: QGp = 0.3d0           ! @ Pop-default
	Double Precision, parameter:: QFp = 1.022252415d0 ! sqrt(1. + 2.*0.15*0.15)
!	Double Precision, parameter:: QFp = 1.022252d0	!	just Pop_C value
    Double Precision,parameter::PQMAX = 7.6149280673d-08 ! hbar*QMAX(Wigner-Seitz cell radius)
    
	Double Precision,parameter::LTAO(4,3) = &
						&  reshape( (/-0.196d-2, -0.250d-2, -0.156d-2,     0.137d-2, &
						&              9.00d5,    5.33d5,    0.0000d0,      -2.91d5,   &
						&              0.0d0,       0.0d0,       9.87682d13,   1.027d14 /),(/4,3/) )
!	Double Precision,parameter::LTAO(4,3) = &
!						&  reshape( (/-0.200d-2, -0.260d-2, -0.160d-2,     0.111d-2, &
!						&              9.01d5,    5.23d5,    0.0000d0,      -2.57d5,   &
!						&              0.0d0,       0.0d0,       9.88d13,   1.02d14 /),(/4,3/) )
!	Double Precision,parameter::LTAO(4,3) = &
!						&  reshape( (/ 0.0d0,       0.0d0, 		-0.156d-2,     0.137d-2, &
!						&              9.00d5,    5.33d5,    0.0000d0,      -2.91d5,   &
!						&              0.0d0,       0.0d0,       9.87682d13,   1.027d14 /),(/4,3/) )
!
!		from Tom_global_variables   !----- 共通変数 -----
	Double Precision dt,Temp,fx
	Double Precision de
	Double Precision :: swk(18,3002)=0.d0
	Double Precision mdos      
	Double Precision gm                 
	Double Precision smh,hhml,hhmt         
	Double Precision alpha
	Double Precision tm(3)                
	Double Precision hm(3)                
	Double Precision qhbar                   
	Double Precision kbTq                  
	Double Precision ef
	integer inum,jtl
	integer iemax
	integer, parameter::PX=5,PY=6,PZ=7,EE=8,TS=9,XXX=1,YYY=2
	integer, parameter::EElost=10
	Double Precision pvx(6), pvy(6), pvz(6)  ! 6 valley p offset
	Double Precision mx(6), my(6), mz(6)     ! mx,my,mz at 6 valleys
!
	Double Precision Doping
	Double Precision Egetlost_scat, Egetlost_drift		!	Energy conservation
	
!
end module Tom2Pop_sub_program_variables
!****************************************************************************************
module Pop_sub_programs
	use Tom2Pop_sub_program_variables
!	use M_random_number
	use N_random_number
	implicit none
	contains
!--------------- 以下に次のsubroutine/function ------------
!	trapz 		! 	数値積分の足し算
!	get_wq 		!	phonon q -> w
!	get_nq 		! 	phonon occupation
!	get_Mq 		! 	(3.8) 積分式の被積分関数
!	fabs_cost  	! 	(3.12) cos(φ)式の絶対値
!	find_root  	! 	|cos(φ)|-1=0の根を探索
!	acrate3    	!	E毎に(3.8) 積分
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
			write(*,*) 'gasdev while cnt=', cnt				!! debug
			write(8,*) 'gasdev while cnt=', cnt				!! debug
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
			Return			! 2018.10.24
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
		beta = sqrt(DOPING/(epsSi*kb*Temp))
		gamma = Ed*(1.+alpha*Ed)
		ks = sqrt(2.d0*mdos*gamma)/hbar
		v0 = hbar*ks/(mdos*(1.d0+2.d0*alpha*Ed))
		xx = hbar*hbar*(2.d0*ks)*(2.d0*ks)/(8.d0*mdos*kb*Temp)
		Gx = (1.d0+0.165651d0*xx**2.+0.0557896d0*xx**4.d0)/ 	&
		&	(1.d0+0.829201d0*xx**2.d0+0.343064d0*xx**4.d0+0.111354d0*xx**6d0)
		bb = 4.d0*ks*ks/(beta*beta*Gx)
		rateDOP = DOPING/(2.d0*pi*hbar*hbar*epsSi*epsSi*v0*4.d0*ks*ks)
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
!	use M_random_number
	use N_random_number
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
	inum_dummy=inum    ! for allocation of Elec(inum,9),iv(inum)
	read(5,*) jtl      ! number of iteration (typ. 1000~10000)
	read(5,*) dt       ! time step [s] (typ. 5E-15 = 5fs)
	SimTime=jtl*dt		! Total Simulation Time [s] 
	read(5,*) Temp      ! temperature [K] (typ.300)
	read(5,*) fx       ! electric field [KV/cm] 
	read(5,*) Doping	!	Dopant Concentration [1/cm3]
!
	write(*,*) inum/1D6,' M particles'
	write(*,*) jtl, ' iteration'
	write(*,*) dt, ' sec'
	write(*,*) Temp, ' K'
	write(*,*) fx, ' KV/cm'
	write(*,*) Doping, ' /cm3'
!
	write(8,*) inum/1D6,' M particles'
	write(8,*) jtl, ' iteration'
	write(8,*) dt, ' sec'
	write(8,*) Temp, ' K'
	write(8,*) fx, ' KV/cm'
	write(8,*) Doping, ' /cm3'
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
	Double Precision ml, mt, eps, Egap
	Double Precision amc      
	Double Precision ei
	Double Precision eee
	integer ie, i
!
!!	Double Precision,parameter::pi = 3.14159265        !π
!
	Double Precision,parameter::q = 1.d0               ! e
	Double Precision,parameter::echarge = 1.d0		!  // electron charge 
!				#define echarge   -1.0             // electron charge 
!				#define ecoulom   -1.60217653d-19  // electron charge in Coulombs
	Double Precision,parameter::hbar  = 6.58211915d-16  ! eV*s 
	Double Precision,parameter::ep0 = 5.52634972d+05    ! e/V/cm    !真空誘電率
	Double Precision,parameter::m0 = 5.68562975d-16    ! eV*s^2/cm^2  !electron mass
!      
!!!			remind LT=   LA=1, TA=2, LO=3, TO=4
!		 aed = -1. for absorption,  +1. for emission
!!	Double Precision,parameter::PQMAX = 7.6149280673d-08 ! hbar*QMAX(Wigner-Seitz cell radius)
!  /************************************************************************/
!  /**************** INTERVALLEY SCATTERING (LO+TO+LA+TA) ******************/
!  /************************************************************************/
	Double Precision, parameter :: Zf=4.d0         ! // f-scat degeneracy of final valley
	Double Precision, parameter :: Zg=1.d0         ! // g-scat degeneracy of final valley
	Double Precision, parameter :: DTAf = 0.6d8   ! // 18 meV
	Double Precision, parameter :: DLAf = 3.5d8   ! // 50 meV
	Double Precision, parameter :: DTOf = 1.5d8   ! // 57 meV
	Double Precision, parameter :: DTAg = 0.7d8   ! // 10 meV
	Double Precision, parameter :: DLAg = 1.d8   ! // 19 meV
	Double Precision, parameter :: DLOg = 6.d8   ! // 63 meV
	Double Precision scatOconst, Nq
	Double Precision ETAf, ELAf, ETOf, ETAg, ELAg, ELOg
!    
!---( 電子の縦有効質量, 横有効質量, 禁制帯幅, 誘電率 )---
!
	Egap = 1.1756 - 8.8131d-05*Temp - 2.6814d-07*Temp*Temp
!	write(*,*) 'Egap_Pop=',Egap					! debug
!	write(8,*) 'Egap_Pop=',Egap					! debug
	ml = 0.9163d0*m0					!// 5.2097e-16 eV*s^2/cm^2
	mt = 0.1880d0*m0*1.1756d0/Egap		!// 1.1220e-16 eV*s^2/cm^2
	eps  = 11.7d0*ep0  ! Si誘電率
!
!---( フォノン散乱の諸パラメータ )---
	kbTq = kb*Temp/echarge
	qhbar   = echarge/hbar     ! h->hbar
!
!---( バンドの非等方性 )---
!
	mdos  = (ml*mt*mt)**(1.d0/3.d0) !相乗平均 effective m 0.328?
	amc  = 3.d0/(1.d0/ml+2.d0/mt)!逆数相加平均 effective m at B edge 0.266?
	tm(1)=sqrt(ml/mdos) !縦 1.67
	tm(2)=sqrt(mt/mdos) !横 0.773
	tm(3)=sqrt(mt/mdos) !横 0.773
	smh  = sqrt(2.d0*mdos*echarge)/hbar    ! h->hbar
	hhml = hbar/ml/echarge*hbar/2.d0 !縦 ! h->hbar
	hhmt = hbar/mt/echarge*hbar/2.d0 !横 ! h->hbar
	hm(1)= hbar/ml !縦 ! h->hbar
	hm(2)= hbar/mt !横 ! h->hbar
	hm(3)= hbar/mt !横 ! h->hbar
!
!---( バンドの非放物線性 )---
!
	alpha = 0.5d0*1.1250281d0/Egap
!	write(*,*) 'alpha_Pop=',alpha				! debug
!	write(8,*) 'alpha_Pop=',alpha				! debug
!
!---( valley offset values )---
!
	do i=1,6
		pvx(i) = 0.d0
        pvy(i) = 0.d0
        pvz(i) = 0.d0
	end do
	pvx(1) = -0.85d0*PQMAX  ! -X axis valley
	mx(1) = ml
	my(1) = mt
	mz(1) = mt
	pvx(2) = +0.85d0*PQMAX  ! +X axis valley
	mx(2) = ml
	my(2) = mt
	mz(2) = mt
	pvy(3) = -0.85d0*PQMAX  ! -Y axis valley
	mx(3) = mt
	my(3) = ml
	mz(3) = mt
	pvy(4) = +0.85d0*PQMAX  ! +Y axis valley
	mx(4) = mt
	my(4) = ml
	mz(4) = mt
	pvz(5) = -0.85d0*PQMAX  ! -Z axis valley
	mx(5) = mt
	my(5) = mt
	mz(5) = ml
	pvz(6) = +0.85d0*PQMAX  ! +Z axis valley
	mx(6) = mt
	my(6) = mt
	mz(6) = ml
!
!---( 散乱レートの計算 )---
!
!					 // Fphon = 2-norm([1 .15 .15]) in constants.h
	ETAf = hbar*get_wq(QFp*QMAX,TA)     ! // ~ 18 meV
!					 // this one is actually on the LO branch here b/c Fphon > 1
	ELAf = hbar*get_wq(QFp*QMAX,LA)     ! // ~ 52 meV
!					  // see how get_wq() is implemented for LA when q > QMAX
	ETOf = hbar*get_wq(QFp*QMAX,TO)     ! // ~ 58 meV
	ETAg = hbar*get_wq(QGp*QMAX,TA)     ! // ~ 10 meV
	ELAg = hbar*get_wq(QGp*QMAX,LA)     ! // ~ 19 meV
	ELOg = hbar*get_wq(QGp*QMAX,LO)     ! // ~ 64 meV
!
	de=0.0005d0   ! energy step  ! Tom(2meV)->Pop(0.5meV)
	iemax=MIN(3002,nint(EMAX/de)+2)	!  Tom(0-2eV;1000div)->Pop(0-1.5eV;3002div)
!
	do ie=1,iemax
!				ei=de*float(ie) ! energy table: de ~ de*iemax
		if (ie == 1) then
			ei = EMIN
		else
			ei=de*float(ie-1) ! energy table: EMIN=1E-10(~0) ~ de*(iemax-1)
		end if
		
!-----( 非有極性光学フォノン散乱 )--- scattering type 1->6; 2->7; 3->8; 4->9; 5->5
!=====[ 光学フォノン散乱 Pop }===== scattering type 6 - 9
! 			 // the value of the overlap integral for inter-valley scattering may
!  			 // be included in the coupling constants (Jacoboni 1983, p. 672)
!
		scatOconst = pi/(2.d0*rho/hbar)
! 			  /************* optical (63 meV LO-g) phonon emission ****************/
		Nq = 1.d0/(exp(ELOg/(kb*Temp))-1.d0)
		swk(6,ie) = Zg*DLOg*DLOg*(Nq+1.d0)*DOSe(ei-ELOg)*scatOconst/ELOg
!					  /************ optical (63 meV LO-g) phonon absorption ***************/
		swk(7,ie)  = Zg*DLOg*DLOg*Nq*DOSe(ei+ELOg)*scatOconst/ELOg
		Nq = 1.d0/(exp(ETOf/(kb*Temp))-1.d0)
!					  /************* optical (59 meV TO-f) phonon emission ****************/
!	DESIGEX					swk(8,ie) = Zf*DTOf*DTOf*(Nq+1.d0)*DOSe(ei-ETOf+DESIGEX)*scatOconst/ETOf
		swk(8,ie) = Zf*DTOf*DTOf*(Nq+1.d0)*DOSe(ei-ETOf)*scatOconst/ETOf
!					  /************ optical (59 meV TO-f) phonon absorption ***************/
!	DESIGEX			 		swk(9,ie) = Zf*DTOf*DTOf*Nq*DOSe(ei+ETOf+DESIGEX)*scatOconst/ETOf
		swk(9,ie) = Zf*DTOf*DTOf*Nq*DOSe(ei+ETOf)*scatOconst/ETOf
!    
!=====[ 音響フォノンValley間散乱 Pop }===== scattering type 10 - 17
!					  /************* acoustic (12 meV TA-g) phonon emission ****************/
		Nq = 1.d0/(exp(ETAg/(kb*Temp))-1.d0) 
		swk(10,ie)  = Zg*DTAg*DTAg*(Nq+1.d0)*DOSe(ei-ETAg)*scatOconst/ETAg
!					  /************ acoustic (12 meV TA-g) phonon absorption ***************/
		swk(11,ie)  = Zg*DTAg*DTAg*Nq*DOSe(ei+ETAg)*scatOconst/ETAg
!					  /************* acoustic (47 meV LA-f) phonon emission ****************/
		Nq = 1.d0/(exp(ELAf/(kb*Temp))-1.d0)
!	DESIGEX					swk(12,ie)  = Zf*DLAf*DLAf*(Nq+1.d0)*DOSe(ei-ELAf+DESIGEX)*scatOconst/ELAf
		swk(12,ie)  = Zf*DLAf*DLAf*(Nq+1.d0)*DOSe(ei-ELAf)*scatOconst/ELAf
!					  /************ acoustic (47 meV LA-f) phonon absorption ***************/
!	DESIGEX 				swk(13,ie)  = Zf*DLAf*DLAf*Nq*DOSe(ei+ELAf+DESIGEX)*scatOconst/ELAf
		swk(13,ie)  = Zf*DLAf*DLAf*Nq*DOSe(ei+ELAf)*scatOconst/ELAf
!					  /************* acoustic (19 meV TA-f) phonon emission ****************/
		Nq = 1.d0/(exp(ETAf/(kb*Temp))-1.d0)
!	DESIGEX			 		swk(14,ie)  = Zf*DTAf*DTAf*(Nq+1.)*DOSe(ei-ETAf+DESIGEX)*scatOconst/ETAf
		swk(14,ie)  = Zf*DTAf*DTAf*(Nq+1.d0)*DOSe(ei-ETAf)*scatOconst/ETAf
!					  /************ acoustic (19 meV TA-f) phonon absorption ***************/
!	DESIGEX    				swk(15,ie)  = Zf*DTAf*DTAf*Nq*DOSe(ei+ETAf+DESIGEX)*scatOconst/ETAf
		swk(15,ie)  = Zf*DTAf*DTAf*Nq*DOSe(ei+ETAf)*scatOconst/ETAf
!					  /************* acoustic (18 meV LA-g) phonon emission ****************/
		Nq = 1.d0/(exp(ELAg/(kb*Temp))-1.d0)
		swk(16,ie)  = Zg*DLAg*DLAg*(Nq+1.)*DOSe(ei-ELAg)*scatOconst/ELAg
!					  /************ acoustic (18 meV LA-g) phonon absorption ***************/
		swk(17,ie)  = Zg*DLAg*DLAg*Nq*DOSe(ei+ELAg)*scatOconst/ELAg
!
!=====[ 音響フォノン散乱 }=====  Pop scattering type #2 - 5
!
		swk(2,ie)= acrate3(ei, LA, EMS)         ! LA intra emission
		swk(3,ie)= acrate3(ei, LA, Absorption)  ! LA intra absorption
		swk(4,ie)= acrate3(ei, TA, EMS)         ! TA intra emission
		swk(5,ie)= acrate3(ei, TA, Absorption)  ! TA intra absorption
!
!=====[ イオン化不純物散乱  }===== Pop scattering type #18
		swk(18,ie)= rateDOP(ei)
!
	end do
!        
!!----- debug scattering rate begin ------
!	do  ie=1,iemax
!		if (ie == 1) then
!			eee = EMIN
!		else
!			eee = de*float(ie-1)
!		end if
!		write (*,*) eee,swk(2,ie),swk(3,ie),swk(4,ie),swk(5,ie),swk(6,ie),swk(7,ie),swk(8,ie),swk(9,ie), &
!				& swk(10,ie),swk(11,ie),swk(12,ie),swk(13,ie),swk(14,ie),swk(15,ie),swk(16,ie),swk(17,ie),swk(18,ie)
!		write (8,*) eee,swk(2,ie),swk(3,ie),swk(4,ie),swk(5,ie),swk(6,ie),swk(7,ie),swk(8,ie),swk(9,ie), &
!				& swk(10,ie),swk(11,ie),swk(12,ie),swk(13,ie),swk(14,ie),swk(15,ie),swk(16,ie),swk(17,ie),swk(18,ie)
!	end do
!!------ debug scattering rate end -------
!
!---( 散乱レートの和の計算 )---
!     
	gm=0.
	do  ie=1,iemax
		do i=2, 18
			swk(1,ie)=swk(1,ie)+swk(i,ie)
		end do  
		if (swk(1,ie) > gm) gm = swk(1,ie)
	end do
	do  ie=1,iemax
		swk(1,ie)=gm-swk(1,ie)
		do i=2, 18
			swk(i,ie)=swk(i-1,ie)+swk(i,ie)
		end do
		do i=1, 18
			swk(i,ie)=swk(i,ie)/gm
		end do
	end do     
!----- debug scattering rate begin ------
	write (*,*) 'total scattering rate=',gm
	write (8,*) 'total scattering rate=',gm
!!      do  ie=1,iemax
!!        if (ie == 1) then
!!				eee = EMIN
!!			else
!!				eee = de*float(ie-1)
!!			end if
!!         write (*,*) eee,swk(1,ie),swk(2,ie),swk(3,ie),swk(4,ie),swk(5,ie),swk(6,ie),swk(7,ie),swk(8,ie),swk(9,ie), &
!!                    &     swk(10,ie),swk(11,ie),swk(12,ie),swk(13,ie),swk(14,ie),swk(15,ie),swk(16,ie),swk(17,ie)
!!         write (8,*) eee,swk(1,ie),swk(2,ie),swk(3,ie),swk(4,ie),swk(5,ie),swk(6,ie),swk(7,ie),swk(8,ie),swk(9,ie), &
!!                    &     swk(10,ie),swk(11,ie),swk(12,ie),swk(13,ie),swk(14,ie),swk(15,ie),swk(16,ie),swk(17,ie)
!!      end do
!!------ debug scattering rate end -------
!
end subroutine param
!========================================================================================
!===( 全粒子の初期条件 )===
!
subroutine initia(t,Elec,iv)
	implicit none
!      
	Double Precision, intent(out)::t
	Double Precision, intent(out)::Elec(:,:)
	integer, intent(out)::iv(:)
	Double Precision rn, ei, ak, cb, sb, fai, sf, r3, v_init
	integer n
	Double Precision kx, ky, kz, gk
	integer MBdist  ! 1:initial=Maxwell-Boltzman distribution
!
	t=0.d0
	call sgrnd(-99)                    !乱数 seed=1
	MBdist=1
!!	MBdist=0
	do  n=1,inum
!--------------------------------------------- begin isotropic Maxwell-Boltzman(Pop)
		if (MBdist == 1) then
			kx = smh*sqrt(kbTq/2.d0)*gasdev()    ! Local kx
			ky = smh*sqrt(kbTq/2.d0)*gasdev()    ! Local ky
			kz = smh*sqrt(kbTq/2.d0)*gasdev()    ! Local kz
!--------------------------------------------- begin isotropic initial p (Tom)
		else
			rn = grnd()
!!					write(*,*)ran1_idum,rn		!!	debud
!!					write(8,*)ran1_idum,rn		!!	debud
			ei  = -kbTq*log(rn)*1.5d0    !--- give E along exp{-E/(3/2kT)} ~39meV @300K
			ak  = smh*sqrt(ei*(1.d0+alpha*ei)) !--- E -> k
			rn = grnd()
			cb  = 1.d0-2.d0*rn ! cosθ
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
		rn = grnd()
		r3=3.*rn
		if (r3 <= 1.0) then
			if (r3 <=0.5) then
				iv(n)=1 !  第１バンド（valley at -X-axis)
			else
				iv(n)=2 !  第１バンド（valley at +X-axis)
			end if
!
			Elec(n,PX)=tm(1)*hbar*kx + pvx(iv(n)) ! 縦 Global px
			Elec(n,PY)=tm(2)*hbar*ky + pvy(iv(n)) ! 横 Global py
			Elec(n,PZ)=tm(3)*hbar*kz + pvz(iv(n)) ! 横 Global pz
			gk=hhml*kx**2+hhmt*(ky**2+kz**2)
		elseif (r3 <= 2.0) then
			if (r3 <=1.5) then
				iv(n)=3 !  第２バンド（valley at -Y-axis)
			else
				iv(n)=4 !  第２バンド（valley at +Y-axis)
			end if
!
			Elec(n,PX)=tm(3)*hbar*kx + pvx(iv(n)) ! 横 Global px
			Elec(n,PY)=tm(1)*hbar*ky + pvy(iv(n)) ! 縦 Global py
			Elec(n,PZ)=tm(2)*hbar*kz + pvz(iv(n)) ! 横 Global pz
			gk=hhml*ky**2+hhmt*(kz**2+kx**2)           
		else
			if (r3 <=2.5) then
				iv(n)=5 !  第３バンド（valley at -Z-axis)
			else
				iv(n)=6 !  第３バンド（valley at +Z-axis)
			end if
!
			Elec(n,PX)=tm(2)*hbar*kx + pvx(iv(n)) ! 横 Global px
			Elec(n,PY)=tm(3)*hbar*ky + pvy(iv(n)) ! 横 Global py
			Elec(n,PZ)=tm(1)*hbar*kz + pvz(iv(n)) ! 縦 Global pz
			gk=hhml*kz**2+hhmt*(kx**2+ky**2)
		end if
		Elec(n,EE)=(sqrt(1.d0+4.d0*alpha*gk)-1.d0)/(2.d0*alpha)  ! Tomizawa (1.7)
!
!         if (n<=100) then   ! debug
!           write(*,*) n,Elec(n,PX),Elec(n,PY),Elec(n,PZ)  ! debug
!           write(8,*) n,Elec(n,PX),Elec(n,PY),Elec(n,PZ)  ! debug
!         end if
!
	end do
end subroutine initia
!========================================================================================
!===( 多粒子モンテカルロ計算：繰り返し )===
!
subroutine monte_carlo(Elec,iv,NELph0,NETph0,NBph,DEph)
	implicit none
	Double Precision, intent(inout)::Elec(:,:)
	integer, intent(inout)::iv(:)
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph
	Double Precision, intent(in)::DEph
	Double Precision t
	integer jt
!            
	do jt=1,jtl
		t=dt*float(jt)
		call emc(t,Elec,iv,NELph0,NETph0,NBph,DEph)
		write(*,'(i10)',advance = 'no') jt   ! iteration no
		write(8,'(i10)',advance = 'no') jt   ! iteration no
		call out(t,Elec,iv)              !  debug
	end do
!
end subroutine monte_carlo
!========================================================================================
!===(時間dtの間の多粒子モンテカルロ計算 )===
!
subroutine emc(t,Elec,iv,NELph0,NETph0,NBph,DEph)
	implicit none
	Double Precision, intent(in)::t
	Double Precision, intent(inout)::Elec(:,:)
	integer, intent(inout)::iv(:)
	Double Precision kx,ky,kz,x,y,tss,eee
	Double Precision tdt, t1, tau, rn
	Double Precision Elost
	integer n, ivv
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph
	Double Precision, intent(in)::DEph
	Double Precision eeeOLD			!	Energy conservation
!
	Egetlost_drift=0.d0			!	Energy conservation 
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
		ivv = iv(n)
!
		t1=t
		do while (tss <= tdt)          ! mean free time > elapsed time -> drift
			tau=tss-t1                   ! time interval before scattering
			eeeOLD=eee										!	Energy conservation
			call drift(tau,kx,ky,kz,eee,x,y,ivv)
			Egetlost_drift=Egetlost_drift+(eee-eeeOLD)		!	Energy conservation
			eeeOLD=eee			!	Energy conservation
			call scat(kx,ky,kz,eee,ivv,NELph0,NETph0,NBph,DEph,Elost)
			Egetlost_scat=Egetlost_scat+(eee-eeeOLD)		!	Energy conservation
			t1=tss
			rn = grnd()
			tss=t1-log(rn)/gm         ! re-new mean free time
		end do ! ---- while end
       !
		tau=tdt-t1
		eeeOLD=eee											!	Energy conservation
		call drift(tau,kx,ky,kz,eee,x,y,ivv)
		Egetlost_drift=Egetlost_drift+(eee-eeeOLD)			!	Energy conservation
!
		Elec(n,PX) = hbar*kx    ! Global px
		Elec(n,PY) = hbar*ky    ! Global py
		Elec(n,PZ) = hbar*kz    ! Global pz
		Elec(n,EE) = eee
		Elec(n,TS) = tss
		Elec(n,XXX) = x
		Elec(n,YYY) = y
		Elec(n,EElost) = Elost
		iv(n)  = ivv
!
	end do
end subroutine emc
!========================================================================================
!===( ドリフト過程の計算 )===
!
subroutine drift(tau,kx,ky,kz,eee,x,y,ivv)  ! Global kx, ky, kz
	implicit none
	Double Precision, intent(inout)::  tau, kx,ky,kz,eee,x,y
	integer, intent(in):: ivv
	Double Precision dkx, skx, sky, skz, cp, sq, gk
!
	dkx=qhbar*fx*tau ! dk=q/h・F・τ
!
	cp  = hm((ivv-1)/2+1)*tau               ! m/h・τ

	gk = eee*(1.d0+alpha*eee)
	sq = sqrt(1.d0+4.d0*alpha*gk)      !--- anisotropic & non-parabolic
	x = x+cp*((kx-pvx(ivv)/hbar)+0.5d0*dkx)/sq   !--- dx=h/m・τ(kx +0.5・dkx)/sq
!
	kx=kx+dkx                   ! Global kx
	skx=(kx-pvx(ivv)/hbar)**2   ! Local skx
	sky=(ky-pvy(ivv)/hbar)**2   ! Local sky
	skz=(kz-pvz(ivv)/hbar)**2   ! Local skz
	if (ivv==1 .OR. ivv==2) then              ! --- after Tomizawa p.3-8
		gk=hhml*skx+hhmt*(sky+skz)
	elseif (ivv==3 .OR. ivv==4) then
		gk=hhml*sky+hhmt*(skz+skx)
	else    ! ivv=5 or ivv=6
		gk=hhml*skz+hhmt*(skx+sky)
	end if
	eee=(sqrt(1.d0+4.d0*alpha*gk)-1.d0)/(2.d0*alpha)            ! Tomizawa (1.7)
!
end subroutine drift
!========================================================================================
!===( 散乱過程の計算 )===
!
subroutine scat(kx,ky,kz,eee,ivv,NELph0,NETph0,NBph,DEph,Elost)  ! Global kx, ky, kz
	implicit none
	Double Precision, intent(inout)::  kx,ky,kz,eee
	integer, intent(inout):: ivv
	Double Precision ki, kf, gk, ei
	Double Precision r1, r2
	integer ie
	Double Precision qM, kspmax, Epmax
	Double Precision aed
	integer LT
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph
	Double Precision, intent(in)::DEph
	Double Precision padjust, Elost
!
!---( 散乱前のエネルギーの計算 )---
!
		Elost = 0.d0
!
!		if (eee < EMIN) then	!
!			padjust = sqrt((EMIN*(1.+alpha*EMIN))/(eee*(1.+alpha*eee)))
!			kx = pvx(ivv)/hbar + (kx - pvx(ivv)/hbar)*padjust
!			ky = pvy(ivv)/hbar + (ky - pvy(ivv)/hbar)*padjust
!			kz = pvz(ivv)/hbar + (kz - pvz(ivv)/hbar)*padjust
!			Elost = Elost + (eee - EMIN)     !				   // add up "gotton" energy
!			eee=EMIN   
!		end if
!
		if (eee > EMAX) then	!
			padjust = sqrt((EMAX*(1.+alpha*EMAX))/(eee*(1.+alpha*eee)))
			kx = pvx(ivv)/hbar + (kx - pvx(ivv)/hbar)*padjust
			ky = pvy(ivv)/hbar + (ky - pvy(ivv)/hbar)*padjust
			kz = pvz(ivv)/hbar + (kz - pvz(ivv)/hbar)*padjust
			Elost = Elost + (eee - EMAX)     !				   // add up "lost" energy
			eee=EMAX
		end if
!
	gk = eee*(1.d0+alpha*eee)
	ki = sqrt(2.d0*mdos*gk)/hbar   ! Local k
!			ei=(sqrt(1.d0+4.d0*alpha*gk)-1.d0)/(2.d0*alpha)  ! Tomizawa (1.7)
	ei = eee
	ie=int(ei/de)+1
	if (ie > iemax) ie=iemax   ! denomination to EMAX ?
!
!---( 散乱機構の選択 )---
!
	r1 = grnd()
	if (r1 <= swk(1,ie)) then       ! self-scattering
!			no operation
	elseif (r1 <= swk(5,ie)) then   !  LA/TA intra emission/absorption
		if (r1 <= swk(2,ie)) then   !  LA intra emission
			LT=LA
			aed=+1.d0
			qM=2.d0*ki
		elseif (r1 <= swk(3,ie)) then   !  LA intra absorption
			LT=LA
			aed=-1.d0
			Epmax = eee + hbar*get_wq(QMAX,LT) 
			kspmax = sqrt(2.*mdos*Epmax*(1.+alpha*Epmax))/hbar 
			qM = ki + kspmax
		elseif (r1 <= swk(4,ie)) then   !  TA intra emission
			LT=TA
			aed=+1.d0  
			qM=2.d0*ki            
		else   !    r1 <= swk(5,ie)     TA intra absorption
			LT=TA
			aed=-1.d0
			Epmax = eee + hbar*get_wq(QMAX,LT)
			kspmax = sqrt(2.d0*mdos*Epmax*(1.d0+alpha*Epmax))/hbar
			qM = ki + kspmax
		end if
!
		call final_state_intra_scatByLATA(kx,ky,kz,ki,eee,ivv,LT,aed,Epmax,kspmax,qM,NELph0,NETph0,NBph,DEph)
!
	elseif (r1 <= swk(6,ie)) then       !  O-g-emis
		LT=LO
		aed=+1.d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
	elseif (r1 <= swk(7,ie)) then   !  O-g-abs
		LT=LO
		aed=-1.d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
	elseif (r1 <= swk(8,ie)) then   !  O-f-emis
		LT=TO
		aed=+1.d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph)
	elseif (r1 <= swk(9,ie)) then   !  O-f-abs
		LT=TO
		aed=-1.d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph)
!
	else if (r1 <= swk(10,ie)) then   ! TA-g emission
		LT=TA
		aed=+1.d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz 
!
	else if (r1 <= swk(11,ie)) then   ! TA-g absorption
		LT=TA
		aed=-1.d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
!
	else if (r1 <= swk(12,ie)) then   ! LA-f emission
		LT=LA
		aed=+1.d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz 	
!
	else if (r1 <= swk(13,ie)) then   ! LA-f absorption
		LT=LA
		aed=-1.d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
!
	else if (r1 <= swk(14,ie)) then   ! TA-f emission
		LT=TA
		aed=+1.d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
!
	else if (r1 <= swk(15,ie)) then   ! TA-f absorption
		LT=TA
		aed=-1.d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
!

	else if (r1 <= swk(16,ie)) then   ! LA-g emission
		LT=LA
		aed=+1.d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz 
!

	else if (r1 <= swk(17,ie)) then   ! LA-g absorption
		LT=LA
		aed=-1.d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph)	 ! Global kx, ky, kz 
!
	else if (r1 <= swk(18,ie)) then   ! ion scattering
		call final_state_impurity_scat(kx,ky,kz,eee,ivv)								 ! Global kx, ky, kz 
	end if
!      
end subroutine scat
!========================================================================================
subroutine final_state_intra_scatByLATA(kx,ky,kz,ki,eee,ivv,LT,aed,Epmax,kspmax,qM,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
! common treatment for intra_ScatByLA&TA:
	implicit none
	Double Precision, intent(inout)::  kx,ky,kz,eee
	Double Precision ki
	integer, intent(inout):: ivv
	integer LT
	Double Precision aed
	Double Precision C, qM, kspmax, Epmax, q0, f0, fx, ephon, acost
	Double Precision qdum, edum, pprime, kprime
	Double Precision cbet, sbet, cgam, sgam, cgamp, cgampp, eta, ceta, seta, phi, sphi, cphi
	integer cnt
	integer BZCLIP
	Double Precision calp, salp, bet
!!	Double Precision,parameter::PQMAX = 7.6149280673e-08 ! hbar*QMAX(Wigner-Seitz cell radius)
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph
	Double Precision, intent(in)::DEph
	integer ie
	Double Precision Pxyz
!      
	C = qM*get_Mq(qM,LT,aed)*(1.d0+2.d0*alpha*(eee-aed*hbar*get_wq(qM,LT))) 
!
	q0 = QMIN + (qM-QMIN)*grnd()
	ephon = hbar*get_wq(q0,LT)
	f0 = C*grnd()
	fx = q0*get_Mq(q0,LT,aed)*(1.d0+2.d0*alpha*(eee-aed*ephon))
	acost = fabs_cost(q0,ki,eee,LT,aed)
	cnt=1									! debug
	do while ((f0 > fx) .OR. (acost > 1.d0))
		q0 = QMIN + (qM-QMIN)*grnd()
		ephon = hbar*get_wq(q0,LT)
		f0 = C*grnd()
		fx = q0*get_Mq(q0,LT,aed)*(1.d0+2.d0*alpha*(eee-aed*ephon))
		acost = fabs_cost(q0,ki,eee,LT,aed)
		cnt=cnt+1								! debug
		if (mod(cnt,100000)==0) then			! debug
			write(*,*) 'cnt=',cnt,'eee=',eee,'acost=',acost,'LT=',LT,'aed=',aed    ! debug
			write(8,*) 'cnt=',cnt,'eee=',eee,'acost=',acost,'LT=',LT,'aed=',aed    ! debug
		end if
	end do
!               
!!!  set the electron energy, consistent w/the momentum exchange above
	eee = eee - aed*ephon
!!!  return the phonon energy and momentum for STATS
!
!		**  phonon counting(+emit & -absorb)
	ie=int(ephon/DEph)+1
	if (ie > NBph) then
		write(*,*) '** WARNING 1 ** too large ephon',ephon
	else if (LT==LA .OR. LT==LO) then
		NELph0(ie)=NELph0(ie)+aed*1.d0
	else if (LT==TA .OR. LT==TO) then
		NETph0(ie)=NETph0(ie)+aed*1.d0
	end if
!
	qdum = q0
	edum = ephon
!        
!--- calculation of final state for LA/TA phonon emission/absorption (type=2-5)
!   /********* non-parabolicity from tomizawa eq. 1.11 page 8 *********/
	pprime = sqrt(2.d0*mdos*eee*(1.d0+alpha*eee)) !
	kprime = pprime/hbar  ! k=p/hbar
!				  /***** momentum conservation for angle beta b/w ki and kprime *****/
	cbet = (ki*ki+kprime*kprime-q0*q0)/(2.d0*ki*kprime) ! low of cosine
!
	if (abs(cbet) > 1.0) then
		write(*,*) 'error abs(cbet) > 1.0'
		write(*,*) 'ki=',ki,' kprime=',kprime,' q0=',q0,' cbet=',cbet    !debug
		write(8,*) 'ki=',ki,' kprime=',kprime,' q0=',q0,' cbet=',cbet    !debug
!		 	 exit(-1)
		stop
	end if
	sbet = sqrt(1.d0-cbet*cbet)
!		          /*** gam = angle b/w ki and Z-axis ***/
	cgam = (hbar*kz-pvz(ivv))*sqrt(mdos/mz(ivv))/(hbar*ki)
	if (abs(cgam) > 1.0) then
		write(*,*) 'error abs(cgam) > 1.0'
		write(*,*) 'cbet=',cbet,' sbet=',sbet,' cgam=',cgam    !debug
		write(8,*) 'cbet=',cbet,' sbet=',sbet,' cgam=',cgam    !debug
!	 			 exit(-1)
		stop
	end if
	sgam = sqrt(1.d0-cgam*cgam)
!			/*** gamp = angle b/w ki and Y-axis ***/
	cgamp = (hbar*ky-pvy(ivv))*sqrt(mdos/my(ivv))/(hbar*ki)
!			/*** gampp = angle b/w ki and X-axis ***/
	cgampp = (hbar*kx-pvx(ivv))*sqrt(mdos/mx(ivv))/(hbar*ki)
!			/*** angle b/w ki projection on X-Y plane and Y-axis ***/
	eta = atan((hbar*kx-pvx(ivv))*my(ivv)/   &
		  &	 ((hbar*ky-pvy(ivv))*mx(ivv)))
	ceta = cos(eta)
	seta = sin(eta)
!
!     begin   once for while
!			/*** rotation phi chosen at random ***/
	phi = 2.d0*pi*grnd()
	sphi = sin(phi) 
	cphi = cos(phi) 
!			/*** compute final state momentum components, like Jacoboni Rev. Mod.
!			 1983, eq. 3.62-3.64, but he has X and Z components mixed up ***/
	kz = ((cbet*cgam   +sbet*cphi*sgam)* &
	   & sqrt(mz(ivv)/mdos)*pprime + pvz(ivv))/hbar
	ky = ((cbet*cgamp  -sbet*cphi*cgam*ceta -sbet*sphi*seta)* &
	   & sqrt(my(ivv)/mdos)*pprime + pvy(ivv))/hbar
	kx = ((cbet*cgampp -sbet*cphi*cgam*seta +sbet*sphi*ceta)* &
	   & sqrt(mx(ivv)/mdos)*pprime + pvx(ivv))/hbar
	cnt=1   !! debug
	BZCLIP=1      ! debug
	if (ivv==1 .OR. ivv==2) then
		Pxyz=abs(kx)*hbar
	else if (ivv==3 .OR. ivv==4) then
		Pxyz=abs(ky)*hbar
	else if (ivv==5 .OR. ivv==6) then
		Pxyz=abs(kz)*hbar
	end if
!		end   once for while
!		do while ((abs(Elec[Pxyz][j]) > PQMAX) && (cnt < 200) && (BZCLIP)) ! <37>{
!			write(*,*) abs(kx),abs(ky),abs(kz), PQMAX/hbar  ! debug
!			write(8,*) abs(kx),abs(ky),abs(kz), PQMAX/hbar  ! debug
!!		do while ((MAX(abs(kx),abs(ky),abs(kz)) > PQMAX/hbar) .AND. (cnt < 200) .AND. (BZCLIP==1))  ! <37>{
	do while ((Pxyz > PQMAX) .AND. (cnt < 200) .AND. (BZCLIP==1))
!			/*** rotation phi chosen at random ***/
		phi = 2.d0*pi*grnd()
		sphi = sin(phi)
		cphi = cos(phi)
!     		/*** compute final state momentum components, like Jacoboni Rev. Mod.
!			1983, eq. 3.62-3.64, but he has X and Z components mixed up ***/
		kz = ((cbet*cgam   +sbet*cphi*sgam)* &
		   & sqrt(mz(ivv)/mdos)*pprime + pvz(ivv))/hbar
		ky = ((cbet*cgamp  -sbet*cphi*cgam*ceta -sbet*sphi*seta)* &
		   & sqrt(my(ivv)/mdos)*pprime + pvy(ivv))/hbar
		kx = ((cbet*cgampp -sbet*cphi*cgam*seta +sbet*sphi*ceta)* &
		   & sqrt(mx(ivv)/mdos)*pprime + pvx(ivv))/hbar
		if (ivv==1 .OR. ivv==2) then
			Pxyz=abs(kx)*hbar
		else if (ivv==3 .OR. ivv==4) then
			Pxyz=abs(ky)*hbar
		else if (ivv==5 .OR. ivv==6) then
			Pxyz=abs(kz)*hbar
		end if
!      		/*** increment counter ***/
		cnt = cnt + 1
	end do    !  }<37> 
!		  // if BZCLIP = 0 this && is always FALSE so it only goes thru ONCE
!		  // if the counter exceeds 200 reps then make it isotropic
	if (cnt > 200) then !  <39>{          // note:  this never happens if BZCLIP = 0
!		  if (cnt > 200) write(*,*) '** WARNING (intra_scatByLATA) **  cnt =',cnt   ! debug
		write(*,*) '** WARNING 4 (intra_scatByLATA) **  cnt =',cnt
!			     /***** compute ISOTROPIC new angles if all else fails above ******/
!       -- begin    once for while
		calp = 1.d0-2.d0*grnd()
		salp = sqrt(1.d0-calp*calp)
		bet = 2.d0*pi*grnd()
!			    /********** compute new isotropic momentum components *************/
		kx = (pvx(ivv) + sqrt(mx(ivv)/mdos)*pprime*salp*cos(bet))/hbar
		ky = (pvy(ivv) + sqrt(my(ivv)/mdos)*pprime*salp*sin(bet))/hbar
		kz = (pvz(ivv) + sqrt(mz(ivv)/mdos)*pprime*calp)/hbar
!       -- end      once for while
		if (ivv==1 .OR. ivv==2) then
			Pxyz=abs(kx)*hbar
		else if (ivv==3 .OR. ivv==4) then
			Pxyz=abs(ky)*hbar
		else if (ivv==5 .OR. ivv==6) then
			Pxyz=abs(kz)*hbar
		end if
!			do while ((fabs(Elec[Pxyz][j]) > PQMAX) && (BZCLIP))     ! <38>{
		cnt=1				!! debug
		do while ((Pxyz > PQMAX) .AND. (BZCLIP==1))   ! <38>{
			calp = 1.d0-2.d0*grnd()
			salp = sqrt(1.d0-calp*calp)
			bet = 2.d0*pi*grnd()
!			  /********** compute new isotropic momentum components *************/
			kx = (pvx(ivv) + sqrt(mx(ivv)/mdos)*pprime*salp*cos(bet))/hbar
			ky = (pvy(ivv) + sqrt(my(ivv)/mdos)*pprime*salp*sin(bet))/hbar
			kz = (pvz(ivv) + sqrt(mz(ivv)/mdos)*pprime*calp)/hbar
			if (ivv==1 .OR. ivv==2) then
				Pxyz=abs(kx)*hbar
			else if (ivv==3 .OR. ivv==4) then
				Pxyz=abs(ky)*hbar
			else if (ivv==5 .OR. ivv==6) then
				Pxyz=abs(kz)*hbar
			end if
			cnt=cnt+1				!! debug
			if (cnt==10000) then				!! debug
				write(*,*) 'isotropic while cnt=', cnt				!! debug
				write(8,*) 'isotropic while cnt=', cnt				!! debug
			end if				!! debug
		end do   ! }<38> 
	end if     ! }<39>
!	}<40>
!--- end calculation of final state for LA/TA phonon emission/absorption (type=2-5)
end subroutine final_state_intra_ScatByLATA
!========================================================================================
subroutine final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
	implicit none
	Double Precision, intent(inout)::  kx,ky,kz,eee
	Double Precision aed
	integer ivv,LT
	integer cnt, REDO, BZCLIP
	Double Precision ephon, ephon0, qphon, qphon0, Ep, bet, calp, salp
	Double Precision pp, pxo, pyo, pzo, dpx, dpy, dpz, r
	Double Precision qdum, edum
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph
	Double Precision, intent(in)::DEph
	integer ie
	integer Pxyz
	integer cnt2				!! debug
	Double Precision, parameter:: ephon_convergence=0.05d0
!
	cnt=0
	REDO=1
!
!			#ifdef POP   // store original momentum components
	if ((aed > 0.) .AND. ((eee - hbar*get_wq(QGp*QMAX,LT)) <= EMIN)) then
!		write (*,*) 'eee - hbar*get_wq(QGp*QMAX,LT) <= EMIN'
		return   !  do nothing OK? 2018.10.10
	end if
	pxo = hbar*kx
	pyo = hbar*ky
	pzo = hbar*kz
!
	call change_equiv_valley(ivv)
!
	dpx = pvx(ivv)-pxo
	dpy = pvy(ivv)-pyo
	dpz = pvz(ivv)-pzo
	qphon0 = abs(2.*PQMAX-sqrt(dpx*dpx + dpy*dpy + dpz*dpz))/hbar
	if (qphon0 < QMIN) write(*,*) '** WARNING 5 inter_scat_g **: qphon0 < QMIN', qphon0/QMAX
!			// original estimate on phonon energy
	ephon = hbar*get_wq(qphon0,LT)
!			?  do ⑧{
	do while ( REDO > 0 )
!			#endif
		ephon0 = ephon
!				// estimate new electron energy
		Ep = eee - aed*ephon0
!			#ifdef POP
		if (Ep > EMIN)  then   !  ⑤{
!			#endif
			BZCLIP=1   !!
!				// compute new momentum in Herring-Vogt space
			pp = sqrt(2.d0*mdos*Ep*(1.d0+alpha*Ep))
!				// choose the set of random angles
			calp = 1.d0-2.d0*grnd()
			salp = sqrt(1.d0-calp*calp)
			bet = 2.d0*pi*grnd()
!				// new momentum components in new valley iv[j]
			kx = (pvx(ivv) + sqrt(mx(ivv)/mdos)*pp*salp*cos(bet))/hbar
			ky = (pvy(ivv) + sqrt(my(ivv)/mdos)*pp*salp*sin(bet))/hbar
			kz = (pvz(ivv) + sqrt(mz(ivv)/mdos)*pp*calp)/hbar
			if (ivv==1 .OR. ivv==2) then
				Pxyz=abs(kx)*hbar
			else if (ivv==3 .OR. ivv==4) then
				Pxyz=abs(ky)*hbar
			else if (ivv==5 .OR. ivv==6) then
				Pxyz=abs(kz)*hbar
			end if
			cnt2=1				!! debug
			do while ((Pxyz > PQMAX) .AND. (BZCLIP==1)) 
!				// choose the set of random angles
				calp = 1.d0-2.d0*grnd()
				salp = sqrt(1.d0-calp*calp)
				bet = 2.d0*pi*grnd()
!				// new momentum components in new valley iv[j]
				kx = (pvx(ivv) + sqrt(mx(ivv)/mdos)*pp*salp*cos(bet))/hbar
				ky = (pvy(ivv) + sqrt(my(ivv)/mdos)*pp*salp*sin(bet))/hbar
				kz = (pvz(ivv) + sqrt(mz(ivv)/mdos)*pp*calp)/hbar
				if (ivv==1 .OR. ivv==2) then
					Pxyz=abs(kx)*hbar
				else if (ivv==3 .OR. ivv==4) then
					Pxyz=abs(ky)*hbar
				else if (ivv==5 .OR. ivv==6) then
					Pxyz=abs(kz)*hbar
				end if
				cnt2=cnt2+1				!! debug
				if (cnt2==10000) then				!! debug
					write(*,*) 'before WARNING 6 while cnt=', cnt2				!! debug
					write(8,*) 'before WARNING 6 while cnt=', cnt2				!! debug
				end if				!! debug
			end do
!	  	  
!				#ifdef POP
!				// difference between final & initial momentum vectors
			dpx = hbar*kx-pxo
			dpy = hbar*ky-pyo
			dpz = hbar*kz-pzo
!				// get qphon = norm of the resultant q-vector reduced to 1st BZ
			qphon = abs(2.d0*PQMAX-sqrt(dpx*dpx + dpy*dpy + dpz*dpz))/hbar
			if (qphon < QMIN) write(*,*) '** WARNING 6 inter_scat_g **: qphon < QMIN', qphon/QMAX
!				// get resulting phonon energy and momentum
			ephon = hbar*get_wq(qphon,LT) 
!				// estimate new electron energy
			Ep = eee - aed*ephon
!				// upper limit for cnt = 200, see right below
			cnt=cnt+1
			if ((Ep > EMIN) .AND. ((abs(ephon-ephon0)/ephon0) < ephon_convergence)) REDO = 0
		else
!				// original estimate on phonon energy
			qphon = QGp*QMAX 
			ephon = hbar*get_wq(qphon,LT) 
		end if
!				// allow it to break out if cnt goes too high
		if ((cnt > 200) .AND. (REDO > 0)) then
			REDO = 0
			qphon = QGp*QMAX
			ephon = hbar*get_wq(qphon,LT)
			if (cnt > 200) write(*,*) '** WARNING 7  inter_scat_g **  cnt =',cnt   ! debug
		end if
!				// empirical note: this happens VERY rarely so it's OK
	end do
!				}⑧ while (REDO)
!				#endif
!  			// set the electron energy, consistent w/the momentum exchange above
	eee = eee - aed*ephon 
!
!			if (aed > 0.) assert(Elec[EE][j] > EMIN)
!			// return the phonon energy and momentum for STATS
!
!		**  phonon counting(+emit & -absorb)
	ie=int(ephon/DEph)+1
	if (ie > NBph) then
		write(*,*) '** WARNING 8 ** too large ephon',ephon
	else if (LT==LA .OR. LT==LO) then
		NELph0(ie)=NELph0(ie)+aed*1.d0
	else if (LT==TA .OR. LT==TO) then
		NETph0(ie)=NETph0(ie)+aed*1.d0
	end if
!
	qdum = qphon
	edum = ephon
!			// consistent within 5% with the electron energy exchange [based on
!			// ephon0] but here i'm keeping qphon + ephon for consistency within
!			// the generated phonon statistics
end subroutine final_state_inter_scat_g
!========================================================================================
subroutine final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
	implicit none
	Double Precision, intent(inout)::  kx,ky,kz,eee
!!	Double Precision, parameter:: EMIN = 1.e-10   ! @ Pop-definition
!!	Double Precision,parameter::PQMAX = 7.6149280673e-08 ! hbar*QMAX(Wigner-Seitz cell radius)
	Double Precision aed
	integer ivv,LT
	integer cnt, REDO, BZCLIP
	Double Precision ephon, ephon0, qphon, Ep, bet, calp, salp
	Double Precision pp, pxo, pyo, pzo, dpx, dpy, dpz, r, EX, E0
	Double Precision qdum, edum
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph
	Double Precision, intent(in)::DEph
	integer ie
	integer Pxyz
	integer cnt2				!! debug
	Double Precision, parameter:: ephon_convergence=0.05d0
!    
	EX=0.d0
	cnt=0
	REDO=1
  
!		eee = E0  = Elec[EE][j]
!		ivv = iv0 = iv[j]

!			#ifdef POP   // store original momentum components and valley
	if ((aed > 0.) .AND. ((eee - hbar*get_wq(QFp*QMAX,LT)) <= EMIN)) then
!			write (*,*) 'eee - hbar*get_wq(QFp*QMAX,LT) <= EMIN'
		return   !  do nothing OK? 2018.10.12
	end if
	pxo = hbar*kx
	pyo = hbar*ky
	pzo = hbar*kz
!
	call change_noneq_valley(ivv)
!			}
!			#ifdef POP
	dpx = pvx(ivv)-pxo
	dpy = pvy(ivv)-pyo
	dpz = pvz(ivv)-pzo
	dpx = PQMAX-abs(dpx)
	dpy = PQMAX-abs(dpy)
	dpz = PQMAX-abs(dpz)
	qphon = sqrt(dpx*dpx + dpy*dpy + dpz*dpz)/hbar
!			#endif
!  
!			/* to do: for f-scat with LA/LO phonons, choose if it's LA or LO based
!			on a random number weighed by the occupation of those branches i.e.
!			Nq(LA) vs Nq(LO) for whatever qphon seems to satisfy the particular
!			transition.  maybe shoot an email to Ravaioli?
!			// new longitudinal axis momentum component
!			Pxyz = PX-1+(int)(iv[j]+1)/2
!			#ifdef POP
	if (qphon < QMIN) write(*,*) '** WARNING 9 inter_scat_f qphon < QMIN **', qphon/QMAX
!			// original estimate on phonon energy
	ephon = hbar*get_wq(qphon,LT)
!			// remember inside get_wq() if qphon > QMAX then phonon is LO not LA
!			do {
	do while ( REDO > 0 )
!			#endif
		ephon0 = ephon
!				// estimate new electron energy
		Ep = eee - aed*ephon0 + EX
!			#ifdef POP
		if (Ep > EMIN) then    !	{
!			#endif
      		BZCLIP=1   !!
!	  			// compute new momentum in Herring-Vogt space
			pp = sqrt(2.d0*mdos*Ep*(1.d0+alpha*Ep))
!				// choose the set of random angles
			calp = 1.d0-2.d0*grnd()
			salp = sqrt(1.d0-calp*calp)
			bet = 2.d0*pi*grnd()
!				// new momentum components in new valley iv[j]
			kx = (pvx(ivv) + sqrt(mx(ivv)/mdos)*pp*salp*cos(bet))/hbar
			ky = (pvy(ivv) + sqrt(my(ivv)/mdos)*pp*salp*sin(bet))/hbar
			kz = (pvz(ivv) + sqrt(mz(ivv)/mdos)*pp*calp)/hbar
			if (ivv==1 .OR. ivv==2) then
				Pxyz=abs(kx)*hbar
			else if (ivv==3 .OR. ivv==4) then
				Pxyz=abs(ky)*hbar
			else if (ivv==5 .OR. ivv==6) then
				Pxyz=abs(kz)*hbar
			end if
			cnt2=1				!! debug
			do while ((Pxyz > PQMAX) .AND. (BZCLIP==1))
!				// choose the set of random angles
				calp = 1.d0-2.d0*grnd()
				salp = sqrt(1.d0-calp*calp)
				bet = 2.d0*pi*grnd()
!				// new momentum components in new valley iv[j]
				kx = (pvx(ivv) + sqrt(mx(ivv)/mdos)*pp*salp*cos(bet))/hbar
				ky = (pvy(ivv) + sqrt(my(ivv)/mdos)*pp*salp*sin(bet))/hbar
				kz = (pvz(ivv) + sqrt(mz(ivv)/mdos)*pp*calp)/hbar
				if (ivv==1 .OR. ivv==2) then
					Pxyz=abs(kx)*hbar
				else if (ivv==3 .OR. ivv==4) then
					Pxyz=abs(ky)*hbar
				else if (ivv==5 .OR. ivv==6) then
					Pxyz=abs(kz)*hbar
				end if
				cnt2=cnt2+1				!! debug
				if (cnt2==10000) then				!! debug
					write(*,*) 'before WARNING 10 while cnt=', cnt2				!! debug
					write(8,*) 'before WARNING 10 while cnt=', cnt2				!! debug
				end if				!! debug
			end do
!
!				#ifdef POP
!				// difference between final & initial momentum vectors
			dpx = hbar*kx-pxo
			dpy = hbar*ky-pyo
			dpz = hbar*kz-pzo
!				// get qphon = norm of the resultant q-vector reduced to 1st BZ
			dpx = PQMAX-abs(dpx)
			dpy = PQMAX-abs(dpy)
			dpz = PQMAX-abs(dpz)
			qphon = sqrt(dpx*dpx + dpy*dpy + dpz*dpz)/hbar	!  // always positive
			if (qphon < QMIN) write(*,*) '** WARNING 10 inter_scat_f qphon < QMIN:', qphon/QMAX
!				// get resulting phonon energy and momentum
			ephon = hbar*get_wq(qphon,LT)
!				// re-estimate new electron energy
			Ep = eee - aed*ephon0 + EX
!				// upper limit for cnt = 200, see right below
			cnt=cnt+1
			if ((Ep > EMIN) .AND. ((abs(ephon-ephon0)/ephon0) < ephon_convergence)) REDO = 0
		else
!				 // go back to basic estimate on phonon energy
			qphon = QFp*QMAX
			ephon = hbar*get_wq(qphon,LT)
		end if
!				// allow it to break out if cnt goes too high
		if ((cnt > 200) .AND. (REDO > 0)) then
			REDO = 0
			qphon = QFp*QMAX
			ephon = hbar*get_wq(qphon,LT)
			if (cnt > 200) write(*,*) '** WARNING 11 inter_scat_f **  cnt =',cnt   ! debug
		end if
!				// empirical note: this happens VERY rarely so it's OK
    end do
!				// reject several unwanted scenarios and tighten ephon guess
!				} while (REDO)
!				#endif
!				// set the electron energy, consistent w/the momentum exchange above
	eee = eee - aed*ephon0 + EX 
!				if (aed > 0.) assert(Elec[EE][j] > EMIN)
!				// return the phonon energy and momentum for STATS
!
!		**  phonon counting(+emit & -absorb)
	ie=int(ephon/DEph)+1
	if (ie > NBph) then
		write(*,*) '** WARNING 12 ** too large ephon',ephon
	else if (LT==LA .OR. LT==LO) then
		NELph0(ie)=NELph0(ie)+aed*1.d0
	else if (LT==TA .OR. LT==TO) then
		NETph0(ie)=NETph0(ie)+aed*1.d0
	end if
!
	qdum = qphon
	edum = ephon
!				// this should be consistent within 5% with the electron energy
!				// exchange [based on ephon0] but here i'm keeping qphon + ephon
!				// for consistency within the generated phonon statistics
!				//if (cnt > 200) write(*,*) 'cnt =',cnt ! kjh: never happen
end subroutine final_state_inter_scat_f
!========================================================================================
subroutine change_equiv_valley(ivv)
	implicit none
	integer ivv
!
	if      (ivv == 1) then
		ivv=2
	else if (ivv == 2) then
		ivv=1
	else if (ivv == 3) then
		ivv=4
	else if (ivv == 4) then
		ivv=3
	else if (ivv == 5) then
		ivv=6
	else if (ivv == 6) then
		ivv=5
	end if  
!
!!  axis = (iv[j] + 1)/2
!!  iv[j] = 4*axis - iv[j] - 1
end subroutine change_equiv_valley
!========================================================================================
subroutine change_noneq_valley(ivv)
	implicit none
	integer ivv
	Double Precision r
!
	r = grnd()
	if      ((ivv == 1) .OR. (ivv == 2))  then
		ivv=3+int(4.d0*r)  ! ivv = 1 or 2 -> {3,4,5,6}
	else if ((ivv == 3) .OR. (ivv == 4))  then    ! ivv = 3 or 4 -> {1,2}{5,6}
		ivv = 1 + int(4.d0*r)
		if (ivv > 2) ivv= ivv+2
	else 
		ivv = 1 + int(4.d0*r)
	end if                                        ! ivv = 5 or 6 -> {1,2,3,4}
end subroutine change_noneq_valley
!========================================================================================
subroutine final_state_impurity_scat(kx,ky,kz,eee,ivv)
	implicit none
	Double Precision, intent(inout)::  kx,ky,kz,eee
	integer ivv
	Double Precision pprime,calp,salp,bet,Pxyz
	integer BZCLIP, cnt
!
	BZCLIP=1   !!
	pprime = sqrt(2.d0*mdos*eee*(1.d0+alpha*eee)) !
!
	calp = 1.d0-2.d0*grnd()
	salp = sqrt(1.d0-calp*calp)
	bet = 2.d0*pi*grnd()
	kx = (pvx(ivv) + sqrt(mx(ivv)/mdos)*pprime*salp*cos(bet))/hbar
	ky = (pvy(ivv) + sqrt(my(ivv)/mdos)*pprime*salp*sin(bet))/hbar
	kz = (pvz(ivv) + sqrt(mz(ivv)/mdos)*pprime*calp)/hbar
!
	if (ivv==1 .OR. ivv==2) then
		Pxyz=abs(kx)*hbar
	else if (ivv==3 .OR. ivv==4) then
		Pxyz=abs(ky)*hbar
	else if (ivv==5 .OR. ivv==6) then
		Pxyz=abs(kz)*hbar
	end if
!
	do while ((Pxyz > PQMAX) .AND. (BZCLIP==1)) 
		calp = 1.d0-2.d0*grnd()
		salp = sqrt(1.d0-calp*calp)
		bet = 2.d0*pi*grnd()
		kx = (pvx(ivv) + sqrt(mx(ivv)/mdos)*pprime*salp*cos(bet))/hbar
		ky = (pvy(ivv) + sqrt(my(ivv)/mdos)*pprime*salp*sin(bet))/hbar
		kz = (pvz(ivv) + sqrt(mz(ivv)/mdos)*pprime*calp)/hbar
!
		if (ivv==1 .OR. ivv==2) then
			Pxyz=abs(kx)*hbar
		else if (ivv==3 .OR. ivv==4) then
			Pxyz=abs(ky)*hbar
		else if (ivv==5 .OR. ivv==6) then
			Pxyz=abs(kz)*hbar
		end if
!
		cnt=cnt+1				!! debug
		if (MOD(cnt,10000)==0) then				!! debug
			write(*,*) 'before WARNING 10 while cnt=', cnt				!! debug
			write(8,*) 'before WARNING 10 while cnt=', cnt				!! debug
		end if				!! debug
	end do
!
end subroutine final_state_impurity_scat
!
!========================================================================================
!===( 出力ルーチン )===
!
subroutine out(t,Elec,iv)
	implicit none
	Double Precision, intent(in)::t
	Double Precision, intent(in)::Elec(:,:)
	integer, intent(in)::iv(:)
	Double Precision ve, s, kx, ky, kz, gk
	Double Precision kk,eee                       ! debug
	Double Precision eeeLost,eeeMin
	integer n, ivv
!
!---( 平均速度と距離とエネルギー保存の計算 )---
!
	ve = 0.d0
	eee = 0.d0   ! debug
	s = 0.d0
	eeeLost=0.d0
	eeeMin=999.d0
	do  n=1,inum
		ivv  = iv(n)
		kx = (Elec(n,PX) - pvx(ivv))/hbar   ! Local kx
!			ky = (Elec(n,PY) - pvy(ivv))/hbar   ! Local ky
!			kz = (Elec(n,PZ) - pvz(ivv))/hbar   ! Local kz
		if (ivv==1 .OR. ivv==2) then                  ! x-direction // E field
			gk = hhml*kx**2 + hhmt*(ky**2+kz**2)
			kk = (kx/tm(1))**2+(ky/tm(2))**2+(kz/tm(3))**2  ! debug
		elseif (ivv==3 .OR. ivv==4) then              ! y-direction ⊥ E field
			gk = hhml*ky**2 + hhmt*(kz**2+kx**2)
			kk = (kx/tm(3))**2+(ky/tm(1))**2+(kz/tm(2))**2  ! debug
		else                               ! z-direction ⊥ E field
			gk = hhml*kz**2 + hhmt*(kx**2+ky**2)
			kk = (kx/tm(2))**2+(ky/tm(3))**2+(kz/tm(1))**2  ! debug
		end if
!
		ve = ve + kx*hm((ivv-1)/2+1)/sqrt(1.d0+4.d0*alpha*gk) 
!			eee = eee+kk                              ! debug
		eee = eee + Elec(n,EE)
		if (Elec(n,EE) < eeeMin) eeeMin=Elec(n,EE)
		s  = s + Elec(n,XXX)
		eeeLost=eeeLost+Elec(n,EElost)
	end do
	ve = ve/float(inum)
!		eee = eee/float(inum)/smh/smh            ! debug kk->E[eV]
	eee = eee/float(inum)
	s  = s /float(inum)
!
!---( 時間 vs 平均速度 vs 距離の出力 )---
	write(*,*) t,s,ve,eee,eeeLost,eeeMin,Egetlost_scat,Egetlost_drift
	write(8,*) t,s,ve,eee,eeeLost,eeeMin,Egetlost_scat,Egetlost_drift
end subroutine out
!  
!========================================================================================
!===( エネルギー分布計算 )===
!
subroutine energy_dist(Elec)
	implicit none
	Double Precision, intent(in)::Elec(:,:)
!		integer, intent(in)::iv(:)
!		Double Precision,allocatable ::eee(:)  ! energy of each particle
!		Double Precision kx, ky, kz, gk
!		Double Precision  kk                       !
	integer n
!
	integer, allocatable ::ehist(:)
	integer:: nhist=301                !ヒストグラム分割数
	integer digitized_energy, iemax
	Double Precision eemax
!
!		allocate(eee(inum))
	allocate(ehist(nhist))
!      
!---( 個々の粒子のエネルギー計算 )---
!
	eemax=0.d0                 ! do not use ::emax=0.d0
	do  n=1,inum
!         ivv  = iv(n)
!         kx = (Elec(n,PX) - pvx(ivv))/hbar   ! Local kx
!         ky = (Elec(n,PY) - pvy(ivv))/hbar   ! Local ky
!         kz = (Elec(n,PZ) - pvz(ivv))/hbar   ! Local kz
!        if (ivv==1 .OR. ivv==2) then                  ! x-direction // E field
!            gk = hhml*kx**2 + hhmt*(ky**2+kz**2)
!            eee(n) =( (kx/tm(1))**2+(ky/tm(2))**2+(kz/tm(3))**2 )/smh/smh  !
!            eee(n) = eee(n)/smh/smh
!         else if (ivv==3 .OR. ivv==4) then              ! y-direction ⊥ E field
!            gk = hhml*ky**2 + hhmt*(kz**2+kx**2)
!            eee(n) =( (kx/tm(3))**2+(ky/tm(1))**2+(kz/tm(2))**2 )/smh/smh  !
!!            eee(n) = eee(n)/smh/smh
!         else                               ! z-direction ⊥ E field
!            gk = hhml*kz**2 + hhmt*(kx**2+ky**2)
!            eee(n) =( (kx/tm(2))**2+(ky/tm(3))**2+(kz/tm(1))**2 )/smh/smh  !
!!            eee(n) = eee(n)/smh/smh
!         end if
!         
!         eee(n)=Elec(n,EE)   ! 2018.10.13 debug
!
		if (Elec(n,EE) > eemax) then
			iemax = n              ! the same particle?
			eemax = Elec(n,EE)
		endif
	end do
!
	do  n=1,nhist
		ehist(n)=0
	end do
!
!			eemax = 0.01d0   ! Debug fine division  10meV/300div -> 0.03meV/div
	do  n=1,inum
		digitized_energy=int(real(nhist-1)*(Elec(n,EE)/eemax))+1
		if (digitized_energy <= nhist) then
			ehist(digitized_energy)=ehist(digitized_energy)+1
		end if
	end do   
!
!---( エネルギーヒストグラムの出力 )---
!
	write(*,*) 'Emax=',eemax, '  particle#', iemax ! debug
	write(8,*) 'Emax=',eemax, '  particle#', iemax ! debug
!
	do n=1,nhist
!			write(*,*) real(n)/real(nhist)*eemax,'  ',ehist(n)
		write(8,*) real(n)/real(nhist)*eemax,'  ',ehist(n)
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
	integer, allocatable::iv(:)
	integer inum_dummy, i
	integer, parameter::NBph=280	!	2018.10.29 280->1000
	Double Precision, parameter::DEph=70.d-3/float(NBph)	!	2.5e-4  ! 70meV/280div
	Double Precision, allocatable ::NELph0(:), NETph0(:)
	Double Precision :: t1, t2
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
	allocate(iv(inum_dummy))
	allocate(NELph0(NBph), NETph0(NBph))
!
!---( 物理定数, 材料定数および諸パラメータ )---
!
	call param
!
!---( 初期条件の設定 )---
!
	call initia(t,Elec,iv)  ! return with t=0.
!		call energy_dist(Elec)                    !  初期エネルギー分布出力
!		call out(t,Elec,iv)                         !  debug
	do i=1,NBph
		NELph0(i)=0.d0
		NETph0(i)=0.d0
	end do
!
!---( 多粒子モンテカルロ計算 )---
!
	call monte_carlo(Elec,iv,NELph0,NETph0,NBph,DEph)
!
!---( 粒子の最終エネルギー分布 )---
	call energy_dist(Elec)
!
!---( フォノンの頻度分布 )---
	do i=1,NBph
		NELph0(i)=NELph0(i)*(1.0d17/inum_dummy)*(1/SimTime)
		NETph0(i)=NETph0(i)*(1.0d17/inum_dummy)*(1/SimTime)
		write (*,*) i, dble(i)*DEph,NELph0(i),NETph0(i)
		write (8,*) i, dble(i)*DEph,NELph0(i),NETph0(i)
	end do
!      
!---( ファイルのクローズ )---
!
	call cpu_time( t2 )
	write(*,*) 'cpu time:', t2-t1, 'seconds.'
	write(8,*) 'cpu time:', t2-t1, 'seconds.'
	deallocate(Elec,iv)
	close(5)
	close(8)
	stop
end program main
