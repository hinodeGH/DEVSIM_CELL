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
!===( 乱数発生関数 )=== 
! 2018.08.16 rnd() -> system random number:call random_number()
! 2018.08.23 rv04  -> Mersenne random number: call grnd()
! 2018.08.25 range mt19937 -> mt19937rv0a:[0 1)  -> (0 1]  for log(grnd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
module random_number
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
		grnd = (dble(y) + 2.0d0 ** 32) / (2.0d0 ** 32) 
	else 
		grnd = dble(y) / (2.0d0 ** 32) 
	endif 
!
	grnd = 1.0d0 - grnd  !  range [0 1)  -> (0 1]  for log(grnd)
!      
end function grnd 
!
!........................................................................................
! 
Double Precision function gasdev() !gaussian random number:Box Muller's method
	implicit none
	integer :: iset=0   ! 0 setting is only once at first function call
	Double Precision gset
	Double Precision  fac,rsq,v1,v2
!
	if  (iset == 0) then
		v1=2.0d0*grnd()-1.0d0
		v2=2.0d0*grnd()-1.0d0
		rsq=v1*v1+v2*v2
		do while (rsq >= 1.0 .OR. rsq == 0.0)
			v1=2.0d0*grnd()-1.0d0
			v2=2.0d0*grnd()-1.0d0
			rsq=v1*v1+v2*v2
		end do
		fac=sqrt(-2.0*log(rsq)/rsq)
		gset=v1*fac
		iset=1
		gasdev = v2*fac
	else
		iset=0
		gasdev = gset
	end if
end function gasdev
!
end module random_number
!========================================================================================
!****************************************************************************************
module Pop_sub_program_variables    !----- 共通変数 -----
	implicit none
	Double Precision,parameter::QMIN = 1.0              ! // nearly zero compared to QMAX
	Double Precision,parameter::QMAX = 1.1569113068e+08 !  2*pi/asi
!    Double Precision,parameter::asi = 5.431e-08         ! // lattice constant, cm
	Double Precision,parameter::Rws = 2.1224e-8         !  asi*pow(3./(16.*pi), 1./3.)
	Double Precision,parameter::alpha = 0.481575201309447
!                                           Egap = 1.1756 - 8.8131e-05*T - 2.6814e-07*T*T;
!		                                     alpha = 0.5*1.1250281/Egap; 
	Double Precision,parameter::mdos = 1.86314779783618E-16      ! pow(ml*mt*mt, 1./3.) eV*s^2/cm^2
!     ml = 0.9163*m0;               // 5.2097e-16 eV*s^2/cm^2
!     mt = 0.1880*m0*1.1756/Egap;   // 1.1220e-16 eV*s^2/cm^2
	integer,parameter:: LA=1
	integer,parameter:: TA=2
    integer,parameter:: LO=3
    integer,parameter:: TO=4
!
    Double Precision,parameter::pi = 3.14159265        !π
    Double Precision,parameter::rho = 1.4516e+12       ! density (eV*s^2/cm^5) (2.329 g/cm^3) 
    Double Precision,parameter::kb = 8.61734318e-05    ! eV/K
    Double Precision,parameter::hbar = 6.58211915e-16  ! eV*s
!
    Double Precision,parameter::DLA = 6.39   ! debug 13.0= 6.39 x ~2
    Double Precision,parameter::DTA = 3.01   ! debug 6.0=3.01 x ~2
    Double Precision,parameter::Temp = 300.0
    Double Precision,parameter::Absorption = -1.0
    Double Precision,parameter::EMS = 1.0
!			  // f- and g-phonon norm, to be multiplied by QMAX
    Double Precision, parameter:: QGp = 0.3           ! @ Pop-default
    Double Precision, parameter:: QFp = 1.022252415 ! sqrt(1. + 2.*0.15*0.15)
    
    Double Precision,parameter::LTAO(4,3) = &
                         &  reshape( (/-0.196e-2, -0.250e-2, -0.156e-2,     0.137e-2, &
                         &              9.00e5,    5.33e5,    0.0000,      -2.91e5,   &
                         &              0.0,       0.0,       9.87682e13,   1.027e14 /),(/4,3/) )
!// Universal Constants from http://physics.nist.gov/cuu/Constants
!!#define pi   3.14159265
!#define pi2  6.28318531
!#define echarge   -1.0             // electron charge 
!#define ecoulom   -1.60217653e-19  // electron charge in Coulombs 
!!#define kb     8.61734318e-05      // eV/K 
!!#define hbar   6.58211915e-16      // eV*s 
!#define m0     5.68562975e-16      // eV*s^2/cm^2
!#define eps0   5.52634972e+05      // e/V/cm permittivity of free space

!// material constants for silicon
!#define asi    5.431e-08        // lattice constant, cm 
!!#define rho    1.4516e+12       // density (eV*s^2/cm^5) (2.329 g/cm^3) 
!#define epssi  6.46582917e+06   // e/V/cm permittivity of Si 11.7*eps0
    
end module Pop_sub_program_variables
!****************************************************************************************
module Pop_sub_programs
	use Pop_sub_program_variables
	use Tom_global_variables
	use random_number
	implicit none
	contains
!--------------- 以下に次のサブルーチン ------------
!    trapz ! 数値積分の足し算
!    get_wq ! phonon q -> w
!    get_nq ! phonon occupation
!    get_Mq ! (3.8) 積分式の被積分関数
!    fabs_cost  ! (3.12) cos(φ)式の絶対値
!    find_root  ! |cos(φ)|-1=0の根を探索
!    acrate3    ! E毎に(3.8) 積分
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
		if (LTtemp == LA) then
			LTtemp = LO                !  LA outside BZ extends into LO
			qTemp = 2.*QMAX - qTemp       ! and flip it back to first BZ
		end if
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
!						/* if (xx < 3.5) return (1./xx -0.5 +xx/12. -xx*xx*xx/720.);
!						else return exp(-xx); */
	get_nq = (1./(exp(xx)-1.))
end function get_nq
!========================================================================================
Double Precision function get_Mq(q, lt, aed)
	implicit none
	Double Precision, intent(inout)::q, aed
	integer, intent(in)::lt
	Double Precision gq
!	
	if (q > QMAX) q = 2.*QMAX-q               !    // returns to the first BZ
!					#ifdef POP
		gq = 3.*(sin(q*Rws)-q*Rws*cos(q*Rws))/((q*Rws)**3)
!					#else
!					gq = 1.0             !  // overlap integral
!					#endif
	get_Mq = (q*q*gq*gq*(get_nq(q,lt)+0.5+0.5*aed)/get_wq(q,lt))
end function get_Mq
!========================================================================================
Double Precision function fabs_cost(q, k, E, lt, aed)
	implicit none
	Double Precision, intent(in)::q, k, E, aed
	integer, intent(in)::lt
	Double Precision wq
!					 |cos(t)| = |+/-q/(2k) + m*w(q)/(hbar*k*q)|   +em/-ab
	wq = get_wq(q,lt)
!					aed = -1 for Absorption and +1 for EMS
	fabs_cost = abs( aed*q/(2.0*k) + mdos*(1.0 + alpha*(2.0*E - aed*hbar*wq))*wq/(hbar*k*q) )
!
end function fabs_cost
!========================================================================================
Double Precision function find_root(q1, q2, k, E, lt, aed)
	implicit none
	Double Precision, intent(inout)::q1, q2
	Double Precision, intent(in):: k, E, aed
	integer, intent(in)::lt
	Double Precision f1, f2, fs, q, eps
!					 find root of |cos(t)| - 1 in interval q1..q2 */
	q=q2
	eps = (q2-q1)/40.               ! 20.
	f1 = fabs_cost(q1,k,E,lt,aed)
	f2 = fabs_cost(q2,k,E,lt,aed)
	fs = (f2-f1)/abs(f2-f1)
	do while (q2-q1 > eps) 
		q = 0.5*(q1+q2)
		if (fs*(fabs_cost(q,k,E,lt,aed)- 1.) > 0.) then
			q2 = q
		else
			q1 = q
		end if
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
	Double Precision, intent(inout)::E, aed
!					 *****  aed = -1. for absorption,  +1. for emission
	integer, intent(in)::lt
	Double Precision,allocatable ::q(:), integ(:)
	Double Precision q1, q2, dq, k, acost, C3, G, DP
	integer i, i1, i2, nloc, N, MULT
!                        /******  ******/
	N=200              !  100
	MULT=20            !  10
!      
	allocate(q(MULT*N))
	allocate(integ(MULT*N))
!						 /* k is the electron wavevector in Herring-Vogt space */
	k = sqrt(2.*mdos*E*(1.+alpha*E))/hbar
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
end module Pop_sub_programs  
!========================================================================================
!****************************************************************************************   
module Tom_global_variables   !----- 共通変数 -----
	implicit none
	Double Precision dt,tem,fx
	Double Precision hwo,de
!	Double Precision :: swk(17,2400)=0.0d0
!	Double Precision :: swk(17,4000)=0.0d0
	Double Precision :: swk(17,3000)=0.0d0
	Double Precision amd      
	Double Precision gm                 
	Double Precision smh,hhml,hhmt         
	Double Precision af,af2,af4          
	Double Precision tm(3)                
	Double Precision hm(3)                
	Double Precision qh                   
	Double Precision bktq                  
	Double Precision ef
	integer inum,jtl
	integer iemax
	integer, parameter::PX=5,PY=6,PZ=7,EE=8,TS=9,XXX=1,YYY=2
	Double Precision pvx(6), pvy(6), pvz(6)  ! 6 valley p offset
	Double Precision mx(6), my(6), mz(6)     ! mx,my,mz at 6 valleys
!
end module Tom_global_variables
!========================================================================================
module Tom_sub_programs
	use Pop_sub_program_variables  !!!!!!
	use Tom_global_variables
	use random_number
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
	read(5,*) tem      ! temperature [K] (typ.300)
	read(5,*) fx       ! electric field [KV/cm] 
!
	write(*,*) inum/1E6,' M particles'
	write(*,*) jtl, ' iteration'
	write(*,*) dt, ' sec'
	write(*,*) tem, ' K'
	write(*,*) fx, ' KV/cm'
!
	write(8,*) inum/1E6,' M particles'
	write(8,*) jtl, ' iteration'
	write(8,*) dt, ' sec'
	write(8,*) tem, ' K'
	write(8,*) fx, ' KV/cm'
!	
	fx=fx*1000.0d0   !  kV/cm -> V/cm
!
end subroutine data_input
!========================================================================================
!===( 散乱レート計算用関数 )===
   Double Precision function DOSe(Ed)
    implicit none
    Double Precision Ed,sgamma,dosconst
!    
    if (Ed > 0.0) then
	  sgamma = sqrt(Ed*(1.+af*Ed)) 
!			// includes 2x for spin, else it would be 4x in the denominator
	  dosconst = 2.0*pi*pi*hbar*hbar*hbar
	  DOSe = ((2.0*amd)**1.5)*sgamma*(1. + af2*Ed)/dosconst
    else
      DOSe=0.0
    end if
  end function DOSe
!========================================================================================
!===( 物理定数, 材料定数および諸パラメータ )===
!
subroutine param
	implicit none      
	Double Precision aml, amt, eg, eps
	Double Precision rou, sv, cl, z2, da, dof, dog, bkt
!!      Double Precision amd, amc
	Double Precision amc      
	Double Precision wo, no, dos, aco, ofe, ofa, oge, oga
	Double Precision ei, sei, sef, eff
	Double Precision eee
	integer ie, i
!
	Double Precision,parameter::pi = 3.14159265        !π
!
	Double Precision,parameter::bk = 8.61734318e-05     ! eV/K
	Double Precision,parameter::q = 1.0d0               ! e
	Double Precision,parameter::hbar  = 6.58211915e-16  ! eV*s 
	Double Precision,parameter::ep0 = 5.52634972e+05    ! e/V/cm    !真空誘電率
	Double Precision,parameter::am0 = 5.68562975e-16    ! eV*s^2/cm^2  !electron mass
!      
!!!      remind   LA=1, TA=2, LO=3, TO=4
	Double Precision::Absorption = -1.0, EMS = 1.0
	Double Precision,parameter::PQMAX = 7.6149280673e-08 ! hbar*QMAX(Wigner-Seitz cell radius)
!  /************************************************************************/
!  /**************** INTERVALLEY SCATTERING (LO+TO+LA+TA) ******************/
!  /************************************************************************/
	Double Precision, parameter :: Zf=4.0         ! // f-scat degeneracy of final valley
	Double Precision, parameter :: Zg=1.0         ! // g-scat degeneracy of final valley
	Double Precision, parameter :: DTAf = 0.6e8   ! // 18 meV
	Double Precision, parameter :: DLAf = 3.5e8   ! // 50 meV
	Double Precision, parameter :: DTOf = 1.5e8   ! // 57 meV
	Double Precision, parameter :: DTAg = 0.7e8   ! // 10 meV
	Double Precision, parameter :: DLAg = 1.0e8   ! // 19 meV
	Double Precision, parameter :: DLOg = 6.0e8   ! // 63 meV
	Double Precision scatOconst, Nq
	Double Precision ETAf, ELAf, ETOf, ETAg, ELAg, ELOg
!    
!---( 電子の縦有効質量, 横有効質量, 禁制帯幅, 誘電率 )---
!
	aml  = 0.916*am0 ! 縦有効質量
	amt  = 0.196*am0 ! 横有効質量     
	eg   = 1.12      ! 禁制帯幅[eV]
	eps  = 11.7*ep0  ! Si誘電率
!
!---( フォノン散乱の諸パラメータ )---
!
!!      rou  = 2329.0 !Si比重
!!      sv   = 9040.0 !Si中の音速
	rou  = 1.4516e+12       ! density (eV*s^2/cm^5) (2.329 g/cm^3) !Si比重
!     PQMAX = 7.6149280673e-08 ! hbar*QMAX(Wigner-Seitz cell radius)
!!      sv   = 9040.0 !Si中の音速
!! 不要？      1.1e+05     // cm/s, optical phonon velocity
!!      cl   = rou*sv*sv    ! 弾性定数
!!  sv cl  不要？
	z2   = 2.0        ! X谷の個数
!!      da   = 9.0*q    ! 谷内音響フォノンの変位ポテンシャル [eV] after TXT A.1.2
!!      dof   =4.0e10*q ! 谷間f光学フォノンの変位ポテンシャル[eV/m] after TXT A.1.2
!       dog  = 6.0e10*q ! 谷間g光学フォノンの変位ポテンシャル[eV/m] after TXT A.1.2
!!      da   = 9.0        ! eV unit 
	dof   =4.0e8      ! eV/cm unit
	dog  = 1.1e9      ! eV/cm unit
	hwo  = 0.06       ! 谷間fg光学フォノンエネルギー[eV]
!
	bkt  = bk*tem
	bktq = bkt/q
	qh   = q/hbar     ! h->hbar
!
!---( バンドの非等方性 )---
!
	amd  = (aml*amt*amt)**(1.0/3.0) !相乗平均 effective m 0.328?
	amc  = 3.0/(1.0/aml+2.0/amt)!逆数相加平均 effective m at B edge 0.266?
	tm(1)=sqrt(aml/amd) !縦 1.67
	tm(2)=sqrt(amt/amd) !横 0.773
	tm(3)=sqrt(amt/amd) !横 0.773
	smh  = sqrt(2.*amd*q)/hbar    ! h->hbar
	hhml = hbar/aml/q*hbar/2.0 !縦 ! h->hbar
	hhmt = hbar/amt/q*hbar/2.0 !横 ! h->hbar
	hm(1)= hbar/aml !縦 ! h->hbar
	hm(2)= hbar/amt !横 ! h->hbar
	hm(3)= hbar/amt !横 ! h->hbar
!
!---( バンドの非放物線性 )---
!
	af   = (1.0-amc/am0)**2/eg !α~0.5
	af2  = 2.0*af
	af4  = 4.0*af
!
	wo     = hwo*q/hbar ! O-phonon angular frequency [J -> rad/s]  h -> hbar ?
	no     = 1.0/(exp(hwo/bktq)-1.0) !probability
!
	dos    = (2.0*amd*q)**(3.0/2.0)/(4.0*pi**2*hbar**3) !DOS-factor  h->hbar
!!      aco    = 2.0*pi*da/q*da*bktq/h*q/cl !音響フォノン
	ofe    = pi*dof/wo*dof/rou/q*(no+1.0) !光学f生成
	ofa    = ofe*no/(1.0+no) !光学f消滅
	oge    = pi*dog/wo*dog/rou/q*(no+1.0) !光学g生成
	oga    = oge*no/(1.0+no) !光学g消滅
!!      write (*,*) wo,no,hwo,bktq,dos,ofe,ofa,oge,oga   ! debug
!!      write (8,*) wo,no,hwo,bktq,dos,ofe,ofa,oge,oga   ! debug
!
!---( valley offset values )---
!
	do i=1,6
		pvx(i) = 0.0d0
        pvy(i) = 0.0d0
        pvz(i) = 0.0d0
	end do
	pvx(1) = -0.85*PQMAX  ! -X axis valley
	mx(1) = aml
	my(1) = amt
	mz(1) = amt
	pvx(2) = +0.85*PQMAX  ! +X axis valley
	mx(2) = aml
	my(2) = amt
	mz(2) = amt
	pvy(3) = -0.85*PQMAX  ! -Y axis valley
	mx(3) = amt
	my(3) = aml
	mz(3) = amt
	pvy(4) = +0.85*PQMAX  ! +Y axis valley
	mx(4) = amt
	my(4) = aml
	mz(4) = amt
	pvz(5) = -0.85*PQMAX  ! -Z axis valley
	mx(5) = amt
	my(5) = amt
	mz(5) = aml
	pvz(6) = +0.85*PQMAX  ! +Z axis valley
	mx(6) = amt
	my(6) = amt
	mz(6) = aml
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
!						write(8,*) ELOg,ETOf,ELAg,ETAg,ELAf,ETAf         ! debug
!
	de=0.0005   ! energy step  ! Tom(2meV)->Pop(0.5meV)
!	iemax=2400  !  Tom(0-2eV;1000div)->Pop(0-1.2eV;2400div)
!	iemax=4000
	iemax=3000
!
	do ie=1,iemax
		ei=de*float(ie) ! initial energy
		sei=sqrt(ei)
!-----( 非有極性光学フォノン散乱 )--- scattering type 1->6; 2->7; 3->8; 4->9; 5->5
!!!        eff=ei-hwo ! phonon emission
!!!        if (eff > 0.0) then
!!!           sef=sqrt(eff*(1.0+af*eff))
!!!           swk(6,ie) = oge*sef*dos*(1.0+2.0*af*eff) !光学g emission
!!!        else
!!!            swk(6,ie) = 0.
!!!        end if
!
!!!        eff=ei+hwo ! phonon absorption
!!!        sef=sqrt(eff*(1.0+af*eff))
!!!        swk(7,ie) = oga*sef*dos*(1.0+2.0*af*eff) !光学g absorption
!
!!!        eff=ei-hwo
!!!        if (eff > 0.0) then
!!!           sef=sqrt(eff*(1.0+af*eff))
!!!           swk(8,ie) = 2.0*z2*ofe*sef*dos*(1.0+2.0*af*eff)!光学f ems
!!!        else
!           swk(3,ie)=swk(2,ie)
!!!           swk(8,ie)=0.
!!!        end if
!
!!!        eff=ei+hwo
!!!           sef=sqrt(eff*(1.+af*eff))
!!!           swk(9,ie) = 2.0*z2*ofa*sef*dos*(1.0+2.0*af*eff) !光学f abs
!=====[ 光学フォノン散乱 Pop }===== scattering type 6 - 9
! 			 // the value of the overlap integral for inter-valley scattering may
!  			 // be included in the coupling constants (Jacoboni 1983, p. 672)
!
		scatOconst = pi/(2.*rho/hbar)
! 			  /************* optical (63 meV LO-g) phonon emission ****************/
		Nq = 1./(exp(ELOg/(kb*Temp))-1.)
		swk(6,ie) = Zg*DLOg*DLOg*(Nq+1.)*DOSe(ei-ELOg)*scatOconst/ELOg
!					  /************ optical (63 meV LO-g) phonon absorption ***************/
		swk(7,ie)  = Zg*DLOg*DLOg*Nq*DOSe(ei+ELOg)*scatOconst/ELOg
!
		Nq = 1./(exp(ETOf/(kb*Temp))-1.)
!					  /************* optical (59 meV TO-f) phonon emission ****************/
!    		swk(8,ie) = Zf*DTOf*DTOf*(Nq+1.)*DOSe(ei-ETOf+DESIGEX)*scatOconst/ETOf
		swk(8,ie) = Zf*DTOf*DTOf*(Nq+1.)*DOSe(ei-ETOf)*scatOconst/ETOf
!					  /************ optical (59 meV TO-f) phonon absorption ***************/
!    		swk(9,ie) = Zf*DTOf*DTOf*Nq*DOSe(ei+ETOf+DESIGEX)*scatOconst/ETOf
		swk(9,ie) = Zf*DTOf*DTOf*Nq*DOSe(ei+ETOf)*scatOconst/ETOf
!    
!=====[ 音響フォノンValley間散乱 Pop }===== scattering type 10 - 17
!					  /************* acoustic (12 meV TA-g) phonon emission ****************/
		Nq = 1./(exp(ETAg/(kb*Temp))-1.) 
		swk(10,ie)  = Zg*DTAg*DTAg*(Nq+1.)*DOSe(ei-ETAg)*scatOconst/ETAg
!					  /************ acoustic (12 meV TA-g) phonon absorption ***************/
		swk(11,ie)  = Zg*DTAg*DTAg*Nq*DOSe(ei+ETAg)*scatOconst/ETAg
!
!					  /************* acoustic (47 meV LA-f) phonon emission ****************/
		Nq = 1./(exp(ELAf/(kb*Temp))-1.)
!    		swk(12,ie)  = Zf*DLAf*DLAf*(Nq+1.)*DOSe(ei-ELAf+DESIGEX)*scatOconst/ELAf
		swk(12,ie)  = Zf*DLAf*DLAf*(Nq+1.)*DOSe(ei-ELAf)*scatOconst/ELAf
!					  /************ acoustic (47 meV LA-f) phonon absorption ***************/
!    		swk(13,ie)  = Zf*DLAf*DLAf*Nq*DOSe(ei+ELAf+DESIGEX)*scatOconst/ELAf
		swk(13,ie)  = Zf*DLAf*DLAf*Nq*DOSe(ei+ELAf)*scatOconst/ELAf
!					  /************* acoustic (19 meV TA-f) phonon emission ****************/
!!   		 Nq = 1./(exp(ETAf/(kb*tempTAf))-1.)
		Nq = 1./(exp(ETAf/(kb*Temp))-1.)
!    		swk(14,ie)  = Zf*DTAf*DTAf*(Nq+1.)*DOSe(ei-ETAf+DESIGEX)*scatOconst/ETAf
		swk(14,ie)  = Zf*DTAf*DTAf*(Nq+1.)*DOSe(ei-ETAf)*scatOconst/ETAf
!					  /************ acoustic (19 meV TA-f) phonon absorption ***************/
!!    		swk(15,ie)  = Zf*DTAf*DTAf*Nq*DOSe(ei+ETAf+DESIGEX)*scatOconst/ETAf
		swk(15,ie)  = Zf*DTAf*DTAf*Nq*DOSe(ei+ETAf)*scatOconst/ETAf
!					  /************* acoustic (18 meV LA-g) phonon emission ****************/
		Nq = 1./(exp(ELAg/(kb*Temp))-1.)
		swk(16,ie)  = Zg*DLAg*DLAg*(Nq+1.)*DOSe(ei-ELAg)*scatOconst/ELAg
!					  /************ acoustic (18 meV LA-g) phonon absorption ***************/
		swk(17,ie)  = Zg*DLAg*DLAg*Nq*DOSe(ei+ELAg)*scatOconst/ELAg
!
!-----( 音響フォノン散乱 )--- scattering type 5->(2) 5
!=====[ 音響フォノン散乱 Pop }===== scattering type 2 - 5
!
		swk(2,ie)= acrate3(ei, LA, EMS)         ! LA intra emission
		swk(3,ie)= acrate3(ei, LA, Absorption)  ! LA intra absorption
		swk(4,ie)= acrate3(ei, TA, EMS)         ! TA intra emission
		swk(5,ie)= acrate3(ei, TA, Absorption)  ! TA intra absorption
!
        end do
!        
!!----- debug scattering rate begin ------
!	do  ie=1,iemax
!		eee=de*float(ie)    
!		write (*,*) eee,swk(2,ie),swk(3,ie),swk(4,ie),swk(5,ie),swk(6,ie),swk(7,ie),swk(8,ie),swk(9,ie), &
!				& swk(10,ie),swk(11,ie),swk(12,ie),swk(13,ie),swk(14,ie),swk(15,ie),swk(16,ie),swk(17,ie)
!		write (8,*) eee,swk(2,ie),swk(3,ie),swk(4,ie),swk(5,ie),swk(6,ie),swk(7,ie),swk(8,ie),swk(9,ie), &
!				& swk(10,ie),swk(11,ie),swk(12,ie),swk(13,ie),swk(14,ie),swk(15,ie),swk(16,ie),swk(17,ie)
!	end do
!!------ debug scattering rate end -------
!
!---( 散乱レートの和の計算 )---
!     
	gm=0.
	do  ie=1,iemax
		do i=2, 17
			swk(1,ie)=swk(1,ie)+swk(i,ie)
		end do  
		if (swk(1,ie) > gm) gm = swk(1,ie)
	end do
	do  ie=1,iemax
		swk(1,ie)=gm-swk(1,ie)
		do i=2, 17
			swk(i,ie)=swk(i-1,ie)+swk(i,ie)
		end do
		do i=1, 17
			swk(i,ie)=swk(i,ie)/gm
		end do
	end do     
!----- debug scattering rate begin ------
	write (*,*) 'total scattering rate=',gm
	 write (8,*) 'total scattering rate=',gm
!!      do  ie=1,iemax
!!        eee=de*float(ie)    
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
	Double Precision,parameter::pi = 3.14159        !π
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
	call sgrnd(1)                    !乱数 seed=1
	MBdist=1
	do  n=1,inum
!--------------------------------------------- begin isotropic Maxwell-Boltzman(Pop)
		if (MBdist == 1) then
			kx = smh*sqrt(bktq/2.0)*gasdev()    ! Local kx
			ky = smh*sqrt(bktq/2.0)*gasdev()    ! Local ky
			kz = smh*sqrt(bktq/2.0)*gasdev()    ! Local kz
!--------------------------------------------- begin isotropic initial p (Tom)
		else
			rn = grnd()
			ei  = -bktq*log(rn)*1.5d0    !--- give E along exp{-E/(3/2kT)} ~39meV @300K
			ak  = smh*sqrt(ei*(1.0+af*ei)) !--- E -> k
			rn = grnd()
			cb  = 1.0-2.0*rn ! cosθ
			sb  = sqrt(1.0-cb*cb) ! sinθ
			rn = grnd()
			fai = 2.0*pi*rn ! Φ
			sf  = sin(fai)
			kx = ak*cb*sf    ! Local kx
			ky = ak*sb*sf    ! Local ky
			kz = ak*cos(fai) ! Local kz
		end if
!--------------------------------------------- end isotropic initial k
		rn = grnd()
		Elec(n,TS) = -log(rn)/gm         ! mean free time along exp(-t/gm)
		Elec(n,XXX) = 0.0 ! x
		Elec(n,YYY) = 0.0 ! y
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
		Elec(n,EE)=(sqrt(1.0+af4*gk)-1.0)/af2  ! Tomizawa (1.7)
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
	Double Precision padjust, Elost
	integer n, ivv
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph
	Double Precision, intent(in)::DEph
!	Double Precision, Parameter :: eee_Min_Limit=1.0e-4 	! 0.4->0.1meV avoid infinit loop at do while in scat
!	Double Precision, Parameter :: eee_Max_Limit=2.0 		! 2eV avoid infinit loop at do while in scat
!	Double Precision, Parameter :: eee_Max_Limit=1.5
!
	Elost = 0.0
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
!		if (eee < eee_Min_Limit) then	   				! ? kjh avoid infinite loop at do while in scat
!			padjust = sqrt((eee_Min_Limit*(1.+af*eee_Min_Limit))/(eee*(1.+af*eee)))
!			kx = pvx(ivv)/hbar + (kx - pvx(ivv)/hbar)*padjust
!			ky = pvy(ivv)/hbar + (ky - pvy(ivv)/hbar)*padjust
!			kz = pvz(ivv)/hbar + (kz - pvz(ivv)/hbar)*padjust
!				Elost = Elost + (eee - eee_Min_Limit)     !				   // add up "got" energy
!			eee=eee_Min_Limit   
!		end if
!		if (eee > eee_Max_Limit) then   				! ? kjh avoid infinite loop at do while in scat
!			padjust = sqrt((eee_Max_Limit*(1.+af*eee_Max_Limit))/(eee*(1.+af*eee)))
!			kx = pvx(ivv)/hbar + (kx - pvx(ivv)/hbar)*padjust
!			ky = pvy(ivv)/hbar + (ky - pvy(ivv)/hbar)*padjust
!			kz = pvz(ivv)/hbar + (kz - pvz(ivv)/hbar)*padjust
!				Elost = Elost + (eee - eee_Max_Limit)     !				   // add up "lost" energy
!			eee=eee_Max_Limit
!		end if
!
		t1=t
		do while (tss <= tdt)          ! mean free time > elapsed time -> drift
			tau=tss-t1                   ! time interval before scattering
			call drift(tau,kx,ky,kz,eee,x,y,ivv)
			call scat(kx,ky,kz,eee,ivv,NELph0,NETph0,NBph,DEph)
			t1=tss
			rn = grnd()
			tss=t1-log(rn)/gm         ! re-new mean free time
		end do ! ---- while end
       !
		tau=tdt-t1
		call drift(tau,kx,ky,kz,eee,x,y,ivv)
!
!!				tentatively skip after drift and/or scat
!		if (eee < eee_Min_Limit) then	   				! ? kjh avoid infinite loop at do while in scat
!			padjust = sqrt((eee_Min_Limit*(1.+af*eee_Min_Limit))/(eee*(1.+af*eee)))
!			kx = pvx(ivv)/hbar + (kx - pvx(ivv)/hbar)*padjust
!			ky = pvy(ivv)/hbar + (ky - pvy(ivv)/hbar)*padjust
!			kz = pvz(ivv)/hbar + (kz - pvz(ivv)/hbar)*padjust
!				Elost = Elost + (eee - eee_Min_Limit)     !				   // add up "got" energy
!			eee=eee_Min_Limit   
!		end if
!		if (eee > eee_Max_Limit) then   				! ? kjh avoid infinite loop at do while in scat
!			padjust = sqrt((eee_Max_Limit*(1.+af*eee_Max_Limit))/(eee*(1.+af*eee)))
!			kx = pvx(ivv)/hbar + (kx - pvx(ivv)/hbar)*padjust
!			ky = pvy(ivv)/hbar + (ky - pvy(ivv)/hbar)*padjust
!			kz = pvz(ivv)/hbar + (kz - pvz(ivv)/hbar)*padjust
!				Elost = Elost + (eee - eee_Max_Limit)     !				   // add up "lost" energy
!			eee=eee_Max_Limit
!		end if
!
		Elec(n,PX) = hbar*kx    ! Global px
		Elec(n,PY) = hbar*ky    ! Global py
		Elec(n,PZ) = hbar*kz    ! Global pz
		Elec(n,EE) = eee
		Elec(n,TS) = tss
		Elec(n,XXX) = x
		Elec(n,YYY) = y
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
	dkx=qh*fx*tau ! dk=q/h・F・τ
!
	cp  = hm((ivv-1)/2+1)*tau               ! m/h・τ

	gk = eee*(1.0+af*eee)
	sq = sqrt(1.0+af4*gk)      !--- anisotropic & non-parabolic
	x = x+cp*((kx-pvx(ivv)/hbar)+0.5*dkx)/sq   !--- dx=h/m・τ(kx +0.5・dkx)/sq
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
	eee=(sqrt(1.0+af4*gk)-1.0)/af2            ! Tomizawa (1.7)
!
end subroutine drift
!========================================================================================
!===( 散乱過程の計算 )===
!
subroutine scat(kx,ky,kz,eee,ivv,NELph0,NETph0,NBph,DEph)  ! Global kx, ky, kz
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
	Double Precision padjust
	Double Precision, Parameter :: eee_Min_Limit=1.0e-4 	! 0.4->0.1meV avoid infinit loop at do while
!	Double Precision, Parameter :: eee_Max_Limit=2.0 		! 2eV avoid infinit loop at do while
	Double Precision, Parameter :: eee_Max_Limit=1.5
!!		Double Precision qdum, edum, pprime, kprime
!!		Double Precision cbet, sbet, cgam, sgam, cgamp, cgampp, eta, ceta, seta, phi, sphi, cphi
!!		integer cnt
!!		integer BZCLIP
!!		Double Precision calp, salp, bet
!!		Double Precision,parameter::PQMAX = 7.6149280673e-08 ! hbar*QMAX(Wigner-Seitz cell radius)
!
!---( 散乱前のエネルギーの計算 )---
!
		if (eee < eee_Min_Limit) then	   				! ? kjh avoid infinite loop at do while in scat
			padjust = sqrt((eee_Min_Limit*(1.+af*eee_Min_Limit))/(eee*(1.+af*eee)))
			kx = pvx(ivv)/hbar + (kx - pvx(ivv)/hbar)*padjust
			ky = pvy(ivv)/hbar + (ky - pvy(ivv)/hbar)*padjust
			kz = pvz(ivv)/hbar + (kz - pvz(ivv)/hbar)*padjust
!				Elost = Elost + (eee - eee_Min_Limit)     !				   // add up "got" energy
			eee=eee_Min_Limit   
		end if
		if (eee > eee_Max_Limit) then   				! ? kjh avoid infinite loop at do while in scat
			padjust = sqrt((eee_Max_Limit*(1.+af*eee_Max_Limit))/(eee*(1.+af*eee)))
			kx = pvx(ivv)/hbar + (kx - pvx(ivv)/hbar)*padjust
			ky = pvy(ivv)/hbar + (ky - pvy(ivv)/hbar)*padjust
			kz = pvz(ivv)/hbar + (kz - pvz(ivv)/hbar)*padjust
!				Elost = Elost + (eee - eee_Max_Limit)     !				   // add up "lost" energy
			eee=eee_Max_Limit
		end if
!
	gk = eee*(1.0+af*eee)
	ki = sqrt(2.0*amd*gk)/hbar   ! Local k
	ei=(sqrt(1.0+af4*gk)-1.0)/af2  ! Tomizawa (1.7)
	ie=int(ei/de)+1
	if (ie > iemax) ie=iemax   ! denomination to eee_Max_Limit ?
!
!---( 散乱機構の選択 )---
!
	r1 = grnd()
	if (r1 <= swk(1,ie)) then       ! self-scattering
!			no operation
	elseif (r1 <= swk(5,ie)) then   !  LA/TA intra emission/absorption
		if (r1 <= swk(2,ie)) then   !  LA intra emission
			LT=LA
			aed=+1.0d0
			qM=2.0*ki
		elseif (r1 <= swk(3,ie)) then   !  LA intra absorption
			LT=LA
			aed=-1.0d0
			Epmax = eee + hbar*get_wq(QMAX,LT) 
			kspmax = sqrt(2.*mdos*Epmax*(1.+af*Epmax))/hbar 
			qM = ki + kspmax
		elseif (r1 <= swk(4,ie)) then   !  TA intra emission
			LT=TA
			aed=+1.0d0  
			qM=2.0*ki            
		else   !    r1 <= swk(5,ie)     TA intra absorption
			LT=TA
			aed=-1.0d0
			Epmax = eee + hbar*get_wq(QMAX,LT)
			kspmax = sqrt(2.*mdos*Epmax*(1.0+af*Epmax))/hbar
			qM = ki + kspmax
		end if
!
		call final_state_intra_scatByLATA(kx,ky,kz,ki,eee,ivv,LT,aed,Epmax,kspmax,qM,NELph0,NETph0,NBph,DEph)
!
	elseif (r1 <= swk(6,ie)) then       !  O-g-emis
		LT=LO
		aed=+1.0d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
	elseif (r1 <= swk(7,ie)) then   !  O-g-abs
		LT=LO
		aed=-1.0d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
	elseif (r1 <= swk(8,ie)) then   !  O-f-emis
		LT=TO
		aed=+1.0d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph)
	elseif (r1 <= swk(9,ie)) then   !  O-f-abs
		LT=TO
		aed=-1.0d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph)
!
	else if (r1 <= swk(10,ie)) then   ! TA-g emission
		LT=TA
		aed=+1.0d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz 
!					  if (STATS) phonon_stats(EMS,ephon,&EphTA,NETph0,NAgrid,Agrid,ipx,ipy);
	else if (r1 <= swk(11,ie)) then   ! TA-g absorption
		LT=TA
		aed=-1.0d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
!					  if (STATS) phonon_stats(ABS,ephon,&EphTA,NETph0,NAgrid,Agrid,ipx,ipy);
	else if (r1 <= swk(12,ie)) then   ! LA-f emission
		LT=LA
		aed=+1.0d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz 	
!
!					  if (STATS) <18>{
!						if (ephon < (hbar*get_wq(QMAX,LA))) <20>{   // can also use q0 < QMAX
!						  phonon_stats(EMS,ephon,&EphLA,NELph0,NAgrid,Agrid,ipx,ipy);
!						}<20> else <19>{
!						  phonon_stats(EMS,ephon,&EphLO,NELph0,NOgrid,Ogrid,ipx,ipy);
!						}<19>
!					  }<18>
!
	else if (r1 <= swk(13,ie)) then   ! LA-f absorption
		LT=LA
		aed=-1.0d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
!
!				#ifndef NOISE
!					  if (STATS) <22>{
!						if (ephon < (hbar*get_wq(QMAX,LA))) <23>{   // can also use q0 < QMAX
!						  phonon_stats(ABS,ephon,&EphLA,NELph0,NAgrid,Agrid,ipx,ipy);
!						}<23> else <24>{
!						  phonon_stats(ABS,ephon,&EphLO,NELph0,NOgrid,Ogrid,ipx,ipy);
!						}<24>
!					  }<22>
!
	else if (r1 <= swk(14,ie)) then   ! TA-f emission
		LT=TA
		aed=+1.0d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
!
!				#ifndef NOISE
!					  if (STATS) phonon_stats(EMS,ephon,&EphTA,NETph0,NAgrid,Agrid,ipx,ipy);
!				#endif
	else if (r1 <= swk(15,ie)) then   ! TA-f absorption
		LT=TA
		aed=-1.0d0
		call final_state_inter_scat_f(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
!
!				#ifndef NOISE
!					  if (STATS) phonon_stats(ABS,ephon,&EphTA,NETph0,NAgrid,Agrid,ipx,ipy);
!				#endif

	else if (r1 <= swk(16,ie)) then   ! LA-g emission
		LT=LA
		aed=+1.0d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz 
!
!				#ifndef NOISE
!					  if (STATS) phonon_stats(EMS,ephon,&EphLA,NELph0,NAgrid,Agrid,ipx,ipy);
!				#endif

	else if (r1 <= swk(17,ie)) then   ! LA-g absorption
		LT=LA
		aed=-1.0d0
		call final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz 
!
!				#ifndef NOISE
!					  if (STATS) phonon_stats(ABS,ephon,&EphLA,NELph0,NAgrid,Agrid,ipx,ipy);
!				#endif
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
	Double Precision,parameter::PQMAX = 7.6149280673e-08 ! hbar*QMAX(Wigner-Seitz cell radius)
	Double Precision, intent(inout)::NELph0(:),NETph0(:)
	integer ,intent(in)::NBph
	Double Precision, intent(in)::DEph
	integer ie
!      
	C = qM*get_Mq(qM,LT,aed)*(1.0+af2*(eee-aed*hbar*get_wq(qM,LT))) 
!
	q0 = QMIN + (qM-QMIN)*grnd()
	ephon = hbar*get_wq(q0,LT)
	f0 = C*grnd()
	fx = q0*get_Mq(q0,LT,aed)*(1.0+af2*(eee-aed*ephon))
	acost = fabs_cost(q0,ki,eee,LT,aed)
	cnt=1									! debug
	do while ((f0 > fx) .OR. (acost > 1.0d0))
		q0 = QMIN + (qM-QMIN)*grnd()
		ephon = hbar*get_wq(q0,LT)
		f0 = C*grnd()
		fx = q0*get_Mq(q0,LT,aed)*(1.0+af2*(eee-aed*ephon))
		acost = fabs_cost(q0,ki,eee,LT,aed)
		cnt=cnt+1								! debug
!		if (mod(cnt,10000)==0) then			! debug
!			write(*,*) 'cnt=',cnt,'eee=',eee,' acost=',acost    ! debug
!			write(8,*) 'cnt=',cnt,'eee=',eee,' acost=',acost    ! debug
!		end if
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
		NELph0(ie)=NELph0(ie)+aed*1.0d0
	else if (LT==TA .OR. LT==TO) then
		NETph0(ie)=NETph0(ie)+aed*1.0d0
	end if
!
	qdum = q0
	edum = ephon
!        
!--- calculation of final state for LA/TA phonon emission/absorption (type=2-5)
!   /********* non-parabolicity from tomizawa eq. 1.11 page 8 *********/
	pprime = sqrt(2.0*amd*eee*(1.+af*eee)) !
!                  this is the longitudinal axis momentum component
!    			Pxyz = PX-1+(int)(iv[j]+1)/2;iv->Pxyz   ! 1,2->5;3,4->6;5,6->7
!				/******** intra-valley acoustic, NON-ISOTROPIC final angle ********/
!!			if (iscat <= 5) then !  <40>
	kprime = pprime/hbar  ! k=p/hbar
!				  /***** momentum conservation for angle beta b/w ki and kprime *****/
	cbet = (ki*ki+kprime*kprime-q0*q0)/(2.0*ki*kprime) ! low of cosine
!	cbet = ((ki/kprime)+(kprime/ki)-(q0/ki)*(q0/kprime))/2.0 !
!	 			 /****************** detect if angle is foobar-ed ******************/
!				   if (abs(cbet) > 1.) then
	if (abs(cbet) > 1.0 .AND. abs(cbet) <= 1.02) then
		write(*,*) '** WARNING 2  |cbet|>1;  cbet=',cbet
			write(*,*) 'ei=',eee+aed*ephon,'e=',eee,' ephon=',ephon,' ki=',ki,' kprime=',kprime,'&
						& q0=',q0,' cbet=',cbet,'QMIN=',QMIN,'qM=',qM,'LT ',LT    !debug
!	 		write(8,*) 'ei=',eee+aed*ephon,'e=',eee,' ephon=',ephon,' ki=',ki,' kprime=',kprime,'&
!						& q0=',q0,' cbet=',cbet,'QMIN=',QMIN,'qM=',qM    !debug
		cbet=abs(cbet)/cbet ! kjh 2% error neglect
	end if
	if (abs(cbet) > 1.02) then
		write(*,*) 'error abs(cbet) > 1.02'
		write(*,*) 'ki=',ki,' kprime=',kprime,' q0=',q0,' cbet=',cbet    !debug
		write(8,*) 'ki=',ki,' kprime=',kprime,' q0=',q0,' cbet=',cbet    !debug
!		 	 exit(-1)
		stop
	end if
	sbet = sqrt(1.-cbet*cbet)
!		          /*** gam = angle b/w ki and Z-axis ***/
	cgam = (hbar*kz-pvz(ivv))*sqrt(amd/mz(ivv))/(hbar*ki)
!			if (abs(cgam) > 1.) then
	if (abs(cgam) > 1.0 .AND. abs(cgam) <= 1.03) then
		write(*,*) '** WARNING 3  |cgam|>1;  cgam=',cgam	   !debug      
		cgam=abs(cgam)/cgam  ! 3% denomination by kjh
	end if
	if (abs(cgam) > 1.03) then
		write(*,*) 'error abs(cgam) > 1.03'
		write(*,*) 'cbet=',cbet,' sbet=',sbet,' cgam=',cgam    !debug
		write(8,*) 'cbet=',cbet,' sbet=',sbet,' cgam=',cgam    !debug
!	 			 exit(-1)
		stop
	end if
	sgam = sqrt(1.-cgam*cgam)
!			/*** gamp = angle b/w ki and Y-axis ***/
	cgamp = (hbar*ky-pvy(ivv))*sqrt(amd/my(ivv))/(hbar*ki)
!			/*** gampp = angle b/w ki and X-axis ***/
	cgampp = (hbar*kx-pvx(ivv))*sqrt(amd/mx(ivv))/(hbar*ki)
!			/*** angle b/w ki projection on X-Y plane and Y-axis ***/
	eta = atan((hbar*kx-pvx(ivv))*my(ivv)/   &
		  &	 ((hbar*ky-pvy(ivv))*mx(ivv)))
	ceta = cos(eta)
	seta = sin(eta)
!
!     begin   once for while
!			/*** rotation phi chosen at random ***/
	phi = 2.*pi*grnd()
	sphi = sin(phi) 
	cphi = cos(phi) 
!			/*** compute final state momentum components, like Jacoboni Rev. Mod.
!			 1983, eq. 3.62-3.64, but he has X and Z components mixed up ***/
	kz = ((cbet*cgam   +sbet*cphi*sgam)* &
	   & sqrt(mz(ivv)/amd)*pprime + pvz(ivv))/hbar
	ky = ((cbet*cgamp  -sbet*cphi*cgam*ceta -sbet*sphi*seta)* &
	   & sqrt(my(ivv)/amd)*pprime + pvy(ivv))/hbar
	kx = ((cbet*cgampp -sbet*cphi*cgam*seta +sbet*sphi*ceta)* &
	   & sqrt(mx(ivv)/amd)*pprime + pvx(ivv))/hbar
	cnt=1   !! debug
	BZCLIP=1      ! debug
!		end   once for while
!		do while ((abs(Elec[Pxyz][j]) > PQMAX) && (cnt < 200) && (BZCLIP)) ! <37>{
!			write(*,*) abs(kx),abs(ky),abs(kz), PQMAX/hbar  ! debug
!			write(8,*) abs(kx),abs(ky),abs(kz), PQMAX/hbar  ! debug
	do while ((MAX(abs(kx),abs(ky),abs(kz)) > PQMAX/hbar) .AND. (cnt < 200) .AND. (BZCLIP==1))  ! <37>{  
!			/*** rotation phi chosen at random ***/
		phi = 2.*pi*grnd()
		sphi = sin(phi)
		cphi = cos(phi)
!     		/*** compute final state momentum components, like Jacoboni Rev. Mod.
!			1983, eq. 3.62-3.64, but he has X and Z components mixed up ***/
		kz = ((cbet*cgam   +sbet*cphi*sgam)* &
		   & sqrt(mz(ivv)/amd)*pprime + pvz(ivv))/hbar
		ky = ((cbet*cgamp  -sbet*cphi*cgam*ceta -sbet*sphi*seta)* &
		   & sqrt(my(ivv)/amd)*pprime + pvy(ivv))/hbar
		kx = ((cbet*cgampp -sbet*cphi*cgam*seta +sbet*sphi*ceta)* &
		   & sqrt(mx(ivv)/amd)*pprime + pvx(ivv))/hbar
!      		/*** increment counter ***/
		cnt = cnt + 1
	end do    !  }<37> 
!		  // if BZCLIP = 0 this && is always FALSE so it only goes thru ONCE
!		  // if the counter exceeds 200 reps then make it isotropic
	if (cnt > 200) then !  <39>{          // note:  this never happens if BZCLIP = 0
!		  if (cnt > 200) write(*,*) '** WARNING (intra_scatByLATA) **  cnt =',cnt   ! debug
		write(*,*) '** WARNING 4 (intra_scatByLATA) **  cnt =',cnt
!			write(*,*) abs(kx),abs(ky),abs(kz), PQMAX/hbar,eee,q0  ! debug
!			write(8,*) abs(kx),abs(ky),abs(kz), PQMAX/hbar,eee,q0  ! debug
!			     /***** compute ISOTROPIC new angles if all else fails above ******/
!       -- begin    once for while
		calp = 1.0-2.0*grnd()
		salp = sqrt(1.0-calp*calp)
		bet = 2.0*pi*grnd()
!			    /********** compute new isotropic momentum components *************/
		kx = (pvx(ivv) + sqrt(mx(ivv)/amd)*pprime*salp*cos(bet))/hbar
		ky = (pvy(ivv) + sqrt(my(ivv)/amd)*pprime*salp*sin(bet))/hbar
		kz = (pvz(ivv) + sqrt(mz(ivv)/amd)*pprime*calp)/hbar
!       -- end      once for while
		cnt=1
!			do while ((fabs(Elec[Pxyz][j]) > PQMAX) && (BZCLIP))     ! <38>{
		do while ((MAX(abs(kx),abs(ky),abs(kz)) > PQMAX/hbar) .AND. (BZCLIP==1))   ! <38>{
			calp = 1.-2.*grnd()
			salp = sqrt(1.-calp*calp)
			bet = 2.*pi*grnd()
!			  /********** compute new isotropic momentum components *************/
			kx = (pvx(ivv) + sqrt(mx(ivv)/amd)*pprime*salp*cos(bet))/hbar
			ky = (pvy(ivv) + sqrt(my(ivv)/amd)*pprime*salp*sin(bet))/hbar
			kz = (pvz(ivv) + sqrt(mz(ivv)/amd)*pprime*calp)/hbar
			cnt=cnt+1
!				write(*,*) 'cnt3= ',cnt,' ivv3= ',ivv   ! debug
!				write(8,*) 'cnt3= ',cnt,' ivv3= ',ivv   ! debug
!				write(*,*) abs(kx),abs(ky),abs(kz), PQMAX/hbar  ! debug
!				write(8,*) abs(kx),abs(ky),abs(kz), PQMAX/hbar  ! debug
		end do   ! }<38> 
	end if     ! }<39>
!	}<40>
!--- end calculation of final state for LA/TA phonon emission/absorption (type=2-5)
end subroutine final_state_intra_ScatByLATA
!========================================================================================
subroutine final_state_inter_scat_g(kx,ky,kz,eee,ivv,LT,aed,NELph0,NETph0,NBph,DEph) ! Global kx, ky, kz
	implicit none
	Double Precision, intent(inout)::  kx,ky,kz,eee
!		Double Precision, parameter:: QGp=0.3  ! @ Pop-default
	Double Precision, parameter:: EMIN = 1.e-10   ! @ Pop-definition
	Double Precision,parameter::PQMAX = 7.6149280673e-08 ! hbar*QMAX(Wigner-Seitz cell radius)
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
!			#else
!			if (LT == LA) ephon = ELAg       ! 210-215*kb~18meV   @Pop-Default
!			else if (LT == TA) ephon = ETAg  ! 140*kb~12meV
!			else ephon = ELOg                ! 700-720*kb~61meV
!			end if
!			if ((aed > 0.) .AND. ((eee - ephon) <= EMIN) then
!				write (*,*) '(eee - ephon) <= EMIN'
!			stop
!			end if
!			qphon = QGp*QMAX;
!			#endif
!
	call change_equiv_valley(ivv)
!			// new longitudinal axis momentum component
!			Pxyz = PX-1+(int)(ivv+1)/2;
!
!			#ifdef POP
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
			pp = sqrt(2.0*amd*Ep*(1.0+af*Ep))
!				// choose the set of random angles
			calp = 1.0-2.0*grnd()
			salp = sqrt(1.0-calp*calp)
			bet = 2.0*pi*grnd()
!				// new momentum components in new valley iv[j]
			kx = (pvx(ivv) + sqrt(mx(ivv)/amd)*pp*salp*cos(bet))/hbar
			ky = (pvy(ivv) + sqrt(my(ivv)/amd)*pp*salp*sin(bet))/hbar
			kz = (pvz(ivv) + sqrt(mz(ivv)/amd)*pp*calp)/hbar
			do while ((MAX(abs(kx),abs(ky),abs(kz)) > PQMAX/hbar) .AND. (BZCLIP==1)) 
!				// choose the set of random angles
				calp = 1.0-2.0*grnd()
				salp = sqrt(1.0-calp*calp)
				bet = 2.0*pi*grnd()
!				// new momentum components in new valley iv[j]
				kx = (pvx(ivv) + sqrt(mx(ivv)/amd)*pp*salp*cos(bet))/hbar
				ky = (pvy(ivv) + sqrt(my(ivv)/amd)*pp*salp*sin(bet))/hbar
				kz = (pvz(ivv) + sqrt(mz(ivv)/amd)*pp*calp)/hbar
			end do
!	  	  
!				#ifdef POP
!				// difference between final & initial momentum vectors
			dpx = hbar*kx-pxo
			dpy = hbar*ky-pyo
			dpz = hbar*kz-pzo
!				// get qphon = norm of the resultant q-vector reduced to 1st BZ
			qphon = abs(2.0*PQMAX-sqrt(dpx*dpx + dpy*dpy + dpz*dpz))/hbar
			if (qphon < QMIN) write(*,*) '** WARNING 6 inter_scat_g **: qphon < QMIN', qphon/QMAX
!				// get resulting phonon energy and momentum
			ephon = hbar*get_wq(qphon,LT) 
!				// estimate new electron energy
			Ep = eee - aed*ephon
!				// upper limit for cnt = 200, see right below
			cnt=cnt+1
			if ((Ep > EMIN) .AND. ((abs(ephon-ephon0)/ephon0) < 0.05)) REDO = 0
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
!				}⑧ while (REDO);
!				#endif
!  			// set the electron energy, consistent w/the momentum exchange above
	eee = eee - aed*ephon 
!
!			if (aed > 0.) assert(Elec[EE][j] > EMIN);
!			// return the phonon energy and momentum for STATS
!
!		**  phonon counting(+emit & -absorb)
	ie=int(ephon/DEph)+1
	if (ie > NBph) then
		write(*,*) '** WARNING 8 ** too large ephon',ephon
	else if (LT==LA .OR. LT==LO) then
		NELph0(ie)=NELph0(ie)+aed*1.0d0
	else if (LT==TA .OR. LT==TO) then
		NETph0(ie)=NETph0(ie)+aed*1.0d0
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
	Double Precision, parameter:: EMIN = 1.e-10   ! @ Pop-definition
	Double Precision,parameter::PQMAX = 7.6149280673e-08 ! hbar*QMAX(Wigner-Seitz cell radius)
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
!    
	EX=0.0
	cnt=0
	REDO=1
  
!		eee = E0  = Elec[EE][j];
!		ivv = iv0 = iv[j];

!			#ifdef POP   // store original momentum components and valley
	if ((aed > 0.) .AND. ((eee - hbar*get_wq(QFp*QMAX,LT)) <= EMIN)) then
!			write (*,*) 'eee - hbar*get_wq(QFp*QMAX,LT) <= EMIN'
		return   !  do nothing OK? 2018.10.12
	end if
	pxo = hbar*kx
	pyo = hbar*ky
	pzo = hbar*kz
!			#else
!			if ((lt == LA) || (lt == LO)) {
!			ephon = ELAf;
!			 } else if (lt == TA) {
!			ephon = ETAf;
!			} else {
!			ephon = ETOf;
!			}
!			qphon = QFp*QMAX;
!			if (aed > 0.) assert((eee - ephon) > EMIN);
!			#endif
!
!			// change valley upon scattering if we're in strained or bulk Si
!			if (SIGEX > 0.0) {
!			#ifdef POP
!			ephon = hbar*get_wq(QFp*QMAX,lt);  // MAKE SURE WHEN THIS lt IS LA OR LO
!			#endif
!			r = ranmc();
!			// DOSe(E-aed*ephon+DESIGEX);
!			if ((iv0+1)/2 == 2)    {       // iv = {3,4} (D2 lower valley)
!			if (r*DOSe(eee-aed*ephon+DESIGEX) > DOSe(eee-aed*ephon-DESIGEX)) {
!			return 0;                  // self-scattering and return
!			} else {
!			change_noneq_valley(j);        // --> {1,2,5,6} continue
!			EX = -DESIGEX;
!			}
!			} else {                 // iv = {1,2,5,6} (D4 upper valley)
!			if (r < 0.5) {
!			iv[j] = 3 + ifloor(2.*ranmc());        // --> {3,4}
!			EX = DESIGEX;
!			} else if (r < (0.5 + DOSe(eee-aed*ephon)/DOSe(eee-aed*ephon+DESIGEX))) {
!			if ((iv0+1)/2 == 1)            // iv = {1,2} --> {5,6}
!			iv[j] = 5 + ifloor(2.*r);
!			else                           // iv = {5,6} --> {1,2}
!			iv[j] = 1 + ifloor(2.*r);
!			} else {
!			return 0;
!			}
!			}
!			// DON'T FORGET TO ADD/SUBTRACT EX = +/-DESIGEX ENERGY BELOW!!!
!			} else {
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
!			Pxyz = PX-1+(int)(iv[j]+1)/2;
!			#ifdef POP
	if (qphon < QMIN) write(*,*) '** WARNING 9 inter_scat_f qphon < QMIN **', qphon/QMAX
!			// original estimate on phonon energy
	ephon = hbar*get_wq(qphon,LT);
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
			pp = sqrt(2.*mdos*Ep*(1.+af*Ep))
!				// choose the set of random angles
			calp = 1.0-2.0*grnd()
			salp = sqrt(1.0-calp*calp)
			bet = 2.0*pi*grnd()
!				// new momentum components in new valley iv[j]
			kx = (pvx(ivv) + sqrt(mx(ivv)/amd)*pp*salp*cos(bet))/hbar
			ky = (pvy(ivv) + sqrt(my(ivv)/amd)*pp*salp*sin(bet))/hbar
			kz = (pvz(ivv) + sqrt(mz(ivv)/amd)*pp*calp)/hbar
			do while ((MAX(abs(kx),abs(ky),abs(kz)) > PQMAX/hbar) .AND. (BZCLIP==1))
!				// choose the set of random angles
				calp = 1.0-2.0*grnd()
				salp = sqrt(1.0-calp*calp)
				bet = 2.0*pi*grnd()
!				// new momentum components in new valley iv[j]
				kx = (pvx(ivv) + sqrt(mx(ivv)/amd)*pp*salp*cos(bet))/hbar
				ky = (pvy(ivv) + sqrt(my(ivv)/amd)*pp*salp*sin(bet))/hbar
				kz = (pvz(ivv) + sqrt(mz(ivv)/amd)*pp*calp)/hbar
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
			ephon = hbar*get_wq(qphon,LT);
!				// re-estimate new electron energy
			Ep = eee - aed*ephon0 + EX
!				// upper limit for cnt = 200, see right below
			cnt=cnt+1
			if ((Ep > EMIN) .AND. ((abs(ephon-ephon0)/ephon0) < 0.05)) REDO = 0
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
!				} while (REDO);
!				#endif
!				// set the electron energy, consistent w/the momentum exchange above
	eee = eee - aed*ephon0 + EX 
!				if (aed > 0.) assert(Elec[EE][j] > EMIN);
!				// return the phonon energy and momentum for STATS
!
!		**  phonon counting(+emit & -absorb)
	ie=int(ephon/DEph)+1
	if (ie > NBph) then
		write(*,*) '** WARNING 12 ** too large ephon',ephon
	else if (LT==LA .OR. LT==LO) then
		NELph0(ie)=NELph0(ie)+aed*1.0d0
	else if (LT==TA .OR. LT==TO) then
		NETph0(ie)=NETph0(ie)+aed*1.0d0
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
!!  axis = (iv[j] + 1)/2;
!!  iv[j] = 4*axis - iv[j] - 1;
end subroutine change_equiv_valley
!========================================================================================
subroutine change_noneq_valley(ivv)
	implicit none
	integer ivv
	Double Precision r
!
	r = grnd()
	if      ((ivv == 1) .OR. (ivv == 2))  then
		ivv=3+int(4.0*r)  ! ivv = 1 or 2 -> {3,4,5,6}
	else if ((ivv == 3) .OR. (ivv == 4))  then    ! ivv = 3 or 4 -> {1,2}{5,6}
		ivv = 1 + int(4.0*r)
		if (ivv > 2) ivv= ivv+2
	else 
		ivv = 1 + int(4.0*r)
	end if                                        ! ivv = 5 or 6 -> {1,2,3,4}
end subroutine change_noneq_valley
!========================================================================================
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
	Double Precision  kk,eee                       ! debug
	integer n, ivv
!
!---( 平均速度と距離の計算 )---
!
	ve = 0.0
	eee = 0.0d0   ! debug
	s = 0.0
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
		ve = ve + kx*hm((ivv-1)/2+1)/sqrt(1.0+af4*gk) 
!			eee = eee+kk                              ! debug
		eee = eee + Elec(n,EE)
		s  = s + Elec(n,XXX)
	end do
	ve = ve/float(inum)
!		eee = eee/float(inum)/smh/smh            ! debug kk->E[eV]
	eee = eee/float(inum)
	s  = s /float(inum)
!
!---( 時間 vs 平均速度 vs 距離の出力 )---
	write(*,*) t,'  ',ve,'  ',s,'  ',eee  ! debug
	write(8,*) t,'  ',ve,'  ',s,'  ',eee  ! debug
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
!			eemax = 0.01   ! Debug fine division  10meV/300div -> 0.03meV/div
	do  n=1,inum
		digitized_energy=int(real(nhist-1)*(Elec(n,EE)/eemax))+1
		if (digitized_energy <= nhist) then
			ehist(digitized_energy)=ehist(digitized_energy)+1
		end if
	end do   
!
!---( エネルギーヒストグラムの出力 )---
!
	write(*,*) 'emax=',eemax, '  particle#', iemax ! debug
	write(8,*) 'emax=',eemax, '  particle#', iemax ! debug
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
	integer, parameter::NBph=280
	Double Precision, parameter::DEph=2.5e-4  ! 70meV/280
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
	allocate(Elec(inum_dummy,17))
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
		NELph0(i)=0.0d0
		NETph0(i)=0.0d0
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
		NELph0(i)=NELph0(i)*(1.0E17/inum_dummy)*(1/SimTime)
		NETph0(i)=NETph0(i)*(1.0E17/inum_dummy)*(1/SimTime)
		write (*,*) i, dble(i)*DEph,NELph0(i),NETph0(i)
		write (8,*) i, dble(i)*DEph,NELph0(i),NETph0(i)
	end do
!      
!---( ファイルのクローズ )---
!
	close(5)
	close(8)
	deallocate(Elec,iv)
	call cpu_time( t2 )
	write(*,*) 'cpu time:', t2-t1, 'seconds.'
	write(8,*) 'cpu time:', t2-t1, 'seconds.'
	stop
end program main
