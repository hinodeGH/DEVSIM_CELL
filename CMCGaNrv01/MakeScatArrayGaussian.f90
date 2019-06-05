!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     Program for Making Sacttering Rate Array
!     for Cellular Monte Carlo(CMC)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! rv00 2019.05.17 プログラミング開始
!
!
!****************************************************************************************
module CMC_sub_program_variables    !----- 共通変数 -----
!
	implicit none
	Double Precision,parameter::aGaN = 5.d-8			!	lattice constant, cm
	Double Precision,parameter::QMAX = 1.256d+08		!  2*pi/aGaN(lattice constant, cm)
	Double Precision,parameter::m0 = 5.68562975d-16    ! eV*s^2/cm^2  !electron mass
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
!
	Double Precision,parameter::LTAO(4,3) = &
						&  reshape( (/	-8.65d-4,	-1.45d-3,	-1.88d-4,	1.98d-4, &
						&				4.91d5,		4.88d5,		-1.45d4,	2.52d4,   &
						&				0.0d0,		0.0d0,		1.38d14,	1.07d14 /),(/4,3/) )
!
	Double Precision,parameter::echarge = 1.d0		!  // electron charge 
!				#define echarge   -1.0             // electron charge 
!				#define ecoulom   -1.60217653d-19  // electron charge in Coulombs
!
!									!----- 共通変数 -----
	Double Precision,SAVE:: Temp
	integer,SAVE:: nkx,nky,nkz		!	# of grid points along each axis; even number?
	Double Precision,SAVE:: kx_max,ky_max,kz_max		!	max k within grid area
	Double Precision,SAVE:: delta_kx,delta_ky,delta_kz	!	k grid pitch
	Double Precision,SAVE:: mdos
	Double Precision,SAVE:: gm
	Double Precision,SAVE:: alpha
	Double Precision,SAVE:: qhbar
	Double Precision,SAVE:: kbTq
!
end module CMC_sub_program_variables
!****************************************************************************************
module CMC_sub_programs
	use CMC_sub_program_variables
	implicit none
	contains
!--------------- 以下に次のsubroutine/function ------------
!	get_wq 		!	phonon q -> w
!	get_nq 		! 	phonon occupation
!	data_input 	! 	総粒子数, 総時間ステップ, 時間増分, 格子温度, 電界
!	param 		! 	物理定数, 材料定数および諸パラメータ;散乱レートの計算→表化
!
!========================================================================================
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
!****************************************************************************************   
!========================================================================================
!
subroutine data_input	!	input data for making CMC transion table
	implicit none
!
	Double Precision FractionOfBZ
!
	read(5,*) Temp				!	格子温度 [K] (typ.300)
	read(5,*) FractionOfBZ		!	Grid area fraction/ BZ (0.1~0.5)
	read(5,*) nkx,nky,nkz		!	how many grid points along each axis; even number?
!
	write(*,*) 'FractionOfBZ: ',FractionOfBZ
	write(*,*) '# of grid points: ',nkx,nky,nkz
	write(8,*) 'FractionOfBZ: ',FractionOfBZ
	write(8,*) '# of grid points: ',nkx,nky,nkz
!
	kx_max = FractionOfBZ*QMAX
	ky_max = FractionOfBZ*QMAX
	kz_max = FractionOfBZ*QMAX
	write(*,*) 'max k within grid area: ',kx_max,ky_max,kz_max
	write(8,*) 'max k within grid area: ',kx_max,ky_max,kz_max
!
end subroutine data_input
!========================================================================================
!===( 物理定数, 材料定数および諸パラメータ )===
!
subroutine param	!	set materials constants and calculate CMC transition table
	implicit none      
	Double Precision Egap
	Double Precision,parameter::eps0 = 5.52634972d+05	!	// e/V/cm permittivity of free space
	Double Precision,parameter::epsGaN = 5.2500322d+06	!	// e/V/cm permittivity of GaN 9.5*eps0
	Double Precision,parameter::epsf =   2.9565971d+06	!	// e/V/cm permittivity of GaN 5.35*eps0
	Double Precision ep
!
	integer ikx, iky, ikz, i_cell, i_cell_max	
	Double Precision, allocatable:: eee(:),kx(:),ky(:),kz(:),kxyz(:)
	Double Precision gk
	Double Precision LargeK,sgmGaussian,kKsgm,ModeSeparator		!	delta function -> Gaussian
!
	integer i_cell_initial, i_cell_final
	integer, allocatable:: mode(:,:)
	Double Precision, allocatable:: wTotal(:,:)
	Double Precision,parameter::QcntrBZ = 0.0d0	!	center of B.Z.
	Double Precision wLAems,wLAabs,wLOems,wLOabs
	Double Precision q_phon_sq,q_phon,e_phon
!!	Double Precision delta_eee
	Double Precision EkConservationLAems,EkConservationLAabs,EkConservationLOems,EkConservationLOabs
	Double Precision Max_sum_rate_to_final,Rate_sum_temp
	Double Precision, allocatable:: Rate_sum_to_final(:)
!!	----	2019.05.24 debug begin
					integer ikx_initial,iky_initial,ikz_initial, ikx_final,iky_final,ikz_final
!!					integer ikx_tmp,iky_tmp,ikz_tmp,OKNG
!!	----	2019.05.24 debug end
!
	Egap = 3.5							!	eV 禁制帯幅
	ep = 1.0/(1.0/epsf-1.0/epsGaN)		!	誘電率
!
	kbTq = kb*Temp/echarge
	qhbar   = echarge/hbar   			!	h->hbar
	mdos = 0.20*m0						!	effective mass
	alpha = (1.d0 - mdos/m0)**2/Egap	!	バンドの非放物線性
!	alpha = 0.	!	放物線!!	----	2019.05.24 debug
!
! ===== subroutine Make_Orthogonal_Grid()
! ===== implicit none
!

	i_cell_max = nkx*nky*nkz
	write(*,*) 'i_cell_max: ',i_cell_max
	allocate( eee(i_cell_max) )
	allocate( kx(i_cell_max) )
	allocate( ky(i_cell_max) )
	allocate( kz(i_cell_max) )
	allocate( kxyz(i_cell_max) )
    delta_kx = 2.*kx_max/nkx
    delta_ky = 2.*ky_max/nky
    delta_kz = 2.*kz_max/nkz
    sgmGaussian = 1.0*(delta_kx*delta_ky*delta_kz)**(1.0/3.0)	!	delta function -> Gaussian
    ModeSeparator = 300.0
    write(*,*) 'delta_kx,delta_ky,delta_kz: ',delta_kx,delta_ky,delta_kz
	write(8,*) 'delta_kx,delta_ky,delta_kz: ',delta_kx,delta_ky,delta_kz
	write(*,*) 'm* ', mdos
	write(8,*) 'm* ', mdos
	write(*,*) 'alpha ', alpha
	write(8,*) 'alpha ', alpha
!
	do ikz=1, nkz	!	ikz loop
		do iky=1, nky	!	iky loop
			do ikx=1, nkx	!	ikx loop
				i_cell = ikx + (iky-1)*nkx + (ikz-1)*nkx*nky	!	linearly aligned cell number
				kx(i_cell) = -kx_max + delta_kx*(ikx-1) + 0.5*delta_kx	!	kx at cell center
				ky(i_cell) = -ky_max + delta_ky*(iky-1) + 0.5*delta_ky	!	ky at cell center
				kz(i_cell) = -kz_max + delta_kz*(ikz-1) + 0.5*delta_kz	!	kz at cell center
				gk = (hbar*hbar/(2.*mdos))*(kx(i_cell)**2 + ky(i_cell)**2 + kz(i_cell)**2)
				eee(i_cell) = 2.*gk/(sqrt(1.+4.*alpha*gk)+1)
				kxyz(i_cell) = Sqrt( kx(i_cell)**2+ky(i_cell)**2+ky(i_cell)**2 )
				write(*,*) i_cell,'ikx iky ikz eee',kx(i_cell),ky(i_cell),kz(i_cell),eee(i_cell)
!!	----	2019.05.24 debug begin
!!				write(8,*) i_cell,'ikx iky ikz eee',kx(i_cell),ky(i_cell),kz(i_cell),eee(i_cell)
!!				ikz_tmp = 1+(i_cell-1)/(nkx*nky)
!!				iky_tmp = 1+( i_cell-1-(ikz_tmp-1)*(nkx*nky) )/nkx
!!				ikx_tmp = i_cell - (ikz_tmp-1)*(nkx*nky) - (iky_tmp-1)*nkx
!!				if (ikx==ikx_tmp .AND. iky==iky_tmp .AND. ikz==ikz_tmp ) then
!!					OKNG=1
!!				else
!!					OKNG=0
!!				end if
!!				write(*,*) i_cell,' ',ikx,' ',iky,' ',ikz,' ',ikx_tmp,' ',iky_tmp,' ',ikz_tmp,OKNG	! debug
!!				write(8,*) i_cell,' ',ikx,' ',iky,' ',ikz,' ',ikx_tmp,' ',iky_tmp,' ',ikz_tmp,OKNG	! debug
!!	----	2019.05.24 debug end
			end do			!	ikx do loop
		end do			!	iky do loop
    end do			!	ikz do loop
! ===== end	subroutine Make_Orthogonal_Grid
!----------------------------------------------------------------------------------------
!!	return	!!	----	2019.05.24 debug
! ===== subroutine Scattering_Rate_Array()
! ===== implicit none
!
	allocate( mode(i_cell_max,i_cell_max) )
	allocate( wTotal(i_cell_max,i_cell_max) )
	do i_cell_initial=1, i_cell_max
		write(*,*) 'i_cell_initial: ',i_cell_initial
		do i_cell_final=1, i_cell_max
			mode(i_cell_initial, i_cell_final) = 0
			if	(i_cell_final == i_cell_initial) then
				wLAems = 0.
				wLAabs = 0.
				wLOems = 0.
				wLOabs = 0.
			else
! ----- initial_center to final_center 
				q_phon_sq = (kx(i_cell_final) - kx(i_cell_initial))**2 + (ky(i_cell_final) - ky(i_cell_initial))**2  &
				&	 + (kz(i_cell_final) - kz(i_cell_initial))**2	!	square of phonon vector
				q_phon = sqrt(q_phon_sq)	!	phonon vector
! ------ LA emission or absorption
				e_phon = hbar*get_wq(q_phon, LA)
! ------ LA emission
!
! ------ Large K
				if ((1.0 - e_phon/eee(i_cell_initial)) > 0. ) then
					LargeK = kxyz(i_cell_initial)*sqrt(1.0 - e_phon/eee(i_cell_initial))
!
					wLAems = (DLA**2*q_phon_sq)/(8.*pi*pi*rho)/get_wq(q_phon, LA)*(get_nq(q_phon, LA)+1.)
					kKsgm = (kxyz(i_cell_final)-LargeK)/sgmGaussian
					if (kKsgm < ModeSeparator) then
						EkConservationLAems =  (kxyz(i_cell_initial)**2/(2.0*eee(i_cell_initial)*LargeK))	&
								& /(sgmGaussian*sqrt(pi))*exp(-kKsgm**2)*delta_kx*delta_ky*delta_kz
!!						write(*,*) wLAems,EkConservationLAems,wLAems*EkConservationLAems
!!						write(8,*) wLAems,EkConservationLAems,wLAems*EkConservationLAems
						mode(i_cell_initial, i_cell_final) = 1
					else
						EkConservationLAems = 0.0
					end if
				else
					EkConservationLAems =  0.0		!	no emission
				end if
!
! ------ LA absorption
!
! ------ Large K
				if ((1.0 + e_phon/eee(i_cell_initial)) > 0. ) then
					LargeK = kxyz(i_cell_initial)*sqrt(1.0 + e_phon/eee(i_cell_initial))
!
					wLAabs = (DLA**2*q_phon_sq)/(8.*pi*pi*rho)/get_wq(q_phon, LA)*get_nq(q_phon, LA)
					kKsgm = (kxyz(i_cell_final)-LargeK)/sgmGaussian
					if (kKsgm < ModeSeparator) then
						EkConservationLAabs =  (kxyz(i_cell_initial)**2/(2.0*eee(i_cell_initial)*LargeK))	&
								& /(sgmGaussian*sqrt(pi))*exp(-kKsgm**2)*delta_kx*delta_ky*delta_kz
						mode(i_cell_initial, i_cell_final) = 2 + mode(i_cell_initial, i_cell_final)
					else
						EkConservationLAabs = 0.0
					end if
				else
					EkConservationLAabs = 0.0		!	no absorption?
				end if
!
! ------ LO emission or absorption				
				e_phon = hbar*get_wq(q_phon, LO)
! ------ LO emission
!
! ------ Large K
				if ((1.0 - e_phon/eee(i_cell_initial)) > 0. ) then
					LargeK = kxyz(i_cell_initial)*sqrt(1.0 - e_phon/eee(i_cell_initial))
!
					wLOems = (1.**2*get_wq(q_phon, LO))/(8.*pi*pi*ep*q_phon_sq)*(get_nq(q_phon,LO)+1.)
					kKsgm = (kxyz(i_cell_final)-LargeK)/sgmGaussian
					if (kKsgm < ModeSeparator) then
						EkConservationLOems =  (kxyz(i_cell_initial)**2/(2.0*eee(i_cell_initial)*LargeK))	&
								& /(sgmGaussian*sqrt(pi))*exp(-kKsgm**2)*delta_kx*delta_ky*delta_kz
						mode(i_cell_initial, i_cell_final) = 4 + mode(i_cell_initial, i_cell_final)
					else
						EkConservationLOems = 0.0
					end if
				else
					EkConservationLOems = 0.0	!	no emission
				end if
! ------ LO absorption
!
!
! ------ Large K
				if ((1.0 + e_phon/eee(i_cell_initial)) > 0. ) then
					LargeK = kxyz(i_cell_initial)*sqrt(1.0 + e_phon/eee(i_cell_initial))
!
					wLOabs = (1.**2*get_wq(q_phon, LO))/(8.*pi*pi*ep*q_phon_sq)*get_nq(q_phon, LO)
					kKsgm = (kxyz(i_cell_final)-LargeK)/sgmGaussian
					if (kKsgm < ModeSeparator) then
						EkConservationLOabs =   (kxyz(i_cell_initial)**2/(2.0*eee(i_cell_initial)*LargeK))	&
								& /(sgmGaussian*sqrt(pi))*exp(-kKsgm**2)*delta_kx*delta_ky*delta_kz
						mode(i_cell_initial, i_cell_final) = 8 + mode(i_cell_initial, i_cell_final)
					else
						EkConservationLOabs = 0.
					end if
				else
					EkConservationLOabs = 0.	!	no absorption
				end if
!!	----	2019.05.24 begin
!!				wTotal(i_cell_initial,i_cell_final) = wLAems*EkConservationLAems + wLAabs*EkConservationLAabs  &
!!												&	+ wLOems*EkConservationLOems + wLOabs*EkConservationLOabs
				wTotal(i_cell_initial,i_cell_final) =  wLOabs*EkConservationLOabs
!!	----	2019.05.24 end
			end if
		end do	!	i_cell_final do loop
	end do	!	i_cell_initial do loop
!
! ===== end subroutine Scattering_Rate_Array
!----------------------------------------------------------------------------------------
!!	----	2019.05.24 begin
!!	return
!!	----	2019.05.24 end
! ===== subroutine Scattering_Rate_to_Probability_Array()
! ===== implicit none
!	--		rate sum of all final state cells, for each initial cell ==> Max_sum_rate_to_final
	allocate ( Rate_sum_to_final(i_cell_max) )
	Max_sum_rate_to_final = 0.
	do i_cell_initial=1, i_cell_max
		write(*,*) 'i_cell_initial: ',i_cell_initial
		Rate_sum_temp = 0.	!	for each i_cell_initial
		do i_cell_final=1, i_cell_max
			Rate_sum_temp = Rate_sum_temp + wTotal(i_cell_initial,i_cell_final)
		end do	!	i_cell_final do loop
		if (Max_sum_rate_to_final < Rate_sum_temp) then
			Max_sum_rate_to_final = Rate_sum_temp
		end if
		Rate_sum_to_final(i_cell_initial) = Rate_sum_temp	!	for each i_cell_initial
		write(*,*) ' i_cell_initial, eee, Rate_sum_to_final: ',		&
					&	i_cell_initial,' ',eee(i_cell_initial),' ',Rate_sum_to_final(i_cell_initial)
		write(8,*) ' i_cell_initial, eee, Rate_sum_to_final: ',		&
					&	i_cell_initial,' ',eee(i_cell_initial),' ',Rate_sum_to_final(i_cell_initial)
	end do	!	i_cell_initial do loop
!!	----	2019.05.24 begin
	return
!!	----	2019.05.24 end
!	--		diagonal components
	do i_cell = 1, i_cell_max
		wTotal(i_cell,i_cell) = Max_sum_rate_to_final - Rate_sum_to_final(i_cell)
	end do	!	i_cell do loop
!	--		rate accumulation for Monte Carlo using uniform random number
	do i_cell_initial=1, i_cell_max
		do i_cell_final=2, i_cell_max
			wTotal(i_cell_initial,i_cell_final) = wTotal(i_cell_initial,i_cell_final - 1) + wTotal(i_cell_initial,i_cell_final)
		end do	!	i_cell_final do loop
	end do	!	i_cell_initial do loop
!	--		Normalize rate to probability
	do i_cell_initial=1, i_cell_max
		do i_cell_final=1, i_cell_max
			wTotal(i_cell_initial,i_cell_final) = wTotal(i_cell_initial,i_cell_final)/Max_sum_rate_to_final
		end do	!	i_cell_final do loop
	end do	!	i_cell_initial do loop
!
! ===== end subroutine Scattering_Rate_to_Probability_Array
!========================================================================================
	deallocate  (eee,kx,ky,kz)
	deallocate (mode,wTotal)
	deallocate (Rate_sum_to_final)
!
end subroutine param
!
end module CMC_sub_programs
!========================================================================================
!****************************************************************************************
program main
!
	use CMC_sub_programs
	implicit none
	Double Precision t1, t2
!
	call cpu_time( t1 )
!---( ファイルのオープン )---
	open(5,file='dataf')
	open(8,file='outf')
!
!---( 計算条件入力　散乱テーブル作成 )---
!
	call data_input
	call param
!
!---( ファイルのクローズ )---
!
	call cpu_time( t2 )
	write(*,*) 'cpu time:', t2-t1, 'seconds.'
	write(8,*) 'cpu time:', t2-t1, 'seconds.'
	close(5)
	close(8)
	stop
end program main
