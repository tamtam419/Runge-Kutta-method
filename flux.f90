program flux_calc
    implicit none
    integer, parameter :: array_number = 1300 !地球中心の距離rのメッシュ数
    integer, parameter :: distance_flux = 10 !地球の半径12800kmを何km毎にとるか
    integer, parameter :: array_crusttype = 64800

    character(100) :: filename !mantleの分布によるファイル名の違い
    character(15) :: sute, kind(array_crusttype) !kindにはcontinentalかoceanicかのどちらかを格納
    integer :: detectortype, distributiontype, secondcontinent
    integer :: i, j, k
    integer :: fi1 = 10, fi2 = 11, fi3 = 12, fi4 = 13, fi5 = 14
    integer :: fi6 = 15, fi7 = 16, fi8 = 17
    !integer :: fi9 = 18, fi10 = 19
    integer :: fo1 = 20
    integer :: ntheta = 180 !theta(緯度方向)の分割数
    integer :: nphi = 360 !phi(経度方向)の分割数
    integer :: htheta = 20, hphi = 20, hr = 8 !crustにおけるさらなる分割数
    integer :: n = 0 !do文中のtheta,phiの繰り返し数,180 x 360 = 64800まで
    integer :: z = 0
    integer :: beforestep = 0 
    integer :: afterstep = 0 !beforestep,afterstepは共にmantleからcrustへの計算の切り替えに使う
    integer :: c = 0, d = 0, e = 0
    integer :: crusttype = 0
    real(8) mh, mtheta, mphi !detectorの角度、標高
    real(8) r, theta, phi, dr, dtheta, dphi, dis !代表点の位置とメッシュの幅
    real(8) :: rep_r = 0.0d0 !
    real(8) ddtheta, ddphi, ddr !crustでの新たなメッシュ幅
    real(8) :: th_tem = 0.0d0, ph_tem = 0.0d0 !計算に必要な一時的な置き場
    real(8) :: Uflux = 0.0d0, Thflux = 0.0d0
    real(8) pi, rho
    real(8) R0(array_crusttype)
    real(8) :: zU = 0.0d0, zTh = 0.0d0
    real(8) :: a(504), b(504) !PREMデータを格納、aが地球中心からの距離、bが密度
	real(8) :: r1(array_crusttype), r2(array_crusttype), r3(array_crusttype)
	real(8) :: r4(array_crusttype), r5(array_crusttype), r6(array_crusttype)
    real(8) :: r7(array_crusttype), r8(array_crusttype), r9(array_crusttype)
	real(8) :: ro1(array_crusttype), ro2(array_crusttype), ro3(array_crusttype)
	real(8) :: ro4(array_crusttype), ro5(array_crusttype), ro6(array_crusttype)
	real(8) :: ro7(array_crusttype), ro8(array_crusttype)
	real(8) :: th(array_crusttype), ph(array_crusttype)	
    real(8) :: meflux(array_number)
    real(8) :: mer(array_number)
    real(8) :: meUflux = 0.0d0, meThflux = 0.0d0
    real(8) :: theta_tem = 0.0d0, mtheta_tem = 0.0d0
	real(8) :: allV = 0.0d0, localV = 0.0d0 !第三大陸、局所的第三大陸の体積
	real(8) :: r2h(array_crusttype)
	real(8) :: meccfluxsum = 0.0d0
	real(8) :: r_pall = 0.0d0
	real(8) :: ocU = 0.0d0, ocTh = 0.0d0 !Oceanic Crustでの係数
	real(8) :: ccU(3), ccTh(3) !Continental Crustでの係数
	real(8) :: maU = 0.0d0, maTh = 0.0d0
	real(8) :: seconU = 0.0d0, seconTh = 0.0d0
    real(8) :: seoceU = 0.0d0, seoceTh = 0.0d0
	real(8) :: maUmass = 0.0d0, maThmass = 0.0d0
	real(8) :: conUmass = 0.0d0, conThmass = 0.0d0
    real(8) :: oceUmass = 0.0d0, oceThmass = 0.0d0
    real(8) :: conseUmass = 0.0d0, conseThmass = 0.0d0
    real(8) :: oceseUmass = 0.0d0, oceseThmass = 0.0d0 
	real(8) :: meocflux = 0.0d0, memaflux = 0.0d0, meccflux(3)
	real(8) :: ocflux(array_number), maflux(array_number), ccflux(array_number)
	real(8) :: Umass = 0.0d0, Thmass = 0.0d0	
	real(8) :: A_U238 = 1.24d7 , A_Th232 = 4.06d6
	real(8) :: pos = 0.59d0 
	real(8) :: n_U238 = 6.0d0, n_Th232 = 4.0d0
    real(8) :: sute1 = 0.0d0, sute2 = 0.0d0

!------------------------------------------------------------------------------------------------
!係数の設定
	pi = 2.0d0 * acos(0.0d0)
	dtheta = pi / dble(ntheta)
	dphi = 2.0d0*pi / dble(nphi)
	ddtheta = dtheta / dble(htheta)
	ddphi = dphi / dble(hphi)
!------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
!初期化、ポイントの座標データ(r,theta,phi)系と配列について
    r = 0.0d0
	theta = 0.0d0
	phi = 0.0d0
	dis = 0.0d0
	mh = 0.0d0
	mtheta = 0.0d0
	mphi = 0.0d0
	meflux = 0.0d0
	mer = 0.0d0
	ccU = 0.0d0
	ccTh = 0.0d0
	meccflux = 0.0d0
	ocflux = 0.0d0
	maflux = 0.0d0
	ccflux = 0.0d0
!------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
!ファイルオープン
	open(fi1, file='crust1_rho.dat' , action = 'read')
	open(fi2, file='crust1_bnds.dat' , action = 'read')
	open(fi3, file='xyz_ro1.dat' , action = 'read')
	open(fi4, file='crust_type.dat' , action = 'read')
    open(fi5, file='PREM500.dat' , action = 'read')
	open(fi6, file='worldmap_con_after.dat' , status = 'replace')
	open(fi7, file='worldmap_oce_after.dat' , status = 'replace')
	open(fi8, file='nplace.dat' , status = 'replace')
!------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
!以下、ファイルから値を出力、各配列に格納、その後再設定

	!PREM
    do i = 1 , 504
        read(fi5,*) a(i) , b(i)
    enddo 

    !CRUST1の密度
    do i = 1 , array_crusttype
        read(fi1,*) ro1(i) , ro2(i) , ro3(i) , ro4(i)&
         , ro5(i) , ro6(i) , ro7(i) , ro8(i) 
    enddo
    
    !CRUST1の密度をg/cm^3 → kg/m^3 に再設定
    do i = 1 , array_crusttype
        ro1(i) = ro1(i) * 1000.0d0
        ro2(i) = ro2(i) * 1000.0d0
        ro3(i) = ro3(i) * 1000.0d0
        ro4(i) = ro4(i) * 1000.0d0
        ro5(i) = ro5(i) * 1000.0d0
        ro6(i) = ro6(i) * 1000.0d0
        ro7(i) = ro7(i) * 1000.0d0
        ro8(i) = ro8(i) * 1000.0d0
    enddo

    !CRUST1の角度
    do i = 1 , array_crusttype
        read(fi3,*) ph(i) , th(i)
    enddo

	!CRUST1の角度をrad、及び緯度・経度を北極からの角度に再設定
    do i = 1 , array_crusttype
        theta_tem = abs(th(i)) * pi / 180.0d0
        if( ph(i) < 0.0d0 ) then
        	ph(i) = ph(i) + 360.0d0
            !ph(i) = abs(ph(i) - 180.0d0)
        endif
        th(i) = abs( ( th(i) - 90.0d0 ) * pi / 180.0d0)
        ph(i) = ph(i) * pi / 180.0d0 
        R0(i) = sqrt( (6378137.0d0*cos( theta_tem ))**2 + (6356752.0d0*sin(theta_tem))**2 ) ! 地球の短径6356752m、長径
    enddo !thetaを北極0の0 ~ 180度。phiを東経を+として、0 ~ 360度に、例えば西経90度は180 + 90 = 270度に。

	!CRUST1の境界面の高さ
	do i = 1 , array_crusttype
		read(fi2,*) r1(i) , r2(i) , r3(i) , r4(i)&
		, r5(i) , r6(i) , r7(i) , r8(i) , r9(i)
	enddo
			
	!CRUST1の境界面の高さを地球の楕円性による、その角度の地球半径を考慮した。さらにkm → mに再設定。
	do i = 1 , array_crusttype
		r1(i) = r1(i) * 1000.0d0 + R0(i)
		r2(i) = r2(i) * 1000.0d0 + R0(i)
		r3(i) = r3(i) * 1000.0d0 + R0(i)
		r4(i) = r4(i) * 1000.0d0 + R0(i)
		r5(i) = r5(i) * 1000.0d0 + R0(i)
		r6(i) = r6(i) * 1000.0d0 + R0(i)
		r7(i) = r7(i) * 1000.0d0 + R0(i)
		r8(i) = r8(i) * 1000.0d0 + R0(i)
		r9(i) = r9(i) * 1000.0d0 + R0(i)
	enddo
		
    !各角度での海底の高さ+200を求めている。これをOceanicかContinentalかの区別に用いる場合もある。
	do i = 1 , array_crusttype
		r2h(i) = r2(i) + 200.0d0
		!r2h(i) = r2(i) + 3000.0d0
	enddo

    !CRUST1の各角度のOceanic、Continentalの区別kind(i)を用いる。
    do i = 1 , array_crusttype
		read(fi4,*) sute1 , sute2 , sute , kind(i)
	enddo
!------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
!第三大陸の体積を求める
    allV = 4.0d0 * pi * ( a(298)**3 - a(274)**3 ) / 3.0d0
    localV = pi * ( a(298)**3 - a(274)**3 ) * ( cos( 7.0d0 * pi / 18.0d0 ) - cos( 13.0d0 * pi / 18.0d0 ) ) / 9.0d0
    !a(298)はマントルの底、a(274)はマントル底から約300kmのところ
!------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
!ディテクター位置とマントル分布のタイプをインプット
!ディテクターの位置KamLANDなら1、Hawaiiなら2、南太平洋海底なら3、南太平洋海面なら4を入力。
    write(*,*) 'Select Detector type : KamLAND = 1, Hawaii = 2, bottom of South Pacific Ocean=3, surface of South Pacific Ocean=4'
    write(*,*) 'Input Detector type : '
    read(*,*) detectortype
!第三大陸の分布、先行研究(第三大陸を考えない場合)と同じなら1、均一第三大陸の場合なら2、局所的第三大陸の場合なら3を入力。
    write(*,*) 'Select Distribution type : Referrence = 1, Uniform Third-continent = 2, Local Third-continent = 3'
    write(*,*) 'Input Distribution type : '
    read(*,*) distributiontype
    write(*,*) 'Select Second Continent type : Nothing Second Continent = 1, Uniform Second-continent = 2, &
    Local Second-continent = 3'
    write(*,*) 'Input Distribution type : '
    read(*,*) secondcontinent
!------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
!インプットされたディテクターから、各値を設定、出力されるfilenameも決定
    if(detectortype == 1) then
        !KamLAND
        mtheta_tem = (36.0d0 + 25.0d0 / 60.0d0 + 36.0d0 / 60.0d0 / 60.0d0) * pi / 180.0d0
        mtheta = (90.0d0 - (36.0d0 + 25.0d0 / 60.0d0 + 36.0d0 / 60.0d0 / 60.0d0)) * pi / 180.0d0
        mphi = (137.0d0 + 18.0d0 / 60.0d0 + 43.0d0 / 60.0d0 / 60.0d0) * pi / 180.0d0
        mh = 358.0d0 + sqrt( (6378137.0d0*cos( mtheta_tem ))**2 + (6356752.0d0*sin(mtheta_tem))**2 )
        filename = 'graph(K).dat'
    else if(detectortype == 2) then
        !ハワイ
        mtheta_tem = (19.72d0) * pi / 180.0d0
        mtheta = ( 90.0d0 - 19.72d0 ) * pi / 180.0d0
        mphi = ( 360.0d0 - 156.32d0 ) * pi / 180.0d0
        mh = -4000.0d0 + sqrt( (6378137.0d0*cos( mtheta_tem ))**2 + (6356752.0d0*sin(mtheta_tem))**2 )
        filename = 'graph(H).dat'
    else if(detectortype == 3) then
        !南太平洋の底
        mtheta_tem = (10.5d0) * pi / 180.0d0
        mtheta = ( 90.0d0 + 10.5d0 ) * pi / 180.0d0
        mphi = ( 360.0d0- 170.5d0 ) * pi / 180.0d0
        mh = -4000.0d0 + sqrt( (6378137.0d0*cos( mtheta_tem ))**2 + (6356752.0d0*sin(mtheta_tem))**2 )
        filename = 'graph(b).dat'
    else if(detectortype == 4) then
        !南太平洋の海面
        mtheta_tem = (10.5d0) * pi / 180.0d0
        mtheta = ( 90.0d0 + 10.5d0 ) * pi / 180.0d0
        mphi = ( 360.0d0- 170.5d0 ) * pi / 180.0d0
        mh = sqrt( (6378137.0d0*cos( mtheta_tem ))**2 + (6356752.0d0*sin(mtheta_tem))**2 )
        filename = 'graph(s).dat'
    endif

    !出力されるファイルの名前をディテクター毎に区別
    open(fo1, file=filename, status='replace')
!------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
!以下、計算を回す
    do k = 1 , ntheta  !ntheta = 180
        do j = 1 , nphi  !nphi = 360
            do i = 0 , 515 !iの値は適当に515、無駄にループを回しているけど計算には影響しない。

                !------------------------------------------------------------------------------------------------
                !doのk,j(theta,phi)により代表点の角度を決定
                n = (k - 1) * 360 + j
                theta = th(n)
                phi = ph(n)
                th_tem = theta
                ph_tem = phi
                !------------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------------
                !crusttypeが出力されるかの確認。正しければrについてのdoを止めた状態で64800個のデータが出力
                if(kind(n) == 'Oceanic') then
                    crusttype = 1
                else if (kind(n) == 'Continental') then
                    crusttype = 2
                endif
                !write(*,*) n, crusttype
                !write(fi8,*) theta , phi , crusttype, kind(n)
                !------------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------------
                !r,dr,rhoの設定、mantleとcrustで場合分けを行い、その起点となるafterstepも求める
                if(i <= 503) then
                    r = a(i)
                    dr = a(i+1) - a(i)
                    rho = b(i)
                endif

				!PREMからcrustへの切り替え
                if( r9(n) >= r .and. r9(n) - r <= dr) then
                    dr = r9(n) - r
                    afterstep = 9
                endif
                !CRUSTの場合
                if(beforestep == 8 ) then
                    !Lower crust bottom to Middle
                    r = r9(n)
                    dr = r8(n) - r9(n)
                    rho = ro8(n)
                    afterstep = 8
                else if(beforestep == 7 ) then
                    !Middle to Upper
                    r = r8(n)
                    dr = r7(n) - r8(n)
                    rho = ro7(n)
                    afterstep = 7
                else if(beforestep == 6 ) then
                    !Upper crust bottom to Lower sediment bottom
                    r = r7(n)
                    dr = r6(n) - r7(n)
                    rho = ro6(n)
                    afterstep = 6
                else if(beforestep == 5 ) then
                    !Lower to Middle
                    r = r6(n)
                    dr = r5(n) - r6(n)
                    rho = ro5(n)
                    afterstep = 5
                else if(beforestep == 4 ) then
                    !Middle to Upper
                    r = r5(n)
                    dr = r4(n) - r5(n)
                    rho = ro4(n)
                    afterstep = 4
                else if(beforestep == 3 ) then
                    !Upper sediment bottom to Ice bottom
                    r = r4(n)
                    dr = r3(n) - r4(n)
                    rho = ro3(n)
                    afterstep = 3
                else if(beforestep == 2 ) then
                    !Ice bottom to Water bottom
                    r = r3(n)
                    dr = r2(n) - r3(n)
                    rho = ro2(n)
                    afterstep = 2
                else if(beforestep == 1 ) then
                    !Water bottom to Surface
                    r = r2(n)
                    dr = r1(n) - r2(n)
                    rho = ro1(n)
                    afterstep = 1
                else if(beforestep == 99 ) then
                    !Crust・Sediment終わり
                    dr = 0.0d0
                    rho = 0.0d0
                endif !crustの場合のr,dr,rho,afterstepについて
                ddr = dr / dble(hr)
                r_pall = r !crustで更に分割して計算する時に使用
                rep_r = sqrt( ( (r+dr)**2 + (r+dr)*r + r**2 ) / 3.0d0 )
                !------------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------------
                !ここまでの値からdisを決定
                dis = SQRT( ( mh*sin(mtheta)*cos(mphi) - (r+dr/2.0d0)*sin(theta)*cos(phi) )**2.0d0 +& 
                ( mh*sin(mtheta)*sin(mphi) - (r+dr/2.0d0)*sin(theta)*sin(phi) )**2.0d0 +&
                ( mh*cos(mtheta) - (r+dr/2.0d0)*cos(theta) )**2.0d0  )
                							
                !dis = SQRT( ( mh*sin(mtheta)*cos(mphi) - r*sin(theta)*cos(phi) )**2.0d0 +& 
				!( mh*sin(mtheta)*sin(mphi) - r*sin(theta)*sin(phi) )**2.0d0 +&
				!( mh*cos(mtheta) - r*cos(theta) )**2.0d0  ) 

				!dis = SQRT( ( mh*sin(mtheta)*cos(mphi) - rep_r*sin(theta)*cos(phi) )**2.0d0 +& 
				!( mh*sin(mtheta)*sin(mphi) - rep_r*sin(theta)*sin(phi) )**2.0d0 +&
				!( mh*cos(mtheta) - rep_r*cos(theta) )**2.0d0  ) 
                !------------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------------
                !MantleでのU,Th係数
                if (i > 273) then !coreからMantle
                    if(distributiontype ==1 ) then
                        !先行研究と同じタイプ
                        if (beforestep == 0 ) then
                            zU = 0.012d-6 * rho
                            zTh = 0.048d-6 * rho
                            maU = 0.012d-6 * rho
                            maTh = 0.048d-6 * rho
                        endif
                    else if(distributiontype == 2 ) then
                        if( i < 299) then
                            !均一第三大陸モデル
                            !zU = 4.813d16 / allV
                            !zTh = 19.252d16 / allV
                            !maU = 4.813d16 / allV
                            !maTh = 19.252d16 / allV
                            zU = 4.797d16 / allV
                            zTh = 19.188d16 / allV
                            maU = 4.797d16 / allV
                            maTh = 19.188d16 / allV
                        endif
                    else if(distributiontype == 3 ) then
                        if( i < 299) then
                            !局所的第三大陸モデル
                            if(phi <= ( 360.0d0 - 140.0d0 ) * pi / 180.0d0 .and. &
                            phi >= 160.0d0 * pi / 180.0d0 ) then
                                if(theta >= 70.5d0  * pi / 180.0d0 .and. &
                                theta <= 130.5d0  * pi / 180.0d0) then
                                    !zU = 4.813d16 / ( 2.0d0 * localV )
                                    !zTh = 19.252d16 / ( 2.0d0 * localV )
                                    !maU = 4.813d16 / ( 2.0d0 * localV )
                                    !maTh = 19.252d16 / ( 2.0d0 * localV )
                                    zU = 4.797d16 / ( 2.0d0 * localV )
                                    zTh = 19.188d16 / ( 2.0d0 * localV )
                                    maU = 4.797d16 / ( 2.0d0 * localV )
                                    maTh = 19.188d16 / ( 2.0d0 * localV )
                                endif
                            else
                                !zU = 4.813d16 / ( 2.0d0 * ( allV - localV ) )
                                !zTh = 19.252d16 / ( 2.0d0 * ( allV - localV ) )
                                !maU = 4.813d16 / ( 2.0d0 * ( allV - localV ) )
                                !maTh = 19.252d16 / ( 2.0d0 * ( allV - localV ) )
                                zU = 4.797d16 / ( 2.0d0 * ( allV - localV ) )
                                zTh = 19.188d16 / ( 2.0d0 * ( allV - localV ) )
                                maU = 4.797d16 / ( 2.0d0 * ( allV - localV ) )
                                maTh = 19.188d16 / ( 2.0d0 * ( allV - localV ) )
                            endif
                        endif
                    endif
                endif !coreからmantleでのU,Th係数について
                !------------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------------
                !第二大陸の計算。U,Thを密度に追加。
                if(secondcontinent == 0) then
                    !第二大陸を考えない場合
                else if(secondcontinent == 1) then
                    if( i >= 449 .and. i <= 468 ) then
                        !均一第二大陸分布
                        zU = zU + 2.728d-3
                        zTh = zTh + 1.055d-2
                        maU = maU + 2.728d-3
                        maTh = maTh + 1.055d-2
                    endif
                else if(secondcontinent ==2) then 
                    if( i >= 449 .and. i <= 468 ) then
                        if(crusttype == 2) then
                            zU = zU + 7.123d-3
                            zTh = zTh + 2.754d-2
                            maU = maU + 7.123d-3
                            maTh = maTh + 2.754d-2
                        endif
                    endif
                endif
                !------------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------------
                !CrustでのU,Thの係数
				if( beforestep <= 8 .and. beforestep >= 1) then !mantleからcrust
					!if( r2h(n) < r1(n) ) then
					!if( kind(n) == 'Oceanic' ) then
                    if( crusttype == 1 ) then 
                        !Oceanic
                        write(fi7,*) phi, theta
                        !if( r <= r6(n) ) then
                        if( beforestep >= 6) then
                            !crust
                            zU = 0.10d-6 * rho
                            zTh = 0.22d-6 * rho
                            ocU = 0.10d-6 * rho
                            ocTh = 0.22d-6 * rho
                            write(*,*) beforestep 
                        !else if( r <= r3(n) ) then
                        else if( beforestep >= 3 ) then
                            !sediment
                            zU = 1.7d-6 * rho
                            zTh = 6.9d-6 * rho
                            seoceU = 1.7d-6 * rho
                            seoceTh = 6.9d-6 * rho
                            write(*,*) beforestep 
                        !else if( r <= r1(n) ) then
                        !	zU = 3200.0d-9
                        !	zTh = 0.02d-9
                        endif !crust to sediment
                    !else if( r2h(n) >= r1(n) ) then
                    !else if( kind(n) == 'Continental' ) then
                    else if(crusttype == 2) then
                        !Continental
                        write(fi6,*) phi, theta
                        !if( r <= r8(n) ) then
                        if( beforestep == 8) then
                            !lower crust
                            zU = 0.2d-6 * rho
                            zTh = 1.2d-6 * rho
                            ccU(1) = 0.2d-6 * rho
                            ccTh(1) = 1.2d-6 * rho
                            write(*,*) beforestep 
                        !else if( r <= r7(n) ) then
                        else if( beforestep == 7 ) then
                            !middle crust
                            zU = 1.6d-6 * rho
                            zTh = 6.1d-6 * rho
                            ccU(2) = 1.6d-6 * rho
                            ccTh(2) = 6.1d-6 * rho
                            write(*,*) beforestep 
                        !else if( r <= r6(n) ) then
                        else if( beforestep == 6 ) then
                            !upper crust
                            zU = 2.8d-6 * rho
                            zTh = 10.7d-6 * rho
                            ccU(3) = 2.8d-6 * rho
                            ccTh(3) = 10.7d-6 * rho
                            write(*,*) beforestep 
                        !else if( r <= r3(n) ) then
                        !else if( beforestep >= 3 .and. beforestep /= 99 ) then                        
                        else if( beforestep >= 1 ) then
                            !sediment
                            zU = 2.8d-6 * rho
                            zTh = 10.7d-6 * rho
                            seconU = 2.8d-6 * rho
                            seconTh = 10.7d-6 * rho
                            write(*,*) beforestep 
                        endif !lower crust to sediment
                    endif !crust type
                endif !mantleからcrust
                !------------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------------
                !特異点での扱い、continentalをoceanicに
                !ここに追加
                !------------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------------
                !flux,massの計算、crustではさらに分割して計算
                if( beforestep <= 9 .and. beforestep >= 1) then 
                    !Crust
                    do c = 1 , htheta 
                        do d = 1 , hphi
                            do e = 0, hr - 1

                            !ddr = dr / dble(hr)
                            !r = r + dble(e-1) * ddr
                            r = r_pall + dble(e) * ddr
                            rep_r = sqrt( ( (r+ddr)**2 + (r+ddr)*r + r**2 ) / 3.0d0 ) 
                            theta = th_tem + dble(c-1) * ddtheta
                            phi = ph_tem + dble(d-1) * ddphi
                            
                            dis = SQRT( ( mh*sin(mtheta)*cos(mphi) - rep_r*sin(theta)*cos(phi) )**2.0d0 +& 
                            ( mh*sin(mtheta)*sin(mphi) - rep_r*sin(theta)*sin(phi) )**2.0d0 +&
                            ( mh*cos(mtheta) - rep_r*cos(theta) )**2.0d0  ) 
                            !dis = SQRT( ( mh*sin(mtheta)*cos(mphi) - r*sin(theta)*cos(phi) )**2.0d0 +& 
                            !( mh*sin(mtheta)*sin(mphi) - r*sin(theta)*sin(phi) )**2.0d0 +&
                            !( mh*cos(mtheta) - r*cos(theta) )**2.0d0  ) 
                            
                            !dis = SQRT( ( mh*sin(mtheta)*cos(mphi) - (r+dr/2.0d0)*sin(theta)*cos(phi) )**2.0d0 +& 
                            !( mh*sin(mtheta)*sin(mphi) - (r+dr/2.0d0)*sin(theta)*sin(phi) )**2.0d0 +&
                            !( mh*cos(mtheta) - (r+dr/2.0d0)*cos(theta) )**2.0d0  ) 

                            !dis = SQRT( ( mh*sin(mtheta)*cos(mphi) - (r+ddr/2.0d0)*sin(theta)*cos(phi) )**2.0d0 +& 
                            !( mh*sin(mtheta)*sin(mphi) - (r+ddr/2.0d0)*sin(theta)*sin(phi) )**2.0d0 +&
                            !( mh*cos(mtheta) - (r+ddr/2.0d0)*cos(theta) )**2.0d0  ) 

                            meccflux(1) = (A_U238 * n_U238 * pos * ccU(1) / (12.0d0 * pi * dis**2 ) &
                            * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) + (A_Th232 * n_Th232 &
                            * pos * ccTh(1) / (12.0d0 * pi * dis**2 ) * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            !meccflux(1) = (A_U238 * n_U238 * pos * ccU(1) / (12.0d0 * pi * dis**2 ) &
                            !* ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) + (A_Th232 * n_Th232 &
                            !* pos * ccTh(1) / (12.0d0 * pi * dis**2 ) * ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            meccflux(2) = (A_U238 * n_U238 * pos * ccU(2) / (12.0d0 * pi * dis**2 ) &
                            * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) + (A_Th232 * n_Th232 &
                            * pos * ccTh(2) / (12.0d0 * pi * dis**2 ) * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            !meccflux(2) = (A_U238 * n_U238 * pos * ccU(2) / (12.0d0 * pi * dis**2 ) &
                            !* ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) + (A_Th232 * n_Th232 &
                            !* pos * ccTh(2) / (12.0d0 * pi * dis**2 ) * ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            meccflux(3) = (A_U238 * n_U238 * pos * ccU(3) / (12.0d0 * pi * dis**2 ) &
                            * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) + (A_Th232 * n_Th232 &
                            * pos * ccTh(3) / (12.0d0 * pi * dis**2 ) * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            !meccflux(3) = (A_U238 * n_U238 * pos * ccU(3) / (12.0d0 * pi * dis**2 ) &
                            !* ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) + (A_Th232 * n_Th232 &
                            !* pos * ccTh(3) / (12.0d0 * pi * dis**2 ) * ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            meccfluxsum = meccflux(1) + meccflux(2) + meccflux(3)
                
                            Uflux = Uflux + (A_U238 * n_U238 * pos * zU / (12.0d0 * pi * dis**2 ) &
                            * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) 
                            !Uflux = Uflux + (A_U238 * n_U238 * pos * zU / (12.0d0 * pi * dis**2 ) &
                            !* ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) 
                            Thflux = Thflux + (A_Th232 * n_Th232 * pos * zTh / (12.0d0 * pi * dis**2 ) &
                            * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            !Thflux = Thflux + (A_Th232 * n_Th232 * pos * zTh / (12.0d0 * pi * dis**2 ) &
                            !* ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)

                            meUflux = (A_U238 * n_U238 * pos * zU / (12.0d0 * pi * dis**2 ) &
                            * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            !meUflux = (A_U238 * n_U238 * pos * zU / (12.0d0 * pi * dis**2 ) &
                            !* ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            meThflux = (A_Th232 * n_Th232 * pos * zTh / (12.0d0 * pi * dis**2 ) &
                            * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            !meThflux = (A_Th232 * n_Th232 * pos * zTh / (12.0d0 * pi * dis**2 ) &
                            !* ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)

                            meocflux = (A_U238 * n_U238 * pos * ocU / (12.0d0 * pi * dis**2 ) &
                            * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) + (A_Th232 * n_Th232 &
                            * pos * ocTh / (12.0d0 * pi * dis**2 ) * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            !meocflux = (A_U238 * n_U238 * pos * ocU / (12.0d0 * pi * dis**2 ) &
                            !* ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) + (A_Th232 * n_Th232 &
                            !* pos * ocTh / (12.0d0 * pi * dis**2 ) * ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            memaflux = (A_U238 * n_U238 * pos * maU / (12.0d0 * pi * dis**2 ) &
                            * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) + (A_Th232 * n_Th232 &
                            * pos * maTh / (12.0d0 * pi * dis**2 ) * ( (r+ddr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)
                            !memaflux = (A_U238 * n_U238 * pos * maU / (12.0d0 * pi * dis**2 ) &
                            !* ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi) + (A_Th232 * n_Th232 &
                            !* pos * maTh / (12.0d0 * pi * dis**2 ) * ( (r+dr)**3 - r**3 ) * sin(theta) * ddtheta * ddphi)

                            maUmass = maUmass + maU * ((r+ddr)**3 - r**3) * sin(theta) * ddtheta * ddphi / 3.0d0
                            maThmass = maThmass + maTh * ((r+ddr)**3 - r**3) * sin(theta) * ddtheta * ddphi / 3.0d0
                            conUmass = conUmass + ( ccU(1) + ccU(2) + ccU(3) ) *&
                             ((r+ddr)**3 - r**3) * sin(theta) * ddtheta * ddphi / 3.0d0
                            conThmass = conThmass + ( ccTh(1) + ccTh(2) + ccTh(3) ) *&
                             ((r+ddr)**3 - r**3) * sin(theta) * ddtheta * ddphi / 3.0d0
                            oceUmass = oceUmass + ( ocU ) *((r+ddr)**3 - r**3) * sin(theta) * ddtheta * ddphi / 3.0d0
                            oceThmass = oceThmass + ( ocTh ) *((r+ddr)**3 - r**3) * sin(theta) * ddtheta * ddphi / 3.0d0
                            conseUmass = conseUmass + seconU *((r+ddr)**3 - r**3) * sin(theta) * ddtheta * ddphi / 3.0d0
                            conseThmass = conseThmass + seconTh *((r+ddr)**3 - r**3) * sin(theta) * ddtheta * ddphi / 3.0d0
                            oceseUmass = oceseUmass + seoceU *((r+ddr)**3 - r**3) * sin(theta) * ddtheta * ddphi / 3.0d0
                            oceseThmass = oceseThmass + seoceTh *((r+ddr)**3 - r**3) * sin(theta) * ddtheta * ddphi / 3.0d0
                            !maUmass = maUmass + maU * r**2.0d0 * sin(theta) * ddr * ddtheta * ddphi
                            !maThmass = maThmass + maTh * r**2.0d0 * sin(theta) * ddr * ddtheta * ddphi
                            !conUmass = conUmass + ( ccU(1) + ccU(2) + ccU(3) + seU ) * r**2.0d0 * sin(theta) * ddr * ddtheta * ddphi
                            !conThmass = conThmass + ( ccTh(1) + ccTh(2) + ccTh(3) + seTh ) * r**2.0d0 * sin(theta) * ddr * ddtheta * ddphi

                            !maUmass = maUmass + maU * rep_r**2.0d0 * sin(theta) * ddr * ddtheta * ddphi
                            !maThmass = maThmass + maTh * rep_r**2.0d0 * sin(theta) * ddr * ddtheta * ddphi
                            !conUmass = conUmass + ( ccU(1) + ccU(2) + ccU(3) + seU ) * rep_r**2.0d0 * sin(theta) * ddr * ddtheta * ddphi
                            !conThmass = conThmass + ( ccTh(1) + ccTh(2) + ccTh(3) + seTh ) * rep_r**2.0d0 * sin(theta) * ddr * ddtheta * ddphi
                            
                            !maUmass = maUmass + maU * r**2.0d0 * sin(theta) * dr * ddtheta * ddphi
                            !maThmass = maThmass + maTh * r**2.0d0 * sin(theta) * dr * ddtheta * ddphi
                            !conUmass = conUmass + ( ccU(1) + ccU(2) + ccU(3) + seU ) * r**2.0d0 * sin(theta) * dr * ddtheta * ddphi
                            !conThmass = conThmass + ( ccTh(1) + ccTh(2) + ccTh(3) + seTh ) * r**2.0d0 * sin(theta) * dr * ddtheta * ddphi

                            !z = int(dis / 1.0d4) + 1
                            z = int(dis / (distance_flux * 1.0d3)) + 1
                            !if(dis < 1.0d5 ) then
                                if(z <= array_number ) then
                                    meflux(z) = meflux(z) + ( meUflux + meThflux )
                                    ocflux(z) = ocflux(z) + meocflux
                                    maflux(z) = maflux(z) + memaflux 
                                    ccflux(z) = ccflux(z) + meccfluxsum 
                                endif
                            !endif	

                            enddo !eについて、つまりrについて
                        enddo !dについて、つまりphiについて
                    enddo !cについて、つまりthetaについて
                else
                    !mantle
                    meccflux(1) = (A_U238 * n_U238 * pos * ccU(1) / (12.0d0 * pi * dis**2 ) &
                    * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi) + (A_Th232 * n_Th232 &
                    * pos * ccTh(1) / (12.0d0 * pi * dis**2 ) * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi)
                    meccflux(2) = abs(A_U238 * n_U238 * pos * ccU(2) / (12.0d0 * pi * dis**2 ) &
                    * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi) + abs(A_Th232 * n_Th232 &
                    * pos * ccTh(2) / (12.0d0 * pi * dis**2 ) * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi)
                    meccflux(3) = (A_U238 * n_U238 * pos * ccU(3) / (12.0d0 * pi * dis**2 ) &
                    * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi) + (A_Th232 * n_Th232 &
                    * pos * ccTh(3) / (12.0d0 * pi * dis**2 ) * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi)
                    meccfluxsum = meccflux(1) + meccflux(2) + meccflux(3)
                
                    Uflux = Uflux + (A_U238 * n_U238 * pos * zU / (12.0d0 * pi * dis**2 ) &
                    * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi) 
                    Thflux = Thflux + (A_Th232 * n_Th232 * pos * zTh / (12.0d0 * pi * dis**2 ) &
                    * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi)

                    meUflux = (A_U238 * n_U238 * pos * zU / (12.0d0 * pi * dis**2 ) &
                    * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi)
                    meThflux = (A_Th232 * n_Th232 * pos * zTh / (12.0d0 * pi * dis**2 ) &
                    * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi)

                    meocflux = (A_U238 * n_U238 * pos * ocU / (12.0d0 * pi * dis**2 ) &
                    * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi) + (A_Th232 * n_Th232 &
                    * pos * ocTh / (12.0d0 * pi * dis**2 ) * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi)
                    memaflux = (A_U238 * n_U238 * pos * maU / (12.0d0 * pi * dis**2 ) &
                    * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi) + (A_Th232 * n_Th232 &
                    * pos * maTh / (12.0d0 * pi * dis**2 ) * ( (r+dr)**3 - r**3 ) * sin(theta) * dtheta * dphi)

                    maUmass = maUmass + maU * ((r+dr)**3 - r**3) * sin(theta) * dtheta * dphi / 3.0d0
                    maThmass = maThmass + maTh * ((r+dr)**3 - r**3) * sin(theta) * dtheta * dphi / 3.0d0
                    conUmass = conUmass + ( ccU(1) + ccU(2) + ccU(3) ) *&
                     ((r+dr)**3 - r**3) * sin(theta) * dtheta * dphi / 3.0d0
                    conThmass = conThmass + ( ccTh(1) + ccTh(2) + ccTh(3) ) *&
                     ((r+dr)**3 - r**3) * sin(theta) * dtheta * dphi / 3.0d0
                    oceUmass = oceUmass + ( ocU ) *&
                     ((r+dr)**3 - r**3) * sin(theta) * dtheta * dphi / 3.0d0
                    oceThmass = oceThmass + ( ocTh ) *&
                      ((r+dr)**3 - r**3) * sin(theta) * dtheta * dphi / 3.0d0
                    conseUmass = conseUmass + seconU *&
                      ((r+dr)**3 - r**3) * sin(theta) * dtheta * dphi / 3.0d0
                    conseThmass = conseThmass + seconTh *&
                      ((r+dr)**3 - r**3) * sin(theta) * dtheta * dphi / 3.0d0
                    oceseUmass = oceseUmass + seoceU *&
                      ((r+dr)**3 - r**3) * sin(theta) * dtheta * dphi / 3.0d0
                    oceseThmass = oceseThmass + seoceTh *&
                      ((r+dr)**3 - r**3) * sin(theta) * dtheta * dphi / 3.0d0
                    !maUmass = maUmass + maU * r**2.0d0 * sin(theta) * dr * dtheta * dphi
                    !maThmass = maThmass + maTh * r**2.0d0 * sin(theta) * dr * dtheta * dphi
                    !conUmass = conUmass + ( ccU(1) + ccU(2) + ccU(3) + seU ) * r**2.0d0 * sin(theta) * dr * dtheta * dphi
                    !conThmass = conThmass + ( ccTh(1) + ccTh(2) + ccTh(3) + seTh ) * r**2.0d0 * sin(theta) * dr * dtheta * dphi

                    !maUmass = maUmass + maU * rep_r**2.0d0 * sin(theta) * dr * dtheta * dphi
                    !maThmass = maThmass + maTh * rep_r**2.0d0 * sin(theta) * dr * dtheta * dphi
                    !conUmass = conUmass + ( ccU(1) + ccU(2) + ccU(3) + seU ) * rep_r**2.0d0 * sin(theta) * dr * dtheta * dphi
                    !conThmass = conThmass + ( ccTh(1) + ccTh(2) + ccTh(3) + seTh ) * rep_r**2.0d0 * sin(theta) * dr * dtheta * dphi

                    !z = int(dis / 1.0d4) + 1
                    z = int(dis / (distance_flux * 1.0d3)) + 1
                    !if(dis < 1.0d5 ) then
                        if(z <= array_number ) then
                            meflux(z) = meflux(z) + ( meUflux + meThflux )
                            ocflux(z) = ocflux(z) + meocflux
                            maflux(z) = maflux(z) + memaflux
                            ccflux(z) = ccflux(z) + meccfluxsum
                        endif							
                    !endif
                endif !flux,massの計算、crustではさらに分割して計算
                !------------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------------
                !繰り返しでの値の初期化
                zU = 0.0d0
                zTh = 0.0d0
                z = 0
                maU = 0.0d0
                maTh = 0.0d0
                ocU = 0.0d0
                ocTh = 0.0d0
                ccU = 0.0d0
                ccTh = 0.0d0
                seconU = 0.0d0
                seconTh = 0.0d0
                seoceU = 0.0d0
                seoceTh = 0.0d0
                meocflux = 0.0d0
                memaflux = 0.0d0
                meccflux = 0.0d0
                meccfluxsum = 0.0d0
                crusttype = 0
                !------------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------------
                !mantleとcrustの場合分けに用いるbeforestep,afterstepの設定。なお、rの繰り返し最期のi = 515で初期化の必要
				!if( r9(n) >= r .and. r9(n) - r <= dr) then
				if( afterstep == 9) then
					beforestep = 8
                else if(afterstep == 8 ) then
					beforestep = 7
				else if(afterstep == 7 ) then
					beforestep = 6
				else if(afterstep == 6 ) then
					beforestep = 5
				else if(afterstep == 5 ) then
		    		beforestep = 4
				else if(afterstep == 4 ) then
					beforestep = 3
				else if(afterstep == 3 ) then
					beforestep = 2
				else if(afterstep == 2 ) then
					beforestep = 1
				else if(afterstep == 1 ) then
					beforestep = 99
				endif

				if(i == 515) then
					beforestep = 0
					afterstep = 0
				endif
                !------------------------------------------------------------------------------------------------
            enddo !iについて、つまりrについて
        enddo !jについて、つまりphiについて
    enddo !kについて、つまりthetaについて

    !------------------------------------------------------------------------------------------------
    !fluxの単位を置き換え
    do i = 1, array_number 
        !mer(i) = dble(i) * 1.0d4
        mer(i) = dble(i) * distance_flux
        meflux(i) = meflux(i) / 1.0d4
        ocflux(i) = ocflux(i) / 1.0d4
        maflux(i) = maflux(i) / 1.0d4
        ccflux(i) = ccflux(i) / 1.0d4
    enddo
    !------------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------------
    !累積フラックスを累積のものにする
    do i = 1, array_number - 1 
        meflux(i+1) = meflux(i+1) + meflux(i)
        ocflux(i+1) = ocflux(i+1) + ocflux(i)
        maflux(i+1) = maflux(i+1) + maflux(i)
        ccflux(i+1) = ccflux(i+1) + ccflux(i)
    enddo
    !------------------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------------------
    !ファイルに出力
    write(fo1,*) 'mesh r ,' , ' flux[/cm^2 s^1] ,' , ' mantle flux[/cm^2 s^1] ,' ,&
	 ' oceanic flux[/cm^2 s^1] ,' , ' continental flux[/cm^2 s^1] '
						
	do i = 1, array_number
		write(fo1,*) mer(i) , meflux(i) , maflux(i) , ocflux(i) , ccflux(i)
	enddo

	!write(fo1,*) ' mantle Umass/E+16 , ' , ' mantle Thmass/E+16 , ' , ' continental crust Umass/E+16 , ' ,&
    ! ' continental crust Thmass/E+16 , ' , ' oceanic crust Umass/E+16 , ' , ' oceanic crust Thmass/E+16 '
	!write(fo1,*) maUmass/1.0d16 , maThmass/1.0d16 , conUmass/1.0d16 , conThmass/1.0d16 , oceUmass/1.0d16 , oceThmass/1.0d16
	!write(fo1,*) ' mantle Umass/4.8d16 , ' , ' mantle Thmass/19.23d16 , ' , ' continental crust Umass/3.16d16 , ' ,&
	! ' continental crust Thmass/12.40d16 , ' , ' oceanic crust Umass/0.04d16 , ' , ' oceanic Thmass/0.09d16 , ' 
	!write(fo1,*) maUmass/4.8d16 , maThmass/19.23d16 , conUmass/3.16d16 , conThmass/12.40d16 , oceUmass/0.04d16 , oceThmass/0.09d16
    !write(fo1,*) ' continental sediment Umass/0.33d16 , ' , ' continental sediment Thmass/1.27d16 '
    !write(fo1,*) oceseUmass/0.33d16, oceseThmass/1.27d16
    !------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------
!ファイルをクローズ
close(fi1)
close(fi2)
close(fi3)
close(fi4)
close(fi5)
close(fi6)
close(fi7)
close(fi8)
close(fo1)
!------------------------------------------------------------------------------------------------
end program flux_calc