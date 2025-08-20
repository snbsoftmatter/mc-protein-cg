! Monte carlo code to simulate a polymer		
	character ::	atom(6000),residue(6000)*4,atomn(6000)
	integer ::		resnum(6000),natom(6000)
	integer,parameter:: N=5501		! Number of particle totalres
	integer,parameter :: N1 = 58            ! Number of particle in beads
	integer,parameter :: N2 = 5443           ! NUmber of particle in surrounding
	integer,parameter :: N3 = 36		!Number of hydrophilic bead  
	integer,parameter :: N4 = 22  !Number of hydrophobic bead
	real,parameter:: L=23.0 	              !box length
	integer	seed,Nstep
	real	theta
	integer nover,naccept,ntrial,iconfig
	real,dimension(N)::Rx,Ry,Rz,Ch
	real,dimension(N1)::phi,psi,vdi,vdai
	integer :: i,j
	real :: junk,junk1
	real :: x1,y1,val1
	real sigma,eps,dens,denb,Rcut,l0,Drmax,Dthmax,beta,V,Rmin
	real Rxijd,Ryijd,Rzijd,Rijdsq	
	real Rxold,Ryold,Rzold,phold,psold,Vold 
	real Rxnew,Rynew,Rznew,phnew,psnew,Vnew,deltaV,ratio,Vend	
	integer,parameter :: maxbin = 45
	integer,parameter:: max_bin=33
	real,dimension(max_bin)::bin1,bin2,bin3,bin4,bine1,bine2,bine3,bine4 
	real step 
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hypphps
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hpbphps
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsipsj		
	real :: vintra,a,b,x,y,p,q
	real :: xmax,xmin,ymax,ymin
	real :: f00,f01,f10,f11,delx,deltx,dely,delty
	real :: vinter,xp,yp,ap,bp
	real :: amax,amin,bmax,bmin
	integer :: a1,a2,a3,a4
	
	common / BLOCK1 / Rx, Ry, Rz
	common / BLOCK2 / hypphps,hpbphps
	common / BLOCK3 / phi,psi
	common / const /  residue
	common / BLOCK4 / vdi,vdai
	common / BLOCK5 / hyhyphiphj,hyhyphipsj,hyhypsiphj,hyhypsipsj
	common / BLOCK6 / hyhpphiphj,hyhpphipsj,hyhppsiphj,hyhppsipsj
	common / BLOCK7 / hphpphiphj,hphpphipsj,hphppsiphj,hphppsipsj
	
	data pi/3.1415927/ 
	sigma = 1.0
	eps=1.0
	dens = N2/(L**3)
	denb = N1/(L**3)
	step =L/(2*max_bin)
	Rcut = 2.5
	l0 = 1.85*sigma
	Drmax = 0.10*sigma ! maximum displacement of the atom
	Dthmax = 0.01
	beta = 1.0
	iconfig = 5 !step after write output cordinate
	iconfig2 = 5 ! Step after which solvent distribution is calculated
	V=0.0
	seed=1
	Rmin = 1.0
	Nstep=0
!......................................................................................

	open(7,file='initial-coordinate.dat') ! Initial coordiante
	open(22,file='model-dihed.dat') ! Initial dihed value for phi psi obtained from md
	
!!.........Intra residue interactions..............................................................
	open(39,file ='phi-psi-free-energy-hydrophilic.dat') ! File of hydrophilic energy profile ph-ps interaction 
	open(79,file ='phi-psi-free-energy-hydrophobic.dat') ! File of hydrophobic energy profile ph-ps interaction	
!........Inter residue interactions (Hydrophilic-hydrophilic)..............................................................
	open(111,file='free-energy-hydrophilic-hydrophilic-phiphj.dat') ! E profile hydrophilic phi-phj
	open(121,file='free-energy-hydrophilic-hydrophilic-phipsj.dat') ! E profile hydrophilic phi-psj
	open(131,file='free-energy-hydrophilic-hydrophilic-psiphj.dat') ! E profile hydrophilic psi-phj
	open(141,file='free-energy-hydrophilic-hydrophilic-psipsj.dat') ! E profile hydrophilic psi-psj

!........Inter residue interactions (Hydrophilic-hydrophobic)..............................................................
	open(211,file='free-energy-hydrophilic-hydrophobic-phiphj.dat') ! E profile hydrophilic hydrophobic phi-phj
	open(221,file='free-energy-hydrophilic-hydrophobic-phipsj.dat') ! E profile hydrophilic hydrophobic phi-psj
	open(231,file='free-energy-hydrophilic-hydrophobic-psiphj.dat') ! E profile hydrophilic hydrophobic psi-phj
	open(241,file='free-energy-hydrophilic-hydrophobic-psipsj.dat') ! E profile hydrophilic hydrophobic psi-psj	
	
!........Inter residue interactions (Hydropobic-hydrophobic)..............................................................
	open(311,file='free-energy-hydrophobic-hydrophobic-phiphj.dat') ! E profile hydrophilic hydrophobic phi-phj
	open(321,file='free-energy-hydrophobic-hydrophobic-phipsj.dat') ! E profile hydrophilic hydrophobic phi-psj
	open(331,file='free-energy-hydrophobic-hydrophobic-psiphj.dat') ! E profile hydrophilic hydrophobic psi-phj
	open(341,file='free-energy-hydrophobic-hydrophobic-psipsj.dat') ! E profile hydrophilic hydrophobic psi-psj	
	
!...........................................................................................................

	open(20,file='energy-step.dat') !potential for each step
	open(30,file="output-coordinate-bead.dat") ! output coordinate file
	open(40,file='acceptance-ratio-step.dat')
	open(55,file='dataset-dihed.data')
	open(50,file='radial-distribution-hydrophilic-bead.dat')
	open(60,file='radial-distribution-hydrophobic-bead.dat')
	open(70,file='radial-distribution-surrounding-bead.dat')
	open(80,file='radial-distribution-bead-bead.dat')
!.........................................................................................

!	Assign inital coordinate
	do i = 1,N
		read(7,*) junk,residue(i),Rx(i),Ry(i),Rz(i),Ch(i)
	end do
	
!	Assign initial dihed angle value
	do i = 1,N1
		read(22,*)junk1,phi(i),psi(i)
	end do
!.......................pair..correlation.................................
	do i = 1,max_bin    ! Initiallization
		bin1(i) = 0.0  ! Hydrophilic particle bead of polymer
		bin2(i) = 0.0  ! HYdrophobic particle bead of polymer
		bin3(i) = 0.0  ! Surrounding LJ particles
		bin4(i) = 0.0   !bead-bead
	end do
!................Intra residue free energy.............................

!	Read free energy profile for hydrophilic ph-ps interaction	

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(39,*) x1,y1,val1
			hypphps(i,j) = val1
		end do
	end do
		
!	Read free energy profile for hydrophobic ph-ps interaction	
	
	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(79,*) x1,y1,val1
			hpbphps(i,j) = val1
		end do
	end do



!................Inter residue free energy.............................

!.....Read E profile hydrophilic phi-phj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(111,*) x1,y1,val1
			hyhyphiphj(i,j) = val1
		end do
	end do

!.....Read E profile hydrophilic phi-psj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(121,*) x1,y1,val1
			hyhyphipsj(i,j) = val1
		end do
	end do

!.....Read E profile hydrophilic psi-phj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(131,*) x1,y1,val1
			hyhypsiphj(i,j) = val1
		end do
	end do

!.....Read E profile hydrophilic psi-psj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(141,*) x1,y1,val1
			hyhypsipsj(i,j) = val1
		end do
	end do

!.....Read E profile hydrophilic hydrophobic phi-phj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(211,*) x1,y1,val1
			hyhpphiphj(i,j) = val1
		end do
	end do

!.....Read E profile hydrophilic hydrophobic phi-psj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(221,*) x1,y1,val1
			hyhpphipsj(i,j) = val1
		end do
	end do
	
!.....Read E profile hydrophilic hydrophobic psi-phj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(231,*) x1,y1,val1
			hyhppsiphj(i,j) = val1
		end do
	end do
	
!.....Read E profile hydrophilic hydrophobic psi-psj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(241,*) x1,y1,val1
			hyhppsipsj(i,j) = val1
		end do
	end do

!.....Read E profile hydrophobic phi-phj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(311,*) x1,y1,val1
			hphpphiphj(i,j) = val1
		end do
	end do

!.....Read E profile hydrophobic phi-psj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(321,*) x1,y1,val1
			hphpphipsj(i,j) = val1
		end do
	end do

!.....Read E profile hydrophobic psi-phj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(331,*) x1,y1,val1
			hphppsiphj(i,j) = val1
		end do
	end do

!.....Read E profile hydrophobic psi-psj...................

	do i = -maxbin,maxbin
		do j = -maxbin,maxbin
			read(341,*) x1,y1,val1
			hphppsipsj(i,j) = val1
		end do
	end do

!....................................................................
!	Check initial overlapping	

	do i = 1,N-1
		do j = i+1,N
		
			Rxijd = Rx(j) - Rx(i)
			Ryijd = Ry(j) - Ry(i)
			Rzijd = Rz(j) - Rz(i)
			Rijdsq = Rxijd**2 + Ryijd**2 + Rzijd**2
			if (Rijdsq .lt. Rmin**2 ) then
				write(*,*)i,j, "Overlap"
			end if
		end do	
	end do


	
!.............................................................MC.................................................
! Monte carlo run of the bead particles
	
	istep1 = 10  !Number of monte carlo steps
	do i = 1,istep1 ! Loop on monte carlos steps
		ntrial = 0
		naccept = 0
		seed = seed + i
		call srand(seed)   ! Random number generator
			do j = 1,N
				if (j.le.N1) then
				phold = phi(j)
				psold = psi(j)
				end if 
				
				Rxold = Rx(j)
				Ryold = Ry(j)
				Rzold = Rz(j)
				
				call energy(phold,psold,Rxold,Ryold,Rzold,Ch,Vold,j,sigma,Rmin) 
				Rxnew = Rxold + (2*rand()-1)*Drmax 
				Rynew = Ryold + (2*rand()-1)*Drmax
				Rznew = Rzold + (2*rand()-1)*Drmax
				if (j .le. N1) then
				phnew = phold + (2*rand()-1)*Dthmax 
				psnew = psold + (2*rand()-1)*Dthmax
				if (phnew .gt. pi) then
					phnew = phnew -(2.0*pi)
				else if (phnew .lt. -pi) then
					phnew = phnew + (2.0*pi)
				end if
				if (psnew .gt. pi) then
					psnew = psnew -(2.0*pi)
				else if (psnew .lt. -pi) then
					psnew = psnew + (2.0*pi)
				end if
				
				end if
				
				call dis(Rxnew,Rynew,Rznew,j,nover,sigma,Rmin)
				
				if(nover .eq. N-1) then
				
				call energy(phnew,psnew,Rxnew,Rynew,Rznew,Ch,Vnew,j,sigma,Rmin) 
					
					deltaV = Vnew - Vold  
					 
					if(beta*deltaV .le. 75.0) then
						if(beta*deltaV .le. 0.0) then
							V = V + deltaV
							Rx(j) = Rxnew
							Ry(j) = Rynew
							Rz(j) = Rznew
							if (j .le. N1) then
							phi(j) = phnew
							psi(j) = psnew							
							end if						
							naccept = naccept + 1
						elseif(exp(-beta*deltaV) .ge. rand()) then
							V = V + deltaV
							Rx(j) = Rxnew
							Ry(j) = Rynew
							Rz(j) = Rznew
							if (j .le. N1) then
							phi(j) = phnew
							psi(j) = psnew							
							end if
							naccept = naccept +1
						end if
					end if
					
				end if

				ntrial = ntrial + 1
			end do
			ratio = real(naccept)/real(ntrial)
			write(40,*)i,ratio

			if(ratio .gt. 0.5) then   
				Drmax=Drmax*1.05
				Dthmax = Dthmax*1.05
				
			else 
				Drmax=Drmax* 0.95
				Dthmax = Dthmax*0.95
				
			end if
			call sum_energy(Vend,sigma,Rmin,Ch) ! Calculate energy after each MC step
1			FORMAT (1(F8.4,'	'))
			write(20,*)i,'	',vend
			if(i .ge. 0) then ! write output coordinate
				if (mod(i,iconfig).eq.0) then
					do k = 1,N
						write(30,101) 'BEAD',k,Rx(k),Ry(k),Rz(k) !output coordinates
					end do
					
					do k = 1,N1
						write(55,119) 'BEAD',k,phi(k),psi(k) !output phi psi
					end do
					
101					format(a4,5x,i5,5x,3(f8.3,5x))
119					format(a4,5x,i5,5x,2(f8.3,5x))
          			end if
          			if (mod(i,iconfig2).eq.0) then
          			 Nstep = Nstep + 1
          			call pair_correlation(i,bine1,bine2,bine3,bine4,dens,step)
					do k = 1, max_bin
						bin1(k) = bin1(k) + bine1(k)
						bin2(k) = bin2(k) + bine2(k)
						bin3(k) = bin3(k) + bine3(k)
						bin4(k) = bin4(k) + bine4(k)
					end do
				end if
          			
			end if 
	enddo
	
	do i=1,max_bin
      	write(50,*) i*step,'     ',bin1(i)/real(Nstep)
      	write(60,*) i*step,'     ',bin2(i)/real(Nstep)
      	write(70,*) i*step,'     ',bin3(i)/real(Nstep)
      	write(80,*) i*step,'     ',bin4(i)/real(Nstep)
      	enddo
	
	close(unit = 20)
	close(unit = 30)
	close(unit = 40)
	close(unit = 55)
	close(unit = 50)
	close(unit = 60)
	close(unit = 70)
	close(unit = 80)
	end
	
!..................Subroutines............................................................
!......potential energy of atom i with all the other atoms in the system..................

	subroutine energy(phj,psj,Rxj,Ryj,Rzj,Chj,Vj,j,sigma,Rmin)	 
	character ::	atom(6000),residue(6000)*4,atomn(6000) 
	integer,parameter::N=5501
	integer,parameter :: N1 = 58            
	integer,parameter :: N2 = 5443  
	integer,parameter :: maxbin = 45         
	real,parameter::L=23.0
	real,dimension(N)::Rx,Ry,Rz,Chj
	real,dimension(N1)::phi,psi,vdi,vdai
	integer i,j
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hypphps
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hpbphps
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsipsj
	real junk,junk1,sigma,Rmin,Rcut,Rcutsq,sigsq
	real fa,fb,V,l0,cons,lambda,dl,eps
	real Rxj,Ryj,Rzj,Rxij,Ryij,Rzij,Rijsq,Rij
	real Rinv,sr2,sr6,sr1,Vij
	real Vb,costheta,Va,phj,psj,Vj
	real xmax,xmin,ymax,ymin,vintra(N1),vinter(N1)	
	common / BLOCK1 / Rx, Ry, Rz
	common / BLOCK2 / hypphps,hpbphps
	common / BLOCK3 / phi,psi
	common / const /  residue
	common / BLOCK4 / vdi,vdai
	common / BLOCK5 / hyhyphiphj,hyhyphipsj,hyhypsiphj,hyhypsipsj
	common / BLOCK6 / hyhpphiphj,hyhpphipsj,hyhppsiphj,hyhppsipsj
	common / BLOCK7 / hphpphiphj,hphpphipsj,hphppsiphj,hphppsipsj
	
	sigma=1.0
	eps=1.0
	Rmin = 1.0
	Rcut = 2.5
	Rcutsq = Rcut*Rcut
	sigsq = sigma*sigma
	fb = 67.0*eps/sigsq
	fa = 35.0*eps
        V = 0.0
        l0 = 1.85*sigma
	cons=0.2
	lambda=0.1*sigma
	dl=1/lambda
	
	do i = 1,N
		if (i .ne. j) then
! Hydrophilic-Hydrophobic and vice versa interaction
			if (residue(j).eq.'HYPL'.and.residue(i).eq.'HYPB'.or.  &
			residue(j).eq.'HYPB'.and.residue(i).eq.'HYPL' .or.     &
			residue(j).eq.'HYPL'.and.residue(i).eq.'HYPL' .or.     &
			residue(j).eq.'HYPB'.and.residue(i).eq.'HYPB') then
				Rxij = Rxj - Rx(i)
				Ryij = Ryj - Ry(i)
				Rzij = Rzj - Rz(i)
				Rxij = Rxij  - L*anint(Rxij/L) 
				Ryij = Ryij  - L*anint(Ryij/L)
				Rzij = Rzij  - L*anint(Rzij/L)
				Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij				
				if(Rijsq .ge. Rmin*Rmin) then
					if(Rijsq .le. Rcutsq) then
						Rij = sqrt(Rijsq/sigsq)
						Rinv=1/Rij
						sr2 = sigsq/Rijsq
						sr6 = sr2*sr2*sr2
						Vij = sr6*sr6
						V = V + 4.0*Vij + cons*((Chj(i)*Chj(j))*Rinv)*exp(-Rij*dl)
					end if
				end if 
! Hydrophilic-solvent interaction
			else if (residue(j).eq.'HYPL'.and. residue(i).eq.'SOLV') then
				Rxij = Rxj - Rx(i)
				Ryij = Ryj - Ry(i)
				Rzij = Rzj - Rz(i)
				Rxij = Rxij  - L*anint(Rxij/L) 
				Ryij = Ryij  - L*anint(Ryij/L)
				Rzij = Rzij  - L*anint(Rzij/L)
				Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij				
				if(Rijsq .ge. Rmin*Rmin) then
					if(Rijsq .le. Rcutsq) then
						sr2 = sigsq/Rijsq
						sr6 = sr2*sr2*sr2
						Vij = sr6 * (sr6 - 1.0)
						Rij = sqrt(Rijsq/sigsq)
						Rinv=1/Rij
						V = V + 4.0*0.5*Vij + cons*((Chj(i)*Chj(j))*Rinv)*exp(-Rij*dl)
					end if
				end if		
! Hydrophobic-solvent interaction
			else if (residue(j).eq.'HYPB'.and. residue(i).eq.'SOLV') then
				Rxij = Rxj - Rx(i)
				Ryij = Ryj - Ry(i)
				Rzij = Rzj - Rz(i)
				Rxij = Rxij  - L*anint(Rxij/L) 
				Ryij = Ryij  - L*anint(Ryij/L)
				Rzij = Rzij  - L*anint(Rzij/L)
				Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij				
				if(Rijsq .ge. Rmin*Rmin) then
					if(Rijsq .le. Rcutsq) then
						sr1=sqrt(Rijsq/sigsq)
						Rinv=1/sr1
						Vij = 13.0*exp(-2.4*sr1)
						V = V + Vij + cons*((Chj(i)*Chj(j))*Rinv)*exp(-sr1*dl)
					end if
				end if
! solvent-hydrophilic interaction
			else if (residue(j).eq.'SOLV'.and.residue(i).eq.'HYPL') then				
				Rxij = Rxj - Rx(i)
				Ryij = Ryj - Ry(i)
				Rzij = Rzj - Rz(i)
				Rxij = Rxij  - L*anint(Rxij/L) 
				Ryij = Ryij  - L*anint(Ryij/L)
				Rzij = Rzij  - L*anint(Rzij/L)
				Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
				if(Rijsq .ge. Rmin*Rmin) then
					if(Rijsq .le. Rcutsq) then
						sr2 = sigsq/Rijsq
						sr6 = sr2*sr2*sr2
						Vij = sr6 * (sr6 - 1.0)
						Rij = sqrt(Rijsq/sigsq)
						Rinv=1/Rij
						V = V + 4.0*0.5*Vij + cons*((Chj(i)*Chj(j))*Rinv)*exp(-Rij*dl)
					end if
				end if		
! solvent-solvent interaction			
			else if (residue(j).eq.'SOLV'.and.residue(i).eq.'SOLV') then				
				Rxij = Rxj - Rx(i)
				Ryij = Ryj - Ry(i)
				Rzij = Rzj - Rz(i)
				Rxij = Rxij  - L*anint(Rxij/L) 
				Ryij = Ryij  - L*anint(Ryij/L)
				Rzij = Rzij  - L*anint(Rzij/L)
				Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
				if(Rijsq .ge. Rmin*Rmin) then
					if(Rijsq .le. Rcutsq) then
						sr2 = sigsq/Rijsq
						sr6 = sr2*sr2*sr2
						Vij = sr6 * (sr6 - 1.0)
						Rij = sqrt(Rijsq/sigsq)
						Rinv=1/Rij
						V = V + 4.0*Vij + cons*((Chj(i)*Chj(j))*Rinv)*exp(-Rij*dl)
					end if
				end if			
! solvent-hydrophobic interaction			
			
			else if (residue(j).eq.'SOLV'.and.residue(i).eq.'HYPB')then
				Rxij = Rxj - Rx(i)
				Ryij = Ryj - Ry(i)
				Rzij = Rzj - Rz(i)
				Rxij = Rxij  - L*anint(Rxij/L) 
				Ryij = Ryij  - L*anint(Ryij/L)
				Rzij = Rzij  - L*anint(Rzij/L)
				Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
				if(Rijsq .ge. Rmin*Rmin) then
					if(Rijsq .le. Rcutsq) then
						sr1=sqrt(Rijsq/sigsq)
						Rinv=1/sr1
						Vij = 13.0*exp(-2.4*sr1)
						V = V + Vij + cons*((Chj(i)*Chj(j))*Rinv)*exp(-sr1*dl)
					end if
				end if
			end if
		end if
	end do
	V = V

!   Harmonic potential due to bonding distance
	if (j.le.N1) then
		if (j.le.N1-1) then		
			Rxij1 = Rxj - Rx(j+1)
	    		Ryij1 = Ryj - Ry(j+1)
        		Rzij1 = Rzj - Rz(j+1)
       	 	Rxij1 = Rxij1  - L*anint(Rxij1/L) 
       	 	Ryij1 = Ryij1  - L*anint(Ryij1/L)
       	 	Rzij1 = Rzij1  - L*anint(Rzij1/L)       	 		
    		    	Rijsq1 = Rxij1*Rxij1 + Ryij1*Ryij1 + Rzij1*Rzij1               	
			 Vb = 0.5*fb*(sqrt(Rijsq1)-l0)**2
		end if
		
  	else
        Vb = 0.0  
	end if 
	
!   Harmonic potential due to bonding angle 	
	
	if (j.le.N1) then		
		if (j.le.N1-2) then
			Rxij1 = Rxj - Rx(j+1)
        		Ryij1 = Ryj - Ry(j+1)
        		Rzij1 = Rzj - Rz(j+1)
        		Rxij2 = Rx(j+1) - Rx(j+2)  
        		Ryij2 = Ry(j+1) - Ry(j+2)
        		Rzij2 = Rz(j+1) - Rz(j+2)
        		Rxij1 = Rxij1  - L*anint(Rxij1/L) 
        		Ryij1 = Ryij1  - L*anint(Ryij1/L)
        		Rzij1 = Rzij1  - L*anint(Rzij1/L)
        		Rxij2 = Rxij2  - L*anint(Rxij2/L) 
        		Ryij2 = Ryij2  - L*anint(Ryij2/L)
        		Rzij2 = Rzij2  - L*anint(Rzij2/L)
        		Rijsq1 = Rxij1*Rxij1 + Ryij1*Ryij1 + Rzij1*Rzij1
        		Rijsq2 = Rxij2*Rxij2 + Ryij2*Ryij2 + Rzij2*Rzij2
!        		costheta = (((Rxij1*Rxij2)+(Ryij1*Ryij2)+(Rzij1*Rzij2))/ &
!                         (sqrt(Rijsq1)*sqrt(Rijsq2)))
			theta = acos(((Rxij1*Rxij2)+(Ryij1*Ryij2)+(Rzij1*Rzij2))/ &
                         (sqrt(Rijsq1)*sqrt(Rijsq2)))
                       Va = 0.5*fa*(theta-1.5)**2
 
		end if

		
  		else
        	Va = 0.0  
	end if  
	
	if (j.le.N1) then
		call intradihedral(phj,psj,j,vintra)
		vdi(j)= vintra(j) 
		call interdihedral(phj,psj,j,vinter)
		vdai(j)= vinter(j)
	end if                    
		
      Vj = V + Vb + Va + vdi(j) + vdai(j)

         return
         end

!......................................overlap checking.....................................
      subroutine dis(Rxj,Ryj,Rzj,j,noverd,sigma,Rmin)      
      character :: 	atom(6000),residue(6000)*4,atomn(6000)
      integer,parameter::N=5501
      integer,parameter::N1=58 
      integer,parameter :: maxbin = 45
      real,dimension(N1)::phi,psi
      real,parameter::L=23.0
      real,dimension(N)::Rx,Ry,Rz
      integer noverd,i,j
      real junk1,Rmin,sigma,Rxj,Ryj,Rzj,Rxij,Ryij,Rzij,dis_R,Rijsq
      common / BLOCK1 / Rx, Ry, Rz
	sigma=1.0
	Rmin = 1.0
	noverd =0
	
	
	do i=1,N
		if (i .ne. j) then
			Rxij = Rxj - Rx(i)
			Ryij = Ryj - Ry(i)
			Rzij = Rzj - Rz(i)
			Rxij = Rxij  - L*anint(Rxij/L) 
			Ryij = Ryij  - L*anint(Ryij/L)
			Rzij = Rzij  - L*anint(Rzij/L)
			Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
			dis_R = sqrt(Rijsq)			
			if (dis_R .ge. Rmin) then
				noverd = noverd + 1
			end if
		end if
	end do
	return
	end
	
!..................................sum of energy..........................................
	subroutine sum_energy(Vd,sigma,Rmin,Chj)	
	character ::	atom(6000),residue(6000)*4,atomn(6000)
	integer,parameter::N=5501
	integer,parameter :: N1 = 58            
	integer,parameter :: N2 = 5443   
	integer,parameter :: maxbin = 45        
	real,parameter::L=23.0
	real,dimension(N)::Rx,Ry,Rz,Chj
	real,dimension(N1)::phi,psi,vdi,vdai
	integer i,j
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hypphps
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hpbphps
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsipsj	
	real junk,junk1,sigma,Rmin,Rcut,Rcutsq
	real sigsq,Rminsq,l0,fb,fa
	real V,Vb,Va,Vc,cons,lambda,dl
	real Rxi,Ryi,Rzi,Rxij,Ryij,Rzij
	real Rijsq,Rij,Rinv,sr2,sr6,sr1,Vij
	real costheta,vintra,vinter,phdi,psdi,Vd		
	common / BLOCK1 / Rx, Ry, Rz
	common / BLOCK2 / hypphps,hpbphps
	common / BLOCK3 / phi,psi
	common / const /  residue
	common / BLOCK4 / vdi,vdai
	common / BLOCK5 / hyhyphiphj,hyhyphipsj,hyhypsiphj,hyhypsipsj
	common / BLOCK6 / hyhpphiphj,hyhpphipsj,hyhppsiphj,hyhppsipsj
	common / BLOCK7 / hphpphiphj,hphpphipsj,hphppsiphj,hphppsipsj
	
	
	sigma=1.0
	eps=1.0
	Rmin = 1.0
	Rcut = 2.5
	Rcutsq = Rcut*Rcut
	sigsq  = sigma*sigma
	Rminsq = Rmin*Rmin
	l0 = 1.85*sigma
	fb = 67.0*eps/sigsq
	fa = 35.0*eps
	V=0.0
	Vb=0.0
	Va=0.0
	Vc=0.0
	cons=0.2
	lambda=0.1*sigma
	dl=1/lambda
	
	
	
	do i = 1,N-1
		if (residue(i) .eq. 'HYPB') then
			Rxi = Rx(i)
			Ryi = Ry(i)
			Rzi = Rz(i)
				do j = i+1,N				        
					if (residue(j).eq.'HYPB'.or.residue(j).eq.'HYPL') then
						Rxij = Rxi - Rx(j)    
						Ryij = Ryi - Ry(j)
						Rzij = Rzi - Rz(j)
						Rxij = Rxij  - L*anint(Rxij/L) 
						Ryij = Ryij  - L*anint(Ryij/L) 
						Rzij = Rzij  - L*anint(Rzij/L)
						Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
						if (Rijsq .ge. Rminsq) then
							if (Rijsq .le. Rcutsq) then
								Rij = sqrt(Rijsq/sigsq)
								Rinv=1/Rij
								sr2 = sigsq/Rijsq
								sr6 = sr2*sr2*sr2		
								Vij = sr6*sr6
								V = V + 4.0*Vij + cons*((Chj(i)*Chj(j))*Rinv)*exp(-Rij*dl)
							end if
						end if
					end if	
				end do
				
				
	        else if (residue(i) .eq. 'HYPB') then
			Rxi = Rx(i)
			Ryi = Ry(i)
			Rzi = Rz(i)
				do j = i+1,N
					if (residue(j).eq.'SOLV') then						
						Rxij = Rxi - Rx(j)   
						Ryij = Ryi - Ry(j)
						Rzij = Rzi - Rz(j)
						Rxij = Rxij  - L*anint(Rxij/L) 
						Ryij = Ryij  - L*anint(Ryij/L) 
						Rzij = Rzij  - L*anint(Rzij/L)
						Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
						if (Rijsq .ge. Rminsq) then
							if (Rijsq .le. Rcutsq) then
								sr1=sqrt(Rijsq/sigsq)
								Rinv=1/sr1
								Vij = 13.0*exp(-2.4*sr1)
								V = V + Vij + cons*((Chj(i)*Chj(j))*Rinv)*exp(-sr1*dl)
							end if
						end if
					end if	
				end do
		else if (residue(i) .eq. 'HYPL') then
			Rxi = Rx(i)
			Ryi = Ry(i)
			Rzi = Rz(i)
				do j = i+1,N
					if (residue(j) .eq. 'HYPB' .or.residue(j).eq. 'HYPL') then
						Rxij = Rxi - Rx(j)    
						Ryij = Ryi - Ry(j)
						Rzij = Rzi - Rz(j)
						Rxij = Rxij  - L*anint(Rxij/L) 
						Ryij = Ryij  - L*anint(Ryij/L) 
						Rzij = Rzij  - L*anint(Rzij/L)
						Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
						if (Rijsq .ge. Rminsq) then
							if (Rijsq .le. Rcutsq) then
								Rij = sqrt(Rijsq/sigsq)
								Rinv=1/Rij
								sr2 = sigsq/Rijsq
								sr6 = sr2*sr2*sr2		
								Vij = sr6*sr6
								V = V + 4.0*Vij + cons*((Chj(i)*Chj(j))*Rinv)*exp(-Rij*dl)
							end if
						end if
					else					      
						Rxij = Rxi - Rx(j)    
						Ryij = Ryi - Ry(j)
						Rzij = Rzi - Rz(j)
						Rxij = Rxij  - L*anint(Rxij/L) 
						Ryij = Ryij  - L*anint(Ryij/L) 
						Rzij = Rzij  - L*anint(Rzij/L)
						Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
						if (Rijsq .ge. Rminsq) then
							if (Rijsq .le. Rcutsq) then
								sr2 = sigsq/Rijsq
								sr6 = sr2*sr2*sr2		
								Vij = sr6*(sr6 - 1.0)
								Rij = sqrt(Rijsq/sigsq)
								Rinv=1/Rij
								V = V + 4.0*0.5*Vij + cons*((Chj(i)*Chj(j))*Rinv)*exp(-Rij*dl)
							end if
						end if
					end if	
				end do
		else if (residue(i) .eq. 'SOLV') then
			Rxi = Rx(i)
			Ryi = Ry(i)
			Rzi = Rz(i)
				do j = i+1,N					
					Rxij = Rxi - Rx(j)    
					Ryij = Ryi - Ry(j)
					Rzij = Rzi - Rz(j)
					Rxij = Rxij  - L*anint(Rxij/L) 
					Ryij = Ryij  - L*anint(Ryij/L) 
					Rzij = Rzij  - L*anint(Rzij/L)
					Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
					if (Rijsq .ge. Rminsq) then
						if (Rijsq .le. Rcutsq) then
							sr2 = sigsq/Rijsq
							sr6 = sr2*sr2*sr2		
							Vij = sr6*(sr6 - 1.0)
							Rij = sqrt(Rijsq/sigsq)
							Rinv=1/Rij
							V = V + 4.0*Vij + cons*((Chj(i)*Chj(j))*Rinv)*exp(-Rij*dl)
						end if
					end if
				end do		
		end if
			
	end do
	V = V
	
!   Harmonic potential due to bonding distance

	do i = 1,N1-1
		Rxi = Rx(i)
		Ryi = Ry(i)
        	Rzi = Rz(i)
	    		Rxij1 = Rxi - Rx(i+1)
	    		Ryij1 = Ryi - Ry(i+1)
        		Rzij1 = Rzi - Rz(i+1)
       		Rxij1 = Rxij1  - L*anint(Rxij1/L) 
       		Ryij1 = Ryij1  - L*anint(Ryij1/L)
       		Rzij1 = Rzij1  - L*anint(Rzij1/L)       	 		
    			Rijsq1 = Rxij1*Rxij1 + Ryij1*Ryij1 + Rzij1*Rzij1    			               	
				Vb = Vb + 0.5*fb*(sqrt(Rijsq1)-l0)**2
	end do
			        
!   Harmonic potential due to bonding angle 

	do i = 1,N1-2
		Rxi = Rx(i)
        	Ryi = Ry(i)
        	Rzi = Rz(i) 
        		Rxij1 = Rxi - Rx(i+1)  
        		Ryij1 = Ryi - Ry(i+1)
        		Rzij1 = Rzi - Rz(i+1)			
			Rxij2 = Rx(i+1) - Rx(i+2)  
        		Ryij2 = Ry(i+1) - Ry(i+2)
        		Rzij2 = Rz(i+1) - Rz(i+2)
        		Rxij1 = Rxij1  - L*anint(Rxij1/L) 
        		Ryij1 = Ryij1  - L*anint(Ryij1/L) 
        		Rzij1 = Rzij1  - L*anint(Rzij1/L)   
        		Rxij2 = Rxij2  - L*anint(Rxij2/L) 
        		Ryij2 = Ryij2  - L*anint(Ryij2/L) 
        		Rzij2 = Rzij2  - L*anint(Rzij2/L)
        		Rijsq1 = Rxij1*Rxij1 + Ryij1*Ryij1 + Rzij1*Rzij1   
			Rijsq2 = Rxij2*Rxij2 + Ryij2*Ryij2 + Rzij2*Rzij2
        		theta = acos(((Rxij1*Rxij2)+(Ryij1*Ryij2)+(Rzij1*Rzij2))/ &
                         (sqrt(Rijsq1)*sqrt(Rijsq2)))
                       Va = Va + 0.5*fa*(theta-1.5)**2
     end do
     
     	do i=1,N1
     		phdi=phi(i)
     		psdi=psi(i)
     			call intradihedral(phdi,psdi,i,vdi)
     			call interdihedral(phdi,psdi,i,vdai)
     			Vc = Vc + vdi(i) + vdai(i)
     	end do

     	Vd = V +Vb+Va +Vc
     	
     	return
     	end



!.................................pair_correlation...................................................
	subroutine pair_correlation(im,hist1,hist2,hist3,hist4,dens,step)	
	character :: 	atom(6000),residue(6000)*4,atomn(6000)
	integer,parameter::N=5501
      	integer,parameter::N1 = 58
      	integer,parameter::N2 = 5443
      	integer,parameter :: N3 = 36		  
	integer,parameter :: N4 = 22        	  
      	real,dimension(N1)::phi,psi
	real,parameter:: L=23.0	
	real,dimension(N)::Rx,Ry,Rz
	integer,parameter :: maxbin = 45
	integer im,bin,Nstep
	integer,parameter:: max_bin=33
	real,dimension(max_bin)::binc1,binc2,binc3,binc4,hist1,hist2,hist3,hist4
	real junk1,Rxi,Ryi,Rzi,Rxij,Ryij,Rzij,Rijsq,disij,dens,denb
	real Rmax,Rmin,step,nideal,nidealb
	
	common / BLOCK1 / Rx, Ry, Rz
	common /const/ residue
	data pi/3.1415927/
	
      	dens = N2/(L**3)
      	denb = N1/(L**3)
	step =L/(2*max_bin)
	
	do i =1,max_bin 
	binc1(i) = 0.0
	binc2(i) = 0.0
	binc3(i) = 0.0
	binc4(i) = 0.0
	enddo
	
	do i = 1,N-1
		Rxi = Rx(i)
		Ryi = Ry(i)
		Rzi = Rz(i)
		
		
		if (residue(i) .eq. 'HYPL') then				
			do j = N1+1,N
				Rxij = Rxi - Rx(j)
				Ryij = Ryi - Ry(j)
				Rzij = Rzi - Rz(j)
				Rxij = Rxij  - L*anint(Rxij/L) 
				Ryij = Ryij  - L*anint(Ryij/L)
				Rzij = Rzij  - L*anint(Rzij/L)
				Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
				 disij = sqrt(Rijsq)
				if (disij .le. L/2.0) then		 
				  bin  = int(disij/step) 
				  binc1(bin) = binc1(bin) + 1
				endif
			enddo
		else if (residue(i) .eq. 'HYPB' ) then
			do j = N1+1,N
				Rxij = Rxi - Rx(j)
				Ryij = Ryi - Ry(j)
				Rzij = Rzi - Rz(j)
				Rxij = Rxij  - L*anint(Rxij/L) 
				Ryij = Ryij  - L*anint(Ryij/L)
				Rzij = Rzij  - L*anint(Rzij/L)
				Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
				 disij = sqrt(Rijsq)
				if (disij .le. L/2.0) then				 
				  bin  = int(disij/step) 
				  binc2(bin) = binc2(bin) + 1
				endif
			enddo
		else if (residue(i) .eq. 'SOLV') then
			do j = i+1,N
				Rxij = Rxi - Rx(j)
				Ryij = Ryi - Ry(j)
				Rzij = Rzi - Rz(j)
				Rxij = Rxij  - L*anint(Rxij/L) 
				Ryij = Ryij  - L*anint(Ryij/L)
				Rzij = Rzij  - L*anint(Rzij/L)
				Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
				 disij = sqrt(Rijsq)
				if (disij .le. L/2.0) then				 
				  bin  = int(disij/step) 
				  binc3(bin) = binc3(bin) + 2
				endif
			enddo
		end if
	end do
			
	do i = 1,N1-1
		Rxi = Rx(i)
		Ryi = Ry(i)
		Rzi = Rz(i)
				
		if (residue(i) .eq. 'HYPL'.or.residue(i).eq. 'HYPB') then				
			do j = i+1,N1
				Rxij = Rxi - Rx(j)
				Ryij = Ryi - Ry(j)
				Rzij = Rzi - Rz(j)
				Rxij = Rxij  - L*anint(Rxij/L) 
				Ryij = Ryij  - L*anint(Ryij/L)
				Rzij = Rzij  - L*anint(Rzij/L)
				Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
				 disij = sqrt(Rijsq)
				if (disij .le. L/2.0) then		 
				  bin  = int(disij/step) 
				  binc4(bin) = binc4(bin) + 2
				endif
			enddo				
		end if
	enddo
	do i =1,max_bin 
		Rmin = real(i-1)*step
		Rmax = Rmin + step
		nideal = 4*pi*dens*(Rmax**3 - Rmin**3)/3
		hist1(i) = binc1(i)/nideal/real(N3)
		hist2(i) = binc2(i)/nideal/real(N4)
		hist3(i) = binc3(i)/nideal/real(N2)
	enddo
	do i =1,max_bin 
		Rmin = real(i-1)*step
		Rmax = Rmin + step
		nidealb = 4*pi*denb*(Rmax**3 - Rmin**3)/3
		hist4(i) = binc4(i)/nidealb/real(N1)
	enddo
	return
	end



 
!................................Intraresidue dihedral coupling................................................      
	
	subroutine intradihedral(x,y,j,v11)	
	character ::	atom(6000),residue(6000)*4
        integer,parameter::N = 5501 
        real,parameter :: binvalue=0.070 
        integer,parameter::N1 = 58
        real,dimension(N)::Rx,Ry,Rz
	real,dimension(N1)::phi,psi,v11		
	real :: junk1
  	integer,parameter :: maxbin = 45
  	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hypphps
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hpbphps
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsipsj
	real x,y,xmax,xmin,ymax,ymin,a,b,p,q,fphps
	integer a1,a2,a3,a4 
  	real f00,f01,f10,f11,delx,deltx,dely,delty 	
  	common / BLOCK1 / Rx, Ry, Rz
  	common / BLOCK2 / hypphps,hpbphps
	common / BLOCK3 / phi,psi
	common / const /  residue
	common / BLOCK5 / hyhyphiphj,hyhyphipsj,hyhypsiphj,hyhypsipsj
	common / BLOCK6 / hyhpphiphj,hyhpphipsj,hyhppsiphj,hyhppsipsj
	common / BLOCK7 / hphpphiphj,hphpphipsj,hphppsiphj,hphppsipsj
  	
  	fphps=0.0
    	a = (x/binvalue)
  	if (a.gt.0.0) then
  	xmax = a*binvalue + binvalue
  	xmin = xmax - binvalue
  	else
  	xmin = a*binvalue - binvalue
  	xmax = xmin + binvalue
  	end if
  	b = (y/binvalue)
  	if (b.gt.0.0) then
  	ymax = b*binvalue + binvalue
  	ymin = ymax - binvalue
  	else
  	ymin = b*binvalue - binvalue
  	ymax = ymin + binvalue
  	end if
  	a1 = xmax/binvalue
	a2 = xmin/binvalue
	a3 = ymax/binvalue
	a4 = ymin/binvalue
	
	
	if (j.le.N1) then	
		if (residue(j).eq.'HYPB') then
			f00 = hpbphps(a2,a4)
			f01 = hpbphps(a2,a3)
			f10 = hpbphps(a1,a4)
			f11 = hpbphps(a1,a3)
			delx = xmax - x
			deltx = xmax-xmin
			dely = ymax - y
			delty = ymax-ymin
			p = delx/deltx
			q = dely/delty
			fphps =   fphps + p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00 
			
		else if (residue(j).eq.'HYPL') then
			f00 = hypphps(a2,a4)
			f01 = hypphps(a2,a3)
			f10 = hypphps(a1,a4)
			f11 = hypphps(a1,a3)
			delx = xmax - x
			deltx = xmax-xmin
			dely = ymax - y
			delty = ymax-ymin
			p = delx/deltx
			q = dely/delty
			fphps =  fphps +  p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00 
	
		end if
	end if
	v11 = fphps	
	return
	end

 	  	
!......... Subroutine to calculate interploation and energy of interresidue dihed   

	subroutine interdihedral(x,y,j,v12)	
	character :: 	atom(6000),residue(6000)*4
	real,parameter::L=23.0
	integer,parameter::N = 5501
	integer,parameter::N1 = 58
	real,parameter :: binvalue=0.070 
	real,dimension(N1)::phi,psi,v12,fenergy
	real :: junk1
	integer,parameter :: maxbin = 45
	real,dimension(N)::Rx,Ry,Rz
	real Rxi,Ryi,Rzi,Rxij,Ryij,Rzij,Rijsq,disij
  	real x,y,xmax,xmin,ymax,ymin,a,b,v,p,q,fphps,fpsph
  	real f00,f01,f10,f11,delx,deltx,dely,delty,fxy  	
  	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hypphps
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hpbphps
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsipsj  		
  	common / BLOCK1 / Rx, Ry, Rz
  	common / BLOCK2 / hypphps,hpbphps
	common / BLOCK3 / phi,psi
	common / const / residue
	common / BLOCK5 / hyhyphiphj,hyhyphipsj,hyhypsiphj,hyhypsipsj
	common / BLOCK6 / hyhpphiphj,hyhpphipsj,hyhppsiphj,hyhppsipsj
	common / BLOCK7 / hphpphiphj,hphpphipsj,hphppsiphj,hphppsipsj
    	
    	
    
    	Rxi = Rx(j)
    	Ryi = Ry(j)
    	Rzi = Rz(j)
    	
    	do i = 1,N1
        	if (j .ne. i) then
            		Rxij = Rxi - Rx(i)
            		Ryij = Ryi - Ry(i)
            		Rzij = Rzi - Rz(i)
            		Rxij = Rxij  - L*anint(Rxij/L)
            		Ryij = Ryij  - L*anint(Ryij/L)
            		Rzij = Rzij  - L*anint(Rzij/L)
            		Rijsq = Rxij*Rxij + Ryij*Ryij + Rzij*Rzij
            		disij = sqrt(Rijsq)
            		if (disij .le. 2.59) then
                		a = phi(i)
                		b = psi(i)   
                		call interpolation(x,y,a,b,j,i,fenergy)   
                		v12 =  fenergy
                	end if
                end if
	end do
	return
	end 
                    
                    
!....................................................................................                    
                                  
	subroutine interpolation(x,y,a,b,j,i,fenergy)	
	character :: atom(6000),residue(6000)*4 
	integer,parameter::N = 5501
	integer,parameter::N1 = 58
        real,dimension(N)::Rx,Ry,Rz
        real,dimension(N1)::phi,psi,fenergy
        real,parameter :: binvalue=0.070 
        integer i,j
        real xp,x,xmax,xmin,yp,y,ymax,ymin
        integer xy1,xy2,xy3,xy4
        real ap,a,amax,amin,bp,b,bmax,bmin
        integer ab1,ab2,ab3,ab4
        real f00,f01,f10,f11,p,q              
        real fphiphj,fphipsj,fpsiphj,fpsipsj
        integer,parameter :: maxbin = 45
        real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hypphps
        real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hpbphps
        real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhyphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhypsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hyhppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphpphipsj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsiphj
	real,dimension(-maxbin:maxbin,-maxbin:maxbin) :: hphppsipsj
	
        common / BLOCK1 / Rx, Ry, Rz
        common / BLOCK2 / hypphps,hpbphps
        common / const / residue
        common / BLOCK5 / hyhyphiphj,hyhyphipsj,hyhypsiphj,hyhypsipsj
        common / BLOCK6 / hyhpphiphj,hyhpphipsj,hyhppsiphj,hyhppsipsj
        common / BLOCK7 / hphpphiphj,hphpphipsj,hphppsiphj,hphppsipsj    
                    
        xp = (x/binvalue)
        if (xp.gt.0.0) then
        	xmax = xp*binvalue + binvalue
                xmin = xmax - binvalue
        else
   	        xmin = xp*binvalue - binvalue
                xmax = xmin + binvalue
        end if
        yp = (y/binvalue)
        if (yp.gt.0.0) then
                ymax = yp*binvalue + binvalue
                ymin = ymax - binvalue
        else
                ymin = yp*binvalue - binvalue
                ymax = ymin + binvalue
        end if
            	xy1 = xmax/binvalue
            	xy2 = xmin/binvalue
            	xy3 = ymax/binvalue
            	xy4 = ymin/binvalue
            	
        ap = (a/binvalue)
        if (ap.gt.0.0) then
        	amax = ap*binvalue + binvalue
                amin = amax - binvalue
        else
   	        amin = ap*binvalue - binvalue
                amax = amin + binvalue
        end if
        bp = (b/binvalue)
        if (bp.gt.0.0) then
                bmax = bp*binvalue + binvalue
                bmin = bmax - binvalue
        else
                bmin = bp*binvalue - binvalue
                bmax = bmin + binvalue
        end if
            	ab1 = amax/binvalue
            	ab2 = amin/binvalue
            	ab3 = bmax/binvalue
            	ab4 = bmin/binvalue
            	
 if (j.le.N1) then
 
     if (i.ne.j) then
     
  	if (residue(j).eq.'HYPL'.and. residue(i).eq.'HYPL') then

	f00 = hyhyphiphj(xy2,ab2)
	f01 = hyhyphiphj(xy2,ab1)
	f10 = hyhyphiphj(xy1,ab2)
	f11 = hyhyphiphj(xy1,ab1)
	delx = xmax - x
	deltx = xmax - xmin
	dely = amax - a
	delty = amax - amin
	p = delx/deltx
	q = dely/delty
	fphiphj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	

	f00 = hyhyphipsj(xy2,ab4)
	f01 = hyhyphiphj(xy2,ab3)
	f10 = hyhyphiphj(xy1,ab4)
	f11 = hyhyphiphj(xy1,ab3)
	delx = xmax - x
	deltx = xmax - xmin
	dely = bmax - b
	delty = bmax - bmin
	p = delx/deltx
	q = dely/delty
	fphipsj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	


	f00 = hyhypsiphj(xy4,ab2)
	f01 = hyhypsiphj(xy4,ab1)
	f10 = hyhypsiphj(xy3,ab2)
	f11 = hyhypsiphj(xy3,ab1)
	delx = ymax - y
	deltx = ymax - ymin
	dely = amax - a
	delty = amax - amin
	p = delx/deltx
	q = dely/delty
	fpsiphj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	


	f00 = hyhypsipsj(xy4,ab4)
	f01 = hyhypsipsj(xy4,ab3)
	f10 = hyhypsipsj(xy3,ab4)
	f11 = hyhypsipsj(xy3,ab3)
	delx = ymax - y
	deltx = ymax - ymin
	dely = bmax - b
	delty = bmax - bmin
	p = delx/deltx
	q = dely/delty
	fpsipsj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	
	fenergy = fphiphj + fphipsj + fpsiphj + fpsipsj
		
	elseif (residue(j).eq.'HYPL'.and.residue(i).eq.'HYPB'.or. &
			residue(j).eq.'HYPB'.and.residue(i).eq.'HYPL') then
  	
  	

	f00 = hyhpphiphj(xy2,ab2)
	f01 = hyhpphiphj(xy2,ab1)
	f10 = hyhpphiphj(xy1,ab2)
	f11 = hyhpphiphj(xy1,ab1)
	delx = xmax - x
	deltx = xmax - xmin
	dely = amax - a
	delty = amax - amin
	p = delx/deltx
	q = dely/delty
	fphiphj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	


	f00 = hyhpphipsj(xy2,ab4)
	f01 = hyhpphiphj(xy2,ab3)
	f10 = hyhpphiphj(xy1,ab4)
	f11 = hyhpphiphj(xy1,ab3)
	delx = xmax - x
	deltx = xmax - xmin
	dely = bmax - b
	delty = bmax - bmin
	p = delx/deltx
	q = dely/delty
	fphipsj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	


	f00 = hyhppsiphj(xy4,ab2)
	f01 = hyhppsiphj(xy4,ab1)
	f10 = hyhppsiphj(xy3,ab2)
	f11 = hyhppsiphj(xy3,ab1)
	delx = ymax - y
	deltx = ymax - ymin
	dely = amax - a
	delty = amax - amin
	p = delx/deltx
	q = dely/delty
	fpsiphj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	


	f00 = hyhppsipsj(xy4,ab4)
	f01 = hyhppsipsj(xy4,ab3)
	f10 = hyhppsipsj(xy3,ab4)
	f11 = hyhppsipsj(xy3,ab3)
	delx = ymax - y
	deltx = ymax - ymin
	dely = bmax - b
	delty = bmax - bmin
	p = delx/deltx
	q = dely/delty
	fpsipsj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	
	fenergy = fphiphj + fphipsj + fpsiphj + fpsipsj
	
	else if (residue(j).eq.'HYPB' .and. residue(i).eq.'HYPB') then

  	

	f00 = hphpphiphj(xy2,ab2)
	f01 = hphpphiphj(xy2,ab1)
	f10 = hphpphiphj(xy1,ab2)
	f11 = hphpphiphj(xy1,ab1)
	delx = xmax - x
	deltx = xmax - xmin
	dely = amax - a
	delty = amax - amin
	p = delx/deltx
	q = dely/delty
	fphiphj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	


	f00 = hphpphipsj(xy2,ab4)
	f01 = hphpphiphj(xy2,ab3)
	f10 = hphpphiphj(xy1,ab4)
	f11 = hphpphiphj(xy1,ab3)
	delx = xmax - x
	deltx = xmax - xmin
	dely = bmax - b
	delty = bmax - bmin
	p = delx/deltx
	q = dely/delty
	fphipsj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	


	f00 = hphppsiphj(xy4,ab2)
	f01 = hphppsiphj(xy4,ab1)
	f10 = hphppsiphj(xy3,ab2)
	f11 = hphppsiphj(xy3,ab1)
	delx = ymax - y
	deltx = ymax - ymin
	dely = amax - a
	delty = amax - amin
	p = delx/deltx
	q = dely/delty
	fpsiphj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	


	f00 = hphppsipsj(xy4,ab4)
	f01 = hphppsipsj(xy4,ab3)
	f10 = hphppsipsj(xy3,ab4)
	f11 = hphppsipsj(xy3,ab3)
	delx = ymax - y
	deltx = ymax - ymin
	dely = bmax - b
	delty = bmax - bmin
	p = delx/deltx
	q = dely/delty
	fpsipsj = p*q*f11 + q*(1-p)*f10 + p*(1-q)*f01 + (1-p-q+p*q)*f00
	
	fenergy = fphiphj + fphipsj + fpsiphj + fpsipsj	 
  	end if
  	
   end if  
  end if
  	
  return
  end

!  	
!  	
!  	
