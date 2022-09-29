	program oscillatore
	
	real, dimension(:), allocatable :: campo
	common/vario/delta
	common/flag/misu,i_corr
	
	call cpu_time(start)
	call ranstart

	open(1, file='input.txt', status='old')		!file coi parametri
	open(2, file='adatiy2c.dat', status='unknown')	!file coi risultati
	open(3, file='ay2toy2c.dat', status='unknown')
	!open(3, file='cammino1.dat', status='unknown')
	!open(4, file='cammino15.dat', status='unknown')
	
	read(1,*) misure         !numero di misure
        read(1,*) i_dec          !decorrelazione fra una misura e l'altra
        read(1,*) P      	 !valore di nlatt * eta costante nella simulazione
        read(1,*) delta		 !parametro dell'algoritmo
        read(1,*) i_term 	 !passi di termalizzazione
        read(1,*) nret		 !numero di reticoli
        read(1,*) misu		 !flag per che misure fare
        read(1,*) i_shift	 !per raggiungere eta minori
        
C===================================================================================
C misu è una variabile per far scegliere che misure fare per non farle tutte insieme
C e sprecare tempo che può assumere valori interi
C misu = 0: viene calcolato <y^2> e <\Delta y^2>   per calcolare poi energia interna
C misu = 1: viene calcolata la correlazione di y   per calcolare poi E_1 - E_0
C misu = 2: viene calcolata la correlazione di y^2 per calcolare poi E_2 - E_0
C===================================================================================

	write(2,*) misure	 !valori che seriviranno per l'analisi
	write(2,*) nret		 !e per i plot
	
	if (misu /= 0) then	 	!se si vuole calcolare la correlazione
		write(2,*) P		!viene scritto su file anche P e lo shift
		write(2,*) i_shift	!così il codice di analisi lo può leggere, ciò
	endif				!dovuto a come veranno stampati su file i dati
C=====================================================================================
C nel caso delle correlzazioni esse vengono scitte su file ad ogni ciclo di k dove
C k va da 1 fino a i_corr dopo aver fatto la somma dei termini fino a nlatt-k
C quindi per ogni eta ci sono misure curve lunghe i_corr, che verrano mediate elemto
C per elemento per ottenere la curva finale, ciò all'interno del codice analisi_corr.f
C=====================================================================================
				
	do l = 1 + i_shift, nret + i_shift
		nlatt = P*l 
		eta = 0.1		 !eta = P/nlatt = 1/l
		
		allocate(campo(nlatt))		 !alloco il cammino
			
		call init(campo, nlatt)		 !inizializzo il cammino
			
		do k=1, i_term
			call metropolis(campo, nlatt, eta)
		enddo
			
		
		do i=1, misure 	!ciclo sulle misure a fissa temperatura
				
			do j=1, i_dec
				call metropolis(campo, nlatt, eta)   !decorrela la matrice
			enddo
					
			call misurazioni(campo, nlatt)	!misurazione delle osservabili
			!if(l == 1) write(3,*) campo
			!if(l == 15 ) write(4,*) campo	
		enddo
		
C		write(4,*) campo	!salvo un cammino per ogni eta
					!evenetualmente potrebbere essere usato
					!per rifar partire una simulazione da dove
					!si era interotto e aumentare la statistica
					
		deallocate(campo)	!dealloco il cammino per poterlo riallocare
					!con una lunghezza diversa 
	enddo
	call ranfinish
	
	call cpu_time(finish)
	print '("tempo di esecuzione = ", f16.8," secondi.")', finish-start
	
	end program oscillatore
C============================================================================
      	
      	subroutine init(campo, nlatt)
      	
      	real, dimension(nlatt) :: campo
      	
      	do i=1, nlatt
    
      		x = 2.0*ran2() - 1.0       ! rand() random fra -1 e 1
                campo(i) = x
                
      	enddo
      	
      	return
      	end
C============================================================================
	subroutine metropolis(campo, nlatt, eta)
	
	real, dimension(nlatt) :: campo
	common/vario/delta
	
	c1 = 1./eta
      	c2 = (1./eta + eta/2.)
	do i = 1, nlatt				!ciclo su tutti i siti
	
		x = ran2()			!numero casuale per il passo
		y = ran2()			!numero casuale per il test
		ip = i + 1
		if(ip > nlatt) ip = 1		!calcolo dei primi vicini
		im = i - 1
		if(im < 1) im = nlatt 
		F = campo(ip) + campo(im) 	!sommo i primi vicini
     		
     		val = campo(i)
     		val_p = val + delta*((2.0*x - 1.0))
     		dS = c2 * val_p**2 - c1 * val_p * F -		!variazione
     & 		     c2 * val**2 + c1 * val * F      		!dell'azione
     		

     		pb = exp(-dS)		!probabilità di accettare la mossa		

     		if(y < pb) then		!test di accettanza
     			campo(i) = val_p
     		endif
	
	enddo
	
	return
	end
	

C============================================================================

	subroutine misurazioni(campo, nlatt)
	common/flag/misu,i_corr
	
	real, dimension(nlatt) :: campo
	
	if (misu == 0) then
		y2 = 0.0
		Dy2 = 0.0
		
		do i = 1, nlatt
			ip = i + 1
			if(ip > nlatt) ip = 1
			
			y2 = y2 + campo(i)**2
			Dy2 = Dy2 + (campo(i) - campo(ip))**2	
		enddo                             
			
		y2 = y2/float(nlatt)	!media sul singolo cammino di y^2 
		Dy2 = Dy2/float(nlatt)	!media sul singolo cammino di Delta y^2 
		
		write(2,*) y2, Dy2	!salvo su file
	endif
	
	
	if (misu == 1) then
		i_corr = nlatt/2 !fino a dove calcolare la correlazione
		!calcolo della correllazione
		do k = 1, i_corr
		
			yc = 0.0
			do i = 1, nlatt-k
				yc = yc + campo(i)*campo(i+k) 
			enddo
			
			yc = yc/float(nlatt-k)
			
			write(2,*) yc 	
		enddo		
	endif
	
	if (misu == 2) then
		i_corr = nlatt/2 !fino a dove calcolare la correlazione
		!calcolo della correlazione
		do k = 1, i_corr
		
			yc = 0.0
			do i = 1, nlatt-k
				yc = yc + (campo(i)*campo(i+k))**2
			enddo
			
			yc = yc/float(nlatt-k)
			write(2,*) yc 	
		enddo
		
		!calcolo del valore medio di y^2 che servira poi per calcolare
		!la funzione di correlazione connessa, calcolo eseguito nel
		!codice analisi_corr.f
		y2 = 0.0
		do j = 1, nlatt
			y2 = y2 + campo(j)**2
		enddo                             
			
		y2 = y2/float(nlatt)
		write(3,*) y2	
	endif
	
	return
        end
	
c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &          ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &          ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     &          rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c      save iv,iy,idum2
c      data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

c=============================================================================
      subroutine ranstart
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      read(23,*) idum
      read(23,*,end=117) idum2
      do i=1,32
         read(23,*) iv(i)
      enddo
      read(23,*) iy
      close(23)
      goto 118                          !!takes account of the first start
 117  if(idum.ge.0) idum = -idum -1     !!
      close(23)
 118  continue                          !!

      return
      end

c=============================================================================
      subroutine ranfinish
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      write(23,*) idum
      write(23,*) idum2
      do i=1,32
         write(23,*) iv(i)
      enddo
      write(23,*) iy
      close(23)

      return
      end
c=============================================================================