!$$$$$$$$$$$$$$$$$$  parametri elettrici  $$$$$$$$$$$$$$$$$$$$$$$$
 
MODULE mod_string
  IMPLICIT  NONE
  SAVE

CONTAINS
  
FUNCTION STRING(inn)
  IMPLICIT NONE
  !  converte l'INTEGER "inn" in una stringa ascii di POS caratteri 
  INTEGER, PARAMETER :: POS= 4
  INTEGER, INTENT(IN) :: inn
  CHARACTER(LEN=POS) :: STRING
  !............................................................Tipi locali
  INTEGER :: cifra, np, mm, num
  
  IF (inn > (10**POS)-1) stop "ERRORE: (inn > (10**3)-1)  in STRING"
  num= inn
  DO np= 1, POS
     mm= pos-np
     cifra= num/(10**mm)            ! divisione fra interi
     STRING(np:np)= ACHAR(48+cifra)
     num= num - cifra*(10**mm)
  END DO
END FUNCTION STRING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION DBSTRING(inn)
  IMPLICIT NONE
  !  converte l'INTEGER "inn" in una stringa ascii di POS caratteri 
  INTEGER, PARAMETER :: POS1= 2 !cifre prima della virgola
  INTEGER, PARAMETER :: POS2= 4 !cifre dopo la virgola
  REAL(8), INTENT(IN) :: inn
  CHARACTER(LEN=(POS1+POS2+1+4)) :: DBSTRING
  !!+1 underscore +4 esponente

  !............................................................Tipi locali
  INTEGER :: cifra, np, mm, num,nn,cont,cont1
  REAL(8) :: dec,app

  IF(inn.ge.0.0d0) THEN
     cont=0 !conta quanti zeri ho prima della virgola
     cont1=0 !conta quanti zeri ho dopo la seconda cifra (>100)
     app=inn
     
     DO WHILE(app.gt.1.0d9)
        app=app/10.0d0
        cont1=cont1+1
     END DO

     num=INT(app) !parte prima della virgola
     
     IF(inn.ne.0.0d0) THEN
        IF(inn.lt.(10.0d0**POS1)) THEN
           DO WHILE(num==0)
              app=app*10.0d0
              num=INT(app)
              cont=cont+1
           END DO
        ELSE
           DO WHILE(num.ge.(10.0d0**POS1))
              app=app/10.0d0
              num=INT(app)
              cont1=cont1+1
           END DO
        END IF

        DO np= 1, POS1
           mm= pos1-np
           cifra= num/(10**mm)            ! divisione fra interi
           DBSTRING(np:np)= ACHAR(48+cifra)
           num= num - cifra*(10**mm)
        END DO
        DBSTRING(POS1+1:POS1+1)='.' 
        
        dec=app-INT(app)
        DO np=POS1+2,POS1+POS2+1
           app=dec*10.0d0
           num=INT(app)
           DBSTRING(np:np)=ACHAR(48+num)
           dec=app-INT(app)
        END DO
        
        IF ((cont==0).and.(cont1==0)) THEN
           DBSTRING(POS1+POS2+2:POS1+POS2+1+4)='E+00'
        ELSE 
           IF(cont.gt.0) THEN
              DBSTRING(POS1+POS2+2:POS1+POS2+1+2)='E-'
!!!trasformo in intero il contatore     
              num= cont
              DO np= 1, 2
                 mm= 2-np
                 cifra= num/(10**mm)            ! divisione fra interi
                 DBSTRING(POS1+POS2+3+np:POS1+POS2+3+np)= ACHAR(48+cifra)
                 num= num - cifra*(10**mm)
              END DO
           END IF
           IF(cont1.gt.0) THEN
              DBSTRING(POS1+POS2+2:POS1+POS2+1+2)='E+'
!!!trasformo in intero il contatore     
              num= cont1
              DO np= 1, 2
                 mm= 2-np
                 cifra= num/(10**mm)            ! divisione fra interi
                 DBSTRING(POS1+POS2+3+np:POS1+POS2+3+np)= ACHAR(48+cifra)
                 num= num - cifra*(10**mm)
              END DO
           END IF
        END IF
     ELSE
        DBSTRING(1:POS1+POS2+1+4)='00.0000E+00'
     END IF !fine numero >=0
  ELSE
     DBSTRING(1:1)='-'
     cont=0 !conta quanti zeri ho prima della virgola
     cont1=0 !conta quanti zeri ho dopo la seconda cifra (>100)
     app=ABS(inn)
     
     DO WHILE(app.gt.1.0d9)
        app=app/10.0d0
        cont1=cont1+1
     END DO
     
     num=INT(app) !parte prima della virgola
     IF(ABS(inn).lt.(10.0d0**(POS1-1))) THEN
        DO WHILE(num==0)
           app=app*10.0d0
           num=INT(app)
           cont=cont+1
        END DO
     ELSE
        DO WHILE(num.ge.(10.0d0**(POS1-1)))
           app=app/10.0d0
           num=INT(app)
           cont1=cont1+1
        END DO
     END IF
     
     DO np= 2, POS1
        mm= pos1-np
        cifra= num/(10**mm)            ! divisione fra interi
        DBSTRING(np:np)= ACHAR(48+cifra)
        num= num - cifra*(10**mm)
     END DO
     DBSTRING(POS1+1:POS1+1)='.' 

     dec=app-INT(app)
     DO np=POS1+2,POS1+POS2+1
        app=dec*10.0d0
        num=INT(app)
        DBSTRING(np:np)=ACHAR(48+num)
        dec=app-INT(app)
     END DO
     
     IF ((cont==0).and.(cont1==0)) THEN
        DBSTRING(POS1+POS2+2:POS1+POS2+1+4)='E+00'
     ELSE 
        IF(cont.gt.0) THEN
           DBSTRING(POS1+POS2+2:POS1+POS2+1+2)='E-'
!!!trasformo in intero il contatore     
           num= cont
           DO np= 1, 2
              mm= 2-np
              cifra= num/(10**mm)            ! divisione fra interi
              DBSTRING(POS1+POS2+3+np:POS1+POS2+3+np)= ACHAR(48+cifra)
              num= num - cifra*(10**mm)
           END DO
        END IF
        IF(cont1.gt.0) THEN
           DBSTRING(POS1+POS2+2:POS1+POS2+1+2)='E+'
!!!trasformo in intero il contatore     
           num= cont1
           DO np= 1, 2
              mm= 2-np
              cifra= num/(10**mm)            ! divisione fra interi
              DBSTRING(POS1+POS2+3+np:POS1+POS2+3+np)= ACHAR(48+cifra)
              num= num - cifra*(10**mm)
           END DO
        END IF
     END IF
  END IF
  
END FUNCTION DBSTRING


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE mod_string
