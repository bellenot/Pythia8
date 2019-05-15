C...Extract parton-level events according to Les Houches Accord.
C...The files produced here can then be used as input 
C...for hadron-level event simulation, also in Pythia8.

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP

C...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA

C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

C...Number of events.
      NEV=100

C...Event kind.
      MSEL=6

C...Files for output.
      MSTP(161)=21
      OPEN(21,FILE='ttsample.init',STATUS='unknown')
      MSTP(162)=22
      OPEN(22,FILE='ttsample.evnt',STATUS='unknown')

C...Initialize.
      CALL PYINIT('CMS','P','PBAR',1960D0)

C...Event loop. List first few events.
      DO 200 IEV=1,NEV
        CALL PYUPEV
        IF(IEV.LE.2) THEN 
          CALL PYLIST(2)
          CALL PYLIST(7)
        ENDIF 
 200  CONTINUE

C...Final statistics.
      CALL PYSTAT(1)
      CALL PYUPIN

      END
