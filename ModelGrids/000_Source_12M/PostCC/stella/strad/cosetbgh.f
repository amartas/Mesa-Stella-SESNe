      SUBROUTINECOSET(METH,NQ,EL,TQ,MAXORD,IDOUB)
      IMPLICITREAL*8(A-H,O-Z)
      COMMON/HNT/HNT(7)
      DIMENSIONPC(6),XX(6)
      DIMENSIONPERST(12,2,3),EL(13),TQ(4)
      DATAPERST/1.D0,1.D0,2.D0,1.D0,.3158D0,.07407D0,.01391D0,2.182D-3,2.945D-4,3.492D-5,3.692D-6,3.524D-7,1.D0,1.D0,.5D0,.16666667D
     *0,.041666667D0,7*1.D0,2.D0,12.D0,24.D0,37.89D0,53.33D0,70.08D0,87.97D0,106.9D0,126.7D0,147.4D0,168.8D0,191.0D0,2.D0,4.5D0,7.33
     *33333D0,10.416667D0,13.7D0,7*1.D0,12.0D0,24.0D0,37.89D0,53.33D0,70.08D0,87.97D0,106.9D0,126.7D0,147.4D0,168.8D0,191.0D0,1.D0,3
     *.D0,6.D0,9.1666667D0,12.5D0,8*1.D0/
      GOTO(09999,09998,09997),METH
      GOTO09996
09999 CONTINUE
      GOTO(09995,09994,09993,09992,09991,09990,09989,09988,09987,09986,09985,09984),NQ
      GOTO09983
09995 CONTINUE
      EL(1)=1.D0
      GOTO09983
09994 CONTINUE
      EL(1)=5.D-1
      EL(3)=5.D-1
      GOTO09983
09993 CONTINUE
      EL(1)=4.166666666666667D-01
      EL(3)=7.5D-1
      EL(4)=1.666666666666667D-01
      GOTO09983
09992 CONTINUE
      EL(1)=3.75D-1
      EL(3)=9.166666666666667D-01
      EL(4)=3.333333333333333D-01
      EL(5)=4.166666666666667D-02
      GOTO09983
09991 CONTINUE
      EL(1)=3.486111111111111D-01
      EL(3)=1.041666666666667D0
      EL(4)=4.861111111111111D-01
      EL(5)=1.041666666666667D-01
      EL(6)=8.333333333333333D-03
      GOTO09983
09990 CONTINUE
      EL(1)=3.29861111111111D-01
      EL(3)=1.141666666666667D0
      EL(4)=6.25D-1
      EL(5)=1.770833333333333D-01
      EL(6)=2.5D-2
      EL(7)=1.388888888888889D-03
      GOTO09983
09989 CONTINUE
      EL(1)=3.155919312169312D-01
      EL(3)=1.225D0
      EL(4)=7.518518518518518D-01
      EL(5)=2.552083333333333D-01
      EL(6)=4.861111111111111D-02
      EL(7)=4.861111111111111D-03
      EL(8)=1.984126984126984D-04
      GOTO09983
09988 CONTINUE
      EL(1)=3.042245370370370D-01
      EL(3)=1.296428571428571D0
      EL(4)=8.685185185185185D-01
      EL(5)=3.357638888888889D-01
      EL(6)=7.777777777777778D-02
      EL(7)=1.064814814814815D-02
      EL(8)=7.936507936507937D-04
      EL(9)=2.48015873015873D-05
      GOTO09983
09987 CONTINUE
      EL(1)=2.948680004409171D-01
      EL(3)=1.358928571428571D0
      EL(4)=9.765542328042328D-01
      EL(5)=4.171875D-1
      EL(6)=1.113541666666667D-01
      EL(7)=1.875D-2
      EL(8)=1.934523809523810D-03
      EL(9)=1.116071428571429D-04
      EL(10)=2.755731922398589D-06
      GOTO09983
09986 CONTINUE
      EL(1)=2.869754464285714D-01
      EL(3)=1.414484126984127D0
      EL(4)=1.077215608465608D0
      EL(5)=4.985670194003527D-01
      EL(6)=0.1484375D0
      EL(7)=2.906057098765432D-02
      EL(8)=3.720238095238095D-03
      EL(9)=2.996858465608466D-04
      EL(10)=1.377865961199295D-05
      EL(11)=2.755731922398589D-07
      GOTO09983
09985 CONTINUE
      EL(1)=2.801895964439367D-01
      EL(3)=1.464484126984127D0
      EL(4)=1.171514550264550D0
      EL(5)=5.793581900352734D-01
      EL(6)=1.883228615520282D-01
      EL(7)=4.143036265432099D-02
      EL(8)=6.211144179894180D-03
      EL(9)=6.252066798941799D-04
      EL(10)=4.041740152851264D-05
      EL(11)=1.515652557319224D-06
      EL(12)=2.505210838544172D-08
      GOTO09983
09984 CONTINUE
      EL(1)=2.742655400315991D-01
      EL(3)=1.509938672438672D0
      EL(4)=1.260271164021164D0
      EL(5)=6.592341820987654D-01
      EL(6)=2.304580026455026D-01
      EL(7)=5.569724610523222D-02
      EL(8)=9.439484126984127D-03
      EL(9)=1.119274966931217D-03
      EL(10)=9.093915343915344D-05
      EL(11)=4.822530864197531D-06
      EL(12)=1.503126503126503D-07
      EL(13)=2.087675698786810D-09
09983 CONTINUE
      DO09982K=1,3
      TQ(K)=PERST(NQ,1,K)
09982 CONTINUE
      GOTO09996
09998 CONTINUE
      GOTO(09979,09978,09977,09976,09975),NQ
      GOTO09974
09979 CONTINUE
      EL(1)=1.D0
      GOTO09974
09978 CONTINUE
      EL(1)=6.666666666666667D-01
      EL(3)=3.333333333333333D-01
      GOTO09974
09977 CONTINUE
      EL(1)=5.454545454545455D-01
      EL(3)=EL(1)
      EL(4)=9.090909090909091D-02
      GOTO09974
09976 CONTINUE
      EL(1)=0.48D0
      EL(3)=0.7D0
      EL(4)=0.2D0
      EL(5)=2.D-2
      GOTO09974
09975 CONTINUE
      EL(1)=4.379562043795620D-01
      EL(3)=8.211678832116788D-01
      EL(4)=3.102189781021898D-01
      EL(5)=5.474452554744526D-02
      EL(6)=3.649635036496350D-03
09974 CONTINUE
      DO09973K=1,3
      TQ(K)=PERST(NQ,2,K)
09973 CONTINUE
      GOTO09996
09997 CONTINUE
      IF(NQ.EQ.1)THEN
      EL(1)=1.D0
      DO09970K=1,2
      TQ(K)=PERST(1,1,K)
09970 CONTINUE
      ELSE
      PC(1)=1.0D0
      XX(1)=1.0D0
      DO09967J=2,NQ
      XX(J)=HNT(J)/HNT(1)+1.D0
09967 CONTINUE
      DO09964INQ=1,NQ
      PC(INQ+1)=0.0D0
      DO09961IB=1,INQ
      I=INQ+2-IB
      PC(I)=PC(I-1)+XX(INQ)*PC(I)
09961 CONTINUE
      PC(1)=XX(INQ)*PC(1)
09964 CONTINUE
      DO09958I=1,NQ+1
      EL(I)=PC(I)/PC(2)
09958 CONTINUE
      TQ(1)=XX(NQ)/PC(1)
      TQ(2)=(XX(NQ)+1.D0)/EL(1)
      ENDIF
      IF(IDOUB.EQ.1.AND.NQ.NE.MAXORD)THEN
      PC(1)=1.D0
      PC(2)=1.D0
      DO09955I=2,NQ
      PC(1)=PC(1)*(HNT(I)+HNT(1))
      PC(2)=PC(2)*HNT(I+1)
09955 CONTINUE
      PC(3)=HNT(1)*PC(1)*(HNT(NQ+1)+HNT(1))/(HNT(2)*PC(2)*HNT(NQ+2))
      TQ(3)=DBLE(NQ+2)/(PC(3)*EL(1))
      ENDIF
09996 CONTINUE
      EL(2)=1.D0
      TQ(4)=.5D0/DBLE(NQ+2)
      RETURN
      END
      SUBROUTINERESCAL(Y,NYDIM,RH,RMAX,RC,LNQ)
      IMPLICITREAL*8(A-H,O-Z)
      COMMON/STCOM1/T,H,HMIN,HMAX,EPS,N,METH,KFLAG,JSTART
      DIMENSIONY(NYDIM,*)
      RH=DMAX1(RH,HMIN/DABS(H))
      RH=DMIN1(RH,HMAX/DABS(H),RMAX)
      R1=1.D0
      DO09952J=2,LNQ
      R1=R1*RH
      DO09949I=1,N
      Y(I,J)=Y(I,J)*R1
09949 CONTINUE
09952 CONTINUE
      H=H*RH
      RC=RC*RH
      RETURN
      END