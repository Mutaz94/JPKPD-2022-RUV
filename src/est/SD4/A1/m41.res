Sun Oct 24 01:58:49 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/SD4/A1/dat41.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2

$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr
$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       24 OCT 2021
Days until program expires : 175
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m41.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1497.31585416023        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5919E+02 -1.2165E+01  9.7743E+00 -4.0359E+00  7.5009E+01  5.0024E+01 -6.7609E+00 -1.4275E+00  2.1356E+00 -3.3301E+01
            -3.2994E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1562.22852959100        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0154E+00  1.0192E+00  9.4423E-01  1.0414E+00  9.1061E-01  9.2875E-01  9.9875E-01  9.5241E-01  9.2860E-01  1.0155E+00
             1.9459E+00
 PARAMETER:  1.1529E-01  1.1902E-01  4.2613E-02  1.4059E-01  6.3646E-03  2.6082E-02  9.8748E-02  5.1239E-02  2.5917E-02  1.1533E-01
             7.6572E-01
 GRADIENT:   6.3340E+01  1.6441E+01  3.6760E+00  1.8123E+01 -2.4499E+01  9.9428E-02  1.0858E+00  6.8027E+00  3.2601E+00  1.3577E+01
             8.1481E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1567.31853019230        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      204
 NPARAMETR:  1.0044E+00  8.2802E-01  6.0580E-01  1.1446E+00  6.5458E-01  9.4551E-01  1.2539E+00  4.0776E-01  8.4369E-01  6.4226E-01
             1.8699E+00
 PARAMETER:  1.0441E-01 -8.8716E-02 -4.0121E-01  2.3501E-01 -3.2376E-01  4.3965E-02  3.2626E-01 -7.9708E-01 -6.9967E-02 -3.4276E-01
             7.2588E-01
 GRADIENT:  -1.0086E+02  1.6122E+01 -1.8808E+01  5.7077E+01  1.7030E+01 -7.9078E+00 -3.1974E+00  1.9472E+00  4.4083E+00 -1.5000E+00
             6.5906E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1579.97480065393        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      383
 NPARAMETR:  1.0450E+00  6.2346E-01  6.0133E-01  1.2353E+00  5.8466E-01  9.5803E-01  1.6182E+00  1.3764E-01  7.8630E-01  7.6822E-01
             1.5058E+00
 PARAMETER:  1.4402E-01 -3.7248E-01 -4.0861E-01  3.1129E-01 -4.3673E-01  5.7126E-02  5.8132E-01 -1.8831E+00 -1.4042E-01 -1.6367E-01
             5.0935E-01
 GRADIENT:   7.1038E+00  1.0818E+01 -8.8986E+00  2.5944E+01  1.2200E+01 -6.4468E-01  5.4649E-01  3.4008E-01 -8.7542E-01  2.6001E-01
            -6.2513E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1583.27307248324        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      558
 NPARAMETR:  1.0309E+00  2.8565E-01  5.1143E-01  1.3653E+00  4.3602E-01  9.5267E-01  2.3510E+00  1.0000E-02  7.7097E-01  8.0370E-01
             1.4800E+00
 PARAMETER:  1.3047E-01 -1.1530E+00 -5.7055E-01  4.1138E-01 -7.3007E-01  5.1511E-02  9.5483E-01 -6.0873E+00 -1.6010E-01 -1.1853E-01
             4.9202E-01
 GRADIENT:  -1.6788E+01  6.8920E+00  1.3086E+01  9.3334E+00 -2.5671E+01 -2.0768E+00 -9.5003E-02  0.0000E+00  3.2361E+00  2.6918E+00
             3.4757E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1584.00716875783        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      733
 NPARAMETR:  1.0362E+00  2.0662E-01  4.8783E-01  1.3745E+00  4.1204E-01  9.5788E-01  2.8009E+00  1.0000E-02  7.5171E-01  7.9762E-01
             1.4566E+00
 PARAMETER:  1.3554E-01 -1.4769E+00 -6.1778E-01  4.1809E-01 -7.8663E-01  5.6964E-02  1.1300E+00 -8.0484E+00 -1.8541E-01 -1.2613E-01
             4.7612E-01
 GRADIENT:   1.7197E+00  2.2993E-01  6.9305E+00 -2.0224E+01 -1.0554E+01  7.5017E-01 -1.2982E+00  0.0000E+00 -1.3634E+00  8.7865E-01
            -8.1770E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1584.36481169806        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  1.0322E+00  1.5472E-01  5.2239E-01  1.4176E+00  4.2895E-01  9.5434E-01  3.5080E+00  1.0000E-02  7.4123E-01  8.1231E-01
             1.4650E+00
 PARAMETER:  1.3168E-01 -1.7661E+00 -5.4933E-01  4.4894E-01 -7.4640E-01  5.3269E-02  1.3550E+00 -9.5655E+00 -1.9945E-01 -1.0788E-01
             4.8185E-01
 GRADIENT:  -4.9762E-01 -2.5962E-01 -4.8212E-01 -5.2601E+00 -6.0581E-03  4.5300E-02  3.8672E-01  0.0000E+00 -2.7387E-02 -4.0309E-01
             6.5324E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1584.37383887623        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1090
 NPARAMETR:  1.0326E+00  1.6045E-01  5.3456E-01  1.4219E+00  4.3713E-01  9.5421E-01  3.4490E+00  1.0000E-02  7.4131E-01  8.1929E-01
             1.4656E+00
 PARAMETER:  1.3212E-01 -1.7298E+00 -5.2631E-01  4.5202E-01 -7.2752E-01  5.3131E-02  1.3381E+00 -9.3054E+00 -1.9934E-01 -9.9311E-02
             4.8225E-01
 GRADIENT:   5.0328E-01 -1.1816E-01  1.0319E-01 -2.0731E+00 -4.2011E-01  6.3086E-03  1.5301E-01  0.0000E+00 -3.0047E-02 -1.6141E-02
             3.2264E-03

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1584.37383887623        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1112
 NPARAMETR:  1.0326E+00  1.6045E-01  5.3456E-01  1.4219E+00  4.3713E-01  9.5421E-01  3.4490E+00  1.0000E-02  7.4131E-01  8.1929E-01
             1.4656E+00
 PARAMETER:  1.3212E-01 -1.7298E+00 -5.2631E-01  4.5202E-01 -7.2752E-01  5.3131E-02  1.3381E+00 -9.3054E+00 -1.9934E-01 -9.9311E-02
             4.8225E-01
 GRADIENT:   5.0328E-01 -1.1816E-01  1.0319E-01 -2.0731E+00 -4.2011E-01  6.3086E-03  1.5301E-01  0.0000E+00 -3.0047E-02 -1.6141E-02
             3.2264E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1112
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0977E-03  3.4013E-02 -2.9867E-04 -2.1143E-02  7.5379E-03
 SE:             2.9644E-02  1.5150E-02  2.1785E-04  2.6708E-02  2.3824E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7046E-01  2.4762E-02  1.7037E-01  4.2858E-01  7.5170E-01

 ETASHRINKSD(%)  6.8963E-01  4.9246E+01  9.9270E+01  1.0523E+01  2.0185E+01
 ETASHRINKVR(%)  1.3745E+00  7.4240E+01  9.9995E+01  1.9939E+01  3.6296E+01
 EBVSHRINKSD(%)  9.6159E-01  6.2063E+01  9.9181E+01  8.1543E+00  1.4974E+01
 EBVSHRINKVR(%)  1.9139E+00  8.5608E+01  9.9993E+01  1.5644E+01  2.7706E+01
 RELATIVEINF(%)  9.6780E+01  5.0016E+00  3.4478E-04  2.8092E+01  3.7407E+00
 EPSSHRINKSD(%)  4.1641E+01
 EPSSHRINKVR(%)  6.5942E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1584.3738388762308     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -849.22301231249264     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1584.374       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.60E-01  5.35E-01  1.42E+00  4.37E-01  9.54E-01  3.45E+00  1.00E-02  7.41E-01  8.19E-01  1.47E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       32.996
Stop Time:
Sun Oct 24 01:58:57 CDT 2021
