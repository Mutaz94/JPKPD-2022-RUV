Sat Oct 23 22:08:57 CDT 2021
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
$DATA ../../../../data/SD3/A2/dat29.csv ignore=@
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
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -684.912999164603        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8331E+02  5.6067E+01  8.7668E+01  2.8446E+01  1.3649E+02  6.1446E+00 -8.0604E+01 -4.1534E+01 -2.7528E+01 -1.2215E+02
            -2.5006E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1595.43125573342        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0225E+00  9.8000E-01  1.0429E+00  1.0450E+00  9.6899E-01  1.0665E+00  1.1627E+00  8.9046E-01  1.1277E+00  9.9634E-01
             2.4749E+00
 PARAMETER:  1.2223E-01  7.9801E-02  1.4198E-01  1.4399E-01  6.8499E-02  1.6441E-01  2.5072E-01 -1.6022E-02  2.2018E-01  9.6334E-02
             1.0062E+00
 GRADIENT:   1.5295E+02 -1.0044E+01 -8.2915E-01 -1.0703E+01  1.6606E+01  1.4800E+01 -3.6217E-02  5.1723E+00  1.4025E+01 -1.0934E+01
            -8.9374E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1602.26658734135        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0010E+00  8.4272E-01  7.8668E-01  1.1171E+00  7.9649E-01  1.0885E+00  1.2483E+00  1.5972E-01  9.1068E-01  1.0352E+00
             2.4952E+00
 PARAMETER:  1.0100E-01 -7.1121E-02 -1.3993E-01  2.1073E-01 -1.2754E-01  1.8480E-01  3.2177E-01 -1.7344E+00  6.4325E-03  1.3458E-01
             1.0144E+00
 GRADIENT:   9.0245E+01 -1.0471E+01 -1.6992E+01  1.2981E+01  3.9155E+01  2.5353E+01 -1.2695E+01  3.4479E-01 -1.8368E+01  4.3263E+00
            -5.6601E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1604.45528643445        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      253
 NPARAMETR:  9.9285E-01  7.8865E-01  4.8878E-01  1.1193E+00  5.9553E-01  1.0676E+00  1.4694E+00  1.0883E-01  8.6986E-01  8.4976E-01
             2.4696E+00
 PARAMETER:  9.2824E-02 -1.3744E-01 -6.1584E-01  2.1273E-01 -4.1831E-01  1.6542E-01  4.8485E-01 -2.1180E+00 -3.9426E-02 -6.2806E-02
             1.0040E+00
 GRADIENT:   1.8677E+00 -1.4161E+01 -5.3151E+01  3.9575E+01  8.2781E+01  4.3779E+00 -6.1565E+00  2.4416E-01 -1.9305E+01  7.7670E+00
            -5.2276E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1613.95135723841        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  9.9272E-01  4.8794E-01  5.4420E-01  1.2979E+00  4.9163E-01  1.0445E+00  2.0798E+00  1.1724E-01  9.2035E-01  7.6406E-01
             2.6614E+00
 PARAMETER:  9.2698E-02 -6.1757E-01 -5.0843E-01  3.6078E-01 -6.1003E-01  1.4358E-01  8.3227E-01 -2.0435E+00  1.6999E-02 -1.6911E-01
             1.0789E+00
 GRADIENT:   4.8780E+00  2.1825E+01  1.4509E+01  3.2717E+01 -2.5828E+01  1.5030E-01  7.2723E+00  2.6826E-01  5.1350E-01  5.5196E-01
             1.8466E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1618.39252238038        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      604
 NPARAMETR:  9.8008E-01  2.3467E-01  5.5232E-01  1.3906E+00  4.4565E-01  1.0367E+00  2.2524E+00  1.5414E-02  9.1786E-01  8.3459E-01
             2.5826E+00
 PARAMETER:  7.9875E-02 -1.3496E+00 -4.9362E-01  4.2970E-01 -7.0822E-01  1.3600E-01  9.1200E-01 -4.0725E+00  1.4295E-02 -8.0815E-02
             1.0488E+00
 GRADIENT:  -3.9054E+00  4.5992E+00  1.6512E+01 -1.0934E+00 -2.5249E+01 -1.7305E-01 -2.5557E+00  4.9757E-03 -1.2120E+00 -1.9974E+00
            -5.5201E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1622.09872271464        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      779
 NPARAMETR:  9.7150E-01  3.3558E-02  6.4549E-01  1.5232E+00  4.7409E-01  1.0305E+00  1.7790E+00  1.0000E-02  8.7655E-01  8.9659E-01
             2.5909E+00
 PARAMETER:  7.1083E-02 -3.2945E+00 -3.3775E-01  5.2081E-01 -6.4637E-01  1.3007E-01  6.7604E-01 -1.0661E+01 -3.1763E-02 -9.1554E-03
             1.0520E+00
 GRADIENT:  -4.0637E+00  8.9355E-01  2.7226E+00  8.0102E+00 -4.4534E+00 -7.3717E-01 -4.3951E-02  0.0000E+00 -3.1741E+00 -2.0977E-01
            -4.9843E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1622.65002648659        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      954
 NPARAMETR:  9.7241E-01  1.0000E-02  6.3341E-01  1.5250E+00  4.6657E-01  1.0314E+00  1.7166E+00  1.0000E-02  8.8150E-01  8.9618E-01
             2.5896E+00
 PARAMETER:  7.2024E-02 -4.6084E+00 -3.5664E-01  5.2202E-01 -6.6235E-01  1.3094E-01  6.4033E-01 -1.5317E+01 -2.6131E-02 -9.6160E-03
             1.0515E+00
 GRADIENT:   4.3791E-01  0.0000E+00 -1.9950E+00 -1.1733E+00  3.4259E+00 -7.7509E-02 -3.1268E-03  0.0000E+00  1.1246E-01  2.4874E-01
             1.1019E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1622.66100815845        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1011
 NPARAMETR:  9.7221E-01  1.0000E-02  6.2677E-01  1.5233E+00  4.6205E-01  1.0316E+00  1.7415E+00  1.0000E-02  8.8199E-01  8.9174E-01
             2.5869E+00
 PARAMETER:  7.1817E-02 -4.6249E+00 -3.6718E-01  5.2086E-01 -6.7209E-01  1.3114E-01  6.5477E-01 -1.5358E+01 -2.5571E-02 -1.4578E-02
             1.0504E+00
 GRADIENT:   3.2530E-03  0.0000E+00 -9.3152E-02  4.5076E-01  7.6182E-02 -2.6863E-02 -3.3188E-03  0.0000E+00 -9.1198E-02  1.1459E-02
             1.0684E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1011
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1400E-03 -5.3871E-04  1.8001E-05 -9.8502E-03 -1.3382E-02
 SE:             2.9244E-02  2.8336E-04  1.6591E-04  2.7886E-02  2.2197E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6890E-01  5.7284E-02  9.1360E-01  7.2392E-01  5.4660E-01

 ETASHRINKSD(%)  2.0292E+00  9.9051E+01  9.9444E+01  6.5782E+00  2.5639E+01
 ETASHRINKVR(%)  4.0173E+00  9.9991E+01  9.9997E+01  1.2724E+01  4.4704E+01
 EBVSHRINKSD(%)  1.9448E+00  9.9203E+01  9.9367E+01  6.1556E+00  2.5348E+01
 EBVSHRINKVR(%)  3.8517E+00  9.9994E+01  9.9996E+01  1.1932E+01  4.4271E+01
 RELATIVEINF(%)  8.6681E+01  4.9581E-04  2.6404E-04  1.0074E+01  3.0115E+00
 EPSSHRINKSD(%)  2.6332E+01
 EPSSHRINKVR(%)  4.5731E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1622.6610081584479     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -703.72247495377519     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1622.661       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.72E-01  1.00E-02  6.27E-01  1.52E+00  4.62E-01  1.03E+00  1.74E+00  1.00E-02  8.82E-01  8.92E-01  2.59E+00
 


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
 #CPUT: Total CPU Time in Seconds,       83.615
Stop Time:
Sat Oct 23 22:09:10 CDT 2021
