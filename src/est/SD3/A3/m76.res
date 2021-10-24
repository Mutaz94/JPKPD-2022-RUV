Sat Oct 23 22:45:26 CDT 2021
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
$DATA ../../../../data/SD3/A3/dat76.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m76.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   180.225811526202        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1008E+02  1.0277E+02  1.6229E+02  5.3876E+01  1.5588E+02 -4.1072E+00 -4.3461E+01 -1.1170E+02 -9.5670E+00 -1.5216E+02
            -4.1243E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1469.77663242108        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.7293E-01  9.9385E-01  9.4747E-01  1.0802E+00  1.0363E+00  9.2307E-01  8.8589E-01  9.5810E-01  7.9298E-01  9.2314E-01
             4.4211E+00
 PARAMETER:  7.2552E-02  9.3833E-02  4.6035E-02  1.7718E-01  1.3570E-01  1.9948E-02 -2.1161E-02  5.7199E-02 -1.3196E-01  2.0030E-02
             1.5864E+00
 GRADIENT:  -4.5439E+01 -1.6413E+01 -2.8697E+01  5.3940E+00  1.5950E+01 -2.2656E+01  4.3376E+00  6.2603E+00  4.5913E+00  1.8106E+01
             2.2893E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1494.87704409904        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.6875E-01  8.3544E-01  5.2713E-01  1.1242E+00  6.0959E-01  9.9121E-01  6.8460E-01  1.2860E-01  9.8567E-01  4.3931E-01
             3.9540E+00
 PARAMETER:  6.8249E-02 -7.9799E-02 -5.4031E-01  2.1711E-01 -3.9497E-01  9.1172E-02 -2.7892E-01 -1.9511E+00  8.5568E-02 -7.2255E-01
             1.4747E+00
 GRADIENT:  -3.8744E+01  2.8126E+01  1.7246E+01  3.7589E+01 -3.3331E+01 -5.9659E+00 -9.0496E-01  1.5936E-01  1.5293E+01  3.3029E+00
             1.5265E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1505.90019455225        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.7238E-01  7.2298E-01  2.7672E-01  1.1297E+00  3.9780E-01  1.0178E+00  1.0412E+00  2.8862E-02  9.1663E-01  3.4791E-01
             3.3919E+00
 PARAMETER:  7.1993E-02 -2.2438E-01 -1.1847E+00  2.2192E-01 -8.2180E-01  1.1760E-01  1.4033E-01 -3.4452E+00  1.2949E-02 -9.5582E-01
             1.3214E+00
 GRADIENT:  -1.8653E+01  5.4642E+01 -4.2198E+00  1.2268E+02 -4.1541E+00 -7.7520E+00  4.1770E+00 -9.0353E-03 -1.7196E+01 -2.8760E+00
             4.9129E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1509.95019353163        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      343
 NPARAMETR:  9.7569E-01  7.5949E-01  2.8596E-01  1.0910E+00  4.1435E-01  1.0251E+00  9.3000E-01  2.8400E-02  9.5245E-01  4.0892E-01
             3.3032E+00
 PARAMETER:  7.5387E-02 -1.7510E-01 -1.1519E+00  1.8711E-01 -7.8104E-01  1.2476E-01  2.7426E-02 -3.4614E+00  5.1281E-02 -7.9425E-01
             1.2949E+00
 GRADIENT:  -3.5642E+01  4.9748E+01  5.7016E+00  6.0123E+01 -3.7131E+01 -7.5882E+00 -1.4595E-01 -8.6097E-03 -1.0402E+01 -3.8177E+00
             1.6434E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1519.76829287376        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  1.0029E+00  4.1875E-01  1.6386E-01  1.0221E+00  2.3579E-01  1.0970E+00  3.4308E-01  1.0000E-02  1.1875E+00  7.5353E-01
             2.8356E+00
 PARAMETER:  1.0292E-01 -7.7048E-01 -1.7088E+00  1.2184E-01 -1.3448E+00  1.9257E-01 -9.6978E-01 -5.7633E+00  2.7183E-01 -1.8299E-01
             1.1422E+00
 GRADIENT:   3.1292E+01  3.2483E+01  2.0836E+01 -4.4752E+00 -5.6047E+01  9.9368E+00  4.4878E-01  0.0000E+00 -3.4341E-01  6.2539E+00
            -4.0178E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1522.29327074545        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      693
 NPARAMETR:  9.8888E-01  3.9585E-01  1.6726E-01  1.0330E+00  2.3883E-01  1.0609E+00  2.4122E-01  1.0000E-02  1.1746E+00  7.1529E-01
             2.9683E+00
 PARAMETER:  8.8814E-02 -8.2673E-01 -1.6882E+00  1.3247E-01 -1.3320E+00  1.5914E-01 -1.3220E+00 -5.6861E+00  2.6092E-01 -2.3506E-01
             1.1880E+00
 GRADIENT:   3.4770E+00 -2.0145E+00 -2.6768E+00 -3.9959E+00  6.6688E+00  3.9602E-01  1.6227E-01  0.0000E+00  1.5847E+00 -2.9688E-01
             1.5789E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1522.35888030694        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      870
 NPARAMETR:  9.8721E-01  3.9516E-01  1.6829E-01  1.0372E+00  2.3896E-01  1.0595E+00  1.4085E-01  1.0000E-02  1.1663E+00  7.1914E-01
             2.9668E+00
 PARAMETER:  8.7129E-02 -8.2846E-01 -1.6820E+00  1.3654E-01 -1.3315E+00  1.5783E-01 -1.8600E+00 -5.6712E+00  2.5386E-01 -2.2971E-01
             1.1875E+00
 GRADIENT:  -2.0784E-02 -2.3113E-01  1.9979E-01  1.0744E-01  2.3622E-01  1.3524E-02  3.7500E-02  0.0000E+00 -1.2052E-02 -1.6128E-02
             2.0519E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1522.37985683972        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1045
 NPARAMETR:  9.8722E-01  3.9399E-01  1.6731E-01  1.0357E+00  2.3802E-01  1.0595E+00  1.0042E-02  1.0000E-02  1.1690E+00  7.2144E-01
             2.9651E+00
 PARAMETER:  8.7135E-02 -8.3142E-01 -1.6879E+00  1.3506E-01 -1.3354E+00  1.5780E-01 -4.5010E+00 -5.6712E+00  2.5616E-01 -2.2651E-01
             1.1869E+00
 GRADIENT:   3.0111E-02  5.7736E-02  2.6894E-01  1.3816E-01 -2.2869E-01 -2.8669E-02  3.0733E-04  0.0000E+00 -4.1834E-02  4.7882E-02
            -5.6515E-02

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1522.38129084027        NO. OF FUNC. EVALS.:  67
 CUMULATIVE NO. OF FUNC. EVALS.:     1112
 NPARAMETR:  9.8763E-01  3.9480E-01  1.6675E-01  1.0352E+00  2.3693E-01  1.0592E+00  1.0000E-02  1.0000E-02  1.1700E+00  7.2281E-01
             2.9626E+00
 PARAMETER:  8.7064E-02 -8.3484E-01 -1.6912E+00  1.3421E-01 -1.3380E+00  1.5798E-01 -5.8199E+00 -5.6712E+00  2.5733E-01 -2.2688E-01
             1.1868E+00
 GRADIENT:  -1.4464E-01 -4.1511E-01  1.1972E-02 -7.7383E-02  9.0337E-01  2.4299E-02  0.0000E+00  0.0000E+00  1.4587E-02 -6.4322E-02
             1.0287E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1112
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2396E-03 -9.9784E-05  1.7776E-04 -1.2081E-02  9.1526E-04
 SE:             2.9022E-02  1.4071E-04  2.1394E-04  2.6111E-02  2.4388E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6593E-01  4.7825E-01  4.0604E-01  6.4361E-01  9.7006E-01

 ETASHRINKSD(%)  2.7734E+00  9.9529E+01  9.9283E+01  1.2524E+01  1.8297E+01
 ETASHRINKVR(%)  5.4699E+00  9.9998E+01  9.9995E+01  2.3480E+01  3.3246E+01
 EBVSHRINKSD(%)  2.4551E+00  9.9529E+01  9.9360E+01  1.0581E+01  1.8572E+01
 EBVSHRINKVR(%)  4.8500E+00  9.9998E+01  9.9996E+01  2.0042E+01  3.3694E+01
 RELATIVEINF(%)  9.5005E+01  3.2700E-04  3.5416E-04  4.0820E+01  2.7965E+00
 EPSSHRINKSD(%)  2.5944E+01
 EPSSHRINKVR(%)  4.5157E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1522.3812908402740     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -603.44275763560131     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1522.381       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.87E-01  3.93E-01  1.67E-01  1.03E+00  2.37E-01  1.06E+00  1.00E-02  1.00E-02  1.17E+00  7.21E-01  2.96E+00
 


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
 #CPUT: Total CPU Time in Seconds,       87.376
Stop Time:
Sat Oct 23 22:45:40 CDT 2021
