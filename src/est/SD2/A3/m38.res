Sat Oct 23 17:59:38 CDT 2021
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
$DATA ../../../../data/SD2/A3/dat38.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      800
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      700
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   764.289369068930        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6910E+02  1.0553E+02  2.4414E+02 -4.7658E+00  1.5465E+02  3.5637E+00 -1.0403E+02 -1.6786E+02 -5.3510E+01 -1.2837E+02
            -6.8372E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1800.86881410116        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.7519E-01  1.0485E+00  8.9113E-01  1.1255E+00  9.8540E-01  9.0898E-01  1.0298E+00  9.5005E-01  9.9767E-01  8.9410E-01
             5.2876E+00
 PARAMETER:  7.4879E-02  1.4739E-01 -1.5260E-02  2.1822E-01  8.5295E-02  4.5653E-03  1.2941E-01  4.8761E-02  9.7672E-02 -1.1937E-02
             1.7654E+00
 GRADIENT:  -8.7324E+01 -2.7924E+01 -3.0767E+01  1.8897E+01  2.8597E+01 -2.0277E+01  6.3488E+00  8.2039E+00  1.6181E+01  2.2914E+01
             5.4673E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1830.03541895551        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.5455E-01  5.7568E-01  2.7624E-01  1.4068E+00  3.5330E-01  9.9740E-01  1.4957E+00  2.6414E-01  1.1212E+00  3.8414E-01
             4.6576E+00
 PARAMETER:  5.3482E-02 -4.5221E-01 -1.1865E+00  4.4133E-01 -9.4044E-01  9.7401E-02  5.0261E-01 -1.2313E+00  2.1436E-01 -8.5676E-01
             1.6385E+00
 GRADIENT:  -1.0029E+02  1.0684E+02 -1.6164E+01  2.1118E+02 -5.1396E+01 -1.0663E+01  1.4640E+01  6.1607E-01  3.2857E+00 -3.7887E-01
             4.3282E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1954.65939358438        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      272
 NPARAMETR:  9.5223E-01  5.5940E-01  3.5024E-01  1.2639E+00  4.0539E-01  1.0277E+00  1.0872E+00  1.4725E-02  1.0806E+00  6.9903E-01
             3.2680E+00
 PARAMETER:  5.1052E-02 -4.8089E-01 -9.4914E-01  3.3421E-01 -8.0290E-01  1.2733E-01  1.8360E-01 -4.1182E+00  1.7750E-01 -2.5806E-01
             1.2842E+00
 GRADIENT:  -5.8949E+01  6.8055E+01 -2.5413E+00  6.8778E+01 -3.8000E+01 -1.0364E+00 -4.6913E+00  2.4042E-03  2.8096E+00 -6.8453E+00
             6.0750E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1976.78451955010        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  9.7431E-01  2.9419E-01  1.8499E-01  1.1030E+00  2.3407E-01  1.0705E+00  6.8156E-01  1.0000E-02  1.2478E+00  8.6959E-01
             3.0181E+00
 PARAMETER:  7.3975E-02 -1.1235E+00 -1.5875E+00  1.9806E-01 -1.3521E+00  1.6814E-01 -2.8338E-01 -9.4429E+00  3.2140E-01 -3.9738E-02
             1.2046E+00
 GRADIENT:  -7.2187E+00  1.6562E+00 -1.1525E+01 -2.0069E-01  1.8274E+01  1.2111E+01  2.0124E-02  0.0000E+00 -1.6155E-01  1.2498E+01
            -1.0336E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1977.64284292407        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      624
 NPARAMETR:  9.7927E-01  2.8435E-01  1.8134E-01  1.0966E+00  2.2755E-01  1.0321E+00  7.1589E-01  1.0000E-02  1.2515E+00  8.1887E-01
             3.0418E+00
 PARAMETER:  7.9055E-02 -1.1575E+00 -1.6074E+00  1.9222E-01 -1.3804E+00  1.3162E-01 -2.3423E-01 -9.6052E+00  3.2435E-01 -9.9833E-02
             1.2125E+00
 GRADIENT:   5.3981E-01  6.7945E-02 -4.0426E-01 -3.5281E-01  1.2211E+00 -4.4121E-01 -5.2583E-01  0.0000E+00  7.9678E-02 -5.7208E-01
            -2.8982E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1977.67126442034        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      805
 NPARAMETR:  9.7901E-01  2.8276E-01  1.7976E-01  1.0948E+00  2.2595E-01  1.0337E+00  7.7544E-01  1.0000E-02  1.2547E+00  8.1536E-01
             3.0379E+00
 PARAMETER:  7.8790E-02 -1.1631E+00 -1.6161E+00  1.9053E-01 -1.3874E+00  1.3318E-01 -1.5432E-01 -9.7311E+00  3.2686E-01 -1.0412E-01
             1.2112E+00
 GRADIENT:   4.6206E-02  5.8004E-01  7.8929E-01 -8.8952E-02 -1.8027E+00  1.9434E-02  2.8720E-02  0.0000E+00  5.0931E-02  9.0901E-02
             3.9521E-02

0ITERATION NO.:   32    OBJECTIVE VALUE:  -1977.67145416766        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      862
 NPARAMETR:  9.7901E-01  2.8269E-01  1.7971E-01  1.0948E+00  2.2600E-01  1.0337E+00  7.7425E-01  1.0000E-02  1.2546E+00  8.1528E-01
             3.0379E+00
 PARAMETER:  7.8784E-02 -1.1634E+00 -1.6164E+00  1.9055E-01 -1.3872E+00  1.3316E-01 -1.5587E-01 -9.7311E+00  3.2680E-01 -1.0423E-01
             1.2112E+00
 GRADIENT:   3.6460E-02  2.1554E-01  2.8595E-01 -8.2566E-03 -7.7557E-01  1.1969E-02  8.8156E-03  0.0000E+00  2.3794E-02  3.1896E-02
            -1.3152E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      862
 NO. OF SIG. DIGITS IN FINAL EST.:  2.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7933E-03  3.2468E-03  1.6407E-04 -1.0796E-02  1.7688E-03
 SE:             2.9086E-02  1.2374E-02  2.1982E-04  2.6778E-02  2.6070E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5084E-01  7.9303E-01  4.5543E-01  6.8681E-01  9.4590E-01

 ETASHRINKSD(%)  2.5566E+00  5.8545E+01  9.9264E+01  1.0290E+01  1.2664E+01
 ETASHRINKVR(%)  5.0479E+00  8.2814E+01  9.9995E+01  1.9521E+01  2.3724E+01
 EBVSHRINKSD(%)  2.3750E+00  5.8846E+01  9.9353E+01  8.4055E+00  1.2995E+01
 EBVSHRINKVR(%)  4.6935E+00  8.3063E+01  9.9996E+01  1.6105E+01  2.4302E+01
 RELATIVEINF(%)  9.5244E+01  3.4235E+00  4.3058E-04  4.0448E+01  4.9953E+00
 EPSSHRINKSD(%)  2.0297E+01
 EPSSHRINKVR(%)  3.6475E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1977.6714541676606     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -691.15750768111889     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1977.671       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  2.83E-01  1.80E-01  1.09E+00  2.26E-01  1.03E+00  7.74E-01  1.00E-02  1.25E+00  8.15E-01  3.04E+00
 


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
 #CPUT: Total CPU Time in Seconds,       56.940
Stop Time:
Sat Oct 23 17:59:50 CDT 2021
