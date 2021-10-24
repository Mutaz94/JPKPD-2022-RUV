Sat Oct 23 14:27:19 CDT 2021
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
$DATA ../../../../data/SD1/A3/dat100.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1450.82135328482        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1502E+02  1.6941E+02  1.7628E+02 -1.7683E+01  2.5817E+02  4.1044E+01 -5.0677E+01 -1.6779E+02 -5.3226E+01 -1.0412E+02
            -4.3338E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2859.55155803032        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.5429E-01  8.5446E-01  7.9118E-01  1.0805E+00  7.6983E-01  8.8791E-01  9.5497E-01  1.0703E+00  9.3700E-01  1.0419E+00
             2.5649E+00
 PARAMETER:  5.3211E-02 -5.7284E-02 -1.3423E-01  1.7743E-01 -1.6159E-01 -1.8883E-02  5.3926E-02  1.6798E-01  3.4932E-02  1.4101E-01
             1.0419E+00
 GRADIENT:   5.7868E+01  2.8697E+00  4.8581E+00 -4.9692E+01 -5.2764E-01 -7.1118E+00  4.2519E+00  1.2468E+01 -3.2195E+00  1.0515E+01
             2.4464E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2868.69215339544        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.5927E-01  5.1734E-01  4.6460E-01  1.2799E+00  4.3275E-01  1.0366E+00  9.7832E-01  8.7285E-01  1.0316E+00  8.5812E-01
             2.4384E+00
 PARAMETER:  5.8419E-02 -5.5905E-01 -6.6657E-01  3.4682E-01 -7.3759E-01  1.3599E-01  7.8081E-02 -3.5989E-02  1.3113E-01 -5.3006E-02
             9.9134E-01
 GRADIENT:   7.5026E+01  1.0202E+02  5.8525E+01  1.5887E+02  2.4440E+01  4.6038E+01 -1.5315E+01  1.8327E+01  1.4539E+01  6.5911E+00
             2.0439E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2873.45912714066        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      314
 NPARAMETR:  9.7521E-01  5.6818E-01  4.8431E-01  1.2179E+00  4.7500E-01  8.9436E-01  1.4320E+00  6.1711E-01  9.8191E-01  6.9865E-01
             2.4975E+00
 PARAMETER:  7.4893E-02 -4.6531E-01 -6.2502E-01  2.9713E-01 -6.4445E-01 -1.1646E-02  4.5909E-01 -3.8271E-01  8.1740E-02 -2.5861E-01
             1.0153E+00
 GRADIENT:   6.1220E+01  5.3351E+01  1.2690E+01  1.4677E+01 -3.4794E+01 -1.1225E+01  1.5221E+01  4.7173E+00  4.5536E+00  2.3635E+00
             2.1936E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2924.08171894112        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      491
 NPARAMETR:  9.1965E-01  3.1049E-01  2.3485E-01  1.1798E+00  2.7224E-01  9.5007E-01  1.4373E+00  5.3501E-01  1.2161E+00  5.3846E-01
             1.8833E+00
 PARAMETER:  1.6233E-02 -1.0696E+00 -1.3488E+00  2.6537E-01 -1.2011E+00  4.8777E-02  4.6278E-01 -5.2547E-01  2.9569E-01 -5.1903E-01
             7.3301E-01
 GRADIENT:  -6.4254E+01  1.6000E+01 -5.6096E+01  5.4364E+01  1.3313E+02  6.6341E+00 -1.7389E+01 -2.4148E+01 -1.3293E+01 -2.4545E+01
            -1.7398E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2953.98898625690        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      671
 NPARAMETR:  9.4407E-01  2.5482E-01  1.8321E-01  1.1108E+00  2.2012E-01  9.2315E-01  1.5053E+00  9.0585E-01  1.4710E+00  5.2623E-01
             1.9304E+00
 PARAMETER:  4.2443E-02 -1.2672E+00 -1.5971E+00  2.0512E-01 -1.4136E+00  2.0035E-02  5.0896E-01  1.1141E-03  4.8595E-01 -5.4201E-01
             7.5774E-01
 GRADIENT:   2.7281E-01  1.3582E+01  1.0716E+01  2.8686E+01 -1.3024E+01  1.5134E+00 -4.6720E+00  1.0800E+00  1.5397E+00 -2.9216E+00
            -8.5134E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2956.92368872743        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  9.4410E-01  2.2264E-01  1.5257E-01  1.0187E+00  1.9479E-01  9.1671E-01  1.4909E+00  7.7992E-01  1.5413E+00  6.4843E-01
             1.9295E+00
 PARAMETER:  4.2475E-02 -1.4022E+00 -1.7802E+00  1.1848E-01 -1.5358E+00  1.3039E-02  4.9937E-01 -1.4857E-01  5.3264E-01 -3.3320E-01
             7.5728E-01
 GRADIENT:  -2.7409E-01  8.9541E-01  1.0525E+00  6.5682E-02 -2.3778E+00  3.3287E-01 -2.9272E-01 -1.7394E-01 -2.1409E-01  1.1235E-01
            -1.4840E+00

0ITERATION NO.:   34    OBJECTIVE VALUE:  -2956.92671428718        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:      980
 NPARAMETR:  9.4424E-01  2.2271E-01  1.5245E-01  1.0185E+00  1.9471E-01  9.1574E-01  1.4929E+00  7.8237E-01  1.5438E+00  6.4774E-01
             1.9309E+00
 PARAMETER:  4.2627E-02 -1.4019E+00 -1.7809E+00  1.1832E-01 -1.5363E+00  1.1979E-02  5.0072E-01 -1.4542E-01  5.3426E-01 -3.3427E-01
             7.5797E-01
 GRADIENT:   9.3582E-02  1.4923E+00  7.3236E-01  6.9043E-02 -2.6466E+00 -5.7020E-02  4.3464E-02  4.2206E-02  2.0486E-01  1.0185E-01
             1.2350E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      980
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4088E-03  1.1055E-02  9.9848E-03 -9.3910E-04  1.3870E-02
 SE:             2.9454E-02  2.5971E-02  1.7920E-02  2.8283E-02  2.2973E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6185E-01  6.7035E-01  5.7740E-01  9.7351E-01  5.4600E-01

 ETASHRINKSD(%)  1.3258E+00  1.2994E+01  3.9966E+01  5.2472E+00  2.3037E+01
 ETASHRINKVR(%)  2.6339E+00  2.4300E+01  6.3959E+01  1.0219E+01  4.0767E+01
 EBVSHRINKSD(%)  1.4392E+00  1.1671E+01  4.0025E+01  4.5843E+00  2.4138E+01
 EBVSHRINKVR(%)  2.8577E+00  2.1980E+01  6.4030E+01  8.9584E+00  4.2450E+01
 RELATIVEINF(%)  9.7132E+01  2.9533E+01  5.4069E+00  5.8064E+01  7.9846E+00
 EPSSHRINKSD(%)  2.1760E+01
 EPSSHRINKVR(%)  3.8784E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2956.9267142871845     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1302.8373545187737     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.20
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2956.927       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.44E-01  2.23E-01  1.52E-01  1.02E+00  1.95E-01  9.16E-01  1.49E+00  7.82E-01  1.54E+00  6.48E-01  1.93E+00
 


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
 #CPUT: Total CPU Time in Seconds,       65.644
Stop Time:
Sat Oct 23 14:27:31 CDT 2021
