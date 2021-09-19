Sat Sep 18 02:45:14 CDT 2021
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
$DATA ../../../../data/int/S1/dat60.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3167.73394204539        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7173E+01 -2.5490E+01  4.8557E+01 -9.0513E+01  3.4709E+01  5.9603E+00 -1.0159E+01 -2.5019E+02 -8.2298E+01 -1.2972E-01
            -1.0983E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3636.92197238962        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.2723E-01  1.0201E+00  1.1705E+00  1.0487E+00  1.0564E+00  1.0214E+00  9.9429E-01  9.7933E-01  1.1188E+00  9.9845E-01
             1.6104E+00
 PARAMETER:  2.4444E-02  1.1993E-01  2.5740E-01  1.4758E-01  1.5488E-01  1.2118E-01  9.4270E-02  7.9113E-02  2.1229E-01  9.8444E-02
             5.7650E-01
 GRADIENT:  -1.5858E+02 -3.7900E+01  1.2010E+01  4.2878E+01  2.5019E+01 -4.2836E+00 -9.1515E-01  5.8115E+00  2.5381E+01  6.2443E+00
             5.5252E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3674.20908915888        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      165
 NPARAMETR:  9.6450E-01  1.0024E+00  8.0612E-01  1.0159E+00  9.6159E-01  9.5874E-01  7.7069E-01  2.5094E-01  1.0419E+00  1.2511E+00
             1.4576E+00
 PARAMETER:  6.3855E-02  1.0237E-01 -1.1552E-01  1.1577E-01  6.0833E-02  5.7862E-02 -1.6047E-01 -1.2825E+00  1.4100E-01  3.2402E-01
             4.7679E-01
 GRADIENT:  -8.1720E+01 -8.3271E-01 -1.7092E+02 -1.2693E+01  8.9587E+01 -2.0795E+01  8.6926E-01 -2.5878E-01 -1.5859E+01  2.9462E+01
             3.9028E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3748.40854692669        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      241
 NPARAMETR:  9.8183E-01  9.8235E-01  9.8816E-01  1.0458E+00  9.7894E-01  9.6970E-01  9.9166E-01  5.4799E-01  1.0293E+00  1.0343E+00
             1.1428E+00
 PARAMETER:  8.1660E-02  8.2196E-02  8.8085E-02  1.4479E-01  7.8711E-02  6.9234E-02  9.1621E-02 -5.0150E-01  1.2889E-01  1.3374E-01
             2.3345E-01
 GRADIENT:  -1.8645E+01 -7.4083E+00 -2.5590E+00  2.8851E+01  4.4345E+01 -1.0796E+01 -1.4345E+00 -9.0122E+00 -1.3423E+01  4.0424E+00
            -6.7501E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3749.16023888847        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      404
 NPARAMETR:  9.8221E-01  9.8857E-01  9.9694E-01  1.0423E+00  9.8519E-01  9.7126E-01  9.9368E-01  5.9535E-01  1.0315E+00  1.0301E+00
             1.1413E+00
 PARAMETER:  8.2049E-02  8.8507E-02  9.6940E-02  1.4143E-01  8.5080E-02  7.0838E-02  9.3657E-02 -4.1861E-01  1.3104E-01  1.2970E-01
             2.3216E-01
 GRADIENT:  -5.3380E+01 -1.2993E+01 -3.2233E+00  1.5723E+01  3.8383E+01 -1.7768E+01 -1.8953E+00 -8.8276E+00 -1.4100E+01  3.4950E+00
            -6.8099E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3750.66188855510        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      582
 NPARAMETR:  9.9457E-01  9.8910E-01  9.8967E-01  1.0420E+00  9.6023E-01  9.9942E-01  1.0103E+00  5.9939E-01  1.0317E+00  1.0044E+00
             1.1412E+00
 PARAMETER:  9.4555E-02  8.9039E-02  8.9614E-02  1.4111E-01  5.9422E-02  9.9416E-02  1.1027E-01 -4.1185E-01  1.3116E-01  1.0435E-01
             2.3211E-01
 GRADIENT:  -2.2657E+01  5.9763E+00  5.5465E+00  7.9042E+00  4.6897E+00 -4.9143E+00 -7.9590E-01 -8.4041E+00 -1.3315E+01  8.4568E-01
            -5.7052E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3750.86052321435        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      758
 NPARAMETR:  1.0049E+00  9.8910E-01  9.7891E-01  1.0420E+00  9.5267E-01  1.0111E+00  1.0175E+00  5.9939E-01  1.0317E+00  9.9608E-01
             1.1412E+00
 PARAMETER:  1.0488E-01  8.9039E-02  7.8683E-02  1.4111E-01  5.1516E-02  1.1102E-01  1.1737E-01 -4.1185E-01  1.3116E-01  9.6076E-02
             2.3211E-01
 GRADIENT:   2.0011E-01  1.2250E+01  2.8391E-02  5.8985E+00 -1.0257E-02 -4.8380E-02 -7.7222E-02 -7.8350E+00 -1.3219E+01 -1.1562E-01
            -5.0908E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3751.43427810083        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      901
 NPARAMETR:  1.0047E+00  9.8769E-01  9.7896E-01  1.0415E+00  9.5268E-01  1.0111E+00  1.0179E+00  6.4036E-01  1.0381E+00  9.9648E-01
             1.1416E+00
 PARAMETER:  1.0471E-01  8.7611E-02  7.8734E-02  1.4063E-01  5.1520E-02  1.1101E-01  1.1777E-01 -3.4572E-01  1.3741E-01  9.6473E-02
             2.3241E-01
 GRADIENT:   3.7456E+01  1.5571E+01 -2.6194E+00  1.6220E+01  3.8879E+00  8.0494E+00  7.6906E-01 -6.5050E+00 -9.9154E+00  1.3607E+00
            -1.3690E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3751.44505982622        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:     1055
 NPARAMETR:  1.0047E+00  9.8767E-01  9.8066E-01  1.0393E+00  9.5266E-01  1.0093E+00  1.0170E+00  6.4034E-01  1.0381E+00  9.9322E-01
             1.1416E+00
 PARAMETER:  1.0472E-01  8.7598E-02  8.0468E-02  1.3855E-01  5.1508E-02  1.0927E-01  1.1686E-01 -3.4576E-01  1.3740E-01  9.3193E-02
             2.3244E-01
 GRADIENT:  -1.0659E-01  8.4317E+00 -1.4256E+00  8.5754E-01 -1.8980E+00 -7.4230E-01 -1.7177E-01 -6.7297E+00 -1.1323E+01  2.7912E-01
            -6.3239E-01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -3751.44755232026        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     1185
 NPARAMETR:  1.0047E+00  9.8767E-01  9.8209E-01  1.0388E+00  9.5266E-01  1.0111E+00  1.0196E+00  6.4034E-01  1.0381E+00  9.9117E-01
             1.1416E+00
 PARAMETER:  1.0472E-01  8.7598E-02  8.1924E-02  1.3810E-01  5.1508E-02  1.1100E-01  1.1937E-01 -3.4576E-01  1.3740E-01  9.1127E-02
             2.3244E-01
 GRADIENT:   8.7267E+00  1.4609E+06  2.9219E+06  2.1156E+06 -4.1295E+00  2.6323E+06  2.4477E+06  8.4488E+05 -2.0545E+02  2.9218E+06
            -6.7133E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1185
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.4106E-04 -2.8393E-02 -1.4259E-02  1.1416E-02 -2.0611E-02
 SE:             2.9880E-02  2.2391E-02  1.4122E-02  2.8686E-02  2.5123E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8822E-01  2.0479E-01  3.1265E-01  6.9067E-01  4.1199E-01

 ETASHRINKSD(%)  1.0000E-10  2.4986E+01  5.2690E+01  3.8968E+00  1.5835E+01
 ETASHRINKVR(%)  1.0000E-10  4.3730E+01  7.7617E+01  7.6418E+00  2.9162E+01
 EBVSHRINKSD(%)  3.2018E-01  2.4915E+01  5.7170E+01  7.5563E+00  1.6218E+01
 EBVSHRINKVR(%)  6.3934E-01  4.3623E+01  8.1656E+01  1.4542E+01  2.9805E+01
 RELATIVEINF(%)  9.9358E+01  2.5178E+01  9.9028E+00  5.7998E+01  2.3889E+01
 EPSSHRINKSD(%)  2.0317E+01
 EPSSHRINKVR(%)  3.6507E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3751.4475523202559     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2097.3581925518451     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.92
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3751.448       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  9.88E-01  9.82E-01  1.04E+00  9.53E-01  1.01E+00  1.02E+00  6.40E-01  1.04E+00  9.91E-01  1.14E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.32E+10
 
 TH 2
+        2.23E+04  7.49E+09
 
 TH 3
+        7.07E+09  3.28E+05  7.57E+09
 
 TH 4
+        1.53E+04 -1.88E+05  2.26E+05  3.55E+09
 
 TH 5
+       -7.29E+09 -3.54E+03 -3.58E+02  5.34E+09  1.61E+10
 
 TH 6
+       -1.26E+01  3.49E+04  3.51E+04  2.40E+04 -3.76E+00  5.80E+09
 
 TH 7
+        5.70E+09 -6.76E+03  6.11E+09 -4.66E+03 -1.18E+00  5.35E+09  4.93E+09
 
 TH 8
+        9.92E+03 -3.34E+09  1.46E+05 -2.30E+09 -1.43E+03 -2.94E+09 -2.71E+09  1.49E+09
 
 TH 9
+       -4.87E+09  5.18E+09 -4.04E+01  3.57E+09  5.38E+09  1.13E+01  1.11E+01 -2.11E+05  3.59E+09
 
 TH10
+       -7.18E-01 -5.53E+02  1.50E+10 -3.68E+02 -1.33E+01  1.31E+10  1.21E+10  3.33E+09 -5.00E+00  1.49E+10
 
 TH11
+       -2.62E+09  2.79E+09 -2.80E+09  1.92E+09  2.89E+09 -2.45E+09  1.10E+01 -6.63E+05  1.76E+05 -2.78E+09  1.04E+09
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       43.826
Stop Time:
Sat Sep 18 02:46:00 CDT 2021
