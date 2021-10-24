Sun Oct 24 00:59:14 CDT 2021
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
$DATA ../../../../data/SD3/D/dat43.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m43.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2134.08542308304        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4011E+02  5.0592E+00 -6.2155E+00  3.2925E+01  1.0404E+01  2.1836E+01  9.6929E+00  8.8335E+00  3.6706E+01  4.4886E+00
             6.5159E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2141.02920693146        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.7768E-01  1.0075E+00  1.0158E+00  1.0180E+00  1.0115E+00  1.0803E+00  9.7518E-01  9.6032E-01  8.6667E-01  9.9714E-01
             9.0256E-01
 PARAMETER:  7.7423E-02  1.0747E-01  1.1571E-01  1.1787E-01  1.1144E-01  1.7725E-01  7.4865E-02  5.9507E-02 -4.3095E-02  9.7137E-02
            -2.5212E-03
 GRADIENT:  -9.0684E+00 -3.4961E+00 -2.1034E+00 -1.2131E+01  9.9985E+00  8.2998E+00 -3.0078E+00  4.5886E+00  2.7822E+00 -3.8289E+00
            -1.4261E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2143.49126840306        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.7618E-01  8.2281E-01  1.0460E+00  1.1514E+00  9.2708E-01  1.0467E+00  1.2024E+00  7.3106E-01  7.4999E-01  9.9305E-01
             8.9848E-01
 PARAMETER:  7.5896E-02 -9.5027E-02  1.4498E-01  2.4100E-01  2.4286E-02  1.4568E-01  2.8435E-01 -2.1326E-01 -1.8769E-01  9.3029E-02
            -7.0479E-03
 GRADIENT:  -9.1910E+00  2.6609E+01  1.8336E+01  2.6992E+01 -1.7917E+01 -3.8310E+00 -5.1761E+00 -1.4769E+00 -6.2261E+00 -1.0954E+00
            -1.8503E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2145.41859154852        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.8124E-01  7.0546E-01  9.4291E-01  1.2079E+00  8.3430E-01  1.0605E+00  1.4292E+00  5.0664E-01  7.2266E-01  9.3061E-01
             9.2207E-01
 PARAMETER:  8.1061E-02 -2.4890E-01  4.1219E-02  2.8886E-01 -8.1158E-02  1.5870E-01  4.5709E-01 -5.7995E-01 -2.2482E-01  2.8087E-02
             1.8868E-02
 GRADIENT:   2.3602E+00  1.6909E+01  7.0143E+00  2.6126E+01 -1.2417E+01  1.6000E+00 -1.0718E+00 -9.3707E-01 -1.9812E+00  2.0545E-02
             3.9470E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2147.02973392233        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  9.7431E-01  4.2950E-01  1.0973E+00  1.3690E+00  8.2968E-01  1.0512E+00  1.9939E+00  7.3841E-01  6.7790E-01  9.3606E-01
             9.1272E-01
 PARAMETER:  7.3973E-02 -7.4514E-01  1.9289E-01  4.1408E-01 -8.6715E-02  1.4990E-01  7.9007E-01 -2.0326E-01 -2.8876E-01  3.3924E-02
             8.6714E-03
 GRADIENT:  -1.4860E+00  2.7581E+00 -4.5969E+00  2.7091E+00  8.4890E+00 -2.6840E-01  1.0721E+00 -7.9284E-01  1.3327E+00 -1.4403E+00
            -5.6567E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2147.79199348299        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  9.7105E-01  2.2351E-01  1.2976E+00  1.5125E+00  8.4005E-01  1.0487E+00  2.8510E+00  1.0092E+00  6.4005E-01  9.6389E-01
             9.2407E-01
 PARAMETER:  7.0625E-02 -1.3983E+00  3.6048E-01  5.1378E-01 -7.4295E-02  1.4758E-01  1.1477E+00  1.0913E-01 -3.4622E-01  6.3220E-02
             2.1037E-02
 GRADIENT:  -5.5165E-02  5.8175E+00  9.5362E+00  2.1263E+01 -1.7743E+01  4.1484E-01 -1.8348E-01  1.1273E+00 -1.3361E+00  1.0278E+00
             5.7113E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2149.18803134507        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  9.6766E-01  7.6341E-02  1.3480E+00  1.5989E+00  8.3333E-01  1.0440E+00  4.9065E+00  1.0139E+00  6.1672E-01  9.8277E-01
             9.1830E-01
 PARAMETER:  6.7125E-02 -2.4725E+00  3.9860E-01  5.6929E-01 -8.2322E-02  1.4304E-01  1.6906E+00  1.1381E-01 -3.8333E-01  8.2620E-02
             1.4768E-02
 GRADIENT:  -1.0986E+00  1.0067E+00  3.0455E+00  3.2847E+01 -6.0559E+00 -3.9804E-01 -1.1938E+00 -2.1211E+00  1.7490E+00  5.5370E-01
            -7.8372E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2149.79713364060        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1234
 NPARAMETR:  9.6792E-01  5.6451E-02  1.3576E+00  1.5986E+00  8.3459E-01  1.0448E+00  5.5566E+00  1.0639E+00  6.0840E-01  9.6467E-01
             9.1858E-01
 PARAMETER:  6.7396E-02 -2.7744E+00  4.0570E-01  5.6910E-01 -8.0814E-02  1.4381E-01  1.8150E+00  1.6193E-01 -3.9693E-01  6.4030E-02
             1.5069E-02
 GRADIENT:   2.2030E-01 -6.3626E-01  1.8411E-01 -1.4065E+00  2.0776E-01  4.4383E-03 -1.2965E+00 -1.9206E-01  1.3220E+00 -2.0201E-01
             1.6042E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -2149.81085177416        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1332
 NPARAMETR:  9.6853E-01  5.5360E-02  1.3643E+00  1.6032E+00  8.3663E-01  1.0452E+00  5.5781E+00  1.0723E+00  6.0783E-01  9.6729E-01
             9.1825E-01
 PARAMETER:  6.7165E-02 -2.7847E+00  4.1067E-01  5.7060E-01 -7.8490E-02  1.4375E-01  1.8244E+00  1.7157E-01 -3.9886E-01  6.7105E-02
             1.4646E-02
 GRADIENT:  -1.2704E+00  4.5354E+01  1.9885E-02 -9.3529E+01 -2.3362E-01 -2.7934E-01  6.6616E+01  9.2261E-02 -1.3874E+00  3.9194E-02
            -7.5114E-02
 NUMSIGDIG:         1.9         2.3         4.1         2.4         2.8         2.3         2.3         1.8         2.4         2.3
                    3.0

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1332
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4694E-03  2.7290E-02 -3.0184E-02 -2.2528E-02 -2.7239E-02
 SE:             2.9957E-02  1.2238E-02  1.7647E-02  2.7465E-02  2.1263E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6088E-01  2.5759E-02  8.7182E-02  4.1208E-01  2.0017E-01

 ETASHRINKSD(%)  1.0000E-10  5.9000E+01  4.0881E+01  7.9891E+00  2.8766E+01
 ETASHRINKVR(%)  1.0000E-10  8.3190E+01  6.5049E+01  1.5340E+01  4.9257E+01
 EBVSHRINKSD(%)  2.5164E-01  7.1356E+01  4.3483E+01  5.6146E+00  2.4759E+01
 EBVSHRINKVR(%)  5.0264E-01  9.1795E+01  6.8058E+01  1.0914E+01  4.3388E+01
 RELATIVEINF(%)  9.9344E+01  3.7043E+00  7.6416E+00  3.9170E+01  1.3869E+01
 EPSSHRINKSD(%)  3.4312E+01
 EPSSHRINKVR(%)  5.6850E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2149.8108517741584     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1230.8723185694857     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2149.811       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  5.59E-02  1.36E+00  1.60E+00  8.37E-01  1.04E+00  5.61E+00  1.07E+00  6.07E-01  9.68E-01  9.18E-01
 


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
 #CPUT: Total CPU Time in Seconds,       48.538
Stop Time:
Sun Oct 24 00:59:24 CDT 2021
