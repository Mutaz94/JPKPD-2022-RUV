Sat Sep 18 13:17:15 CDT 2021
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
$DATA ../../../../data/spa/S2/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1669.30938282706        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -4.1907E+01 -7.7054E+01  4.0450E+01 -1.3521E+02 -2.6906E+01  4.8750E+00 -2.6454E+01 -1.3472E+01  2.7915E+00 -1.8591E+01
             7.5669E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1688.42657482760        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0603E+00  1.1629E+00  8.4470E-01  9.7793E-01  1.0455E+00  9.7073E-01  1.3239E+00  1.0940E+00  7.9749E-01  1.1580E+00
             9.4360E-01
 PARAMETER:  1.5855E-01  2.5089E-01 -6.8780E-02  7.7686E-02  1.4445E-01  7.0294E-02  3.8060E-01  1.8982E-01 -1.2628E-01  2.4669E-01
             4.1944E-02
 GRADIENT:   1.1741E+02  1.6785E+00 -1.3114E+01  1.1490E+01  9.2237E+00 -6.2299E+00  2.5043E-01  5.7817E+00  6.7142E+00  1.0244E+01
            -9.1932E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1689.14325343065        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0467E+00  9.9877E-01  9.6796E-01  1.0668E+00  1.0120E+00  9.7411E-01  1.6546E+00  1.0441E+00  6.3838E-01  1.1456E+00
             9.9688E-01
 PARAMETER:  1.4567E-01  9.8765E-02  6.7438E-02  1.6470E-01  1.1196E-01  7.3773E-02  6.0354E-01  1.4319E-01 -3.4882E-01  2.3589E-01
             9.6870E-02
 GRADIENT:   7.8396E+01  2.4422E+00  1.6574E+01 -3.2597E+01 -2.0248E+01 -3.1898E+00  1.3889E+01 -9.1454E-01 -2.4360E+00  5.9799E+00
             1.0857E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1690.38088781883        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0221E+00  1.0372E+00  7.5806E-01  1.0401E+00  9.1147E-01  9.8508E-01  1.5001E+00  7.3977E-01  7.0016E-01  1.0326E+00
             9.5602E-01
 PARAMETER:  1.2181E-01  1.3650E-01 -1.7699E-01  1.3932E-01  7.3042E-03  8.4966E-02  5.0551E-01 -2.0142E-01 -2.5645E-01  1.3209E-01
             5.5022E-02
 GRADIENT:   1.1775E+01  2.6469E+00 -5.0672E+00  4.2697E+00  2.3893E+00  3.9730E-02  4.0995E+00  2.2191E+00 -4.3403E-01  4.2248E+00
            -2.9043E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1690.43924312028        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0196E+00  1.0373E+00  6.7592E-01  1.0301E+00  8.5026E-01  9.8613E-01  1.4613E+00  5.6798E-01  7.1951E-01  9.5561E-01
             9.5849E-01
 PARAMETER:  1.1939E-01  1.3661E-01 -2.9167E-01  1.2968E-01 -6.2213E-02  8.6035E-02  4.7933E-01 -4.6567E-01 -2.2919E-01  5.4591E-02
             5.7608E-02
 GRADIENT:   2.2969E+00  7.9416E-01 -3.9960E+00  3.5278E+00  1.4920E+00 -1.9259E-01  1.4315E+00  1.7091E+00  2.7205E-01  2.2087E+00
            -1.1374E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1690.44121748569        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0192E+00  1.0363E+00  6.4100E-01  1.0254E+00  8.2348E-01  9.8666E-01  1.4487E+00  4.8072E-01  7.2530E-01  9.2287E-01
             9.5926E-01
 PARAMETER:  1.1901E-01  1.3563E-01 -3.4473E-01  1.2507E-01 -9.4215E-02  8.6574E-02  4.7064E-01 -6.3248E-01 -2.2116E-01  1.9738E-02
             5.8410E-02
 GRADIENT:   1.7117E-01  1.0076E-01 -3.1607E+00  2.4972E+00  1.0156E+00 -1.9801E-01  6.8170E-01  1.3461E+00  3.5156E-01  1.4276E+00
            -5.4935E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1690.44287625077        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      463
 NPARAMETR:  1.0191E+00  1.0351E+00  6.1808E-01  1.0223E+00  8.0550E-01  9.8701E-01  1.4418E+00  4.1563E-01  7.2842E-01  9.0152E-01
             9.5948E-01
 PARAMETER:  1.1896E-01  1.3452E-01 -3.8114E-01  1.2206E-01 -1.1629E-01  8.6923E-02  4.6589E-01 -7.7795E-01 -2.1688E-01 -3.6778E-03
             5.8641E-02
 GRADIENT:  -5.3192E+01 -4.8333E+00 -3.8679E+00 -7.9104E+00 -5.4253E-01 -4.9630E+00 -2.3911E+00  1.0337E+00 -4.7508E-01  9.0164E-01
            -3.5825E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1691.42344697395        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  1.0423E+00  9.1319E-01  7.0193E-01  1.1065E+00  8.1031E-01  9.9430E-01  1.6459E+00  4.5387E-01  6.8632E-01  9.5107E-01
             9.6094E-01
 PARAMETER:  1.4138E-01  9.1914E-03 -2.5392E-01  2.0122E-01 -1.1034E-01  9.4284E-02  5.9832E-01 -6.8994E-01 -2.7641E-01  4.9835E-02
             6.0152E-02
 GRADIENT:   9.6943E-01  1.7147E+00  1.8840E+00  8.2556E-01 -6.6160E-01 -6.4234E-02  1.1672E-01 -3.1754E-01 -7.0079E-01 -6.9660E-01
            -7.6793E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1691.43978626108        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      816
 NPARAMETR:  1.0413E+00  8.7777E-01  7.3184E-01  1.1287E+00  8.1540E-01  9.9351E-01  1.6983E+00  5.3816E-01  6.7711E-01  9.5950E-01
             9.6200E-01
 PARAMETER:  1.4051E-01 -3.0371E-02 -2.1220E-01  2.2107E-01 -1.0407E-01  9.3489E-02  6.2961E-01 -5.1961E-01 -2.8992E-01  5.8654E-02
             6.1258E-02
 GRADIENT:   3.5071E-01  3.5795E-01  6.2923E-01 -4.6413E-01 -1.1826E+00 -8.9944E-02 -1.0219E-01  2.9065E-02 -2.8461E-01  2.3562E-02
            -1.2522E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1691.44102227793        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      978
 NPARAMETR:  1.0410E+00  8.7184E-01  7.3995E-01  1.1329E+00  8.1924E-01  9.9359E-01  1.7087E+00  5.5110E-01  6.7672E-01  9.6403E-01
             9.6232E-01
 PARAMETER:  1.4021E-01 -3.7147E-02 -2.0117E-01  2.2477E-01 -9.9377E-02  9.3567E-02  6.3574E-01 -4.9583E-01 -2.9050E-01  6.3368E-02
             6.1590E-02
 GRADIENT:   7.6250E-03  1.1347E-03 -9.1744E-03  1.1979E-02  3.4616E-03  1.0763E-03 -2.2067E-03  1.4225E-03 -1.4197E-03  3.4744E-03
            -2.2595E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      978
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.9983E-04  1.2104E-02 -2.7321E-02 -1.8133E-02 -1.1382E-02
 SE:             2.9861E-02  2.3460E-02  1.0749E-02  2.1659E-02  2.2003E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8665E-01  6.0590E-01  1.1028E-02  4.0248E-01  6.0497E-01

 ETASHRINKSD(%)  1.0000E-10  2.1405E+01  6.3990E+01  2.7439E+01  2.6287E+01
 ETASHRINKVR(%)  1.0000E-10  3.8228E+01  8.7033E+01  4.7348E+01  4.5663E+01
 EBVSHRINKSD(%)  4.1881E-01  2.0543E+01  6.7984E+01  2.8179E+01  2.3769E+01
 EBVSHRINKVR(%)  8.3586E-01  3.6867E+01  8.9750E+01  4.8417E+01  4.1889E+01
 RELATIVEINF(%)  9.8707E+01  8.3614E+00  1.3695E+00  6.2720E+00  7.1632E+00
 EPSSHRINKSD(%)  4.4965E+01
 EPSSHRINKVR(%)  6.9711E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1691.4410222779297     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -956.29019571419155     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.70
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1691.441       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  8.72E-01  7.40E-01  1.13E+00  8.19E-01  9.94E-01  1.71E+00  5.51E-01  6.77E-01  9.64E-01  9.62E-01
 


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
+        1.03E+03
 
 TH 2
+       -8.59E+00  3.29E+02
 
 TH 3
+        1.94E+01  1.27E+02  7.58E+02
 
 TH 4
+       -7.74E+00  3.28E+02 -4.08E+02  9.92E+02
 
 TH 5
+       -3.30E+00 -2.12E+02 -7.59E+02  3.96E+02  1.04E+03
 
 TH 6
+        2.56E-01 -2.44E+00  1.95E+00 -9.61E-01 -3.84E+00  1.99E+02
 
 TH 7
+        1.97E+00  2.89E+01 -1.15E+01 -2.11E+01  5.19E+00 -1.84E-02  3.14E+01
 
 TH 8
+       -7.29E-01 -1.11E+01 -7.19E+01  8.19E+00  1.33E+01 -5.20E-01  2.46E+00  2.23E+01
 
 TH 9
+        1.04E+00 -2.14E+01 -4.23E+01 -3.59E+00  2.41E+01 -8.96E-01  1.84E+01  1.00E+01  1.34E+02
 
 TH10
+       -3.69E+00 -3.78E+00 -3.87E+01 -1.96E+01 -6.34E+01  1.19E+00  1.43E+00  2.00E+01  1.67E+01  7.59E+01
 
 TH11
+       -6.76E+00 -9.06E+00 -2.93E+01 -3.29E+00  3.48E-01  3.72E+00  2.43E+00  6.49E+00  1.45E+01  1.49E+01  2.25E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.547
Stop Time:
Sat Sep 18 13:17:32 CDT 2021
