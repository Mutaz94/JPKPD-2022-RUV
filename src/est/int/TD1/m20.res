Wed Sep 29 06:06:45 CDT 2021
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
$DATA ../../../../data/int/TD1/dat20.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3256.13533424710        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9873E+02  4.9444E+01  1.1533E+02  9.1760E+01  1.5202E+02  3.4668E+01 -9.1466E+00 -2.4968E+02 -3.2917E+01 -4.7156E+00
            -8.7199E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3623.57707603497        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      129
 NPARAMETR:  7.6983E-01  9.5521E-01  9.1642E-01  9.4581E-01  9.0146E-01  9.3838E-01  1.0098E+00  1.1596E+00  1.0106E+00  9.8358E-01
             1.5301E+00
 PARAMETER: -1.6158E-01  5.4172E-02  1.2714E-02  4.4283E-02 -3.7390E-03  3.6401E-02  1.0978E-01  2.4806E-01  1.1056E-01  8.3442E-02
             5.2536E-01
 GRADIENT:  -5.8604E+02 -7.1584E+01 -7.2782E+00 -1.0866E+02 -2.6578E+01 -1.6889E+02  1.8358E+00  2.3681E+01 -4.0116E-01  2.0570E+01
             6.1834E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3643.36111800925        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  7.8661E-01  8.3248E-01  6.8490E-01  1.0101E+00  7.2877E-01  9.2269E-01  1.3097E+00  6.4817E-01  1.1279E+00  6.6524E-01
             1.4850E+00
 PARAMETER: -1.4003E-01 -8.3350E-02 -2.7848E-01  1.1009E-01 -2.1640E-01  1.9535E-02  3.6977E-01 -3.3361E-01  2.2035E-01 -3.0760E-01
             4.9538E-01
 GRADIENT:  -5.4922E+02  3.7098E+00 -8.2265E+01 -6.8606E+01 -3.0219E+01 -1.5602E+02  2.1289E+01  8.2263E+00  2.6846E+01  7.0529E+00
             5.6721E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3810.39982762300        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      483
 NPARAMETR:  9.0708E-01  9.8478E-01  8.9321E-01  9.7799E-01  8.9691E-01  9.7379E-01  9.0471E-01  4.2155E-01  1.0227E+00  1.0548E+00
             1.1473E+00
 PARAMETER:  2.4704E-03  8.4660E-02 -1.2931E-02  7.7740E-02 -8.7965E-03  7.3437E-02 -1.3796E-04 -7.6381E-01  1.2246E-01  1.5335E-01
             2.3737E-01
 GRADIENT:  -1.5921E+02 -3.6126E+00  2.5395E+01 -2.7254E+01 -2.4036E+01 -2.7163E+01 -1.4424E+01 -1.0013E+01 -1.5109E+00  2.1791E+01
             1.6763E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3830.55838308960        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  9.2306E-01  9.6741E-01  9.1176E-01  9.8939E-01  9.0413E-01  9.7662E-01  1.0941E+00  7.8402E-01  1.0005E+00  9.8058E-01
             1.0960E+00
 PARAMETER:  1.9938E-02  6.6868E-02  7.6165E-03  8.9328E-02 -7.7935E-04  7.6343E-02  1.8989E-01 -1.4331E-01  1.0054E-01  8.0393E-02
             1.9163E-01
 GRADIENT:  -1.1725E+02 -1.3121E+01  2.0575E+01 -1.4567E+01 -8.6951E+00 -2.1161E+01  4.6531E+00 -1.2747E+01 -2.1565E-01  1.7138E+01
             1.4338E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3831.44598022608        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      833
 NPARAMETR:  9.2390E-01  9.6746E-01  9.1169E-01  9.8946E-01  9.0393E-01  1.0285E+00  1.0949E+00  7.9132E-01  1.0003E+00  9.7930E-01
             1.0952E+00
 PARAMETER:  2.0847E-02  6.6921E-02  7.5439E-03  8.9404E-02 -1.0031E-03  1.2809E-01  1.9067E-01 -1.3406E-01  1.0034E-01  7.9086E-02
             1.9096E-01
 GRADIENT:  -1.0359E+02 -1.2887E+01  2.0180E+01 -1.4711E+01 -9.1161E+00  8.8310E-01  4.7262E+00 -1.2490E+01 -1.7863E-01  1.7009E+01
             1.4303E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3832.34518355446        NO. OF FUNC. EVALS.: 202
 CUMULATIVE NO. OF FUNC. EVALS.:     1035
 NPARAMETR:  9.2911E-01  9.6754E-01  9.1153E-01  9.8954E-01  9.0448E-01  1.0281E+00  1.0906E+00  7.9367E-01  1.0003E+00  9.6803E-01
             1.0945E+00
 PARAMETER:  2.6474E-02  6.6997E-02  7.3702E-03  8.9485E-02 -3.9955E-04  1.2774E-01  1.8676E-01 -1.3109E-01  1.0033E-01  6.7504E-02
             1.9032E-01
 GRADIENT:  -9.1911E+01 -1.2644E+01  1.9760E+01 -1.4581E+01 -7.3283E+00  1.8419E+00  3.8343E+00 -1.2518E+01 -3.4158E-01  1.4843E+01
             1.4171E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3832.58322525681        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1231
 NPARAMETR:  9.2952E-01  9.6871E-01  9.1133E-01  9.8951E-01  9.0470E-01  1.0257E+00  1.0634E+00  7.9337E-01  1.0000E+00  9.5790E-01
             1.0946E+00
 PARAMETER:  2.6913E-02  6.8206E-02  7.1513E-03  8.9453E-02 -1.5372E-04  1.2533E-01  1.6149E-01 -1.3146E-01  1.0005E-01  5.6989E-02
             1.9038E-01
 GRADIENT:  -9.1477E+01 -1.2043E+01  1.9966E+01 -1.4705E+01 -6.5755E+00  9.3361E-01 -1.0071E-01 -1.2636E+01 -9.9071E-01  1.2172E+01
             1.4042E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3832.66808570379        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1420
 NPARAMETR:  9.2952E-01  9.6929E-01  9.1133E-01  9.8951E-01  9.0470E-01  1.0255E+00  1.0634E+00  7.9442E-01  1.0000E+00  9.5286E-01
             1.0946E+00
 PARAMETER:  2.6913E-02  6.8813E-02  7.1522E-03  8.9453E-02 -1.5278E-04  1.2521E-01  1.6149E-01 -1.3015E-01  1.0005E-01  5.1716E-02
             1.9039E-01
 GRADIENT:  -8.8380E+01 -2.0776E+01 -2.7589E+05 -7.7046E+01 -6.5844E+00 -2.2039E+05 -1.5455E+00  2.1189E+05  2.7579E+05 -2.7595E+05
            -1.4517E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1420
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3377E-02 -1.8984E-02 -2.4459E-02  1.9380E-02 -2.2321E-02
 SE:             2.9567E-02  2.3259E-02  1.7545E-02  2.7762E-02  2.3229E-02
 N:                     100         100         100         100         100

 P VAL.:         1.4235E-01  4.1439E-01  1.6330E-01  4.8513E-01  3.3659E-01

 ETASHRINKSD(%)  9.4665E-01  2.2081E+01  4.1222E+01  6.9942E+00  2.2180E+01
 ETASHRINKVR(%)  1.8843E+00  3.9286E+01  6.5451E+01  1.3499E+01  3.9440E+01
 EBVSHRINKSD(%)  2.9818E-01  2.1824E+01  4.8047E+01  8.1234E+00  1.8430E+01
 EBVSHRINKVR(%)  5.9547E-01  3.8886E+01  7.3009E+01  1.5587E+01  3.3463E+01
 RELATIVEINF(%)  9.9402E+01  3.0578E+01  1.6900E+01  5.8993E+01  2.6413E+01
 EPSSHRINKSD(%)  2.6046E+01
 EPSSHRINKVR(%)  4.5308E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3832.6680857037868     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2178.5787259353760     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.78
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3832.668       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.30E-01  9.69E-01  9.11E-01  9.90E-01  9.05E-01  1.03E+00  1.06E+00  7.94E-01  1.00E+00  9.53E-01  1.09E+00
 


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
+        1.60E+08
 
 TH 2
+       -2.09E+00  1.47E+08
 
 TH 3
+       -8.15E+07  7.81E+07  8.30E+07
 
 TH 4
+       -2.30E+00  2.34E+02 -7.65E+07  1.41E+08
 
 TH 5
+       -8.20E+07  7.87E+07  8.37E+07 -7.71E+07  1.69E+08
 
 TH 6
+        9.31E+01 -1.48E+00 -5.90E+07  4.48E-01 -1.52E+00  4.18E+07
 
 TH 7
+        4.32E+07  3.26E+01 -4.41E+07  1.10E+01 -4.44E+07  3.13E+07  4.68E+07
 
 TH 8
+       -7.18E+07  6.89E+07  7.33E+07 -6.74E+07 -5.62E+00  2.08E+03 -3.89E+07  6.45E+07
 
 TH 9
+       -6.76E+07 -5.18E+00  7.57E+07  5.23E+01  7.19E+00  2.16E+03 -4.02E+07 -6.67E+07  6.89E+07
 
 TH10
+        7.79E+07 -7.47E+07 -7.95E+07  7.32E+07 -8.00E+07 -2.27E+03  3.84E+07  1.47E+01  1.52E+01  7.60E+07
 
 TH11
+       -3.55E+07 -2.69E+01  3.62E+07 -2.73E+00 -8.29E+00 -1.04E+03 -1.92E+07  4.15E+04 -3.30E+07  2.28E+00  3.18E+07
 
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
 #CPUT: Total CPU Time in Seconds,       53.990
Stop Time:
Wed Sep 29 06:07:41 CDT 2021
