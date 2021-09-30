Wed Sep 29 12:07:37 CDT 2021
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
$DATA ../../../../data/spa/A1/dat42.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m42.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1346.23188919429        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0377E+02  8.3472E+01 -2.8155E+00  1.3319E+02  9.6921E+01  5.2144E+01 -1.6007E+01  2.7935E+00 -2.6961E+01 -3.4242E+01
            -6.1839E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1520.67024110231        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0797E+00  9.1340E-01  1.0305E+00  1.0431E+00  9.0721E-01  8.6559E-01  9.9498E-01  9.0260E-01  1.0928E+00  8.8751E-01
             2.0854E+00
 PARAMETER:  1.7667E-01  9.4198E-03  1.3005E-01  1.4219E-01  2.6175E-03 -4.4341E-02  9.4972E-02 -2.4738E-03  1.8873E-01 -1.9332E-02
             8.3495E-01
 GRADIENT:   2.7112E+02  2.7381E+01  7.7850E+00  3.9516E+01 -1.1491E+01 -3.9589E+01  4.2003E+00  6.6115E+00  9.5406E+00  6.2991E+00
             2.5498E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1527.19555430188        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0712E+00  8.5823E-01  6.7389E-01  1.0652E+00  7.2489E-01  9.0146E-01  9.6930E-01  4.6461E-01  1.0427E+00  5.7990E-01
             2.0527E+00
 PARAMETER:  1.6876E-01 -5.2886E-02 -2.9469E-01  1.6312E-01 -2.2174E-01 -3.7444E-03  6.8821E-02 -6.6656E-01  1.4180E-01 -4.4489E-01
             8.1915E-01
 GRADIENT:   2.2325E+02  1.8607E+01 -2.2277E+01  8.2248E+01  5.2789E+01 -2.3270E+01 -4.0834E+00  2.4095E+00  2.2735E+00 -3.4331E+00
             1.1306E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1531.54749124281        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      272
 NPARAMETR:  1.0398E+00  7.4368E-01  5.5321E-01  1.0883E+00  6.1419E-01  9.3292E-01  1.1705E+00  1.7095E-01  1.0109E+00  6.6174E-01
             1.9126E+00
 PARAMETER:  1.3904E-01 -1.9614E-01 -4.9201E-01  1.8462E-01 -3.8746E-01  3.0563E-02  2.5740E-01 -1.6664E+00  1.1082E-01 -3.1289E-01
             7.4844E-01
 GRADIENT:  -4.5666E+01 -7.1130E+00 -4.5476E+01  2.6577E+01  6.3087E+01 -2.1145E+01 -3.6055E-01  6.7701E-01  7.1666E+00  8.8693E+00
             3.2641E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1539.08171836771        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  1.0555E+00  5.0779E-01  4.7786E-01  1.1725E+00  4.7612E-01  9.9728E-01  1.6573E+00  5.2483E-02  8.9477E-01  5.5393E-01
             1.8730E+00
 PARAMETER:  1.5398E-01 -5.7768E-01 -6.3844E-01  2.5914E-01 -6.4210E-01  9.7272E-02  6.0516E-01 -2.8473E+00 -1.1186E-02 -4.9071E-01
             7.2755E-01
 GRADIENT:  -7.6375E+00  1.1449E+01  3.3740E+00 -2.3370E+00 -7.0376E+00  4.8897E+00  4.0977E+00  6.1089E-02 -9.2074E-03  3.0986E+00
             3.1940E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1540.67215156778        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      623
 NPARAMETR:  1.0602E+00  3.5437E-01  4.9526E-01  1.2490E+00  4.4904E-01  9.8268E-01  1.8280E+00  2.4544E-02  8.8410E-01  5.9101E-01
             1.8587E+00
 PARAMETER:  1.5843E-01 -9.3743E-01 -6.0268E-01  3.2234E-01 -7.0065E-01  8.2524E-02  7.0323E-01 -3.6073E+00 -2.3186E-02 -4.2592E-01
             7.1985E-01
 GRADIENT:   8.4168E+00  1.8374E+00  9.4211E+00 -4.5007E+00 -1.4861E+01  4.3837E-01 -4.8939E+00  1.0667E-02 -1.1213E+00 -1.7160E+00
            -1.6137E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1541.92629678165        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      799
 NPARAMETR:  1.0482E+00  2.1119E-01  5.7178E-01  1.3477E+00  4.6937E-01  9.7405E-01  2.7353E+00  1.3999E-02  8.6860E-01  6.6241E-01
             1.8606E+00
 PARAMETER:  1.4712E-01 -1.4550E+00 -4.5899E-01  3.9840E-01 -6.5636E-01  7.3711E-02  1.1062E+00 -4.1687E+00 -4.0870E-02 -3.1187E-01
             7.2088E-01
 GRADIENT:  -2.1510E+00  2.9564E+00  9.0140E+00  3.5235E+00 -1.2528E+01 -3.9458E-01  1.3295E+00  3.6895E-03  1.1642E+00  5.3126E-01
             6.2080E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1542.19840250451        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      976
 NPARAMETR:  1.0473E+00  1.5763E-01  5.4362E-01  1.3616E+00  4.4829E-01  9.7512E-01  3.2396E+00  1.0000E-02  8.5256E-01  6.4154E-01
             1.8468E+00
 PARAMETER:  1.4622E-01 -1.7475E+00 -5.0950E-01  4.0869E-01 -7.0232E-01  7.4807E-02  1.2755E+00 -4.7695E+00 -5.9510E-02 -3.4388E-01
             7.1345E-01
 GRADIENT:  -6.3095E-01 -1.3606E+00  2.3256E-01  1.2328E+01 -1.2945E+00  1.4705E-01 -4.2652E+00  0.0000E+00 -1.2152E+00  1.0048E+00
            -4.0005E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1542.52583974001        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1152
 NPARAMETR:  1.0475E+00  1.1304E-01  5.0125E-01  1.3595E+00  4.1784E-01  9.7587E-01  3.9762E+00  1.0000E-02  8.5051E-01  6.1420E-01
             1.8305E+00
 PARAMETER:  1.4637E-01 -2.0800E+00 -5.9064E-01  4.0710E-01 -7.7265E-01  7.5574E-02  1.4803E+00 -5.5556E+00 -6.1917E-02 -3.8743E-01
             7.0461E-01
 GRADIENT:   6.8979E+00  2.0100E+01 -1.7505E+01 -1.3815E+01  2.0400E+01  1.0080E+00  3.0191E+01  0.0000E+00 -1.9031E+01 -9.1620E+00
            -9.2615E+00

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1542.59791637883        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1286
 NPARAMETR:  1.0458E+00  1.0675E-01  5.0124E-01  1.3599E+00  4.1589E-01  9.7433E-01  4.0944E+00  1.0000E-02  8.5256E-01  6.1653E-01
             1.8394E+00
 PARAMETER:  1.4533E-01 -2.1308E+00 -5.9260E-01  4.0613E-01 -7.7527E-01  7.4995E-02  1.5143E+00 -5.6643E+00 -5.9866E-02 -3.8501E-01
             7.0694E-01
 GRADIENT:   4.1133E+02  2.3521E+01 -9.6157E+01 -1.3443E+02  6.6129E+01  3.5664E-01  6.3717E+01  0.0000E+00 -2.9974E+02 -1.5384E+02
            -8.4378E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1286
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.3984E-04  3.4301E-02 -1.6657E-04 -2.0279E-02  4.6952E-03
 SE:             2.9474E-02  1.4273E-02  2.3722E-04  2.7139E-02  2.0940E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7997E-01  1.6251E-02  4.8256E-01  4.5492E-01  8.2259E-01

 ETASHRINKSD(%)  1.2590E+00  5.2184E+01  9.9205E+01  9.0820E+00  2.9848E+01
 ETASHRINKVR(%)  2.5021E+00  7.7136E+01  9.9994E+01  1.7339E+01  5.0787E+01
 EBVSHRINKSD(%)  1.3839E+00  6.5344E+01  9.9112E+01  7.3181E+00  2.6046E+01
 EBVSHRINKVR(%)  2.7486E+00  8.7989E+01  9.9992E+01  1.4101E+01  4.5308E+01
 RELATIVEINF(%)  9.6199E+01  5.4890E+00  3.5912E-04  3.5642E+01  2.5094E+00
 EPSSHRINKSD(%)  3.7808E+01
 EPSSHRINKVR(%)  6.1321E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1542.5979163788295     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -807.44708981509132     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.54
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.32
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1542.598       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.07E-01  5.00E-01  1.36E+00  4.17E-01  9.75E-01  4.11E+00  1.00E-02  8.52E-01  6.16E-01  1.83E+00
 


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
+        6.53E+04
 
 TH 2
+       -4.79E+00  5.17E+04
 
 TH 3
+        4.11E+01 -1.80E+04  2.02E+04
 
 TH 4
+       -2.83E+00 -9.28E+03 -5.15E+02  5.51E+03
 
 TH 5
+        3.31E+01  1.61E+04 -5.42E+03  2.55E+01  2.34E+04
 
 TH 6
+       -6.26E+01 -5.40E+01  2.63E+01  2.99E+00 -1.75E+01  2.81E+02
 
 TH 7
+       -3.28E-01  9.23E+02 -7.07E+02 -3.67E+02  6.58E+02 -9.58E-01  7.56E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.40E+02 -4.92E+02 -1.57E+02 -1.27E+02  9.30E+01  7.53E+01 -2.82E+00  0.00E+00  2.04E+05
 
 TH10
+        8.63E+01 -2.63E+02 -1.65E+02 -5.53E+01  5.81E+01  2.27E+01 -8.75E+00  0.00E+00 -4.84E+02  2.64E+04
 
 TH11
+        3.81E+00 -4.19E+01 -6.53E+01 -2.70E+01  3.57E+01  5.39E+00 -5.53E-01  0.00E+00 -6.83E+01 -9.09E+00  9.48E+02
 
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
 #CPUT: Total CPU Time in Seconds,       23.924
Stop Time:
Wed Sep 29 12:08:03 CDT 2021
