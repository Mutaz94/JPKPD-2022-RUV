Wed Sep 29 23:55:40 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat5.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m5.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -548.535182477476        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4577E+02  7.5341E+01  1.4608E+02  2.9048E+01  1.1290E+02  2.3987E+01 -7.0134E+01 -8.3456E+01 -6.6466E+00 -1.5114E+02
            -2.6680E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1436.06387602714        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.6155E-01  1.0694E+00  9.2053E-01  1.1517E+00  1.0186E+00  8.1806E-01  1.0168E+00  9.5574E-01  8.3489E-01  1.0507E+00
             5.2860E+00
 PARAMETER:  6.0791E-02  1.6709E-01  1.7197E-02  2.4124E-01  1.1840E-01 -1.0083E-01  1.1668E-01  5.4728E-02 -8.0452E-02  1.4941E-01
             1.7651E+00
 GRADIENT:  -7.8551E+01  3.3086E+01 -2.1089E+01  7.4844E+01 -1.2642E+01 -3.1676E+01  1.0658E+01  6.6654E+00  1.8558E+01  2.1663E+01
             3.9978E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1560.87417747699        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.2536E-01  8.4443E-01  2.0261E+00  1.2388E+00  1.5315E+00  8.4437E-01  1.5699E+00  1.5721E-02  9.0772E-01  7.2570E-01
             3.4424E+00
 PARAMETER:  2.2431E-02 -6.9097E-02  8.0609E-01  3.1414E-01  5.2624E-01 -6.9166E-02  5.5104E-01 -4.0528E+00  3.1828E-03 -2.2061E-01
             1.3362E+00
 GRADIENT:  -3.4382E+01  1.0017E+01 -8.5978E+00  3.5126E+01  2.2564E+01 -3.0166E+01  2.5363E+01  3.3979E-04  3.3765E+01  4.0815E+00
             1.2728E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1580.76482110078        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  9.3105E-01  8.7630E-01  4.5055E-01  1.0843E+00  5.8906E-01  9.3671E-01  1.4458E+00  1.1380E-02  7.1578E-01  2.9952E-01
             2.9317E+00
 PARAMETER:  2.8555E-02 -3.2051E-02 -6.9728E-01  1.8092E-01 -4.2922E-01  3.4614E-02  4.6866E-01 -4.3759E+00 -2.3438E-01 -1.1056E+00
             1.1756E+00
 GRADIENT:  -1.1112E+01  2.4445E+01 -1.0610E+00  3.5117E+01  9.9564E+00  1.7425E-01  1.3640E+01  8.8628E-04 -4.1355E+00  1.9791E+00
             2.0682E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1582.74118228959        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      392
 NPARAMETR:  9.6293E-01  8.3136E-01  3.8333E-01  1.0859E+00  5.2673E-01  9.5792E-01  1.2973E+00  1.0000E-02  7.8269E-01  1.9397E-01
             2.8596E+00
 PARAMETER:  6.2224E-02 -8.4693E-02 -8.5886E-01  1.8244E-01 -5.4107E-01  5.7009E-02  3.6027E-01 -4.7964E+00 -1.4502E-01 -1.5401E+00
             1.1507E+00
 GRADIENT:   3.0192E+01  5.4166E+00 -1.0082E+01  2.5120E+01  1.4495E+01  4.1858E+00 -3.3750E+00  0.0000E+00 -4.5761E+00  2.2806E-01
            -8.7127E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1584.94975878743        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      567
 NPARAMETR:  9.4790E-01  6.0799E-01  3.1851E-01  1.1340E+00  4.0488E-01  9.4106E-01  1.5952E+00  1.0000E-02  8.3453E-01  7.4786E-02
             2.7920E+00
 PARAMETER:  4.6491E-02 -3.9759E-01 -1.0441E+00  2.2579E-01 -8.0416E-01  3.9254E-02  5.6702E-01 -5.3141E+00 -8.0892E-02 -2.4931E+00
             1.1268E+00
 GRADIENT:  -9.3577E-01 -3.2171E+00 -2.2329E+00 -5.5152E+00  8.0253E+00 -2.1815E+00 -6.0835E-01  0.0000E+00 -2.5109E-02 -2.3720E-01
            -4.6452E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1587.20338725228        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  9.3216E-01  5.8237E-01  3.0895E-01  1.1285E+00  3.8144E-01  9.3890E-01  1.5783E+00  1.0000E-02  8.3979E-01  4.7049E-01
             2.7423E+00
 PARAMETER:  2.9748E-02 -4.4065E-01 -1.0746E+00  2.2087E-01 -8.6379E-01  3.6958E-02  5.5634E-01 -5.3603E+00 -7.4605E-02 -6.5399E-01
             1.1088E+00
 GRADIENT:   1.7537E-01  2.1233E+01  1.4796E+01 -3.5964E+00  4.1284E+01 -1.1789E+00  7.9275E+00  0.0000E+00 -1.1495E+01 -1.4184E+00
             3.7218E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1588.21561492915        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  9.0903E-01  5.4974E-01  2.9960E-01  1.1147E+00  3.5715E-01  9.3563E-01  1.4564E+00  1.0000E-02  8.9510E-01  6.0307E-01
             2.5799E+00
 PARAMETER:  4.6206E-03 -4.9831E-01 -1.1053E+00  2.0855E-01 -9.2959E-01  3.3465E-02  4.7599E-01 -5.3603E+00 -1.0821E-02 -4.0572E-01
             1.0478E+00
 GRADIENT:  -5.0373E+01  3.4638E+01  5.5518E+01 -3.4106E+01  6.2966E+00 -5.9122E+00  1.7595E+00  0.0000E+00 -1.4349E+01  1.0635E+00
             1.0540E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1588.30255780083        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      873
 NPARAMETR:  9.0054E-01  5.1719E-01  2.8506E-01  1.1093E+00  3.3709E-01  9.3811E-01  1.3996E+00  1.0000E-02  9.3493E-01  6.0924E-01
             2.5051E+00
 PARAMETER: -4.7639E-03 -5.5935E-01 -1.1551E+00  2.0373E-01 -9.8741E-01  3.6110E-02  4.3618E-01 -5.3603E+00  3.2722E-02 -3.9554E-01
             1.0183E+00
 GRADIENT:  -6.7682E+01  3.6722E+01  7.2817E+01 -4.0888E+01 -1.6114E-01 -7.0823E+00 -1.4637E+00  0.0000E+00 -1.1082E+01 -1.0105E+00
            -9.0754E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1592.51350170157        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  9.3609E-01  4.8093E-01  2.8591E-01  1.1527E+00  3.3667E-01  9.5441E-01  1.3242E+00  1.0000E-02  1.0176E+00  6.3120E-01
             2.4368E+00
 PARAMETER:  3.3959E-02 -6.3203E-01 -1.1521E+00  2.4209E-01 -9.8866E-01  5.3337E-02  3.8082E-01 -5.3603E+00  1.1744E-01 -3.6013E-01
             9.9067E-01
 GRADIENT:  -1.9771E+01  3.5292E+00  2.3272E+01 -8.5413E+00  4.0016E+00 -7.9880E-01 -6.7506E+00  0.0000E+00  8.4918E+00 -3.5304E+00
            -4.3639E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1593.32141092174        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1229
 NPARAMETR:  9.3339E-01  4.2693E-01  3.0234E-01  1.1813E+00  3.3538E-01  9.4846E-01  1.5904E+00  1.0000E-02  9.6668E-01  6.1150E-01
             2.4816E+00
 PARAMETER:  3.1063E-02 -7.5115E-01 -1.0962E+00  2.6661E-01 -9.9249E-01  4.7089E-02  5.6399E-01 -5.3603E+00  6.6110E-02 -3.9184E-01
             1.0089E+00
 GRADIENT:  -2.5886E+01  1.2088E+00  2.9857E+01 -8.4706E+00 -6.6497E+00 -1.6193E+00 -4.1930E+00  0.0000E+00  4.8201E+00 -4.5048E+00
            -3.2665E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1597.64049281479        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1407
 NPARAMETR:  9.4599E-01  3.6211E-01  2.0971E-01  1.1158E+00  2.5901E-01  9.5687E-01  1.4578E+00  1.0000E-02  1.0602E+00  6.1404E-01
             2.4987E+00
 PARAMETER:  4.4475E-02 -9.1582E-01 -1.4620E+00  2.0961E-01 -1.2509E+00  5.5914E-02  4.7692E-01 -5.3603E+00  1.5850E-01 -3.8769E-01
             1.0158E+00
 GRADIENT:   1.6601E+00 -3.8666E+00 -7.0103E+00 -2.6717E+00  1.0933E+01 -5.2673E-01  3.7250E-01  0.0000E+00  6.3162E-02 -2.4535E+00
             1.7874E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1597.76751527214        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1582
 NPARAMETR:  9.4570E-01  3.6689E-01  2.0596E-01  1.1114E+00  2.5574E-01  9.5889E-01  1.3302E+00  1.0000E-02  1.0716E+00  6.5020E-01
             2.4860E+00
 PARAMETER:  4.4171E-02 -9.0268E-01 -1.4801E+00  2.0558E-01 -1.2636E+00  5.8025E-02  3.8531E-01 -5.3603E+00  1.6912E-01 -3.3048E-01
             1.0107E+00
 GRADIENT:   4.2620E-02 -1.0580E-01 -2.0509E-01 -7.1454E-02  2.7695E-01  2.3790E-02 -2.5335E-02  0.0000E+00 -2.0689E-01 -7.9458E-03
            -1.1166E-01

0ITERATION NO.:   61    OBJECTIVE VALUE:  -1597.76751527214        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1604
 NPARAMETR:  9.4570E-01  3.6689E-01  2.0596E-01  1.1114E+00  2.5574E-01  9.5889E-01  1.3302E+00  1.0000E-02  1.0716E+00  6.5020E-01
             2.4860E+00
 PARAMETER:  4.4171E-02 -9.0268E-01 -1.4801E+00  2.0558E-01 -1.2636E+00  5.8025E-02  3.8531E-01 -5.3603E+00  1.6912E-01 -3.3048E-01
             1.0107E+00
 GRADIENT:   4.2620E-02 -1.0580E-01 -2.0509E-01 -7.1454E-02  2.7695E-01  2.3790E-02 -2.5335E-02  0.0000E+00 -2.0689E-01 -7.9458E-03
            -1.1166E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1604
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2520E-03  1.8276E-02 -6.7545E-05 -1.1986E-02  6.9290E-03
 SE:             2.9144E-02  1.5554E-02  2.2011E-04  2.6863E-02  2.1649E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6573E-01  2.4001E-01  7.5894E-01  6.5546E-01  7.4893E-01

 ETASHRINKSD(%)  2.3630E+00  4.7891E+01  9.9263E+01  1.0006E+01  2.7473E+01
 ETASHRINKVR(%)  4.6702E+00  7.2846E+01  9.9995E+01  1.9011E+01  4.7398E+01
 EBVSHRINKSD(%)  2.2107E+00  4.9055E+01  9.9282E+01  8.8262E+00  2.6872E+01
 EBVSHRINKVR(%)  4.3726E+00  7.4046E+01  9.9995E+01  1.6873E+01  4.6523E+01
 RELATIVEINF(%)  9.5543E+01  5.1756E+00  3.7908E-04  4.5843E+01  2.3230E+00
 EPSSHRINKSD(%)  2.8770E+01
 EPSSHRINKVR(%)  4.9262E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1597.7675152721390     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -678.82898206746631     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.66
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.73
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1597.768       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.46E-01  3.67E-01  2.06E-01  1.11E+00  2.56E-01  9.59E-01  1.33E+00  1.00E-02  1.07E+00  6.50E-01  2.49E+00
 


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
+        1.30E+03
 
 TH 2
+       -3.27E+01  1.22E+03
 
 TH 3
+       -3.93E+01  2.18E+03  1.35E+04
 
 TH 4
+       -1.35E+01  1.25E+02 -7.17E+02  6.54E+02
 
 TH 5
+        1.28E+02 -4.06E+03 -1.53E+04 -1.63E+02  2.18E+04
 
 TH 6
+        5.62E+00 -1.23E+01  4.69E+01 -1.19E+01 -1.22E+01  1.96E+02
 
 TH 7
+        2.00E+00  4.46E+01 -4.12E+01 -1.77E+00 -4.69E+01 -1.35E-01  1.04E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.09E+01 -1.95E+01  2.02E+02 -1.18E+01  2.05E+01 -1.65E+00  3.71E+00  0.00E+00  1.14E+02
 
 TH10
+       -5.97E+00  3.89E+01 -1.37E+02  6.30E+00  1.52E+02  2.54E+00  2.24E+01  0.00E+00  1.76E+00  1.34E+02
 
 TH11
+       -1.92E+01 -1.63E-01 -1.05E+02 -7.36E+00  7.31E+01  3.63E+00  3.17E+00  0.00E+00  8.41E+00  1.49E+01  7.03E+01
 
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
 #CPUT: Total CPU Time in Seconds,       31.456
Stop Time:
Wed Sep 29 23:56:32 CDT 2021
