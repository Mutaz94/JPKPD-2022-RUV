Sat Sep 25 12:20:06 CDT 2021
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
$DATA ../../../../data/spa/S2/dat50.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1697.60371914728        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1009E+01 -3.7426E+01 -3.6298E+01 -4.3644E+01  1.4737E-01  3.1320E+00  1.4917E+01  1.8704E+01 -2.7508E+00  3.8335E+01
             8.8430E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1706.02582821089        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0141E+00  9.4320E-01  1.2332E+00  1.0496E+00  1.0531E+00  9.8275E-01  7.6205E-01  8.0144E-01  1.0800E+00  7.1910E-01
             1.0727E+00
 PARAMETER:  1.1400E-01  4.1527E-02  3.0958E-01  1.4844E-01  1.5170E-01  8.2603E-02 -1.7174E-01 -1.2134E-01  1.7700E-01 -2.2976E-01
             1.7014E-01
 GRADIENT:   8.6467E+01 -1.6931E+01  2.0217E+01 -2.6382E+01  1.6442E+01 -4.5380E+00  5.0966E+00 -8.9866E-01  1.1092E+01 -1.9504E+01
             1.7634E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1709.58474882563        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0089E+00  9.5374E-01  1.1521E+00  1.0457E+00  1.0622E+00  1.0166E+00  4.2060E-01  4.2033E-01  1.1039E+00  9.5603E-01
             9.8504E-01
 PARAMETER:  1.0886E-01  5.2637E-02  2.4155E-01  1.4470E-01  1.6033E-01  1.1645E-01 -7.6607E-01 -7.6673E-01  1.9885E-01  5.5030E-02
             8.4924E-02
 GRADIENT:   7.8164E+01 -1.9176E+01 -1.1966E+01 -1.0083E+01  3.4304E+01  1.0492E+01 -1.2014E+00  7.6734E-01  1.9374E+00  5.7137E+00
            -1.3358E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1711.48185093673        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      252
 NPARAMETR:  9.8302E-01  8.7927E-01  1.0714E+00  1.1012E+00  9.7534E-01  9.9129E-01  6.1547E-01  3.5184E-01  1.0252E+00  8.7076E-01
             1.0154E+00
 PARAMETER:  8.2876E-02 -2.8667E-02  1.6900E-01  1.9637E-01  7.5028E-02  9.1251E-02 -3.8537E-01 -9.4457E-01  1.2484E-01 -3.8384E-02
             1.1531E-01
 GRADIENT:  -2.6439E+01 -8.6335E+00 -1.1605E+01 -6.7924E+00  1.6157E+01 -4.2475E+00 -3.3763E-01  7.2039E-01 -1.0326E+00  4.4705E+00
             2.8157E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1712.22022005752        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  9.9124E-01  6.9410E-01  1.0620E+00  1.2187E+00  8.8941E-01  1.0004E+00  7.3781E-01  1.7047E-01  9.3119E-01  8.3829E-01
             1.0117E+00
 PARAMETER:  9.1205E-02 -2.6514E-01  1.6011E-01  2.9782E-01 -1.7200E-02  1.0040E-01 -2.0406E-01 -1.6692E+00  2.8705E-02 -7.6387E-02
             1.1162E-01
 GRADIENT:  -4.4419E+00  6.5032E-01 -1.1523E+00 -1.6336E+00 -1.6895E+00  2.8938E-01 -7.5910E-02  1.7562E-01 -1.2956E+00  1.9283E+00
             6.0128E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1712.39972230652        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  9.9047E-01  5.0466E-01  1.0656E+00  1.3329E+00  8.2785E-01  9.9721E-01  8.5979E-01  3.0534E-02  8.5961E-01  8.2849E-01
             1.0073E+00
 PARAMETER:  9.0421E-02 -5.8388E-01  1.6352E-01  3.8733E-01 -8.8920E-02  9.7203E-02 -5.1067E-02 -3.3889E+00 -5.1276E-02 -8.8152E-02
             1.0726E-01
 GRADIENT:  -9.9821E-01  5.6494E-01  1.3761E+00 -5.7419E-01 -2.3642E+00  1.3380E-02 -5.6550E-02  4.9467E-03 -7.3349E-01 -3.2173E-02
            -6.3413E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1712.41009745993        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      780
 NPARAMETR:  9.9001E-01  4.3901E-01  1.0641E+00  1.3721E+00  8.0802E-01  9.9606E-01  9.3353E-01  1.1774E-02  8.3815E-01  8.2692E-01
             1.0069E+00
 PARAMETER:  8.9964E-02 -7.2324E-01  1.6217E-01  4.1636E-01 -1.1317E-01  9.6054E-02  3.1217E-02 -4.3418E+00 -7.6560E-02 -9.0052E-02
             1.0692E-01
 GRADIENT:   1.7572E-01  3.9352E-01  7.1808E-01  1.0943E+00 -9.8861E-01 -4.2703E-02  2.4734E-02  7.3535E-04 -1.0710E-01  2.6738E-02
            -1.0395E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1712.41077588250        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      960             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8991E-01  4.3523E-01  1.0629E+00  1.3736E+00  8.0661E-01  9.9613E-01  9.3223E-01  1.0000E-02  8.3740E-01  8.2650E-01
             1.0070E+00
 PARAMETER:  8.9860E-02 -7.3188E-01  1.6099E-01  4.1744E-01 -1.1491E-01  9.6126E-02  2.9824E-02 -4.6666E+00 -7.7458E-02 -9.0554E-02
             1.0694E-01
 GRADIENT:   4.0038E+01  4.9343E+00  4.6545E-01  5.0308E+01  8.9566E-01  3.9948E+00  5.6319E-02  0.0000E+00  7.5701E-01  4.2957E-02
             9.5111E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1712.41079432008        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1136
 NPARAMETR:  9.8988E-01  4.3515E-01  1.0630E+00  1.3737E+00  8.0663E-01  9.9611E-01  9.2939E-01  1.0000E-02  8.3748E-01  8.2670E-01
             1.0069E+00
 PARAMETER:  8.9824E-02 -7.3207E-01  1.6106E-01  4.1751E-01 -1.1488E-01  9.6103E-02  2.6769E-02 -4.6666E+00 -7.7353E-02 -9.0312E-02
             1.0688E-01
 GRADIENT:  -6.8478E-03 -5.3577E-03 -5.5267E-03 -5.5346E-02 -2.3713E-02 -5.9639E-04  1.4090E-03  0.0000E+00  5.6177E-03  1.6048E-03
            -9.1135E-03

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1712.41079467674        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     1199
 NPARAMETR:  9.8990E-01  4.3512E-01  1.0630E+00  1.3737E+00  8.0667E-01  9.9612E-01  9.2875E-01  1.0000E-02  8.3748E-01  8.2674E-01
             1.0069E+00
 PARAMETER:  8.9827E-02 -7.3212E-01  1.6110E-01  4.1754E-01 -1.1487E-01  9.6104E-02  2.5985E-02 -4.6666E+00 -7.7346E-02 -9.0277E-02
             1.0690E-01
 GRADIENT:  -1.4405E-02  8.8014E-04 -1.4096E-02  4.7060E-02 -1.7482E-02 -1.3464E-03 -6.4483E-04  0.0000E+00  1.3730E-03 -5.2611E-04
            -2.2920E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1199
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.4846E-04 -9.4003E-03 -2.3006E-04 -2.7425E-03 -1.9196E-02
 SE:             2.9829E-02  7.6201E-03  2.1549E-04  2.8521E-02  2.5147E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9603E-01  2.1734E-01  2.8571E-01  9.2340E-01  4.4525E-01

 ETASHRINKSD(%)  6.9541E-02  7.4472E+01  9.9278E+01  4.4500E+00  1.5754E+01
 ETASHRINKVR(%)  1.3903E-01  9.3483E+01  9.9995E+01  8.7020E+00  2.9025E+01
 EBVSHRINKSD(%)  4.2456E-01  7.5017E+01  9.9270E+01  4.4131E+00  1.4291E+01
 EBVSHRINKVR(%)  8.4732E-01  9.3758E+01  9.9995E+01  8.6315E+00  2.6540E+01
 RELATIVEINF(%)  9.7139E+01  2.0257E-01  6.4285E-04  4.7293E+00  3.3806E+00
 EPSSHRINKSD(%)  4.1894E+01
 EPSSHRINKVR(%)  6.6236E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1712.4107946767440     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -977.25996811300581     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.91
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1712.411       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.90E-01  4.35E-01  1.06E+00  1.37E+00  8.07E-01  9.96E-01  9.29E-01  1.00E-02  8.37E-01  8.27E-01  1.01E+00
 


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
+        1.14E+03
 
 TH 2
+       -2.44E+01  4.36E+02
 
 TH 3
+        1.06E+01  2.26E+02  5.36E+02
 
 TH 4
+       -9.30E+00  4.45E+02 -5.79E+01  7.67E+02
 
 TH 5
+        1.84E+00 -5.36E+02 -9.10E+02 -1.05E+01  1.85E+03
 
 TH 6
+        4.90E-01 -2.69E+00  5.80E-01 -2.87E+00 -1.44E+00  1.99E+02
 
 TH 7
+        1.99E-01 -2.49E+00  5.49E-01 -3.22E+00 -1.66E+00  1.16E+00  2.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.80E+00 -6.34E+01  6.09E+00  7.83E+00 -1.09E+01 -7.15E-01  1.43E+01  0.00E+00  2.40E+02
 
 TH10
+       -2.08E+00  2.34E+01 -2.28E+01 -4.63E+00 -8.16E+01 -3.26E-01  8.38E+00  0.00E+00 -6.06E+00  1.58E+02
 
 TH11
+       -8.59E+00 -1.75E+01 -4.16E+01 -8.90E+00  1.94E+01 -1.22E+00  3.58E+00  0.00E+00  9.12E+00  3.72E+01  2.29E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.182
Stop Time:
Sat Sep 25 12:20:26 CDT 2021
