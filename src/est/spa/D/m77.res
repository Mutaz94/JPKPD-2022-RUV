Sat Sep 18 15:37:56 CDT 2021
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
$DATA ../../../../data/spa/D/dat77.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   20034.2654539617        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7546E+02  3.2096E+02 -1.1400E+01  2.1423E+02  1.9732E+02 -2.4607E+03 -9.6927E+02 -7.9218E+01 -1.5143E+03 -5.4185E+02
            -3.7707E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -570.680195916634        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4139E+00  1.0668E+00  8.7373E-01  1.5355E+00  1.3573E+00  1.8755E+00  1.0807E+00  9.7969E-01  9.7601E-01  9.9183E-01
             1.4819E+01
 PARAMETER:  4.4632E-01  1.6470E-01 -3.4989E-02  5.2883E-01  4.0547E-01  7.2887E-01  1.7764E-01  7.9483E-02  7.5716E-02  9.1799E-02
             2.7959E+00
 GRADIENT:   9.5221E+00  2.1436E+01  3.6875E+00  2.3909E+01 -8.6888E+00  2.9740E+01 -1.2571E+00  4.3493E+00  4.4215E+00  1.2535E+00
             5.3224E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -585.913944410884        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2817E+00  6.1261E-01  9.2894E-01  1.6584E+00  2.4033E+00  1.6615E+00  2.9652E+00  1.1377E-01  7.7659E-01  1.4691E+00
             1.4206E+01
 PARAMETER:  3.4818E-01 -3.9002E-01  2.6288E-02  6.0588E-01  9.7685E-01  6.0772E-01  1.1869E+00 -2.0736E+00 -1.5284E-01  4.8468E-01
             2.7537E+00
 GRADIENT:  -1.2506E+01  1.5922E+01  5.5674E+00  2.0326E+01 -6.4202E+00  7.3618E-01  3.8871E+00  4.1926E-02  6.7588E+00  1.3825E+00
             3.6354E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -622.889128108554        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.1689E+00  1.4483E-01  4.0948E-01  1.4819E+00  8.2918E+00  1.5713E+00  3.5292E+00  1.0000E-02  1.1675E-01  6.1406E+00
             1.4250E+01
 PARAMETER:  2.5607E-01 -1.8322E+00 -7.9287E-01  4.9330E-01  2.2153E+00  5.5193E-01  1.3611E+00 -6.5815E+00 -2.0477E+00  1.9149E+00
             2.7567E+00
 GRADIENT:   1.5809E+01  3.0634E+00  6.5282E+00  1.5764E+00  1.5028E+01 -1.0154E+01  1.2984E+00  0.0000E+00  7.1193E-01 -1.8409E+00
             2.5591E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -644.501917285369        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0985E+00  1.7657E-01  2.5748E-01  1.3673E+00  5.1343E+00  1.8646E+00  1.1199E+00  1.0000E-02  1.4903E-01  1.6829E+00
             1.4293E+01
 PARAMETER:  1.9399E-01 -1.6341E+00 -1.2568E+00  4.1282E-01  1.7359E+00  7.2306E-01  2.1326E-01 -7.1661E+00 -1.8036E+00  6.2049E-01
             2.7598E+00
 GRADIENT:   2.4835E+00  1.8715E+01  5.1791E+00  3.6452E+01 -2.9582E+01  1.1358E+01  4.3590E+00  0.0000E+00  7.4594E-01  8.5048E+00
             1.4462E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -671.823850275861        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  7.9381E-01  1.2320E-01  1.2111E-01  1.0490E+00  1.2864E+01  1.9394E+00  1.2193E-01  1.0000E-02  7.6926E-02  6.9850E-01
             1.3463E+01
 PARAMETER: -1.3091E-01 -1.9939E+00 -2.0111E+00  1.4784E-01  2.6545E+00  7.6239E-01 -2.0043E+00 -8.8336E+00 -2.4649E+00 -2.5882E-01
             2.6999E+00
 GRADIENT:  -7.5433E+01  3.0902E+01 -4.1511E+01  1.0836E+02 -3.2409E+00 -7.3630E+00  3.8309E-01  0.0000E+00  4.6491E-01  1.6514E-01
             2.8002E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -705.285858655631        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  5.8869E-01  3.5739E-02  3.6922E-02  3.9116E-01  1.7703E+02  1.5444E+00  1.3238E-01  1.0000E-02  1.0000E-02  9.6384E+00
             1.2040E+01
 PARAMETER: -4.2985E-01 -3.2315E+00 -3.1989E+00 -8.3863E-01  5.2763E+00  5.3467E-01 -1.9221E+00 -9.0922E+00 -9.5780E+00  2.3658E+00
             2.5883E+00
 GRADIENT:   2.7884E+01  8.0951E+00  1.2873E+01 -2.5061E+01 -6.4196E-02 -3.1559E+01  1.5320E-01  0.0000E+00  0.0000E+00  2.1253E-02
            -1.9993E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -709.705445975769        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      583
 NPARAMETR:  4.0258E-01  1.3507E-02  1.5323E-02  1.9318E-01  1.0206E+03  1.7107E+00  1.0873E-01  1.0000E-02  1.0000E-02  4.3368E+01
             1.2265E+01
 PARAMETER: -8.0985E-01 -4.2045E+00 -4.0784E+00 -1.5441E+00  7.0281E+00  6.3691E-01 -2.1189E+00 -9.0698E+00 -1.5776E+01  3.8697E+00
             2.6067E+00
 GRADIENT:   1.4331E+00  5.2786E+00  6.5099E+00 -1.5530E+01 -5.4809E-03  9.0266E+00  1.0713E-02  0.0000E+00  0.0000E+00  1.3503E-03
             6.0955E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -710.569367359094        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      759
 NPARAMETR:  4.0627E-01  1.0000E-02  1.5990E-02  2.0030E-01  6.1960E+02  1.6780E+00  2.3261E-01  1.0000E-02  1.0000E-02  4.4426E+01
             1.2182E+01
 PARAMETER: -8.0073E-01 -4.5731E+00 -4.0358E+00 -1.5079E+00  6.5291E+00  6.1758E-01 -1.3584E+00 -9.1472E+00 -1.6073E+01  3.8938E+00
             2.6000E+00
 GRADIENT:   7.9998E-02  0.0000E+00 -4.3675E-01 -1.4639E-01 -1.9994E-04  9.5477E-01  9.0687E-04  0.0000E+00  0.0000E+00  5.1643E-05
             4.3441E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -710.579630620900        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      943
 NPARAMETR:  4.1300E-01  1.0000E-02  1.6665E-02  2.0714E-01  5.4235E+02  1.6755E+00  2.4401E-01  1.0000E-02  1.0000E-02  4.1458E+01
             1.2181E+01
 PARAMETER: -7.8430E-01 -4.5559E+00 -3.9944E+00 -1.4743E+00  6.3959E+00  6.1613E-01 -1.3106E+00 -9.1544E+00 -1.5833E+01  3.8247E+00
             2.5999E+00
 GRADIENT:   5.2626E-03  0.0000E+00  5.0624E-02 -6.9637E-02  3.8391E-05  6.8446E-03  6.6598E-04  0.0000E+00  0.0000E+00  7.5059E-06
             1.9097E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -710.579989883430        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1124
 NPARAMETR:  4.1312E-01  1.0000E-02  1.6662E-02  2.0713E-01  1.2016E+02  1.6755E+00  6.5056E-02  1.0000E-02  1.0000E-02  1.8713E+01
             1.2182E+01
 PARAMETER: -7.8401E-01 -4.5559E+00 -3.9946E+00 -1.4744E+00  4.8888E+00  6.1609E-01 -2.6325E+00 -9.1544E+00 -1.5833E+01  3.0292E+00
             2.6000E+00
 GRADIENT:   2.7409E-02  0.0000E+00  5.9204E-02 -1.1532E-01 -1.4036E-04 -1.3424E-02  5.0627E-05  0.0000E+00  0.0000E+00  7.1558E-05
             8.9173E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -710.580021809064        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1299
 NPARAMETR:  4.1308E-01  1.0000E-02  1.6663E-02  2.0715E-01  1.3882E+02  1.6756E+00  3.8516E-02  1.0000E-02  1.0000E-02  1.0681E+01
             1.2181E+01
 PARAMETER: -7.8412E-01 -4.5559E+00 -3.9946E+00 -1.4743E+00  5.0332E+00  6.1616E-01 -3.1567E+00 -9.1544E+00 -1.5833E+01  2.4685E+00
             2.5998E+00
 GRADIENT:  -9.1891E-03  0.0000E+00 -1.1427E-02 -2.2989E-03 -6.9670E-06  3.4705E-03  1.7467E-05  0.0000E+00  0.0000E+00 -1.6741E-05
            -4.8675E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -710.580031986829        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1479
 NPARAMETR:  4.1307E-01  1.0000E-02  1.6663E-02  2.0715E-01  1.3839E+02  1.6757E+00  1.2659E-02  1.0000E-02  1.0000E-02  1.4649E+01
             1.2180E+01
 PARAMETER: -7.8414E-01 -4.5559E+00 -3.9945E+00 -1.4743E+00  5.0301E+00  6.1625E-01 -4.2694E+00 -9.1544E+00 -1.5833E+01  2.7844E+00
             2.5998E+00
 GRADIENT:  -2.6795E-02  0.0000E+00 -8.7267E-03  1.7512E-03 -1.5931E-05  3.0289E-02  1.8884E-06  0.0000E+00  0.0000E+00 -1.0650E-05
            -2.0362E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -710.580035963975        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1641
 NPARAMETR:  4.1307E-01  1.0000E-02  1.6663E-02  2.0715E-01  1.4787E+02  1.6755E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7209E+01
             1.2181E+01
 PARAMETER: -7.8413E-01 -4.5559E+00 -3.9945E+00 -1.4743E+00  5.0963E+00  6.1614E-01 -4.6496E+00 -9.1544E+00 -1.5833E+01  2.9454E+00
             2.5998E+00
 GRADIENT:  -7.0804E-03  0.0000E+00 -9.9873E-03 -3.8868E-03 -1.0048E-06 -1.7223E-03  0.0000E+00  0.0000E+00  0.0000E+00 -5.2628E-08
            -3.1293E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1641
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5012E-03  4.0977E-06  9.4481E-05 -1.7540E-04  2.0153E-05
 SE:             2.8977E-02  1.6178E-06  2.7975E-04  3.3051E-04  1.1537E-04
 N:                     100         100         100         100         100

 P VAL.:         9.5868E-01  1.1311E-02  7.3557E-01  5.9563E-01  8.6133E-01

 ETASHRINKSD(%)  2.9220E+00  9.9995E+01  9.9063E+01  9.8893E+01  9.9614E+01
 ETASHRINKVR(%)  5.7586E+00  1.0000E+02  9.9991E+01  9.9988E+01  9.9999E+01
 EBVSHRINKSD(%)  3.1942E+00  9.9991E+01  9.8972E+01  9.8780E+01  9.9622E+01
 EBVSHRINKVR(%)  6.2863E+00  1.0000E+02  9.9989E+01  9.9985E+01  9.9999E+01
 RELATIVEINF(%)  7.2052E-01  1.7105E-07  2.0040E-05  2.1142E-05  1.7836E-05
 EPSSHRINKSD(%)  6.4589E+00
 EPSSHRINKVR(%)  1.2501E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -710.58003596397464     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       24.570790599763541     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.72
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -710.580       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.13E-01  1.00E-02  1.67E-02  2.07E-01  1.48E+02  1.68E+00  1.00E-02  1.00E-02  1.00E-02  1.72E+01  1.22E+01
 


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
+        2.14E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.45E+04  0.00E+00  2.78E+06
 
 TH 4
+       -6.67E+02  0.00E+00 -2.59E+05  2.63E+04
 
 TH 5
+        1.84E-03  0.00E+00 -4.85E-02  3.45E-03  1.37E-08
 
 TH 6
+       -6.48E+00  0.00E+00  6.84E+02 -8.23E+01  6.70E-06  5.62E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        4.25E-04  0.00E+00  1.11E-03 -8.61E-04 -1.08E-07  1.53E-05  0.00E+00  0.00E+00  0.00E+00  3.03E-06
 
 TH11
+       -2.01E+01  0.00E+00  4.27E+02 -2.47E+01 -2.26E-05  4.51E-01  0.00E+00  0.00E+00  0.00E+00 -1.31E-07  2.56E+00
 
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
 #CPUT: Total CPU Time in Seconds,       27.519
Stop Time:
Sat Sep 18 15:38:25 CDT 2021
