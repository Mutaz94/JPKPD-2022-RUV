Thu Sep 30 03:10:07 CDT 2021
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
$DATA ../../../../data/spa1/D/dat50.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:   14575.1187588498        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2020E+02  2.4545E+02 -6.0033E+01  1.0706E+02  2.5699E+02 -1.3612E+03 -5.2558E+02 -3.5307E+01 -1.0736E+03 -5.8315E+02
            -2.9508E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -696.336372543862        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2034E+00  1.0223E+00  9.9788E-01  1.7167E+00  1.2297E+00  2.0368E+00  1.2852E+00  9.4214E-01  1.6398E+00  1.2353E+00
             1.3877E+01
 PARAMETER:  2.8515E-01  1.2202E-01  9.7882E-02  6.4042E-01  3.0674E-01  8.1139E-01  3.5092E-01  4.0398E-02  5.9456E-01  3.1131E-01
             2.7302E+00
 GRADIENT:  -6.9887E+01  1.9690E+01 -9.5848E+00  2.6560E+01 -9.1231E+00  3.8321E+01  1.1705E+00  5.3954E+00  1.4645E+00  4.2883E+00
             2.8466E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -728.141103549654        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.2208E+00  6.8251E-01  2.3816E+00  2.4773E+00  5.3319E+00  2.0467E+00  6.4307E+00  4.5347E-01  2.6476E+00  8.3421E+00
             1.1479E+01
 PARAMETER:  2.9952E-01 -2.8198E-01  9.6779E-01  1.0072E+00  1.7737E+00  8.1624E-01  1.9611E+00 -6.9082E-01  1.0737E+00  2.2213E+00
             2.5405E+00
 GRADIENT:  -3.3658E+01  1.3697E+01  1.0170E+01  7.2700E+01 -1.4828E+00 -1.2124E+01  1.2219E+01 -3.3791E-01  6.0781E+01  5.4719E+00
             2.3285E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -795.152824135796        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.1778E+00  4.9815E-01  7.5668E-01  1.5212E+00  4.3083E+00  1.9393E+00  2.4185E+00  2.1724E-02  1.3804E+00  5.4809E+00
             9.4073E+00
 PARAMETER:  2.6367E-01 -5.9684E-01 -1.7881E-01  5.1947E-01  1.5605E+00  7.6231E-01  9.8315E-01 -3.7294E+00  4.2235E-01  1.8013E+00
             2.3415E+00
 GRADIENT:   2.9611E+01  3.0006E+01  2.4942E+00 -1.8759E+01 -2.2408E+01  2.3921E+01  6.3632E+00 -4.9568E-03 -5.0667E+01  2.0868E+00
             8.8980E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -849.989781306259        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  9.4232E-01  1.4987E-01  2.0253E-01  1.2041E+00  1.9257E+01  1.5381E+00  4.0446E-01  1.0000E-02  1.3393E+00  6.8558E+00
             7.3212E+00
 PARAMETER:  4.0595E-02 -1.7980E+00 -1.4968E+00  2.8575E-01  3.0579E+00  5.3057E-01 -8.0521E-01 -5.5743E+00  3.9217E-01  2.0251E+00
             2.0908E+00
 GRADIENT:   7.8104E+01  4.4211E+01 -5.3778E+01  1.2812E+02 -3.9564E+00 -5.6768E+01  1.3676E+00  0.0000E+00  2.7544E+00  1.8506E+00
            -1.5687E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -886.892366600829        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      396
 NPARAMETR:  5.5914E-01  3.5829E-02  5.0982E-02  5.8413E-01  4.3141E+01  1.3693E+00  1.0000E-02  1.0000E-02  1.0201E+00  5.5185E+00
             8.0584E+00
 PARAMETER: -4.8135E-01 -3.2290E+00 -2.8763E+00 -4.3763E-01  3.8645E+00  4.1430E-01 -4.9379E+00 -9.3634E+00  1.1986E-01  1.8081E+00
             2.1867E+00
 GRADIENT:  -5.2921E+01  6.1527E+00 -1.4867E+02  2.5651E+02 -5.1791E-02 -5.6137E+01  0.0000E+00  0.0000E+00 -2.4572E+01  5.3710E-02
            -3.4387E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -920.502497578925        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      572
 NPARAMETR:  4.9411E-01  1.6267E-02  3.2699E-02  3.5499E-01  5.7879E+01  1.5077E+00  1.0000E-02  1.0000E-02  8.9119E-01  4.5154E+00
             8.1522E+00
 PARAMETER: -6.0500E-01 -4.0186E+00 -3.3204E+00 -9.3567E-01  4.1583E+00  5.1056E-01 -7.4920E+00 -1.0697E+01 -1.5192E-02  1.6075E+00
             2.1983E+00
 GRADIENT:  -5.1042E+00  4.8486E-01  2.9691E+00  1.4534E+00  1.9289E-02  1.1863E+00  0.0000E+00  0.0000E+00 -2.3648E+00 -7.5183E-06
             1.0028E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -921.663413833343        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      758             RESET HESSIAN, TYPE I
 NPARAMETR:  4.5340E-01  1.1807E-02  2.5764E-02  2.9589E-01  6.4647E+01  1.4768E+00  1.0000E-02  1.0000E-02  8.9541E-01  4.6342E+00
             8.1243E+00
 PARAMETER: -6.9098E-01 -4.3391E+00 -3.5588E+00 -1.1178E+00  4.2689E+00  4.8987E-01 -8.6154E+00 -1.1119E+01 -1.0472E-02  1.6335E+00
             2.1949E+00
 GRADIENT:   7.8724E+01  1.4753E-01  1.0740E+02  3.5647E+01  2.0196E-02  1.0307E+01  0.0000E+00  0.0000E+00  7.1042E-01  1.8963E-04
             2.0878E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -921.676038992503        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:      907
 NPARAMETR:  4.5029E-01  1.0812E-02  2.5490E-02  2.9508E-01  5.8821E+01  1.4749E+00  1.0000E-02  1.0000E-02  8.9290E-01  4.7634E+00
             8.1188E+00
 PARAMETER: -6.9787E-01 -4.4271E+00 -3.5695E+00 -1.1205E+00  4.1745E+00  4.8858E-01 -8.6154E+00 -1.1119E+01 -1.3283E-02  1.6610E+00
             2.1942E+00
 GRADIENT:  -1.3003E+00  1.6518E-02 -7.5224E+00  8.3125E+00  2.4261E-02 -6.0978E-01  0.0000E+00  0.0000E+00 -4.2243E-01 -5.7915E-05
            -2.0473E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -921.736552366717        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1066             RESET HESSIAN, TYPE I
 NPARAMETR:  4.5051E-01  1.0000E-02  2.4974E-02  2.9015E-01  4.3843E+01  1.4769E+00  1.0000E-02  1.0000E-02  8.8824E-01  6.3151E+00
             8.1174E+00
 PARAMETER: -6.9737E-01 -4.5280E+00 -3.5899E+00 -1.1374E+00  3.8806E+00  4.8994E-01 -8.6154E+00 -1.1119E+01 -1.8508E-02  1.9430E+00
             2.1940E+00
 GRADIENT:   8.5937E+01 -2.4773E-02  1.0138E+02  4.4188E+01  4.3199E-02  1.0805E+01  0.0000E+00  0.0000E+00 -8.7396E-01  9.3470E-04
             1.8243E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -921.763158816564        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1203
 NPARAMETR:  4.4798E-01  1.0000E-02  2.4962E-02  2.8992E-01  4.1486E+01  1.4750E+00  1.0000E-02  1.0000E-02  8.8992E-01  6.3918E+00
             8.1194E+00
 PARAMETER: -7.0301E-01 -4.5070E+00 -3.5904E+00 -1.1381E+00  3.8254E+00  4.8864E-01 -8.6154E+00 -1.1119E+01 -1.6628E-02  1.9550E+00
             2.1943E+00
 GRADIENT:   1.6636E+00  2.0429E-03 -6.7972E+00  5.4538E+00  3.6769E-02 -1.0519E-01  0.0000E+00  0.0000E+00 -7.6772E-01 -8.8998E-05
            -2.3899E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -921.775004611320        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1389
 NPARAMETR:  4.4745E-01  1.0000E-02  2.5009E-02  2.8949E-01  3.4870E+01  1.4749E+00  1.0000E-02  1.0000E-02  8.9096E-01  9.0903E+00
             8.1311E+00
 PARAMETER: -7.0418E-01 -4.5118E+00 -3.5885E+00 -1.1396E+00  3.6516E+00  4.8861E-01 -8.6154E+00 -1.1119E+01 -1.5454E-02  2.3072E+00
             2.1957E+00
 GRADIENT:  -5.7733E-01  5.2176E-03 -1.4534E+00 -3.1274E-01  3.6873E-02 -7.5158E-02  0.0000E+00  0.0000E+00 -2.4587E-01  5.9661E-04
            -1.6787E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -921.837143543982        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1573
 NPARAMETR:  4.4823E-01  1.0000E-02  2.4997E-02  2.8939E-01  1.2112E+01  1.4756E+00  1.0000E-02  1.0000E-02  8.9072E-01  1.3958E-01
             8.1279E+00
 PARAMETER: -7.0245E-01 -4.5635E+00 -3.5890E+00 -1.1400E+00  2.5942E+00  4.8906E-01 -8.6154E+00 -1.1119E+01 -1.5722E-02 -1.8691E+00
             2.1953E+00
 GRADIENT:  -1.2053E+00  0.0000E+00 -9.0341E-02 -2.0216E+00 -2.2603E-03 -2.9903E-02  0.0000E+00  0.0000E+00 -1.0390E-01  2.7748E-06
            -2.5173E-01

0ITERATION NO.:   63    OBJECTIVE VALUE:  -921.838801046070        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1665
 NPARAMETR:  4.4829E-01  1.0000E-02  2.4990E-02  2.8959E-01  1.2075E+01  1.4752E+00  1.0000E-02  1.0000E-02  8.9142E-01  1.4809E-01
             8.1282E+00
 PARAMETER: -7.0232E-01 -4.5643E+00 -3.5893E+00 -1.1393E+00  2.5911E+00  4.8881E-01 -8.6154E+00 -1.1119E+01 -1.4940E-02 -1.8099E+00
             2.1953E+00
 GRADIENT:  -1.0110E+00  0.0000E+00 -1.4093E+00 -4.2281E-01 -2.0209E-03 -1.3373E-01  0.0000E+00  0.0000E+00 -7.8482E-03  3.2767E-06
            -3.0141E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1665
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2464E-03  1.2053E-06  7.5240E-05 -1.9509E-02 -3.9451E-06
 SE:             2.9151E-02  1.0770E-06  2.4677E-04  2.4978E-02  1.5432E-05
 N:                     100         100         100         100         100

 P VAL.:         9.1133E-01  2.6310E-01  7.6044E-01  4.3478E-01  7.9823E-01

 ETASHRINKSD(%)  2.3415E+00  9.9996E+01  9.9173E+01  1.6320E+01  9.9948E+01
 ETASHRINKVR(%)  4.6282E+00  1.0000E+02  9.9993E+01  2.9976E+01  1.0000E+02
 EBVSHRINKSD(%)  2.4083E+00  9.9995E+01  9.9234E+01  1.6832E+01  9.9948E+01
 EBVSHRINKVR(%)  4.7585E+00  1.0000E+02  9.9994E+01  3.0831E+01  1.0000E+02
 RELATIVEINF(%)  3.3307E+00  2.6504E-08  3.9111E-05  4.5067E-01  1.8241E-06
 EPSSHRINKSD(%)  1.2674E+01
 EPSSHRINKVR(%)  2.3741E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -921.83880104606953     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2.9002678413968397     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.03
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.35
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -921.839       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.48E-01  1.00E-02  2.50E-02  2.90E-01  1.21E+01  1.48E+00  1.00E-02  1.00E-02  8.91E-01  1.48E-01  8.13E+00
 


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
+        2.39E+03
 
 TH 2
+        0.00E+00  1.91E+03
 
 TH 3
+       -1.57E+04  0.00E+00  9.89E+05
 
 TH 4
+       -7.57E+01  0.00E+00 -9.72E+04  1.09E+04
 
 TH 5
+        3.68E-01  0.00E+00 -5.76E+00  4.48E-01  3.56E-03
 
 TH 6
+        1.72E+00  0.00E+00 -1.29E+02 -2.50E+01  4.08E-03  7.99E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.43E-01  0.00E+00  9.72E+02 -1.21E+02 -1.54E-02 -2.14E+00  0.00E+00  0.00E+00  1.17E+02
 
 TH10
+       -4.37E-03  0.00E+00  1.48E-02  1.81E-03 -4.24E-05  5.60E-04  0.00E+00  0.00E+00  4.53E-03  7.42E-03
 
 TH11
+       -2.49E+01  0.00E+00  3.46E+02 -2.17E+01 -4.64E-03  8.60E-01  0.00E+00  0.00E+00  4.42E+00 -4.76E-05  7.65E+00
 
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
 #CPUT: Total CPU Time in Seconds,       39.451
Stop Time:
Thu Sep 30 03:10:48 CDT 2021
